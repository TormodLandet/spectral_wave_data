.. _amp2-free-surface-potential:

Free-surface-potential files (``amp=2``)
=========================================

Shape classes 1 and 2 can be written with ``amp=2``.  In this mode the stored
``c`` array contains the velocity potential **at the free surface** rather than
the linear-theory propagated value used by ``amp=1``.  This allows a nonlinear
wave generator (such as a High-Order Spectral Method, HOSM) to hand off its
fully nonlinear surface potential to SpectralWaveData, which then evaluates
kinematics at arbitrary depths without going through the source model again.

When you open a shape-1 or shape-2 file that has ``amp=2``, SpectralWaveData
automatically selects the H2-operator implementation.  No extra flag is needed:
the choice is made at construction time based on the file header.

.. note::

   ``amp=2`` is only supported for shape 1 and shape 2 (long-crested waves).
   Multi-directional shape classes 4 and 5 do not currently support ``amp=2``.

.. note::

   The physical grid size ``nx = 2n`` (where ``n`` is the highest spectral
   index stored in the file) must be a **power of two**.  Attempting to open a
   file where ``nx`` is not a power of two raises an error.


What happens when you call ``update_time``
------------------------------------------

At each call to ``update_time(t)`` the library determines which four SWD time
steps surround the requested time and ensures all four are present in an
internal cache.  For each step that is not yet cached, the library:

1. Reads ``h(0:n)`` (elevation) and ``c(0:n)`` (surface potential) from file.
2. Reconstructs the spatial fields ``η(x)`` and ``ψ(x)`` via an inverse FFT.
3. Runs the H2 operator to transfer ``ψ`` from the free surface down to a
   fixed reference level ``z_\text{ref}`` on ``nlayers = 20`` sigma-layers
   (see :ref:`amp2-h2-operator` below).
4. Stores the resulting spectral sigma-layer coefficients in the cache.

The four cached steps are then used together with second-order finite
differences to reconstruct time derivatives, and the result is handed to the
same sigma-coordinate evaluator used by :doc:`shape 7 <shape_7>`.

The reference level ``z_\text{ref}`` is determined **once** at construction
time by scanning all time steps in the file, finding the global minimum wave
trough, and setting

.. math::

   z_\text{ref} = \text{round}(1.2\, \eta_{\min})

This scan reads only the elevation array ``h`` for each step (one inverse FFT
per step) and is fast relative to H2 generation.  A fixed ``z_\text{ref}``
ensures the sigma layers span consistent physical depths across the entire
simulation.


Supported functions
--------------------

The ``amp=2`` implementation delegates all kinematic evaluation to the same
sigma-coordinate engine used by :doc:`shape_7`.  The set of implemented and
stub functions is therefore identical to the one described in
:ref:`the shape-7 function table <shape7-function-table>` — please refer to
that page for details.  In particular, ``grad_phi_2nd``, ``acc_euler``, and
``acc_particle`` currently return zeros, and ``convergence`` / ``strip`` raise
an error.


Performance guidance
---------------------

**Sequential stepping is cheap.**  When ``update_time`` is called in
monotonically increasing order with small steps relative to the SWD time
spacing ``dt``, the cache is advanced by at most one new H2 evaluation per
call.

**Seeking and rewinds are expensive.**  Jumping to a time that falls outside
the current four-step window requires up to four H2 evaluations to refill the
cache.  Each H2 evaluation involves :math:`O(M \cdot n_\text{step})` forward
and backward FFTs on an extended grid of size :math:`8n_x`.  For typical
parameters (:math:`n_x = 1024`, ``M_kin = 5``, ``nstep = 20``) this amounts
to hundreds of FFTs per H2 call.  If your application queries a single
isolated time instant (rather than advancing through a simulation),
the latency at each call will be noticeable.  If your application
advances time continuously — the common case in CFD — the cost is
amortised over many small solver time steps.


Thread safety
--------------

SpectralWaveData objects are **not safe to share between threads**.  Each
object maintains an internal cache of spectral data for the current time
window; concurrent reads and writes to that cache from multiple threads
produce undefined behaviour.  This is the same constraint that applies to all
shape classes (they all maintain a four-step window of file data).

The safe pattern — which is typical for CFD applications — is to create one
SpectralWaveData object per boundary patch or per MPI rank and never share an
object between threads.  Distinct objects are fully independent and can be used
concurrently without locking.


Environment variables
----------------------

The following environment variables may be set **before** opening the file to
control the H2 operator behaviour.  All variables are optional; the defaults
give sensible production settings.

``SWD_NUM_H2_STEPS`` (integer, default ``20``)
   Number of H2 iteration steps from the free surface down to ``z_\text{ref}``.
   The total number of sigma layers stored internally is ``nstep + 1``.
   Acceptable values are 2–100.  Larger values give higher accuracy at the cost
   of proportionally more FFT work per time step.  Primarily useful for
   sensitivity studies.

   Example::

      export SWD_NUM_H2_STEPS=10

``SWD_H2_ZREF`` (real, default: compute from ``eta_min``)
   When **positive** (or not set), the library scans all (or a window of)
   time steps at construction time, finds the global minimum wave trough
   :math:`\eta_{\min}`, and sets :math:`z_\text{ref} = 1.2\,\eta_{\min}`.
   When **negative**, the given value is used directly as :math:`z_\text{ref}`
   and the scan is skipped entirely — useful when opening many SWD files from
   the same wave field where the reference level is already known.

   For shape-2 files the library validates that :math:`d + z_\text{ref} > 0`;
   an explicit negative value that violates this constraint raises an error.

   Example (skip prepass, use :math:`z_\text{ref} = -5\,\text{m}`)::

      export SWD_H2_ZREF=-5.0

``SWD_WINDOW_TMIN`` / ``SWD_WINDOW_TMAX`` (real, user time)
   When ``SWD_H2_ZREF`` is not set (or is positive), the prepass scan is
   restricted to the time window ``[TMIN, TMAX]`` in user time.  This can
   substantially reduce construction time for very long simulation files when
   only a short segment will actually be queried.  The times are given in the
   same coordinate as the ``time`` argument to ``update_time``.

   .. note::

      ``SWD_WINDOW_TMIN``/``TMAX`` affect only the prepass scan.  They do not
      restrict which time steps you may later request via ``update_time``.
      Seeking outside the scanned window may result in a slightly suboptimal
      choice of ``z_\text{ref}`` if the deepest trough lies outside the window.

   Example (restrict to the first 60 s)::

      export SWD_WINDOW_TMIN=0.0
      export SWD_WINDOW_TMAX=60.0


.. _amp2-h2-operator:

Appendix: the H2 operator
---------------------------

The H2 operator propagates the free-surface velocity potential
:math:`\psi(x,t) = \phi(x,\eta(x,t),t)` downward from the instantaneous free
surface to the flat plane :math:`z = z_\text{ref}` in many small steps.  It is
a pseudospectral method based on a Taylor expansion of the velocity potential
around a flat surface:

.. math::

   \phi(x, z_0 + \Delta z, t)
   = \sum_{m=0}^{M-1} \frac{(\Delta z)^m}{m!}
     \left(\frac{\partial^m \phi}{\partial z^m}\right)_{z=z_0}

where the vertical derivatives are evaluated spectrally using the relation
:math:`\partial^n \phi / \partial z^n = \mathcal{F}^{-1}\{|k|^n \hat{\phi}\}`
(modified by the depth factor :math:`\tanh(|k|h)` for finite-depth shape 2
files).

The full downward transfer from :math:`z=\eta` to :math:`z=z_\text{ref}` is
split into ``nstep`` incremental steps (default 20) to improve accuracy.  At
each step the potential is first transferred to :math:`z=0` using the
H-operator (the Taylor expansion above), then lifted back to the next
intermediate height using a second application of the same operator — hence
the name *H2*.

After all ``nstep`` steps, the library produces ``nstep + 1`` sigma layers with
exact positions

.. math::

   \sigma_m = \frac{m-1}{n_\text{step}}, \quad m = 1,\ldots,n_\text{step}+1

(so :math:`\sigma_1 = 0` at :math:`z_\text{ref}`, and
:math:`\sigma_{n_\text{step}+1} = 1` at the free surface).
The free-surface layer (:math:`\sigma = 1`) stores the original input potential
:math:`\psi` directly; all other layers contain the H2-computed result for
that :math:`\sigma`-level.

Non-linear products of the form :math:`\eta^m \cdot k^n\hat{\phi}` produce
aliasing in the spatial domain.  To suppress this, all intermediate
computations are performed on an **extended grid** of size :math:`8n_x`,
and the result is projected back to the original :math:`n_x`-point grid after
each step by spectral truncation.

The nonlinearity order is fixed at ``M_kin = 5``.
The FFT operations use an internal pocketfft backend (BSD-3-Clause licence)
compiled into the library.  All arithmetic is done in double precision
regardless of the single-precision storage in the SWD file.

.. note::

   **Performance and grid size.**  PocketFFT handles any grid size, but
   :math:`n_x = 2n` runs fastest when :math:`n_x` is a product of small
   primes (ideally a power of two).  Grids from HOSM simulations are
   practically always powers of two, so this is rarely a concern.
   For the extended grid :math:`8n_x` the same rule applies.

**References:** The H2 operator is described in
West *et al.* (1987) *J. Geophys. Res.* **92** 11803–11824 and
Dommermuth & Yue (1987) *J. Fluid Mech.* **184** 267–288.
