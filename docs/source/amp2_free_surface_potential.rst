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


Environment variable override
------------------------------

The number of sigma layers used for the H2 expansion can be overridden by
setting the environment variable ``SWD_NUM_H2_LAYERS`` before opening the
file::

   export SWD_NUM_H2_LAYERS=10

Acceptable values are integers in the range 2–100.  The default is 20.  This
variable is primarily intended for testing and sensitivity studies; the default
of 20 is recommended for production use.


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
split into ``nstep = 20`` incremental steps to improve accuracy.  At each step
the potential is first transferred to :math:`z=0` using the H-operator (the
Taylor expansion above), then lifted back to the next intermediate height using
a second application of the same operator — hence the name *H2*.

Non-linear products of the form :math:`\eta^m \cdot k^n\hat{\phi}` produce
aliasing in the spatial domain.  To suppress this, all intermediate
computations are performed on an **extended grid** of size :math:`8n_x`,
and the result is projected back to the original :math:`n_x`-point grid after
each step by spectral truncation.

The nonlinearity order is fixed at ``M_kin = 5``.  The FFT operations use an
internal pocketfft backend (BSD-3-Clause licence) compiled into the library.
All arithmetic is done in double precision regardless of the single-precision
storage in the SWD file.

**References:** The H2 operator is described in
West *et al.* (1987) *J. Geophys. Res.* **92** 11803–11824 and
Dommermuth & Yue (1987) *J. Fluid Mech.* **184** 267–288.
