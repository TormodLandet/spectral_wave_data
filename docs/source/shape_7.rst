Shape class 7
-------------

.. note::

  Shape class 7 is **experimental**. The file format and API may change in
  future versions.

This shape class describes long-crested waves in a **wave-following sigma-coordinate**
representation.  Instead of evaluating the velocity potential with a single vertical basis function
(as in shape 2), the potential is stored on :math:`N_\sigma` horizontal layers that deform with the
free surface.

The sigma coordinate is defined as

.. math::

   \sigma(x, z, t) = \frac{z - z_\text{ref}}{\zeta(x,t) - z_\text{ref}},
   \qquad \sigma \in [0, 1]

where :math:`z_\text{ref}` is a fixed reference level chosen safely below all wave troughs
(typically :math:`z_\text{ref} \approx -3H_s`).  The bottom layer at :math:`\sigma = 0` lies on the
flat plane :math:`z = z_\text{ref}`; the surface layer at :math:`\sigma = 1` follows the free
surface :math:`z = \zeta(x,t)`.

On each layer :math:`\sigma_m` the velocity potential is stored spectrally:

.. math::

   \phi_m(x, t) = \sum_{j=0}^n \mathcal{Re}\bigl\{c_{j,m}(t)\,X_j(x)\bigr\},
   \qquad X_j(x) = e^{-ik_jx},\quad k_j = j\Delta k

where :math:`c_{j,m}(t)` are the spectral amplitudes on layer :math:`m`. The surface elevation is
stored identically to shape 2:

.. math::

   \zeta(x, t) = \sum_{j=0}^n \mathcal{Re}\bigl\{h_j(t)\,X_j(x)\bigr\}

The layered representation avoids the issues with the exponential behaviour in z for large k where
:math:`k_j z_\text{crest}` is large and can cause issues when working in limited precision, both
when computing :math:`e^{k_j z}\rightarrow\inf` and when computing :math:`e^{-k_j z}\rightarrow 0`.


Kinematics
^^^^^^^^^^

**Layered region** :math:`z_\text{ref} \le z \le \zeta(x,t)`:

The sigma coordinate at the query point is

.. math::

   \sigma = \frac{z - z_\text{ref}}{\zeta(x,t) - z_\text{ref}}

and :math:`H = \zeta - z_\text{ref}` is the local column height. The potential is interpolated
linearly (impl=1) between the two adjacent :math:`\sigma`-layers bracketing :math:`\sigma`:

.. math::

   \phi(x, z, t) \approx \phi_m + s\,(\phi_{m+1} - \phi_m),
   \qquad s = \frac{\sigma - \sigma_m}{\sigma_{m+1} - \sigma_m}

The particle velocities are then obtained by differentiating with respect to the physical
coordinates, accounting for the curvilinear nature of the :math:`\sigma` coordinate:

.. math::

   u(x, z, t) = u_\sigma - \frac{\sigma\,\zeta_x}{H}\,\phi_\sigma

.. math::

   w(x, z, t) = \frac{\phi_\sigma}{H}

where

.. math::

   u_\sigma = \sum_{j=1}^n k_j\,\mathcal{Im}\bigl\{c_{j,m^*}(t)\,X_j(x)\bigr\}
   \quad\text{(interpolated between layers)}

.. math::

   \phi_\sigma = \frac{\phi_{m+1} - \phi_m}{\sigma_{m+1} - \sigma_m}

.. math::

   \zeta_x = \sum_{j=1}^n k_j\,\mathcal{Im}\bigl\{h_j(t)\,X_j(x)\bigr\}

**Below** :math:`z_\text{ref}`:

A shape-2 style extrapolation is applied using the bottom-layer coefficients :math:`c_{j,1}(t)` (the
layer at :math:`\sigma = 0`).  Setting :math:`z' = z - z_\text{ref}` (negative) the velocity
potential is

.. math::

   \phi(x, z, t) = \sum_{j=0}^n \mathcal{Re}\bigl\{c_{j,1}(t)\,X_j(x)\bigr\}\,Z_j(z')

with the same finite/infinite-depth :math:`Z_j` as shape 2, evaluated at :math:`z'` relative to
:math:`z_\text{ref}`:

.. math::

   Z_j(z') = \frac{\cosh k_j(z' + d_\text{eff})}{\cosh k_j d_\text{eff}},
   \qquad d_\text{eff} = d + z_\text{ref}

or :math:`Z_j(z') = e^{k_j z'}` for infinite depth.

**Above** :math:`\zeta(x,t)`:

The surface value (:math:`\sigma = 1`) is returned.


Pressure
^^^^^^^^

The fully nonlinear Bernoulli pressure is

.. math::

   p = -\rho\,\frac{\partial\phi}{\partial t}
       - \frac{1}{2}\rho\,\bar{\nabla}\phi\cdot\bar{\nabla}\phi
       - \rho g\,\bar{z}

where :math:`\partial\phi/\partial t` is the Euler time derivative of the interpolated potential
and :math:`\bar{\nabla}\phi` is the particle velocity vector.


Parameters
^^^^^^^^^^

Shape class 7 requires the following additional parameters beyond the common
SWD header:

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - :math:`n`
     - int
     - Number of spectral components (highest wavenumber index)
   * - :math:`\Delta k`
     - float
     - Spacing of wave numbers
   * - :math:`d`
     - float
     - Water depth (:math:`d > 0`) or :math:`-1` for infinite depth
   * - :math:`z_\text{ref}`
     - float
     - Constant :math:`z`-position of the lowest :math:`\sigma`-layer
       (:math:`\sigma=0`).  Should be below all wave troughs.
   * - :math:`N_\sigma` (``nsig``)
     - int
     - Number of :math:`\sigma`-layers (:math:`N_\sigma \ge 2`)
   * - :math:`\sigma_1, \dots, \sigma_{N_\sigma}`
     - float
     - Tabulated layer positions.  Must satisfy :math:`\sigma_1 = 0` and
       :math:`\sigma_{N_\sigma} = 1`.


Implementation notes
^^^^^^^^^^^^^^^^^^^^

**Sum-first strategy (impl=1).**  Because the piecewise-linear
interpolation weights are purely geometric (independent of :math:`k_j`), the
summation over wave numbers can be performed *before* the layer
interpolation.  This reduces the key inner loop from the
:math:`(N_x, N_\sigma, N_k)` tensor contraction to two :math:`(N_x, N_\sigma)` matrix
multiplies:

.. math::

   \phi_m(x_i) = \sum_{j=1}^n \mathcal{Re}\{c_{j,m}\,X_{j,i}\}
   = \mathcal{Re}\bigl\{\mathbf{X}_i \cdot \mathbf{c}_m\bigr\}

.. math::

   u_m(x_i) = \sum_{j=1}^n k_j\,\mathcal{Im}\{c_{j,m}\,X_{j,i}\}
   = \mathcal{Im}\bigl\{(\mathbf{k} \odot \mathbf{X}_i) \cdot \mathbf{c}_m\bigr\}

followed by scalar linear interpolation in :math:`\sigma`.

**Temporal interpolation.**  No time derivatives :math:`\dot{h}_j` or
:math:`\dot{c}_{j,m}` are stored in the SWD file.  The Fortran API
reconstructs them numerically using 2nd-order finite differences from the
four-step temporal window, exactly as for the other shape classes but without
the stored derivative data.

**Fortran acceleration.**  The inner spectral summation loop iterates over
:math:`j = 1, \dots, n_{sum}` and all :math:`N_\sigma` layers simultaneously,
exploiting the recursive evaluation :math:`X_j = \kappa_1 \cdot X_{j-1}`
with :math:`\kappa_1 = e^{-i \Delta k\, x_{swd}}`.

See :doc:`swd_format` for the binary file layout.
