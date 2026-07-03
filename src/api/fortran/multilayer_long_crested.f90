module multilayer_long_crested_def

! Shared sigma-coordinate kinematics and temporal-window application for
! multilayer long-crested wave shapes.
!
! The ml_state type holds all quantities needed for kinematic evaluation at
! the current time step.  It is filled by ml_apply_window after the owning
! shape class has assembled the four-step spectral window.
!
! Used by:
!   - spectral_wave_data_shape_7_impl_1  (file-backed sigma-coordinate waves)
!   - spectral_wave_data_shape_1_or_2_impl_7  (amp=2, H2-operator generated)

use, intrinsic :: iso_c_binding, only: c_float
use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal
use spectral_interpolation_def, only: spectral_interpolation

implicit none
private

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

public :: ml_state         ! State type (current-time kinematics data)
public :: ml_init          ! Allocate and initialise ml_state
public :: ml_close         ! Deallocate ml_state
public :: ml_apply_window  ! Fill h_cur/ht_cur/c_cur/ct_cur from 4-step window

! Kinematics  (match the spectral_wave_data abstract interface signatures)
public :: ml_elev, ml_elev_t, ml_grad_elev, ml_grad_elev_2nd
public :: ml_phi, ml_stream, ml_phi_t, ml_grad_phi, ml_pressure

!##############################################################################
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

!------------------------------------------------------------------------------
! State type
!------------------------------------------------------------------------------
type :: ml_state
    integer  :: n              ! Highest spectral index (n+1 components, j=0..n)
    integer  :: nsumx          ! Number of components used in summation (<=n)
    real(wp) :: dk             ! Wavenumber spacing
    real(wp) :: d              ! Water depth (d>0) or -1 for infinite depth
    real(wp) :: zref           ! z-pos of bottom sigma-layer (sigma=0); <= 0
    real(wp) :: tanhdkd_eff    ! tanh(dk*(d+zref)) for below-zref extrapolation
    integer  :: nlayers        ! Number of sigma-layers
    real(wp), allocatable :: sig(:)    ! (1:nlayers) sigma positions
                                       ! sig(1)=0, sig(nlayers)=1, strictly increasing
    real(wp) :: cbeta          ! cos(beta)
    real(wp) :: sbeta          ! sin(beta)
    real(wp) :: x0             ! SWD x-origin
    real(wp) :: grav           ! Acceleration of gravity [m/s^2]
    real(wp) :: rho            ! Water density [kg/m^3]
    logical  :: dc_bias        ! .true. = include k=0 contribution
    complex(wp), allocatable :: h_cur(:)    ! (0:n)  elevation spectral coeff
    complex(wp), allocatable :: ht_cur(:)   ! (0:n)  d/dt elevation
    complex(wp), allocatable :: c_cur(:,:)  ! (0:n, 1:nlayers) potential
    complex(wp), allocatable :: ct_cur(:,:) ! (0:n, 1:nlayers) d/dt potential
end type ml_state

contains

!==============================================================================
! Construction / destruction
!==============================================================================

subroutine ml_init(st, n, nsumx, dk, d, zref, tanhdkd_eff, nlayers, sig, &
                   cbeta, sbeta, x0, grav, rho, dc_bias)
! Allocate and initialise an ml_state from scalar parameters.
type(ml_state), intent(out) :: st
integer,  intent(in) :: n, nsumx, nlayers
real(wp), intent(in) :: dk, d, zref, tanhdkd_eff
real(wp), intent(in) :: sig(nlayers)
real(wp), intent(in) :: cbeta, sbeta, x0, grav, rho
logical,  intent(in) :: dc_bias
!
st%n           = n
st%nsumx       = nsumx
st%dk          = dk
st%d           = d
st%zref        = zref
st%tanhdkd_eff = tanhdkd_eff
st%nlayers     = nlayers
allocate(st%sig(nlayers))
st%sig         = sig
st%cbeta       = cbeta
st%sbeta       = sbeta
st%x0          = x0
st%grav        = grav
st%rho         = rho
st%dc_bias     = dc_bias
allocate(st%h_cur(0:n), st%ht_cur(0:n))
allocate(st%c_cur(0:n, nlayers), st%ct_cur(0:n, nlayers))
st%h_cur  = cmplx(0.0_wp, 0.0_wp, wp)
st%ht_cur = cmplx(0.0_wp, 0.0_wp, wp)
st%c_cur  = cmplx(0.0_wp, 0.0_wp, wp)
st%ct_cur = cmplx(0.0_wp, 0.0_wp, wp)
!
end subroutine ml_init

!==============================================================================

subroutine ml_close(st)
! Deallocate an ml_state.
type(ml_state), intent(inout) :: st
!
if (allocated(st%sig))    deallocate(st%sig)
if (allocated(st%h_cur))  deallocate(st%h_cur)
if (allocated(st%ht_cur)) deallocate(st%ht_cur)
if (allocated(st%c_cur))  deallocate(st%c_cur)
if (allocated(st%ct_cur)) deallocate(st%ct_cur)
st%n       = 0
st%nsumx   = 0
st%nlayers = 0
!
end subroutine ml_close

!==============================================================================
! Window application
!==============================================================================

subroutine ml_apply_window(st, h_win, c_win, tpol, i1, i2, i3, i4, delta, dt)
! Fill st%h_cur / st%ht_cur / st%c_cur / st%ct_cur from a four-step spectral
! window using 2nd-order finite-difference time derivatives + tpol interpolation.
!
!   i1 = column index for t_(i-1)
!   i2 = column index for t_i
!   i3 = column index for t_(i+1)
!   i4 = column index for t_(i+2)
!   delta = (tswd - t_i) / dt   in [0, 1)
!   dt    = SWD time step
type(ml_state),               intent(inout) :: st
complex(c_float),             intent(in)    :: h_win(0:, :)      ! (0:n, 4)
complex(c_float),             intent(in)    :: c_win(0:, :, :)   ! (0:n, 1:nlayers, 4)
type(spectral_interpolation), intent(in)    :: tpol
integer,  intent(in) :: i1, i2, i3, i4
real(wp), intent(in) :: delta, dt
!
integer  :: i, m
real(wp) :: dt2
complex(wp) :: fval, dfval
complex(wp), parameter :: czero = cmplx(0.0_wp, 0.0_wp, wp)
!
dt2 = 2.0_wp * dt
!
! Elevation h_cur / ht_cur
do concurrent (i = 0 : st%nsumx)
    call tpol%scheme(delta,                                                      &
        cmplx(h_win(i,i1), kind=wp), cmplx(h_win(i,i2), kind=wp),              &
        cmplx(h_win(i,i3), kind=wp), cmplx(h_win(i,i4), kind=wp),              &
        (-3.0_wp*cmplx(h_win(i,i1),kind=wp) + 4.0_wp*cmplx(h_win(i,i2),kind=wp) &
         - cmplx(h_win(i,i3),kind=wp)) / dt2,                                   &
        (cmplx(h_win(i,i3),kind=wp) - cmplx(h_win(i,i1),kind=wp)) / dt2,       &
        (cmplx(h_win(i,i4),kind=wp) - cmplx(h_win(i,i2),kind=wp)) / dt2,       &
        (cmplx(h_win(i,i2),kind=wp) - 4.0_wp*cmplx(h_win(i,i3),kind=wp)        &
         + 3.0_wp*cmplx(h_win(i,i4),kind=wp)) / dt2,                           &
        st%h_cur(i), st%ht_cur(i))
end do
!
! Potential c_cur / ct_cur (all sigma layers)
do m = 1, st%nlayers
    do concurrent (i = 0 : st%nsumx)
        call tpol%scheme(delta,                                                      &
            cmplx(c_win(i,m,i1),kind=wp), cmplx(c_win(i,m,i2),kind=wp),           &
            cmplx(c_win(i,m,i3),kind=wp), cmplx(c_win(i,m,i4),kind=wp),           &
            (-3.0_wp*cmplx(c_win(i,m,i1),kind=wp) + 4.0_wp*cmplx(c_win(i,m,i2),kind=wp) &
             - cmplx(c_win(i,m,i3),kind=wp)) / dt2,                               &
            (cmplx(c_win(i,m,i3),kind=wp) - cmplx(c_win(i,m,i1),kind=wp)) / dt2, &
            (cmplx(c_win(i,m,i4),kind=wp) - cmplx(c_win(i,m,i2),kind=wp)) / dt2, &
            (cmplx(c_win(i,m,i2),kind=wp) - 4.0_wp*cmplx(c_win(i,m,i3),kind=wp)  &
             + 3.0_wp*cmplx(c_win(i,m,i4),kind=wp)) / dt2,                        &
            st%c_cur(i,m), st%ct_cur(i,m))
    end do
end do
!
! DC-bias suppression
if (.not. st%dc_bias) then
    st%h_cur(0)    = czero
    st%ht_cur(0)   = czero
    st%c_cur(0,:)  = czero
    st%ct_cur(0,:) = czero
end if
!
end subroutine ml_apply_window

!==============================================================================
! Private helpers
!==============================================================================

subroutine calc_elev_and_slope(st, xswd, eta, zeta_x, zeta_t)
! Compute wave elevation, its x-derivative and its Euler d/dt at xswd.
type(ml_state), intent(in)  :: st
real(wp),       intent(in)  :: xswd
real(wp),       intent(out) :: eta, zeta_x, zeta_t
!
integer  :: j
real(wp) :: kval
complex(wp) :: kappa1, Xfun
!
kappa1 = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
Xfun   = 1.0_wp
eta    = st%h_cur(0) % re
zeta_x = 0.0_wp
zeta_t = st%ht_cur(0) % re
kval   = 0.0_wp
do j = 1, st%nsumx
    Xfun   = kappa1 * Xfun
    kval   = kval + st%dk
    eta    = eta    + real(st%h_cur(j)  * Xfun)
    zeta_x = zeta_x + kval * aimag(st%h_cur(j) * Xfun)
    zeta_t = zeta_t + real(st%ht_cur(j) * Xfun)
end do
!
end subroutine calc_elev_and_slope

!==============================================================================

subroutine depth_Zfun(st, kval, z_prime, Zfun, Zfun_z, Rfun)
! Vertical basis function (and its z-derivative) for below-zref shape-2
! extrapolation.  z_prime = z - zref (<= 0).
! Rfun carries the Rfun recursion state; caller must initialise to 0.0
! before the j=1..nsumx loop and pass it unchanged between calls.
type(ml_state), intent(in)    :: st
real(wp),       intent(in)    :: kval, z_prime
real(wp),       intent(inout) :: Rfun
real(wp),       intent(out)   :: Zfun, Zfun_z
!
real(wp) :: kappa2, kappa3, Sfun, Tfun, Ufun, Vfun
real(wp), parameter :: Rfun_eps = 100.0_wp * epsilon(1.0_wp)
!
kappa2 = exp( kval * z_prime)
kappa3 = 1.0_wp / kappa2
Sfun   = kappa2
Tfun   = kappa3
if (st%d <= 0.0_wp) then
    ! Infinite depth
    Zfun   = Sfun
    Zfun_z = kval * Sfun
else
    Rfun = (Rfun + st%tanhdkd_eff) / (1.0_wp + st%tanhdkd_eff * Rfun)
    if (1.0_wp - Rfun < Rfun_eps) then
        Zfun   = Sfun
        Zfun_z = kval * Sfun
    else
        Ufun   = (1.0_wp + Rfun) * 0.5_wp
        Vfun   = 1.0_wp - Ufun
        Zfun   = Ufun * Sfun + Vfun * Tfun
        Zfun_z = kval * (Ufun * Sfun - Vfun * Tfun)
    end if
end if
!
end subroutine depth_Zfun

!==============================================================================

subroutine interp_sigma(st, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, c_arr)
! Piecewise-linear sigma-layer interpolation of phi and u_x at sigma_eval.
! c_arr is either st%c_cur or st%ct_cur: assumed shape (0:, 1:).
type(ml_state), intent(in)  :: st
real(wp),       intent(in)  :: xswd, sigma_eval
real(wp),       intent(out) :: phi_val, u_val, dphi_dsigma
complex(wp),    intent(in)  :: c_arr(0:, 1:)
!
integer  :: j, m, m1, m2
real(wp) :: kval, s, dsigma
real(wp) :: phi_m(st%nlayers), u_m(st%nlayers)
complex(wp) :: kappa1, Xfun, cx
!
kappa1 = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
phi_m  = 0.0_wp
u_m    = 0.0_wp
kval   = 0.0_wp
Xfun   = 1.0_wp
! j=0 DC contribution (dc_bias zeroing already applied to c_arr)
do m = 1, st%nlayers
    phi_m(m) = phi_m(m) + c_arr(0,m) % re
end do
do j = 1, st%nsumx
    kval = kval + st%dk
    Xfun = kappa1 * Xfun
    do m = 1, st%nlayers
        cx       = c_arr(j, m) * Xfun
        phi_m(m) = phi_m(m) + cx % re
        u_m(m)   = u_m(m)   + kval * cx % im
    end do
end do
! Find bracket [m1, m2] containing sigma_eval
m1 = st%nlayers - 1
do m = 1, st%nlayers - 1
    if (sigma_eval <= st%sig(m+1)) then
        m1 = m
        exit
    end if
end do
m2 = m1 + 1
dsigma      = st%sig(m2) - st%sig(m1)
s           = (sigma_eval - st%sig(m1)) / dsigma
phi_val     = phi_m(m1) + s * (phi_m(m2) - phi_m(m1))
u_val       = u_m(m1)   + s * (u_m(m2)   - u_m(m1))
dphi_dsigma = (phi_m(m2) - phi_m(m1)) / dsigma
!
end subroutine interp_sigma

!==============================================================================
! Kinematic functions
!==============================================================================

function ml_phi(st, x, y, z) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer  :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z
real(wp) :: sigma_eval, phi_val, u_val, dphi_dsigma
complex(wp) :: kappa1, Xfun, cx
!
xswd = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
!
if (real(z,wp) >= st%zref) then
    call calc_elev_and_slope(st, xswd, eta, zeta_x, zeta_t)
    sigma_eval = min(max((real(z,wp) - st%zref) / (eta - st%zref), 0.0_wp), 1.0_wp)
    call interp_sigma(st, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, st%c_cur)
    res = phi_val
else
    ! Below zref: shape-2 style extrapolation from bottom layer
    z_prime = real(z,wp) - st%zref
    kappa1 = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
    Xfun   = 1.0_wp
    res    = st%c_cur(0,1) % re
    Rfun   = 0.0_wp
    kval   = 0.0_wp
    do j = 1, st%nsumx
        kval = kval + st%dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(st, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx  = st%c_cur(j,1) * Xfun
        res = res + cx % re * Zfun
    end do
end if
!
end function ml_phi

!==============================================================================

function ml_stream(st, x, y, z) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
res = 0.0_knd
!
end function ml_stream

!==============================================================================

function ml_phi_t(st, x, y, z) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer  :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z
real(wp) :: sigma_eval, phi_t_sig, phi_val, u_val, dphi_dsigma, H
complex(wp) :: kappa1, Xfun
!
xswd = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
!
if (real(z,wp) >= st%zref) then
    call calc_elev_and_slope(st, xswd, eta, zeta_x, zeta_t)
    H = eta - st%zref
    sigma_eval = min(max((real(z,wp) - st%zref) / H, 0.0_wp), 1.0_wp)
    ! Phi_t at fixed sigma from ct_cur
    call interp_sigma(st, xswd, sigma_eval, phi_t_sig, u_val, dphi_dsigma, st%ct_cur)
    ! phi_sigma from c_cur (for chain-rule correction at fixed physical z)
    call interp_sigma(st, xswd, sigma_eval, phi_val,   u_val, dphi_dsigma, st%c_cur)
    ! Euler time derivative: Phi_t|sigma - sigma * zeta_t / H * dphi_dsigma
    res = phi_t_sig - sigma_eval * zeta_t / H * dphi_dsigma
else
    z_prime = real(z,wp) - st%zref
    kappa1  = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
    Xfun    = 1.0_wp
    res     = st%ct_cur(0,1) % re
    Rfun    = 0.0_wp
    kval    = 0.0_wp
    do j = 1, st%nsumx
        kval = kval + st%dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(st, kval, z_prime, Zfun, Zfun_z, Rfun)
        res  = res + real(st%ct_cur(j,1) * Xfun) * Zfun
    end do
end if
!
end function ml_phi_t

!==============================================================================

function ml_grad_phi(st, x, y, z) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y, z
real(knd)             :: res(3)
!
integer  :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z
real(wp) :: sigma_eval, phi_val, u_val, dphi_dsigma, H
complex(wp) :: kappa1, Xfun, cx
!
xswd = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
!
if (real(z,wp) >= st%zref) then
    call calc_elev_and_slope(st, xswd, eta, zeta_x, zeta_t)
    H = eta - st%zref
    sigma_eval = min(max((real(z,wp) - st%zref) / H, 0.0_wp), 1.0_wp)
    call interp_sigma(st, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, st%c_cur)
    res(1) = (u_val - sigma_eval * zeta_x / H * dphi_dsigma) * st%cbeta
    res(2) = (u_val - sigma_eval * zeta_x / H * dphi_dsigma) * st%sbeta
    res(3) = dphi_dsigma / H
else
    z_prime = real(z,wp) - st%zref
    kappa1  = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
    Xfun    = 1.0_wp
    Rfun    = 0.0_wp
    kval    = 0.0_wp
    res     = 0.0_knd
    do j = 1, st%nsumx
        kval = kval + st%dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(st, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx      = st%c_cur(j,1) * Xfun
        res(1)  = res(1) + kval * cx % im * Zfun
        res(3)  = res(3) + cx % re * Zfun_z
    end do
    res(2) = res(1) * st%sbeta
    res(1) = res(1) * st%cbeta
end if
!
end function ml_grad_phi

!==============================================================================

function ml_elev(st, x, y) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y
real(knd)             :: res
!
integer  :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd   = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
kappa1 = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
Xfun   = 1.0_wp
res    = st%h_cur(0) % re
do j = 1, st%nsumx
    Xfun = kappa1 * Xfun
    res  = res + real(st%h_cur(j) * Xfun)
end do
!
end function ml_elev

!==============================================================================

function ml_elev_t(st, x, y) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y
real(knd)             :: res
!
integer  :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd   = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
kappa1 = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
Xfun   = 1.0_wp
res    = st%ht_cur(0) % re
do j = 1, st%nsumx
    Xfun = kappa1 * Xfun
    res  = res + real(st%ht_cur(j) * Xfun)
end do
!
end function ml_elev_t

!==============================================================================

function ml_grad_elev(st, x, y) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y
real(knd)             :: res(3)
!
integer  :: j
real(wp) :: xswd, elev_x_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd       = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
kappa1     = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
Xfun       = 1.0_wp
elev_x_swd = 0.0_wp
kval       = 0.0_wp
do j = 1, st%nsumx
    Xfun       = kappa1 * Xfun
    kval       = kval + st%dk
    elev_x_swd = elev_x_swd + kval * aimag(st%h_cur(j) * Xfun)
end do
res(1) = elev_x_swd * st%cbeta
res(2) = elev_x_swd * st%sbeta
res(3) = 0.0_knd
!
end function ml_grad_elev

!==============================================================================

function ml_grad_elev_2nd(st, x, y) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y
real(knd)             :: res(3)
!
integer  :: j
real(wp) :: xswd, elev_xx_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd        = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
kappa1      = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
Xfun        = 1.0_wp
elev_xx_swd = 0.0_wp
kval        = 0.0_wp
do j = 1, st%nsumx
    Xfun        = kappa1 * Xfun
    kval        = kval + st%dk
    elev_xx_swd = elev_xx_swd - kval * kval * real(st%h_cur(j) * Xfun)
end do
res(1) = elev_xx_swd * st%cbeta * st%cbeta
res(2) = elev_xx_swd * st%sbeta * st%cbeta
res(3) = elev_xx_swd * st%sbeta * st%sbeta
!
end function ml_grad_elev_2nd

!==============================================================================

function ml_pressure(st, x, y, z) result(res)
type(ml_state), intent(in) :: st
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer  :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z
real(wp) :: sigma_eval, phi_val, u_val, dphi_dsigma, H
real(wp) :: phi_xswd, phi_z, phi_t_val, phi_sig
complex(wp) :: kappa1, Xfun, cx
!
xswd = st%x0 + real(x,wp) * st%cbeta + real(y,wp) * st%sbeta
!
if (real(z,wp) >= st%zref) then
    call calc_elev_and_slope(st, xswd, eta, zeta_x, zeta_t)
    H = eta - st%zref
    sigma_eval = min(max((real(z,wp) - st%zref) / H, 0.0_wp), 1.0_wp)
    ! Spatial gradients from c_cur
    call interp_sigma(st, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, st%c_cur)
    phi_xswd = u_val - sigma_eval * zeta_x / H * dphi_dsigma
    phi_z    = dphi_dsigma / H
    phi_sig  = dphi_dsigma  ! save dphi/dsigma for chain-rule correction
    ! Time derivative at fixed sigma from ct_cur
    call interp_sigma(st, xswd, sigma_eval, phi_t_val, u_val, dphi_dsigma, st%ct_cur)
    ! Euler chain-rule: phi_t|z = Phi_t|sigma - sigma * zeta_t / H * phi_sigma
    phi_t_val = phi_t_val - sigma_eval * zeta_t / H * phi_sig
else
    z_prime  = real(z,wp) - st%zref
    kappa1   = exp(cmplx(0.0_wp, -st%dk * xswd, kind=wp))
    Xfun     = 1.0_wp
    Rfun     = 0.0_wp
    kval     = 0.0_wp
    phi_xswd = 0.0_wp
    phi_z    = 0.0_wp
    phi_t_val = st%ct_cur(0,1) % re
    do j = 1, st%nsumx
        kval = kval + st%dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(st, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx        = st%c_cur(j,1) * Xfun
        phi_xswd  = phi_xswd  + kval * cx % im * Zfun
        phi_z     = phi_z     + cx % re * Zfun_z
        phi_t_val = phi_t_val + real(st%ct_cur(j,1) * Xfun) * Zfun
    end do
end if
!
res = (-phi_t_val - 0.5_wp*(phi_xswd**2 + phi_z**2) - real(z,wp)*st%grav) * st%rho
!
end function ml_pressure

!==============================================================================

end module multilayer_long_crested_def
