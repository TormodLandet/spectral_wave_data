module hosm_h2_operator_def

! Port of the Python HosmH2Operator class (1-D, long-crested only).
!
! Computes the velocity potential at multiple sigma levels by iteratively
! applying the H2-operator (Taylor expansion in surface elevation) starting
! from the HOSM free-surface potential psi.
!
! All arithmetic is in float64/complex128 regardless of kind_swd_internal.
!
! The implementation is based on a Python reference implementation from the
! DNV-internal dnvgl-waves, "Waves", library which is not open source.
! The H2-operator implementation in Waves was *heavily* inspired by (copied
! from) the H2-operator Python code kindly provided by Christopher Lawrence
! with some minor extensions by Tormod Landet (both also from DNV).

use, intrinsic :: iso_c_binding, only: c_double
use swd_fft_def, only: swd_fft_plan, fft_init, fft_destroy, &
                       fft_r2c, fft_c2r, fft_dealias, &
                       fft_resample_up, fft_resample_down, &
                       fft_swd_to_real, fft_real_to_swd

implicit none
private

integer, parameter :: dp = c_double

!##############################################################################
!                    P U B L I C    Q U A N T I T I E S
!##############################################################################

public :: h2op_state
public :: h2op_init, h2op_close
public :: h2op_calc_potential
public :: h2op_convert_to_swd_layers

!##############################################################################

type :: h2op_state
    !-- grid parameters (set in h2op_init)
    integer  :: M_kin           ! nonlinearity order
    integer  :: nstep           ! number of H2 steps from surface to zref
    integer  :: nx              ! original grid size (power of 2)
    integer  :: nxE             ! extended grid size (= 8*nx, power of 2)
    integer  :: n_swd           ! highest SWD spectral index (= nx/2)
    real(dp) :: dk              ! wavenumber spacing of original grid
    real(dp) :: h_depth         ! water depth (>0) or -1 for infinite
    real(dp) :: zref            ! z-position of bottom sigma layer (<= 0)
    !-- precomputed spectral arrays (0:nxE/2)
    real(dp), allocatable :: kE(:)      ! wavenumbers: kE(j) = j*dk
    real(dp), allocatable :: thkhE(:)   ! tanh(kE * dist) or 1.0 for infinite depth
    !-- per-object FFT plans (thread-safe: never shared across objects)
    type(swd_fft_plan) :: plan_nx   ! for original grid
    type(swd_fft_plan) :: plan_nxE  ! for extended grid
    !-- result of last h2op_calc_potential (0:nxE-1, 0:nstep-1)
    real(dp), allocatable :: phiE(:,:)
    logical :: has_result = .false.
end type h2op_state

contains

!==============================================================================
! h2op_init
!==============================================================================

subroutine h2op_init(op, M_kin, nstep, nx, dk, h_depth, zref, err_msg)
! Initialise the H2 operator.  nx must be a power of two.
type(h2op_state), intent(out) :: op
integer,          intent(in)  :: M_kin, nstep, nx
real(dp),         intent(in)  :: dk, h_depth, zref
character(len=*), intent(out) :: err_msg

integer  :: j, nxE
real(dp) :: dist

err_msg = ''

! Validate nx is a power of two
if (nx <= 0 .or. iand(nx, nx-1) /= 0) then
    write(err_msg,'(a,i0)') 'h2op_init: nx must be a power of two, got nx=', nx
    return
end if

op%M_kin   = M_kin
op%nstep   = nstep
op%nx      = nx
op%nxE     = 8 * nx        ! oversampling factor hardcoded to 8
op%n_swd   = nx / 2
op%dk      = dk
op%h_depth = h_depth
op%zref    = zref
nxE = op%nxE

! Allocate spectral arrays for extended grid
allocate(op%kE(0:nxE/2))
allocate(op%thkhE(0:nxE/2))
do j = 0, nxE/2
    op%kE(j) = real(j, dp) * dk
end do

! thkhE: the even-order H-operator factor tanh(k * dist)
! dist = h_depth + zref  (distance from sea floor to the zref level)
if (h_depth < 0.0_dp) then
    ! infinite depth: thkhE = 1 for k>0, 0 for k=0
    op%thkhE(0) = 0.0_dp
    op%thkhE(1:nxE/2) = 1.0_dp
else
    dist = h_depth + zref   ! both h_depth>0 and zref<0, so dist<h_depth
    do j = 0, nxE/2
        op%thkhE(j) = tanh(op%kE(j) * dist)
    end do
end if

! Create FFT plans
call fft_init(op%plan_nx, nx, err_msg)
if (err_msg /= '') return
call fft_init(op%plan_nxE, nxE, err_msg)
if (err_msg /= '') return

allocate(op%phiE(0:nxE-1, 0:nstep-1))
op%phiE       = 0.0_dp
op%has_result = .false.

end subroutine h2op_init

!==============================================================================
! h2op_close
!==============================================================================

subroutine h2op_close(op)
type(h2op_state), intent(inout) :: op
call fft_destroy(op%plan_nx)
call fft_destroy(op%plan_nxE)
if (allocated(op%kE))   deallocate(op%kE)
if (allocated(op%thkhE)) deallocate(op%thkhE)
if (allocated(op%phiE)) deallocate(op%phiE)
op%has_result = .false.
end subroutine h2op_close

!==============================================================================
! h2op_calc_potential
!==============================================================================

subroutine h2op_calc_potential(op, eta_nx, psi_nx, err_msg)
! Compute the velocity potential at all nstep sigma layers by applying the
! H2-operator iteratively from the free surface (z=eta) down to z=zref.
!
! Inputs (on the original nx grid):
!   eta_nx(0:nx-1) : surface elevation
!   psi_nx(0:nx-1) : free-surface velocity potential psi = phi(x, z=eta, t)
!
! Output stored in op%phiE(0:nxE-1, 0:nstep-1):
!   Column ilayer holds the total potential (summed over M_kin orders) at the
!   sigma-layer position sigma = 1 - (ilayer+1)/nstep on the extended grid.
!   ilayer=0 is just below the free surface; ilayer=nstep-1 is at z=zref.

type(h2op_state), intent(inout) :: op
real(dp),         intent(in)    :: eta_nx(0:), psi_nx(0:)
character(len=*), intent(out)   :: err_msg

integer  :: ilayer, um, un, nxE, nhalf_E
real(dp) :: sigma_val, this_sigma, prev_sigma

! Working arrays on extended grid
real(dp),    allocatable :: etaE(:)         ! (0:nxE-1) surface elevation
real(dp),    allocatable :: psiE(:)         ! (0:nxE-1) surface potential
real(dp),    allocatable :: phi0E(:,:)      ! (0:nxE-1, 0:M_kin-1) H-op result
real(dp),    allocatable :: phi2E(:,:)      ! (0:nxE-1, 0:M_kin-1) H2-op result
real(dp),    allocatable :: layer_h(:)      ! (0:nxE-1) previous layer height
real(dp),    allocatable :: layer_psi(:)    ! (0:nxE-1) previous layer potential
real(dp),    allocatable :: this_h(:)       ! (0:nxE-1) current layer height
real(dp),    allocatable :: rho_work(:)     ! (0:nxE-1) scratch
real(dp),    allocatable :: etaE_pow(:)     ! (0:nxE-1) etaE^(un+1)
real(dp),    allocatable :: height_total(:) ! (0:nxE-1) etaE - zref

! Spectral workspace (0:nxE/2)
complex(dp), allocatable :: spec_phi(:)
complex(dp), allocatable :: spec_rho(:)
real(dp),    allocatable :: k_power(:)      ! (0:nxE/2) k^(un+1)
integer  :: facn

err_msg = ''
nxE     = op%nxE
nhalf_E = nxE / 2

allocate(etaE(0:nxE-1), psiE(0:nxE-1))
allocate(phi0E(0:nxE-1, 0:op%M_kin-1))
allocate(phi2E(0:nxE-1, 0:op%M_kin-1))
allocate(layer_h(0:nxE-1), layer_psi(0:nxE-1))
allocate(this_h(0:nxE-1), rho_work(0:nxE-1))
allocate(etaE_pow(0:nxE-1), height_total(0:nxE-1))
allocate(spec_phi(0:nhalf_E), spec_rho(0:nhalf_E))
allocate(k_power(0:nhalf_E))

! Upsample eta and psi from nx to nxE
call fft_resample_up(op%plan_nx, op%plan_nxE, eta_nx, op%nx, nxE, etaE, err_msg)
if (err_msg /= '') return
call fft_resample_up(op%plan_nx, op%plan_nxE, psi_nx, op%nx, nxE, psiE, err_msg)
if (err_msg /= '') return

! Total height above zref = etaE - zref
height_total = etaE - op%zref

! Initialise prev layer as free surface
layer_h(:)   = height_total(:)   ! previous height = full wave column
layer_psi(:) = psiE(:)           ! previous potential = surface potential

! Iterate from free surface down to zref
do ilayer = 0, op%nstep - 1
    this_sigma = 1.0_dp - real(ilayer + 1, dp) / real(op%nstep, dp)
    this_h(:)  = height_total(:) * this_sigma   ! fraction of total height

    if (ilayer < op%nstep - 1) then
        ! H2-operator: transfer from layer_h to this_h
        call perform_h2_operator(op, layer_psi, layer_h, this_h, &
                                 phi0E, phi2E, rho_work, etaE_pow, &
                                 spec_phi, spec_rho, k_power, err_msg)
        if (err_msg /= '') return
    else
        ! Final step: H-operator to transfer to the flat reference level
        phi0E = 0.0_dp
        phi0E(:, 0) = layer_psi(:)
        call perform_h_operator(op, phi0E, layer_h, &
                                rho_work, etaE_pow, &
                                spec_phi, spec_rho, k_power, err_msg)
        if (err_msg /= '') return
        phi2E = phi0E  ! H-op result is the final phi2E
    end if

    ! Dealias: zero modes above nx/2 for each M_kin layer
    do um = 0, op%M_kin - 1
        call fft_dealias(op%plan_nxE, phi2E(:, um), nxE, op%nx, err_msg)
        if (err_msg /= '') return
    end do

    ! Store sum over M_kin orders as the potential at this sigma level
    op%phiE(:, ilayer) = 0.0_dp
    do um = 0, op%M_kin - 1
        op%phiE(:, ilayer) = op%phiE(:, ilayer) + phi2E(:, um)
    end do

    ! Update for next iteration: layer position and potential
    layer_h(:)   = this_h(:)
    layer_psi(:) = op%phiE(:, ilayer)
end do

op%has_result = .true.

deallocate(etaE, psiE, phi0E, phi2E, layer_h, layer_psi, this_h)
deallocate(rho_work, etaE_pow, height_total)
deallocate(spec_phi, spec_rho, k_power)

end subroutine h2op_calc_potential

!==============================================================================
! h2op_convert_to_swd_layers
!==============================================================================

subroutine h2op_convert_to_swd_layers(op, eta_nx, psi_nx, h_swd, c_swd, err_msg)
! Convert H2 results to SWD shape-7 spectral coefficients.
!
! Produces nlayers = nstep+1 layers with exact sigma positions
!   sigma(m) = (m-1) / nstep,  m = 1..nstep+1
! so sigma runs from 0.0 (zref) in steps of 1/nstep up to 1.0 (free surface).
!
! Layer mapping:
!   m = 1..nstep : op%phiE(:, nstep-m) resampled from nxE to nx grid
!   m = nstep+1  : psi_nx (free-surface potential from SWD file, no resampling needed)
!
! Outputs:
!   h_swd(0:n_swd)              : elevation spectral coefficients (from eta_nx)
!   c_swd(0:n_swd, 1:nstep+1)   : potential spectral coefficients per sigma layer
type(h2op_state), intent(in)  :: op
real(dp),         intent(in)  :: eta_nx(0:), psi_nx(0:)
complex(dp),      intent(out) :: h_swd(0:)     ! (0:n_swd)
complex(dp),      intent(out) :: c_swd(0:, 1:) ! (0:n_swd, 1:nstep+1)
character(len=*), intent(out) :: err_msg

integer  :: m, nxE, n_swd
real(dp), allocatable :: phi_nx(:), phi_nxE_slice(:)

err_msg = ''
if (.not. op%has_result) then
    err_msg = 'h2op_convert_to_swd_layers: calc_potential has not been called'
    return
end if

nxE   = op%nxE
n_swd = op%n_swd
allocate(phi_nx(0:op%nx-1), phi_nxE_slice(0:nxE-1))

! Elevation spectral coefficients
call fft_real_to_swd(op%plan_nx, eta_nx, op%nx, n_swd, h_swd, err_msg)
if (err_msg /= '') return

! Layers m=1..nstep: phiE(:, nstep-m) is the H2-computed potential at
! sigma=(m-1)/nstep.  nstep-m runs from nstep-1 (sigma=0) down to 0 (sigma=0.95).
do m = 1, op%nstep
    phi_nxE_slice = op%phiE(:, op%nstep - m)
    call fft_resample_down(op%plan_nxE, op%plan_nx, &
                           phi_nxE_slice, nxE, op%nx, phi_nx, err_msg)
    if (err_msg /= '') return
    call fft_real_to_swd(op%plan_nx, phi_nx, op%nx, n_swd, c_swd(:, m), err_msg)
    if (err_msg /= '') return
end do

! Layer m=nstep+1: free-surface potential psi at sigma=1 (no resampling)
call fft_real_to_swd(op%plan_nx, psi_nx, op%nx, n_swd, c_swd(:, op%nstep + 1), err_msg)

deallocate(phi_nx, phi_nxE_slice)

end subroutine h2op_convert_to_swd_layers

!==============================================================================
! Private: H-operator (transfer potential from etaE to z=0 reference)
!==============================================================================

subroutine perform_h_operator(op, phi0E, etaE, &
                               rho_work, etaE_pow, &
                               spec_phi, spec_rho, k_power, err_msg)
! Apply the H-operator in-place on phi0E.
! phi0E(:,0) is the input surface potential; phi0E(:,1:M_kin-1) start as zero.
! After the call, phi0E(:,um) holds the um-th order correction.
!
! The Taylor expansion is: phi(z=0) = sum_{m=0}^{M-1} eta^m/m! * d^m phi/dz^m|_{z=eta}
! Vertical derivatives are applied spectrally: d^n phi/dz^n <-> |k|^n * phi-hat,
! with tanh(k*(h+zref)) applied for even-order terms (see amp2_free_surface_potential.rst).
! Note: the +psiE_M[um] term from the Python reference is always zero here because
! the calling convention sets psiE_M[um]=0 for all um>0.
type(h2op_state), intent(in)    :: op
real(dp),         intent(inout) :: phi0E(0:, 0:)   ! (nxE, 0:M_kin-1)
real(dp),         intent(in)    :: etaE(0:)         ! surface elevation
real(dp),         intent(inout) :: rho_work(0:)
real(dp),         intent(inout) :: etaE_pow(0:)
complex(dp),      intent(inout) :: spec_phi(0:), spec_rho(0:)
real(dp),         intent(inout) :: k_power(0:)
character(len=*), intent(out)   :: err_msg
!
integer :: um, un, j, nxE, facn
!
err_msg = ''
nxE = op%nxE
phi0E(:, 1:op%M_kin-1) = 0.0_dp   ! higher-order terms start zero
!
do um = 1, op%M_kin - 1
    k_power(:) = op%kE(:)           ! E = kE^1 at start of un loop
    etaE_pow(:) = etaE(:)           ! etaE^(un+1) = etaE^1 for un=0
    facn = 1
    do un = 0, um - 1
        ! rho = irfft( k_power * [thkhE *] rfft(phi0E(:,um-un-1)) )
        call fft_r2c(op%plan_nxE, phi0E(:, um-un-1), spec_phi, err_msg)
        if (err_msg /= '') return
        if (mod(un, 2) == 0) then
            do j = 0, nxE/2
                spec_rho(j) = cmplx(k_power(j) * op%thkhE(j), 0.0_dp, dp) * spec_phi(j)
            end do
        else
            do j = 0, nxE/2
                spec_rho(j) = cmplx(k_power(j), 0.0_dp, dp) * spec_phi(j)
            end do
        end if
        call fft_c2r(op%plan_nxE, spec_rho, rho_work, err_msg)
        if (err_msg /= '') return
        phi0E(:, um) = phi0E(:, um) - etaE_pow(:) * rho_work(:) / real(facn, dp)
        ! Prepare for next un iteration
        if (un < um - 1) then
            facn = facn * (un + 2)
            k_power(:)  = op%kE(:) * k_power(:)   ! kE^(un+2)
            etaE_pow(:) = etaE(:)  * etaE_pow(:)  ! etaE^(un+2)
        end if
    end do
end do
!
end subroutine perform_h_operator

!==============================================================================
! Private: H2-operator (transfer potential from eta to eta2)
!==============================================================================

subroutine perform_h2_operator(op, psi_prev, eta_prev, eta_next, &
                                phi0E, phi2E, rho_work, etaE_pow, &
                                spec_phi, spec_rho, k_power, err_msg)
! Apply the H2-operator:
!   1. H-operator from eta_prev to z=0 -> phi0E
!   2. Taylor expansion from z=0 to eta_next -> phi2E
! See amp2_free_surface_potential.rst for the mathematical description.
type(h2op_state), intent(in)    :: op
real(dp),         intent(in)    :: psi_prev(0:)   ! potential at previous layer
real(dp),         intent(in)    :: eta_prev(0:)   ! height of previous layer above zref
real(dp),         intent(in)    :: eta_next(0:)   ! height of current layer above zref
real(dp),         intent(inout) :: phi0E(0:, 0:)  ! (nxE, M_kin) scratch/output
real(dp),         intent(inout) :: phi2E(0:, 0:)  ! (nxE, M_kin) output
real(dp),         intent(inout) :: rho_work(0:)
real(dp),         intent(inout) :: etaE_pow(0:)
complex(dp),      intent(inout) :: spec_phi(0:), spec_rho(0:)
real(dp),         intent(inout) :: k_power(0:)
character(len=*), intent(out)   :: err_msg

integer :: um, un, j, nxE, facn

err_msg = ''
nxE = op%nxE

! Step 1: H-operator from eta_prev to z=0
phi0E(:, 0) = psi_prev(:)
call perform_h_operator(op, phi0E, eta_prev, &
                        rho_work, etaE_pow, spec_phi, spec_rho, k_power, err_msg)
if (err_msg /= '') return

! Step 2: Taylor expansion from z=0 to eta_next
phi2E(:, :) = phi0E(:, :)   ! initialise from phi0E

do um = 1, op%M_kin - 1
    k_power(:)  = op%kE(:)     ! kE^1
    etaE_pow(:) = eta_next(:)  ! eta_next^1
    facn = 1
    do un = 0, um - 1
        call fft_r2c(op%plan_nxE, phi0E(:, um-un-1), spec_phi, err_msg)
        if (err_msg /= '') return
        if (mod(un, 2) == 0) then
            do j = 0, nxE/2
                spec_rho(j) = cmplx(k_power(j) * op%thkhE(j), 0.0_dp, dp) * spec_phi(j)
            end do
        else
            do j = 0, nxE/2
                spec_rho(j) = cmplx(k_power(j), 0.0_dp, dp) * spec_phi(j)
            end do
        end if
        call fft_c2r(op%plan_nxE, spec_rho, rho_work, err_msg)
        if (err_msg /= '') return
        phi2E(:, um) = phi2E(:, um) + etaE_pow(:) * rho_work(:) / real(facn, dp)
        if (un < um - 1) then
            facn = facn * (un + 2)
            k_power(:)  = op%kE(:)    * k_power(:)
            etaE_pow(:) = eta_next(:) * etaE_pow(:)
        end if
    end do
end do

end subroutine perform_h2_operator

!==============================================================================

end module hosm_h2_operator_def
