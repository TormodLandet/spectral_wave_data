module swd_fft_lib
!
! ============================================================
! LAYER 2 (middle) of the SWD FFT stack:
!
!   swd_fft.f90                  (Layer 1: SWD-specific stateful facade)
!     -> swd_fft_lib.f90         *** THIS FILE ***
!        -> swd_fft_backend.f90  (Layer 3: thin PocketFFT or FFTW wrapper)
!
! Responsibilities:
!   - Own the swd_fft_plan type and its lifecycle (fft_init / fft_destroy).
!   - Provide the 1-D high-level API consumed by the H2 operator and the
!     write-side helpers:
!       fft_r2c, fft_c2r
!       fft_dealias
!       fft_resample_up, fft_resample_down
!       fft_swd_to_real, fft_real_to_swd
!   - Provide irfft2: the 2-D (and 1-D) real inverse FFT function consumed by
!     Layer 1 (Odin's swd_fft) for grid-based field evaluation.
!     Zero-padding / truncation (spectral resampling) is implemented here —
!     this is the only place in the stack where that logic lives.
!   - Dispatch all raw transforms to swd_fft_backend.  No direct calls to
!     PocketFFT or FFTW appear in this module.
!
! Consumers:
!   H2 operator:      hosm_h2_operator.f90
!                     spectral_wave_data_shape_1_or_2_impl_7.f90
!   Odin's facade:    swd_fft.f90  (via `use swd_fft_lib, only: irfft2`)
!
! Conventions (matching numpy.fft):
!   fft_r2c / irfft2 r2c side : scale = 1    (rfft  equivalent)
!   fft_c2r                   : scale = 1/n  (irfft equivalent)
!   irfft2 c2r side           : scale = 1    (unnormalized; caller pre-scales coeffs)
!
! SWD spectral coefficient convention for shape 1/2:
!   The file stores h_j such that
!     eta(x_m) = real(h_0) + sum_{j=1}^{n-1} real(h_j * exp(-i j dk x_m))
!              + real(h_n) * cos(n dk x_m)
!   with
!     h_0 = mean(eta)          [DC term — assumed zero for ocean wave data]
!     h_j = (2/N) * conj(rfft(eta)[j])   for j=1..n-1
!     h_n = rfft(eta)[n] / N             (Nyquist; real-valued)
!   where N = 2*n is the physical grid size.
! ============================================================

use swd_fft_backend, only: backend_plan1d_alloc, backend_plan1d_free, &
                            backend_r2c_1d, backend_c2r_1d, &
                            backend_plan2d_alloc, backend_plan2d_free, backend_c2r_2d

use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_null_ptr

implicit none
private

integer, parameter :: dp = c_double   ! always float64, regardless of kind_swd_internal

!##############################################################################
!                    P U B L I C    Q U A N T I T I E S
!##############################################################################

public :: swd_fft_plan
public :: fft_init, fft_destroy
public :: fft_r2c, fft_c2r
public :: fft_dealias
public :: fft_resample_up, fft_resample_down
public :: fft_swd_to_real, fft_real_to_swd
public :: irfft2

!##############################################################################
!                  E N D    P U B L I C    Q U A N T I T I E S
!##############################################################################

!------------------------------------------------------------------------------
! Plan type — wraps the opaque backend handle plus n for convenience.
!------------------------------------------------------------------------------
type :: swd_fft_plan
    type(c_ptr) :: handle = c_null_ptr
    integer     :: n      = 0
end type swd_fft_plan

contains

!==============================================================================
! fft_init / fft_destroy
!==============================================================================

subroutine fft_init(plan, n, err_msg)
! Create a plan for 1-D real FFTs of length n.
type(swd_fft_plan), intent(out) :: plan
integer,            intent(in)  :: n
character(len=*),   intent(out) :: err_msg
err_msg = ''
call backend_plan1d_alloc(plan%handle, n, err_msg)
if (err_msg == '') plan%n = n
end subroutine fft_init

!==============================================================================

subroutine fft_destroy(plan)
! Free the plan.  Safe to call on an already-destroyed plan.
type(swd_fft_plan), intent(inout) :: plan
call backend_plan1d_free(plan%handle)   ! no-op for null handle
plan%n = 0
end subroutine fft_destroy

!==============================================================================
! fft_r2c / fft_c2r
!==============================================================================

subroutine fft_r2c(plan, real_in, cmplx_out, err_msg)
! Forward real-to-complex FFT.  real_in(0:n-1) -> cmplx_out(0:n/2).
! Scale = 1  (matches numpy.fft.rfft).
type(swd_fft_plan), intent(in)  :: plan
real(dp),           intent(in)  :: real_in(0:)
complex(dp),        intent(out) :: cmplx_out(0:)
character(len=*),   intent(out) :: err_msg
err_msg = ''
call backend_r2c_1d(plan%handle, real_in, cmplx_out, err_msg)
end subroutine fft_r2c

!==============================================================================

subroutine fft_c2r(plan, cmplx_in, real_out, err_msg)
! Backward complex-to-real FFT.  cmplx_in(0:n/2) -> real_out(0:n-1).
! Scale = 1/n  (matches numpy.fft.irfft).
type(swd_fft_plan), intent(in)  :: plan
complex(dp),        intent(in)  :: cmplx_in(0:)
real(dp),           intent(out) :: real_out(0:)
character(len=*),   intent(out) :: err_msg
err_msg = ''
call backend_c2r_1d(plan%handle, cmplx_in, real_out, err_msg)
end subroutine fft_c2r

!==============================================================================
! Spectral operations
!==============================================================================

subroutine fft_dealias(plan_hi, signal, n_hi, n_lo, err_msg)
! Dealias a real signal on the extended grid by zeroing spectral modes
! above n_lo/2.  signal(:) is modified in-place.
! n_hi = plan_hi%n, n_lo = nx (original grid size, n_lo < n_hi).
type(swd_fft_plan), intent(in)    :: plan_hi
real(dp),           intent(inout) :: signal(0:)
integer,            intent(in)    :: n_hi, n_lo
character(len=*),   intent(out)   :: err_msg
!
complex(dp), allocatable :: spec(:)
integer :: nhalf_hi, nhalf_lo
!
err_msg = ''
nhalf_hi = n_hi / 2
nhalf_lo = n_lo / 2
allocate(spec(0:nhalf_hi))
!
call fft_r2c(plan_hi, signal, spec, err_msg)
if (err_msg /= '') return
!
! Zero modes above n_lo/2
spec(nhalf_lo+1 : nhalf_hi) = cmplx(0.0_dp, 0.0_dp, dp)
! Ensure the Nyquist mode of the lower grid is real (for even n_lo)
spec(nhalf_lo) = cmplx(real(spec(nhalf_lo), dp), 0.0_dp, dp)
!
call fft_c2r(plan_hi, spec, signal, err_msg)
!
deallocate(spec)
end subroutine fft_dealias

!==============================================================================

subroutine fft_resample_up(plan_lo, plan_hi, sig_lo, n_lo, n_hi, sig_hi, err_msg)
! Upsample a real signal from n_lo to n_hi grid points.
! Amplitude-preserving: sig_hi samples the same continuous signal as sig_lo.
! Requires n_hi > n_lo, both powers of two.
type(swd_fft_plan), intent(in)  :: plan_lo, plan_hi
real(dp),           intent(in)  :: sig_lo(0:)
integer,            intent(in)  :: n_lo, n_hi
real(dp),           intent(out) :: sig_hi(0:)
character(len=*),   intent(out) :: err_msg
!
complex(dp), allocatable :: spec_lo(:), spec_hi(:)
real(dp) :: scale
integer  :: nhalf_lo, nhalf_hi
!
err_msg = ''
nhalf_lo = n_lo / 2
nhalf_hi = n_hi / 2
allocate(spec_lo(0:nhalf_lo), spec_hi(0:nhalf_hi))
!
! Forward on small grid
call fft_r2c(plan_lo, sig_lo, spec_lo, err_msg)
if (err_msg /= '') return
!
! Scale up: DFT amplitude of larger grid = (n_hi/n_lo) * DFT of smaller grid
scale = real(n_hi, dp) / real(n_lo, dp)
spec_hi(0:nhalf_lo) = spec_lo(0:nhalf_lo) * scale
! Nyquist of small grid: must be purely real after splitting
spec_hi(nhalf_lo) = cmplx(real(spec_hi(nhalf_lo), dp), 0.0_dp, dp)
spec_hi(nhalf_lo+1 : nhalf_hi) = cmplx(0.0_dp, 0.0_dp, dp)
!
! Inverse on large grid
call fft_c2r(plan_hi, spec_hi, sig_hi, err_msg)
!
deallocate(spec_lo, spec_hi)
end subroutine fft_resample_up

!==============================================================================

subroutine fft_resample_down(plan_hi, plan_lo, sig_hi, n_hi, n_lo, sig_lo, err_msg)
! Downsample a real signal from n_hi to n_lo grid points.
! Amplitude-preserving: keeps only frequency content up to n_lo/2 * dk.
! Requires n_hi > n_lo, both powers of two.
type(swd_fft_plan), intent(in)  :: plan_hi, plan_lo
real(dp),           intent(in)  :: sig_hi(0:)
integer,            intent(in)  :: n_hi, n_lo
real(dp),           intent(out) :: sig_lo(0:)
character(len=*),   intent(out) :: err_msg
!
complex(dp), allocatable :: spec_hi(:), spec_lo(:)
real(dp) :: scale
integer  :: nhalf_lo, nhalf_hi
!
err_msg = ''
nhalf_lo = n_lo / 2
nhalf_hi = n_hi / 2
allocate(spec_hi(0:nhalf_hi), spec_lo(0:nhalf_lo))
!
! Forward on large grid
call fft_r2c(plan_hi, sig_hi, spec_hi, err_msg)
if (err_msg /= '') return
!
! Truncate and scale down
scale = real(n_lo, dp) / real(n_hi, dp)
spec_lo(0:nhalf_lo) = spec_hi(0:nhalf_lo) * scale
spec_lo(nhalf_lo)   = cmplx(real(spec_lo(nhalf_lo), dp), 0.0_dp, dp)
!
! Inverse on small grid
call fft_c2r(plan_lo, spec_lo, sig_lo, err_msg)
!
deallocate(spec_hi, spec_lo)
end subroutine fft_resample_down

!==============================================================================
! SWD spectral coefficient <-> real-space conversions
!==============================================================================

subroutine fft_swd_to_real(plan, h_swd, n_swd, n_real, real_out, err_msg)
! Convert SWD spectral coefficients h_swd(0:n_swd) to a real-space array
! real_out(0:n_real-1), where n_real = 2*n_swd (even, power of two required).
!
! SWD convention:
!   h(0)     = mean(eta)  [DC term — assumed zero for ocean wave data; handled
!                           correctly for non-zero values via N*h(0) scaling]
!   h(j)     = 2/N * conj(rfft(eta)[j])  for j=1..n_swd-1
!   h(n_swd) = rfft(eta)[n_swd] / N      [Nyquist, real-only]
! where N = n_real = 2*n_swd.
!
! Inverse: rfft(eta)[j]     = (N/2) * conj(h(j))  for j=1..n_swd-1
!          rfft(eta)[0]      = N * h(0)
!          rfft(eta)[n_swd] = N * h(n_swd)
type(swd_fft_plan), intent(in)  :: plan
complex(dp),        intent(in)  :: h_swd(0:)
integer,            intent(in)  :: n_swd, n_real
real(dp),           intent(out) :: real_out(0:)
character(len=*),   intent(out) :: err_msg
!
complex(dp), allocatable :: H(:)  ! rfft coefficients for irfft
real(dp) :: half_n
integer  :: j
!
err_msg = ''
allocate(H(0:n_swd))
half_n = real(n_real, dp) * 0.5_dp
!
! Build rfft input from SWD coefficients
H(0) = cmplx(real(n_real, dp) * real(h_swd(0), dp), 0.0_dp, dp)
do j = 1, n_swd - 1
    H(j) = half_n * conjg(h_swd(j))
end do
H(n_swd) = cmplx(real(n_real, dp) * real(h_swd(n_swd), dp), 0.0_dp, dp)
!
call fft_c2r(plan, H, real_out, err_msg)
!
deallocate(H)
end subroutine fft_swd_to_real

!==============================================================================

subroutine fft_real_to_swd(plan, real_in, n_real, n_swd, h_swd, err_msg)
! Convert a real-space array real_in(0:n_real-1) to SWD spectral coefficients
! h_swd(0:n_swd), where n_swd = n_real/2.
!
! This is the inverse of fft_swd_to_real.
type(swd_fft_plan), intent(in)  :: plan
real(dp),           intent(in)  :: real_in(0:)
integer,            intent(in)  :: n_real, n_swd
complex(dp),        intent(out) :: h_swd(0:)
character(len=*),   intent(out) :: err_msg
!
complex(dp), allocatable :: H(:)
real(dp) :: inv_n, inv_half_n
integer  :: j
!
err_msg = ''
allocate(H(0:n_swd))
inv_n      = 1.0_dp / real(n_real, dp)
inv_half_n = 2.0_dp * inv_n   ! = 2/n_real
!
call fft_r2c(plan, real_in, H, err_msg)
if (err_msg /= '') return
!
h_swd(0) = cmplx(real(H(0), dp) * inv_n, 0.0_dp, dp)
do j = 1, n_swd - 1
    h_swd(j) = inv_half_n * conjg(H(j))
end do
h_swd(n_swd) = cmplx(real(H(n_swd), dp) * inv_n, 0.0_dp, dp)
!
deallocate(H)
end subroutine fft_real_to_swd

!==============================================================================
! irfft2 — 2-D (and 1-D) real inverse FFT with spectral resampling
!==============================================================================

function irfft2(fh, nx_in, ny_in, nx_out, ny_out) result(fE)
! 2-D real inverse FFT of the complex half-spectrum fh(nx_in/2+1, ny_in),
! with optional zero-padding or truncation to produce an (nx_out, ny_out) result.
!
! Scale = 1 (unnormalized, matching FFTW's c2r convention).
! The caller is responsible for pre-scaling spectral coefficients.
!
! For ny_in = ny_out = 1 this degenerates to a 1-D unnormalized irfft.
!
! Note: the DC component fh(1,1) (zero wavenumber) is assumed to be zero
! in ocean wave applications (no mean water-level shift in the wave field).
! The code handles non-zero DC correctly, but it is not tested.
!
! Based on Odin Gramstad's irfft2 from swd_fft_fftw3.f90; ported to use
! the backend-neutral swd_fft_backend interface so both PocketFFT and FFTW
! are supported via a single copy of this logic.
complex(dp), intent(in) :: fh(nx_in/2 + 1, ny_in)
integer,     intent(in) :: nx_in, ny_in, nx_out, ny_out
real(dp), allocatable   :: fE(:, :)
!
complex(dp), allocatable :: fhE(:, :)
complex(dp), allocatable :: fhInt(:, :)
type(c_ptr) :: plan2d
character(len=200) :: err_msg
integer :: nx_int, ny_int
!
allocate(fE(nx_out, ny_out))
!
! Allocate a 2-D plan for the output grid.
! NOTE (temporary): this stateless facade function currently creates and frees
! the plan on every call.  Once Odin's stateful swd_fft object (Layer 1) is
! merged it will own and reuse the plan across calls; the plan-based backend
! interface (backend_plan2d_alloc/backend_c2r_2d/backend_plan2d_free) is already
! in place for that.
call backend_plan2d_alloc(plan2d, nx_out, ny_out, err_msg)
if (err_msg /= '') error stop 'irfft2: backend_plan2d_alloc failed'
!
! Fast path: same size — no resampling needed
if (nx_in == nx_out .and. ny_in == ny_out) then
    call backend_c2r_2d(plan2d, fh, fE, err_msg)
    if (err_msg /= '') error stop 'irfft2: backend_c2r_2d failed (same-size path)'
    call backend_plan2d_free(plan2d)
    return
end if
!
allocate(fhE(nx_out/2 + 1, ny_out))
!
if (nx_in <= nx_out .and. ny_in <= ny_out) then
    ! Upsample (zero-pad) in all dimensions
    call swd_zeropad(fh, nx_in, ny_in, fhE, nx_out, ny_out)
elseif (nx_in >= nx_out .and. ny_in >= ny_out) then
    ! Downsample (truncate) in all dimensions
    call swd_truncate(fh, nx_in, ny_in, fhE, nx_out, ny_out)
else
    ! Mixed: upsample in one dimension, downsample in the other.
    ! Go through an intermediate grid that is large enough in both dimensions.
    nx_int = max(nx_out, nx_in)
    ny_int = max(ny_out, ny_in)
    allocate(fhInt(nx_int/2 + 1, ny_int))
    call swd_zeropad(fh, nx_in, ny_in, fhInt, nx_int, ny_int)
    call swd_truncate(fhInt, nx_int, ny_int, fhE, nx_out, ny_out)
    deallocate(fhInt)
end if
!
call backend_c2r_2d(plan2d, fhE, fE, err_msg)
if (err_msg /= '') error stop 'irfft2: backend_c2r_2d failed'
!
call backend_plan2d_free(plan2d)
deallocate(fhE)
end function irfft2

!==============================================================================
! Private spectral resampling helpers
! Based on Odin Gramstad's zeropad / truncate from swd_fft_fftw3.f90; ported
! to subroutine form and placed here so the logic is shared by both backends.
! Nyquist-frequency splitting/combining is handled correctly (verified against
! scipy.signal.resample in Odin's original test suite).
! DC is assumed zero for ocean wave data; it is copied correctly regardless.
!==============================================================================

subroutine swd_zeropad(fh, nx, ny, fhE, nxE, nyE)
! Zero-pad fh(nx/2+1, ny) into fhE(nxE/2+1, nyE) for upsampling (nxE>=nx, nyE>=ny).
! Corresponds to evaluating the same continuous signal on a finer spatial grid.
complex(dp), intent(in)  :: fh(nx/2+1, ny)
integer,     intent(in)  :: nx, ny
complex(dp), intent(out) :: fhE(nxE/2+1, nyE)
integer,     intent(in)  :: nxE, nyE
!
integer :: nxh, nyh, nxhE
!
nxh  = nx  / 2 + 1
nxhE = nxE / 2 + 1
nyh  = ny  / 2 + 1
!
fhE = cmplx(0.0_dp, 0.0_dp, dp)
!
if (ny == 1 .and. nyE == 1) then
    ! 1-D case
    fhE(1:nxh, 1) = fh(1:nxh, 1)
    if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
        ! Even nx: the Nyquist mode must be split when moving to a larger grid
        ! to avoid double-counting at both +nx/2 and -nx/2.
        fhE(nxh, 1) = 0.5_dp * fhE(nxh, 1)
    end if
else
    ! 2-D case: copy the low-frequency block and the negative-ky block
    fhE(1:nxh, 1:nyh) = fh(1:nxh, 1:nyh)
    if (mod(ny, 2) == 0) then
        fhE(1:nxh, nyE - nyh + 3 : nyE) = fh(1:nxh, nyh + 1 : ny)
    else
        fhE(1:nxh, nyE - nyh + 2 : nyE) = fh(1:nxh, nyh + 1 : ny)
    end if
    ! Nyquist in x: split across the zero-padding boundary
    if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
        fhE(nxh, :) = 0.5_dp * fhE(nxh, :)
    end if
    ! Nyquist in y: split into both Nyquist slots (even ny only)
    if (mod(ny, 2) == 0 .and. ny /= nyE) then
        fhE(:, nyh) = 0.5_dp * fhE(:, nyh)
        fhE(:, nyE - nyh + 2) = fhE(:, nyh)
    end if
end if
end subroutine swd_zeropad

!==============================================================================

subroutine swd_truncate(fhE, nxE, nyE, fh, nx, ny)
! Truncate fhE(nxE/2+1, nyE) into fh(nx/2+1, ny) for downsampling (nxE>=nx, nyE>=ny).
! Keeps only the frequency content within the bandwidth of the coarser grid.
complex(dp), intent(in)  :: fhE(nxE/2+1, nyE)
integer,     intent(in)  :: nxE, nyE
complex(dp), intent(out) :: fh(nx/2+1, ny)
integer,     intent(in)  :: nx, ny
!
integer :: nxh, nyh, nxhE
!
nxh  = nx  / 2 + 1
nxhE = nxE / 2 + 1
nyh  = ny  / 2 + 1
!
fh = cmplx(0.0_dp, 0.0_dp, dp)
!
if (ny == 1 .and. nyE == 1) then
    ! 1-D case
    fh(1:nxh, 1) = fhE(1:nxh, 1)
    if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
        ! Re-combine the split Nyquist terms from the larger grid
        fh(nxh, 1) = 2.0_dp * real(fh(nxh, 1), dp)
    end if
else
    ! 2-D case
    fh(1:nxh, 1:nyh) = fhE(1:nxh, 1:nyh)
    if (mod(ny, 2) == 0) then
        fh(1:nxh, nyh + 1 : ny) = fhE(1:nxh, nyE - nyh + 3 : nyE)
    else
        fh(1:nxh, nyh + 1 : ny) = fhE(1:nxh, nyE - nyh + 2 : nyE)
    end if
    ! Re-combine Nyquist in y (even ny only).
    ! zeropad() split the Nyquist y-column into two halved copies; reconstruct
    ! the original by doubling.  The column is real-valued by Hermitian symmetry.
    if (mod(ny, 2) == 0 .and. ny /= nyE) then
        fh(1:nxh, nyh) = 2.0_dp * real(fh(1:nxh, nyh), dp)
    end if
    ! Re-combine Nyquist in x (even nx only)
    if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
        fh(nxh, 1) = 2.0_dp * real(fh(nxh, 1), dp)
        fh(nxh, 2:ny) = fh(nxh, 2:ny) + conjg(fh(nxh, ny:2:-1))
    end if
end if
end subroutine swd_truncate

!==============================================================================

end module swd_fft_lib
