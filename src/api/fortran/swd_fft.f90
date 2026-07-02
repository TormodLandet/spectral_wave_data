module swd_fft_def

! Fortran facade for the pocketfft C++ library (and maybe FFTW3/MKL in the future)
!
! Provides 1-D real-FFT operations on float64/complex128 arrays via
! iso_c_binding calls to pocketfft_c_api.cpp.
!
! All arrays use Fortran-0-based indexing: real(0:n-1), complex(0:n/2).
!
! Conventions (matching numpy.fft):
!   fft_r2c : scale = 1    (rfft  equivalent)
!   fft_c2r : scale = 1/n  (irfft equivalent)
!
! SWD spectral coefficient convention for shape 1/2:
!   The file stores h_j such that
!     eta(x_m) = real(h_0) + sum_{j=1}^{n-1} real(h_j * exp(-i j dk x_m))
!              + real(h_n) * cos(n dk x_m)
!   with
!     h_0 = mean(eta)
!     h_j = (2/N_real) * conj(rfft(eta)[j])   for j=1..n-1
!     h_n = rfft(eta)[n] / N_real              (Nyquist; real-valued)
!   where N_real = 2*n is the physical grid size.

use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, &
                                       c_associated

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

!##############################################################################
!                  E N D    P U B L I C    Q U A N T I T I E S
!##############################################################################

!------------------------------------------------------------------------------
! Plan type — wraps the opaque C pointer plus n for convenience.
!------------------------------------------------------------------------------
type :: swd_fft_plan
    type(c_ptr) :: handle = c_null_ptr
    integer     :: n      = 0
end type swd_fft_plan

!------------------------------------------------------------------------------
! C interface bindings (private)
!------------------------------------------------------------------------------
interface
    integer(c_int) function c_rfft_plan_create(n, plan_out) &
            bind(c, name='swd_rfft_plan_create')
        import c_int, c_ptr
        integer(c_int), value :: n
        type(c_ptr), intent(out) :: plan_out
    end function

    subroutine c_rfft_plan_destroy(plan) &
            bind(c, name='swd_rfft_plan_destroy')
        import c_ptr
        type(c_ptr), value :: plan
    end subroutine

    ! out: complex(c_double)(0:n/2)  — passed as void* on C side
    integer(c_int) function c_rfft_forward(plan, real_in, cmplx_out) &
            bind(c, name='swd_rfft_forward')
        import c_int, c_ptr, c_double
        type(c_ptr),    value       :: plan
        real(c_double), intent(in)  :: real_in(*)
        complex(c_double), intent(out) :: cmplx_out(*)
    end function

    ! in: complex(c_double)(0:n/2)  — passed as const void* on C side
    integer(c_int) function c_rfft_backward(plan, cmplx_in, real_out) &
            bind(c, name='swd_rfft_backward')
        import c_int, c_ptr, c_double
        type(c_ptr),    value      :: plan
        complex(c_double), intent(in)  :: cmplx_in(*)
        real(c_double), intent(out) :: real_out(*)
    end function
end interface

contains

!==============================================================================
! fft_init / fft_destroy
!==============================================================================

subroutine fft_init(plan, n, err_msg)
! Create a plan for 1-D real FFTs of length n.
type(swd_fft_plan), intent(out) :: plan
integer,            intent(in)  :: n
character(len=*),   intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_plan_create(int(n, c_int), plan%handle)
if (ios /= 0) then
    write(err_msg,'(a,i0)') 'swd_fft: failed to create FFT plan for n=', n
    return
end if
plan%n = n
end subroutine fft_init

!==============================================================================

subroutine fft_destroy(plan)
type(swd_fft_plan), intent(inout) :: plan
if (c_associated(plan%handle)) then
    call c_rfft_plan_destroy(plan%handle)
    plan%handle = c_null_ptr
    plan%n      = 0
end if
end subroutine fft_destroy

!==============================================================================
! fft_r2c / fft_c2r
!==============================================================================

subroutine fft_r2c(plan, real_in, cmplx_out, err_msg)
! Forward real-to-complex FFT.  real_in(0:n-1) -> cmplx_out(0:n/2).
! Scale = 1 (matches numpy.fft.rfft).
type(swd_fft_plan), intent(in)  :: plan
real(dp),           intent(in)  :: real_in(0:)
complex(dp),        intent(out) :: cmplx_out(0:)
character(len=*),   intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_forward(plan%handle, real_in, cmplx_out)
if (ios /= 0) err_msg = 'swd_fft: rfft_forward failed'
end subroutine fft_r2c

!==============================================================================

subroutine fft_c2r(plan, cmplx_in, real_out, err_msg)
! Backward complex-to-real FFT.  cmplx_in(0:n/2) -> real_out(0:n-1).
! Scale = 1/n (matches numpy.fft.irfft).
type(swd_fft_plan), intent(in)  :: plan
complex(dp),        intent(in)  :: cmplx_in(0:)
real(dp),           intent(out) :: real_out(0:)
character(len=*),   intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_backward(plan%handle, cmplx_in, real_out)
if (ios /= 0) err_msg = 'swd_fft: rfft_backward failed'
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
!   h(0)     = mean(eta)  [real-only stored]
!   h(j)     = 2/N * conj(rfft(eta)[j])  for j=1..n_swd-1
!   h(n_swd) = rfft(eta)[n_swd] / N      [Nyquist, real-only]
! where N = n_real = 2*n_swd.
!
! Inverse: rfft(eta)[j] = (N/2) * conj(h(j))  for j=1..n_swd-1
!          rfft(eta)[0]  = N * h(0)
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

end module swd_fft_def
