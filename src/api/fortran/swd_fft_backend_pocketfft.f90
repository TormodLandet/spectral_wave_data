module swd_fft_backend
!
! ============================================================
! LAYER 3 (bottom) of the SWD FFT stack:
!
!   swd_fft.f90                  (Layer 1: SWD-specific stateful facade)
!     -> swd_fft_lib.f90         (Layer 2: shared library facade)
!        -> swd_fft_backend.f90   *** THIS FILE ***
!
! This is the PocketFFT backend.
! Concrete implementation lives in:
!   src/thirdparty/pocketfft/swd_pocketfft_c_api.cpp  (thin C wrapper)
!   src/thirdparty/pocketfft/pocketfft_hdronly.h      (vendored C++17, header-only)
!
! Responsibilities:
!   - Expose raw 1-D and 2-D real FFT primitives with a backend-neutral interface.
!   - 1-D operations are plan-based; the plan handle is an opaque C pointer.
!       backend_plan_alloc  / backend_plan_free
!       backend_r2c_1d  (scale = 1,   matches numpy rfft)
!       backend_c2r_1d  (scale = 1/n, matches numpy irfft)
!   - 2-D c2r is plan-less; the backend creates/caches plans internally.
!       backend_c2r_2d  (scale = 1, unnormalized — FFTW convention)
!       For ny = 1 this degenerates to a 1-D unnormalized transform.
!   - All transforms use Fortran column-major (contiguous in the first index).
!   - No SWD-specific knowledge; all spectral coefficient conventions handled above.
!
! Selected by: -DSWD_FFT_BACKEND=POCKETFFT  (the default)
! The mutually exclusive alternative is swd_fft_backend_fftw.f90.
! ============================================================

use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, &
                                       c_associated

implicit none
private

!##############################################################################
!                    P U B L I C    Q U A N T I T I E S
!##############################################################################

public :: backend_plan_alloc
public :: backend_plan_free
public :: backend_r2c_1d
public :: backend_c2r_1d
public :: backend_c2r_2d

!##############################################################################
!                  E N D    P U B L I C    Q U A N T I T I E S
!##############################################################################

!------------------------------------------------------------------------------
! C bindings to swd_pocketfft_c_api.cpp
!------------------------------------------------------------------------------
interface

    integer(c_int) function c_rfft_plan_create(n, plan_out) &
            bind(c, name='swd_rfft_plan_create')
        import c_int, c_ptr
        integer(c_int), value          :: n
        type(c_ptr),    intent(out)    :: plan_out
    end function

    subroutine c_rfft_plan_destroy(plan) &
            bind(c, name='swd_rfft_plan_destroy')
        import c_ptr
        type(c_ptr), value :: plan
    end subroutine

    ! 1-D r2c: real[n] -> complex[n/2+1], scale = 1
    integer(c_int) function c_rfft_forward(plan, real_in, cmplx_out) &
            bind(c, name='swd_rfft_forward')
        import c_int, c_ptr, c_double
        type(c_ptr),       value          :: plan
        real(c_double),    intent(in)     :: real_in(*)
        complex(c_double), intent(out)    :: cmplx_out(*)
    end function

    ! 1-D c2r: complex[n/2+1] -> real[n], scale = 1/n
    integer(c_int) function c_rfft_backward(plan, cmplx_in, real_out) &
            bind(c, name='swd_rfft_backward')
        import c_int, c_ptr, c_double
        type(c_ptr),       value          :: plan
        complex(c_double), intent(in)     :: cmplx_in(*)
        real(c_double),    intent(out)    :: real_out(*)
    end function

    ! 2-D c2r: complex[nx/2+1, ny] -> real[nx, ny], scale = 1 (unnormalized)
    ! Fortran column-major layout; for ny=1 degenerates to 1-D.
    integer(c_int) function c_rfft2_c2r(nx, ny, cmplx_in, real_out) &
            bind(c, name='swd_rfft2_c2r')
        import c_int, c_double
        integer(c_int),    value          :: nx, ny
        complex(c_double), intent(in)     :: cmplx_in(*)
        real(c_double),    intent(out)    :: real_out(*)
    end function

end interface

contains

!==============================================================================

subroutine backend_plan_alloc(handle, n, err_msg)
! Allocate a 1-D FFT plan for transforms of length n.
type(c_ptr),      intent(out) :: handle
integer,          intent(in)  :: n
character(len=*), intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_plan_create(int(n, c_int), handle)
if (ios /= 0) then
    write(err_msg, '(a,i0)') 'swd_fft_backend(pocketfft): plan_create failed for n=', n
    handle = c_null_ptr
end if
end subroutine backend_plan_alloc

!==============================================================================

subroutine backend_plan_free(handle)
! Free a plan previously allocated by backend_plan_alloc.  No-op for null handle.
type(c_ptr), intent(inout) :: handle
if (c_associated(handle)) then
    call c_rfft_plan_destroy(handle)
    handle = c_null_ptr
end if
end subroutine backend_plan_free

!==============================================================================

subroutine backend_r2c_1d(handle, real_in, cmplx_out, err_msg)
! Forward 1-D real-to-complex FFT.
! real_in(0:n-1) -> cmplx_out(0:n/2).  Scale = 1  (numpy.fft.rfft convention).
type(c_ptr),       intent(in)  :: handle
real(c_double),    intent(in)  :: real_in(0:)
complex(c_double), intent(out) :: cmplx_out(0:)
character(len=*),  intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_forward(handle, real_in, cmplx_out)
if (ios /= 0) err_msg = 'swd_fft_backend(pocketfft): rfft_forward failed'
end subroutine backend_r2c_1d

!==============================================================================

subroutine backend_c2r_1d(handle, cmplx_in, real_out, err_msg)
! Backward 1-D complex-to-real FFT.
! cmplx_in(0:n/2) -> real_out(0:n-1).  Scale = 1/n  (numpy.fft.irfft convention).
type(c_ptr),       intent(in)  :: handle
complex(c_double), intent(in)  :: cmplx_in(0:)
real(c_double),    intent(out) :: real_out(0:)
character(len=*),  intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft_backward(handle, cmplx_in, real_out)
if (ios /= 0) err_msg = 'swd_fft_backend(pocketfft): rfft_backward failed'
end subroutine backend_c2r_1d

!==============================================================================

subroutine backend_c2r_2d(nx, ny, cmplx_in, real_out, err_msg)
! 2-D complex-to-real inverse FFT.  Scale = 1 (unnormalized, FFTW convention).
! Input:  cmplx_in(nx/2+1, ny)  — Fortran column-major.
! Output: real_out(nx, ny)      — Fortran column-major.
! For ny = 1 this is equivalent to a 1-D unnormalized irfft.
! The caller is responsible for pre-scaling coefficients as needed.
integer,           intent(in)  :: nx, ny
complex(c_double), intent(in)  :: cmplx_in(*)
real(c_double),    intent(out) :: real_out(*)
character(len=*),  intent(out) :: err_msg
integer :: ios
err_msg = ''
ios = c_rfft2_c2r(int(nx, c_int), int(ny, c_int), cmplx_in, real_out)
if (ios /= 0) err_msg = 'swd_fft_backend(pocketfft): rfft2_c2r failed'
end subroutine backend_c2r_2d

!==============================================================================

end module swd_fft_backend
