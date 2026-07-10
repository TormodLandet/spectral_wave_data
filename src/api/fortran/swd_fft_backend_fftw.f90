module swd_fft_backend
!
! =============================================================================
! LAYER 3 (bottom) of the SWD FFT stack:
!
!   swd_fft.f90                  (Layer 1: SWD-specific stateful facade)
!     -> swd_fft_lib.f90         (Layer 2: shared library facade)
!        -> swd_fft_backend.f90   *** THIS FILE ***
!
! This is the FFTW3 backend.  It IS the concrete implementation: it calls the
! system FFTW3 library (or MKL's impl.) directly through the standard Fortran
! interface `include 'fftw3.f03'`.  There is no separate C shim like there is
! for the PocketFFT implemantation since FFTW3 has a Fortran 2003 interface.
!
! Build / linking (see src/api/fortran/Cmake/swd_fft.cmake):
!   - Selected by:  -DSWD_FFT_BACKEND=FFTW
!   - Requires:     libfftw3 (double precision) + fftw3.f03 on the include path
!                   or the FFTW3-compatible implementation in Intel MKL.
!   - NOTE: linking against "normal" FFTW makes the resulting SWD library GPL,
!           not MIT.  The distributed Python wheels therefore use the PocketFFT
!           backend (swd_fft_backend_pocketfft.f90).  This file exists mainly so
!           that GPL FFTW can be used in-house if wanted or with MKL or FFTW3
!           if an appropriate licence is held (FFTW3 has a commercial offering
!           for those who do not want to use a GPL-licenced library).
!
! Responsibilities (identical public interface to the PocketFFT backend):
!   - Both 1-D and 2-D operations are plan-based; a plan handle is an opaque C
!     pointer, owned and freed by the caller (no global/module state here).
!       backend_plan1d_alloc / backend_plan1d_free  (wraps FFTW r2c + c2r plans
!                                                    and scratch work buffers)
!       backend_r2c_1d  (scale = 1,   matches numpy rfft)
!       backend_c2r_1d  (scale = 1/n, matches numpy irfft — FFTW is unnormalized,
!                        so this backend divides by n explicitly)
!       backend_plan2d_alloc / backend_plan2d_free  (wraps the FFTW 2-D c2r plan
!                                                    and its work buffers)
!       backend_c2r_2d  (scale = 1, unnormalized — the FFTW convention; the
!                        caller pre-scales coefficients)
!       For ny = 1 the 2-D transform degenerates to a 1-D unnormalized transform.
!   - All transforms use Fortran column-major (contiguous in the first index).
!   - No SWD-specific knowledge; spectral coefficient conventions handled on
!     layers above this backend layer.
!
! Thread-safety: like the rest of SWD (and the PocketFFT backend), this backend
! assumes single-threaded use per process.  FFTW plan creation is not
! thread-safe; plan *execution* is, but we do not rely on that here.  There is
! no module-level mutable state — every plan lives in a caller-owned handle.
!
! The mutually exclusive alternative is swd_fft_backend_pocketfft.f90.
! =============================================================================

! The FFTW3 Fortran interface (fftw3.f03) references many iso_c_binding kinds
! (c_int32_t, c_intptr_t, c_size_t, c_funptr, ...), so import the whole module
! rather than a restricted `only:` list.
use, intrinsic :: iso_c_binding

implicit none
private

! FFTW3 standard Fortran interface (plan constructors, executors, constants).
include 'fftw3.f03'

!##############################################################################
!                    P U B L I C    Q U A N T I T I E S
!##############################################################################

public :: backend_plan1d_alloc
public :: backend_plan1d_free
public :: backend_r2c_1d
public :: backend_c2r_1d
public :: backend_plan2d_alloc
public :: backend_plan2d_free
public :: backend_c2r_2d

!##############################################################################
!                  E N D    P U B L I C    Q U A N T I T I E S
!##############################################################################

!------------------------------------------------------------------------------
! 1-D plan state.  A c_ptr to one of these is handed back as the opaque handle.
! We keep persistent scratch buffers so the plans are created once (with a fixed
! alignment) and merely re-used via new-array-free execute on the same buffers.
!------------------------------------------------------------------------------
type :: fftw_state_1d
    integer                        :: n        = 0
    type(c_ptr)                    :: plan_r2c = c_null_ptr
    type(c_ptr)                    :: plan_c2r = c_null_ptr
    real(c_double),    allocatable :: rbuf(:)   ! (n)      real scratch
    complex(c_double), allocatable :: cbuf(:)   ! (n/2+1)  complex scratch
end type fftw_state_1d

!------------------------------------------------------------------------------
! 2-D c2r plan state.  A c_ptr to one of these is the opaque 2-D plan handle.
! Persistent work buffers keep the plan valid (fixed arrays/alignment) and give
! us a safe place to copy the caller's input into — FFTW's multi-dimensional
! c2r transform DESTROYS its input array, so we must never transform the
! caller's (intent(in)) data in place.
!------------------------------------------------------------------------------
type :: fftw_state_2d
    integer                        :: nx   = 0
    integer                        :: ny   = 0
    type(c_ptr)                    :: plan = c_null_ptr
    complex(c_double), allocatable :: cin(:,:)   ! (nx/2+1, ny) complex scratch
    real(c_double),    allocatable :: rout(:,:)  ! (nx, ny)     real scratch
end type fftw_state_2d

contains

!==============================================================================
! 1-D plan lifecycle
!==============================================================================

subroutine backend_plan1d_alloc(handle, n, err_msg)
! Allocate a 1-D FFT plan (r2c + c2r) for transforms of length n.
type(c_ptr),      intent(out) :: handle
integer,          intent(in)  :: n
character(len=*), intent(out) :: err_msg
!
type(fftw_state_1d), pointer :: st
integer :: ios
!
err_msg = ''
handle  = c_null_ptr
!
if (n <= 0) then
    write(err_msg, '(a,i0)') 'swd_fft_backend(fftw): invalid n=', n
    return
end if
!
allocate(st, stat=ios)
if (ios /= 0) then
    err_msg = 'swd_fft_backend(fftw): failed to allocate plan state'
    return
end if
!
st%n = n
allocate(st%rbuf(n), st%cbuf(n/2 + 1), stat=ios)
if (ios /= 0) then
    err_msg = 'swd_fft_backend(fftw): failed to allocate scratch buffers'
    deallocate(st)
    return
end if
!
! FFTW_ESTIMATE does not touch the buffers at plan time, so the (uninitialised)
! scratch arrays are fine here.  FFTW_UNALIGNED lets us execute on these
! allocatable buffers without assuming a particular SIMD alignment.
st%plan_r2c = fftw_plan_dft_r2c_1d(int(n, c_int), st%rbuf, st%cbuf, &
                                   ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
st%plan_c2r = fftw_plan_dft_c2r_1d(int(n, c_int), st%cbuf, st%rbuf, &
                                   ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
!
if (.not. c_associated(st%plan_r2c) .or. .not. c_associated(st%plan_c2r)) then
    err_msg = 'swd_fft_backend(fftw): FFTW failed to create 1-D plan(s)'
    if (c_associated(st%plan_r2c)) call fftw_destroy_plan(st%plan_r2c)
    if (c_associated(st%plan_c2r)) call fftw_destroy_plan(st%plan_c2r)
    deallocate(st)
    return
end if
!
handle = c_loc(st)
end subroutine backend_plan1d_alloc

!==============================================================================

subroutine backend_plan1d_free(handle)
! Free a 1-D plan previously allocated by backend_plan1d_alloc.  No-op if null.
type(c_ptr), intent(inout) :: handle
!
type(fftw_state_1d), pointer :: st
!
if (.not. c_associated(handle)) return
call c_f_pointer(handle, st)
if (c_associated(st%plan_r2c)) call fftw_destroy_plan(st%plan_r2c)
if (c_associated(st%plan_c2r)) call fftw_destroy_plan(st%plan_c2r)
if (allocated(st%rbuf)) deallocate(st%rbuf)
if (allocated(st%cbuf)) deallocate(st%cbuf)
deallocate(st)
handle = c_null_ptr
end subroutine backend_plan1d_free

!==============================================================================
! 1-D transforms
!==============================================================================

subroutine backend_r2c_1d(handle, real_in, cmplx_out, err_msg)
! Forward 1-D real-to-complex FFT.
! real_in(0:n-1) -> cmplx_out(0:n/2).  Scale = 1  (numpy.fft.rfft convention).
type(c_ptr),       intent(in)  :: handle
real(c_double),    intent(in)  :: real_in(0:)
complex(c_double), intent(out) :: cmplx_out(0:)
character(len=*),  intent(out) :: err_msg
!
type(fftw_state_1d), pointer :: st
!
err_msg = ''
if (.not. c_associated(handle)) then
    err_msg = 'swd_fft_backend(fftw): rfft_forward called with null plan'
    return
end if
call c_f_pointer(handle, st)
!
st%rbuf(:) = real_in(0:st%n - 1)
call fftw_execute_dft_r2c(st%plan_r2c, st%rbuf, st%cbuf)
cmplx_out(0:st%n/2) = st%cbuf(:)
end subroutine backend_r2c_1d

!==============================================================================

subroutine backend_c2r_1d(handle, cmplx_in, real_out, err_msg)
! Backward 1-D complex-to-real FFT.
! cmplx_in(0:n/2) -> real_out(0:n-1).  Scale = 1/n  (numpy.fft.irfft convention).
! FFTW's c2r is unnormalized, so we divide by n here to match the PocketFFT
! backend / numpy.  FFTW also destroys the complex input, hence the copy into
! the scratch buffer (which keeps cmplx_in intent(in) intact).
type(c_ptr),       intent(in)  :: handle
complex(c_double), intent(in)  :: cmplx_in(0:)
real(c_double),    intent(out) :: real_out(0:)
character(len=*),  intent(out) :: err_msg
!
type(fftw_state_1d), pointer :: st
real(c_double) :: inv_n
!
err_msg = ''
if (.not. c_associated(handle)) then
    err_msg = 'swd_fft_backend(fftw): rfft_backward called with null plan'
    return
end if
call c_f_pointer(handle, st)
!
st%cbuf(:) = cmplx_in(0:st%n/2)
call fftw_execute_dft_c2r(st%plan_c2r, st%cbuf, st%rbuf)
inv_n = 1.0_c_double / real(st%n, c_double)
real_out(0:st%n - 1) = st%rbuf(:) * inv_n
end subroutine backend_c2r_1d

!==============================================================================
! 2-D transform (plan-less interface, internally cached FFTW plans)
!==============================================================================

subroutine backend_c2r_2d(handle, cmplx_in, real_out, err_msg)
! 2-D complex-to-real inverse FFT.  Scale = 1 (unnormalized, FFTW convention).
! Input:  cmplx_in(nx/2+1, ny)  — Fortran column-major (contiguous in x).
! Output: real_out(nx, ny)      — Fortran column-major.
! nx, ny are recorded in the plan created by backend_plan2d_alloc.
! For ny = 1 this is equivalent to a 1-D unnormalized irfft.
! The caller is responsible for pre-scaling coefficients as needed.
type(c_ptr),       intent(in)  :: handle
complex(c_double), intent(in)  :: cmplx_in(*)
real(c_double),    intent(out) :: real_out(*)
character(len=*),  intent(out) :: err_msg
!
type(fftw_state_2d), pointer :: st
integer :: nxh, ncplx, nreal
!
err_msg = ''
if (.not. c_associated(handle)) then
    err_msg = 'swd_fft_backend(fftw): c2r_2d called with null plan'
    return
end if
call c_f_pointer(handle, st)
!
nxh   = st%nx/2 + 1
ncplx = nxh   * st%ny
nreal = st%nx * st%ny
!
! Copy the caller's spectrum into the scratch input buffer.  FFTW's 2-D c2r
! overwrites its input, so we must not transform cmplx_in directly.
st%cin = reshape(cmplx_in(1:ncplx), [nxh, st%ny])
!
call fftw_execute_dft_c2r(st%plan, st%cin, st%rout)
!
real_out(1:nreal) = reshape(st%rout, [nreal])
end subroutine backend_c2r_2d

!==============================================================================

subroutine backend_plan2d_alloc(handle, nx, ny, err_msg)
! Allocate a 2-D c2r plan (+ work buffers) for size (nx, ny).
type(c_ptr),      intent(out) :: handle
integer,          intent(in)  :: nx, ny
character(len=*), intent(out) :: err_msg
!
type(fftw_state_2d), pointer :: st
integer :: ios
!
err_msg = ''
handle  = c_null_ptr
!
if (nx <= 0 .or. ny <= 0) then
    write(err_msg, '(a,i0,a,i0)') &
        'swd_fft_backend(fftw): invalid 2-D size nx=', nx, ' ny=', ny
    return
end if
!
allocate(st, stat=ios)
if (ios /= 0) then
    err_msg = 'swd_fft_backend(fftw): failed to allocate 2-D plan state'
    return
end if
!
st%nx = nx
st%ny = ny
allocate(st%cin(nx/2 + 1, ny), st%rout(nx, ny), stat=ios)
if (ios /= 0) then
    err_msg = 'swd_fft_backend(fftw): failed to allocate 2-D work buffers'
    deallocate(st)
    return
end if
!
! Column-major Fortran arrays cin(nx/2+1, ny) and rout(nx, ny) are, in C
! row-major terms, [ny][nx/2+1] complex and [ny][nx] real.  FFTW's c2r treats
! the LAST (contiguous) axis as the real-transform axis, so the logical FFTW
! dimensions are (n0, n1) = (ny, nx): a c2c along y and a c2r along x.
st%plan = fftw_plan_dft_c2r_2d(int(ny, c_int), int(nx, c_int),  &
                               st%cin, st%rout,                  &
                               ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
if (.not. c_associated(st%plan)) then
    err_msg = 'swd_fft_backend(fftw): FFTW failed to create 2-D c2r plan'
    deallocate(st%cin, st%rout)
    deallocate(st)
    return
end if
!
handle = c_loc(st)
end subroutine backend_plan2d_alloc

!==============================================================================

subroutine backend_plan2d_free(handle)
! Free a 2-D plan previously allocated by backend_plan2d_alloc.  No-op if null.
type(c_ptr), intent(inout) :: handle
!
type(fftw_state_2d), pointer :: st
!
if (.not. c_associated(handle)) return
call c_f_pointer(handle, st)
if (c_associated(st%plan)) call fftw_destroy_plan(st%plan)
if (allocated(st%cin))  deallocate(st%cin)
if (allocated(st%rout)) deallocate(st%rout)
deallocate(st)
handle = c_null_ptr
end subroutine backend_plan2d_free

!==============================================================================

end module swd_fft_backend
