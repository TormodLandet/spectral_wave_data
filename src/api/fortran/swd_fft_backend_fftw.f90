module swd_fft_backend
!
! ============================================================
! LAYER 3 (bottom) of the SWD FFT stack:
!
!   swd_fft.f90                  (Layer 1: SWD-specific stateful facade)
!     -> swd_fft_lib.f90         (Layer 2: shared library facade)
!        -> swd_fft_backend.f90   *** THIS FILE ***
!
! This is the FFTW3 backend.  It IS the concrete implementation: it calls the
! system (or vendored) FFTW3 library directly through the standard Fortran
! interface `include 'fftw3.f03'`.  There is no separate C shim.
!
! Build / linking (see src/api/fortran/Cmake/swd_fft.cmake):
!   - Selected by:  -DSWD_FFT_BACKEND=FFTW
!   - Requires:     libfftw3 (double precision) + fftw3.f03 on the include path
!                   (RHEL: dnf install fftw-devel, Debian: apt install libfftw3-dev)
!   - NOTE: linking against GPL FFTW makes the resulting SWD library GPL, not
!           MIT.  The distributed Python wheels therefore use the PocketFFT
!           backend (swd_fft_backend_pocketfft.f90).  This file exists mainly so
!           that FFTW/MKL can be used in-house and to prevent bit-rot.
!
! Responsibilities (identical public interface to the PocketFFT backend):
!   - 1-D operations are plan-based; the plan handle is an opaque C pointer
!     wrapping an FFTW r2c plan, an FFTW c2r plan and scratch work buffers.
!       backend_plan_alloc  / backend_plan_free
!       backend_r2c_1d  (scale = 1,   matches numpy rfft)
!       backend_c2r_1d  (scale = 1/n, matches numpy irfft — FFTW is unnormalized,
!                        so this backend divides by n explicitly)
!   - 2-D c2r is plan-less at the interface; FFTW plans are cached internally in
!     a module-level cache keyed by (nx, ny) and reused across calls for speed.
!       backend_c2r_2d  (scale = 1, unnormalized — the FFTW convention; the
!                        caller pre-scales coefficients)
!       For ny = 1 this degenerates to a 1-D unnormalized transform.
!   - All transforms use Fortran column-major (contiguous in the first index).
!   - No SWD-specific knowledge; spectral coefficient conventions handled above.
!
! Thread-safety: like the rest of SWD (and the PocketFFT backend), this backend
! assumes single-threaded use per process.  FFTW plan creation and the module
! plan cache below are NOT thread-safe; plan *execution* is, but we do not rely
! on that here.
!
! The mutually exclusive alternative is swd_fft_backend_pocketfft.f90.
! ============================================================

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

public :: backend_plan_alloc
public :: backend_plan_free
public :: backend_r2c_1d
public :: backend_c2r_1d
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
! 2-D c2r plan cache.  FFTW needs an explicit plan per (nx, ny); we build it
! lazily on first use and reuse it on subsequent calls with the same size.
! Persistent work buffers keep the plan valid (fixed arrays/alignment) and give
! us a safe place to copy the caller's input into — FFTW's multi-dimensional
! c2r transform DESTROYS its input array, so we must never transform the
! caller's (intent(in)) data in place.
!------------------------------------------------------------------------------
type :: fftw_cache_2d
    integer                        :: nx   = 0
    integer                        :: ny   = 0
    type(c_ptr)                    :: plan = c_null_ptr
    complex(c_double), allocatable :: cin(:,:)   ! (nx/2+1, ny) complex scratch
    real(c_double),    allocatable :: rout(:,:)  ! (nx, ny)     real scratch
end type fftw_cache_2d

type(fftw_cache_2d), allocatable, save :: cache2d(:)
integer,                          save :: n_cache2d = 0

contains

!==============================================================================
! 1-D plan lifecycle
!==============================================================================

subroutine backend_plan_alloc(handle, n, err_msg)
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
end subroutine backend_plan_alloc

!==============================================================================

subroutine backend_plan_free(handle)
! Free a plan previously allocated by backend_plan_alloc.  No-op for null handle.
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
end subroutine backend_plan_free

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

subroutine backend_c2r_2d(nx, ny, cmplx_in, real_out, err_msg)
! 2-D complex-to-real inverse FFT.  Scale = 1 (unnormalized, FFTW convention).
! Input:  cmplx_in(nx/2+1, ny)  — Fortran column-major (contiguous in x).
! Output: real_out(nx, ny)      — Fortran column-major.
! For ny = 1 this is equivalent to a 1-D unnormalized irfft.
! The caller is responsible for pre-scaling coefficients as needed.
integer,           intent(in)  :: nx, ny
complex(c_double), intent(in)  :: cmplx_in(*)
real(c_double),    intent(out) :: real_out(*)
character(len=*),  intent(out) :: err_msg
!
integer :: idx, nxh, ncplx, nreal
!
err_msg = ''
if (nx <= 0 .or. ny <= 0) then
    write(err_msg, '(a,i0,a,i0)') 'swd_fft_backend(fftw): invalid 2-D size nx=', &
        nx, ' ny=', ny
    return
end if
!
nxh   = nx/2 + 1
ncplx = nxh * ny
nreal = nx  * ny
!
call get_cache2d_entry(nx, ny, idx, err_msg)
if (err_msg /= '') return
!
! Copy the caller's spectrum into the scratch input buffer.  FFTW's 2-D c2r
! overwrites its input, so we must not transform cmplx_in directly.
cache2d(idx)%cin = reshape(cmplx_in(1:ncplx), [nxh, ny])
!
call fftw_execute_dft_c2r(cache2d(idx)%plan, cache2d(idx)%cin, cache2d(idx)%rout)
!
real_out(1:nreal) = reshape(cache2d(idx)%rout, [nreal])
end subroutine backend_c2r_2d

!==============================================================================

subroutine get_cache2d_entry(nx, ny, idx, err_msg)
! Return the index of the cache entry for (nx, ny), creating it (and its FFTW
! plan + work buffers) on first use.
integer,          intent(in)  :: nx, ny
integer,          intent(out) :: idx
character(len=*), intent(out) :: err_msg
!
type(fftw_cache_2d), allocatable :: tmp(:)
integer :: i, ios, nxh
!
err_msg = ''
idx     = 0
!
! Look for an existing plan of this size.
do i = 1, n_cache2d
    if (cache2d(i)%nx == nx .and. cache2d(i)%ny == ny) then
        idx = i
        return
    end if
end do
!
! Grow the cache array if necessary.
if (.not. allocated(cache2d)) then
    allocate(cache2d(4), stat=ios)
    if (ios /= 0) then
        err_msg = 'swd_fft_backend(fftw): failed to allocate 2-D plan cache'
        return
    end if
else if (n_cache2d == size(cache2d)) then
    allocate(tmp(2*size(cache2d)), stat=ios)
    if (ios /= 0) then
        err_msg = 'swd_fft_backend(fftw): failed to grow 2-D plan cache'
        return
    end if
    tmp(1:n_cache2d) = cache2d(1:n_cache2d)
    call move_alloc(tmp, cache2d)
end if
!
n_cache2d = n_cache2d + 1
idx       = n_cache2d
nxh       = nx/2 + 1
!
cache2d(idx)%nx = nx
cache2d(idx)%ny = ny
allocate(cache2d(idx)%cin(nxh, ny), cache2d(idx)%rout(nx, ny), stat=ios)
if (ios /= 0) then
    err_msg = 'swd_fft_backend(fftw): failed to allocate 2-D work buffers'
    n_cache2d = n_cache2d - 1
    return
end if
!
! Column-major Fortran arrays cin(nx/2+1, ny) and rout(nx, ny) are, in C
! row-major terms, [ny][nx/2+1] complex and [ny][nx] real.  FFTW's c2r treats
! the LAST (contiguous) axis as the real-transform axis, so the logical FFTW
! dimensions are (n0, n1) = (ny, nx): a c2c along y and a c2r along x.
cache2d(idx)%plan = fftw_plan_dft_c2r_2d(int(ny, c_int), int(nx, c_int),  &
                                         cache2d(idx)%cin, cache2d(idx)%rout, &
                                         ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
if (.not. c_associated(cache2d(idx)%plan)) then
    err_msg = 'swd_fft_backend(fftw): FFTW failed to create 2-D c2r plan'
    deallocate(cache2d(idx)%cin, cache2d(idx)%rout)
    n_cache2d = n_cache2d - 1
    return
end if
end subroutine get_cache2d_entry

!==============================================================================

end module swd_fft_backend
