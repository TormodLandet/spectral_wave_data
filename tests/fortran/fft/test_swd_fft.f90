program test_swd_fft
! Unit tests for the SWD FFT stack (swd_fft_lib + the compiled swd_fft_backend).
!
! These tests are backend-independent: they pass with either
! -DSWD_FFT_BACKEND=POCKETFFT or -DSWD_FFT_BACKEND=FFTW, and the 2-D tests use a
! brute-force DFT as ground truth so they also verify that the two backends
! agree with each other (and with the documented axis/scale conventions).
!
! Tests:
!   1. fft_r2c / fft_c2r round-trip recovers the original signal
!   2. fft_resample_up followed by fft_resample_down recovers the original signal
!   3. fft_swd_to_real / fft_real_to_swd round-trip recovers the coefficients
!   4. fft_dealias zeroes spectral modes above the cutoff
!   5. backend_c2r_2d matches a brute-force unnormalized inverse DFT (nx /= ny)
!   6. backend_c2r_2d with ny = 1 matches the brute-force 1-D case
!   7. irfft2 up-sampling (zero-pad) reproduces the field at coincident points
!   8. irfft2 down-sampling (truncate) reproduces a band-limited field
!
! Exit code: 0 = all passed, 1 = at least one failure.

use, intrinsic :: iso_c_binding, only: c_double
use swd_fft_lib,     only: swd_fft_plan, fft_init, fft_destroy, &
                           fft_r2c, fft_c2r, fft_dealias, &
                           fft_resample_up, fft_resample_down, &
                           fft_swd_to_real, fft_real_to_swd, irfft2
use swd_fft_backend, only: backend_c2r_2d

implicit none

integer, parameter :: dp = c_double
real(dp), parameter :: twopi = 2.0_dp * 3.14159265358979323846_dp
integer :: failures
character(len=200) :: err_msg

failures = 0

call test_r2c_c2r_roundtrip(failures, err_msg)
call test_resample_roundtrip(failures, err_msg)
call test_swd_real_roundtrip(failures, err_msg)
call test_dealias(failures, err_msg)
call test_backend_c2r_2d(failures, err_msg)
call test_backend_c2r_2d_ny1(failures, err_msg)
call test_irfft2_zeropad(failures, err_msg)
call test_irfft2_truncate(failures, err_msg)

if (failures == 0) then
    write(*,'(a)') 'ALL FFT UNIT TESTS PASSED'
    stop 0
else
    write(*,'(a,i0,a)') failures, ' FFT UNIT TEST(S) FAILED'
    stop 1
end if

contains

!------------------------------------------------------------------------------

subroutine check(label, cond, failures)
character(len=*), intent(in)  :: label
logical,          intent(in)  :: cond
integer,          intent(inout) :: failures
if (.not. cond) then
    write(*,'(a,a)') 'FAIL: ', trim(label)
    failures = failures + 1
else
    write(*,'(a,a)') 'pass: ', trim(label)
end if
end subroutine check

!------------------------------------------------------------------------------
! Test 1: fft_r2c -> fft_c2r round-trip
!------------------------------------------------------------------------------
subroutine test_r2c_c2r_roundtrip(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer,  parameter :: n = 64
real(dp) :: x_orig(0:n-1), x_back(0:n-1)
complex(dp) :: spec(0:n/2)
type(swd_fft_plan) :: plan
integer :: i
real(dp) :: max_err

call fft_init(plan, n, err_msg)
call check('fft_init n=64', err_msg == '', failures)
if (err_msg /= '') return

! Construct a signal with a few harmonics
do i = 0, n-1
    x_orig(i) = 1.5_dp * cos(2*3.14159265_dp*real(3*i,dp)/n) &
               + 0.7_dp * sin(2*3.14159265_dp*real(7*i,dp)/n)
end do

call fft_r2c(plan, x_orig, spec, err_msg)
call check('fft_r2c', err_msg == '', failures)
call fft_c2r(plan, spec, x_back, err_msg)
call check('fft_c2r', err_msg == '', failures)

max_err = maxval(abs(x_back - x_orig))
call check('r2c->c2r max error < 1e-12', max_err < 1.0e-12_dp, failures)

call fft_destroy(plan)
end subroutine test_r2c_c2r_roundtrip

!------------------------------------------------------------------------------
! Test 2: fft_resample_up -> fft_resample_down round-trip
!------------------------------------------------------------------------------
subroutine test_resample_roundtrip(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer,  parameter :: n_lo = 32, n_hi = 128
real(dp) :: sig_lo(0:n_lo-1), sig_hi(0:n_hi-1), sig_back(0:n_lo-1)
type(swd_fft_plan) :: plan_lo, plan_hi
integer :: i
real(dp) :: max_err, pi

pi = 3.14159265358979323846_dp
call fft_init(plan_lo, n_lo, err_msg)
call check('resample fft_init lo', err_msg == '', failures)
call fft_init(plan_hi, n_hi, err_msg)
call check('resample fft_init hi', err_msg == '', failures)
if (err_msg /= '') return

! Signal band-limited to n_lo/2 modes
do i = 0, n_lo-1
    sig_lo(i) = cos(2*pi*real(2*i,dp)/n_lo) + 0.5_dp*sin(2*pi*real(5*i,dp)/n_lo)
end do

call fft_resample_up(plan_lo, plan_hi, sig_lo, n_lo, n_hi, sig_hi, err_msg)
call check('fft_resample_up', err_msg == '', failures)
call fft_resample_down(plan_hi, plan_lo, sig_hi, n_hi, n_lo, sig_back, err_msg)
call check('fft_resample_down', err_msg == '', failures)

max_err = maxval(abs(sig_back - sig_lo))
call check('resample round-trip max error < 1e-12', max_err < 1.0e-12_dp, failures)

call fft_destroy(plan_lo)
call fft_destroy(plan_hi)
end subroutine test_resample_roundtrip

!------------------------------------------------------------------------------
! Test 3: fft_swd_to_real -> fft_real_to_swd round-trip
!------------------------------------------------------------------------------
subroutine test_swd_real_roundtrip(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: n_swd = 16, n_real = 32
complex(dp) :: h_orig(0:n_swd), h_back(0:n_swd)
real(dp) :: eta(0:n_real-1)
type(swd_fft_plan) :: plan
integer :: j
real(dp) :: max_err

call fft_init(plan, n_real, err_msg)
call check('swd round-trip fft_init', err_msg == '', failures)
if (err_msg /= '') return

! Construct synthetic SWD coefficients (real h_0, complex interior, real h_n)
h_orig(0)     = cmplx(0.1_dp,  0.0_dp, dp)
h_orig(n_swd) = cmplx(0.02_dp, 0.0_dp, dp)
do j = 1, n_swd-1
    h_orig(j) = cmplx(0.3_dp/j, -0.15_dp/j, dp)
end do

call fft_swd_to_real(plan, h_orig, n_swd, n_real, eta, err_msg)
call check('fft_swd_to_real', err_msg == '', failures)
call fft_real_to_swd(plan, eta, n_real, n_swd, h_back, err_msg)
call check('fft_real_to_swd', err_msg == '', failures)

max_err = maxval(abs(h_back - h_orig))
call check('swd->real->swd max error < 1e-12', max_err < 1.0e-12_dp, failures)

call fft_destroy(plan)
end subroutine test_swd_real_roundtrip

!------------------------------------------------------------------------------
! Test 4: fft_dealias zeroes modes above n_lo/2
!------------------------------------------------------------------------------
subroutine test_dealias(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: n_hi = 64, n_lo = 16
real(dp) :: sig(0:n_hi-1)
complex(dp) :: spec(0:n_hi/2)
type(swd_fft_plan) :: plan_hi, plan_lo
integer :: i
real(dp) :: energy_high, pi

pi = 3.14159265358979323846_dp
call fft_init(plan_hi, n_hi, err_msg)
call check('dealias fft_init hi', err_msg == '', failures)
call fft_init(plan_lo, n_lo, err_msg)
call check('dealias fft_init lo', err_msg == '', failures)
if (err_msg /= '') return

! Signal with energy above n_lo/2
do i = 0, n_hi-1
    sig(i) = sin(2*pi*real(20*i,dp)/n_hi) + sin(2*pi*real(30*i,dp)/n_hi)
end do

call fft_dealias(plan_hi, sig, n_hi, n_lo, err_msg)
call check('fft_dealias', err_msg == '', failures)

! After dealiasing, the forward transform should have near-zero power above mode n_lo/2
call fft_r2c(plan_hi, sig, spec, err_msg)
energy_high = 0.0_dp
do i = n_lo/2 + 1, n_hi/2
    energy_high = energy_high + abs(spec(i))**2
end do
call check('dealias zeroes high modes (energy < 1e-20)', energy_high < 1.0e-20_dp, failures)

call fft_destroy(plan_hi)
call fft_destroy(plan_lo)
end subroutine test_dealias

!------------------------------------------------------------------------------
! Ground-truth helper: unnormalized forward 2-D DFT of a real field, returning
! the half spectrum chat(0:nx/2, 0:ny-1) in the convention used by the backend:
!   real axis = x = first (column-major) index, reduced to nx/2+1; y = full.
!   chat(kx,ky) = sum_{px,py} f(px,py) * exp(-2*pi*i (kx*px/nx + ky*py/ny))
! The matching unnormalized inverse (backend_c2r_2d) therefore returns nx*ny*f.
!------------------------------------------------------------------------------
subroutine forward_half_2d(f, nx, ny, chat)
integer,     intent(in)  :: nx, ny
real(dp),    intent(in)  :: f(0:nx-1, 0:ny-1)
complex(dp), intent(out) :: chat(0:nx/2, 0:ny-1)
integer :: kx, ky, px, py
real(dp) :: ang
complex(dp) :: acc
do ky = 0, ny-1
    do kx = 0, nx/2
        acc = cmplx(0.0_dp, 0.0_dp, dp)
        do py = 0, ny-1
            do px = 0, nx-1
                ang = -twopi * (real(kx,dp)*px/nx + real(ky,dp)*py/ny)
                acc = acc + f(px,py) * cmplx(cos(ang), sin(ang), dp)
            end do
        end do
        chat(kx,ky) = acc
    end do
end do
end subroutine forward_half_2d

!------------------------------------------------------------------------------
! Test 5: backend_c2r_2d matches the brute-force unnormalized inverse DFT.
! Uses nx /= ny so that any accidental x/y transpose is caught.
!------------------------------------------------------------------------------
subroutine test_backend_c2r_2d(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: nx = 8, ny = 6
real(dp)    :: f(0:nx-1, 0:ny-1), g(nx, ny)
complex(dp) :: chat(0:nx/2, 0:ny-1)
integer  :: px, py
real(dp) :: max_err

err_msg = ''
do py = 0, ny-1
    do px = 0, nx-1
        f(px,py) = 1.3_dp + 0.7_dp*cos(twopi*(1.0_dp*px/nx))            &
                 + 0.5_dp*sin(twopi*(2.0_dp*py/ny))                     &
                 + 0.4_dp*cos(twopi*(1.0_dp*px/nx + 1.0_dp*py/ny))      &
                 + 0.2_dp*cos(twopi*(3.0_dp*px/nx - 2.0_dp*py/ny))
    end do
end do

call forward_half_2d(f, nx, ny, chat)
call backend_c2r_2d(nx, ny, chat, g, err_msg)
call check('backend_c2r_2d no error', err_msg == '', failures)

max_err = 0.0_dp
do py = 0, ny-1
    do px = 0, nx-1
        max_err = max(max_err, abs(g(px+1,py+1) - real(nx*ny,dp)*f(px,py)))
    end do
end do
call check('backend_c2r_2d matches brute-force DFT (err < 1e-8)', &
           max_err < 1.0e-8_dp, failures)
end subroutine test_backend_c2r_2d

!------------------------------------------------------------------------------
! Test 6: backend_c2r_2d degenerates correctly to a 1-D transform for ny = 1.
!------------------------------------------------------------------------------
subroutine test_backend_c2r_2d_ny1(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: nx = 8, ny = 1
real(dp)    :: f(0:nx-1, 0:ny-1), g(nx, ny)
complex(dp) :: chat(0:nx/2, 0:ny-1)
integer  :: px
real(dp) :: max_err

err_msg = ''
do px = 0, nx-1
    f(px,0) = 0.9_dp + 0.6_dp*cos(twopi*(1.0_dp*px/nx))    &
            + 0.3_dp*sin(twopi*(3.0_dp*px/nx))
end do

call forward_half_2d(f, nx, ny, chat)
call backend_c2r_2d(nx, ny, chat, g, err_msg)
call check('backend_c2r_2d(ny=1) no error', err_msg == '', failures)

max_err = 0.0_dp
do px = 0, nx-1
    max_err = max(max_err, abs(g(px+1,1) - real(nx*ny,dp)*f(px,0)))
end do
call check('backend_c2r_2d(ny=1) matches brute-force DFT (err < 1e-8)', &
           max_err < 1.0e-8_dp, failures)
end subroutine test_backend_c2r_2d_ny1

!------------------------------------------------------------------------------
! Test 7: irfft2 up-sampling (zero-pad).  Doubling nx and ny must reproduce the
! same continuous field; the coarse-grid samples coincide with every other
! fine-grid sample.  ny is even, exercising the Nyquist-splitting path.
!------------------------------------------------------------------------------
subroutine test_irfft2_zeropad(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: nxc = 8, nyc = 6, nxf = 16, nyf = 12
real(dp)    :: f(0:nxc-1, 0:nyc-1), base(nxc, nyc)
complex(dp) :: chat(0:nxc/2, 0:nyc-1)
real(dp), allocatable :: fine(:,:)
integer  :: px, py
real(dp) :: max_err

err_msg = ''
do py = 0, nyc-1
    do px = 0, nxc-1
        f(px,py) = 1.1_dp + 0.6_dp*cos(twopi*(1.0_dp*px/nxc))          &
                 + 0.4_dp*sin(twopi*(2.0_dp*py/nyc))                   &
                 + 0.3_dp*cos(twopi*(2.0_dp*px/nxc + 1.0_dp*py/nyc))
    end do
end do

call forward_half_2d(f, nxc, nyc, chat)
call backend_c2r_2d(nxc, nyc, chat, base, err_msg)
call check('irfft2 zeropad: base transform no error', err_msg == '', failures)

fine = irfft2(chat, nxc, nyc, nxf, nyf)
call check('irfft2 zeropad: output shape', &
           size(fine,1) == nxf .and. size(fine,2) == nyf, failures)

! Coarse point (px,py) coincides with fine point (2*px, 2*py).
max_err = 0.0_dp
do py = 0, nyc-1
    do px = 0, nxc-1
        max_err = max(max_err, abs(fine(2*px+1, 2*py+1) - base(px+1, py+1)))
    end do
end do
call check('irfft2 zeropad matches at coincident points (err < 1e-8)', &
           max_err < 1.0e-8_dp, failures)
end subroutine test_irfft2_zeropad

!------------------------------------------------------------------------------
! Test 8: irfft2 down-sampling (truncate).  A field band-limited below both the
! coarse and fine Nyquist frequencies is sampled losslessly on the coarse grid.
! nx and ny are both even on the fine grid, exercising the Nyquist-recombine
! path in swd_truncate (the line fixed during the facade refactor).
!------------------------------------------------------------------------------
subroutine test_irfft2_truncate(failures, err_msg)
integer,          intent(inout) :: failures
character(len=*), intent(out)   :: err_msg

integer, parameter :: nxb = 16, nyb = 12, nxs = 8, nys = 6
real(dp)    :: fb(0:nxb-1, 0:nyb-1), base(nxb, nyb)
complex(dp) :: chat(0:nxb/2, 0:nyb-1)
real(dp), allocatable :: small(:,:)
integer  :: px, py
real(dp) :: max_err

err_msg = ''
! Band-limited to |kx| <= 3 < nxs/2 and |ky| <= 2 < nys/2 so truncation is exact.
do py = 0, nyb-1
    do px = 0, nxb-1
        fb(px,py) = 1.0_dp + 0.7_dp*cos(twopi*(1.0_dp*px/nxb))         &
                  + 0.5_dp*sin(twopi*(3.0_dp*px/nxb))                  &
                  + 0.4_dp*cos(twopi*(2.0_dp*py/nyb))                  &
                  + 0.3_dp*cos(twopi*(3.0_dp*px/nxb + 2.0_dp*py/nyb))
    end do
end do

call forward_half_2d(fb, nxb, nyb, chat)
call backend_c2r_2d(nxb, nyb, chat, base, err_msg)
call check('irfft2 truncate: base transform no error', err_msg == '', failures)

small = irfft2(chat, nxb, nyb, nxs, nys)
call check('irfft2 truncate: output shape', &
           size(small,1) == nxs .and. size(small,2) == nys, failures)

! Small point (px,py) coincides with big point (2*px, 2*py).
max_err = 0.0_dp
do py = 0, nys-1
    do px = 0, nxs-1
        max_err = max(max_err, abs(small(px+1, py+1) - base(2*px+1, 2*py+1)))
    end do
end do
call check('irfft2 truncate matches at coincident points (err < 1e-8)', &
           max_err < 1.0e-8_dp, failures)
end subroutine test_irfft2_truncate

end program test_swd_fft
