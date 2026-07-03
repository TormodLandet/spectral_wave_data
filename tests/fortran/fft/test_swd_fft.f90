program test_swd_fft
! Unit tests for the swd_fft_def module.
!
! Tests:
!   1. fft_r2c / fft_c2r round-trip recovers the original signal
!   2. fft_resample_up followed by fft_resample_down recovers the original signal
!   3. fft_swd_to_real / fft_real_to_swd round-trip recovers the coefficients
!   4. fft_dealias zeroes spectral modes above the cutoff
!
! Exit code: 0 = all passed, 1 = at least one failure.

use, intrinsic :: iso_c_binding, only: c_double
use swd_fft_def, only: swd_fft_plan, fft_init, fft_destroy, &
                       fft_r2c, fft_c2r, fft_dealias, &
                       fft_resample_up, fft_resample_down, &
                       fft_swd_to_real, fft_real_to_swd

implicit none

integer, parameter :: dp = c_double
integer :: failures
character(len=200) :: err_msg

failures = 0

call test_r2c_c2r_roundtrip(failures, err_msg)
call test_resample_roundtrip(failures, err_msg)
call test_swd_real_roundtrip(failures, err_msg)
call test_dealias(failures, err_msg)

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

end program test_swd_fft
