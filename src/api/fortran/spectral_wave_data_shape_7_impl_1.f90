module spectral_wave_data_shape_7_impl_1_def

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention, &
                             swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use spectral_interpolation_def, only: spectral_interpolation
use swd_version, only: version

implicit none
private

! This module provides an implementation of shape class 7 of the
! spectral-wave-data API.
!
! Shape 7: long-crested waves with a sigma-coordinate multi-layer potential.
!
! The velocity potential is stored on nlayers sigma-layers that deform with the
! free surface.  No time-derivative arrays are stored in the file; they are
! computed numerically from the four-step temporal window.
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_shape_7_impl_1
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type, extends(spectral_wave_data) :: spectral_wave_data_shape_7_impl_1
    integer            :: n              ! Highest spectral index (n+1 components, j=0..n)
    integer            :: nsumx          ! Number of components included in summation (<=n)
    real(wp)           :: dk             ! Constant spacing of wave numbers
    real(wp)           :: d              ! Water depth (or -1 for infinite depth)
    real(wp)           :: zref           ! z-position of the bottom sigma-layer (sigma=0)
    real(wp)           :: tanhdkd_eff    ! tanh(dk*(d+zref)) for below-zref extrapolation
    integer            :: nlayers        ! Number of sigma-layers
    real(wp), allocatable :: sig(:)      ! (1:nlayers) sigma-layer positions; sig(1)=0, sig(nlayers)=1
    integer            :: icur           ! Column index of most recently read data (1:4)
    integer            :: istp           ! Most recent SWD step in the four-step window
    integer            :: ipt(4,4)       ! Circular-buffer column mapping
    complex(c_float), allocatable :: h_win(:,:)   ! (0:n, 4)  elevation window
    complex(c_float), allocatable :: c_win(:,:,:) ! (0:n, nlayers, 4)  potential window
    complex(wp), allocatable :: h_cur(:)          ! (0:n) interpolated elevation
    complex(wp), allocatable :: ht_cur(:)         ! (0:n) interpolated d/dt elevation
    complex(wp), allocatable :: c_cur(:,:)        ! (0:n, nlayers) interpolated potential per layer
    complex(wp), allocatable :: ct_cur(:,:)       ! (0:n, nlayers) interpolated d/dt potential per layer
    type(spectral_interpolation) :: tpol          ! Temporal interpolation scheme
contains
    procedure :: close              ! Destructor
    procedure :: update_time        ! Obtain spectral data for current time
    procedure :: phi                ! Velocity potential at (x,y,z)
    procedure :: stream             ! Stream function (not supported, returns 0)
    procedure :: phi_t              ! Euler time derivative of potential
    procedure :: grad_phi           ! Particle velocity (gradient of potential)
    procedure :: grad_phi_2nd       ! Second order spatial gradients (stub)
    procedure :: acc_euler          ! Euler acceleration
    procedure :: acc_particle       ! Particle acceleration
    procedure :: elev               ! Surface elevation
    procedure :: elev_t             ! d/dt of surface elevation
    procedure :: grad_elev          ! Gradient of surface elevation
    procedure :: grad_elev_2nd      ! Second order gradients of elevation
    procedure :: pressure           ! Fully nonlinear Bernoulli pressure
    procedure :: bathymetry         ! Local water depth
    procedure :: bathymetry_nvec    ! Normal vector of sea floor
    procedure :: convergence        ! Convergence CSV (stub)
    procedure :: strip              ! Strip to time window (stub)
    procedure :: get_int            ! Extract integer parameter
    procedure :: get_logical        ! Extract logical parameter
    procedure :: get_real           ! Extract real parameter
    procedure :: get_chr            ! Extract char parameter
end type spectral_wave_data_shape_7_impl_1

interface spectral_wave_data_shape_7_impl_1
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_7_impl_1) :: self
!
logical :: opened
!
inquire(unit=self % unit, opened=opened)
if (opened) close(self % unit)
if (allocated(self % cid))    deallocate(self % cid)
if (allocated(self % sig))    deallocate(self % sig)
if (allocated(self % h_win))  deallocate(self % h_win)
if (allocated(self % c_win))  deallocate(self % c_win)
if (allocated(self % h_cur))  deallocate(self % h_cur)
if (allocated(self % ht_cur)) deallocate(self % ht_cur)
if (allocated(self % c_cur))  deallocate(self % c_cur)
if (allocated(self % ct_cur)) deallocate(self % ct_cur)
!
self % file = '0'
self % unit = 0
self % n = 0
self % nsumx = 0
self % nlayers = 0
!
end subroutine close

!==============================================================================

function constructor(file, x0, y0, t0, beta, rho, nsumx, ipol, norder, &
                     dc_bias) result(self)
character(len=*),    intent(in)  :: file
real(knd),           intent(in)  :: x0, y0
real(knd),           intent(in)  :: t0
real(knd),           intent(in)  :: beta
real(knd), optional, intent(in)  :: rho
integer,   optional, intent(in)  :: nsumx
integer,   optional, intent(in)  :: ipol
integer,   optional, intent(in)  :: norder
logical,   optional, intent(in)  :: dc_bias
type(spectral_wave_data_shape_7_impl_1) :: self
!
integer :: i, ios, err_id, m
integer(int64) :: ipos1, ipos2
integer(c_int) :: fmt, shp, amp, n, order, nid, nsteps, nstrip, nlayers_c
real(c_float)  :: d, dk, dt, grav, lscale, magic, zref_c, sig_c
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::constructor'
character(len=250) :: err_msg(5)
complex(wp) :: fval, dfval
real(wp) :: dt_tpol
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
!
call self % error % clear()
!
if (present(rho)) then
    self % rho = rho
else
    self % rho = 1025.0_wp
end if

if (t0 < 0.0_wp) then
    err_msg(1) = 'The temporal seed t0 should be zero or positive.'
    write(err_msg(2),'(a,f0.8)') 't0 = ', t0
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if
self % t0 = t0
self % x0 = x0
self % y0 = y0
self % file = file

call swd_validate_binary_convention(self % file, err_id, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    call self % error % set_id_msg(err_proc, err_id, err_msg(1:2))
    return
end if

call open_swd_file(newunit=self % unit, file=self % file, status='old', &
                   as_little_endian=.true., iostat=ios)
if (ios /= 0) then
    err_msg(1) = 'Not able to open SWD file:'
    err_msg(2) = self % file
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:2))
    return
end if

! Read common SWD header
read(self % unit, end=98, err=99) magic
read(self % unit, end=98, err=99) fmt
if (fmt /= 100) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a)') 'fmt=', fmt, ' is an unknown swd format parameter.'
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
    return
end if
self % fmt = fmt

read(self % unit, end=98, err=99) shp
self % shp = shp

read(self % unit, end=98, err=99) amp
self % amp = amp
if (amp /= 1) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a)') 'amp=', amp, &
        ' is not supported for shape 7. Only amp=1 is valid.'
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
    return
end if

read(self % unit, end=98, err=99) cprog
self % prog = trim(cprog)

read(self % unit, end=98, err=99) cdate
self % date = trim(cdate)

read(self % unit, end=98, err=99) nid
allocate(character(len=int(nid)) :: cid)
allocate(character(len=int(nid)) :: self % cid)
read(self % unit, end=98, err=99) cid
self % cid = cid(:nid)

read(self % unit, end=98, err=99) grav
self % grav = grav

read(self % unit, end=98, err=99) lscale
self % lscale = lscale

read(self % unit, end=98, err=99) nstrip
self % nstrip = nstrip

read(self % unit, end=98, err=99) nsteps
self % nsteps = nsteps

read(self % unit, end=98, err=99) dt
self % dt = dt

read(self % unit, end=98, err=99) order
self % order = order

! Read shape-7 specific header
read(self % unit, end=98, err=99) n
self % n = n

read(self % unit, end=98, err=99) dk
self % dk = dk

read(self % unit, end=98, err=99) d
self % d = d

read(self % unit, end=98, err=99) zref_c
self % zref = zref_c
if (self % zref > 0.0_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.5)') 'zref = ', self % zref
    err_msg(3) = 'zref must be <= 0 (should lie safely below all wave troughs)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if
if (self % d > 0.0_wp .and. (self % d + self % zref) <= 0.0_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.5,a,f0.5)') 'd = ', self % d, ', zref = ', self % zref
    err_msg(3) = 'Effective depth d_eff = d + zref must be > 0 for finite-depth shape 7'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if

read(self % unit, end=98, err=99) nlayers_c
self % nlayers = nlayers_c
if (self % nlayers < 2) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,i0)') 'nlayers = ', self % nlayers
    err_msg(3) = 'nlayers must be >= 2 for shape 7'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if

allocate(self % sig(self % nlayers), stat=i)
if (i /= 0) then
    err_msg(1) = 'Not able to allocate sigma array'
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:1))
    return
end if
do m = 1, self % nlayers
    read(self % unit, end=98, err=99) sig_c
    self % sig(m) = sig_c
end do

! Validate sigma-layer positions
if (abs(self % sig(1)) > 1.0e-6_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.8)') 'sig(1) = ', self % sig(1)
    err_msg(3) = 'sig(1) must be 0.0 for shape 7 (bottom layer at z = zref)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if
if (abs(self % sig(self % nlayers) - 1.0_wp) > 1.0e-6_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a,f0.8)') &
        'sig(nlayers=', self % nlayers, ') = ', self % sig(self % nlayers)
    err_msg(3) = 'sig(nlayers) must be 1.0 for shape 7 (surface layer)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if
do m = 1, self % nlayers - 1
    if (self % sig(m+1) <= self % sig(m)) then
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a,i0,a,i0)') &
            'Non-monotone sigma layers at indices ', m, ' and ', m+1
        err_msg(3) = 'sigma-layer positions must be strictly increasing for shape 7'
        call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
        return
    end if
end do

! Set norder (ignored for shape 7; kinematics use sigma-coordinates)
self % norder = self % order
if (present(norder)) then
    if (norder /= 0) self % norder = norder
end if

if (present(dc_bias)) then
    self % dc_bias = dc_bias
else
    self % dc_bias = .false.
end if

if (present(nsumx)) then
    if (nsumx < 0) then
        self % nsumx = self % n
    else
        self % nsumx = nsumx
        if (nsumx > self % n) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a,i0)') 'n in SWD file = ', self % n
            write(err_msg(3),'(a,i0)') 'Requested nsumx = ', nsumx
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
            return
        end if
    end if
else
    self % nsumx = self % n
end if

! Temporal interpolation setup
if (self % nsteps == 1) then
    dt_tpol = 1.0_wp
else
    dt_tpol = self % dt
end if
if (present(ipol)) then
    call self % tpol % construct(ischeme=ipol, delta_t=dt_tpol, ierr=i)
else
    call self % tpol % construct(ischeme=0, delta_t=dt_tpol, ierr=i)
end if
self % ipol = self % tpol % ischeme
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    err_msg(2) = 'ipol is out of bounds.'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if

self % sbeta = sin(beta * pi / 180.0_wp)
self % cbeta = cos(beta * pi / 180.0_wp)
self % tmax = self % dt * (self % nsteps - 1) - self % t0
if (self % tmax < self % dt) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    err_msg(2) = 'Constructor parameter t0 is too large.'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if

! Effective depth for below-zref shape-2 extrapolation
if (self % d > 0.0_wp) then
    ! Finite depth: use effective depth d_eff = d + zref
    self % tanhdkd_eff = tanh(self % dk * (self % d + self % zref))
else
    ! Infinite depth: extrapolation uses exp(k*z') only
    self % tanhdkd_eff = 1.0_wp  ! Indicates deep-water mode
end if

! Allocate window arrays
allocate(self % h_win(0:self % n, 4),              &
         self % c_win(0:self % n, self % nlayers, 4), &
         self % h_cur(0:self % n),                 &
         self % ht_cur(0:self % n),                &
         self % c_cur(0:self % n, self % nlayers),  &
         self % ct_cur(0:self % n, self % nlayers), stat=i)
if (i /= 0) then
    err_msg(1) = 'Not able to allocate spectral window arrays.'
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:1))
    return
end if

! Initialise to zero
self % h_win = czero_c
self % c_win = czero_c

! Read first timestep into column 2
associate(h => self % h_win, c => self % c_win)
    inquire(self % unit, pos=self % ipos0)
    read(self % unit, end=98, err=99) h(:, 2)
    do m = 1, self % nlayers
        read(self % unit, end=98, err=99) c(:, m, 2)
    end do
    ipos1 = self % ipos0
    inquire(self % unit, pos=ipos2)
    ! Size of one complete timestep block in file (storage units)
    self % size_complex = (ipos2 - ipos1) / ((1 + self % nlayers) * (self % n + 1))
    self % size_step    = (ipos2 - ipos1)
end associate

self % istp = 1
self % icur = 3
self % ipt(:,1) = [1,2,3,4]
self % ipt(:,2) = [2,3,4,1]
self % ipt(:,3) = [3,4,1,2]
self % ipt(:,4) = [4,1,2,3]

return

98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
return

99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
!
end function constructor

!==============================================================================

subroutine update_time(self, time)
! Advance the temporal window to cover the requested simulation time.
! Because the SWD file does not store time derivatives, we compute them
! numerically from the four time-step window using 2nd-order finite differences.
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd), intent(in) :: time
!
integer(int64) :: ipos
integer :: istp_max, i, j, m, i1, i2, i3, i4, istp_min, ios, imove
real(wp) :: delta, teps, dt2
complex(wp) :: fval, dfval
complex(c_float) :: cdum
complex(wp), parameter :: czero = cmplx(0.0_wp, 0.0_wp, wp)
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::update_time'
character(len=250) :: err_msg(6)
!
teps = spacing(time) * 10.0_wp
if (time > self % tmax + teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Requested time is too large!'
    write(err_msg(3),'(a,f0.5)') 'User time = ', time
    write(err_msg(4),'(a,f0.5)') 'Max user time = ', self % tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if
self % tswd = self % t0 + min(time, self % tmax - teps)
if (self % tswd < -teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'time corresponds to negative swd time!'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
else if (self % tswd < 0.0_wp) then
    self % tswd = 0.0_wp
end if

if (self % nsteps == 1) then
    istp_min = 1
    delta = 0.0_wp
    istp_max = 1
    self % icur = 1
else
    istp_min = int((self % tswd - teps) / self % dt)
    delta = self % tswd / self % dt - istp_min
    istp_max = istp_min + 3
end if

dt2 = 2.0_wp * self % dt  ! denominator for central differences

associate(h => self % h_win, c => self % c_win, ic => self % icur, ip => self % ipt)

    imove = istp_max - self % istp

    if (imove < 0 .or. imove > 4) then
        ! Jump — refill buffer from scratch
        if (istp_min == 0) then
            self % istp = 0
            ic = 2
        else
            self % istp = istp_min - 1
            ic = 1
        end if
        ipos = self % ipos0 + int(self % istp, int64) * self % size_step
        read(self % unit, pos=ipos - self % size_complex, iostat=ios) cdum
        if (ios /= 0) then
            err_msg(1) = 'User time beyond what is available from SWD file:'
            err_msg(2) = self % file
            err_msg(3) = 'The file has less content than expected.'
            write(err_msg(4),'(a,f0.5)') 'Requested swd-time = ', self % tswd
            call self % error % set_id_msg(err_proc, 1003, err_msg(1:4))
            return
        end if
    end if

    do j = 1, istp_max - self % istp
        if (self % istp == self % nsteps) then
            ! Right-padding for end of file — use spectral_interpolation pad_right.
            ! Because we have no stored derivatives, we estimate the derivative from
            ! columns i1, i2, i3 using a one-sided 2nd-order formula before padding.
            ic = ic + 1
            if (ic > 4) ic = 1
            i1 = ip(1,ic)
            i2 = ip(2,ic)
            i3 = ip(3,ic)
            i4 = ip(4,ic)
            do concurrent (i = 0 : self % nsumx)
                ! Elevation: estimate derivatives then pad_right
                call self % tpol % pad_right(                                  &
                    cmplx(h(i,i1), kind=wp), cmplx(h(i,i2), kind=wp),        &
                    cmplx(h(i,i3), kind=wp),                                  &
                    (cmplx(h(i,i2), kind=wp) - cmplx(h(i,i1), kind=wp)) / self % dt, &
                    (cmplx(h(i,i3), kind=wp) - cmplx(h(i,i1), kind=wp)) / dt2,      &
                    (cmplx(h(i,i3), kind=wp) - cmplx(h(i,i2), kind=wp)) / self % dt, &
                    fval, dfval)
                h(i,i4) = fval
            end do
            do m = 1, self % nlayers
                do concurrent (i = 0 : self % nsumx)
                    call self % tpol % pad_right(                                  &
                        cmplx(c(i,m,i1), kind=wp), cmplx(c(i,m,i2), kind=wp),   &
                        cmplx(c(i,m,i3), kind=wp),                               &
                        (cmplx(c(i,m,i2),kind=wp) - cmplx(c(i,m,i1),kind=wp)) / self % dt, &
                        (cmplx(c(i,m,i3),kind=wp) - cmplx(c(i,m,i1),kind=wp)) / dt2,      &
                        (cmplx(c(i,m,i3),kind=wp) - cmplx(c(i,m,i2),kind=wp)) / self % dt, &
                        fval, dfval)
                    c(i,m,i4) = fval
                end do
            end do
            self % istp = istp_max
        else
            ! Read next step from file
            read(self % unit, end=98, err=99) h(:, ic)
            do m = 1, self % nlayers
                read(self % unit, end=98, err=99) c(:, m, ic)
            end do
            self % istp = self % istp + 1
            ic = ic + 1
            if (ic > 4) ic = 1
        end if
    end do

    if (istp_min == 0) then
        ! Left-padding when tswd < dt_swd: estimate derivatives from cols 2,3,4
        do concurrent (i = 0 : self % nsumx)
            call self % tpol % pad_left(                                      &
                cmplx(h(i,2), kind=wp), cmplx(h(i,3), kind=wp),             &
                cmplx(h(i,4), kind=wp),                                      &
                (cmplx(h(i,3), kind=wp) - cmplx(h(i,2), kind=wp)) / self % dt, &
                (cmplx(h(i,4), kind=wp) - cmplx(h(i,2), kind=wp)) / dt2,       &
                (cmplx(h(i,4), kind=wp) - cmplx(h(i,3), kind=wp)) / self % dt, &
                fval, dfval)
            h(i,1) = fval
        end do
        do m = 1, self % nlayers
            do concurrent (i = 0 : self % nsumx)
                call self % tpol % pad_left(                                      &
                    cmplx(c(i,m,2), kind=wp), cmplx(c(i,m,3), kind=wp),         &
                    cmplx(c(i,m,4), kind=wp),                                    &
                    (cmplx(c(i,m,3),kind=wp) - cmplx(c(i,m,2),kind=wp)) / self % dt, &
                    (cmplx(c(i,m,4),kind=wp) - cmplx(c(i,m,2),kind=wp)) / dt2,       &
                    (cmplx(c(i,m,4),kind=wp) - cmplx(c(i,m,3),kind=wp)) / self % dt, &
                    fval, dfval)
                c(i,m,1) = fval
            end do
        end do
    end if

    ! Now compute central-difference time derivatives from the full 4-step window
    ! and apply the temporal interpolation scheme for the current time within the window.
    i1 = ip(1,ic); i2 = ip(2,ic); i3 = ip(3,ic); i4 = ip(4,ic)
    do concurrent (i = 0 : self % nsumx)
        ! Numerical time derivatives at the 4 window points:
        !   df_im1 = (-3*f_im1 + 4*f_i - f_ip1) / (2*dt)  [2nd-order forward]
        !   df_i   = (f_ip1 - f_im1) / (2*dt)              [2nd-order central]
        !   df_ip1 = (f_ip2 - f_i)   / (2*dt)              [2nd-order central]
        !   df_ip2 = (f_i - 4*f_ip1 + 3*f_ip2) / (2*dt)   [2nd-order backward]
        call self % tpol % scheme(delta,                                       &
            cmplx(h(i,i1), kind=wp), cmplx(h(i,i2), kind=wp),               &
            cmplx(h(i,i3), kind=wp), cmplx(h(i,i4), kind=wp),               &
            (-3.0_wp*cmplx(h(i,i1),kind=wp) + 4.0_wp*cmplx(h(i,i2),kind=wp) &
             - cmplx(h(i,i3),kind=wp)) / dt2,                                &
            (cmplx(h(i,i3),kind=wp) - cmplx(h(i,i1),kind=wp)) / dt2,       &
            (cmplx(h(i,i4),kind=wp) - cmplx(h(i,i2),kind=wp)) / dt2,       &
            (cmplx(h(i,i2),kind=wp) - 4.0_wp*cmplx(h(i,i3),kind=wp)        &
             + 3.0_wp*cmplx(h(i,i4),kind=wp)) / dt2,                        &
            self % h_cur(i), self % ht_cur(i))
    end do
    do m = 1, self % nlayers
        do concurrent (i = 0 : self % nsumx)
            call self % tpol % scheme(delta,                                   &
                cmplx(c(i,m,i1),kind=wp), cmplx(c(i,m,i2),kind=wp),         &
                cmplx(c(i,m,i3),kind=wp), cmplx(c(i,m,i4),kind=wp),         &
                (-3.0_wp*cmplx(c(i,m,i1),kind=wp) + 4.0_wp*cmplx(c(i,m,i2),kind=wp) &
                 - cmplx(c(i,m,i3),kind=wp)) / dt2,                          &
                (cmplx(c(i,m,i3),kind=wp) - cmplx(c(i,m,i1),kind=wp)) / dt2,&
                (cmplx(c(i,m,i4),kind=wp) - cmplx(c(i,m,i2),kind=wp)) / dt2,&
                (cmplx(c(i,m,i2),kind=wp) - 4.0_wp*cmplx(c(i,m,i3),kind=wp) &
                 + 3.0_wp*cmplx(c(i,m,i4),kind=wp)) / dt2,                  &
                self % c_cur(i,m), self % ct_cur(i,m))
        end do
    end do
end associate

if (.not. self % dc_bias) then
    self % h_cur(0)    = czero
    self % ht_cur(0)   = czero
    self % c_cur(0,:)  = czero
    self % ct_cur(0,:) = czero
end if

return

98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
return

99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
!
end subroutine update_time

!==============================================================================
! Helper: compute wave elevation and free-surface slope at (xswd).
! Returns eta and zeta_x.
!==============================================================================
subroutine calc_elev_and_slope(self, xswd, eta, zeta_x, zeta_t)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(wp), intent(in) :: xswd
real(wp), intent(out) :: eta, zeta_x, zeta_t
!
integer :: j
real(wp) :: kval
complex(wp) :: kappa1, Xfun
!
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
eta = self % h_cur(0) % re
zeta_x = 0.0_wp
zeta_t = self % ht_cur(0) % re
kval = 0.0_wp
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    kval = kval + self % dk
    eta = eta + real(self % h_cur(j) * Xfun)
    zeta_x = zeta_x + kval * aimag(self % h_cur(j) * Xfun)
    zeta_t = zeta_t + real(self % ht_cur(j) * Xfun)
end do
end subroutine calc_elev_and_slope

!==============================================================================
! Helper: depth factor for below-zref shape-2 extrapolation.
! Uses Rfun recursion for tanh(k*(z'+d_eff)).
! Returns Zfun and its z-derivative Zfun_z.
!==============================================================================
subroutine depth_Zfun(self, kval, z_prime, Zfun, Zfun_z, Rfun)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(wp), intent(in) :: kval, z_prime
real(wp), intent(inout) :: Rfun   ! Carries state between successive j calls
real(wp), intent(out) :: Zfun, Zfun_z
!
real(wp) :: kappa2, kappa3, Sfun, Tfun, Ufun, Vfun
real(wp), parameter :: Rfun_eps = 100.0_wp * epsilon(1.0_wp)
!
kappa2 = exp(kval * z_prime)
kappa3 = 1.0_wp / kappa2
Sfun = kappa2
Tfun = kappa3
if (self % d <= 0.0_wp) then
    ! Infinite depth: only exp(k z')
    Zfun   = Sfun
    Zfun_z = kval * Sfun
else
    Rfun = (Rfun + self % tanhdkd_eff) / (1.0_wp + self % tanhdkd_eff * Rfun)
    if (1.0_wp - Rfun < Rfun_eps) then
        Zfun   = Sfun
        Zfun_z = kval * Sfun
    else
        Ufun = (1.0_wp + Rfun) * 0.5_wp
        Vfun = 1.0_wp - Ufun
        Zfun   = Ufun * Sfun + Vfun * Tfun
        Zfun_z = kval * (Ufun * Sfun - Vfun * Tfun)
    end if
end if
end subroutine depth_Zfun

!==============================================================================
! Helper: layer interpolation of phi and u (horizontal velocity) at sigma_eval.
! Uses piecewise-linear interpolation between adjacent layers.
! Returns phi_sigma_eval, u_sigma_eval, and dphi_dsigma.
!==============================================================================
subroutine interp_sigma(self, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, &
                         c_arr, ct_arr_or_null)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(wp),    intent(in) :: xswd
real(wp),    intent(in) :: sigma_eval          ! sigma in [0,1]
real(wp),    intent(out) :: phi_val            ! phi at sigma_eval
real(wp),    intent(out) :: u_val              ! d(phi)/d(xswd) at sigma_eval
real(wp),    intent(out) :: dphi_dsigma        ! d(phi)/d(sigma) at sigma_eval
complex(wp), intent(in) :: c_arr(0:, 1:)       ! c_cur or ct_cur, shape (0:n, nlayers)
logical,     intent(in) :: ct_arr_or_null      ! .true. = ignore k factor in u_val
!
integer :: j, m1, m2, m
real(wp) :: kval, phi_m(size(self%sig)), u_m(size(self%sig))
real(wp) :: s, dsigma
complex(wp) :: kappa1, Xfun, cx
!
! Sum over spectral components for each layer
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
phi_m = 0.0_wp
u_m   = 0.0_wp
kval  = 0.0_wp
Xfun  = 1.0_wp
if (.not. ct_arr_or_null) then
    ! j=0 DC contribution
    do m = 1, self % nlayers
        phi_m(m) = phi_m(m) + c_arr(0,m) % re
    end do
end if
do j = 1, self % nsumx
    kval = kval + self % dk
    Xfun = kappa1 * Xfun
    do m = 1, self % nlayers
        cx = c_arr(j, m) * Xfun
        phi_m(m) = phi_m(m) + cx % re
        u_m(m)   = u_m(m)   + kval * cx % im
    end do
end do

! Find the bracket [m1, m2] containing sigma_eval
m1 = self % nlayers - 1
do m = 1, self % nlayers - 1
    if (sigma_eval <= self % sig(m+1)) then
        m1 = m
        exit
    end if
end do
m2 = m1 + 1

dsigma = self % sig(m2) - self % sig(m1)
s = (sigma_eval - self % sig(m1)) / dsigma

phi_val      = phi_m(m1) + s * (phi_m(m2) - phi_m(m1))
u_val        = u_m(m1)   + s * (u_m(m2)   - u_m(m1))
dphi_dsigma  = (phi_m(m2) - phi_m(m1)) / dsigma
!
end subroutine interp_sigma

!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer :: j, m
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z, sigma_eval
real(wp) :: phi_val, u_val, dphi_dsigma
complex(wp) :: kappa1, Xfun, cx
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta

if (real(z,wp) >= self % zref) then
    call calc_elev_and_slope(self, xswd, eta, zeta_x, zeta_t)
    sigma_eval = min(max((real(z,wp) - self % zref) / (eta - self % zref), 0.0_wp), 1.0_wp)
    call interp_sigma(self, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, &
                      self % c_cur, .false.)
    res = phi_val
else
    ! Below zref: shape-2 style with bottom layer
    z_prime = real(z,wp) - self % zref
    kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
    Xfun = 1.0_wp
    res = self % c_cur(0,1) % re
    Rfun = 0.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(self, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx = self % c_cur(j,1) * Xfun
        res = res + cx % re * Zfun
    end do
end if
!
end function phi

!==============================================================================

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
! stream() is not supported for shape 7. Returns zero.
res = 0.0_knd
!
end function stream

!==============================================================================

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z, sigma_eval
real(wp) :: phi_t_sig, phi_val, u_val, dphi_dsigma, H
complex(wp) :: kappa1, Xfun, cx
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta

if (real(z,wp) >= self % zref) then
    call calc_elev_and_slope(self, xswd, eta, zeta_x, zeta_t)
    H = eta - self % zref
    sigma_eval = min(max((real(z,wp) - self % zref) / H, 0.0_wp), 1.0_wp)
    ! Phi_t at fixed sigma from ct_cur
    call interp_sigma(self, xswd, sigma_eval, phi_t_sig, u_val, dphi_dsigma, &
                      self % ct_cur, .false.)
    ! phi_sigma from c_cur (for chain-rule correction at fixed physical z)
    call interp_sigma(self, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, &
                      self % c_cur, .false.)
    ! Euler time derivative at fixed z: Phi_t|sigma - sigma * zeta_t / H * phi_sigma
    res = phi_t_sig - sigma_eval * zeta_t / H * dphi_dsigma
else
    z_prime = real(z,wp) - self % zref
    kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
    Xfun = 1.0_wp
    res = self % ct_cur(0,1) % re
    Rfun = 0.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(self, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx = self % ct_cur(j,1) * Xfun
        res = res + cx % re * Zfun
    end do
end if
!
end function phi_t

!==============================================================================

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res(3)
!
integer :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z, sigma_eval
real(wp) :: phi_val, u_val, dphi_dsigma, H
complex(wp) :: kappa1, Xfun, cx
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta

if (real(z,wp) >= self % zref) then
    call calc_elev_and_slope(self, xswd, eta, zeta_x, zeta_t)
    H = eta - self % zref
    sigma_eval = min(max((real(z,wp) - self % zref) / H, 0.0_wp), 1.0_wp)
    call interp_sigma(self, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, &
                      self % c_cur, .false.)
    ! phi_x = u_interp - (sigma * zeta_x / H) * dphi_dsigma
    ! phi_z = dphi_dsigma / H
    res(1) = (u_val - sigma_eval * zeta_x / H * dphi_dsigma) * self % cbeta
    res(2) = (u_val - sigma_eval * zeta_x / H * dphi_dsigma) * self % sbeta
    res(3) = dphi_dsigma / H
else
    z_prime = real(z,wp) - self % zref
    kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
    Xfun = 1.0_wp
    Rfun = 0.0_wp
    kval = 0.0_wp
    res = 0.0_knd
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(self, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx = self % c_cur(j,1) * Xfun
        res(1) = res(1) + kval * cx % im * Zfun
        res(3) = res(3) + cx % re * Zfun_z
    end do
    res(2) = res(1) * self % sbeta
    res(1) = res(1) * self % cbeta
end if
!
end function grad_phi

!==============================================================================

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res(6)
!
res = 0.0_knd  ! Not yet implemented for sigma-coordinate shape
!
end function grad_phi_2nd

!==============================================================================

function acc_euler(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res(3)
!
! acc_euler() requires chain-rule corrections in the moving sigma-coordinate
! frame and is not yet implemented for shape 7. Returns zero.
res = 0.0_knd
!
end function acc_euler

!==============================================================================

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res(3)
!
! acc_particle() depends on acc_euler() and grad_phi_2nd(), both of which
! are not yet implemented for shape 7. Returns zero.
res = 0.0_knd
!
end function acc_particle

!==============================================================================

function elev(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res
!
integer :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
res = self % h_cur(0) % re
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    res = res + real(self % h_cur(j) * Xfun)
end do
!
end function elev

!==============================================================================

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res
!
integer :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
res = self % ht_cur(0) % re
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    res = res + real(self % ht_cur(j) * Xfun)
end do
!
end function elev_t

!==============================================================================

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res(3)
!
integer :: j
real(wp) :: xswd, elev_x_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
elev_x_swd = 0.0_wp
kval = 0.0_wp
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    kval = kval + self % dk
    elev_x_swd = elev_x_swd + kval * aimag(self % h_cur(j) * Xfun)
end do
res(1) = elev_x_swd * self % cbeta
res(2) = elev_x_swd * self % sbeta
res(3) = 0.0_knd
!
end function grad_elev

!==============================================================================

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res(3)
!
integer :: j
real(wp) :: xswd, elev_xx_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
elev_xx_swd = 0.0_wp
kval = 0.0_wp
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    kval = kval + self % dk
    elev_xx_swd = elev_xx_swd - kval * kval * real(self % h_cur(j) * Xfun)
end do
res(1) = elev_xx_swd * self % cbeta * self % cbeta
res(2) = elev_xx_swd * self % sbeta * self % cbeta
res(3) = elev_xx_swd * self % sbeta * self % sbeta
!
end function grad_elev_2nd

!==============================================================================

function pressure(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd)             :: res
!
integer :: j
real(wp) :: xswd, eta, zeta_x, zeta_t, z_prime, kval, Rfun, Zfun, Zfun_z, sigma_eval
real(wp) :: phi_val, u_val, dphi_dsigma, H
real(wp) :: phi_xswd, phi_z, phi_t_val, phi_sig
complex(wp) :: kappa1, Xfun, cx
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta

if (real(z,wp) >= self % zref) then
    call calc_elev_and_slope(self, xswd, eta, zeta_x, zeta_t)
    H = eta - self % zref
    sigma_eval = min(max((real(z,wp) - self % zref) / H, 0.0_wp), 1.0_wp)
    ! Spatial gradients from c_cur
    call interp_sigma(self, xswd, sigma_eval, phi_val, u_val, dphi_dsigma, &
                      self % c_cur, .false.)
    phi_xswd = u_val - sigma_eval * zeta_x / H * dphi_dsigma
    phi_z    = dphi_dsigma / H
    phi_sig  = dphi_dsigma   ! sigma-derivative of phi (for phi_t chain-rule)
    ! Phi_t at fixed sigma from ct_cur
    call interp_sigma(self, xswd, sigma_eval, phi_t_val, u_val, dphi_dsigma, &
                      self % ct_cur, .false.)
    ! Apply chain-rule: phi_t|z = Phi_t|sigma - sigma * zeta_t / H * phi_sigma
    phi_t_val = phi_t_val - sigma_eval * zeta_t / H * phi_sig
else
    z_prime = real(z,wp) - self % zref
    kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
    Xfun = 1.0_wp
    Rfun = 0.0_wp
    kval = 0.0_wp
    phi_xswd = 0.0_wp
    phi_z    = 0.0_wp
    phi_t_val = self % ct_cur(0,1) % re
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        call depth_Zfun(self, kval, z_prime, Zfun, Zfun_z, Rfun)
        cx = self % c_cur(j,1) * Xfun
        phi_xswd = phi_xswd + kval * cx % im * Zfun
        phi_z    = phi_z    + cx % re * Zfun_z
        phi_t_val = phi_t_val + real(self % ct_cur(j,1) * Xfun) * Zfun
    end do
end if
!
res = (-phi_t_val - 0.5_wp*(phi_xswd**2 + phi_z**2) - real(z,wp)*self % grav) * self % rho
!
end function pressure

!==============================================================================

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd),        intent(in) :: x, y, z
character(len=*), intent(in) :: csv
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::convergence'
character(len=250) :: err_msg(1)
!
err_msg(1) = 'convergence() not implemented for shape 7'
call self % error % set_id_msg(err_proc, 1004, err_msg)
!
end subroutine convergence

!==============================================================================

subroutine strip(self, tmin, tmax, file_swd)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd),        intent(in) :: tmin, tmax
character(len=*), intent(in) :: file_swd
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::strip'
character(len=250) :: err_msg(1)
!
err_msg(1) = 'strip() not implemented for shape 7'
call self % error % set_id_msg(err_proc, 1004, err_msg)
!
end subroutine strip

!==============================================================================

function bathymetry(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res
!
res = self % d
!
end function bathymetry

!==============================================================================

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd)             :: res(3)
!
res(1) = 0.0_knd
res(2) = 0.0_knd
res(3) = 1.0_knd
!
end function bathymetry_nvec

!==============================================================================

function get_int(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
integer                      :: res
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_int'
character(len=250) :: err_msg(1)
!
select case(name)
case('fmt')
    res = self % fmt
case('shp')
    res = self % shp
case('amp')
    res = self % amp
case('nid')
    res = len(self % cid)
case('nstrip')
    res = self % nstrip
case('nsteps')
    res = self % nsteps
case('order')
    res = self % order
case('norder')
    res = self % norder
case('n')
    res = self % n
case('ipol')
    res = self % ipol
case('impl')
    res = 1
case('nsumx')
    res = self % nsumx
case('nsumy')
    res = -1
case('nlayers')
    res = self % nlayers
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_int

!==============================================================================

function get_logical(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
logical                      :: res
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_logical'
character(len=250) :: err_msg(1)
!
select case(name)
case('dc_bias')
    res = self % dc_bias
case default
    res = .false.
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_logical

!==============================================================================

function get_real(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
real(knd)                    :: res
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_real'
character(len=250) :: err_msg(1)
!
select case(name)
case('t0')
    res = self % t0
case('x0')
    res = self % x0
case('y0')
    res = self % y0
case('beta')
    res = atan2(self % sbeta, self % cbeta) * 180.0_wp / pi
    if (res < 0.0_wp) res = res + 360.0_wp
case('rho')
    res = self % rho
case('magic')
    res = swd_magic_number
case('grav')
    res = self % grav
case('lscale')
    res = self % lscale
case('dt')
    res = self % dt
case('dk')
    res = self % dk
case('d')
    res = self % d
case('zref')
    res = self % zref
case('tmax')
    res = self % tmax
case('lmin')
    res = 2.0_knd * pi / (self % dk * self % nsumx)
case('lmax')
    res = 2.0_knd * pi / self % dk
case('sizex')
    res = 2.0_knd * pi / self % dk
case('sizey')
    res = 0.0_knd
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_real

!==============================================================================

function get_chr(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
character(len=:), allocatable :: res
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_chr'
character(len=250) :: err_msg(1)
!
select case(name)
case('file', 'file_swd')
    res = self % file
case('version')
    res = version
case('class')
    res = 'spectral_wave_data_shape_7_impl_1'
case('prog')
    res = self % prog
case('date')
    res = self % date
case('cid')
    res = self % cid
case default
    res = 'unknown name specified'
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_chr

!==============================================================================

end module spectral_wave_data_shape_7_impl_1_def
