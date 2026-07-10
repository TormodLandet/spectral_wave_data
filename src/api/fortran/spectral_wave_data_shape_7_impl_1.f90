module spectral_wave_data_shape_7_impl_1_def

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention, &
                             swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use spectral_interpolation_def, only: spectral_interpolation
use multilayer_long_crested_def, only: multilayer_state, multilayer_init, multilayer_close, multilayer_apply_window, &
    multilayer_phi, multilayer_stream, multilayer_phi_t, multilayer_grad_phi, multilayer_elev, multilayer_elev_t, &
    multilayer_grad_elev, multilayer_grad_elev_2nd, multilayer_pressure
use swd_version, only: version

implicit none
private

public :: spectral_wave_data_shape_7_impl_1

type, extends(spectral_wave_data) :: spectral_wave_data_shape_7_impl_1
    type(multilayer_state) :: st                   ! Shared kinematics state
    real(wp)           :: d                ! Water depth (or -1 for infinite depth)
    real(wp)           :: zref             ! z-pos of bottom sigma-layer
    integer            :: icur             ! Column index of most recently read data (1:4)
    integer            :: istp             ! Most recent SWD step in the four-step window
    integer            :: ipt(4,4)         ! Circular-buffer column mapping
    complex(c_float), allocatable :: h_win(:,:)   ! (0:n, 4)
    complex(c_float), allocatable :: c_win(:,:,:) ! (0:n, nlayers, 4)
    type(spectral_interpolation) :: tpol
contains
    procedure :: close
    procedure :: update_time
    procedure :: phi
    procedure :: stream
    procedure :: phi_t
    procedure :: grad_phi
    procedure :: grad_phi_2nd
    procedure :: acc_euler
    procedure :: acc_particle
    procedure :: elev
    procedure :: elev_t
    procedure :: grad_elev
    procedure :: grad_elev_2nd
    procedure :: pressure
    procedure :: bathymetry
    procedure :: bathymetry_nvec
    procedure :: convergence
    procedure :: strip
    procedure :: get_int
    procedure :: get_logical
    procedure :: get_real
    procedure :: get_chr
    procedure :: elev_fft           ! Surface elevation on a regular grid using FFT (not implemented)
    procedure :: grad_phi_fft       ! Grad phi on a regular grid using FFT (not implemented)
end type spectral_wave_data_shape_7_impl_1

interface spectral_wave_data_shape_7_impl_1
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_7_impl_1) :: self
logical :: opened
inquire(unit=self % unit, opened=opened)
if (opened) close(self % unit)
if (allocated(self % cid))    deallocate(self % cid)
if (allocated(self % h_win))  deallocate(self % h_win)
if (allocated(self % c_win))  deallocate(self % c_win)
call multilayer_close(self % st)
self % file = '0'
self % unit = 0
end subroutine close

!==============================================================================

function constructor(file, x0, y0, t0, beta, rho, nsumx, ipol, norder, &
                     dc_bias) result(self)
character(len=*),    intent(in)  :: file
real(knd),           intent(in)  :: x0, y0, t0, beta
real(knd), optional, intent(in)  :: rho
integer,   optional, intent(in)  :: nsumx, ipol, norder
logical,   optional, intent(in)  :: dc_bias
type(spectral_wave_data_shape_7_impl_1) :: self
!
integer :: i, ios, err_id, m
integer(int64) :: ipos1, ipos2
integer(c_int) :: fmt, shp, amp, n_c, order, nid, nsteps, nstrip, nlayers_c
real(c_float)  :: d_c, dk_c, dt_c, grav_c, lscale_c, magic, zref_c, sig_c
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::constructor'
character(len=250) :: err_msg(5)
complex(wp) :: fval, dfval
real(wp) :: dt_tpol, tanhdkd_eff
integer  :: n_f, nsumx_f, nlayers_f
real(wp) :: dk_f, d_f, zref_f
real(wp), allocatable :: sig_wp(:)
logical :: dc_bias_f
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
!
call self % error % clear()
!
if (present(rho)) then
    self % rho = rho
else
    self % rho = 1025.0_wp
end if
!
if (t0 < 0.0_wp) then
    err_msg(1) = 'The temporal seed t0 should be zero or positive.'
    write(err_msg(2),'(a,f0.8)') 't0 = ', t0
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if
self % t0 = t0; self % x0 = x0; self % y0 = y0; self % file = file
!
call swd_validate_binary_convention(self % file, err_id, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    call self % error % set_id_msg(err_proc, err_id, err_msg(1:2))
    return
end if
!
call open_swd_file(newunit=self % unit, file=self % file, status='old', &
                   as_little_endian=.true., iostat=ios)
if (ios /= 0) then
    err_msg(1) = 'Not able to open SWD file:'; err_msg(2) = self % file
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:2)); return
end if
!
read(self % unit, end=98, err=99) magic
read(self % unit, end=98, err=99) fmt
if (fmt /= 100) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a)') 'fmt=', fmt, ' is an unknown swd format parameter.'
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:2)); return
end if
self % fmt = fmt
read(self % unit, end=98, err=99) shp;       self % shp   = shp
read(self % unit, end=98, err=99) amp;       self % amp   = amp
if (amp /= 1) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a)') 'amp=', amp, &
        ' is not supported for shape 7. Only amp=1 is valid.'
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:2)); return
end if
read(self % unit, end=98, err=99) cprog;     self % prog  = trim(cprog)
read(self % unit, end=98, err=99) cdate;     self % date  = trim(cdate)
read(self % unit, end=98, err=99) nid
allocate(character(len=int(nid)) :: cid, self % cid)
read(self % unit, end=98, err=99) cid;       self % cid   = cid(:nid)
read(self % unit, end=98, err=99) grav_c;    self % grav  = grav_c
read(self % unit, end=98, err=99) lscale_c;  self % lscale = lscale_c
read(self % unit, end=98, err=99) nstrip;    self % nstrip = nstrip
read(self % unit, end=98, err=99) nsteps;    self % nsteps = nsteps
read(self % unit, end=98, err=99) dt_c;      self % dt    = dt_c
read(self % unit, end=98, err=99) order;     self % order = order
read(self % unit, end=98, err=99) n_c
read(self % unit, end=98, err=99) dk_c
read(self % unit, end=98, err=99) d_c
read(self % unit, end=98, err=99) zref_c
n_f    = int(n_c)
dk_f   = real(dk_c, wp)
d_f    = real(d_c,  wp)
zref_f = real(zref_c, wp)
self % d    = d_f
self % zref = zref_f
if (zref_f > 0.0_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.5)') 'zref = ', zref_f
    err_msg(3) = 'zref must be <= 0 (should lie safely below all wave troughs)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
end if
if (d_f > 0.0_wp .and. (d_f + zref_f) <= 0.0_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.5,a,f0.5)') 'd = ', d_f, ', zref = ', zref_f
    err_msg(3) = 'Effective depth d_eff = d + zref must be > 0 for finite-depth shape 7'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
end if
!
read(self % unit, end=98, err=99) nlayers_c
nlayers_f = int(nlayers_c)
if (nlayers_f < 2) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,i0)') 'nlayers = ', nlayers_f
    err_msg(3) = 'nlayers must be >= 2 for shape 7'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
end if
allocate(sig_wp(nlayers_f))
do m = 1, nlayers_f
    read(self % unit, end=98, err=99) sig_c; sig_wp(m) = real(sig_c, wp)
end do
if (abs(sig_wp(1)) > 1.0e-6_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,f0.8)') 'sig(1) = ', sig_wp(1)
    err_msg(3) = 'sig(1) must be 0.0 for shape 7 (bottom layer at z = zref)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
end if
if (abs(sig_wp(nlayers_f) - 1.0_wp) > 1.0e-6_wp) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a,f0.8)') 'sig(nlayers=', nlayers_f, ') = ', sig_wp(nlayers_f)
    err_msg(3) = 'sig(nlayers) must be 1.0 for shape 7 (surface layer)'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
end if
do m = 1, nlayers_f - 1
    if (sig_wp(m+1) <= sig_wp(m)) then
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a,i0,a,i0)') &
            'Non-monotone sigma layers at indices ', m, ' and ', m+1
        err_msg(3) = 'sigma-layer positions must be strictly increasing for shape 7'
        call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
    end if
end do
!
self % norder = self % order
if (present(norder)) then
    if (norder /= 0) self % norder = norder
end if
!
if (present(dc_bias)) then
    dc_bias_f = dc_bias
else
    dc_bias_f = .false.
end if
self % dc_bias = dc_bias_f
!
if (present(nsumx)) then
    if (nsumx < 0) then
        nsumx_f = n_f
    else
        if (nsumx > n_f) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a,i0)') 'n in SWD file = ', n_f
            write(err_msg(3),'(a,i0)') 'Requested nsumx = ', nsumx
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:3)); return
        end if
        nsumx_f = nsumx
    end if
else
    nsumx_f = n_f
end if
!
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
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2)); return
end if
!
self % sbeta = sin(beta*pi/180.0_wp)
self % cbeta = cos(beta*pi/180.0_wp)
self % tmax  = self % dt * (self % nsteps - 1) - self % t0
if (self % tmax < self % dt) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    err_msg(2) = 'Constructor parameter t0 is too large.'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2)); return
end if
!
if (d_f > 0.0_wp) then
    tanhdkd_eff = tanh(dk_f * (d_f + zref_f))
else
    tanhdkd_eff = 1.0_wp
end if
!
call multilayer_init(self % st, n_f, nsumx_f, dk_f, d_f, zref_f, tanhdkd_eff, &
             nlayers_f, sig_wp, self % cbeta, self % sbeta, self % x0,  &
             self % grav, self % rho, dc_bias_f)
deallocate(sig_wp)
!
allocate(self % h_win(0:self % st % n, 4),                      &
         self % c_win(0:self % st % n, self % st % nlayers, 4), &
         stat=i)
if (i /= 0) then
    err_msg(1) = 'Not able to allocate spectral window arrays.'
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:1)); return
end if
self % h_win = czero_c
self % c_win = czero_c
!
associate(h => self % h_win, c => self % c_win)
    inquire(self % unit, pos=self % ipos0)
    read(self % unit, end=98, err=99) h(:, 2)
    do m = 1, self % st % nlayers
        read(self % unit, end=98, err=99) c(:, m, 2)
    end do
    ipos1 = self % ipos0; inquire(self % unit, pos=ipos2)
    self % size_complex = (ipos2 - ipos1) / &
                          ((1 + self % st % nlayers) * (self % st % n + 1))
    self % size_step    = (ipos2 - ipos1)
end associate
!
self % istp = 1; self % icur = 3
self % ipt(:,1) = [1,2,3,4]; self % ipt(:,2) = [2,3,4,1]
self % ipt(:,3) = [3,4,1,2]; self % ipt(:,4) = [4,1,2,3]
return
98 continue
err_msg(1) = 'End of file when reading data from file:'; err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2)); return
99 continue
err_msg(1) = 'Error when reading data from file:'; err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
end function constructor

!==============================================================================

subroutine update_time(self, time)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd), intent(in) :: time
!
integer(int64) :: ipos
integer :: istp_max, i, j, m, i1, i2, i3, i4, istp_min, ios, imove
real(wp) :: delta, teps, dt2
complex(wp) :: fval, dfval
complex(c_float) :: cdum
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
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4)); return
end if
self % tswd = self % t0 + min(time, self % tmax - teps)
if (self % tswd < -teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'time corresponds to negative swd time!'
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2)); return
else if (self % tswd < 0.0_wp) then
    self % tswd = 0.0_wp
end if
!
if (self % nsteps == 1) then
    istp_min = 1; delta = 0.0_wp; istp_max = 1; self % icur = 1
else
    istp_min = int((self % tswd - teps) / self % dt)
    delta    = self % tswd / self % dt - istp_min
    istp_max = istp_min + 3
end if
dt2 = 2.0_wp * self % dt
!
associate(h => self % h_win, c => self % c_win, ic => self % icur, ip => self % ipt)
    imove = istp_max - self % istp
    if (imove < 0 .or. imove > 4) then
        if (istp_min == 0) then
            self % istp = 0; ic = 2
        else
            self % istp = istp_min - 1; ic = 1
        end if
        ipos = self % ipos0 + int(self % istp, int64) * self % size_step
        read(self % unit, pos=ipos - self % size_complex, iostat=ios) cdum
        if (ios /= 0) then
            err_msg(1) = 'User time beyond what is available from SWD file:'
            err_msg(2) = self % file
            err_msg(3) = 'The file has less content than expected.'
            write(err_msg(4),'(a,f0.5)') 'Requested swd-time = ', self % tswd
            call self % error % set_id_msg(err_proc, 1003, err_msg(1:4)); return
        end if
    end if
    do j = 1, istp_max - self % istp
        if (self % istp == self % nsteps) then
            ic = ic + 1; if (ic > 4) ic = 1
            i1 = ip(1,ic); i2 = ip(2,ic); i3 = ip(3,ic); i4 = ip(4,ic)
            do concurrent (i = 0 : self % st % nsumx)
                call self % tpol % pad_right(                                         &
                    cmplx(h(i,i1),kind=wp), cmplx(h(i,i2),kind=wp), cmplx(h(i,i3),kind=wp), &
                    (cmplx(h(i,i2),kind=wp)-cmplx(h(i,i1),kind=wp))/self%dt,         &
                    (cmplx(h(i,i3),kind=wp)-cmplx(h(i,i1),kind=wp))/dt2,             &
                    (cmplx(h(i,i3),kind=wp)-cmplx(h(i,i2),kind=wp))/self%dt,         &
                    fval, dfval)
                h(i,i4) = fval
            end do
            do m = 1, self % st % nlayers
                do concurrent (i = 0 : self % st % nsumx)
                    call self % tpol % pad_right(                                              &
                        cmplx(c(i,m,i1),kind=wp), cmplx(c(i,m,i2),kind=wp), cmplx(c(i,m,i3),kind=wp), &
                        (cmplx(c(i,m,i2),kind=wp)-cmplx(c(i,m,i1),kind=wp))/self%dt,         &
                        (cmplx(c(i,m,i3),kind=wp)-cmplx(c(i,m,i1),kind=wp))/dt2,             &
                        (cmplx(c(i,m,i3),kind=wp)-cmplx(c(i,m,i2),kind=wp))/self%dt,         &
                        fval, dfval)
                    c(i,m,i4) = fval
                end do
            end do
            self % istp = istp_max
        else
            read(self % unit, end=98, err=99) h(:, ic)
            do m = 1, self % st % nlayers
                read(self % unit, end=98, err=99) c(:, m, ic)
            end do
            self % istp = self % istp + 1; ic = ic + 1; if (ic > 4) ic = 1
        end if
    end do
    if (istp_min == 0) then
        do concurrent (i = 0 : self % st % nsumx)
            call self % tpol % pad_left(                                           &
                cmplx(h(i,2),kind=wp), cmplx(h(i,3),kind=wp), cmplx(h(i,4),kind=wp), &
                (cmplx(h(i,3),kind=wp)-cmplx(h(i,2),kind=wp))/self%dt,            &
                (cmplx(h(i,4),kind=wp)-cmplx(h(i,2),kind=wp))/dt2,               &
                (cmplx(h(i,4),kind=wp)-cmplx(h(i,3),kind=wp))/self%dt,            &
                fval, dfval)
            h(i,1) = fval
        end do
        do m = 1, self % st % nlayers
            do concurrent (i = 0 : self % st % nsumx)
                call self % tpol % pad_left(                                                &
                    cmplx(c(i,m,2),kind=wp), cmplx(c(i,m,3),kind=wp), cmplx(c(i,m,4),kind=wp), &
                    (cmplx(c(i,m,3),kind=wp)-cmplx(c(i,m,2),kind=wp))/self%dt,            &
                    (cmplx(c(i,m,4),kind=wp)-cmplx(c(i,m,2),kind=wp))/dt2,                &
                    (cmplx(c(i,m,4),kind=wp)-cmplx(c(i,m,3),kind=wp))/self%dt,            &
                    fval, dfval)
                c(i,m,1) = fval
            end do
        end do
    end if
    i1 = ip(1,ic); i2 = ip(2,ic); i3 = ip(3,ic); i4 = ip(4,ic)
    call multilayer_apply_window(self % st, h, c, self % tpol, i1, i2, i3, i4, &
                         delta, self % dt)
end associate
return
98 continue
err_msg(1) = 'End of file when reading data from file:'; err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2)); return
99 continue
err_msg(1) = 'Error when reading data from file:'; err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
end subroutine update_time

!==============================================================================
! Kinematic TBPs — all delegate to ml_* from multilayer_long_crested
!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = multilayer_phi(self % st, x, y, z)
end function phi

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = multilayer_stream(self % st, x, y, z)
end function stream

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = multilayer_phi_t(self % st, x, y, z)
end function phi_t

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = multilayer_grad_phi(self % st, x, y, z)
end function grad_phi

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(6)
res = 0.0_knd
end function grad_phi_2nd

function acc_euler(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = 0.0_knd
end function acc_euler

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = 0.0_knd
end function acc_particle

function elev(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
res = multilayer_elev(self % st, x, y)
end function elev

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
res = multilayer_elev_t(self % st, x, y)
end function elev_t

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res = multilayer_grad_elev(self % st, x, y)
end function grad_elev

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res = multilayer_grad_elev_2nd(self % st, x, y)
end function grad_elev_2nd

function pressure(self, x, y, z) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = multilayer_pressure(self % st, x, y, z)
end function pressure

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd), intent(in) :: x, y, z
character(len=*), intent(in) :: csv
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::convergence'
character(len=250) :: err_msg(1)
err_msg(1) = 'convergence() not implemented for shape 7'
call self % error % set_id_msg(err_proc, 1004, err_msg)
end subroutine convergence

subroutine strip(self, tmin, tmax, file_swd)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
real(knd), intent(in) :: tmin, tmax
character(len=*), intent(in) :: file_swd
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::strip'
character(len=250) :: err_msg(1)
err_msg(1) = 'strip() not implemented for shape 7'
call self % error % set_id_msg(err_proc, 1004, err_msg)
end subroutine strip

function bathymetry(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
res = self % d
end function bathymetry

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res(1) = 0.0_knd; res(2) = 0.0_knd; res(3) = 1.0_knd
end function bathymetry_nvec

function get_int(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
integer :: res
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_int'
character(len=250) :: err_msg(1)
select case(name)
case('fmt');     res = self % fmt
case('shp');     res = self % shp
case('amp');     res = self % amp
case('nid');     res = len(self % cid)
case('nstrip');  res = self % nstrip
case('nsteps');  res = self % nsteps
case('order');   res = self % order
case('norder');  res = self % norder
case('n');       res = self % st % n
case('ipol');    res = self % ipol
case('impl');    res = 1
case('nsumx');   res = self % st % nsumx
case('nsumy');   res = -1
case('nlayers'); res = self % st % nlayers
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
end function get_int

function get_logical(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
logical :: res
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_logical'
character(len=250) :: err_msg(1)
select case(name)
case('dc_bias'); res = self % dc_bias
case default
    res = .false.
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
end function get_logical

function get_real(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
real(knd) :: res
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_real'
character(len=250) :: err_msg(1)
select case(name)
case('t0');     res = self % t0
case('x0');     res = self % x0
case('y0');     res = self % y0
case('beta')
    res = atan2(self % sbeta, self % cbeta) * 180.0_wp / pi
    if (res < 0.0_wp) res = res + 360.0_wp
case('rho');    res = self % rho
case('magic');  res = swd_magic_number
case('grav');   res = self % grav
case('lscale'); res = self % lscale
case('dt');     res = self % dt
case('dk');     res = self % st % dk
case('d');      res = self % d
case('zref');   res = self % zref
case('tmax');   res = self % tmax
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
end function get_real

function get_chr(self, name) result(res)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self
character(len=*), intent(in) :: name
character(len=:), allocatable :: res
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::get_chr'
character(len=250) :: err_msg(1)
select case(name)
case('prog');    res = self % prog
case('date');    res = self % date
case('cid');     res = self % cid
case('file');    res = self % file
case('version'); res = version
case default
    res = 'Unknown parameter'
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
end function get_chr

!==============================================================================

function elev_fft(self, nx_fft_in, ny_fft_in) result(elev)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self ! Actual class
integer, optional, intent(in) :: nx_fft_in, ny_fft_in
real(knd), allocatable :: elev(:, :)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::elev_fft'
character(len=:), allocatable :: err_msg(:)

allocate(elev(1,1))
elev = huge(elev)

err_msg = ["not implemented"]
call self % error % set_id_msg(err_proc, 1004, err_msg)

end function elev_fft

!==============================================================================

function grad_phi_fft(self, z, nx_fft_in, ny_fft_in) result(grad_phi)
class(spectral_wave_data_shape_7_impl_1), intent(inout) :: self ! Actual class
real(wp), intent(in) :: z
integer, optional, intent(in) :: nx_fft_in, ny_fft_in
real(knd), allocatable :: grad_phi(:, :, :)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_7_impl_1::grad_phi_fft'
character(len=:), allocatable :: err_msg(:)

allocate(grad_phi(1,1,1))
grad_phi = huge(grad_phi)

err_msg = ["not implemented"]
call self % error % set_id_msg(err_proc, 1004, err_msg)

end function grad_phi_fft

!==============================================================================

end module spectral_wave_data_shape_7_impl_1_def
