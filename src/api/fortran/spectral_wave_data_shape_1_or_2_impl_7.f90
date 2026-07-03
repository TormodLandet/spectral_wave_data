module spectral_wave_data_shape_1_or_2_impl_7_def

! Lazy H2-operator implementation for SWD shape 1 and 2 with amp=2.
!
! When a shape 1 (deep water) or shape 2 (finite depth) SWD file is written
! with amp=2, the c-array contains the velocity potential evaluated on the
! free surface, not evaluated at z=0.
!
! This implementation runs the HOSM H2-operator on the fly to transfer the
! surface potential down to a fixed sigma-coordinate grid, then reuses the
! shape-7 sigma-coordinate kinematics to answer queries.
!
! Constraints:
!   - nx must be a power of two  (nx = 2*n where n is SWD's spectral count)
!   - long-crested (1-D) only
!
! Reference: docs/source/amp2_free_surface_potential.rst

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float, c_double

use kind_values,            only: knd => kind_swd_interface, wp => kind_swd_internal
use open_swd_file_def,      only: open_swd_file, swd_validate_binary_convention, &
                                  swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use spectral_interpolation_def, only: spectral_interpolation
use multilayer_long_crested_def, only: ml_state, ml_init, ml_close, ml_apply_window, &
    ml_phi, ml_stream, ml_phi_t, ml_grad_phi, ml_elev, ml_elev_t, &
    ml_grad_elev, ml_grad_elev_2nd, ml_pressure
use swd_fft_def,            only: swd_fft_plan, fft_init, fft_destroy, &
                                  fft_swd_to_real
use hosm_h2_operator_def,   only: h2op_state, h2op_init, h2op_close, &
                                  h2op_calc_potential, h2op_convert_to_swd_layers
use swd_version,            only: version

implicit none
private

public :: spectral_wave_data_shape_1_or_2_impl_7

!==============================================================================
! Number of H2 sigma layers: read from environment, default 20.
!==============================================================================
integer, parameter :: nlayers_default = 20
integer, parameter :: nlayers_min     = 2
integer, parameter :: nlayers_max     = 100

! H2 nonlinearity order and step count (hardcoded for now)
integer, parameter :: M_kin_h2  = 5
integer, parameter :: nstep_h2  = 20

!==============================================================================
! The new shape class
!==============================================================================

type, extends(spectral_wave_data) :: spectral_wave_data_shape_1_or_2_impl_7
    !-- SWD file metadata
    integer  :: n        ! highest spectral index in file (n_swd = n)
    integer  :: nx       ! physical grid size = 2*n (power of two)
    real(wp) :: dk_val   ! wavenumber spacing (stored as wp for interface compat)
    real(wp) :: d        ! water depth (> 0 for shape 2, -1 for shape 1)
    real(wp) :: zref     ! fixed reference z-position of bottom sigma layer
    integer  :: nsumx    ! number of spectral components used in summation
    !-- temporal window (same 4-column scheme as shape 7)
    integer  :: icur              ! column of most recently generated data (1:4)
    integer  :: istp              ! most recent H2-computed SWD step in window
    integer  :: ipt(4,4)          ! circular-buffer column mapping
    complex(c_float), allocatable :: h_win(:,:)    ! (0:n, 4) elevation window
    complex(c_float), allocatable :: c_win(:,:,:)  ! (0:n, nlayers, 4) potential window
    integer  :: nlayers           ! number of sigma layers used
    type(spectral_interpolation) :: tpol
    !-- shared kinematics state (filled at each update_time call)
    type(ml_state) :: st
    !-- H2 operator
    type(h2op_state) :: h2op
    !-- FFT plan for n_swd<->nx conversion (plan_nx lives also in h2op, but
    !   we need it before h2op is ready, so keep a reference copy here)
    type(swd_fft_plan) :: plan_nx_prepass
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
end type spectral_wave_data_shape_1_or_2_impl_7

interface spectral_wave_data_shape_1_or_2_impl_7
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================
! Constructor
!==============================================================================

function constructor(file, x0, y0, t0, beta, rho, nsumx, ipol, norder, &
                     dc_bias) result(self)
character(len=*),    intent(in)  :: file
real(knd),           intent(in)  :: x0, y0, t0, beta
real(knd), optional, intent(in)  :: rho
integer,   optional, intent(in)  :: nsumx
integer,   optional, intent(in)  :: ipol
integer,   optional, intent(in)  :: norder
logical,   optional, intent(in)  :: dc_bias
type(spectral_wave_data_shape_1_or_2_impl_7) :: self

integer :: i, ios, err_id, nlayers_h2
integer(int64) :: ipos1, ipos2
integer(c_int) :: fmt, shp, amp, n, order, nid, nsteps, nstrip
real(c_float)  :: d_c, dk_c, dt_c, grav_c, lscale_c, magic_c
real(c_double) :: dk64, zref64, d64
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
character(len=*), parameter :: err_proc = &
    'spectral_wave_data_shape_1_or_2_impl_7::constructor'
character(len=250) :: err_msg(6)
real(wp) :: dt_tpol, eta_min_global
character(len=20) :: env_val
integer  :: env_status
complex(c_float), allocatable :: h_step(:)
real(c_double), allocatable :: eta_nx(:)
real(wp), allocatable :: sig(:)
real(wp) :: tanhdkd_eff_wp

call self%error%clear()

if (present(rho)) then
    self%rho = rho
else
    self%rho = 1025.0_wp
end if
if (t0 < 0.0_wp) then
    err_msg(1) = 'The temporal seed t0 should be zero or positive.'
    write(err_msg(2),'(a,f0.8)') 't0 = ', t0
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if
self%t0 = t0
self%x0 = x0
self%y0 = y0
self%file = file

call swd_validate_binary_convention(self%file, err_id, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self%file)
    call self%error%set_id_msg(err_proc, err_id, err_msg(1:2))
    return
end if

call open_swd_file(newunit=self%unit, file=self%file, status='old', &
                   as_little_endian=.true., iostat=ios)
if (ios /= 0) then
    err_msg(1) = 'Not able to open SWD file:'
    err_msg(2) = self%file
    call self%error%set_id_msg(err_proc, 1001, err_msg(1:2))
    return
end if

! --- Read header ---
read(self%unit, end=98, err=99) magic_c
read(self%unit, end=98, err=99) fmt
if (fmt /= 100) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self%file)
    write(err_msg(2),'(a,i0)') 'Unknown fmt=', fmt
    call self%error%set_id_msg(err_proc, 1003, err_msg(1:2))
    return
end if
self%fmt = fmt
read(self%unit, end=98, err=99) shp
self%shp = shp
read(self%unit, end=98, err=99) amp
self%amp = amp
read(self%unit, end=98, err=99) cprog;  self%prog = trim(cprog)
read(self%unit, end=98, err=99) cdate;  self%date = trim(cdate)
read(self%unit, end=98, err=99) nid
allocate(character(len=int(nid)) :: cid)
allocate(character(len=int(nid)) :: self%cid)
read(self%unit, end=98, err=99) cid;  self%cid = cid(:nid)
read(self%unit, end=98, err=99) grav_c;   self%grav   = grav_c
read(self%unit, end=98, err=99) lscale_c; self%lscale = lscale_c
read(self%unit, end=98, err=99) nstrip;   self%nstrip = nstrip
read(self%unit, end=98, err=99) nsteps;   self%nsteps = nsteps
read(self%unit, end=98, err=99) dt_c;     self%dt     = dt_c
read(self%unit, end=98, err=99) order;    self%order  = order
read(self%unit, end=98, err=99) n;        self%n      = n
read(self%unit, end=98, err=99) dk_c

! Shape 2 has an extra depth parameter
if (shp == 2) then
    read(self%unit, end=98, err=99) d_c
    self%d = d_c
else
    self%d = -1.0_wp   ! infinite depth (shape 1)
end if

! --- Validate ---
if (iand(int(2*n), int(2*n)-1) /= 0) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self%file)
    write(err_msg(2),'(a,i0,a)') 'nx = 2*n = ', 2*n, ' is not a power of two.'
    err_msg(3) = 'amp=2 with lazy H2 requires nx to be a power of two.'
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if

self%nx    = 2 * int(n)
dk64       = real(dk_c, c_double)
self%dk_val = real(dk64, wp)

! norder (for API compatibility, unused in sigma-coordinate evaluation)
self%norder = self%order
if (present(norder)) then
    if (norder /= 0) self%norder = norder
end if
if (present(dc_bias)) then
    self%dc_bias = dc_bias
else
    self%dc_bias = .false.
end if
if (present(nsumx)) then
    if (nsumx < 0) then
        self%nsumx = int(n)
    else
        self%nsumx = min(nsumx, int(n))
    end if
else
    self%nsumx = int(n)
end if

! Number of H2 sigma layers (from environment or default)
call get_environment_variable('SWD_NUM_H2_LAYERS', env_val, status=env_status)
if (env_status == 0) then
    read(env_val, *, iostat=ios) nlayers_h2
    if (ios /= 0 .or. nlayers_h2 < nlayers_min .or. nlayers_h2 > nlayers_max) then
        nlayers_h2 = nlayers_default
    end if
else
    nlayers_h2 = nlayers_default
end if
self%nlayers = nlayers_h2

! --- Temporal interpolation ---
dt_tpol = self%dt
if (self%nsteps == 1) dt_tpol = 1.0_wp
if (present(ipol)) then
    call self%tpol%construct(ischeme=ipol, delta_t=dt_tpol, ierr=i)
else
    call self%tpol%construct(ischeme=0, delta_t=dt_tpol, ierr=i)
end if
self%ipol = self%tpol%ischeme
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self%file)
    err_msg(2) = 'ipol is out of bounds.'
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if

self%sbeta = sin(beta*pi/180.0_wp)
self%cbeta = cos(beta*pi/180.0_wp)
self%tmax  = self%dt * (self%nsteps - 1) - self%t0
if (self%tmax < self%dt) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self%file)
    err_msg(2) = 'Constructor parameter t0 is too large.'
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if

! --- FFT plan for nx (used in eta-prepass and swd-coefficient conversion) ---
call fft_init(self%plan_nx_prepass, self%nx, err_msg(1))
if (err_msg(1) /= '') then
    call self%error%set_id_msg(err_proc, 1005, err_msg(1:1))
    return
end if

! --- Eta-prepass: scan all timesteps to find global eta_min for fixed zref ---
inquire(self%unit, pos=self%ipos0)

! Allocate scratch and measure step size from the first timestep
allocate(h_step(0:int(n)))
allocate(eta_nx(0:self%nx-1))
ipos1 = self%ipos0
read(self%unit, pos=ipos1, end=98, err=99) h_step(:)         ! h(0:n)
read(self%unit, end=98, err=99) h_step(:)                    ! ht(0:n)
read(self%unit, end=98, err=99) h_step(:)                    ! c(0:n)
read(self%unit, end=98, err=99) h_step(:)                    ! ct(0:n)
inquire(self%unit, pos=ipos2)
self%size_complex = (ipos2 - ipos1) / (4 * (int(n) + 1))
self%size_step    = ipos2 - ipos1

! Now scan all timesteps for eta_min
eta_min_global = 0.0_wp

do i = 0, int(nsteps) - 1
    ! Seek to h(0:n) of step i
    read(self%unit, pos = self%ipos0 + int(i, int64)*self%size_step, &
         end=98, err=99) h_step(:)
    ! Convert to real-space eta on nx grid
    call fft_swd_to_real(self%plan_nx_prepass, &
                         cmplx(h_step(:), kind=c_double), &
                         int(n), self%nx, eta_nx, err_msg(1))
    if (err_msg(1) /= '') then
        call self%error%set_id_msg(err_proc, 1005, err_msg(1:1))
        return
    end if
    if (i == 0 .or. minval(eta_nx) < eta_min_global) then
        eta_min_global = minval(eta_nx)
    end if
end do

! Compute zref: as 1.2 * eta_min
! Note: the Python code used for reference rounds the zref value for nicer debug 
! prints, but any floating point value is equally good for us here
d64   = real(self%d, c_double)
zref64 = 1.2_c_double * eta_min_global

! Clamp: zref must satisfy d + zref > 0 for finite depth
if (d64 > 0.0_c_double) then
    if (zref64 <= -d64) then
        zref64 = -d64 + 1.0e-4_c_double * d64   ! seabed clamp: just above -d
    end if
end if
! zref must be <= 0 (shape-7 convention)
if (zref64 > 0.0_c_double) zref64 = 0.0_c_double

self%zref = real(zref64, wp)

deallocate(eta_nx)

! --- Initialise H2 operator ---
call h2op_init(self%h2op, M_kin=M_kin_h2, nstep=nstep_h2, nx=self%nx, &
               dk=real(dk64, c_double), h_depth=d64, zref=real(zref64, c_double), &
               err_msg=err_msg(1))
if (err_msg(1) /= '') then
    call self%error%set_id_msg(err_proc, 1005, err_msg(1:1))
    return
end if

! --- Allocate 4-column window arrays (same layout as shape 7) ---
allocate(self%h_win(0:self%n, 4), stat=i)
allocate(self%c_win(0:self%n, self%nlayers, 4), stat=ios)
i = i + ios
if (i /= 0) then
    err_msg(1) = 'Not able to allocate spectral window arrays.'
    call self%error%set_id_msg(err_proc, 1005, err_msg(1:1))
    return
end if
self%h_win = cmplx(0.0_c_float, 0.0_c_float, c_float)
self%c_win = cmplx(0.0_c_float, 0.0_c_float, c_float)

! --- Build sigma layer positions (uniform 0..1) ---
allocate(sig(nlayers_h2))
do i = 1, nlayers_h2
    sig(i) = real(i - 1, wp) / real(nlayers_h2 - 1, wp)
end do
sig(1)          = 0.0_wp
sig(nlayers_h2) = 1.0_wp

! tanhdkd_eff for ml_state (shape-2 below-zref extrapolation)
if (self%d > 0.0_wp) then
    tanhdkd_eff_wp = tanh(self%dk_val * (self%d + self%zref))
else
    tanhdkd_eff_wp = 1.0_wp
end if

! --- Initialise the shared ml_state ---
call ml_init(self%st, &
    n           = self%n,          &
    nsumx       = self%nsumx,      &
    dk          = self%dk_val,     &
    d           = self%d,          &
    zref        = self%zref,       &
    tanhdkd_eff = tanhdkd_eff_wp,  &
    nlayers     = nlayers_h2,      &
    sig         = sig,             &
    cbeta       = self%cbeta,      &
    sbeta       = self%sbeta,      &
    x0          = self%x0,         &
    grav        = self%grav,       &
    rho         = self%rho,        &
    dc_bias     = self%dc_bias)

deallocate(sig, h_step)

! --- Circular buffer initialisation (same as shape 7) ---
self%istp = 0
self%icur = 3
self%ipt(:,1) = [1,2,3,4]
self%ipt(:,2) = [2,3,4,1]
self%ipt(:,3) = [3,4,1,2]
self%ipt(:,4) = [4,1,2,3]

return
98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self%file
call self%error%set_id_msg(err_proc, 1003, err_msg(1:2))
return
99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self%file
call self%error%set_id_msg(err_proc, 1003, err_msg(1:2))

end function constructor

!==============================================================================
! close
!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_1_or_2_impl_7) :: self
logical :: opened
inquire(unit=self%unit, opened=opened)
if (opened) close(self%unit)
call ml_close(self%st)
call h2op_close(self%h2op)
call fft_destroy(self%plan_nx_prepass)
if (allocated(self%cid))   deallocate(self%cid)
if (allocated(self%h_win)) deallocate(self%h_win)
if (allocated(self%c_win)) deallocate(self%c_win)
self%file = '0'
self%unit = 0
self%n    = 0
end subroutine close

!==============================================================================
! update_time — advance the H2 temporal window to cover `time`
!==============================================================================

subroutine update_time(self, time)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
real(knd), intent(in) :: time

integer(int64) :: ipos
integer  :: istp_max, istp_min, j, imove, ios
real(wp) :: delta, teps, dt2
complex(wp) :: fval, dfval
character(len=*), parameter :: err_proc = &
    'spectral_wave_data_shape_1_or_2_impl_7::update_time'
character(len=250) :: err_msg(6)
character(len=250) :: h2_err
complex(c_float), allocatable :: h_file(:), c_file(:)  ! scratch for file reads

teps = spacing(time) * 10.0_wp
if (time > self%tmax + teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self%file)
    err_msg(2) = 'Requested time is too large!'
    write(err_msg(3),'(a,f0.5)') 'User time = ', time
    write(err_msg(4),'(a,f0.5)') 'Max user time = ', self%tmax
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if
self%tswd = self%t0 + min(time, self%tmax - teps)
if (self%tswd < -teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self%file)
    err_msg(2) = 'time corresponds to negative swd time!'
    call self%error%set_id_msg(err_proc, 1004, err_msg(1:2))
    return
else if (self%tswd < 0.0_wp) then
    self%tswd = 0.0_wp
end if

if (self%nsteps == 1) then
    istp_min = 1
    delta    = 0.0_wp
    istp_max = 1
    self%icur = 1
else
    istp_min = int((self%tswd - teps) / self%dt)
    delta    = self%tswd / self%dt - istp_min
    istp_max = istp_min + 3
end if

dt2 = 2.0_wp * self%dt

allocate(h_file(0:self%n), c_file(0:self%n))

associate(h => self%h_win, c => self%c_win, ic => self%icur, ip => self%ipt)
    imove = istp_max - self%istp
    if (imove < 0 .or. imove > 4) then
        ! Jump or rewind: refill buffer from scratch
        if (istp_min == 0) then
            self%istp = 0
            ic = 2
        else
            self%istp = istp_min - 1
            ic = 1
        end if
    end if

    do j = 1, istp_max - self%istp
        if (self%istp >= self%nsteps) then
            ! Right-padding at end of time domain
            ic = ic + 1
            if (ic > 4) ic = 1
            call do_pad_right(self, h, c, ic, ip, dt2, fval, dfval)
            self%istp = istp_max
        else
            ! Generate H2 for SWD step self%istp+1 into CURRENT ic, then advance
            ! (mirrors shape-7/shape-1 pattern: read-into-current-ic, then ic+1)
            call generate_h2_column(self, self%istp + 1, ic, h_file, c_file, &
                                    h2_err)
            if (h2_err /= '') then
                err_msg(1) = h2_err
                call self%error%set_id_msg(err_proc, 1005, err_msg(1:1))
                deallocate(h_file, c_file)
                return
            end if
            self%istp = self%istp + 1
            ic = ic + 1
            if (ic > 4) ic = 1
        end if
    end do

    if (istp_min == 0) then
        ! Left-padding when tswd < dt_swd
        call do_pad_left(self, h, c, ip(:,ic), dt2, fval, dfval)
    end if

    ! Apply the 4-step temporal window to fill ml_state h_cur/ht_cur/c_cur/ct_cur
    call ml_apply_window(self%st, &
        h_win=h, c_win=c, tpol=self%tpol, &
        i1=ip(1,ic), i2=ip(2,ic), i3=ip(3,ic), i4=ip(4,ic), &
        delta=delta, dt=real(self%dt, wp))
end associate

deallocate(h_file, c_file)

end subroutine update_time

!==============================================================================
! generate_h2_column — internal: run H2 for one SWD step and store in window
!==============================================================================

subroutine generate_h2_column(self, istp_target, ic_col, h_file, c_file, err_msg)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
integer,          intent(in)    :: istp_target, ic_col
complex(c_float), intent(inout) :: h_file(0:), c_file(0:)   ! scratch
character(len=*), intent(out)   :: err_msg

real(c_double), allocatable :: eta_nx(:), psi_nx(:)
complex(c_double), allocatable :: h_swd_layer(:), c_swd_layer(:,:)
real(c_double), allocatable :: sigma_layer(:)
integer  :: ios
integer(int64) :: ipos_step
character(len=*), parameter :: err_proc = 'generate_h2_column'

err_msg = ''

! Seek to the h-array of the target step (0-based step index = istp_target-1)
ipos_step = self%ipos0 + int(istp_target - 1, int64) * self%size_step

! Read h(0:n) and c(0:n) from file (skip ht and ct)
read(self%unit, pos=ipos_step, end=98, err=99) h_file(:)   ! h(0:n)
! skip ht: advance by one h-array worth of bytes
ipos_step = ipos_step + int(self%size_complex, int64) * int(self%n + 1, int64)
read(self%unit, pos=ipos_step + int(self%size_complex, int64) * int(self%n + 1, int64), &
     end=98, err=99) c_file(:)   ! c(0:n)   [skipping ht]

! Convert h-spectral to eta on nx grid
allocate(eta_nx(0:self%nx-1), psi_nx(0:self%nx-1))
call fft_swd_to_real(self%h2op%plan_nx, &
                     cmplx(h_file(:), kind=c_double), &
                     self%n, self%nx, eta_nx, err_msg)
if (err_msg /= '') return

! Convert c-spectral to psi on nx grid
call fft_swd_to_real(self%h2op%plan_nx, &
                     cmplx(c_file(:), kind=c_double), &
                     self%n, self%nx, psi_nx, err_msg)
if (err_msg /= '') return

! Run the H2 operator
call h2op_calc_potential(self%h2op, eta_nx, psi_nx, err_msg)
if (err_msg /= '') return

! Convert H2 results to SWD shape-7 spectral coefficients
allocate(h_swd_layer(0:self%n))
allocate(c_swd_layer(0:self%n, 1:self%nlayers))
allocate(sigma_layer(1:self%nlayers))

call h2op_convert_to_swd_layers(self%h2op, eta_nx, psi_nx, &
    self%nlayers, h_swd_layer, c_swd_layer, sigma_layer, err_msg)
if (err_msg /= '') return

! Store in window (downcast to c_float for storage)
self%h_win(:, ic_col) = cmplx(h_swd_layer(:), kind=c_float)
self%c_win(:, :, ic_col) = cmplx(c_swd_layer(:,:), kind=c_float)

deallocate(eta_nx, psi_nx, h_swd_layer, c_swd_layer, sigma_layer)
return
98 continue
err_msg = 'End of file when reading data for H2 generation from: ' // trim(self%file)
return
99 continue
err_msg = 'Error when reading data for H2 generation from: ' // trim(self%file)
end subroutine generate_h2_column

!==============================================================================
! Padding helpers (same finite-difference logic as shape 7)
!==============================================================================

subroutine do_pad_right(self, h, c, ic, ip, dt2, fval, dfval)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
complex(c_float), intent(inout) :: h(0:, 1:), c(0:, 1:, 1:)
integer, intent(in) :: ic, ip(4,4)
real(wp), intent(in) :: dt2
complex(wp), intent(inout) :: fval, dfval
integer :: i, m, i1, i2, i3, i4
i1 = ip(1,ic); i2 = ip(2,ic); i3 = ip(3,ic); i4 = ip(4,ic)
do concurrent (i = 0 : self%nsumx)
    call self%tpol%pad_right( &
        cmplx(h(i,i1),kind=wp), cmplx(h(i,i2),kind=wp), cmplx(h(i,i3),kind=wp), &
        (cmplx(h(i,i2),kind=wp)-cmplx(h(i,i1),kind=wp))/self%dt, &
        (cmplx(h(i,i3),kind=wp)-cmplx(h(i,i1),kind=wp))/dt2, &
        (cmplx(h(i,i3),kind=wp)-cmplx(h(i,i2),kind=wp))/self%dt, &
        fval, dfval)
    h(i, i4) = fval
end do
do m = 1, self%nlayers
    do concurrent (i = 0 : self%nsumx)
        call self%tpol%pad_right( &
            cmplx(c(i,m,i1),kind=wp), cmplx(c(i,m,i2),kind=wp), cmplx(c(i,m,i3),kind=wp), &
            (cmplx(c(i,m,i2),kind=wp)-cmplx(c(i,m,i1),kind=wp))/self%dt, &
            (cmplx(c(i,m,i3),kind=wp)-cmplx(c(i,m,i1),kind=wp))/dt2, &
            (cmplx(c(i,m,i3),kind=wp)-cmplx(c(i,m,i2),kind=wp))/self%dt, &
            fval, dfval)
        c(i, m, i4) = fval
    end do
end do
end subroutine do_pad_right

!==============================================================================

subroutine do_pad_left(self, h, c, icols, dt2, fval, dfval)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
complex(c_float), intent(inout) :: h(0:, 1:), c(0:, 1:, 1:)
integer, intent(in) :: icols(4)  ! [i1,i2,i3,i4] = ipt(:,ic)
real(wp), intent(in) :: dt2
complex(wp), intent(inout) :: fval, dfval
integer :: i, m
integer :: i1, i2, i3, i4
i1=icols(1); i2=icols(2); i3=icols(3); i4=icols(4)
do concurrent (i = 0 : self%nsumx)
    call self%tpol%pad_left( &
        cmplx(h(i,i2),kind=wp), cmplx(h(i,i3),kind=wp), cmplx(h(i,i4),kind=wp), &
        (cmplx(h(i,i3),kind=wp)-cmplx(h(i,i2),kind=wp))/self%dt, &
        (cmplx(h(i,i4),kind=wp)-cmplx(h(i,i2),kind=wp))/dt2, &
        (cmplx(h(i,i4),kind=wp)-cmplx(h(i,i3),kind=wp))/self%dt, &
        fval, dfval)
    h(i,i1) = fval
end do
do m = 1, self%nlayers
    do concurrent (i = 0 : self%nsumx)
        call self%tpol%pad_left( &
            cmplx(c(i,m,i2),kind=wp), cmplx(c(i,m,i3),kind=wp), cmplx(c(i,m,i4),kind=wp), &
            (cmplx(c(i,m,i3),kind=wp)-cmplx(c(i,m,i2),kind=wp))/self%dt, &
            (cmplx(c(i,m,i4),kind=wp)-cmplx(c(i,m,i2),kind=wp))/dt2, &
            (cmplx(c(i,m,i4),kind=wp)-cmplx(c(i,m,i3),kind=wp))/self%dt, &
            fval, dfval)
        c(i, m, i1) = fval
    end do
end do
end subroutine do_pad_left

!==============================================================================
! Kinematics — all delegate to ml_* (same as shape 7)
!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = ml_phi(self%st, x, y, z)
end function phi

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = ml_stream(self%st, x, y, z)
end function stream

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = ml_phi_t(self%st, x, y, z)
end function phi_t

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = ml_grad_phi(self%st, x, y, z)
end function grad_phi

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(6)
res = 0.0_knd  ! not yet implemented
end function grad_phi_2nd

function acc_euler(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = 0.0_knd  ! not yet implemented
end function acc_euler

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res(3)
res = 0.0_knd  ! not yet implemented
end function acc_particle

function elev(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
res = ml_elev(self%st, x, y)
end function elev

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
res = ml_elev_t(self%st, x, y)
end function elev_t

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res = ml_grad_elev(self%st, x, y)
end function grad_elev

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res = ml_grad_elev_2nd(self%st, x, y)
end function grad_elev_2nd

function pressure(self, x, y, z) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y, z
real(knd) :: res
res = ml_pressure(self%st, x, y, z)
end function pressure

function bathymetry(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res
if (self%d > 0.0_wp) then
    res = real(self%d, knd)
else
    res = -1.0_knd
end if
end function bathymetry

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(in) :: self
real(knd), intent(in) :: x, y
real(knd) :: res(3)
res = [0.0_knd, 0.0_knd, 1.0_knd]
end function bathymetry_nvec

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
real(knd), intent(in) :: x, y, z
character(len=*), intent(in) :: csv
character(len=*), parameter :: err_proc = &
    'spectral_wave_data_shape_1_or_2_impl_7::convergence'
character(len=100) :: err_msg(1)
err_msg(1) = 'convergence() not supported for amp=2 H2 implementation'
call self%error%set_id_msg(err_proc, 1003, err_msg)
end subroutine convergence

subroutine strip(self, tmin, tmax, file_swd)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
real(knd), intent(in) :: tmin, tmax
character(len=*), intent(in) :: file_swd
character(len=*), parameter :: err_proc = &
    'spectral_wave_data_shape_1_or_2_impl_7::strip'
character(len=100) :: err_msg(1)
err_msg(1) = 'strip() not supported for amp=2 H2 implementation'
call self%error%set_id_msg(err_proc, 1003, err_msg)
end subroutine strip

function get_int(self, name) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
character(len=*), intent(in) :: name
integer :: res
select case(name)
case('fmt');    res = self%fmt
case('shp');    res = self%shp
case('amp');    res = self%amp
case('nstrip'); res = self%nstrip
case('nsteps'); res = self%nsteps
case('n');      res = self%n
case('order');  res = self%order
case('norder'); res = self%norder
case('ipol');   res = self%ipol
case('nsumx');  res = self%nsumx
case('nlayers'); res = self%nlayers
case default
    res = 0
end select
end function get_int

function get_logical(self, name) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
character(len=*), intent(in) :: name
logical :: res
select case(name)
case('dc_bias'); res = self%dc_bias
case default;    res = .false.
end select
end function get_logical

function get_real(self, name) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
character(len=*), intent(in) :: name
real(knd) :: res
select case(name)
case('dt');     res = real(self%dt, knd)
case('t0');     res = real(self%t0, knd)
case('x0');     res = real(self%x0, knd)
case('y0');     res = real(self%y0, knd)
case('tmax');   res = real(self%tmax, knd)
case('grav');   res = real(self%grav, knd)
case('lscale'); res = real(self%lscale, knd)
case('rho');    res = real(self%rho, knd)
case('dk');     res = real(self%dk_val, knd)
case('d');      res = real(self%d, knd)
case('zref');   res = real(self%zref, knd)
case('cbeta');  res = real(self%cbeta, knd)
case('sbeta');  res = real(self%sbeta, knd)
case default;   res = 0.0_knd
end select
end function get_real

function get_chr(self, name) result(res)
class(spectral_wave_data_shape_1_or_2_impl_7), intent(inout) :: self
character(len=*), intent(in) :: name
character(len=:), allocatable :: res
select case(name)
case('prog');    res = trim(self%prog)
case('date');    res = trim(self%date)
case('cid');     res = trim(self%cid)
case('version'); res = version
case default;    res = ''
end select
end function get_chr

!==============================================================================

end module spectral_wave_data_shape_1_or_2_impl_7_def
