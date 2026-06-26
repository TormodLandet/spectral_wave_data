module swd_write_shape_7_def

use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float, c_null_char

implicit none
private

! Please adjust wp to the working precision applied in the wave generator
integer, parameter :: wp = kind(1.0d0) ! (double precision)


! This module provides a class for writing SWD files for shape 7.
!
! Shape 7: long-crested waves with a sigma-coordinate multi-layer potential.
! The velocity potential is stored on nlayers horizontal sigma-layers that
! deform with the free surface.  No time-derivative arrays (ht, ct) are
! stored; the Fortran reader reconstructs them numerically.
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_write_shape_7
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type swd_write_shape_7
    integer :: lu        ! Unit number associated with swd file
    integer :: nsteps    ! Expected number of time steps
    integer :: it        ! Current time step count
    integer :: n         ! Highest spectral index (n_swd in theory; n+1 components)
    integer :: nlayers   ! Number of sigma-layers
contains
    procedure :: init    ! Open the swd file and write the header section
    procedure :: update  ! Add spectral data at current time step
    procedure :: close   ! Close the file
end type swd_write_shape_7

contains

!==============================================================================

subroutine init(self, file, prog, cid, grav, lscale, amp, &
                n, order, dk, dt, nsteps, d, zref, nlayers, sig)
! Open a new SWD file and write the common header plus the shape-7 header.
class(swd_write_shape_7), intent(out) :: self   ! Object to create
character(len=*), intent(in) :: file         ! Requested name for new SWD file
character(len=*), intent(in) :: prog         ! Name of wave generator (incl. version)
character(len=*), intent(in) :: cid          ! Text describing input parameters
real(wp),         intent(in) :: grav         ! Acceleration of gravity [m/s^2]
real(wp),         intent(in) :: lscale       ! Length units per metre
integer,          intent(in) :: amp          ! Amplitude content flag
integer,          intent(in) :: n            ! Highest spectral index
integer,          intent(in) :: order        ! Expansion order (< 0 for fully nonlinear)
real(wp),         intent(in) :: dk           ! Wave-number spacing [1/m]
real(wp),         intent(in) :: dt           ! Time step [s]
integer,          intent(in) :: nsteps       ! Total number of time steps to write
real(wp),         intent(in) :: d            ! Water depth [m], or -1 for infinite depth
real(wp),         intent(in) :: zref         ! z-position of the bottom sigma-layer (sigma=0)
integer,          intent(in) :: nlayers      ! Number of sigma-layers
real(wp),         intent(in) :: sig(nlayers) ! Sigma-layer positions sig(1)=0.0, sig(nlayers)=1.0
!
integer :: nid
integer, parameter :: fmt = 100
integer, parameter :: shp = 7
integer, parameter :: nstrip = 0
integer :: m
! C-compatible character variables
character(kind=c_char, len=:), allocatable :: ccid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
!
open(newunit=self % lu, file=file, access='stream', form='unformatted', &
     status='replace', action='write', convert='little_endian')
! The CONVERT option is an Intel, GFortran, HP and IBM extension.
! Modify this parameter for other compilers if your compiler complains.
!
self % nsteps = nsteps
self % it = 0
self % n = n
self % nlayers = nlayers
!
nid = len_trim(cid) + 1  ! one extra for the C null character
allocate(character(len=nid) :: ccid)
ccid = trim(cid) // c_null_char
cdate = timestamp() // c_null_char
cprog = trim(prog) // c_null_char
!
! ---- Common SWD header ----
write(self % lu) 37.0221_c_float        ! magic number
write(self % lu) int(fmt,    c_int)     ! file format version
write(self % lu) int(shp,    c_int)     ! shape class = 7
write(self % lu) int(amp,    c_int)     ! amplitude content flag
write(self % lu) cprog                  ! wave generator name (30 chars)
write(self % lu) cdate                  ! creation timestamp (20 chars)
write(self % lu) int(nid,    c_int)     ! length of cid string (incl. null)
write(self % lu) ccid(:nid)             ! input description string
write(self % lu) real(grav,  c_float)   ! gravity
write(self % lu) real(lscale,c_float)   ! length scale
write(self % lu) int(nstrip, c_int)     ! nstrip (always 0 for shape 7)
write(self % lu) int(nsteps, c_int)     ! number of time steps
write(self % lu) real(dt,    c_float)   ! time step
write(self % lu) int(order,  c_int)     ! expansion order
! ---- Shape-7 header ----
write(self % lu) int(n,      c_int)     ! highest spectral index
write(self % lu) real(dk,    c_float)   ! wave-number spacing
write(self % lu) real(d,     c_float)   ! water depth (or -1)
write(self % lu) real(zref,  c_float)   ! z-position of bottom sigma-layer
write(self % lu) int(nlayers, c_int)    ! number of sigma-layers
do m = 1, nlayers
    write(self % lu) real(sig(m), c_float) ! sigma position for layer m
end do

contains

  function timestamp() result(timestr)
    ! Return a string with the current time in the form YYYY-MM-DD hh-mm-ss
    character (len = 19) :: timestr
    integer :: y, mo, d_val, h, mn, s
    integer :: values(8)
    call date_and_time(values=values)
    y    = values(1)
    mo   = values(2)
    d_val= values(3)
    h    = values(5)
    mn   = values(6)
    s    = values(7)
    write(timestr, '(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
         y,'-',mo,'-',d_val,' ',h,':',mn,':',s
  end function timestamp

end subroutine init

!==============================================================================

subroutine update(self, h, c_layers)
! Write one time step: elevation amplitudes h(0:n) followed by the potential
! amplitudes c_layers(0:n, 1:nlayers) for all sigma-layers.
! No time derivatives are stored.
class(swd_write_shape_7), intent(inout) :: self     ! Object to update
complex(wp), intent(in) :: h(0:)                    ! h(0:n)  elevation spectral amplitudes
complex(wp), intent(in) :: c_layers(0:, 1:)         ! c_layers(0:n, 1:nlayers)  potential per layer
!
integer :: j, m
!
self % it = self % it + 1
!
! Write elevation amplitudes
do j = 0, self % n
    write(self % lu) cmplx(h(j), kind=c_float)
end do
!
! Write potential amplitudes for each sigma-layer
do m = 1, self % nlayers
    do j = 0, self % n
        write(self % lu) cmplx(c_layers(j, m), kind=c_float)
    end do
end do
!
end subroutine update

!==============================================================================

subroutine close(self)
! Validate the number of written time steps and close the file.
class(swd_write_shape_7), intent(inout) :: self
!
if (self % it /= self % nsteps) then
    print*, "WARNING: from swd_write_shape_7 :: close"
    print*, "Specified number of time steps = ", self % nsteps
    print*, "Number of provided time steps  = ", self % it
end if
!
close(self % lu)
!
end subroutine close

!==============================================================================

end module swd_write_shape_7_def
