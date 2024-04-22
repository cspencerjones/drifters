module MOM_particles_mod

use constants_mod,    only: radius, pi, omega, HLF
use MOM_grid,         only : ocean_grid_type
use MOM_time_manager, only : time_type, get_date, operator(-)
use MOM_variables,    only : thermo_var_ptrs

use fms_mod, only: field_exist, get_global_att_value
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_chksum
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP

use mpp_mod, only: mpp_gather
use fms_mod, only: clock_flag_default
use fms_io_mod, only: get_instance_filename
use mpp_domains_mod, only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_parameter_mod, only: CGRID_NE, BGRID_NE, CORNER, AGRID
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST


use diag_manager_mod, only: send_data

use MOM_particles_framework, only: particles_framework_init
use MOM_particles_framework, only: particles_gridded, xyt, particle, particles, buffer
use MOM_particles_framework, only: verbose, really_debug,debug,use_roundoff_fix
use MOM_particles_framework, only: find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell
use MOM_particles_framework, only: bilin,yearday,count_parts,parts_chksum
use MOM_particles_framework, only: linlinx,linliny
use MOM_particles_framework, only: checksum_gridded,add_new_part_to_list
use MOM_particles_framework, only: send_parts_to_other_pes,move_trajectory,move_all_trajectories
use MOM_particles_framework, only: record_posn,check_position,print_part,print_parts,print_fld
use MOM_particles_framework, only: add_new_part_to_list,delete_particle_from_list,destroy_particle
use MOM_particles_framework, only: grd_chksum2,grd_chksum3
use MOM_particles_framework, only: offset_part_dates
use MOM_particles_framework, only: count_parts_in_list,list_chksum
use MOM_particles_framework, only: monitor_a_part,move_part_between_cells, update_halo_particles
use MOM_particles_framework, only: is_point_within_xi_yj_bounds
use MOM_particles_framework, only: find_layer, find_depth


use MOM_particles_io,        only: particles_io_init,write_restart,write_trajectory
use MOM_particles_io,        only: read_restart_parts, particles_io_init



implicit none ; private

public particles_init !, particles_end, particles_run, particles_stock_pe, particles
public particles_end, particles_run, particles
public particles_save_restart
public particles_to_z_space, particles_to_k_space

real, parameter :: pi_180=pi/180.  !< Converts degrees to radians
real, parameter :: r180_pi=180./pi !< Converts radians to degrees
real, parameter :: Rearth=6360000. !< Radius of earth (m)


#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

! ##############################################################################
subroutine particles_init(parts, Grid, Time, dt, u, v, h)

 type(particles), pointer, intent(out) :: parts
 type(ocean_grid_type), target, intent(in) :: Grid !< Grid type from parent model
 type(time_type), intent(in) :: Time !< Time type from parent model
 real, intent(in)            :: dt !< particle timestep in seconds
 real, dimension(:,:,:),intent(in)      :: u, v !< Horizontal velocity fields
 real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers

 integer :: io_layout(2)
 integer :: stdlogunit, stderrunit

 ! Get the stderr and stdlog unit numbers
 stderrunit=stderr()
 stdlogunit=stdlog()
 write(stdlogunit,*) "particles: "//trim(version)

 call particles_framework_init(parts, Grid, Time, dt)
 call particles_to_z_space(parts,h)
 call mpp_clock_begin(parts%clock_ior)
 call particles_io_init(parts,Grid%Domain%io_layout)
 call read_restart_parts(parts,Time, u, v,h)
 call parts_chksum(parts, 'read_restart_particles')
 call mpp_clock_end(parts%clock_ior)



 if (really_debug) call print_parts(stderrunit,parts,'particles_init, initial status')

end subroutine particles_init


subroutine interp_flds(grd, i, j, k, xi, yj, uo, vo)
! Arguments
 type(particles_gridded), pointer :: grd
 integer, intent(in) :: i, j
 real, intent(in) :: k
 real, intent(in) :: xi, yj
 real, intent(out) :: uo, vo
 ! Local variables
 real :: cos_rot, sin_rot
 real :: hxm, hxp
 integer :: kint
 real  ::  xiu,yjv
 integer  ::  iu, jv
 kint = ceiling(k)

 cos_rot=bilin(grd, grd%cos, i, j, xi, yj) ! If true, uses the inverted bilin function
 sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

 yjv=yj+0.5
 if (yjv>1) then
    yjv=yjv-1.
    jv=j+1
 else
    jv=j
 endif
 uo=linlinx(grd, grd%uo(grd%isd:,:,kint), i+1, j, xi, yj)

 xiu = xi+0.5
 if (xiu>1) then
     xiu= xiu-1.
     iu=i+1
 else
    iu=i
 endif
 vo=linliny(grd, grd%vo(:,grd%jsd:,kint), i, j+1, xi, yj)

 ! Rotate vectors from local grid to lat/lon coordinates
 call rotate(uo, vo, cos_rot, sin_rot)


contains


 subroutine rotate(u, v, cos_rot, sin_rot)
   ! Arguments
   real, intent(inout) :: u, v
   real, intent(in) :: cos_rot, sin_rot
   ! Local variables
   real :: u_old, v_old

   u_old=u
   v_old=v
   u=cos_rot*u_old+sin_rot*v_old
   v=cos_rot*v_old-sin_rot*u_old

 end subroutine rotate

end subroutine interp_flds

! ##############################################################################



!> The main driver the steps updates particles
subroutine particles_run(parts, time, uo, vo, ho, tv, use_uh)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:,:),intent(in) :: uo !< Ocean zonal velocity (m/s)
  real, dimension(:,:,:),intent(in) :: vo !< Ocean meridional velocity (m/s)
  real, dimension(:,:,:),intent(in) :: ho !< Ocean layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs), intent(in) :: tv !< structure containing pointers to available thermodynamic fields
  logical :: use_uh !<use uh rather than u

  ! Local variables
  integer :: iyr, imon, iday, ihr, imin, isec, k
  type(particles_gridded), pointer :: grd
  logical :: lerr, sample_traj, write_traj, lverbose
  real :: grdd_u_particle, grdd_v_particle
  integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
  real :: mask
  real, dimension(:,:), allocatable :: uC_tmp, vC_tmp, uA_tmp, vA_tmp
  real, dimension(:,:,:),allocatable :: h_upoints, h_vpoints
  real, dimension(:,:), allocatable :: iCount
  integer :: stderrunit


  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd
!  grd%u_particle(:,:)=0.
!  grd%v_particle(:,:)=0.

  !Initializing _on_ocean_fields
!  grd%Uvel_on_ocean(:,:,:)=0. ;   grd%Vvel_on_ocean(:,:,:)=0.

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  parts%current_year=iyr
  parts%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if ( (parts%traj_sample_hrs>0)  .and. (.not. parts%ignore_traj) ) then
    if (mod(60*60*24*iday+ 60*60*ihr + 60*imin + isec ,60*60*parts%traj_sample_hrs).eq.0) &
      sample_traj=.true.
  endif
  write_traj=.false.
  if ((parts%traj_write_hrs>0) .and. (.not. parts%ignore_traj))  then
    if (mod(60*60*24*(iday-1)+ 60*60*ihr + 60*imin + isec, 60*60*parts%traj_write_hrs).eq.0) &
      write_traj=.true.
  endif
  lverbose=.false.
  if (parts%verbose_hrs>0) then
     if (mod(24*iday+ihr+(imin/60.),float(parts%verbose_hrs)).eq.0) lverbose=verbose
  endif

  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamonds: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', parts%current_year, parts%current_yearday


 ! SPENCER: here is where we need to pass all ocean velocities

  if (use_uh) then
    h_upoints=0.5*(ho(grd%isd:grd%ied,grd%jsd:grd%jed,1:grd%ke)+ho(grd%isd-1:grd%ied-1,grd%jsd:grd%jed,1:grd%ke))
    h_vpoints=0.5*(ho(grd%isd:grd%ied,grd%jsd:grd%jed,1:grd%ke)+ho(grd%isd:grd%ied,grd%jsd-1:grd%jed-1,1:grd%ke))
    do k=1,grd%ke
      grd%uo(:,:,k) = uo(grd%isd:grd%ied,grd%jsd:grd%jed,k) /h_upoints(:,:,k) / parts%dt / grd%dx
      grd%vo(:,:,k) = vo(grd%isd:grd%ied,grd%jsd:grd%jed,k) /h_vpoints(:,:,k) / parts%dt / grd%dx
    enddo
  else
    grd%uo(:,:,:) = uo(:,:,:)
    grd%vo(:,:,:) = vo(:,:,:)
  endif
  do k=2,grd%ke
      grd%hdepth(grd%isd:grd%ied,grd%jsd:grd%jed,k) = grd%hdepth(grd%isd:grd%ied,grd%jsd:grd%jed,k-1)+ho(grd%isd:grd%ied,grd%jsd:grd%jed,k)
  enddo

  ! Make sure that gridded values agree with mask  (to get ride of NaN values)
  do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
    ! Initializing all gridded values to zero
!    if (grd%msk(i,j).lt. 0.5) then
!      grd%uo(i,j,:) = 0.0 ;  grd%vo(i,j,:) = 0.0
!    endif
!    if (grd%uo(i,j) .ne. grd%uo(i,j)) grd%uo(i,j)=0.
!    if (grd%vo(i,j) .ne. grd%vo(i,j)) grd%vo(i,j)=0.
  enddo; enddo

  if (debug) call parts_chksum(parts, 'run parts (top)')
  if (debug) call checksum_gridded(parts%grd, 'top of s/r run')

  !Move to k-space if not already

  call particles_to_k_space(parts,ho)

  if (parts%initial_traj) then 
       call record_posn(parts, ho, tv%T)
       parts%initial_traj=.False.
  endif

  call evolve_particles(parts)
  if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after evolve()     ')
  call move_part_between_cells(parts)  !Markpoint6
  if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after move_lists() ')
  if (debug) call parts_chksum(parts, 'run parts (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(parts%grd, 's/r run after evolve')
  call send_parts_to_other_pes(parts)
  if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after send_parts() ')

  ! For each part, record
  ! sample_traj = .true.
  if (sample_traj) call record_posn(parts, ho, tv%T)
  if (write_traj) then
    call move_all_trajectories(parts)
    call write_trajectory(parts%trajectories, parts%save_short_traj)
  endif

  ! Dump particles to screen
  if (really_debug) call print_parts(stderrunit,parts,'particles_run, status')
  if (debug) call parts_chksum(parts, 'run parts (bot)')
end subroutine particles_run

! ##############################################################################
!> Checks whether all particles are in k-space and if not, moves them to k-space
subroutine particles_to_k_space(parts,h)
   ! Arguments
   type(particles), pointer :: parts !< Container for all types and memory
   real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers 
   !Local variables
   type(particles_gridded), pointer :: grd
   type(particle), pointer :: part
   integer :: grdi, grdj

   integer :: stdlogunit, stderrunit

   ! Get the stderr and stdlog unit numbers                                     
   stderrunit=stderr()


   ! For convenience
   grd=>parts%grd

   do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    part=>parts%list(grdi,grdj)%first
    do while (associated(part)) ! loop over all parts 
    call find_layer(grd, part%depth, h, part%k, part%ine,part%jne, part%xi,part%yj, part%k_space)
    part=>part%next
    enddo
   enddo ; enddo

end subroutine particles_to_k_space

! ##############################################################################

!> Checks whether all particles are in k-space and if not, moves them to k-space

subroutine particles_to_z_space(parts,h)
   ! Arguments                                                                 
   type(particles), pointer :: parts !< Container for all types and memory
   real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers 
   !Local variables
   type(particles_gridded), pointer :: grd
   type(particle), pointer :: part
   integer :: grdi, grdj


   integer :: stdlogunit, stderrunit

   ! Get the stderr and stdlog unit numbers
   stderrunit=stderr()


   ! For convenience 
   grd=>parts%grd

   do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    part=>parts%list(grdi,grdj)%first
    do while (associated(part)) ! loop over all parts 
    call find_depth(grd, part%k, h, part%depth, part%ine,part%jne, part%xi,part%yj, part%k_space)
    part=>part%next
    enddo
   enddo ; enddo

end subroutine particles_to_z_space



! ##############################################################################
!> Evolves particles forward by updating velocity and position with a time-stepping scheme
subroutine evolve_particles(parts)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  ! Local variables
  type(particles_gridded), pointer :: grd
  type(particle), pointer :: part
  real :: uveln, vveln, lonn, latn
  real :: xi, yj
  integer :: i, j
  integer :: grdi, grdj
  integer :: stderrunit
  logical :: bounced, Runge_not_Verlet

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  Runge_not_Verlet=parts%Runge_not_Verlet


  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    part=>parts%list(grdi,grdj)%first
    do while (associated(part)) ! loop over all parts
      !if (part%static_part .lt. 0.5) then  !Only allow non-static particles to evolve !LUYU: not needed, because all the particles are non-static.

        !Checking it everything is ok:
        if (.not. is_point_in_cell(parts%grd, part%lon, part%lat, part%ine, part%jne) ) then
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
          enddo
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
          enddo
          call print_part(stderrunit, part, 'evolve_particle, part is not in proper starting cell')
          write(stderrunit,'(a,i3,2(i4,3f8.2))') 'evolve_particle: pe,lon/lat(i,j)=', mpp_pe(), &
                   part%ine,part%lon,grd%lon(part%ine-1,part%jne-1),grd%lon(part%ine,part%jne), &
                   part%jne,part%lat,grd%lat(part%ine-1,part%jne-1),grd%lat(part%ine,part%jne)
          if (debug) call error_mesg('diamonds, evolve_particle','part is in wrong starting cell!',FATAL)
        endif
        if (debug) call check_position(grd, part, 'evolve_particle (top)')
!       Interpolate gridded velocity fields to part and generate uvel and vvel
        call interp_flds(grd,part%ine,part%jne,part%k,part%xi,part%yj,part%uvel, part%vvel)
          !Time stepping schemes:       
        call Runge_Kutta_stepping(parts,part, uveln, vveln,lonn, latn, i, j, xi, yj)
        
        part%uvel=uveln
        part%vvel=vveln

        if (Runge_not_Verlet) then
          part%lon=lonn  ;   part%lat=latn
          part%ine=i     ;   part%jne=j
          part%xi=xi     ;   part%yj=yj
        else
          !if (.not. interactive_particles_on) call update_verlet_position(parts,part)
        endif

        !call interp_flds(grd, i, j, xi, yj, part%uo, part%vo, part%ui, part%vi, part%ua, part%va, part%ssh_x, part%ssh_y, part%sst)
        !if (debug) call print_part(stderr(), part, 'evolve_particle, final posn.')
        if (debug) call check_position(grd, part, 'evolve_particle (bot)')
      !endif
      part=>part%next
    enddo ! loop over all parts
  enddo ; enddo



end subroutine evolve_particles

!> Calculate explicit and implicit accelerations, new velocity, and new position, using the fourth order Runge-Kutta  method
subroutine Runge_Kutta_stepping(parts, part, uveln, vveln, lonn, latn, i, j, xi, yj)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(particle), pointer, intent(inout) :: part !< particle
  real, intent(out) :: uveln !< New zonal velocity (m/s)
  real, intent(out) :: vveln !< New meridional velocity (m/s)
  real, intent(out) :: lonn !< New longitude (degree E)
  real, intent(out) :: latn !< New latitude (degree N)
  integer, intent(out) :: i !< New i-index of containing cell
  integer, intent(out) :: j !< New i-index of containing cell
  real, intent(out) :: xi !< New non-dimensional x-position
  real, intent(out) :: yj !< New non-dimensional y-position
  ! Local variables
  type(particles_gridded), pointer :: grd
  real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1
  real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2
  real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3
  real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4
  real :: x1, xdot1, xddot1, y1, ydot1, yddot1, xddot1n, yddot1n
  real :: x2, xdot2, xddot2, y2, ydot2, yddot2, xddot2n, yddot2n
  real :: x3, xdot3, xddot3, y3, ydot3, yddot3, xddot3n, yddot3n
  real :: x4, xdot4, xddot4, y4, ydot4, yddot4, xddot4n, yddot4n
  real :: xn, xdotn, xddotn, yn, ydotn, yddotn, xddotnn, yddotnn
  real :: dt, dt_2, dt_6, dydl
  real :: reg_dldx
  integer :: i1,j1,i2,j2,i3,j3,i4,j4
  integer :: stderrunit
  logical :: bounced, on_tangential_plane, error_flag
  ! 4th order Runge-Kutta to solve:
  !    d/dt X = V,  d/dt V = A
  ! with I.C.'s:
  !    X=X1 and V=V1
  !
  !  A1 = A(X1)
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3; A4=A(X4)
  !
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6


  ! Get the stderr unit number
  stderrunit = stderr()
  ! For convenience
  grd=>parts%grd
  ! Common constants
  dt=parts%dt
  dt_2=0.5*dt
  dt_6=dt/6.

  i=part%ine
  j=part%jne
  xi=part%xi
  yj=part%yj
  bounced=.false.
  on_tangential_plane=.false.
  if ((part%lat>89.) .and. (parts%grd%grid_is_latlon)) on_tangential_plane=.true.
  i1=i;j1=j


  reg_dldx=min(abs(grd%lon(i+1,j)-grd%lon(i,j)),abs(grd%lon(i,j)-grd%lon(i-1,j)))/grd%dx(i,j)

  ! A1 = A(X1)
  lon1=part%lon; lat1=part%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)

  call  convert_from_meters_to_grid(lat1,parts%grd%grid_is_latlon,parts%grd%grid_is_regular,dxdl1,dydl,reg_dldx)
  !dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  !dydl=r180_pi/Rearth
  uvel1=part%uvel; vvel1=part%vvel
  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel1,vvel1,xdot1,ydot1)
  u1=uvel1*dxdl1; v1=vvel1*dydl


  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  if (on_tangential_plane) then
    x2=x1+dt_2*xdot1; y2=y1+dt_2*ydot1
    xdot2=xdot1+dt_2*xddot1; ydot2=ydot1+dt_2*yddot1
    call rotpos_from_tang(x2,y2,lon2,lat2)
    call rotvec_from_tang(lon2,xdot2,ydot2,uvel2,vvel2)
  else
    lon2=lon1+dt_2*u1; lat2=lat1+dt_2*v1
    uvel2=uvel1; vvel2=vvel1
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag, part%id)
  i2=i; j2=j
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  call  convert_from_meters_to_grid(lat2,parts%grd%grid_is_latlon,parts%grd%grid_is_regular,dxdl2,dydl,reg_dldx)
  !dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
  if (.not.on_tangential_plane) then
    call interp_flds(grd, i, j, part%k, xi, yj, uvel2, vvel2)
  endif
  u2=uvel2*dxdl2; v2=vvel2*dydl

  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  if (on_tangential_plane) then
    x3=x1+dt_2*xdot2; y3=y1+dt_2*ydot2
    xdot3=xdot1+dt_2*xddot2; ydot3=ydot1+dt_2*yddot2
    call rotpos_from_tang(x3,y3,lon3,lat3)
    call rotvec_from_tang(lon3,xdot3,ydot3,uvel3,vvel3)
  else
    lon3=lon1+dt_2*u2; lat3=lat1+dt_2*v2
    uvel3=uvel1; vvel3=vvel1
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, part%id)
  i3=i; j3=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  call  convert_from_meters_to_grid(lat3,parts%grd%grid_is_latlon,parts%grd%grid_is_regular,dxdl3,dydl,reg_dldx)
  !dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
  if (.not.on_tangential_plane) then
    call interp_flds(grd, i, j, part%k, xi, yj, uvel3, vvel3)
  endif
  u3=uvel3*dxdl3; v3=vvel3*dydl

  !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
  if (on_tangential_plane) then
    x4=x1+dt*xdot3; y4=y1+dt*ydot3
    xdot4=xdot1+dt*xddot3; ydot4=ydot1+dt*yddot3
    call rotpos_from_tang(x4,y4,lon4,lat4)
    call rotvec_from_tang(lon4,xdot4,ydot4,uvel4,vvel4)
  else
    lon4=lon1+dt*u3; lat4=lat1+dt*v3
    uvel4=uvel1; vvel4=vvel1
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag, part%id)
  i4=i; j4=j
  call  convert_from_meters_to_grid(lat4,parts%grd%grid_is_latlon,parts%grd%grid_is_regular,dxdl4,dydl,reg_dldx)
  !dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
  if (.not.on_tangential_plane) then
    call interp_flds(grd, i, j, part%k, xi, yj, uvel4, vvel4)
  endif
  u4=uvel4*dxdl4; v4=vvel4*dydl

  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  if (on_tangential_plane) then
    xn=x1+dt_6*( (xdot1+xdot4)+2.*(xdot2+xdot3) )
    yn=y1+dt_6*( (ydot1+ydot4)+2.*(ydot2+ydot3) )
    xdotn=xdot1+dt_6*( (xddot1+xddot4)+2.*(xddot2+xddot3) )
    ydotn=ydot1+dt_6*( (yddot1+yddot4)+2.*(yddot2+yddot3) )
    xddotn=( (xddot1n+xddot4n)+2.*(xddot2n+xddot3n) )/6.  !Alon
    yddotn=( (yddot1n+yddot4n)+2.*(yddot2n+yddot3n) )/6.  !Alon
    call rotpos_from_tang(xn,yn,lonn,latn)
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
  else
    lonn=part%lon+dt_6*( (u1+u4)+2.*(u2+u3) )
    latn=part%lat+dt_6*( (v1+v4)+2.*(v2+v3) )
    uveln=part%uvel
    vveln=part%vvel
  endif

  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag, part%id)
end subroutine Runge_Kutta_stepping

! ###############################################################################
!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_to_tang(lon, lat, x, y, id_in)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: lat !< Latitude (degree N)
  real, intent(out) :: x !< x-coordinate in tangent plane
  real, intent(out) :: y !< y-coordinate in tangent plane
  integer(kind=8), intent(in), optional :: id_in !< part identifier
  ! Local variables
  real :: r,colat,clon,slon
  integer :: stderrunit, id

  stderrunit = stderr()
  id=0
  if (present(id_in)) then
        id=id_in
  endif

  if (lat>90.) then
      write(stderrunit,*) 'drifters, rotpos_to_tang: lat>90 already!',lat, lon, id
      call error_mesg('drifters, rotpos_to_tang','Something went very wrong!',FATAL)
  endif
  if (lat==90.) then
    write(stderrunit,*) 'drifters, rotpos_to_tang: lat==90 already!',lat, lon
    call error_mesg('drifters, rotpos_to_tang','Something went wrong!',FATAL)
  endif

  colat=90.-lat
  r=Rearth*(colat*pi_180)
  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  x=r*clon
  y=r*slon

end subroutine rotpos_to_tang

! ###################################################################################
!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x !< x-coordinate in tangent plane
  real, intent(in) :: y !< y-coordinate in tangent plane
  real, intent(out) :: lon !< Longitude (degree E)
  real, intent(out) :: lat !< Latitude (degree N)
  ! Local variables
  real :: r

  r=sqrt(x**2+y**2)
  lat=90.-(r180_pi*r/Rearth)
  lon=r180_pi*acos(x/r)*sign(1.,y)

end subroutine rotpos_from_tang

! ###################################################################################
!> Calculate velocity oriented in geographic coordinates from tangent plane velocity
subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(in) :: ydot !< y-component of velocity in tangent plane (m/s)
  real, intent(out) :: uvel !< Zonal velocity (m/s)
  real, intent(out) :: vvel !< Meridional velocity (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  uvel=-slon*xdot+clon*ydot
  vvel=-clon*xdot-slon*ydot

end subroutine rotvec_from_tang

! ###################################################################################
!> Calculates tangent plane velocity from velocity in velocity oriented in geographic coordinates
subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: uvel !< Zonal velocity (m/s)
  real, intent(in) :: vvel !< Meridional velocity (m/s)
  real, intent(out) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(out) :: ydot !< y-component of velocity in tangent plane (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  xdot=-slon*uvel-clon*vvel
  ydot=clon*uvel-slon*vvel

end subroutine rotvec_to_tang

! ####################################################################################
!> Returns metric converting distance in meters to grid distance
subroutine  convert_from_meters_to_grid(lat_ref,grid_is_latlon,grid_is_regular,dlon_dx,dlat_dy,reg_dlon_dx)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  logical, intent(in) :: grid_is_regular !<True if grid is cartesian
  real, intent(out) :: dlon_dx !< Metric dlon/dx
  real, intent(out) :: dlat_dy !< Metric dlat/dy
  real, intent(in), optional :: reg_dlon_dx !<dx used if the grid is cartesian

  if (grid_is_latlon) then
    if (grid_is_regular) then
       dlon_dx=reg_dlon_dx
    else
       dlon_dx=(180./pi)/(Rearth*cos((lat_ref)*(pi/180.)))
    endif
    dlat_dy=(180./pi)/Rearth
  else
    dlon_dx=1.
    dlat_dy=1.
  endif

end subroutine convert_from_meters_to_grid

! ###################################################################################
!> Returns metric converting grid distances to meters
subroutine convert_from_grid_to_meters(lat_ref, grid_is_latlon, grid_is_regular, dx_dlon, dy_dlat,reg_dlon_dx)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  logical, intent(in) :: grid_is_regular !<True if grid is cartesian
  real, intent(out) :: dx_dlon !< Metric dx/dlon
  real, intent(out) :: dy_dlat !< Metric dy/dlat
  real, intent(in), optional :: reg_dlon_dx !<dx used if the grid is cartesian

  if (grid_is_latlon) then
    if (grid_is_regular) then
       dx_dlon=1./reg_dlon_dx
    else
       dx_dlon=(pi/180.)*Rearth*cos((lat_ref)*(pi/180.))
    endif
    dy_dlat=(pi/180.)*Rearth
  else
    dx_dlon=1.
    dy_dlat=1.
  endif

end subroutine convert_from_grid_to_meters

! ###################################################################################
!> Moves part's cell indexes,(i,j), checking for collisional with coasts
subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error, id)
  ! Arguments
  type(particles_gridded), pointer :: grd !< Container for gridded fields
  real, intent(inout) :: lon !< Longitude (degree E)
  real, intent(inout) :: lat !< Latitude (degree N)
  real, intent(inout) :: uvel !< Zonal velocity (m/s)
  real, intent(inout) :: vvel !< Meridional velocity (m/s)
  real, intent(inout) :: xi !< Non-dimension x-position within cell
  real, intent(inout) :: yj !< Non-dimension y-position within cell
  integer, intent(inout) :: i !< i-index of cell
  integer, intent(inout) :: j !< j-index of cell
  logical, intent(out) :: bounced !< True if part collided with coast
  logical, intent(out) :: error !< True if adjustments could not be made consistently
  integer(kind=8), intent(in) :: id !< part identifier
  ! Local variables
  logical lret, lpos
  real, parameter :: posn_eps=0.05
  integer :: icount, i0, j0, inm, jnm
  real :: xi0, yj0, lon0, lat0
  integer :: stderrunit
  logical :: point_in_cell_using_xi_yj

  ! Get the stderr unit number
  stderrunit = stderr()

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
!  print *, 'Alon:', lon, lat, i, j, xi, yj, lret
  xi0=xi; yj0=yj ! original xi,yj


  !Removing this while debuggin
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    point_in_cell_using_xi_yj=is_point_within_xi_yj_bounds(xi,yj)
    if (.not. point_in_cell_using_xi_yj) then

      if (lret) then
        write(stderrunit,*) 'drifters, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'drifters, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'drifters, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'drifters, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'drifters, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'drifters, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'drifters, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,explain=.true.)
        write(stderrunit,*) 'drifters, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
        write(stderrunit,*) 'drifters, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','particle is_point_in_cell=True but xi, yi are out of cell',FATAL)
        error=.true.; return
     endif
    else
      if (.not.lret) then
        write(stderrunit,*) 'drifters, adjust: WARNING!!! lret=F but |xi,yj|<1',mpp_pe()
        write(stderrunit,*) 'drifters, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'drifters, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'drifters, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'drifters, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'drifters, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'drifters, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,  explain=.true.)
        write(stderrunit,*) 'drifters, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'drifters, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','particle is_point_in_cell=False but xi, yi are out of cell',FATAL)
        error=.true.; return
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  endif ! debug

  if (lret) return ! part was already in cell

  ! Find inm, jnm (as if adjusting i,j) based on xi,yj
  ! ignoring the mand mask.
  ! NOTE:  This search appears to have *NO* active role
  ! in the algorithm other than to flag a warning.
  icount=0
  inm=i0; jnm=j0 ! original i,j
  do while (debug .and. .not.lret .and. icount<4)
    icount=icount+1
    if (xi.lt.0.) then
      if (inm>grd%isd) then
        inm=inm-1
      endif
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
!    elseif (yj.ge.1.) then   !Alon:maybe it should be .ge.
      if (jnm<grd%jed) then
        jnm=jnm+1
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj) ! Update xi and yj
  enddo
  if (abs(inm-i0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'drifters, adjust: inm,i0,inm-i0=',inm,i0,inm-i0
   !stop 'Moved too far in i without mask!'
  endif
  if (abs(jnm-j0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'drifters, adjust: jnm,i0,jnm-j0=',jnm,j0,inm-j0
   !stop 'Moved too far in j without mask!'
  endif

  ! Adjust i,j based on xi,yj while bouncing off of masked land cells
  icount=0
  lret=pos_within_cell(grd, lon, lat, i0, j0, xi, yj)
  do while ( .not.lret.and. icount<4 )
    icount=icount+1
    if (xi.lt.0.) then
      if (i>grd%isd) then
!        if (grd%msk(i-1,j)>0.) then
          if (i>grd%isd+1) i=i-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'drifters, adjust: bouncing part from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
!      endif
    elseif (xi.ge.1.) then    !Alon!!!!
!    elseif (xi.gt.1.) then
      if (i<grd%ied) then
!        if (grd%msk(i+1,j)>0.) then
          if (i<grd%ied) i=i+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'drifters, adjust: bouncing part from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
!        endif
      endif
    endif
    if (yj.lt.0.) then
      if (j>grd%jsd) then
!        if (grd%msk(i,j-1)>0.) then
          if (j>grd%jsd+1) j=j-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'drifters, adjust: bouncing part from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
!      endif
    elseif (yj.ge.1.) then     !Alon.
!    elseif (yj.gt.1.) then
      if (j<grd%jed) then
!        if (grd%msk(i,j+1)>0.) then
          if (j<grd%jed) j=j+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'drifters, adjust: bouncing part from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
!      endif
    endif
    if (bounced) then
      if (xi>=1.) xi=1.-posn_eps   !Alon.
!      if (xi>1.) xi=1.-posn_eps   !
      if (xi<0.) xi=posn_eps
      if (yj>=1.) yj=1.-posn_eps  !Alon.
!      if (yj>1.) yj=1.-posn_eps
      if (yj<0.) yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
    endif
    if (debug) then
      if (grd%msk(i,j)==0.) stop 'drifters, adjust: part is in land! This should not happen...'
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  enddo
 !if (debug) then
 !  if (abs(i-i0)>2) then
 !    stop 'drifters, adjust: Moved too far in i!'
 !  endif
 !  if (abs(j-j0)>2) then
 !    stop 'drifters, adjust: Moved too far in j!'
 !  endif
 !endif

  if (.not.bounced.and.lret) return
 ! if (.not.bounced.and.lret.and.grd%msk(i,j)>0.) return ! Landed in ocean without bouncing so all is well

  if (.not.bounced.and..not.lret) then ! This implies the part traveled many cells without getting far enough
    if (debug) then
      write(stderrunit,*) 'drifters, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'drifters, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'drifters, adjust: i0,j0=',i0,j0
      write(stderrunit,*) 'drifters, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'drifters, adjust: xi,yj=',xi,yj
      write(stderrunit,*) 'drifters, adjust: i,j=',i,j
      write(stderrunit,*) 'drifters, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'drifters, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
      write(stderrunit,*) 'drifters, adjust: lret=',lret
    endif

    if (abs(i-i0)+abs(j-j0)==0) then
      if (use_roundoff_fix) then
        ! This is a special case due to round off where is_point_in_cell()
        ! returns false but xi and yj are between 0 and 1.
        ! It occurs very rarely but often enough to have brought down
        ! ESM2G four times since the spin-up began. (as of 8/10/2010)
        ! This temporary fix arbitrarily moves the part toward the
        ! center of the current cell.
        xi=(xi-0.5)*(1.-posn_eps)+0.5
        yj=(yj-0.5)*(1.-posn_eps)+0.5
      endif
      call error_mesg('drifters, adjust', 'part did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
      write(stderrunit,*) 'drifters, adjust: The adjusting particle is: ', id,  mpp_pe()
      write(stderrunit,*) 'drifters, adjust: The adjusting lon,lat,u,v: ', lon, lat, uvel, vvel
      write(stderrunit,*) 'drifters, adjust: The adjusting xi,ji: ', xi, yj
      lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj,explain=.true.)
    else
      call error_mesg('drifters, adjust', 'part iterated many times without bouncing!', WARNING)
    endif
  endif
!  if (xi>1.) xi=1.-posn_eps    !Alon
  if (xi>=1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
!  if (yj>1.) yj=1.-posn_eps
  if (yj<=0.) yj=posn_eps        !Alon
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(0,*) 'i0, j0,=', i0,j0
    write(0,*) 'xi0, yj0,=', xi0,yj0
    write(0,*) 'grd%msk(i0, j0)=', grd%msk(i0,j0)
    write(0,*) 'lon0, lat0,=', lon0,lat0
    write(0,*) 'i,j,lon, lat,grd%msk(i,j)=', i,j,lon,lat,grd%msk(i,j)
    write(stderrunit,*) 'drifters, adjust: Should not get here! part is not in cell after adjustment', id, mpp_pe()
    if (debug) error=.true.
  endif
 end subroutine adjust_index_and_ground

! ###################################################################################

subroutine particles_save_restart(parts, h, directory, time, stamp)
! Arguments
type(particles), pointer :: parts
character(len=*), intent(in) :: directory !< The directory where the restart files are to be written
type(time_type),          intent(in)    :: time       !< The current model time
logical, optional, intent(in)    :: stamp !< If present and true, add time-stamp
real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers



  if (.not.associated(parts)) return
  call mpp_clock_begin(parts%clock_iow)
  call parts_chksum(parts, 'write_restart parts')
  call write_restart(parts, h, directory, time,stamp)
  call mpp_clock_end(parts%clock_iow)

end subroutine particles_save_restart

! ##############################################################################

subroutine particles_end(parts,h,temp,salt)
! Arguments
type(particles), pointer :: parts
real,dimension(:,:,:),intent(in) :: h
real,dimension(:,:,:),optional,intent(in) :: temp, salt

! Local variables
type(particle), pointer :: this, next

  if (.not.associated(parts)) return

  ! particles_save_restart() is called directly by SIS1 and SIS2 so
  ! we do not need to call it a second time. If particles were controlled
  ! by the coupler then the particles would need to take responsibility for
  ! the restarts at the end of the run.
  !if (present(temp) .and. present(salt)) then
  !  call particles_save_restart(parts,h,temp,salt)
  !else
  !  call particles_save_restart(parts,h)
  !endif

  call mpp_clock_begin(parts%clock_ini)
  ! Delete parts and structures
  call move_all_trajectories(parts, delete_parts=.true.)

  !call write_trajectory(parts%trajectories)

  deallocate(parts%grd%lon)
  deallocate(parts%grd%lat)
  deallocate(parts%grd%lonc)
  deallocate(parts%grd%latc)
  !deallocate(parts%grd%dx)
  !deallocate(parts%grd%dy)
  !deallocate(parts%grd%area)
  deallocate(parts%grd%msk)
  deallocate(parts%grd%cos)
  deallocate(parts%grd%sin)
  deallocate(parts%grd%ocean_depth)
! ! deallocate(parts%grd%domain)
! ! deallocate(parts%grd)
  call dealloc_buffer(parts%obuffer_n)
  call dealloc_buffer(parts%obuffer_s)
  call dealloc_buffer(parts%obuffer_e)
  call dealloc_buffer(parts%obuffer_w)
  call dealloc_buffer(parts%ibuffer_n)
  call dealloc_buffer(parts%ibuffer_s)
  call dealloc_buffer(parts%ibuffer_e)
  call dealloc_buffer(parts%ibuffer_w)
  call dealloc_buffer(parts%ibuffer_io)
  call dealloc_buffer(parts%obuffer_io)
  call mpp_clock_end(parts%clock_ini)
  !deallocate(parts)

  if (mpp_pe()==mpp_root_pe()) write(*,'(a,i8)') 'drifters: particles_end complete',mpp_pe()

  contains

  subroutine dealloc_buffer(buff)
  ! Arguments
  type(buffer), pointer :: buff
  ! Local variables
    if (associated(buff)) then
      if (associated(buff%data)) deallocate(buff%data)
      deallocate(buff)
    endif
  end subroutine dealloc_buffer

end subroutine particles_end

!> Returns true of a point is in or on the rectangle with opposite corners A and B
logical function point_in_interval(Ax, Ay, Bx, By, px, py)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: px !< x-position of point
  real, intent(in) :: py !< y-position of point
  point_in_interval=.False.
  if ((px <= max(Ax,Bx)) .and. (px >= min(Ax,Bx))) then
    if ((py <= max(Ay,By)) .and. (py >= min(Ay,By))) then
      point_in_interval=.True.
    endif
  endif
end function point_in_interval

!> Returns true if point q is on a line through points A and B
logical function point_is_on_the_line(Ax, Ay, Bx, By, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of point A
  real, intent(in) :: Ay !< y-position of point A
  real, intent(in) :: Bx !< x-position of point B
  real, intent(in) :: By !< y-position of point B
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: tol, dxc,dyc,dxl,dyl,cross
  !tol=1.e-12
  tol=0.0
  dxc = qx - Ax
  dyc = qy - Ay
  dxl = Bx - Ax
  dyl = By - Ay
  cross = dxc * dyl - dyc * dxl
  if (abs(cross)<=tol) then
    point_is_on_the_line=.True.
  else
   point_is_on_the_line=.False.
  endif
end function point_is_on_the_line

!> Returns True if a point q is inside a triangle ABC
!!
!! This function decides whether a point (qx,qy) is inside the triangle ABC.
!! There is also the option to include the boundary of the triangle.
logical function point_in_triangle(Ax, Ay, Bx, By, Cx, Cy, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: l0,l1,l2,p0,p1,p2
  real :: v0x,v1x,v2x,v0y,v1y,v2y,dot00,dot01,dot02,dot11,dot12

  point_in_triangle = .False.
  if ((Ax==qx .and. Ay==qy) .or. (Bx==qx .and. By==qy) .or. (Cx==qx .and. Cy==qy)) then ! Exclude the pathelogical case
      point_in_triangle = .False.
  else
    if (((point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) .or. (point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy))) .or. (point_is_on_the_line(Bx,By,Cx,Cy,qx,qy)))) then
      point_in_triangle = .False.
    else
      ! Compute point in triangle using Barycentric coordinates (the same as sum_sign_dot_prod routines)
      l0=(qx-Ax)*(By-Ay)-(qy-Ay)*(Bx-Ax)
      l1=(qx-Bx)*(Cy-By)-(qy-By)*(Cx-Bx)
      l2=(qx-Cx)*(Ay-Cy)-(qy-Cy)*(Ax-Cx)

      p0=sign(1., l0); if (l0==0.)  p0=0.
      p1=sign(1., l1); if (l1==0.)  p1=0.
      p2=sign(1., l2); if (l2==0.)  p2=0.

      if ( (abs(p0)+abs(p2))+(abs(p1)) == abs((p0+p2)+(p1)) )  point_in_triangle = .True.
    endif
  endif
end function point_in_triangle

!> Returns area of a triangle
real function Area_of_triangle(Ax, Ay, Bx, By, Cx, Cy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  Area_of_triangle    =   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) )
end function Area_of_triangle



end module MOM_particles_mod
