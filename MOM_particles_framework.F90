!> This is the core module for MOM drifter diagnostics

module MOM_particles_framework

! This file is part of MOM6. See LICENSE.md for the license.

use constants_mod, only: radius, pi, omega, HLF
use MOM_grid, only : ocean_grid_type
use MOM_time_manager, only : get_time
use MOM_file_parser, only : get_param

use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_pe, mpp_root_pe, mpp_chksum
use mpp_mod, only: COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_mod, only: COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8
use mpp_mod, only: COMM_TAG_9, COMM_TAG_10
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id, input_nml_file
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use mpp_domains_mod, only: domain2D
use mpp_parameter_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE, CORNER, AGRID,CENTER
use mpp_domains_mod, only: mpp_update_domains, mpp_define_domains
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use mpp_domains_mod, only: mpp_define_io_domain
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: clock_flag_default
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init

implicit none ; private

integer :: buffer_width=20 ! size of buffer dimension for comms
integer :: buffer_width_traj=17  
logical :: folded_north_on_pe = .false. !< If true, indicates the presence of the tri-polar grid
logical :: verbose=.false. !< Be verbose to stderr
logical :: debug=.false. !< Turn on debugging
logical :: really_debug=.false. !< Turn on debugging
logical :: parallel_reprod=.true. !< Reproduce across different PE decompositions
logical :: use_slow_find=.true. !< Use really slow (but robust) find_cell for reading restarts
logical :: ignore_ij_restart=.false. !< Read i,j location from restart if available (needed to use restarts on different grids)
logical :: use_roundoff_fix=.true. !< Use a "fix" for the round-off discrepancy between is_point_in_cell() and pos_within_cell()
logical :: old_bug_bilin=.true. !< If true, uses the inverted bilinear function (use False to get correct answer)
character(len=10) :: restart_input_dir = 'INPUT/' !< Directory to look for restart files
integer, parameter :: delta_buf=25 !< Size by which to increment buffers
real, parameter :: pi_180=pi/180. !< Converts degrees to radians
logical :: fix_restart_dates=.true. !< After a restart, check that parts were created before the current model date
logical :: do_unit_tests=.false. !< Conduct some unit tests
logical :: force_all_pes_traj=.false. !< Force all pes write trajectory files regardless of io_layout


!Public params !Niki: write a subroutine to expose these
public buffer_width,buffer_width_traj
public verbose, really_debug, debug, restart_input_dir,old_bug_bilin,use_roundoff_fix
public ignore_ij_restart, use_slow_find
public force_all_pes_traj

!Public types
public particles_gridded, xyt, particle, particles, buffer!, bond

!Public subs
public particles_framework_init
public send_parts_to_other_pes
public update_halo_particles
public pack_traj_into_buffer2, unpack_traj_from_buffer2
public increase_ibuffer
public add_new_part_to_list, count_out_of_order, check_for_duplicates
public insert_part_into_list, create_particle, delete_particle_from_list, destroy_particle
public print_fld,print_part, print_parts,record_posn, push_posn, append_posn, check_position
public move_trajectory, move_all_trajectories
public find_cell, find_cell_by_search, count_parts, is_point_in_cell, pos_within_cell
public find_layer, find_depth
public bilin, yearday, parts_chksum, list_chksum, count_parts_in_list
public find_u, find_v, jacob2D
public checksum_gridded
public grd_chksum2,grd_chksum3
public fix_restart_dates, offset_part_dates
public move_part_between_cells
public find_individual_particle
public monitor_a_part
public is_point_within_xi_yj_bounds
public test_check_for_duplicate_ids_in_list
public check_for_duplicates_in_parallel
public split_id, id_from_2_ints, generate_id

!> Container for gridded fields
type :: particles_gridded
  type(domain2D), pointer :: domain !< MPP parallel domain
  integer :: halo !< Nominal halo width
  integer :: isc !< Start i-index of computational domain
  integer :: iec !< End i-index of computational domain
  integer :: jsc !< Start j-index of computational domain
  integer :: jec !< End j-index of computational domain
  integer :: isd !< Start i-index of data domain
  integer :: ied !< End i-index of data domain
  integer :: jsd !< Start j-index of data domain
  integer :: jed !< End j-index of data domain
  integer :: isg !< Start i-index of global domain
  integer :: ieg !< End i-index of global domain
  integer :: jsg !< Start j-index of global domain
  integer :: jeg !< End j-index of global domain
  integer :: ke         !< The number of layers in the vertical.
  integer :: is_offset=0 !< add to i to recover global i-index
  integer :: js_offset=0 !< add to j to recover global j-index
  integer :: my_pe !< MPI PE index
  integer :: pe_N !< MPI PE index of PE to the north
  integer :: pe_S !< MPI PE index of PE to the south
  integer :: pe_E !< MPI PE index of PE to the east
  integer :: pe_W !< MPI PE index of PE to the west
  logical :: grid_is_latlon !< Flag to say whether the coordinate is in lat-lon degrees, or meters
  logical :: grid_is_regular !< Flag to say whether point in cell can be found assuming regular Cartesian grid
  logical :: no_TS
  real :: Lx !< Length of the domain in x direction
  real, dimension(:,:), allocatable :: lon !< Longitude of cell corners (degree E)
  real, dimension(:,:), allocatable :: lat !< Latitude of cell corners (degree N)
  real, dimension(:,:), allocatable :: lonc !< Longitude of cell centers (degree E)
  real, dimension(:,:), allocatable :: latc !< Latitude of cell centers (degree N)
  real, dimension(:,:), allocatable :: dx !< Length of cell edge (m)
  real, dimension(:,:), allocatable :: dy !< Length of cell edge (m)
  real, dimension(:,:), allocatable :: area !< Area of cell (m^2)
  real, dimension(:,:), allocatable :: msk !< Ocean-land mask (1=ocean)
  real, dimension(:,:), allocatable :: cos !< Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), allocatable :: sin !< Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), allocatable :: ocean_depth !< Depth of ocean (m)
  real, dimension(:,:,:), allocatable :: uo !< Ocean zonal flow (m/s)
  real, dimension(:,:,:), allocatable :: vo !< Ocean meridional flow (m/s)
  real, dimension(:,:,:), allocatable :: hdepth !< Cumulative thickness from ocen surf
  real, dimension(:,:), allocatable :: tmp !< Temporary work space
  real, dimension(:,:), allocatable :: tmpc !< Temporary work space
  real, dimension(:,:), allocatable :: parity_x !< X component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  real, dimension(:,:), allocatable :: parity_y !< Y component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  integer, dimension(:,:), allocatable :: particle_counter_grd !< Counts particles created for naming purposes
  !>@{
  !! Diagnostic handle
  integer :: id_uo=-1, id_vo=-1, id_unused=-1
  integer :: id_count=-1, id_chksum=-1
  !>@}

end type particles_gridded

!>xyt is a data structure containing particle position and velocity fields.
type :: xyt
  real :: lon, lat, day      !< Current position (degrees) and day
  real :: lat_old, lon_old   !< Previous position (degrees)
  real :: uvel, vvel         !< Current velocity components (m/s)
  real :: uvel_old, vvel_old !< Previous velocity components (m/s)
  real :: theta
  integer :: year, particle_num  !< Current year and particle number
  integer(kind=8) :: id = -1 !< Particle Identifier
  real :: k !<Current vertical level i which the particle resides
  real :: depth !<Current depth of the particle
  logical :: k_space      !<Logical indicating whether particle is in k (vs z)
  type(xyt), pointer :: next=>null()  !< Pointer to the next position in the list
end type xyt

!>particle types are data structures describing a tracked particle
type :: particle
  type(particle), pointer :: prev=>null(), next=>null()
  ! State variables (specific to the particle, needed for restarts)
  real :: lon, lat, depth, uvel, vvel !< position (degrees) and zonal and meridional velocities (m/s)
  real :: theta !<temperature at the location of the particle
  real :: lon_old, lat_old, uvel_old, vvel_old  !< previous position (degrees) and zonal
                                                !< and meridional velocities (m/s)
  real :: start_lon, start_lat, start_d, start_day       !< origination position (degrees) and day
  integer :: start_year                         !< origination year
  real :: halo_part  !< equal to zero for particles on the computational domain, and 1 for particles on the halo
  integer(kind=8) :: id,drifter_num             !< particle identifier
  integer :: ine, jne                           !< nearest index in NE direction (for convenience)
  real :: k                 !<vertical level of particle
  logical :: k_space        !<flag for whether depth is stored in kspace (vs z)
  real :: xi, yj                                !< non-dimensional coords within current cell (0..1)
  real :: uo, vo                                !< zonal and meridional ocean velocities experienced
  real :: hdepth                                !<depth from surface at bottom of layer
                                                !< by the particle (m/s)
  type(xyt), pointer :: trajectory=>null()
end type particle

!>A buffer structure for message passing
type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

!> A wrapper for the particle linked list (since an array of pointers is not allowed)
type :: linked_list
  type(particle), pointer :: first=>null() !< Pointer to the beginning of a linked list of parts
end type linked_list

!> A grand data structure for the particles in the local MOM domain
type :: particles !; private
  type(particles_gridded) :: grd !< Container with all gridded data
  type(linked_list), dimension(:,:), allocatable :: list !< Linked list of particles
  type(xyt), pointer :: trajectories=>null() !< A linked list for detached segments of trajectories
  real :: dt !< Time-step between particle calls
  integer :: current_year !< Current year (years)
  real :: current_yearday !< Current year-day, 1.00-365.99, (days)
  integer :: traj_sample_hrs !< Period between sampling for trajectories (hours)
  integer :: traj_write_hrs !< Period between writing of trajectories (hours)
  integer :: verbose_hrs !< Period between terminal status reports (hours)
  !>@{
  !! Handles for clocks
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com, clock_ini, clock_ior, clock_iow, clock_dia
  integer :: clock_trw, clock_trp
  !>@}
  logical :: restarted=.false. !< Indicate whether we read state from a restart or not
  logical :: Runge_not_Verlet=.True. !< True=Runge-Kutta, False=Verlet.
  logical :: xystagger=.False.
  logical :: ignore_missing_restart_parts=.False. !< True allows the model to ignore particles missing in the restart.
  logical :: halo_debugging=.False. !< Use for debugging halos (remove when its working)
  logical :: save_short_traj=.false. !< True saves only lon,lat,time,id in particle_trajectory.nc
  logical :: ignore_traj=.False. !< If true, then model does not write trajectory data at all
  logical :: initial_traj=.True. !< If true, then model will write trajectory data before starting the run
  logical :: use_new_predictive_corrective =.False. !< Flag to use Bob's predictive corrective particle scheme- Added by Alon
  integer(kind=8) :: debug_particle_with_id = -1 !< If positive, monitors a part with this id
  integer :: nstep=1 !< Number of timestepped (used for staggering particle motion in x and y
  type(buffer), pointer :: obuffer_n=>null() !< Buffer for outgoing parts to the north
  type(buffer), pointer :: ibuffer_n=>null() !< Buffer for incoming parts from the north
  type(buffer), pointer :: obuffer_s=>null() !< Buffer for outgoing parts to the south
  type(buffer), pointer :: ibuffer_s=>null() !< Buffer for incoming parts from the south
  type(buffer), pointer :: obuffer_e=>null() !< Buffer for outgoing parts to the east
  type(buffer), pointer :: ibuffer_e=>null() !< Buffer for incoming parts from the east
  type(buffer), pointer :: obuffer_w=>null() !< Buffer for outgoing parts to the west
  type(buffer), pointer :: ibuffer_w=>null() !< Buffer for incoming parts from the west
  type(buffer), pointer :: obuffer_io=>null() !< Buffer for outgoing parts during i/o
  type(buffer), pointer :: ibuffer_io=>null() !< Buffer for incoming parts during i/o
end type particles


#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

!> Set a value in the buffer at position (counter,n) after incrementing counter
interface push_buffer_value
  module procedure push_buffer_rvalue, push_buffer_ivalue
end interface

!> Get a value in the buffer at position (counter,n) after incrementing counter
interface pull_buffer_value
  module procedure pull_buffer_rvalue, pull_buffer_ivalue
end interface

contains


! ##############################################################################

subroutine particles_framework_init(parts, Grid, Time, dt)
  type(particles), pointer :: parts !< Particles to be allocated
  type(ocean_grid_type), target, intent(in) :: Grid !< MOM6 grid
  real, intent(in) :: dt !< Time step (s)
  type(time_type), intent(in) :: Time !< Model time

  ! Namelist parameters (and defaults)
  integer :: halo=4 ! Width of halo region
  integer :: traj_sample_hrs=24 ! Period between sampling of position for trajectory storage
  integer :: traj_write_hrs=24 ! Period between writing sampled trajectories to disk
  integer :: verbose_hrs=24 ! Period between verbose messages
  real :: Lx=360. ! Length of domain in x direction, used for periodicity (use a huge number for non-periodic)
  logical :: Runge_not_Verlet=.True. ! True=Runge Kutta, False=Verlet.
  logical :: xystagger=.True. !If true then stagger timestepping, with particle first advected in the xy direction and then in the y direction
  logical :: grid_is_latlon=.True. ! True means that the grid is specified in lat lon, and uses to radius of the earth to convert to distance
  logical :: grid_is_regular=.False. ! Flag to say whether point in cell can be found assuming regular Cartesian grid
  logical :: ignore_missing_restart_parts=.False. ! True Allows the model to ignore particles missing in the restart.
  logical :: no_TS=.False.
  logical :: halo_debugging=.False. ! Use for debugging halos (remove when its working)
  logical :: save_short_traj=.false. ! True saves only lon,lat,time,id in particle_trajectory.nc
  logical :: ignore_traj=.False. ! If true, then model does not traj trajectory data at all
  logical :: initial_traj=.True. !< If true, then model will write trajectory data before starting the run
  logical :: use_new_predictive_corrective =.False. ! Flag to use Bob's predictive corrective particle scheme- Added by Alon
  logical :: do_unit_tests=.false. ! Conduct some unit tests
  logical :: input_freq_distribution=.false. ! Flag to show if input distribution is freq or mass dist (=1 if input is a freq dist, =0 to use an input mass dist)
  logical :: read_old_restarts=.false. ! Legacy option that does nothing
  integer(kind=8) :: debug_particle_with_id = -1 ! If positive, monitors a part with this id
  integer :: generate_days=-1 ! If positive, is the period in days between generation of new particles on a grid. If 0, generate once. Negative do nothing.
  real :: generate_lons(3) ! Start,end and delta longitude, if generating particles
  real :: generate_lats(3) ! Start,end and delta longitude, if generating particles
  real :: generate_d(3) ! Start, end and delta depth if positive. Start end and delta k if negative

  namelist /particles_nml/ verbose, halo,  traj_sample_hrs, traj_write_hrs, save_short_traj,  &
         verbose_hrs,  &
         debug, really_debug, ignore_missing_restart_parts, &
         parallel_reprod, use_slow_find, ignore_ij_restart, use_new_predictive_corrective, halo_debugging, &
         fix_restart_dates, use_roundoff_fix, Runge_not_Verlet, &
         restart_input_dir, old_bug_bilin,do_unit_tests, force_all_pes_traj, &
         grid_is_latlon,Lx, no_TS, &
         grid_is_regular, &
         generate_days, generate_lons, generate_lats, generate_d, &
         ignore_traj, initial_traj, debug_particle_with_id, read_old_restarts

  ! Local variables
  integer :: ierr, iunit, i, j, id_class, is, ie, js, je, np
  integer :: iyr, imon, iday, ihr, imin, isec
  type(particles_gridded), pointer :: grd
  real :: lon_mod, big_number
  logical :: lerr, lgenerate
  integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "particles_framework: "//trim(version)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=particles_nml, iostat=ierr)
#else
  iunit = open_namelist_file()
  read  (iunit, particles_nml,iostat=ierr)
  call close_file (iunit)
#endif
  ierr = check_nml_error(ierr,'particles_nml')

  if (really_debug) debug=.true. ! One implies the other...

  write (stdlogunit, particles_nml)

  ! Allocate memory
  allocate(parts)

  grd=>parts%grd
  grd%domain => Grid%domain%mpp_domain

! Clocks
  parts%clock=mpp_clock_id( 'Particles', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  parts%clock_com=mpp_clock_id( 'Particles-communication', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  parts%clock_ini=mpp_clock_id( 'Particles-initialization', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  parts%clock_ior=mpp_clock_id( 'Particles-I/O read', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  parts%clock_iow=mpp_clock_id( 'Particles-I/O write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )


  call mpp_clock_begin(parts%clock)
  call mpp_clock_begin(parts%clock_ini)

  grd%isg = Grid%isg; grd%ieg = Grid%ieg
  grd%jsg = Grid%jsg; grd%jeg = Grid%jeg
  grd%isc = Grid%isc; grd%iec = Grid%iec
  grd%jsc = Grid%jsc; grd%jec = Grid%jec
  grd%isd = Grid%isd; grd%ied = Grid%ied
  grd%jsd = Grid%jsd; grd%jed = Grid%jed
  grd%ke = Grid%ke
  grd%is_offset = Grid%idg_offset
  grd%js_offset = Grid%jdg_offset

  call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
  call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
  call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
  call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)


  folded_north_on_pe = ((Grid%Domain%y_flags == FOLD_NORTH_EDGE) .and. (grd%jec + grd%js_offset == grd%jeg))
 ! Allocate grid of pointers
  allocate( parts%list(grd%isd:grd%ied, grd%jsd:grd%jed) )
  do j = grd%jsd,grd%jed ; do i = grd%isd,grd%ied
    parts%list(i,j)%first => null()
  enddo ; enddo

  big_number=1.0E15
  allocate( grd%lon(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=big_number
  allocate( grd%lat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=big_number
  allocate( grd%lonc(grd%isd:grd%ied, grd%jsd:grd%jed) );grd%lon(:,:)=big_number
  allocate( grd%latc(grd%isd:grd%ied, grd%jsd:grd%jed) );grd%lat(:,:)=big_number
  allocate( grd%dx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dx(:,:)=0.
  allocate( grd%dy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dy(:,:)=0.
  allocate( grd%area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%area(:,:)=0.
  allocate( grd%msk(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%msk(:,:)=0.
  allocate( grd%cos(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cos(:,:)=1.
  allocate( grd%sin(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sin(:,:)=0.
  allocate( grd%ocean_depth(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ocean_depth(:,:)=0.
  allocate( grd%uo(grd%isd:grd%ied, grd%jsd:grd%jed,grd%ke) ); grd%uo(:,:,:)=0.
  allocate( grd%vo(grd%isd:grd%ied, grd%jsd:grd%jed,grd%ke) ); grd%vo(:,:,:)=0.
  allocate( grd%tmp(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%tmp(:,:)=0.
  allocate( grd%tmpc(grd%isc:grd%iec, grd%jsc:grd%jec) ); grd%tmpc(:,:)=0.
  allocate( grd%parity_x(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_x(:,:)=1.
  allocate( grd%parity_y(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_y(:,:)=1.
  allocate( grd%particle_counter_grd(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%particle_counter_grd(:,:)=0


  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
  grd%lon(is:ie,js:je)=Grid%geolonBu(is:ie,js:je)
  grd%lat(is:ie,js:je)=Grid%geolatBu(is:ie,js:je)
  grd%area(is:ie,js:je)=Grid%areaBu(is:ie,js:je) !sis2 has *(4.*pi*radius*radius)
  grd%ocean_depth(is:ie,js:je) = Grid%bathyT(is:ie,js:je)
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
!  grd%dx(is:ie,js:je)=Grid%dxCu(is:ie,js:)!je
  grd%dx(:,:)=Grid%dxCu(:,:)
!  grd%dy(is:ie,js:je)=Grid%dyCv(is:ie,js:)!je
  grd%dy(:,:)=Grid%dyCv(:,:)!je
  grd%msk(is:ie,js:je)=Grid%mask2dBu(is:ie,js:je)
  grd%cos(is:ie,js:je)=Grid%cos_rot(is:ie,js:je)
  grd%sin(is:ie,js:je)=Grid%sin_rot(is:ie,js:je)

  call mpp_update_domains(grd%lon, grd%domain)
  call mpp_update_domains(grd%lat, grd%domain)
  call mpp_update_domains(grd%dy, grd%dx, grd%domain, gridtype=CGRID_NE, flags=SCALAR_PAIR)
  call mpp_update_domains(grd%area, grd%domain)
  call mpp_update_domains(grd%msk, grd%domain)
  call mpp_update_domains(grd%cos, grd%domain, position=CENTER)
  call mpp_update_domains(grd%sin, grd%domain, position=CENTER)
  call mpp_update_domains(grd%ocean_depth, grd%domain)
  call mpp_update_domains(grd%parity_x, grd%parity_y, grd%domain, gridtype=AGRID) ! If either parity_x/y is -ve, we need rotation of vectors



  ! Sanitize lon and lat in the southern halo
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=grd%lon(i,j+1)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo

  ! fix halos on edge of the domain
  !1) South
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i,j+1)-grd%lon(i,j+2)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo
  !2) North
  do j=grd%jec+1,grd%jed; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i,j-1)-grd%lon(i,j-2)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j-1)-grd%lat(i,j-2)
  enddo; enddo
  !3) West
  do i=grd%isc-1,grd%isd,-1; do j=grd%jsd,grd%jed
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i+1,j)-grd%lon(i+2,j)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i+1,j)-grd%lat(i+2,j)
  enddo; enddo
  !4) East
  do i=grd%iec+1,grd%ied; do j=grd%jsd,grd%jed
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i-1,j)-grd%lon(i-2,j)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i-1,j)-grd%lat(i-2,j)
  enddo; enddo


  if ((Lx.gt.1E15 ) .and. (mpp_pe().eq.mpp_root_pe())) then
          call error_mesg('particles, framework', 'Model does not enjoy the domain being larger than 1E15. Not sure why. Probably to do with floating point precision.', WARNING)
  endif
  if ((.not. grid_is_latlon) .and. (Lx.eq.360.)) then
    if (mpp_pe().eq.mpp_root_pe())  then
            call error_mesg('particles, framework', 'Since the lat/lon grid is off, the x-direction is being set as non-periodic. Set Lx not equal to 360 override.', WARNING)
    endif
    Lx=-1.
  endif


 !The fix to reproduce across PE layout change, from AJA
  if (Lx>0.) then
    j=grd%jsc; do i=grd%isc+1,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i-1,j),Lx)
      if (abs(grd%lon(i,j)-lon_mod)>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo
    j=grd%jsc; do i=grd%isc-1,grd%isd,-1
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i+1,j) ,Lx)
      if (abs(grd%lon(i,j)-  lon_mod )>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo
    do j=grd%jsc+1,grd%jed; do i=grd%isd,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i,j-1) ,Lx)
      if (abs(grd%lon(i,j)-(lon_mod ))>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo; enddo
    do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i,j+1) ,Lx)
      if (abs(grd%lon(i,j)- lon_mod )>(Lx/2.)) &
        grd%lon(i,j)=  lon_mod
    enddo; enddo
  endif

  is=grd%isd; ie=grd%ied; js=grd%jsd; je=grd%jed
  grd%lon(is:ie,js:je)=Grid%geolonBu(is:ie,js:je)
  grd%lat(is:ie,js:je)=Grid%geolatBu(is:ie,js:je)

  ! WE SHOULD JUST COPY geolonC,geolatC instead (MJH)
  ! lonc, latc used for searches
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
    grd%lonc(i,j)=0.25*( (grd%lon(i,j)+grd%lon(i-1,j-1)) &
                        +(grd%lon(i-1,j)+grd%lon(i,j-1)) )
    grd%latc(i,j)=0.25*( (grd%lat(i,j)+grd%lat(i-1,j-1)) &
                        +(grd%lat(i-1,j)+grd%lat(i,j-1)) )
  enddo; enddo

  if (debug) then
    write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'particles, particles_init: (',mpp_pe(),') [ij][se]c=', &
         grd%isc,grd%iec,grd%jsc,grd%jec, &
         ' [lon|lat][min|max]=', minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)
  endif
! Final check for NaN's in the latlon grid:
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
    if (grd%lat(i,j) .ne. grd%lat(i,j)) then
      write(stderrunit,*) 'Lat not defined properly', mpp_pe(),i,j,grd%lat(i,j)
      call error_mesg('particles,grid defining', 'Latitude contains NaNs', FATAL)
    endif
    if (grd%lon(i,j) .ne. grd%lon(i,j)) then
      write(stderrunit,*) 'Lon not defined properly', mpp_pe(),i,j,grd%lon(i,j)
      call error_mesg('particles, grid defining', 'Longatudes contains NaNs', FATAL)
    endif
  enddo; enddo

  if (ignore_traj) buffer_width_traj=0 ! If this is true, then all traj files should be ignored


 ! Parameters
  parts%dt=dt
  parts%traj_sample_hrs=traj_sample_hrs
  parts%traj_write_hrs=traj_write_hrs
  parts%save_short_traj=save_short_traj
  parts%ignore_traj=ignore_traj
  parts%initial_traj=initial_traj
  parts%verbose_hrs=verbose_hrs
  parts%grd%halo=halo
  parts%grd%Lx=Lx
  parts%grd%no_TS=no_TS
  parts%grd%grid_is_latlon=grid_is_latlon
  parts%grd%grid_is_regular=grid_is_regular
  parts%Runge_not_Verlet=Runge_not_Verlet
  parts%ignore_missing_restart_parts=ignore_missing_restart_parts
  parts%use_new_predictive_corrective=use_new_predictive_corrective  !Alon
  parts%debug_particle_with_id=debug_particle_with_id



  if (do_unit_tests) then
   if (unit_tests(parts)) call error_mesg('particles, particles_init', 'Unit tests failed!', FATAL)
  endif

  ! Generate a grid of particles if requested
  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  if (3600*ihr + 60*imin +isec == 0) then ! Make sure we are on a day boundary
    lgenerate = .false.
    if (generate_days>0) then
      if ( mod(365*(iyr-1)+iday-1, generate_days)==0 ) lgenerate = .true.
    elseif (generate_days==0) then ! .and. 365*(iyr-1)+iday-1==0) then
      lgenerate = .true.
    endif
    if (lgenerate) then 
       call generate_grid_of_particles(parts, &
                          generate_lons(1), generate_lons(2), generate_lons(3), &
                          generate_lats(1), generate_lats(2), generate_lats(3), &
                          generate_d(1), generate_d(2), generate_d(3))
    endif
  endif

  call mpp_clock_end(parts%clock_ini)
  call mpp_clock_end(parts%clock)

end subroutine particles_framework_init

!> Generate particles on a grid
subroutine generate_grid_of_particles(parts, lon_start, lon_end, dlon, lat_start, lat_end, dlat,d_start,d_end,dd)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
real, intent(in) :: lon_start !< Start longitude of grid of particles
real, intent(in) :: lon_end !< End longitude of grid of particles
real, intent(in) :: dlon !< Separation longitude of particles on grid
real, intent(in) :: lat_start !< Start latitude of grid of particles
real, intent(in) :: lat_end !< End latitude of grid of particles
real, intent(in) :: dlat !< Separation latitude of particles on grid
real, intent(in) :: d_start !< Start depth of grid of particles
real, intent(in) :: d_end !< End depth of grid of particles
real, intent(in) :: dd !< Separation depth of particles on grid

! Local variables
type(particles_gridded), pointer :: grd => null()
type(particle) :: localpart
integer :: i, j, d, ie, je, de
integer :: num
real :: lat_min, lat_max
logical :: lres


  grd=>parts%grd
  lat_min = minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat_max = maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )

  ! Adjust local grid range to match generating grid
  lat_min = max( int((lat_min - lat_start)/dlat)*dlat + lat_start, lat_start)
  lat_max = min( int((lat_max - lat_start)/dlat)*dlat + lat_start, lat_end)

  ie = int( (lon_end-lon_start)/dlon - 0.5 )
  je = int( (lat_end-lat_start)/dlat - 0.5 )
  de = int( (abs(d_end)-abs(d_start))/abs(dd) - 0.5 )
  !CSJ hard-coded this for now
  !localpart%k_space=.true.
  !localpart%k=0.5
  num = 0
  do d = 0,de
     if (d_end<0) then
        localpart%k_space=.true.
        localpart%k=abs(d_start) + abs(dd*float(d))
     elseif (d_end>0) then
        localpart%k_space=.false.
        localpart%depth=abs(d_start) + abs(dd*float(d))
     endif
     do j =  0,je
        localpart%lat = lat_start + dlat*float(j)
!        if (localpart%lat >= lat_min .and. localpart%lat <= lat_max) then
         do i = 0,ie
           localpart%lon = lon_start + dlon*float(i)
           lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
           num=num+1
           if (lres) then
             if (grd%msk(localpart%ine,localpart%jne)>-1.) then
               localpart%drifter_num = num
               localpart%id = generate_id(grd, localpart%ine, localpart%jne)
               lres=pos_within_cell(grd, localpart%lon, localpart%lat,localpart%ine,localpart%jne, localpart%xi, localpart%yj)
               call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart)
             endif
           endif
         enddo
!       endif
     enddo
   enddo

  call parts_chksum(parts, 'after generated particles')

end subroutine generate_grid_of_particles

! ##############################################################################

!> Adjust part dates to allow use of restarts from later dates
subroutine offset_part_dates(parts,Time)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
type(time_type), intent(in) :: Time !< Model time
! Local variables
type(particle), pointer :: this
integer :: iyr, imon, iday, ihr, imin, isec, yr_offset
real :: latest_start_year, part_start_year
real :: current_time_val
integer :: grdi, grdj

  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  latest_start_year=iyr-999999.

  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      part_start_year=float(this%start_year)+this%start_day/367.
      if (part_start_year>latest_start_year) latest_start_year=part_start_year
      this=>this%next
    enddo
  enddo ; enddo
  call mpp_max(latest_start_year)

  current_time_val=float(iyr)+yearday(iyr, imon, iday, ihr, imin, isec)/367.
  if (latest_start_year<=current_time_val) return ! No conflicts!

  yr_offset=int(latest_start_year+1.)-iyr
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,i8,a)') &
    'particles: parts found with creation dates after model date! Adjusting part dates by ',yr_offset,' years'
  call parts_chksum(parts, 'before adjusting start dates')
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      this%start_year=this%start_year-yr_offset
      this=>this%next
    enddo
  enddo ; enddo
  call parts_chksum(parts, 'after adjusting start dates')

end subroutine offset_part_dates

! ###############################################################################################
!> Moves particles between lists if they have moved from cell to cell
subroutine move_part_between_cells(parts)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
! Local variables
type(particles_gridded), pointer :: grd => null()
type(particle), pointer :: moving_part => null(), this => null()
integer :: grdi, grdj
logical :: quick
! For convenience
grd=>parts%grd

do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))

      if ((this%ine.ne.grdi) .or. (this%jne.ne.grdj))  then
        moving_part=>this
        this=>this%next

        !Removing the particle from the old list
        if (associated(moving_part%prev)) then
          moving_part%prev%next=>moving_part%next
        else
          parts%list(grdi,grdj)%first=>moving_part%next
        endif
        if (associated(moving_part%next)) moving_part%next%prev=>moving_part%prev

        !Inserting the particle into the new list
        call insert_part_into_list(parts%list(moving_part%ine,moving_part%jne)%first,moving_part)

        !Clear moving_part
        moving_part=>null()

      else
        this=>this%next
      endif
    enddo
enddo ; enddo

end subroutine move_part_between_cells

! #############################################################################
!> Populates the halo lists with parts from neighbor processers
subroutine update_halo_particles(parts)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
! Local variables
type(particle), pointer :: kick_the_bucket, this
integer :: nparts_to_send_e, nparts_to_send_w
integer :: nparts_to_send_n, nparts_to_send_s
integer :: nparts_rcvd_from_e, nparts_rcvd_from_w
integer :: nparts_rcvd_from_n, nparts_rcvd_from_s
type(particles_gridded), pointer :: grd
integer :: i, nparts_start, nparts_end
integer :: stderrunit
integer :: grdi, grdj
integer :: halo_width
integer :: temp1, temp2
real :: current_halo_status
logical :: halo_debugging

  halo_width=parts%grd%halo
  halo_debugging=parts%halo_debugging

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  ! For debugging, MP1
  if (halo_debugging) then
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>parts%list(grdi,grdj)%first
      do while (associated(this))
          write(stderrunit,*) 'A', this%id, mpp_pe(), this%halo_part, grdi, grdj
        this=>this%next
      enddo
    enddo; enddo
    ! Use when debugging:
  endif

  ! Step 1: Clear the current halos
  call mpp_sync_self()
  do grdj = grd%jsd,grd%jsc-1 ;  do grdi = grd%isd,grd%ied
    call delete_all_parts_in_list(parts, grdj, grdi)
  enddo ; enddo

  do grdj = grd%jec+1,grd%jed ;  do grdi = grd%isd,grd%ied
    call delete_all_parts_in_list(parts,grdj,grdi)
  enddo ; enddo

  do grdj = grd%jsd,grd%jed ;    do grdi = grd%isd,grd%isc-1
    call delete_all_parts_in_list(parts,grdj,grdi)
  enddo ; enddo

  do grdj = grd%jsd,grd%jed ;    do grdi = grd%iec+1,grd%ied
    call delete_all_parts_in_list(parts,grdj,grdi)
  enddo ; enddo

  call mpp_sync_self()

  ! For debugging
  if (halo_debugging) then
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>parts%list(grdi,grdj)%first
        do while (associated(this))
        write(stderrunit,*) 'B', this%id, mpp_pe(), this%halo_part, grdi, grdj
      this=>this%next
      enddo
    enddo; enddo
  endif
  if (debug) then
    nparts_start=count_parts(parts)
  endif

  call mpp_sync_self()

  ! Step 2: Updating the halos  - This code is mostly copied from send_to_other_pes

  ! Find number of parts that headed east/west
  nparts_to_send_e=0
  nparts_to_send_w=0
  ! parts on eastern side of the processor
  do grdj = grd%jsc,grd%jec ; do grdi = grd%iec-halo_width+2,grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nparts_to_send_e=nparts_to_send_e+1
      current_halo_status=kick_the_bucket%halo_part
      kick_the_bucket%halo_part=1.
      call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_e, nparts_to_send_e)
      kick_the_bucket%halo_part=current_halo_status
    enddo
  enddo; enddo

  ! parts on the western side of the processor
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%isc+halo_width-1
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nparts_to_send_w=nparts_to_send_w+1
      current_halo_status=kick_the_bucket%halo_part
      kick_the_bucket%halo_part=1.
      call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_w, nparts_to_send_w)
      kick_the_bucket%halo_part=current_halo_status
    enddo
  enddo; enddo

  ! Send parts east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nparts_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nparts_to_send_e.gt.0) then
      call mpp_send(parts%obuffer_e%data, nparts_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send parts west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nparts_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nparts_to_send_w.gt.0) then
      call mpp_send(parts%obuffer_w%data, nparts_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive parts from west
  if (grd%pe_W.ne.NULL_PE) then
    nparts_rcvd_from_w=-999
    call mpp_recv(nparts_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nparts_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_w.gt.0) then
      call increase_ibuffer(parts%ibuffer_w, nparts_rcvd_from_w,buffer_width)
      call mpp_recv(parts%ibuffer_w%data, nparts_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nparts_rcvd_from_w
        call unpack_part_from_buffer2(parts, parts%ibuffer_w, i, grd )
      enddo
    endif
  else
    nparts_rcvd_from_w=0
  endif

  ! Receive parts from east
  if (grd%pe_E.ne.NULL_PE) then
    nparts_rcvd_from_e=-999
    call mpp_recv(nparts_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nparts_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_e.gt.0) then
      call increase_ibuffer(parts%ibuffer_e, nparts_rcvd_from_e,buffer_width)
      call mpp_recv(parts%ibuffer_e%data, nparts_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nparts_rcvd_from_e
        call unpack_part_from_buffer2(parts, parts%ibuffer_e, i, grd )
      enddo
    endif
  else
    nparts_rcvd_from_e=0
  endif

  ! Find number of parts that headed north/south
  nparts_to_send_n=0
  nparts_to_send_s=0

  ! parts on north side of the processor
  do grdj = grd%jec-halo_width+2,grd%jec ; do grdi = grd%isd,grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nparts_to_send_n=nparts_to_send_n+1
      current_halo_status=kick_the_bucket%halo_part
      kick_the_bucket%halo_part=1.
      call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_n, nparts_to_send_n )
      kick_the_bucket%halo_part=current_halo_status
    enddo
  enddo; enddo

  ! parts on south side of the processor
  do grdj = grd%jsc,grd%jsc+halo_width-1 ; do grdi = grd%isd,grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nparts_to_send_s=nparts_to_send_s+1
      current_halo_status=kick_the_bucket%halo_part
      kick_the_bucket%halo_part=1.
      call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_s, nparts_to_send_s )
      kick_the_bucket%halo_part=current_halo_status
    enddo
  enddo; enddo

 ! Send parts north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
      call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else
      call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nparts_to_send_n.gt.0) then
      if(folded_north_on_pe) then
        call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
        call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
      endif
    endif
  endif

  ! Send parts south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nparts_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nparts_to_send_s.gt.0) then
      call mpp_send(parts%obuffer_s%data, nparts_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif

  ! Receive parts from south
  if (grd%pe_S.ne.NULL_PE) then
    nparts_rcvd_from_s=-999
    call mpp_recv(nparts_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nparts_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_s.gt.0) then
      call increase_ibuffer(parts%ibuffer_s, nparts_rcvd_from_s,buffer_width)
      call mpp_recv(parts%ibuffer_s%data, nparts_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nparts_rcvd_from_s
        call unpack_part_from_buffer2(parts, parts%ibuffer_s, i, grd  )
      enddo
    endif
  else
    nparts_rcvd_from_s=0
  endif

  ! Receive parts from north
  if (grd%pe_N.ne.NULL_PE) then
    nparts_rcvd_from_n=-999
    if(folded_north_on_pe) then
      call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
      call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nparts_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_n.gt.0) then
      call increase_ibuffer(parts%ibuffer_n, nparts_rcvd_from_n,buffer_width)
      if(folded_north_on_pe) then
        call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
        call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nparts_rcvd_from_n
        call unpack_part_from_buffer2(parts, parts%ibuffer_n, i, grd )
      enddo
    endif
  else
    nparts_rcvd_from_n=0
  endif

  ! For debugging
  if (halo_debugging) then
    call mpp_sync_self()
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>parts%list(grdi,grdj)%first
      do while (associated(this))
        write(stderrunit,*)  'C', this%id, mpp_pe(), this%halo_part,  grdi, grdj
        this=>this%next
      enddo
    enddo; enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (debug) then
    nparts_end=count_parts(parts)
    i=nparts_rcvd_from_n+nparts_rcvd_from_s+nparts_rcvd_from_e+nparts_rcvd_from_w &
     -nparts_to_send_n-nparts_to_send_s-nparts_to_send_e-nparts_to_send_w
    if (nparts_end-(nparts_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_end=',nparts_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_start=',nparts_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: error=',nparts_end-(nparts_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_to_send_n=',nparts_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_to_send_s=',nparts_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_to_send_e=',nparts_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_to_send_w=',nparts_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_rcvd_from_n=',nparts_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_rcvd_from_s=',nparts_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_rcvd_from_e=',nparts_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, update_halos: nparts_rcvd_from_w=',nparts_rcvd_from_w,' on PE',mpp_pe()
    endif
  endif
  if (debug) then
    i=0
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      this=>parts%list(grdi,grdj)%first
      do while (associated(this))
        call check_position(grd, this, 'exchange (bot)')
        if (this%ine.lt.parts%grd%isc .or. &
            this%ine.gt.parts%grd%iec .or. &
            this%jne.lt.parts%grd%jsc .or. &
            this%jne.gt.parts%grd%jec) i=i+1
        this=>this%next
      enddo ! while
    enddo ; enddo
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'particles, update_halos: # of parts outside computational domain = ',i
      call error_mesg('particles, update_halos:', 'there are parts still in halos!', FATAL)
    endif ! root_pe
  endif ! debug
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

end subroutine update_halo_particles


! #############################################################################
!> Destroys all parts in a list
subroutine delete_all_parts_in_list(parts, grdj, grdi)
  type(particles), pointer :: parts !< Container for all types and memory
  integer :: grdi !< i-index of list
  integer :: grdj !< j-index of list
  ! Local variables
  type(particle), pointer :: kick_the_bucket, this
  this=>parts%list(grdi,grdj)%first
  do while (associated(this))
    kick_the_bucket=>this
    this=>this%next
    call destroy_particle(kick_the_bucket)
  enddo
  parts%list(grdi,grdj)%first=>null()
end  subroutine delete_all_parts_in_list


! #############################################################################
!> Send parts in halo lists to other processors
subroutine send_parts_to_other_pes(parts)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
! Local variables
type(particle), pointer :: kick_the_bucket, this
integer :: nparts_to_send_e, nparts_to_send_w
integer :: nparts_to_send_n, nparts_to_send_s
integer :: nparts_rcvd_from_e, nparts_rcvd_from_w
integer :: nparts_rcvd_from_n, nparts_rcvd_from_s
type(particles_gridded), pointer :: grd
integer :: i, nparts_start, nparts_end
integer :: stderrunit
integer :: grdi, grdj

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  if (debug) then
    nparts_start=count_parts(parts, with_halos=.true.)
  endif

  ! Find number of parts that headed east/west
  nparts_to_send_e=0
  nparts_to_send_w=0
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      if (this%halo_part .lt. 0.5) then
        if (this%ine.gt.parts%grd%iec) then
          kick_the_bucket=>this
          this=>this%next
          nparts_to_send_e=nparts_to_send_e+1
          call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_e, nparts_to_send_e  )
          call move_trajectory(parts, kick_the_bucket)
          call delete_particle_from_list(parts%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%ine.lt.parts%grd%isc) then
          kick_the_bucket=>this
          this=>this%next
          nparts_to_send_w=nparts_to_send_w+1
          call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_w, nparts_to_send_w  )
          call move_trajectory(parts, kick_the_bucket)
          call delete_particle_from_list(parts%list(grdi,grdj)%first,kick_the_bucket)
        else
          this=>this%next
        endif
      else
         this=>this%next
      endif
    enddo
  enddo ; enddo

  ! Send parts east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nparts_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nparts_to_send_e.gt.0) then
      call mpp_send(parts%obuffer_e%data, nparts_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send parts west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nparts_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nparts_to_send_w.gt.0) then
      call mpp_send(parts%obuffer_w%data, nparts_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive parts from west
  if (grd%pe_W.ne.NULL_PE) then
    nparts_rcvd_from_w=-999
    call mpp_recv(nparts_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nparts_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_w.gt.0) then
      call increase_ibuffer(parts%ibuffer_w, nparts_rcvd_from_w,buffer_width)
      call mpp_recv(parts%ibuffer_w%data, nparts_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nparts_rcvd_from_w
        call unpack_part_from_buffer2(parts, parts%ibuffer_w, i, grd )
      enddo
    endif
  else
    nparts_rcvd_from_w=0
  endif

  ! Receive parts from east
  if (grd%pe_E.ne.NULL_PE) then
    nparts_rcvd_from_e=-999
    call mpp_recv(nparts_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nparts_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_e.gt.0) then
      call increase_ibuffer(parts%ibuffer_e, nparts_rcvd_from_e,buffer_width)
      call mpp_recv(parts%ibuffer_e%data, nparts_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nparts_rcvd_from_e
        call unpack_part_from_buffer2(parts, parts%ibuffer_e, i, grd)
      enddo
    endif
  else
    nparts_rcvd_from_e=0
  endif

  ! Find number of parts that headed north/south
  ! (note: this block should technically go ahead of the E/W recv block above
  !  to handle arbitrary orientation of PEs. But for simplicity, it is
  !  here to accomodate diagonal transfer of parts between PEs -AJA)
  nparts_to_send_n=0
  nparts_to_send_s=0
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      if (this%halo_part .lt. 0.5) then
        if (this%jne.gt.parts%grd%jec) then
          kick_the_bucket=>this
          this=>this%next
          nparts_to_send_n=nparts_to_send_n+1
          call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_n, nparts_to_send_n)
          call move_trajectory(parts, kick_the_bucket)
          call delete_particle_from_list(parts%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%jne.lt.parts%grd%jsc) then
          kick_the_bucket=>this
          this=>this%next
          nparts_to_send_s=nparts_to_send_s+1
          call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_s, nparts_to_send_s)
          call move_trajectory(parts, kick_the_bucket)
          call delete_particle_from_list(parts%list(grdi,grdj)%first,kick_the_bucket)
        else
          this=>this%next
        endif
      else
          this=>this%next
      endif
    enddo
  enddo ; enddo

  ! Send parts north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
       call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nparts_to_send_n.gt.0) then
       if(folded_north_on_pe) then
          call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
       else
          call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
       endif
    endif
  endif

  ! Send parts south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nparts_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nparts_to_send_s.gt.0) then
      call mpp_send(parts%obuffer_s%data, nparts_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif

  ! Receive parts from south
  if (grd%pe_S.ne.NULL_PE) then
    nparts_rcvd_from_s=-999
    call mpp_recv(nparts_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nparts_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_s.gt.0) then
      call increase_ibuffer(parts%ibuffer_s, nparts_rcvd_from_s,buffer_width)
      call mpp_recv(parts%ibuffer_s%data, nparts_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nparts_rcvd_from_s
        call unpack_part_from_buffer2(parts, parts%ibuffer_s, i, grd )
      enddo
    endif
  else
    nparts_rcvd_from_s=0
  endif

  ! Receive parts from north
  if (grd%pe_N.ne.NULL_PE) then
    nparts_rcvd_from_n=-999
    if(folded_north_on_pe) then
       call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nparts_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_n.gt.0) then
      call increase_ibuffer(parts%ibuffer_n, nparts_rcvd_from_n,buffer_width)
      if(folded_north_on_pe) then
         call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
         call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nparts_rcvd_from_n
        call unpack_part_from_buffer2(parts, parts%ibuffer_n, i, grd)
      enddo
    endif
  else
    nparts_rcvd_from_n=0
  endif

  if (debug) then
    nparts_end=count_parts(parts, with_halos=.true.)
    i=nparts_rcvd_from_n+nparts_rcvd_from_s+nparts_rcvd_from_e+nparts_rcvd_from_w &
     -nparts_to_send_n-nparts_to_send_s-nparts_to_send_e-nparts_to_send_w
    if (nparts_end-(nparts_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_end=',nparts_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_start=',nparts_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: error=',nparts_end-(nparts_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_to_send_n=',nparts_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_to_send_s=',nparts_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_to_send_e=',nparts_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_to_send_w=',nparts_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_rcvd_from_n=',nparts_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_rcvd_from_s=',nparts_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_rcvd_from_e=',nparts_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'particles, send_parts_to_other_pes: nparts_rcvd_from_w=',nparts_rcvd_from_w,' on PE',mpp_pe()
      call error_mesg('particles, send_parts_to_other_pes:', 'We lost some parts!', FATAL)
    endif
  endif

  if (debug) then
    i=0
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      this=>parts%list(grdi,grdj)%first
      do while (associated(this))
        call check_position(grd, this, 'exchange (bot)', grdi, grdj)
        if (this%ine.lt.parts%grd%isc .or. &
            this%ine.gt.parts%grd%iec .or. &
            this%jne.lt.parts%grd%jsc .or. &
            this%jne.gt.parts%grd%jec) i=i+1
        this=>this%next
      enddo ! while
    enddo ; enddo
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'particles, send_parts_to_other_pes: # of parts outside computational domain = ',i
      call error_mesg('particles, send_parts_to_other_pes:', 'there are parts still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

  call mpp_sync_self()

end subroutine send_parts_to_other_pes

!> Pack a part into a buffer
subroutine pack_part_into_buffer2(part, buff, n)
! Arguments
type(particle), pointer :: part !< particle to pack into buffer
type(buffer), pointer :: buff !< Buffer to pack part into
integer, intent(in) :: n !< Position in buffer to place part
!integer, intent(in) :: max_bonds  ! Change this later
! Local variables
integer :: counter, k, max_bonds, id_cnt, id_ij
integer :: kspace_int

  max_bonds=0

  if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width)
  if (n>buff%size) call increase_ibuffer(buff,n,buffer_width)


 if (part%k_space) then
   kspace_int=1
 else
   kspace_int=0
 endif


  counter = 0
  call push_buffer_value(buff%data(:,n), counter, part%lon)
  call push_buffer_value(buff%data(:,n), counter, part%lat)
  call push_buffer_value(buff%data(:,n), counter, part%depth)
  call push_buffer_value(buff%data(:,n), counter, part%uvel)
  call push_buffer_value(buff%data(:,n), counter, part%vvel)
  call push_buffer_value(buff%data(:,n), counter, part%xi)
  call push_buffer_value(buff%data(:,n), counter, part%yj)
  call push_buffer_value(buff%data(:,n), counter, part%theta)
  call push_buffer_value(buff%data(:,n), counter, part%start_lon)
  call push_buffer_value(buff%data(:,n), counter, part%start_lat)
  call push_buffer_value(buff%data(:,n), counter, part%start_d)
  call push_buffer_value(buff%data(:,n), counter, part%start_year)
  call push_buffer_value(buff%data(:,n), counter, part%start_day)
  call push_buffer_value(buff%data(:,n), counter, INT(part%id))
  call push_buffer_value(buff%data(:,n), counter, INT(part%drifter_num))
  call push_buffer_value(buff%data(:,n), counter, part%ine)
  call push_buffer_value(buff%data(:,n), counter, part%jne)
  call push_buffer_value(buff%data(:,n), counter, part%k)
  call push_buffer_value(buff%data(:,n), counter, kspace_int)
  call push_buffer_value(buff%data(:,n), counter, part%halo_part)



end subroutine pack_part_into_buffer2

!> Set a real value in the buffer at position (counter) after incrementing counter
subroutine push_buffer_rvalue(vbuf, counter, val)
  real, dimension(:), intent(inout) :: vbuf    !< Buffer vector to push value into
  integer,            intent(inout) :: counter !< Position to increment
  real,               intent(in)    :: val     !< Value to place in buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in push_buffer_rvalue'
  vbuf(counter) = val

end subroutine push_buffer_rvalue

!> Set an integer value in the buffer at position (counter) after incrementing counter
subroutine push_buffer_ivalue(vbuf, counter, val)
  real, dimension(:), intent(inout) :: vbuf    !< Buffer vector to push value into
  integer,            intent(inout) :: counter !< Position to increment
  integer,            intent(in)    :: val     !< Value to place in buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in push_buffer_ivalue'
  vbuf(counter) = float(val)

end subroutine push_buffer_ivalue

!> Get a real value from the buffer at position (counter) after incrementing counter
subroutine pull_buffer_rvalue(vbuf, counter, val)
  real, dimension(:), intent(in)    :: vbuf    !< Buffer vector to pull value from
  integer,            intent(inout) :: counter !< Position to increment
  real,               intent(out)   :: val     !< Value to get from buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in pull_buffer_rvalue'
  val = vbuf(counter)

end subroutine pull_buffer_rvalue

!> Get an integer value from the buffer at position (counter) after incrementing counter
subroutine pull_buffer_ivalue(vbuf, counter, val)
  real, dimension(:), intent(in)    :: vbuf    !< Buffer vector to pull value from
  integer,            intent(inout) :: counter !< Position to increment
  integer,            intent(out)   :: val     !< Value to get from buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in pull_buffer_ivalue'
  val = nint(vbuf(counter))

end subroutine pull_buffer_ivalue

  subroutine increase_buffer(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=delta
    else
      new_size=old%size+delta
    endif
    allocate(new)
    allocate(new%data(buffer_width,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new
   !write(stderr(),*) 'particles, increase_buffer',mpp_pe(),' increased to',new_size

  end subroutine increase_buffer

  !> Unpacks a part entry from a buffer to a new part
subroutine unpack_part_from_buffer2(parts, buff, n, grd, force_append)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
type(buffer), pointer :: buff !< Buffer from which to unpack part
integer, intent(in) :: n !< Position in buffer to unpack
type(particles_gridded), pointer :: grd !< Container for gridded fields
logical, optional :: force_append !< <undocumented>

! Local variables
!real :: lon, lat, uvel, vvel, xi, yj
!real :: start_lon, start_lat, start_day, start_mass
!integer :: ine, jne, start_year
logical :: lres
type(particle) :: localpart
type(particle), pointer :: this
integer :: other_part_ine, other_part_jne
integer :: counter, k, max_bonds, id_cnt, id_ij,tmp
integer(kind=8) :: id
integer :: kspace_int
integer :: stderrunit
logical :: force_app
logical :: quick

  ! Get the stderr unit number
  stderrunit = stderr()

  quick=.false.
  max_bonds=0

  force_app = .false.
  if(present(force_append)) force_app = force_append



  counter = 0
  call pull_buffer_value(buff%data(:,n), counter, localpart%lon)
  call pull_buffer_value(buff%data(:,n), counter, localpart%lat)
  call pull_buffer_value(buff%data(:,n), counter, localpart%depth)
  call pull_buffer_value(buff%data(:,n), counter, localpart%uvel)
  call pull_buffer_value(buff%data(:,n), counter, localpart%vvel)
  call pull_buffer_value(buff%data(:,n), counter, localpart%xi)
  call pull_buffer_value(buff%data(:,n), counter, localpart%yj)
  call pull_buffer_value(buff%data(:,n), counter, localpart%theta)
  call pull_buffer_value(buff%data(:,n), counter, localpart%start_lon)
  call pull_buffer_value(buff%data(:,n), counter, localpart%start_lat)
  call pull_buffer_value(buff%data(:,n), counter, localpart%start_d)
  call pull_buffer_value(buff%data(:,n), counter, localpart%start_year)
  call pull_buffer_value(buff%data(:,n), counter, localpart%start_day)
  call pull_buffer_value(buff%data(:,n), counter, tmp)
  localpart%id=INT(tmp,kind=8)
  call pull_buffer_value(buff%data(:,n), counter, tmp)
  localpart%drifter_num=INT(tmp,kind=8)
  call pull_buffer_value(buff%data(:,n), counter, localpart%ine)
  call pull_buffer_value(buff%data(:,n), counter, localpart%jne)
  call pull_buffer_value(buff%data(:,n), counter, localpart%k)
  call pull_buffer_value(buff%data(:,n), counter, kspace_int)
  call pull_buffer_value(buff%data(:,n), counter, localpart%halo_part)

 if (kspace_int>0.5) then
   localpart%k_space=.true.
 else
   localpart%k_space=.false.
 endif


  !These quantities no longer need to be passed between processors
  localpart%uvel_old=localpart%uvel
  localpart%vvel_old=localpart%vvel
  localpart%lon_old=localpart%lon
  localpart%lat_old=localpart%lat

  !added by CSJ in order to get a re-entrant channel (more permanent solution needed)
!  if (localpart%lon>22.0) then
       ! write(stderrunit,'("moving part=(",i3,a,2f8.2)')&
       !  & mpp_pe(),' from lon,lat=',localpart%lon,localpart%lat
!       localpart%lon=localpart%lon-22.0
       ! write(stderrunit,*) localpart%lon,localpart%lat
!  elseif (localpart%lon<0.0) then
!       localpart%lon=localpart%lon+22.0
!  endif


  ! force_app=.true.
  if(force_app) then !force append with origin ine,jne (for I/O)
    call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart, quick, this)
  else
    lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    if (lres) then
      lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
      call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart,quick,this)
    else
      lres=find_cell_wide(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
      if (lres) then
        lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
        call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart,quick,this)
      else
        write(stderrunit,'("particles, unpack_part_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
         & mpp_pe(),') Failed to find i,j=',localpart%ine,localpart%jne,' for lon,lat=',localpart%lon,localpart%lat
        write(stderrunit,*) localpart%lon,localpart%lat
        write(stderrunit,*) localpart%uvel,localpart%vvel
        write(stderrunit,*) localpart%uvel_old,localpart%vvel_old
        write(stderrunit,*) localpart%lon_old,localpart%lat_old
        write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
        write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
        write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
        write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
        write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
        write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
        write(stderrunit,*) lres
        call error_mesg('particles, unpack_part_from_buffer', 'can not find a cell to place part in!', FATAL)
      endif
    endif
  endif

  this=>null()

end subroutine unpack_part_from_buffer2

!> Increase size of buffer
!!
!! This routine checks if the buffer size is smaller than nparts
!! If it is, the buffer size is increased by delta_buf.
!! The buffer increases by more than 1 so that the buffer does not have to increase every time.
subroutine increase_ibuffer(old, num_parts, width)
! Arguments
type(buffer), pointer :: old !< Buffer to expand
integer, intent(in) :: num_parts !< Number of parts
integer, intent(in) :: width !< Width of buffer (first dimension)
! Local variables
type(buffer), pointer :: new
integer :: new_size, old_size

  if (.not.associated(old)) then
    new_size=num_parts+delta_buf
    old_size=0
  else
    old_size=old%size
    if (num_parts<old%size) then
      new_size=old%size
    else
      new_size=num_parts+delta_buf
    endif
  endif

  if (old_size.ne.new_size) then
    allocate(new)
    allocate(new%data(width,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new
  endif

end subroutine increase_ibuffer

  subroutine increase_ibuffer_traj(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size, old_size

    if (.not.associated(old)) then
      new_size=delta+delta_buf
      old_size=0
    else
      old_size=old%size
      if (delta<old%size) then
        new_size=old%size+delta
      else
        new_size=delta+delta_buf
      endif
    endif

    if (old_size.ne.new_size) then
      allocate(new)
      allocate(new%data(buffer_width_traj,new_size))
      new%size=new_size
      if (associated(old)) then
        new%data(:,1:old%size)=old%data(:,1:old%size)
        deallocate(old%data)
        deallocate(old)
      endif
      old=>new
    endif

  end subroutine increase_ibuffer_traj

  subroutine increase_buffer_traj(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=delta
    else
      new_size=old%size+delta
    endif
    allocate(new)
    allocate(new%data(buffer_width_traj,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new

  end subroutine increase_buffer_traj

  subroutine pack_traj_into_buffer2(traj, buff, n)
  ! Arguments
  type(xyt), pointer :: traj
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables
  integer :: cnt, ij

    if (.not.associated(buff)) call increase_buffer_traj(buff,delta_buf)
    if (n>buff%size) call increase_buffer_traj(buff,delta_buf)

    buff%data(1,n)=traj%lon
    buff%data(2,n)=traj%lat
    buff%data(3,n)=traj%day
    buff%data(4,n)=traj%year
    call split_id(traj%id, cnt, ij)
    buff%data(5,n)=float(cnt)
    buff%data(6,n)=float(ij)
    buff%data(7,n)=traj%uvel
    buff%data(8,n)=traj%vvel
    buff%data(9,n)=traj%theta
    buff%data(10,n)=traj%uvel_old !Alon
    buff%data(11,n)=traj%vvel_old !Alon
    buff%data(12,n)=traj%lon_old !Alon
    buff%data(13,n)=traj%lat_old !Alon
    buff%data(14,n)=traj%particle_num !Alon
    buff%data(15,n)=traj%k
    buff%data(16,n)=traj%k_space
    buff%data(17,n)=traj%depth

  end subroutine pack_traj_into_buffer2

  subroutine unpack_traj_from_buffer2(first, buff, n)
  ! Arguments
  type(xyt), pointer :: first
  type(buffer), pointer :: buff
  integer, intent(in) :: n
 ! Local variables
  type(xyt) :: traj
  integer :: stderrunit, cnt, ij
  ! Get the stderr unit number
  stderrunit = stderr()

    traj%lon=buff%data(1,n)
    traj%lat=buff%data(2,n)
!    traj%year=nint(buff%data(3,n))
    traj%day=buff%data(3,n)
    traj%year=buff%data(4,n) 
    cnt = nint(buff%data(5,n))
    ij = nint(buff%data(6,n))
    traj%id = id_from_2_ints(cnt, ij)
    traj%uvel=buff%data(7,n)
    traj%vvel=buff%data(8,n)
    traj%theta=buff%data(9,n)
    traj%uvel_old=buff%data(10,n) !Alon
    traj%vvel_old=buff%data(11,n) !Alon
    traj%lon_old=buff%data(12,n) !Alon
    traj%lat_old=buff%data(13,n) !Alon
    traj%particle_num=buff%data(14,n)
    traj%k=buff%data(15,n) 
    traj%k_space=buff%data(16,n)
    traj%depth=buff%data(17,n)

    call append_posn(first, traj)

  end subroutine unpack_traj_from_buffer2


! ##############################################################################
!> Add a new part to a list by copying values
!!
!! The input part are a part with set values whose memory is assumed to be
!! temporary. This routine allocates memory for a new part and copies the
!! the input values into it. The memory for the new part is pointed to
!! by newpart_return (if present).
subroutine add_new_part_to_list(first, partvals, quick, newpart_return)
! Arguments
type(particle), pointer :: first !< List of particles
type(particle), intent(in) :: partvals !< part values to copy
type(particle), intent(out), pointer, optional :: newpart_return !< New part
logical, intent(in), optional :: quick !< If true, use the quick insertion algorithm
! Local variables
type(particle), pointer :: new=>null()

  new=>null()
  call create_particle(new, partvals)

  if (present(newpart_return)) then
    newpart_return=>new
    !newpart_return=>null()
  endif

  if (present(quick)) then
    if(quick) then
      call insert_part_into_list(first, new, quick=.true.)
    else
      call insert_part_into_list(first, new)
    endif
  else
    call insert_part_into_list(first, new)
  endif

  !Clear new
  new=>null()

end subroutine add_new_part_to_list

! ##############################################################################
!> Scans all lists and checks that parts are in a sorted order in each list
subroutine count_out_of_order(parts,label)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
character(len=*) :: label !< Label to add to messages
! Local variables
type(particle), pointer :: this, next
integer :: i, icnt1, icnt2, icnt3
integer :: grdi, grdj

  icnt1=0; icnt3=0
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    next=>null()
    if (associated(this)) then
      if (associated(this%next)) next=>this%next
    endif
    do while (associated(next))
      if (.not. inorder(this,next)) icnt1=icnt1+1
      if (inorder(this,next).and.inorder(next,this)) icnt3=icnt3+1
      this=>next
      next=>next%next
    enddo
  enddo;enddo
  call mpp_sum(icnt1)

  i=0;
  icnt2=0
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      i=1
      if (this%ine<parts%grd%isc .or. &
          this%ine>parts%grd%iec .or. &
          this%jne<parts%grd%jsc .or. &
          this%jne>parts%grd%jec) icnt2=icnt2+1
      this=>this%next
      if (i>1.and..not.associated(this%prev)) then
        call error_mesg('particles, count_out_of_order', 'Pointer %prev is unassociated. This should not happen!', FATAL)
      endif
    enddo
  enddo; enddo
  call mpp_sum(icnt2)

  if ((debug.or.icnt1.ne.0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,3(x,a,i6),x,a)') 'particles, count_out_of_order:', &
      '# out of order=', icnt1,'# in halo=',icnt2,'# identicals=',icnt3,label
  endif

  call check_for_duplicates(parts,label)

end subroutine count_out_of_order


! ##############################################################################
!> Scans all lists and checks for duplicate identifiers between lists
subroutine check_for_duplicates(parts,label)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
character(len=*) :: label !< Label to add to message
! Local variables
type(particle), pointer :: this1, next1, this2, next2
integer :: icnt_id, icnt_same
integer :: grdi, grdj
integer :: grdi_inner, grdj_inner

  icnt_id=0
  icnt_same=0
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this1=>parts%list(grdi,grdj)%first
    do while (associated(this1))
      do grdj_inner = grdj,parts%grd%jec ; do grdi_inner = parts%grd%isc,parts%grd%iec
        if ( .not.(grdj_inner==grdj .and. grdi_inner<grdi) ) then
          if (grdj_inner==grdj .and. grdi_inner==grdi ) then
            this2=>this1%next
          else
            this2=>parts%list(grdi_inner,grdj_inner)%first
          endif
          do while (associated(this2))
            if (sameid(this1,this2)) icnt_id=icnt_id+1
            if (samepart(this1,this2)) icnt_same=icnt_same+1
            this2=>this2%next
          enddo
        endif
      enddo ; enddo
      this1=>this1%next
    enddo
  enddo ; enddo
  call mpp_sum(icnt_id)
  call mpp_sum(icnt_same)

  if ((debug.or.icnt_id>0.or.icnt_same>0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,2(x,a,i9),x,a)') 'particles, check_for_duplicates:', &
      '# with same id=', icnt_id,'# identical parts=',icnt_same,label
  endif

end subroutine check_for_duplicates

!> Generate an iceberg id from a counter at calving-location and the calving location itself.
!! Note that this updates grd%iceberg_counter_grd.
!!
!! \todo If we initialized grd%iceberg_counter_grd to 0 when the model is first run we could move
!! this increment line to before the id generation and then the counter would be an actual count.
integer(kind=8) function generate_id(grd, i, j)
  type(particles_gridded), pointer    :: grd !< Container for gridded fields
  integer,                intent(in) :: i   !< i-index of calving location
  integer,                intent(in) :: j   !< j-index of calving location
  ! Local variables
  integer :: ij ! Hash of i,j

  ! Increment counter in calving cell
  grd%particle_counter_grd(i,j) = grd%particle_counter_grd(i,j) + 1
  ! ij is unique number for each grid cell (32-bit integers allow for ~1/100th degree global resolution)
  ij = ij_component_of_id(grd, i, j)
  ! Generate a 64-bit id
  generate_id = id_from_2_ints( grd%particle_counter_grd(i,j), ij )

end function generate_id

!> Calculate the location-derived component of a particle id which is a hash of the i,j-indexes for the cell
integer function ij_component_of_id(grd, i, j)
  type(particles_gridded), pointer    :: grd !< Container for gridded fields
  integer,                intent(in) :: i   !< i-index of calving location
  integer,                intent(in) :: j   !< j-index of calving location
  ! Local variables
  integer :: ij ! Hash of i,j
  integer :: iNg ! Zonal size of the global grid
  integer :: ig ! Global i-index
  integer :: jg ! Global j-index

  ! Using the current grid shape maximizes the numbers of IDs that can be represented
  ! allowing up to 30-minute uniform global resolution, or potentially finer if non-uniform. 
  iNg = grd%ieg - grd%isg + 1

  ig = i + grd%is_offset ! Calculate global i-index
  jg = j + grd%js_offset ! Calculate global j-index

  ! ij_component_of_id is unique number for each grid cell (32-bit integers allow for ~1/100th degree global resolution)
  ij_component_of_id = ig + ( iNg * ( jg - 1 ) )

end function ij_component_of_id

! ##############################################################################
!> Prints a particular part's vitals
!!
!! All lists are scanned and if a part has the identifier equal to
!! debug_particle_with_id then the state of that part is printed.
subroutine monitor_a_part(parts, label)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
character(len=*) :: label !< Label to add to message
! Local variables
type(particle), pointer :: this
integer :: grdi, grdj
integer :: stderrunit

  if (parts%debug_particle_with_id<0) return
  stderrunit=stderr() ! Get the stderr unit number

  do grdj = parts%grd%jsd,parts%grd%jed ; do grdi = parts%grd%isd,parts%grd%ied
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      if (this%id == parts%debug_particle_with_id) then
        call print_part(stderrunit, this, 'MONITOR: '//label, grdi, grdj)
      endif
      this=>this%next
    enddo
  enddo ; enddo

end subroutine monitor_a_part

! ##############################################################################
!> Inserts a part into a list
subroutine insert_part_into_list(first, newpart, quick)
! Arguments
type(particle), pointer :: first !< List of parts
type(particle), pointer :: newpart !< New part to insert
logical, intent(in), optional :: quick !< If true, use the quick insertion algorithm
                                       !! \todo Delete arguments since the code does not appear to use it.
! Local variables
type(particle), pointer :: this, prev
logical :: quickly

  quickly = .false.

  if (associated(first)) then
    if (.not. parallel_reprod .or. quickly) then
      newpart%next=>first
      newpart%prev=>null()
      first%prev=>newpart
      first=>newpart
    else
      if (inorder(newpart,first)) then
        ! Insert at front of list
        newpart%next=>first
        newpart%prev=>null()
        first%prev=>newpart
        first=>newpart
      else
        this=>first
        prev=>null()
        do while( associated(this) )
          if (inorder(newpart,this) ) then
            exit
          endif
          prev=>this
          this=>this%next
        enddo
        prev%next=>newpart
        newpart%prev=>prev
        if (associated(this)) this%prev=>newpart
        newpart%next=>this
      endif
    endif
  else
    ! list is empty so create it
    first=>newpart
    first%next=>null()
    first%prev=>null()
  endif

end subroutine insert_part_into_list


! ##############################################################################

logical function inorder(part1, part2)  !MP Alon - Change to include particle_num
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  if (part1%start_year<part2%start_year) then ! want newer first
    inorder=.true.
    return
  else if (part1%start_year>part2%start_year) then
    inorder=.false.
    return
  endif
  if (part1%start_day<part2%start_day) then ! want newer first
    inorder=.true.
    return
  else if (part1%start_day>part2%start_day) then
    inorder=.false.
    return
  endif
  if (part1%start_lon<part2%start_lon) then ! want eastward first
    inorder=.true.
    return
  else if (part1%start_lon>part2%start_lon) then
    inorder=.false.
    return
  endif
  if (part1%start_lat<part2%start_lat) then ! want southern first
    inorder=.true.
    return
  else if (part1%start_lat>part2%start_lat) then
    inorder=.false.
    return
  endif
  if (part1%start_d<part2%start_d) then ! want top first
    inorder=.true.
  else if (part1%start_d>part2%start_d) then
    inorder=.false.
    return
  endif
  inorder=.true. ! passing the above tests mean the parts 1 and 2 are identical?
end function inorder

! ##############################################################################

  real function time_hash(part)!  Alon: Think about removing this.
  ! Arguments
  type(particle), pointer :: part
    time_hash=part%start_day+366.*float(part%start_year)
  end function time_hash

! ##############################################################################

  real function pos_hash(part)
  ! Arguments
  type(particle), pointer :: part
    pos_hash=part%start_lon+360.*(part%start_lat+90.)
  end function pos_hash

! ##############################################################################

logical function sameid(part1, part2) !  Alon: MP updat this.
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  sameid=.false.
  if (part1%start_year.ne.part2%start_year) return
  if (part1%start_day.ne.part2%start_day) return
  if (part1%start_lon.ne.part2%start_lon) return
  if (part1%start_lat.ne.part2%start_lat) return
  if (part1%start_d.ne.part2%start_d) return
  sameid=.true. ! passing the above tests means that parts 1 and 2 have the same id
end function sameid

! ##############################################################################

logical function samepart(part1, part2)
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  samepart=.false.
  if (.not. sameid(part1, part2)) return
  if (part1%lon.ne.part2%lon) return
  if (part1%lat.ne.part2%lat) return
  if (part1%uvel.ne.part2%uvel) return
  if (part1%vvel.ne.part2%vvel) return
  if (part1%uvel_old.ne.part2%uvel_old) return  !Alon
  if (part1%vvel_old.ne.part2%vvel_old) return  !Alon
  if (part1%lon_old.ne.part2%lon_old) return  !Alon
  if (part1%lat_old.ne.part2%lat_old) return  !Alon
  samepart=.true. ! passing the above tests mean that parts 1 and 2 are identical
end function samepart

! ##############################################################################

real function yearday(iyr, imon, iday, ihr, imin, isec)
! Arguments
integer, intent(in) :: iyr, imon, iday, ihr, imin, isec

! Local variables
integer :: sec
integer :: yeardayint
type(time_type) :: Time_test

   call get_time((set_date(iyr, imon, iday, ihr, imin, isec) - set_date(iyr, 1, 1, 0, 0, 0)), sec, yeardayint)
   yearday = float(yeardayint) + (float(ihr) + (float(imin) + float(isec)/60.)/60.)/24.

end function yearday

! ##############################################################################

subroutine create_particle(part, partvals)
! Arguments
type(particle), pointer :: part
type(particle), intent(in) :: partvals
! Local variables
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  ! write(stderrunit,'(a8,i3)') 'partvalk', (partvals%k)

  if (associated(part)) then
    write(stderrunit,*) 'particles, create_particle: part already associated!!!!',mpp_pe()
    call error_mesg('particles, create_particle', 'part already associated. This should not happen!', FATAL)
  endif
  allocate(part)
  part=partvals
  part%prev=>null()
  part%next=>null()

end subroutine create_particle


! ##############################################################################

subroutine delete_particle_from_list(first, part)
! Arguments
type(particle), pointer :: first, part
! Local variables

  ! Connect neighbors to each other
  if (associated(part%prev)) then
    part%prev%next=>part%next
  else
    first=>part%next
  endif
  if (associated(part%next)) part%next%prev=>part%prev

  ! Bye-bye part
  call destroy_particle(part)

end subroutine delete_particle_from_list

! ##############################################################################

subroutine destroy_particle(part)
! Arguments
type(particle), pointer :: part
! Local variables

  ! Bye-bye part
  deallocate(part)

end subroutine destroy_particle

! ##############################################################################
!> Print the state of a particular part
subroutine print_part(iochan, part, label, il, jl)
! Arguments
integer, intent(in) :: iochan !< Standard channel to use (usually stdout or stderr)
type(particle), pointer :: part !< part to print
character(len=*) :: label !< Label to use in messages
integer, optional, intent(in) :: il !< i-index of cell part should be in
integer, optional, intent(in) :: jl !< j-index of cell part should be in
! Local variables

  write(iochan,'("particles, print_part: ",2a,i5,a,i12,a,2f10.4,i5,f7.2,f5.1)') &
    label, 'pe=(', mpp_pe(), ') #=', part%id, ' start lon,lat,yr,day=', &
    part%start_lon, part%start_lat, part%start_year, part%start_day, part%halo_part
  if (present(il).and.present(jl)) then
    write(iochan,'("particles, print_part: ",2a,i5,a,i12,a,2i5)') &
      label, 'pe=(', mpp_pe(), ') #=', part%id, ' List i,j=',il,jl
  endif
  write(iochan,'("particles, print_part: ",2a,i5,a,i12,a,2i5,a,2l2)') &
    label, 'pe=(', mpp_pe(), ') #=', part%id, &
    ' i,j=', part%ine, part%jne, &
    ' p,n=', associated(part%prev), associated(part%next)
  write(iochan,'("particles, print_part: ",2a,i5,a,i12,3(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', part%id, &
    ' xi,yj=', part%xi, part%yj, &
    ' lon,lat=', part%lon, part%lat, &
    ' lon_old,lat_old=', part%lon_old, part%lat_old
  write(iochan,'("particles, print_part: ",2a,i5,a,i12,2(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', part%id, &
    ' u,v=', part%uvel, part%vvel, &
    ' uvel_old,vvel_old=', part%uvel_old, part%vvel_old
  write(iochan,'("particles, print_part: ",2a,i5,a,i12,(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', part%id, &
    ' depth,k=', part%depth, part%k
!  write(iochan,'("particles, print_part: ",2a,i5,a,i12,2(a,2f14.8))') &
!    label, 'pe=(', mpp_pe(), ') #=', part%id, &
!    ' axn,ayn=', part%axn, part%ayn, &
!    ' bxn,byn=', part%bxn, part%byn
end subroutine print_part


! ##############################################################################
!> Print the state of all parts
subroutine print_parts(iochan, parts, label)
! Arguments
integer, intent(in) :: iochan !< Standard channel to use (usually stdout or stderr)
type(particles), pointer :: parts !< Container for all types and memory
character(len=*) :: label !< Label to use in messages
! Local variables
integer :: nparts, nnparts
type(particle), pointer :: this
integer :: grdi, grdj

  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while(associated(this))
      call print_part(iochan, this, label)
      this=>this%next
    enddo
  enddo ; enddo
  nparts=count_parts(parts)
  nnparts=nparts
  call mpp_sum(nnparts)
  if (nparts.gt.0) write(iochan,'("particles, ",a," there are",i5," parts out of",i6," on PE ",i4)') label, nparts, nnparts, mpp_pe()

end subroutine print_parts


! ##############################################################################
!> Returns number of parts across all lists
integer function count_parts(parts, with_halos)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
logical, optional :: with_halos !< If true, include halo lists
! Local variables
integer :: grdi, grdj, is, ie, js, je
logical :: include_halos

  include_halos = .false.
  if (present(with_halos)) include_halos = with_halos
  if (include_halos) then
   is = parts%grd%isd ; ie = parts%grd%ied ; js = parts%grd%jsd ; je = parts%grd%jed
  else
   is = parts%grd%isc ; ie = parts%grd%iec ; js = parts%grd%jsc ; je = parts%grd%jec
  endif

  count_parts=0
  do grdj = js,je ; do grdi = is,ie
    count_parts=count_parts+count_parts_in_list(parts%list(grdi,grdj)%first)
  enddo ; enddo

end function count_parts


!> Returns number of parts in a list
integer function count_parts_in_list(first)
! Arguments
type(particle), pointer :: first !< List of parts
! Local variables
type(particle), pointer :: this

  count_parts_in_list=0
  this=>first
  do while(associated(this))
    count_parts_in_list=count_parts_in_list+1
    this=>this%next
  enddo

end function count_parts_in_list

! ##############################################################################
!> Add a record to the trajectory of each part
subroutine record_posn(parts, h, thetao)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
real,dimension(:,:,:), intent(in) :: h
real,dimension(:,:,:), optional, intent(in) :: thetao
! Local variables
type(xyt) :: posn
logical :: kspace_copy
type(particle), pointer :: this
integer :: grdi, grdj
type(particles_gridded), pointer :: grd
integer :: stderrunit

  stderrunit=stderr()

  grd=>parts%grd
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      posn%lon=this%lon
      posn%lat=this%lat
      posn%k=this%k
      kspace_copy = this%k_space
      call find_depth(grd,this%k,h,this%depth,this%ine,this%jne,this%xi,this%yj,kspace_copy)
      if (grd%no_TS) then
         posn%theta = -999.0 
      else 
          posn%theta = bilin(grd,thetao(:,:,floor(this%k)), this%ine, this%jne, this%xi, this%yj) 
      endif
!      write(stderrunit,'(a,i3,a,i4,3f12.4)') &
!                     'particles, theta: pe=(',mpp_pe(),') k, xi, xj, theta', &
!                     floor(this%k), this%lon, this%lat, posn%theta  
      posn%depth=this%depth
      posn%year=parts%current_year
      posn%day=parts%current_yearday
      posn%id=this%id
      posn%particle_num=this%drifter_num
      if (.not. parts%save_short_traj) then !Not totally sure that this is correct
        posn%uvel=this%uvel
        posn%vvel=this%vvel
!        posn%axn=this%axn
!        posn%ayn=this%ayn
!        posn%bxn=this%bxn
!        posn%byn=this%byn
        posn%uvel_old=this%uvel_old
        posn%vvel_old=this%vvel_old
        posn%lon_old=this%lon_old
        posn%lat_old=this%lat_old
      !  !posn%halo_part=this%halo_part
      endif

      call push_posn(this%trajectory, posn)

      this=>this%next
    enddo
  enddo ; enddo

end subroutine record_posn

! ##############################################################################
!> Add trajectory values as a new record in a trajectory
subroutine push_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory !< Trajectory list
type(xyt) :: posn_vals !< Values to add
! Local variables
type(xyt), pointer :: new_posn

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>trajectory
  trajectory=>new_posn

end subroutine push_posn

subroutine append_posn(trajectory, posn_vals)
! This routine appends a new position leaf to the end of the given trajectory
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn,next,last

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>null()
  if(.NOT. associated(trajectory)) then
     trajectory=>new_posn
  else
     ! Find end of the trajectory and point it to the  new leaf
     next=>trajectory
     do while (associated(next))
        last=>next
        next=>next%next
     enddo
     last%next=>new_posn
  endif
end subroutine append_posn

! ##############################################################################
!> Disconnect a trajectory from a part and add it to a list of trajectory segments
subroutine move_trajectory(parts, part)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
type(particle), pointer :: part !< part containing trajectory
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals

  if (parts%ignore_traj) return

  ! If the trajectory is empty, ignore it
  if (.not.associated(part%trajectory)) return

  ! Find end of part trajectory and point it to start of existing trajectories
  next=>part%trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>parts%trajectories

  parts%trajectories=>part%trajectory
  part%trajectory=>null()

end subroutine move_trajectory
! ##############################################################################
!> Scan all parts in a list and disconnect trajectories and more to the list of trajectory segments
!! \todo The argument delete_parts should be removed.
subroutine move_all_trajectories(parts, delete_parts)
! Arguments
type(particles),    pointer    :: parts !< Container for all types and memory
logical, optional, intent(in) :: delete_parts !< If true, delete parts after disconnecting its trajectory
! Local variables
type(particle), pointer :: this, next
logical :: delete_parts_after_moving_traj
integer :: grdi, grdj

  if (parts%ignore_traj) return

  delete_parts_after_moving_traj = .false.
  if (present(delete_parts)) delete_parts_after_moving_traj = delete_parts
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      next=>this%next
      call move_trajectory(parts, this)
   !  if (delete_parts_after_moving_traj) call destroy_particle(this)
      this=>next
    enddo
  enddo ; enddo

end subroutine move_all_trajectories

! ##############################################################################

logical function find_cell_by_search(grd, x, y, i, j)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
integer, intent(inout) :: i, j
! Local variables
integer :: is,ie,js,je,di,dj,io,jo,icnt
real :: d0,d1,d2,d3,d4,d5,d6,d7,d8,dmin
logical :: explain=.false.
real :: Lx

911 continue

  Lx=grd%Lx
  find_cell_by_search=.false.
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec

  ! Start at nearest corner
  d1=dcost(x,y,grd%lonc(is+1,js+1),grd%latc(is+1,js+1),Lx)
  d2=dcost(x,y,grd%lonc(ie-1,js+1),grd%latc(ie-1,js+1),Lx)
  d3=dcost(x,y,grd%lonc(ie-1,je-1),grd%latc(ie-1,je-1),Lx)
  d4=dcost(x,y,grd%lonc(is+1,je-1),grd%latc(is+1,je-1),Lx)
  dmin=min(d1,d2,d3,d4)
  if (d1==dmin) then; i=is+1; j=js+1
  elseif (d2==dmin) then; i=ie-1; j=js+1
  elseif (d3==dmin) then; i=ie-1; j=je-1
  elseif (d4==dmin) then; i=is+1; j=je-1
  else
    call error_mesg('particles, find_cell_by_search:', 'This should never EVER happen! (1)', FATAL)
  endif

  if (explain) then
    write(0,'(i3,a,2i4,f9.4)') mpp_pe(),'Initial corner i-is,j-js,cost=',i-is,j-js,dmin
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d4,d3
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d1,d2
  endif

  if (is_point_in_cell(grd, x, y, i, j)) then
    find_cell_by_search=.true.
    return
  endif

  do icnt=1, 1*(ie-is+je-js)
    io=i; jo=j

    d0=dcost(x,y,grd%lonc(io,jo),grd%latc(io,jo),Lx)
    d1=dcost(x,y,grd%lonc(io,jo+1),grd%latc(io,jo+1),Lx)
    d2=dcost(x,y,grd%lonc(io-1,jo+1),grd%latc(io-1,jo+1),Lx)
    d3=dcost(x,y,grd%lonc(io-1,jo),grd%latc(io-1,jo),Lx)
    d4=dcost(x,y,grd%lonc(io-1,jo-1),grd%latc(io-1,jo-1),Lx)
    d5=dcost(x,y,grd%lonc(io,jo-1),grd%latc(io,jo-1),Lx)
    d6=dcost(x,y,grd%lonc(io+1,jo-1),grd%latc(io+1,jo-1),Lx)
    d7=dcost(x,y,grd%lonc(io+1,jo),grd%latc(io+1,jo),Lx)
    d8=dcost(x,y,grd%lonc(io+1,jo+1),grd%latc(io+1,jo+1),Lx)

  ! dmin=min(d0,d1,d3,d5,d7)
    dmin=min(d0,d1,d2,d3,d4,d5,d6,d7,d8)
    if (d0==dmin) then; di=0; dj=0
    elseif (d2==dmin) then; di=-1; dj=1
    elseif (d4==dmin) then; di=-1; dj=-1
    elseif (d6==dmin) then; di=1; dj=-1
    elseif (d8==dmin) then; di=1; dj=1
    elseif (d1==dmin) then; di=0; dj=1
    elseif (d3==dmin) then; di=-1; dj=0
    elseif (d5==dmin) then; di=0; dj=-1
    elseif (d7==dmin) then; di=1; dj=0
    else
      call error_mesg('particles, find_cell_by_search:', 'This should never EVER happen!', FATAL)
    endif

    i=min(ie, max(is, io+di))
    j=min(je, max(js, jo+dj))

    if (explain) then
      write(0,'(i3,a,2i4,f9.5,a,2i4,a,2i4)') mpp_pe(),'Current position i,j,cost=',i,j,dmin,' di,dj=',di,dj,' old io,jo=',io,jo
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    endif

    if (is_point_in_cell(grd, x, y, i, j)) then
      find_cell_by_search=.true.
      return
    endif

    if ((i==io.and.j==jo) &
        .and. .not.find_better_min(grd, x, y, 3, i, j) &
       ) then
      ! Stagnated
      find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 1, i, j)
      if (.not. find_cell_by_search) find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 3, i, j)
      i=min(ie, max(is, i))
      j=min(je, max(js, j))
      if (is_point_in_cell(grd, x, y, i, j)) then
        find_cell_by_search=.true.
      else
  !     find_cell_by_search=find_cell(grd, x, y, io, jo)
  !     if (find_cell_by_search) then
  !       if (explain) then
  !         call print_fld(grd, grd%lat, 'Lat')
  !         call print_fld(grd, grd%lon, 'Lat')
  !         do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
  !           grd%tmp(i,j) = dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
  !         enddo; enddo
  !         call print_fld(grd, grd%tmp, 'Cost')
  !         stop 'Avoid recursing'
  !       endif
  !       write(0,'(i3,a,2i5,a,2i3,a,2f8.3)') mpp_pe(),'particles, find_cell_by_search: false negative io,jo=',io,jo,' di,dj=',di,dj,' targ=',x,y
  !       explain=.true.; goto 911
  !     endif
      endif
      return
    endif

  enddo

  find_cell_by_search=find_cell(grd, x, y, i, j)
  if (find_cell_by_search) then
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    write(0,'(i3,a,2f9.5)') mpp_pe(),'x,y ',x,y
    write(0,'(i3,a,4i5)') mpp_pe(),'io,jo ',io,jo,di,dj
    write(0,'(i3,a,2i5,a,2i3)') mpp_pe(),'particles, find_cell_by_search: false negative 2 i,j=',i-is,j-js,' di,dj=',di,dj
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'particles, find_cell_by_search: false negative 2 io,jo=',io,jo
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'particles, find_cell_by_search: false negative 2 i,j=',i,j,' targ=',x,y
    return
  endif
  find_cell_by_search=.false.

  contains

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  real function dcost(x1, y1, x2, y2,Lx)
  ! Arguments
  real, intent(in) :: x1 !< x-coordinate or point 1
  real, intent(in) :: x2 !< x-coordinate or point 2
  real, intent(in) :: y1 !< y-coordinate or point 1
  real, intent(in) :: y2 !< y-coordinate or point 2
  real, intent(in) :: Lx !< Periodicity of x-coordinate
  ! Local variables
  real :: x1m
  x1m=apply_modulo_around_point(x1,x2,Lx)
  !x1m=modulo(x1-(x2-180.),360.)+(x2-180.)
  ! dcost=(x2-x1)**2+(y2-y1)**2
    dcost=(x2-x1m)**2+(y2-y1)**2
  end function dcost

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  logical function find_better_min(grd, x, y, w, oi, oj)
  ! Arguments
  type(particles_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye
  real :: dmin, dcst
  real :: Lx

  Lx=grd%Lx

  xs=max(grd%isc, oi-w)
  xe=min(grd%iec, oi+w)
  ys=max(grd%jsc, oj-w)
  ye=min(grd%jec, oj+w)

  find_better_min=.false.
  dmin=dcost(x,y,grd%lonc(oi,oj),grd%latc(oi,oj),Lx)
  do j=ys,ye; do i=xs,xe
      dcst=dcost(x,y,grd%lonc(i,j),grd%latc(i,j),Lx)
      if (dcst<dmin) then
        find_better_min=.true.
        dmin=dcst
        oi=i;oj=j
      endif
  enddo; enddo

  end function find_better_min

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  logical function find_cell_loc(grd, x, y, is, ie, js, je, w, oi, oj)
  ! Arguments
  type(particles_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: is, ie, js, je, w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye

    xs=max(is, oi-w)
    xe=min(ie, oi+w)
    ys=max(js, oj-w)
    ye=min(je, oj+w)

    find_cell_loc=.false.
    do j=ys,ye; do i=xs,xe
        if (is_point_in_cell(grd, x, y, i, j)) then
          oi=i; oj=j; find_cell_loc=.true.
          return
        endif
    enddo; enddo

  end function find_cell_loc

end function find_cell_by_search

! ##############################################################################

logical function find_cell(grd, x, y, oi, oj)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell=.false.; oi=-999; oj=-999

  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
        return
      endif
  enddo; enddo

end function find_cell

! ##############################################################################

logical function find_cell_wide(grd, x, y, oi, oj)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell_wide=.false.; oi=-999; oj=-999

  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

! ##############################################################################

!> Returns the i,j of cell containing an particle with the given identifier
subroutine find_individual_particle(parts, id, ine, jne, part_found, search_data_domain)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
integer(kind=8), intent(in) :: id !< part identifier
integer, intent(out) :: ine !< i-index of cell containing part
integer, intent(out) :: jne !< j-index of cell containing part
logical, intent(in) :: search_data_domain !< If true, search halos too
real, intent(out) :: part_found !< Returns 1.0 if part is found, 0. otherwise
! Local variables
type(particle), pointer :: this
type(particles_gridded), pointer :: grd
integer :: grdi, grdj
integer :: ilim1, ilim2, jlim1, jlim2

part_found=0.0
ine=999
jne=999
  ! For convenience
    grd=>parts%grd

    if (search_data_domain) then
        ilim1 = grd%isd  ; ilim2=grd%ied  ; jlim1 = grd%jsd  ; jlim2=grd%jed
    else
        ilim1 = grd%isc  ; ilim2=grd%iec  ; jlim1 = grd%jsc  ; jlim2=grd%jec
    endif


    do grdj = jlim1, jlim2 ; do grdi = ilim1, ilim2
    !do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    !do grdj = parts%grd%jsd,parts%grd%jed ; do grdi = parts%grd%isd,parts%grd%ied
      this=>parts%list(grdi,grdj)%first
      do while (associated(this))
        if (id .eq. this%id) then
          ine=this%ine
          jne=this%jne
          part_found=1.0
          !print *, 'found this one'
          return
        endif
        this=>this%next
      enddo
    enddo ; enddo
end subroutine  find_individual_particle

! ##############################################################################

logical function is_point_in_cell(grd, x, y, i, j, explain)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
logical, intent(in), optional :: explain
! Local variables
real :: xlo, xhi, ylo, yhi
integer :: stderrunit
real :: Lx
real :: tol

  ! Get the stderr unit number
  stderrunit=stderr()
  Lx=grd%Lx

  ! Safety check index bounds
  if (i-1.lt.grd%isd.or.i.gt.grd%ied.or.j-1.lt.grd%jsd.or.j.gt.grd%jed) then
    write(stderrunit,'(a,i3,(a,3i4))') &
                     'particles, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('particles, is_point_in_cell', 'test is off the PE!', FATAL)
  endif

  is_point_in_cell=.false.

  ! Test crude bounds
  xlo=min( apply_modulo_around_point(grd%lon(i-1,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i-1,j  )   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j  )   ,x,  Lx) )
  xhi=max( apply_modulo_around_point(grd%lon(i-1,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i-1,j  )   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j  )   ,x,  Lx) )

  ! The modolo function inside sum_sign_dot_prod leads to a roundoff.
  !Adding adding a tolorance to the crude bounds avoids excluding the cell which
  !would be correct after roundoff. This is a bit of a hack.
  tol=0.1
  if (x.lt.(xlo-tol) .or. x.gt.(xhi+tol)) return

  ylo=min( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  yhi=max( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  if (y.lt.ylo .or. y.gt.yhi) return

  if ((grd%lat(i,j).gt.89.999).and. (grd%grid_is_latlon))   then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, Lx,explain=explain)
  elseif ((grd%lat(i-1,j).gt.89.999) .and. (grd%grid_is_latlon))  then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i  ,j  ),grd%lat(i-1,j  ), &
                                        grd%lon(i-1,j-1),grd%lat(i-1,j  ), &
                                        x, y,Lx, explain=explain)
  elseif ((grd%lat(i-1,j-1).gt.89.999) .and. (grd%grid_is_latlon))  then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j  ),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y,Lx, explain=explain)
  elseif ((grd%lat(i,j-1).gt.89.999) .and. (grd%grid_is_latlon)) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i-1,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, Lx,explain=explain)
  else
  is_point_in_cell=sum_sign_dot_prod4(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y,Lx, explain=explain)
  endif

end function is_point_in_cell
! ##############################################################################

logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y, Lx, explain)

  ! Get the stderr unit numberrguments
real, intent(in) :: x0 !< Longitude of first corner
real, intent(in) :: y0 !< Latitude of first corner
real, intent(in) :: x1 !< Longitude of second corner
real, intent(in) :: y1 !< Latitude of second corner
real, intent(in) :: x2 !< Longitude of third corner
real, intent(in) :: y2 !< Latitude of third corner
real, intent(in) :: x3 !< Longitude of fourth corner
real, intent(in) :: y3 !< Latitude of fourth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in), optional :: explain !< If true, print debugging
! Local variables
real :: p0,p1,p2,p3,xx
real :: l0,l1,l2,l3
real :: xx0,xx1,xx2,xx3
integer :: stderrunit
  stderrunit=stderr()

  sum_sign_dot_prod4=.false.
  xx= apply_modulo_around_point(x,x0,Lx)
  xx0= apply_modulo_around_point(x0,x0,Lx)
  xx1= apply_modulo_around_point(x1,x0,Lx)
  xx2= apply_modulo_around_point(x2,x0,Lx)
  xx3= apply_modulo_around_point(x3,x0,Lx)

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y0-y3)-(y-y3)*(xx0-xx3)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    sum_sign_dot_prod4=.true.
  endif


  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: x=',mpp_pe(),':', &
                           x0,x1,x2,x3, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: y=',mpp_pe(),':', &
                           y0,y1,y2,y3, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: l=',mpp_pe(),':', &
                           l0,l1,l2,l3
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: p=',mpp_pe(),':', &
                           p0,p1,p2,p3, abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) )
   endif
  endif

end function sum_sign_dot_prod4

! ##############################################################################

logical function sum_sign_dot_prod5(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, Lx, explain)
! Arguments
real, intent(in) :: x0 !< Longitude of first corner
real, intent(in) :: y0 !< Latitude of first corner
real, intent(in) :: x1 !< Longitude of second corner
real, intent(in) :: y1 !< Latitude of second corner
real, intent(in) :: x2 !< Longitude of third corner
real, intent(in) :: y2 !< Latitude of third corner
real, intent(in) :: x3 !< Longitude of fourth corner
real, intent(in) :: y3 !< Latitude of fourth corner
real, intent(in) :: x4 !< Longitude of fifth corner
real, intent(in) :: y4 !< Latitude of fifth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in), optional :: explain !< If true, print debugging
! Local variables
real :: p0,p1,p2,p3,p4,xx
real :: l0,l1,l2,l3,l4
real :: xx0,xx1,xx2,xx3,xx4
integer :: stderrunit
  stderrunit=stderr()

  sum_sign_dot_prod5=.false.
  xx= apply_modulo_around_point(x,x0,Lx)
  xx0= apply_modulo_around_point(x0,x0,Lx)
  xx1= apply_modulo_around_point(x1,x0,Lx)
  xx2= apply_modulo_around_point(x2,x0,Lx)
  xx3= apply_modulo_around_point(x3,x0,Lx)
  xx4= apply_modulo_around_point(x4,x0,Lx)

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y4-y3)-(y-y3)*(xx4-xx3)
  l4=(xx-xx4)*(y0-y4)-(y-y4)*(xx0-xx4)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.
  p4=sign(1., l4); if (l4.eq.0.) p4=0.

  if ( ((abs(p0)+abs(p2))+(abs(p1)+abs(p3)))+abs(p4) - abs(((p0+p2)+(p1+p3))+p4) .lt. 0.5 ) then
    sum_sign_dot_prod5=.true.
  endif

  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: x=',mpp_pe(),':', &
                           x0,x1,x2,x3,x4, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3,xx4, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: y=',mpp_pe(),':', &
                           y0,y1,y2,y3,y4, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod5: l=',mpp_pe(),':', &
                           l0,l1,l2,l3,l4
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: p=',mpp_pe(),':', &
                           p0,p1,p2,p3,p4
   endif
  endif

end function sum_sign_dot_prod5

! ##############################################################################

logical function pos_within_cell(grd, x, y, i, j, xi, yj, explain)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
real, intent(out) :: xi, yj
logical, intent(in), optional :: explain
! Local variables
real :: x1,y1,x2,y2,x3,y3,x4,y4,xx,yy,fac
integer :: stderrunit
real :: Lx, dx,dy
real :: Delta_x
logical :: is_point_in_cell_using_xi_yj

  ! Get the stderr unit number
  stderrunit=stderr()
  Lx=grd%Lx

  pos_within_cell=.false.; xi=-999.; yj=-999.
  if (i-1<grd%isd) return
  if (j-1<grd%jsd) return
  if (i>grd%ied) return
  if (j>grd%jed) return

  x1=grd%lon(i-1,j-1)
  y1=grd%lat(i-1,j-1)
  x2=grd%lon(i  ,j-1)
  y2=grd%lat(i  ,j-1)
  x3=grd%lon(i  ,j  )
  y3=grd%lat(i  ,j  )
  x4=grd%lon(i-1,j  )
  y4=grd%lat(i-1,j  )

  if (present(explain)) then
    if(explain) then
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',x
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',y
    endif
  endif

  if ((.not. grd%grid_is_latlon) .and. (grd%grid_is_regular))  then
    dx=abs((grd%lon(i  ,j  )-grd%lon(i-1  ,j  )))
    dy=abs((grd%lat(i  ,j  )-grd%lat(i  ,j-1  )))
    x1=grd%lon(i  ,j  )-(dx/2)
    y1=grd%lat(i  ,j  )-(dy/2)

    Delta_x= apply_modulo_around_point(x,x1,Lx)-x1

    xi=((Delta_x)/dx)+0.5
    !xi=((x-x1)/dx)+0.5
    yj=((y-y1)/dy)+0.5

  elseif ((max(y1,y2,y3,y4)<89.999) .or.(.not. grd%grid_is_latlon)) then
      ! This returns non-dimensional position xi,yj for quad cells (not at a
      ! pole)
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, Lx,explain=explain)
  else
  if (debug) write(stderrunit,*) 'particles, pos_within_cell: working in tangential plane!'
    xx=(90.-y)*cos(x*pi_180)
    yy=(90.-y)*sin(x*pi_180)
    x1=(90.-y1)*cos(grd%lon(i-1,j-1)*pi_180)
    y1=(90.-y1)*sin(grd%lon(i-1,j-1)*pi_180)
    x2=(90.-y2)*cos(grd%lon(i  ,j-1)*pi_180)
    y2=(90.-y2)*sin(grd%lon(i  ,j-1)*pi_180)
    x3=(90.-y3)*cos(grd%lon(i  ,j  )*pi_180)
    y3=(90.-y3)*sin(grd%lon(i  ,j  )*pi_180)
    x4=(90.-y4)*cos(grd%lon(i-1,j  )*pi_180)
    y4=(90.-y4)*sin(grd%lon(i-1,j  )*pi_180)
    if (present(explain)) then
      if(explain) then
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',xx
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',yy
      endif
    endif
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, xx, yy, xi, yj, Lx, explain=explain)
    if (is_point_in_cell(grd, x, y, i, j)) then
      if (abs(xi-0.5)>0.5.or.abs(yj-0.5)>0.5) then
        ! Scale internal coordinates to be consistent with is_point_in_cell()
        ! Note: this is intended to fix the inconsistency between the tangent plane
        ! and lat-lon calculations
        fac=2.*max( abs(xi-0.5), abs(yj-0.5) ); fac=max(1., fac)
        xi=0.5+(xi-0.5)/fac
        yj=0.5+(yj-0.5)/fac
        if (debug) call error_mesg('particles, pos_within_cell', 'in cell so scaling internal coordinates!', WARNING)
      endif
    else
      if (abs(xi-0.5)<=0.5.and.abs(yj-0.5)<=0.5) then
        if (debug) call error_mesg('particles, pos_within_cell', 'out of cell but coordinates <=0.5!', WARNING)
      endif
    endif
  endif

  if (present(explain)) then
        if(explain) write(stderrunit,'(a,2f12.6)') 'pos_within_cell: xi,yj=',xi,yj
  endif

 !if (.not. is_point_in_cell(grd, x, y, i, j) ) then
 !   write(stderrunit,'(a,i3,a,8f8.2,a)') 'particles, pos_within_cell: (',mpp_pe(),') ', &
 !                   x1, y1, x2, y2, x3, y3, x4, y4, ' NOT IN CELL!'
 !endif

  if (xi.ge.0. .and. xi.le.1. .and. yj.ge.0. .and. yj.le.1.) then
    pos_within_cell=is_point_in_cell(grd, x, y, i, j, explain=explain)
    if (.not. pos_within_cell .and. verbose) then
      if (debug) call error_mesg('particles, pos_within_cell', 'pos_within_cell is in cell BUT is_point_in_cell disagrees!', WARNING)
    endif
   !pos_within_cell=.true. ! commenting this out makes pos_within_cell agree with is_point_in_cell
  endif

  contains

  subroutine calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, Lx, explain)
  ! Arguments
  real, intent(in) :: x1 !< Longitude of first corner
  real, intent(in) :: y1 !< Latitude of first corner
  real, intent(in) :: x2 !< Longitude of second corner
  real, intent(in) :: y2 !< Latitude of second corner
  real, intent(in) :: x3 !< Longitude of third corner
  real, intent(in) :: y3 !< Latitude of third corner
  real, intent(in) :: x4 !< Longitude of fourth corner
  real, intent(in) :: y4 !< Latitude of fourth corner
  real, intent(in) :: x !< Longitude of point
  real, intent(in) :: y !< Latitude of point
  real, intent(out) :: xi !< Non-dimensional x-position within cell
  real, intent(out) :: yj !< Non-dimensional y-position within cell
  real, intent(in) :: Lx !< Length of domain in zonal direction
  logical, intent(in), optional :: explain !< If true, print debugging
  ! Local variables
  real :: alpha, beta, gamma, delta, epsilon, kappa, a, b, c, d, dx, dy, yy1, yy2
  logical :: expl=.false.
  integer :: stderrunit
  ! Get the stderr unit number
  stderrunit=stderr()

  expl=.false.
  if (present(explain)) then
     if(explain) expl=.true.
  endif
  alpha=x2-x1
  delta=y2-y1
  beta=x4-x1
  epsilon=y4-y1
  gamma=(x3-x1)-(alpha+beta)
  kappa=(y3-y1)-(delta+epsilon)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs alpha,beta,gamma',alpha,beta,gamma,delta,epsilon,kappa
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs delta,epsilon,kappa',alpha,beta,gamma,delta,epsilon,kappa

  a=(kappa*beta-gamma*epsilon)
  dx= apply_modulo_around_point(x,x1,Lx)-x1

  dy=y-y1
  b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
  c=(alpha*dy-delta*dx)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs dx,dy=',dx,dy
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs A,B,C=',a,b,c
  if (abs(a)>1.e-12) then
    d=0.25*(b**2)-a*c
    if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs D=',d
    if (d.ge.0.) then
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for b/2a, sqrt(d) = ',-0.5*b/a,sqrt(d)/a
      yy1=-(0.5*b+sqrt(d))/a
      yy2=-(0.5*b-sqrt(d))/a
      if (abs(yy1-0.5).lt.abs(yy2-0.5)) then; yj=yy1; else; yj=yy2; endif
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for y = ',yy1,yy2,yj
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3)') 'calc_xiyj: b<0 in quadratic root solver!!!!',mpp_pe()
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,d,dx,dy
      call error_mesg('particles, calc_xiyj', 'We have complex roots. The grid must be very distorted!', FATAL)
    endif
  else
    if (b.ne.0.) then
      yj=-c/b
    else
      yj=0.
    endif
  endif

  a=(alpha+gamma*yj)
  b=(delta+kappa*yj)
  if (a.ne.0.) then
    xi=(dx-beta*yj)/a
  elseif (b.ne.0.) then
    xi=(dy-epsilon*yj)/b
  else
    c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
    if (c.ne.0.) then
      xi=(epsilon*dx-beta*dy)/c
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3,1p2e12.4)') 'calc_xiyj: coeffs a,b',mpp_pe(),a,b
      call error_mesg('particles, calc_xiyj', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
    endif
  endif
  if (expl) write(stderrunit,'(a,2e12.4)') 'calc_xiyj: xi,yj=',xi,yj

  end subroutine calc_xiyj

end function pos_within_cell


! ##############################################################################

real function find_layer1D(grd, depth,hdepth)
!Arguments
type(particles_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: depth
real, dimension(:),intent(in) :: hdepth
!real, intent(out) :: k

!Local                                                                          
integer :: klev
real :: rdepth
integer :: stderrunit


  ! Get the stderr unit number                                                                                
  stderrunit=stderr()




if (depth.lt.hdepth(1)) then
   find_layer1D=depth/hdepth(1)
   return
endif



do klev=2,grd%ke
   if (depth.lt.hdepth(klev)) then
      find_layer1D=klev - 1 + (depth-hdepth(klev-1))/(hdepth(klev)-hdepth(klev-1))
      return
   endif
enddo


write(stderrunit,'(a)')"particles: depth specified is deeper than deepest level"
write(stderrunit,*) 'depths', hdepth(1), hdepth(2)

end function find_layer1D

! ##############################################################################

!>Finds what layer the particle is in
subroutine find_layer(grd, depth,hdepth,k,ine,jne,xi,yj,k_space)
!Arguments
type(particles_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: depth
real, dimension(:,:,:),intent(in) :: hdepth
real, intent(out) :: k
integer, intent(in) :: ine
integer, intent(in) :: jne
real, intent(in) :: xi
real, intent(in) :: yj
logical, intent(inout) :: k_space


!Local
integer :: klev
real :: rdepth
integer :: stderrunit
real, dimension(grd%ke) :: hdepth1D
real, dimension(grd%ke) :: hdepth1Dcum                                                

  ! Get the stderr unit number
  stderrunit=stderr()

if (k_space)then
   return
endif


do klev=1,grd%ke
!    hdepth1D(klev) = bilin(grd, hdepth(:,:,klev), ine, jne, xi, yj)
     hdepth1D(klev) = hdepth(ine,jne,klev)
    if (klev.eq.1)then
        hdepth1Dcum(klev) =hdepth1D(klev)
    else
        hdepth1Dcum(klev)=hdepth1Dcum(klev-1)+hdepth1D(klev)
    endif 
enddo 

   k = find_layer1D(grd, depth,hdepth1Dcum)

k_space=.true.


end subroutine find_layer

! ##############################################################################
real function find_depth1D(grd,k,h)
!Arguments
type(particles_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: k
real, dimension(:),intent(in) :: h
!real, intent(out) :: depth

!Local
integer :: klev
real :: rdepth
integer :: kint
integer :: stderrunit


  ! Get the stderr unit number

  stderrunit=stderr()
  
kint=floor(k)

!if (kint.eq.0)then
!   depth= h(1)*k
!else
!   depth = h(kint)+(h(kint+1)-h(kint))*(k-kint)
!endif
find_depth1D=0
do klev=0,kint
   if (klev.eq.kint)then
     find_depth1D = find_depth1D+h(klev+1)*(k-kint)
   else
     find_depth1D=find_depth1D+h(klev+1)
   endif
enddo

end function find_depth1D

! ############################################################################## 
!>Finds what depth the particle is at  
subroutine find_depth(grd,k,h,depth,ine,jne,xi,yj,k_space)
!Arguments
type(particles_gridded), pointer :: grd !< Container for gridded fields 
real, intent(in) :: k
real, dimension(:,:,:),intent(in) :: h
real, intent(out) :: depth
integer, intent(in) :: ine
integer, intent(in) :: jne
real, intent(in) :: xi
real, intent(in) :: yj
logical, intent(inout) :: k_space

!Local
real, dimension(grd%ke) :: hdepth1D
integer :: klev

if (.not.k_space)then
   return
endif


do klev=1,grd%ke
!    hdepth1D(klev) = bilin(grd, h(:,:,klev), ine, jne, xi, yj)
    hdepth1D(klev) = h(ine,jne,klev)
enddo

k_space=.false.

depth = find_depth1D(grd,k,hdepth1D)

end subroutine find_depth

! ##############################################################################


!> Checks that a part's position metrics are consistent
subroutine check_position(grd, part, label, il, jl)
! Arguments
type(particles_gridded), pointer :: grd !< Container for gridded fields
type(particle), pointer :: part !< part to check
character(len=*) :: label !< Label to add to messages
integer, optional, intent(in) :: il !< i-index of cell part should be in
integer, optional, intent(in) :: jl !< j-index of cell part should be in
! Local variables
real :: xi, yj
logical :: lret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  lret=pos_within_cell(grd, part%lon, part%lat, part%ine, part%jne, xi, yj)
  if (xi.ne.part%xi.or.yj.ne.part%yj) then
    write(stderrunit,'("particles: check_position (",i4,") b%x,x,-=",3(es12.4,x),a)') mpp_pe(),part%xi,xi,part%xi-xi,label
    write(stderrunit,'("particles: check_position (",i4,") b%y,y,-=",3(es12.4,x),a)') mpp_pe(),part%yj,yj,part%yj-yj,label
    call print_part(stderrunit, part, 'check_position', il, jl)
    call error_mesg('particles, check_position, '//trim(label),'part has inconsistent xi,yj!',FATAL)
  endif
  !if (grd%msk(part%ine, part%jne)==0.) then
  !  call print_part(stderrunit, part, 'check_position, '//trim(label), il, jl)
  !  call error_mesg('particles, check_position, '//trim(label),'part is in a land cell!',FATAL)
  !endif

end subroutine check_position




subroutine checksum_gridded(grd, label)
! Arguments
type(particles_gridded) :: grd
character(len=*) :: label
! Local variables

  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'particles: checksumming gridded data @ ',trim(label)

  ! static
!  call grd_chksum2(grd, grd%lon, 'lon')
!  call grd_chksum2(grd, grd%lat, 'lat')
!  call grd_chksum2(grd, grd%lonc, 'lonc')
!  call grd_chksum2(grd, grd%latc, 'latc')
!  call grd_chksum2(grd, grd%dx, 'dx')
!  call grd_chksum2(grd, grd%dy, 'dy')
!  call grd_chksum2(grd, grd%msk, 'msk')
!  call grd_chksum2(grd, grd%cos, 'cos')
!  call grd_chksum2(grd, grd%sin, 'sin')
!  call grd_chksum2(grd, grd%ocean_depth, 'depth')

end subroutine checksum_gridded

! ##############################################################################

subroutine grd_chksum3(grd, fld, txt)
! Arguments
type(particles_gridded) :: grd
real, dimension(:,:,:), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, k, halo, icount, io, jo
real :: mean, rms, SD, minv, maxv
real, dimension(lbound(fld,1):ubound(fld,1), lbound(fld,2):ubound(fld,2), lbound(fld,3):ubound(fld,3)) :: tmp

  halo=grd%halo
  mean=0.
  rms=0.
  sd=0.
  icount=0
  i=lbound(fld,1)+halo
  j=lbound(fld,2)+halo
  k=lbound(fld,3)
  minv=fld(i,j,k)
  maxv=fld(i,j,k)
  tmp(:,:,:)=0.
  io=grd%isd-lbound(fld,1)
  jo=grd%jsd-lbound(fld,2)
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        icount=icount+1
        mean=mean+fld(i,j,k)
        rms=rms+fld(i,j,k)**2
        minv=min(minv,fld(i,j,k))
        maxv=max(maxv,fld(i,j,k))
        tmp(i,j,k)=fld(i,j,k)*float(i+io+2*(j+jo)+3*(k-1))
      enddo
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        sd=sd+(fld(i,j,k)-mean)**2
      enddo
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  j=mpp_chksum( tmp(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("particles, grd_chksum3: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  j=mpp_chksum( tmp(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("particles, grd_chksum3* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#endif

end subroutine grd_chksum3

! ##############################################################################

subroutine grd_chksum2(grd, fld, txt)
! Arguments
type(particles_gridded) :: grd
real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, icount
real :: mean, rms, SD, minv, maxv

  grd%tmp(:,:)=0.

  mean=0.
  rms=0.
  sd=0.
  icount=0
  minv=fld(grd%isc,grd%jsc)
  maxv=fld(grd%isc,grd%jsc)
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      icount=icount+1
      mean=mean+fld(i,j)
      rms=rms+fld(i,j)**2
      minv=min(minv,fld(i,j))
      maxv=max(maxv,fld(i,j))
      grd%tmp(i,j)=fld(i,j)*float(i+2*j)
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      sd=sd+(fld(i,j)-mean)**2
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(grd%isc:grd%iec,grd%jsc:grd%jec) )
  j=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("particles, grd_chksum2: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i8)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(grd%isd:grd%ied,grd%jsd:grd%jed) )
  j=mpp_chksum( grd%tmp(grd%isd:grd%ied,grd%jsd:grd%jed) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("particles, grd_chksum2* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#endif

end subroutine grd_chksum2

! ##############################################################################
!> Calculates checksums for all parts
subroutine parts_chksum(parts, txt, ignore_halo_violation)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
character(len=*), intent(in) :: txt !< Label to use in messages
logical, optional :: ignore_halo_violation
! Local variables
integer :: i, nparts, ichk1, ichk2, ichk3, ichk4, ichk5, ipart
real, allocatable :: fld(:,:), fld2(:,:)
integer, allocatable :: icnt(:,:)
type(particle), pointer :: this
type(particles_gridded), pointer :: grd
logical :: check_halo
integer :: grdi, grdj

! For convenience
  grd=>parts%grd

  nparts=count_parts(parts)
  call mpp_max(nparts)
  nparts = max(nparts, 1)
  allocate( fld( nparts, 2 ) ) !Changed from 11 to 19 by Alon
  allocate( fld2( nparts, 2 ) ) !Changed from 11 to 19 by Alon
  allocate( icnt( grd%isd:grd%ied, grd%jsd:grd%jed ) )
  fld(:,:)=0.
  fld2(:,:)=0.
  icnt(:,:)=0
  grd%tmp(:,:)=0.

  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    this=>parts%list(grdi,grdj)%first
    i=0; ichk5=0
    do while(associated(this))
      i=i+1
      ipart=part_chksum(this)
      fld(i,1) = this%lon
      fld(i,2) = this%lat
!      fld(i,3) = this%uvel
!      fld(i,4) = this%vvel
!      fld(i,5) = this%uvel_old !added by Alon
!      fld(i,6) = this%vvel_old !added by Alon
!      fld(i,7) = this%lon_old !added by Alon
!      fld(i,8) = this%lat_old !added by Alon
!      fld(i,17) = time_hash(this) !Changed from 9 to 17 by Alon
!      fld(i,18) = pos_hash(this) !Changed from 10 to 18 by Alon
!      fld(i,19) = float(ipart) !Changed from 11 to 19 by Alon
      icnt(this%ine,this%jne)=icnt(this%ine,this%jne)+1
      fld2(i,:) = fld(i,:)*float( icnt(this%ine,this%jne) ) !*float( i )
      grd%tmp(this%ine,this%jne)=grd%tmp(this%ine,this%jne)+time_hash(this)*pos_hash(this)!+log(this%mass) !LUYU: we assume zero mass for now.
      ichk5=ichk5+ipart
      this=>this%next
    enddo
  enddo ; enddo

  ichk1=mpp_chksum( fld )
  ichk2=mpp_chksum( fld2 )
  ichk3=mpp_chksum( grd%tmp )
  ichk4=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( ichk5 )
  nparts=count_parts(parts)

  if (nparts.ne.sum(icnt(:,:))) then
    write(*,'("particles, parts_chksum: ",2(a,i8))') &
      '# parts =', nparts, ' sum(icnt) =',sum(icnt(:,:))
    call error_mesg('particles, parts_chksum:', 'mismatch in part count!', FATAL)
  endif

  check_halo=.true.
  if (present(ignore_halo_violation)) then
    if (ignore_halo_violation) check_halo=.false.
  endif
  if (check_halo.and.nparts.ne.sum(icnt(grd%isc:grd%iec, grd%jsc:grd%jec))) then
    write(*,'("particles, parts_chksum: ",2(a,i8))') &
      '# parts =', nparts, ' sum(icnt(comp_dom)) =',sum(icnt(:,:))
    call error_mesg('particles, parts_chksum:', 'mismatch in part count on computational domain!', FATAL)
  endif

  call mpp_sum(nparts)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("particles, parts_chksum: ",a18,6(x,a,"=",i22))') &
      txt, 'chksum', ichk1, 'chksum2', ichk2, 'chksum3', ichk3, 'chksum4', ichk4, 'chksum5', ichk5, '#', nparts

  grd%tmp(:,:)=real(icnt(:,:))
!  call grd_chksum2(grd,grd%tmp,'# of parts/cell')

  deallocate( fld )
  deallocate( fld2 )
  deallocate( icnt )

  if (debug) call count_out_of_order(parts,txt)

end subroutine parts_chksum

! ##############################################################################
!> Checksum a list of parts
integer function list_chksum(first)
! Arguments
type(particle), pointer :: first !< List of parts
! Local variables
integer :: i
type(particle), pointer :: this

  this=>first
  i=0; list_chksum=0
  do while(associated(this))
    i=i+1
    list_chksum=list_chksum+part_chksum(this)*i
    this=>this%next
  enddo

end function list_chksum
! ##############################################################################

integer function part_chksum(part )
! Arguments
type(particle), pointer :: part
! Local variables
real :: rtmp(14) !Changed from 28 to 34 by Alon
integer :: itmp(19), i8=0, ichk1, ichk2, ichk3 !Changed from 28 to 34 by Alon
integer :: i

  rtmp(:)=0.
  rtmp(1)=part%lon
  rtmp(2)=part%lat
  rtmp(3)=part%uvel
  rtmp(4)=part%vvel
  rtmp(5)=part%start_lon
  rtmp(6)=part%start_lat
  rtmp(7)=part%start_d
  rtmp(8)=part%start_day
  rtmp(9)=part%xi
  rtmp(10)=part%yj
  rtmp(11)=part%uvel_old !Added by Alon
  rtmp(12)=part%vvel_old !Added by Alon
  rtmp(13)=part%lat_old !Added by Alon
  rtmp(14)=part%lon_old !Added by Alon

  itmp(1:14)=transfer(rtmp,i8) !Changed from 28 to 36 by Alon
  itmp(15)=part%start_year !Changed from 29 to 37 by Alon
  itmp(16)=part%ine !Changed from 30 to 38 by Alon
  itmp(17)=part%jne !Changed from 31 to 39 by Alon
  call split_id(part%id, itmp(18), itmp(19))

  ichk1=0; ichk2=0; ichk3=0
  do i=1,19 !Changd from 28 to 37 by Alon
   ichk1=ichk1+itmp(i)
   ichk2=ichk2+itmp(i)*i
   ichk3=ichk3+itmp(i)*i*i
  enddo
  part_chksum=ichk1+ichk2+ichk3

end function part_chksum

! ##############################################################################

real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
! Local variables


    bilin=(fld(i,j  )*xi+fld(i-1,j  )*(1.-xi))*yj &
         +(fld(i,j-1)*xi+fld(i-1,j-1)*(1.-xi))*(1.-yj)

end function bilin

! #############################################################################

real function linlinx(grd,fld,x,y,i,j,xi,yj)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
real, intent(in) :: fld(grd%isd:grd%ied+1,grd%jsd:grd%jed+1), xi, yj
integer, intent(in) :: i, j
real :: linlinxnum, linlinxden
! Local variables
logical :: explain=.false.
integer :: stderrunit
  stderrunit=stderr()

    linlinx = fld(i,j  )*xi + fld(i-1,j)*(1-xi)

end function linlinx

real function find_u(grd,fld,x, y, i, j, xi, yj, dx_dlon, dy_dlat)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
real, intent(in) :: fld(grd%isd:grd%ied+1,grd%jsd:grd%jed+1), xi, yj
!real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
real, intent(in) :: dx_dlon, dy_dlat
real :: Ul, Vl, u, v
real :: linlinxnum, jacden
! Local variables
real :: x0 !< Longitude of first corner
real :: y0 !< Latitude of first corner
real :: x1 !< Longitude of second corner
real :: y1 !< Latitude of second corner
real :: x2 !< Longitude of third corner
real :: y2 !< Latitude of third corner
real :: x3 !< Longitude of fourth corner
real :: y3 !< Latitude of fourth corner
real :: xx, xx0,xx1,xx2,xx3
real :: l0, l1, l2, l3
real :: dxidt, detadt
real :: Lx
integer :: stderrunit
  stderrunit=stderr()

x0 = grd%lon(i-1,j-1)
x1 = grd%lon(i  ,j-1)
x2 = grd%lon(i  ,j  )
x3 = grd%lon(i-1,j  )

y0 = grd%lat(i-1,j-1)
y1 = grd%lat(i  ,j-1)
y2 = grd%lat(i  ,j  )
y3 = grd%lat(i-1,j  )

Lx = grd%Lx

xx= apply_modulo_around_point(x,x0,Lx)
xx0= apply_modulo_around_point(x0,x0,Lx)
xx1= apply_modulo_around_point(x1,x0,Lx)
xx2= apply_modulo_around_point(x2,x0,Lx)
xx3= apply_modulo_around_point(x3,x0,Lx)

Ul = linlinx(grd,fld,x, y, i, j, xi, yj)
Vl = linliny(grd,fld,x, y, i, j, xi, yj)
u = ((-(1-yj)*Ul - (1-xi)*Vl)*x0 &
      + ((1-yj)*Ul - xi*Vl)*x1   &
      + (yj*Ul + xi*Vl)*x2       &
      + (-yj * Ul + (1-xi)*Vl)*x3)

v = ((-(1-yj)*Ul - (1-xi)*Vl)*y0 &
      + ((1-yj)*Ul - xi*Vl)*y1   &
      + (yj*Ul + xi*Vl)*y2       &
      + (-yj * Ul + (1-xi)*Vl)*y3)

jacden = jacob2D(xi, yj, grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y,grd%Lx, l0, l1, l2, l3, &
                                        explain=.False.)



dxidt = u/(jacden*dy_dlat)
detadt = v/(jacden*dy_dlat)


find_u  = dxidt*l0 + detadt*l1

!write(stderrunit,*) 'find_u: u, v, dxidt, detadt, find_u, l0, l1', u, v, dxidt, detadt, find_u, l0, l1
!write(stderrunit,*) 'find_u: u, Ul, Vl', u, Ul, Vl
end function find_u

real function linliny(grd,fld,x, y, i, j, xi, yj)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
real, intent(in) :: fld(grd%isd:grd%ied+1,grd%jsd:grd%jed+1), xi, yj
!real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
real :: linlinynum, linlinyden
! Local variables
integer :: stderrunit
  stderrunit=stderr()

    linliny = fld(i,j  )*yj + fld(i,j-1)*(1-yj)

end function linliny

real function find_v(grd,fld,x, y, i, j, xi, yj, dx_dlon, dy_dlat)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
real, intent(in) :: fld(grd%isd:grd%ied+1,grd%jsd:grd%jed+1), xi, yj
!real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
real, intent(in) :: dx_dlon, dy_dlat
real :: Ul, Vl, u, v
real :: linlinxnum, jacden
! Local variables
real :: x0 !< Longitude of first corner
real :: y0 !< Latitude of first corner
real :: x1 !< Longitude of second corner
real :: y1 !< Latitude of second corner
real :: x2 !< Longitude of third corner
real :: y2 !< Latitude of third corner
real :: x3 !< Longitude of fourth corner
real :: y3 !< Latitude of fourth corner
real :: xx, xx0,xx1,xx2,xx3
real :: l0, l1, l2, l3
real :: dxidt, detadt
real :: Lx
integer :: stderrunit
  stderrunit=stderr()

x0 = grd%lon(i-1,j-1)
x1 = grd%lon(i  ,j-1)
x2 = grd%lon(i  ,j  )
x3 = grd%lon(i-1,j  )

y0 = grd%lat(i-1,j-1)
y1 = grd%lat(i  ,j-1)
y2 = grd%lat(i  ,j  )
y3 = grd%lat(i-1,j  )

Lx = grd%Lx

xx= apply_modulo_around_point(x,x0,Lx)
xx0= apply_modulo_around_point(x0,x0,Lx)
xx1= apply_modulo_around_point(x1,x0,Lx)
xx2= apply_modulo_around_point(x2,x0,Lx)
xx3= apply_modulo_around_point(x3,x0,Lx)

Ul = linlinx(grd,fld,x, y, i, j, xi, yj)
Vl = linliny(grd,fld,x, y, i, j, xi, yj)
v = ((-(1-yj)*Ul - (1-xi)*Vl)*y0 &
      + ((1-yj)*Ul - xi*Vl)*y1   &
      + (yj*Ul + xi*Vl)*y2       &
      + (-yj * Ul + (1-xi)*Vl)*y3)

jacden = jacob2D(xi, yj, grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y,grd%Lx, l0, l1, l2, l3, &
                                       explain=.False.)

dxidt = u/(jacden*dx_dlon)
detadt = v/(jacden*dx_dlon)


find_v  = dxidt*l2 + detadt*l3

!write(stderrunit,*) 'find_v: u, v, dxidt, detadt, find_v', u, v, dxidt, detadt, find_v

end function find_v

real function jacob2D(xi, yj, x0, y0, x1, y1, x2, y2, x3, y3, x, y, Lx, l0, l1,l2, l3, explain)
! Get the stderr unit numberrguments
real, intent(in) :: xi
real, intent(in) :: yj
real, intent(in) :: x0 !< Longitude of first corner
real, intent(in) :: y0 !< Latitude of first corner
real, intent(in) :: x1 !< Longitude of second corner
real, intent(in) :: y1 !< Latitude of second corner
real, intent(in) :: x2 !< Longitude of third corner
real, intent(in) :: y2 !< Latitude of third corner
real, intent(in) :: x3 !< Longitude of fourth corner
real, intent(in) :: y3 !< Latitude of fourth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in) :: explain
real, intent(out) :: l0,l1,l2,l3
! Local variables
real :: xx
real :: xx0,xx1,xx2,xx3
integer :: stderrunit
  stderrunit=stderr()


  xx= apply_modulo_around_point(x,x0,Lx)
  xx0= apply_modulo_around_point(x0,x0,Lx)
  xx1= apply_modulo_around_point(x1,x0,Lx)
  xx2= apply_modulo_around_point(x2,x0,Lx)
  xx3= apply_modulo_around_point(x3,x0,Lx)

   l0 = xx0*(yj-1) + xx1*(1-yj) + xx2*yj-xx3*yj
   l1 = xx0*(xi-1)- xx1*xi + xx2*xi + xx3*(1-xi)
   l2 = y0*(yj-1) + y1*(1-yj) + y2*yj-y3*yj
   l3 = y0*(xi-1) - y1*xi + y2 *xi +y3*(1-xi)
   
  jacob2D = l0*l3 - l2*l3
end function jacob2D

! ##############################################################################

subroutine print_fld(grd, fld, label)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed)
character(len=*) :: label
! Local variables
integer :: i, j
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  write(stderrunit,'("pe=",i3,x,a8,32i10)') mpp_pe(),label,(i,i=grd%isd,grd%ied)
  do j=grd%jed,grd%jsd,-1
    write(stderrunit,'("pe=",i3,x,i8,32es10.2)') mpp_pe(),j,(fld(i,j),i=grd%isd,grd%ied)
  enddo

end subroutine print_fld

! ##############################################################################

!> Invoke some unit tests
logical function unit_tests(parts)
  type(particles), pointer:: parts !< Container for all types and memory
  ! Local variables
  type(particles_gridded), pointer :: grd
  integer :: stderrunit,i,j,c1,c2
  integer :: k1
  integer(kind=8) :: id
  real, dimension(parts%grd%ke) :: hdummy 
  real, dimension(parts%grd%ke) :: hddummy

  ! This function returns True is a unit test fails
  unit_tests=.false.
  ! For convenience
  grd=>parts%grd
 
 hdummy(1) = 100
 hddummy(1) = 100
  do k1 =2,grd%ke
      hdummy(k1) = 100
      hddummy(k1) = hddummy(k1-1)+100
  enddo

  stderrunit=stderr()

  i=grd%isc; j=grd%jsc
  call localTest( unit_tests, bilin(grd, grd%lon, i, j, 0., 1.), grd%lon(i-1,j) )
  call localTest( unit_tests, bilin(grd, grd%lon, i, j, 1., 1.), grd%lon(i,j) )
  call localTest( unit_tests, bilin(grd, grd%lat, i, j, 1., 0.), grd%lat(i,j-1) )
  call localTest( unit_tests, bilin(grd, grd%lat, i, j, 1., 1.), grd%lat(i,j) )
  if (grd%ke .gt. 4) then
      call localTest( unit_tests, find_layer1D(grd, 350.,hddummy), 3.5)
      call localTest( unit_tests, find_depth1D(grd,3.5,hdummy), 350.)
  endif

  ! Test 64-bit ID conversion
  i = 1440*1080 ; c1 = 2**30 + 2**4 + 1
  id = id_from_2_ints(c1, i)

end function unit_tests

! ##############################################################################

!> Returns true if non-dimensional position xi,yj is in unit interval
!!
!! Includes South and East boundaries, and excludes North and West.
!! \todo Double check definition of is_point_within_xi_yj_bounds()
logical function is_point_within_xi_yj_bounds(xi,yj)
! Arguments
real, intent(in) :: xi !< Non-dimensional x-position
real, intent(in) :: yj !< Non-dimensional y-position
! Local variables
!Includes South and East boundaries, and excludes North and West  (double check this is the way that is needed)
  is_point_within_xi_yj_bounds=.False.
  if ((xi .ge. 0 )  .and.  (xi .lt. 1)) then
    if ((yj .ge. 0 )  .and.  (yj .lt. 1)) then
      is_point_within_xi_yj_bounds=.True.
    endif
  endif
end function is_point_within_xi_yj_bounds

! ##############################################################################
!> Modulo value of x in an interval [y-(Lx/2)  y+(Lx/2)]
!!
!! Gives the modulo value of x in an interval [y-(Lx/2)  y+(Lx/2)]  , modulo Lx
!! If Lx<=0, then it returns x without applying modulo arithmetic.
real function apply_modulo_around_point(x, y, Lx)
! Arguments
real, intent(in) :: x !< Value to apply modulo arithmetic to
real, intent(in) :: y !< Center of modulo range
real, intent(in) :: Lx !< Modulo width
!Local_variables
real ::Lx_2

  if (Lx>0.) then
    Lx_2=Lx/2.
    apply_modulo_around_point=modulo(x-(y-Lx_2),Lx)+(y-Lx_2)
  else
    apply_modulo_around_point=x
  endif

end function apply_modulo_around_point

! #################################################################################
!> Combine a counter and i,j hash into an id
integer(kind=8) function id_from_2_ints(counter, ijhash)
  integer, intent(in) :: counter !< The counter value assigned at calving
  integer, intent(in) :: ijhash  !< A hash of i,j calving location

  id_from_2_ints = int(counter,8) * (int(2,8)**32) + int(ijhash,8)

end function id_from_2_ints

! #################################################################################

!> Split an particle ID into two parts
subroutine split_id(id, counter, ijhash)
  integer(kind=8), intent(in)  :: id      !< A unique id assigned when a part is created
  integer,         intent(out) :: counter !< The counter value assigned at calving
  integer,         intent(out) :: ijhash  !< A hash of i,j calving location
  ! Local variables
  integer(kind=8) :: i8

  counter = ishft(id,-32)
  !counter = i8
  ijhash = int(id,4)

end subroutine split_id

! #################################################################################
!> Checks answer to right answer and prints results if different
subroutine localTest(unit_test, answer, rightAnswer)
  logical, intent(inout) :: unit_test !< Set to true answer is wrong
  real, intent(in) :: answer !< Calculated answer
  real, intent(in) :: rightAnswer !< Correct answer
  ! Local variables
  integer :: stderrunit
  stderrunit=stderr()
  if (answer==rightAnswer) return
  unit_test=.true.
  write(stderrunit,*) 'a=',answer,'b=',rightAnswer
end subroutine localTest

!> Check for duplicates of particles on and across processors and issue an error
!! if any are detected
subroutine check_for_duplicates_in_parallel(parts)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  ! Local variables
  type(particles_gridded), pointer :: grd
  type(particle), pointer :: this
  integer :: stderrunit, i, j, k, nparts, nparts_total
  integer(kind=8), dimension(:), allocatable :: ids ! List of ids of all parts on this processor

  stderrunit=stderr()
  grd=>parts%grd
  nparts = count_parts(parts)
  nparts_total = nparts
  call mpp_sum(nparts_total) ! Total number of parts
  if (nparts_total==0) return ! Skip the rest of the test

  k = 0
  if (nparts>0) then
    allocate(ids(nparts))
    do j = grd%jsc,grd%jec ; do i = grd%isc,grd%iec
      this=>parts%list(i,j)%first
      do while (associated(this))
        k = k + 1
        ids(k) = this%id
        this=>this%next
      enddo
    enddo ; enddo
  endif
  if (k /= nparts) then
    write(stderrunit,*) 'counted parts=',k,'count_parts()=',nparts
    call error_mesg('particles, check_for_duplicates:', 'Mismatch between concatenation of lists and count_parts()!', FATAL)
  endif
  k = check_for_duplicate_ids_in_list(nparts, ids, verbose=.true.)
  if (k /= 0) call error_mesg('particles, check_for_duplicates:', 'Duplicate part detected across PEs!', FATAL)
  if (nparts>0) deallocate(ids)
end subroutine check_for_duplicates_in_parallel

! #################################################################################
!> Returns error count of duplicates of integer values in a distributed list
integer function check_for_duplicate_ids_in_list(nparts, ids, verbose)
  ! Arguments
  integer,                       intent(in)    :: nparts !< Length of ids
  integer(kind=8), dimension(:), intent(inout) :: ids !< List of ids
  logical,                       intent(in)    :: verbose !< True if messages should be written
  ! Local variables
  integer :: stderrunit, i, j, k, nparts_total, ii
  integer(kind=8) :: lowest_id, nonexistent_id, id, lid
  logical :: have_part

  stderrunit=stderr()
  nparts_total = nparts
  call mpp_sum(nparts_total) ! Total number of parts

  ! Establish lowest id or 0 across all PEs
  lowest_id = 0
  if (nparts>0) lowest_id = minval(ids)
  call mpp_min(lowest_id)
  id = lowest_id
  nonexistent_id = lowest_id - 1
  if (nonexistent_id >= lowest_id) then
    write(stderrunit,*) 'Underflow in particle ids!',nonexistent_id,lowest_id,mpp_pe()
  endif
  ! Sort the list "ids" (largest first)
  do j = 1, nparts-1
    do i = j+1, nparts
      if (ids(j) < ids(i)) then
        ! Swap
        id = ids(i)
        ids(i) = ids(j)
        ids(j) = id
      endif
    enddo
  enddo
  ! Check for duplicates on processor
  check_for_duplicate_ids_in_list = 0
  do k = 1, nparts-1
    if (ids(k) == ids(k+1)) then
      if (verbose) write(stderrunit,*) 'Duplicated part on PE with id=',ids(k),'pe=',mpp_pe()
      check_for_duplicate_ids_in_list = check_for_duplicate_ids_in_list + 1
    endif
  enddo
  ! Check for duplicates across processor
  j = 1 ! Pointer to first part in my list
  do k = 1, nparts_total
    ! Set i to first id in my list
    if (j <= nparts) then
      id = ids(j)
      have_part = .true.
    else
      id = nonexistent_id
      have_part = .false.
    endif
    lid = id
    call mpp_max(lid)
    if (have_part .and. id == lid) then
      ii = 1 ! This part is mine
      j = j + 1
    else
      ii = 0 ! This part is not mine
    endif
    call mpp_sum(ii)
    if (ii > 1) then
      if (verbose) write(stderrunit,*) 'Duplicated part across PEs with id=',id,lid,' seen',ii,' times pe=',mpp_pe(),k,j,nparts
      check_for_duplicate_ids_in_list = check_for_duplicate_ids_in_list + 1
    elseif (ii == 0) then
      if (verbose) write(stderrunit,*) 'part not accounted for on all PEs with id=',id,lid,' seen',ii,' times pe=',mpp_pe(),k,j,nparts
    endif
  enddo

end function check_for_duplicate_ids_in_list

! #################################################################################

!> Unit test for check_for_duplicate_ids_in_list()
subroutine test_check_for_duplicate_ids_in_list()
  ! Local variables
  integer :: k
  integer(kind=8), dimension(:), allocatable :: ids
  integer :: error_count

  allocate(ids(5))
  do k = 1,5
    ids(k) = k + 5*mpp_pe()
  enddo
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count /= 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('particles, test_check_for_duplicate_ids_in_list:', 'Unit test for clean list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = ids(4)
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('particles, test_check_for_duplicate_ids_in_list:', 'Unit test for dirty list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = 7 + 5*mpp_pe()
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('particles, test_check_for_duplicate_ids_in_list:', 'Unit test for a really dirty list failed!', FATAL)
  endif
  deallocate(ids)
end subroutine test_check_for_duplicate_ids_in_list

end module MOM_particles_framework
