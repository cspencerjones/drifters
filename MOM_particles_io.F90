module MOM_particles_io

use constants_mod, only: pi, omega, HLF

use mpp_domains_mod, only: domain2D
use mpp_domains_mod, only: mpp_domain_is_tile_root_pe,mpp_get_domain_tile_root_pe, mpp_define_io_domain
use mpp_domains_mod, only: mpp_get_tile_pelist,mpp_get_tile_npes,mpp_get_io_domain,mpp_get_tile_id

use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_gather, mpp_chksum
use mpp_mod, only: COMM_TAG_11, COMM_TAG_12, COMM_TAG_13, COMM_TAG_14

use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING, NOTE
use fms_mod, only: field_exist, file_exist, read_data, write_data

use fms_io_mod, only: get_instance_filename
use fms_io_mod, only : save_restart, restart_file_type, free_restart_type, set_meta_global
use fms_io_mod, only : register_restart_axis, register_restart_field, set_domain, nullify_domain
use fms_io_mod, only : read_unlimited_axis =>read_compressed, field_exist, get_field_size

use mpp_mod,    only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod,    only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use fms_mod,    only : clock_flag_default

use time_manager_mod, only: time_type, get_time, set_date, operator(-)
use MOM_time_manager,  only : time_type_to_real, real_to_time, get_ticks_per_second
use MOM_time_manager,  only : get_date

use MOM_grid, only : ocean_grid_type

use MOM_particles_framework, only: particles_gridded, xyt, particle, particles, buffer
use MOM_particles_framework, only: pack_traj_into_buffer2,unpack_traj_from_buffer2
use MOM_particles_framework, only: find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell,append_posn
use MOM_particles_framework, only: find_layer, find_depth
!use particles_framework, only: count_bonds, form_a_bond
use MOM_particles_framework, only: find_individual_particle
use MOM_particles_framework, only: push_posn
use MOM_particles_framework, only: add_new_part_to_list,destroy_particle
use MOM_particles_framework, only: increase_ibuffer,grd_chksum2,grd_chksum3
use MOM_particles_framework, only: bilin
!params !Niki: write a subroutine to get these
use MOM_particles_framework, only: buffer_width, buffer_width_traj
use MOM_particles_framework, only: verbose, really_debug, debug, restart_input_dir
use MOM_particles_framework, only: ignore_ij_restart, use_slow_find
use MOM_particles_framework, only: force_all_pes_traj
use MOM_particles_framework, only: check_for_duplicates_in_parallel
use MOM_particles_framework, only: split_id, id_from_2_ints, generate_id

implicit none ; private

include 'netcdf.inc'

public particles_io_init
public read_restart_parts, write_restart,write_trajectory

!Local Vars
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
!I/O vars
type(domain2d), pointer, save :: io_domain=>NULL()
integer, save :: io_tile_id(1), io_tile_root_pe, io_npes
integer, allocatable,save :: io_tile_pelist(:)
logical :: is_io_tile_root_pe = .true.

integer :: clock_trw,clock_trp

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

!> Initialize parallel i/o
subroutine particles_io_init(parts, io_layout)
type(particles), pointer :: parts !< particles container
integer, intent(in) :: io_layout(2) !< Decomposition of i/o processors

integer :: np
integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "particles_framework: "//trim(version)

  !I/O layout init
  io_tile_id=-1

  io_domain => mpp_get_io_domain(parts%grd%domain)

  if (.not. associated(io_domain)) then
    call mpp_define_IO_domain(parts%grd%domain, io_layout)
    io_domain => mpp_get_io_domain(parts%grd%domain)
  endif


  io_tile_id = mpp_get_tile_id(io_domain)
  is_io_tile_root_pe = mpp_domain_is_tile_root_pe(io_domain)
  io_tile_root_pe = mpp_get_domain_tile_root_pe(io_domain)
  np=mpp_get_tile_npes(io_domain)
  allocate(io_tile_pelist(np))
  call mpp_get_tile_pelist(io_domain,io_tile_pelist)
  io_npes = io_layout(1)*io_layout(2)

  clock_trw=mpp_clock_id( 'particles-traj write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  clock_trp=mpp_clock_id( 'particles-traj prepare', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

end subroutine particles_io_init

! ##############################################################################

!> Write an particle restart file
subroutine write_restart(parts, h, directory, time, stamp)
! Arguments
type(particles), pointer :: parts !< particles container
real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers
character(len=*), intent(in) :: directory
    !< The directory where the restart files are to be written 
type(time_type),          intent(in)    :: time       !< The current model time
logical,  intent(in) , optional   :: stamp !< If present and true, add time-stamp
! Local variables
!type(bond), pointer :: current_bond
integer :: i,j,id
character(len=35) :: filename
character(len=35) :: filename_bonds
type(particle), pointer :: this=>NULL()
integer :: stderrunit
!I/O vars
type(restart_file_type) :: parts_restart
type(restart_file_type) :: parts_bond_restart
integer :: nparts, nbonds
integer :: n_static_parts
logical :: check_bond_quality
character(len=40) :: restart_time
type(particles_gridded), pointer :: grd
real, allocatable, dimension(:) :: lon,          &
                                   lat,          &
                                   depth,        &
                                   uvel,         &
                                   vvel,         &
                                   start_lon,    &
                                   start_lat,    &
                                   start_d


integer, allocatable, dimension(:) :: ine,              &
                                      jne,              &
                                      drifter_num,      &
                                      id_cnt,           &
                                      id_ij,            &
                                      first_id_cnt,     &
                                      other_id_cnt,     &
                                      first_id_ij,      &
                                      other_id_ij,      &
                                      first_part_jne,         &
                                      first_part_ine,         &
                                      other_part_jne,         &
                                      other_part_ine


integer :: grdi, grdj

! Get the stderr unit number
 stderrunit=stderr()


  ! For convenience
  grd=>parts%grd

  !First add the parts on the io_tile_root_pe (if any) to the I/O list
  nparts = 0
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while (associated(this))
      nparts = nparts +1
      this=>this%next
    enddo
  enddo ; enddo


   allocate(lon(nparts))
   allocate(lat(nparts))
   allocate(depth(nparts))
   allocate(uvel(nparts))
   allocate(vvel(nparts))
   allocate(start_lon(nparts))
   allocate(start_lat(nparts))
   allocate(start_d(nparts))

   allocate(ine(nparts))
   allocate(jne(nparts))
   allocate(drifter_num(nparts))
   allocate(id_cnt(nparts))
   allocate(id_ij(nparts))

 
  call get_instance_filename("drifters.res.nc", filename)

 
  call set_domain(parts%grd%domain)

  call register_restart_axis(parts_restart,filename,'i',nparts)

  call set_meta_global(parts_restart,'file_format_major_version',ival=(/file_format_major_version/))
  call set_meta_global(parts_restart,'file_format_minor_version',ival=(/file_format_minor_version/))
  call set_meta_global(parts_restart,'time_axis',ival=(/0/))

  !Now start writing in the io_tile_root_pe if there are any parts in the I/O list

  ! Define Variables
  id = register_restart_field(parts_restart,filename,'lon',lon,longname='longitude',units='degrees_E')
  id = register_restart_field(parts_restart,filename,'lat',lat,longname='latitude',units='degrees_N')
  id = register_restart_field(parts_restart,filename,'depth',depth,longname='depth below surface',units='m')
  id = register_restart_field(parts_restart,filename,'uvel',uvel,longname='zonal velocity',units='m/s')
  id = register_restart_field(parts_restart,filename,'vvel',vvel,longname='meridional velocity',units='m/s')
  id = register_restart_field(parts_restart,filename,'ine',ine,longname='i index',units='none')
  id = register_restart_field(parts_restart,filename,'jne',jne,longname='j index',units='none')
  id = register_restart_field(parts_restart,filename,'start_lon',start_lon, &
                                            longname='longitude of starting location',units='degrees_E')
  id = register_restart_field(parts_restart,filename,'start_lat',start_lat, &
                                            longname='latitude of starting location',units='degrees_N')
  id = register_restart_field(parts_restart,filename,'start_d',start_d, &
                                            longname='depth of starting location',units='m')
  id = register_restart_field(parts_restart,filename,'drifter_num',drifter_num, &
                                            longname='identification of the drifter', units='dimensionless')
  id = register_restart_field(parts_restart,filename,'id_cnt',id_cnt, &
                                            longname='counter component of particle id', units='dimensionless')
  id = register_restart_field(parts_restart,filename,'id_ij',id_ij, &
                                            longname='position component of particle id', units='dimensionless')

  ! Write variables

  i = 0
  do grdj = parts%grd%jsc,parts%grd%jec ; do grdi = parts%grd%isc,parts%grd%iec
    this=>parts%list(grdi,grdj)%first
    do while(associated(this))
      i = i + 1
      call find_depth(grd,this%k,h,this%depth,this%ine,this%jne,this%xi,this%yj,this%k_space)
      lon(i) = this%lon; lat(i) = this%lat; depth(i) = this%depth
      uvel(i) = this%uvel; vvel(i) = this%vvel
      ine(i) = this%ine; jne(i) = this%jne
      start_lon(i) = this%start_lon; start_lat(i) = this%start_lat
      start_d(i) = this%depth
      id_cnt(i) = this%id; drifter_num(i) = this%drifter_num !; id_ij(i) = this%id
      call split_id(this%id, id_cnt(i), id_ij(i))
      this=>this%next
    enddo
  enddo ; enddo
  if (present(stamp)) then
    restart_time = convert_date_to_string(time)
    call save_restart(parts_restart, restart_time)
  else
    call save_restart(parts_restart)
  endif
  if (really_debug) print *, 'Finish save_restart.' !for debugging
  call free_restart_type(parts_restart)

  deallocate(              &
             lon,          &
             lat,          &
             depth,        &
             uvel,         &
             vvel,         &
             start_lon,    &
             start_lat,    &
             start_d) 


  deallocate(           &
             ine,       &
             jne,       &
             drifter_num,       &
             id_cnt,    &
             id_ij) 

  call nullify_domain()

end subroutine write_restart



!> This function converts a date into a string, valid with ticks and for dates
!up to year 99,999,999
character(len=40) function convert_date_to_string(date)
  type(time_type), intent(in) :: date !< The date to be translated into a string.
  
  ! Local variables
  character(len=40) :: sub_string
  real    :: real_secs
  integer :: yrs, mons, days, hours, mins, secs, ticks, ticks_per_sec
  character(len=40)    :: date_string

  call get_date(date, yrs, mons, days, hours, mins, secs, ticks)
  write (date_string, '(i8.4)') yrs
  write (sub_string, '("-", i2.2, "-", I2.2, ".", i2.2, ":", i2.2, ":")') &
         mons, days, hours, mins
  date_string = trim(adjustl(date_string)) // trim(sub_string)
  if (ticks > 0) then
    ticks_per_sec = get_ticks_per_second()
    real_secs = secs + ticks/ticks_per_sec
    if (ticks_per_sec <= 100) then
      write (sub_string, '(F7.3)') real_secs
    else
      write (sub_string, '(F10.6)') real_secs
    endif
  else
    write (sub_string, '(i2.2)') secs
  endif
  convert_date_to_string = trim(date_string) // trim(adjustl(sub_string))

end function convert_date_to_string

! ##############################################################################

subroutine read_restart_parts(parts,Time, u, v, h)
! Arguments
type(particles), pointer :: parts
type(time_type), intent(in) :: Time
real, dimension(:,:,:),intent(in) :: u, v
real, dimension(:,:,:),intent(in) :: h

!Local variables
integer :: n, siz(4), nparts_in_file, nparts_read
logical :: lres, found_restart, found, replace_drifter_num
logical :: explain
logical :: multiPErestart  ! Not needed with new restart read; currently kept for compatibility
real :: lon0, lon1, lat0, lat1
real :: pos_is_good, pos_is_good_all_pe
character(len=33) :: filename, filename_base
type(particles_gridded), pointer :: grd=>NULL()

type(particle) :: localpart
integer :: stderrunit, i, j, k, cnt, ij

real, allocatable,dimension(:) :: lon, &
                                  lat, &
                                  depth,  &
                                  drifter_num,  &
                                  id, &
                                  start_lon, &
                                  start_lat, &
                                  start_d
integer, allocatable, dimension(:) :: id_cnt, &
                                      id_ij
  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>parts%grd

  if (.not. associated(grd)) print *,'parts%grd not associated!'

  if (allocated(grd%uo)) deallocate(grd%uo)
  if (allocated(grd%vo)) deallocate(grd%vo)
  if (allocated(grd%hdepth)) deallocate(grd%hdepth)

  allocate(grd%uo(grd%isd:grd%ied,grd%jsd:grd%jed,grd%ke))
  allocate(grd%vo(grd%isd:grd%ied,grd%jsd:grd%jed,grd%ke))
  allocate(grd%hdepth(grd%isd:grd%ied,grd%jsd:grd%jed,grd%ke))

  do k=1,grd%ke
    do j=grd%jsd,grd%jed
      do i=grd%isd,grd%ied
         grd%uo(i,j,k) = u(i,j,k)
       enddo
    enddo
  enddo

  do k=1,grd%ke
    do j=grd%jsd,grd%jed
      do i=grd%isd,grd%ied
         grd%vo(i,j,k) = v(i,j,k)
       enddo
    enddo
  enddo

  do k=1,grd%ke
    do j=grd%jsd,grd%jed
      do i=grd%isd,grd%ied
         if (k.eq.1)then
             grd%hdepth(i,j,k) = h(i,j,k)
         else
             grd%hdepth(i,j,k) = h(i,j,k)+grd%hdepth(i,j,k-1)
         endif
       enddo
    enddo
  enddo

  ! Zero out nparts_in_file
  nparts_in_file = 0

  filename_base=trim(restart_input_dir)//'drifters.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))

  if (found_restart) then
    filename = filename_base
    call get_field_size(filename,'i',siz, field_found=found, domain=grd%domain)

    nparts_in_file = siz(1)
    replace_drifter_num = field_exist(filename, 'drifter_num') ! True if using a 32-bit drifter_num in restart file
    allocate(lon(nparts_in_file))
    allocate(lat(nparts_in_file))
    allocate(start_lon(nparts_in_file))
    allocate(start_lat(nparts_in_file))
    allocate(start_d(nparts_in_file))
    allocate(depth(nparts_in_file))
    if (replace_drifter_num) then
      allocate(id(nparts_in_file))
      allocate(drifter_num(nparts_in_file))
    else
      allocate(id_cnt(nparts_in_file))
      allocate(id_ij(nparts_in_file))
    endif

    call read_unlimited_axis(filename,'lon',lon,domain=grd%domain)
    call read_unlimited_axis(filename,'lat',lat,domain=grd%domain)
    call read_unlimited_axis(filename,'start_lon',start_lon,domain=grd%domain)
    call read_unlimited_axis(filename,'start_lat',start_lat,domain=grd%domain)
    call read_unlimited_axis(filename,'start_d',start_d,domain=grd%domain)
    call read_unlimited_axis(filename,'depth',depth,domain=grd%domain)
    if (replace_drifter_num) then
      call read_unlimited_axis(filename,'drifter_num',id,domain=grd%domain)
      call read_unlimited_axis(filename,'drifter_num',drifter_num,domain=grd%domain)
    else
      call read_int_vector(filename, 'id_cnt', id_cnt, grd%domain)
      call read_int_vector(filename, 'id_ij', id_ij, grd%domain)
    endif
  end if ! found_restart ln 569

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )

  do n=1, nparts_in_file
    localpart%lon=lon(n)
    localpart%lat=lat(n)
    localpart%depth=depth(n)

    if (use_slow_find) then
      lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    else
      lres=find_cell_by_search(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    endif

    if (really_debug) then
      write(stderrunit,'(a,i8,a,3f9.4,a,i8)') 'particles, read_restart_part: part ',n,' is at ',localpart%lon,localpart%lat,localpart%depth,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'particles, read_restart_parts: lres = ',lres
    endif

    if (lres) then ! True if the particle resides on the current processors computational grid
      localpart%start_lon=start_lon(n)
      localpart%start_lat=start_lat(n)
      localpart%start_d=start_d(n)
      if (replace_drifter_num) then
        localpart%id = generate_id(grd, localpart%ine, localpart%jne)
        localpart%drifter_num = drifter_num(n)
      else
        localpart%id = id_from_2_ints(id_cnt(n), id_ij(n))
      endif
      localpart%halo_part=0.
      lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
      localpart%k_space=.false.
      call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart)
    endif
  end do ! ln 650

  if (found_restart) then
    deallocate(lon,          &
               lat,          &
               depth,        &
               start_lon,    &
               start_lat,    &
               start_d)
    if (replace_drifter_num) then
      deallocate(id)
      deallocate(drifter_num)
    else
      deallocate(id_cnt)
      deallocate(id_ij)
    endif
  endif

  call check_for_duplicates_in_parallel(parts)

end subroutine read_restart_parts

!> Read a vector of integers from file and use a default value if variable is missing
subroutine read_int_vector(filename, varname, values, domain, value_if_not_in_file)
  character(len=*),  intent(in)  :: filename !< Name of file to read from
  character(len=*),  intent(in)  :: varname !< Name of variable to read
  integer,           intent(out) :: values(:) !< Returned vector of integers
  type(domain2D),    intent(in)  :: domain !< Parallel decomposition
  integer, optional, intent(in)  :: value_if_not_in_file !< Value to use if variable is not in file

  if (present(value_if_not_in_file).and..not.field_exist(filename, varname)) then
    values(:)=value_if_not_in_file
  else
    call read_unlimited_axis(filename,varname,values,domain=domain)
  endif
end subroutine read_int_vector

! ##############################################################################
!> Write a trajectory-based diagnostics file
subroutine write_trajectory(trajectory, save_short_traj)
! Arguments
type(xyt), pointer :: trajectory !< An particle trajectory
logical, intent(in) :: save_short_traj !< If true, record less data
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid, idcntid, idijid, drnumid
integer :: kid, depthid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid, sssid
integer :: cnid, hiid
integer :: mid, did, wid, lid, mbid, hdid
character(len=37) :: filename
character(len=7) :: pe_name
type(xyt), pointer :: this, next
integer :: stderrunit, cnt, ij
!I/O vars
type(xyt), pointer :: traj4io=>null()
integer :: ntrajs_sent_io,ntrajs_rcvd_io
integer :: from_pe,np
type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
logical :: io_is_in_append_mode

  ! Get the stderr unit number
  stderrunit=stderr()
  traj4io=>null()
  obuffer_io=>null()
  ibuffer_io=>null()

  !Assemble the list of trajectories from all pes in this I/O tile
  call mpp_clock_begin(clock_trp)

  !First add the trajs on the io_tile_root_pe (if any) to the I/O list
  if(is_io_tile_root_pe .OR. force_all_pes_traj ) then
     if(associated(trajectory)) then
        this=>trajectory
        do while (associated(this))
           call append_posn(traj4io, this)
           this=>this%next
        enddo
        trajectory => null()
     endif
  endif

  if(.NOT. force_all_pes_traj ) then

  !Now gather and append the parts from all pes in the io_tile to the list on corresponding io_tile_root_pe
  ntrajs_sent_io =0
  ntrajs_rcvd_io =0
  if(is_io_tile_root_pe) then
     !Receive trajs from all pes in this I/O tile !FRAGILE!SCARY!
     do np=2,size(io_tile_pelist) ! Note: np starts from 2 to exclude self
        from_pe=io_tile_pelist(np)
        call mpp_recv(ntrajs_rcvd_io, glen=1, from_pe=from_pe, tag=COMM_TAG_11)
        if (ntrajs_rcvd_io .gt. 0) then
           call increase_ibuffer(ibuffer_io, ntrajs_rcvd_io,buffer_width_traj)
           call mpp_recv(ibuffer_io%data, ntrajs_rcvd_io*buffer_width_traj,from_pe=from_pe, tag=COMM_TAG_12)
           do i=1, ntrajs_rcvd_io
                 call unpack_traj_from_buffer2(traj4io, ibuffer_io, i)
           enddo
       endif
     enddo
  else
     ! Pack and send trajectories to the root PE for this I/O tile
     do while (associated(trajectory))
       ntrajs_sent_io = ntrajs_sent_io +1
       call pack_traj_into_buffer2(trajectory, obuffer_io, ntrajs_sent_io)
       this => trajectory ! Need to keep pointer in order to free up the links memory
       trajectory => trajectory%next ! This will eventually result in trajectory => null()
       deallocate(this) ! Delete the link from memory
     enddo

     call mpp_send(ntrajs_sent_io, plen=1, to_pe=io_tile_root_pe, tag=COMM_TAG_11)
     if (ntrajs_sent_io .gt. 0) then
        call mpp_send(obuffer_io%data, ntrajs_sent_io*buffer_width_traj, to_pe=io_tile_root_pe, tag=COMM_TAG_12)
     endif
  endif

  endif !.NOT. force_all_pes_traj

  call mpp_clock_end(clock_trp)

  !Now start writing in the io_tile_root_pe if there are any parts in the I/O list
  call mpp_clock_begin(clock_trw)

  if((force_all_pes_traj .OR. is_io_tile_root_pe) .AND. associated(traj4io)) then

    call get_instance_filename("drifter_trajectories.nc", filename)
    if(io_tile_id(1) .ge. 0 .AND. .NOT. force_all_pes_traj) then !io_tile_root_pes write
       if(io_npes .gt. 1) then !attach tile_id  to filename only if there is more than one I/O pe
          if (io_tile_id(1)<10000) then
             write(filename,'(A,".",I4.4)') trim(filename), io_tile_id(1)
          else
             write(filename,'(A,".",I6.6)') trim(filename), io_tile_id(1)
          endif
       endif
    else !All pes write, attach pe# to filename
       if (mpp_npes()<10000) then
          write(filename,'(A,".",I4.4)') trim(filename), mpp_pe()
       else
          write(filename,'(A,".",I6.6)') trim(filename), mpp_pe()
       endif
    endif
    io_is_in_append_mode = .false.
    iret = nf_create(filename, NF_NOCLOBBER, ncid)
    if (iret .ne. NF_NOERR) then
      iret = nf_open(filename, NF_WRITE, ncid)
      io_is_in_append_mode = .true.
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'particles, write_trajectory: nf_open failed'
    endif
    if (verbose) then
      if (io_is_in_append_mode) then
        write(*,'(2a)') 'particles, write_trajectory: appending to ',filename
      else
        write(*,'(2a)') 'particles, write_trajectory: creating ',filename
      endif
    endif

    if (io_is_in_append_mode) then
      iret = nf_inq_dimid(ncid, 'i', i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'particles, write_trajectory: nf_inq_dimid i failed'
      lonid = inq_varid(ncid, 'lon')
      latid = inq_varid(ncid, 'lat')
      kid = inq_varid(ncid, 'k')
      depthid = inq_varid(ncid, 'depth')
      yearid = inq_varid(ncid, 'year')
      dayid = inq_varid(ncid, 'day')
      drnumid = inq_varid(ncid, 'drifter_num')
      idcntid = inq_varid(ncid, 'id_cnt')
      idijid = inq_varid(ncid, 'id_ij')
      if (.not.save_short_traj) then
        uvelid = inq_varid(ncid, 'uvel')
        vvelid = inq_varid(ncid, 'vvel')
      endif
    else
      ! Dimensions
      iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'particles, write_trajectory: nf_def_dim i failed'

      ! Variables
      lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
      latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
      kid = def_var(ncid,'k', NF_DOUBLE, i_dim)
      depthid = def_var(ncid,'depth', NF_DOUBLE, i_dim)
      yearid = def_var(ncid, 'year', NF_INT, i_dim)
      dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
      drnumid = def_var(ncid, 'drifter_num', NF_INT, i_dim)
      idcntid = def_var(ncid, 'id_cnt', NF_INT, i_dim)
      idijid = def_var(ncid, 'id_ij', NF_INT, i_dim)
      if (.not. save_short_traj) then
        uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
        vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
      endif
      ! Attributes
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_major_version', NF_INT, 1, 0)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_minor_version', NF_INT, 1, 1)
      call put_att(ncid, lonid, 'long_name', 'longitude')
      call put_att(ncid, lonid, 'units', 'degrees_E')
      call put_att(ncid, latid, 'long_name', 'latitude')
      call put_att(ncid, latid, 'units', 'degrees_N')
      call put_att(ncid, kid, 'long_name', 'k')
      call put_att(ncid, kid, 'units', 'layer number')
      call put_att(ncid, depthid, 'long_name', 'depth')
      call put_att(ncid, depthid, 'units', 'm')
      call put_att(ncid, yearid, 'long_name', 'year')
      call put_att(ncid, yearid, 'units', 'years')
      call put_att(ncid, dayid, 'long_name', 'year day')
      call put_att(ncid, dayid, 'units', 'days')
      call put_att(ncid, drnumid, 'long_name', 'identification of the drifter')
      call put_att(ncid, drnumid, 'units', 'dimensionless')
      call put_att(ncid, idcntid, 'long_name', 'counter component of particle id')
      call put_att(ncid, idcntid, 'units', 'dimensionless')
      call put_att(ncid, idijid, 'long_name', 'position component of particle id')
      call put_att(ncid, idijid, 'units', 'dimensionless')

      if (.not. save_short_traj) then
        call put_att(ncid, uvelid, 'long_name', 'zonal spped')
        call put_att(ncid, uvelid, 'units', 'm/s')
        call put_att(ncid, vvelid, 'long_name', 'meridional spped')
        call put_att(ncid, vvelid, 'units', 'm/s')
      endif
    endif

    ! End define mode
    iret = nf_enddef(ncid)

    ! Write variables
    this=>traj4io
    if (io_is_in_append_mode) then
      iret = nf_inq_dimlen(ncid, i_dim, i)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'particles, write_trajectory: nf_inq_dimlen i failed'
    else
      i = 0
    endif
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon)
      call put_double(ncid, latid, i, this%lat)
      call put_double(ncid, kid, i, this%k)
      call put_double(ncid, depthid, i, this%depth)
      call put_int(ncid, yearid, i, this%year)
      call put_double(ncid, dayid, i, this%day)
      call put_int(ncid, drnumid, i, this%particle_num)
      call split_id(this%id, cnt, ij)
      call put_int(ncid, idcntid, i, cnt)
      call put_int(ncid, idijid, i, ij)
      if (.not. save_short_traj) then
        call put_double(ncid, uvelid, i, this%uvel)
        call put_double(ncid, vvelid, i, this%vvel)
      endif
      next=>this%next
      deallocate(this)
      this=>next
    enddo

    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'particles, write_trajectory: nf_close failed',mpp_pe(),filename

  endif !(is_io_tile_root_pe .AND. associated(traj4io))
  call mpp_clock_end(clock_trw)


this=>trajectory


end subroutine write_trajectory


! ##############################################################################

integer function inq_var(ncid, var, unsafe)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
logical, optional, intent(in) :: unsafe
! Local variables
integer :: iret
integer :: stderrunit
logical :: unsafely=.false.

if(present(unsafe)) unsafely=unsafe
  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
    if (.not. unsafely) then
      write(stderrunit,*) 'particles, inq_var: nf_inq_varid ',var,' failed'
      call error_mesg('particles, inq_var', 'netcdf function returned a failure!', FATAL)
    else
      inq_var=-1
    endif
  endif

end function inq_var

! ##############################################################################

integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid, ntype, idim
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, def_var: nf_def_var failed for ',trim(var)
    call error_mesg('particles, def_var', 'netcdf function returned a failure!', FATAL)
  endif

end function def_var

! ##############################################################################

integer function inq_varid(ncid, var)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_inq_varid(ncid, var, inq_varid)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, inq_varid: nf_inq_varid failed for ',trim(var)
    call error_mesg('particles, inq_varid', 'netcdf function returned a failure!', FATAL)
  endif

end function inq_varid

! ##############################################################################

subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid, id
character(len=*), intent(in) :: att, attval
! Local variables
integer :: vallen, iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  vallen=len_trim(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('particles, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

! ##############################################################################

real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, get_double: nf_get_var1_double failed reading'
    call error_mesg('particles, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

! ##############################################################################

integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, get_int: nf_get_var1_int failed reading'
    call error_mesg('particles, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

! ##############################################################################
!> Write a real to a netcdf file
subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index of position to write
real, intent(in) :: val !< Value to write
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, put_double: nf_put_vara_double failed writing'
    call error_mesg('particless, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

! ##############################################################################

!> Write an integer to a netcdf file
subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index of position to write
integer, intent(in) :: val !< Value to write
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'particles, put_int: nf_put_vara_int failed writing'
    call error_mesg('particles, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int


! ##############################################################################

logical function find_restart_file(filename, actual_file, multiPErestart, tile_id)
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: actual_file
  logical, intent(out) :: multiPErestart
  integer, intent(in) :: tile_id

  character(len=6) :: pe_name

  find_restart_file = .false.

  ! If running as ensemble, add the ensemble id string to the filename
  call get_instance_filename(filename, actual_file)

  ! Prefer combined restart files.
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) return

  ! Uncombined restart
  if(tile_id .ge. 0) then
    write(actual_file,'(A,".",I4.4)') trim(actual_file), tile_id
  else
  if (mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()
  endif
  actual_file=trim(actual_file)//trim(pe_name)
  endif
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) then
     multiPErestart=.true.
     return
  endif

  ! No file found, Reset all return parameters
  find_restart_file=.false.
  actual_file = ''
  multiPErestart=.false.

end function find_restart_file


!######################################################################################

end module MOM_particles_io
