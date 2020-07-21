subroutine read_fnc_mesh(pfFile,coeff_dx,offset,ncells,selection, &
                          cell_volumes,cell_coords,vertices_coords, sel_ncells)

  use netcdf
  use pf_params, only: pf_read_debug
  implicit none

  integer, intent(in) :: ncells
  integer, intent(in) :: sel_ncells
  real*8, intent(in) :: coeff_dx
  real*8, dimension(3), intent(in) :: offset
  integer*4, dimension(sel_ncells), intent(in) :: selection
  character(len=256),intent(in) :: pfFile
  real*8, dimension(ncells), intent(out) :: cell_volumes
  !f2py depend(selection) :: sel_ncells=len(selection)
  real*8, dimension(ncells,3), intent(out) :: cell_coords
  real*8, dimension(ncells*8,3), intent(out) :: vertices_coords

  ! Local variables
  integer :: ncid,ncerr,measid
  integer :: rank,k,nc,idx,idim
  integer, dimension(5) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  character(len=256) :: dim_name
  real*8 :: dx
  real*8, allocatable, dimension(:) :: buffer
  real*8, allocatable, dimension(:) :: voxel_scales

  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'

  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inq_varid(ncid, "coords", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  rank = 2

  ncerr = nf90_inquire_variable(ncid, measid, dimids = meas_dim_ids(1:rank))
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  allocate(idims(rank))

  do k=1,rank
    ncerr = nf90_inquire_dimension(ncid, meas_dim_ids(k), len = idims(k), name = dim_name)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo

  if ( idims(1)/=3 .OR. idims(2)/=ncells ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    write(*,*) 'idims(1)=',idims(1),3
    write(*,*) 'idims(2)=',idims(2),ncells
    stop
  endif

  allocate(start(rank))
  allocate(count(rank))
  allocate(buffer(ncells))

  start=(/1,1/)
  count=(/1,ncells/)

  ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  do k = 1,ncells
    cell_coords(k,1) = buffer(k) + offset(1)
  enddo

  start=(/2,1/)
  count=(/1,ncells/)

  ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  do k = 1,ncells
    cell_coords(k,2) = buffer(k) + offset(2)
  enddo

  start=(/3,1/)
  count=(/1,ncells/)

  ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  do k = 1,ncells
    cell_coords(k,3) = buffer(k) + offset(3)
  enddo

  deallocate(start)
  deallocate(idims)
  deallocate(count)
  deallocate(buffer)

  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inq_varid(ncid, "voxel_scales", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  rank = 1

  ncerr = nf90_inquire_variable(ncid, measid, dimids = meas_dim_ids(1:rank))
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  allocate(idims(rank))

  do k=1,rank
    ncerr = nf90_inquire_dimension(ncid, meas_dim_ids(k), len = idims(k), name = dim_name)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension', &
                              k,':',trim(dim_name),':',idims(k)
  enddo

  if ( idims(1)/=ncells ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    stop
  endif

  allocate(start(rank))
  allocate(count(rank))
  allocate(voxel_scales(ncells))

  start=(/1/)
  count=(/ncells/)

  ncerr = nf90_get_var(ncid, measid, voxel_scales,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ! Close file
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  deallocate(start)
  deallocate(idims)
  deallocate(count)

  do nc = 1, sel_ncells
    idx = selection(nc)+1
    dx = 2**(1.0 + voxel_scales(nc))
    cell_volumes(nc) = (dx * coeff_dx)**3

    vertices_coords(1+(nc-1)*8,1) = cell_coords(idx,1)
    vertices_coords(2+(nc-1)*8,1) = cell_coords(idx,1) + dx
    vertices_coords(3+(nc-1)*8,1) = cell_coords(idx,1) + dx
    vertices_coords(4+(nc-1)*8,1) = cell_coords(idx,1)
    vertices_coords(5+(nc-1)*8,1) = cell_coords(idx,1)
    vertices_coords(6+(nc-1)*8,1) = cell_coords(idx,1) + dx
    vertices_coords(7+(nc-1)*8,1) = cell_coords(idx,1) + dx
    vertices_coords(8+(nc-1)*8,1) = cell_coords(idx,1)

    vertices_coords(1+(nc-1)*8,2) = cell_coords(idx,2)
    vertices_coords(2+(nc-1)*8,2) = cell_coords(idx,2)
    vertices_coords(3+(nc-1)*8,2) = cell_coords(idx,2) + dx
    vertices_coords(4+(nc-1)*8,2) = cell_coords(idx,2) + dx
    vertices_coords(5+(nc-1)*8,2) = cell_coords(idx,2)
    vertices_coords(6+(nc-1)*8,2) = cell_coords(idx,2)
    vertices_coords(7+(nc-1)*8,2) = cell_coords(idx,2) + dx
    vertices_coords(8+(nc-1)*8,2) = cell_coords(idx,2) + dx

    vertices_coords(1+(nc-1)*8,3) = cell_coords(idx,3)
    vertices_coords(2+(nc-1)*8,3) = cell_coords(idx,3)
    vertices_coords(3+(nc-1)*8,3) = cell_coords(idx,3)
    vertices_coords(4+(nc-1)*8,3) = cell_coords(idx,3)
    vertices_coords(5+(nc-1)*8,3) = cell_coords(idx,3) + dx
    vertices_coords(6+(nc-1)*8,3) = cell_coords(idx,3) + dx
    vertices_coords(7+(nc-1)*8,3) = cell_coords(idx,3) + dx
    vertices_coords(8+(nc-1)*8,3) = cell_coords(idx,3) + dx

  end do

  do idim=1,3
    do nc=1,ncells
      cell_coords(nc,idim) = cell_coords(nc,idim) * coeff_dx
    enddo
  enddo


  deallocate(voxel_scales)

end subroutine read_fnc_mesh

subroutine read_fnc_frame(pfFile,frame_number,ncells,nnodes,selection,cell_volumes, &
  connectivity,read_indices,scale_types,data_node,sel_ncells,nvars)

  use netcdf
  use pf_params
  implicit none

  integer, intent(in) :: frame_number
  integer, intent(in) :: ncells
  integer, intent(in) :: nnodes
  integer, intent(in) :: sel_ncells
  integer, intent(in) :: nvars
  integer, dimension(nvars), intent(in) :: read_indices
  integer, dimension(nvars), intent(in) :: scale_types
  integer*4, dimension(sel_ncells), intent(in) :: selection
  real*8, dimension(sel_ncells), intent(in) :: cell_volumes
  integer*4, dimension(sel_ncells,8), intent(in) :: connectivity
  character(len=256),intent(in) :: pfFile
  real*8, dimension(nvars,nnodes), intent(out) :: data_node

  ! Local variables
  integer :: ncid,ncerr,measid
  integer :: rank,k,idx,iv,ivar
  integer*4 :: no, glo_cell,nc
  integer, dimension(5) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  real*8 :: eps, scale_nv, cell_weight, ivol,min_val,max_val,mean_val,mean2_val,std_val
  character(len=256) :: dim_name
  real*8, allocatable, dimension(:) :: tmp, tmp_s, volume_node
  real*8, allocatable, dimension(:,:) :: buffer
  real*8, allocatable, dimension(:,:) :: data_cell

  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inq_varid(ncid, "measurements", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inquire_variable(ncid, measid, ndims = rank)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  if (rank/=3) then
    write(*,*) "Expecting array of rank 3 for measurements"
    stop
  endif

  ncerr = nf90_inquire_variable(ncid, measid, dimids = meas_dim_ids(1:rank))
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  allocate(idims(rank))

  do k=1,rank
    ncerr = nf90_inquire_dimension(ncid, meas_dim_ids(k), len = idims(k), name = dim_name)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo

  if (idims(1)/=ncells) then
    write(*,'(A,1X,I9)') "Expecting Number of cells:",ncells
    stop
  endif

  allocate(start(rank))
  allocate(count(rank))

  allocate(buffer(idims(1),idims(2)))

  start=(/1,1,frame_number+1/)
  count=(/idims(1),idims(2),1/)

  ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  deallocate(start)
  deallocate(idims)
  deallocate(count)

  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  if (debug_scaling) write(*,*) 'allocate data cell'

  allocate(tmp(sel_ncells))
  allocate(tmp_s(sel_ncells))
  allocate(data_cell(sel_ncells,nvars))

  do ivar = 1,nvars

    if (debug_scaling) write(*,*) 'scaling ivar=',ivar

    ! C to F index
    idx = read_indices(ivar) + 1

    do nc = 1,sel_ncells
      glo_cell = selection(nc) + 1
      tmp(nc) = buffer(glo_cell,idx)
    enddo

    if (debug_scaling) write(*,*) 'tmp is ready'

    call scale_var(scale_types(ivar),tmp,tmp_s,sel_ncells)
    do nc = 1,sel_ncells
      data_cell(nc,ivar) = tmp_s(nc)
    enddo

    if (pf_read_debug) write(*,*) 'cell data is scaled'

  enddo

  if (debug_scaling) write(*,*) 'scatter to vertices'

  allocate(volume_node(nnodes))

  if (debug_scaling) write(*,*) 'Initialize data_node and volume_node'
  data_node(:,:) = 0.0d0
  volume_node(:) = 0.0d0

  scale_nv = 1.0d0/8.0d0
  eps = 1.0e5
  do iv = 1,8
    do nc = 1,sel_ncells
      no = connectivity(nc,iv) + 1
      cell_weight = cell_volumes(nc)*scale_nv
      do ivar = 1,nvars
        data_node(ivar,no) = data_node(ivar,no) + data_cell(nc,ivar)*cell_weight
      enddo
      volume_node(no) = volume_node(no) + cell_weight
    enddo
  enddo

  eps = minval(cell_volumes)*0.01
  if (debug_scaling) write(*,*) 'gather to node'
  if (debug_scaling) write(*,*) 'eps=',eps

  do no=1,nnodes
    ivol =  1.0d0 / max(volume_node(no),eps)
    do ivar = 1,nvars
      data_node(ivar,no) = data_node(ivar,no) * ivol
    enddo
  enddo

  if (debug_scaling) then
    do ivar = 1,nvars
      min_val =1.0d8
      max_val =-1.0d8
      mean_val = 0.0d0
      mean2_val = 0.0d0
      do no=1,nnodes
        min_val = min(data_node(ivar,no),min_val)
        max_val = max(data_node(ivar,no),max_val)
        mean_val = mean_val + data_node(ivar,no)
        mean2_val = mean2_val + data_node(ivar,no)**2
      enddo
      mean_val = mean_val/nnodes
      mean2_val = mean2_val/nnodes
      std_val = sqrt(mean2_val-mean_val*mean_val)
      write(*,'(A,I1,1X,A,E16.8,1X,A,E16.8,A)') 'ivar=',ivar,'(min=',min_val,' max=',max_val,')'
      write(*,'(A,I1,1X,A,E16.8,1X,A,E16.8,A)') 'ivar=',ivar,'(mean=',mean_val,'std=',std_val,')'
    enddo
  endif

  deallocate(data_cell)
  deallocate(buffer)
  deallocate(tmp)
  deallocate(tmp_s)


end subroutine read_fnc_frame

subroutine find_closest_cell(pt_coords,cell_coords,iclosest,dist_pt,ncells)

  implicit none

  integer :: ncells
  real*8, dimension(3), intent(in) :: pt_coords
  real*8, dimension(ncells,3), intent(in) :: cell_coords
  integer, intent(out) :: iclosest
  real*8, intent(out) :: dist_pt


  integer :: icell, i
  real*8 :: min_dist2, dist2

  min_dist2 = 1.0d20

  icell = 1
  do i = 1,ncells
    dist2 = (pt_coords(1)-cell_coords(i,1))**2                                             &
          + (pt_coords(2)-cell_coords(i,2))**2                                             &
          + (pt_coords(3)-cell_coords(i,3))**2
    if (min_dist2>dist2) then
      min_dist2 = dist2
      icell = i
    endif
  enddo

  dist_pt = sqrt(min_dist2)
  iclosest = icell - 1  ! F to C index

end subroutine find_closest_cell

subroutine get_cell_data(icell,nframes,pfFile,read_indices,scale_types,temporal_data,nvars)

  use netcdf
  use pf_params
  implicit none

  integer, intent(in) :: icell, nvars, nframes
  integer, dimension(nvars), intent(in) :: read_indices
  integer, dimension(nvars), intent(in) :: scale_types
  character(len=256),intent(in) :: pfFile
  real*8, dimension(nframes,nvars), intent(out) :: temporal_data

  ! Local variables
  integer :: ncid,ncerr,measid
  integer :: rank,k,idx,ivar,nt
  integer, dimension(5) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  character(len=256) :: dim_name
  real*8, allocatable, dimension(:) :: tmp, tmp_s
  real*8, allocatable, dimension(:,:) :: buffer

  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inq_varid(ncid, "measurements", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  ncerr = nf90_inquire_variable(ncid, measid, ndims = rank)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  if (rank/=3) then
    write(*,*) "Expecting array of rank 3 for measurements"
    stop
  endif

  ncerr = nf90_inquire_variable(ncid, measid, dimids = meas_dim_ids(1:rank))
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  allocate(idims(rank))

  do k=1,rank
    ncerr = nf90_inquire_dimension(ncid, meas_dim_ids(k), len = idims(k), name = dim_name)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo

  if (idims(3)/=nframes) then
    write(*,'(A,1X,I9)') "Expecting Number of frames:",nframes
    stop
  endif

  allocate(start(rank))
  allocate(count(rank))

  allocate(buffer(idims(2),idims(3)))
  idx = icell+1 ! C to F
  start=(/idx,1,1/)
  count=(/1,idims(2),idims(3)/)

  ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  deallocate(start)
  deallocate(idims)
  deallocate(count)

  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop

  allocate(tmp(nframes))
  allocate(tmp_s(nframes))

  do ivar = 1,nvars

    if (debug_scaling) write(*,*) 'scaling ivar=',ivar

    ! C to F index
    idx = read_indices(ivar) + 1

    do nt = 1,nframes
      tmp(nt) = buffer(idx,nt)
    enddo

    if (debug_scaling) write(*,*) 'tmp is ready'

    call scale_var(scale_types(ivar),tmp,tmp_s,nframes)

    do nt = 1,nframes
      temporal_data(nt,ivar) = tmp_s(nt)
    enddo

    if (debug_scaling) write(*,*) 'cell data is scaled'

  enddo

  deallocate(tmp)
  deallocate(tmp_s)

end subroutine get_cell_data
