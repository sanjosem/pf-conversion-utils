subroutine read_fnc_mesh(pfFile,coeff_dx,offset,ncells,selection, &
                          cell_volumes,cell_coords,vertices_coords)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: ncells
  integer, intent(in) :: sel_ncells
  real*8, intent(in) :: coeff_dx
  real*8, dimension(3), intent(in) :: offset
  integer*4, dimension(sel_ncells), intent(in) :: selection
  character(len=256),intent(in) :: pfFile 
  real*8, dimension(ncells), intent(out) :: cell_volumes
  real*8, dimension(ncells,3), intent(out) :: cell_coords
  real*8, dimension(ncells*8,3), intent(out) :: vertices_coords

  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nc,idx,idim
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  character(len=256) :: dim_name
  real*8 :: dx
  real*8, allocatable, dimension(:) :: buffer
  real*8, allocatable, dimension(:) :: voxel_scales
  
  pf_read_debug = .TRUE.
  
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
  
    vertices_coords(1+nc,1) = cell_coords(idx,1)
    vertices_coords(2+nc,1) = cell_coords(idx,1) + dx
    vertices_coords(3+nc,1) = cell_coords(idx,1) + dx
    vertices_coords(4+nc,1) = cell_coords(idx,1)
    vertices_coords(5+nc,1) = cell_coords(idx,1)
    vertices_coords(6+nc,1) = cell_coords(idx,1) + dx
    vertices_coords(7+nc,1) = cell_coords(idx,1) + dx
    vertices_coords(8+nc,1) = cell_coords(idx,1)
  
    vertices_coords(1+nc,2) = cell_coords(idx,2)
    vertices_coords(2+nc,2) = cell_coords(idx,2)
    vertices_coords(3+nc,2) = cell_coords(idx,2) + dx
    vertices_coords(4+nc,2) = cell_coords(idx,2) + dx
    vertices_coords(5+nc,2) = cell_coords(idx,2)
    vertices_coords(6+nc,2) = cell_coords(idx,2)
    vertices_coords(7+nc,2) = cell_coords(idx,2) + dx
    vertices_coords(8+nc,2) = cell_coords(idx,2) + dx
  
    vertices_coords(1+nc,3) = cell_coords(idx,3)
    vertices_coords(2+nc,3) = cell_coords(idx,3)
    vertices_coords(3+nc,3) = cell_coords(idx,3)
    vertices_coords(4+nc,3) = cell_coords(idx,3)
    vertices_coords(5+nc,3) = cell_coords(idx,3) + dx
    vertices_coords(6+nc,3) = cell_coords(idx,3) + dx
    vertices_coords(7+nc,3) = cell_coords(idx,3) + dx
    vertices_coords(8+nc,3) = cell_coords(idx,3) + dx
  
  end do
  
  do idim=1,3
    do nc=1,ncells
      cell_coords(nc,idim) = cell_coords(nc,idim) * coeff_dx
    enddo
  enddo
  
  
  deallocate(voxel_scales)

end subroutine read_fnc_mesh

subroutine read_fnc_frame(pfFile,frame_number,ncells,selection,connectivity,data_node)
  
  use netcdf
  use pf_params
  implicit none
  
  integer, intent(in) :: frame_number
  integer, intent(in) :: ncells
  integer, intent(in) :: sel_ncells
  real*8, intent(in) :: coeff_dx
  real*8, dimension(3), intent(in) :: offset
  integer*4, dimension(sel_ncells), intent(in) :: selection
  real*8, dimension(sel_ncells), intent(in) :: cell_volumes
  real*8, dimension(sel_ncells,8), intent(in) :: connectivity
  character(len=256),intent(in) :: pfFile 
  real*8, dimension(nvars,nnodes), intent(out) :: data_node
  
  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nc,idx,idim
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  real*8 :: eps
  real*8, allocatable, dimension(:,:) :: buffer
  integer, allocatable, dimension(:) :: retrieve_index
  character(len=64), allocatable, dimension(:) :: names
  
  pf_read_debug = .TRUE.
  
  ncerr = nf90_open(filename, NF90_NOWRITE, ncid)
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
    WRITE(message,'(A,1X,I9)') "Expecting Number of cells:",ncells
    write(*,*) message
    stop
  endif
    
  allocate(start(rank))
  allocate(count(rank))
  
  allocate(buffer(idims(1),idims(2)))
  allocate(tmp(sel_ncells))
  allocate(data_cell(sel_ncells,nvars))

  start=(/1,1,frame_number/)
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

  do ivar = 1,nvars
    retrieve_index(ivar) = var_indices(ivar)+1
    names(ivar) = var_names(ivar)
    if (TRIM(names(ivar))=='static_pressure') ipress = retrieve_index(ivar) 
    if (TRIM(names(ivar))=='density') irho = retrieve_index(ivar) 
  enddo
  
  do ivar = 1,nvars
    if (retrieve_index(ivar)<=0) then
      if (TRIM(names(ivar))=='static_pressure') then
        retrieve_index(ivar) = irho
        names(ivar) = 'pressure_from_density'
      elseif (TRIM(names(ivar))=='density') then
        retrieve_index(ivar) = ipress
        names(ivar) = 'density_from_pressure'
      endif
    endif
  enddo
  
  do ivar = 1,nvars

    idx = retrieve_index(ivar)
    
    do nc = 1,sel_ncells
      gloc_cell = selection(nc) + 1
      tmp(nc) = buffer(gloc_cell,idx)
    enddo
    
    call scale_var(names(ivar),tmp,data_cell(:,ivar))
      
  enddo
  
  deallocate(retrieve_index)
  deallocate(names)
  
  scale_nv = 1.0d0/8.0d0
  eps = 1.0e5
  do iv = 1,8
    do nc = 1,sel_ncells
      no = connectivity(nc,iv) + 1
      cell_weight = volume_cell(nc)*scale_nv
      eps = min(eps,0.1d0*cell_weight)
      do ivar = 1,nvar
        data_node(ivar,no) = data_node(ivar,no) + data_cell(nc,ivar)*cell_weight
      enddo
      volume_node(no) = volume_node(no) + cell_weight
    enddo
  enddo
  
  do no=1,nnodes
    volume_node(no) = 1.0d0 / max(volume_node(no),eps)
  enddo
  
  do no=1,nnodes
    do ivar = 1,nvar
      data_node(ivar,no) = data_node(ivar,no) * volume_node(no)
    enddo
  enddo
  
  deallocate(data_cell)
  deallocate(buffer)
  deallocate(tmp)
  

end subroutine read_fnc_frame
