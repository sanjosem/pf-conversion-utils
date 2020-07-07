subroutine read_snc_coordinates(pfFile,coeff_dx,offset,nnodes,coords)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: nnodes
  real*8, intent(in) :: coeff_dx
  real*8, dimension(3), intent(in) :: offset
  character(len=256),intent(in) :: pfFile 
  real*8, dimension(nnodes,3), intent(out) :: coords

  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nf
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  character(len=256) :: dim_name
  real*8, dimension(3) :: minv,maxv
  real*8, allocatable, dimension(:) :: buffer
  
  pf_read_debug = .TRUE.
  
  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'
  
  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "vertex_coords", measid)
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
  
  if ( idims(1)/=3 .OR. idims(2)/=nnodes ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    write(*,*) 'idims(1)=',idims(1),3
    write(*,*) 'idims(2)=',idims(2),nnodes
    stop
  endif
  
  allocate(start(rank))
  allocate(count(rank))
  allocate(buffer(nnodes))
  
  do k = 1,3
    start=(/k,1/)
    count=(/1,nnodes/)
    
    ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    
    do nf = 1,nnodes
      coords(nf,k) = (buffer(nf) + offset(k)) * coeff_dx
    enddo
    
  enddo
  
  deallocate(start)
  deallocate(idims)
  deallocate(count)
  deallocate(buffer)
  
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  minv = minval(coords,dim=1)
  maxv = maxval(coords,dim=1)
  
  write(*,'(A)') 'Bounding box in meters:'
  write(*,'(A,1X,E10.3,1X,E10.3)') '  x: ',minv(1),maxv(1)
  write(*,'(A,1X,E10.3,1X,E10.3)') '  y: ',minv(2),maxv(2)
  write(*,'(A,1X,E10.3,1X,E10.3)') '  z: ',minv(3),maxv(3)

end subroutine read_snc_coordinates

subroutine compute_vertex_nb(pfFile,nfaces,nvertices,min_nv,max_nv,vert_per_face)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: nfaces, nvertices 
  integer, intent(out) :: min_nv,max_nv
  character(len=256),intent(in) :: pfFile 
  integer, dimension(nfaces), intent(out) :: vert_per_face

  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nf
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims
  character(len=256) :: dim_name
  integer, allocatable, dimension(:) :: buffer
  
  pf_read_debug = .TRUE.
  
  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'
  
  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "first_vertex_refs", measid)
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
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo
  
  if ( idims(1)/=nfaces ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    write(*,*) 'idims(1)=',idims(1),nfaces
    stop
  endif
  
  allocate(buffer(nfaces))
        
  ncerr = nf90_get_var(ncid, measid, buffer)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  deallocate(idims)
  
  do nf=1,nfaces-1
    vert_per_face(nf) = buffer(nf+1) - buffer(nf) ! C numbering
  enddo
  vert_per_face(nfaces) = nvertices - buffer(nfaces)

  deallocate(buffer)
  
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  min_nv = minval(vert_per_face)
  max_nv = maxval(vert_per_face)
  
  write(*,'(A)') 'Reading connectivity:'
  write(*,'(A,1X,I10)') '  Total number of faces: ',nfaces
  write(*,'(A,1X,I3,1X,I3)') '  Min/Max vertices per faces: ',min_nv, max_nv

end subroutine compute_vertex_nb


subroutine create_initial_connectivity(pfFile,vert_per_face,max_nv,nfaces,face_to_node)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: max_nv 
  integer, intent(in) :: nfaces
!f2py optional , depend(vert_per_face) :: nfaces=len(vert_per_face)
  character(len=256),intent(in) :: pfFile 
  integer, dimension(nfaces), intent(in) :: vert_per_face
  integer, dimension(nfaces,max_nv), intent(out) :: face_to_node

  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nf,nv,nvert, idv
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims
  character(len=256) :: dim_name
  integer, allocatable, dimension(:) :: vertex_to_node, first_vertex
  
  pf_read_debug = .TRUE.
  
  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'
  
  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "first_vertex_refs", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  allocate(first_vertex(nfaces))
        
  ncerr = nf90_get_var(ncid, measid, first_vertex)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "vertex_refs", measid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  rank = 1
  allocate(idims(rank))
  
  ncerr = nf90_inquire_variable(ncid, measid, dimids = meas_dim_ids(1:rank))
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
    
  do k=1,rank
    ncerr = nf90_inquire_dimension(ncid, meas_dim_ids(k), len = idims(k), name = dim_name)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo
  
  allocate(vertex_to_node(idims(1)))
        
  ncerr = nf90_get_var(ncid, measid, vertex_to_node)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  deallocate(idims)
  
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
    
  do nf=1,nfaces
    nvert = min(vert_per_face(nf),max_nv)
    do nv = 1,nvert
      idv = first_vertex(nf)+nv ! C to F numbering 
      face_to_node(nf,nv) = vertex_to_node(idv) ! keep C numbering
    enddo
  enddo
  
  deallocate(vertex_to_node)
  deallocate(first_vertex)

end subroutine create_initial_connectivity


subroutine compute_face_node_weight(pfFile,coeff_dx,vert_per_face,face_to_node,nnodes,nfaces, &
  max_nv,face_weight, node_weight,face_area)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: nfaces
  integer, intent(in) :: nnodes
  integer, intent(in) :: max_nv
  real*8, intent(in) :: coeff_dx
  character(len=256),intent(in) :: pfFile 
  integer, dimension(nfaces), intent(in) :: vert_per_face
  integer, dimension(nfaces,max_nv), intent(in) :: face_to_node
  real*8, dimension(nfaces), intent(out) :: face_weight, face_area
  real*8, dimension(nnodes), intent(out) :: node_weight

  ! Local variables
  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nf,nv,glo,no
  real*8 :: scale_nv, scale_area
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims
  character(len=256) :: dim_name
  
  pf_read_debug = .TRUE.
  
  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'
  
  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "surfel_area", measid)
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
    if (pf_read_debug) write(*,'(A,1X,I3,1X,A,1X,A,1X,A,I9)') 'dimension',k,':',trim(dim_name),':',idims(k)
  enddo
  
  if ( idims(1)/=nfaces ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    write(*,*) 'idims(1)=',idims(1),nfaces
    stop
  endif
        
  ncerr = nf90_get_var(ncid, measid, face_area)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  deallocate(idims)
  
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  scale_area = coeff_dx**2
  do nf=1,nfaces
    face_area(nf) = face_area(nf) * scale_area
  enddo
  
  do nf=1,nfaces
    scale_nv = 1.0d0/vert_per_face(nf)
    face_weight(nf) = face_area(nf) * scale_nv
  enddo
  
  do no=1,nnodes
    node_weight(no) = 0.0d0
  enddo
  
  do nf=1,nfaces
    do nv=1,vert_per_face(nf)
      glo = face_to_node(nf,nv)+1 ! C to F numbering
      node_weight(glo) = node_weight(glo) + face_weight(nf) 
    enddo
  enddo
  
end subroutine compute_face_node_weight

subroutine read_face_norm(pfFile,nfaces,face_norm)
  
  use netcdf
  implicit none
  
  integer, intent(in) :: nfaces
  character(len=256),intent(in) :: pfFile 
  real*8, dimension(nfaces,3), intent(out) :: face_norm

  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,nf
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  character(len=256) :: dim_name
  real*8, allocatable, dimension(:) :: buffer
  
  pf_read_debug = .TRUE.
  
  if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'
  
  ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
  ncerr = nf90_inq_varid(ncid, "surfel_normal", measid)
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
  
  if ( idims(1)/=3 .OR. idims(2)/=nfaces ) then
    write(*,'(A)') 'Wrong dimensions! stoping'
    write(*,*) 'idims(1)=',idims(1),3
    write(*,*) 'idims(2)=',idims(2),nfaces
    stop
  endif
  
  allocate(start(rank))
  allocate(count(rank))
  allocate(buffer(nfaces))
  
  do k = 1,3
    start=(/k,1/)
    count=(/1,nfaces/)
    
    ncerr = nf90_get_var(ncid, measid, buffer,start = start, count=count )
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop
    
    do nf = 1,nfaces
      face_norm(nf,k) = buffer(nf)
    enddo
    
  enddo
  
  deallocate(start)
  deallocate(idims)
  deallocate(count)
  deallocate(buffer)
  
  ncerr = nf90_close(ncid)
  ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
  if (ncerr /= NF90_NOERR) stop
  
end subroutine read_face_norm

subroutine read_snc_frame(pfFile,frame_number,nnodes,face_selection,node_selection, &
                face_weight, face_to_node,vert_per_face,read_indices,scale_types, &
                data_sel_node,sel_nfaces,sel_nnodes,nvars,nfaces,max_nv)
  
  use netcdf
  use pf_params
  implicit none
  
  integer, intent(in) :: frame_number
  integer, intent(in) :: nnodes
  integer, intent(in) :: nfaces
!f2py optional , depend(vert_per_face) :: nfaces=len(vert_per_face)
  integer, intent(in) :: max_nv
!f2py optional , depend(face_to_node) :: max_nv=size(face_to_node,dim=2)
  integer, intent(in) :: sel_nfaces
!f2py optional , depend(face_selection) :: sel_nfaces=len(face_selection)
  integer, intent(in) :: sel_nnodes
!f2py optional , depend(node_selection) :: sel_nnodes=len(node_selection)
  integer, intent(in) :: nvars
!f2py optional , depend(read_indices) :: nvars=len(read_indices)
  integer*4, dimension(sel_nfaces), intent(in) :: face_selection
  integer*4, dimension(sel_nnodes), intent(in) :: node_selection
  integer*4, dimension(nfaces,max_nv), intent(in) :: face_to_node
  integer*4, dimension(nfaces), intent(in) :: vert_per_face
  real*8, dimension(nfaces), intent(in) :: face_weight
  integer, dimension(nvars), intent(in) :: read_indices
  integer, dimension(nvars), intent(in) :: scale_types
  character(len=256),intent(in) :: pfFile 
  real*8, dimension(nvars,sel_nnodes), intent(out) :: data_sel_node
  
  ! Local variables

  logical :: pf_read_debug
  
  integer :: ncid,ncerr,measid
  integer :: rank,k,idx,iv,ivar,nvert
  integer*4 :: no, glo_no,nf,glo_face
  integer, dimension(NF90_MAX_VAR_DIMS) :: meas_dim_ids
  integer, dimension(:), allocatable :: idims,start,count
  real*8 :: eps, iweight,min_val,max_val,mean_val,mean2_val,std_val
  character(len=256) :: dim_name
  real*8, allocatable, dimension(:) :: tmp, tmp_s, node_weight
  real*8, allocatable, dimension(:,:) :: buffer
  real*8, allocatable, dimension(:,:) :: data_face, data_node
  
  
  pf_read_debug = .TRUE.
  
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
  
  if (idims(1)/=nfaces) then
    write(*,'(A,1X,I9)') "Expecting Number of faces:",nfaces
    stop
  endif
    
  allocate(start(rank))
  allocate(count(rank))
  
  allocate(buffer(idims(1),idims(2)))
  
  start=(/1,1,frame_number+1/) ! C to F numbering
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
  
  if (pf_read_debug) write(*,*) 'allocate data face'
  
  allocate(tmp(sel_nfaces))
  allocate(tmp_s(sel_nfaces))
  allocate(data_face(sel_nfaces,nvars))
  
  do ivar = 1,nvars
    
    if (pf_read_debug) write(*,*) 'scaling ivar=',ivar
    
    ! C to F index
    idx = read_indices(ivar) + 1 
    
    if (pf_read_debug) write(*,*) 'corresponding idx=',idx
    
    do nf = 1,sel_nfaces
      glo_face = node_selection(nf) + 1  ! C to F numbering
      tmp(nf) = buffer(glo_face,idx)
    enddo
    
    if (pf_read_debug) write(*,*) 'tmp is ready'
    
    call scale_var(scale_types(ivar),tmp,tmp_s,sel_nfaces)
    do nf = 1,sel_nfaces
      data_face(nf,ivar) = tmp_s(nf)
    enddo
    
    if (pf_read_debug) write(*,*) 'face data is scaled'
      
  enddo
  
  if (pf_read_debug) write(*,*) 'scatter to vertices'
  
  allocate(node_weight(nnodes))
  allocate(data_node(nvars,nnodes))
  
  if (pf_read_debug) write(*,*) 'Initialize data_node and node_weight'
  data_node(:,:) = 0.0d0
  node_weight(:) = 0.0d0
  
  eps = 1.0e5
  do nf = 1,sel_nfaces
    glo_face = face_selection(nf) + 1 ! C to F numbering
    nvert = vert_per_face(glo_face)
    do iv = 1, nvert
      no = face_to_node(glo_face,iv) + 1 ! C to F numbering
      do ivar = 1,nvars
        data_node(ivar,no) = data_node(ivar,no) + data_face(nf,ivar)*face_weight(glo_face)
      enddo
      node_weight(no) = node_weight(no) + face_weight(glo_face)
    enddo
  enddo
  
  eps = minval(face_weight)*0.1
  if (pf_read_debug) write(*,*) 'gather to node'
  if (pf_read_debug) write(*,*) 'eps=',eps
  
  do no=1,sel_nnodes
    glo_no = node_selection(no) + 1 ! C to F numbering
    iweight =  1.0d0 / max(node_weight(glo_no),eps)
    do ivar = 1,nvars
      data_sel_node(ivar,no) = data_node(ivar,glo_no) * iweight
    enddo
  enddo
  
  if (pf_read_debug) then
    do ivar = 1,nvars
      min_val =1.0d8
      max_val =-1.0d8
      mean_val = 0.0d0
      mean2_val = 0.0d0
      do no=1,sel_nnodes
        min_val = min(data_sel_node(ivar,no),min_val)
        max_val = max(data_sel_node(ivar,no),max_val)
        mean_val = mean_val + data_sel_node(ivar,no)
        mean2_val = mean2_val + data_sel_node(ivar,no)**2
      enddo
      mean_val = mean_val/nnodes
      mean2_val = mean2_val/nnodes
      std_val = sqrt(mean2_val-mean_val*mean_val)
      write(*,'(A,I1,1X,A,E16.8,1X,A,E16.8,A)') 'ivar=',ivar,'(min=',min_val,' max=',max_val,')'
      write(*,'(A,I1,1X,A,E16.8,1X,A,E16.8,A)') 'ivar=',ivar,'(mean=',mean_val,'std=',std_val,')'
    enddo
  endif
  
  deallocate(data_face)
  deallocate(data_node)
  deallocate(node_weight)
  deallocate(buffer)
  deallocate(tmp)
  deallocate(tmp_s)
  
end subroutine read_snc_frame
