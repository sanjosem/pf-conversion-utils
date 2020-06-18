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


subroutine compute_face_node_weight(pfFile,coeff_dx,vert_per_face,face_to_node,nnodes,nfaces,max_nv,face_weight, &
                        node_weight,face_area)
  
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
