module pf_params

  implicit none

  real*8 :: coeff_dx
  real*8 :: timestep
  real*8 :: coeff_vel
  real*8 :: coeff_rho
  real*8 :: coeff_press
  real*8 :: offset_press
  real*8 :: weight_p2r
  real*8 :: weight_r2p

  integer*4, parameter :: type_length = 0
  integer*4, parameter :: type_time = 1
  integer*4, parameter :: type_pressure = 2
  integer*4, parameter :: type_pressure_from_density = 20
  integer*4, parameter :: type_density = 3
  integer*4, parameter :: type_density_from_pressure = 30
  integer*4, parameter :: type_velocity = 4

  logical, parameter :: debug_scaling = .FALSE.
  logical, parameter :: pf_read_debug = .FALSE.

contains

  subroutine scale_var(scale_type,var_lat,var_SI,n)

    integer, intent(in) :: n
    !f2py depend(var_lat) :: n=len(var_lat)
    integer, intent(in) :: scale_type
    real*8, dimension(n), intent(in) :: var_lat
    real*8, dimension(n), intent(out) :: var_SI

    integer:: i

    if (debug_scaling) write(*,*) 'Scaling type:', scale_type

    select case (scale_type)
    case(type_length)
        do i=1,n
          var_SI(i) = var_lat(i) * coeff_dx
        end do
    case(type_time)
        do i=1,n
          var_SI(i) = var_lat(i) * timestep
        end do
      case(type_pressure)
        do i=1,n
          var_SI(i) = (var_lat(i) + offset_press) * coeff_press
        end do
      case(type_pressure_from_density)
        do i=1,n
          var_SI(i) = (var_lat(i) * weight_r2p + offset_press) * coeff_press
        end do
      case(type_density)
        do i=1,n
          var_SI(i) = var_lat(i) * coeff_rho
        end do
      case(type_density_from_pressure)
        if (debug_scaling) write(*,*) weight_p2r
        if (debug_scaling) write(*,*) coeff_rho
        do i=1,n
          var_SI(i) = var_lat(i) * weight_p2r * coeff_rho
        end do
      case(type_velocity)
        do i=1,n
          var_SI(i) = var_lat(i) * coeff_vel
        end do
    end select

  end subroutine scale_var

  subroutine read_time(pfFile,dt,ntime,avg_period,sampling_period,time_center)
    use netcdf

    implicit none

    integer, intent(in) :: ntime
    real*8, intent(in) :: dt
    character(len=256),intent(in) :: pfFile
    real*8, intent(out) :: avg_period,sampling_period
    real*8, dimension(ntime), intent(out) :: time_center

    ! Local variables
    integer :: ncid,ncerr,measid
    integer :: rank,k,nt
    integer, dimension(5) :: meas_dim_ids
    integer, dimension(:), allocatable :: idims
    character(len=256) :: dim_name
    real*8 :: minv,maxv, check
    real*8, dimension(:), pointer :: buffer, work
    real*8, dimension(:), allocatable :: tstart, tend

    ! store in module
    timestep = dt

    if (pf_read_debug) write(*,'(A,1X,A,1X,A)') 'Opening file',TRIM(pfFile),'for reading'

    ncerr = nf90_open(pfFile, NF90_NOWRITE, ncid)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop

    ncerr = nf90_inq_varid(ncid, "start_time", measid)
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

    if ( idims(1)/=ntime ) then
      write(*,'(A)') 'Wrong dimensions! stoping'
      write(*,*) 'idims(1)=',idims(1),ntime
      stop
    endif

    allocate(buffer(ntime))
    allocate(tstart(ntime))
    allocate(tend(ntime))

    ncerr = nf90_get_var(ncid, measid, buffer)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop

    call scale_var(type_time,buffer,tstart,ntime)

    ncerr = nf90_inq_varid(ncid, "end_time", measid)
    if (ncerr /= NF90_NOERR) stop

    ncerr = nf90_get_var(ncid, measid, buffer)
    if (ncerr /= NF90_NOERR) stop

    call scale_var(type_time,buffer,tend,ntime)

    do nt=1,ntime
      time_center(nt) = 0.5*(tend(nt) + tstart(nt))
    enddo

    write(*,'(A,I6,A,E16.6,1X,A)') '  -> First frame (',0,'): ',time_center(1),' s'
    if (ntime>1) then
      write(*,'(A,I6,A,E16.6,1X,A)') '  -> Last frame (',ntime,'): ',time_center(ntime),' s'
    endif

    deallocate(idims)

    ncerr = nf90_close(ncid)
    ! if (ncerr /= NF90_NOERR) call handle_error(ncerr)
    if (ncerr /= NF90_NOERR) stop


    do nt=1,ntime
      buffer(nt) = (tend(nt) - tstart(nt))
    enddo

    minv = minval(buffer)
    maxv = maxval(buffer)
    check = abs(maxv-minv)/dt
    if (check>8) then
      write(*,'(A)') '  ! Averaging period is not constant'
      write(*,'(A,E16.6,1X,E16.6)') '  ! Min/max period [s]: ',minv,maxv
      avg_period = 0.0d0
    else
      avg_period = sum(buffer)/real(ntime,8)
    endif

    if (ntime>1) then

      do nt=1,ntime-1
        buffer(nt) = time_center(nt+1) - time_center(nt)
      enddo

      work => buffer(1:(ntime-1))

      minv = minval(work)
      maxv = maxval(work)
      check = abs(maxv-minv)/dt
      if (check>4) then
        write(*,'(A)') '  ! Sampling period is not constant'
        write(*,'(A,E16.6,1X,E16.6)') '  ! Min/max sampling [s]: ',minv,maxv
        sampling_period = 0.0d0
      else
        sampling_period = sum(work)/real(ntime-1,8)
      endif

    else
      sampling_period = 0.0d0
    endif

    deallocate(buffer)

  end subroutine read_time

end module pf_params
