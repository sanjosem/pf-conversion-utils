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
    !f2py optional , depend(var_lat) :: n=len(var_lat)
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

end module pf_params
