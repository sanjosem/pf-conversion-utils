module pf_params
  
  implicit none
  
  integer :: nvars
  real*8 :: coeff_dx
  real*8 :: coeff_dt
  real*8 :: coeff_vel
  real*8 :: coeff_rho
  real*8 :: coeff_press
  real*8 :: offset_press
  real*8 :: weight_p2r
  real*8 :: weight_r2p
  integer, allocatable, dimension(:) :: var_indices
  character(len=64), allocatable, dimension(:) :: var_names

contains
  
  subroutine scale_var(scale_type,var_lat,var_SI)
    
    integer, intent(in) :: n
    character(len=64), intent(in) :: scale_type
    real*8, dimension(n), intent(in) :: var_lat
    real*8, dimension(n), intent(out) :: var_SI
    
    integer:: i
    
    select case (TRIM(scale_type))
    case('length')
        do i=1,n
          var_SI(n) = var_lat(n) * coeff_dx
        end do 
    case('time')
        do i=1,n
          var_SI(n) = var_lat(n) * coeff_dt
        end do 
      case('static_pressure')
        do i=1,n
          var_SI(n) = (var_lat(n) + offset_press) * coeff_press
        end do 
      case('pressure_from_density')
        do i=1,n
          var_SI(n) = (var_lat(n) * weight_r2p + offset_press) * coeff_press
        end do 
      case('density')
        do i=1,n
          var_SI(n) = var_lat(n) * coeff_rho
        end do 
      case('density_from_pressure')
        do i=1,n
          var_SI(n) = var_lat(n) * weight_p2r * coeff_rho
        end do 
      case('x_velocity','y_velocity','z_velocity')
        do i=1,n
          var_SI(n) = var_lat(n) * coeff_vel
        end do 
    end select
    
  end subroutine scale_var

end module pf_params
