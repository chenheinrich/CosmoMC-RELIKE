! =============================================================     
! Header file for the custom function custom_xe(param_array, z) 
! =============================================================
interface			      

      
  function custom_xe(transformed_parameters, z)
    use settings
    use constants
    use ReionizationTanh
    implicit none

    real(mcp) custom_xe 
    real(mcp), dimension(:), intent(in) :: transformed_parameters
    real(mcp), intent(in) :: z
  end function custom_xe

  subroutine get_zre_from_tau(tau, zre)
    ! Getting z_re from optical_depth 
    use settings
    use constants
    use ReionizationTanh
    implicit none

    real(mcp), intent(in):: tau
    real(mcp), intent(out) :: zre
  end subroutine get_zre_from_tau

  subroutine get_transformed_parameters(parameters, transformed_parameters)
    use settings
    implicit none

    real(mcp), dimension(:), intent(in) :: parameters
    real(mcp), dimension(:), intent(out) :: transformed_parameters
  end subroutine get_transformed_parameters

end interface