! ===============================================================     
! Custom function custom_xe(z) must follow the following format
! and returns ionization fraction at z for 6 < z < 30. 
! For z < 6, we assume fully ionized hydrogen and singly ionized
! helium, and a tanh transition at 3.5 for second helium ioniza-
! tion.
! ===============================================================
!function custom_xe(z)
!    use settings
!    implicit none
!    real(mcp) custom_xe
!    real(mcp), intent(in) :: z
!    custom_xe = <put your custom xe(z) function here>
!    return
!end function custom_xe

function custom_xe(transformed_parameters, z) result (xe)
    use settings
    use ReionizationTanh
    implicit none

    real(mcp), dimension(:), intent(in) :: transformed_parameters
    real(mcp), intent(in) :: z
    real(mcp) :: xe

    real(mcp) :: zre
    logical, parameter :: INCLUDE_HELIUM_FALSE = .false.
    logical, parameter :: USE_OPTICAL_DEPTH_FALSE = .false.

    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist 

    zre = transformed_parameters(1)
    call reion_init(USE_OPTICAL_DEPTH_FALSE, zre, INCLUDE_HELIUM_FALSE, Reion, ReionHist)
    xe = Reionization_xe(1._mcp/(1._mcp + z))

end function custom_xe

subroutine get_transformed_parameters(parameters, transformed_parameters)
    use settings
    implicit none

    real(mcp), dimension(:), intent(in) :: parameters
    real(mcp), dimension(:), intent(out) :: transformed_parameters
    real(mcp) :: tau, zre

    tau = parameters(1)
    call get_zre_from_tau(tau, zre)

    transformed_parameters(1) = zre

end subroutine get_transformed_parameters

subroutine reion_init(use_optical_depth_in, var_value, include_helium_fullreion_in, Reion, ReionHist)
    use settings
    use ReionizationTanh
    use MassiveNu
    implicit none

    real(mcp), intent(in):: var_value
    logical, intent(in) :: use_optical_depth_in, include_helium_fullreion_in  

    Type(ReionizationParams), intent(out) :: Reion
    Type(ReionizationHistory), intent(out) :: ReionHist 

    ! Values from fiducial cosmology (Planck 2015 best-fit tanh)
    real(mcp) :: YHe = 0.2453368, akthom = 3.867559814364194E-007, tau0 = 14172.4718396201 
    integer :: FeedbackLevel = 0

    call Reionization_SetDefParams_KDE(Reion, include_helium_fullreion_in) 
    Reion%use_optical_depth = use_optical_depth_in

    if (Reion%use_optical_depth .eqv. .true.) then
        Reion%optical_depth = var_value
    else 
        Reion%redshift = var_value
    end if

    call ThermalNuBackground%Init() ! needed for getting dtauda in ReionizationTanh module

    call Reionization_Init(Reion, ReionHist, YHe, akthom, tau0, FeedbackLevel)

end subroutine reion_init


subroutine get_zre_from_tau(tau, zre)

    use settings
    use ReionizationTanh
    implicit none

    real(mcp), intent(in):: tau
    real(mcp), intent(out) :: zre

    logical, parameter :: INCLUDE_HELIUM_TRUE = .true.
    logical, parameter :: USE_OPTICAL_DEPTH_TRUE = .true.

    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist 
    
    call reion_init(USE_OPTICAL_DEPTH_TRUE, tau, INCLUDE_HELIUM_TRUE , Reion, ReionHist)
    
    zre = Reion%redshift 

end subroutine get_zre_from_tau	
