! ===============================================================     
! Custom function custom_xe(z) must follow the following format
! and returns ionization fraction at z without helium ionization
! ===============================================================
!function custom_xe(z)
!    use settings
!    implicit none
!    real(mcp) custom_xe
!    real(mcp), intent(in) :: z
!    custom_xe = 2.0_mcp
!    return
!end function custom_xe

function custom_xe(transformed_parameters, z) result (xe)
    use settings
    use KdeReionizationTanh
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
    use KdeReionizationTanh
    implicit none

    real(mcp), intent(in):: var_value
    logical, intent(in) :: use_optical_depth_in, include_helium_fullreion_in  

    Type(ReionizationParams), intent(out) :: Reion
    Type(ReionizationHistory), intent(out) :: ReionHist 

    ! from tanh ML from fiducial cosmology, needed by Reionization_Init 
    real(mcp) :: YHe = 0.2453368, akthom = 3.867559814364194E-007, tau0 = 14172.4718396201 
    integer :: FeedbackLevel = 0

    call Reionization_SetDefParams_KDE(Reion, include_helium_fullreion_in) 
    Reion%use_optical_depth = use_optical_depth_in

    if (Reion%use_optical_depth == .true.) then
        Reion%optical_depth = var_value
    else 
        Reion%redshift = var_value
    end if

    call init_massive_nu(.true.) ! needed for getting dt/da
    call Reionization_Init(Reion, ReionHist, YHe, akthom, tau0, FeedbackLevel)

end subroutine reion_init


subroutine get_zre_from_tau(tau, zre)

    use settings
    use KdeReionizationTanh
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

function custom_xe2(transformed_parameters, z) result (custom_xe)

    use settings
    use constants
    use Reionization 
    implicit none

    real(mcp) custom_xe 
    real(mcp), dimension(:), intent(in) :: transformed_parameters
    real(mcp), intent(in) :: z

    real(mcp) :: Reion_redshift, YHe
    real(mcp) :: delta_redshift = 0.5_mcp, Reionization_zexp = 1.5_mcp
    real(mcp) :: a, xe, tgh, xod, xstart, fraction, WindowVarMid, WindowVarDelta, fHe

    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist 

    real(mcp) :: optical_depth
    real(mcp) :: xe2
    logical :: include_helium_fullreion_custom_xe = .false. ! want false for getting mjs (PC components)
    
    Reion_redshift = transformed_parameters(1) 

    ! TODO not shared, dangerous (need to make uniform)
    YHe = 0.2453368
    fHe = YHe/(mass_ratio_He_H*(1.0_mcp-YHe))
    fraction = 1._mcp + fHe
    xstart = 0._mcp ! TODO could use recomb xe if needed

    a = 1._mcp/(1._mcp + z)
	WindowVarMid = (1._mcp + Reion_redshift) ** Reionization_zexp
    WindowVarDelta = Reionization_zexp * &
      (1._mcp + Reion_redshift) ** (Reionization_zexp - 1._mcp) *delta_redshift

    xod = (WindowVarMid - 1._mcp / a ** Reionization_zexp) / WindowVarDelta
    if (xod > 100) then
        tgh=1.0_mcp
    else
        tgh=tanh(xod)
    end if

    custom_xe = (fraction - xstart) * (tgh + 1._mcp) / 2._mcp + xstart

    !TODO shared with and defined in reionization.f90 for now: 
    ! helium_fullreion_redshift, helium_fullreion_deltaredshift
    ! include_helium_fullreion
    ! cannot just set include_helium_fullreion = .false. here
    ! because that would turn off helium for zre calculation (we want true there)
    if(include_helium_fullreion_custom_xe .and. z<helium_fullreion_redshiftstart) then
        xod = (1+helium_fullreion_redshift - 1._mcp/a)/helium_fullreion_deltaredshift
            if (xod > 100) then
            tgh=1.d0
            else
            tgh=tanh(xod)
            end if
        custom_xe =  custom_xe + fHe*(tgh+1._mcp)/2._mcp
    end if

end function custom_xe2