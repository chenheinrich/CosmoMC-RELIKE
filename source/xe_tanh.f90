module ReionizationTanh
  use Precision
  use MiscUtils
  use MpiUtils, only : MpiStop
  use cosmology
  implicit none  

!This module puts smooth tanh reionization of specified mid-point (z_{re}) and width
!The tanh function is in the variable (1+z)**Rionization_zexp

!Rionization_zexp=1.5 has the property that for the same z_{re} 
!the optical depth agrees with infinitely sharp model for matter domination 
!So tau and zre can be mapped into each other easily (for any symmetric window)
!However for generality the module maps tau into z_{re} using a binary search
!so could be easily modified for other monatonic parameterizations.

!AL March 2008
!AL July 2008 - added trap for setting optical depth without use_optical_depth

!See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

  character(LEN=*), parameter :: Reionization_Name = 'CAMB_reionization'

  real(dl), parameter :: Reionization_DefFraction = -1._dl 
  !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other
  
  real(dl) :: Reionization_AccuracyBoost = 1.3_dl
  !VM ENDS
  real(dl) :: Rionization_zexp = 1.5_dl
  logical  :: include_helium_fullreion = .true.
  real(dl) :: helium_fullreion_redshift  = 3.5_dl
  real(dl) :: helium_fullreion_deltaredshift  = 0.5
  real(dl) :: helium_fullreion_redshiftstart  = 5._dl
  integer, private, parameter :: reion_max_nbasis = 100
	! model x_e(t) = \bar{x}_e + \sum_j m_j S_j(z)
  Type BasisFineGrid
		integer :: nz, nbasis
		real(dl), dimension(:), allocatable :: z
	  real(dl), dimension(:,:), allocatable :: sj
    real(dl) :: zmin, zmax, dz, sigma
	end Type BasisFineGrid	
	Type ReonizationBasis
	  character(LEN=1024) :: file_name
		logical :: init_done = .false.
    logical :: input_file_format_mortoson = .false.
		integer :: nz, nbasis
    real(dl) :: xe_fiducial
		real(dl), dimension(:), allocatable :: z
	  real(dl), dimension(:,:), allocatable :: sj
    real(dl) :: zmin, zmax, dz_min
    Type(BasisFineGrid) :: fine
		contains
		procedure,public :: init=>reonizationbasis_init
    procedure,private :: eval_xe=>reonizationbasis_eval_xe
		procedure,public :: eval_basis=>reonizationbasis_eval_basis
		procedure,private :: eval_fiducial=>reonizationbasis_eval_fiducial
    procedure,private :: setup_basis=>reonizationbasis_setup_basis
	end Type ReonizationBasis

  type ReionizationParams
    logical :: Reionization
    logical :: use_optical_depth
    real(dl) :: redshift, delta_redshift, fraction
    real(dl) :: optical_depth
    integer :: nbasis
    real(dl) :: mj(reion_max_nbasis)
    real(dl) :: AccBoost
    real(dl) :: tau
    logical  :: done_tau = .false.

  end type ReionizationParams

  type ReionizationHistory
    !These two are used by main code to bound region where xe changing 
    real(dl) :: tau_start, tau_complete
    !This is set from main code          
    real(dl) :: akthom, fHe
    !The rest are internal to this module.
    real(dl) :: WindowVarMid, WindowVarDelta
  end type ReionizationHistory
  
  real(dl) :: Reionization_maxz = 50._dl 
  
  real(dl), private, parameter :: Reionization_tol = 1d-5

  Type(ReionizationParams), private, pointer ::  ThisReion
  Type(ReionizationHistory), private, pointer :: ThisReionHist

  Type(ReonizationBasis), save :: xe_basis
  logical :: use_basis	

contains
  ! ----------------------------------------------------------------

  !Background evolution copied from CosmoMC for fiducial cosmology
  ! (Planck 2015 best-fit tanh dz = 0.1 model)
function dtauda_fixed_cosmology(a)
    !get d tau / d a
    use precision
    !use ModelParams
    use MassiveNu
    !use LambdaGeneral
    implicit none
    real(dl) dtauda_fixed_cosmology
    real(dl), intent(IN) :: a
    real(dl) rhonu,grhoa2, a2, pnu
    integer nu_i
    
    ! Got these by printing from equations.f90 
    real(dl) :: grhok = 0.000000000000000E+000
    real(dl) :: grhoc = 3.995382819431314E-008
    real(dl) :: grhob = 7.424529124136492E-009
    real(dl) :: grhog = 8.260148100755502E-012
    real(dl) :: w_lam = -1.00000000000000     
    real(dl) :: grhornomass = 3.809408986356058E-012
    real(dl) :: grhov = 1.035798870387910E-007

    real(dl) :: Num_Nu_massive = 1
    real(dl) :: nu_mass_eigenstates = 1

    real(dl) :: nu_masses = 353.712881263422     
    real(dl) :: grhormass = 1.903307442173650E-012
    
    a2=a**2
    !  8*pi*G*rho*a**4.
    grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
    !if (is_cosmological_constant) then
    grhoa2=grhoa2+grhov*a2**2
    !else
    !    grhoa2=grhoa2+ grho_de(a)
    !end if

	if (Num_Nu_massive /= 0) then
		!Get massive neutrino density relative to massless
        !call Nu_rho(a*nu_masses,rhonu)
        call ThermalNuBackground%rho_P(a*nu_masses, rhonu, pnu)
        grhoa2=grhoa2+rhonu*grhormass
	end if
    dtauda_fixed_cosmology=sqrt(3/grhoa2)

end function dtauda_fixed_cosmology

  subroutine reonizationbasis_init(this,file_name,xef,smooth_sigma) 
  	character(LEN=1024), intent(in) :: file_name
    real(dl), intent(in) :: xef,smooth_sigma
  	class(ReonizationBasis) :: this
  	if (this%init_done .eqv. .false.) then
      this%file_name = file_name
      call this%setup_basis(smooth_sigma)
      this%xe_fiducial = xef
    	this%init_done = .true.
    end if  
  end subroutine reonizationbasis_init
  
	subroutine reonizationbasis_setup_basis(this,smooth_sigma) 
		integer :: status, i, j
    real(dl), intent(in) :: smooth_sigma
		class(ReonizationBasis) :: this			
		open(unit=10,file=this%file_name,status='old', form ='formatted',&
     iostat=status,action='read')		
		if(status .ne. 0) then
			write (*,*) 'Error (reionization): open basis file failed'
			stop
		endif
    ! we assume that the first two entries are nz and nbasis 
		read(10,*,iostat=status) this%nz,this%nbasis
		if(status .ne. 0) then
			write (*,*) 'Error (reionization): read basis file failed'
			stop
		endif
		! We need enough points for interpolation
		if((this%nz .gt. 5) .neqv. .true.) then
			write (*,*) 'Error (reionization): insufficient number redshift bins'
			stop
		endif
		! check if nbasis > 0
    if((this%nbasis .gt. 0) .neqv. .true.) then
			write (*,*) 'Error (reionization): insufficient number of basis'
			stop
		endif
    ! memory allocation of the basis and fiducial xe
    allocate(this%z(this%nz),this%sj(this%nz,this%nbasis),STAT=status)
		if(status .ne. 0) then
			write (*,*) 'Error (reionization): memory allocation failed'
			stop
		endif
		if(this%input_file_format_mortoson .eqv. .false.) then
      ! read file
			do i=1,this%nz,1 		  
				! we assume first column gives z then 
        ! the subsequent columns provide the value of the basis components at z
				read(10,*,iostat=status) this%z(i),(this%sj(i,j),j=1,this%nbasis,1)
        if(status .ne. 0) then
				  write (*,*) 'Error (reionization): read basis file failed'
				  stop
			 	endif
			end do
    else 
      read(10,*,iostat=status) (this%z(i),i=1,this%nz,1)
      if(status .ne. 0) then
				write (*,*) 'Error (reionization): read basis file failed'
				stop
		 	endif
      ! THEN LOOP OVER ALL BASIS
      do i=1,this%nbasis,1 		
				read(10,*,iostat=status) (this%sj(j,i),j=1,this%nz,1)
        if(status .ne. 0) then
					write (*,*) 'Error (reionization): read basis file failed'
					stop
			 	endif
      end do  
    end if
		close(unit=10)
		! check if z vector is given in asceding order 
		do i=1,this%nz-1,1 		  
			if( (this%z(i)<this%z(i+1)) .neqv. .true.) then
				write (*,*) 'Error (reionization): z must be given in ascending order'
				stop
		 	endif
		end do
    ! finally get some aux variables
  	this%zmin = 2.0*this%z(1)-this%z(2)
    this%zmax = 2.0*this%z(this%nz)-this%z(this%nz-1)	
    this%dz_min = this%z(2)-this%z(1)
		do i=2,this%nz,1 		  
			if(this%dz_min .gt. this%z(i)-this%z(i-1)) then
			  this%dz_min = this%z(i)-this%z(i-1)
		 	endif
		end do		
		! ----- SETUP FINE GRID -------
    this%fine%zmin = 0.0_dl
    this%fine%zmax = Reionization_maxz-1
		this%fine%dz   = this%dz_min/7.0_dl 
		this%fine%nz   = int((this%fine%zmax-this%fine%zmin)/this%fine%dz)+1
		this%fine%zmax = this%fine%zmin+(this%fine%nz-1)*this%fine%dz
		allocate(this%fine%z(this%fine%nz), & 
			this%fine%sj(this%fine%nz,this%nbasis),STAT=status)
		if(status .ne. 0) then
			write (*,*) 'Error (reionization): memory allocation failed'
			stop
		endif		
		!$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i=1,this%fine%nz,1
			this%fine%z(i)=this%fine%zmin+(i-1)*this%fine%dz
			call this%eval_basis(this%fine%z(i),this%fine%sj(i,1:this%nbasis))
    end do	
		!$OMP END PARALLEL DO 
		! SETUP SIGMA (LN(Z))
		this%fine%sigma = smooth_sigma
	end subroutine reonizationbasis_setup_basis 	   
  ! ----------------------------------------------------------------
  subroutine reonizationbasis_eval_basis(this,z,eval)
		real(dl), intent(in) :: z
		real(dl), dimension(:), intent(inout) :: eval
		integer  :: i, status, lo, hi, mid
  	real(dl) :: dz
    class(ReonizationBasis) :: this
    ! Following Morthoson - we do linear interpolation here
		! Following Morthoson - WE ADD AN EXTRA BIN z(1)-dz,z(1),...z(n),z(n)+dz
    ! for transition between first/last value to zero
    if(z .lt. this%zmin) then        
			do i=1,this%nbasis,1
				eval(i) = 0.0_dl
			end do 
    else if(z .gt. this%zmax) then  
			do i=1,this%nbasis,1
				eval(i) = 0.0_dl
			end do
    else 
	  	if(z .lt. this%z(1)) then        
  			dz = this%z(2)-this%z(1)
        do i=1,this%nbasis,1
  				eval(i)=(this%sj(1,i)-0.0_dl)/dz
          eval(i)=eval(i)*(z-this%zmin)+0.0_dl
  			end do 
	  	else if(z .gt. this%z(this%nz)) then        
  			dz = this%z(this%nz)-this%z(this%nz-1)
        do i=1,this%nbasis,1
  				eval(i)=(0.0_dl-this%sj(this%nz,i))/dz
          eval(i)=eval(i)*(z-this%z(this%nz))+this%sj(this%nz,i)
  			end do 
      else 
        lo = 1
        hi = this%nz

        do while(hi-lo .gt. 1)
          mid = INT(0.5*(hi+lo))
          if(this%z(mid) .gt. z) then
            hi = mid 
          else 
            lo = mid
          endif
        end do
        dz = this%z(hi)-this%z(lo)

				do i=1,this%nbasis,1
  			  eval(i) = (this%sj(hi,i)-this%sj(lo,i))/dz
          eval(i) = eval(i)*(z-this%z(lo))+this%sj(lo,i)
  			end do
      end if
    endif

  end subroutine reonizationbasis_eval_basis
  ! ----------------------------------------------------------------
	subroutine reonizationbasis_eval_fiducial(this,z,eval,xlowz,z_recom,xe_recom)
    real(dl), dimension(:), intent(in) :: z_recom,xe_recom
    real(dl), intent(in) :: xlowz, z
		real(dl), intent(out) :: eval
		integer :: i, status, lo, hi, mid
  	real(dl) :: dz,xod,th, recom_tmp
    class(ReonizationBasis) :: this
    ! IDEA WE WILL CREATE AN EXTRA BIN
    ! z(1)-dz, z(1), ... z(n), z(n)+dz
    ! for transition between first/last value to zero
    ! Following Mortoson - we do linear interpolation here
    if(z .lt. this%zmin) then        
      eval = xlowz
    else if(z .gt. this%zmax) then  
      lo = 1
      hi = size(xe_recom)
      do while(hi-lo .gt. 1)
        mid = INT(0.5*(hi+lo))
        if(z_recom(mid) .gt. z) then !CH: .lt. -> .gt. for increasing z_recom
          hi = mid 
        else 
          lo = mid
        endif
      end do
			recom_tmp=(xe_recom(hi)-xe_recom(lo))/(z_recom(hi)-z_recom(lo))
			recom_tmp=recom_tmp*(z-z_recom(lo))+xe_recom(lo)
			eval=recom_tmp
    else 
	  	if(z .lt. this%z(1)) then
  			dz = this%z(2)-this%z(1)
        eval=(this%xe_fiducial-xlowz)/dz
        eval=eval*(z-(this%z(1)-dz))+xlowz
	  	else if(z .gt. this%z(this%nz)) then          
	      lo = 1
	      hi = size(xe_recom)
	      do while(hi-lo .gt. 1)
	        mid = INT(0.5*(hi+lo))
	        if(z_recom(mid) .gt. this%zmax) then !CH: changed z_recom decreasing -> increasing
	          hi = mid 
	        else 
	          lo = mid
	        endif
	      end do
				recom_tmp=(xe_recom(hi)-xe_recom(lo))/(z_recom(hi)-z_recom(lo))
				recom_tmp=recom_tmp*(this%zmax-z_recom(lo))+xe_recom(lo)
				dz = this%z(this%nz)-this%z(this%nz-1)
        eval=(recom_tmp-this%xe_fiducial)/dz
        eval=eval*(z-this%z(this%nz))+this%xe_fiducial
      else 
        eval=this%xe_fiducial
      end if
    endif
		! INCLUDE HELIUM EFFECTS
    if(include_helium_fullreion .and. z<helium_fullreion_redshiftstart) then
    	!Effect of Helium becoming fully ionized is small so details not important
    	xod=(helium_fullreion_redshift-z)/helium_fullreion_deltaredshift
      if (xod > 100) then
      	th=1.d0
      else
      	th=tanh(xod)
      end if
			eval=eval+ThisReionHist%fHe*(th+1._dl)/2._dl
		end if
  end subroutine reonizationbasis_eval_fiducial
  ! ----------------------------------------------------------------
	subroutine reonizationbasis_eval_xe(this,z,xlowz,z_recom,xe_recom,mj,eval)
    class(ReonizationBasis) :: this
    real(dl), dimension(:), intent(in) :: mj,z_recom,xe_recom
    real(dl), intent(in) :: xlowz, z
    real(dl), intent(out) :: eval
    real(dl) :: xe_fiducial,xe_tmp,sigma2,gauss,norm
    integer :: i 
    ! gaussian smoothing in the fine grid
		sigma2 = 2*(this%fine%sigma)**2
    eval = 0.0
    norm = 0.0     
    ! TRAPEZOIDAL RULE
		call this%eval_fiducial(this%fine%z(1),xe_fiducial,xlowz, &
			z_recom,xe_recom)
    xe_tmp = xe_fiducial + &
			dot_product(mj(1:this%nbasis),this%fine%sj(1,1:this%nbasis))
		gauss = exp(-(log((1.0_dl+this%fine%z(1))/(1.0_dl+z)))**2/sigma2)/(1.0_dl&
    	+this%fine%z(1))
		eval = eval + xe_tmp*gauss
    norm = norm + gauss
    
		do i=2,this%fine%nz-1,1
			call this%eval_fiducial(this%fine%z(i),xe_fiducial,xlowz, & 
				z_recom,xe_recom)
	    xe_tmp=xe_fiducial + &
				dot_product(mj(1:this%nbasis),this%fine%sj(i,1:this%nbasis))
			gauss = exp(-(log((1.0_dl+this%fine%z(i))/(1.0_dl+z)))**2/sigma2)/(1.0_dl&
				+this%fine%z(i)) 
      eval = eval + 2.0_dl*xe_tmp*gauss
      norm = norm + 2.0_dl*gauss
    end do
		
		call this%eval_fiducial(this%fine%z(this%fine%nz),xe_fiducial, & 
			xlowz,z_recom,xe_recom)
    xe_tmp = xe_fiducial + &
			dot_product(mj(1:this%nbasis),this%fine%sj(this%fine%nz,1:this%nbasis))
		gauss = exp(-(log((1+this%fine%z(this%fine%nz))/(1+z)))**2/sigma2)/(1.0_dl& 
    	+this%fine%z(this%fine%nz))
		eval = eval + xe_tmp*gauss
    norm = norm + gauss    
		! FINAL RESULT
    eval = 0.5_dl*this%fine%dz*eval 
    norm = 0.5_dl*this%fine%dz*norm 
		eval = eval/norm
  end subroutine reonizationbasis_eval_xe     
  ! ----------------------------------------------------------------
  function Reionization_xe2(a,z_recom,xe_recom)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    real(dl), intent(in) :: a
    real(dl), dimension(:), intent(in) :: z_recom,xe_recom
    real(dl) :: Reionization_xe2
    real(dl) :: xstart
  	real(dl) :: tmp
  	real(dl) :: z
	  
    z = 1._dl/a-1._dl
    call xe_basis%eval_xe(z,ThisReion%fraction,z_recom,xe_recom, &
        ThisReion%mj, Reionization_xe2)

  end function Reionization_xe2
  !VM ENDS
  ! ----------------------------------------------------------------
  function Reionization_xe(a, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    real(dl), intent(in) :: a
    real(dl), intent(in), optional :: tau, xe_recomb
    !CH: tau doesn't seem to be used...?
    !CH BEGINS
    logical :: write_xez = .false.
    !CH ENDS
    real(dl) Reionization_xe
    real(dl) tgh, xod
    real(dl) xstart
    if (present(xe_recomb)) then
    xstart = xe_recomb
    else
    xstart = 0._dl
    end if  
		!VM BEGINS
		if(use_basis .eqv. .false.) then
		!VM ENDS
      xod = (ThisReionHist%WindowVarMid - 1._dl/a**Rionization_zexp)/ThisReionHist%WindowVarDelta
      if (xod > 100) then
        tgh=1.d0
      else
        tgh=tanh(xod)
      end if
      Reionization_xe =(ThisReion%fraction-xstart)*(tgh+1._dl)/2._dl+xstart
      !CH BEGINS: output xe(z) without helium if not using basis (for mj projection)
      ! file unit 54 specified in modules.f90 in subroutine inithermo
      if (.false.) then 
        if (write_xez .eqv. .true.) then
          if(use_basis .eqv. .false.) then
          write(54,*) (1._dl/a)-1._dl, Reionization_xe
          end if
        end if 
      end if
      !CH ENDS
      if (include_helium_fullreion .and. a>(1/(1+ helium_fullreion_redshiftstart))) then
        !Effect of Helium becoming fully ionized at z <~ 3.5 is very small so details not important
        xod = (1+helium_fullreion_redshift - 1._dl/a)/helium_fullreion_deltaredshift
        if (xod > 100) then
           tgh=1.d0
        else
           tgh=tanh(xod)
        end if
        Reionization_xe =  Reionization_xe + ThisReionHist%fHe*(tgh+1._dl)/2._dl
      end if
		!VM BEGINS
		else
			write(*,*) 'Error (reionization): illegal call with use_basis=true'
			stop
		end if
		!VM ENDS
  end function Reionization_xe  
  ! ----------------------------------------------------------------
  function Reionization_timesteps(ReionHist)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    Type(ReionizationHistory) :: ReionHist
    integer Reionization_timesteps
    Reionization_timesteps = 50 
  end  function Reionization_timesteps 
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  subroutine Reionization_SetParamsForZre(Reion,ReionHist)
    !use AMLutils
    use MiscUtils
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist 
	  
    Reion%delta_redshift = 0.015 * (1._dl + Reion%redshift)
    ReionHist%WindowVarMid = (1._dl+Reion%redshift)**Rionization_zexp
    ReionHist%WindowVarDelta = Rionization_zexp*&
      (1._dl+Reion%redshift)**(Rionization_zexp-1._dl)*Reion%delta_redshift
  end subroutine Reionization_SetParamsForZre
  ! ----------------------------------------------------------------
  subroutine Reionization_Init(Reion,ReionHist,Yhe,akthom,tau0,FeedbackLevel)
    use constants
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist
    real(dl), intent(in) :: akthom, tau0, Yhe 
    integer, intent(in) :: FeedbackLevel
    !VM BEGINS
    real(dl) :: astart, zstart, zend, aend, epsrel
		!VM ENDS
    ReionHist%akthom = akthom  
    ReionHist%fHe =  YHe/(mass_ratio_He_H*(1.d0-YHe))
    ReionHist%tau_start=tau0
    ReionHist%tau_complete=tau0
    ThisReion => Reion
    ThisReionHist => ReionHist

 	if (Reion%Reionization .and. (use_basis .neqv. .true.)) then
 	  
      if (Reion%optical_depth /= 0._dl .and. .not. Reion%use_optical_depth) &
        write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'
      
      if (Reion%use_optical_depth.and.Reion%optical_depth<0.001 &
        .or. .not.Reion%use_optical_depth .and. Reion%Redshift<0.001) then
             Reion%Reionization = .false.
      end if  

    end if    

    if (Reion%Reionization .and. (use_basis .neqv. .true.)) then
		
      if (Reion%fraction==Reionization_DefFraction) &
         Reion%fraction = 1._dl + ReionHist%fHe  !H + singly ionized He

      if (Reion%use_optical_depth) then
        call Reionization_SetFromOptDepth(Reion,ReionHist)
        if (FeedbackLevel > 0)  then
          write(*,'("Reion redshift       =  ",f6.3)') Reion%redshift
        end if
      end if
      
      call Reionization_SetParamsForZre(ThisReion,ThisReionHist) 
      
      !this is a check, agrees very well in default parameterization
      if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') &
            Reionization_GetOptDepth(Reion, ReionHist)          
      !Get relevant times       
      astart=1.d0/(1.d0+Reion%redshift + Reion%delta_redshift*8)

      ReionHist%tau_start = max(0.05_dl, rombint(dtauda_fixed_cosmology,0._dl,astart,1d-3))
          !Time when a very small reionization fraction (assuming tanh fitting)
      ReionHist%tau_complete = min(tau0,ReionHist%tau_start+ rombint(dtauda_fixed_cosmology,& 
       astart,1.d0/(1.d0+max(0.d0,Reion%redshift-Reion%delta_redshift*8)),1d-3))

    else if (Reion%Reionization .and. use_basis) then
      if(Reion%fraction==Reionization_DefFraction) then
        Reion%fraction = 1._dl + ReionHist%fHe  !H + singly ionized He
      end if
      if(log((1._dl+Reionization_maxz)/(1._dl+xe_basis%zmax)) & 
        <16.0_dl*xe_basis%fine%sigma) then
       Reionization_maxz=(1+xe_basis%zmax)*exp(16.0*xe_basis%fine%sigma)-1
      end if     
      ! eta begins -------
      epsrel = 1e-4_dl
      ! eta START -------
      !CH BEGINS (changed for zmax = 50)
      !zstart = (1+xe_basis%zmax)*exp(8.0*xe_basis%fine%sigma)-1
      zstart = 55.0_dl
      !CH ENDS (changed for zmax = 50)
      astart = 1.d0/(1.d0+zstart)  
      ReionHist%tau_start =  max(0.05_dl,rombint(dtauda_fixed_cosmology,0.0_dl,astart,epsrel))
      ! eta COMPLETE --- not including helium - why even std camb does that??   
      !CH temp begins
      !zend = (1+xe_basis%zmin)*exp(-8.0*xe_basis%fine%sigma)-1  
      zend = 5.5_dl
      !CH temp ends 
      zend = max(zend, 0.5)
      aend = 1.d0/(1.d0+zend)
      ReionHist%tau_complete = min(tau0,rombint(dtauda_fixed_cosmology,0.0_dl,aend,epsrel))
    
      ! eta ends -------		
    end if   
    !VM ENDS   
  end subroutine Reionization_Init
  
  subroutine Reionization_SetDefParams(Reion)
    Type(ReionizationParams) :: Reion
    Reion%Reionization = .true.
    Reion%use_optical_depth = .false.
    Reion%optical_depth = 0._dl
    Reion%redshift = 10
    Reion%fraction = Reionization_DefFraction
    Reion%delta_redshift = 0.5_dl
		!VM BEGINS
		use_basis = .false.
		!VM ENDS
  end subroutine Reionization_SetDefParams

  subroutine Reionization_SetDefParams_KDE(Reion, include_helium_fullreion_in)
    Type(ReionizationParams) :: Reion
    logical, intent(in) :: include_helium_fullreion_in
    Reion%Reionization = .true.
    Reion%use_optical_depth = .true.
    Reion%optical_depth = 0._dl
    Reion%fraction = Reionization_DefFraction 
    Reion%delta_redshift = 0.5_dl
    use_basis = .false.
    include_helium_fullreion = include_helium_fullreion_in
  end subroutine Reionization_SetDefParams_KDE

  subroutine Reionization_Validate(Reion, OK)
    Type(ReionizationParams),intent(in) :: Reion
    logical, intent(inout) :: OK
    !VM BEGINS
    !if (Reion%Reionization) then
    if (Reion%Reionization .and. (use_basis .neqv. .true.)) then
    !VM ENDS
      if (Reion%use_optical_depth) then
        if (Reion%optical_depth<0 .or. Reion%optical_depth > 0.9  .or. &
           include_helium_fullreion .and. Reion%optical_depth<0.01) then
         OK = .false.
         write(*,*) 'Optical depth is strange. You have:', Reion%optical_depth 
        end if
      else
        if(Reion%redshift < 0 .or. &
          Reion%Redshift+Reion%delta_redshift*3 > Reionization_maxz .or. &
          include_helium_fullreion .and. &
          Reion%redshift < helium_fullreion_redshift) then
          OK = .false.
          write(*,*) 'Reionization redshift strange. You have: ',Reion%Redshift
        end if
      end if
      if(Reion%fraction/= Reionization_DefFraction .and. &
        (Reion%fraction < 0 .or. Reion%fraction > 1.5)) then
        OK = .false.
        write(*,*) 'Reionization fraction strange. You have: ',Reion%fraction
      end if
      if (Reion%delta_redshift > 3 .or. Reion%delta_redshift<0.1 ) then
        !Very narrow windows likely to cause problems in interpolation etc.
        !Very broad likely to conflic with quasar data at z=6
        OK = .false.
        write(*,*) &
        'Reionization delta_redshift is strange. You have:',Reion%delta_redshift
      end if 
    end if
  end  subroutine Reionization_Validate
  
  function Reionization_doptdepth_dz(z)
   real(dl) :: Reionization_doptdepth_dz
   real(dl), intent(in) :: z
   real(dl) a
   a = 1._dl/(1._dl+z)
   
   Reionization_doptdepth_dz = &
       Reionization_xe(a)*ThisReionHist%akthom*dtauda_fixed_cosmology(a)
   !VM ENDS
  end function Reionization_doptdepth_dz


 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rombint2(f,a,b,tol, maxit, minsteps)
        use precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.

! Modified by AL to specify max iterations and minimum number of steps
! (min steps useful to stop wrong results on periodic or sharp functions)
        implicit none
        integer, parameter :: MAXITER=20,MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint2
        real(dl), intent(in) :: a,b,tol
        integer, intent(in):: maxit,minsteps
     
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
      
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
        do
          i=i+1
          if (i > maxit.or.(i > 5.and.abs(error) < tol) .and. nint > minsteps) exit
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
          do k=1,nint
            g0=g0+f(a+(k+k-1)*h)
          end do
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
          do j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
          end do  
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        end do
  
        rombint2=g0
        if (i > maxit .and. abs(error) > tol)  then
          write(*,*) 'Warning: Rombint2 failed to converge; '
          write (*,*)'integral, error, tol:', rombint2,error, tol
        end if
        
    end function rombint2
  
  function Reionization_GetOptDepth(Reion, ReionHist) 
    Type(ReionizationParams), target :: Reion
    Type(ReionizationHistory), target :: ReionHist
    real(dl) Reionization_GetOptDepth     

    ThisReion => Reion
    ThisReionHist => ReionHist 

  	if(use_basis .eqv. .false.) then	
      Reionization_GetOptDepth =  &
        rombint2(Reionization_doptdepth_dz,0.d0,Reionization_maxz,&
        Reionization_tol, 20, nint(Reionization_maxz/Reion%delta_redshift*5))
    end if

  end function Reionization_GetOptDepth
  
  subroutine Reionization_zreFromOptDepth(Reion, ReionHist)
    !General routine to find zre parameter given optical depth
    !Not used for Rionization_zexp = 1.5
    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist
    real(dl) try_b, try_t
    real(dl) tau
    integer i

  	if(use_basis .eqv. .false.) then	

        try_b = 0
        try_t = Reionization_maxz
        i=0
        do 
        i=i+1  
        Reion%redshift = (try_t + try_b)/2
        call Reionization_SetParamsForZre(Reion,ReionHist)
        tau = Reionization_GetOptDepth(Reion, ReionHist)
        if (tau > Reion%optical_depth) then
            try_t = Reion%redshift
        else
            try_b = Reion%redshift
        end if
        if (abs(try_b - try_t) < 2e-3/Reionization_AccuracyBoost) exit
        
        
        !HACK need to put in proper stop command
        if (i>100) then
            call mpiStop('Reionization_zreFromOptDepth: failed to converge')
            !stop
        end if
        end do
        if (abs(tau - Reion%optical_depth) > 0.002) then
            write(*,*) 'Reionization_zreFromOptDepth: Did not converge to optical depth'
            write(*,*) 'tau =',tau, 'optical_depth = ', Reion%optical_depth
            write(*,*) try_t, try_b
            !stop
        !HACK need to put in proper stop command
            call mpiStop('Reionization_zreFromOptDepth: Did not converge to optical depth')
        end if

    end if
  end subroutine Reionization_zreFromOptDepth 

  subroutine Reionization_SetFromOptDepth(Reion, ReionHist)
    Type(ReionizationParams) :: Reion
    Type(ReionizationHistory) :: ReionHist
    ! This subroutine calculates the redshift of reionization
    ! This implementation is approximate but quite accurate and fast d
    real(dl) dz, optd
    real(dl) z, tmp, tmpHe
    integer na    

    Reion%redshift = 0
    if (Reion%Reionization .and. Reion%optical_depth /= 0) then
      !Do binary search to find zre from z
      !This is general method
      call Reionization_zreFromOptDepth(Reion, ReionHist)
      if (.false.) then
        !Use equivalence with sharp for special case
        optd=0
        na=1
        dz=1._dl/2000/Reionization_AccuracyBoost
        tmp = dz*Reion%fraction*ThisReionHist%akthom
        tmpHe = dz*(Reion%fraction+ReionHist%fHe)*ThisReionHist%akthom
        z=0
        do while (optd < Reion%optical_depth)
          z=na*dz
          if(include_helium_fullreion .and. z < helium_fullreion_redshift) then
          optd=optd+ tmpHe*dtauda_fixed_cosmology(1._dl/(1._dl+z))
          else
          optd=optd+tmp*dtauda_fixed_cosmology(1._dl/(1._dl+z))
          end if
          na=na+1
        end do
      end if
    else
      Reion%Reionization = .false.
    end if
  end  subroutine Reionization_SetFromOptDepth 

end module ReionizationTanh

