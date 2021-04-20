module KernelDensityEstimate
  use KdeReionizationTanh
  use settings 
  implicit none

  include "custom_xe_tanh.h"

  Type RelikeKdeChains
    integer :: nrow, ncol
    real(mcp), dimension(:,:), allocatable :: data
    contains
    procedure, public :: setup_chains
  end Type RelikeKdeChains
  
  Type RelikeKdeParams 
    integer :: nbasis = 5
    real(mcp) :: mjs(5)

    real(mcp) :: xef = 0.15
    real(mcp) :: smooth_sigma = 0.0

    ! Could be input by user, otherwise default below adopted
    real(mcp) :: fcov
    character(LEN=1024) :: fn_pc, fn_invcov, fn_invcov_spec
    character(LEN=1024) :: fn_chain, fn_chain_spec
    character(LEN=1024) :: fn_mean
    character(LEN=1024) :: kde_or_gaussian

    !default used if input not specified 
    real(mcp) :: default_fcov = 0.14_mcp
    character(LEN=1024) :: default_fn_pc = 'kde_data/pl18_zmax30/pc.dat'
    character(LEN=1024) :: default_fn_invcov= "kde_data/pl18_zmax30/invcov.dat"
    character(LEN=1024) :: default_fn_invcov_spec= "kde_data/pl18_zmax30/invcov_spec.dat"

    character(LEN=1024) :: default_fn_mean= "kde_data/pl18_zmax30/mean.dat"
    character(LEN=1024) :: default_fn_chain = 'kde_data/pl18_zmax30/chains/chains.txt'  
    character(LEN=1024) :: default_fn_chain_spec = 'kde_data/pl18_zmax30/chains/spec.txt'

    character(LEN=1024) :: default_kde_or_gaussian = 'kde'
    
  end Type RelikeKdeParams

  Type RelikeKde
    logical :: init_done = .false.
    real(mcp), dimension(:,:), allocatable :: invcov
    real(mcp), dimension(:), allocatable :: mean
    Type(ReonizationBasis) :: kde_xe_basis
    Type(RelikeKdeParams) :: kde_params
    contains
    procedure, public :: init=>kde_init
    procedure, private :: get_mjs
    procedure, private :: kde_impl
    procedure, private :: gauss_impl
    procedure, private :: setup_invcov
    procedure, private :: setup_mean
  end Type RelikeKde
    
  Type(RelikeKde), private, pointer :: ThisKde
  Type(RelikeKdeChains), pointer, save :: ThisKdeChains 
    
contains

  subroutine kde_read_params(this, Ini)

    class(RelikeKdeParams) :: this
    Type(TSettingIni) :: Ini

    logical :: AllowBlanck = .true.
    
    this%fcov = Ini%Read_Double('kde_fcov', this%default_fcov)
    write (*,*) '... You have set fcov (used for KDE only): ', this%fcov

    if (this%fcov <= 0.0) then
      write (*,*) 'Error (KernelDensityEstimate.f90): kde_fcov must be specified at input and > 0.'
      stop
    end if

    this%fn_pc = Ini%Read_String_Default('fn_pc', \
      this%default_fn_pc, AllowBlanck)

    this%fn_invcov = Ini%Read_String_Default('fn_invcov', \
      this%default_fn_invcov, AllowBlanck)
    this%fn_invcov_spec = Ini%Read_String_Default('fn_invcov_spec', \
      this%default_fn_invcov_spec, AllowBlanck)
  
    this%fn_chain = Ini%Read_String_Default('fn_chain', \
      this%default_fn_chain, AllowBlanck)
    this%fn_chain_spec = Ini%Read_String_Default('fn_chain_spec', \
      this%default_fn_chain_spec, AllowBlanck)

    this%fn_mean = Ini%Read_String_Default('fn_mean', \
      this%default_fn_mean, AllowBlanck)

    this%kde_or_gaussian = Ini%Read_String_Default('kde_or_gaussian', \
      this%default_kde_or_gaussian, AllowBlanck)

  end subroutine kde_read_params

  subroutine kde_init(this)  !called after kde_read_params
    real(mcp), dimension(:), allocatable :: mjs
    integer :: status
    class(RelikeKde) :: this
    integer :: i,j
    
    if (this%init_done .eqv. .false.) then
      print *, '... Setting up PCs and likelihood'
      ! Setup reionization PC class, load inverse PC chain mean and inverse covariance.
      call this%kde_xe_basis%init(this%kde_params%fn_pc, this%kde_params%xef, this%kde_params%smooth_sigma) 
      call this%setup_invcov() 
      call this%setup_mean() 
      this%init_done = .true.
    end if  

  end subroutine kde_init

  function xe_fiducial(basis, z)
    !Evaluates fiducial xe (ie - xe(z) when all PC components are zero
    !xe_fiducial = 1+fHe for z <= 6, xe_fiducial = 0 for z>=30
    !xe_fiducial = xef = 0.15 for 6.25 <= z <= 29.75
    !interpolated linearly for the first and last bin [6, 6.25] and [29.75, 30].
    use settings
    real(mcp), intent(in) :: z
    class(ReonizationBasis), intent(in) :: basis
    real(mcp) :: xe_fiducial
    real(mcp) :: zmin, zmax, xef, xe_lowz
    integer :: nz
    
    nz = basis%nz
    zmin = basis%zmin 
    zmax = basis%zmax 
    xef = basis%xe_fiducial

    call get_default_xe_lowz(xe_lowz)

    if (z < zmin) then
        xe_fiducial = xe_lowz
    else if (z > zmax) then
        xe_fiducial = 0.0_mcp
    else
        if (z < basis%z(1)) then
            xe_fiducial = (z - zmin) * (xef - xe_lowz) / (basis%z(2)-basis%z(1)) + xe_lowz
        else if (z > basis%z(nz)) then
            xe_fiducial = &
              (z - basis%z(nz)) * (0.0_mcp - xef) / (basis%z(nz) - basis%z(nz-1)) + xef
        else
            xe_fiducial = xef
        end if
    end if

  end function xe_fiducial

  subroutine get_default_xe_lowz(xe_lowz)
    use settings
    use constants
    real(mcp), intent(out) :: xe_lowz
    real(mcp) :: YHe, fHe
    YHe = 0.2453368
    fHe = YHe/(mass_ratio_He_H*(1.0_mcp-YHe))
    xe_lowz = 1._mcp + fHe
  end subroutine get_default_xe_lowz

  subroutine get_mjs(this, transformed_parameters, mjs) 
    ! Returns mjs as the array of projected principal component amplitude
    ! for the xe(z) model described by f; projection done in redshift range
    ! [6, 30], via m_j = int_{zmin}^{zmax} (xe(z)-x_fiducial(z)) * b_j(z) where
    ! b_j(z) is the jth principal component function.
    use settings
    real(mcp), dimension(:), intent(out), allocatable :: mjs
    real(mcp), dimension(:), intent(in) :: transformed_parameters
    
    include "custom_xe_tanh.h"

    real(mcp), dimension(:,:), allocatable :: integrand
    real(mcp), dimension(:), allocatable :: xe, xe_fid_array

    real(mcp), dimension(:), allocatable :: tmp
    real(mcp) :: z, zmin, zmax
    real(mcp) :: integral, step_size
    integer :: status, i, j, n_simpson = 1000
    integer :: nmjs, nz
    class(RelikeKde) :: this

    nmjs = this%kde_xe_basis%nbasis

    ! If n is odd we add +1 to make it even
    if ((n_simpson/2)*2.ne.n_simpson) n_simpson = n_simpson + 1
    nz = n_simpson + 1

    allocate(mjs(nmjs), xe(nz), xe_fid_array(nz), integrand(nz, nmjs), STAT=status) !might want it as a derived parameter
		if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate): memory allocation failed'
			stop
		endif 

    allocate(tmp(nmjs), STAT=status) !might want it as a derived parameter
		
    ! Compute z grid for simpson integration, and evaluate basis functions there
    zmin = this%kde_xe_basis%zmin
    zmax = this%kde_xe_basis%zmax
    step_size = (zmax-zmin)/dfloat(n_simpson)

    ! Make xe table 
    do i = 1, nz
      z = zmin + (dfloat(i)-1) * step_size
      xe(i) = custom_xe(transformed_parameters, z) 
      xe_fid_array(i) = xe_fiducial(this%kde_xe_basis, z)
      call this%kde_xe_basis%eval_basis(z, integrand(i, 1:nmjs))
    end do

    ! Calculate integrand and performs integral for the projection
    do j = 1, nmjs
      do i = 1, nz
        integrand(i, j) = integrand(i, j) * (xe(i) - xe_fid_array(i))
      end do
      call simpson(integrand(:,j), step_size, n_simpson, mjs(j))
    end do
    mjs = mjs/(zmax - zmin) 

  end subroutine get_mjs
  
  subroutine Kde_CalcDerivedParams(parameters, derived_parameters)

    real(mcp), dimension(:), intent(inout) :: parameters
    real(mcp), dimension(:), allocatable :: transformed_parameters 
    real(mcp), dimension(:), allocatable :: mjs
    real(mcp), dimension(:), allocatable, intent(out) :: derived_parameters
    logical :: do_print_derived = .false.
    integer :: status, row
    integer :: j, num_params

    allocate(mjs(5), transformed_parameters(1), derived_parameters(6), STAT=status) !might want to add mjs as derived parameters
		if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate): memory allocation failed'
			stop
		endif

    ! Get mjs for xe(z) with parameters in param_array 
    call get_transformed_parameters(parameters, transformed_parameters)

    call ThisKde%get_mjs(transformed_parameters, mjs)

    num_params = 0 

    do j = 1, 1
      num_params = num_params + 1
      derived_parameters(num_params) = transformed_parameters(j)
    end do

    do j = 1, 5
      num_params = num_params + 1
      derived_parameters(num_params) = mjs(j)
    end do
    
    if (do_print_derived .eqv. .true.) then
      do j = 1, 6
        print *, 'derived_parameters(j) = ', derived_parameters(j)
      end do
    end if

  end subroutine Kde_CalcDerivedParams

  subroutine RelikeKde_OneModel(parameters, likelihood, derived_parameters) 
    
    real(mcp), dimension(:), intent(inout) :: parameters
    real(mcp), intent(out) :: likelihood
    real(mcp), dimension(:), allocatable :: mjs
    real(mcp), dimension(:), allocatable, intent(out) :: derived_parameters
    integer :: status, j

    allocate(mjs(5), STAT=status) 
		if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate): memory allocation failed'
			stop
		endif

    call Kde_CalcDerivedParams(parameters, derived_parameters)

    do j = 1, 5
      mjs(j) = derived_parameters(1+j)
    end do

    if (ThisKde%kde_params%kde_or_gaussian == 'kde') then
      call ThisKde%kde_impl(mjs, likelihood)
    else if (ThisKde%kde_params%kde_or_gaussian == 'gaussian') then
      call ThisKde%gauss_impl(mjs, likelihood)
    end if
    
  end subroutine RelikeKde_OneModel

  subroutine gauss_impl(this, mjs, result)

    real(mcp), dimension(:), intent(in) :: mjs
    real(mcp), intent(out) :: result
    class(RelikeKde) :: this

    integer :: sz, row, col
    integer :: status
    integer nmjs
    real(mcp) prod, weight, norm
    real(mcp), dimension(:), allocatable :: dmjs

    allocate(dmjs(5), STAT=status) 
		
    dmjs = mjs - this%mean
    prod = dot_product(dmjs, matmul(this%invcov * this%kde_params%fcov, dmjs)) 
    result = exp(-0.5_mcp * prod) 

  end subroutine gauss_impl

  subroutine kde_impl(this, mjs, result)

    real(mcp), dimension(:), intent(in) :: mjs
    real(mcp), intent(out) :: result
    class(RelikeKde) :: this

    integer :: sz, row, col
    integer :: status
    integer nmjs
    real(mcp) prod, weight, norm
    real(mcp), dimension(:), allocatable :: dmjs
  
    allocate(dmjs(5), STAT=status) 
		
    result = 0.0_mcp
    do row = 1, ThisKdeChains%nrow 
      dmjs = mjs - ThisKdeChains%data(row, 2:2+nmjs)
      prod = dot_product(dmjs, matmul(this%invcov, dmjs)) 
      weight = ThisKdeChains%data(row, 1)
      result = result + exp(-0.5_mcp * prod) * weight 
    end do

  end subroutine kde_impl

  subroutine setup_invcov(this)
    class(RelikeKde) :: this
    real(mcp), dimension(:,:), allocatable :: invcov
    integer :: status
    integer :: nrow, ncol
    integer :: row, col, i, j
    
    open(unit = 100, file = this%kde_params%fn_invcov_spec)
    read(100,*) nrow, ncol
    close(100)
  
    allocate(this%invcov(nrow, ncol), STAT=status)
    if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate.f90): memory allocation failed.'
			stop
		endif

    open(unit = 101, file = this%kde_params%fn_invcov)   
    read(101, *, iostat=status) ((this%invcov(row, col), col = 1, ncol), row = 1, nrow)
    if(status .ne. 0) then
      write (*,*) 'Error (KernelDensityEstimate.f90): reading inverse covariance file failed.'
      stop
    endif
    close(101)   
    
    this%invcov = this%invcov/this%kde_params%fcov

  end subroutine setup_invcov

  subroutine setup_mean(this)
    class(RelikeKde) :: this
    real(mcp), dimension(:), allocatable :: mean
    integer :: status
    integer :: nrow=5, ncol
    integer :: row, col, i, j
    
    allocate(this%mean(nrow), STAT=status)
    if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate.f90): memory allocation failed.'
			stop
		endif

    open(unit = 101, file = this%kde_params%fn_mean)   
    read(101, *, iostat=status) (this%mean(row), row = 1, nrow)
    if(status .ne. 0) then
      write (*,*) 'Error (KernelDensityEstimate.f90): reading PC mean file failed.'
      stop
    endif
    close(101)   
    
  end subroutine setup_mean

  subroutine setup_chains(this, fn_chain, fn_spec)
    class(RelikeKdeChains) :: this
    character(LEN=1024), intent(in) :: fn_chain
    character(LEN=1024), intent(in) :: fn_spec
    integer :: nrow, ncol
    integer :: row, col
    integer :: i, j
    real(mcp), dimension(:,:), allocatable :: data
    integer :: status

    open(unit = 100, file=fn_spec)
    read(100,*) nrow, ncol
    close(100)
  
    this%nrow = nrow
    this%ncol = ncol
    allocate(this%data(this%nrow, this%ncol),STAT=status)
    if(status .ne. 0) then
			write (*,*) 'Error (KernelDensityEstimate.f90): memory allocation failed.'
			stop
		endif

    print *, '... Reading chain file'  
    open(unit = 101, file=fn_chain)    
    read(101, *, iostat=status) ((this%data(row, col), col = 1, ncol), row = 1, nrow)
    if(status .ne. 0) then
      write (*,*) 'Error (KernelDensityEstimate.f90): reading chain file failed.'
      stop
    endif
    close(101)   
    
  end subroutine setup_chains 

  subroutine RelikeKde_Init(Kde, KdeChains) ! TODO: might need more input params
    
    use constants
    Type(RelikeKde), target :: Kde
    Type(RelikeKdeChains), target :: KdeChains
    
    ThisKde => Kde
    ThisKdeChains => KdeChains

    call ThisKde%init()

    ! No need to setup chains if using Gaussian approximation
    if (Kde%kde_params%kde_or_gaussian == 'kde') then
      call ThisKdeChains%setup_chains(Kde%kde_params%fn_chain, Kde%kde_params%fn_chain_spec)
    end if
    
  end subroutine RelikeKde_Init

  subroutine simpson(table, step_size, n, integral) 
    !==========================================================
    ! Integration of f(x) on [a,b] 
    ! where a = table(1) and b = table(n+1)
    ! Method: Simpson rule for n intervals  
    !----------------------------------------------------------
    ! IN:
    ! table       - table of f(x) to integrate (supplied by a user)
    ! step_size   - step size dx
    ! n           - number of intervals (not number of points)
    ! OUT:
    ! integral    - Result of integration
    !==========================================================
    implicit none

    real(mcp), intent(out) :: integral
    real(mcp), dimension(:), intent(in) :: table
    real(mcp), intent(in) :: step_size
    integer, intent(in) :: n 

    real(mcp) :: s
    integer i 

    s = 0.0
    do i = 3, n-1, 2
      s = s + 2.0*table(i) + 4.0*table(i+1)
    end do
    
    integral = (s + table(1) + table(n+1) + 4.0*table(2)) * step_size/3.0

    return
  end subroutine simpson

end module KernelDensityEstimate

