module shared

  use nrtype

  
  implicit none

  !! (CODATA)
  real(dp), parameter :: fermi = 197.3269631_dp                             ! hbar*c (MeV fm) 
  real(dp), parameter :: m_Iron =  55.934_dp*931.494_dp/fermi               ! m_iron (fm) 
  real(dp), parameter :: m_Cascade = 1321.71_dp/fermi                       ! m_Xi (fm) 
  real(dp), parameter :: m_reduced = (m_Cascade*m_Iron)/(m_Cascade+m_Iron)  ! Reduced mass (fm)

  real(dp), parameter :: alpha = 0.0072973525376_dp ! fine structure constant
  real(dp), parameter :: Z = 26.0_dp                ! Proton Number of Iron
  real(dp), parameter :: Anum = 56.0_dp             ! Nucleon Number of Iron

  !!! Coupling Constants from QMC model, fitted to the binding energy of Lambda hyperons
  real(dp), parameter :: Gsig = 8.65_dp        !!! Coupling constant to sigma meson (fm^2)
  real(dp), parameter :: Gom = 5.64_dp         !!! Coupling Constant to omega meson (fm^2)
  real(dp), parameter :: Grho = 4.71_dp        !!! Coupling Constant to rho meson (fm^2)
  !!! Parameters from the QMC paper, taken from Guichon et. al. 2018, Quark-Meson-Coupling model for finite nuclei, nuclear matter and beyond
  real(dp), parameter :: w1 = 0.3400_dp        !!! Fitting to reproduce cascade mass
  real(dp), parameter :: w2 = 0.3688_dp        !!! Second constant to reproduce cascade mass
  real(dp), parameter :: wom = 1.0_dp/3.0_dp   !!! Constant modifying coupling of omega meson to cascade
  real(dp), parameter :: d = 0.176_dp          !!! Scalar Polarisability (fm)
  
  
  real(dp), parameter :: msig = 504.0_dp/fermi ! Mass of sigma meson (fm)
  real(dp), parameter :: mom = 787.0_dp/fermi  ! Mass of omega meson (fm)
  real(dp), parameter :: mrho = 770.0_dp/fermi ! Mass of rho meson (fm)

  !!! Charge normalisation, such that int(Charge_dens)=Z
  real(dp), parameter :: ChargeNorm = 29.022146894129477_dp
  
!!! DEFINING GAMMA AND EULER_GAMMA_FUNC(1+2Gamma), for use in the point Coulomg solutions.
!!! The number 1-6 corresponds to the n=1, ..., 6 solution
!!! The euler_gamma_func is calculated externally using wolfram alpha
!!! If you want to use this for a different nuclei, you will need to recalculated all these values!
  
  real(dp), parameter :: gamma1 = sqrt(1.0_dp-(Z*alpha)**2.0_dp)
  real(dp), parameter :: euler_gamma_func1 = 1.934572274986814_dp
  
  real(dp), parameter :: gamma2 = sqrt(4.0_dp-(Z*alpha)**2.0_dp)
  real(dp), parameter :: euler_gamma_func2 = 23.3575454077_dp
  
  real(dp), parameter :: gamma3 = 2.99399433611_dp ! gamma = SQRT[k^2-z^2*alpha^2]
  real(dp), parameter :: euler_gamma_func3 = 703.99147588_dp! Calculated fully

  real(dp), parameter :: gamma4 = 0.9818360783_dp ! gamma = SQRT[k^2-z^2*alpha^2]
  real(dp), parameter :: euler_gamma_func4 = 39550.4425533_dp! Calculated fully

  real(dp), parameter :: gamma5 = 4.996398911683_dp ! gamma = SQRT[k^2-z^2*alpha^2]
  real(dp), parameter :: euler_gamma_func5 = 3.56786274744e6_dp! Calculated fully

  real(dp), parameter :: gamma6 =5.996999423434_dp ! gamma = SQRT[k^2-z^2*alpha^2]
  real(dp), parameter :: euler_gamma_func6 =4.71795908500e8_dp! Calculated fully


  !!! Constants for the Woods-Saxon Density
  real(dp), parameter :: R_c = 1.1_dp    ! (fm)
  real(dp), parameter :: a0  = 0.6_dp    ! (fm)
  real(dp), parameter :: rho0 = 0.15_dp  ! (fm^-3)
  
  
  integer(I8B), parameter :: n_steps = 150*1e5_I8B ! number of steps (150 fm)
  integer(I8B), parameter :: match   = 1500*1e3_I8B!150*1e3_I8B  ! match point number (1.5 fm)

  integer, parameter :: max_attempts = 1e4  ! max number of attempts to solve

  real(dp), parameter :: h_step  = 1e-5_dp!  ! step size = 1/100000 fm

  real(dp), parameter :: yErrReq = 1.0e-6_dp      ! required precision (!) 
  !!! Error for the integral solvers
  real(dp), parameter :: abs_err = 1.0e-10_dp 
  real(dp), parameter :: rel_err = 0.0_dp

  !!! Number of states to be solved for
  integer, parameter :: num_states = 6

  !!! Number of nodes for a given solution
  integer, dimension(num_states), parameter :: expected_nodes = (/ 0, 0, 0, 0, 0, 0 /)

  !!! The eigenvalue of kappa for each state
  integer, dimension(num_states), parameter :: kappa = (/ -1, -2, -3, -4, -5, -6 /)

  real(dp), dimension(6), parameter :: n_state = (/ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp /)
  
  !! Globals: !!

  integer  :: state
  real(dp), dimension(:,:), allocatable :: F_state, G_state

 
end module shared
