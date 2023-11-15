module funcs

  use shared
  use myIntegrate
  
  implicit none


contains  
    
  !!!!! Point-Coulomb Function
  function Coulomb(rad)

    use shared  

    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: Coulomb

    Coulomb = -alpha*Z/rad

  end function Coulomb



  !!!!! Function with the integrand for the charge density
  function ChargeNormInt(rad)

    use shared
    use myIntegrate

    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: ChargeNormInt

    ChargeNormInt = 4.0_dp*pi_d*rad**2*rho2(rad)

  end function ChargeNormInt


  !!!!! Defining the charge density, properly normalised here
  function ChargeRho(rad)

    use shared
    use myIntegrate

    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: ChargeRho
    real(dp) :: a1, b1, Norm
    
    !!! ChargeNorm definined in shared.f90
    !!! The normalisation was calculated by integrating the above function to infty, then the result was set as a parameter.
    ChargeRho = Z*rho2(rad)/ChargeNorm
    
  end function Chargerho


  !!! Defining the integrand for the first r' integral of the finite Coulomb potential
  function Int1(x)

    use shared

    real(dp), intent(in) :: x
    real(dp)             :: Int1

    Int1 = x*chargerho(x)

  end function Int1


!!! Defining the integrand for the second q integral of the finite Coulomb potential
!!! We note this integrand is defined as an integral
!!! Here the integral interfaces to the integral function defined in myIntegrate
  !!! Using this it integrates a function of the type sin(qx)*f(x), and thus provides a way of defining an integral for a parameter q, unique to our problem.
  function Int2(q)

    use shared
    use myIntegrate

    real(dp), intent(in) :: q
    real(dp)             :: Int2, a, b, FiniteCoulomb, omega, x

!!! integral limits.
!!! We set the lower limit to a small value, but do not wish to integrate to values of b too large, as this leads to instability in the sin term.
    
    a = 1.0e-3_dp
    b = 10.0_dp


    Int2 = integral(Int1, a, b, q, abs_err, rel_err)/q**2
    
  end function Int2
    
!!! Function to calculate the finite Coulomb potential, V(r).
  function FiniteCoulomb(rad)

    use shared
    use myIntegrate

    real(dp), intent(in) ::rad

    real(dp) :: a, b, FiniteCoulomb, omega

!!! integral limits.
!!! We set the lower limit to a small value, but do not wish to integrate to values of b too large, as this leads to instability in the sin term.
    a = 1.0e-3_dp
    b = 10.0_dp!huge(1.0_dp)
    
    FiniteCoulomb = -8.0_dp*alpha*integral(Int2, a, b, rad, abs_err, rel_err)/rad
    
  end function FiniteCoulomb


!!! Defining the nucleon density, which is here a Woods-Saxon form.
function rho(rad)

    use shared

    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: rho
    real(dp) :: c1,c2

    rho = rho0/(1.0_dp+exp((rad-R_c*Anum**(1.0_dp/3.0_dp))/a0))
    
  end function rho


  !!! Defining the proton density
  function rho2(rad)

    use shared
    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: rho2
 
    rho2 = rho0/(1.0_dp + exp((rad-R_c*Z**(1.0_dp/3.0_dp))/a0))

  end function rho2

!!! Defining the term for the Isospin deoendent interaction
!!! Here this is effectively rho_p-rho_n
  
  function RhoIso(rad)

    use shared
    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: RhoIso

    RhoIso = rho2(rad) - rho0/(1.0_dp + exp((rad-R_c*(Anum-Z)**(1.0_dp/3.0_dp))/a0))
  end function RhoIso
  

!!! One of the derivative terms as defined in the potential used in the Schrodinger Equation case.
  function diffrho(rad)

    use shared

    implicit none

    real(dp), intent(in) :: rad
    real(dp)             :: diffrho, arg

      
    arg = exp((rad-R_c*(Anum)**(1.0_dp/3.0_dp))/a0)
    diffrho = rho0*(-2.0_dp*a0*(1.0_dp+arg) - rad*(1.0_dp-arg))*arg/(a0*2.0_dp*rad*(1.0_dp+arg)**3.0_dp)

  end function diffrho


!!! One of the derivative terms as defined in the potential used in the Schrodinger Equation case.
  function diff2rho(rad)

    use shared

    implicit none

    real(dp), intent(in) :: rad
    real(dp)             :: diff2rho, arg

      
    arg = exp((rad-R_c*(Anum)**(1.0_dp/3.0_dp))/a0)
    
    diff2rho =2.0_dp*rho0**2.0_dp*(-2.0_dp*a0*(1.0_dp+arg)-rad*(1.0_dp-arg))*arg/(a0**2.0_dp*rad*(1.0_dp+arg)**4.0_dp)
    
  end function diff2rho
  
!!! One of the derivative terms as defined in the potential used in the Schrodinger Equation case.
  function diffrhosquare(rad)

    use shared

    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: diffrhosquare, arg

    arg = exp((rad-R_c*(Anum)**(1.0_dp/3.0_dp))/a0)

    diffrhosquare = -2.0_dp*rho0*arg/(a0*(1.0_dp+arg)**3.0_dp)
    
  end function diffrhosquare


!!! The laplacian of the Woods-Saxon nucleon density
  function LaplaceRho(rad)

    use shared
    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: LaplaceRho, arg

    arg = exp((rad-R_c*(Anum)**(1.0_dp/3.0_dp))/a0)

    LaplaceRho = rho0*(2.0_dp*arg**2.0_dp-arg*(1.0_dp+arg))/(a0**2.0_dp*(1.0_dp+arg)**3.0_dp) -2.0_dp/(rad*a0)*arg*(1.0_dp+arg)/(1.0_dp+arg)**3.0_dp
  end function LaplaceRho


!!! The nuclear potential as used in the Schrodinger equation Case  
  function NuclearPot(rad)

    use shared
    implicit none
    real(dp), intent(in) :: rad
    real(dp)             :: OM1, T1, T2, T3, b, drho_drho, d2rho, drho2_drho,  NuclearPot

    
    OM1 = wom*Gom
    T1 = -w1*Gsig
    T2 = d*(w2/2.0_dp+w1)*Gsig**2.0_dp
    T3 = -d**2.0_dp*(w1+w2)*Gsig**3.0_dp

    drho_drho = -d*Gsig**2.0_dp*(w1+w2)*diffrho(rad)**2.0_dp/(msig**2.0_dp)
    
    drho2_drho = Gsig**3.0_dp*d**2.0_dp*(w1+w2)*diffrho(rad)*diffrhosquare(rad)/(msig**2.0_dp)
    
    d2rho = Gsig**3.0_dp*d**2.0_dp*(w1+2.0_dp*w2)*rho(rad)*diff2rho(rad)/(msig**2.0_dp)
    
    b = Grho/4.0_dp

    NuclearPot = OM1*rho(rad) + T1*rho(rad) + T2*rho(rad)**2.0_dp + T3*rho(rad)**3.0_dp - b*rhoiso(rad) + drho_drho + drho2_drho + d2rho

  end function NuclearPot
  

!!! The potential due to the scalar meson interaction
  function ScalarPot(rad)

    use shared
    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: ScalarPot

    ScalarPot = w1*( Gsig*rho(rad)*(1.0_dp-d*Gsig*rho(rad)+d**2.0_dp*Gsig**2.0_dp*rho(rad)**2.0))-w1*Gsig/(msig**2.0_dp)*laplacerho(rad)

  end function ScalarPot


!!! The potential due to the vector meson interactions
  function VectorPot(rad)

    use shared
    implicit none

    real(dp), intent(in) :: rad
    real(dp) :: VectorPot

    vectorPot = wom*(Gom*rho(rad)) - Grho/4.0_dp*rhoiso(rad)-wom*Gom/mom**2.0_dp*laplacerho(rad)
    
  end function VectorPot
  
  
    
!!!! Below we have the solutions to the Point Coulomb equations, for states n=1, ..., 6, for both the F and G states.


  !!! EXACT FUNCTIONS n=1

  function exact_F1(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F1

    real(dp), parameter :: gamma = gamma1
    real(dp), parameter :: euler_func = euler_gamma_func1
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((2.0_dp*(1.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F1 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F1

  

  function exact_G1(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G1
    
    real(dp), parameter :: gamma = gamma1
    real(dp), parameter :: euler_func = euler_gamma_func1
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((2.0_dp*(1.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G1 = A*(1.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G1



  !!!! EXACT FUNCTIONS n=2


 function exact_F2(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F2

    real(dp), parameter :: gamma = gamma2
    real(dp), parameter :: euler_func = euler_gamma_func2
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((4.0_dp*(2.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F2 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F2

  

  function exact_G2(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G2
    
    real(dp), parameter :: gamma = gamma2
    real(dp), parameter :: euler_func = euler_gamma_func2
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((4.0_dp*(2.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G2 = A*(2.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G2
  
!!! Exact Functions n=3

  function exact_F3(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F3

    real(dp), parameter :: gamma = gamma3
    real(dp), parameter :: euler_func = euler_gamma_func3
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((6.0_dp*(3.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F3 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F3

  

  function exact_G3(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G3
    
    real(dp), parameter :: gamma = gamma3
    real(dp), parameter :: euler_func = euler_gamma_func3
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((6.0_dp*(3.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G3 = A*(3.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G3
 

!!! EXACT FUNCTION N=4

  function exact_F4(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F4

    real(dp), parameter :: gamma = gamma4
    real(dp), parameter :: euler_func = euler_gamma_func4
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((8.0_dp*(4.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F4 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F4

  

  function exact_G4(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G4
    
    real(dp), parameter :: gamma = gamma4
    real(dp), parameter :: euler_func = euler_gamma_func4
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((8.0_dp*(4.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G4 = A*(4.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G4
  

!!! EXACT FUNCTIONS n=5

  function exact_F5(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F5

    real(dp), parameter :: gamma = gamma5
    real(dp), parameter :: euler_func = euler_gamma_func5
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((10.0_dp*(5.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F5 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F5

  

  function exact_G5(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G5
    
    real(dp), parameter :: gamma = gamma5
    real(dp), parameter :: euler_func = euler_gamma_func5
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((10.0_dp*(5.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G5 = A*(5.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G5


!!! Exact functions n=6

  function exact_F6(r)

    use shared
    
    implicit none

    real(dp),  intent(in) :: r
    real(dp) :: exact_F6

    real(dp), parameter :: gamma = gamma6
    real(dp), parameter :: euler_func = euler_gamma_func6
    real(dp) :: A, C
    
    
    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((12.0_dp*(6.0_dp+gamma))) * sqrt(C/(euler_func))
 
    exact_F6 = A*Z*alpha*r*(2.0_dp*C*r)**(gamma)*exp(-C*r)/sqrt(4.0_dp*pi_d)

  end function exact_F6

  

  function exact_G6(r)

    use shared

    implicit none

    real(dp), intent(in) :: r
    real(dp) :: exact_G6
    
    real(dp), parameter :: gamma = gamma6
    real(dp), parameter :: euler_func = euler_gamma_func6
    
    real(dp) :: A, C


    C = Z*alpha*m_reduced
    A = 1.0_dp/sqrt((12.0_dp*(6.0_dp+gamma))) * sqrt(C/euler_func)

    exact_G6 = A*(6.0_dp+gamma)*(2.0_dp*C*r)**(gamma)*exp(-C*r)*r/sqrt(4.0_dp*pi_d)    

  end function exact_G6
  
end module funcs
