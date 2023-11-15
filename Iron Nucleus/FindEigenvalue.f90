program FindEigenvalue

  use shared
  use funcs
  use myIntegrate
  
  implicit none

  !! Define Quantities: !!

  !! CODATA(98)

  real(dp) ::radius

  real(dp) :: eigenvalue(num_states), test_eigenvalue
  real(dp) :: error_est(num_states)
  integer  :: nodes(num_states)
  integer  :: current_sign, new_sign
  
  integer  :: r_it, attempt, guess ! iterators

  real(dp) :: scalefactor, sum_integral, Fdiff, deltaEV
  real(dp) :: Fdiff_up, Fdiff_down, Fderiv
  real(dp) :: exact_solution(6)

  integer :: accept_ev

  real(dp) :: V_Coulomb(n_steps), Rval(n_steps)
  
  real(dp) :: proportional_factor, EVshift

  allocate(F_state(n_steps,num_states))
  allocate(G_state(n_steps,num_states))

  
  proportional_factor = h_step


  !!! Code commented out to calculate the finite Coulomb potential and save it to a file called FiniteCoulombval.dat
  
  !open(10, file='FiniteCoulombVal.dat', status='unknown') 
  !do r_it = 1, n_steps
  !   write(10, *) r_it*h_step, FiniteCoulomb(r_it*h_step)
  !end do
  !close(10)
 
!!! Reading in Coulomb Potential:
  open(10, file = 'FiniteCoulombVal.dat', status = 'old')
  do r_it = 1, n_steps
     read(10, *) Rval(r_it), V_Coulomb(r_it)
  end do
  close(10)

  
  !! Solver: !!
  !! Iteratively solve Dirac equation to find E_e = E - m
  
  !!! The exact solutions to the point Coulomb potetnail
  exact_solution(1) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp))         - m_reduced  ! E - m 
  exact_solution(2) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp/4.0_dp))  - m_reduced  ! E - m
  exact_solution(3) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp/9.0_dp))  - m_reduced  ! E - m
  exact_solution(4) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp/16.0_dp)) - m_reduced  ! E - m
  exact_solution(5) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp/25.0_dp)) - m_reduced  ! E - m
  exact_solution(6) = m_reduced*sqrt(1.0_dp-((Z*alpha)**2.0_dp/36.0_dp)) - m_reduced  ! E - m

  !!! Writing the exact solutions to the terminal
  do state = 1, num_states 
     write(*,'(A,I1,A,F13.7,A,F18.12,A)')"state ",state," exact = ",exact_solution(state)," = ",exact_solution(state)*fermi," MeV"
  end do

  !!! Setting the proportional factor initially to the step size
  proportional_factor = h_step


!!! Initialising the eigenvalues to the solution
  eigenvalue(:)   = (1.0_dp-proportional_factor)*exact_solution(:) ! Setting initial guess
!!! Initialising the F and G states
  F_state(:,:)    = 0.0_dp
  G_state(:,:)    = 0.0_dp



!!! Solving for the 5g and 6h states. If you wish to solve for other solutions you might need to adjust the match point in shared.f90
!!! This is due to the wavefunction spreading out as you go to higher angular momentum states
  do state = 5, 6! num_states

     do attempt = 1, max_attempts

        write(*,'(a,F18.12,a)')"eigenvalue = ",eigenvalue(state)*fermi," MeV"
        write(*,'(a,F18.12,a)')"exact soln = ",exact_solution(state)*fermi," MeV"
        ! calculate F and G
        test_eigenvalue = eigenvalue(state)
        call Runge_up(test_eigenvalue, V_Coulomb)
        call Runge_down(test_eigenvalue, V_Coulomb)
        
        ! scale the solutions to make the 
        ! match point smooth
        scalefactor = G_state(match+1,state)/G_state(match,state)
        G_state(:match,state) = G_state(:match,state)*scalefactor
        F_state(:match,state) = F_state(:match,state)*scalefactor
        
        ! calculate integral |F|^2 + |G|^2 (extended Simpson's Rule)
        sum_integral = 0.0_dp
        sum_integral = sum_integral + h_step*(abs(F_state(1,state))**2 + abs(G_state(1,state))**2)
        sum_integral = sum_integral + h_step*(abs(F_state(n_steps,state))**2 + abs(G_state(n_steps,state))**2)
        
        do r_it = 2, n_steps-2, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 2.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        do r_it = 3, n_steps-1, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 4.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        sum_integral = sum_integral/3.0_dp
        ! normalize
        F_state(:,state) = F_state(:,state)/sqrt(abs(sum_integral))
        G_state(:,state) = G_state(:,state)/sqrt(abs(sum_integral))
        ! calculate discontinuity
        Fdiff = (F_state(match+1,state) - F_state(match,state))/F_state(match,state)

        !!! Calculate the energy shift
        EVshift = -G_state(match,state)*Fdiff*F_state(match,state)
                
        print *,"EV shift = ",EVshift*fermi," MeV"


        !!! Determining when the energy shift is small enough, otherwise add EVshift to eigenvalue and repeat
        if ( abs(EVshift) < yErrReq ) then
           print *,"EV shift is small enough now... done."
           write(*,'(a,e12.6,a,e12.6)')"EVshift = ",abs(Evshift)," vs yErrReq = ",yErrReq
           go to 53
        end if

        eigenvalue(state) = eigenvalue(state) + EVshift


        !!! Recording the wavefunction for plotting
        open(61,file='Wavefunctions.dat',status='unknown')
        do r_it = 1, n_steps,1000
           write(61,*)r_it*h_step,F_state(r_it,state),G_state(r_it,state)
        end do
        close(61)
        !
        if ( attempt == max_attempts ) STOP ' too many attempts '
        !
     end do
    



     
     !!! Code to determine when solution to high precision
53   continue

     do guess = 1, 200

        write(*,'(a,F18.12,a)')"eigenvalue = ",eigenvalue(state)*fermi," MeV"
        write(*,'(a,F18.12,a)')"exact soln = ",exact_solution(state)*fermi," MeV"

        F_state(:,state) = 0.0_dp
        G_state(:,state) = 0.0_dp

        test_eigenvalue = eigenvalue(state)
     
        call Runge_up(test_eigenvalue, V_Coulomb)
        call Runge_down(test_eigenvalue, V_Coulomb)
    

        ! scale the solutions to make the
        ! match point smooth
        scalefactor = G_state(match+1,state)/G_state(match,state)
        G_state(:match,state) = G_state(:match,state)*scalefactor
        F_state(:match,state) = F_state(:match,state)*scalefactor
        ! calculate integral |F|^2 + |G|^2 (extended Simpson's Rule)
        sum_integral = 0.0_dp
        sum_integral = sum_integral + h_step*(abs(F_state(1,state))**2 + abs(G_state(1,state))**2)
        sum_integral = sum_integral + h_step*(abs(F_state(n_steps,state))**2 + abs(G_state(n_steps,state))**2)
        do r_it = 2, n_steps-2, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 2.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        do r_it = 3, n_steps-1, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 4.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        sum_integral = sum_integral/3.0_dp
        ! normalize
        F_state(:,state) = F_state(:,state)/sqrt(abs(sum_integral))
        G_state(:,state) = G_state(:,state)/sqrt(abs(sum_integral))
        ! calculate discontinuity
        Fdiff = (F_state(match+1,state) - F_state(match,state))/F_state(match,state)

        F_state(:,state) = 0.0_dp
        G_state(:,state) = 0.0_dp

        test_eigenvalue = (1.0_dp+proportional_factor)*eigenvalue(state)
        call Runge_up(test_eigenvalue, V_Coulomb)
        call Runge_down(test_eigenvalue, V_Coulomb)
 

        ! scale the solutions to make the 
        ! match point smooth
        scalefactor = G_state(match+1,state)/G_state(match,state)
        G_state(:match,state) = G_state(:match,state)*scalefactor
        F_state(:match,state) = F_state(:match,state)*scalefactor
        ! calculate integral |F|^2 + |G|^2 (extended Simpson's Rule)
        sum_integral = 0.0_dp
        sum_integral = sum_integral + h_step*(abs(F_state(1,state))**2 + abs(G_state(1,state))**2)
        sum_integral = sum_integral + h_step*(abs(F_state(n_steps,state))**2 + abs(G_state(n_steps,state))**2)
        do r_it = 2, n_steps-2, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 2.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        do r_it = 3, n_steps-1, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 4.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        sum_integral = sum_integral/3.0_dp
        ! normalize
        F_state(:,state) = F_state(:,state)/sqrt(abs(sum_integral))
        G_state(:,state) = G_state(:,state)/sqrt(abs(sum_integral))
        ! calculate discontinuity
        Fdiff_up = (F_state(match+1,state) - F_state(match,state))/F_state(match,state)

        F_state(:,state) = 0.0_dp
        G_state(:,state) = 0.0_dp

        test_eigenvalue = (1.0_dp-proportional_factor)*eigenvalue(state)
        call Runge_up(test_eigenvalue, V_Coulomb)
        call Runge_down(test_eigenvalue, V_Coulomb)
     

        ! scale the solutions to make the 
        ! match point smooth
        scalefactor = G_state(match+1,state)/G_state(match,state)
        G_state(:match,state) = G_state(:match,state)*scalefactor
        F_state(:match,state) = F_state(:match,state)*scalefactor
        ! calculate integral |F|^2 + |G|^2 (extended Simpson's Rule)
        sum_integral = 0.0_dp
        sum_integral = sum_integral + h_step*(abs(F_state(1,state))**2 + abs(G_state(1,state))**2)
        sum_integral = sum_integral + h_step*(abs(F_state(n_steps,state))**2 + abs(G_state(n_steps,state))**2)
        do r_it = 2, n_steps-2, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 2.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        do r_it = 3, n_steps-1, 2
           !
           radius = r_it*h_step
           sum_integral = sum_integral &
                & + 4.0_dp*h_step*(abs(F_state(r_it,state))**2 + abs(G_state(r_it,state))**2)
           !
        end do
        sum_integral = sum_integral/3.0_dp
        ! normalize
        F_state(:,state) = F_state(:,state)/sqrt(abs(sum_integral))
        G_state(:,state) = G_state(:,state)/sqrt(abs(sum_integral))
        ! calculate discontinuity
        Fdiff_down = (F_state(match+1,state) - F_state(match,state))/F_state(match,state)

        write(*,*)"central = ",Fdiff
        write(*,*)"up      = ",Fdiff_up
        write(*,*)"down    = ",Fdiff_down


        !!! Determining where the solution is in a minimum
        if ( ( abs(Fdiff_down) < abs(Fdiff) ) .and. ( abs(Fdiff_up) < abs(Fdiff) ) ) then
           !
           if ( abs(Fdiff_up) .le. abs(Fdiff_down) ) then
              eigenvalue(state) = (1.0_dp+proportional_factor)*eigenvalue(state)
           end if
           if ( abs(Fdiff_down) .le. abs(Fdiff_up) ) then
              eigenvalue(state) = (1.0_dp-proportional_factor)*eigenvalue(state)
           end if
        else
           if ( abs(Fdiff_up) < abs(Fdiff) ) then
              eigenvalue(state) = (1.0_dp+proportional_factor)*eigenvalue(state)
           end if
           if ( abs(Fdiff_down) < abs(Fdiff) ) then
              eigenvalue(state) = (1.0_dp-proportional_factor)*eigenvalue(state)
           end if
           if ( ( abs(Fdiff_down) > abs(Fdiff) ) .and. ( abs(Fdiff_up) > abs(Fdiff) ) ) then
              ! "Solution is in a minimum"
              proportional_factor = 0.1_dp*proportional_factor
           end if
        end if
        
        if ( abs(proportional_factor) < yErrReq ) then
           print *,"proportional factor is small enough now... done"
           go to 51
        end if

        if ( guess == 200 ) STOP "can't find EV"
        !
     end do ! guess

!!! Included here is code to ensure the correct number of nodes are present in the solution
     
51   continue
     !
     write(*,*)"THE EVALUE OF STATE ",state," IS ",eigenvalue(state)*fermi," MeV"
     write(*,*)"THE EXACT SOL OF STATE  ",state," IS ",exact_solution(state)*fermi," MeV"

     !!! Comment out the below line if you want to include checking for the correct number of nodes
     go to 52

     !! COUNT NODES:
     if ( m_reduced - eigenvalue(state) .gt. 0.0_dp ) then
        nodes(state) = 0
        current_sign = int(G_state(10,state)/abs(G_state(10,state)))
        do r_it = 11, n_steps
           if ( abs(G_state(r_it,state)) > 0.0_dp ) new_sign = int(G_state(r_it,state)/abs(G_state(r_it,state)))
           if ( new_sign /= current_sign ) then 
              nodes(state) = nodes(state) + 1
              current_sign = new_sign
           end if
        end do
     else
        nodes(state) = 99
     end if
     !
     write(*,*)" This state appears to have ",nodes(state)," nodes"
     !!
     if ( ( nodes(state) == expected_nodes(state) ) ) then
        accept_ev = 1
     else
        accept_ev = 0
     end if
     !
     if ( accept_ev == 1 ) then
        write(*,*)"<< CORRECT NUMBER OF NODES! >>"
        go to 52
        !
     else if ( accept_ev == 0 ) then
        ! Incorrect number of nodes - change starting point
        write(*,*)"INCORRECT NUMBER OF NODES!"
        eigenvalue(state) = eigenvalue(state)*(1.05_dp)
     end if
     !
     !
52   continue
     !
     write(*,*)"THE EIGENVALUE OF STATE ",state," IS ",eigenvalue(state)*fermi," MeV"
     !
     open(61,file='Wavefunctions.dat',status='unknown')
     do r_it = 1, n_steps, 1000
        write(61,*)r_it*h_step,F_state(r_it,state),G_state(r_it,state)
     end do
     close(61)
     !
  end do ! state

  write(*,*) "The Transition Energy from (6,5) to (5,4) is", (eigenvalue(5)-eigenvalue(6))*fermi, "MeV"
  
end program FindEigenvalue



!!! This subroutine contains the matrix that defines the derivatives for the F and G functions
!!! The key parameter is the NuclearPotCheck, which differentiates between different solutions
!!! If set to 0, the solver will solve for the Point Coulomb Case
!!! If set to 1, the solver will solve for just the finite Coulomb Case
!!! If set to 2, the solver will solve for the finite Coulomb + vector potential - scalar potential, as definied in funcs.f90
!!! If set to 3, the solver will solve for the nuclear_pot + finite Coulomb, as defined in funcs.f90
subroutine myDerivs(r, params, VC, dydr)

  use shared
  use funcs

  implicit none

  real(dp), intent(in)  :: r, params(3), VC
  real(dp), intent(out) :: dydr(3)
  !
  real(dp) :: F_local, G_local, lambda
  !
  integer, parameter :: NuclearPotCheck = 1
  F_local = params(1)
  G_local = params(2)
  lambda  = params(3)
  !
  !!! Solving for Point Coulomb
  if (NuclearPotCheck == 0) then
     dydr(1) = F_local*kappa(state)/r - (lambda - Coulomb(r))*G_local
     dydr(2) = -G_local*kappa(state)/r + (lambda + 2.0_dp*m_reduced -Coulomb(r))*F_local
  end if
  !!! Solving for Finite Coulomb
  if (NuclearPotCheck == 1) then
     dydr(1) =  F_local*kappa(state)/r - (lambda - VC)*G_local 
     dydr(2) = -G_local*kappa(state)/r + (lambda + 2.0_dp*m_reduced - VC)*F_local
  end if
  !!! Solving for Finite Coulomb + Vector Pot - Scalar Pot
  if (NuclearPotCheck == 2) then
     dydr(1) = F_local*kappa(state)/r - (lambda+scalarPot(r) - VC- VectorPot(r))*G_local
     dydr(2) = -G_local*kappa(state)/r + (lambda + 2.0_dp*m_reduced - scalarPot(r) - VC-vectorpot(r))*F_local
  end if
  !!! Solving for Finite Coulomb + Nuclear Pot
  if (NuclearPotCheck == 3) then
     dydr(1) =  F_local*kappa(state)/r - (lambda - VC - NuclearPot(r))*G_local 
     dydr(2) = -G_local*kappa(state)/r + (lambda + 2.0_dp*m_reduced - VC-nuclearpot(r))*F_local
  end if
  
  !!! Final term is related to the value for Lambda, which is set to zero here.
  dydr(3) = 0.0_dp
  !
end subroutine myDerivs


!!! Code to solve for the left side solution to the wavefunction
subroutine Runge_up(eigen, V_Coul)

  use shared
  use funcs

  implicit none

  real(dp), intent(in) :: eigen, V_Coul(n_steps)

  real(dp) :: tempF(n_steps),  tempG(n_steps), tempL(n_steps)
  real(dp), dimension(3) :: KF1, KF2, KF3, KF4, temp_params, tempdydr

  real(dp) :: temp_rad, newrad
  real(dp) :: newF, newG, newL

  real(dp) :: V0

  integer :: ki

  tempL(1) = eigen

  !!! Initialising the solver using the exact solutions to the point Coulomb case
  do ki = 1, 5
     !
     temp_rad = h_step*ki
     V0 = Coulomb(temp_rad)
     !

      if ( state == 1 ) F_state(ki, state) = exact_F1(temp_rad)
      if ( state == 2 ) F_state(ki, state) = exact_F2(temp_rad)
      if ( state == 3 ) F_state(ki, state) = exact_F3(temp_rad)
      if ( state == 4 ) F_state(ki, state) = exact_F4(temp_rad)
      if ( state == 5 ) F_state(ki, state) = exact_F5(temp_rad)
      if ( state == 6 ) F_state(ki, state) = exact_F6(temp_rad)

     tempF(ki) = F_state(ki,state)


     if ( state == 1 ) G_state(ki, state) = exact_G1(temp_rad)
     if ( state == 2 ) G_state(ki, state) = exact_G2(temp_rad)
     if ( state == 3 ) G_state(ki, state) = exact_G3(temp_rad)
     if ( state == 4 ) G_state(ki, state) = exact_G4(temp_rad)
     if ( state == 5 ) G_state(ki, state) = exact_G5(temp_rad)
     if ( state == 6 ) G_state(ki, state) = exact_G6(temp_rad)

     tempL(ki) = eigen
     !

  end do

  !!! Continuing to integrate out, using the derivatives to iterate. This is done using the Runge Kutta method
  do ki = 5, match-1

     temp_params = (/ tempF(ki), tempG(ki), tempL(ki) /)

     call myDerivs(temp_rad,temp_params, V_Coul(ki), tempdydr)

     KF1 = tempdydr

     newrad = temp_rad + h_step/2.0_dp
     newF = tempF(ki) + h_step*KF1(1)/2.0_dp
     newG = tempG(ki) + h_step*KF1(2)/2.0_dp
     newL = tempL(ki) + h_step*KF1(3)/2.0_dp

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad,temp_params, V_Coul(ki), tempdydr)

     KF2 = tempdydr

     newrad = temp_rad + h_step/2.0_dp
     newF = tempF(ki) + h_step*KF2(1)/2.0_dp
     newG = tempG(ki) + h_step*KF2(2)/2.0_dp
     newL = tempL(ki) + h_step*KF2(3)/2.0_dp

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad,temp_params, V_Coul(ki), tempdydr)

     KF3 = tempdydr

     newrad = temp_rad + h_step
     newF = tempF(ki) + h_step*KF3(1)
     newG = tempG(ki) + h_step*KF3(2)
     newL = tempL(ki) + h_step*KF3(3)

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad,temp_params, V_Coul(ki), tempdydr)

     KF4 = tempdydr

     tempF(ki+1) = tempF(ki) + h_step*(KF1(1)+2.0_dp*KF2(1)+2.0_dp*KF3(1)+KF4(1))/6.0_dp
     tempG(ki+1) = tempG(ki) + h_step*(KF1(2)+2.0_dp*KF2(2)+2.0_dp*KF3(2)+KF4(2))/6.0_dp
     tempL(ki+1) = tempL(ki) + h_step*(KF1(3)+2.0_dp*KF2(3)+2.0_dp*KF3(3)+KF4(3))/6.0_dp

     F_state(ki+1,state) = tempF(ki+1)
     G_state(ki+1,state) = tempG(ki+1)

     temp_rad = temp_rad + h_step

  end do

end subroutine Runge_up


!!! Code to solve for the right side solution to the wavefunction
subroutine Runge_down(eigen, V_Coul)

  use shared
  use funcs

  implicit none

  real(dp), intent(in) :: eigen, V_Coul(n_steps)

  real(dp) :: tempF(n_steps),  tempG(n_steps), tempL(n_steps)
  real(dp), dimension(3) :: KF1, KF2, KF3, KF4, temp_params, tempdydr


  real(dp) :: temp_rad, newrad
  real(dp) :: newF, newG, newL

  integer :: ki

  !!! Initialising the solutions using the exact solutions to the Point Coulomb Case

  do ki = n_steps, n_steps-5, -1
     !
     temp_rad = h_step*ki
     !
     if ( state == 1 ) F_state(ki,state) = exact_F1(temp_rad)
     if ( state == 2 ) F_state(ki,state) = exact_F2(temp_rad)
     if ( state == 3 ) F_state(ki,state) = exact_F3(temp_rad)
     if ( state == 4 ) F_state(ki,state) = exact_F4(temp_rad)
     if ( state == 5 ) F_state(ki,state) = exact_F5(temp_rad)
     if ( state == 6 ) F_state(ki,state) = exact_F6(temp_rad)

     tempF(ki) = F_state(ki,state)

     if ( state == 1 ) G_state(ki,state) = exact_G1(temp_rad)
     if ( state == 2 ) G_state(ki,state) = exact_G2(temp_rad)
     if ( state == 3 ) G_state(ki,state) = exact_G3(temp_rad)
     if ( state == 4 ) G_state(ki,state) = exact_G4(temp_rad)
     if ( state == 5 ) G_state(ki,state) = exact_G5(temp_rad)
     if ( state == 6 ) G_state(ki,state) = exact_G6(temp_rad)

     tempG(ki) = G_state(ki,state)

     tempL(ki) = eigen
     !
  end do

  !!! Iterating inwards towards the match point, using the derivatives to iterate
  do ki = n_steps-5, match+2, -1

     temp_params = (/ tempF(ki), tempG(ki), tempL(ki) /)

     call myDerivs(temp_rad, temp_params, V_coul(ki), tempdydr)

     KF1 = tempdydr

     newrad = temp_rad - h_step/2.0_dp
     newF = tempF(ki) - h_step*KF1(1)/2.0_dp
     newG = tempG(ki) - h_step*KF1(2)/2.0_dp
     newL = tempL(ki) - h_step*KF1(3)/2.0_dp

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad, temp_params,V_Coul(ki), tempdydr)

     KF2 = tempdydr

     newrad = temp_rad - h_step/2.0_dp
     newF = tempF(ki) - h_step*KF2(1)/2.0_dp
     newG = tempG(ki) - h_step*KF2(2)/2.0_dp
     newL = tempL(ki) - h_step*KF2(3)/2.0_dp

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad, temp_params, V_Coul(ki), tempdydr)

     KF3 = tempdydr

     newrad = temp_rad - h_step
     newF = tempF(ki) - h_step*KF3(1)
     newG = tempG(ki) - h_step*KF3(2)
     newL = tempL(ki) - h_step*KF3(3)

     temp_params = (/ newF, newG, newL /)

     call myDerivs(newrad, temp_params, V_Coul(ki), tempdydr)

     KF4 = tempdydr

     tempF(ki-1) = tempF(ki) - h_step*(KF1(1)+2.0_dp*KF2(1)+2.0_dp*KF3(1)+KF4(1))/6.0_dp
     tempG(ki-1) = tempG(ki) - h_step*(KF1(2)+2.0_dp*KF2(2)+2.0_dp*KF3(2)+KF4(2))/6.0_dp
     tempL(ki-1) = tempL(ki) - h_step*(KF1(3)+2.0_dp*KF2(3)+2.0_dp*KF3(3)+KF4(3))/6.0_dp

     F_state(ki-1,state) = tempF(ki-1)
     G_state(ki-1,state) = tempG(ki-1)

     temp_rad = temp_rad - h_step

  end do
  
end subroutine Runge_down
