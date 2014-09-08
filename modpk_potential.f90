MODULE potential
  USE modpkparams
  use internals, only : pi
  use modpk_qsf
  use modpk_numerics
  use modpk_errorhandling, only : raise
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot, getH, getdHdalpha, getEps, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, geteta, zpower, getH_with_t, stability_check_on_H, &
       getEps_with_t, d3Vdphi3, &
       logP_of_observ, fisher_rao_metric, guess_EOI_field

  public :: norm
  public :: bundle, field_bundle

  type :: bundle
    real(dp) :: N=0e0_dp
    real(dp) :: exp_scalar=0e0_dp
    real(dp) :: dlogThetadN=0e0_dp
    contains
      procedure, public :: calc_exp_scalar =>bundle_exp_scalar
  end type bundle

  type(bundle) :: field_bundle


CONTAINS


  recursive function pot(phi) result(V_potential)
    !
    !     Returns V(phi) given phi, phi can be an array for multifield,
    !     The code implement multifield potential in the form of V = \sum V(phi_i),
    !     More complicated form of potentials can be customized
    !
    real(dp) :: V_potential
    real(dp), intent(in) :: phi(:)

    real(dp) :: m2_V(size(phi)) !! m2_V is the diagonal mass square matrix
    real(dp) :: c1_V(size(phi)) !! some constants for a more flexible potential construction in case 9
    real(dp) :: lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    !real(dp), dimension(size(phi)/2) :: lambda_waterfall, mass_waterfall, &
    !  mass_infl, couple_water_infl
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(num_inflaton)

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix
    integer :: i,j, temp_choice

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope

    real(dp) :: m_light2, M_heavy2, param0, param_closest, dist,&
      phi_light, phi_light0

    real(dp) :: p_exp

    real(dp) :: V0
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: ii, jj

    select case(potential_choice)

    case(1)
      ! m_i^2 phi_i^2 --- N-quadratic
       m2_V = 10.e0_dp**(vparams(1,:))
       V_potential = 0.5e0_dp*sum(m2_V*phi*phi)

    case(2)
      ! N-flation (axions)
       lambda = 10.e0_dp**vparams(1,:)
       finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
       V_potential = sum(lambda**4*(1.e0_dp+cos(finv*phi)))

    case(3)
       lambda = 10.e0_dp**vparams(1,:)
       V_potential = sum(0.25e0_dp*lambda*phi**4)

    case(4)
       lambda = 10.e0_dp**vparams(1,:)
       V_potential = sum(lambda*phi)

    case(5)
       lambda = 10.e0_dp**vparams(1,:)
       V_potential = sum(lambda*1.5e0_dp*phi**(2.e0_dp/3.e0_dp))

    case(6)
       lambda = 10.e0_dp**vparams(1,:)
       mu = 10.e0_dp**vparams(2,:)
       V_potential = sum(lambda**4 - mu*phi**4/4.e0_dp)

    case(7) !product of exponentials
       lambda = vparams(2, :)
       V_potential = vparams(1,1)*M_Pl**4 * exp(dot_product(lambda, phi/M_Pl))

    case(8)
       !Canonical two-field hybrid
       if (size(phi) /= 2) then
         print*, "Potential_choice =", Potential_choice
         print*, "Number of fields =", size(phi)
         call raise%fatal_cosmo(&
             "This potential requires two fields.  &
             Set num_inflaton=2.", &
             __FILE__, __LINE__)
       end if

       lambda_hybrid = vparams(1,1)
       mass_hybrid = vparams(1,2)
       mu_hybrid =vparams(1,3)
       nu_hybrid =vparams(1,4)
       V_potential = (lambda_hybrid**4)*((1.0_dp - phi(1)**2/mass_hybrid**2)**2 +&
         phi(2)**2/mu_hybrid**2 +&
         phi(1)**2*phi(2)**2/nu_hybrid**4)

    case(9)
      m2_V = vparams(1,:)
      c1_V = vparams(2,:)
      V_potential = 0.5e0_dp*sum(m2_V*phi*phi) + sum(c1_V)

    case(10)
      print*, "MODPK: Potential_choice=", Potential_choice
      call raise%fatal_code(&
        "This potential choice is broken for some reason.",&
        __FILE__, __LINE__)

    case(11)
      !N-quadratic w/one quartic intxn term
      !phi_i^2 + phi_{lightest}^2*phi_i^2

      m2_V = 10.0e0_dp**(vparams(1,:))
      lambda4 = 10.0e0_dp**vparams(2,:)
      phi_light_index = minloc(m2_V,1)

      V_potential = 0.5e0_dp*sum(m2_V*phi*phi)+ &
        (1.0e0_dp/24.0e0_dp)*(phi(phi_light_index)**2)*sum(lambda4*phi*phi)

    case(12)
      !Mass matrix with diagonal terms = m_i^2
      !Off-diagonal terms = \eps
      m2_V = 10.e0_dp**(vparams(1,:))
      lambda2 = 10.e0_dp**(vparams(2,1))
      V_potential = 0e0_dp

      do i=1, num_inflaton
        do j=1, num_inflaton

          if (i==j) then
            m2_matrix(i,j) = m2_V(i)
            V_potential = V_potential +0.5e0_dp* m2_matrix(i,j)*phi(i)*phi(j)
          else
            m2_matrix(i,j) = lambda2
            V_potential = V_potential + m2_matrix(i,j)*phi(i)*phi(j)
          end if

        end do
      end do

    case(13)
      !Quasi--single-field/turning trajectory
      !phi(1) = phi_light

      !vparams(1,:)  = Veff vparams
      !vparams(2,:)  = \lambda_i
      !vparams(3,:)  = \alpha_i
      !vparams(4+,:) = turn_funct params

      !Only need vparams(X,2:num_inflaton)
      lambda4  = 10.0e0_dp**vparams(2,:)
      alpha2 = vparams(3,:)

      !Run through potential with "effective" single-field potential_choice
      temp_choice = potential_choice
      potential_choice = effective_V_choice
      V_potential = pot((/phi(1)/))
      potential_choice = temp_choice

      !Choice of function to use set in parameters_multimodecode.txt
#define PHI_I phi(1), i, turning_choice(i-1)
#define DELTAPHI (phi(i) - turning_function(PHI_I))
#define EXPTERM exp( (-0.5e0_dp/alpha2(i))*(phi(i) - turning_function(PHI_I))**2)

      do i=2, num_inflaton
        V_potential = V_potential - lambda4(i)*&
          (EXPTERM  - 1.0e0_dp)
      end do

    case(14)

      !Multifield step potential
      m2_V = 10.e0_dp**(vparams(1,:))
      location_phi = vparams(2,:)
      step_size = vparams(3,:)
      step_slope = vparams(4,:)

      !Check for /0 error
      if (any(abs(step_slope) < 1e-15)) then
        print*, "step_slope=", step_slope
        call raise%fatal_code(&
          "This is a division by zero error.  &
          Set the slope in tanh(phi/step_slope) greater than about 1e-15.",&
          __FILE__, __LINE__)
      end if

      V_potential=0e0_dp

      do i=1,num_inflaton
        V_potential = V_potential + &
            0.5*m2_V(i)*(phi(i)**2) * &
            (1.0e0_dp + step_size(i)* &
            tanh( (phi(i)-location_phi(i))/step_slope(i)))
      end do

    case(15)

      !"Numerical" QSF trajectory
      !(1/2)m_L^2 phi_L^2 + (1/2)M_heavy^2 D^2
      !for D^2 = MIN( sum_i (phi_i - funct_i(param))^2; wrt param)
      !Where param_min is found numerically

      m_light2 = 10.0e0_dp**vparams(1,1)
      M_heavy2 = 10.0e0_dp**vparams(1,2)

      param0 = vparams(2,1)

      !Find the parameter that gives the minimum distance
      !between phi and the parametrized curve
      if (allocated(qsf_runref%phi)) then
        if ( size(qsf_runref%phi)/=num_inflaton) then
          deallocate(qsf_runref%phi)
          allocate(qsf_runref%phi(num_inflaton))
        end if
      else
          allocate(qsf_runref%phi(num_inflaton))
      end if
      qsf_runref%phi = phi

      !Closest approach to parameterized curve
      param_closest = qsf_runref%min_dist()
      dist = distance(param_closest)

      !Check not getting too far off line
      if (dist >5e-1_dp) then
      !if (dist >1e0_dp) then
        print*, "QSF: dist =", dist
        print*, "QSF: param_closest =", param_closest
        print*, "QSF: phi =", phi
        print*, "QSF: turning_function_parametric =", &
          turning_function_parametric(param_closest)

        call raise%fatal_cosmo(&
          "The trajectory is significantly deviating &
           from the parametric curve. &
           Check that you have actually found the closest parameter &
           or that you're not taking steps that are too large.", &
        __FILE__, __LINE__)

      end if

      !Get the integrated distance this closest point is up the curve
      phi_light = qsf_runref%phi_light(param_closest)

      !DEBUG
      if (param_closest < param0 .or. &
        phi_light > qsf_runref%phi_light_vs_param(size(qsf_runref%phi_light_vs_param,1),1)) then
        print*, "param_closest = ", param_closest
        print*, "phi = ", phi
        print*, "phi_light = ", phi_light
        print*, "qsf_runref%param = ", qsf_runref%param
        print*, "dist = ", dist
        print*, "dist2 = ", distance(4.0e0_dp)
        stop
      end if

      !Reset param guess for next time through
      qsf_runref%param = param_closest

      V_potential = 0.5e0_dp*m_light2*phi_light**2 &
        + 0.5e0_dp*M_heavy2*dist**2

    case(16)
      ! (1/p) lambda_i |phi_i|^p --- N-monomial


      p_exp = vparams(2,1)
      m2_V = vparams(1,:)
      V_potential = (1.0e0_dp/p_exp)*sum(m2_V*abs(phi)**p_exp)

    case(17)
      ! Generalized axions

      V0 = vparams(1,1)
      A_i = vparams(2,:)
      B_ij = vparams(3:num_inflaton+2,:)

      V_potential = V0
      do ii=1,size(A_i)
        V_potential = V_potential + A_i(ii)*cos(2.0e0_dp*pi*phi(ii))
      end do

      do ii=1,size(A_i); do jj=1, size(A_i)
        V_potential = V_potential + B_ij(ii,jj)*&
          cos(2.0e0_dp*pi*phi(ii) - 2.0e0_dp*pi*phi(jj))
      end do; end do


    case default
      print*, "potential_choice =", potential_choice
      call raise%fatal_code(&
          "Need to set the potential V(phi) for this potential choice.",&
          __FILE__, __LINE__)

    end select

  END FUNCTION pot


  !Needs "result" because array-valued and recursive.
  recursive function dVdphi(phi) result(first_deriv)
    !
    !     Returns dV/dPhi given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: first_deriv(size(phi))

    real(dp) :: dphi(size(phi)), phiplus(size(phi))
    real(dp) :: m2_V(size(phi)),  c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i,j, temp_choice

    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(size(phi))

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope

    real(dp), dimension(size(phi)) :: stepsize
    real(dp), dimension(:), allocatable :: numderiv

    real(dp) :: param_closest
    real(dp) :: m_light2, M_heavy2
    real(dp) :: phi_light, dphi_light
    real(dp), dimension(num_inflaton) :: turnfunct, dturnfunct, &
      dparam_dphi

    real(dp) :: p_exp
    integer :: ii

    real(dp) :: V0, twopi
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: alpha, beta

    if (vnderivs) then
       ! MULTIFIELD
       do i=1, num_inflaton
          phiplus = phi
          phiplus(i) = phi(i) + 0.5e0_dp*phi(i)*findiffdphi**(1.0e0_dp/3.0e0_dp)
          dphi = phiplus - phi
          first_deriv(i) = (pot(phi+dphi)-pot(phi-dphi))/(2.0e0_dp*dphi(i))
          if (first_deriv(i).eq.0.e0_dp) then

            call raise%fatal_code(&
             'first_deriv(i)=0, possibly signaling a problem with &
             accuracy of numerical derivatives. &
             Try using vnderivs=F if possible.',&
             __FILE__, __LINE__)
          end if
       end do
       ! END MULTIFIELD
    else
       select case(potential_choice)
       !MULTIFIELD
       case(1)
          m2_V = 10.e0_dp**(vparams(1,:))
          first_deriv = m2_V*phi
       case(2)
          lambda = 10.e0_dp**vparams(1,:)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
          first_deriv = -lambda**4*finv*sin(finv*phi)
       case(3)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda*phi**3
       case(4)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda
       case(5)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda*phi**(-1./3.)
       case(6)
          mu = 10.e0_dp**vparams(2,:)
          first_deriv = -mu*phi**3
       case(7)
          lambda = vparams(2,:)
          first_deriv = lambda*vparams(1,1)*M_Pl**3 * exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          first_deriv(1) = (lambda_hybrid**4)*(4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)&
            *phi(1)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)
          first_deriv(2) = (lambda_hybrid**4)*(2.0_dp*phi(2)/mu_hybrid**2 +&
            2.0_dp*phi(1)**2*phi(2)/nu_hybrid**4)
       case(9)!Saddle type things
          m2_V = (vparams(1,:))
          first_deriv = m2_V*phi
       case(10)
         print*, "MODPK: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(11)
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda4 = 10.e0_dp**vparams(2,:)
         phi_light_index = minloc(m2_V,1)

         !for i /= lightest
         first_deriv = m2_V*phi + (1.0e0_dp/12.0e0_dp)*(phi(phi_light_index)**2)*lambda4*phi

         !i=lightest
         first_deriv(phi_light_index) = m2_V(phi_light_index)*phi(phi_light_index) + &
           (1.0e0_dp/6.0e0_dp)*(phi(phi_light_index)**3)*lambda4(phi_light_index)

       case(12)

         !Mass matrix with diagonal terms = m_i^2
         !Off-diagonal terms = \eps
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda2 = 10.e0_dp**(vparams(2,1))
         first_deriv = 0e0_dp

         do i=1, num_inflaton
           do j=1, num_inflaton
             if (i==j) then
               m2_matrix(i,j) = m2_V(i)
             else
               m2_matrix(i,j) = lambda2
             end if

             first_deriv(i) = first_deriv(i) + m2_matrix(i,j)*phi(j)
           end do
         end do

       case(13)

         lambda4  = 10.0e0_dp**vparams(2,:)
         alpha2 = vparams(3,:)

#define HEAVY 2:num_inflaton
         !deriv wrt phi_light
         temp_choice = potential_choice
         potential_choice = effective_V_choice
         first_deriv = dVdphi((/phi(1)/))
         first_deriv(HEAVY) = 0e0_dp
         potential_choice = temp_choice

         do i=2,num_inflaton
          first_deriv(1) = first_deriv(1) &
            - (lambda4(i)/alpha2(i)) * &
            (phi(i) - turning_function(PHI_I)) * &
            dturndphi(PHI_I) * &
            EXPTERM
         end do

         !deriv wrt heavy fields
         do i=2, num_inflaton
           first_deriv(i) = (lambda4(i)/alpha2(i)) *&
             (phi(i) -  turning_function(PHI_I)) * &
             EXPTERM
         end do

       case(14)

         !Multifield step potential
         m2_V = 10.e0_dp**(vparams(1,:))
         location_phi = vparams(2,:)
         step_size = vparams(3,:)
         step_slope = vparams(4,:)


         do i=1,num_inflaton

          first_deriv(i) = m2_V(i)*phi(i)* &
            (1.0e0_dp + step_size(i)* &
            tanh( (phi(i)-location_phi(i))/step_slope(i))) + &
            0.5e0_dp*m2_V(i)*phi(i)**2*step_size(i)* &
            (MySech((phi(i)-location_phi(i))/step_slope(i)))**(2)/&
            step_slope(i)

         end do

       case(15)
         !Numerical QSF

         select case(turning_choice(1))
         case (:3)
           m_light2 = 10.0e0_dp**vparams(1,1)
           M_heavy2 = 10.0e0_dp**vparams(1,2)

           qsf_runref%phi = phi

           !Find point of closest approach
           param_closest = qsf_runref%min_dist()

           !Reset param guess for next time through
           qsf_runref%param = param_closest

           !Get the integrated distance this closest point is up the curve
           !and the derivs of the integrated distance
           phi_light = qsf_runref%phi_light(param_closest)
           dphi_light = qsf_runref%dphi_light_dparam(param_closest)

           !Get turn funct and its derivs
           turnfunct = turning_function_parametric(param_closest)
           dturnfunct = dturndparam(param_closest)

           dparam_dphi = dparam_closest_dphi(param_closest)

           !Build dVdphi
           first_deriv = m_light2*phi_light*dphi_light*dparam_dphi &
             + M_heavy2*(phi - turnfunct) &
             - M_heavy2*dparam_dphi*sum( dturnfunct*(phi - turnfunct))

         case default

           !Take derivs of potential directly via central diff
           stepsize = 1.0e-8_dp
           call num_first_deriv(pot, phi, stepsize, numderiv)
           first_deriv = numderiv

         end select

       case(16)
         ! (1/p) lambda_i |phi_i|^p --- N-monomial

         p_exp = vparams(2,1)
         m2_V = vparams(1,:)
         first_deriv = m2_V*abs(phi)**(p_exp-1.0e0_dp)*sign(1.0e0_dp,phi)

         !Regularize around phi=0
         do ii=1,size(phi)
           if (abs(phi(ii))<1e-5_dp) &
             first_deriv(ii)=0e0_dp
         end do

       case(17)
         ! Generalized axions

         V0 = vparams(1,1)
         A_i = vparams(2,:)
         B_ij = vparams(3:num_inflaton+2,:)

         twopi = 2.0e0_dp*pi

         first_deriv = -twopi*sin(twopi*phi)

         do alpha=1, size(first_deriv)
           do beta=1, size(first_deriv)
             if (alpha==beta) cycle
             first_deriv(alpha) = first_deriv(alpha) + &
               twopi*(B_ij(alpha,beta)+B_ij(beta,alpha))* &
               sin(twopi*phi(beta) - twopi*phi(alpha))
           end do
         end do


       !END MULTIFIELD
       case default

         print*, "potential_choice =", potential_choice
         call raise%fatal_code(&
           "Need to set first derivative for this potential choice &
           or use numerical derivatives (vnderivs=.true.)",&
           __FILE__, __LINE__)

       end select

    end if

  END FUNCTION dVdphi


  !Needs "result" because array-valued and recursive.
  recursive function d2Vdphi2(phi) result(second_deriv)
    !
    !     Returns d^2V/dPhi^2 given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: second_deriv(size(phi),size(phi))

    real(dp) :: m2_V(size(phi)), c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i,j, temp_choice

    real(dp) :: dphi,phiplus
    !  mass_infl, couple_water_infl
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(size(phi))

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope
    real(dp), dimension(size(phi)) :: stepsize
    real(dp), dimension(:,:), allocatable :: numderiv
    real(dp), dimension(:), allocatable :: numderiv1

    real(dp) :: param_closest
    real(dp) :: m_light2, M_heavy2
    real(dp) :: phi_light, dphi_light, d2phi_light
    real(dp), dimension(num_inflaton) :: turnfunct, dturnfunct, d2turnfunct, &
      dparam_dphi
    real(dp), dimension(num_inflaton, num_inflaton) :: d2param_dphi2, delta
    integer :: ii, jj, ll

    real(dp) :: p_exp

    real(dp) :: V0, twopi, fourpi_sq
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: alpha, beta

    if (vnderivs) then
       !MULTIFIELD
       if (size(phi) .ne. 1) then
         call raise%fatal_code(&
          'The 2nd order numerical derivative has not &
          been implemented for more than one field.',&
          __FILE__, __LINE__)
       end if
       phiplus = phi(1) + 0.2e0_dp*phi(1)*findiffdphi**(1.e0_dp/4.e0_dp)
       dphi = phiplus - phi(1)
       second_deriv = (pot(phi+2.e0_dp*dphi)+pot(phi-2.e0_dp*dphi)- &
         2.e0_dp*pot(phi))/(4.e0_dp*dphi*dphi)
    else
       second_deriv(:,:) = 0e0_dp

       select case(potential_choice)
       case(1)
          m2_V = 10.e0_dp**(vparams(1,:))
          forall (i = 1:size(phi)) second_deriv(i,i) = m2_V(i)
       case(2)
          lambda = 10.e0_dp**vparams(1,:)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
          forall (i = 1:size(phi)) second_deriv(i,i) = -lambda(i)**4*finv(i)*finv(i)*cos(finv(i)*phi(i))
       case(3)
          lambda = 10.e0_dp**vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = 3.e0_dp*lambda(i)*phi(i)**2
       case(4)
          forall (i = 1:size(phi)) second_deriv(i,i) = 0.e0_dp
       case(5)
          lambda = 10.e0_dp**vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = -lambda(i)/3.e0_dp*phi(i)**(-4./3.)
       case(6)
          mu = 10.e0_dp**vparams(2,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = -3.e0_dp*mu(i)*phi(i)**2
       case(7)
          lambda = vparams(2, :)
          forall (i=1:size(phi), j=1:size(phi)) second_deriv(i,j) = &
               lambda(i)*lambda(j)*vparams(1,1)*M_Pl**2 *exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          second_deriv(1,1) = (lambda_hybrid**4)*(4.0_dp*(2.0_dp*phi(1)/mass_hybrid**2)*phi(1)/mass_hybrid**2 +&
          4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)

          second_deriv(2,1) = (lambda_hybrid**4)*(4.0_dp*phi(1)*phi(2)/nu_hybrid**4)

          second_deriv(1,2) = second_deriv(2,1)

          second_deriv(2,2) = (lambda_hybrid**4)*(2.0_dp/mu_hybrid**2 +&
            2.0_dp*phi(1)**2/nu_hybrid**4)
       case(9)
          m2_V = vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = m2_V(i)
       case(10)
         print*, "MODPK: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(11)
         M2_V = 10.e0_dp**(vparams(1,:))
         Lambda4 = 10.e0_dp**vparams(2,:)
         Phi_light_index = minloc(m2_V,1)

         do i=1,num_inflaton
           do j=1, num_inflaton

             if (i/=Phi_light_index .and. j/= Phi_light_index) then
               if (i/=j) then
                 second_deriv(i,j) = 0e0_dp
               else
                 second_deriv(i,j) = m2_V(i) + &
                   (1.0e0_dp/12.0e0_dp)*Lambda4(i)*phi(Phi_light_index)**2
               end if

             else if (i==Phi_light_index .and. j/=Phi_light_index) then

               second_deriv(i,j) = 0e0_dp

             else if (i/=Phi_light_index .and. j==Phi_light_index) then

               second_deriv(i,j) = (1.0e0_dp/6.0e0_dp)*Lambda4(i)* &
                 phi(Phi_light_index)*phi(i)

             else if (i==Phi_light_index .and. j==Phi_light_index) then

               second_deriv(i,j) = m2_V(Phi_light_index) + &
                 0.5e0_dp*Lambda4(Phi_light_index)*phi(Phi_light_index)**2

             end if

           end do

         end do

       case(12)

         !Mass matrix with diagonal terms = m_i^2
         !Off-diagonal terms = \eps
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda2 = 10.e0_dp**(vparams(2,1))
         second_deriv = 0e0_dp

         do i=1, num_inflaton
           do j=1, num_inflaton
             if (i==j) then
               m2_matrix(i,j) = m2_V(i)
             else
               m2_matrix(i,j) = lambda2
             end if

           end do
         end do

         second_deriv = m2_matrix

       case(13)
         second_deriv = 0e0_dp

         lambda4  = 10.0e0_dp**vparams(2,:)
         alpha2 = vparams(3,:)

         !d2V/dphi_light2
         temp_choice = potential_choice
         potential_choice = effective_V_choice
         second_deriv = d2Vdphi2((/phi(1)/))
         second_deriv(1,HEAVY) = 0e0_dp
         second_deriv(HEAVY,:) = 0e0_dp
         potential_choice = temp_choice

         do i=2,num_inflaton
           second_deriv(1,1) = second_deriv(1,1) &
             - (lambda4(i)/alpha2(i)) * &
             ( -1e0_dp*dturndphi(PHI_I)**2 + &
             (phi(i) - turning_function(PHI_I)) * d2turndphi2(PHI_I) + &
             (1.0e0_dp/alpha2(i)) * &
             ((phi(i) - turning_function(PHI_I))**2) *&
             (dturndphi(PHI_I)**2)) * &
             EXPTERM
         end do

         !d2V/dphi_light dphi_heavy
         do i=2,num_inflaton
           second_deriv(1,i) = (-lambda4(i)/alpha2(i)) * dturndphi(PHI_I) *&
             (1.0e0_dp - ((phi(i) - turning_function(PHI_I))**2)/alpha2(i)) * &
             EXPTERM
         end do
         second_deriv(:,1)=second_deriv(1,:)

         !d2V/dphi_heavy dphi_heavy
         do i=2, num_inflaton; do j=2, num_inflaton
             if (i /= j) then
               second_deriv(i,j) = 0e0_dp
             else
               second_deriv(i,j) = (lambda4(i)/alpha2(i)) * &
                 (1.0e0_dp - ((phi(i) - turning_function(PHI_I))**2)/alpha2(i)) * &
                 EXPTERM
             end if
         end do; end do

       case(14)

         !Multifield step potential
         m2_V = 10.e0_dp**(vparams(1,:))
         location_phi = vparams(2,:)
         step_size = vparams(3,:)
         step_slope = vparams(4,:)

         second_deriv=0e0_dp

         do i=1,num_inflaton

           second_deriv(i,i) = m2_V(i)*&
              (1.0e0_dp + step_size(i)* &
              tanh( (phi(i)-location_phi(i))/step_slope(i))) +&
              2.0e0_dp*m2_V(i)*phi(i)*step_size(i)* &
              (MySech( (phi(i)-location_phi(i))/step_slope(i)))**(2)/&
              step_slope(i) - &
              m2_V(i)*phi(i)**2 * (step_size(i)/step_slope(i)**2) *&
              tanh((phi(i)-location_phi(i))/step_slope(i)) *&
              (MySech((phi(i)-location_phi(i))/step_slope(i)))**(2)

         end do

       case(15)

         select case(turning_choice(1))
         case (:3)

           !Numerical QSF
           m_light2 = 10.0e0_dp**vparams(1,1)
           M_heavy2 = 10.0e0_dp**vparams(1,2)

           qsf_runref%phi = phi

           !Find point of closest approach
           param_closest = qsf_runref%min_dist()

           !Reset param guess for next time through
           qsf_runref%param = param_closest

           !Get the integrated distance this closest point is up the curve
           !and the derivs of the integrated distance
           phi_light = qsf_runref%phi_light(param_closest)
           dphi_light = qsf_runref%dphi_light_dparam(param_closest)
           d2phi_light = qsf_runref%d2phi_light_dparam2(param_closest)

           !Get turn funct and its derivs
           turnfunct = turning_function_parametric(param_closest)
           dturnfunct = dturndparam(param_closest)
           d2turnfunct = d2turndparam2(param_closest)


           dparam_dphi = dparam_closest_dphi(param_closest)
           d2param_dphi2 = d2param_closest_dphi2(param_closest)

           !Identity matrix
           delta=0e0_dp
           do ii=1, num_inflaton
             delta(ii,ii) = 1.0e0_dp
           end do

           !Build d2Vdphi2
           second_deriv = 0e0_dp
           do ii=1,num_inflaton; do jj=1,num_inflaton

             do ll=1, num_inflaton
               second_deriv(ii,jj) = second_deriv(ii,jj) &
                 - M_heavy2*d2turnfunct(ll)*dparam_dphi(jj)* &
                   (phi(ll)-turnfunct(ll))*dparam_dphi(ii) &
                 - M_heavy2*dturnfunct(ll)*&
                   (delta(ll,jj)-dturnfunct(ll)*dparam_dphi(jj))*dparam_dphi(ii) &
                 - M_heavy2*dturnfunct(ll)*&
                   (phi(ll)-turnfunct(ll))*d2param_dphi2(ii,jj)
             end do

             second_deriv(ii,jj) = second_deriv(ii,jj) &
               + m_light2*(dphi_light**2)*dparam_dphi(ii)*dparam_dphi(jj) &
               + m_light2*phi_light*d2phi_light*dparam_dphi(ii)*dparam_dphi(jj) &
               + m_light2*phi_light*dphi_light*d2param_dphi2(ii,jj) &
               + M_heavy2*(delta(ii,jj) - dturnfunct(ii)*dparam_dphi(jj))

           end do; end do

         case default

           !Take second derivs of potential directly via central diff
           stepsize = 1.0e-4_dp
           call num_second_deriv(pot, phi, stepsize, numderiv)
           second_deriv = numderiv

         end select

       case(16)
         ! (1/p) lambda_i |phi_i|^p --- N-monomial

         p_exp = vparams(2,1)
         m2_V = vparams(1,:)
         second_deriv = 0e0_dp
         do ii=1,size(phi)
           second_deriv(ii,ii) =(p_exp-1.0e0_dp)*m2_V(ii)*abs(phi(ii))**(p_exp-2.0e0_dp)
         end do

         !Regularize around phi=0
         do ii=1,size(phi)
           if (abs(phi(ii))<1e-5) &
             second_deriv(ii,ii)=0e0_dp
         end do

       case(17)
         ! Generalized axions

         V0 = vparams(1,1)
         A_i = vparams(2,:)
         B_ij = vparams(3:num_inflaton+2,:)

         twopi = 2.0e0_dp*pi
         fourpi_sq = 4.0e0_dp*pi**2

         !Identity matrix
         delta=0e0_dp
         do ii=1, num_inflaton
           delta(ii,ii) = 1.0e0_dp
         end do

         !Optimize me...
         do alpha=1, size(A_i)
           do beta=1, size(A_i)
             second_deriv(alpha,beta) = fourpi_sq*(B_ij(alpha,beta)+B_ij(beta,alpha))*&
               cos(twopi*phi(beta) - twopi*phi(alpha)) -&
               fourpi_sq*delta(alpha,beta)*A_i(alpha)*&
               cos(twopi*phi(alpha))

             if (alpha==beta) then
               do ii=1, size(A_i)
                 second_deriv(alpha,beta) = second_deriv(alpha,beta) + &
                   fourpi_sq*(B_ij(alpha,ii) + B_ij(ii,alpha))*&
                   cos(twopi*phi(ii) - twopi*phi(alpha))
               end do
             end if

           end do
         end do

       case default

         print*, "potential_choice =", potential_choice
         call raise%fatal_code(&
           "Need to set second_deriv for this potential choice.",&
            __FILE__, __LINE__)

       end select
       !END MULTIFIELD
    end if

  END FUNCTION d2Vdphi2

#undef HEAVY
#undef PHI_I
#undef DELTAPHI
#undef EXPTERM

  !Only used to get the 3rd order SR parameters for the SR approximation
  !of the scalar running, alpha_s, and if you wish to give an analytical
  !Jacobian for the mode equation evolution.
  function d3Vdphi3(phi) result(third_deriv)
    real(dp), intent(in) :: phi(:)
    real(dp) :: third_deriv(size(phi),size(phi),size(phi))
    real(dp) :: m2_V(size(phi))
    integer :: ii

    real(dp) :: p_exp

    third_deriv = 0e0_dp

    select case(potential_choice)
    case(1)
      m2_V = 10.e0_dp**(vparams(1,:))
      third_deriv=0e0_dp
      do ii=1,size(phi)
        third_deriv(ii,ii,ii)=m2_V(ii)
      end do
    case(15)
      !DEBUG
      print*, "testing --- fake d3Vdphi3 for numerical qsf..."
      third_deriv = 0e0_dp

    case (16)
      ! (1/p) lambda_i |phi_i|^p --- N-monomial

      p_exp = vparams(2,1)
      m2_V = 10.e0_dp**(vparams(1,:))

      third_deriv = 0e0_dp
      do ii=1,size(phi)
        third_deriv(ii,ii,ii) = (p_exp-2.0e0_dp)*(p_exp-1.0e0_dp)*&
          m2_V(ii)*abs(phi(ii))**(p_exp-3.0e0_dp)
      end do

    case default

      print*, "potential_choice =", potential_choice
      call raise%fatal_code(&
        "Need to set third derivative for this potential choice.",&
        __FILE__, __LINE__)

    end select

  end function d3Vdphi3


  FUNCTION initialphi(phi0)
    !
    !     Sets initial value of phi (depending on potential, may use
    !     either the user-specified value phi0 or something derived
    !     from the potential parameters)
    !
    real(dp), INTENT(IN) :: phi0(:)
    real(dp) :: initialphi(size(phi0))

    real(dp) :: phii(size(phi0))
    real(dp) :: Ninit, finv, lambda, mu, phesq
    real(dp) :: x1, x2


    Ninit = 70.e0_dp

    if (size(phi0) .gt. 1) then
       phii = phi0 ! MULTIFIELD
    else  !SINGLE FIELD
       select case(potential_choice)
       case(1)
          phii = 2.e0_dp*sqrt(Ninit+0.5e0_dp)
       case(2)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,1))
          phii = 2.e0_dp/finv*asin(exp(-0.5e0_dp*Ninit*finv*finv)/ &
               sqrt(1.e0_dp+0.5e0_dp*finv*finv))
       case(3)
          phii = sqrt(8.e0_dp*(Ninit+1.e0_dp))
       case(4)
          phii = sqrt(2.e0_dp*Ninit+0.5e0_dp)
       case(5)
          phii = sqrt(4.e0_dp/3.e0_dp*Ninit+2.e0_dp/9.e0_dp)
       case(6)
          lambda = 10.e0_dp**vparams(1,1)
          mu = 10.e0_dp**vparams(2,1)
          x1 = lambda**4/mu
          phesq = ((sqrt(2.e0_dp)*x1)**(-4./3.)+1.e0_dp/(4.e0_dp*x1))**(-0.5)
          if (vparams(1,1)<-3.e0_dp) then
             phii = sqrt(phesq/(2.e0_dp**1.5*Ninit/sqrt(phesq)+1.e0_dp))
          else
             x2 = 4.e0_dp*Ninit + 2.e0_dp*x1/phesq + 0.5e0_dp*phesq
             phii = sqrt(x2)*sqrt(1.e0_dp-sqrt(1.e0_dp-4.e0_dp*x1/x2/x2))
          end if
       case default
          phii = phi0
       end select
    end if

    initialphi = phii

  END FUNCTION initialphi

  FUNCTION getEps(phi,dphi)
    !
    !     Returns epsilon given phi and dphi/dalpha, alpha=ln(a)
    !     for single field, slowroll parameter epsilon_H = 2 M_pl^2 [(dH/dphi)/H]^2
    !     for canonical multi-field,
    !
    real(dp) :: getEps
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    !MULTIFIELD
    getEps = 0.5e0_dp*(M_Pl)**2 * dot_product(dphi,dphi)

    if (getEps >=3.0e0_dp) then
      print*, "MODPK: epsilon =", getEps

      call raise%fatal_cosmo(&
        "Epsilon is >3.0 in SUBROUTINE getEps. &
         This means H is complex.  &
         (This is not universe I'd like to live in.) &
         This error might arise if there is a large separation &
         in scales (stiff problem) and the integrator walks &
         to a bad position in parameter space. &
         Try reducing the integration stepsize.", &
        __FILE__, __LINE__)

    end if

    !END MULTIFIELD

  END FUNCTION getEps


  FUNCTION getH(phi,dphi)
    !
    !     Returns H given phi and dphi/dalpha
    !
    real(dp) :: getH
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    getH=pot(phi)/3.0e0_dp/M_Pl**2 / &
      (1.0e0_dp - dot_product(dphi, dphi)/6.0e0_dp/M_Pl**2)

    if (getH < 0.0e0_dp) then
      call raise%fatal_cosmo(&
        "H is complex. &
        Try smaller stepsize in integrator.",&
        __FILE__,__LINE__)
    else
      getH = sqrt(getH)
    end if

  END FUNCTION getH

  !For when using t-integrator
  FUNCTION getH_with_t(phi,dphidt)
    !
    !     Returns H given phi and dphi/dt
    !
    real(dp) :: getH_with_t
    real(dp), INTENT(IN) :: phi(:), dphidt(:)

    ! MULTIFIELD
    getH_with_t= (pot(phi) + 0.5e0_dp*dot_product(dphidt, dphidt))/ &
      3.0e0_dp/M_pl**2

    if (getH_with_t < 0.0e0_dp) then
      call raise%fatal_cosmo(&
        "H is complex.  &
        Try smaller stepsize in integrator.",&
        __FILE__,__LINE__)
    else
      getH_with_t = sqrt(getH_with_t)
    end if

    ! MULTIFIELD

  END FUNCTION getH_with_t

  FUNCTION getEps_with_t(phi,dphi)
    !
    !     Returns epsilon given phi and dphi/dt
    !
    real(dp) :: getEps_with_t, hubble
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    !MULTIFIELD
    hubble =getH_with_t(phi,dphi)
    getEps_with_t = 0.5e0_dp*(M_Pl)**2 * &
      dot_product(dphi,dphi)/hubble**2

    if (getEps_with_t >3.0e0_dp) then
      print*, "MODPK: epsilon =", getEps_with_t

      call raise%fatal_cosmo(&
        "Epsilon is >3.0 in SUBROUTINE getEps_with_t. &
        This means H is complex.  &
        (This is not universe I'd like to live in.) &
        This error might arise if there is a large separation &
        in scales (stiff problem) and the integrator walks &
        to a bad position in parameter space. &
        Try reducing the integration stepsize.", &
        __FILE__, __LINE__)
    end if
    !END MULTIFIELD


  END FUNCTION getEps_with_t


  FUNCTION getdHdalpha(phi,dphi)
    !
    !     Returns dH/dalpha given phi and dphi/dalpha
    !
    real(dp) :: getdHdalpha
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    ! MULTIFIELD
    getdHdalpha = -dot_product(dphi, dphi) * getH(phi,dphi)/2.0e0_dp/M_Pl**2
    ! END MULTIFIELD
  END FUNCTION getdHdalpha


  FUNCTION getdepsdalpha(phi,dphi)
    !
    !    Returns depsilon/dalpha given phi and dphi/dalpha
    !    Gets this by differentiating Peiris et al Eq A3 (2003)
    !
    real(dp) :: getdepsdalpha, H, dHdalpha, eps
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    H=getH(phi,dphi)
    dHdalpha=getdHdalpha(phi,dphi)
    eps=getEps(phi,dphi)

    ! MULTIFIELD
    getdepsdalpha=6.0e0_dp*dHdalpha/H*(1.0e0_dp-eps/3.0e0_dp) &
      -dot_product(dVdphi(phi), dphi)/(H**2*M_Pl**2)
    ! END MULTIFIELD

  END FUNCTION getdepsdalpha

  FUNCTION geteta(phi, dphi)
    !
    !    Return the eta parameter eta = deps/dalpha / eps
    !
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    real(dp) :: geteta, eps

    eps = getEps(phi, dphi)
    geteta = getdepsdalpha(phi, dphi) / eps

  END FUNCTION geteta

  ! Powerspectrum for psi ptbs mode matrix, outputting full IJ matrix,
  !adiabatic P(k), and isocurv P(k); includes cross correlations
  subroutine powerspectrum(psi, dpsi, phi, dphi, scalefactor, power_spectrum, using_q)
    use internals
    use modpk_observables

    type(power_spectra), intent(inout) :: power_spectrum
    real(dp), dimension(:), intent(in) :: dphi, phi

    logical, optional, intent(in) :: using_q
    logical :: use_q

    ! Hacked matrix to vect
    complex(dp), dimension(:), intent(in) :: psi, dpsi
    real(dp), intent(in) :: scalefactor

    complex(dp), dimension(size(dphi),size(dphi)) :: power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: d_power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: cross_matrix
    real(dp) :: power_adiab
    real(dp) :: power_isocurv
    real(dp) :: power_pnad
    real(dp) :: power_cross
    real(dp) :: power_entropy

    real(dp) :: AAprod, ABprod, BAprod, BBprod
    type(KahanSum) :: pnad_sumAA, pnad_sumAB, pnad_sumBA, pnad_sumBB

    real(dp) :: H, Pdot

    !Pre-factors for Pnad P(k)
    real(dp), dimension(num_inflaton) :: A_vect, B_vect
    real(dp), dimension(num_inflaton) :: A_vect_adiab, B_vect_adiab

    ! size(dphi)=num_inflaton
    complex(kind=dp), dimension(size(dphi),size(dphi)) :: ptb_matrix,&
      dptb_matrix
    real(dp), dimension(size(dphi)-1) :: pk_iso_vect

    !proj along adiab dir
    real(dp), dimension(size(dphi)) :: omega_z
    !proj along isocurv dirs
    real(dp), dimension(size(dphi)-1,size(dphi)) :: s_iso, d_s_iso
    real(dp), parameter :: div_zero_tol=1e-20_dp

    real(dp) :: zeta2, phi_dot_0_scaled, prefactor
    integer :: numb_infl
    integer :: i, j, ll, kk

    real(dp) :: prod_exponent, phi_adiab, d_phi_adiab, hubble
    real(dp), dimension(size(s_iso)) :: phi_iso, d_phi_iso
    real(dp), dimension(size(phi)) :: d_omega_z, Vprime

    real(dp) :: power_total
    real(dp) :: power_pressure, power_press_cross, power_press_cross2, &
      power_press_adiab

    !Variable passed to powerspectrum
    !Psi_ij =a q_ij is default
    !q_ij optional
    if (present(using_q)) then
      use_q = using_q
    else
      use_q =.false.
    end if

    numb_infl=size(dphi)

    ! Convert hacked vector to matrix
    forall (i=1:numb_infl, j=1:numb_infl)
      ptb_matrix(i,j) = psi((i-1)*numb_infl+j)
      dptb_matrix(i,j) = dpsi((i-1)*numb_infl+j)
    end forall

    ! Make projection vector along adiabatic and isocurv directions
    ! NB: phi_dot_0_scaled = sqrt(2*epsilon) = phi_dot_0/H
    phi_dot_0_scaled = sqrt(dot_product(dphi,dphi))
    omega_z = dphi/phi_dot_0_scaled

    !If there's a major hierarchy in scales for the vector components, then
    !there can be a spurious component of the isocurvature direction
    !oriented along the adiabatic one.  When you're approaching the
    !adiabatic limit this can be the dominant contribution, so we force the
    !smallest term to be zero whenever adding it doesn't affect the value of
    !the norm.
    !call renormalize_remove_smallest(omega_z)


    !Cosmology params to make phi_adiab and d_phi_adiab
    hubble = getH(phi,dphi)
    Vprime = dVdphi(phi)
    !domega/dalpha
    d_omega_z = (1e0_dp/(hubble**2)*phi_dot_0_scaled)* &
      (omega_z*sum(omega_z*Vprime) - Vprime )

    phi_adiab = sum(omega_z*phi)
    !d phi_adiab/dalpha
    d_phi_adiab = sum(omega_z*dphi + d_omega_z*phi)


    !Build power matrices for fields, vels, and cross correls
    !All other quantities are background, so only power matrices dependent on
    !choice of \psi vs q ptb variable

    !NB: These are different for \psi_ij and q_ij variables
    !NB: Default is for \psi_ij, making the q_ij eqns a little more complicated
    !NB: \psi_ij = a*q_ij
    !NB: d\psi_ij = a*(q_ij + dq_ij)

    power_matrix=0e0_dp !in delta_phi
    d_power_matrix=0e0_dp !in psi
    cross_matrix=0e0_dp !in psi

    if (use_q) then
      !Don't divide out by scalefact
      prefactor= (k**3/2.0e0_dp/(pi**2))/(2e0_dp*k)
    else
      prefactor= (k**3/2.0e0_dp/(pi**2)/scalefactor**2)/(2e0_dp*k)
    end if

    do i=1,numb_infl; do j=1, numb_infl; do ll=1,numb_infl

      power_matrix(i,j) =power_matrix(i,j)+ &
        ptb_matrix(i,ll)*conjg(ptb_matrix(j,ll))*prefactor

      if (.not. use_q) then
        cross_matrix(i,j) =cross_matrix(i,j)+ &
          ptb_matrix(i,ll)*conjg(dptb_matrix(j,ll))*prefactor
      else
        cross_matrix(i,j) =cross_matrix(i,j)+ &
          ptb_matrix(i,ll)*conjg(ptb_matrix(j,ll) + dptb_matrix(j,ll))*prefactor
      end if

      if (.not. use_q) then
        d_power_matrix(i,j) =d_power_matrix(i,j)+&
          dptb_matrix(i,ll)*conjg(dptb_matrix(j,ll))*prefactor
      else
        d_power_matrix(i,j) =d_power_matrix(i,j)+&
          (ptb_matrix(i,ll)+dptb_matrix(i,ll))*&
          conjg(ptb_matrix(j,ll)+dptb_matrix(j,ll))*prefactor
      end if

    end do; end do; end do

    !------------------------------------------------------------
    !Adiabatic power spectrum
    power_adiab = dot_product(matmul(power_matrix,omega_z),omega_z)
    power_adiab = (1e0_dp/phi_dot_0_scaled**2)*power_adiab
    !------------------------------------------------------------

    !Isocurvature power spectra
    if (numb_infl>1) then

      !------------------------------------------------------------
      !Build the isocurv unit vects
      call build_isocurv_basis()

      !DEBUG
      !call build_isocurv_basisALTERNATIVE()

      !Vector of iso-spectra
      !Project power_matrix onto directs perpend to adiab direction
      pk_iso_vect = 0e0_dp
      do i=1,size(s_iso,1); do j=1, numb_infl; do ll=1,numb_infl
        pk_iso_vect(i) = pk_iso_vect(i) + &
          s_iso(i,j)*s_iso(i,ll)*power_matrix(j,ll)
      end do; end do; end do

      pk_iso_vect = pk_iso_vect*(1e0_dp/phi_dot_0_scaled**2)

      power_isocurv = 0e0_dp
      power_isocurv = sum(pk_iso_vect)

      !------------------------------------------------------------
      !Cross spectra between adiab and "projected" isocurvature
      power_cross = 0e0_dp
      do i =1, numb_infl; do j=1, numb_infl; do ll=1, size(s_iso,1)
        power_cross = power_cross + &
          omega_z(i)*s_iso(ll,j)*&
          (power_matrix(i,j) + power_matrix(j,i))
      end do; end do; end do
      power_cross = power_cross/(phi_dot_0_scaled**2)

      !------------------------------------------------------------

      !P(k) of total non-adiab pressure ptbs
      !dP_nad_i(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*dPsi_ij)*\hat{a}_j
      !<P_nad P_nad*>
      power_pnad =0e0_dp

      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)

      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp

      !DEBUG
      !call pnad_sumAA%clear()
      !call pnad_sumAB%clear()
      !call pnad_sumBA%clear()
      !call pnad_sumBB%clear()

      do i=1,numb_infl; do j=1,numb_infl

        AAprod = AAprod +A_vect(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect(j)*d_power_matrix(i,j)

        !DEBUG
        !call pnad_sumAA%add(real(A_vect(i)*A_vect(j)*power_matrix(i,j)))
        !call pnad_sumAB%add(real(A_vect(i)*B_vect(j)*cross_matrix(i,j)))
        !call pnad_sumBA%add(real(B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))))
        !call pnad_sumBB%add(real(B_vect(i)*B_vect(j)*d_power_matrix(i,j)))

      end do; end do
      power_pnad = (AAprod + BBprod) + (ABprod + BAprod)

      !Total pressure spectrum
      !<P P*>
      A_vect = get_A_vect_Ptotal(phi,dphi)
      B_vect = get_B_vect_Ptotal(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect(j)*d_power_matrix(i,j)
      end do; end do
      power_pressure = (AAprod + BBprod) + (ABprod + BAprod)

      !Total-adiab pressure cross-spectrum
      !<P_nad P_ad*>
      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect(i)*A_vect_adiab(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect_adiab(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect_adiab(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect_adiab(j)*d_power_matrix(i,j)
      end do; end do
      power_press_cross = (AAprod + BBprod) + (ABprod + BAprod)

      !<P_ad P_nad*>
      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect_adiab(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect_adiab(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect_adiab(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect_adiab(i)*B_vect(j)*d_power_matrix(i,j)
      end do; end do
      power_press_cross2 = (AAprod + BBprod) + (ABprod + BAprod)

      !<P_ad P_ad*>
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect_adiab(i)*A_vect_adiab(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect_adiab(i)*B_vect_adiab(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect_adiab(i)*A_vect_adiab(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect_adiab(i)*B_vect_adiab(j)*d_power_matrix(i,j)
      end do; end do
      power_press_adiab = (AAprod + BBprod) + (ABprod + BAprod)

      !The values (AA + BB) --> -(AB+BA) as approaches adiab limit.
      !Taking diff of "large" numbs means large error in the difference
      !Check if power_pnad is smaller than DP accuracy and set to zero

      prod_exponent = abs(real(AAprod)) + abs(real(BBprod)) +&
        abs(real(ABprod)) + abs(real(BAprod))

      if (power_pnad>1e-60_dp) then

        prod_exponent = log10(prod_exponent/power_pnad)

        if( prod_exponent >14e0_dp) then
          !Reached DP accuracy, all subseq results are numer error
          power_pnad = 0e0_dp
        end if

      else
        power_pnad=0e0_dp
      endif


      !Rescale so spectrum for S=(H/Pdot)dP_nad - total entropy ptb
      H=getH(phi,dphi)
      Pdot=getPdot(phi,dphi)
      power_entropy = ((H/Pdot)**2)*power_pnad

    else
      power_spectrum%isocurv =  0e0_dp
      power_spectrum%pnad    =  0e0_dp
      power_spectrum%entropy =  0e0_dp
      power_spectrum%cross_ad_iso =  0e0_dp
    end if

    power_spectrum%phi_ij  =  power_matrix
    power_spectrum%adiab   =  power_adiab
    power_spectrum%isocurv =  power_isocurv
    power_spectrum%pnad    =  power_pnad
    power_spectrum%entropy =  power_entropy
    power_spectrum%pressure =  power_pressure
    power_spectrum%press_ad =  power_press_adiab
    power_spectrum%cross_ad_iso =  power_cross


    contains

    subroutine build_isocurv_basisALTERNATIVE()

      !real(dp), dimension(num_inflaton), intent(in) :: omega_z
      real(dp), dimension(num_inflaton) :: Vprime
      integer :: i, j
      real(dp), dimension(num_inflaton,num_inflaton) :: perp
      real(dp) :: s_1_norm

      Vprime = dVdphi(phi)

      !Align first isocurv vector perpendicular to omega_z
      !NB: only works with trivial field-space metric
      !NB: need to upgrade delta_ij --> G_ij
      do i=1, num_inflaton; do j=1, num_inflaton
          if (i==j) then
            perp(i,j) = 1.0e0_dp - omega_z(i)*omega_z(j)
          else
            perp(i,j) = - omega_z(i)*omega_z(j)
          end if
      end do; end do

      s_1_norm = 0e0_dp
      s_iso = 0e0_dp
      do i=1,num_inflaton
        do j=1,num_inflaton
          s_iso(1,i) = s_iso(1,i) - &
            perp(i,j)*Vprime(j)

          s_1_norm = s_1_norm + &
            perp(i,j)*Vprime(i)*Vprime(j)
        end do
      end do
      s_iso(1,:) = s_iso(1,:)/sqrt(s_1_norm)

    end subroutine build_isocurv_basisALTERNATIVE


    !Basis perpend to adiab direction.  Gram-Schmidt
    subroutine build_isocurv_basis()

      real(dp), dimension(size(dphi),size(dphi)) :: spanning
      !field dirs
      real(dp), dimension(size(dphi),size(dphi)) :: phi_hat
      real(dp), dimension(size(dphi)-1,size(dphi)) :: field_basis
      real(dp), dimension(size(dphi)) :: phi_hat_temp
      real(dp) :: check_dot
      real(dp), dimension(size(dphi)-1) :: normalization, dnormalization
      integer :: i, j, adiab_index

      !Build field space vects
      do i=1,size(phi_hat,1); do j=1,size(phi_hat,2)
        if (i==j) then
          phi_hat(i,j) = 1e0_dp
        else
          phi_hat(i,j) = 0e0_dp
        end if
      end do; end do

      !Find field direction closest aligned to omega_z
      !This check makes sure get linearly indep field space vects
      check_dot=1e-30_dp
      adiab_index=0
      do i=1,size(phi_hat,1)
        if (abs(dot_product(phi_hat(i,:),omega_z(:)))>check_dot) then
          check_dot = abs(dot_product(phi_hat(i,:),omega_z(:)))
          adiab_index=i
        end if
      end do

      if (adiab_index ==0) then
        call raise%fatal_cosmo(&
        "It appears that no field space directions have projection &
        along the adiab direction.",&
        __FILE__, __LINE__)
      end if

      !Set the "most adiab" dir to 1
      if (adiab_index /= 1) then
        phi_hat_temp = phi_hat(1,:)
        phi_hat(1,:) = phi_hat(adiab_index,:)
        phi_hat(adiab_index,:) = phi_hat_temp
      end if

      spanning = 0e0_dp
      spanning(1,:) = omega_z/norm(omega_z)

      s_iso = 0e0_dp

      do i=2,size(spanning,1)

        spanning(i,:) = phi_hat(i,:)
        do j=1, i-1
          spanning(i,:) = spanning(i,:) - &
            projection(spanning(i,:),spanning(j,:))
        end do

        if (norm(spanning(i,:)) > div_zero_tol) then
          s_iso(i-1,:) = spanning(i,:)/norm(spanning(i,:))
        else
          print*, "MODPK: i=",i
          call raise%fatal_cosmo(&
            "The component spanning(i,:) has zero norm.",&
            __FILE__, __LINE__)
        end if

        !If there's a major hierarchy in scales for the vector components, then
        !there can be a spurious component of the isocurvature direction
        !oriented along the adiabatic one.  When you're approaching the
        !adiabatic limit this can be the dominant contribution, so we force the
        !smallest term to be zero whenever adding it doesn't affect the value of
        !the norm.
        !call renormalize_remove_smallest(s_iso(i-1,:))


        if (abs(dot_product(omega_z,s_iso(i-1,:)))>1e-12 .or.&
          isnan(abs(dot_product(omega_z,s_iso(i-1,:))))) then

          write(*,*), "omega_z.s_iso =",dot_product(omega_z,s_iso(i-1,:))," for i=",i-1
          call raise%fatal_cosmo(&
            "The isocurvature projection has a large adiabatic component.",&
            __FILE__, __LINE__)

        end if

      end do

    end subroutine build_isocurv_basis

    !Takes a vector and sets to zero those elements that are so small they're
    !not affecting the normalization

    !If there's a major hierarchy in scales for the vector components, then
    !there can be a spurious component of the isocurvature direction
    !oriented along the adiabatic one.  When you're approaching the
    !adiabatic limit this can be the dominant contribution, so we force the
    !smallest term to be zero whenever adding it doesn't affect the value of
    !the norm.
    subroutine renormalize_remove_smallest(vect)

      real(dp), dimension(:), intent(inout) :: vect
      real(dp), dimension(size(vect)) :: vect_temp
      real(dp) :: check
      integer :: smallest_index
      integer :: j

      check= 1e10_dp
      vect_temp = vect
      smallest_index= 1
      do j=1,size(vect)
        if (abs(vect(j)) < check) then
          check= abs(vect(j))
          smallest_index= j
        end if
      end do
      vect_temp(smallest_index) = 0e0_dp

      !1e-14 for double-precision round-off error
      if (abs(norm(vect)-norm(vect_temp))/abs(norm(vect)) < 1e-14) then
        vect = vect_temp/norm(vect_temp)
      end if

    end subroutine renormalize_remove_smallest

    !Functions for calculating the total non-adiabatic pressure perturbation
    !matrices

    !dP_nad(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*dPsi_ij)*\hat{a}_j
    function get_A_vect(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, firstterm, Vdot, c2, eps, gamm
      real(dp), dimension(size(phi)) :: secondterm, thirdterm, Vprime

      H=getH(phi,dphi)
      H2=H**2
      Vprime=dVdphi(phi)
      !c2 = getcs2(phi,dphi)
      Vdot = H*sum(Vprime*dphi)
      eps = geteps(phi,dphi)
      !gamm=(1.0e0_dp+c2)/(1.0e0_dp-c2)

      firstterm = (2.0e0_dp/3.0e0_dp/H2/sum(dphi*dphi))
      secondterm=(-3.0e0_dp*H2*sum(dphi*dphi) - (Vdot/H))*Vprime
      thirdterm = Vdot*H*dphi*(1.0e0_dp + 0.5e0_dp*sum(dphi*dphi))
      A=firstterm*(secondterm+thirdterm)

    end function get_A_vect


    function get_B_vect(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: firstterm, H, H2, Vdot, eps, c2
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      c2 = getcs2(phi,dphi)

      B = (1.0e0_dp-c2)*H2*dphi

    end function get_B_vect

    function get_A_vect_Ptotal(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, eps
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      eps = geteps(phi,dphi)
      Vprime =dVdphi(phi)

      A = -H2*dphi - eps*H2*dphi - Vprime

    end function get_A_vect_Ptotal

    function get_B_vect_Ptotal(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: H

      H=getH(phi,dphi)

      B = (H**2)*dphi

    end function get_B_vect_Ptotal

    function get_A_vect_Padiab(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, eps, cs2
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      eps = geteps(phi,dphi)
      Vprime =dVdphi(phi)
      cs2 = getcs2(phi,dphi)

      A = cs2*(-H2*dphi - eps*H2*dphi + Vprime)

    end function get_A_vect_Padiab

    function get_B_vect_Padiab(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: H, cs2

      H=getH(phi,dphi)
      cs2 = getcs2(phi,dphi)

      B = (H**2)*dphi*cs2

    end function get_B_vect_Padiab


    function getPdot(phi,dphi) result(Pdot)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: Pdot, H, H2
      real(dp), dimension(size(phi)) :: Vprime

      Vprime=dVdphi(phi)

      H=getH(phi,dphi)
      H2=H**2

      Pdot = -H*sum(3.0e0_dp*H2*dphi*dphi + 2e0_dp*Vprime*dphi)

    end function getPdot


    function getrhodot(phi,dphi) result(rhodot)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: rhodot, H, H2

      H=getH(phi,dphi)

      rhodot = -3.0e0_dp*H**3*sum(dphi*dphi)

    end function getrhodot


    function getcs2(phi,dphi) result(cs2)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: cs2

      cs2 = getpdot(phi,dphi)/getrhodot(phi,dphi)

    end function getcs2


  end subroutine powerspectrum


  !Projection of v orthogonally onto line spanned by u
  pure function projection(v,u) result(proj)

    real(dp), dimension(:), intent(in) :: v, u
    real(dp), dimension(size(v)) :: proj

    proj = u*(dot_product(v,u)/dot_product(u,u))

  end function projection


  !Norm of v
  pure function norm(v)

    real(dp), dimension(:), intent(in) :: v
    real(dp) :: norm

    norm = sqrt(dot_product(v,v))

  end function norm


    !END MULTIFIELD


  pure FUNCTION tensorpower(v, a)
    USE internals
    real(dp) :: tensorpower
    real(dp), INTENT(IN) :: a
    COMPLEX(KIND=DP), intent(in) :: v

    !MULTIFIELD: Calculates P_h(k) given v, a
    tensorpower = abs(v)**2/(2.0e0_dp*k) / a**2 * (k**3)*4.0e0_dp/(PI**2)/(M_Pl**2)
    !END MULTIFIELD

  END FUNCTION tensorpower

  pure FUNCTION zpower(u_zeta, dsigma, a)
    USE internals
    real(dp) :: zpower
    real(dp), INTENT(IN) :: dsigma
    real(dp), INTENT(IN) :: a
    COMPLEX(KIND=DP), INTENT(IN) :: u_zeta

    zpower = abs(u_zeta**2)/dsigma**2/a**2 /(2.0e0_dp*k) * (k**3)/(2.0e0_dp*PI**2)

  END FUNCTION zpower


  function trace_d2logVdphi2(phi) result(trace)
    !
    !     Returns trace of d^2V/dPhi^2 given phi
    !     Used in calculating bundle exp_scalar
    !

    real(dp), intent(in) :: phi(:)
    real(dp) :: trace, V_ab(size(phi),size(phi)), V_a(size(phi)), V
    integer :: i

    V_ab = d2Vdphi2(phi)
    V_a = dVdphi(phi)
    V = pot(phi)

    trace = 0e0_dp
    do i=1,size(phi)
      trace = trace + V_ab(i,i)/V - (V_a(i)/V)**2
    end do

  end function trace_d2logVdphi2

  !Calc bundle exp_scalar by integrating tr(d2Vdphi2) to efold by Riemann sum
  subroutine bundle_exp_scalar(self,phi, efold)

    class(bundle) :: self
    real(dp), intent(in) :: phi(:), efold
    real(dp) :: dN

    dN = efold - self%N

    self%dlogThetadN= self%dlogThetadN - &
      dN*trace_d2logVdphi2(phi)

    self%exp_scalar=exp(self%dlogThetadN)

    self%N=efold

  end subroutine bundle_exp_scalar


  pure FUNCTION MySech(x)
    real(dp), intent(in)  :: x
    real(dp) :: MySech

    IF(ABS(x).GT.40.0e0_dp) THEN
       MySech=0.0e0_dp
    ELSE
       MySech=1.0e0_dp/COSH(x)
    END IF
  END FUNCTION MySech

  !Checks to see if H^2=V/(3-eps) is stable, ie, if H>0 and if V~0, then it will
  !return FALSE for unstable.
  subroutine stability_check_on_H(stable,phi,dphi,using_t)
    logical, intent(inout) :: stable
    real(dp), dimension(:), intent(in) :: phi, dphi
    logical, intent(in) :: using_t
    real(dp) :: eps, V

    !DEBUG
    !overriding stability check...
    !stable = .true.
    !return

    V=pot(phi)

    if (using_t) then
      eps = getEps_with_t(phi, dphi)
    else
      eps = getEps(phi, dphi)
    end if

    stable = (3.0e0_dp -eps >1.0e-6)


  end subroutine stability_check_on_H

  !From a histogram-estimate of a PDF, return the log(P) for a
  !vector-valued observable, given the model parameters, which are
  !implicit in the construction of the histogram
  function logP_of_observ(observ, hist) result(logP)
    implicit none

    real(dp), dimension(:,:), intent(in) :: hist
    real(dp), dimension(:), intent(in) :: observ
    real(dp) :: logP

    real(dp), parameter :: logzero=-1e30_dp
    real(dp), parameter :: tol_pdf = 1e-3_dp
    integer(dp) :: ndimns

    real(dp), dimension(size(hist,2)-1) :: binsize
    integer :: bin_position
    integer :: ii, jj

    ndimns=size(hist,2)-1

    !Find binsize from histogram
    do jj=1, ndimns
      do ii=2,size(hist,1)
        !Hist bins organized in increasing, repeated "chunks"
        if (hist(ii,jj) /= hist(1,jj)) then
          binsize(jj)=hist(ii,jj)-hist(1,jj)
          exit
        end if
      end do
    end do

    !Find which bin observ is in
    !Could probably be improved with locate alg...
    bin_position=0
    do ii=1,size(hist,1)
      if (all(observ .ge. hist(ii,1:ndimns)) .and.&
        all(observ .le. hist(ii,1:ndimns)+binsize)) then
        bin_position = ii
        exit
      end if
    end do

    if (bin_position==0) then
      !If not in any bin, then set zero probability
      logP = logzero
      !DEBUG
      print*, "setting logP=logzero A"
      print*, "observ = ", observ
      stop
    else if (log(hist(bin_position,ndimns+1))<logzero) then
      logP = logzero
      !DEBUG
      print*, "setting logP=logzero B"
      print*, "observ = ", observ
      stop
    else
      logP = log( hist(bin_position,ndimns+1))
    end if

  end function logP_of_observ

  !Given parameters, returns the Fisher-Rao metric at that point in the
  !parameter manifold.  Distance measured by F-R metric is the Fisher
  !information
  !g_mu_nu = <dlogP/dparam_mu dlogP/dparam_nu>
  function fisher_rao_metric(param, stepsize) result(g_mu_nu)
    implicit none

    real(dp), dimension(:), intent(in) :: param
    real(dp), dimension(:), intent(in) :: stepsize
    real(dp), dimension(size(param),size(param)) :: g_mu_nu

    real(dp), dimension(size(param)) :: observ

    real(dp), dimension(size(param)) :: dlogP
    real(dp), dimension(size(param)) :: int_limit_min, int_limit_max

    real(dp), dimension(size(param)) :: param_plus, param_minus
    real(dp), dimension(size(param)) :: param_plus2, param_minus2
    real(dp), dimension(size(param)) :: param_plus3, param_minus3

    integer(dp) :: ii, jj, kk
    integer(dp) :: ndimns, nbins_int(size(param))
    real(dp), dimension(size(param)) :: binsize

    real(dp), dimension(:,:), allocatable :: int_grid

    real(dp), dimension(:,:), allocatable :: dataset
    real(dp), dimension(:,:), allocatable :: hist_param

    type :: hist_array
      real(dp), dimension(:,:), allocatable :: hist
    end type
    type(hist_array), dimension(size(param)) :: hist_plus, hist_minus
    type(hist_array), dimension(size(param)) :: hist_plus2, hist_minus2
    type(hist_array), dimension(size(param)) :: hist_plus3, hist_minus3

    real(dp) :: dV, P_o, dlogP_i, dlogP_j
    logical :: quadrature

    ndimns = size(param)

    call obtain_derivs()

    !Get run data for param
    call inflation_sample(param, dataset)
    !Build histograms from data
    hist_param = histogram_Nd(dataset, method=2, norm=1)


    !Average over observables for
    !g_mu_nu = <dlogP/dparam_mu dlogP/dparam_nu>

    quadrature = .false.
    if (quadrature) then

      call int_dlogP_quadrature()

    else
      !Integrate using the dataset from inflation_sample(param),
      !since we have assumed this is a representative sample in order
      !to build PDFs

      call int_dlogP_MCMC()

    end if


    contains

      subroutine int_limits()

        !Find the integration limits
        do ii=1,ndimns
          int_limit_min(ii) = minval(hist_param(:,ii))
          int_limit_max(ii) = maxval(hist_param(:,ii))

          do jj=1,size(hist_plus)
            int_limit_min(ii) = min(int_limit_min(ii) ,&
              minval(hist_plus(jj)%hist(:,ii)), &
              minval(hist_minus(jj)%hist(:,ii)))

            int_limit_max(ii) = max(int_limit_max(ii) ,&
              maxval(hist_plus(jj)%hist(:,ii)), &
              maxval(hist_minus(jj)%hist(:,ii)))
          end do
        end do

      end subroutine int_limits


      subroutine obtain_derivs()

        !Calc derivatives dlogP/dparam_mu
        do ii=1,size(param)
          !Get run data for param_plus and param_minus
          !Build histograms from data
          param_plus = param
          param_plus(ii) = param_plus(ii) + stepsize(ii)
          param_minus = param
          param_minus(ii) = param_minus(ii) - stepsize(ii)

          param_plus2 = param
          param_plus2(ii) = param_plus2(ii) + 2.0e0_dp*stepsize(ii)
          param_minus2 = param
          param_minus2(ii) = param_minus2(ii) -2.0e0_dp* stepsize(ii)

          param_plus3 = param
          param_plus3(ii) = param_plus3(ii) + 3.0e0_dp*stepsize(ii)
          param_minus3 = param
          param_minus3(ii) = param_minus3(ii) -3.0e0_dp* stepsize(ii)

          call inflation_sample(param_plus, dataset)
          hist_plus(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

          call inflation_sample(param_minus, dataset)

          hist_minus(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

          call inflation_sample(param_plus2, dataset)
          hist_plus2(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

          call inflation_sample(param_minus2, dataset)
          hist_minus2(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

          call inflation_sample(param_plus3, dataset)
          hist_plus3(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

          call inflation_sample(param_minus3, dataset)
          hist_minus3(ii)%hist = histogram_Nd(dataset, method=2, norm=1)
          deallocate(dataset)

        end do
      end subroutine obtain_derivs

      subroutine int_dlogP_quadrature()

        !NB: Currently only integrates in quadrature.
        !For higher dimensions, implement some other technique
        !like Metropolis-Hastings, etc.

        !Set the integration limits
        call int_limits()

        !Estimate how many integration bins are relevant
        !Could be improved drastically...
        do jj=1,ndimns
          do ii=2,size(hist_param)
            if (hist_param(ii,jj)>hist_param(1,jj)) then
              binsize(jj) = (maxval(hist_param(:,jj)) - minval(hist_param(:,jj)))
              nbins_int(jj) = &
                binsize(jj) /(hist_param(ii,jj)-hist_param(1,jj))
              exit
            end if
          end do
        end do
        !DEBUG
        nbins_int = 20.0*nbins_int

        !Switch to observable space binsize
        binsize = (int_limit_max-int_limit_min)/nbins_int

        !Make an integration grid
        allocate(int_grid(product(nbins_int),ndimns))
        call make_grid(nbins_int, int_limit_max, int_limit_min, int_grid)

        !Perform the integration
        g_mu_nu =0e0_dp
        do ii=1,size(g_mu_nu,1); do jj=1,size(g_mu_nu,2)
          if (jj<ii) then
            g_mu_nu(ii,jj) = g_mu_nu(jj,ii)
            cycle
          end if
          do kk=1,size(int_grid,1)

            observ = int_grid(kk,:)+binsize/2.0

            dV=product(binsize)
            P_o = exp(logP_of_observ(observ,hist_param))
            dlogP_i = (0.5e0_dp/stepsize(ii))*&
              (logP_of_observ(observ,hist_plus(ii)%hist) &
                - logP_of_observ(observ,hist_minus(ii)%hist))
            dlogP_j = (0.5e0_dp/stepsize(jj))*&
              (logP_of_observ(observ,hist_plus(jj)%hist) &
                - logP_of_observ(observ,hist_minus(jj)%hist))

            g_mu_nu(ii,jj) = g_mu_nu(ii,jj) + &
              dV * P_o * dlogP_i * dlogP_j

          end do

        end do; end do

      end subroutine int_dlogP_quadrature

      subroutine int_dlogP_MCMC()

        print*, "logP", logP_of_observ(observ,hist_param)

        print*, "logP", logP_of_observ(observ,hist_plus(1)%hist)

        print*, "dlogP O(2)", (1.0/stepsize(1))*(logP_of_observ(observ,hist_plus(1)%hist) - logP_of_observ(observ,hist_param))
        print*, "dlogP O(4)", (1.0/stepsize(1))*(&
          (1.0e0_dp/12.0e0_dp)*logP_of_observ(observ,hist_minus2(1)%hist)&
          -(2.0e0_dp/3.0e0_dp)*logP_of_observ(observ,hist_minus(1)%hist)&
          +(2.0e0_dp/3.0e0_dp)*logP_of_observ(observ,hist_plus(1)%hist)&
          -(1.0e0_dp/12.0e0_dp)*logP_of_observ(observ,hist_plus2(1)%hist))
        stop

        g_mu_nu = 0e0_dp
        do ii=1,size(g_mu_nu,1)
          do jj=1,size(g_mu_nu,2)
            do kk=1, size(dataset,1)

              observ = dataset(kk,:)

              !Get the derivs at this observ
                  dlogP_i = (0.5e0_dp/stepsize(ii))*&
                    (logP_of_observ(observ,hist_plus(ii)%hist) &
                      - logP_of_observ(observ,hist_minus(ii)%hist))
                  dlogP_j = (0.5e0_dp/stepsize(jj))*&
                    (logP_of_observ(observ,hist_plus(jj)%hist) &
                      - logP_of_observ(observ,hist_minus(jj)%hist))

              !Integrate
              g_mu_nu(ii,jj) = g_mu_nu(ii,jj) + &
                dlogP_i*dlogP_j

            end do
          end do
        end do

        !Rescale
        g_mu_nu = g_mu_nu/real(size(dataset,1),kind=dp)


        end subroutine int_dlogP_MCMC

  end function fisher_rao_metric

  !Given a parametrization of the inflationary model, obtain an
  !MCMC sample of points in observable space from which we can
  !build PDFs
  subroutine inflation_sample(param, dataset)
    use modpk_rng, only : init_random_seed_serial, normal_scalar
    implicit none

    real(dp), dimension(:), intent(in) :: param
    real(dp), dimension(:,:), allocatable, intent(out) :: dataset
    integer :: i, j

    print*, "std=", param(1)


!DEBUG
allocate(dataset(1000,1))
call init_random_seed_serial()
do i=1,size(dataset,1)
  do j=1,size(dataset,2)
    dataset(i,j) = normal_scalar(0.e0_dp,abs(param(1)))
end do; end do
dataset=dataset


  end subroutine inflation_sample

  !Guess the end of inflation point of phi based off of the value of phi at
  !the point where the pivot scale leaves the horizon and the heaviest field's
  !EOI value.  Based off N-1 constants of motion in the SR limit.
  function guess_EOI_field(phi_horizon, phi_end) result(guess)

    real(dp), dimension(:) :: phi_horizon, phi_end
    integer, dimension(1) :: heavy_index
    real(dp) :: m2_V(size(phi_horizon)), p_exp
    real(dp), dimension(size(phi_horizon)) :: guess
    real(dp) :: phi_end_heavy, phi_piv_heavy
    real(dp), dimension(1) :: temp_end, temp_piv, temp_hvy
    integer :: ii
    real(dp) :: expon, m_hvy

    select case(potential_choice)

    case(16)
      ! (1/p) lambda_i |phi_i|^p --- N-monomial

      p_exp = vparams(2,1)
      m2_V = vparams(1,:)

      !DEBUG
      !phi_end = 0e0_dp


      !Find heaviest field index
      heavy_index = maxloc(m2_V)
      temp_end = phi_end(heavy_index)
      temp_piv = phi_horizon(heavy_index)
      temp_hvy = m2_V(heavy_index)
      phi_end_heavy = abs(temp_end(1))
      phi_piv_heavy = abs(temp_piv(1))
      m_hvy = temp_hvy(1)

      print*, p_exp

      !Find guess
      if (abs(p_exp-2.0e0_dp) < 1.0e-6_dp) then
        guess = 0e0_dp
        do ii=1, size(phi_horizon)
          guess(ii) = abs(phi_horizon(ii))*(phi_end_heavy/phi_piv_heavy)&
            **(m2_V(ii)/m_hvy)
        end do

      else

        expon = 2.0e0_dp - p_exp
        do ii=1, size(phi_horizon)
          guess(ii) = (abs(phi_horizon(ii))**expon - &
            (m2_V(ii)/m_hvy)*(phi_piv_heavy**expon + phi_end_heavy**expon)) &
            **(1.0e0_dp/expon)
        end do

      end if

      !print*, "guess:"
      !print*, guess

      !print*, "real:"
      !print*, phi_end


      print*, "ratio:"
      do ii=1, size(phi_end)
        write(400,'(4e12.3)'), m2_V(ii)/m_hvy, abs(phi_end(ii)), abs(guess(ii)), abs(phi_horizon(ii))
        write(*,'(4e12.3)'), m2_V(ii)/m_hvy, abs(phi_end(ii)), abs(guess(ii)), abs(phi_horizon(ii))
      end do

      stop

    case default
      print*, "MODPK: potential choice =", potential_choice, &
        "not implemented in guess_EOI_field"
      stop
    end select

  end function guess_EOI_field




END MODULE potential


!Module that lets us compare the mode evolution to the slow-roll expectation for
!sum-separable potentials, using the results of Battefeld-Easther astro-ph/0610296
module modpk_deltaN_SR
  use modpkparams, only : dp, vparams
  use modpk_observables, only : power_spectra
  use internals, only : pi
  use potential, only : pot, dVdphi, d2Vdphi2, d3Vdphi3
  implicit none

  !The horizon crossing approximation
  logical :: HC_approx=.false.

  contains

    !Evalutes adiabatic power spectrum under SR approximation at field positions
    !phi_end, given that the mode of interest crossed the horizon at phi_pivot
    !Assumes a massless, uncorrelated mode subhorizon
    function PR_SR(phi_pivot,phi_end, spectrum) result(PR)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) ::PR
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp) :: H_piv, V_piv, P_dphi
      real(dp), dimension(size(phi_pivot)) :: eps_i, u_i
      type(power_spectra), intent(in), optional :: spectrum
      integer :: ii, jj

      V_piv = pot(phi_pivot)
      H_piv = sqrt(V_piv/3.0e0_dp)
      dN = dNdphi_SR(phi_pivot,phi_end)

      eps_i = eps_SR(phi_pivot)
      u_i = (V_i_sum_sep(phi_pivot)+Z_i_BE(phi_end))/V_piv

      P_dphi = (H_piv/2.0e0_dp/pi)**2

      PR = sum(dN*dN)*P_dphi

      if (present(spectrum)) then
        PR=0e0_dp
        do ii=1,size(dN); do jj=1,size(dN)
          PR = PR+ &
            dN(ii)*dN(jj)*spectrum%phi_ij(ii,jj)
        end do; end do
      end if

    end function PR_SR

    function r_SR(phi_pivot,phi_end, spectrum) result(r)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: r
      real(dp), dimension(size(phi_pivot)) :: dN
      type(power_spectra), intent(in), optional :: spectrum
      real(dp) :: P_tens, P_scal, H, V

      integer :: ii, jj

      dN = dNdphi_SR(phi_pivot,phi_end)

      r = 8.0e0_dp/sum(dN*dN)

      if (present(spectrum)) then
        V = pot(phi_pivot)
        H = sqrt(V/3.0e0_dp)
        P_tens = 8.0e0_dp*(H/2.0e0_dp/pi)**2

        P_scal=0e0_dp
        do ii=1,size(dN); do jj=1,size(dN)
          P_scal = P_scal + &
            dN(ii)*dN(jj)*spectrum%phi_ij(ii,jj)
        end do; end do

        r = P_tens/P_scal

      end if


    end function r_SR

    function nt_SR(phi_pivot) result(nt)
      real(dp), dimension(:), intent(in) :: phi_pivot
      real(dp) :: nt, eps(size(phi_pivot))

      eps = eps_SR(phi_pivot)

      nt = -2.0e0_dp*sum(eps)
      !nt = -2.0e0_dp*sum(eps)/(1.0e0_dp - sum(eps))

    end function nt_SR

    !Assumes independent, GRFs at horizon exit
    function fNL_SR(phi_pivot,phi_end) result(fnl)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: fnl
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2N
      integer :: ii, jj

      dN = dNdphi_SR(phi_pivot,phi_end)
      d2N = d2Ndphi2_SR(phi_pivot,phi_end)

      fnl=0e0_dp
      do ii=1,size(dN); do jj=1,size(dN)
        fnl = fnl + dN(ii)*dN(jj)*d2N(ii,jj)
      end do; end do
      fnl = fnl*(-5.0e0_dp/6.0e0_dp)/(sum(dN*dN))**2

    end function fNL_SR

    !Assumes independent, GRFs at horizon exit
    !From Eq 41 in astro-ph/0611075
    function tauNL_SR(phi_pivot,phi_end) result(taunl)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: taunl
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2N
      integer :: aa, bb, cc

      dN = dNdphi_SR(phi_pivot,phi_end)
      d2N = d2Ndphi2_SR(phi_pivot,phi_end)

      taunl=0e0_dp
      do aa=1,size(dN); do bb=1,size(dN); do cc=1,size(dN)
        taunl = taunl + d2N(aa,bb)*d2N(aa,cc)*dN(bb)*dN(cc)
      end do; end do; end do
      taunl = taunl*(1.0e0_dp)/(sum(dN*dN))**3

    end function tauNL_SR

    function ns_SR(phi_pivot,phi_end) result(ns)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: ns, eps_piv, V
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2V
      real(dp), dimension(size(phi_pivot)) :: dV, dN
      integer :: ii, jj
      real(dp), dimension(size(phi_pivot)) :: eps_i, eta_i, u_i

      !dV = dVdphi(phi_pivot)
      d2V = d2Vdphi2(phi_pivot)
      dN = dNdphi_SR(phi_pivot,phi_end)
      !eps_i = eps_SR(phi_pivot)
      !eta_i = eta_SR(phi_pivot)
      eps_piv = sum(eps_SR(phi_pivot))
      V = pot(phi_pivot)

      !u_i = (V_i_sum_sep(phi_pivot)+Z_i_BE(phi_end))/V

      !ns = 1.0e0_dp - 2.0e0_dp*eps_piv &
      !  - (4.0e0_dp/sum(u_i**2/eps_i))*&
      !  (1.0e0_dp - sum(eta_i*u_i**2/2.0e0_dp/eps_i))

      ns = 1.0e0_dp &
        - 2.0e0_dp*eps_piv &
        - (2.0e0_dp/sum(dN*dN))

      do ii=1,size(phi_pivot); do jj=1,size(phi_pivot)
        ns = ns +&
          2.0e0_dp*d2V(ii,jj)*dN(ii)*dN(jj)/V/sum(dN*dN)
      end do; end do

    end function ns_SR

    !Running of ns
    !Formula as in Eq 6.14 1203.3792
    function alpha_s_SR(phi_pivot,phi_end) result(alpha_s)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: alpha_s


      real(dp), dimension(size(phi_pivot)) :: V_i, Z_i, u_i
      real(dp) :: V, eps_piv
      real(dp), dimension(size(phi_pivot)) :: eps_i_piv, eta_i_piv, &
        xi_i_piv
      real(dp) :: sum_ui_over_epsi
      real(dp) :: term1, term2, term3, term4, term5, term6

      V = pot(phi_pivot)
      eps_i_piv = eps_SR(phi_pivot)
      eps_piv = sum(eps_i_piv)
      xi_i_piv = xi_SR(phi_pivot)
      eta_i_piv = eta_SR(phi_pivot)

      V_i = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      u_i = (V_i + Z_i)/V

      sum_ui_over_epsi = sum(u_i**2/eps_i_piv)

      term1 = -8.0e0_dp * eps_piv**2
      term2 = 4.0e0_dp * sum(eps_i_piv*eta_i_piv)

      term3 = (-16.0e0_dp/sum_ui_over_epsi**2)*&
        (1.0e0_dp - sum(eta_i_piv*u_i**2/2.0e0_dp/eps_i_piv))**2

      term4 = (-8.0e0_dp/sum_ui_over_epsi)*&
        sum(eta_i_piv*u_i*(1.0e0_dp - &
          eta_i_piv*u_i**2/2.0e0_dp/eps_i_piv))

      term5 = (4.0e0_dp*eps_piv/sum_ui_over_epsi)*&
        sum(eta_i_piv*u_i**2/eps_i_piv)

      term6 = (-2.0e0_dp/sum_ui_over_epsi)*&
        sum(xi_i_piv*u_i**2/eps_i_piv)

      alpha_s = term1 + term2 + term3 + term4 + term5 + term6


    end function alpha_s_SR

    !Deriv of N wrt phi on horiz cross surface
    !Battefeld-Easther eq 29
    function dNdphi_SR(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_pivot)) :: dNdphi_SR
      real(dp), dimension(size(phi_pivot)) :: eps, V_i, Z_i
      real(dp) :: V

      eps = eps_SR(phi_pivot)
      V_i = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      V = pot(phi_pivot)

      dNdphi_SR = ((1.0e0_dp/sqrt(2.0e0_dp*eps))/V)*(V_i + Z_i)

    end function dNdphi_SR

    !2nd Deriv of N wrt phi on horiz cross surface
    !Battefeld-Easther eq 32
    function d2Ndphi2_SR(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2Ndphi2_SR
      real(dp), dimension(size(phi_end),size(phi_end)) :: delta

      real(dp), dimension(size(phi_end)) :: eta_piv, eps_piv, V_i_piv, Z_i
      real(dp) :: V_piv
      real(dp), dimension(size(phi_end),size(phi_end)) :: dZ_ij

      integer :: ll, kk

      eta_piv = eta_SR(phi_pivot)
      eps_piv = eps_SR(phi_pivot)
      V_i_piv = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      V_piv = pot(phi_pivot)
      dZ_ij = dZdphi_ij_BE(phi_pivot,phi_end)

      !Identity matrix
      delta=0e0_dp
      do ll=1, size(phi_end)
        delta(ll,ll) = 1.0e0_dp
      end do

      do ll=1, size(phi_end); do kk=1,size(phi_end)
        d2Ndphi2_SR(ll,kk) = delta(kk,ll)*&
            (1.0e0_dp - &
              (eta_piv(ll)/2.0e0_dp/eps_piv(ll))*&
              (V_i_piv(ll)+Z_i(ll))/V_piv&
             ) +&
             (1.0e0_dp/sqrt(2.0e0_dp*eps_piv(ll))/V_piv)*dZ_ij(ll,kk)
      end do; end do

    end function d2Ndphi2_SR

    function eps_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: eps_SR
      real(dp), dimension(size(phi)) :: dV
      real(dp) :: V

      dV = dVdphi(phi)
      V = pot(phi)
      eps_SR = 0.5e0_dp*( dV**2/V**2)

    end function eps_SR

    function eta_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: eta_SR
      real(dp), dimension(size(phi),size(phi)) :: d2V
      real(dp) :: V
      integer :: i

      d2V = d2Vdphi2(phi)
      V = pot(phi)

      do i=1, size(eta_SR)
        eta_SR(i) = d2V(i,i)**2/V
      end do

    end function eta_SR

    function xi_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: xi_SR
      real(dp), dimension(size(phi),size(phi),size(phi)) :: d3V
      real(dp) :: V
      real(dp), dimension(size(phi)) :: dV
      integer :: i

      d3V = d3Vdphi3(phi)
      V = pot(phi)
      dV = dVdphi(phi)

      do i=1, size(xi_SR)
        xi_SR(i) = dV(i)*d3V(i,i,i)**2/V**2
      end do

    end function xi_SR

    !The function Z_i from Battefeld-Easther that encodes all details from the
    !end of inflation surface.  Eq. 31
    function Z_i_BE(phi_end)
      real(dp), dimension(:), intent(in) :: phi_end
      real(dp), dimension(size(phi_end)) :: Z_i_BE

      real(dp), dimension(size(phi_end)) :: eps_end
      real(dp) :: eps, V

      if (HC_approx) then
        Z_i_BE=0e0_dp
        return
      end if

      eps_end = eps_SR(phi_end)
      eps = sum(eps_end)
      V = pot(phi_end)

      Z_i_BE = V*eps_end/eps - V_i_sum_sep(phi_end)

    end function Z_i_BE

    !Deriv of Z_i wrt fields at horiz cross
    !Battefeld-Easther Eq. 33
    function dZdphi_ij_BE(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_end),size(phi_end)) :: dZdphi_ij_BE

      real(dp), dimension(size(phi_end)) :: eps_end, eps_piv, eta_end
      real(dp) :: eps_t_end, V_end, V_piv
      real(dp), dimension(size(phi_end),size(phi_end)) :: delta

      integer :: jj, ll, kk

      if (HC_approx) then
        dZdphi_ij_BE=0e0_dp
        return
      end if

      eps_end = eps_SR(phi_end)
      eps_piv = eps_SR(phi_pivot)
      eps_t_end = sum(eps_end)
      V_end = pot(phi_end)
      V_piv = pot(phi_pivot)
      eta_end = eta_SR(phi_end)

      !Identity matrix
      delta=0e0_dp
      do ll=1, size(phi_end)
        delta(ll,ll) = 1.0e0_dp
      end do

      !Summation over jj
      dZdphi_ij_BE = 0e0_dp
      do ll=1, size(phi_end); do kk=1,size(phi_end)
        do jj=1,size(phi_end)
          dZdphi_ij_BE(ll,kk) = dZdphi_ij_BE(ll,kk) + &
            (-V_end**2/V_piv)*sqrt(2.0e0_dp/eps_piv(kk))*&
            (&
              eps_end(jj)*&
              ((eps_end(ll)/eps_t_end) - delta(ll,jj))*&
              ((eps_end(kk)/eps_t_end) - delta(kk,jj))*&
              (1.0e0_dp - (eta_end(jj)/eps_t_end))&
            )
        end do
      end do; end do


    end function dZdphi_ij_BE

    !For a sum-separable potential V=\sum_i V_i.  This returns only the V_i part
    function V_i_sum_sep(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: V_i_sum_sep

      real(dp), dimension(:,:), allocatable :: vparams_temp
      integer :: vrows, jj
      real(dp) :: V

      !The idea: make temp copy of vparams; change vparams as if it had only the
      !one field; get V; restore vparams
      !NB: vparams(vrows,num_inflaton)

      vrows = size(vparams,1)

      allocate(vparams_temp(size(vparams,1),size(vparams,2)))
      vparams_temp = vparams

      do jj=1,size(phi)
        deallocate(vparams)
        allocate(vparams(vrows,1))
        vparams(:,1) = vparams_temp(:,jj)
        V_i_sum_sep(jj) = pot((/phi(jj)/))
      end do

      deallocate(vparams)
      allocate(vparams(size(vparams_temp,1),size(vparams_temp,2)))
      vparams = vparams_temp
      deallocate(vparams_temp)

    end function V_i_sum_sep



end module modpk_deltaN_SR
