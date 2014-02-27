MODULE potential
  USE modpkparams
  use internals, only : pi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot, getH, getdHdalpha, getEps, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, geteta, zpower, getH_with_t, stability_check_on_H, getEps_with_t,&
       effective_V_choice, turning_choice

  public :: norm
  public :: bundle, field_bundle

  type :: bundle
    real(dp) :: N=0e0_dp
    real(dp) :: exp_scalar=0e0_dp
    real(dp) :: dlogThetadN=0e0_dp
    contains
      procedure :: calc_exp_scalar =>bundle_exp_scalar
  end type bundle

  type(bundle) :: field_bundle

  integer :: effective_V_choice
  !For turning trajs, set choice of turning function here.
  !Need one choice for every heavy direction (assumes one light direction)
  !Ref functions turning_function and derivs
  integer, dimension(:), allocatable :: turning_choice


CONTAINS


  recursive function pot(phi) result(V_potential)
    !
    !     Returns V(phi) given phi, phi can be an array for multifield, 
    !     The code implement multifield potential in the form of V = \sum V(phi_i), 
    !     More complicated form of potentials can be customized 
    !
    !real(dp) :: pot
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
    real(dp) :: lambda4(num_inflaton), alpha2(num_inflaton)

    real(dp) :: lambda2
    real(dp), dimension(num_inflaton,num_inflaton) :: m2_matrix
    integer :: i,j, temp_choice

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
         print*, "Potential_choice", Potential_choice, "requires two fields."
         print*, "Number of fields =", size(phi)
         stop
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
      ! V_potential = 0.5*M2**2*cos(0.5*theta2)*(sin(0.5*theta2)*(phi(1)-c2)+&
      !       cos(0.5*theta2)phi(2)+tan((1/pi)*theta2*atan(s2*(cos(0.5*theta2)*(phi(1)-c2)-&
      !       sin(0.5*theta2)phi(2))))(-cos(0.5*theta2)*(phi(1)-c2)+sin(0.5*theta2)phi(2)))**2
      theta2 = Pi/10.0
      c2 = 0.0
      M2 = 2.0E-4
      s2 = 100.0*sqrt(3.0)
      mphi1 = 1E-7
      phi1shift = 100.0

      potlarge = (mphi1**2*(-phi1shift + phi(1))**2)/2.
      potsmall = (M2**2*Cos(theta2/2.)*(Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) +&
                 (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
                 Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2)/2.
      V_potential=potlarge+potsmall
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
      V_potential = pot(phi)
      potential_choice = temp_choice

      !Choice of function to use set in parameters_multimodecode.txt
#define PHI_I phi(1), i, turning_choice(i-1)
#define DELTAPHI (phi(i) - turning_function(PHI_I))
#define EXPTERM exp( (-0.5e0_dp/alpha2(i))*(phi(i) - turning_function(PHI_I))**2)

      !Check to see if too far outside valley in any massive field direction
      !do i=2,num_inflaton
      !  if (abs(DELTAPHI**3/DELTAPHI**2)<0.5e0_dp) then
      !    print*, "ERROR: Trajectory too far outside the valley for"
      !    print*, "potential_choice=",potential_choice
      !    print*, "to be applicable."
      !    stop
      !  end if
      !end do

      !print*, "philight", phi(1)

      do i=2, num_inflaton
        V_potential = V_potential - lambda4(i)*&
          (EXPTERM  - 1.0e0_dp)
      end do

    case default
       write(*,*) 'MODPK: Need to set pot(phi) in modpk_potential.f90 for potential_choice =',potential_choice
       STOP
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
    real(dp) :: lambda4(num_inflaton), alpha2(num_inflaton)

    real(dp) :: lambda2
    real(dp), dimension(num_inflaton,num_inflaton) :: m2_matrix

    if (vnderivs) then
       ! MULTIFIELD
       do i=1, num_inflaton
          phiplus = phi
          phiplus(i) = phi(i) + 0.5e0_dp*phi(i)*findiffdphi**(1.0e0_dp/3.0e0_dp)
          dphi = phiplus - phi
          first_deriv(i) = (pot(phi+dphi)-pot(phi-dphi))/(2.0e0_dp*dphi(i))
          if (first_deriv(i).eq.0.e0_dp) then
             write(*,*) 'MODPK: For i=', i
             write(*,*) 'MODPK: first_deriv(i)=0, possibly signaling a problem with accuracy of numerical derivatives.'
             write(*,*) 'MODPK: Try using vnderivs=F if possible.'
             STOP
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

          first_deriv(1) = (lambda_hybrid**4)*(4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)*phi(1)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)
          first_deriv(2) = (lambda_hybrid**4)*(2.0_dp*phi(2)/mu_hybrid**2 +&
            2.0_dp*phi(1)**2*phi(2)/nu_hybrid**4)
       case(9)!Saddle type things
          m2_V = (vparams(1,:))
          first_deriv = m2_V*phi
       case(10)!Modified Langlois model
      theta2 = Pi/10.0
      c2 = 0.0
      M2 = 2.0E-4
      s2 = 100.0*sqrt(3.0)
      mphi1 = 1E-7
      phi1shift = 100.0
          first_deriv = (/& 
                     mphi1**2*(-phi1shift + phi(1)) + M2**2*Cos(theta2/2.)*&
         (Sin(theta2/2.) + (s2*theta2*Cos(theta2/2.)*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) - &
           Cos(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
              Pi))*(Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) + &
           (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
            Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))&
                     ,&
                    M2**2*Cos(theta2/2.)*(Cos(theta2/2.) - &
                   (s2*theta2*(1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**&
                    2*Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
           (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) + &
          Sin(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
             Pi))*(Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) + &
          (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
           Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))&
                   /)

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
         first_deriv = dVdphi(phi)
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


       !END MULTIFIELD
       case default
          write(*,*) 'MODPK: Need to set first_deriv in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
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
    real(dp) :: lambda4(num_inflaton), alpha2(num_inflaton)

    real(dp) :: lambda2
    real(dp), dimension(num_inflaton,num_inflaton) :: m2_matrix

    if (vnderivs) then
       !MULTIFIELD
       if (size(phi) .ne. 1) then
          write(*,*), 'MODPK: num_inflaton =', num_inflaton
          write(*,*), 'MODPK: 2nd order numerical derivative for multifield not implemented !'
          STOP
       end if
       phiplus = phi(1) + 0.2e0_dp*phi(1)*findiffdphi**(1.e0_dp/4.e0_dp)
       dphi = phiplus - phi(1)
       second_deriv = (pot(phi+2.e0_dp*dphi)+pot(phi-2.e0_dp*dphi)- &
         2.e0_dp*pot(phi))/(4.e0_dp*dphi*dphi)
    else
       second_deriv(:,:) = 0
       
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
      theta2 = Pi/10.0
      c2 = 0.0
      M2 = 2.0E-4
      s2 = 100.0*sqrt(3.0)
      mphi1 = 1E-7
      phi1shift = 100.0
      second_deriv(1,1) =  mphi1**2 + M2**2*Cos(theta2/2.)*(Sin(theta2/2.) +& 
            (s2*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                    Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
               (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
             (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) -& 
            Cos(theta2/2.)*Tan((theta2*Atan(s2*&
                   (Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2 +&
        M2**2*Cos(theta2/2.)*(Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) +& 
           (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
            Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))*&
         ((-2*s2**3*theta2*Cos(theta2/2.)**2*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              (Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2) -& 
           (2*s2*theta2*Cos(theta2/2.)**2*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2)/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) +&
           (2*s2**2*theta2**2*Cos(theta2/2.)**2*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
              Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))/&
           (Pi**2*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2))

       second_deriv(1,2) = M2**2*Cos(theta2/2.)*(Sin(theta2/2.) + &
           (s2*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) - &
           Cos(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
              Pi))*(Cos(theta2/2.) - (s2*theta2*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) + &
           Sin(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
              Pi)) + M2**2*Cos(theta2/2.)*&
         (Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) + &
           (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
            Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))*&
         ((2*s2**3*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2) + &
           (2*s2*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) - &
           (2*s2**2*theta2**2*Cos(theta2/2.)*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
              Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))/&
            (Pi**2*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2))

        second_deriv(2,1) = M2**2*Cos(theta2/2.)*(Sin(theta2/2.) + &
           (s2*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) - &
           Cos(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
              Pi))*(Cos(theta2/2.) - (s2*theta2*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) + &
           Sin(theta2/2.)*Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
              Pi)) + M2**2*Cos(theta2/2.)*&
         (Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) + &
           (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
            Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))*&
         ((2*s2**3*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2) + &
           (2*s2*theta2*Cos(theta2/2.)*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) - &
           (2*s2**2*theta2**2*Cos(theta2/2.)*&
              (1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
              Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))/&
            (Pi**2*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2))

        second_deriv(2,2) = M2**2*Cos(theta2/2.)*(Cos(theta2/2.) -& 
            (s2*theta2*(1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
                  Pi))**2*Sin(theta2/2.)*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
             (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) +& 
            Sin(theta2/2.)*Tan((theta2*Atan(s2*&
                   (Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2 + &
        M2**2*Cos(theta2/2.)*(Cos(theta2/2.)*phi(2) + (-c2 + phi(1))*Sin(theta2/2.) + &
           (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
            Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))*&
         ((-2*s2**3*theta2*(1.0/Cos((theta2*Atan(s2*&
         (Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*Sin(theta2/2.)**2*&
              (Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))*&
              (-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.)))/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2) -& 
       (2*s2*theta2*(1.0/Cos((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/&
                 Pi))**2*Sin(theta2/2.)**2)/&
            (Pi*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)) + &
           (2*s2**2*theta2**2*(1.0/Cos((theta2*&
                   Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))**2*&
              Sin(theta2/2.)**2*(-(Cos(theta2/2.)*(-c2 + phi(1))) + phi(2)*Sin(theta2/2.))*&
              Tan((theta2*Atan(s2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))))/Pi))/&
            (Pi**2*(1 + s2**2*(Cos(theta2/2.)*(-c2 + phi(1)) - phi(2)*Sin(theta2/2.))**2)**2))

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
         second_deriv = d2Vdphi2(phi)
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

       case default
          write(*,*) 'MODPK: Need to set second_deriv in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
       end select
       !END MULTIFIELD
    end if

  END FUNCTION d2Vdphi2

  !Function that parameterizes the turn for quasi--single-field trajectories
  !funct_i(phi_light)
  !Function parameters passed here via vparams(4+,:)
  real(dp) function turning_function(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    !heavy_field_index => which heavy field? for picking (light,heavy_i) direction
    !  This is the "i" in funct_i => need heavy_field_index>1

    !turning_choice => pick a function for turn in (light,heavy_i) direction

    !Parameters for turning_function are set in:
    !vparams(j+4,heavy_field_index) = j^th param for turn in heavy_field_index
    !   direction
    

    if (heavy_field_index <2) then
      print*, "Set heavy_field_index>1. heavy_field_index=", heavy_field_index
      stop
    end if

    select case(turning_choice)
    case(1)
      !No turn, effectively single-field 
      turning_function = 0e0_dp
    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      if (size(vparams,1) <6) then
        print*, "Not enough vparams to set turning_choice=", turning_choice
        stop
      end if

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

      turning_function = sqrt( (focal_length**2)*(sin(asympt_angle)**2) + &
        ((phi - offset_phi)**2) * (tan(asympt_angle)**2) ) &
        -focal_length*sin(asympt_angle)

    case default
       write(*,*) 'MODPK: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select

  end function turning_function

  !Function that parameterizes the turn for quasi--single-field trajectories
  real(dp) function dturndphi(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    select case(turning_choice)
    case(1)
      dturndphi = 0e0_dp
    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

#define funct turning_function(phi,heavy_field_index,turning_choice)
      dturndphi = (tan(asympt_angle)**2)* &
        (phi - offset_phi)/(funct+focal_length*sin(asympt_angle))
#undef funct

    case default
       write(*,*) 'MODPK: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select


  end function dturndphi

  !Function that parameterizes the turn for quasi--single-field trajectories
  real(dp) function d2turndphi2(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    select case(turning_choice)
    case(1)
      d2turndphi2 = 0e0_dp

    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

#define funct turning_function(phi,heavy_field_index,turning_choice)
#define dfunct dturndphi(phi,heavy_field_index,turning_choice)

      d2turndphi2 = ((tan(asympt_angle)**2)/&
        (funct+focal_length*sin(asympt_angle))) * &
        (1.0e0_dp  - (phi-offset_phi)*&
        (dfunct/ (funct &
        +focal_length*sin(asympt_angle))))

#undef funct
#undef dfunct


    case default
       write(*,*) 'MODPK: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select


  end function d2turndphi2


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
      print*, "ERROR: epsilon =", getEps, ">=3.0"
      print*, "ERROR: in getEps"
    end if

    !END MULTIFIELD

  END FUNCTION getEps


  FUNCTION getH(phi,dphi)
    !
    !     Returns H given phi and dphi/dalpha
    !
    real(dp) :: getH
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    ! MULTIFIELD
    getH=SQRT(pot(phi)/3.0e0_dp/M_Pl**2 / &
      (1.0e0_dp - dot_product(dphi, dphi)/6.0e0_dp/M_Pl**2))
    ! MULTIFIELD
  END FUNCTION getH

  !For when using t-integrator
  FUNCTION getH_with_t(phi,dphidt)
    !
    !     Returns H given phi and dphi/dt
    !
    real(dp) :: getH_with_t
    real(dp), INTENT(IN) :: phi(:), dphidt(:)

    ! MULTIFIELD
    getH_with_t= sqrt((pot(phi) + 0.5e0_dp*dot_product(dphidt, dphidt))/ &
      3.0e0_dp/M_pl**2)
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
      print*, "ERROR: epsilon =", getEps_with_t, ">3.0"
      print*, "ERROR: in getEps_with_t"
      stop
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
  subroutine powerspectrum(psi, dpsi, phi, dphi, a, power_spectrum, using_q)
    use internals
    use powersp

    type(power_spectra), intent(inout) :: power_spectrum
    real(dp), dimension(:), intent(in) :: dphi, phi

    logical, optional, intent(in) :: using_q
    logical :: use_q

    ! Hacked matrix to vect
    complex(dp), dimension(:), intent(in) :: psi, dpsi
    real(dp), intent(in) :: a

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
    do i=1,numb_infl; do j=1, numb_infl
      ptb_matrix(i,j) = psi((i-1)*numb_infl+j)
      dptb_matrix(i,j) = dpsi((i-1)*numb_infl+j)
    end do; end do

    ! Make projection vector along adiabatic and isocurv directions
    ! NB: phi_dot_0_scaled = sqrt(2*epsilon) = phi_dot_0/H
    phi_dot_0_scaled = sqrt(dot_product(dphi,dphi))
    omega_z = dphi/phi_dot_0_scaled

    !DEBUG
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
    power_matrix=0e0_dp
    d_power_matrix=0e0_dp
    cross_matrix=0e0_dp

    if (use_q) then
      !Don't divide out by scalefact
      prefactor= (k**3/2.0e0_dp/(pi**2))/(2e0_dp*k)
    else
      prefactor= (k**3/2.0e0_dp/(pi**2)/a**2)/(2e0_dp*k)
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
        !DEBUG
        !print*, "-------------------------"
        !print*, "s_iso",s_iso(i,j)
        !print*, "s_iso",s_iso(i,ll)
        !print*, "power_matrix",power_matrix(j,ll)
        !print*, "pk_iso_vect",pk_iso_vect(i)
        !print*, "addition",s_iso(i,j)*s_iso(i,ll)*power_matrix(j,ll)
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
      power_pnad =0e0_dp

      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)

      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      call pnad_sumAA%clear()
      call pnad_sumAB%clear()
      call pnad_sumBA%clear()
      call pnad_sumBB%clear()

      do i=1,numb_infl; do j=1,numb_infl

        AAprod = AAprod +A_vect(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect(j)*d_power_matrix(i,j)

        !DEBUG
        call pnad_sumAA%add(real(A_vect(i)*A_vect(j)*power_matrix(i,j)))
        call pnad_sumAB%add(real(A_vect(i)*B_vect(j)*cross_matrix(i,j)))
        call pnad_sumBA%add(real(B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))))
        call pnad_sumBB%add(real(B_vect(i)*B_vect(j)*d_power_matrix(i,j)))

      end do; end do
      power_pnad = (AAprod + BBprod) + (ABprod + BAprod)

      !The values (AA + BB) --> -(AB+BA) as approaches adiab limit.
      !Taking diff of "large" numbs means large error in the difference
      !Check that power_pnad is smaller than DP accuracy and set to zero

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

    power_spectrum%matrix  =  power_matrix
    power_spectrum%adiab   =  power_adiab
    power_spectrum%isocurv =  power_isocurv
    power_spectrum%pnad    =  power_pnad
    power_spectrum%entropy =  power_entropy
    power_spectrum%cross_ad_iso =  power_cross


    contains

    !subroutine build_isocurv_basisALTERNATIVE(omega)
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

      !DEBUG
      print*, "from build_isocurv_basisALTERNATIVE"
      print*, "s_iso", s_iso(1,:)
      print*, "omega_z", omega_z
      print*, "s1.omega_z",dot_product(s_iso(1,:),omega_z)
      print*, "first comp", s_iso(1,1)*omega_z(1)
      print*, "second comp", s_iso(1,2)*omega_z(2)
      print*, "norm(omega)", norm(omega_z)
      print*, "norm(s1)", norm(s_iso(1,:))
      !stop

       
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
        print*, "It appears no field space directions have projection"
        print*, "along the adiab direction."
        stop
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
          print*, "spanning(",i,",:) has zero norm."
          stop
        end if

        !DEBUG
        !If there's a major hierarchy in scales for the vector components, then
        !there can be a spurious component of the isocurvature direction
        !oriented along the adiabatic one.  When you're approaching the
        !adiabatic limit this can be the dominant contribution, so we force the
        !smallest term to be zero whenever adding it doesn't affect the value of
        !the norm.
        !call renormalize_remove_smallest(s_iso(i-1,:))


        if (abs(dot_product(omega_z,s_iso(i-1,:)))>1e-12 .or.&
          isnan(abs(dot_product(omega_z,s_iso(i-1,:))))) then

          print*, "Isocurv projection has large adiab component."
          write(*,*), "omega_z.s_iso =",dot_product(omega_z,s_iso(i-1,:))," for i=",i-1
          stop
        end if

      end do


      !Get the rate of change of isocurv directions
      !DEBUG
      !if (size(s_iso)<2) return

      !d_s_iso = 0e0_dp
      !field_basis = phi_hat(2:size(phi_hat),:)


      !!Build the normalization vector and its deriv
      !!Use to build the ds_iso/dalpha vectors
      !normalization = 1e0_dp
      !dnormalization = 0e0_dp
      !d_s_iso = 0e0_dp
      !do i=1,size(s_iso,1)
      !  if (i==1) then
      !    normalization(i) = normalization(i) - (sum(field_basis(i,:)*omega_z))**2

      !    dnormalization(i) = -2e0_dp* &
      !      sum(field_basis(i,:)*omega_z)*sum(field_basis(i,:)*d_omega_z)

      !    d_s_iso(i,:) = (-1e0_dp/normalization)*&
      !      (sum(field_basis(i,:)*d_omega_z)*omega_z +&
      !      sum(field_basis(i,:)*omega_z)*d_omega_z + &
      !      s_iso(i,:)*dnormalization(i))
      !  else
      !    do j=1, i-1

      !      normalization(i) = normalization(i) - &
      !        (sum(field_basis(i,:)*s_iso(j,:)))**2

      !      dnormalization(i) = dnormalization(i) - 2e0_dp* &
      !        sum(field_basis(i,:)*s_iso(j,:))*sum(field_basis(i,:)*d_s_iso(j,:))
      !    end do

      !    do j=1, i-1
      !      d_s_iso(i,:) = d_s_iso(i,:) + &
      !        sum(field_basis(i,:)*d_s_iso(j,:))*s_iso(j,:) + &
      !        sum(field_basis(i,:)*s_iso(j,:))*d_s_iso(j,:)
      !    end do

      !    d_s_iso(i,:) = (-1e0_dp/normalization)*(d_s_iso(i,:) + s_iso(i,:)*dnormalization(i))

      !  end if

      !end do


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

      !DEBUG
      !print*, "vect", vect
      !print*, "vect_temp", vect_temp
      !print*, "diff norms", abs(norm(vect)-norm(vect_temp))

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
  subroutine bundle_exp_scalar(this,phi, efold)

    class(bundle) :: this
    real(dp), intent(in) :: phi(:), efold
    real(dp) :: dN

    dN = efold - this%N

    this%dlogThetadN= this%dlogThetadN - &
      dN*trace_d2logVdphi2(phi)

    this%exp_scalar=exp(this%dlogThetadN)

    this%N=efold

  end subroutine bundle_exp_scalar


  FUNCTION MySech(x)
    real(dp)  :: x,MySech

    IF(ABS(x).GT.40.0) THEN
       MySech=0.0
    ELSE
       MySech=1.0/COSH(x)
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

    !DEBUG
    !print*, "stable???", stable, "eps=", eps, "V", V, "using_t", using_t

    !DEBUG
    !print*, "SETTING UNSTABLE SO USE_T"
    !stable=.false.

  end subroutine stability_check_on_H

#undef HEAVY
#undef PHI_I
#undef DELTAPHI
#undef EXPTERM

END MODULE potential
