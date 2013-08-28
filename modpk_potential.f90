MODULE potential
  USE modpkparams
  use internals, only : pi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot, getH, getHdot, getEps, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, geteta, zpower

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


CONTAINS

  FUNCTION MySech(x)
    real(dp)  :: x,MySech

    IF(ABS(x).GT.40.0) THEN
       MySech=0.0
    ELSE
       MySech=1.0/COSH(x)
    END IF
    RETURN
  END FUNCTION MySech

  FUNCTION pot(phi)
    !
    !     Returns V(phi) given phi, phi can be an array for multifield, 
    !     The code implement multifield potential in the form of V = \sum V(phi_i), 
    !     More complicated form of potentials can be customized 
    !
    real(dp) :: pot
    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: m2_V(size(phi)) !! m2_V is the diagonal mass square matrix
    real(dp) :: c1_V(size(phi)) !! some constants for a more flexible potential construction in case 9
    real(dp) :: lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    !real(dp), dimension(size(phi)/2) :: lambda_waterfall, mass_waterfall, &
    !  mass_infl, couple_water_infl
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid
    
    select case(potential_choice)
    case(1) 
    ! MULTIFIELD
       m2_V = 10.d0**(vparams(1,:))
       pot = 0.5*sum(m2_V*phi*phi)
    case(2)
       lambda = 10.d0**vparams(1,:)
       finv = 1.d0/(10.d0**vparams(2,:))
       pot = sum(lambda**4*(1.d0+cos(finv*phi)))
    case(3)
       lambda = 10.d0**vparams(1,:)
       pot = sum(0.25*lambda*phi**4)
    case(4)
       lambda = 10.d0**vparams(1,:)
       pot = sum(lambda*phi)
    case(5)
       lambda = 10.d0**vparams(1,:)
       pot = sum(lambda*1.5d0*phi**(2./3.))
    case(6)
       lambda = 10.d0**vparams(1,:)
       mu = 10.d0**vparams(2,:)
       pot = sum(lambda**4 - mu*phi**4/4.d0)
    case(7) !product of exponentials
       lambda = vparams(2, :)
       pot = vparams(1,1)*M_Pl**4 * exp(dot_product(lambda, phi/M_Pl)) 
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
       pot = (lambda_hybrid**4)*((1.0_dp - phi(1)**2/mass_hybrid**2)**2 +&
         phi(2)**2/mu_hybrid**2 +&
         phi(1)**2*phi(2)**2/nu_hybrid**4)
     case(9)
       m2_V = vparams(1,:)
       c1_V = vparams(2,:)
       pot = 0.5*sum(m2_V*phi*phi) + sum(c1_V)
     case(10)
      ! pot = 0.5*M2**2*cos(0.5*theta2)*(sin(0.5*theta2)*(phi(1)-c2)+&
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
      pot=potlarge+potsmall
    case default
       write(*,*) 'MODPK: Need to set pot(phi) in modpk_potential.f90 for potential_choice =',potential_choice
       STOP
    end select

    RETURN
  END FUNCTION pot


  FUNCTION dVdphi(phi)
    !
    !     Returns dV/dPhi given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: dVdphi(size(phi))
    real(dp) :: dphi(size(phi)), phiplus(size(phi))
    real(dp) :: m2_V(size(phi)),  c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i

    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    if (vnderivs) then
       ! MULTIFIELD
       do i=1, num_inflaton
          phiplus = phi
          phiplus(i) = phi(i) + 0.5*phi(i)*findiffdphi**(1./3.)
          dphi = phiplus - phi
          dVdphi(i) = (pot(phi+dphi)-pot(phi-dphi))/(2.*dphi(i))
          if (dVdphi(i).eq.0.d0) then
             write(*,*) 'MODPK: For i=', i
             write(*,*) 'MODPK: dVdphi(i)=0, possibly signaling a problem with accuracy of numerical derivatives.'
             write(*,*) 'MODPK: Try using vnderivs=F if possible.'
             STOP
          end if
       end do
       ! END MULTIFIELD
    else
       select case(potential_choice)
       !MULTIFIELD
       case(1)
          m2_V = 10.d0**(vparams(1,:))
          dVdphi = m2_V*phi
       case(2)
          lambda = 10.d0**vparams(1,:)
          finv = 1.d0/(10.d0**vparams(2,:))
          dVdphi = -lambda**4*finv*sin(finv*phi)
       case(3)
          lambda = 10.d0**vparams(1,:)
          dVdphi = lambda*phi**3
       case(4)
          lambda = 10.d0**vparams(1,:)
          dVdphi = lambda
       case(5)
          lambda = 10.d0**vparams(1,:)
          dVdphi = lambda*phi**(-1./3.)
       case(6)
          mu = 10.d0**vparams(2,:)
          dVdphi = -mu*phi**3
       case(7) 
          lambda = vparams(2,:)
          dVdphi = lambda*vparams(1,1)*M_Pl**3 * exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          dVdphi(1) = (lambda_hybrid**4)*(4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)*phi(1)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)
          dVdphi(2) = (lambda_hybrid**4)*(2.0_dp*phi(2)/mu_hybrid**2 +&
            2.0_dp*phi(1)**2*phi(2)/nu_hybrid**4)
       case(9)!Saddle type things
          m2_V = (vparams(1,:))
          dVdphi = m2_V*phi
       case(10)!Modified Langlois model
      theta2 = Pi/10.0
      c2 = 0.0
      M2 = 2.0E-4
      s2 = 100.0*sqrt(3.0)
      mphi1 = 1E-7
      phi1shift = 100.0
          dVdphi = (/& 
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

       !END MULTIFIELD
       case default
          write(*,*) 'MODPK: Need to set dVdphi in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
       end select

    end if

    RETURN
  END FUNCTION dVdphi


  FUNCTION d2Vdphi2(phi)
    !
    !     Returns d^2V/dPhi^2 given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: d2VdPhi2(size(phi),size(phi))
    real(dp) :: m2_V(size(phi)), c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i, j

    real(dp) :: dphi,phiplus
    !  mass_infl, couple_water_infl
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    if (vnderivs) then
       !MULTIFIELD
       if (size(phi) .ne. 1) then
          write(*,*), 'MODPK: num_inflaton =', num_inflaton
          write(*,*), 'MODPK: 2nd order numerical derivative for multifield not implemented !'
          STOP
       end if
       phiplus = phi(1) + 0.2*phi(1)*findiffdphi**(1./4.)
       dphi = phiplus - phi(1)
       d2Vdphi2 = (pot(phi+2.*dphi)+pot(phi-2.*dphi)-2.*pot(phi))/(4.*dphi*dphi)
    else
       d2Vdphi2(:,:) = 0
       
       select case(potential_choice)
       case(1)          
          m2_V = 10.d0**(vparams(1,:))
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = m2_V(i)   
       case(2)
          lambda = 10.d0**vparams(1,:)
          finv = 1.d0/(10.d0**vparams(2,:))
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = -lambda(i)**4*finv(i)*finv(i)*cos(finv(i)*phi(i))
       case(3)
          lambda = 10.d0**vparams(1,:)
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = 3.d0*lambda(i)*phi(i)**2
       case(4)
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = 0.d0
       case(5)
          lambda = 10.d0**vparams(1,:)
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = -lambda(i)/3.d0*phi(i)**(-4./3.)
       case(6)
          mu = 10.d0**vparams(2,:)
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = -3.d0*mu(i)*phi(i)**2
       case(7) 
          lambda = vparams(2, :)
          forall (i=1:size(phi), j=1:size(phi)) d2Vdphi2(i,j) = &
               lambda(i)*lambda(j)*vparams(1,1)*M_Pl**2 *exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          d2Vdphi2(1,1) = (lambda_hybrid**4)*(4.0_dp*(2.0_dp*phi(1)/mass_hybrid**2)*phi(1)/mass_hybrid**2 +&
          4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)

          d2Vdphi2(2,1) = (lambda_hybrid**4)*(4.0_dp*phi(1)*phi(2)/nu_hybrid**4)

          d2Vdphi2(1,2) = d2Vdphi2(2,1)

          d2Vdphi2(2,2) = (lambda_hybrid**4)*(2.0_dp/mu_hybrid**2 +&
            2.0_dp*phi(1)**2/nu_hybrid**4)
       case(9)          
          m2_V = vparams(1,:)
          forall (i = 1:size(phi)) d2Vdphi2(i,i) = m2_V(i)
       case(10)
      theta2 = Pi/10.0
      c2 = 0.0
      M2 = 2.0E-4
      s2 = 100.0*sqrt(3.0)
      mphi1 = 1E-7
      phi1shift = 100.0
      d2Vdphi2(1,1) =  mphi1**2 + M2**2*Cos(theta2/2.)*(Sin(theta2/2.) +& 
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

       d2Vdphi2(1,2) = M2**2*Cos(theta2/2.)*(Sin(theta2/2.) + &
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

        d2Vdphi2(2,1) = M2**2*Cos(theta2/2.)*(Sin(theta2/2.) + &
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

        d2Vdphi2(2,2) = M2**2*Cos(theta2/2.)*(Cos(theta2/2.) -& 
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

       case default
          write(*,*) 'MODPK: Need to set d2Vdphi2 in modpk_potential.f90 or use numerical derivatives (vnderivs=T)'
          STOP
       end select
       !END MULTIFIELD
    end if

    RETURN
  END FUNCTION d2Vdphi2


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


    Ninit = 70.d0

    if (size(phi0) .gt. 1) then 
       phii = phi0 ! MULTIFIELD
    else  !SINGLE FIELD
       select case(potential_choice)
       case(1)
          phii = 2.d0*sqrt(Ninit+0.5d0)
       case(2)
          finv = 1.d0/(10.d0**vparams(2,1))
          phii = 2.d0/finv*asin(exp(-0.5d0*Ninit*finv*finv)/ &
               sqrt(1.d0+0.5d0*finv*finv))
       case(3)
          phii = sqrt(8.d0*(Ninit+1.d0))
       case(4)
          phii = sqrt(2.d0*Ninit+0.5d0)
       case(5)
          phii = sqrt(4.d0/3.d0*Ninit+2.d0/9.d0)
       case(6)
          lambda = 10.d0**vparams(1,1)
          mu = 10.d0**vparams(2,1)
          x1 = lambda**4/mu
          phesq = ((sqrt(2.d0)*x1)**(-4./3.)+1.d0/(4.d0*x1))**(-0.5)
          if (vparams(1,1)<-3.d0) then
             phii = sqrt(phesq/(2.d0**1.5*Ninit/sqrt(phesq)+1.d0))
          else
             x2 = 4.d0*Ninit + 2.d0*x1/phesq + 0.5d0*phesq
             phii = sqrt(x2)*sqrt(1.d0-sqrt(1.d0-4.d0*x1/x2/x2))
          end if
       case default
          phii = phi0
       end select
    end if

    initialphi = phii

    RETURN
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
    getEps = 0.5d0*(M_Pl)**2 * dot_product(dphi,dphi)
    !!getEps = 2.d0*(M_Pl**2)*(((getHdot(phi,dphi)/dphi)/getH(phi,dphi))**2)
    !END MULTIFIELD
    RETURN

  END FUNCTION getEps


  FUNCTION getH(phi,dphi)
    !
    !     Returns H given phi and dphi/dalpha
    !
    real(dp) :: getH
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    ! MULTIFIELD
    getH=SQRT(pot(phi)/3./M_Pl**2 / (1.0 - dot_product(dphi, dphi)/6.0/M_Pl**2))
    ! MULTIFIELD
    RETURN
  END FUNCTION getH


  FUNCTION getHdot(phi,dphi)
    !
    !     Returns dH/dalpha given phi and dphi/dalpha
    !
    real(dp) :: getHdot
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    ! MULTIFIELD
    getHdot = -dot_product(dphi, dphi) * getH(phi,dphi)/2./M_Pl**2
    ! END MULTIFIELD
    RETURN
  END FUNCTION getHdot


  FUNCTION getdepsdalpha(phi,dphi)
    !
    !    Returns depsilon/dalpha given phi and dphi/dalpha
    !    Gets this by differentiating Peiris et al Eq A3 (2003)
    !
    real(dp) :: getdepsdalpha, H, dHdalpha, eps
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    
    H=getH(phi,dphi)
    dHdalpha=getHdot(phi,dphi)
    eps=getEps(phi,dphi)

    ! MULTIFIELD
    getdepsdalpha=6.*dHdalpha/H*(1.-eps/3.)-dot_product(dVdphi(phi), dphi)/(H**2*M_Pl**2)
    ! END MULTIFIELD

    RETURN
  END FUNCTION getdepsdalpha

  FUNCTION geteta(phi, dphi)
    !
    !    Return the eta parameter eta = deps/dalpha / eps
    !
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    real(dp) :: geteta, eps

    eps = getEps(phi, dphi)
    geteta = getdepsdalpha(phi, dphi) / eps

    RETURN
  END FUNCTION geteta

  ! Powerspectrum for psi ptbs mode matrix, outputting full IJ matrix,
  !adiabatic P(k), and isocurv P(k); includes cross correlations
  !subroutine powerspectrum(psi, dpsi, phi, dphi, a, power_matrix, power_adiab, power_isocurv)
  subroutine powerspectrum(psi, dpsi, phi, dphi, a, power_spectrum)
    use internals
    use powersp

    type(power_spectra), intent(inout) :: power_spectrum
    real(dp), dimension(:), intent(in) :: dphi, phi
    ! Hacked matrix to vect
    complex(dp), dimension(:), intent(in) :: psi, dpsi
    real(dp), intent(in) :: a

    complex(dp), dimension(size(dphi),size(dphi)) :: power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: d_power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: cross_matrix
    real(dp) :: power_adiab
    real(dp) :: power_isocurv
    real(dp) :: power_pnad
    real(dp) :: power_entropy

    real(dp) :: AAprod, ABprod, BAprod, BBprod

    real(dp) :: H, Pdot

    !Pre-factors for Pnad P(k)
    real(dp), dimension(num_inflaton) :: A_vect, B_vect

    ! size(dphi)=num_inflaton
    complex(kind=dp), dimension(size(dphi),size(dphi)) :: psi_matrix,&
      dpsi_matrix
    real(dp), dimension(size(dphi)-1) :: pk_iso_vect

    !proj along adiab dir
    real(dp), dimension(size(dphi)) :: omega_z
    !proj along isocurv dirs
    real(dp), dimension(size(dphi)-1,size(dphi)) :: s_iso
    real(dp), parameter :: div_zero_tol=1e-20_dp

    real(dp) :: zeta2, phi_dot_0_scaled, prefactor
    integer :: numb_infl
    integer :: i, j, ll

    real(dp) :: prod_exponent

    numb_infl=size(dphi)

    ! Convert hacked vector to matrix
    do i=1,numb_infl; do j=1, numb_infl
      psi_matrix(i,j) = psi((i-1)*numb_infl+j)
      dpsi_matrix(i,j) = dpsi((i-1)*numb_infl+j)
    end do; end do

    ! Make projection vector along adiabatic and isocurv directions
    ! NB: phi_dot_0_scaled = sqrt(2*epsilon) = phi_dot_0/H
    phi_dot_0_scaled = sqrt(dot_product(dphi,dphi))
    omega_z = dphi/phi_dot_0_scaled

    !Build power matrices for fields, vels, and cross correls
    power_matrix=0e0_dp
    d_power_matrix=0e0_dp
    cross_matrix=0e0_dp

    prefactor= (k**3/2.0e0_dp/(pi**2)/a**2)/(2e0_dp*k)
    do i=1,numb_infl;do j=1, numb_infl; do ll=1,numb_infl

      power_matrix(i,j) =power_matrix(i,j)+ &
        psi_matrix(i,ll)*conjg(psi_matrix(j,ll))*prefactor

      d_power_matrix(i,j) =d_power_matrix(i,j)+&
        dpsi_matrix(i,ll)*conjg(dpsi_matrix(j,ll))*prefactor

      cross_matrix(i,j) =cross_matrix(i,j)+ &
        psi_matrix(i,ll)*conjg(dpsi_matrix(j,ll))*prefactor

    end do; end do; end do

    !Adiabatic power spectrum
    power_adiab = 0e0_dp
    do i=1,numb_infl; do j=1,numb_infl
      power_adiab = power_adiab + omega_z(i)*omega_z(j)*power_matrix(i,j)
    end do; end do
    power_adiab = (1e0_dp/phi_dot_0_scaled**2)*power_adiab

    !Isocurvature power spectra
    if (numb_infl>1) then

      !Build the isocurv unit vects
      call build_isocurv_basis()

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

      !P(k) of total non-adiab pressure ptbs
      !dP_nad_i(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*Psi_ij)*\hat{a}_j
      power_pnad =0e0_dp

      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)

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
    end if

    !DEBUG [JF]
    !write(21, *) Log(a)-Log(a_init), power_adiab, power_isocurv, power_pnad,&
    !power_entropy

    !write(29, *)Log(a)-Log(a_init),power_isocurv
    !write(28, *) Log(a)-Log(a_init),power_pnad
    !write(27, *) Log(a)-Log(a_init),power_entropy

    power_spectrum%matrix  =  power_matrix
    power_spectrum%adiab   =  power_adiab
    power_spectrum%isocurv =  power_isocurv
    power_spectrum%pnad    =  power_pnad
    power_spectrum%entropy =  power_entropy


    contains


    !Basis perpend to adiab direction.  Gram-Schmidt
    subroutine build_isocurv_basis()

      real(dp), dimension(size(dphi),size(dphi)) :: spanning
      !field dirs
      real(dp), dimension(size(dphi),size(dphi)) :: phi_hat
      real(dp), dimension(size(dphi)) :: phi_hat_temp
      real(dp) :: check_dot
      integer :: i, j, adiab_index

      !Build field space vects
      do i=1,size(phi_hat,1); do j=1,size(phi_hat,2)
        if (i==j) then
          phi_hat(i,j) = 0e0_dp
        else
          phi_hat(i,j) = 1e0_dp
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
        print*, "It appears all field space directions are aligned with"
        print*, "the adiab direction."
        stop
      end if

      !Set the "most adiab" dir to 1
      if (adiab_index /= 1) then
        phi_hat_temp = phi_hat(1,:)
        phi_hat(1,:) = phi_hat(adiab_index,:)
        phi_hat(adiab_index,:) = phi_hat_temp
      end if

      spanning = 0e0_dp
      spanning(1,:)=omega_z/norm(omega_z)

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
          print*, "spanning(",i,",:) has zero norm..."
          stop
        end if

        if (abs(dot_product(omega_z,s_iso(i-1,:)))>1e-12) then
          print*, "Isocurv projection has large adiab component."
          write(*,*), "omega_z.s_iso =",dot_product(omega_z,s_iso(i-1,:))," for i=",i-1
          stop
        end if
      end do


    end subroutine build_isocurv_basis

    !Functions for calculating the total non-adiabatic pressure perturbation
    !matrices

    !dP_nad(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*Psi_ij)*\hat{a}_j
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
    tensorpower = abs(v)**2/(2*k) / a**2 * (k**3)*4./(PI**2)/(M_Pl**2)
    !END MULTIFIELD

  END FUNCTION tensorpower

  pure FUNCTION zpower(u_zeta, dsigma, a)
    USE internals
    real(dp) :: zpower
    real(dp), INTENT(IN) :: dsigma
    real(dp), INTENT(IN) :: a
    COMPLEX(KIND=DP), INTENT(IN) :: u_zeta

    zpower = abs(u_zeta**2)/dsigma**2/a**2 /(2*k) * (k**3)/(2*PI**2)

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



END MODULE potential
