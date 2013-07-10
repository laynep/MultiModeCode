MODULE potential
  USE modpkparams
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot, getH, getHdot, getEps, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, geteta, zpower

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
    real(dp) :: lambda(size(phi)), finv(size(phi)), mu(size(phi))
    
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
    ! END MULTIFIELD
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
    real(dp) :: m2_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    integer :: i

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
    real(dp) :: m2_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    integer :: i, j

    real(dp) :: dphi,phiplus

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

  ![ LP: ] Powerspectrum for psi ptbs mode matrix, outputting full IJ matrix,
  !adiabatic P(k), and isocurv P(k); includes cross correlations
  !MULTIFIELD: Calculates the adiabatic perturbation P_R(k) given psi, dphi/dalpha, a
  subroutine powerspectrum(psi, dphi, a, power_matrix, power_adiab, power_isocurv)
    use internals
    real(dp), dimension(:,:), intent(out) :: power_matrix
    real(dp), intent(out) :: power_adiab
    real(dp), intent(out), optional :: power_isocurv

    real(dp), dimension(:), intent(in) :: dphi
    ![ LP: ] Hacked matrix to vect
    complex(kind=dp), dimension(:), intent(in) :: psi
    real(dp), intent(in) :: a

    ![ LP: ] size(dphi)=num_inflaton
    complex(kind=dp), dimension(size(dphi),size(dphi)) :: psi_matrix
    real(dp), dimension(size(dphi)) :: omega_z![ LP: ] proj along adiab dir
    real(dp), dimension(size(dphi)) :: s_iso![ LP: ] proj along isoc dir
    real(dp) :: zeta2, phi_dot_0_scaled
    integer :: numb_infl
    integer :: i, j, ll

    numb_infl=size(dphi)

    ![ LP: ] Convert hacked vector to matrix
    do i=1,numb_infl; do j=1, numb_infl
      psi_matrix(i,j) = psi((i-1)*numb_infl+j)
    end do; end do

    ![ LP: ] Make projection vector along adiabatic and isocurv directions
    ![ LP: ] NB: phi_dot_0_scaled = sqrt(2*epsilon) = phi_dot_0/H
    phi_dot_0_scaled = sqrt(dot_product(dphi,dphi))
    omega_z = dphi/phi_dot_0_scaled

    if (present(power_isocurv)) then
      power_isocurv=0e0_dp
      print*, "isocurv power not working yet"
      stop
      !s_iso = XXX
    end if

    power_matrix=0e0_dp
    do i=1,numb_infl;do j=1, numb_infl; do ll=1,numb_infl
      power_matrix(i,j) =power_matrix(i,j)+ (k**3/2.0e0_dp/(pi**2)/a**2)*psi_matrix(i,ll)*conjg(psi_matrix(j,ll))/(2e0_dp*k)
    end do; end do; end do

    !DEBUG [JF]
    !print*, "power_matrix", power_matrix
    write(3, *), Log(a)-Log(a_init), power_matrix

    power_adiab = 0e0_dp
    !if (present(power_isocurv)) power_isocurv = 0e0_dp
    ![ LP: ] Project power spectra
    do i=1,numb_infl; do j=1,numb_infl
      power_adiab = power_adiab + omega_z(i)*omega_z(j)*power_matrix(i,j)
      !if (present(power_isocurv)) power_isocurv = &
      !power_isocurv + s_iso(i)*s_iso(j)*power_matrix(i,j)
    end do; end do
    power_adiab = (1e0_dp/phi_dot_0_scaled**2)*power_adiab
    !if (present(power_isocurv)) power_isocurv = (1e0_dp/phi_dot_0_scaled**2)*power_isocurv

    !DEBUG [JF]
    !print*, "powerspectrum", power_adiab
    !END MULTIFIELD

    !DEBUG [JF]
     write(20, *), Log(a)-Log(a_init), power_adiab

  end subroutine powerspectrum


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

END MODULE potential
