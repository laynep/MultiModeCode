    MODULE access_modpk
  USE camb_interface
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: potinit, evolve, total_efold

! Number of k values for computing spline of P(k).
! Set to lower numbers for smoother potentials.
  INTEGER*4, PARAMETER, PUBLIC :: pkspline_n = 500

  DOUBLE PRECISION, PARAMETER, PUBLIC :: pkspline_kmin = log(1.d-5), pkspline_kmax = log(5.d0)
  DOUBLE PRECISION, PUBLIC :: pkspline_k(pkspline_n), pkspline_p(pkspline_n), &
	pkspline_p2der(pkspline_n), pkspline_pt(pkspline_n), &
	pkspline_pt2der(pkspline_n)

CONTAINS

  SUBROUTINE potinit
    USE modpkparams
    USE powersp
    USE background_evolution, ONLY : backgrnd
    USE potential, ONLY : initialphi
    IMPLICIT NONE
    DOUBLE PRECISION :: k,klo,khi

    k_start = 1.d2
    eval_ps = 5.d2
    useq_ps = 1.d2
    !
    !     Solve the background equations
    !
    pk_bad = 0
    phi_init = initialphi(phi_init0)
    CALL backgrnd
    RETURN
  END SUBROUTINE potinit

  SUBROUTINE total_efold
    USE modpkparams
    USE background_evolution, ONLY : backgrnd_efold
    USE potential, ONLY : initialphi
    IMPLICIT NONE
 
    pk_bad = 0
    phi_init = initialphi(phi_init0)
    CALL backgrnd_efold
    RETURN
  END SUBROUTINE total_efold

  SUBROUTINE evolve(kin, pow_adiabatic, pow_isocurvature, powt, powz)
    USE modpk_odeint
    USE ode_path
    USE modpkparams
    USE internals
    USE powersp
    USE potential, ONLY: pot,powerspectrum, dVdphi, getH, getHdot
    USE modpk_utils, ONLY : locate, polint, derivs, qderivs, rkqs_c, array_polint
    IMPLICIT NONE

    !!INTEGER*4, PARAMETER :: NVAR = 10

    INTEGER*4 :: i,j
    DOUBLE PRECISION :: accuracy,h1,hmin,x1,x2 
    COMPLEX(KIND=DP), DIMENSION(2*num_inflaton + 2*(num_inflaton**2)+4) :: y 
    double precision :: identity(num_inflaton**2)
    DOUBLE PRECISION, INTENT(IN) :: kin
    DOUBLE PRECISION, INTENT(OUT) :: pow_adiabatic, pow_isocurvature, powt, powz

    ![ LP: ] Upgrade pow --> pow(I,J)
    double precision, dimension(:,:) :: power_matrix

    DOUBLE PRECISION :: dum, ah, alpha_ik, dalpha, dh
    DOUBLE PRECISION, DIMENSION(num_inflaton) :: p_ik,delphi

    !     Set initial conditions
    !     x = alpha     e-folds
    ![ LP: ] Background
    !     y(1:n) = phi                 dydx(1:n)=dphi/dalpha
    !     y(n+1:2n) = dphi/dalpha      dydx(n+1:2n)=d^2phi/dalpha^2

    ![ LP: ] Mode matrix ptb, psi_IJ
    ![ LP: ] NB: the psi portion of the y-vector is 1:n=psi_1(1:n) and
    ![ LP: ] NB: 2n:3n=psi_2(1:n), etc.
    !     y(2n+1:2n+n**2) = psi               dydx(2n+1:3n)=dpsi/dalpha
    !     y(2n+n**2+1:2n+2n**2) = dpsi/dalpha       dydx(3n+1:4n)=d^2psi/dalpha^2

    ![ LP: ] Tensors
    !     y(2n+2n**2+1) = v                  dydx(4n+1)=dv/dalpha
    !     y(2n+2n**2+2) = dv/dalpha          dydx(4n+2)=d^2v/dalpha^2

    ! --- u_zeta is the adiabatic mode ignoring coupings to other modes, used to compare with the full zeta perturbation
    !! --- full_zeta - u_zeta gives the super-horizon evolution 
    !     y(2n+2n**2+3) = u_zeta     dydx(4n+4)=d^2u_zeta/dalpha^2  
    !     y(2n+2n**2+4) = du_zeta/dalpha     dydx(4n+4)=d^2u_zeta/dalpha^2  


    ![ LP: ] Make the powerspectrum array.
    if (allocated(pow_ptb_ij)) then
      deallocate(pow_ptb_ij)
    end if
    allocate(pow_ptb_ij(num_inflaton,num_inflaton))

    k=kin*Mpc2Mpl

    ah=LOG(k/k_start)     !! start where k = 100 aH, deep in the horizon, ah = log(aH)
    i= locate(aharr(1:nactual_bg), ah)
    
    IF(i.eq.0.) THEN
       PRINT*,'MODPK: The background solution worked, but the k you requested is outside'
       PRINT*,'MODPK: the bounds of the background you solved for. Please reconsider'
       PRINT*,'MODPK: your phi_init and N_pivot combo.'
       PRINT*,'MODPK: QUITTING'
       print*, ah, aharr(1)
       STOP
    END IF

    ! nactual_bg here is set by the background evolution
    j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
    
    !MULTIFIELD
    CALL array_polint(aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
    !END MULTIFIELD
    CALL polint(aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
    CALL polint(aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)
    a_ik=EXP(alpha_ik)*a_init
    x1=alpha_ik

    IF(x1.le.0.) THEN
       PRINT*,'MODPK: The phi_init you specified is too small to give'
       PRINT*,'MODPK: sufficient efolds of inflation. We cannot self-consistently'
       PRINT*,'MODPK: solve this for you. Please adjust phi_init and try again.'
       PRINT*,'MODPK: QUITTING'
       STOP
    END IF
    IF((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1 .OR. dalpha .GT. 0.1 .OR. dh .GT. 0.1) THEN
       PRINT*,'MODPK: The interpolation in SUBROUTINE EVOLVE has suspiciously large'
       PRINT*,'MODPK: errors. Your model smells fishy.'
       PRINT*,'MODPK: QUITTING'
       STOP
    ENDIF

    ![ LP: ] Set the initial conditions.
    ![ LP: ] Make an identity vector analog of identity matrix
    call make_identity(identity)


    ![ LP: ] Background - from previous evolution
    y(1:num_inflaton) = cmplx(p_ik)             !phi(x1)
    y(num_inflaton+1:2*num_inflaton) = cmplx(-dVdphi(p_ik)/3./h_ik/h_ik)  !dphi/dalpha(x1) slowroll approx
    ![ LP: ] mode matrix - diagonalize
    y(2*num_inflaton+1:2*num_inflaton+num_inflaton**2) = (1.d0, 0)*identity  !cmplx(1/sqrt(2*k))
    y(2*num_inflaton+num_inflaton**2+1:2*num_inflaton+2*num_inflaton**2) = &
           cmplx(0., -k/exp(ah))*identity
    ![ LP: ] tensors - ???
    y(2*num_inflaton+2*num_inflaton**2+1) = (1.d0, 0) !cmplx(1/sqrt(2*k))
    y(2*num_inflaton+2*num_inflaton**2+2) = cmplx(0., -k/exp(ah))  
    ![ LP: ] u_zeta - ???
    y(2*num_inflaton+2*num_inflaton**2+3) = (1.d0, 0) !cmplx(1/sqrt(2*k))
    y(2*num_inflaton+2*num_inflaton**2+4) = cmplx(0., -k/exp(ah))  

    !     Call the integrator
    !
    ode_underflow = .false.
    ode_ps_output = .true.
    ode_infl_end = .true.
    save_steps = .true. 
    pk_bad = 0

    !MULTIFIELD, need to evolve towards the end of inflation
    x2 = lna(nactual_bg) + 5.d0
    !END MULTIFIELD

    h1=0.05 !guessed start stepsize
    accuracy=1.0d-9 !4.0d-2 !2!6 !has a big impact on the speed of the code
    hmin=0.0 !minimum stepsize

    CALL odeint(y, x1, x2, accuracy, h1, hmin, derivs, qderivs, rkqs_c)
    nactual_mode = kount  ! update nactual after evolving the modes
    
    IF(.NOT. ode_underflow) THEN 
       !power_matrix = pow_ptb_ij
       powt = powt_ik
       powz = powz_ik
       pow_adiabatic = pow_adiab_ik
       pow_isocurvature = pow_isocurv_ik
    ELSE
       pow=0.
       powt=0.
       pk_bad=1
    ENDIF


    RETURN

    contains

      ![ LP: ] A "packed" identity vector analog of identity matrix.
      subroutine make_identity(identityvector)

        double precision, dimension(:), intent(out) :: identityvector
        integer :: i, j

        identityvector=0d0

        do i=1,num_inflaton; do j=1, num_inflaton
          if (i==j) then
            identityvector((i-1)*num_inflaton+j)=1d0
          end if
        end do; end do

      end subroutine make_identity

  END SUBROUTINE evolve

END MODULE access_modpk
