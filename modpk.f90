MODULE access_modpk
  use modpkparams, only : dp
  USE camb_interface
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: potinit, evolve, total_efold

! Number of k values for computing spline of P(k).
! Set to lower numbers for smoother potentials.
  INTEGER*4, PARAMETER, PUBLIC :: pkspline_n = 500

  real(dp), PARAMETER, PUBLIC :: pkspline_kmin = log(1.e-5_dp), pkspline_kmax = log(5.e0_dp)
  real(dp), PUBLIC :: pkspline_k(pkspline_n), pkspline_p(pkspline_n), &
	pkspline_p2der(pkspline_n), pkspline_pt(pkspline_n), &
	pkspline_pt2der(pkspline_n)

CONTAINS

  SUBROUTINE potinit
    USE modpkparams
    USE powersp
    USE background_evolution, ONLY : backgrnd
    USE potential, ONLY : initialphi
    use modpk_icsampling, only : bad_ic, sampling_techn, reg_samp
    IMPLICIT NONE
    real(dp) :: k,klo,khi

    k_start = 1.d2
    eval_ps = 5.d2
    useq_ps = 1.d2
    !
    !     Solve the background equations
    !
    pk_bad = 0
    phi_init = initialphi(phi_init0)

    !NB: For eqen sampling, dphi_init set in trial_background
    CALL backgrnd

    !When scanning ICs, let some backgrnd errors override
    if (sampling_techn/=reg_samp) then
      if (pk_bad == bad_ic) then
        print*, "--------------- BAD IC; RESTARTING -------------------"
        return
      end if
    end if

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

  SUBROUTINE evolve(kin, powerspectrum_out)
    USE modpk_odeint
    USE ode_path
    USE modpkparams
    USE internals
    USE powersp
    USE potential, ONLY: pot,powerspectrum, dVdphi, getH, getdHdalpha, field_bundle
    USE modpk_utils, ONLY : locate, polint, derivs, qderivs, rkqs_c, array_polint
    use modpk_icsampling, only : bad_ic, sampling_techn, reg_samp
    IMPLICIT NONE

    !!INTEGER*4, PARAMETER :: NVAR = 10
    type(power_spectra), intent(out) :: powerspectrum_out

    integer*4 :: i,j
    real(dp) :: accuracy,h1,hmin,x1,x2 
    complex(kind=dp), dimension(2*num_inflaton + 2*(num_inflaton**2)+4) :: y
    real(dp), INTENT(IN) :: kin
    real(dp) :: pow_isocurvature
    real(dp) :: pow_adiabatic,  powt, powz
    real(dp), dimension(:,:), allocatable :: power_matrix
    real(dp) :: dum, ah, alpha_ik, dalpha, dh
    real(dp), DIMENSION(num_inflaton) :: p_ik,delphi
    real(dp), DIMENSION(num_inflaton) :: dp_ik

    character(36) :: e2_fmt = '(a25, es12.4, a3, es11.4, a1)'

    ! the psi portion of the y-vector is 1:n=psi_1(1:n) and
    ! 2n:3n=psi_2(1:n), etc.

    !     x = alpha     e-folds
    ! Background
    !     y(1:n) = phi                 dydx(1:n)=dphi/dalpha
    !     y(n+1:2n) = dphi/dalpha      dydx(n+1:2n)=d^2phi/dalpha^2

    ! Mode matrix ptb, psi_IJ
    !     y(2n+1:2n+n**2) = psi               dydx(2n+1:3n)=dpsi/dalpha
    !     y(2n+n**2+1:2n+2n**2) = dpsi/dalpha       dydx(3n+1:4n)=d^2psi/dalpha^2

    ! Tensors
    !     y(2n+2n**2+1) = v                  dydx(4n+1)=dv/dalpha
    !     y(2n+2n**2+2) = dv/dalpha          dydx(4n+2)=d^2v/dalpha^2

    ! --- u_zeta is the adiabatic mode ignoring coupings to other modes, used to compare with the full zeta perturbation
    !! --- full_zeta - u_zeta gives the super-horizon evolution 
    !     y(2n+2n**2+3) = u_zeta     dydx(4n+4)=d^2u_zeta/dalpha^2  
    !     y(2n+2n**2+4) = du_zeta/dalpha     dydx(4n+4)=d^2u_zeta/dalpha^2  

    ! Set aliases for indices for above
    index_ptb_y = 2*num_inflaton+1
    index_ptb_vel_y = 2*num_inflaton+1+num_inflaton**2
    index_tensor_y = 2*num_inflaton+2*num_inflaton**2+1
    index_uzeta_y = index_tensor_y + 2

    ! Make the powerspectrum array.
    if (allocated(powerspectrum_out%matrix)) then
      deallocate(powerspectrum_out%matrix)
    end if
    if (allocated(power_internal%matrix)) then
      deallocate(power_internal%matrix)
    end if
    allocate(power_internal%matrix(num_inflaton,num_inflaton))
    allocate(powerspectrum_out%matrix(num_inflaton,num_inflaton))


    !Evaluation scale
    k=kin*Mpc2Mpl
    powerspectrum_out%k=k


    !! start where k = 100 aH, deep in the horizon, ah = log(aH), k_start=1d2
    ah=LOG(k/k_start)
    i= locate(aharr(1:nactual_bg), ah)

    IF(i.eq.0.) THEN
       PRINT*,'MODPK: The background solution worked, but the k you requested', k,' is outside'
       PRINT*,'MODPK: the bounds of the background you solved for. Please reconsider'
       PRINT*,'MODPK: your phi_init and N_pivot combo.'

       !Override the stop.
       if (sampling_techn/=reg_samp) then
         pk_bad=bad_ic
         return
       end if


       PRINT*,'MODPK: QUITTING'
       write(*,e2_fmt) "log(k/k_start):", ah
       write(*,e2_fmt) "aharr(1):", aharr(1)
       STOP
    END IF

    ! nactual_bg here is set by the background evolution
    j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)

    !MULTIFIELD
    CALL array_polint(aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
    CALL array_polint(aharr(j:j+4), dphiarr(:,j:j+4), ah, dp_ik, delphi)
    !END MULTIFIELD

    CALL polint(aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
    CALL polint(aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)

    a_ik=EXP(alpha_ik)*a_init

    x1=alpha_ik


    IF(x1.le.0.) THEN
       PRINT*,'MODPK: The phi_init you specified is too small to give'
       PRINT*,'MODPK: sufficient efolds of inflation. We cannot self-consistently'
       PRINT*,'MODPK: solve this for you. Please adjust phi_init and try again.'

       !Override the stop.
       if (sampling_techn/=reg_samp) then
         pk_bad=bad_ic
         return
       end if

       PRINT*,'MODPK: QUITTING'
       STOP
    END IF
    IF((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1 .OR. dalpha .GT. 0.1 .OR. dh .GT. 0.1) THEN
       PRINT*,'MODPK: The interpolation in SUBROUTINE EVOLVE has suspiciously large'
       PRINT*,'MODPK: errors. Your model smells fishy.'
       PRINT*,'MODPK: QUITTING'
       STOP
    ENDIF

    call set_ic(y)

    power_internal = powerspectrum_out

    !     Call the integrator
    !
    ode_underflow = .false.
    ode_ps_output = .true.
    ode_infl_end = .true.

    !DEBUG
    !save_steps = .true.
    save_steps = .false.

    pk_bad = 0

    !MULTIFIELD, need to evolve towards the end of inflation
    x2 = lna(nactual_bg) + 5.d0
    !END MULTIFIELD


    h1=0.1 !guessed start stepsize
    accuracy=1.0d-8 !has a big impact on the speed of the code
    hmin=1e-30_dp !minimum stepsize

    CALL odeint(y, x1, x2, accuracy, h1, hmin, derivs, qderivs, rkqs_c)
    nactual_mode = kount  ! update nactual after evolving the modes

    IF(.NOT. ode_underflow) THEN
      powerspectrum_out = power_internal
      powerspectrum_out%bundle_exp_scalar=field_bundle%exp_scalar
    ELSE
      powerspectrum_out%adiab= 0e0_dp
      powerspectrum_out%isocurv=0e0_dp
      powerspectrum_out%tensor=0e0_dp
      pk_bad=1
    ENDIF


    RETURN

    contains

      ! A "packed" identity vector analog of identity matrix.
      subroutine make_identity(identityvector)

        real(dp), dimension(:), intent(out) :: identityvector
        integer :: i, j

        identityvector=0d0

        do i=1,num_inflaton; do j=1, num_inflaton
          if (i==j) then
            identityvector((i-1)*num_inflaton+j)=1d0
          end if
        end do; end do

      end subroutine make_identity

      subroutine set_ic(y)

        ! Note that ah=log(aH) and overall scaled by sqrt(2k)

        complex(dp), dimension(:), intent(out) :: y
        real(dp), dimension(num_inflaton**2) :: identity

        ! Identity vector analog of identity matrix
        call make_identity(identity)

        ! Background - from previous evolution
        y(1:num_inflaton) = cmplx(p_ik)             !phi(x1)
        y(num_inflaton+1:2*num_inflaton) = cmplx(dp_ik) !Not in exact SR

        ! mode matrix - diagonalize, Bunch-Davies
        y(index_ptb_y:index_ptb_vel_y-1) = (1.d0, 0)*identity  !cmplx(1/sqrt(2*k))
        y(index_ptb_vel_y:index_tensor_y-1) = cmplx(0., -k/exp(ah))*identity

        ! tensors
        y(index_tensor_y) = (1.d0, 0) !cmplx(1/sqrt(2*k))
        y(index_tensor_y+1) = cmplx(0., -k/exp(ah))

        ! u_zeta
        y(index_uzeta_y) = (1.d0, 0) !cmplx(1/sqrt(2*k))
        y(index_uzeta_y+1) = cmplx(0., -k/exp(ah))

      end subroutine set_ic

  END SUBROUTINE evolve

END MODULE access_modpk
