MODULE access_modpk
  use modpkparams, only : dp
  USE camb_interface
  use modpk_output, only : out_opt
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
    USE modpk_observables
    USE background_evolution, ONLY : backgrnd
    USE potential, ONLY : initialphi
    use modpk_icsampling, only : bad_ic, sampling_techn, reg_samp
    IMPLICIT NONE
    real(dp) :: k,klo,khi


    !
    !     Solve the background equations
    !
    pk_bad = 0
    phi_init = initialphi(phi_init0)

    !NB: For eqen sampling, dphi_init set in trial_background
    CALL backgrnd

    !When scanning ICs, let some backgrnd errors be overridden
    if (sampling_techn/=reg_samp .and. pk_bad == bad_ic) then
      if (out_opt%modpkoutput) then
        print*, "--------------- BAD IC; RESTARTING -------------------"
      end if
      return
    end if

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
    USE modpk_observables
    USE potential, ONLY: pot,powerspectrum, dVdphi, getH, getdHdalpha, field_bundle, getEps, &
      pot, d2Vdphi2
    USE modpk_utils, ONLY : derivs, qderivs, rkqs_c
    use modpk_numerics, only : locate, polint, array_polint
    use modpk_icsampling, only : bad_ic, sampling_techn, reg_samp
    use modpk_qsf
    IMPLICIT NONE

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
    if (allocated(powerspectrum_out%phi_ij)) then
      deallocate(powerspectrum_out%phi_ij)
    end if
    if (allocated(power_internal%phi_ij)) then
      deallocate(power_internal%phi_ij)
    end if
    allocate(power_internal%phi_ij(num_inflaton,num_inflaton))
    allocate(powerspectrum_out%phi_ij(num_inflaton,num_inflaton))

    !Evaluation scale
    k=kin*Mpc2Mpl
    powerspectrum_out%k=k

    !When to start evaluating P(k), k<aH/eval_ps
    if (num_inflaton==1) then
      eval_ps = 5.0e2_dp
    else
      eval_ps = 1.0e0_dp
    end if
    useq_ps = 1.0e2_dp !When switch variables to q=\delta \phi (k<aH/useq_ps)

    !How far inside the horizon to set the modes' (Bunch-Davies) IC; k = k_start*aH
    call set_consistent_BD_scale(k_start)

    !! start where k = k_start* aH, deep in the horizon, ah = log(aH)
    ah=LOG(k/k_start)
    i= locate(log_aharr(1:nactual_bg), ah)

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
       write(*,e2_fmt) "log_aharr(1):", log_aharr(1)
       STOP
    END IF

    ! nactual_bg here is set by the background evolution
    j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)

    !MULTIFIELD
    CALL array_polint(log_aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
    CALL array_polint(log_aharr(j:j+4), dphiarr(:,j:j+4), ah, dp_ik, delphi)
    !END MULTIFIELD

    CALL polint(log_aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
    CALL polint(log_aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)

    if (sampling_techn == qsf_parametric) &
      call get_param_guess(ah)

    a_ik=exp(alpha_ik)*a_init

    x1=alpha_ik

    IF(x1.le.0.) THEN
       PRINT*,'MODPK: The phi_init you specified is too small to give'
       PRINT*,'MODPK: sufficient efolds of inflation. We cannot self-consistently'
       PRINT*,'MODPK: solve this for you. Please adjust phi_init and try again.'
       print*, "alpha=",x1

       print*, "dalpha", dalpha

       PRINT*,'MODPK: QUITTING'
       STOP
    END IF
    IF((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1 .OR. dalpha .GT. 0.1 .OR. dh .GT. 0.1) THEN
       PRINT*,'MODPK: The interpolation in SUBROUTINE EVOLVE has suspiciously large'
       PRINT*,'MODPK: errors. Your model smells fishy.'
       PRINT*,'MODPK: QUITTING'
       if ((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1)&
         print*, "MODPK: Error in dphi interpolation.",&
         (sqrt(dot_product(delphi,delphi))/num_inflaton)
       if (dalpha .GT. 0.1) print*, "MODPK: Error in alpha interpolation.", dalpha
       if (dh > 0.1) print*, "MODPK: Error in Hubble interpolation", dh

       STOP
    ENDIF

    call set_ic(y)

    power_internal = powerspectrum_out

    !     Call the integrator
    !
    ode_underflow = .false.
    ode_ps_output = .true.
    ode_infl_end = .true.

    save_steps = .false.

    pk_bad = 0

    !MULTIFIELD, need to evolve towards the end of inflation
    x2 = lna(nactual_bg) + 5.e0_dp
    !END MULTIFIELD


    !h1=0.1 !guessed start stepsize
    h1=1e-5 !guessed start stepsize

    !Some fast-roll cases need high accuracy; activate conditionally in odeint_c
    if (use_high_accuracy) then
      !has a big impact on the speed of the code
      accuracy=1.0e-7_dp
    else
      accuracy=1.0e-6_dp
    end if

    hmin=1e-30_dp !minimum stepsize

    CALL odeint(y, x1, x2, accuracy, h1, hmin, derivs, qderivs, rkqs_c)
    nactual_mode = kount  ! update nactual after evolving the modes

    if(.not. ode_underflow) then
      powerspectrum_out = power_internal
      powerspectrum_out%bundle_exp_scalar=field_bundle%exp_scalar
    else
      powerspectrum_out%adiab= 0e0_dp
      powerspectrum_out%isocurv=0e0_dp
      powerspectrum_out%tensor=0e0_dp
      pk_bad=1
    endif


    contains

      ! A "packed" identity vector analog of identity matrix.
      subroutine make_identity(identityvector)

        real(dp), dimension(:), intent(out) :: identityvector
        integer :: i, j

        identityvector=0e0_dp

        do i=1,num_inflaton; do j=1, num_inflaton
          if (i==j) then
            identityvector((i-1)*num_inflaton+j)=1e0_dp
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
        y(index_ptb_y:index_ptb_vel_y-1) = (1.e0_dp, 0)*identity  !cmplx(1/sqrt(2*k))
        y(index_ptb_vel_y:index_tensor_y-1) = cmplx(0., -k/exp(ah))*identity

        ! tensors
        y(index_tensor_y) = (1.e0_dp, 0) !cmplx(1/sqrt(2*k))
        y(index_tensor_y+1) = cmplx(0., -k/exp(ah))

        ! u_zeta
        y(index_uzeta_y) = (1.e0_dp, 0) !cmplx(1/sqrt(2*k))
        y(index_uzeta_y+1) = cmplx(0., -k/exp(ah))

      end subroutine set_ic

      !Find the scale at which we can set the Bunch-Davies initial state
      !self-consistently
      !Note that there might be a correction due to massive modes m_heavy>H
      subroutine set_consistent_BD_scale(k_start)
        real(dp), intent(out) :: k_start

        real(dp) :: ah
        real(dp) :: horiz_fract
        integer :: ah_index

        real(dp) :: tol

        real(dp) :: alpha_ik, dalpha
        real(dp) :: h_ik, dh, a_ik
        real(dp) :: eps, V, dV(num_inflaton), d2V(num_inflaton,num_inflaton)
        real(dp), dimension(num_inflaton) :: p_ik, dp_ik, delphi
        real(dp) :: check1
        real(dp), dimension(num_inflaton,num_inflaton) :: check2, check3, check4

        logical :: bd_consistent

        bd_consistent = .false.

        horiz_fract=1e0_dp

        do while (.not. bd_consistent)

          horiz_fract = horiz_fract*2.0e0_dp

          ah=LOG(k/horiz_fract)
          ah_index= locate(log_aharr(1:nactual_bg), ah)

          if (ah_index==0) then
            !The background isn't able to set this IC
            !Set the start scale so far inside, that it fails
            !when it returns from this function.
            k_start = 1e20_dp
            return
          end if

          j=min(max(ah_index-(4-1)/2,1),nactual_bg+1-4)
          call array_polint(log_aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
          call array_polint(log_aharr(j:j+4), dphiarr(:,j:j+4), ah, dp_ik, delphi)

          call polint(log_aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
          call polint(log_aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)

          if (sampling_techn == qsf_parametric) &
            call get_param_guess(ah)

          a_ik = exp(alpha_ik)*a_init
          eps = getEps(p_ik, dp_ik)

          V = pot(p_ik)
          dV = dVdphi(p_ik)
          d2V = d2Vdphi2(p_ik)

          !Check the corrections to the conformal-time mode equations.
          !When k is large, the IC should be a plane wave.
          !If the mass-matrix is diagonal and relevant, then could use Hankel
          !function solution
          check1 = abs((eps-2.0e0_dp)/(horiz_fract**2))
          check2 = abs( d2V/ (horiz_fract**2 * h_ik**2))

          do i=1, num_inflaton; do j=1, num_inflaton
            check3(i,j) = abs( (dp_ik(i)*dV(j) + &
              dp_ik(j)*dV(j))/(h_ik**2*horiz_fract**2))
          end do; end do

          do i=1, num_inflaton; do j=1, num_inflaton
            check4(i,j) = abs( (3.0e0_dp - eps)*dp_ik(i)*dp_ik(j)/horiz_fract**2)
          end do; end do

          tol = 1e-5_dp
          if   (check1 < tol  .and. &
            all(check2 < tol) .and. &
            all(check3 < tol) .and. &
            all(check4 < tol)) then

            bd_consistent = .true.

          end if

          k_start = horiz_fract

        end do

      end subroutine set_consistent_BD_scale

      !For numerical QSF trajectories, need to restart the guess for Newtonian
      !optimization
      subroutine get_param_guess(ah0)
        implicit none

        real(dp), intent(in) :: ah0


        !Get a guess for the initial param (should be pretty good)
        call polint(log_aharr(j:j+4),param_arr(j:j+4), &
          ah0, qsf_runref%param, dh)
        call qsf_runref%get_param(param=qsf_runref%param)

      end subroutine get_param_guess

  end subroutine evolve

end module access_modpk
