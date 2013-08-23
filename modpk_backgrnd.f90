MODULE background_evolution
  USE camb_interface
  USE modpkparams
  USE modpk_odeint
  USE ode_path
  USE powersp
  USE potential, ONLY : pot,getH, getHdot, dVdphi, getEps, d2Vdphi2
  USE modpk_utils, ONLY : locate, polint, bderivs, rkqs_r, array_polint
  IMPLICIT NONE
  PUBLIC :: backgrnd

CONTAINS

  SUBROUTINE backgrnd
    use modpk_icsampling, only : save_iso_N, N_iso_ref, phi_iso_N, &
      dphi_iso_N, sampling_techn, eqen_samp, bad_ic

    INTEGER*4 :: i,j, rescl_count

    real(dp) :: phi_init_trial(size(phi_init))
    real(dp) :: alpha_e,dalpha,V_end,dv,dh
    real(dp) :: a_end_inst
    real(dp) :: V_i, ph, alpha_pivot, bb(size(phi_init))
    real(dp) :: Np_last

    CHARACTER(16) :: fmt = '(a25,es10.3)'

    !MULTIFIELD
    CHARACTER(19) :: array_fmt
    CHARACTER(len=5) :: ci

    real(dp) :: alpha_iso_N

    write(ci, '(I5)'), size(phi_init)
    ci = adjustl(ci)
    array_fmt = '(a25,'//trim(ci)//'es10.3)'
    !END MULTIFIELD

    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dx
    !     y(2)=dphi/dx          dydx(1)=d^2phi/dx^2

    if (.not. allocated(phi_init)) then
      PRINT*, 'MODPK: Please initialize phi_ini as an array!'
      stop
    end if

    if (modpkoutput) write(*, array_fmt) 'phi_init =', phi_init
    if (modpkoutput .and. sampling_techn==eqen_samp) &
      write(*, array_fmt) 'dphi_init =', dphi_init0

    if (size(phi_init) .eq. 1) then !! for single field

       phidot_sign = -(dVdphi(phi_init))/ABS(dVdphi(phi_init))
       !if < 0 field rolls from large to small
       !if > 0, the field rolls from small to large

       !Check whether phi_infl_end is in the correct direction from phi_init
       IF(.NOT.(slowroll_infl_end)) THEN
          IF (phidot_sign(1) .GT.0 .AND. phi_init(1) .gt. phi_infl_end(1)) THEN
             PRINT*, 'MODPK: Initial phi is smaller than final phi.'
             PRINT*, 'MODPK: Please check your initial conditions'
             PRINT*, 'MODPK: QUITTING'
             STOP
          ENDIF
          IF (phidot_sign(1) .LT.0 .AND. phi_init(1) .lt. phi_infl_end(1)) THEN
             PRINT*, 'MODPK: Initial phi is larger than final phi.'
             PRINT*, 'MODPK: Please check your initial conditions'
             PRINT*, 'MODPK: QUITTING'
             STOP
          ENDIF

       ENDIF

       !Try the background evolution for a user-defined phi_init and rescale if necessary
       rescl_count=0
       rescale_factor = 0.1d0

       DO
          phi_init_trial=phi_init
          CALL trial_background(phi_init_trial, alpha_e, V_end)
          !! no rescaling implemented for multifield, so the code will exit here for multifield
          IF ((pk_bad/=0) .OR. phi_init_trial(1).EQ.phi_init(1)) EXIT
          rescl_count=rescl_count+1
          IF (rescl_count .EQ. 50) THEN
             pk_bad=2
             PRINT*,'MODPK: phi_init rescaling did not work after 50 tries.'
             EXIT
          END IF
          if (modpkoutput) then
             PRINT*,'MODPK: phi_init was inconsistent. Rescaling:'
             WRITE(*,fmt) ' phi_init =', phi_init_trial
          end if
          phi_init = phi_init_trial
       END DO
    !MULTIFIELD
    else
       CALL trial_background(phi_init, alpha_e, V_end)
       if (pk_bad .ne. 0) then
          print*, 'MODPK: pk_bad = ', pk_bad
          stop
       end if
    end if
    !END MULTIFIELD

    IF(pk_bad==0) THEN
       !Matching condition
       V_i=pot(phi_init)

       !MULTIFIELD
       if ((alpha_e - lna(1)) .lt. N_pivot) then
          if (modpkoutput) write(*, *) 'MODPK: Not enough efolds obtained. Please adjust your initial value'
          print*, "MODPK: alpha_e - lna(1) =", alpha_e - lna(1),"< ",N_pivot
          pk_bad = 6
          return

       !DEBUG
       !else if ((alpha_e - lna(1)) .gt. 2*N_pivot) then
       !else if ((alpha_e - lna(1)) .gt. 20*N_pivot) then
       else if ((alpha_e - lna(1)) .gt. Nefold_max) then

          if (modpkoutput) write(*, *) 'MODPK: Too many efolds obtained. Please rescale your initial value'
          pk_bad = 6
          !return
       end if
       !MULTIFIELD

       IF (instreheat) THEN
          a_init = EXP(-71.1d0 - alpha_e + LOG(V_i/V_end)/4.d0 + LOG((M_Pl**4)/V_i)/4.d0)
          Np_last = 0.5d0*N_pivot
          do while (abs(N_pivot-Np_last)>0.01d0)  ! determine N_pivot iteratively
             Np_last = N_pivot
             alpha_pivot = alpha_e-N_pivot

             i=locate(lna(1:nactual_bg),alpha_pivot)
             j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
             CALL array_polint(lna(j:j+4), phiarr(:, j:j+4), alpha_pivot, phi_pivot, bb)
             CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), alpha_pivot, dphi_pivot, bb)
             N_pivot = -71.1d0-log(V_end/(M_Pl**4))/4.d0-log(k_pivot*Mpc2Mpl/H_pivot)
          end do
          a_end = EXP(alpha_e)*a_init
       ELSE
          alpha_pivot = alpha_e-N_pivot

          !DEBUG
          if (save_iso_N) then
            alpha_iso_N = alpha_e - N_iso_ref
            if (alpha_iso_N<0e0_dp) then
              print*, "alpha_iso_N=",alpha_iso_N,"<0"
              pk_bad=bad_ic
              return
            end if
          end if

          i=locate(lna(1:nactual_bg),alpha_pivot)
          j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
          CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
          CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), alpha_pivot, phi_pivot, bb)
          CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), alpha_pivot, dphi_pivot, bb)

          !DEBUG
          if (save_iso_N) then
            if (.not. allocated(phi_iso_N)) then
              print*, "phi_iso_N not allocated..."
              stop
            endif

            i=locate(lna(1:nactual_bg),alpha_iso_N)
            j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
            CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), &
              alpha_iso_N, phi_iso_N, bb)
            CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), &
              alpha_iso_N, dphi_iso_N, bb)
            !DEBUG
            !print*, "fields 0", phi_init0
            !print*, "fields N", phi_iso_N
            !print*, "fields piv", phi_pivot
            !print*, "vels 0", dphi_init0
            !print*, "vels N", dphi_iso_N
            !print*, "vels piv", dphi_pivot
            !stop
          end if

          a_init=k_pivot*Mpc2Mpl/H_pivot/EXP(alpha_pivot)
          a_end=EXP(alpha_e)*a_init
          a_end_inst=EXP(-71.1d0+LOG(V_i/V_end)/4.d0+LOG((M_Pl**4)/V_i)/4.d0)
          IF (a_end .GT. a_end_inst) THEN
             PRINT*,'MODPK: inflation ends too late with this N_pivot.'
             pk_bad=3
             RETURN
          ENDIF
       END IF

       a_pivot = EXP(alpha_pivot)*a_init
       N_tot = alpha_e - lna(1)

       if (modpkoutput) then
          WRITE(*, fmt) ' H_pivot =', H_pivot
          WRITE(*, fmt) ' a_end = ', a_end
          WRITE(*, array_fmt) ' phi_pivot =', phi_pivot
          WRITE(*, fmt) ' a_pivot =', a_pivot
          WRITE(*, fmt) ' N_pivot =', N_pivot
          WRITE(*, fmt) ' N_tot =', N_tot
       end if

       DO i=1,nactual_bg
          aharr(i)=LOG(a_init*EXP(lna(i))*hubarr(i))
       END DO
    ENDIF

    RETURN
  END SUBROUTINE backgrnd

  SUBROUTINE trial_background(phi_init_trial, alpha_e, V_end)
    use modpk_icsampling, only : sampling_techn, eqen_samp

    INTEGER*4 :: i,j
    INTEGER*4, PARAMETER :: BNVAR=2

    !! MULTIFIELD
    real(dp), DIMENSION(:), INTENT(INOUT) :: phi_init_trial
    real(dp), DIMENSION(:) :: y(BNVAR*size(phi_init_trial))
    !! END MUTLTIFIELD

    real(dp) :: accuracy, h1, hmin, x1, x2
    real(dp) :: alpha_e, dalpha, V_end, dv, ep
    real(dp) :: ph, alpha_pivot, aa(size(phi_init_trial)), bb(size(phi_init_trial))
    real(dp) :: vv(nsteps) !!epsarr(nsteps),

    real(dp), DIMENSION(:) :: Vp(num_inflaton)
    real(dp) :: Vz, dotphi, thetaN, grad_V

    CHARACTER(16) :: fmt = '(a25,es10.3)'

    !MULTIFIELD
    CHARACTER(19) :: array_fmt
    CHARACTER(len=5) :: ci

    write(ci, '(I5)'), size(phi_init_trial)
    ci = adjustl(ci)
    array_fmt = '(a25, '//trim(ci)//'es10.3)'
    !END MULTIFIELD

    !MULTIFIELD
    h_init=SQRT(pot(phi_init_trial)/(6.d0*M_Pl**2) * &
         (1.+SQRT(1.+2./3.* M_Pl**2 * dot_product(dVdphi(phi_init_trial), dVdphi(phi_init_trial)) / pot(phi_init_trial)**2.)))
    !END MULTIFIELD

    x1=0.0 !starting value
    x2=Nefold_max !ending value

    !MULTIFIELD
    !Set the ICS
    y(1 : size(y)/2) = phi_init_trial  !phi(x1)

    if (sampling_techn==eqen_samp) then
      !Set w/equal energy, not necess close to SR
      y(size(y)/2+1 : (size(y))) = dphi_init0

    else
      !dphi/dalpha(x1) slowroll approx
      y(size(y)/2+1 : (size(y))) = &
        -dVdphi(phi_init_trial)/3./h_init/h_init
      ![ JF ] DEBUG
      print*, "-------"
      Print*, "dy",  y(size(y)/2+1 : (size(y))) 
      Print*, "potential", pot(phi_init_trial)
      Print*, "first deriv of V", dVdphi(phi_init_trial)
      Print*, "second deriv of V", d2Vdphi2(phi_init_trial)
      !DEBUG
      y(size(y)/2+1 : (size(y))) = 0e0_dp
      !stop
    end if
    !END MULTIFIELD

    !Call the integrator
    ode_underflow = .FALSE.
    ode_ps_output = .FALSE.
    ode_infl_end = .TRUE.
    save_steps = .true.
    pk_bad = 0

    !MULTIFIELD
    IF(getEps(y(1:size(y)/2),y(size(y)/2+1:size(y))) .GT. 1.) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF
    if (modpkoutput) write(*,'(a25, L2)') 'slowroll start =', slowroll_start
    !END MULTIFIELD

    !guessed start stepsize
    if (potential_choice.eq.6) then
       h1 = 0.001
    else
       h1 = 0.1
    end if
    dxsav=1.d-7
    accuracy=1.0d-7
    hmin=0.0 !minimum stepsize
    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs_r)

    IF(.NOT. ode_underflow) THEN
       lna(1:kount)=xp(1:kount)
       phiarr(:,1:kount)=yp(1:size(y)/2, 1:kount)
       dphiarr(:,1:kount)=yp(size(y)/2+1:size(y),1:kount)
       !MULTIFIELD
       DO i=1,kount
          vv(i) = pot(phiarr(:,i))
          hubarr(i) = getH(phiarr(:, i),dphiarr(:,i))
          epsarr(i) = getEps(phiarr(:,i),dphiarr(:,i))
          sigma_arr(i) = sqrt(dot_product((phiarr(:,i)-phi_init),(phiarr(:,i)-phi_init)))

          !! compute d\theta/dN_e
          dotphi = sqrt(dot_product(dphiarr(:,i), dphiarr(:,i)))
          Vp = dVdphi(phiarr(:,i))
          Vz = dot_product(Vp, phiarr(:,i))/dotphi
          grad_V = sqrt(dot_product(Vp, Vp))
          dtheta_dN = sqrt((grad_V + Vz)*(grad_V - Vz))/(dotphi*hubarr(i)**2)
       END DO
       !END MULTIFIELD
       !
       !     Determine the parameters needed for converting k(Mpc^-1) to K
       !
       nactual_bg=kount
       IF(slowroll_infl_end) THEN
          ep=1.d0
          i=locate(epsarr(1:kount),ep)
          j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
          CALL polint(epsarr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
          CALL polint(epsarr(j:j+4), vv(j:j+4), ep, V_end, dv)
          !MULTIFIELD
          CALL array_polint(epsarr(j:j+4), phiarr(:, j:j+4), ep, phi_infl_end, bb)
       ELSE
          if (size(phi_init) .eq. 1) then
             ep = abs(phi_init(1) - phi_infl_end(1))
             i = locate(sigma_arr(1:kount), ep)
             j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL polint(sigma_arr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
             V_end = pot(phi_infl_end)
          else
             ep = delsigma
             i = locate(sigma_arr(1:kount), ep)
             j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL polint(sigma_arr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
             CALL polint(sigma_arr(j:j+4), vv(j:j+4), ep, V_end, dv)
             CALL array_polint(sigma_arr(j:j+4), phiarr(:, j:j+4), ep, phi_infl_end, bb)
          end if
          !END MULTIFIELD
       ENDIF

       if (modpkoutput) WRITE(*,array_fmt) ' phi_end =', phi_infl_end

       IF (instreheat) THEN
          !Set a plausible pivot (ignoring logarithmic V_k and V_end terms) for
          !the smallest k likely to be requested by CAMB.
          N_pivot=-LOG10(1.d-6)+60.
       ENDIF

       !MULTIFIELD
       !
       !Rescaling the initial value only implemented for single field cases.
       !For multi-field, one may need to specify whether to rescale a sigle compoent field,
       !or to homogeneously rescale all the fields. It is better to work on a case by case basis.
       !

       IF (size(phi_init) .eq. 1) THEN
          IF(alpha_e.LT.(N_pivot+20.)) THEN
             IF ((potential_choice.eq.6).and.(vparams(1,1)<-2.d0)) THEN
                phi_init_trial = phi_init*0.9d0
             ELSE
                phi_init_trial = phi_init+(phi_init-phi_infl_end)*rescale_factor
             ENDIF
             RETURN
          END IF

          IF (alpha_e.GT.N_pivot+55) THEN
             ep = alpha_e-(N_pivot+50)
             i = locate(lna(1:kount),ep)
             j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), ep, aa, bb)
             phi_init_trial = aa
             RETURN
          END IF
       ENDIF
    ELSE
       pk_bad=4
    ENDIF

    !END MULTIFIELD
    RETURN

  END SUBROUTINE trial_background

  SUBROUTINE backgrnd_efold

    INTEGER*4, PARAMETER :: BNVAR=2
    real(dp), DIMENSION(:) :: y(BNVAR*num_inflaton)
    real(dp) :: accuracy, h1, hmin, x1, x2


    h_init=SQRT(pot(phi_init)/(6.d0*M_Pl**2) * &
         (1.+SQRT(1.+2./3.* M_Pl**2 * dot_product(dVdphi(phi_init), dVdphi(phi_init)) / pot(phi_init)**2.)))

    x1=0.0 !starting value
    x2=Nefold_max !ending value

    y(1 : size(y)/2) = phi_init  !phi(x1)
    y(size(y)/2+1 : (size(y))) = -dVdphi(phi_init)/3./h_init/h_init !dphi/dalpha(x1) slowroll approx

    !Call the integrator
    ode_underflow = .FALSE.
    ode_ps_output = .FALSE.
    ode_infl_end = .TRUE.
    save_steps = .TRUE.
    pk_bad = 0

    IF(getEps(y(1:size(y)/2),y(size(y)/2+1:size(y))) .GT. 1.) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF
    if (modpkoutput) write(*,'(a25, L2)') 'slowroll start =', slowroll_start

    !guessed start stepsize
    if (potential_choice.eq.6) then
       h1 = 0.001
    else
       h1 = 0.1
    end if
    dxsav=1.d-7
    accuracy=1.0d-6
    hmin=0.0 !minimum stepsize
    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs_r)

    IF(.NOT. ode_underflow) THEN
       lna(1:kount)=xp(1:kount)
       nactual_bg=kount
    ELSE
       pk_bad = 4
    END IF

    N_tot = lna(nactual_bg) - lna(1)

  END SUBROUTINE backgrnd_efold


END MODULE BACKGROUND_EVOLUTION
