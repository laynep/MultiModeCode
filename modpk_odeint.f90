MODULE modpk_odeint
  use modpkparams, only : dp
  use camb_interface, only : pk_bad
  use modpk_icsampling, only : ic_sampling, bad_ic, ic_flags
  use dvode_f90_m, only : vode_opts, set_normal_opts, dvode_f90, get_stats, &
    set_intermediate_opts
  use modpk_io, only : out_opt
  use csv_file, only : csv_write
  implicit none

  interface odeint
     module procedure odeint_r
     module procedure odeint_c
  end interface

  public :: odeint

contains

  subroutine odeint_r(ystart,x1,x2,eps,h1,hmin,derivs,rkqs_r)
    use ode_path
    use internals
    use modpk_observables
    use modpkparams
    use potential
    use modpk_utils, only : reallocate_rv, reallocate_rm, bderivs_dvode, &
      jacobian_background_DVODE
    use modpk_qsf

    implicit none
    real(dp), DIMENSION(:), INTENT(INOUT) :: ystart
    real(dp), INTENT(IN) :: x1,x2,eps,h1,hmin
    !MULTIFIELD
    real(dp), DIMENSION(num_inflaton) :: phi, dphi
    !END MULTIFIELD

    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         use modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         real(dp), DIMENSION(:), INTENT(IN) :: y
         real(dp), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs

       SUBROUTINE rkqs_r(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
         use modpkparams
         IMPLICIT NONE
         real(dp), DIMENSION(:), INTENT(INOUT) :: y
         real(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
         real(dp), INTENT(INOUT) :: x
         real(dp), INTENT(IN) :: htry,eps
         real(dp), INTENT(OUT) :: hdid,hnext
         INTERFACE
            SUBROUTINE derivs(x,y,dydx)
              use modpkparams
              IMPLICIT NONE
              real(dp), INTENT(IN) :: x
              real(dp), DIMENSION(:), INTENT(IN) :: y
              real(dp), DIMENSION(:), INTENT(OUT) :: dydx
            END SUBROUTINE derivs
         END INTERFACE
       END SUBROUTINE rkqs_r
    END INTERFACE

    real(dp), PARAMETER :: TINY=1.0e-30_dp
    INTEGER, PARAMETER :: MAXSTP=nsteps
    INTEGER*4 :: nstp,i
    real(dp) :: h,hdid,hnext,x,xsav
    real(dp), DIMENSION(SIZE(ystart)) :: dydx, y, yscal
    real(dp) :: z, scalefac

    real(dp) :: infl_efolds, infl_efolds_start
    logical :: infl_checking
    logical :: leave

    !For DVODE integrator
    integer :: neq, istats(31)
    integer :: itask, istate
    real(dp) :: rtol, rstats(22), nefold_out, dN_step
    real(dp), dimension(:), allocatable :: atol
    type (vode_opts) :: ode_integrator_opt

    !Inits for checking whether in
    !extended inflation period (for IC scan)
    infl_checking = .false.
    infl_efolds_start = 0e0_dp

    ode_underflow=.FALSE.
    infl_ended=.FALSE.
    x=x1
    h=SIGN(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    y(:)=ystart(:)
    NULLIFY(xp,yp)
    IF (save_steps) THEN
       xsav=x-2.0e0_dp*dxsav
       ALLOCATE(xp(256))
       if (ic_sampling == ic_flags%qsf_parametric) ALLOCATE(param_p(256))
       ALLOCATE(yp(SIZE(ystart),SIZE(xp)))
    END IF

    if (tech_opt%use_dvode_integrator) then
      !Options for first call to dvode integrator
      call initialize_dvode()
    end if

    DO nstp=1,MAXSTP

       if (any(isnan(y))) then
         print*, "MODPK in odeint_r"
         print*, "MODPK: y has a NaN value."
         print*, "E-fold",x
         print*, "nstp",nstp
         print*, "y", y
         stop
       end if

       !Calc bundle exp_scalar by integrating tr(d2Vdphi2)
       if (nstp==1) then
         field_bundle%N=0e0_dp
         field_bundle%dlogThetadN=0e0_dp
         field_bundle%exp_scalar=1e0_dp
       end if
       call field_bundle%calc_exp_scalar(y(1:num_inflaton),x)

       !Record the background trajectory
       if (out_opt%save_traj) call print_traj()

       CALL derivs(x,y,dydx)
       !If get bad deriv, then override this error when IC sampling
       if (pk_bad==bad_ic) return

       IF (save_steps .AND. (ABS(x-xsav) > ABS(dxsav))) &
            CALL save_a_step

      if (tech_opt%use_dvode_integrator) then

        if (tech_opt%use_analytical_jacobian) then
          call dvode_f90(bderivs_dvode,neq,y,x,nefold_out, &
            itask,istate,ode_integrator_opt,J_FCN=jacobian_background_DVODE)
        else
          call dvode_f90(bderivs_dvode,neq,y,x,nefold_out, &
            itask,istate,ode_integrator_opt)
        end if
        call get_stats(rstats,istats)

        if (istate<0) then
          print*, "MODPK in dvode_f90 istate=", istate
          stop
        end if

      else

        yscal(:)=ABS(y(:))+ABS(h*dydx(:))+TINY

        IF ((x+h-x2)*(x+h-x1) > 0.0) h = x2 - x

        CALL rkqs_r(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

        IF (hdid == h) THEN
           nok=nok+1
        ELSE
           nbad=nbad+1
        END IF

      end if

      call check_for_eternal_inflation()

      !MULTIFIELD
      phi = y(1:num_inflaton)
      dphi = y(num_inflaton+1 : 2*num_inflaton)

      call check_inflation_started_properly()

      !END MULTIFIELD

      call check_inflation_ended_properly(leave)
      if (leave) then
        !Record the background trajectory
        if (out_opt%save_traj) call print_traj()
        return
      end if

      if (abs(x2-x)<1e-10) then
        print*, "reached end of N-integration...."
        print*, "y=",y
        print*, "Efolds=",x
        stop
      end if

      IF (ode_underflow) RETURN

      if ( .not. tech_opt%use_dvode_integrator) then
        IF (ABS(hnext) < hmin) THEN
           write(*,*) 'stepsize smaller than minimum in odeint'
           STOP
        END IF
      end if

      !Set up next N-step
      if (tech_opt%use_dvode_integrator) then
        if (itask/=2) nefold_out = min(x + dN_step, x2)

        !Increase accuracy requirements when not in SR
        if (tech_opt%accuracy_setting>0) then
          if (getEps(phi,dphi)>0.2e0_dp) then
            rtol=1e-12_dp
            atol=1e-12_dp
            istate=3
          end if
        end if

      else
        h=hnext
      end if

    END DO

    !It got to the end without going through enough inflation to even be called
    !"slowroll_start"
    if (getEps(y(1:num_inflaton),y(num_inflaton+1:2*num_inflaton))>1.0e0_dp) then
      pk_bad = bad_ic
      print*, "MODPK: N-integration finished with eps>1.0 and"
      print*, "MODPK: without inflating or only transient periods of inflation."
      return

    else

      PRINT*, 'too many steps in odeint_r', nstp, MAXSTP
      PRINT*, "E-fold", x
      print*, "Step size", h
      print*, "epsilon=", getEps(y(1:num_inflaton),y(num_inflaton+1:2*num_inflaton))
      print*, "V=", pot(y(1:num_inflaton))
      PRINT*, "y=", y
      ode_underflow=.TRUE.
    end if

  CONTAINS

    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
    SUBROUTINE save_a_step
      kount=kount+1
      IF (kount > SIZE(xp)) THEN
         xp=>reallocate_rv(xp,2*SIZE(xp))
         if (ic_sampling==ic_flags%qsf_parametric) &
           param_p=>reallocate_rv(param_p,2*SIZE(xp))
         yp=>reallocate_rm(yp,SIZE(yp,1), SIZE(xp))
      END IF
      xp(kount)=x
      yp(:,kount)=y(:)
      xsav=x

      !For numerical QSF trajs
      if (ic_sampling==ic_flags%qsf_parametric) then
        param_p(kount) = qsf_runref%param
      end if

    END SUBROUTINE save_a_step

    subroutine initialize_dvode()

      neq = size(y)

      if (allocated(atol)) deallocate(atol)
      allocate(atol(neq))

      !Relative tolerance
      !Absolute tolerance
      if (tech_opt%accuracy_setting==2) then
        rtol = 1.e-10_dp
        atol = 1.0e-14_dp
      else if (tech_opt%accuracy_setting==1) then
        rtol = 1.e-6_dp
        atol = 1.0e-6_dp
      else if (tech_opt%accuracy_setting==0) then
        rtol = 1.e-5_dp
        atol = 1.0e-5_dp
      else if (tech_opt%accuracy_setting==-1) then
        rtol = tech_opt%dvode_rtol_back
        atol = tech_opt%dvode_atol_back(1:neq)
      else
        print*, "MODPK: initialize_dvode"
        print*, "MODPK: accuracy_setting =", tech_opt%accuracy_setting
        print*, "MODPK: is not supported."
        stop
      end if

      istate = 1 !Set =1 for 1st call to integrator

      itask  = 1 !Indicates normal usage, see dvode_f90_m.f90 for other values
      !itask  = 2 !Take only one time-step and output

      if (itask /=2) then
        !Integrate until nefold_out
        !dN_step = sign(0.001e0_dp,x2-x1)
        !dN_step = sign(0.01e0_dp,x2-x1)
        dN_step = sign(0.001e0_dp,x2-x1)
        nefold_out = x + dN_step
      else
        !Take only one step
        nefold_out = Nefold_max
      end if

      !Force initial step-size guess very small
      if (tech_opt%use_analytical_jacobian) then
        ode_integrator_opt = set_intermediate_opts(dense_j=.true.,&
          abserr_vector=atol,&
          relerr=rtol,&
          user_supplied_jacobian=.true., &
          mxstep=50000,&
          H0=1e-9_dp)
      else
        ode_integrator_opt = set_intermediate_opts(dense_j=.true.,&
          abserr_vector=atol,&
          relerr=rtol,&
          user_supplied_jacobian=.false., &
          mxstep=50000,&
          H0=1e-9_dp)
      end if

    end subroutine initialize_dvode

    subroutine print_traj()
      integer :: ii
      character(1024) :: cname
      logical :: adv

      !Write the column header
      if (out_opt%first_trajout) then

        !First column
        call csv_write(&
          out_opt%trajout,&
          'N', &
          advance=.false.)

        !Next num_inflaton columns
        do ii=1,num_inflaton
          write(cname, "(A3,I4.4)") "phi", ii
          call csv_write(&
            out_opt%trajout,&
            trim(cname), &
            advance=.false.)
        end do

        do ii=1,num_inflaton
          write(cname, "(A4,I4.4)") "dphi", ii
          call csv_write(&
            out_opt%trajout,&
            trim(cname), &
            advance=.false.)
        end do

        call csv_write(&
          out_opt%trajout,&
          (/'V', 'eps','H','eta'/), &
          advance=.false.)

        do ii=1,num_inflaton
          write(cname, "(A2,I4.4)") "dV", ii
          call csv_write(&
            out_opt%trajout,&
            trim(cname), &
            advance=.false.)
        end do

        !Advance here
        do ii=1,2*num_inflaton
          if (ii==2*num_inflaton) then
            adv=.true.
          else
            adv=.false.
          end if
          write(cname, "(A3,I4.4)") "d2V", ii
          call csv_write(&
            out_opt%trajout,&
            trim(cname), &
            advance=adv)
        end do

        out_opt%first_trajout = .false.
      end if

      !Write the trajectory
      write(out_opt%trajout,'(100E18.10)'),&
        x, &
        y(:), &
        pot(y(1:num_inflaton)),&
        getEps(y(1:num_inflaton),y(num_inflaton+1:2*num_inflaton)), &
        getH(y(1:num_inflaton),y(num_inflaton+1:2*num_inflaton)), &
        geteta(y(1:num_inflaton),y(num_inflaton+1:2*num_inflaton)), &
        dVdphi(y(1:num_inflaton)), &
        d2Vdphi2(y(1:num_inflaton))
    end subroutine print_traj

    subroutine check_for_eternal_inflation

       IF ((x-x2)*(x2-x1) > 0.0e0_dp) THEN

          PRINT*,'MODPK: This could be a model for which inflation does not end.'
          PRINT*,'MODPK: Either adjust phi_init or use slowroll_infl_end for a potential'
          PRINT*,'MODPK: for which inflation does not end by breakdown of slowroll.'
          PRINT*,'MODPK: QUITTING in odeint_r'
          WRITE(*, *) 'x, x1, x2 :', x, x1, x2
          WRITE(*,*) 'vparams: ', (vparams(i,:),i=1,size(vparams,1))
          IF (.NOT.instreheat) WRITE(*,*) 'N_pivot: ', N_pivot
          STOP
       END IF

     end subroutine check_for_eternal_inflation

     subroutine check_inflation_started_properly()

       IF(getEps(phi,dphi) .LT. 1 .AND. .NOT.(slowroll_start)) then

         if (ic_sampling==ic_flags%slowroll_samp .or. &
           ic_sampling==ic_flags%iso_N .or.&
           ic_sampling==ic_flags%reg_samp .or. &
           ic_sampling==ic_flags%qsf_random .or. &
           ic_sampling==ic_flags%qsf_parametric) then

           slowroll_start=.true.
         else
           !If scan ICs, say inflating iff eps<1 for "extended" period,
           !3 efolds - protects against transient inflation epochs, i.e.,
           !at traj turn-around or chance starting with dphi=0

           if (.not. infl_checking) then
             infl_checking = .true.
             infl_efolds_start = x
           else

             infl_efolds = x - infl_efolds_start
             if (infl_efolds > 3.0) then
               slowroll_start=.true.
             end if

           end if
         end if
       else if (infl_checking) then
         infl_checking=.false.
       endif

     end subroutine check_inflation_started_properly

     subroutine check_inflation_ended_properly(leave)
       logical, intent(inout) :: leave

       leave = .false.

       IF(ode_infl_end) THEN
          IF (slowroll_infl_end) THEN
             IF(getEps(phi, dphi) .GT. 1 .AND. slowroll_start) THEN
                infl_ended = .TRUE.
                ystart(:) = y(:)
                IF (save_steps) CALL save_a_step

                leave = .true.
                RETURN
             ENDIF
          ELSE
             IF(getEps(phi, dphi) .GT. 1 .AND. slowroll_start) THEN
                PRINT*,'MODPK: You asked for a no-slowroll-breakdown model, but inflation'
                PRINT*,'MODPK: already ended via slowroll violation before your phi_end was'
                PRINT*,'MODPK: reached. Please take another look at your inputs.'
                PRINT*,'MODPK: QUITTING'
                PRINT*,'EPSILON =', getEps(phi, dphi)
                STOP
             ENDIF

             !MULTIFIELD
             IF (size(phi) .eq. 1) THEN
                IF (phidot_sign(1).GT.0..AND.(phi(1).GT.(phi_infl_end(1)+0.1))) THEN
                   infl_ended=.TRUE.
                   ystart(:)=y(:)
                   IF (save_steps) CALL save_a_step

                   leave=.true.
                   RETURN
                ENDIF
                IF (phidot_sign(1).LT.0..AND.(phi(1).LT.(phi_infl_end(1)-0.1))) THEN
                   infl_ended=.TRUE.
                   ystart(:)=y(:)
                   IF (save_steps) CALL save_a_step

                   leave=.true.
                   RETURN
                ENDIF
             ELSE
               ! for multifield, determine the total field distance travelled

               if(check_stop_when_not_slowroll_infl_end(phi,dphi)) then
                 infl_ended = .true.
                 ystart(:) = y(:)
                 IF (save_steps) CALL save_a_step

                 leave=.true.
                 RETURN
               end if

             END IF
             !END MULTIFIELD
          ENDIF
       ENDIF

     end subroutine check_inflation_ended_properly

  END SUBROUTINE odeint_r

  ! Only called for ptb mode eqns
  SUBROUTINE odeint_c(ystart, x1, x2, eps, h1, hmin, derivs, qderivs, rkqs_c)
    USE ode_path
    USE internals
    USE modpk_observables
    USE modpkparams
    USE potential, only : tensorpower, getH, getEps, zpower,&
      powerspectrum
    use modpk_qsf, only : turning_choice
    use modpk_utils, only : reallocate_rv, reallocate_rm, mode_derivs_dvode, &
      qderivs_dvode, jacobian_psi_modes_DVODE

    IMPLICIT NONE
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(INOUT) :: ystart
    real(dp), INTENT(IN) :: x1,x2,eps,h1,hmin
    real(dp) :: eps_adjust

    real(dp), DIMENSION(num_inflaton) :: phi, delphi  ! the classical field phi and dphi are real
    COMPLEX(KIND=DP), DIMENSION(size(ystart)) :: ytmp
    LOGICAL :: use_q, compute_zpower

    INTERFACE
       SUBROUTINE derivs(x, y, dydx)
         USE modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs

       SUBROUTINE qderivs(x, y, dydx)
         USE modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE qderivs

       SUBROUTINE rkqs_c(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
         USE modpkparams
         IMPLICIT NONE
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(INOUT) :: y
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
         real(dp), INTENT(INOUT) :: x
         real(dp), INTENT(IN) :: htry,eps
         real(dp), INTENT(OUT) :: hdid,hnext
         INTERFACE
            SUBROUTINE derivs(x, y, dydx)
              USE modpkparams
              IMPLICIT NONE
              real(dp), INTENT(IN) :: x
              COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
              COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: dydx
            END SUBROUTINE derivs
         END INTERFACE
       END SUBROUTINE rkqs_c
    END INTERFACE

    real(dp), PARAMETER :: TINY=1.0e-40_dp
    INTEGER, PARAMETER :: MAXSTP=nsteps
    INTEGER*4 :: nstp,i
    real(dp) :: h,hdid,hnext,x,xsav
    COMPLEX(KIND=DP), DIMENSION(size(ystart)) :: dydx, y, yscal
    real(dp) :: scalefac, hubble, a_switch, dotphi

    complex(dp), dimension(num_inflaton**2) :: psi, dpsi
    complex(dp), dimension(num_inflaton**2) :: qij, dqij

    real(dp) :: nk_sum, Nprime(num_inflaton), Nprimeprime(num_inflaton,num_inflaton)
    integer :: ii, jj, kk

    !For DVODE integrator
    real(dp), dimension(size(ystart)*2) :: yreal, dyrealdx

    integer :: neq, istats(31)
    integer :: itask, istate
    real(dp) :: rtol, rstats(22), nefold_out, dN_step
    real(dp), dimension(:), allocatable :: atol, atol_real, atol_compl
    type (vode_opts) :: ode_integrator_opt

    character(1024) :: cname


    !DEBUG
    !for two-knot QSF trajectory
    real(dp) :: a_turn1, N_turn1, H_turn1, phi_test
    real(dp) :: a_turn2, N_turn2, H_turn2
    phi_test = real(ystart(1))
    a_turn1 = 0e0_dp
    N_turn1=0e0_dp
    H_turn1=0e0_dp
    a_turn2=0e0_dp
    N_turn2=0e0_dp
    H_turn2 = 0e0_dp

    ode_underflow=.FALSE.
    infl_ended=.FALSE.
    x=x1
    h=SIGN(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    y(:)=ystart(:)
    NULLIFY(xp,yp)

    IF (save_steps) THEN
       xsav=x-2.e0_dp*dxsav
       ALLOCATE(xp(256))
       !MULTIFIELD
       ALLOCATE(yp(2*SIZE(ystart),SIZE(xp)))  !! store real and imiganary seperately
       !END MULTIFIELD
    END IF

    use_q = .false.
    compute_zpower = .true.
    eps_adjust = eps

    if (tech_opt%use_dvode_integrator) then
      !Options for first call to dvode integrator
      call initialize_dvode_MODES()
    end if

   DO nstp=1,MAXSTP

       if (any(isnan(real(y))) .or. any(isnan(aimag(y)))) then
         print*, "MODPK in odeint_c"
         print*, "MODPK: y has a NaN value."
         print*, "E-fold",x
         print*, "nstp",nstp
         print*, "y", y
         stop
       end if

       IF (use_q) THEN
         ! super-h use Q
         CALL qderivs(x, y, dydx)
       ELSE
         ! sub-h use psi
         CALL derivs(x, y, dydx)
       END IF

       IF (save_steps .AND. (ABS(x-xsav) > ABS(dxsav))) &
         CALL save_a_step

       if (tech_opt%use_dvode_integrator) then

         !Cmplx --> real
         yreal(1:neq/2) = real(y)
         yreal(neq/2+1:neq) = aimag(y)

         if (use_q) then
           call dvode_f90(qderivs_dvode,neq,yreal,x,nefold_out, &
             itask,istate,ode_integrator_opt)
         else
           call dvode_f90(mode_derivs_dvode,neq,yreal,x,nefold_out, &
             itask,istate,ode_integrator_opt)
         end if
         call get_stats(rstats,istats)

         if (istate<0) then
           print*, "MODPK in dvode_f90 istate=", istate
           stop
         end if

         !Set complex y from real y's
         y = cmplx(yreal(1:neq/2),yreal(neq/2+1:neq), kind=dp)

       else

         ! for yscal, evaluate real and imaginary parts separately,
         ! and then assemble them into complex format

         !"Trick" to get constant fractional errors except very near
         !zero-crossings. (Numerical Recipes)
         yscal(:)=cmplx(ABS(real(y(:),kind=dp))+ABS(h*real(dydx(:),kind=dp))+TINY, &
           ABS(real(y(:)*(0,-1),kind=dp))+ABS(h*real(dydx(:)*(0,-1),kind=dp))+TINY)


         IF ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x

         IF (use_q) THEN
            CALL rkqs_c(y,dydx,x,h,eps_adjust,yscal,hdid,hnext,qderivs)
         ELSE
            CALL rkqs_c(y,dydx,x,h,eps_adjust,yscal,hdid,hnext,derivs)
         END IF

         IF (hdid == h) THEN
            nok=nok+1
         ELSE
            nbad=nbad+1
         END IF

       end if

       !MULTIFIELD
       phi = real(y(1:num_inflaton),kind=dp)
       delphi = real(y(num_inflaton+1 : 2*num_inflaton),kind=dp)
       dotphi = sqrt(dot_product(delphi, delphi))

       if (out_opt%modes) call print_modes()

       scalefac = a_init*exp(x)

       !Increase accuracy requirements when not in SR
       if (tech_opt%accuracy_setting>0) then
         if (getEps(phi,delphi)>0.2e0_dp) then
           eps_adjust=1e-10_dp
           if (tech_opt%accuracy_setting==2) then
             if (getEps(phi,delphi)>0.9e0_dp) then
               eps_adjust=1e-16_dp
             end if
           end if
         end if
       else
         eps_adjust = eps
       end if
       !END MULTIFIELD

       IF(getEps(phi,delphi) .LT. 1 .AND. .NOT.(slowroll_start)) slowroll_start=.true.

       IF(ode_ps_output) THEN

         !DEBUG
         if (potential_choice == 13 .and. &
           all(turning_choice == 5)) then
           call find_turning_scales()
         end if

         ! if k<aH/eval_ps, then k<<aH
         if(k .lt. a_init*exp(x)*getH(phi, delphi)/eval_ps) &
           call evaluate_powerspectra()

       END IF

       call check_inflation_ended_properly_MODES()

       call check_for_eternal_inflation_MODES()

       IF(ode_infl_end) THEN
          IF (infl_ended) THEN
            IF (use_q) THEN
               ytmp(:) = y(:)
               ! bckgrd
               ystart(1:2*num_inflaton) = y(1:2*num_inflaton)

               ! ptbs
               ystart(index_ptb_y:index_ptb_vel_y-1) = &
                 ytmp(index_ptb_y:index_ptb_vel_y-1)*scalefac/a_switch
               ystart(index_ptb_vel_y:index_tensor_y-1) = &
                 ytmp(index_ptb_vel_y:index_tensor_y-1)&
                 *scalefac/a_switch + ystart(index_ptb_y:index_ptb_vel_y-1)

               ! tensors
               ystart(index_tensor_y) =&
                 ytmp(index_tensor_y)*scalefac/a_switch
               ystart(index_tensor_y+1) =&
                 ytmp(index_tensor_y+1)*scalefac/a_switch&
                 + ystart(index_tensor_y)
            ELSE
               ystart(:) = y(:)
            END IF
            IF (save_steps) CALL save_a_step

            !For outputting field values
            if (out_opt%fields_end_infl) then
              !Make header column

              if (out_opt%first_fields_end_out) then

                !First column
                call csv_write(out_opt%fields_end_out,&
                  'k',&
                  advance=.false.)

                !Next num_inflaton columns
                do ii=1,num_inflaton
                  write(cname, "(A3,I4.4)") "phi_piv", ii
                  if (ii==num_inflaton) then
                    call csv_write(&
                      out_opt%fields_end_out,&
                      trim(cname), &
                      advance=.true.)
                  else
                    call csv_write(&
                      out_opt%fields_end_out,&
                      trim(cname), &
                      advance=.false.)
                  end if
                end do

                out_opt%first_fields_end_out = .false.

              end if

              !Write data
              call csv_write(out_opt%fields_end_out,&
                (/k, phi/),&
                advance=.true.)
            end if

            !!For outputting field values
            !if (out_opt%fields_end_infl) &
            !  write(out_opt%fields_end_out,'(500E28.20)') k, phi

            RETURN
          END IF
       ENDIF

       !switch to the Q variable for super-horizon evolution
       !only apply the switch on y(1:4*num_inflaton+2)
       IF (k .LT. a_init*exp(x)*getH(phi, delphi)/useq_ps &
         .and. (.not. use_q)) THEN
         call switch_to_qvar()
       end if

       IF (ode_underflow) RETURN

       if (.not. tech_opt%use_dvode_integrator) then
         IF (ABS(hnext) < hmin) THEN
            WRITE(*,*) 'stepsize smaller than minimum in odeint'
            stop
         end if
       end if

       !Set up next N-step
       if (tech_opt%use_dvode_integrator) then
         if (itask/=2) nefold_out = min(x + dN_step, x2)
       else
         h=hnext
       end if

    end do

    print*,'too many steps in odeint_c'
    print*,'N =', x
    print*, 'stepsize, h =', h
    print*, 'background, y =', y(1:num_inflaton)
    print*, 'accuracy =', eps_adjust, eps
    print*, "epsilon", getEps(phi,delphi)
    ode_underflow=.TRUE.

    stop

  contains

    subroutine print_modes()

      character(1024) :: cname
      integer :: ii

      !Make column headers
      if (out_opt%first_modeout) then
        out_opt%first_modeout = .false.

        !First column --- reals
        call csv_write(out_opt%modeout(1),&
          'N',&
          advance=.false.)
        call csv_write(out_opt%modeout(3),&
          'N',&
          advance=.false.)

        !Next num_inflaton columns
        do ii=index_ptb_y,index_ptb_vel_y-1
          write(cname, "(A8,I4.4)") "Re[mode]", ii
          if (ii==index_ptb_vel_y-1) then
            call csv_write(&
              out_opt%modeout(1),&
              trim(cname), &
              advance=.true.)
          else
            call csv_write(&
              out_opt%modeout(1),&
              trim(cname), &
              advance=.false.)
          end if

          write(cname, "(A9,I4.4)") "Re[qmode]", ii
          if (ii==index_ptb_vel_y-1) then
            call csv_write(&
              out_opt%modeout(3),&
              trim(cname), &
              advance=.true.)
          else
            call csv_write(&
              out_opt%modeout(3),&
              trim(cname), &
              advance=.false.)
          end if
        end do

        !First column --- imags
        call csv_write(out_opt%modeout(2),&
          'N',&
          advance=.false.)
        call csv_write(out_opt%modeout(4),&
          'N',&
          advance=.false.)

        !Next num_inflaton columns
        do ii=index_ptb_y,index_ptb_vel_y-1
          write(cname, "(A8,I4.4)") "Im[mode]", ii
          if (ii==index_ptb_vel_y-1) then
            call csv_write(&
              out_opt%modeout(2),&
              trim(cname), &
              advance=.true.)
          else
            call csv_write(&
              out_opt%modeout(2),&
              trim(cname), &
              advance=.false.)
          end if

          write(cname, "(A9,I4.4)") "Im[qmode]", ii
          if (ii==index_ptb_vel_y-1) then
            call csv_write(&
              out_opt%modeout(4),&
              trim(cname), &
              advance=.true.)
          else
            call csv_write(&
              out_opt%modeout(4),&
              trim(cname), &
              advance=.false.)
          end if
        end do

      end if

      if (.not. use_q) then
        write(out_opt%modeout(1),'(100E30.22)') &
          x - (n_tot - N_pivot), &
          real(y(index_ptb_y:index_ptb_vel_y-1))/sqrt(2*k)
        write(out_opt%modeout(2),'(100E30.22)') &
          x - (n_tot - N_pivot),&
          aimag(y(index_ptb_y:index_ptb_vel_y-1))/sqrt(2*k)
      else
        write(out_opt%modeout(3),'(100E30.22)') &
          x - (n_tot - N_pivot), &
          real(y(index_ptb_y:index_ptb_vel_y-1))/sqrt(2*k)
        write(out_opt%modeout(4),'(100E30.22)') &
          x - (n_tot - N_pivot),&
          aimag(y(index_ptb_y:index_ptb_vel_y-1))/sqrt(2*k)
      end if

    end subroutine print_modes

    SUBROUTINE save_a_step
      USE modpkparams
      IMPLICIT NONE

      COMPLEX(KIND=DP), DIMENSION(size(ystart)) :: ytmp
      COMPLEX(KIND=DP), DIMENSION(num_inflaton**2) :: ptb_tmp, dptb_tmp

      kount=kount+1
      IF (kount > SIZE(xp)) THEN
         xp=>reallocate_rv(xp,2*SIZE(xp))
         yp=>reallocate_rm(yp,SIZE(yp,1), SIZE(xp))
      END IF
      xp(kount) = x

      IF (use_q) THEN  ! convert from (a_switch*Q) to v
         ytmp(:) = y(:)
         ptb_tmp =ytmp(index_ptb_y:index_ptb_vel_y-1)
         dptb_tmp =ytmp(index_ptb_vel_y:index_tensor_y-1)

         ytmp(index_ptb_y:index_ptb_vel_y-1) = ptb_tmp*scalefac/a_switch
         ytmp(index_ptb_vel_y:index_tensor_y-1) = dptb_tmp*scalefac/a_switch + ptb_tmp

         ytmp(index_tensor_y) = &
           ytmp(index_tensor_y)*scalefac/a_switch
         ytmp(index_tensor_y+1) = &
           ytmp(index_tensor_y+1)*scalefac/a_switch + ytmp(index_tensor_y)
      END IF

      yp(1:size(yp,1)/2,kount) = real(ytmp(:),kind=dp)
      yp(size(yp,1)/2+1:size(yp,1),kount) = real(ytmp(:)*(0,-1),kind=dp)
      xsav=x

    END SUBROUTINE save_a_step
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.

    subroutine initialize_dvode_MODES()

      neq = 2*size(y) !Need reals for dvode, so 2* bc y is complex

      if (allocated(atol)) deallocate(atol)
      if (allocated(atol_real)) deallocate(atol_real)
      if (allocated(atol_compl)) deallocate(atol_compl)
      allocate(atol(neq))
      allocate(atol_real(neq/2))
      allocate(atol_compl(neq/2))

      if (tech_opt%accuracy_setting==-1) then

        rtol = tech_opt%dvode_rtol_modes
        atol(1:neq/2) = tech_opt%dvode_atol_modes_real(1:neq/2)
        atol(neq/2+1:neq) = tech_opt%dvode_atol_modes_imag(1:neq/2)

      else

        !DEBUG
        !print*, "testing accuracy settings........................"

        !Set rtol to 10^-(m+1) where m = # decimal places that are important
        !rtol = 1e-5_dp

        !Set atol_i where |y_i| is insignificant
        !atol(1:num_inflaton) = 1e-8_dp
        !atol(num_inflaton+1:2*num_inflaton) = 1e-6_dp
        !atol(neq/2+1:neq/2+num_inflaton) = 1e2_dp !Set to zero exactly in derivs
        !atol(neq/2+num_inflaton+1:neq/2+2*num_inflaton) = 1e2_dp !Set to zero exactly in derivs

        !atol(index_ptb_y:index_ptb_vel_y-1) = 1e-5_dp
        !atol(index_ptb_vel_y:index_tensor_y-1) = 1e-5_dp
        !atol(neq/2+index_ptb_y:neq/2+index_ptb_vel_y-1) = 1e-5_dp
        !atol(neq/2+index_ptb_vel_y:neq/2+index_tensor_y-1) = 1e-5_dp

        !atol(index_tensor_y:index_tensor_y+1) = 1.0e-8_dp
        !atol(index_uzeta_y:index_uzeta_y+1) = 1.0e-3_dp
        !atol(neq/2+index_tensor_y:neq/2+index_tensor_y+1) = 1e0_dp
        !atol(neq/2+index_uzeta_y:neq/2+index_uzeta_y+1) = 1e0_dp

        !atol((neq/2)+1:neq)=atol(1:neq/2)

        !DEBUG
        rtol = 1e-6_dp
        atol = 1e-14_dp

        print*, "testing accuracy..."



      end if

      itask = 1 !Indicates normal usage, see dvode_f90_m.f90 for other values
      !itask = 2
      istate = 1 !Set =1 for 1st call to integrator

      if (itask /=2) then
        !Integrate until nefold_out
        !dN_step = sign(0.01e0_dp,x2-x1)
        !dN_step = sign(0.005e0_dp,x2-x1)
        dN_step = sign(0.001e0_dp,x2-x1)
        nefold_out = x + dN_step
      else
        !Take only one step
        nefold_out = Nefold_max
      end if

      if (tech_opt%use_analytical_jacobian) then
        ode_integrator_opt = set_intermediate_opts(dense_j=.true., abserr_vector=atol,&
          relerr=rtol, user_supplied_jacobian=.true., mxstep=1000, &
          mxhnil=1)
      else
        ode_integrator_opt = set_intermediate_opts(dense_j=.true., abserr_vector=atol,&
          relerr=rtol, user_supplied_jacobian=.false., mxstep=10000, &
          mxhnil=1)
      end if

    end subroutine initialize_dvode_MODES

    subroutine check_for_eternal_inflation_MODES()

       IF ((x-x2)*(x2-x1) >= 0.0) THEN
          PRINT*,'MODPK: This could be a model for which inflation does not end.'
          PRINT*,'MODPK: Either adjust phi_init or use slowroll_infl_end for a potential'
          PRINT*,'MODPK: for which inflation does not end by breakdown of slowroll.'
          PRINT*,'MODPK: QUITTING in odeint_c'
          WRITE(*, *) 'vparams: ', (vparams(i,:),i=1, size(vparams,1))
          WRITE(*, *) 'x1, x, x2 :', x1, x, x2
          WRITE(*, *) 'phi_back :', real(y(1:num_inflaton))
          WRITE(*, *) 'epsilon :', geteps(real(y(1:num_inflaton)),&
            real(y(num_inflaton+1:2*num_inflaton)))
          IF (.NOT.instreheat) WRITE(*,*) 'N_pivot: ', N_pivot
          STOP
       END IF

    end subroutine check_for_eternal_inflation_MODES

    subroutine evaluate_powerspectra()

      !MULTIFIELD
      IF (use_q) THEN

        !Y's are in \bar{Q}=Q/a_switch
        qij = y(index_ptb_y:index_ptb_vel_y-1)/a_switch
        dqij = y(index_ptb_vel_y:index_tensor_y-1)/a_switch

        ! with isocurv calculation
        call powerspectrum(qij, dqij, phi, delphi, &
          scalefac, power_internal, using_q=.true.)


        power_internal%tensor=tensorpower(y(index_tensor_y) &
           *scalefac/a_switch, scalefac)
      ELSE

        psi = y(index_ptb_y:index_ptb_vel_y-1)
        dpsi = y(index_ptb_vel_y:index_tensor_y-1)

        call powerspectrum(psi, dpsi, phi, delphi, &
          scalefac, power_internal)

        power_internal%tensor=tensorpower(y(index_tensor_y), scalefac)
      END IF

      !Record spectrum
      if (out_opt%spectra) then
        call write_spectra()
      end if

      if (compute_zpower) then  !! compute only once upon horizon exit
         power_internal%powz = zpower(y(index_uzeta_y), dotphi, scalefac)
         compute_zpower = .false.
      end if

      ! for single field, end mode evolution when the mode is frozen out of the horizon
      ! for multifield, need to evolve the modes until the end of inflation to
      ! include superhorizon evolution
      IF (num_inflaton .EQ. 1) infl_ended = .TRUE.
      !END MULTIFIELD

    end subroutine evaluate_powerspectra

    subroutine write_spectra()

      !Write the column header
      if (out_opt%first_spectraout) then

        call csv_write(&
          out_opt%spectraout,&
          (/'N', 'P_ad', 'P_iso', 'P_ent', 'P_nad', 'P_tens'/), &
          advance=.true.)

        out_opt%first_spectraout = .false.
      end if

      call csv_write(out_opt%spectraout,&
        (/ x - (n_tot - N_pivot), &
          power_internal%adiab, &
          power_internal%isocurv,&
          power_internal%entropy, &
          power_internal%pnad, &
          power_internal%tensor/),&
        advance=.true.)

    end subroutine write_spectra

    subroutine check_inflation_ended_properly_MODES()

       IF(ode_infl_end) THEN
          IF (slowroll_infl_end) THEN
             IF(getEps(phi, delphi) .GT. 1 .AND. slowroll_start) infl_ended=.TRUE.
          ELSE
             IF(getEps(phi, delphi) .GT. 1 .AND. slowroll_start) THEN
                PRINT*,'MODPK: You asked for a no-slowroll-breakdown model, but inflation'
                PRINT*,'MODPK: already ended via slowroll violation before your phi_end was'
                PRINT*,'MODPK: reached. Please take another look at your inputs.'
                PRINT*,'MODPK: QUITTING'
                PRINT*,'EPSILON =', getEps(phi, delphi), 'phi =', phi
                STOP
             ENDIF

             !MULTIFIELD
             IF (SIZE(phi) .EQ. 1) THEN
                IF (phidot_sign(1).GT.0..AND.(phi(1).GT.(phi_infl_end(1)+0.1))) infl_ended=.TRUE.
                IF (phidot_sign(1).LT.0..AND.(phi(1).LT.(phi_infl_end(1)-0.1))) infl_ended=.TRUE.
             ELSE
               if (check_stop_when_not_slowroll_infl_end(phi,delphi)) infl_ended = .TRUE.
             END IF
             !END MULTIFIELD
          ENDIF
        end if

    end subroutine check_inflation_ended_properly_MODES

    subroutine switch_to_qvar()
          use_q = .TRUE.
          a_switch = scalefac
          !set intial condition in (Q*a_switch)
          ytmp(:) = y(:)

          y(index_ptb_y:index_ptb_vel_y-1) = ytmp(index_ptb_y:index_ptb_vel_y-1)
          y(index_ptb_vel_y:index_tensor_y-1) = &
            ytmp(index_ptb_vel_y:index_tensor_y-1) &
            - y(index_ptb_y:index_ptb_vel_y-1)

          y(index_tensor_y) = ytmp(index_tensor_y)
          y(index_tensor_y+1) = ytmp(index_tensor_y+1) &
            - y(index_tensor_y)

          if (tech_opt%use_dvode_integrator) then
            !Reset istate to let integrator know it's a new variable
            istate=1

!DEBUG
!print*, "Resetting absolute error tolerances"
!atol(index_ptb_y:index_ptb_vel_y-1) = 1e-4_dp
!atol(neq/2+index_ptb_y:neq/2+index_ptb_vel_y-1) = 1e-4_dp

!atol = 1e-14_dp

          end if


    end subroutine switch_to_qvar

    !DEBUG
    !Capture the first and second turn positions in N
    !for two-knot QSF trajectory
    subroutine find_turning_scales()
      real(dp) :: phi_knot

      phi_knot = vparams(4,2) + vparams(6,2)/ &
        tan(vparams(5,2))
      if ( phi(1) <= phi_knot .and. &
        phi_test > phi_knot) then

        a_turn1 = a_init*exp(x)
        N_turn1 = x
        H_turn1 = getH(phi,delphi)

        !DEBUG
        print*, "a_* H_* at turn1 =", a_turn1*H_turn1
        print*, "a_* m_H at turn1 =", a_turn1*sqrt(10.0e0_dp**vparams(2,2))
        print*, "k_* at turn1 =", a_turn1*H_turn1/Mpc2Mpl
        print*, "k_m at turn1 =", a_turn1*sqrt(10.0e0_dp**vparams(2,2))/Mpc2Mpl
        print*, "N at turn1 =", N_turn1
        print*, "H at turn1 =", H_turn1
        print*, "m_H/H at turn1 =",sqrt(10.0e0_dp**vparams(2,2))/ H_turn1
        print*, "m_L/H at turn1 =",sqrt(10.0e0_dp**vparams(1,2))/ H_turn1
        print*, "phi at turn1 =", phi(1)
        print*, "----------------"

      end if

      phi_knot = vparams(4,2) - vparams(6,2)/ &
        tan(vparams(5,2))
      if ( phi(1) <= phi_knot .and. &
        phi_test > phi_knot) then

        a_turn2 = a_init*exp(x)
        N_turn2 = x
        H_turn2 = getH(phi,delphi)

        !DEBUG
        print*, "a_* H_* at turn2 =", a_turn2*H_turn2
        print*, "a_* m_H at turn2 =", a_turn2*sqrt(10.0e0_dp**vparams(2,2))
        print*, "k_* at turn2 =", a_turn2*H_turn2/Mpc2Mpl
        print*, "k_m at turn2 =", a_turn2*sqrt(10.0e0_dp**vparams(2,2))/Mpc2Mpl
        print*, "N at turn2 =", N_turn2
        print*, "H at turn2 =", H_turn2
        print*, "m_H/H at turn2 =",sqrt(10.0e0_dp**vparams(2,2))/ H_turn2
        print*, "m_L/H at turn2 =",sqrt(10.0e0_dp**vparams(1,2))/ H_turn2
        print*, "phi at turn2 =", phi(1)
        print*, "----------------"

      end if

      phi_test = phi(1)

    end subroutine find_turning_scales

  END SUBROUTINE odeint_c


  !Uses t as integration variable, copy of odeint_r --- only to be used with
  !background
  SUBROUTINE odeint_with_t(ystart,x1,x2,eps,h1,hmin,derivs,rkqs_r)
    USE ode_path
    USE internals
    USE modpk_observables
    USE modpkparams
    USE potential
    USE modpk_utils, ONLY : reallocate_rv, reallocate_rm, use_t

    IMPLICIT NONE
    real(dp), DIMENSION(:), INTENT(INOUT) :: ystart
    real(dp), INTENT(IN) :: x1,x2,eps,h1,hmin
    !MULTIFIELD
    real(dp), DIMENSION(num_inflaton) :: p, delp
    !END MULTIFIELD

    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         use modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         real(dp), DIMENSION(:), INTENT(IN) :: y
         real(dp), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs

       SUBROUTINE rkqs_r(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
         use modpkparams
         IMPLICIT NONE
         real(dp), DIMENSION(:), INTENT(INOUT) :: y
         real(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
         real(dp), INTENT(INOUT) :: x
         real(dp), INTENT(IN) :: htry,eps
         real(dp), INTENT(OUT) :: hdid,hnext
         INTERFACE
            SUBROUTINE derivs(x,y,dydx)
              use modpkparams
              IMPLICIT NONE
              real(dp), INTENT(IN) :: x
              real(dp), DIMENSION(:), INTENT(IN) :: y
              real(dp), DIMENSION(:), INTENT(OUT) :: dydx
            END SUBROUTINE derivs
         END INTERFACE
       END SUBROUTINE rkqs_r
    END INTERFACE

    real(dp), PARAMETER :: TINY=1.0d-30
    INTEGER, PARAMETER :: MAXSTP=nsteps
    INTEGER*4 :: nstp,i
    real(dp) :: h,hdid,hnext,x,xsav
    real(dp), DIMENSION(SIZE(ystart)) :: dydx, y, yscal
    real(dp) :: z, scalefac

    real(dp) :: infl_efolds, infl_efolds_start
    logical :: infl_checking, stability
    real(dp) :: Nefolds


    !Don't save_steps, bc interfere with N-integrator
    !DEBUG
    save_steps = .false.
    stability = .false.

    !Inits for checking whether in
    !extended inflation period (for IC scan)
    infl_checking = .false.
    infl_efolds_start = 0e0_dp
    infl_efolds = 0e0_dp

    ode_underflow=.FALSE.
    infl_ended=.FALSE.

    x=x1
    h=SIGN(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    y(:)=ystart(:)
    NULLIFY(xp,yp)
    IF (save_steps) THEN
       xsav=x-2.e0_dp*dxsav
       ALLOCATE(xp(256))
       ALLOCATE(yp(SIZE(ystart),SIZE(xp)))
    END IF

    DO nstp=1,MAXSTP

       if (any(isnan(y))) then
         print*, "MODPK in odeint_r"
         print*, "MODPK: y has a NaN value."
         print*, "E-fold",x
         print*, "nstp",nstp
         print*, "y", y
         stop
       end if

       !use_t = .true.

       CALL derivs(x,y,dydx)

       !If get bad deriv, then override this error when IC sampling
       !if (pk_bad==bad_ic) then
       !  return
       !end if

       yscal(:)=ABS(y(:))+ABS(h*dydx(:))+TINY

       IF (save_steps .AND. (ABS(x-xsav) > ABS(dxsav))) &
            CALL save_a_step
       IF ((x+h-x2)*(x+h-x1) > 0.0) h = x2 - x

       CALL rkqs_r(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

       IF (hdid == h) THEN
          nok=nok+1
       ELSE
          nbad=nbad+1
       END IF

       !MULTIFIELD
       p = y(1:num_inflaton)
       delp = y(num_inflaton+1 : 2*num_inflaton)

       Nefolds = y(2*num_inflaton+1)


       IF(getEps_with_t(p,delp) .LT. 1 .AND. .NOT.(slowroll_start)) then
         if (ic_sampling==ic_flags%reg_samp) then
           slowroll_start=.true.
         else
           !If scan ICs, say inflating iff eps<1 for "extended" period,
           !3 efolds - protects against transient inflation epochs, i.e.,
           !at traj turn-around or chance starting with dphi=0

           if (.not. infl_checking) then
             infl_checking = .true.
             infl_efolds_start = Nefolds
           else

             infl_efolds = Nefolds - infl_efolds_start
             if (infl_efolds > 3.0) then
               slowroll_start=.true.
             end if

           end if
         end if
       else if (infl_checking) then
         infl_checking=.false.
       endif
       !END MULTIFIELD

       !Check if H is stable now, so switch to integrate w/N
       call stability_check_on_H(stability,p,delp,using_t=.true.)
       if (stability) then
         !DEBUG
         print*, "Was using t, but H now stable!!"
         stop

         ystart = y
         use_t=.false.
         return
       end if


      if (abs(x2-x)<1e-10) then
        use_t=.false.
        print*, "reached end of t-integration...."

    !DEBUG
    print*, "N", Nefolds
    print*, "eps", getEps_with_t(p, delp)
    print*, "t", x
    print*, "t_end", x2
    print*, "nstp", nstp
    print*, "N last", y(2*num_inflaton+1)
    ystart = y
        return
      end if


       !IF(ode_infl_end) THEN
       !   IF (slowroll_infl_end) THEN
       !      IF(getEps_with_t(p, delp) .GT. 1 .AND. slowroll_start) THEN
       !         infl_ended = .TRUE.
       !         ystart(:) = y(:)
       !         IF (save_steps) CALL save_a_step
       !         RETURN
       !      ENDIF
       !   ELSE
       !      IF(getEps_with_t(p, delp) .GT. 1 .AND. slowroll_start) THEN
       !         PRINT*,'MODPK: You asked for a no-slowroll-breakdown model, but inflation'
       !         PRINT*,'MODPK: already ended via slowroll violation before your phi_end was'
       !         PRINT*,'MODPK: reached. Please take another look at your inputs.'
       !         PRINT*,'MODPK: QUITTING'
       !         PRINT*,'EPSILON =', getEps_with_t(p, delp)
       !         STOP
       !      ENDIF

       !   ENDIF
       !ENDIF

       IF (ode_underflow) RETURN
       IF (ABS(hnext) < hmin) THEN
          write(*,*) 'stepsize smaller than minimum in odeint'
          STOP
       END IF
       h=hnext
    END DO

    PRINT*,'too many steps in odeint_with_t', nstp
    PRINT*,"t=", x, "Step size", h
    print*, "N", Nefolds
    print*, "eps", getEps_with_t(p, delp)
    print*, "t", x
    print*, "t_end", x2
    print*, "nstp", nstp

    ode_underflow=.TRUE.

  CONTAINS

    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
    SUBROUTINE save_a_step
      kount=kount+1
      IF (kount > SIZE(xp)) THEN
         xp=>reallocate_rv(xp,2*SIZE(xp))
         yp=>reallocate_rm(yp,SIZE(yp,1), SIZE(xp))
      END IF
      xp(kount)=x
      yp(:,kount)=y(:)
      xsav=x
    END SUBROUTINE save_a_step



  END SUBROUTINE odeint_with_t

  !When not requiring inflation to end, the stopping requirements will likely
  !vary with the potential that you're using, so here's the function that you
  !will need to edit.
  logical function check_stop_when_not_slowroll_infl_end(phi, dphi) &
    result(stopping)
    use modpkparams, only : potential_choice, phi_init, delsigma
    real(dp), dimension(:), intent(in) :: phi, dphi

    stopping = .false.

    select case(potential_choice)
    case(13)

      ![ LP: ] Complete ad-hoc bullshit
      if (phi(1) < 244.95-0.4) stopping=.true.

    case default
      ! for multifield, determine the total field distance travelled
      if (sqrt(dot_product(phi-phi_init, phi-phi_init)) .gt. (delsigma+0.1)) then
        stopping = .true.
      end if

    end select

  end function check_stop_when_not_slowroll_infl_end

END MODULE modpk_odeint
