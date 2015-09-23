MODULE modpk_utils
  !Module that contains the main subroutines that define the equations we need to
  !solve for the background and the modes, as functions of different variables.
  !Also contains some utility functions for general use and the Runge-Kutta
  !methods.
  use modpkparams
  use modpk_sampling, only : ic_sampling, ic_flags
  use modpk_errorhandling, only : raise, run_outcome
  use modpk_reheat, only : reheater
  use csv_file, only : csv_write
  implicit none

!Define some macros for global use
#include 'modpk_macros.f90'

  interface rkck
     module procedure rkck_r
     module procedure rkck_c
  end interface

  !If true, then switch to using Q variable
  logical, private :: using_q_superh=.false.
  !If true, then switch to using cosmic time; for pre-inflation integration
  logical, private :: using_cosmic_time=.false.

  logical :: use_t

  !Constraints to place on the ODE solutions, for use with the
  !DVODE integrator
  type :: ode_constraints
    integer, dimension(:), allocatable :: vect_indices
    logical, dimension(:), allocatable :: indices_ready
    integer :: num_constraints
    real(dp), dimension(:), allocatable :: lower_bound, upper_bound

    contains

      procedure, public :: init => constraint_initializer
      procedure, public :: set_eps_limits => constraint_set_eps_limits

  end type ode_constraints

  type(ode_constraints) :: dvode_constraints

CONTAINS

  !Wrapper for using the bderivs with the dvode_f90_m integrator
  subroutine bderivs_dvode(neq, t, y, ydot)
    implicit none
    integer, intent (in) :: neq
    real(dp), intent (in) :: t
    real(dp), intent (in) :: y(neq)
    real(dp), intent (out) :: ydot(neq)

    call bderivs(t,y,ydot)

  end subroutine bderivs_dvode

  !Wrapper for using the mode derivs with the dvode_f90_m integrator
  subroutine mode_derivs_dvode(neq, t, y, ydot)
    implicit none
    integer, intent (in) :: neq
    real(dp), intent (in) :: t
    real(dp), intent (in) :: y(neq)
    real(dp), intent (out) :: ydot(neq)

    complex(dp), dimension(neq/2) :: y_comp, ydot_comp


    !Take the real values and pack into complex

    y_comp = cmplx(y(1:neq/2), y(neq/2+1:neq),kind=dp)

    call derivs(t,y_comp,ydot_comp)

    !Take the complex values and unpack into reals
    ydot(1:neq/2) = real(ydot_comp,kind=dp)
    ydot(neq/2+1:neq) = aimag(ydot_comp)

  end subroutine mode_derivs_dvode

  !Wrapper for using the qderivs with the dvode_f90_m integrator
  subroutine qderivs_dvode(neq, t, y, ydot)
    implicit none
    integer, intent (in) :: neq
    real(dp), intent (in) :: t
    real(dp), intent (in) :: y(neq)
    real(dp), intent (out) :: ydot(neq)

    complex(dp), dimension(neq/2) :: y_comp, ydot_comp


    !Take the real values and pack into complex

    y_comp = cmplx(y(1:neq/2), y(neq/2+1:neq),kind=dp)

    call qderivs(t,y_comp,ydot_comp)

    !Take the complex values and unpack into reals
    ydot(1:neq/2) = real(ydot_comp,kind=dp)
    ydot(neq/2+1:neq) = aimag(ydot_comp)

  end subroutine qderivs_dvode


  !Background derivatives y'=f(y)
  SUBROUTINE bderivs(x,y,yprime)
    USE modpkparams
    USE potential, ONLY: pot,dVdphi,d2Vdphi2,getH,getdHdalpha,getEps, &
      getEps_with_t, getH_with_t
    USE camb_interface, ONLY : pk_bad
    use modpk_deltaN, only : V_i_sum_sep
    real(dp), INTENT(IN) :: x
    real(dp), DIMENSION(:), INTENT(IN) :: y
    real(dp), DIMENSION(:), INTENT(OUT) :: yprime

    !MULTIFIELD
    !real(dp), DIMENSION(size(y)/2) :: phi, delphi
    real(dp), DIMENSION(num_inflaton) :: phi, delphi
    real(dp) :: hubble,dhubble, eps
    real(dp), dimension(num_inflaton) :: rho_radn, rho_fields
    !END MULTIFIEND

    integer :: i
    !
    !     x=alpha
    !     y(1:num_inflaton)=phi
    !     dydx(1:num_inflaton)=dphi/dalpha
    !     y(num_inflaton+1:2*num_inflaton)=dphi/dalpha
    !     dydx(num_inflaton+1:2*num_inflaton)=d^2phi/dalpha^2
    !
    !     Can also evolve a radiation fluid
    !     If using_t, then final portion of y is e-folds

    !MULTIFIELD
    phi = y(IND_FIELDS)
    delphi = y(IND_VEL)
    !END MULTIFIELD

    if (use_t) then
      eps = getEps_with_t(phi,delphi)
    else
      eps = getEps(phi,delphi)
    end if

    !Instability check since H^2=V/(3-eps) is not
    !numerically stable as V~0 for
    !H>0, since requires eps->3
    IF(eps .ge. 3.0e0_dp) THEN
       write(*,*) 'MODECODE: Pot=', pot(phi)
       write(*,*) 'MODECODE: Eps=',eps
       write(*,*) 'MODECODE: Using t?', use_t
       if (use_t) then
         write(*,*) 'MODECODE: t=',x
       else
         write(*,*) 'MODECODE: E-fold=',x
       end if
       write(*,*) "MODECODE: Phi=",phi
       write(*,*) "MODECODE: Dphi=",delphi
       write(*,*) 'MODECODE: vparams= ', (vparams(i,:),i=1,size(vparams,1))
       if (.not.instreheat) write(*,*) 'MODECODE: N_pivot: ', N_pivot

       !Can sometimes get here when IC sampling
       !close to the point where V=0, since H^2=V/(3-eps)
       !Might override this error
       call raise%fatal_cosmo(&
         'H is imaginary in bderivs.  &
         Check if V~=0, since makes H unstable.  &
         Or check if start step-size too large.  &
         You might be able to override this error if &
         you know how to auto-correct it.',&
         __FILE__, __LINE__)

       !In the case of the hilltop potential, the integrator
       !in a trial step can go here very occasionally because
       !the trial step is too large and it has come too close to V=0.
       !We will stop it going this way, and the code will find the
       !correct epsilon=1 point which has, by definition, to be
       !before this problematic region is reached.
       IF(potential_choice.eq.6) THEN
          yprime(1)=0.0e0_dp
          yprime(2)=0.0e0_dp
          RETURN
       ENDIF
    END IF

    !MULTIFIELD
    if (use_t) then
      !Derivs in cosmic time
      hubble = getH_with_t(phi,delphi)

      yprime(IND_FIELDS) = delphi
      yprime(IND_VEL) = -3.0e0_dp*hubble*delphi - dVdphi(phi)

      !E-folds
      yprime(IND_EFOLDS) = hubble

      !Auxiliary constraints
      if (tech_opt%use_dvode_integrator .and. &
        tech_opt%use_ode_constraints) then
       call raise%fatal_cosmo(&
         'Auxiliary constraints not yet implemented with t-derivs.', &
         __FILE__, __LINE__)
      end if

    else if (reheater%evolving_gamma) then
      !Include a post-inflationary radiation field
      !Derivs in cosmic time

      rho_radn = y(IND_RADN)
      hubble = reheater%getH_with_radn(phi, delphi, sum(rho_radn))
      dhubble = reheater%getdH_with_radn(phi, delphi, sum(rho_radn))

      rho_fields = 0.5e0_dp*hubble**2*delphi**2 + V_i_sum_sep(phi)

      if (any(rho_radn<0) .or. any (rho_fields<0)) then
        print*, "MODECODE: rho_radn:", rho_radn
        print*, "MODECODE: rho_fields:", rho_fields
        call raise%fatal_cosmo(&
          'The individual energy densities are negative.',&
          __FILE__, __LINE__)
      end if

      !Fields
      yprime(IND_FIELDS) = delphi
      yprime(IND_VEL) = &
        -(3.0e0_dp+ reheater%Gamma_i/hubble + dhubble/hubble)*delphi &
        - dVdphi(phi)/hubble**2

      !E-folds
      yprime(IND_EFOLDS) = hubble

      !Radiation
      yprime(IND_RADN) = &
        -4.0e0_dp*rho_radn &
        + reheater%Gamma_i*rho_fields/hubble

      !Auxiliary constraints
      if (tech_opt%use_dvode_integrator .and. &
        tech_opt%use_ode_constraints) then
        yprime(IND_CONST_EPS_RADN) = &
          sum(yprime(IND_FIELDS)*yprime(IND_VEL))
      end if

    else

      !Normal

      !Derivs wrt e-folds
      hubble=getH(phi,delphi)
      dhubble=getdHdalpha(phi,delphi)

      yprime(IND_FIELDS) = delphi
      yprime(IND_VEL) = &
        -((3.0e0_dp+dhubble/hubble)*delphi+&
        dVdphi(phi)/hubble/hubble)

      !Auxiliary constraints
      if (tech_opt%use_dvode_integrator .and. &
        tech_opt%use_ode_constraints) then
        yprime(IND_CONST_EPS_BACK) = &
          sum(yprime(IND_FIELDS)*yprime(IND_VEL))
      end if

    end if


    !END MULTIFIELD

  END SUBROUTINE bderivs


  ! Full y (back+mode matrix+tensor) derivatives y'=f(y)
  SUBROUTINE derivs(x,y,yprime)
    USE modpkparams
    USE internals
    USE potential, ONLY: pot, dVdphi, d2Vdphi2, getH, getdHdalpha, getEps, getEta
    USE camb_interface, ONLY : pk_bad
    IMPLICIT NONE
    real(dp), INTENT(IN) :: x
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: yprime

    ! background quantity
    real(dp) :: hubble, dhubble, scale_factor, epsilon_, dotphi, eta
    real(dp) :: thetaN2, Vzz, Vz, grad_V  !! thetaN2 = (d\theta/dNe)^2
    real(dp), dimension(num_inflaton) :: phi, delphi, Vp
    real(dp), dimension(num_inflaton, num_inflaton) :: Cab, d2V

    complex(dp), dimension(num_inflaton*num_inflaton) :: psi, dpsi ! scalar ptb mode matrix
    complex(dp) :: v_tensor, dv_tensor                             ! tensor perturbations
    complex(dp) :: u_zeta, du_zeta

    integer :: i, j

    !     x = alpha
    !     y(1:n) = phi                 dydx(1:n)=dphi/dalpha
    !     y(n+1:2n) = dphi/dalpha      dydx(n+1:2n)=d^2phi/dalpha^2
    !     y(2n+1:2n+n**2) = psi               dydx(2n+1:3n)=dpsi/dalpha
    !     y(2n+n**2+1:2n+2n**2) = dpsi/dalpha       dydx(3n+1:4n)=d^2psi/dalpha^2
    !     y(2n+2n**2+1) = v                  dydx(4n+1)=dv/dalpha
    !     y(2n+2n**2+2) = dv/dalpha          dydx(4n+2)=d^2v/dalpha^2
    !     y(2n+2n**2+3) = u_zeta     dydx(4n+4)=d^2u_zeta/dalpha^2
    !     y(2n+2n**2+4) = du_zeta/dalpha     dydx(4n+4)=d^2u_zeta/dalpha^2

    phi = real(y(IND_FIELDS),kind=dp)
    delphi = real(y(IND_VEL),kind=dp)
    dotphi = sqrt(dot_product(delphi, delphi))

    !Aliases to potential derivatives
    epsilon_ = getEps(phi, delphi)
    eta = geteta(phi, delphi)
    d2V = d2Vdphi2(phi)
    Vp = dVdphi(phi)
    Vzz = dot_product(delphi, matmul(d2V, delphi))/dotphi**2
    Vz = dot_product(Vp, delphi)/dotphi
    grad_V = sqrt(dot_product(Vp, Vp))

    IF(dot_product(delphi, delphi) .GT. 6.e0_dp) THEN
       write(*,*) 'MODECODE: Using t?', use_t
       if (use_t) then
         write(*,*) 'MODECODE: t=',x
       else
         write(*,*) 'MODECODE: E-fold=',x
       end if
       write(*,*) 'MODECODE: vparams= ', (vparams(i,:),i=1,size(vparams,1))
       if (.not.instreheat) write(*,*) 'MODECODE: N_pivot: ', N_pivot

       !Can sometimes get here when IC sampling
       !close to the point where V=0, since H^2=V/(3-eps)
       !Might override this error
       call raise%fatal_cosmo(&
         'H is imaginary in derivs.  &
         Check if V~=0, since makes H unstable.  &
         Or check if start step-size too large.  &
         You might be able to override this error if &
         you know how to auto-correct it.',&
         __FILE__, __LINE__)


    END IF

    hubble=getH(phi,delphi)
    dhubble=getdHdalpha(phi,delphi)
    scale_factor=a_init*exp(x)

    ! Alias y's into real variable names
    ! NB: if using_q_superh, psi(i,j)-->q_ptb(i,j)
    psi = y(IND_MODES)
    dpsi = y(IND_MODES_VEL)

    v_tensor  = y(IND_TENSOR)
    dv_tensor  = y(IND_TENSOR_VEL)
    u_zeta = y(IND_UZETA)
    du_zeta = y(IND_UZETA_VEL)

    ! Build the mass matrix, Cab
    call build_mass_matrix(Cab)

    ! -----------------------------
    ! Set the RHS of y'(x) = f(y(x))
    ! -----------------------------

    ! background
    yprime(IND_FIELDS) = cmplx(delphi,kind=dp)
    yprime(IND_VEL) =&
      cmplx(-((3.0e0_dp+dhubble/hubble)*delphi+dVdphi(phi)/hubble/hubble),kind=dp)

    ! ptb matrix
    yprime(IND_MODES) = dpsi

    if (using_q_superh) then
      yprime(IND_MODES_VEL) = -(3.0e0_dp - epsilon_)*dpsi &
        - (k/scale_factor/hubble)**2*psi &
        - dot(Cab, psi)/hubble**2
    else
      yprime(IND_MODES_VEL) = -(1.0e0_dp - epsilon_)*dpsi &
        - (k/scale_factor/hubble)**2*psi &
        + (2.0e0_dp - epsilon_)*psi - dot(Cab, psi)/hubble**2
    end if

    ! tensors
    yprime(IND_TENSOR) = dv_tensor
    if (using_q_superh) then
      yprime(IND_TENSOR_VEL) = -(3.0e0_dp - epsilon_)*dv_tensor - &
        (k/scale_factor/hubble)**2*v_tensor
    else
      yprime(IND_TENSOR_VEL) = -(1.0e0_dp - epsilon_)*dv_tensor - &
        (k/scale_factor/hubble)**2*v_tensor + (2.0e0_dp - epsilon_)*v_tensor
    end if

    ! adiabatic ptb
    yprime(IND_UZETA) = du_zeta
    thetaN2 = (grad_V + Vz)*(grad_V - Vz)/(dotphi*hubble**2)**2
    yprime(IND_UZETA_VEL) = -(1.0e0_dp - epsilon_)*du_zeta -&
      (k/scale_factor/hubble)**2*u_zeta &
      + (2.0e0_dp + 5.0e0_dp*epsilon_ - 2.0e0_dp*epsilon_**2 + &
      2.0e0_dp*epsilon_*eta + thetaN2 - Vzz/hubble**2)*u_zeta

    if (tech_opt%use_dvode_integrator .and. &
      tech_opt%use_ode_constraints) then
      yprime(IND_CONST_EPS_MODES) = &
        cmplx(sum(yprime(IND_FIELDS)*yprime(IND_VEL)),kind=dp)
    end if

    contains

      subroutine build_mass_matrix(mass_matrix)

        real(dp), dimension(num_inflaton, num_inflaton), intent(out) :: mass_matrix

        if (potential_choice .eq. 7) then
          ! for exponential potential 7, Cab is exactly zero
          ! set this in order to prevent numerical error
          mass_matrix = 0e0_dp
        else
           forall (i=IND_FIELDS, j=IND_FIELDS) &
                mass_matrix(i,j) = d2V(i,j) +  &
                (delphi(i)*Vp(j) + delphi(j)*Vp(i)) &
                + (3e0_dp-epsilon_)*hubble**2 * delphi(i)*delphi(j)
        end if

      end subroutine build_mass_matrix

      ! Only works for square matrices
      pure function dot(matrixA,hacked_vector) result(outvect)

        real(dp), dimension(:,:), intent(in) :: matrixA
        complex(dp), dimension(size(matrixA,1),size(matrixA,2)) :: matrixB, &
          matrixC
        complex(dp), dimension(:), intent(in) :: hacked_vector
        integer :: i, j, k, n
        complex(dp), dimension(size(hacked_vector)) :: outvect

        n=size(matrixA,1)

        outvect=0e0_dp
        matrixB=0e0_dp
        matrixC=0e0_dp

        matrixB=convert_hacked_vector_to_matrix(hacked_vector)

        matrixC = matmul(matrixA, matrixB)

        outvect = convert_matrix_to_hacked_vector(matrixC)

      end function dot

  END SUBROUTINE derivs

  ! Full y (back+mode matrix(in Q at horizon cross)+tensor) derivatives y'=f(y)
  subroutine qderivs(x,y,yprime)
    use modpkparams
    implicit none
    real(dp), intent(in) :: x
    complex(kind=dp), dimension(:), intent(in) :: y
    complex(kind=dp), dimension(:), intent(out) :: yprime

    using_q_superh=.true.

    call derivs(x,y,yprime)

    using_q_superh=.false.

  end subroutine qderivs


  SUBROUTINE rkqs_r(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    USE ode_path
    use camb_interface, only : pk_bad
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
    INTEGER*4 :: ndum
    real(dp) :: errmax,h,htemp,xnew
    real(dp), DIMENSION(size(y)) :: yerr,ytemp
    real(dp), PARAMETER :: SAFETY=0.9e0_dp,PGROW=-0.2e0_dp,PSHRNK=-0.25e0_dp,&
         ERRCON=1.89e-4_dp

    if (size(y)==size(dydx) .and. size(dydx)==size(yscal)) then
       ndum = size(y)
    else
      call raise%fatal_code(&
       'Wrong array sizes in rkqs', __FILE__, __LINE__)
    end if
    h=htry
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)

       !For possible error overrides in derivs
       if (pk_bad /= run_outcome%success) return

       errmax=maxval(abs(yerr(:)/yscal(:)))/eps
       if (errmax <= 1.0e0_dp) exit
       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1e0_dp*abs(h)),h)
       xnew=x+h
       !Errors if use xnew==x, as one might expect
       !if (xnew == x) then
       if (abs(h)<1e-30_dp) then
         call raise%warning('Stepsize underflow in rkqs_r',&
           __FILE__, __LINE__)
          ode_underflow = .true.
          return
       endif
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0e0_dp*h
    end if
    hdid=h
    x=x+h
    y(:)=ytemp(:)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END SUBROUTINE rkqs_r

  SUBROUTINE rkqs_c(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    USE ode_path
    USE modpkparams
    IMPLICIT NONE
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(INOUT) :: y
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: yscal, dydx
    real(dp), INTENT(INOUT) :: x
    real(dp), INTENT(IN) :: htry,eps
    real(dp), INTENT(OUT) :: hdid,hnext
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         USE modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE

    INTEGER*4 :: ndum
    real(dp) :: errmax,h,htemp,xnew
    COMPLEX(KIND=DP), DIMENSION(size(y)) :: yerr,ytemp
    real(dp), PARAMETER :: SAFETY=0.9e0_dp,PGROW=-0.2e0_dp,PSHRNK=-0.25e0_dp,&
         ERRCON=1.89e-4_dp

    real(dp) :: yerr_r(2*size(y)), yscal_r(2*size(yscal))

    if (size(y)==size(dydx) .and. size(dydx)==size(yscal)) then
       ndum = size(y)
    else
      call raise%fatal_code(&
       'Wrong array sizes in rkqs', __FILE__, __LINE__)
    end if
    h=htry
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)

       ! in doing error estimation and step size rescaling, we switch to real components
       yerr_r(1 : size(yerr)) = real(yerr,kind=dp)
       yerr_r(size(yerr)+1 : 2*size(yerr)) = real(yerr*(0,-1),kind=dp)

       yscal_r(1 : size(yscal)) = real(yscal,kind=dp)
       yscal_r(size(yscal)+1 : 2*size(yscal)) = real(yscal*(0,-1),kind=dp)

       errmax=maxval(abs(yerr_r(:)/yscal_r(:)))/eps

       if (errmax <= 1.0e0_dp) exit

       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1e0_dp*abs(h)),h)
       xnew=x+h

       !Errors if use xnew==x, as one might suspect
       !if (xnew == x) then
       if (abs(h)<1e-30_dp) then
         call raise%warning('Stepsize underflow in rkqs_c',&
           __FILE__, __LINE__)
          ode_underflow = .true.
          return
       endif
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0e0_dp*h
    end if
    hdid=h
    x=x+h
    y(:)=ytemp(:)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END SUBROUTINE rkqs_c

  SUBROUTINE rkck_r(y,dydx,x,h,yout,yerr,derivs)
    IMPLICIT NONE
    real(dp), DIMENSION(:), INTENT(IN) :: y,dydx
    real(dp), INTENT(IN) :: x,h
    real(dp), DIMENSION(:), INTENT(OUT) :: yout,yerr
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         use modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         real(dp), DIMENSION(:), INTENT(IN) :: y
         real(dp), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    INTEGER*4 :: ndum
    real(dp), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    real(dp), PARAMETER :: A2=0.2e0_dp,A3=0.3e0_dp,A4=0.6e0_dp,A5=1.0e0_dp,&
         A6=0.875e0_dp,B21=0.2e0_dp,B31=3.0e0_dp/40.0e0_dp,B32=9.0e0_dp/40.0e0_dp,&
         B41=0.3e0_dp,B42=-0.9e0_dp,B43=1.2e0_dp,B51=-11.0e0_dp/54.0e0_dp,&
         B52=2.5e0_dp,B53=-70.0e0_dp/27.0e0_dp,B54=35.0e0_dp/27.0e0_dp,&
         B61=1631.0e0_dp/55296.0e0_dp,B62=175.0e0_dp/512.0e0_dp,&
         B63=575.0e0_dp/13824.0e0_dp,B64=44275.0e0_dp/110592.0e0_dp,&
         B65=253.0e0_dp/4096.0e0_dp,C1=37.0e0_dp/378.0e0_dp,&
         C3=250.0e0_dp/621.0e0_dp,C4=125.0e0_dp/594.0e0_dp,&
         C6=512.0e0_dp/1771.0e0_dp,DC1=C1-2825.0e0_dp/27648.0e0_dp,&
         DC3=C3-18575.0e0_dp/48384.0e0_dp,DC4=C4-13525.0e0_dp/55296.0e0_dp,&
         DC5=-277.0e0_dp/14336.0e0_dp,DC6=C6-0.25e0_dp

    if (size(y)==size(dydx) .and. size(dydx)==size(yout) .and. size(yout)==size(yerr)) then
       ndum = size(y)
    else
      call raise%fatal_code(&
       'Wrong array sizes in rkck', __FILE__, __LINE__)
    end if

    ytemp=y+B21*h*dydx
    call derivs(x+A2*h,ytemp,ak2)
    ytemp=y+h*(B31*dydx+B32*ak2)
    call derivs(x+A3*h,ytemp,ak3)
    ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
    call derivs(x+A4*h,ytemp,ak4)
    ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
    call derivs(x+A5*h,ytemp,ak5)
    ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call derivs(x+A6*h,ytemp,ak6)
    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END SUBROUTINE rkck_r

  SUBROUTINE rkck_c(y,dydx,x,h,yout,yerr,derivs)
    USE modpkparams
    IMPLICIT NONE
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y, dydx
    real(dp), INTENT(IN) :: x,h
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: yout,yerr
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         USE modpkparams
         IMPLICIT NONE
         real(dp), INTENT(IN) :: x
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
         COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    INTEGER*4 :: ndum
    COMPLEX(KIND=DP), DIMENSION(size(y)) :: ytemp, ak2,ak3,ak4,ak5,ak6
    real(dp), PARAMETER :: A2=0.2e0_dp,A3=0.3e0_dp,A4=0.6e0_dp,A5=1.0e0_dp,&
         A6=0.875e0_dp,B21=0.2e0_dp,B31=3.0e0_dp/40.0e0_dp,B32=9.0e0_dp/40.0e0_dp,&
         B41=0.3e0_dp,B42=-0.9e0_dp,B43=1.2e0_dp,B51=-11.0e0_dp/54.0e0_dp,&
         B52=2.5e0_dp,B53=-70.0e0_dp/27.0e0_dp,B54=35.0e0_dp/27.0e0_dp,&
         B61=1631.0e0_dp/55296.0e0_dp,B62=175.0e0_dp/512.0e0_dp,&
         B63=575.0e0_dp/13824.0e0_dp,B64=44275.0e0_dp/110592.0e0_dp,&
         B65=253.0e0_dp/4096.0e0_dp,C1=37.0e0_dp/378.0e0_dp,&
         C3=250.0e0_dp/621.0e0_dp,C4=125.0e0_dp/594.0e0_dp,&
         C6=512.0e0_dp/1771.0e0_dp,DC1=C1-2825.0e0_dp/27648.0e0_dp,&
         DC3=C3-18575.0e0_dp/48384.0e0_dp,DC4=C4-13525.0e0_dp/55296.0e0_dp,&
         DC5=-277.0e0_dp/14336.0e0_dp,DC6=C6-0.25e0_dp

    if (size(y)==size(dydx) .and. size(dydx)==size(yout) .and. size(yout)==size(yerr)) then
       ndum = size(y)
    else
      call raise%fatal_code(&
       'Wrong array sizes in rkqs', __FILE__, __LINE__)
    end if

    ytemp=y+B21*h*dydx
    call derivs(x+A2*h,ytemp,ak2)
    ytemp=y+h*(B31*dydx+B32*ak2)
    call derivs(x+A3*h,ytemp,ak3)
    ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
    call derivs(x+A4*h,ytemp,ak4)
    ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
    call derivs(x+A5*h,ytemp,ak5)
    ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call derivs(x+A6*h,ytemp,ak6)
    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END SUBROUTINE rkck_c



  FUNCTION reallocate_rv(p,n)
    real(dp), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER*4, INTENT(IN) :: n
    INTEGER*4 :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)  !! allocate memeory of size n at new address to be returned
    if (ierr /= 0) then
      call raise%fatal_code(&
       'reallocate_rv: problem in attempt to allocate memory', __FILE__, __LINE__)
    end if
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))  !! copy old array contents to the new memory block
    deallocate(p)  !! free the old meomory block
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END FUNCTION reallocate_rv

  FUNCTION reallocate_rm(p,n,m)
    real(dp), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER*4, INTENT(IN) :: n,m
    INTEGER*4 :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) then
      call raise%fatal_code(&
       'reallocate_rm: problem in attempt to allocate memory', __FILE__, __LINE__)
    end if
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))= p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END FUNCTION reallocate_rm


  ! Take a "hacked" vector, i.e., a matrix B that is in the
  !form of a vector where the rows of B(i,1:N) --> first N terms in B(1:N), etc
  !Returns the matrix form. NB: the vector must be a representation of a square
  !matrix.

  pure function convert_hacked_vector_to_matrix(matrix_as_vector) &
    result(real_matrix)

    complex(dp), dimension(:), intent(in) :: matrix_as_vector
    complex(dp), dimension(:,:), allocatable :: real_matrix
    integer :: n, i, j

    n=int(sqrt(real(size(matrix_as_vector))))
    allocate(real_matrix(n,n))

    forall (i=1:n, j=1:n) &
      real_matrix(i,j) = matrix_as_vector((i-1)*n+j)

  end function convert_hacked_vector_to_matrix

  ! Take a matrix and "hack" it into a vector, i.e.,
  !a vector where the rows of B(i,1:N) --> first N terms in B(1:N), etc
  pure function convert_matrix_to_hacked_vector(matrix) &
    result(hacked_vector)

    complex(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), allocatable :: hacked_vector
    integer :: n, i, j

    n=size(matrix,1)*size(matrix,2)
    allocate(hacked_vector(n))

    forall (i=1:size(matrix,1), j=1:size(matrix,2)) &
      hacked_vector((i-1)*size(matrix,2)+j) = matrix(i,j)

  end function convert_matrix_to_hacked_vector

  !For ODE system y_i'(t) = f_i[y_j(t)], this returns the Jacobian df_i/dy_j,
  !with df(i)/dy(j) loaded into PD(i,j).
  !Valid for background evolution.  Used with the DVODE integrator when invoking
  !stiff solving methods.
  subroutine jacobian_background_DVODE(neq, t, y, ml, mu, pd, nrowpd)
    use modpkparams
    use potential, only : getEps, pot, dVdphi, d2Vdphi2
    implicit none

    integer, intent(in) :: neq, nrowpd, ml, mu
    real(dp), intent(in) :: t, y(neq)
    real(dp), intent(out) :: pd(nrowpd,neq)

    real(dp), dimension(num_inflaton) :: phi, dphi, dV
    real(dp) :: eps, V_pot
    real(dp), dimension(num_inflaton, num_inflaton) :: d2V
    real(dp), dimension(2*num_inflaton, 2*num_inflaton) :: delta

    integer :: ii, jj

    !f(i) = y(i+num_inflaton)  for   1 <= i <= num_inflaton
    !f(i) = -(3-eps)*y(i-num_inflaton) - H**-2 * dVdphi(i-num_inflaton)
    !                 for   num_inflaton + 1 <= i <= 2*num_inflaton

    phi = y(IND_FIELDS)
    dphi = y(IND_VEL)

    eps = getEps(phi,dphi)
    V_pot = pot(phi)
    dV = dVdphi(phi)
    d2V = d2Vdphi2(phi)

    !Bc 2nd order ODE
    pd = 0e0_dp
    do ii=1,num_inflaton
      pd(ii,ii+num_inflaton) = 1.0e0_dp
    end do

    !Derivs wrt phi
    forall ( ii=IND_FIELDS, jj=IND_FIELDS)&
      pd(ii+num_inflaton,jj) = -(3.0e0_dp - eps)*&
        (d2V(ii,jj)/V_pot - (dV(ii)*dV(jj)/V_pot**2))

    !Derivs wrt dphi
    delta = 0e0_dp
    do ii=1,2*num_inflaton
      delta(ii,ii) = 1e0_dp
    end do

    forall ( ii=IND_FIELDS, jj=IND_FIELDS)&
      pd(ii+num_inflaton,jj+num_inflaton) = &
        -(3.0e0_dp-eps)*delta(ii,jj) &
        + dphi(ii)*dphi(jj) &
        + dV(ii)*dphi(jj)/V_pot

  end subroutine jacobian_background_DVODE

  !For ODE system y_i'(t) = f_i[y_j(t)], this returns the Jacobian df_i/dy_j,
  !with df(i)/dy(j) loaded into PD(i,j).
  !Valid for mode evolution with psi variable, ie, not when use_q.
  !Used with the DVODE integrator when invoking stiff solving methods.
  subroutine jacobian_psi_modes_DVODE(neq, t, y, ml, mu, pd, nrowpd)
    use modpkparams
    use internals
    USE potential, ONLY: pot, dVdphi, d2Vdphi2, d3Vdphi3, &
      getH, getdHdalpha, getEps, getEta
    implicit none

    integer, intent(in) :: neq, nrowpd, ml, mu
    real(dp), intent(in) :: t, y(neq)
    real(dp), intent(out) :: pd(nrowpd,neq)

    complex(dp), dimension(neq/2) :: y_comp, ydot_comp

    real(dp), dimension(2*num_inflaton,2*num_inflaton) :: pd_back
    real(dp), dimension(2*num_inflaton) :: y_back

    real(dp), dimension(num_inflaton) :: phi, dphi
    complex(dp), dimension(num_inflaton*num_inflaton) :: psi, dpsi
    real(dp), dimension(nrowpd,neq) :: delta
    integer :: ii, jj, kk, ll

    real(dp) :: hubble, dhubble, scale_factor, epsilon_, dotphi, eta
    real(dp) :: Vzz, Vz, grad_V, V_pot
    real(dp), dimension(num_inflaton) :: dV
    real(dp), dimension(num_inflaton, num_inflaton) :: Cab, d2V
    real(dp), dimension(num_inflaton, num_inflaton, num_inflaton) :: dCdphi, &
      dCddelphi, d3V


    y_back = y(1:2*num_inflaton)
    call jacobian_background_DVODE(NEQ=size(pd_back,2), t=t, y=y_back, &
      ml=ml, mu=mu, pd=pd_back, nrowpd=size(pd_back,1))
    pd(1:size(pd_back,1),1:size(pd_back,2)) = pd_back

    phi = real(y(IND_FIELDS),kind=dp)
    dphi = real(y(IND_VEL),kind=dp)
    psi = y(IND_MODES)
    dpsi = y(IND_MODES_VEL)

    delta=0e0_dp
    do ii=1,nrowpd
      delta(ii,ii) = 1e0_dp
    end do

    !Aliases to potential derivatives
    epsilon_ = getEps(phi, dphi)
    V_pot = pot(phi)
    d3V = d3Vdphi3(phi)
    d2V = d2Vdphi2(phi)
    dV = dVdphi(phi)
    Vzz = dot_product(dphi, matmul(d2V, dphi))/dotphi**2
    Vz = dot_product(dV, dphi)/dotphi
    grad_V = sqrt(dot_product(dV, dV))
    hubble=getH(phi,dphi)
    dhubble=getdHdalpha(phi,dphi)
    scale_factor=a_init*EXP(t)

    ! Build the mass matrix, Cab and its derivatives
    forall (ii=IND_FIELDS, jj=IND_FIELDS) &
         Cab(ii,jj) = d2V(ii,jj) +  &
         (dphi(ii)*dV(jj) + dphi(jj)*dV(ii)) &
         + (3e0_dp-epsilon_)*hubble**2 * dphi(ii)*dphi(jj)

    forall (ii=IND_FIELDS, ll=IND_FIELDS, kk=IND_FIELDS)
        dCdphi(ii,ll,kk) = (1.0e0_dp/hubble**2)*&
          (d3V(ii,ll,kk) - (dV(kk)/V_pot)*d2V(ii,ll) &
          - (dV(kk)/V_pot)*(dphi(ii)*dV(ll)+dphi(ll)*dV(ii)) &
          + dphi(ii)*d2V(ll,kk) + dphi(ll)*d2V(ii,kk))

        dCddelphi(ii,ll,kk) = (1.0e0_dp/V_pot)*&
          (-dphi(kk)*d2V(ii,ll) &
          -dphi(kk)*(dphi(ii)*dV(ll) + dphi(ll)*dV(ii)) &
          +(3.0e0_dp - epsilon_)*(delta(ii,kk)*dV(ll) + delta(ll,kk)*dV(ii))&
          -dphi(kk)*dphi(ii)*dphi(ll) &
          +(3.0e0_dp - epsilon_)*delta(ii,kk)*dphi(ll)&
          +(3.0e0_dp - epsilon_)*dphi(ii)*delta(ll,kk))

    end forall


    if (using_q_superh) then
      !DEBUG
      print*, "using_q_superh not yet implemented in jacobian_psi_modes_DVODE"
      stop
    else
      !DEBUG
      print*, "Cab", Cab
      print*, "dCdphi", dCdphi
      print*, "dCddelphi", dCddelphi
      print*, "back jacobian:"
      print*, pd_back
      stop
    end if


  end subroutine jacobian_psi_modes_DVODE

  subroutine constraint_initializer(self, num_constraints)
    class(ode_constraints) :: self
    integer, intent(in) :: num_constraints

    !Only for use with DVODE integrator
    if (.not. tech_opt%use_dvode_integrator) then
      call raise%fatal_code(&
        "Need to use DVODE integrator if &
        you require auxiliary constraints.",&
        __FILE__,__LINE__)
    end if

    !Initialize
    self%num_constraints = num_constraints
    if (allocated(self%vect_indices)) deallocate(self%vect_indices)
    if (allocated(self%indices_ready)) deallocate(self%indices_ready)
    if (allocated(self%lower_bound)) deallocate(self%lower_bound)
    if (allocated(self%upper_bound)) deallocate(self%upper_bound)
    allocate(self%vect_indices(self%num_constraints))
    allocate(self%lower_bound(self%num_constraints))
    allocate(self%upper_bound(self%num_constraints))
    allocate(self%indices_ready(self%num_constraints))

    self%indices_ready = .false.
    self%vect_indices = -1

  end subroutine constraint_initializer

  subroutine constraint_set_eps_limits(self, &
      evolve_modes, evolve_radn_back)
    class(ode_constraints) :: self
    logical, intent(in) :: evolve_modes, evolve_radn_back

    integer :: ii

    !Consistency
    if (.not. allocated(self%vect_indices) .or. &
        .not. allocated(self%lower_bound) .or. &
        .not. allocated(self%upper_bound)) then
      call raise%fatal_code(&
        "Initialize DVODE constraints &
        prior to using this subroutine.",&
        __FILE__,__LINE__)
    end if

    !Find next available position in vect_indices
    do ii=1,size(self%vect_indices)
      if (.not. self%indices_ready(ii)) then
        !Set the constraint on this index
        self%indices_ready(ii) = .true.

        if (evolve_modes) then
          self%vect_indices(ii) = IND_CONST_EPS_MODES
        else if (evolve_radn_back) then
          self%vect_indices(ii) = IND_CONST_EPS_RADN
        else
          self%vect_indices(ii) = IND_CONST_EPS_BACK
        end if

        !Set upper and lower limits for epsilon
        self%lower_bound(ii) = 0.0e0_dp
        self%upper_bound(ii) = 3.0e0_dp

        exit
      end if

      !If get to end of vect_indices, then they're already all set
      if (ii==size(self%vect_indices)) then
        call raise%fatal_code(&
          "All constraint indices used prior to setting &
          epsilon<3 constraint.",&
          __FILE__,__LINE__)
      end if

    end do


  end subroutine constraint_set_eps_limits


end module modpk_utils
