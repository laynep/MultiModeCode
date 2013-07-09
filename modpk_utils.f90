MODULE modpk_utils
  use modpkparams
  IMPLICIT NONE

  INTERFACE rkck
     module procedure rkck_r
     module procedure rkck_c
  END INTERFACE

  ![ LP: ] Local variable; if true, then switch to using Q variable
  logical, private :: using_q_superh=.false.

CONTAINS

  ![ LP: ] Background derivatives y'=f(y)
  SUBROUTINE bderivs(x,y,yprime)
    USE modpkparams
    USE potential, ONLY: pot,dVdphi,d2Vdphi2,getH,getHdot
    USE camb_interface, ONLY : pk_bad
    real(dp), INTENT(IN) :: x
    real(dp), DIMENSION(:), INTENT(IN) :: y
    real(dp), DIMENSION(:), INTENT(OUT) :: yprime

    !MULTIFIELD
    real(dp), DIMENSION(size(y)/2) :: p, delp
    real(dp) :: hubble,dhubble
    !END MULTIFIEND

    integer :: i
    !
    !     x=alpha
    !     y(FIRST HALF)=phi              dydx(FIRST HALF)=dphi/dalpha
    !     y(SECOND HALF)=dphi/dalpha      dydx(SECOND HALF)=d^2phi/dalpha^2
    !

    !MULTIFIELD
    p = y(1 : size(y)/2)
    delp = y(size(y)/2+1 : size(y))
    !END MULTIFIELD

    IF(dot_product(delp,delp) .gt. 6.0d0) THEN
       WRITE(*,*) 'MODPK: H is imaginary:',x,p,delp
       !In the case of the hilltop potential, the integrator
       !in a trial step can go here very occasionally because
       !the trial step is too large and it has come too close to V=0.
       !We will stop it going this way, and the code will find the
       !correct epsilon=1 point which has, by definition, to be
       !before this problematic region is reached.
       IF(potential_choice.eq.6) THEN
          yprime(1)=0.0d0
          yprime(2)=0.0d0
          RETURN
       ENDIF
       WRITE(*,*) 'MODPK: QUITTING'
       write(*,*) 'vparams: ', (vparams(i,:),i=1,size(vparams,1))
       if (.not.instreheat) write(*,*) 'N_pivot: ', N_pivot
       STOP
    END IF

    hubble=getH(p,delp)
    dhubble=getHdot(p,delp)

    !MULTIFIELD
    yprime(1 : size(y)/2) = delp
    yprime(size(y)/2+1 : size(y)) = -((3.0+dhubble/hubble)*delp+dVdphi(p)/hubble/hubble)
    !END MULTIFIELD 

    RETURN

  END SUBROUTINE bderivs


  ![ LP: ] Full y (back+mode matrix+tensor) derivatives y'=f(y)
  SUBROUTINE derivs(x,y,yprime)
    USE modpkparams
    USE internals
    USE potential, ONLY: pot, dVdphi, d2Vdphi2, getH, getHdot, getEps, getEta
    IMPLICIT NONE
    real(dp), INTENT(IN) :: x
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(IN) :: y
    COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: yprime

    ! background quantity
    real(dp) :: hubble, dhubble, a, epsilon, dotphi, eta
    real(dp) :: thetaN2, Vzz, Vz, grad_V  !! thetaN2 = (d\theta/dNe)^2
    real(dp), DIMENSION(num_inflaton) :: phi, delphi, Vp
    real(dp), DIMENSION(num_inflaton, num_inflaton) :: Cab, Vpp

    COMPLEX(DP), DIMENSION(num_inflaton,num_inflaton) :: psi, dpsi ! scalar ptb mode matrix
    COMPLEX(DP) :: v, dv                                        ! tensor perturbations
    COMPLEX(DP) :: u_zeta, du_zeta

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

    phi = dble(y(1:num_inflaton))
    delphi = dble(y(num_inflaton+1:2*num_inflaton))
    dotphi = sqrt(dot_product(delphi, delphi))

    ![ LP: ] Aliases to potential derivatives
    epsilon = getEps(phi, delphi)
    eta = geteta(phi, delphi)
    Vpp = d2Vdphi2(phi)
    Vp = dVdphi(phi)
    Vzz = dot_product(delphi, matmul(Vpp, delphi))/dotphi**2
    Vz = dot_product(Vp, delphi)/dotphi
    grad_V = sqrt(dot_product(Vp, Vp))

    IF(dot_product(delphi, delphi) .GT. 6.d0) THEN
       WRITE(*,*) 'MODPK: H is imaginary:',x,phi,delphi
       WRITE(*,*) 'MODPK: QUITTING'
       write(*,*) 'vparams: ', (vparams(i, :), i=1, size(vparams,1))
       if (.not.instreheat) write(*,*) 'N_pivot: ', N_pivot
       STOP
    END IF

    hubble=getH(phi,delphi)
    dhubble=getHdot(phi,delphi)
    a=EXP(x)

    ![ LP: ] Alias y's into real variable names
    ![ LP: ] NB: if using_q_superh, psi(i,j)-->q_ptb(i,j)
    !DEBUG
    psi = convert_hacked_vector_to_matrix(y(index_ptb_y:index_ptb_vel_y-1))
    dpsi = convert_hacked_vector_to_matrix(y(index_ptb_vel_y:index_tensor_y-1))
    !psi = y(index_ptb_y:index_ptb_vel_y-1)
    !dpsi = y(index_ptb_vel_y:index_tensor_y-1)

    v  = y(index_tensor_y)
    dv  = y(index_tensor_y+1)
    u_zeta = y(index_uzeta_y)
    du_zeta = y(index_uzeta_y+1)

    ![ LP: ] Build the mass matrix, Cab
    call build_mass_matrix(Cab)

    ![ LP: ] -----------------------------
    ![ LP: ] Set the RHS of y'(x) = f(y(x))
    ![ LP: ] -----------------------------

    ![ LP: ] background
    yprime(1:num_inflaton) = cmplx(delphi)
    yprime(num_inflaton+1:2*num_inflaton) =&
      cmplx(-((3.0e0_dp+dhubble/hubble)*delphi+dVdphi(phi)/hubble/hubble))

    ![ LP: ] ptb matrix
    yprime(index_ptb_y:index_ptb_vel_y-1) = &
      convert_matrix_to_hacked_vector(dpsi)
    if (using_q_superh) then
      yprime(index_ptb_vel_y:index_tensor_y-1) = convert_matrix_to_hacked_vector(&
        -(3.0e0_dp - epsilon)*dpsi - (k/a_init/a/hubble)**2*psi -&
        matmul(Cab, psi)/hubble**2)
    else
      yprime(index_ptb_vel_y:index_tensor_y-1) = convert_matrix_to_hacked_vector(&
        -(1.0e0_dp - epsilon)*dpsi - (k/a_init/a/hubble)**2*psi &
        + (2.0e0_dp - epsilon)*psi - matmul(Cab, psi)/hubble**2)
    end if

    ![ LP: ] tensors
    yprime(index_tensor_y) = dv
    if (using_q_superh) then
      yprime(index_tensor_y+1) = -(3.0e0_dp - epsilon)*dv - (k/a_init/a/hubble)**2*v
    else
      yprime(index_tensor_y+1) = -(1.0e0_dp - epsilon)*dv - &
        (k/a_init/a/hubble)**2*v + (2.0e0_dp - epsilon)*v
    end if

    ![ LP: ] adiabatic ptb
    yprime(index_uzeta_y) = du_zeta
    thetaN2 = (grad_V + Vz)*(grad_V - Vz)/(dotphi*hubble**2)**2 
    yprime(index_uzeta_y+1) = -(1.0e0_dp - epsilon)*du_zeta -&
      (k/a_init/a/hubble)**2*u_zeta &
      + (2.0e0_dp + 5.0e0_dp*epsilon - 2.0e0_dp*epsilon**2 + &
      2.0e0_dp*epsilon*eta + thetaN2 - Vzz/hubble**2)*u_zeta

    RETURN

    contains

      subroutine build_mass_matrix(mass_matrix)

        !DEBUG
        integer :: ii, jj

        real(dp), dimension(num_inflaton, num_inflaton), intent(out) :: mass_matrix

        if (potential_choice .eq. 7) then
          ! for exponential potential 7, Cab is excatly zero, need to set this in order to prevent numerical error
          mass_matrix = 0e0_dp
        else
           forall (i=1:num_inflaton, j=1:num_inflaton) & 
                mass_matrix(i,j) = Vpp(i,j) +  &
                (delphi(i)*Vp(j) + delphi(j)*Vp(i)) &
                + (3e0_dp-epsilon)*hubble**2 * delphi(i)*delphi(j)
        end if
        !DEBUG
        !write(314,*),log(a),((mass_matrix(ii,jj),ii=1,size(mass_matrix,1)),jj=1,size(mass_matrix,2))

      end subroutine build_mass_matrix

      !DEBUG
      !pure function dot(matrixA,matrixB)

      !  real(dp), dimension(:,:), intent(in) :: matrixA
      !  complex(dp), dimension(:,:), intent(in) :: matrixB
      !  integer :: i, j, k
      !  complex(dp), dimension(size(matrixA,1),size(matrixB,2)) :: dot

      !  dot=0e0_dp

      !  do i=1,size(matrixA,1);do j=1,size(matrixB,2)
      !    do k=1,size(matrixA,2)
      !      dot(i,j)=dot(i,j)+ matrixA(i,k)*matrixB(k,j)
      !    end do
      !  end do; end do



      !end function dot

  END SUBROUTINE derivs

  ![ LP: ] Full y (back+mode matrix(in Q at horizon cross)+tensor) derivatives y'=f(y)
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
    real(dp), PARAMETER :: SAFETY=0.9d0,PGROW=-0.2d0,PSHRNK=-0.25d0,&
         ERRCON=1.89e-4
    if (size(y)==size(dydx) .and. size(dydx)==size(yscal)) then
       ndum = size(y)
    else
       write(*,*) 'Wrong array sizes in rkqs'
       stop
    end if
    h=htry
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)
       errmax=maxval(abs(yerr(:)/yscal(:)))/eps
       if (errmax <= 1.0) exit
       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1d0*abs(h)),h)
       xnew=x+h
       if (xnew == x) then 
          print*, 'stepsize underflow in rkqs'
          ode_underflow = .true.
          return
       endif
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0d0*h
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
    real(dp), PARAMETER :: SAFETY=0.9d0,PGROW=-0.2d0,PSHRNK=-0.25d0,&
         ERRCON=1.89e-4

    real(dp) :: yerr_r(2*size(y)), yscal_r(2*size(yscal))

    if (size(y)==size(dydx) .and. size(dydx)==size(yscal)) then
       ndum = size(y)
    else
       write(*,*) 'Wrong array sizes in rkqs'
       stop
    end if
    h=htry
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)

       ! in doing error estimation and step size rescaling, we switch to real components
       yerr_r(1 : size(yerr)) = dble(yerr)
       yerr_r(size(yerr)+1 : 2*size(yerr)) = dble(yerr*(0,-1)) 

       yscal_r(1 : size(yscal)) = dble(yscal)
       yscal_r(size(yscal)+1 : 2*size(yscal)) = dble(yscal*(0,-1))  

       errmax=maxval(abs(yerr_r(:)/yscal_r(:)))/eps

       if (errmax <= 1.0) exit

       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1d0*abs(h)),h)
       xnew=x+h
       if (xnew == x) then 
          print*, 'stepsize underflow in rkqs'
          ode_underflow = .true.
          return
       endif
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0d0*h
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
    real(dp), PARAMETER :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,&
         A6=0.875d0,B21=0.2d0,B31=3.0d0/40.0d0,B32=9.0d0/40.0d0,&
         B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.0d0/54.0d0,&
         B52=2.5d0,B53=-70.0d0/27.0d0,B54=35.0d0/27.0d0,&
         B61=1631.0d0/55296.0d0,B62=175.0d0/512.0d0,&
         B63=575.0d0/13824.0d0,B64=44275.0d0/110592.0d0,&
         B65=253.0d0/4096.0d0,C1=37.0d0/378.0d0,&
         C3=250.0d0/621.0d0,C4=125.0d0/594.0d0,&
         C6=512.0d0/1771.0d0,DC1=C1-2825.0d0/27648.0d0,&
         DC3=C3-18575.0d0/48384.0d0,DC4=C4-13525.0d0/55296.0d0,&
         DC5=-277.0d0/14336.0d0,DC6=C6-0.25d0
    if (size(y)==size(dydx) .and. size(dydx)==size(yout) .and. size(yout)==size(yerr)) then
       ndum = size(y)
    else
       write(*,*) 'Wrong array sizes in rkck'
       stop
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
    real(dp), PARAMETER :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,&
         A6=0.875d0,B21=0.2d0,B31=3.0d0/40.0d0,B32=9.0d0/40.0d0,&
         B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.0d0/54.0d0,&
         B52=2.5d0,B53=-70.0d0/27.0d0,B54=35.0d0/27.0d0,&
         B61=1631.0d0/55296.0d0,B62=175.0d0/512.0d0,&
         B63=575.0d0/13824.0d0,B64=44275.0d0/110592.0d0,&
         B65=253.0d0/4096.0d0,C1=37.0d0/378.0d0,&
         C3=250.0d0/621.0d0,C4=125.0d0/594.0d0,&
         C6=512.0d0/1771.0d0,DC1=C1-2825.0d0/27648.0d0,&
         DC3=C3-18575.0d0/48384.0d0,DC4=C4-13525.0d0/55296.0d0,&
         DC5=-277.0d0/14336.0d0,DC6=C6-0.25d0
    if (size(y)==size(dydx) .and. size(dydx)==size(yout) .and. size(yout)==size(yerr)) then
       ndum = size(y)
    else
       write(*,*) 'Wrong array sizes in rkck'
       stop
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


  FUNCTION locate(xx,x)
    IMPLICIT NONE
    real(dp), DIMENSION(:), INTENT(IN) :: xx
    real(dp), INTENT(IN) :: x
    INTEGER*4 :: locate
    INTEGER*4 :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END FUNCTION locate

  ![ LP: ] Polynomial interpolation
  ![ LP: ] Given array XA and YA (of same length) and given a value X, returns
  !value Y such that if P(x) is a polynomial st P(XA_i)=YA_i, then Y=P(X) ---
  !and an error estimate DY
  SUBROUTINE polint(xa,ya,x,y,dy)
    IMPLICIT NONE
    real(dp), DIMENSION(:), INTENT(IN) :: xa,ya
    real(dp), INTENT(IN) :: x
    real(dp), INTENT(OUT) :: y,dy
    INTEGER*4 :: m,n,ns
    INTEGER*4, DIMENSION(1) :: imin
    real(dp), DIMENSION(size(xa)) :: c,d,den,ho,absho
    if (size(xa)==size(ya)) then
       n=size(xa)
    else
       write(*,*) 'Wrong array sizes in polint'
       stop
    end if
    c=ya
    d=ya
    ho=xa-x
    absho=abs(ho)
    imin=minloc(absho(:))
    ns=imin(1)
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0)) then
          write(*,*) 'polint: calculation failure'
          stop
       end if
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END SUBROUTINE polint

  !MULTIFIELD 
  SUBROUTINE array_polint(xa, ya, x, y, dy)
    IMPLICIT NONE

    real(dp), INTENT(IN) :: xa(:), ya(:,:)
    real(dp), INTENT(IN) :: x
    real(dp), INTENT(OUT) :: y(size(ya(:,1))),dy(size(ya(:,1)))
    INTEGER :: i

    do i = 1, size(ya(:,1))
       call polint(xa, ya(i,:), x, y(i), dy(i))
    end do

  END SUBROUTINE array_polint
  !END MULTIFIELD

  FUNCTION reallocate_rv(p,n)
    real(dp), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER*4, INTENT(IN) :: n
    INTEGER*4 :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)  !! allocate memeory of size n at new address to be returned
    if (ierr /= 0) then
       write(*,*) 'reallocate_rv: problem in attempt to allocate memory'
       stop
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
       write(*,*) 'reallocate_rm: problem in attempt to allocate memory'
       stop
    end if
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))= p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
  END FUNCTION reallocate_rm


  ![ LP: ] Take a "hacked" vector, i.e., a matrix B that is in the
  !form of a vector where the rows of B(i,1:N) --> first N terms in B(1:N), etc
  !Returns the matrix form. NB: the vector must be a representation of a square
  !matrix.

  function convert_hacked_vector_to_matrix(matrix_as_vector) &
    result(real_matrix)

    complex(dp), dimension(:), intent(in) :: matrix_as_vector
    complex(dp), dimension(:,:), allocatable :: real_matrix
    integer :: n, i, j

    n=int(sqrt(real(size(matrix_as_vector))))
    allocate(real_matrix(n,n))

    do i=1,n; do j=1,n
      real_matrix(i,j) = matrix_as_vector((i-1)*n+j)
    end do; end do

  end function convert_hacked_vector_to_matrix

  ![ LP: ] Take a matrix and "hack" it into a vector, i.e.,
  !a vector where the rows of B(i,1:N) --> first N terms in B(1:N), etc
  function convert_matrix_to_hacked_vector(matrix) &
    result(hacked_vector)

    complex(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), allocatable :: hacked_vector
    integer :: n, i, j

    n=size(matrix,1)*size(matrix,2)
    allocate(hacked_vector(n))

    do i=1,size(matrix,1); do j=1,size(matrix,2)
      hacked_vector((i-1)*size(matrix,2)+j) = matrix(i,j)
    end do; end do

  end function convert_matrix_to_hacked_vector

END MODULE modpk_utils
