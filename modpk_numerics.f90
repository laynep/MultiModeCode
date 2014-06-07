!Some numerical routines
module modpk_numerics
  use modpkparams, only : dp
  implicit none

  contains

    !Newton's method for finding zeros to f(x)=0.
    !Modified from: http://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_newton.html
    function zero_finder(f, fp, x0, iters, debugging) result(x)

      ! Estimate the zero of f(x) using Newton's method. 
      ! Input:
      !   f:  the function to find a root of
      !   fp: function returning the derivative f'
      !   x0: the initial guess
      !   debug: logical, prints iterations if debug=.true.
      ! Returns:
      !   the estimate x satisfying f(x)=0 (assumes Newton converged!) 
      !   the number of iterations iters

      implicit none

      interface
        function f(x)
          use modpkparams
          implicit none
          real(dp), intent(in) :: x
          real(dp) :: f
        end function f

        function fp(x)
          use modpkparams
          implicit none
          real(dp), intent(in) :: x
          real(dp) :: fp
        end function fp
      end interface

      !real(dp), intent(out) :: x
      real(dp) :: x

      real(dp), intent(in) :: x0
      logical, intent(in), optional :: debugging
      logical :: debug
      integer, intent(out), optional :: iters
      integer, parameter :: maxiter = 200
      real(dp), parameter :: tol = 1.e-12_dp
      real(dp), parameter :: dtol = 1.e-5_dp

      ! Declare any local variables:
      real(dp) :: deltax, fx, fxprime
      integer :: k

      ! initial guess
      x = x0

      if (.not. present(debugging)) then
        debug=.false.
      else
        debug=debugging
      end if

      if (debug) then
          print 11, x
 11       format('Initial guess: x = ', e22.15)
      endif

      ! Newton iteration to find a zero of f(x) 

      do k=1,maxiter

          ! evaluate function and its derivative:
          fx = f(x)
          fxprime = fp(x)

          if (abs(fx) < tol) then
              exit  ! jump out of do loop
          endif

          !if (abs(fxprime) < tol) then
          !  print*, "ERROR: fxprime =", fxprime,"<", dtol
          !  print*, "ERROR: in fx/fxprime"
          !  stop
          !end if

          ! compute Newton increment x:
          deltax = fx/fxprime

          ! update x:
          x = x - deltax

          if (debug) then
              print 12, k,x
 12           format('After', i3, ' iterations, x = ', e22.15)
          endif

      enddo


      if (k > maxiter) then
          ! might not have converged

          fx = f(x)
          if (abs(fx) > tol) then
              print *, '*** Warning: has not yet converged'
              endif
          endif 

      ! number of iterations taken:
      if (present(iters)) iters = k-1

    end function zero_finder

    subroutine num_first_deriv(f, x, h, df)
      implicit none

      ! Estimate the (partial) derivs of f(x) using central diff
      ! Input:
      !   f:  the function to find derivs of
      !   x:  the point to differentiate around
      !   h:  the step-size for finite diff
      ! Returns:
      !   the estimate df of each (partial) deriv at x

      interface
        function f(x)
          use modpkparams
          implicit none
          real(dp), dimension(:), intent(in) :: x
          real(dp) :: f
        end function f
      end interface

      real(dp), dimension(:), intent(in) :: x, h
      real(dp), dimension(:), allocatable, intent(out) :: df

      integer :: i
      real(dp), dimension(size(x)) :: x_step

      allocate(df(size(x)))

      do i=1, size(x)
        x_step = 0e0_dp
        x_step(i) = h(i)
        df(i) = (0.5e0_dp/h(i))*&
          (f(x + x_step) - f(x - x_step))
      end do

    end subroutine num_first_deriv

    subroutine num_second_deriv(f, x, h, d2f)
      ! Estimate the second (partial) derivs of f(x) using central diff
      ! Input:
      !   f:  the function to find derivs of
      !   x:  the point to differentiate around
      !   h:  the step-size for finite diff (assumes same step each dimn)
      ! Returns:
      !   the estimate df of each (partial) deriv at x

      interface
        function f(x)
          use modpkparams
          implicit none
          real(dp), dimension(:), intent(in) :: x
          real(dp) :: f
        end function f
      end interface

      real(dp), dimension(:), intent(in) :: x, h
      real(dp), dimension(:,:), allocatable, intent(out) :: d2f
      real(dp), dimension(size(x)) :: x_step, y_step

      integer :: i, j

      allocate(d2f(size(x),size(x)))

      !Non-mixed derivs
      do i=1, size(x)
        x_step = 0e0_dp
        x_step(i) = h(i)
        d2f(i,i) = (1.0e0_dp/h(i)**2)*&
          (f(x + x_step) - 2e0_dp*f(x) + f(x - x_step))
      end do

      !Mixed derivs
      do i=1, size(x); do j=1, size(x)
        if (i==j) cycle

        x_step = 0e0_dp
        y_step = 0e0_dp
        x_step(i) = h(i)
        y_step(j) = h(j)
        d2f(i,j) = (0.25e0_dp/h(i)/h(j))*&
          (f(x + x_step + y_step) - f(x + x_step - y_step) &
          -f(x - x_step + y_step) + f(x - x_step - y_step))
      end do; end do

    end subroutine num_second_deriv

  pure FUNCTION locate(xx,x)
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

  !Find first entry in array that is (1) less than x and (2) the next entry is greater
  !than x
  pure integer function stupid_locate(array,x)
    real(dp), DIMENSION(:), intent(in) :: array
    real(dp), intent(in) :: x

    real(dp) :: check_point

    integer :: i

    check_point = 0
    do i=1,size(array)-1
      if (array(i) .le. x .and. array(i+1) > x) then
        check_point = i
        exit
      end if
    end do

    stupid_locate=check_point

  end function stupid_locate

  ! Polynomial interpolation
  ! Given array XA and YA (of same length) and given a value X, returns
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

  !Integrates a one-dimensional real function from x0 to xend using
  !the trapezoid rule with given number of steps
  function integrate_1d(funct, x0, xend, nsteps) result(area)
    implicit none

    real(dp) :: area
    real(dp), intent(in) :: x0, xend
    integer, intent(in) :: nsteps
    interface
      function funct(x)
        use modpkparams
        real(dp), intent(in) :: x
      end function funct
    end interface

    real(dp) :: xrange, dx, x_a, x_b
    integer :: ii

    xrange = xend - x0
    dx = xrange/real(nsteps,dp)

    area = 0e0_dp
    x_a = x0
    do ii=1,nsteps
      x_b = x_a + dx
      area = area +&
        (dx)*(funct(x_a) + funct(x_b))/2.0e0_dp
      x_a = x_b
    end do

  end function integrate_1d

end module modpk_numerics
