!Some numerical routines
module modpk_numerics
  use modpkparams, only : dp, num_inflaton
  implicit none

  interface num_first_deriv
  	module procedure num_first_deriv
  	module procedure num_first_deriv_vectorfunct
  end interface

  interface num_second_deriv
  	module procedure num_second_deriv
  	module procedure num_second_deriv_vectorfunct
  end interface

  contains

    !**********************************************************
    !Heapsorts an array based on the first column only.

    !Adapted from Numerical Recipes pg 231.

    pure subroutine heapsort(table)
      implicit none

    	real(dp), dimension(:,:), intent(inout) :: table
    	integer :: n, l, ir, i, j, i_1, i_2
    	real(dp), dimension(size(table,2)) :: rra	!row temporary placeholder.

      if (size(table,1)==1) return

    	rra=0_dp
    	n=size(table,1)
    	l = (n/2)+1	!note the integer division.
    	ir = n
    do1:	do
    		if(l > 1) then
    			l = l-1
    			call vect_eq_tablerow_d(rra,l,table)
    		else
    			call vect_eq_tablerow_d(rra,ir,table)
    			call row_equal_d(ir,1,table)
    			ir = ir -1
    			if(ir==1)then
    				do i_1=1,size(table,2)
    					table(1,i_1) = rra(i_1)
    				end do
    				return
    			end if
    		end if
    		i = l
    		j = l+l
    do2:		do while(j <= ir)
    			if(j < ir) then
    				if(table(j,1) < table(j+1,1)) then
    					j = j+1
    				end if
    			end if
    			if(rra(1) < table(j,1)) then
    				call row_equal_d(i,j,table)
    				i = j
    				j =j+j
    			else
    				j = ir + 1
    			end if
    		end do do2
    		do i_2=1,size(table,2)
    			table(i,i_2) = rra(i_2)
    		end do
    	end do do1


    end subroutine heapsort

    !this subroutine makes a vector (of rank equal to the number of columns in table) equal to the ith row in a table.
    pure subroutine vect_eq_tablerow_d(vect,i,table)
    implicit none

    	real(dp), dimension(:), intent(inout) :: vect
    	real(dp), dimension(:,:), intent(in) :: table
    	integer, intent(in) :: i

    	vect(:) = table(i,:)

    end subroutine vect_eq_tablerow_d



    !this subroutine changes the ith row of a table to equal the jth row.
    pure subroutine row_equal_d(i,j,table)
    implicit none

    	real(dp), dimension(:,:), intent(inout) :: table
    	integer, intent(in) :: i,j

    	table(i,:) = table(j,:)

    end subroutine row_equal_d

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
      integer, parameter :: maxiter = 1000
      real(dp), parameter :: tol = 1.e-6_dp !Need this to be pretty conservative

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
              print *, '*** Warning: has not yet converged', abs(fx) - tol
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

    subroutine num_first_deriv_vectorfunct(f, x, h, df)
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
          real(dp), intent(in) :: x
          real(dp), dimension(num_inflaton) :: f
        end function f
      end interface

      real(dp), intent(in) :: x, h
      real(dp), dimension(:), allocatable, intent(out) :: df

      integer :: i

      allocate(df(num_inflaton))

      df = (0.5e0_dp/h)*(f(x + h) - f(x - h))


    end subroutine num_first_deriv_vectorfunct

    subroutine num_second_deriv_vectorfunct(f, x, h, d2f)
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
          real(dp), intent(in) :: x
          real(dp), dimension(num_inflaton) :: f
        end function f
      end interface

      real(dp), intent(in) :: x, h
      real(dp), dimension(:,:), allocatable, intent(out) :: d2f

      integer :: i, j

      allocate(d2f(num_inflaton,1))

      d2f(:,1) = (1.0e0_dp/h**2)*&
        (f(x + h) - 2e0_dp*f(x) + f(x - h))

    end subroutine num_second_deriv_vectorfunct

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

    !Given array XX and a value X, returns value XOUT st X is between XX(XOUT)
    !and XX(XOUT+1)
    subroutine hunt(xx,x,jlo)
      implicit none
      real(dp), dimension(:) :: xx
      logical ascnd
      integer :: n, jlo, inc, jhi, jm
      real(dp) :: x

      n = size(xx)

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
    end subroutine

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
    !Scalar interpolating function with vector argument
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
      dx = xrange/dble(nsteps)

      area = 0e0_dp
      x_a = x0
      do ii=1,nsteps
        x_b = x_a + dx
        area = area +&
          (dx)*(funct(x_a) + funct(x_b))/2.0e0_dp
        x_a = x_b
      end do

    end function integrate_1d

    !Build a histogram over an N-dimensional real dataset
    pure function histogram_Nd(dataset, real_binsize, method, norm) &
        result(hist)
      implicit none

      real(dp), dimension(:,:), intent(in) :: dataset
      real(dp), dimension(:), intent(in), optional :: real_binsize
      integer, intent(in), optional :: method
      real(dp), dimension(:,:), allocatable :: hist
      integer, intent(in), optional :: norm

      integer(dp) :: ndimns

      real(dp), dimension(size(dataset,2)) :: binsize
      integer(dp), dimension(size(dataset,2)) :: numb_bins
      real(dp), dimension(size(dataset,2)) :: data_max, data_min
      real(dp), dimension(size(dataset,2)) :: bin, bin_max
      integer(dp) :: ndata
      integer :: ii, jj, kk


      integer, dimension(:), allocatable :: bincount

      ndimns = size(dataset,2)

      !Get the binsize
      if (.not. present(real_binsize)) then
        binsize = determine_binsize(dataset, method)
      else
        binsize = real_binsize
      end if


      !Data stats
      ndata = size(dataset,1)
      do ii=1, size(data_max)
        data_max(ii) = maxval(dataset(:,ii))
        data_min(ii) = minval(dataset(:,ii))
      end do

      numb_bins = floor((data_max-data_min)/binsize)

      !NB: Some binsize techn give estimates of binsize, not bin #
      !Minor rescale of binsize to make sure fits dataset perfectly
      binsize = (data_max-data_min)/dble(numb_bins)


      !Make the histogram array
      allocate(hist(product(numb_bins),1+ndimns))
      hist =0e0_dp

      !Load bin positions
      call make_grid(numb_bins, data_max, data_min, hist)

      !Count the number of datapoints in each bin
      allocate(bincount(product(numb_bins)))
      do ii=1,size(hist,1)
        bin = hist(ii,1:ndimns)
        bin_max = bin+binsize
        bincount(ii)=0
        do jj=1,size(dataset,1)
          if ( (all(dataset(jj,:) .ge. bin) &
            .and. (all(dataset(jj,:) .le. bin_max)))) then
            bincount(ii) = bincount(ii) +1
          end if
        end do
      end do

      !Load hist
      if (present(norm)) then
        if (norm==1) then
          !Normalize to PDF (int P dx =1)
          hist(:,ndimns+1) = dble(bincount) / dble(sum(bincount)) /&
            product(binsize)
        else if(norm==2) then
          !Normalize to PMF (sum P_i =1)
          hist(:,ndimns+1) = dble(bincount) / dble(sum(bincount))
        else
          hist(:,ndimns+1) = bincount(:)
        end if
      else
        hist(:,ndimns+1) = bincount(:)
      end if


    end function histogram_Nd

    !Function to determine the binsize for an N-dimensional histogram
    pure function determine_binsize(dataset, method) result(binsize)
      implicit none

      !Method:
      !  1 = Freedman-Diaconis
      !  2 = Sturges
      !  3 = Scott normal reference rule
      !  Default = Square-root

      real(dp), dimension(:,:), intent(in) :: dataset
      real(dp), dimension(size(dataset,2)) :: binsize
      integer, intent(in) :: method

      real(dp) :: ndata, IQR
      real(dp), dimension(size(dataset,2)) :: data_max, data_min, &
        numb_bins
      integer :: n_sub

      real(dp), dimension(size(dataset,1),1) :: data_vect

      integer :: ii

      !Data stats
      ndata = size(dataset,1)
      do ii=1, size(data_max)
        data_max(ii) = maxval(dataset(:,ii))
        data_min(ii) = minval(dataset(:,ii))
      end do

      !Find binsize using named method
      select case(method)
      case(1)
        !Freedman-Diaconis binning

        !Determine quartiles
        do ii=1, size(dataset,2)
          n_sub = int(ndata/4)
          data_vect(:,1) = dataset(:,ii) !Not ideal for large datasets
          call heapsort(data_vect)
          IQR = data_vect(3*n_sub,1) - data_vect(n_sub,1)
          binsize(ii) = 2.0e0_dp*IQR*ndata**(-1.0e0_dp/3.0e0_dp)
        end do

      case(2)
        !Sturges binning
        numb_bins = ceiling(log(ndata)/log(2.0e0_dp) + 1)

        binsize = (data_max-data_min)/numb_bins

      case(3)
        !Scott normal reference binning
        binsize = (3.5e0_dp/ndata**(1.0e0_dp/3.0e0_dp))*&
          stand_dev(dataset)

      case default
        !Square root binning
        numb_bins = ceiling(sqrt(ndata))

        binsize = (data_max-data_min)/numb_bins

      end select


    end function determine_binsize

    !Makes a grid of points between data_max and data_min where there are
    !(numb_bins) of bins of size binsize
    pure subroutine make_grid(numb_bins,  &
        data_max, data_min, grid)
      implicit none

      integer(dp), dimension(:), intent(in) :: numb_bins
      real(dp), dimension(size(numb_bins)) :: binsize
      real(dp), dimension(:,:), intent(out) :: grid
      integer :: ii, jj, kk
      integer :: ndimns

      integer :: chunk, nchunks, nvects

      real(dp), dimension(:), allocatable :: vect
      integer :: base
      real(dp), dimension(size(numb_bins)), intent(in) :: data_max, data_min

      ndimns = size(numb_bins)

      binsize = (data_max-data_min)/dble(numb_bins)


      do ii=ndimns,1,-1
        chunk=1
        if (ii<ndimns) then
          do jj=ii, ndimns-1
            chunk =  chunk*numb_bins(jj+1)
          end do
        end if

        if (allocated(vect)) deallocate(vect)
        allocate(vect(numb_bins(ii)))
        do jj=0, numb_bins(ii)-1
          vect(jj+1) = dble(jj)
        end do
        nvects = product(numb_bins)/numb_bins(ii)/chunk

        base = 0
        do kk=1, nvects
          do jj=1,numb_bins(ii)
            grid((jj-1)*chunk+1+base:jj*chunk+base ,ii) = &
              vect(jj)*binsize(ii) + data_min(ii)
          end do
          base = base + numb_bins(ii)*chunk
        end do

      end do

    end subroutine

    !Simple arithmetic mean
    pure function mean(dataset, weight)
      implicit none

      real(dp), dimension(:,:), intent(in) :: dataset
      real(dp), dimension(size(dataset,1)), intent(in), optional :: weight
      real(dp), dimension(size(dataset,2)) :: mean

      integer :: i

      if (present(weight)) then
        do i=1,size(mean)
          mean(i) = sum(weight(:)*dataset(:,i))
        end do
      else
        do i=1,size(mean)
          mean(i) = sum(dataset(:,i))/real(size(dataset,1),dp)
        end do
      end if

    end function mean

    !Standard deviation
    pure function stand_dev(dataset, weight)
      implicit none

      real(dp), dimension(:,:), intent(in) :: dataset
      real(dp), dimension(size(dataset,1)), intent(in), optional :: weight
      real(dp), dimension(size(dataset,2)) :: stand_dev

      integer :: i

      if (present(weight)) then
        stand_dev = sqrt(mean(dataset**2, weight) - mean(dataset,weight)**2)
      else
        stand_dev = sqrt(mean(dataset**2) - mean(dataset)**2)
      end if

    end function stand_dev



end module modpk_numerics
