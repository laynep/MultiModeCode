
module modpk_rng
  !Module that contains routines for random number generation.
  !Some of this has been taken (where indicated) from
  !RosettaCode.org and from other sites on the net.
  use modpkparams, only : dp
  implicit none

interface init_random_seed
	module procedure init_random_seed_parallel
	module procedure init_random_seed_serial
end interface

interface shuffle
	module procedure shuffle_dp_dimn
	module procedure shuffle_int_dimn
	module procedure shuffle_dp_1
	module procedure shuffle_int_1
end interface

interface normal
  module procedure normal_array
  module procedure normal_scalar
end interface

!Default to public.


contains

!Generate an array of n numbers sampled from a Gaussian normal distribution with mean and standard deviation given.
!Adapted from RosettaCode.org.
function normal_array(n,mean,std) result(normal)
  implicit none

  integer, intent(in) :: n
  integer :: i
  real(dp) :: normal(n), temp
  real(dp), intent(in) :: mean, std
  real(dp), parameter :: pi = 4.0*ATAN(1.0)

  !Get uniform distribution.
  call random_number(normal)

  !Now convert to normal distribution
  do i = 1, n-1, 2
    temp = std * SQRT(-2.0*LOG(normal(i))) * COS(2*pi*normal(i+1)) + mean
    normal(i+1) = std * SQRT(-2.0*LOG(normal(i))) * SIN(2*pi*normal(i+1)) + mean
    normal(i) = temp
  end do

end function normal_array

!Generate a scalar sampled from a Gaussian normal distribution with mean and standard deviation given.
!Adapted from RosettaCode.org.
function normal_scalar(mean,std) result(scalar)
  implicit none

  integer :: i
  integer, parameter :: n=4
  real(dp) :: normal(n), temp, scalar
  real(dp), intent(in) :: mean, std
  real(dp), parameter :: pi = 4.0*ATAN(1.0)

  !Get uniform distribution.
  call random_number(normal)

  !Now convert to normal distribution
  do i = 1, n-1, 2
    temp = std * SQRT(-2.0*LOG(normal(i))) * COS(2*pi*normal(i+1)) + mean
    normal(i+1) = std * SQRT(-2.0*LOG(normal(i))) * SIN(2*pi*normal(i+1)) + mean
    normal(i) = temp
  end do

  scalar = normal(1)

end function normal_scalar

!Generate a random seed based on the clock time.
!Adapted from GNU site.
subroutine init_random_seed_serial()
implicit none
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed

	call random_seed(size = n)
	allocate(seed(n))

	call system_clock(count=clock)

	seed = clock + 37* (/ (i - 1, i = 1, n) /)
	call random_seed(PUT = seed)

	deallocate(seed)
end subroutine init_random_seed_serial


!Generate a random seed based on the clock time --- For use in parallel.
subroutine init_random_seed_parallel(rank)
implicit none
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed
	integer, intent(in) :: rank

	call random_seed(size = n)
	allocate(seed(n))

	call system_clock(count=clock)

	seed = clock + 37*rank* (/ (i - 1, i = 1, n) /)
	call random_seed(PUT = seed)

	deallocate(seed)
end subroutine init_random_seed_parallel


!Subroutine to shuffle an array by Knuth Shuffle.
!Inspired by RosettaCode.org.
subroutine shuffle_int_1(a)
implicit none
	integer, intent(inout) :: a(:)
	integer :: i, randpos, temp
	real :: r

	!Count backwards
	do i = size(a), 2, -1
		!Get a random number.
		call random_number(r)
		randpos = int(r * i) + 1
		!Exchange the rows.
		temp = a(randpos)
		a(randpos) = a(i)
		a(i) = temp
	end do

end subroutine shuffle_int_1

subroutine shuffle_dp_dimn(a)
  implicit none

  real(dp), dimension(:,:), intent(inout) :: a
	real(dp), dimension(size(a,2)) :: temp
	integer :: i, hold
	real(dp) :: rand

	!Count backwards.
	do i=size(a,1), 1, -1
		!Generate a random int from 1-i.
		call random_number(rand)
		hold=int(rand*i)+1
		!Swap the ith row with the rand row.
		temp(:)=a(hold,:)
		a(hold,:)=a(i,:)
		a(i,:) = temp(:)
	end do

end subroutine shuffle_dp_dimn

subroutine shuffle_dp_1(a)
  implicit none

	real(dp), dimension(:), intent(inout) :: a
	real(dp) :: temp
	integer :: i, hold
	real(dp) :: rand

	!Count backwards.
	do i=size(a,1), 1, -1
		!Generate a random int from 1-i.
		call random_number(rand)
		hold=int(rand*i)+1
		!Swap the ith row with the rand row.
		temp=a(hold)
		a(hold)=a(i)
		a(i) = temp
	end do

end subroutine shuffle_dp_1

subroutine shuffle_int_dimn(a)
  implicit none

	integer, dimension(:,:), intent(inout) :: a
	integer, dimension(size(a,2)) :: temp
	integer :: i, hold
	real(dp) :: rand

	!Count backwards.
	do i=size(a,1), 1, -1
		!Generate a random int from 1-i.
		call random_number(rand)
		hold=int(rand*i)+1
		!Swap the ith row with the rand row.
		temp(:)=a(hold,:)
		a(hold,:)=a(i,:)
		a(i,:) = temp(:)
	end do

end subroutine shuffle_int_dimn

!Subroutine that takes two sets, shuffles them together n times, then divides them along the columns, returning same sized arrays that are a mix of both sets.
subroutine shuffle_cut(setA, setB, n)
  implicit none

	real(dp), dimension(:,:), intent(inout) :: setA, setB
	integer, optional, intent(in) :: n
	real(dp), dimension((size(setA,1)+size(setB,1)),size(setA,2)) :: work
	real(dp) ::  rand
	integer :: i, iend

	!Load the work function.
	do i=1,size(work,1)
		if (i .le. size(setA,1)) then
			work(i,:)=setA(i,:)
		else
			work(i,:)=setB(i-size(setA,1),:)
		end if
	end do

	!Shuffle the set.
	if (present(n)) then
		iend=n
	else
		iend=2
	end if
	do i=1,iend
		call shuffle(work)
	end do

	!Unload work.
	do i=1,size(work,1)
		if (i .le. size(setA,1)) then
			setA(i,:)=work(i,:)
		else
			setB(i-size(setA,1),:)=work(i,:)
		end if
	end do

end subroutine shuffle_cut


!Function to make a cluster of "n" points centered at "center" with given std "sigma."
function make_blob(n, center, sigma)
  implicit none

	real(dp), dimension(:), intent(in) :: center
	real(dp), optional, intent(in) :: sigma
	real(dp):: std
	integer, intent(in) :: n
	real(dp), dimension(:,:), allocatable :: make_blob
	integer :: i,j

	!Optional argument for standard deviation.
	if (present(sigma)) then
		std = sigma
	else
		std = 1e0_dp
	end if

	call init_random_seed()

	allocate(make_blob(n,size(center)))

	do j=1,size(center)
		make_blob(:,j) = normal(n,center(j),std)
	end do

end function make_blob

function rand_sign() result(sign)

  real(dp) :: sign
  real(dp) :: rand

  call random_number(rand)

  if (rand >0.5) then
    sign = -1.0e0_dp
  else
    sign = 1.0e0_dp
  end if

end function rand_sign


end module modpk_rng
