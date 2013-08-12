!Implements various IC sampling techniques
module modpk_icsampling
  use potential
  use modpkparams
  implicit none

  integer, parameter :: test_samp=1, eqen_samp=2

  contains

    !Grab a new IC according to the sampling_techn variable
    subroutine get_ic(sampling_techn, energy_scale, numb_samples)

      integer, intent(in) :: sampling_techn, numb_samples
      real(dp), intent(in), optional :: energy_scale

      real(dp), dimension(num_inflaton) :: y_background
      real(dp), dimension(:), allocatable :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max
      integer :: u, i

      real(dp), dimension(:,:), allocatable :: sample

      real(dp) :: penalty_fact

      namelist /priors/ phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max, &
        penalty_fact

      if (sampling_techn == eqen_samp) then

        allocate(phi0_priors_min(num_inflaton))
        allocate(phi0_priors_max(num_inflaton))
        allocate(dphi0_priors_min(num_inflaton))
        allocate(dphi0_priors_max(num_inflaton))

        if (slowroll_start) then
          print*, "Setting slowroll_start = .false."
          slowroll_start = .false.
        end if

        !Read phi0 priors from file
	      open(newunit=u, file="parameters_multimodecode.txt", &
          status="old", delim = "apostrophe")
        read(unit=u, nml=priors)
        close(u)


        !Not ready for primetime
        !call penalized_constrained_montecarlo(sample,energy_scale, &
        !  numb_samples, penalty_fact,&
        !  phi0_priors_min, phi0_priors_max, &
        !  dphi0_priors_min, dphi0_priors_max)

        call alt_eqen_ic(y_background, energy_scale, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)

      else
        print*, "ERROR: Sampling technique hasn't been implemented."
        stop
      end if

    end subroutine get_ic

    !MCMC-style sampling that penalizes points for clustering together by giving
    !an inverse weighting to points that, if charged, would experience a large
    !force perpendicular to the "proposal" direction, i.e., points that are
    !trying to clump unnecessarily.
    subroutine penalized_constrained_montecarlo(sample,energy_scale, &
          numb_samples, penalty_fact, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max, &
          !optionals
          delE_tol_in, coulomb_samples, prop_variance)

      use modpk_rng, only : normal

      real(dp), dimension(:), intent(in) :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max
      real(dp), intent(in) :: energy_scale, penalty_fact
      integer, intent(in) :: numb_samples
      real(dp), dimension(:,:), allocatable, intent(out) :: sample

      integer, intent(in), optional :: coulomb_samples
      real(dp), intent(in), optional :: delE_tol_in
      real(dp), intent(in), dimension(2*num_inflaton), optional :: prop_variance

      logical, dimension(size(sample,1)) :: off_shell_sample
      real(dp) :: delE_tol, rand
      real(dp), dimension(2*num_inflaton) :: variance, proposal
      integer :: i, live_point, rand_index
      real(dp), dimension(:,:), allocatable :: charged_sample
      real(dp), dimension(2*num_inflaton) :: live_vect, penalty_vect, radius
      real(dp) :: old_measure, new_measure, accept_ratio


      call build_unconstrained_sample(sample,numb_samples, energy_scale,&
        phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max)

      !Energy off-set tolerance
      if (present(delE_tol_in)) then
        delE_tol = delE_tol_in
      else
        delE_tol = 1e-10_dp
      end if

      !Number of points to use for the "Coulomb" penalty
      if (present(coulomb_samples)) then
        allocate(charged_sample(coulomb_samples,2*num_inflaton))
      else
        if (size(sample,1) > 10) then
          allocate(charged_sample(10,2*num_inflaton))
        else
          allocate(charged_sample(1,2*num_inflaton))
        end if
      end if

      !Proposal Gaussian's variance
      if (present(prop_variance)) then
        variance = prop_variance
      else
        variance = 1.0e0_dp
      end if

      !Do the sampling...

      off_shell_sample = sample(:,1) < delE_tol
      do while (any(off_shell_sample))

        !Find point in sample w/smallest energy_measure
        !DEBUG
        live_point = 1
        live_vect = sample(live_point,2:2*num_inflaton+1)

        !Get a proposal vector from Gaussian w/zero mean and var=prop_variance
        do i=1,size(proposal)
          proposal(i) = normal(0e0_dp,sqrt(variance(i)))
        end do

        !Get M=charged_sample number of points at random from the sample
        do i=1, size(charged_sample,1)
          do
            call random_number(rand)
            rand_index = ceiling(size(sample,1)*rand)
            if (rand_index/=live_point) exit
          end do
          charged_sample(i,:) = sample(rand_index,2:size(sample,2))
        end do

        !Calculate "penalty" vector on live point from charge points
        !Inverse square on phase space distance
        !Force penalty vector to have norm(penalty)/norm(proposal)=penalty_fact
        penalty_vect = 0e0_dp
        do i=1,size(charged_sample,1)
          radius(:) = live_vect(:) - charged_sample(i,:)
          penalty_vect = penalty_vect+ (1e0_dp/dot_product(radius,radius))
        end do
        penalty_vect = penalty_vect* &
          penalty_fact*norm(proposal)/norm(penalty_vect)

        !DEBUG
        !print*, "penalty_vect", penalty_vect
        !print*, "proposal", proposal
        !print*, "live", live_vect
        !print*, "diff", proposal - penalty_vect
        !print*, "norm", norm(penalty_vect)/norm(proposal)
        !print*, "penalty_fact", penalty_fact
        !stop


        !Add the penalty and proposal vectors and take trial step
        proposal = proposal + penalty_vect
        proposal = proposal + live_vect

        !Calculate energy measure and acceptance ratio
        old_measure = energy_measure(live_vect,energy_scale)
        new_measure = energy_measure(proposal,energy_scale)
!DEBUG
print*, "proposal", proposal
print*, "live", live_vect
print*, "old", old_measure
print*, "new", new_measure
stop
        if (old_measure<1e-10_dp .and. new_measure<1e-10_dp) then
          !Just random walk if the measures are too low
          accept_ratio = 1.0e0_dp
        else
          accept_ratio=min(1.0e0_dp,new_measure/old_measure)
        end if

        !Accept if rand > min(1,accept_ratio)
        call random_number(rand)
        if (accept_ratio>=rand) then
!DEBUG
print*, "ACCEPT"
          sample(live_point,2:2*num_inflaton+1) = proposal
          sample(live_point,1) = new_measure
          !Check for off_shell tolerance
          if (new_measure < delE_tol) then
            off_shell_sample(live_point)=.false.
          end if
        end if

!DEBUG
print*, new_measure

      end do


    end subroutine penalized_constrained_montecarlo

    subroutine build_unconstrained_sample(sample,numb_samples,&
          energy_scale, phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)

      real(dp), dimension(:), intent(in) :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max
      integer :: u, i, j
      integer, intent(in) :: numb_samples
      real(dp), intent(in) :: energy_scale

      real(dp), dimension(:,:), allocatable, intent(out) :: sample
      integer, parameter :: excess = 1
      real(dp) :: rand

      !Array w/each row a point in phase space + an energy
      !sample(:,1)=energy_measure(:)
      allocate(sample(numb_samples*excess,1+2*num_inflaton))

      sample=1e0_dp
      do i=1,size(sample,1)
        do j=2,size(sample,2)
          call random_number(rand)
          if (j>num_inflaton+1) then
            sample(i,j) = rand*(dphi0_priors_max(j-1-num_inflaton)- &
              dphi0_priors_min(j-1-num_inflaton))+&
              dphi0_priors_min(j-1-num_inflaton)
          else
            sample(i,j) = rand*(phi0_priors_max(j-1)- &
              phi0_priors_min(j-1))+&
              phi0_priors_min(j-1)
          end if
        end do
      end do

      do i=1,size(sample,1)
        sample(i,1) = energy_measure(sample(i,2:2*num_inflaton+1), &
          energy_scale)
      end do


    end subroutine build_unconstrained_sample

    function energy_measure(vect,energy_scale) result(delE)

      real(dp), dimension(2*num_inflaton), intent(in) :: vect
      real(dp), dimension(num_inflaton):: fields, veloc
      real(dp), intent(in) :: energy_scale
      real(dp) :: delE
      real(dp) :: current_energy

      fields = vect(1:num_inflaton)
      veloc = vect(num_inflaton+1:2*num_inflaton)

      !Potential
      current_energy = pot(fields)

      !V + kinetic
      current_energy = current_energy + &
        0.5e0_dp*sum(veloc*veloc)

      delE = abs(energy_scale - current_energy)

      delE = exp(-delE)

    end function energy_measure







    !Require the background ICs to start with the same energy.
    !Alternates selecting the fields and vels 
    !Assumes canonical kinetic term
    subroutine alt_eqen_ic(y, energy_scale, &
        phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max)

      real(dp), dimension(num_inflaton*2), intent(out) :: y
      real(dp), dimension(num_inflaton), intent(in) :: phi0_priors_min, &
        phi0_priors_max, dphi0_priors_min, dphi0_priors_max
      real(dp), intent(in) :: energy_scale

      real(dp), dimension(num_inflaton*2) :: delta_prior

      real(dp), dimension(num_inflaton) :: phi, dphi

	    real(dp) :: rand, rho_not_alloc
	    integer :: param_constr, a, b, i, ll
      integer, parameter :: maxtry=1000


      delta_prior(1:num_inflaton) = phi0_priors_max-phi0_priors_min
      delta_prior(num_inflaton+1:2*num_inflaton) = &
        dphi0_priors_max-dphi0_priors_min

	    do ll=1,maxtry

        if(ll==maxtry) then
          print*, "Couldn't find an IC with energy=energy_scale in the max"
          print*, "number of tries."
          stop
        end if

        !Choose y from unif prior
        do i=1,size(y)
          call random_number(rand)
          y(i) = rand*delta_prior(i)
          if (i > num_inflaton) then
            y(i) = y(i) + dphi0_priors_min(i)
          else
            y(i) = y(i) + phi0_priors_min(i)
          end if
        end do

        !Randomly pick a dimn to set by energy constraint
        call random_number(rand)
        param_constr = ceiling(rand*size(y))

        !Set y(constraint) to minim energy value
        y(param_constr) = min_fixed_energy(param_constr)


        phi = y(1:num_inflaton)
        dphi = y(num_inflaton+1:2*num_inflaton)


	    	!Energy density not allocated
	    	rho_not_alloc = energy_scale**4 - pot(phi) -&
          0.5e0_dp*sum(dphi*dphi)

!DEBUG
!print*, "constraining", param_constr
!print*, "energy remaining", rho_not_alloc
!print*, "y", y

	    	if (rho_not_alloc<0) then
!DEBUG
print*, "Off shell ------------- cycling"
print*, "E diff", rho_not_alloc
print*, "KE", 0.5e0_dp*sum(dphi*dphi)
print*, "PE", pot(phi)
print*, "energy_scale", energy_scale
print*, "fields", phi
print*, "vels", dphi
          cycle
        else
          exit
        end if

	    end do
		
	    !Set y(constraint) to fix total energy
      call set_y_by_energy_constraint(y, param_constr, energy_scale)


      !DEBUG
      !print*, "stopping at end of alt_eqen_ic"
      !stop

    end subroutine alt_eqen_ic

    !Find min of potential w/all but one param constrained by method of
    !gradient descent.

    !---Model dependent---
    !Can put in an analytical soln to this problem here, rather than
    !relying on gradient descent
    function min_fixed_energy(constraint) result(minim)

      integer, intent(in) :: constraint
      real(dp) :: minim

      if (constraint>num_inflaton) then
        !Kinetic var; assume canonical KE
        minim=0e0_dp
      else

        !DEBUG
        if (potential_choice /= 1) then
          print*, "Just setting the IC to 0 -----"
        end if
        minim = 0e0_dp

      end if

    end function min_fixed_energy

    !Given only one dof in Y(constraint) and the fact that the energy
    !constraint equation has one or more real solutions, solve for the
    !remaining Y(constraint).
    subroutine set_y_by_energy_constraint(y, constraint, energy_scale)

      !use modpk_utils, only : laguerre_root_finder
      use modpk_rng, only : rand_sign

      real(dp), dimension(:), intent(inout) :: y
      real(dp), intent(in) :: energy_scale
      integer, intent(in) :: constraint

      real(dp), dimension(num_inflaton) :: phi, dphi
      integer :: poly_degree
      complex(dp) :: root
      real(dp) :: rand
      real(dp) :: m2_V(num_inflaton), energy_remaining

      !NB: y(constraint) already set to its minim value

      !Force y(constraint)=0 in case this isn't the min
      y(constraint) = 0e0_dp

      phi = y(1:num_inflaton)
      dphi = y(num_inflaton+1:2*num_inflaton)

      !w/y(constraint)=0
      energy_remaining = energy_scale**4 - pot(phi) - &
        0.5e0_dp*sum(dphi*dphi)

      if (constraint>num_inflaton) then

        !Constraining a velocity
        !Assume canonical KE
        !Choose sqrt's sign randomly
        y(constraint) = rand_sign()*sqrt(2e0_dp* &
          (energy_remaining))

        !return

      else

        !Constraining a field and have to solve the potential
        !Try to do this by hand.
        if (potential_choice==1) then
          poly_degree = 2

          m2_V = 10.d0**(vparams(1,:))

          y(constraint) = rand_sign()*sqrt((2e0_dp/m2_V(constraint))*&
            energy_remaining)


        else
          !DEBUG
          print*, "Specify the degree of the polynomial that makes"
          print*, "the energy constraint equation or put the"
          print*, "algebraic constraint in by-hand."
          stop
        end if

        !root = laguerre_root_finder(poly_degree)

      end if


          !DEBUG
          print*, "y(constraint)", y(constraint), "constraint", constraint
          if (constraint<num_inflaton) print*, "const mass", m2_V(constraint)
          print*, "KE", .5e0_dp*sum(dphi*dphi)
          print*, "PE", .5e0_dp*sum(m2_V*phi*phi)
          print*, "sign", rand_sign()
          !stop


    end subroutine set_y_by_energy_constraint


end module modpk_icsampling