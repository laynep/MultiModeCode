!For use in modpk_sampling.f90

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

    !Takes a uniform sample of implicit surface defined in the "constrained"
    !dimensions. Then moves the sample by electro-static charge until in
    !equilibrium.  Then spawns new points and repeats until full sample built.
    !Inspired largely from
    ! http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.27.2922&rep=rep1&type=pdf
    subroutine implicit_surface_sampling(sample,energy_scale, &
      numb_samples,&
      phi0_priors_min, phi0_priors_max, &
      dphi0_priors_min, dphi0_priors_max)

      real(dp), dimension(:,:), allocatable :: sample
      real(dp), intent(in) :: energy_scale
      integer, intent(in) :: numb_samples
      real(dp), dimension(num_inflaton), intent(in) :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max

      real(dp), dimension(:,:), allocatable :: temp_sample
      real(dp), dimension(num_inflaton) :: phi
      real(dp), dimension(2*num_inflaton) :: xvecti, xvectj, radius,&
        deltax
      real(dp) :: alpha, sigma, energy, tol

      integer :: i,j, firstvel
      logical, dimension(:), allocatable :: in_equil

      !Make the arrays
      if (allocated(sample)) then
        deallocate(sample)
      endif
      allocate(sample(numb_samples,2*num_inflaton))
      if (numb_samples<20) then
        allocate(temp_sample(numb_samples,2*num_inflaton))
      else
        allocate(temp_sample(numb_samples/10,2*num_inflaton))
      end if
      allocate(in_equil(size(temp_sample,1)))

      !Seed the IC surface
      !Use the first velocity (in time) as the constrained variable
      firstvel = num_inflaton+1
      do i=1,size(temp_sample,1)
        call unconstrained_ic(temp_sample(i,:), &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)
        phi = temp_sample(i,1:num_inflaton)
        temp_sample(i,firstvel) = constrain_first_vel(temp_sample(i,:),&
          sgn=1.0e0_dp)
      end do

      !Walk temporary points around surface by electrostatic-like pot
      !Uses grad desc
      in_equil = .false.
      tol = 1e-3_dp
      alpha=0.1e0_dp
      !DEBUG
      !Supposed to be sqrt(surface area/sample size)
      sigma=0.3e0_dp*sqrt(1.0e0_dp/size(temp_sample,1))
      do while (.not. all(in_equil))

        deltax=0e0_dp

        do i=1,size(temp_sample,1)
          xvecti = temp_sample(i,:)

          !DEBUG
          !O(N^2)
          do j=1,size(temp_sample,1)-1
            xvectj = temp_sample(j,:)
            radius = xvecti-xvectj
            energy = alpha*&
              exp(-1e0_dp*abs(sum(radius*radius))/(2.0e0_dp*sigma**2))

            !DEBUG
            deltax = deltax + radius*energy

            !Check if in equil
            if (sum(deltax*deltax)<tol) then
              in_equil(i) =.true.
            end if

          end do

          !Take step on constraint surface
          temp_sample(i,:)=temp_sample(i,:) + deltax
          temp_sample(i,firstvel) = constrain_first_vel(temp_sample(i,:),&
            temp_sample(i,firstvel)/abs(temp_sample(i,firstvel)),&
            randvel=.false.)

        end do


      end do


      contains

        function constrain_first_vel(y,sgn,randvel) result(constrain_y)

          real(dp), dimension(2*num_inflaton), intent(in) :: y
          real(dp), dimension(num_inflaton-1) :: vely
          real(dp) :: sgn
          logical, optional :: randvel
          real(dp) :: constrain_y
          real(dp) :: rand

          vely = y(num_inflaton+2:2*num_inflaton)

          constrain_y = sgn*sqrt(2.0e0_dp*(energy_scale**4-pot(phi))&
            - sum(vely*vely))

          !Pick first vel's sign at random
          if (present(randvel) .and. randvel) then
            call random_number(rand)
            rand = 2.0e0_dp*rand
            if (ceiling(rand)==2) then
              constrain_y=-1.0e0_dp*constrain_y
            endif
          end if

        end function constrain_first_vel

    end subroutine implicit_surface_sampling

!Some routines for post-processing the IC data
module modpk_postprocessing
  use modpkparams, only : dp, num_inflaton
  use modpk_sampling, only : ic_and_observables
  implicit none

  contains

    function calc_clumping_penalty(ics, sigma) result(E)

      type(ic_and_observables), dimension(:), intent(in) :: ics
      real(dp), intent(in) :: sigma
      !real(dp),dimension(size(ics),size(ics)) :: Eij
      real(dp),dimension(size(ics)) :: E
      real(dp),dimension(2*num_inflaton) :: xi,xj,rad

      integer :: i, j

      E=0e0_dp
      do i=1,size(ics)
        xi = ics(i)%ic(:)
        do j=1, size(ics)
          xj = ics(j)%ic(:)
          rad = xi-xj
          E(i) = E(i) + exp(-1e0_dp/2.0e0_dp/sigma**2/sum(rad*rad))
        end do
      end do

    end function calc_clumping_penalty

end module modpk_postprocessing
