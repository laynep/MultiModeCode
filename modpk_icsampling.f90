!Implements various IC sampling techniques:
!reg_samp = Take phi_init0 from param file and set dphi in SR
!eqen_samp = Sample an iso-E surface by alternating which dimn is set by the
!energy constraint
!slowroll_samp = Sample phi0, set dphi in SR

module modpk_icsampling
  use potential, only : norm, pot, getH
  use modpkparams, only : dp, slowroll_start, num_inflaton, &
    potential_choice, vparams
  use internals, only : pi
  implicit none

  integer, parameter :: reg_samp=1, eqen_samp=2, slowroll_samp=3, &
    fromfile_samp=4, mass_loop_samp=5
  integer, parameter :: bad_ic=6
  integer :: sampling_techn
  real(dp) :: penalty_fact
  logical :: save_iso_N=.false.
  real(dp) :: N_iso_ref
  real(dp), dimension(:), allocatable :: phi_iso_N, dphi_iso_N

  !Type to save the ICs and observs. Add new observables here
  type :: ic_and_observables
    real(dp), dimension(:), allocatable :: ic
    real(dp) :: As, ns, r, nt, alpha_s
    real(dp) :: n_iso, n_pnad,n_ent
    real(dp) :: A_iso, A_pnad, A_ent, A_bundle
    contains
      procedure :: printout => ic_print_observables
      procedure :: load_observables => ic_load_observables
  end type ic_and_observables

  integer, private :: icfile


  contains


    !Grab a new IC according to the sampling_techn variable
    subroutine get_ic(phi0, dphi0,sampling_techn, &
        priors_min, priors_max, &
         numb_samples,energy_scale)

      integer, intent(in) :: sampling_techn, numb_samples
      real(dp), intent(in), optional :: energy_scale
      real(dp), dimension(num_inflaton), intent(out) :: phi0, dphi0
      real(dp), dimension(2,num_inflaton), intent(in) :: priors_min, &
        priors_max

      real(dp), dimension(2*num_inflaton) :: y_background
      real(dp), dimension(num_inflaton) :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max

      real(dp), dimension(:,:), allocatable :: sample
      real(dp) :: H, rho
      integer :: i,j

      phi0_priors_max=priors_max(1,:)
      phi0_priors_min=priors_min(1,:)
      dphi0_priors_max=priors_max(2,:)
      dphi0_priors_min=priors_min(2,:)

      if (sampling_techn == eqen_samp) then

        !Tell solver starting out of SR, even if in SR
        !in case it's only transient, e.g., starting with dphi=0
        !but not on "attractor"
        if (slowroll_start) slowroll_start = .false.

        !Not ready for primetime
        !call implicit_surface_sampling(sample,energy_scale, &
        !  numb_samples,&
        !  phi0_priors_min, phi0_priors_max, &
        !  dphi0_priors_min, dphi0_priors_max)

        !Not ready for primetime
        !call penalized_constrained_montecarlo(sample,energy_scale, &
        !  numb_samples, penalty_fact,&
        !  phi0_priors_min, phi0_priors_max, &
        !  dphi0_priors_min, dphi0_priors_max)

        call eqen_ic(y_background, energy_scale, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max, &
          velconst=.true.)

      else if (sampling_techn == slowroll_samp) then

        !Force solver to start in SR
        slowroll_start = .true.

        !Will override what vels it chooses later
        call unconstrained_ic(y_background, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)

      else if (sampling_techn == fromfile_samp) then

        !Get IC from a file
#ifdef MPI
!DEBUG
print*, "sampling_techn==fromfile_samp doesn't work with MPI"
stop
#endif

!DEBUG
print*, "sampling_techn==fromfile_samp not yet working..."
stop
        open(newunit=icfile, file="initial_conditions.txt",&
          form="unformatted",status="old")
        !DEBUG
        print*, "icfile", icfile
        do j=1,10
          read(icfile), (y_background(i),i=1,size(y_background))
        print*, "size", size(y_background)
        end do
        !DEBUG
        print*, "INITIAL COND", y_background
        stop

      else if (sampling_techn==mass_loop_samp) then
        !Get new vparams
        call mass_spectrum_prior(vparams)

        !Set IC with zero velocity
        call equal_displacement_ic(y_background, energy_scale)

      else
        print*, "ERROR: Sampling technique hasn't been implemented."
        stop
      end if

      !Load initial vals from sample
      phi0 = y_background(1:num_inflaton)
      dphi0 = y_background(num_inflaton+1:2*num_inflaton)

      !NB: Convert dphidt-->dphidN
      rho=0.5e0_dp*sum(dphi0*dphi0)+pot(phi0)
      H=sqrt(rho/3.0e0_dp)
      dphi0 = H*y_background(num_inflaton+1:2*num_inflaton)

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
    subroutine eqen_ic(y, energy_scale, &
        phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max,velconst)

      real(dp), dimension(num_inflaton*2), intent(out) :: y
      real(dp), dimension(num_inflaton), intent(in) :: phi0_priors_min, &
        phi0_priors_max, dphi0_priors_min, dphi0_priors_max
      real(dp), intent(in) :: energy_scale
      logical, intent(in) :: velconst

      real(dp), dimension(num_inflaton) :: phi, dphi

	    real(dp) :: rand, rho_not_alloc, KE

	    integer :: param_constr, i, ll, j
      integer, parameter :: maxtry=200


	    do ll=1,maxtry

        call unconstrained_ic(y, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)

        !Pick a dimn to set by energy constraint
        call random_number(rand)
        if (velconst) then
          !Parameterize the closed surface in velocity space:
          !KE = (1/2)sum(dotphi*dotphi)
          !by spherical coords.
          !Use prior on angles since model indep assuming canon KE
          !Gives "uniform area" sample on vels

	        !Energy density not allocated
          phi = y(1:num_inflaton)
          dphi=y(num_inflaton+1:2*num_inflaton)
	        KE = energy_scale**4 - pot(phi)

          call constrain_via_thetas(dphi,sqrt(2.0e0_dp*KE))

          y(num_inflaton+1:2*num_inflaton) = dphi

    	    !Energy density not allocated
    	    rho_not_alloc = energy_scale**4 - pot(phi) -&
            0.5e0_dp*sum(dphi*dphi)


	    	  if (abs(rho_not_alloc)<1e-25) then
            print*, "IC off shell ------------- cycling", ll
            if (ll==maxtry) then
              print*, "    Energy overrun =", rho_not_alloc
              print*, "    E**4 =", energy_scale**4
              print*, "    KE =", 0.5e0_dp*sum(dphi*dphi)
              print*, "    PE =", pot(phi)
              exit
            end if
            cycle
          else
            return
          end if


        else
          param_constr = ceiling(rand*size(y))

          !Set y(constraint) to minim energy value
          y(param_constr) = min_fixed_energy(param_constr)

          phi = y(1:num_inflaton)
          dphi = y(num_inflaton+1:2*num_inflaton)

	    	  !Energy density not allocated
	    	  rho_not_alloc = energy_scale**4 - pot(phi) -&
            0.5e0_dp*sum(dphi*dphi)

	    	  if (abs(rho_not_alloc)<1e-25_dp) then
            print*, "IC off shell ------------- cycling", ll
            if (ll==maxtry) then
              print*, "    Energy overrun =", rho_not_alloc
              print*, "    E**4 =", energy_scale**4
              print*, "    KE =", 0.5e0_dp*sum(dphi*dphi)
              print*, "    PE =", pot(phi)
              exit
            end if
            cycle
          else

	          !Set y(constraint) to fix total energy
            call set_y_by_energy_constraint(y, param_constr, energy_scale)
            phi = y(1:num_inflaton)
            dphi = y(num_inflaton+1:2*num_inflaton)

            !DEBUG
            return

            !if (any(phi<phi0_priors_min) .or. any(dphi<dphi0_priors_min) &
            !  .or. any(phi>phi0_priors_max) &
            !  .or. any(dphi>dphi0_priors_max)) then
            !  if(any(phi<phi0_priors_min)) then
            !    print*,"phi<phi0_priors_min"
            !  endif
            !  if(any(dphi<dphi0_priors_min)) then
            !    print*,"dphi<dphi0_priors_min"
            !  endif
            !  if(any(phi>phi0_priors_max)) then
            !    print*,"phi>phi0_priors_max"
            !  endif
            !  if(any(dphi>dphi0_priors_max)) then
            !    print*,"dphi>dphi0_priors_max"
            !    !DEBUG
            !    !do i=1, size(phi)
            !    !  if (dphi(i)>dphi0_priors_max(i)) then
            !    !    print*, "vals", dphi(i), dphi0_priors_max(i)
            !    !  end if
            !    !end do
            !  endif
            !  !DEBUG
            !  !stop
            !  cycle
            !else
            !  return
            !end if

          end if
        end if

	    end do

      print*, "Couldn't find an IC with energy=energy_scale in the max"
      print*, "number of tries."
      stop

    end subroutine eqen_ic

    !Uniform sampling of fields, unconstrained
    !Also grabs momenta uniformly
    subroutine unconstrained_ic(y, &
        phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max)

      real(dp), dimension(num_inflaton*2), intent(out) :: y
      real(dp), dimension(num_inflaton), intent(in) :: phi0_priors_min, &
        phi0_priors_max, dphi0_priors_min, dphi0_priors_max

      real(dp), dimension(num_inflaton*2) :: delta_prior


	    real(dp) :: rand
	    integer :: i


      delta_prior(1:num_inflaton) = phi0_priors_max-phi0_priors_min
      delta_prior(num_inflaton+1:2*num_inflaton) = &
        dphi0_priors_max-dphi0_priors_min


      !Choose y from unif prior
      do i=1,size(y)
        call random_number(rand)
        y(i) = rand*delta_prior(i)
        if (i > num_inflaton) then
          y(i) = y(i) + dphi0_priors_min(i-num_inflaton)
        else
          y(i) = y(i) + phi0_priors_min(i)
        end if
      end do


    end subroutine unconstrained_ic


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
      real(dp) :: rand, H
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


      else

        !Constraining a field and have to solve the potential
        !Try to do this by hand.
        if (potential_choice==1) then

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

    end subroutine set_y_by_energy_constraint


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


    !Print the cosmo observables to file
    subroutine ic_print_observables(this, outunit)

      class(ic_and_observables) :: this
      integer, intent(in) :: outunit

      if (num_inflaton*2 +12 > 120000) then
        print*, "Don't be silly."
        print*, "Too many fields to print out properly."
        print*, "Fix formatting."
        stop
      end if

      write(outunit, '(120000E18.10)') this%ic(:), this%As, this%ns,&
        this%r, this%nt, this%alpha_s, this%A_iso, this%A_pnad,&
        this%A_ent, this%A_bundle, this%n_iso, this%n_pnad, this%n_ent

    end subroutine ic_print_observables

    !Load the cosmo observables into the ic_and_observables type
    subroutine ic_load_observables(this, phi0,dphi0, As, ns, r, nt,&
      alpha_s, A_iso, A_pnad, A_ent, A_bundle, &
      n_iso, n_pnad, n_ent)

      class(ic_and_observables) :: this
      real(dp), intent(in) :: As, ns, r, nt, alpha_s
      real(dp), intent(in) :: A_iso, A_pnad, A_ent, A_bundle
      real(dp), intent(in) :: n_iso, n_pnad, n_ent
      real(dp), dimension(:), intent(in) :: phi0, dphi0

      integer :: ind

      if (.not. allocated(this%ic)) then
        allocate(this%ic(size(phi0)+size(dphi0)))
      else if (size(phi0) + size(dphi0) /= size(this%ic)) then
        deallocate(this%ic)
        allocate(this%ic(size(phi0) + size(dphi0)))
      end if

      !num_inflaton
      ind = size(phi0)

      this%ic(1:ind) = phi0
      this%ic(ind+1:ind+size(dphi0)) = dphi0
      this%As = As
      this%ns = ns
      this%r = r
      this%nt = nt
      this%alpha_s = alpha_s
      this%A_iso= A_iso
      this%A_pnad=A_pnad
      this%A_ent=A_ent
      this%A_bundle=A_bundle
      this%n_iso= n_iso
      this%n_pnad=n_pnad
      this%n_ent=n_ent

    end subroutine ic_load_observables

    subroutine constrain_via_thetas(y, radius)

      real(dp), dimension(:), intent(out) :: y
      real(dp), dimension(size(y)) :: dphi
      real(dp), allocatable :: theta(:)
      real(dp), intent(in) :: radius
      real(dp) :: rand
      integer :: n, i, j

      n = size(y)

      if (n>1) then
        allocate(theta(n-1))
      else
        print*, "Eqen sampling doesn't work for num_inflaton=1"
        stop
      end if

      !Get the angles with unif priors
      call random_number(rand)
      theta(size(theta)) = rand*2.0e0_dp*pi
      if (size(theta)>1) then
        do i=1,size(theta)-1
          call random_number(rand)
          theta(i) =rand*pi
        end do
      end if

      !Change to Cartesian, i.e., field-space coords
      dphi = radius
      do i=1,n-1
        dphi(i) = dphi(i)*cos(theta(i))

        if (i>1) then
          do j=1,i-1
            dphi(i) = dphi(i)*sin(theta(j))
          end do
        end if

      end do

      do j=1,size(theta)
        dphi(n) = dphi(n)*sin(theta(j))
      end do

      y = dphi

    end subroutine constrain_via_thetas

    !For a sum-separable potential only
    !Set IC for each field at an equal displacement from the local minimum
    !as measured by energy
    subroutine equal_displacement_ic(y,energy_scale)

      real(dp), dimension(num_inflaton*2), intent(out) :: y
      real(dp), intent(in) :: energy_scale
      real(dp) :: m2(num_inflaton), E4

      !N-quadratic, m^2*phi^2
      if (potential_choice==1) then
        !Energy displacement for every field IC
        E4 = energy_scale**4/num_inflaton
        !Masses
        m2 = 10.d0**(vparams(1,:))

        y(1:num_inflaton) = sqrt(2.0e0_dp*E4/m2(:))

        y(num_inflaton+1:2*num_inflaton) = 0e0_dp
      else
        write(*,*) "Error in equal_displacement_ic:"
        write(*,*) "potential_choice=",potential_choice,"not supported."
        stop
      end if

    end subroutine equal_displacement_ic

    !Sets the masses in N-quadratic inflation according to hep-th/0512102
    !Easther-McAllister
    !........
    !Not yet working.
    !subroutine mass_spectrum_nflation(vpnew)

    !Chooses masses from flat/log priors
    subroutine mass_spectrum_prior(vpnew, mass_ratio)

      real(dp), intent(inout), dimension(1,num_inflaton) :: vpnew
      real(dp), optional :: mass_ratio
      real(dp) :: ratio_max, rand
      integer :: i, prior
      integer, parameter :: unif_param=0, log_param=1

      if (size(vparams,1) /=1 .or. size(vparams,2) /= num_inflaton) then
        write(*,*), "Error: vparams of wrong order."
        write(*,*), "size(vparams,2) /= num_inflaton or"
        write(*,*), "size(vparams,1) /= 1."
        stop
      end if

      !m2 set at random (unif on m^2) w/max prior range set by max mass
      !ratio; ratio_max=m2_lightest/m2_heaviest
      if (present(mass_ratio)) then
        ratio_max = mass_ratio
      else
        ratio_max = 10.0e0_dp
      end if

      !prior = log_param
      prior = unif_param

      !vpnew(1,1) = -10.091514981121351e0_dp
      vpnew(1,1) = minval(vparams(1,:))
      do i=2,num_inflaton
        call random_number(rand)
        if (prior == unif_param) then
          rand = rand*(ratio_max-1.0e0_dp) + 1.0e0_dp
          vpnew(1,i) = vpnew(1,1)+2.0e0_dp*log10(rand)
        else if (prior == log_param) then
          rand = rand*2.0e0_dp*log10(ratio_max)
          vpnew(1,i) = vpnew(1,1)+rand
        else
          write(*,*) "Prior not supported. prior=", prior
          stop
        end if
      end do

      print*, "mass ratio",sqrt((10e0**vpnew(1,:))/(10e0**vpnew(1,1)))

    end subroutine mass_spectrum_prior

end module modpk_icsampling
