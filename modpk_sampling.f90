module modpk_sampling
  !Module that implements various sampling techniques for the initial conditions
  !and model parameters.  These routines should generally be called prior to
  !doing any integration.  Monte Carlo methodology.
  use potential, only : norm, pot, getH, getEps_with_t, getEps, getH_with_t
  use modpk_qsf
  use modpkparams, only : dp, slowroll_start, num_inflaton, &
    potential_choice, vparams
  use internals, only : pi
  use modpk_errorhandling, only : raise
  use modpk_io, only : out_opt
  implicit none

  integer :: ic_sampling
  type :: ic_samp_flags
    integer :: reg_samp = 1
    integer :: eqen_samp = 2
    integer :: slowroll_samp = 3
    integer :: fromfile_samp = 4
    integer :: parameter_loop_samp = 5
    integer :: iso_N = 6
    integer :: param_unif_prior = 7
    integer :: qsf_random = 8
    integer :: qsf_parametric = 9
    integer :: fisher_inf = 10
    integer :: single_axion = 11
  end type
  type(ic_samp_flags) :: ic_flags

  logical :: allow_superplanckian

  integer :: inflaton_sampling
  type :: inflaton_samp_flags
    integer :: unif = 1
    integer :: logunif = 2
  end type
  type(inflaton_samp_flags) :: inflaton_flags

  integer :: param_sampling
  type :: param_samp_flags
    integer :: reg_constant = 1
    integer :: unif_prior = 2
    integer :: log_prior = 3
    integer :: num_QSF = 4
    integer :: MarPast_dist = 5
  end type
  type(param_samp_flags) :: param_flags
  real(dp), dimension(:,:), allocatable :: vp_prior_max, vp_prior_min
  real(dp), dimension(:), allocatable :: prior_auxparams_min, prior_auxparams_max
  logical :: use_first_priorval


  logical :: save_iso_N=.false., set_Nisoref_by_Npivot=.false.
  real(dp) :: N_iso_ref
  real(dp), dimension(:), allocatable :: phi_iso_N, dphi_iso_N

  !For choosing an equal-energy prior
  integer, parameter, private :: equal_area_prior=1, unif_prior=2
  integer, parameter, private :: eqen_prior = equal_area_prior


  contains


    !Grab a new IC according to the ic_sampling variable
    subroutine get_ic(phi0, dphi0, &
        priors_min, priors_max, &
         numb_samples,energy_scale)

      use modpk_rng, only : normal
      use modpk_numerics, only : heapsort

      integer, intent(in) :: numb_samples
      real(dp), intent(in), optional :: energy_scale
      real(dp), dimension(num_inflaton), intent(out) :: phi0, dphi0
      real(dp), dimension(2,num_inflaton), intent(in) :: priors_min, &
        priors_max

      real(dp), dimension(2*num_inflaton) :: y_background
      real(dp), dimension(num_inflaton) :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max

      real(dp), dimension(:,:), allocatable :: sample
      real(dp) :: H
      integer :: i,j, kk, ii

      logical :: with_velocity, with_eqen
      real(dp) :: beta, rand
      real(dp) :: finv(num_inflaton)
      real(dp) :: lambda(num_inflaton)

      !DEBUG
      !Hard-coded parameter priors for potential_choice==11
      real(dp), dimension(size(vparams,1),size(vparams,2)) :: &
        vparams_priors_min, vparams_priors_max

      real(dp), dimension(:,:), allocatable :: knot_temp

      real(dp) :: param0, phi_light0

      real(dp) :: ic_radius, p_exp


      phi0_priors_max=priors_max(1,:)
      phi0_priors_min=priors_min(1,:)
      dphi0_priors_max=priors_max(2,:)
      dphi0_priors_min=priors_min(2,:)

      !-----------------------------------------
      if (ic_sampling == ic_flags%eqen_samp) then

        !Tell solver starting out of SR, even if in SR
        !in case it's only transient, e.g., starting with dphi=0
        !but not on "attractor"
        if (slowroll_start) slowroll_start = .false.


        call eqen_ic(y_background, energy_scale, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max, &
          velconst=.true.)

      !-----------------------------------------
      else if (ic_sampling == ic_flags%iso_N) then

        !Start in SR with ICs sampled N_iso_ref before inflation ends
        !Only works for N-quadratic

        if (potential_choice /= 1 .and. potential_choice /= 16) then
          print*, "MODECODE: potential_choice=", potential_choice
          call raise%fatal_code(&
            "Can't implicitly define the iso_N surface for this potential choice.",&
            __FILE__, __LINE__)
        end if

        !Can't also record iso-N, since using N_iso_ref as override
        save_iso_N = .false.


        !Setting N_iso_ref<0 defaults to N_iso_ref = N_pivot, which may vary
        if (set_Nisoref_by_Npivot) N_iso_ref = N_pivot

        y_background = 0e0_dp

        if (potential_choice ==16) then
          p_exp = vparams(2,1)
        else
          p_exp = 2.0e0_dp
        end if

        ic_radius = sqrt(2.0e0_dp * N_iso_ref * p_exp)

        call sample_nsphere(y_background(1:num_inflaton),ic_radius)


      !-----------------------------------------
      else if (ic_sampling == ic_flags%slowroll_samp) then

        !Force solver to start in SR
        slowroll_start = .true.

        !Will override what vels it chooses later
        call unconstrained_ic(y_background, &
          phi0_priors_min, phi0_priors_max, &
          dphi0_priors_min, dphi0_priors_max)

      !-----------------------------------------
      else if (ic_sampling==ic_flags%parameter_loop_samp) then

        if (potential_choice>2 .and. potential_choice/=11 .and. &
          potential_choice/=12 .and. potential_choice /=16) then

          print*, "MODECODE: potential_choice=", potential_choice
          call raise%fatal_code(&
            "The sampling technique parameter_loop_samp &
            doesn't work for this potential choice.",&
            __FILE__,__LINE__)

        end if

        !Get new vparams
        if (potential_choice==1) then
          !N-quadratic -- m^2 phi^2

          call n_quadratic_mass_looping(vparams)

          if (numb_samples > 1) then
            call random_number(rand)
            beta = rand*(0.3e0_dp) + 0.35e0_dp
          else
            beta = 0.5e0_dp
          end if

          call mass_spectrum_nflation(vparams,beta)

          !DEBUG
          print*, "vparams from mass_spectrum_nflation"
          print*, vparams(1,:)
          stop

        else if (potential_choice==2) then
          !N-flation (axions-cosine)
          call n_flation_looping(vparams, energy_scale)
          !call n_flation_random_prior(vparams)

        !DEBUG
        !else if (potential_choice==11) then
        else if (potential_choice==11 .or. potential_choice==12) then
          !N-quadratic w/intxn for lightest field
          !m^2 phi^2 + phi_lightest^2 phi_i^2
          call n_quadratic_mass_intxn_looping(vparams)

        end if


        !Set IC with either (V=0 and KE) or (KE=0 and V)
        with_velocity = .false.
        with_eqen = .true.
        if (num_inflaton==1) then

          y_background(1:num_inflaton) = 20.0e0_dp
          y_background(num_inflaton+1:2*num_inflaton) = 0e0_dp

        else if (with_eqen) then
          !Set IC with equal energy over pre-defined field range
          if(out_opt%modpkoutput) print*, "Setting IC with equal energy over range"

          call eqen_ic(y_background, energy_scale, &
            phi0_priors_min, phi0_priors_max, &
            dphi0_priors_min, dphi0_priors_max, &
            velconst=.true.)

        else if (with_velocity) then
          !Set IC with zero potential
          if(out_opt%modpkoutput) print*, "Setting IC with zero potential"
          call zero_potential_ic(y_background, energy_scale)
        else
          !Set IC with zero velocity
          if(out_opt%modpkoutput) print*, "Setting IC with zero velocity"
          call equal_displacement_ic(y_background, energy_scale)
        end if


      !-----------------------------------------
      else if (ic_sampling==ic_flags%param_unif_prior) then

        if (potential_choice /= 11) then
          print*, "MODECODE: potential_choice=", potential_choice
          call raise%fatal_code(&
            "The sampling technique param_unif_prior &
            doesn't work for this potential choice.",&
            __FILE__,__LINE__)
        end if

        if(out_opt%modpkoutput) print*, "BAD MODECODE: Parameter priors hard-coded in..."

        !Some inspiration for these priors from 1112.0326

        !Masses: m^2 = 10**vparams(1,:)
        vparams_priors_min(1,:) = -13.5e0_dp
        vparams_priors_max(1,:) = -8.0e0_dp

        !Intxn:  lambda^4 = 10**vparams(2,:)
        vparams_priors_min(2,:) = -20.0e0_dp
        vparams_priors_max(2,:) = -10.0e0_dp

        do i=1,size(vparams,1); do j=1,size(vparams,2)
          call random_number(rand)
          vparams(i,j) = (vparams_priors_max(i,j)-vparams_priors_min(i,j))*rand + &
            vparams_priors_min(i,j)
        end do; end do



        !Get IC
        if (num_inflaton==1) then

          y_background(1:num_inflaton) = 25.0e0_dp
          y_background(num_inflaton+1:2*num_inflaton) = 0e0_dp

        else
          !Set IC with equal energy over pre-defined field range
          if(out_opt%modpkoutput) print*, "Setting IC with equal energy over range"

          call eqen_ic(y_background, energy_scale, &
            phi0_priors_min, phi0_priors_max, &
            dphi0_priors_min, dphi0_priors_max, &
            velconst=.true.)
        end if

      !-----------------------------------------
      else if (ic_sampling==ic_flags%qsf_random) then

        turning_choice = 4

        !Set IC from parameters file
        !Default heavy fields to 0
        y_background = 0e0_dp
        y_background(1) = phi0(1)

        ![ LP: ] Indices for knots:
        !knot_positions(HEAVYFIELD, KNOT#, (LIGHT_POS, HEAVY_POS))
        if (allocated(knot_positions)) deallocate(knot_positions)
        if (num_inflaton>1 .and. maxval(number_knots_qsfrandom)>0) then
          allocate(knot_positions(num_inflaton-1, maxval(number_knots_qsfrandom), 2))
        else
          print*, "MODECODE: Num_inflaton=", num_inflaton
          print*, "MODECODE: number_knots_qsfrandom=", number_knots_qsfrandom

          call raise%fatal_cosmo(&
            "Incorrect specifications for knot_positions.",&
            __FILE__, __LINE__)

        end if
        knot_positions = 0e0_dp

        do i=1,size(knot_positions,1) !# heavy fields
          do j=1, size(knot_positions,2) !# knots

            if (j > number_knots_qsfrandom(i)) then
              knot_positions(i,j,:) = 0e0_dp
            else
              !Knot light field position
              call random_number(rand)
              !knot_positions(i,j,1) = phi0(1)*rand

              !DEBUG
              if (custom_knot_range) then
                if (knot_range_min(i) .ge. knot_range_max(i)) then
                  print*, "MODECODE: knot_range_min =", knot_range_min
                  print*, "MODECODE: knot_range_max =", knot_range_max

                  call raise%fatal_code(&
                    "The knot_range_min > knot_range_max.",&
                    __FILE__, __LINE__)
                end if

                knot_positions(i,j,1) = (knot_range_max(i)-knot_range_min(i))*rand &
                  + knot_range_min(i)
              else
                knot_positions(i,j,1) = phi0(1)*rand
              end if

              !Knot heavy field position
              knot_positions(i,j,2) = normal(0.0e0_dp,stand_dev_qsfrandom(i))

            end if

          end do
        end do

        !Sort the knots by light field position
        do i=1, size(knot_positions,1)
          if (allocated(knot_temp)) deallocate(knot_temp)
          allocate(knot_temp(number_knots_qsfrandom(i),2))
          knot_temp = knot_positions(i,1:number_knots_qsfrandom(i),:)
          call heapsort(knot_temp)
          knot_positions(i,1:number_knots_qsfrandom(i),:) = knot_temp

        end do

      !-----------------------------------------
      else if (ic_sampling == ic_flags%qsf_parametric) then
        !"Numerical" QSF trajectory, needs to set initial position based off the
        !parametric curve.

        param0 = vparams(2,1)
        phi_light0 = vparams(3,1)

        !Integrate through traj and set-up interpolation if first time
        call qsf_runref%initialize_traj(phi_light0,param0)

        !Find the initial param that coincides with setting
        !IC in minimum of valley (dist=0) with phi_light=phi_light0
        call qsf_runref%get_param(phi_light=phi_light0)

        !Set initial condition from this parameter in the valley
        y_background(1:num_inflaton) = turning_function_parametric(qsf_runref%param)
        y_background(num_inflaton+1:2*num_inflaton) = 0e0_dp

      else if (ic_sampling == ic_flags%single_axion) then
        !Kinetic-dominated initial condition for single-field axion
        !Sets initial field value with uniform prior over range
        !-pi/f < phi < pi/f

        if (.not.(potential_choice == 2 .or. potential_choice==18)) then
          print*, "MODECODE: potential_choice=", potential_choice
          call raise%fatal_code(&
            "The sampling technique single_axion &
            doesn't work for this potential choice.",&
            __FILE__,__LINE__)
        end if

        if (num_inflaton>1) then
          call raise%fatal_code(&
            "This IC sampler only works for &
            one field.",&
            __FILE__,__LINE__)
        end if

        !Axion params
        lambda = 10.e0_dp**vparams(1,:)
        finv = 1.e0_dp/(10.e0_dp**vparams(2,:))


        !Set field IC with uniform random prior

        !Something weird with the clock in init_random_seed_serial?
        call random_number(rand); call random_number(rand)
        y_background(1:num_inflaton) = rand*&
            (2.0e0_dp*pi/finv)-(pi/finv)

        !Set velocity much higher than the "axion scale"=lambda
        !y_background(num_inflaton+1:2*num_inflaton) = &
        !    1.0e6_dp*2.0e0_dp*lambda**2

        !Set KE at SSB scale f^4
        y_background(num_inflaton+1:2*num_inflaton) = &
            !finv**(-2)
            1.0e0_dp


        !y_background(num_inflaton+1:2*num_inflaton) = sqrt(2.0e0_dp)
        !y_background(num_inflaton+1:2*num_inflaton) = 0e0_dp


        y_background(1:num_inflaton) = (0.95*pi/2.0e0_dp)/finv
        y_background(num_inflaton+1:2*num_inflaton) = 0e0_dp



      !-----------------------------------------
      else
        call raise%fatal_code(&
            "The sampling technique hasn't been implemented.",&
            __FILE__, __LINE__)
      end if

      !Load initial vals from sample
      phi0 = y_background(1:num_inflaton)
      dphidt_init0 = y_background(num_inflaton+1:2*num_inflaton)

      !Convert dphidt-->dphidN
      H=getH_with_t(phi0, dphidt_init0)
      dphi0 = (1.0e0/H)*y_background(num_inflaton+1:2*num_inflaton)

      !Check if any fields are super-Planckian, if that's not allowed
      if (.not. allow_superplanckian .and. any(abs(phi0)>1.0e0_dp)) then
        call raise%fatal_cosmo(&
            "A field value is super-Planckian, which you have not allowed.",&
            __FILE__, __LINE__)
      end if

    end subroutine get_ic



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
      integer, parameter :: maxtry=1000, maxtry_unif=10*maxtry

      !For unif_prior
      real(dp), dimension(num_inflaton) :: phi0_min, &
        phi0_max, dphi0_min, dphi0_max

      if(out_opt%modpkoutput) &
        print*, "MODECODE: IC with equal energy and eqen_prior=", eqen_prior

      if (eqen_prior==equal_area_prior) then
        !Use the pseudo--"equal-area" prior of Easther-Price (1304.4244)

	      do ll=1,maxtry

          y = 0e0_dp

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

            call sample_nsphere(dphi,sqrt(2.0e0_dp*KE))

            y(num_inflaton+1:2*num_inflaton) = dphi

    	      !Energy density not allocated
    	      rho_not_alloc = energy_scale**4 - pot(phi) -&
              0.5e0_dp*sum(dphi*dphi)

	      	  if (abs(rho_not_alloc)<1e-25 .or. isnan(rho_not_alloc)) then
              if(mod(maxtry,ll)==20) then
                if (mod(ll,200)==0 .and. out_opt%modpkoutput)&
                  print*, "IC off shell ------------- cycling", ll
              end if
              if (ll==maxtry) then
                if (out_opt%modpkoutput) then
                  print*, "MODECODE: Energy overrun =", rho_not_alloc
                  print*, "MODECODE: E**4 =", energy_scale**4
                  print*, "MODECODE: KE =", 0.5e0_dp*sum(dphi*dphi)
                  print*, "MODECODE: PE =", pot(phi)
                  call raise%warning(&
                    'Energy overrun when setting initial conditions &
                    with equal energy prior',&
                    __FILE__, __LINE__)
                end if
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

	      	  if (rho_not_alloc<1e-25_dp) then
              if(out_opt%modpkoutput) print*, "IC off shell ------------- cycling", ll
              if (ll==maxtry) then
                if (out_opt%modpkoutput) then
                  print*, "MODECODE: Energy overrun =", rho_not_alloc
                  print*, "MODECODE: E**4 =", energy_scale**4
                  print*, "MODECODE: KE =", 0.5e0_dp*sum(dphi*dphi)
                  print*, "MODECODE: PE =", pot(phi)
                  call raise%warning(&
                    "Energy overrun when setting initial conditions with equal &
                    energy prior",&
                    __FILE__, __LINE__)
                end if
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


      else if (eqen_prior==unif_prior) then
        !Define the constraint surface with the first velocity
        !v_1^2 = 2E_0^4 - 2V - v_2^2 - v_3^2 - ...
        !Then does a uniform sampling over the remaining fields and velocities.

        !Useful because we can then re-weight the ICs in post-processing
        !\Integral d(IC_constrained) --> \Integral P(IC_constraint) d(IC_constrained)

        !Uses the prior ranges in parameters_multimodecode.txt

	      do ll=1,maxtry_unif

          y = 0e0_dp

          call unconstrained_ic(y, &
            phi0_priors_min, phi0_priors_max, &
            dphi0_priors_min, dphi0_priors_max)

          phi = y(1:num_inflaton)
          dphi = y(num_inflaton+1:2*num_inflaton)

          !Constrain first velocity
          dphi(1) = 0e0_dp

	    	  !Energy density not allocated
	    	  rho_not_alloc = energy_scale**4 - pot(phi) -&
            0.5e0_dp*sum(dphi*dphi)

	      	if (rho_not_alloc<1e-25_dp) then
            if(mod(ll,200)==0) then
              if (mod(ll,200)==0 .and. out_opt%modpkoutput)&
                print*, "IC off shell ------------- cycling", ll
            end if
            if (ll==maxtry) then
              if (out_opt%modpkoutput) then
                print*, "MODECODE: Energy overrun =", rho_not_alloc
                print*, "MODECODE: E**4 =", energy_scale**4
                print*, "MODECODE: KE =", 0.5e0_dp*sum(dphi*dphi)
                print*, "MODECODE: PE =", pot(phi)
                call raise%warning(&
                  'Energy overrun when setting initial conditions &
                  with equal energy prior',&
                  __FILE__, __LINE__)
              end if
              exit
            end if
            cycle
          else
	          !Set y(constraint) to fix total energy
            call set_y_by_energy_constraint(y, num_inflaton+1, energy_scale)
            phi = y(1:num_inflaton)
            dphi = y(num_inflaton+1:2*num_inflaton)

            return
          end if

        end do



      else
        print*, "MODECODE: eqen_prior = ", eqen_prior

        call raise%fatal_code(&
            "The prior for the equal energy sampler is not recognized.",&
            __FILE__, __LINE__)

      end if

      call raise%fatal_code(&
          "Couldn't find an IC with energy=energy_scale &
          in the max number of tries.",&
          __FILE__, __LINE__)

    end subroutine eqen_ic

    !Uniform sampling of fields, unconstrained
    !Also grabs momenta uniformly
    subroutine unconstrained_ic(y, &
        phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max)

      real(dp), dimension(num_inflaton*2), intent(inout) :: y
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

          m2_V = 10.e0_dp**(vparams(1,:))

          y(constraint) = rand_sign()*sqrt((2e0_dp/m2_V(constraint))*&
            energy_remaining)


        else
          !DEBUG
          call raise%fatal_code(&
          "Specify the degree of the polynomial that makes &
          the energy constraint equation or put the &
          algebraic constraint in by-hand.", &
          __FILE__, __LINE__)
        end if

      end if

    end subroutine set_y_by_energy_constraint

    subroutine sample_nsphere(y, radius)
      use modpk_rng, only : normal_array

      real(dp), dimension(:), intent(out) :: y
      real(dp), intent(in) :: radius

      !Sample Gaussian with zero mean and fixed std
      y = normal_array(size(y), 0e0_dp, 1e0_dp)
      y = (radius/norm(y))*y

    end subroutine sample_nsphere


    !For a sum-separable potential only
    !Set IC for each field at an equal displacement from the local minimum
    !as measured by energy
    subroutine equal_displacement_ic(y,energy_scale)

      real(dp), dimension(num_inflaton*2), intent(out) :: y
      real(dp), intent(in) :: energy_scale
      real(dp) :: E4
      real(dp), dimension(num_inflaton) :: m2, l4, f

      !Energy displacement for every field IC
      E4 = energy_scale**4/num_inflaton

      !N-quadratic, m^2*phi^2
      if (potential_choice==1) then
        !Masses
        m2 = 10.e0_dp**(vparams(1,:))

        y(1:num_inflaton) = sqrt(2.0e0_dp*E4/m2(:))

      else if (potential_choice==2) then
        !N-flation axions

        l4 = 10.e0_dp**vparams(1,:)
        f = 10.e0_dp**vparams(2,:)

        y(1:num_inflaton) = acos(E4/l4 - 1.0e0_dp)*f
        if (any(isnan(y))) then
          call raise%fatal_code(&
            "The variable y has a NaN in equal_displacement_ic.",&
            __FILE__, __LINE__)
        end if

      !N-quadratic with intxn w/lightest field, m^2*phi^2+phi_light^2 phi_i^2
      !Just displace same as m^2phi^2
      else if (potential_choice==11 .or. potential_choice==12) then
        !Masses
        m2 = 10.e0_dp**(vparams(1,:))

        y(1:num_inflaton) = sqrt(2.0e0_dp*E4/m2(:))
      else
        write(*,*) "MODECODE: potential_choice=",potential_choice

        call raise%fatal_code(&
            "This potential choice isn't supported.",&
            __FILE__, __LINE__)
      end if

      y(num_inflaton+1:2*num_inflaton) = 0e0_dp

    end subroutine equal_displacement_ic

    !For a sum-separable potential only
    !Set IC for each field at an equal displacement from the local minimum
    !as measured by energy
    subroutine zero_potential_ic(y,energy_scale)

      real(dp), dimension(num_inflaton*2), intent(inout) :: y
      real(dp), intent(in) :: energy_scale
      real(dp) :: KE
      !real(dp), dimension(num_inflaton) :: phi, dphi

      y(1:num_inflaton) = 0.1e-10

      KE = energy_scale**4 - pot(y(1:num_inflaton))

      call sample_nsphere(y(num_inflaton+1:2*num_inflaton),&
        sqrt(2.0e0_dp*KE))

    end subroutine zero_potential_ic


    !Chooses masses from flat/log priors
    subroutine n_quadratic_mass_looping(vpnew, mass_ratio)

      real(dp), intent(inout), dimension(1,num_inflaton) :: vpnew
      real(dp), optional :: mass_ratio
      real(dp) :: ratio_max, rand
      integer :: i, prior
      integer, parameter :: unif_param=0, log_param=1

      !if (size(vparams,1) /=1 .or. size(vparams,2) /= num_inflaton) then
      !  write(*,*), "Error: vparams of wrong order."
      !  write(*,*), "size(vparams,2) /= num_inflaton or"
      !  write(*,*), "size(vparams,1) /= 1."
      !  stop
      !end if

      !m2 set at random (unif on m^2) w/max prior range set by max mass
      !ratio; ratio_max=m2_lightest/m2_heaviest
      if (present(mass_ratio)) then
        ratio_max = mass_ratio
      else
        ratio_max = 9.0e0_dp
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

          write(*,*) "MODECODE: prior=", prior
          call raise%fatal_code(&
            "This prior is not supported.",&
            __FILE__, __LINE__)

        end if
      end do

      print*, "mass ratio",sqrt((10e0**vpnew(1,:))/(10e0**vpnew(1,1)))

    end subroutine n_quadratic_mass_looping

    !Chooses masses from flat/log priors
    subroutine n_quadratic_mass_intxn_looping(vpnew)

      real(dp), intent(inout), dimension(2,num_inflaton) :: vpnew
      real(dp) :: rand, intxn_min, intxn_max
      real(dp) :: log_intxn_min, log_intxn_max
      integer :: i, prior
      integer, parameter :: unif_param=0, log_param=1

      !Set masses as in N-quadratic
      call n_quadratic_mass_looping(vpnew(1,:))

      !DEBUG
      if (potential_choice==12) then
        print*, "-------------------------------------------"
        print*, "Not varying interaction term in mass matrix"
        print*, "for potential_choice==", potential_choice
        print*, "-------------------------------------------"
        return
      end if

      !Set intxn term over interval

      prior = log_param

      do i=1, size(vpnew,2)
        call random_number(rand)
        if (prior==unif_param) then
          intxn_max = 1e-12_dp
          intxn_min = 0e0_dp
          vpnew(2,i) = log10(rand*(intxn_max-intxn_min)+intxn_min)

        else if (prior==log_param) then
          log_intxn_max = 2e0_dp + minval(vpnew(1,:))
          log_intxn_min = -4e0_dp +minval(vpnew(1,:))

          vpnew(2,i) = rand*(log_intxn_max-log_intxn_min)+log_intxn_min

        else

          write(*,*) "MODECODE: prior=", prior
          call raise%fatal_code(&
              "The prior is not recognized.",&
              __FILE__, __LINE__)

        end if

      end do

    end subroutine n_quadratic_mass_intxn_looping


    !Place priors on the terms in expansion of the cos term into
    ! 0.5*m^2phi^2 - (1/24)*alpha^2*phi^4
    subroutine n_flation_looping(vpnew, energy_scale, mass_ratio)

      real(dp), intent(inout), dimension(2,num_inflaton) :: vpnew
      real(dp), intent(in) :: energy_scale
      real(dp), optional :: mass_ratio
      real(dp), dimension(num_inflaton) :: logmasses2, alpha4, masses2, alpha_max4,&
        lambda,f
      real(dp) :: ratio_max, rand
      integer :: i, prior
      integer, parameter :: unif_param=0, log_param=1

      if (size(vparams,1) /=2 .or. size(vparams,2) /= num_inflaton) then
        call raise%fatal_code(&
            "The vparams are of the wrong order. &
            size(vparams,2) /= num_inflaton or &
            size(vparams,1) /= 2.", &
            __FILE__, __LINE__)
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

      !Get masses**2
      logmasses2 = 0e0_dp
      !logmasses2(1) = minval(vparams(1,:))
      logmasses2(1) = -10.091514981121351e0_dp
      do i=2,num_inflaton
        call random_number(rand)
        if (prior == unif_param) then
          rand = rand*(ratio_max-1.0e0_dp) + 1.0e0_dp
          logmasses2(i) = logmasses2(1)+2.0e0_dp*log10(rand)
        else if (prior == log_param) then
          rand = rand*2.0e0_dp*log10(ratio_max)
          logmasses2(i) = logmasses2(1)+rand
        else

          write(*,*) "MODECODE: prior=", prior
          call raise%fatal_code(&
            "This prior is not supported.",&
            __FILE__, __LINE__)
        end if
      end do
      masses2 = 10.0e0_dp**logmasses2

      !Get quartic expansion term
      do i=1,num_inflaton
        call random_number(rand)
        alpha4(i) = (masses2(i)*rand)**2
      end do

      !Convert to vparams for lambda and f
      lambda(:) = (masses2(:)**0.5)/(alpha4(:)**0.25)
      f(:) = (masses2(:)**0.5)/(alpha4(:)**0.5)

      vpnew(1,:) = log10(lambda(:))
      vpnew(2,:) = log10(f(:))

      !Set max val of alpha st quadratic energy = quartic energy initially
      !alpha_max4= (6.0e0_dp*num_inflaton*masses2**2)/energy_scale**4

      !alpha4=0e0_dp

      !do i=1,num_inflaton
      !  call random_number(rand)
      !  if (prior == unif_param) then
      !    alpha4(i) = rand*alpha_max4(i)
      !  else if (prior == log_param) then
      !    rand = rand*log10(alpha_max4(i))-20.0e0_dp
      !    alpha4(i) = 10.0e0_dp**(rand)
      !  else
      !    write(*,*) "Prior not supported. prior=", prior
      !    stop
      !  end if
      !end do


      !DEBUG
      print*, "mass ratio",sqrt((10e0**logmasses2(:))/(10e0**logmasses2(1)))
      print*, "alpha4",alpha4
      print*, "lambda", lambda
      print*, "f", f
      print*, "masses2", masses2


    end subroutine n_flation_looping


    subroutine n_flation_random_prior(vpnew)

      real(dp), intent(inout), dimension(2,num_inflaton) :: vpnew
      real(dp) :: rand
      real(dp), dimension(num_inflaton) :: ranges_l, ranges_f
      integer :: i

      ranges_l = (/-5.0e0_dp, -1.96e0_dp  /)
      ranges_f = (/ 0.5e0_dp, 1.3e0_dp  /)

      do i=1,num_inflaton
        call random_number(rand)
        vpnew(1,i)=rand*(ranges_l(2)-ranges_l(1))+ ranges_l(1)
        call random_number(rand)
        vpnew(2,i)=rand*(ranges_f(2)-ranges_f(1))+ ranges_f(1)
      end do


    end subroutine n_flation_random_prior

    !Sets the masses in N-quadratic inflation according to hep-th/0512102
    !Easther-McAllister
    !........

    SUBROUTINE mass_spectrum_nflation(msqd, beta)
      implicit none

      !real(dp), Parameter ::  PI=3.1415926535897932385

      integer argnum

      Integer loop1
      real(dp) h1,start,finish,deltaM,delta

      real(dp) ::  mass4

      real(dp), intent(in) :: beta
      real(dp), dimension(:,:), intent(out) :: msqd
      integer  :: nval

      nval = size(msqd,2)

      mass4 = 1e0_dp

      if(beta.ne.0e0_dp) then ! beta = 0 all masses identical

        ! set mass terms
        msqd(1,1) = (1e0_dp-sqrt(beta))**2

        deltaM = ((1e0_dp+sqrt(beta))**2 -(1e0_dp-sqrt(beta))**2)/real(nval-1)

        delta = 1e0_dp/real(nval-1)

        do loop1 = 2,nval-1
        start = msqd(1,loop1-1)
        finish = bisect(start,delta,start,start+deltaM)
        msqd(1,loop1) = finish
        enddo

        msqd(1,nval) = (1e0_dp+sqrt(beta))**2


        do loop1 =1,nval
        msqd(1,loop1) = mass4**2 * msqd(1,loop1)
        enddo


      else
        do loop1=1,nval
        msqd(1,loop1) = mass4**2
        enddo
      endif


      !Setting mass scale by hand.
      !See Eq 4.6 in arXiv: hep-th/0512102
      !msqd = msqd*(1e-5_dp**2)

      !For Planck-normalization on SR, single-field attractor solution
      msqd = msqd*(4.30e-11_dp)



      !Convert to vparams
      msqd = log10(msqd)


    contains

      pure real(dp) function partial(x,y)
        implicit none


        real(dp), intent(in) :: x, y
        real(dp) s

        call qtrap(x,y,s)
        partial = s

      end function

      !--------------------------------------------------------------------
      !   marchenko-pastur distribution, with modification to exclude delta
      !   -fn at origin. (Multiplied by beta)

      pure function pdf(x)
        implicit none

        real(dp) :: pdf


        real(dp) sb
        real(dp), intent(in) :: x

        sb = sqrt(beta)


        if(x.lt.(1e0_dp-sb)**2) then
          pdf = 0e0_dp
          return
        endif

        if(x.gt.(1e0_dp+sb)**2) then
          pdf = 0e0_dp
          return
        endif
        pdf = sqrt((x-(1e0_dp-sb)**2)*((1e0_dp+sb)**2-x))/(2e0_dp *pi* x)/beta
        return
      end function

      !---------------------------------------------------------
      !
      !     a quick and dirty bisection rountine.
      !

      pure real(dp) function bisect(base,target,lostart,histart)
        implicit none


        real(dp) hi,lo
        real(dp), intent(in) :: target, base,lostart,histart

        real(dp) midpt,midval,hival,loval

        lo = lostart
        hi = histart

        hival = partial(base,hi)-target

        loval = partial(base,lo)-target


        do while(hival.lt.0e0_dp)
        hi = hi + (hi-lo)
        hival = partial(base,hi) -target
        enddo


        do while( (hi-lo).gt.10d-10)
        midpt = lo+(hi-lo)/2e0_dp
        midval = partial(base,midpt)-target
        if(loval*midval.ge.0e0_dp) then ! loval and midval have same sign
          loval = midval
          lo =midpt
        else
          hival=midval
          hi = midpt
        endif
        enddo

        bisect = lo + (hi-lo)/2e0_dp
        return
      end function

      pure SUBROUTINE qtrap(a,b,s)
        implicit none
        INTEGER JMAX
        real(dp), intent(in) :: a,b
        real(dp), intent(out) :: s
        real(dp) :: EPS
        PARAMETER (EPS=1.d-8, JMAX=200)
        INTEGER j
        real(dp) olds
        olds=-1.e30_dp

        do j=1,JMAX
        call trapzd(a,b,s,j)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        if (s.eq.0.e0_dp.and.olds.eq.0.e0_dp.and.j.gt.6) return
        olds=s
        end do

      END SUBROUTINE

      pure SUBROUTINE trapzd(a,b,s,n)
        implicit none
        INTEGER, intent(in) :: n
        real(dp), intent(in) :: a,b
        real(dp), intent(out) :: s
        INTEGER it,j
        real(dp) del,sum,tnm,x


        if (n.eq.1) then
          s=0.5e0_dp*(b-a)*(pdf(a)+pdf(b))
        else

          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5e0_dp*del
          sum=0.e0_dp

          do j=1,it
          sum=sum+pdf(x)
          x=x+del
          end do

          s=0.5e0_dp*(s+(b-a)*sum/tnm)
        endif
        return
      END SUBROUTINE


    end SUBROUTINE mass_spectrum_nflation

    subroutine get_new_N_pivot(Npiv, prior_min, prior_max)

      real(dp), intent(inout) :: Npiv
      real(dp), intent(in) :: prior_max, prior_min

      real(dp) :: rand

      call random_number(rand)
      Npiv = (prior_max - prior_min)*rand + prior_min

    end subroutine get_new_N_pivot

    !Subroutine to setup the vparams based off the chosen sampling technique.
    subroutine get_vparams()
      use modpk_rng, only : shuffle

      real(dp) :: rand
      integer :: ii, jj
      real(dp), dimension(:,:), allocatable :: rectang_RM, square_RM
      real(dp), dimension(:), allocatable :: eigvals
      integer :: mp_M, mp_N, flag_eigvals
      real(dp) :: rand_range, sigma, beta, beta_min, beta_max, sigma_min, sigma_max
      real(dp), dimension(:), allocatable :: work_eigvals

      !Don't need to marginalize over vparams
      if (param_sampling == param_flags%num_QSF .or. &
          param_sampling == param_flags%reg_constant) then
        return
      end if

      !Set up prior; if use_first_priorval, then only have to specify the first
      !column in each prior array.
      if (use_first_priorval) then
        do ii=1,size(vparams,2)
          vp_prior_min(:,ii) = vp_prior_min(:,1)
          vp_prior_max(:,ii) = vp_prior_max(:,1)
        end do
      end if

      !Uniform prior
      if (param_sampling == param_flags%unif_prior) then

        do ii=1,size(vparams,1); do jj=1,size(vparams,2)
          call random_number(rand)
          vparams(ii,jj) = (vp_prior_max(ii,jj)-vp_prior_min(ii,jj))*rand &
            +vp_prior_min(ii,jj)
        end do; end do

      !"Log" prior => Uniform prior on x for vparams= 10^x
      else if (param_sampling == param_flags%log_prior) then

        do ii=1,size(vparams,1); do jj=1,size(vparams,2)
          call random_number(rand)
          vparams(ii,jj) = 10.0e0_dp**&
            ((vp_prior_max(ii,jj)-vp_prior_min(ii,jj))*rand+vp_prior_min(ii,jj))
        end do; end do

#ifdef MKL_EXIST
      !Sample masses from the Marcenko-Pastur distribution from random matrix arguments
      else if (param_sampling == param_flags%MarPast_dist) then
        !other params: 1 = beta, 2 = sigma

        !Build a big random matrix

        !mp_M ~ # rows
        !mp_N ~ # cols
        !beta = M/N

        !Get beta from uniform prior
        beta_min = prior_auxparams_min(1)
        beta_max = prior_auxparams_max(1)
        call random_number(rand)
        beta = rand*(beta_max-beta_min)+beta_min


        !Save this for output
        auxparams(1) = beta


        mp_M = min(1000,num_inflaton)
        mp_N = mp_M/beta


        if (allocated(rectang_RM)) deallocate(rectang_RM)
        if (allocated(square_RM)) deallocate(square_RM)
        if (allocated(eigvals)) deallocate(eigvals)

        allocate(rectang_RM(mp_M,mp_N))
        allocate(square_RM(mp_M,mp_M))
        allocate(eigvals(mp_M))

        !Get sigma from uniform prior
        sigma_min = prior_auxparams_min(2)
        sigma_max = prior_auxparams_max(2)
        call random_number(rand)
        sigma = rand*(sigma_max - sigma_min) + sigma_min

        !Save this for output
        auxparams(2) = sigma

        !sigma is stand dev of dist for elements of RM
        rand_range=sqrt(3.0e0_dp)*sigma
        do jj=1,mp_N; do ii=1,mp_M
          call random_number(rand)
          rand = rand*2.0e0_dp*abs(rand_range)-rand_range
          rectang_RM(ii,jj) = rand

        end do; end do


        !Matrix multiplication rectang_RM.rectang_RM^T
        call dgemm( 'n', 't', &
          mp_M, mp_M, mp_N, &
          1.0e0_dp, &
          rectang_RM, mp_M, &
          rectang_RM, mp_M, &
          0e0_dp, &
          square_RM, mp_M)

        square_RM = square_RM/mp_N

        !Find the eigenvalues
        if (allocated(work_eigvals)) deallocate(work_eigvals)
        allocate(work_eigvals(max(1,(3*mp_M-1)*10)))
        call dsyev('N','U', &
          mp_M, square_RM, mp_M, &
          eigvals, &
          work_eigvals, size(work_eigvals), &
          flag_eigvals)

        !Mix 'em up
        call shuffle(eigvals)

        !Set masses to eigenvalues
        if (potential_choice /= 1) then
          print*, "MODECODE: potential_choice=", potential_choice

          call raise%fatal_code(&
              "This potential is not yet supported for the MP distribution.",&
              __FILE__, __LINE__)
        else

          if (size(eigvals) < size(vparams,2)) then

            print*, "MODECODE: size(eigvals) < num_inflaton"

            call raise%fatal_code(&
                "Need to increase the size of the matrix (mp_M) that determines &
                the eigenvalues for the MP distribution.",&
                __FILE__, __LINE__)

          else
            vparams(1,1:size(vparams,2)) = log10(eigvals(1:size(vparams,2)))
          end if
        end if

        deallocate(rectang_RM)
        deallocate(square_RM)
        deallocate(eigvals)

#endif
      else

        print*, "MODECODE: param_sampling=", param_sampling

        call raise%fatal_code(&
            "This choice of param_sampling is not supported.",&
            __FILE__, __LINE__)

      end if

    end subroutine get_vparams


end module modpk_sampling
