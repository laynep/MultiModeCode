module modpk_icsampling
  !Module that implements various sampling techniques for the initial conditions
  !and model parameters.  These routines should generally be called prior to
  !doing any integration.  Monte Carlo methodology.
  use potential, only : norm, pot, getH
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
    integer :: iso_N = 6
  end type
  type(ic_samp_flags) :: ic_flags

  integer :: param_sampling
  type :: param_samp_flags
    integer :: reg_constant = 1
    integer :: unif_prior = 2
    integer :: log_prior = 3
  end type
  type(param_samp_flags) :: param_flags
  real(dp), dimension(:,:), allocatable :: vp_prior_max, vp_prior_min
  logical :: use_first_priorval


  logical :: save_iso_N=.false.
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
      real(dp) :: H, rho
      integer :: i,j, kk

      logical :: with_velocity, with_eqen
      real(dp) :: beta, rand

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
      else
        call raise%fatal_code(&
            "The sampling technique hasn't been implemented.",&
            __FILE__, __LINE__)
      end if

      !Load initial vals from sample
      phi0 = y_background(1:num_inflaton)
      dphi0 = y_background(num_inflaton+1:2*num_inflaton)

      !Convert dphidt-->dphidN
      rho=0.5e0_dp*sum(dphi0*dphi0)+pot(phi0)
      H=sqrt(rho/3.0e0_dp)

      dphi0 = (1.0e0/H)*y_background(num_inflaton+1:2*num_inflaton)

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

              return

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

        if (potential_choice /= 1) then
          print*, "Just setting the IC to 0 -----"
          print*, "Fix me if this isn't correct.",  __FILE__, __LINE__
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

      real(dp) :: rand
      integer :: ii, jj

      !Don't need to marginalize over vparams
      if (param_sampling == param_flags%reg_constant) then
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

      else

        print*, "MODECODE: param_sampling=", param_sampling

        call raise%fatal_code(&
            "This choice of param_sampling is not supported.",&
            __FILE__, __LINE__)

      end if

    end subroutine get_vparams


end module modpk_icsampling
