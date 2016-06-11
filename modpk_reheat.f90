!A module that implements perturbative reheating after the end of
!inflation.

module modpk_reheat
  use modpkparams, only : dp, slowroll_infl_end, vparams, num_inflaton, Mpc2Mpl, &
    k_pivot, N_pivot, lna, nactual_bg, potential_choice, &
    auxparams
  use modpk_observables, only : power_spectra
  use internals, only : pi, k
  use potential, only : getH, geteps, getkineticenergy, &
    getw, dVdphi, getH_with_t, pot, d2Vdphi2
  use modpk_errorhandling, only : raise
  use csv_file, only : csv_write
  use modpk_numerics, only : array_polint
  use modpk_observables, only : observables
  use modpk_sampling, only : prior_auxparams_min, prior_auxparams_max

  implicit none

  type :: reheat_model_type
    !Use this module?
    logical :: use_reheat = .false.

    !Reheating options, only perturbative implemented currently
    integer :: reheat_model
    integer :: perturbative = 1

    !Sampling option for reheating decay parameters Gamma
    integer :: gamma_sampler
    integer :: uniform = 1
    integer :: fixed_mass = 2

  end type reheat_model_type
  type (reheat_model_type) :: reheat_opts

  type :: oscillation_counter
    real(dp), dimension(:), allocatable :: last_position
    integer, dimension(:), allocatable :: counter
    logical :: init_count=.false.
    integer :: total_oscillations = 5
    contains
      procedure, public :: init => osc_count_initializer
      procedure, public :: count_oscillations => &
        osc_count_count_oscillations

  end type oscillation_counter
  type(oscillation_counter) :: osc_count

  !Would bind options here, but can't use allocatable in namelist
  type :: reheat_state

    logical :: reheating_phase = .false.
    logical :: inflation_ended = .false.
    logical :: evolving_fluids = .false.

    complex(dp), dimension(:, :), allocatable :: q_horizcross
    real(dp), dimension(:, :), allocatable ::  eta_horizcross
    real(dp), dimension( :), allocatable :: phi_infl_end
    real(dp), dimension( :), allocatable :: dphi_infl_end
    real(dp) :: h_horizcross, eps_horizcross, efolds_end
    real(dp) :: h_end

    real(dp), dimension(:,:), allocatable :: c_ij_avg, c_ij_min, c_ij_max
    integer :: count_avg, count_moving_avg

    real(dp), dimension(:), allocatable :: dN, W_i
    real(dp), dimension(:), allocatable :: Gamma_i
    real(dp), dimension(:,:), allocatable :: r_ij, &
      rho_matter_decay, rho_radn_decay
    real(dp) :: hubble_prev
    real(dp), dimension(:), allocatable :: rho_matter_prev, rho_radn_prev
    logical, dimension(:), allocatable :: fields_decayed
    integer, dimension(:), allocatable :: fields_decay_order

    real(dp), dimension(:), allocatable :: y_radn_0

    real(dp) :: decay_factor_H

    type(observables) :: observs

    !DEBUG
    real(dp) :: Omega_phi, Omega_chi
    logical :: int_with_t

    type(power_spectra) :: pk, pk_hc

    contains

      procedure, public :: ending_check => reheat_ending_check
      procedure, public :: did_reheat_start => reheat_did_reheat_start
      procedure, public :: save_horizcross => reheat_save_horizcross
      procedure, public :: init => reheat_initializer
      procedure, public :: get_Gamma => reheat_get_Gamma
      procedure, public :: getH_with_radn => reheat_getH_with_radn
      procedure, public :: getdH_with_radn => reheat_getdH_with_radn
      procedure, public :: getr_ij => reheat_getr_ij
      procedure, public :: getW_i => reheat_getW_i
      procedure, public :: get_powerspectrum => reheat_get_powerspectrum

  end type reheat_state

  type (reheat_state) :: reheater

  contains

    function reheat_ending_check(self,phi,dphidN,q_modes,dqdN_modes,efolds) result(stopping)

      use potential, only : pot

      class(reheat_state) :: self

      real(dp), intent(in), dimension(:) :: phi, dphidN
      real(dp), dimension(size(phi)) :: test_phi
      complex(dp), intent(in), dimension(:), optional :: q_modes, dqdN_modes
      real(dp), intent(in), optional :: efolds

      real(dp) :: dot_phi0
      real(dp), dimension(size(dphidN)) :: omega
      integer :: ii, jj, kk
      complex(kind=dp), dimension(size(dphidN),size(dphidN)) :: ptb_matrix,&
        dptb_matrix

      logical :: stopping

      real(dp), dimension(size(phi)) :: dV, rho_i
      complex(dp), dimension(size(phi),size(phi)) :: xi, test_xi
      real(dp) :: lapse, hubble
      complex(dp), dimension(size(phi)) :: temp_vect
      complex(dp) :: q_horizcross_sr
      real(dp) :: pk_zeta, pk_zeta2

      complex(dp), dimension(size(phi)) :: xi_diag, zeta
      complex(dp), dimension(size(phi),size(phi)) :: C_ij
      complex(dp), dimension(size(phi),size(phi)) :: delta_rho
      complex(dp), parameter :: imag=(0.0e0_dp,1.0e0_dp)

      stopping = .false.

      !Consistency checks
      if (reheat_opts%use_reheat) then
        slowroll_infl_end = .false.
      end if

      if (.not. self%reheating_phase) then
        !Haven't tried calculating average yet
        self%count_avg = 0
        self%count_moving_avg = 0

      else if (present(q_modes) .and. present(dqdN_modes)) then

        !Convert hacked vector to matrix for modes
        !Builds the perturbation spectrum during reheating
        !ptb_matrix is \q_ij from 1410.0685
        !u_i = a \delta \phi_i = a \q_ij \hat a^j
        forall (ii=1:size(phi), jj=1:size(phi))
          ptb_matrix(ii,jj) = q_modes((ii-1)*size(phi)+jj)
          dptb_matrix(ii,jj) = dqdN_modes((ii-1)*size(phi)+jj)
        end forall

        !Convert to \zeta_i = -(H/\dot \rho_i)* \delta \rho_i
        !\zeta_i = \xi_ij \hat a^j
        !\zeta = (1/\dot \rho)*\sum_i \dot \rho_i \zeta_i
        hubble = getH(phi,dphidN)
        dV = dVdphi(phi)

        forall (ii=1:size(phi))
          temp_vect(ii) = sum(dphidN(:)*ptb_matrix(:,ii))
        end forall

        !forall (ii=1:size(phi), jj=1:size(phi))
        !  xi(ii,jj) =(1.0e0_dp/3.0e0_dp)*&
        !    (dptb_matrix(ii,jj)/dphidN(ii) &
        !    -0.5e0_dp*temp_vect(jj) &
        !    + dV(ii)*ptb_matrix(ii,jj)/hubble**2/dphidN(ii)**2)
        !end forall

        !Connection to Ewan and Joel's \zeta_i = C_ij \delta \phi^j
        !Assume q_ij is diagonal at horizon crossing...

        q_horizcross_sr = - imag*self%h_horizcross/&
            sqrt(2.0e0_dp*k**3)

        do ii=1,size(phi); do jj=1,size(phi)

          delta_rho(ii,jj) = ((hubble**2)*dphidN(ii)*dptb_matrix(ii,jj)&
            -(hubble**2)*dphidN(ii)**2*0.5e0_dp*temp_vect(jj)&
            +dV(ii)*ptb_matrix(ii,jj))

          !Ignores background pressure (w=0)
          test_phi = 0.0e0_dp
          test_phi(ii) = phi(ii)
          rho_i(ii) = 0.5e0_dp*hubble**2*dphidN(ii)**2 + pot(test_phi)

          C_ij(ii,jj) = delta_rho(ii,jj)/3.0e0_dp/rho_i(ii) &
            !/self%q_horizcross(jj,jj)
            /q_horizcross_sr

          !C_ij(ii,jj) = sqrt(C_ij(ii,jj)*conjg(C_ij(ii,jj)))

          !Sign chosen to match expectation from Meyers-Tarrant for
          !two-field N-flation
          C_ij(ii,jj) = (-real(C_ij(ii,jj))/abs(real(C_ij(ii,jj))))*&
            sqrt(C_ij(ii,jj)*conjg(C_ij(ii,jj)))

        end do; end do


        !Calculate the average of the C_ij during reheating period
        !Necessary for the Meyers-Tarrant matching
        !call average_cij(real(c_ij,dp))

        !To check
        !pk_zeta=0e0_dp
        !do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
        !  pk_zeta = pk_zeta +&
        !    dphidN(ii)**2*dphidN(jj)**2*xi(ii,kk)*conjg(xi(jj,kk))
        !end do; end do; end do
        !pk_zeta = pk_zeta*(k**3/2.0e0_dp/pi**2)/sum(dphidN**2)**2

        !DEBUG
        !print*, "this is pk_zeta", pk_zeta
        !print*, efolds, ',', pk_zeta, ',', (1.0/8.0/pi**2)*self%h_horizcross**2/self%eps_horizcross

        !To check
        !pk_zeta=0e0_dp
        !do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
        !  pk_zeta = pk_zeta +&
        !    dphidN(ii)**2*dphidN(jj)**2*C_ij(ii,kk)*conjg(C_ij(jj,kk))
        !    !dphidN(ii)**2*dphidN(jj)**2*real(C_ij(ii,kk))*real(C_ij(jj,kk))
        !end do; end do; end do
        !pk_zeta = pk_zeta/sum(dphidN**2)**2
        !pk_zeta = pk_zeta*((self%h_horizcross/pi)**2)/2.0e0_dp

        !!DEBUG
        !print*, "this is pk_zeta2", pk_zeta
        !stop


      end if

      select case(reheat_opts%reheat_model)
      case default

        call raise%fatal_code(&
          "Please specify your reheating model.",&
          __FILE__,__LINE__)

      case(1)

        !*********
        !IMPORTANT
        !*********
        !Hasn't been tested for cases where not initially in slow roll.  Will likely break!
        !Hasn't been tested for situations where need multiple samples.
        !*********
        !IMPORTANT
        !*********

        !Count the number of times each field has oscillated around its minimum
        !after inflation ends with epsilon>1

        !Don't stop if inflation hasn't ended
        if (geteps(phi,dphidN) < 1.0e0_dp .and. &
          .not. osc_count%init_count) then
          stopping = .false.
          return

        else

          !Start counting oscillations around the minimum phi_i=0
          if (.not. osc_count%init_count) then

            !Warning
            call raise%warning("Using reheat testing module.  &
              Be vewy, vewy careful!")

            call osc_count%init(phi)

          end if

          !Check to see if fields have passed zero
          call osc_count%count_oscillations(phi)

          !Calculate the average of the C_ij during reheating period
          !Necessary for the Meyers-Tarrant matching
          if (present(q_modes) .and.  all(osc_count%counter>0)) then
            call average_cij(real(c_ij,dp))
          end if

          !Stop after certain number of oscillations for every field
          if (minval(osc_count%counter) >= osc_count%total_oscillations) then
            stopping = .true.
            !Reset the counter
            osc_count%init_count=.false.
          end if

        end if

      end select


    contains

      !Calculate the average of the c_ij
      subroutine average_cij(c_ij)
        real(dp), dimension(:,:), intent(in) :: c_ij
        integer :: ii
        real(dp) :: rand

        if (.not. allocated(self%c_ij_avg)) &
          allocate(self%c_ij_avg(size(c_ij,1),size(c_ij,2)))

        !The stupid way seems to work best...
        !Avg = (max + min)/2
        do ii=1,size(c_ij,1); do jj=1,size(c_ij,2)
          if (c_ij(ii,jj) > self%c_ij_max(ii,jj)) &
            self%c_ij_max(ii,jj) = c_ij(ii,jj)

          if (c_ij(ii,jj) < self%c_ij_min(ii,jj)) &
            self%c_ij_min(ii,jj) = c_ij(ii,jj)
        end do; end do

        self%c_ij_avg = 0.5e0_dp*(self%c_ij_max + self%c_ij_min)

        !DEBUG
        !call csv_write(&
        !  6,&
        !  (/efolds, c_ij, self%c_ij_avg /), &
        !  advance=.true.)

      end subroutine average_cij


    end function reheat_ending_check

    !Check if inflation has ended and if we can start checking for
    !reheating physics.
    !Save some things for future use
    subroutine reheat_did_reheat_start(self, phi,dphidN,efolds,slowroll_start)

      class(reheat_state) :: self
      real(dp), dimension(:), intent(in) :: phi, dphidN
      real(dp), intent(in) :: efolds
      logical, intent(in) :: slowroll_start

      integer :: i, j

      if (.not. slowroll_start) then
        !If inflation hasn't started, reheating hasn't started
        self%reheating_phase = .false.
        self%inflation_ended = .false.
      else if( slowroll_start .and. &
        .not. self%reheating_phase .and. &
        getEps(phi, dphidN) .ge. 1) then
        !If you think you're inflating (slowroll_start) and
        !if not already in reheating (reheating_phase) and
        !if inflation has just ended (eps>1) then
        !you've just started reheating.

        !Save some values at end of inflation
        self%inflation_ended = .true.

        !Indicate started reheating
        self%reheating_phase = .true.

        self%efolds_end = efolds

        if (allocated(self%phi_infl_end)) &
          deallocate(self%phi_infl_end)
        allocate(self%phi_infl_end(size(phi)))
        self%phi_infl_end = phi

        if (allocated(self%dphi_infl_end)) &
          deallocate(self%dphi_infl_end)
        allocate(self%dphi_infl_end(size(phi)))
        self%dphi_infl_end = dphidN

        self%H_end = getH(phi,dphidN)


      else if (.not. self%inflation_ended) then
        !Otherwise just keep inflating
        self%reheating_phase = .false.
        self%inflation_ended = .false.
      end if

      if (self%inflation_ended .and. .not. self%reheating_phase) then
        call raise%fatal_cosmo(&
          'Inflation has ended, but the flag to indicate reheating &
          has started was not invoked.', &
          __FILE__, __LINE__)
      end if


    end subroutine reheat_did_reheat_start

    subroutine reheat_save_horizcross(self, phi, dphidN, q_modes)
      class(reheat_state) :: self
      real(dp), dimension(:), intent(in) :: phi, dphidN
      complex(dp), dimension(:) :: q_modes
      integer :: ii, jj

      if (allocated(self%q_horizcross)) &
        deallocate(self%q_horizcross)
      allocate(self%q_horizcross(num_inflaton,num_inflaton))
      if (allocated(self%eta_horizcross)) &
        deallocate(self%eta_horizcross)
      allocate(self%eta_horizcross(num_inflaton,num_inflaton))


      !ptb_matrix is \q_ij from 1410.0685
      !u_i = a \delta \phi_i = a \q_ij \hat a^j
      forall (ii=1:size(phi), jj=1:size(phi))
        self%q_horizcross(ii,jj) = q_modes((ii-1)*num_inflaton+jj)
      end forall

      self%h_horizcross = getH(phi,dphidN)
      self%eps_horizcross = getEps(phi,dphidN)
      self%eta_horizcross = d2Vdphi2(phi)/pot(phi)

    end subroutine reheat_save_horizcross

    subroutine reheat_initializer(self, reset_Gamma)
      class(reheat_state) :: self

      logical, intent(in) :: reset_Gamma

      if (allocated(self%q_horizcross))&
        deallocate(self%q_horizcross)
      if (allocated(self%phi_infl_end))&
        deallocate(self%phi_infl_end)
      if (allocated(self%dphi_infl_end))&
        deallocate(self%dphi_infl_end)
      if (allocated(self%c_ij_avg))&
        deallocate(self%c_ij_avg)

      if (allocated(self%dN))&
        deallocate(self%dN)
      if (allocated(self%W_i))&
        deallocate(self%W_i)

      if (reset_Gamma) then
        if (allocated(self%Gamma_i))&
          deallocate(self%Gamma_i)
      end if

      if (allocated(self%r_ij))&
        deallocate(self%r_ij)
      if (allocated(self%fields_decayed))&
        deallocate(self%fields_decayed)
      if (allocated(self%fields_decay_order))&
        deallocate(self%fields_decay_order)
      if (allocated(self%rho_matter_decay))&
        deallocate(self%rho_matter_decay)
      if (allocated(self%rho_radn_decay))&
        deallocate(self%rho_radn_decay)
      if (allocated(self%rho_matter_prev))&
        deallocate(self%rho_matter_prev)
      if (allocated(self%rho_radn_prev))&
        deallocate(self%rho_radn_prev)

      self%hubble_prev = 0.0e0_dp

      !Set min and max to ridiculously big/small values, respectively
      if (allocated(self%c_ij_min))&
        deallocate(self%c_ij_min)
      allocate(self%c_ij_min(num_inflaton, num_inflaton))
      self%c_ij_min= huge(1.0e0_dp)

      if (allocated(self%c_ij_max))&
        deallocate(self%c_ij_max)
      allocate(self%c_ij_max(num_inflaton, num_inflaton))
      self%c_ij_max=-huge(1.0e0_dp)

      if (allocated(self%q_horizcross)) &
        deallocate(self%q_horizcross)
      if (allocated(self%eta_horizcross)) &
        deallocate(self%eta_horizcross)

      self%h_horizcross=0.0e0_dp
      self%eps_horizcross=0.0e0_dp
      self%efolds_end=0.0e0_dp
      self%H_end=0.0e0_dp
      self%reheating_phase = .false.
      self%inflation_ended = .false.
      self%evolving_fluids = .false.

      self%count_avg = 0

      call self%pk%init()
      call self%pk_hc%init()

      !Say a field has decayed into radiation whenever
      !Gamma_i >= 3*H*decay_factor_H
      self%decay_factor_H = 1.0e0_dp/3.0e0_dp

      self%Omega_phi = 0.0_dp
      self%Omega_chi = 0.0_dp

      self%int_with_t = .false.

    end subroutine reheat_initializer


    !Sample the reheating decay parameters \Gamma_i uniformly
    subroutine reheat_get_Gamma(self)
      class(reheat_state) :: self

      real(dp), dimension(:), allocatable :: rand
      real(dp) :: R_ratio, R_min, R_max
      integer :: ii

      real(dp) :: A_var, B_var
      real(dp) :: m2_V(num_inflaton)

      !Don't have any Gamma_i yet
      if (.not. allocated(self%Gamma_i)) then
        allocate(self%Gamma_i(num_inflaton))
      else
        !Don't get new Gamma_i without deallocating first
        !ie, need to call init
        return
      end if

      if (reheat_opts%gamma_sampler == reheat_opts%fixed_mass) then
        !The decay rates are fixed wrt the masses with the relationship
        !\Gamma_i ~ A m_i^B

        if (size(auxparams)<2 &
          .or. size(prior_auxparams_max)<2 &
          .or. size(prior_auxparams_min)<2) then
          call raise%fatal_code(&
            'Not enough auxparams to sample reheat parameters.', &
            __FILE__, __LINE__)
        end if

        allocate(rand(2))
        call random_number(rand)
        rand = rand*(prior_auxparams_max-prior_auxparams_min)&
          +prior_auxparams_min
        A_var = 10.0**rand(1)
        B_var = rand(2)

        m2_V = 10.e0_dp**(vparams(1,:))


        self%Gamma_i(:) = A_var * self%H_end * sqrt(m2_V/maxval(m2_V))**B_var

        print*, "This is A and B", A_var, B_var
        print*, "This is Gamma_i", self%Gamma_i
        print*, "This is H_end", self%H_end

      else if (reheat_opts%gamma_sampler == reheat_opts%uniform) then

        allocate(rand(num_inflaton))
        call random_number(rand)
        call random_number(rand)


        if (potential_choice==21 .or. potential_choice==1) then

          !DEBUG
          if (num_inflaton>2) then
            print*, "testing only works for n=2"
            stop
          end if

          !R_min = 1e-4_dp
          !R_max = 1e-2_dp

          R_min = 1e-3_dp
          R_max = 1e2_dp

          call random_number(R_ratio)
          R_ratio = (log10(R_max) - log10(R_min))*R_ratio + log10(R_min)
          R_ratio = 10.0**R_ratio


          !R_ratio = 1.01
          !R_ratio = 1.0e-4
          R_ratio = 0.10545737E+02


          !DEBUG
          !print*, "this is R", R_ratio

          !self%Gamma_i(1) = 1e-14_dp
          !self%Gamma_i(1) = 5e-1 * R_min * self%H_end

          if (R_ratio<1.0) then
            self%Gamma_i(1) = 5e-1 * R_ratio * self%H_end
          else
            self%Gamma_i(1) = 1e-3 * self%H_end
          end if

          !DEBUG
          !self%Gamma_i(1) = 5e-1 *  self%H_end
          self%Gamma_i(2) = self%Gamma_i(1)/R_ratio

        end if

        !DEBUG
        !print*, "this is gamma", self%Gamma_i

      else
        call raise%fatal_code(&
          'Sampling technique for gamma should be set to 1 &
          to use this method.', &
          __FILE__, __LINE__)
      end if



    end subroutine reheat_get_Gamma

    function reheat_getH_with_radn(self, phi, dphi, rho_radn, use_t) &
        result(hubble)
      use potential
      class(reheat_state) :: self

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), intent(in) :: rho_radn
      logical, intent(in), optional :: use_t

      logical :: use_t_local

      real(dp) :: hubble

      if (present(use_t)) then
        use_t_local = use_t
      else
        use_t_local = .false.
      end if

      if (use_t_local) then

        !Cosmic time: dphi = dphidt
        hubble = sqrt((sum(0.5_dp*dphi**2) + pot(phi) + rho_radn)/3.0)

      else

        !E-folds: dphi = dphidN
        hubble = sqrt((rho_radn + pot(phi))/(3.0e0_dp - 0.5e0_dp*sum(dphi**2)))
      end if

    end function reheat_getH_with_radn

    function reheat_getdH_with_radn(self, phi, dphidN, rho_radn) &
        result(dhubble)
      use potential
      use modpk_deltaN, only :V_i_sum_sep
      class(reheat_state) :: self

      real(dp), dimension(:), intent(in) :: phi, dphidN
      real(dp), intent(in) :: rho_radn

      real(dp) :: dhubble, eps, hubble, V
      real(dp) :: denom, numer
      real(dp), dimension(size(phi)) :: rho_fields, dV

      !E-folds

      hubble = self%getH_with_radn(phi,dphidN,rho_radn)
      eps = geteps(phi,dphidN)
      rho_fields = 0.5e0_dp*hubble**2*dphidN**2 + V_i_sum_sep(phi)
      V = pot(phi)
      dV = dVdphi(phi)

      denom = 2.0e0_dp*hubble*(3.0e0_dp - eps)**2 &
        +(2.0e0_dp*eps/hubble)*(V + rho_radn)
      numer = (3.0e0_dp - eps)*&
        (sum(dV*dphidN) - 4.0e0_dp*rho_radn +&
          sum(self%Gamma_i*rho_fields)) &
        -&
        (V+rho_radn)*(&
        sum(dphidN**2 * (3.0e0_dp + self%Gamma_i/hubble)) &
        + sum(dphidN*dV/hubble**2))

      dhubble = numer/denom


    end function reheat_getdH_with_radn

    subroutine reheat_getr_ij(self, rho_fields, rho_matter, rho_radn, &
        all_fluid_descr)
      class (reheat_state) :: self

      real(dp), dimension(:), intent(in) :: rho_fields, rho_radn, rho_matter
      logical, intent(in) :: all_fluid_descr

      real(dp), dimension(size(rho_fields)) :: rho_error, real_rho_matter
      logical, dimension(size(rho_fields)) :: save_field_rij
      real(dp) :: hubble, rho_decayed, rho_notdecayed
      integer :: ii, jj

      real(dp), dimension(2) :: hubble_vect
      real(dp), dimension(size(rho_fields),2) :: rho_mat

      !Using the sudden decay approximation where scalar fields
      !are taken to have decayed if H<\Gamma_i, even if they
      !have some energy left in the scalar field energy

      if (.not. self%evolving_fluids) then
        call raise%fatal_code(&
          "Need to be evolving with radiation fluid &
          to use this subroutine.",&
          __FILE__,__LINE__)
      end if

      call initialize()

      if (all_fluid_descr) then
        !Ignore the rho_fields and use only rho_matter
        real_rho_matter = rho_matter
      else
        !Ignore the rho_matter and use only rho_fields
        real_rho_matter = rho_fields
      end if
      hubble = sqrt( ( sum(real_rho_matter) + sum(rho_radn))/3.0e0_dp)

      !Find which fields have decayed
      save_field_rij = .false.
      do jj=1, size(rho_fields)

        if (3.0e0_dp*hubble*self%decay_factor_H &
          <= self%Gamma_i(jj)) then
          if (.not. self%fields_decayed(jj)) then
            !This is first time step after decaying
            !Save everything
            save_field_rij(jj) = .true.

            !Save energy density in fields and radn
            !at decay time for this field as a vector
            !NB: Use polynomial interpolation on this time step and
            !previous to get more accurate value
            hubble_vect(1) = &
              3.0e0_dp*self%hubble_prev*self%decay_factor_H
            hubble_vect(2) = &
              3.0e0_dp*hubble*self%decay_factor_H
            rho_mat(:,1) = self%rho_matter_prev
            rho_mat(:,2) = real_rho_matter
            call array_polint(hubble_vect, &
              rho_mat, &
              self%Gamma_i(jj),&
              self%rho_matter_decay(:,jj), &
              rho_error)
            if (any(rho_error >1e-3)) then
              call raise%fatal_code(&
                "Interpolation error.",&
                __FILE__,__LINE__)
            end if

            rho_mat(:,1) = self%rho_radn_prev
            rho_mat(:,2) = rho_radn
            call array_polint(hubble_vect, &
              rho_mat, &
              self%Gamma_i(jj),&
              self%rho_radn_decay(:,jj), &
              rho_error)
            if (any(rho_error >1e-3)) then
              call raise%fatal_code(&
                "Interpolation error.",&
                __FILE__,__LINE__)
            end if

          end if


        else
          self%fields_decayed(jj) = .false.
          save_field_rij(jj) = .false.

        end if
      end do

      !Save new energy densities
      self%rho_matter_prev = real_rho_matter
      self%rho_radn_prev = rho_radn
      self%hubble_prev = hubble


      !Calculate r_ij only at the time step where the field decays
      !into radiation, ie, whenever save_field_rij is .true.

      !Evaluated at decay time t_jj, ii indexes energy components at
      !this time
      do jj=1, size(rho_fields)

        if (save_field_rij(jj)) then

          !Find energy density decayed and still in fields
          !Add radiation energy (decayed) and
          !matter energy (notdecayed) separately
          rho_decayed = 0.0e0_dp
          rho_notdecayed = 0.0e0_dp
          do ii=1,size(rho_fields)
            if (self%fields_decayed(ii)) then
              rho_decayed = rho_decayed + &
                self%rho_matter_decay(ii,jj) + &
                self%rho_radn_decay(ii,jj)
            else
              rho_notdecayed = rho_notdecayed + &
                self%rho_matter_decay(ii,jj) + &
                self%rho_radn_decay(ii,jj)
            end if
          end do


          self%r_ij(:,jj) = 3.0e0_dp* &
            (self%rho_matter_decay(:,jj) +self%rho_radn_decay(:,jj))/&
            (4.0e0_dp*rho_decayed + 3.0e0_dp*rho_notdecayed)

          !Label this field as having decayed into radiation
          !for next time
          self%fields_decayed(jj) = .true.
          self%fields_decay_order(jj) = count(self%fields_decayed)

          !DEBUG
          print*, "This field decayed:", jj

        end if

      end do

      contains

        subroutine initialize()

          !Start counting which fields have decayed
          if (.not. allocated(self%fields_decayed)) then
            allocate(self%fields_decayed(size(rho_fields)))
            self%fields_decayed = .false.
          end if
          if (.not. allocated(self%fields_decay_order)) then
            allocate(self%fields_decay_order(size(rho_fields)))
            self%fields_decay_order = -1
          end if
          if (.not. allocated(self%rho_matter_decay)) then
            allocate(self%rho_matter_decay(size(rho_fields),size(rho_fields)))
            self%rho_matter_decay = 0.0e0_dp
          end if
          if (.not. allocated(self%rho_radn_decay)) then
            allocate(self%rho_radn_decay(size(rho_fields),size(rho_fields)))
            self%rho_radn_decay = 0.0e0_dp
          end if
          if (.not. allocated(self%r_ij)) then
            allocate(self%r_ij(size(rho_fields),size(rho_fields)))
            self%r_ij = 0.0e0_dp
          end if
          if (.not. allocated(self%rho_matter_prev)) then
            allocate(self%rho_matter_prev(size(rho_fields)))
            self%rho_matter_prev = rho_fields
          end if
          if (.not. allocated(self%rho_radn_prev)) then
            allocate(self%rho_radn_prev(size(rho_fields)))
            self%rho_radn_prev = rho_radn
          end if
        end subroutine initialize


    end subroutine reheat_getr_ij

    subroutine reheat_getW_i(self)
      class(reheat_state) :: self

      real(dp) :: sum_term, product_term
      real(dp), dimension(num_inflaton,num_inflaton) :: r_proper
      real(dp), dimension(num_inflaton) :: W_i_temp
      real(dp), dimension(0:num_inflaton-1) :: a_vect

      integer :: ii, jj, kk
      real(dp) :: summer

      if (.not. self%evolving_fluids) then
        call raise%fatal_code(&
          "Need to be evolving with radiation fluid &
          to use this subroutine.",&
          __FILE__,__LINE__)
      end if

      if (.not. allocated(self%r_ij)) then
        call raise%fatal_code(&
          "Need to evaluate the r_ij &
          before the W_i.",&
          __FILE__,__LINE__)
      end if

      if (.not. all(self%fields_decayed)) then
        call raise%fatal_code(&
          "Not all fields have decayed to radiation &
          before calculating the W_i.",&
          __FILE__,__LINE__)
      end if

      !Initialize
      if (.not. allocated(self%W_i)) then
        allocate(self%W_i(num_inflaton))
        self%W_i = 0.0e0_dp
      end if

      !The definition of the W_i relies on ordering the fields
      !according to when they decay to radiation.
      !We want to identify the "right" W_i by this criterion, then
      !reorder them back so that they match the same field basis
      !we are using in the rest of the code.

      !Define the "properly ordered" r_ij
      do ii=1,size(self%fields_decay_order)
        do jj=1,size(self%fields_decay_order)
          r_proper(&
            self%fields_decay_order(ii), &
            self%fields_decay_order(jj)) = &
            self%r_ij(ii,jj)
        end do
      end do

      !Updated 12/17/15
      call define_a_vect(a_vect)
      do ii=1, num_inflaton
        summer = 0.0_dp
        do jj=0, num_inflaton-1
          summer = summer + &
            a_vect(jj)*r_proper(ii,num_inflaton-jj)
        end do

        self%W_i(ii) = summer

      end do

      !Reorder the W_i so that it aligns with the original field
      !space basis
      W_i_temp = self%W_i
      do ii=1, num_inflaton
        self%W_i(self%fields_decay_order(ii)) = W_i_temp(ii)
      end do

      contains

        subroutine define_a_vect(a_vect)

          real(dp), dimension(0:num_inflaton-1), intent(out) :: a_vect
          integer :: jj, kk
          real(dp) :: summer

          a_vect(0) = 1.0_dp

          do jj=1, num_inflaton-1
            summer = 0.0_dp
            do kk=0, jj-1
              summer = summer + &
                a_vect(kk)*r_proper(num_inflaton-jj,num_inflaton-kk)
            end do
            a_vect(jj) = summer/3.0_dp
          end do

        end subroutine define_a_vect

    end subroutine reheat_getW_i

    !Calculate the power spectrum in the sudden decay approximation for
    !perturbative reheating
    subroutine reheat_get_powerspectrum(self)
      class(reheat_state) :: self

      real(dp), dimension(num_inflaton) :: dN
      real(dp)  :: pk_horizcross, sum_term
      integer :: jj, ii
      integer :: jj_decay, ii_decay

      !DEBUG
      real(dp) :: phi_piv_test, chi_piv_test, X_test
      real(dp) :: dV_test(num_inflaton)
      real(dp) :: v_debug, q_debug, p_debug, R_debug

      if (.not. self%evolving_fluids) then
        call raise%fatal_code(&
          "Need to be evolving with radiation fluid &
          to use this subroutine.",&
          __FILE__,__LINE__)
      end if

      if (.not. allocated(self%c_ij_avg)) then
        call raise%fatal_code(&
          "Need to evaluate the C_ij &
          before the power spectrum.",&
          __FILE__,__LINE__)
      end if

      !DEBUG
      !if (self%Gamma_i(1)/self%Gamma_i(2)<1.0_dp) then
      !  !print*, "this is r_ij", self%r_ij(1,1), self%r_ij(1,2), self%r_ij(2,1), self%r_ij(2,2)
      !  v_debug = 0.60_dp
      !  q_debug = 0.63_dp*log(self%Omega_phi)/log(1.0_dp-self%Omega_phi)
      !  p_debug = ((4.0_dp*self%Omega_phi)/(3.0_dp+self%Omega_phi))**(-1.0_dp/v_debug) - q_debug
      !  R_debug = self%Gamma_i(1)/self%Gamma_i(2)

      !  print*, "this is r_ij prediction for R<1:",&
      !  1.0_dp - (p_debug + q_debug/R_debug)**(-v_debug),&
      !  self%Omega_phi,&
      !  abs(self%Omega_phi-self%Omega_chi),&
      !  self%Omega_chi

      !else
      !  !print*, "this is r_ij", self%r_ij(1,1), self%r_ij(1,2), self%r_ij(2,1), self%r_ij(2,2)

      !  v_debug = 1.6666_dp
      !  q_debug = (4.0*(1.0-self%Omega_phi)/(4.0-self%Omega_phi))**(-1.0/v_debug) - 1.0
      !  R_debug = self%Gamma_i(1)/self%Gamma_i(2)

      !  print*, "this is r_ij prediction for R>1:", &
      !    self%Omega_chi, &
      !    abs(self%Omega_phi-self%Omega_chi), &
      !    self%Omega_phi, &
      !    1.0_dp - (1.0 + q_debug*sqrt(R_debug))**(-v_debug)


      !end if

      !DEBUG
      !print*, "OVERRIDING r_ij!"
      !self%r_ij(1,1) = self%r_ij(2,1)
      !self%r_ij(1,2) = self%r_ij(1,1)
      !self%r_ij(2,1) = self%r_ij(1,2)
      !self%r_ij(2,2) = self%r_ij(2,2)

      !Get the W_i vector
      call self%getW_i()

      !Calculate derivatives of e-foldings
      !in "our" basis
      dN = 0.0e0_dp
      do jj=1,num_inflaton
        dN(:) = dN(:) + &
          self%W_i(jj)*self%c_ij_avg(jj,:)
      end do



      !DEBUG
      !phi_piv_test = 0.499002
      !phi_piv_test = 0.1
      !phi_piv_test=1.60830E-01
      !chi_piv_test = 1.5006E+01
      !C_{\phi\phi} = 2 / (3 \phi_*),
      !C_{\phi\chi} \approx \chi_* / (2Mp^2),
      !C_{\chi\phi} = 0,                      
      !C_{\chi\chi} \approx \chi_* / (2Mp^2),



      !DEBUG
      !print*, "this is dN:", dN

      !print*, "OVERRIDING dN"
      !dN(1) = 0.75074301E+01
      !dN(2) = 0.18677625E+03
      !print*, "this is expected:", "???", chi_piv_test/2.0
      !print*, "these are c's", self%c_ij_avg(1,1), self%c_ij_avg(1,2), self%c_ij_avg(2,1), self%c_ij_avg(2,2)
      !print*, "this is c_ij prediction:", 7.50, -4.8, 7.50, 476.0
      !print*, "this is W_i", self%W_i
      !print*, "this is W_i summed",self%W_i(1) + self%W_i(2)

      !Load the power spectrum attributes
      self%pk%k = self%pk_hc%k

      !Assuming SR
      !NB: pk_hc%adiab is actually the *field* power spectrum ~ H^2/4 pi^2
      self%pk%adiab = self%pk_hc%adiab*sum(dN**2)
      self%pk%tensor = 8.0e0_dp*self%pk_hc%adiab

      sum_term = 0.0e0_dp
      do ii=1,num_inflaton; do jj=1,num_inflaton
        sum_term = sum_term + &
          self%eta_horizcross(ii,jj)*dN(ii)*dN(jj)
      end do; end do
      print*, "this is R, r, ns:", self%Gamma_i(1)/self%Gamma_i(2), &
        self%pk%tensor/self%pk%adiab, &
      !print*, "this is r, ns:", 8.0/sum(dN**2), &
        1.0e0_dp&
        -2.0e0_dp*self%eps_horizcross &
        -(2.0e0_dp/sum(dN**2))*&
        (1.0e0_dp -  sum_term)


      call self%observs%set_zero()
      self%observs%As = self%pk%adiab
      self%observs%ns = 1.0e0_dp&
        -2.0e0_dp*self%eps_horizcross &
        -(2.0e0_dp/sum(dN**2))*&
        (1.0e0_dp -  sum_term)
      self%observs%r  = self%pk%tensor/self%pk%adiab


      !No isocurvature
      self%pk%isocurv = 0.0e0_dp
      self%pk%cross_ad_iso = 0.0e0_dp
      self%pk%pnad = 0.0e0_dp
      self%pk%pressure = 0.0e0_dp
      self%pk%press_ad = 0.0e0_dp
      self%pk%entropy = 0.0e0_dp
      self%pk%bundle_exp_scalar = 0.0e0_dp


      !DEBUG
      !print*, "Stopping here"
      !stop

    end subroutine reheat_get_powerspectrum

    subroutine osc_count_initializer(self, phi)
      class(oscillation_counter) :: self
      real(dp), dimension(:), intent(in) :: phi

      !Start counting oscillations
      osc_count%init_count = .true.

      !Make the counting arrays
      if (allocated(osc_count%last_position)) deallocate(osc_count%last_position)
      if (allocated(osc_count%counter)) deallocate(osc_count%counter)
      allocate(osc_count%counter(size(phi)))
      allocate(osc_count%last_position(size(phi)))

      osc_count%counter = 0
      osc_count%last_position = phi

    end subroutine osc_count_initializer

    subroutine osc_count_count_oscillations(self, phi)
      class(oscillation_counter) :: self
      real(dp), dimension(:), intent(in) :: phi

      integer :: ii
      integer :: sign_phi, sign_last

      do ii=1,size(phi)
        sign_phi = int(phi(ii)/abs(phi(ii)))
        sign_last = int(osc_count%last_position(ii)/abs(osc_count%last_position(ii)))
        if (sign_phi /= sign_last ) then
          osc_count%counter(ii) =osc_count%counter(ii)+1
        end if
      end do

      !Save latest field-space position
      osc_count%last_position = phi

    end subroutine osc_count_count_oscillations

end module modpk_reheat
