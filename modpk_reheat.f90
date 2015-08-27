!A module that implements perturbative reheating after the end of
!inflation.

module modpk_reheat
  use modpkparams, only : dp, slowroll_infl_end, vparams, num_inflaton, Mpc2Mpl, &
    k_pivot, N_pivot
  use internals, only : pi, k
  use potential, only : getH, geteps, getkineticenergy, getw, dVdphi
  use modpk_errorhandling, only : raise
  use csv_file, only : csv_write

  implicit none

  logical :: use_reheat=.false.
  type :: reheat_model_flags
    integer :: perturbative = 1
  end type reheat_model_flags
  integer :: reheat_model

  type :: oscillation_counter
    real(dp), dimension(:), allocatable :: last_position
    integer, dimension(:), allocatable :: counter
    logical :: init_count=.false.
  end type oscillation_counter
  type(oscillation_counter) :: osc_count

  type :: reheat_state
    complex(dp), dimension(:, :), allocatable :: q_horizcross
    real(dp), dimension( :), allocatable :: phi_infl_end
    real(dp) :: h_horizcross, eps_horizcross, efolds_end
    logical :: reheating_phase, inflation_ended

    real(dp), dimension(:,:), allocatable :: c_ij_avg, c_ij_min, c_ij_max
    real(dp), dimension(:,:,:), allocatable :: c_ij_moving_avg
    integer :: count_avg, count_moving_avg
    real(dp) :: moving_avg_tol = 1.0e-2_dp

    contains

      procedure, public :: ending_check => reheat_ending_check
      procedure, public :: did_reheat_start => reheat_did_reheat_start
      procedure, public :: save_horizcross => reheat_save_horizcross
      procedure, public :: init => reheat_initializer

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
      integer :: sign_phi, sign_last

      integer, parameter :: total_oscillations = 5

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

      call consistency_checks()

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
        call average_cij(real(c_ij,dp))

        !To check
        !pk_zeta=0e0_dp
        !do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
        !  pk_zeta = pk_zeta +&
        !    dphidN(ii)**2*dphidN(jj)**2*xi(ii,kk)*conjg(xi(jj,kk))
        !end do; end do; end do
        !pk_zeta = pk_zeta*(k**3/2.0e0_dp/pi**2)/sum(dphidN**2)**2

        !!DEBUG
        !print*, "this is pk_zeta"
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

        !DEBUG
        !print*, "this is pk_zeta2", pk_zeta





      end if


      select case(reheat_model)
      case default

        call raise%fatal_code("Please specify your reheating model.",__FILE__,__LINE__)

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
              Not safe!  &
              Be vewy, vewy careful.")

            !Start counting oscillations
            osc_count%init_count = .true.

            !Save inflation ending position in order to normalize modes
            !slowroll_infl_end = .true.

            !Make the counting arrays
            if (allocated(osc_count%last_position)) deallocate(osc_count%last_position)
            if (allocated(osc_count%counter)) deallocate(osc_count%counter)
            allocate(osc_count%counter(size(phi)))
            allocate(osc_count%last_position(size(phi)))

            osc_count%counter = 0
            osc_count%last_position = phi

          end if

          !Check to see if fields have passed zero
          do ii=1,size(phi)
            sign_phi = int(phi(ii)/abs(phi(ii)))
            sign_last = int(osc_count%last_position(ii)/abs(osc_count%last_position(ii)))
            if (sign_phi /= sign_last ) then
              osc_count%counter(ii) =osc_count%counter(ii)+1
              !print*, "counting new one", ii, osc_count%counter(ii)
            end if
          end do

          !Save latest field-space position
          osc_count%last_position = phi

          !Stop after certain number of oscillations for every field
          if (minval(osc_count%counter) >= total_oscillations) then
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
        if (.not. allocated(self%c_ij_moving_avg)) &
          allocate(self%c_ij_moving_avg(1000,size(c_ij,1),size(c_ij,2)))

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
        print*, "this is c_ij_avg", self%c_ij_avg(1,2), c_ij(1,2)


        !if (self%count_avg == 0) then
        !  self%c_ij_avg(:,:) = c_ij(:,:)
        !else
        !  self%c_ij_avg(:,:) = (c_ij(:,:) + &
        !    self%count_avg*self%c_ij_avg(:,:))&
        !    /(self%count_avg + 1)
        !end if

        !call random_number(rand)
        !if (rand < self%moving_avg_tol) then
        !  !Reset the count for individual average
        !  self%count_avg = 0

        !  !Perform binning on average
        !  self%count_moving_avg = self%count_moving_avg + 1

        !  self%c_ij_moving_avg(self%count_moving_avg,:,:) = &
        !    self%c_ij_avg(:,:)

        !  !self%c_ij_moving_avg(self%count_moving_avg,:,:) = &
        !  !  (self%c_ij_avg(:,:) + &
        !  !  (self%count_moving_avg-1)*self%c_ij_moving_avg(:,:))&
        !  !  /(self%count_moving_avg)
        !else
        !  self%count_avg = self%count_avg + 1
        !end if

        !if (self%count_moving_avg ==0) &
        !  self%c_ij_moving_avg(1,:,:)=self%c_ij_avg(:,:)

        !call csv_write(&
        !  6,&
        !  (/efolds, real(C_ij(1,1)), real(C_ij(1,2)), &
        !  real(C_ij(2,1)), real(C_ij(2,2)) /)   , &
        !  advance=.false.)

        !call csv_write(&
        !  6,&
        !  dphidN   , &
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
        self%efolds_end = efolds
        if (allocated(self%phi_infl_end)) &
          deallocate(self%phi_infl_end)
        allocate(self%phi_infl_end(size(phi)))
        self%phi_infl_end = phi

        self%reheating_phase = .true.

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


      !ptb_matrix is \q_ij from 1410.0685
      !u_i = a \delta \phi_i = a \q_ij \hat a^j
      forall (ii=1:size(phi), jj=1:size(phi))
        self%q_horizcross(ii,jj) = q_modes((ii-1)*num_inflaton+jj)
      end forall

      self%h_horizcross = getH(phi,dphidN)
      self%eps_horizcross = getEps(phi,dphidN)

    end subroutine reheat_save_horizcross

    subroutine reheat_initializer(self)
      class(reheat_state) :: self

      if (allocated(self%q_horizcross))&
        deallocate(self%q_horizcross)
      if (allocated(self%phi_infl_end))&
        deallocate(self%phi_infl_end)
      if (allocated(self%c_ij_avg))&
        deallocate(self%c_ij_avg)
      if (allocated(self%c_ij_moving_avg))&
        deallocate(self%c_ij_moving_avg)

      !Set min and max to ridiculously big/small values
      if (allocated(self%c_ij_min))&
        deallocate(self%c_ij_min)
      allocate(self%c_ij_min(num_inflaton, num_inflaton))
      self%c_ij_min= huge(1.0e0_dp)

      if (allocated(self%c_ij_max))&
        deallocate(self%c_ij_max)
      allocate(self%c_ij_max(num_inflaton, num_inflaton))
      self%c_ij_max=-huge(1.0e0_dp)

      self%h_horizcross=0.0e0_dp
      self%eps_horizcross=0.0e0_dp
      self%efolds_end=0.0e0_dp
      self%reheating_phase = .false.
      self%inflation_ended = .false.


      self%count_avg = 0
      self%count_moving_avg = 0

    end subroutine reheat_initializer

    subroutine consistency_checks()

      if (use_reheat) then
        slowroll_infl_end = .false.
      end if

    end subroutine consistency_checks

end module modpk_reheat
