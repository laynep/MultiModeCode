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
    logical :: reheating_phase
  end type reheat_state
  type (reheat_state) :: reheat_saver

  contains

    function reheat_ending_check(phi,dphidN,q_modes,dqdN_modes,efolds) result(stopping)

      use potential, only : pot

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
      !real(dp), dimension(size(phi),size(phi)) :: C_ij
      complex(dp), parameter :: imag=(0.0e0_dp,1.0e0_dp)

      stopping = .false.

      call consistency_checks()

      !Convert hacked vector to matrix for modes
      !Builds the perturbation spectrum during reheating
      if (present(q_modes) .and. present(dqdN_modes)) then

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

        forall (ii=1:size(phi), jj=1:size(phi))
          xi(ii,jj) =(1.0e0_dp/3.0e0_dp)*&
            (dptb_matrix(ii,jj)/dphidN(ii) &
            -0.5e0_dp*temp_vect(jj) &
            + dV(ii)*ptb_matrix(ii,jj)/hubble**2/dphidN(ii)**2)
        end forall

        !Connection to Ewan and Joel's \zeta_i = C_ij \delta \phi^j
        !Assume q_ij is diagonal at horizon crossing...

        q_horizcross_sr = - imag*reheat_saver%h_horizcross/&
            sqrt(2.0e0_dp*k**3)

        !forall (ii=1:size(phi), jj=1:size(phi))
        do ii=1,size(phi); do jj=1,size(phi)

          delta_rho(ii,jj) = ((hubble**2)*dphidN(ii)*dptb_matrix(ii,jj)&
            -(hubble**2)*dphidN(ii)**2*0.5e0_dp*temp_vect(jj)&
            +dV(ii)*ptb_matrix(ii,jj))

          !Ignores background pressure (w=0)
          test_phi = 0.0e0_dp
          test_phi(ii) = phi(ii)
          rho_i(ii) = 0.5e0_dp*hubble**2*dphidN(ii)**2 + pot(test_phi)

          C_ij(ii,jj) = delta_rho(ii,jj)/3.0e0_dp/rho_i(ii) &
            !/reheat_saver%q_horizcross(jj,jj)
            /q_horizcross_sr

          C_ij(ii,jj) = sqrt(C_ij(ii,jj)*conjg(C_ij(ii,jj)))

        !end forall
        end do; end do


        !DEBUG
        call csv_write(&
          6,&
          (/efolds, real(C_ij(1,1)), real(C_ij(1,2)), real(C_ij(2,1)), real(C_ij(2,2)) /)   , &
          advance=.true.)

        !To check
        !pk_zeta=0e0_dp
        !do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
        !  pk_zeta = pk_zeta +&
        !    dphidN(ii)**2*dphidN(jj)**2*xi(ii,kk)*conjg(xi(jj,kk))
        !end do; end do; end do
        !pk_zeta = pk_zeta*(k**3/2.0e0_dp/pi**2)/sum(dphidN**2)**2

        !!DEBUG
        !print*, "this is pk_zeta"
        !print*, efolds, ',', pk_zeta, ',', (1.0/8.0/pi**2)*reheat_saver%h_horizcross**2/reheat_saver%eps_horizcross

        !To check
        !pk_zeta=0e0_dp
        !do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
        !  pk_zeta = pk_zeta +&
        !    dphidN(ii)**2*dphidN(jj)**2*C_ij(ii,kk)*conjg(C_ij(jj,kk))
        !    !dphidN(ii)**2*dphidN(jj)**2*real(C_ij(ii,kk))*real(C_ij(jj,kk))
        !end do; end do; end do
        !pk_zeta = pk_zeta/sum(dphidN**2)**2
        !pk_zeta = pk_zeta*((reheat_saver%h_horizcross/pi)**2)/2.0e0_dp

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



    end function reheat_ending_check

    subroutine consistency_checks()

      if (use_reheat) then
        slowroll_infl_end = .false.
      end if

    end subroutine consistency_checks

end module modpk_reheat
