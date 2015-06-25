!A module that implements perturbative reheating after the end of
!inflation.

module modpk_reheat
  use modpkparams, only : dp, slowroll_infl_end, vparams, num_inflaton, Mpc2Mpl, k_pivot
  use internals, only : pi, k
  use potential, only : getH, geteps, getkineticenergy, getw, dVdphi
  use modpk_errorhandling, only : raise

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

  contains

    function reheat_ending_check(phi,dphidN,q_modes,dqdN_modes,efolds) result(stopping)

      real(dp), intent(in), dimension(:) :: phi, dphidN
      complex(dp), intent(in), dimension(:), optional :: q_modes, dqdN_modes
      real(dp), intent(in), optional :: efolds

      real(dp) :: dot_phi0
      real(dp), dimension(size(dphidN)) :: omega
      integer :: ii, jj, kk
      complex(kind=dp), dimension(size(dphidN),size(dphidN)) :: ptb_matrix,&
        dptb_matrix

      logical :: stopping
      integer :: sign_phi, sign_last

      real(dp), dimension(size(phi)) :: dV, rho_i
      complex(dp), dimension(size(phi),size(phi)) :: xi
      real(dp) :: lapse, hubble
      complex(dp), dimension(size(phi)) :: temp_vect
      real(dp) :: pk_zeta
      complex(dp), dimension(size(phi)) :: xi_diag
      complex(dp), dimension(size(phi),size(phi)) :: C_ij
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
            - 0.5e0_dp*temp_vect(jj) &
            + dV(ii)*ptb_matrix(ii,jj)/hubble**2/dphidN(ii)**2)
        end forall

        !To check
        pk_zeta=0e0_dp
        do ii=1,size(phi); do jj=1,size(phi); do kk=1,size(phi)
          pk_zeta = pk_zeta +&
            dphidN(ii)**2*dphidN(jj)**2*xi(ii,kk)*conjg(xi(jj,kk))
        end do; end do; end do
        pk_zeta = pk_zeta*(k**3/2.0e0_dp/pi**2)/sum(dphidN**2)**2

        !Connection to Ewan and Joel's \zeta_i = C_ij \delta \phi^j
        C_ij=(imag*exp(-imag))*xi*sqrt(2.0e0_dp*k**3)/hubble

        !DEBUG
        print*, "from here"
        !if (present(efolds)) print*, efolds, C_ij


        !print*, ptb_matrix/(hubble*exp(imag)/sqrt(2.0e0_dp*k**3))
        !print*, ptb_matrix*((-imag*exp(-imag))*sqrt(2.0e0_dp*k**3)/hubble)

        C_ij=ptb_matrix*((imag*exp(-imag))*sqrt(2.0e0_dp*k**3)/hubble)*(1.0e0_dp+imag)
        do ii=1,size(phi)
          print*, "c", C_ij(ii,ii)
          print*, "doesn't work for imag portion..."
        end do



        !print*, "this is pk_zeta", pk_zeta
        print*, "this is k", k/Mpc2Mpl, k_pivot
        stop

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
              Not entirely safe!  &
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

          !Stop after at least 5 oscillations for every field
          if (minval(osc_count%counter) >= 5) then
            stopping = .true.
            !Reset the counter
            osc_count%init_count=.false.
          end if

        end if





      end select


      !print*, "this is rho:", &
      !print*,((1.0e0_dp/vparams(2,1))*vparams(1,ii)*abs(phi(ii))**(vparams(2,1)) + &
        !0.5e0_dp*dphidN(ii)**2*getH(phi,dphidN)**2,ii=1,size(phi))

      !DEBUG
      !dot_phi0 = sqrt(dot_product(dphidN,dphidN))
      !omega = dphidN/ sqrt(dot_product(dphidN,dphidN))
      !print*, "omega", omega
      !print*, "-H/dot_phi", -getH(phi,dphidN)/dot_phi0
      !stop

      !print*,sum((1.0e0_dp/vparams(2,1))*vparams(1,:)*abs(phi(:))**(vparams(2,1)) + &
      !  0.5e0_dp*dphidN(:)**2*getH(phi,dphidN)**2)

      !print*, "eps", geteps(phi,dphidN)
      !print*, "w", getw(phi,dphidN)
      !print*, "m/H",sqrt(vparams(1,:))/getH(phi,dphidN)


      !print*, geteps(phi,dphidN), getw(phi,dphidN)

      !print*, "phi0", phi
      !print*, "dphi0", dphidN

      !if (geteps(phi,dphidN)>1.0) stop




    end function reheat_ending_check

    subroutine consistency_checks()

      if (use_reheat) then
        slowroll_infl_end = .false.
      end if

    end subroutine consistency_checks

end module modpk_reheat
