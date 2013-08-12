program test_mmodpk
  use modpkparams
  use potential
  use background_evolution
  use modpk_utils
  use camb_interface
  use ode_path
  use access_modpk
  use internals

  implicit none

  !Writing fmts
  character(16), parameter:: e_fmt = '(a25, es12.4)'
  character(36), parameter:: e2_fmt = '(a25, es12.4, a3, es11.4, a1)'
  character(16), parameter:: i_fmt = '(a25,I3)'
  character(16) :: array_fmt
  character(len=2) :: ci

  !Run-specific input params
  integer :: i, vparam_arrays

  !Cosmology
  real(dp) :: dlnk, As, ns, nt, r

  !Sampling parameters for ICs
  integer, parameter :: test_samp=1, eqen_samp=2
  integer :: sampling_techn
  real(dp) :: energy_scale

  integer :: u

  !For run-time alloc w/out re-compile
  namelist /init/ num_inflaton, potential_choice, &
    modpkoutput, slowroll_infl_end, instreheat, vparam_arrays

  !namelist /ic_sampling/ sampling_techn, energy_scale

  namelist /params/ phi_init0, vparams, &
    N_pivot, k_pivot, dlnk


  !------------------------------------------------


  !Read initializing params from file (num_inflaton)
	open(newunit=u, file="parameters_multimodecode.txt", &
    status="old", delim = "apostrophe")
  read(unit=u, nml=init)


  call allocate_vars()

  !Read other params from file
  !	read(unit=u, nml=ic_sampling)
	read(unit=u, nml=params)
	close(unit=u)

  call output_initial_data()

  call calculate_pk_observables(k_pivot,dlnk,As,ns,r,nt)

  contains

    subroutine calculate_pk_observables(k_pivot,dlnk,As,ns,r,nt)

      real(dp), intent(in) :: k_pivot,dlnk
      real(dp), intent(out) :: As,ns,r,nt
      real(dp) :: epsilon, eta
      real(dp) :: ps0, pt0, ps1, pt1, ps2, pt2, x1, x2
      real(dp) :: ps0_iso,ps1_iso,ps2_iso
      real(dp) :: pz0, pz1, pz2

      !Initialize potential and calc background
      call potinit

      !DEBUG
      call DEBUG_writing_etc()

      call evolve(k_pivot, ps0, pt0, pz0,ps0_iso)
      call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1, ps1_iso)
      call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2, ps2_iso)

      epsilon = getEps(phi_pivot, dphi_pivot)
      eta = geteta(phi_pivot, dphi_pivot)

      As = ps0
      ns = 1.d0+log(ps2/ps1)/dlnk/2.d0
      r=pt0/ps0
      nt=log(pt2/pt1)/dlnk/2.d0

      call output_observables((/ps0,ps1,ps2/),(/pt0,pt1,pt2/), &
        (/pz0,pz1,pz2/),(/ps0_iso,ps1_iso,ps2_iso/), &
        ns,r,nt, epsilon,eta)

    end subroutine calculate_pk_observables


    subroutine allocate_vars()

      !Model dependent
      allocate(vparams(vparam_arrays,num_inflaton))

      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(phi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))

    end subroutine allocate_vars

    subroutine output_observables(As,At,Az,A_iso,ns,r,nt, epsilon,eta)

      real(dp), intent(in) :: ns,r,nt, epsilon,eta
      real(dp), dimension(:), intent(in) :: As, At, Az, A_iso

      write(*, i_fmt) "Number of Inflaton =", num_inflaton
      write(*, i_fmt) "Potential Choice =", potential_choice
      write(*, e_fmt) "log10(m^2) =", vparams(1,1)
      write(*, e_fmt) "N_pivot =", N_pivot
      write(*, e2_fmt) "phi_pivot =", phi_pivot(1), '(', sqrt(4*N_pivot + phi_infl_end(1)**2), ')'
      ! [JF] The commented out option below just includes the ind of inflation field coordinates which are negliable in the SR.
      write(*, e2_fmt) "N_tot =", N_tot,'(', 0.25*dot_product(phi_init, phi_init), ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure) 
      write(*, e2_fmt) "Ps =", As(1), '(', N_pivot*H_pivot**2/(4*PI**2), ')'
      write(*, *), As(1), Az(1)
      write(*, *), As(2), Az(2)
      write(*, *), As(3), Az(3)
      write(*, e2_fmt), "Isocurvature P =", A_iso(1)
      write(*, e2_fmt), "Isocurvature P =", A_iso(2)
      write(*, e2_fmt), "Isocurvature P =", A_iso(3)
      write(*, e2_fmt) "Pt/Ps =", r, '(', 16*epsilon, ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, e2_fmt) "n_s =", ns, '(', 1-2*epsilon-1/(N_pivot),')'
      write(*, e2_fmt) "n_t =", nt, '(', -2*epsilon, ')'

    end subroutine output_observables

    subroutine output_initial_data()
      write(ci, '(I2)'), num_inflaton
      ci = adjustl(ci)
      array_fmt = '(a25,'//trim(ci)//'es10.3)'
      write(*, *) 'Testing two field with V(phi) = 1/2 m_I^2 phi_I^2+1/2 m_J^2 phi_J^2'
      write(*, *), "vparams(1,:) =", vparams(1,:)
    end subroutine output_initial_data


    subroutine DEBUG_writing_etc()
      PRINT*, "Writing background solution to phiarr.txt"
      open(1, file='phiarr.txt')
      i = 1
      do while ( lna(i+1) > 1e-8 )
      !PRINT*, lna(i), phiarr(:,i)
      write(1, *), lna(i), phiarr(:,i)
      i = i + 1
      end do
      close(1)
      PRINT*, "Done writing"
      ![DEBUG] [JF]
      PRINT*, "Writing powerspectrum solution to pow.txt"
      open (unit = 20, file = "pow.txt", status = 'replace')
      PRINT*, "Writing field correlation solution to powmatrix.txt"
      open (unit = 3, file = "powmatrix.txt", status = 'replace')
    end subroutine DEBUG_writing_etc



end program test_mmodpk
