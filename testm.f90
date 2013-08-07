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

  CHARACTER(16) :: e_fmt = '(a25, es12.4)'
  CHARACTER(36) :: e2_fmt = '(a25, es12.4, a3, es11.4, a1)'
  CHARACTER(16) :: i_fmt = '(a25,I3)'
  CHARACTER(16) :: array_fmt
  CHARACTER(len=2) :: ci
  integer :: i

  real(dp) :: ps0, pt0, ps1, pt1, ps2, pt2, dlnk, x1, x2
  real(dp) :: ps0_iso,ps1_iso,ps2_iso
  real(dp) :: pz0, pz1, pz2
  real(dp) :: epsilon, eta

  integer :: u

  !For run-time alloc w/out re-compile
  namelist /init/ num_inflaton, potential_choice, &
    modpkoutput, slowroll_infl_end, instreheat

  namelist /params/ phi_init0, vparams, &
    N_pivot, k_pivot, dlnk



  !Read initializing params from file (num_inflaton)
	open(newunit=u, file="parameters_multimodecode.txt", &
    status="old", delim = "apostrophe")
  read(unit=u, nml=init)

  call allocate_vars()

  !Read other params from file
	read(unit=u, nml=params)
	close(unit=u)

  call output_initial_data()

  call potinit

  !DEBUG
  call DEBUG_writing_etc()

  call evolve(k_pivot, ps0, pt0, pz0,ps0_iso)
  call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1, ps1_iso)
  call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2, ps2_iso)

  epsilon = getEps(phi_pivot, dphi_pivot)
  eta = geteta(phi_pivot, dphi_pivot)

  call output_observables()

  contains

    subroutine allocate_vars()

      allocate(vparams(1,num_inflaton)) !Origionally this line was allocate(vparams(1,1))
      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(phi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))

    end subroutine allocate_vars

    subroutine output_observables()

      write(*, i_fmt) "Number of Inflaton =", num_inflaton
      write(*, i_fmt) "Potential Choice =", potential_choice
      write(*, e_fmt) "log10(m^2) =", vparams(1,1)
      write(*, e_fmt) "N_pivot =", N_pivot
      !write(*, e2_fmt) "phi_pivot =", phi_pivot(1), '(', sqrt(4*N_pivot + phi_infl_end(1)**2), ')'
      ! [JF] Need to look up this formatting confusion next time you have a book or the internet....
      write(*, *),  "phi_pivot =", phi_pivot ! [JF] This line should be replaced with the previous line.
      ! [JF] The commented out option below just includes the ind of inflation field coordinates which are negliable in the SR.
      write(*, e2_fmt) "N_tot =", N_tot,'(', 0.25*dot_product(phi_init, phi_init), ')' !'(', 0.25*(dot_product(phi_init, phi_init) - dot_product(phi_infl_end, phi_infl_end)), ')'

      !write(*, e2_fmt) "Ps =", ps0, '(', H_pivot**2/(8*PI**2*epsilon), ')'
      write(*, e2_fmt) "Ps =", ps0, '(', N_pivot*H_pivot**2/(4*PI**2), ')' ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, *), ps0, pz0
      write(*, *), ps1, pz1
      write(*, *), ps2, pz2
      write(*, *), "Isocurvature P =", ps0_iso, ps1_iso, ps2_iso
      write(*, e2_fmt) "Pt/Ps =", pt0/ps0, '(', 16*epsilon, ')'

      write(*, e2_fmt) "n_s =", 1.d0+log(ps2/ps1)/dlnk/2.d0, '(', 1-2*epsilon-eta,')'
      write(*, e2_fmt) "n_s =", 1.d0+log(ps2/ps1)/dlnk/2.d0, '(', 1-2*epsilon-1/(N_pivot),')' ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, e2_fmt) "n_t =", log(pt2/pt1)/dlnk/2.d0, '(', -2*epsilon, ')'

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
