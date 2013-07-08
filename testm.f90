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
  real(dp) :: pz0, pz1, pz2
  real(dp) :: epsilon, eta
  real(dp), PARAMETER :: As = 2.d-9
  real(dp), ALLOCATABLE :: lambda(:)

  !! code initialization
  num_inflaton = 2

  allocate(vparams(1,num_inflaton)) !Origionally this line was allocate(vparams(1,1))
  allocate(phi_init0(num_inflaton))
  allocate(phi_init(num_inflaton))
  allocate(phidot_sign(num_inflaton))
  allocate(phiarr(num_inflaton, 1:nsteps))
  allocate(dphiarr(num_inflaton, 1:nsteps))
  allocate(phi_infl_end(num_inflaton))
  allocate(phi_pivot(num_inflaton))
  allocate(dphi_pivot(num_inflaton))

  write(ci, '(I2)'), num_inflaton
  ci = adjustl(ci)
  array_fmt = '(a25,'//trim(ci)//'es10.3)'
  !! end code initialization

  write(*, *) 'Testing two field with V(phi) = 1/2 m_I^2 phi_I^2+1/2 m_J^2 phi_J^2'

  ! \sum_i m_i^2 \phi_i^2

  modpkoutput = .true.
  slowroll_infl_end = .true.
  instreheat = .false.
  N_pivot = 55
  k_pivot = 0.002

  potential_choice = 1
  phi_init0 =(/ 10.31001, 12.93651/)!(/ 11.31001, 11.31001/) !(/ 10., 10., 10., 10./) !(/ 15., 15./) !Origionally this was 20*M_Pl
  !delsigma = 18*M_Pl
  ! phi_infl_end = phi_init0 - delsigma
  vparams(1,:) = (/ -10.422895047377294, -10.422895047377294 + log10(81.0) /) !(/ -12., -12., -12., -12. /)  ! (/ -12., -12. /) !Origionally this was log(96*PI**2*As/(phi_infl_end(1)**2+4*N_pivot)**2)/log(10.d0)
  !vparams(1,:) = (/ -12, -12 /)

  write(*, *), "vparams(1,:) =", vparams(1,:)

  call potinit

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
       open (unit = 2, file = "pow.txt", status = 'replace')
       PRINT*, "Writing field correlation solution to powmatrix.txt"
       open (unit = 3, file = "powmatrix.txt", status = 'replace')
  call evolve(k_pivot, ps0, pt0, pz0)
    ![DEBUG] [JF]
     close(2)
     close(3)
     PRINT*, "finished writing powerspec ralated things"
  dlnk = 0.05
  call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1)
  call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2)

  epsilon = getEps(phi_pivot, dphi_pivot)
  eta = geteta(phi_pivot, dphi_pivot)

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
  write(*, e2_fmt) "Pt/Ps =", pt0/ps0, '(', 16*epsilon, ')'

  write(*, e2_fmt) "n_s =", 1.d0+log(ps2/ps1)/dlnk/2.d0, '(', 1-2*epsilon-eta,')'
  write(*, e2_fmt) "n_s =", 1.d0+log(ps2/ps1)/dlnk/2.d0, '(', 1-2*epsilon-1/(N_pivot),')' ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
  write(*, e2_fmt) "n_t =", log(pt2/pt1)/dlnk/2.d0, '(', -2*epsilon, ')'
  

  !return


!!$  ! \exp(\sum_i \lambda_i \phi_i)
!!$  print*
!!$  write(*, *) 'Testing single field with V(phi) = V_0 exp(lambda phi)'
!!$  
!!$  modpkoutput = .true.
!!$  slowroll_infl_end = .false.
!!$  instreheat = .false.
!!$  N_pivot = 50
!!$  k_pivot = 0.001
!!$
!!$  potential_choice = 7
!!$  deallocate(vparams)
!!$  allocate(vparams(2, num_inflaton))
!!$  allocate(lambda(num_inflaton))
!!$
!!$  phi_init0 = 0
!!$  delsigma = 15*M_Pl 
!!$  vparams(2,:) = (/ -sqrt(0.01), -sqrt(0.02), -sqrt(0.02), -sqrt(0.01) /)
!!$  lambda = vparams(2,:)  
!!$
!!$  phi_infl_end = phi_init0 - delsigma*lambda/sqrt(dot_product(lambda,lambda))
!!$
!!$  vparams(1,1) = 12*PI**2*exp(-dot_product(lambda, (phi_infl_end+lambda*N_pivot)))*As*dot_product(lambda,lambda)
!!$
!!$  write(*, i_fmt) "Number of Inflaton =", num_inflaton
!!$  write(*, i_fmt) "Potential Choice =", potential_choice
!!$
!!$  call potinit
!!$  call evolve(k_pivot, ps0, pt0, pz0)
!!$  dlnk = 0.1
!!$  call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1)
!!$  call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2)
!!$
!!$  epsilon = getEps(phi_pivot, dphi_pivot)
!!$  eta = geteta(phi_pivot, dphi_pivot)
!!$  delsigma = sqrt(dot_product(phiarr(:,1)-phiarr(:,nactual_bg), phiarr(:,1)-phiarr(:,nactual_bg)))
!!$
!!$  write(*, e_fmt) "N_pivot =", N_pivot
!!$  write(*, e2_fmt) "N_tot =", lna(nactual_bg) - lna(1), '(', &
!!$       delsigma/sqrt(dot_product(lambda, lambda)), ')'
!!$  write(*, e2_fmt) "H_pivot =", H_pivot, '(', sqrt(vparams(1,1)/3)*exp(0.5*dot_product(lambda, phi_pivot)), ')'
!!$  write(*, array_fmt) "phi_pivot =", phi_pivot
!!$  write(*, array_fmt) "(phi_pivot) =", phi_infl_end + lambda*N_pivot
!!$  write(*, array_fmt) "dphi_pivot =", dphi_pivot
!!$  write(*, array_fmt) "(dphi_pivot) =", -lambda 
!!$  write(*, e2_fmt) "Ps =", ps0, '(', H_pivot**2 /(8*PI**2*0.5*dot_product(lambda, lambda)), ')'
!!$  write(*, e2_fmt) "Pt/Ps =", pt0/ps0, '(', 16*epsilon, ')'
!!$
!!$  write(*, *), ps0, pz0
!!$  write(*, *), ps1, pz1
!!$  write(*, *), ps2, pz2  
!!$
!!$  write(*, e2_fmt) "n_s =", 1.d0+log(ps2/ps1)/dlnk/2.d0, '(', 1-2*epsilon-eta,')'
!!$  write(*, e2_fmt) "n_t =", log(pt2/pt1)/dlnk/2.d0, '(', -2*epsilon, ')'
!!$
!!$

!!$  open(1, file='result.txt')
!!$  do i=1, nactual_mode
!!$     write(1, '(f10.1, 12e18.6)'), xp(i), yp(:,i)*sqrt(2*k_pivot*Mpc2Mpl)
!!$  end do
!!$  close(1)


end program test_mmodpk
