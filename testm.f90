program test_mmodpk
  use modpkparams
  use potential
  use background_evolution
  use modpk_utils
  use camb_interface
  use ode_path
  use access_modpk
  use internals
  use modpk_icsampling
  use modpk_rng, only : init_random_seed

#ifdef MPI
  use mpi
#endif

  implicit none

  !Writing fmts
  character(16), parameter:: e_fmt = '(a25, es12.4)'
  character(36), parameter:: e2_fmt = '(a25, es12.4, a3, es11.4, a1)'
  character(16), parameter:: i_fmt = '(a25,I3)'
  character(16) :: array_fmt
  character(len=2) :: ci

  !Run-specific input params
  integer :: i, vparam_rows

  !Parallel variables
  integer :: numtasks, rank

  !Cosmology
  real(dp) :: dlnk, As, ns, nt, r

  !Sampling parameters for ICs
  integer :: numb_samples
  real(dp) :: energy_scale
  real(dp), dimension(:,:), allocatable :: priors_min, priors_max

  type(ic_and_observables), dimension(:), allocatable :: ic_output
  type(ic_and_observables), dimension(:), allocatable :: ic_output_iso_N

  integer :: u

  !For run-time alloc w/out re-compile
  namelist /init/ num_inflaton, potential_choice, &
    modpkoutput, slowroll_infl_end, instreheat, vparam_rows

  namelist /ic_sampling/ sampling_techn, energy_scale, numb_samples, &
    save_iso_N, N_iso_ref

  namelist /params/ phi_init0, vparams, &
    N_pivot, k_pivot, dlnk


  !------------------------------------------------


  !Read initializing params from file
	open(newunit=u, file="parameters_multimodecode.txt", &
    status="old", delim = "apostrophe")
  read(unit=u, nml=init)

  call allocate_vars()

  !Read other params from file
	read(unit=u, nml=ic_sampling)
	read(unit=u, nml=params)
	close(unit=u)

  call output_initial_data()

  if (sampling_techn==reg_samp) then

    call calculate_pk_observables(k_pivot,dlnk,As,ns,r,nt)

  else if (sampling_techn == eqen_samp .or. &
    sampling_techn == slowroll_samp) then

#ifdef MPI
    call mpi_parallelize()
	  !Set random seed for each process.
	  call init_random_seed(rank)
#else
	  call init_random_seed()
#endif

    !Initialize the sampler
    call init_sampler(priors_min, priors_max)

    do i=1,numb_samples

      call calculate_pk_observables_per_IC(k_pivot,dlnk,As,ns,r,nt)

      !Load & print output array
      !Save in ic_output in case want to post-process.
      !Comment-out if don't want to keep and just do write(1,*)
      call ic_output(i)%load_observables(phi_init0, dphi_init0,As,ns,r,nt)
      call ic_output(i)%printout(1)
      if (save_iso_N) then
        call ic_output_iso_N(i)%load_observables(phi_iso_N, dphi_iso_N, &
          As,ns,r,nt)
        call ic_output_iso_N(i)%printout(2)
      end if

    end do

  else
    print*, "ERROR: sampling technique",sampling_techn,"not implemented."
    stop
  end if

  contains

    subroutine calculate_pk_observables(k_pivot,dlnk,As,ns,r,nt)

      real(dp), intent(in) :: k_pivot,dlnk
      real(dp), intent(out) :: As,ns,r,nt
      real(dp) :: epsilon, eta
      real(dp) :: ps0, pt0, ps1, pt1, ps2, pt2, x1, x2
      real(dp) :: ps0_iso,ps1_iso,ps2_iso
      real(dp) :: pz0, pz1, pz2
      real(dp), dimension(:,:), allocatable :: pk_arr, pk_iso_arr
      logical :: calc_full_pk

      !Initialize potential and calc background
      call potinit

      !DEBUG
      call DEBUG_writing_etc()

      call evolve(k_pivot, ps0, pt0, pz0,ps0_iso)
      call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1, ps1_iso)
      call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2, ps2_iso)

      !Get full spectrum for adiab and isocurv at equal intvs in lnk
      call get_full_pk(pk_arr,pk_iso_arr,dlnk,calc_full_pk)

      epsilon = getEps(phi_pivot, dphi_pivot)
      eta = geteta(phi_pivot, dphi_pivot)

      As = ps0
      ns = 1.d0+log(ps2/ps1)/dlnk/2.d0
      r=pt0/ps0
      nt=log(pt2/pt1)/dlnk/2.d0

      call output_observables(pk_arr,pk_iso_arr, &
        (/ps0,ps1,ps2/),(/pt0,pt1,pt2/), &
        (/pz0,pz1,pz2/),(/ps0_iso,ps1_iso,ps2_iso/), &
        ns,r,nt, epsilon,eta,calc_full_pk)


    end subroutine calculate_pk_observables


    subroutine get_full_pk(pk_arr,pk_iso_arr,dlnk,calc_full_pk)

      real(dp), dimension(:,:), allocatable, intent(out) :: pk_arr
      real(dp), allocatable, optional, intent(out) :: pk_iso_arr(:,:)
      real(dp), intent(in) :: dlnk

      real(dp) :: kmin, kmax, incr
      logical, intent(inout) :: calc_full_pk
      real(dp) :: p_scalar, p_tensor, p_zeta, p_iso, k_input
      integer :: i, steps, u

      namelist /full_pk/ kmin, kmax, steps, calc_full_pk

	    open(newunit=u, file="parameters_multimodecode.txt", &
        status="old", delim = "apostrophe")
      read(unit=u, nml=full_pk)
      close(u)

      !If don't want full spectrum, return
      if (.not. calc_full_pk) return

      !Make the output arrays
      if (allocated(pk_arr)) deallocate(pk_arr)
      if (allocated(pk_iso_arr) .and. present(pk_iso_arr)) &
        deallocate(pk_iso_arr)
      allocate(pk_arr(steps, 2))
      if (present(pk_iso_arr)) allocate(pk_iso_arr(steps, 2))
      pk_arr=0e0_dp
      if (present(pk_iso_arr)) pk_iso_arr=0e0_dp

      k_input=kmin
      incr=(kmax/kmin)**(1/real(steps-1))
      do i=1,steps
        k_input = kmin*incr**(i-1)
        call evolve(k_input, p_scalar, p_tensor, p_zeta, p_iso)

        pk_arr(i,:)=(/k_input,p_scalar/)
        if (present(pk_iso_arr)) pk_iso_arr(i,:)=(/k_input,p_iso/)

      end do


    end subroutine get_full_pk


    subroutine allocate_vars()

      !Model dependent
      allocate(vparams(vparam_rows,num_inflaton))
      if (sampling_techn/=reg_samp) then
        allocate(priors_max(2,num_inflaton))
        allocate(priors_min(2,num_inflaton))
        allocate(dphi_init0(num_inflaton))
        allocate(dphi_init(num_inflaton))
      end if

      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(phi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))


    end subroutine allocate_vars

    subroutine output_observables(pk_arr, pk_iso_arr,&
        As,At,Az,A_iso,ns,r,nt, epsilon,eta,calc_full_pk)

      real(dp), dimension(:,:), intent(in) :: pk_arr
      real(dp), dimension(:,:), intent(in), optional :: pk_iso_arr
      real(dp), intent(in) :: ns,r,nt, epsilon,eta
      real(dp), dimension(:), intent(in) :: As, At, Az, A_iso

      logical :: calc_full_pk

      integer :: i

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
      write(*, e2_fmt), "Bundle Width =", field_bundle%width
      write(*, e2_fmt) "Pt/Ps =", r, '(', 16*epsilon, ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, e2_fmt) "n_s =", ns, '(', 1-2*epsilon-1/(N_pivot),')'
      write(*, e2_fmt) "n_t =", nt, '(', -2*epsilon, ')'

      if (calc_full_pk) then
        do i=1,size(pk_arr,1)
          write(101,*) pk_arr(i,:)
          if(present(pk_iso_arr)) write(102,*) pk_iso_arr(i,:)
        end do
        write(*,*) "Adiab P(k) written to file #101"
        if (present(pk_iso_arr)) write(*,*) "Iso-P(k) written to file #102"
      end if

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

    !Calculate observables, but grab a new IC each time called
    subroutine calculate_pk_observables_per_IC(k_pivot,dlnk,As,ns,r,nt)

      real(dp), intent(in) :: k_pivot,dlnk
      real(dp), intent(out) :: As,ns,r,nt
      real(dp) :: epsilon, eta
      real(dp) :: ps0, pt0, ps1, pt1, ps2, pt2, x1, x2
      real(dp) :: ps0_iso,ps1_iso,ps2_iso
      real(dp) :: pz0, pz1, pz2
      real(dp), dimension(:,:), allocatable :: pk_arr, pk_iso_arr
      logical :: calc_full_pk

      call get_ic(phi_init0, dphi_init0, sampling_techn, &
        priors_min, priors_max, &
         numb_samples,energy_scale)

      !Initialize potential and calc background
      call potinit

      !DEBUG
      !call DEBUG_writing_etc()

      call evolve(k_pivot, ps0, pt0, pz0,ps0_iso)
      call evolve(k_pivot*exp(-dlnk), ps1, pt1, pz1, ps1_iso)
      call evolve(k_pivot*exp(dlnk), ps2, pt2, pz2, ps2_iso)

      !Get full spectrum for adiab and isocurv at equal intvs in lnk
      call get_full_pk(pk_arr,pk_iso_arr,dlnk,calc_full_pk)

      epsilon = getEps(phi_pivot, dphi_pivot)
      eta = geteta(phi_pivot, dphi_pivot)

      As = ps0
      ns = 1.d0+log(ps2/ps1)/dlnk/2.d0
      r=pt0/ps0
      nt=log(pt2/pt1)/dlnk/2.d0

      call output_observables(pk_arr,pk_iso_arr, &
        (/ps0,ps1,ps2/),(/pt0,pt1,pt2/), &
        (/pz0,pz1,pz2/),(/ps0_iso,ps1_iso,ps2_iso/), &
        ns,r,nt, epsilon,eta,calc_full_pk)


    end subroutine calculate_pk_observables_per_IC

    subroutine init_sampler(priors_min, priors_max)


      real(dp), dimension(:,:), intent(out) :: priors_min, &
        priors_max

      real(dp), dimension(:), allocatable :: phi0_priors_min, &
        dphi0_priors_min, phi0_priors_max, dphi0_priors_max

      integer :: u, i

      namelist /priors/ phi0_priors_min, phi0_priors_max, &
        dphi0_priors_min, dphi0_priors_max, &
        penalty_fact

      if (allocated(phi0_priors_max)) then
        print*, "ERROR: Priors allocated before initialization."
        stop
      else
        allocate(phi0_priors_max(num_inflaton))
        allocate(dphi0_priors_max(num_inflaton))
        allocate(phi0_priors_min(num_inflaton))
        allocate(dphi0_priors_min(num_inflaton))
        phi0_priors_min=0e0_dp
        phi0_priors_max=0e0_dp
        dphi0_priors_min=0e0_dp
        dphi0_priors_max=0e0_dp
      end if

      if (save_iso_N) then
        allocate(phi_iso_N(num_inflaton))
        allocate(dphi_iso_N(num_inflaton))
      end if

      !Read phi0 priors from file
	    open(newunit=u, file="parameters_multimodecode.txt", &
        status="old", delim = "apostrophe")
      read(unit=u, nml=priors)
      close(u)

      !Make ouput array(s)
      allocate(ic_output(numb_samples))
      if (save_iso_N) allocate(ic_output_iso_N(numb_samples))
      do i=1,size(ic_output)
        allocate(ic_output(i)%ic(2*num_inflaton))
        if (save_iso_N) allocate(ic_output_iso_N(i)%ic(2*num_inflaton))
      end do

      priors_max(1,:) = phi0_priors_max
      priors_max(2,:) = dphi0_priors_max
      priors_min(1,:) = phi0_priors_min
      priors_min(2,:) = dphi0_priors_min

    end subroutine init_sampler

!NOT WORKING YET
#ifdef MPI
    subroutine mpi_parallelize()

      integer :: ierr, rc

	    call mpi_init(ierr)
	    	if(ierr .ne. 0) then
	    		print*,"Error parallelizing."
	    		call mpi_abort(mpi_comm_world, rc, ierr)
	    		stop
	    	end if

	    !Obtains info on processors.
	    call mpi_comm_rank(mpi_comm_world, rank, ierr)
	    call mpi_comm_size(mpi_comm_world, numtasks, ierr)
	    print*,'Number of tasks=',numtasks,' My rank=',rank
	    call mpi_barrier(mpi_comm_world,ierr)


    end subroutine mpi_parallelize
#endif


end program test_mmodpk
