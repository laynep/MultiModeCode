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
  character(36), parameter:: e2_fmt = '(a25, es17.9, a3, es11.4, a1)'
  character(16), parameter:: i_fmt = '(a25,I3)'
  character(16) :: array_fmt
  character(len=2) :: ci

  !Run-specific input params
  integer :: i, vparam_rows

  !Parallel variables
  integer :: numtasks, rank

  !Cosmology
  real(dp) :: dlnk, As, ns, nt, r, alpha_s
  real(dp) :: A_iso, A_pnad, A_ent, A_bundle
  real(dp) :: n_iso, n_pnad, n_ent


  !Sampling parameters for ICs
  integer :: numb_samples
  integer :: outsamp, outsamp_N_iso
  real(dp) :: energy_scale
  real(dp), dimension(:,:), allocatable :: priors_min, priors_max
  logical :: output_badic

  type(ic_and_observables), dimension(:), allocatable :: ic_output
  type(ic_and_observables), dimension(:), allocatable :: ic_output_iso_N

  integer :: u
	character(len=15) :: isoNname, sampname

  !For run-time alloc w/out re-compile
  namelist /init/ num_inflaton, potential_choice, &
    modpkoutput, slowroll_infl_end, instreheat, vparam_rows

  namelist /ic_sampling/ sampling_techn, energy_scale, numb_samples, &
    save_iso_N, N_iso_ref,output_badic

  namelist /params/ phi_init0, dphi_init0, vparams, &
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

    call calculate_pk_observables(k_pivot,dlnk)

  !Eqen sampling
  else if (sampling_techn == eqen_samp .or. &
    !Set vels in SR and fields on iso-N surface for N-quad
    sampling_techn == iso_N .or.&
    !Set vels in SR
    sampling_techn == slowroll_samp .or.&
    !Grab IC from file
    sampling_techn == fromfile_samp .or. &
    !Loop over different vparams for given num_inflaton
    sampling_techn == parameter_loop_samp) then

#ifdef MPI

    call mpi_parallelize()
	  !Set random seed for each process.
	  call init_random_seed(rank)

    call open_output_files()

#else

    !Open output files
    open(newunit=outsamp,file="ic_eqen.txt")
    open(newunit=outsamp_N_iso,file="ic_isoN.txt")

	  !Set random seed
	  call init_random_seed()

#endif

    !Initialize the sampler
    call init_sampler(priors_min, priors_max)

    do i=1,numb_samples

        if (modpkoutput) write(*,*) "---------------------------------------------"
        if (modpkoutput) write(*,*) "Sample numb", i, "of", numb_samples
        if (modpkoutput) write(*,*) "---------------------------------------------"

      call calculate_pk_observables(k_pivot,dlnk)

    end do

  else
    print*, "ERROR: sampling technique",sampling_techn,"not implemented."
    stop
  end if

#ifdef MPI
	!Halts processors here.
	call mpi_barrier(mpi_comm_world,i)
#endif

  contains

    subroutine open_output_files()

      !Open output files for each rank
      outsamp=300+rank
      outsamp_N_iso=300+rank+numtasks

		  write(sampname,'(a,i4.4,a)')'ic_eqen',outsamp,'.txt'
		  write(isoNname,'(a,i4.4,a)')'ic_isoN',outsamp_N_iso,'.txt'

		  open(unit=outsamp,status='new',file=sampname)
      if (save_iso_N) open(unit=outsamp_N_iso,status='new',file=isoNname)

    end subroutine open_output_files


    subroutine get_full_pk(pk_arr,pk_iso_arr,calc_full_pk)

      real(dp), dimension(:,:), allocatable, intent(out) :: pk_arr
      real(dp), allocatable, optional, intent(out) :: pk_iso_arr(:,:)

      real(dp) :: kmin, kmax, incr
      logical, intent(inout) :: calc_full_pk
      real(dp) :: p_scalar, p_tensor, p_zeta, p_iso, k_input
      integer :: i, steps, u

      type(power_spectra) :: pk

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
        call evolve(k_input, pk)

        pk_arr(i,:)=(/k_input,pk%adiab/)
        if (present(pk_iso_arr)) pk_iso_arr(i,:)=(/k_input,pk%isocurv/)

      end do


    end subroutine get_full_pk


    subroutine allocate_vars()

      !Model dependent
      allocate(vparams(vparam_rows,num_inflaton))
      allocate(priors_max(2,num_inflaton))
      allocate(priors_min(2,num_inflaton))

      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(phi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))
      allocate(dphi_init0(num_inflaton))
      allocate(dphi_init(num_inflaton))


    end subroutine allocate_vars

    subroutine output_observables(pk_arr, pk_iso_arr,&
        As,At,Az,A_iso,A_pnad,A_ent,A_cross,ns,r,nt, alpha_s,&
        eps,eta,calc_full_pk)

      real(dp), dimension(:,:), intent(in) :: pk_arr
      real(dp), dimension(:,:), intent(in), optional :: pk_iso_arr
      real(dp), intent(in) :: r,nt, eps,eta, alpha_s
      real(dp), dimension(:), intent(in) :: As, At, Az, A_iso, &
        A_pnad,A_ent, ns
      real(dp), intent(in) :: A_cross

      logical :: calc_full_pk

      integer :: i

      real(dp) :: phi_piv_pred, N_tot_pred, Ps_pred, r_pred, ns_pred, nt_pred,&
      alphas_pred

      !Predictions for N-quad
      phi_piv_pred = sqrt(4*N_pivot + phi_infl_end(1)**2)
      N_tot_pred = 0.25*dot_product(phi_init, phi_init)
      Ps_pred = N_pivot*H_pivot**2/(4*PI**2)
      r_pred = 16*eps
      ns_pred = 1-2*eps-1/(N_pivot)
      nt_pred = -2*eps
      alphas_pred = 8.0*eps*(2.0*eta - 3.0*eps)

      write(*, i_fmt) "Number of Inflaton =", num_inflaton
      write(*, i_fmt) "Potential Choice =", potential_choice
      !write(*, e_fmt) "log10(m^2) =", vparams(1,1)
      write(*, e_fmt) "N_pivot =", N_pivot
      write(*, e2_fmt) "phi_pivot =", phi_pivot(1), '(', phi_piv_pred , ')'
      ! [JF] The commented out option below just includes the ind of inflation field coordinates which are negliable in the SR.
      write(*, e2_fmt) "N_tot =", N_tot,'(', N_tot_pred , ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure) 
      write(*, e2_fmt) "Ps =", As(1), '(', Ps_pred , ')'
      !write(*, *), As(1), Az(1)
      !write(*, *), As(2), Az(2)
      !write(*, *), As(3), Az(3)
      write(*, e2_fmt), "Isocurvature P =", A_iso(1)
      write(*, e2_fmt), "Pnad P =", A_pnad(1)
      write(*, e2_fmt), "Entropy P =", A_ent(1)
      write(*, e2_fmt), "(TEST) Cross Ad-Iso P =", A_cross
      write(*, e2_fmt), "Bundle Expand Scalar =", field_bundle%exp_scalar
      write(*, e2_fmt) "r = Pt/Ps =", r, '(', r_pred, ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, e2_fmt) "n_s =", ns(1), '(', ns_pred ,')'
      if (size(ns)>1) then
        write(*, e2_fmt) "n_iso =", ns(2)
        write(*, e2_fmt) "n_pnad =", ns(3)
        write(*, e2_fmt) "n_ent =", ns(4)
      end if
      write(*, e2_fmt) "n_t =", nt, '(', nt_pred , ')'
      write(*, e2_fmt) "alpha_s =", alpha_s, '(', alphas_pred , ')'


      if (calc_full_pk) then
        do i=1,size(pk_arr,1)
          write(101,*) pk_arr(i,:)
          if(present(pk_iso_arr)) write(102,*) pk_iso_arr(i,:)
        end do
        write(*,*) "Adiab P(k) written to file #101"
        if (present(pk_iso_arr)) write(*,*) "Iso-P(k) written to file #102"
      end if

!DEBUG
!if (ns(1)<0.92 .or. ns(1)>0.97) then
!if (A_iso(1)>1e-12) then
if (abs(ns(1)-0.963)<1e-5 .and. abs(alpha_s)>5e-3) then
  print*, 'Outlier'
  print*, ns(1), alpha_s
  stop
end if

    end subroutine output_observables

    subroutine output_initial_data()
      integer :: i

      write(ci, '(I2)'), num_inflaton
      ci = adjustl(ci)
      array_fmt = '(a25,'//trim(ci)//'es10.3)'
      !write(*, *) 'Testing two field with V(phi) = 1/2 m_I^2 phi_I^2+1/2 m_J^2 phi_J^2'
      if (size(vparams,1)>1) then
        do i=1, size(vparams,1)
          write(*, '(A8,I1,A5,100E12.3)'), "vparams(",i,",:) =", vparams(i,:)
        end do
      else
        write(*, *), "vparams(1,:) =", vparams(1,:)
      end if
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
      PRINT*, "Writing powerspectrum solution to pow.txt"
      open (unit = 20, file = "pow.txt", status = 'replace')
      PRINT*, "Writing field correlation solution to powmatrix.txt"
      open (unit = 3, file = "powmatrix.txt", status = 'replace')
    end subroutine DEBUG_writing_etc

    !Calculate observables, optionally grab a new IC each time called
    subroutine calculate_pk_observables(k_pivot,dlnk)

      real(dp), intent(in) :: k_pivot,dlnk
      real(dp) :: As,ns,r,nt, alpha_s
      real(dp) :: runofrun
      real(dp) :: A_iso, A_pnad, A_ent, A_bundle, A_cross
      real(dp) :: n_iso, n_pnad, n_ent
      real(dp) :: epsilon, eta
      real(dp) :: ps0, pt0, ps1, pt1, ps2, pt2, x1, x2
      real(dp) :: ps0_iso,ps1_iso,ps2_iso
      real(dp) :: pz0, pz1, pz2
      real(dp) :: pnad0, pnad1, pnad2
      real(dp) :: pent0, pent1, pent2
      real(dp), dimension(:,:), allocatable :: pk_arr, pk_iso_arr
      logical :: calc_full_pk, leave

      type(power_spectra) :: pk0, pk1, pk2, pk3, pk4

      pk_bad=0
      leave = .false.

      if (sampling_techn/=reg_samp) then
        call get_ic(phi_init0, dphi_init0, sampling_techn, &
          priors_min, priors_max, &
          numb_samples,energy_scale)
      end if

      !Initialize potential and calc background
      call potinit

      call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
        A_iso, A_pnad, A_ent, A_bundle, &
        n_iso, n_pnad, n_ent, &
        leave)
      if (leave) return

      call evolve(k_pivot, pk0)
        call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
          A_iso, A_pnad, A_ent, A_bundle, &
          n_iso, n_pnad, n_ent, &
          leave)
        if (leave) return
!DEBUG
!print*, "Not evaluating second and third evolve routines"
      call evolve(k_pivot*exp(-dlnk), pk1)
        call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
          A_iso, A_pnad, A_ent, A_bundle, &
          n_iso, n_pnad, n_ent, &
          leave)
        if (leave) return
      call evolve(k_pivot*exp(dlnk), pk2)
        call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
          A_iso, A_pnad, A_ent, A_bundle, &
          n_iso, n_pnad, n_ent, &
          leave)
        if (leave) return

      !!Uncomment here and below for alpha_s from 5-pt stencil
      !!or running of running
      !call evolve(k_pivot*exp(-2.0e0_dp*dlnk), pk3)
      !  call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
      !    A_iso, A_pnad, A_ent, A_bundle, &
      !    n_iso, n_pnad, n_ent, &
      !    leave)
      !  if (leave) return
      !call evolve(k_pivot*exp(2.0e0_dp*dlnk), pk4)
      !  call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
      !    A_iso, A_pnad, A_ent, A_bundle, &
      !    n_iso, n_pnad, n_ent, &
      !    leave)
      !  if (leave) return

      ps0= pk0%adiab
      ps1= pk1%adiab
      ps2= pk2%adiab
      pt0= pk0%tensor
      pt1= pk1%tensor
      pt2= pk2%tensor
      ps0_iso=  pk0%isocurv
      ps1_iso=  pk1%isocurv
      ps2_iso=  pk2%isocurv
      pz0= pk0%powz
      pz1= pk1%powz
      pz2= pk2%powz
      pnad0=pk0%pnad
      pnad1=pk1%pnad
      pnad2=pk2%pnad
      pent0=pk0%entropy
      pent1=pk1%entropy
      pent2=pk2%entropy

      A_iso=ps0_iso
      A_pnad=pnad0
      A_ent=pent0
      A_cross = pk0%cross_ad_iso

      A_bundle=field_bundle%exp_scalar


      epsilon = getEps(phi_pivot, dphi_pivot)
      eta = geteta(phi_pivot, dphi_pivot)

      As = ps0
      ns = 1.d0+log(ps2/ps1)/dlnk/2.d0
      r=pt0/ps0
      nt=log(pt2/pt1)/dlnk/2.d0

      alpha_s = log(ps2*ps1/ps0**2)/dlnk**2


      !alpha_s from 5-pt stencil
      !alpha_s = (1.0e0_dp/12.0e0_dp/dlnk**2)*&
      !  (-log(pk4%adiab) + 16.0e0_dp*log(pk2%adiab) - &
      !  30.0e0_dp*log(pk0%adiab) + 16.0e0_dp*log(pk1%adiab) - &
      !  log(pk3%adiab))

      !runofrun = (1.0e0_dp/2.0e0_dp/dlnk**3)*&
      !  (log(pk4%adiab) -2* log(pk2%adiab) + 2*log(pk1%adiab) -log(pk3%adiab))

      !print*, "running of running =", runofrun

      n_iso=log(ps2_iso/ps1_iso)/dlnk/2.d0
      n_pnad=log(pnad2/pnad1)/dlnk/2.d0
      n_ent=log(pent2/pent1)/dlnk/2.d0


      !Get full spectrum for adiab and isocurv at equal intvs in lnk
      call get_full_pk(pk_arr,pk_iso_arr,calc_full_pk)

      if (modpkoutput) then
        call output_observables(pk_arr,pk_iso_arr, &
          (/ps0,ps1,ps2/),(/pt0,pt1,pt2/), &
          (/pz0,pz1,pz2/),(/ps0_iso,ps1_iso,ps2_iso/), &
          (/pnad0,pnad1,pnad2/),(/pent0,pent1,pent2/),A_cross,&
          (/ns,n_iso,n_pnad,n_ent/),r,nt,alpha_s, &
          epsilon,eta,calc_full_pk)
      end if

      if (sampling_techn/=reg_samp) then
        !Load & print output array
        !Save in ic_output in case want to post-process.
        !Comment-out if don't want to keep and just do write(1,*)
        call ic_output(i)%load_observables(phi_init0, dphi_init0,As,ns,r,nt,&
        alpha_s, A_iso, A_pnad, A_ent, A_bundle, n_iso, n_pnad, n_ent,&
        A_cross)
        if (output_badic .or. pk_bad/=bad_ic) then
          call ic_output(i)%printout(outsamp)
        endif

        if (save_iso_N) then
          call ic_output_iso_N(i)%load_observables(phi_iso_N, dphi_iso_N, &
            As,ns,r,nt, alpha_s,&
            A_iso, A_pnad, A_ent, A_bundle, n_iso, n_pnad, n_ent,&
            A_cross)
          if (output_badic .or. pk_bad/=bad_ic) then
            call ic_output_iso_N(i)%printout(outsamp_N_iso)
          end if
        end if
      end if


    end subroutine calculate_pk_observables


    subroutine test_bad(pk_bad,As,ns,r,nt,alpha_s,&
      A_iso, A_pnad, A_ent, A_bundle, &
      n_iso, n_pnad, n_ent, &
      leave)

      integer :: pk_bad
      logical :: leave
      real(dp) ::As,ns,r,nt, alpha_s
      real(dp) :: A_iso, A_pnad, A_ent, A_bundle
      real(dp) :: n_iso, n_pnad, n_ent

      !If pk_bad==bad_ic, then restart IC
      !If pk_bad==4, then ode_underflow
      if (pk_bad==bad_ic .or. pk_bad==4) then

        !Set all observs to 0
        As =0e0_dp
        ns=0e0_dp
        r=0e0_dp
        nt=0e0_dp
        alpha_s=0e0_dp
        A_iso=0e0_dp
        A_pnad=0e0_dp
        A_ent=0e0_dp
        A_bundle=0e0_dp
        n_iso=0e0_dp
        n_pnad=0e0_dp
        n_ent=0e0_dp
        leave = .true.

      end if


    end subroutine test_bad

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
