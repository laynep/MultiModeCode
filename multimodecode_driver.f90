program multimodecode
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
  use modpk_output, only : out_opt
  use modpk_deltaN_SR

  implicit none


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
  integer :: out_adiab, out_isoc
  real(dp) :: energy_scale
  real(dp), dimension(:,:), allocatable :: priors_min, priors_max

  !Other sampling params
  real(dp) :: N_pivot_prior_min, N_pivot_prior_max
  logical :: varying_N_pivot
  logical :: more_potential_params

  type(ic_and_observables), dimension(:), allocatable :: ic_output, ic_output_iso_N

  integer :: u

  !For run-time alloc w/out re-compile
  namelist /init/ num_inflaton, potential_choice, &
    slowroll_infl_end, instreheat, vparam_rows, &
    more_potential_params, use_deltaN_SR, evaluate_modes

  namelist /ic_sampling/ sampling_techn, energy_scale, numb_samples, &
    save_iso_N, N_iso_ref, varying_N_pivot

  namelist /params/ phi_init0, dphi_init0, vparams, &
    N_pivot, k_pivot, dlnk

  namelist /more_params/ effective_V_choice, turning_choice, &
    number_knots_qsfrandom, stand_dev_qsfrandom, &
    knot_range_min, knot_range_max, custom_knot_range

  namelist /print_out/ out_opt

  !------------------------------------------------


  !Read initializing params from file
	open(newunit=u, file="parameters_multimodecode.txt", &
    status="old", delim = "apostrophe")
  read(unit=u, nml=init)

  call allocate_vars()

  !Read other params from file
	read(unit=u, nml=ic_sampling)
	read(unit=u, nml=params)
	read(unit=u, nml=print_out)
	close(unit=u)

  call output_initial_data()

  if (sampling_techn==reg_samp) then

    call out_opt%open_files()

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
    sampling_techn == parameter_loop_samp .or. &
    sampling_techn == param_unif_prior .or. &
    sampling_techn == qsf_random  &
    ) then

    call out_opt%open_files(ICs=.true.)


	  !Set random seed
	  call init_random_seed()

    !Initialize the sampler
    call init_sampler(priors_min, priors_max)

    do i=1,numb_samples

      if (out_opt%modpkoutput) write(*,*) "---------------------------------------------"
      if (out_opt%modpkoutput) write(*,*) "Sample numb", i, "of", numb_samples
      if (out_opt%modpkoutput) write(*,*) "---------------------------------------------"

      call calculate_pk_observables(k_pivot,dlnk)

    end do

  else
    print*, "ERROR: sampling technique",sampling_techn,"not implemented."
    stop
  end if

  contains

    subroutine get_full_pk(pk_arr,pk_iso_arr,calc_full_pk)

      real(dp), dimension(:,:), allocatable, intent(out) :: pk_arr
      real(dp), allocatable, optional, intent(out) :: pk_iso_arr(:,:)

      real(dp) :: kmin, kmax, incr
      logical, intent(inout) :: calc_full_pk
      real(dp) :: p_scalar, p_tensor, p_zeta, p_iso
      real(dp), dimension(:), allocatable :: k_input
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

      !Make the arrays for k values to sample
      allocate(k_input(steps))
      !k_input=kmin
      incr=(kmax/kmin)**(1/real(steps-1))
      do i=1,steps
        k_input(i) = kmin*incr**(i-1)
      end do

      do i=1,steps
        call evolve(k_input(i), pk)

        pk_arr(i,:)=(/k_input(i),pk%adiab/)
        if (present(pk_iso_arr)) pk_iso_arr(i,:)=(/k_input(i),pk%isocurv/)

      end do

    end subroutine get_full_pk


    subroutine allocate_vars()

      !Prepare extra params if necessary
      if (more_potential_params) then
        allocate(turning_choice(num_inflaton-1))
        allocate(number_knots_qsfrandom(num_inflaton-1))
        allocate(stand_dev_qsfrandom(num_inflaton-1))
        allocate(knot_range_min(num_inflaton-1))
        allocate(knot_range_max(num_inflaton-1))

        read(unit=u, nml=more_params)

      end if

      !Model dependent
      if (potential_choice==8) then
        allocate(vparams(1,4))
      else
        allocate(vparams(vparam_rows,num_inflaton))
      end if
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
      alphas_pred = 8.0*eps*(2.0*eta - 3.0*eps) !This prediction isn't as good.

      !DEBUG
      print*, "----------------"
      print*, "---FIX OUTPUT---"
      print*, "----------------"
      do i=1,size(vparams,1)
        print*, "vparams", vparams(i,:)
      end do

      write(*, out_opt%i_fmt) "Number of Inflaton =", num_inflaton
      write(*, out_opt%i_fmt) "Potential Choice =", potential_choice
      write(*, out_opt%e_fmt) "N_pivot =", N_pivot
      write(*, out_opt%e2_fmt) "phi_pivot =", phi_pivot(1), '(', phi_piv_pred , ')'
      ! [JF] The commented out option below just includes the ind of inflation field coordinates which are negliable in the SR.
      write(*, out_opt%e2_fmt) "N_tot =", N_tot,'(', N_tot_pred , ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, out_opt%e2_fmt) "Ps =", As(1), '(', Ps_pred , ')'
      write(*, out_opt%e2_fmt), "Isocurvature P =", A_iso(1)
      write(*, out_opt%e2_fmt), "Pnad P =", A_pnad(1)
      write(*, out_opt%e2_fmt), "Entropy P =", A_ent(1)
      write(*, out_opt%e2_fmt), "Cross Ad-Iso P =", A_cross
      write(*, out_opt%e2_fmt), "Bundle Expand Scalar =", field_bundle%exp_scalar
      write(*, out_opt%e2_fmt) "r = Pt/Ps =", r, '(', r_pred, ')'

      ! [JF] This SR expression should hold for an arbitrary number of fields but I should check more carefully (holds for 2 for sure)
      write(*, out_opt%e2_fmt) "n_s =", ns(1), '(', ns_pred ,')'
      if (size(ns)>1) then
        write(*, out_opt%e2_fmt) "n_iso =", ns(2)
        write(*, out_opt%e2_fmt) "n_pnad =", ns(3)
        write(*, out_opt%e2_fmt) "n_ent =", ns(4)
      end if
      write(*, out_opt%e2_fmt) "n_t =", nt, '(', nt_pred , ')'
      write(*, out_opt%e2_fmt) "alpha_s =", alpha_s, '(', alphas_pred , ')'


      if (calc_full_pk) then

        !Open output files
        open(newunit=out_adiab,file="out_pk_adiab.txt")
        open(newunit=out_isoc, file="out_pk_isocurv.txt")

        do i=1,size(pk_arr,1)
          write(out_adiab,*) pk_arr(i,:)
          if(present(pk_iso_arr)) write(out_isoc,*) pk_iso_arr(i,:)
        end do
        write(*,*) "Adiab P(k) written to out_pk_adiab.txt"
        if (present(pk_iso_arr)) write(*,*) "Iso-P(k) written to out_pk_isocurv.txt"
      end if

    end subroutine output_observables

    subroutine output_initial_data()
      integer :: i

      call out_opt%formatting(num_inflaton)

      if (size(vparams,1)>1) then
        do i=1, size(vparams,1)
          if (out_opt%modpkoutput) &
            write(*, '(A8,I1,A5,100E12.3)'), "vparams(",i,",:) =", vparams(i,:)
        end do
      else
        if (out_opt%modpkoutput) &
          write(*, *), "vparams(1,:) =", vparams(1,:)
      end if
    end subroutine output_initial_data


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

        if (varying_N_pivot) then
          save_iso_N = .false.
          call get_new_N_pivot(N_pivot, N_pivot_prior_min, N_pivot_prior_max)
        end if

      end if

      !Initialize potential and calc background
      call potinit

      call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
        A_iso, A_pnad, A_ent, A_bundle, &
        n_iso, n_pnad, n_ent, &
        leave)
      if (leave) return

      if (use_deltaN_SR) then
        call calculate_SR_observables()
      end if


      if (.not. evaluate_modes) return

      !Evaluate the mode functions
      call evolve(k_pivot, pk0)
        call test_bad(pk_bad,As,ns,r,nt,alpha_s,&
          A_iso, A_pnad, A_ent, A_bundle, &
          n_iso, n_pnad, n_ent, &
          leave)
        if (leave) return
!DEBUG
!print*, "Not evaluating second and third evolve routines"
!stop
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
      ns = 1.e0_dp+log(ps2/ps1)/dlnk/2.e0_dp
      r=pt0/ps0
      nt=log(pt2/pt1)/dlnk/2.e0_dp

      alpha_s = log(ps2*ps1/ps0**2)/dlnk**2


      !alpha_s from 5-pt stencil
      !alpha_s = (1.0e0_dp/12.0e0_dp/dlnk**2)*&
      !  (-log(pk4%adiab) + 16.0e0_dp*log(pk2%adiab) - &
      !  30.0e0_dp*log(pk0%adiab) + 16.0e0_dp*log(pk1%adiab) - &
      !  log(pk3%adiab))

      !runofrun = (1.0e0_dp/2.0e0_dp/dlnk**3)*&
      !  (log(pk4%adiab) -2* log(pk2%adiab) + 2*log(pk1%adiab) -log(pk3%adiab))

      !print*, "running of running =", runofrun

      n_iso=log(ps2_iso/ps1_iso)/dlnk/2.e0_dp
      n_pnad=log(pnad2/pnad1)/dlnk/2.e0_dp
      n_ent=log(pent2/pent1)/dlnk/2.e0_dp


      !Get full spectrum for adiab and isocurv at equal intvs in lnk
      call get_full_pk(pk_arr,pk_iso_arr,calc_full_pk)

      if (out_opt%modpkoutput) then
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
        if (out_opt%output_badic .or. pk_bad/=bad_ic) then
          call ic_output(i)%printout(out_opt%outsamp)
        endif

        if (save_iso_N) then
          call ic_output_iso_N(i)%load_observables(phi_iso_N, dphi_iso_N, &
            As,ns,r,nt, alpha_s,&
            A_iso, A_pnad, A_ent, A_bundle, n_iso, n_pnad, n_ent,&
            A_cross)
          if (out_opt%output_badic .or. pk_bad/=bad_ic) then
            call ic_output_iso_N(i)%printout(out_opt%outsamp_N_iso)
          end if
        end if
      end if


    end subroutine calculate_pk_observables

    !Calculate observables for the power spectrum, as well as fNL, using the
    !delta-N formalism in slow-roll
    subroutine calculate_SR_observables()
      integer :: j, i
      real(dp) :: ah, alpha_ik, dalpha, N_end, del_N, Npiv_renorm
      real(dp), dimension(num_inflaton) :: phi_pivot, phi_end, del_phi

      !Find field values at end of inflation
      !Note that eps=1 perhaps twice, so take the last one.
      CALL array_polint(epsarr(nactual_bg-4:nactual_bg), phiarr(:,nactual_bg-4:nactual_bg),&
        1.0e0_dp,  phi_end, del_phi)
      CALL polint(epsarr(nactual_bg-4:nactual_bg), lna(nactual_bg-4:nactual_bg),&
        1.0e0_dp,  N_end, del_N)

      !Find field values at horizon crossing
      Npiv_renorm = N_end - N_pivot

      i= locate(lna(1:nactual_bg), Npiv_renorm)
      j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
      CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), Npiv_renorm, phi_pivot, del_phi)

      print*, "testing SR approx"
      print*, "P_R", PR_SR(phi_pivot,phi_end)
      print*, "ns", ns_SR(phi_pivot,phi_end)
      print*, "nt", nt_SR(phi_pivot)
      print*, "r", r_SR(phi_pivot,phi_end)
      print*, "fnl", fnl_SR(phi_pivot,phi_end), (-5.0/6.0)/N_pivot



    end subroutine calculate_SR_observables


    subroutine test_bad(pk_bad,As,ns,r,nt,alpha_s,&
      A_iso, A_pnad, A_ent, A_bundle, &
      n_iso, n_pnad, n_ent, &
      leave)

      integer,  intent(in)     :: pk_bad
      logical,  intent(inout)  :: leave
      real(dp), intent(inout)  :: As,ns,r,nt, alpha_s
      real(dp), intent(inout)  :: A_iso, A_pnad, A_ent, A_bundle
      real(dp), intent(inout)  :: n_iso, n_pnad, n_ent

      !If pk_bad==bad_ic, then restart IC
      !If pk_bad==4, then ode_underflow

      if (pk_bad==bad_ic .or. pk_bad==4) then
        As = 0e0_dp
        ns = 0e0_dp
        r = 0e0_dp
        nt = 0e0_dp
        alpha_s = 0e0_dp
        A_iso = 0e0_dp
        A_pnad = 0e0_dp
        A_ent = 0e0_dp
        A_bundle = 0e0_dp
        n_iso = 0e0_dp
        n_pnad = 0e0_dp
        n_ent = 0e0_dp
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
        N_pivot_prior_min, N_pivot_prior_max

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


end program multimodecode
