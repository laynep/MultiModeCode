MODULE camb_interface
  IMPLICIT NONE
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
END MODULE camb_interface

MODULE modpkparams
  IMPLICIT NONE

  !Double precision.
  integer, parameter :: dp = selected_real_kind(15, 307)

  !Quad precision; remember to compile with -r16
  !INTEGER, parameter :: DP = selected_real_kind(33, 4931)

  LOGICAL :: use_modpk, vnderivs, instreheat

! increase max_vparams to use more potential parameters
  INTEGER*4, parameter :: max_vparams = 9

  INTEGER :: potential_choice

  INTEGER*4 :: nactual_bg, nactual_mode
  INTEGER, PARAMETER :: nsteps=1e5
  real(dp), PARAMETER :: M_Pl=1.0e0_dp
  real(dp), PARAMETER :: Mpc2Mpl=2.6245e-57_dp
  real(dp) :: k_pivot, N_pivot, N_tot, H_pivot
  real(dp) :: a_end, a_pivot
  real(dp) :: a_init
  real(dp) :: h_init, rescale_factor

  real(dp), ALLOCATABLE :: phidot_sign(:)
  real(dp) :: Nefold_max=100000.e0_dp
  real(dp) :: lna(nsteps)
  real(dp) :: hubarr(nsteps), log_aharr(nsteps), epsarr(nsteps), dtheta_dN(nsteps)
  LOGICAL :: slowroll_infl_end
  LOGICAL :: slowroll_start=.false.

  !MULTIFIELD
  INTEGER :: num_inflaton
  real(dp), DIMENSION(:,:), ALLOCATABLE :: vparams
  real(dp), ALLOCATABLE :: phi_init0(:), phi_init(:)
  real(dp), ALLOCATABLE :: dphi_init0(:), dphi_init(:)
  real(dp), ALLOCATABLE:: phi_pivot(:), dphi_pivot(:), phi_infl_end(:)

  real(dp), ALLOCATABLE :: phiarr(:,:), dphiarr(:,:) !The first index is the multifield index
  real(dp) :: sigma_arr(nsteps)
  real(dp) :: delsigma = M_Pl      !specify total field distance travelled before inflation ends
  !END MULTIFIELD

  real(dp) :: findiffdphi

  real(dp) :: modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r

  logical :: use_deltaN_SR
  logical :: evaluate_modes

END MODULE modpkparams

!Behold the beauty that is Fortran IO.
module modpk_output
  implicit none

  type :: print_options

    !Things to print
    logical :: modpkoutput
    logical :: output_reduced
    logical :: save_traj
    logical :: output_badic
    logical :: fields_horiz
    logical :: fields_end_infl
    logical :: spectra
    logical :: modes

    !File unit numbers for output
    integer :: trajout
    integer :: spectraout
    integer :: fields_h_out
    integer :: fields_end_out
    integer :: outsamp
    integer :: outsamp_N_iso
    integer, dimension(4) :: modeout

    !Writing fmts
    character(16) :: e_fmt = '(a25, es12.4)'
    character(36) :: e2_fmt = '(a25, es17.9, a3, es11.4, a1)'
    character(16) :: i_fmt = '(a25,I3)'
    character(16) :: array_fmt
    character(len=2) :: ci

    contains
      procedure :: open_files => output_file_open
      procedure :: formatting => make_formatting

  end type

  type(print_options) :: out_opt

  contains

    !Open output files
    subroutine output_file_open(this,ICs)
      class(print_options) :: this
      logical, intent(in), optional :: ICs

      if (this%save_traj) &
        open(newunit=this%trajout, &
          file="out_trajectory.txt")
      if (this%spectra) &
        open(newunit=this%spectraout, &
          file="out_powerspectra.txt")
      if (this%fields_horiz) &
        open(newunit=this%fields_h_out, &
          file="out_fields_horizon_cross.txt")
      if (this%fields_end_infl) &
        open(newunit=this%fields_end_out, &
          file="out_fields_infl_end.txt")
      if (this%modes) then
        open(newunit=this%modeout(1), &
          file="out_modes_1.txt")
        open(newunit=this%modeout(2), &
          file="out_modes_2.txt")
        open(newunit=this%modeout(3), &
          file="out_modes_3.txt")
        open(newunit=this%modeout(4), &
          file="out_modes_4.txt")
      end if

      if (present(ICs) .and. ICs) then
        open(newunit=this%outsamp,&
          file="out_ic_eqen.txt")
        open(newunit=this%outsamp_N_iso,&
          file="out_ic_isoN.txt")
      end if

    end subroutine output_file_open

    subroutine make_formatting(this, num_inflaton)
      class(print_options) :: this
      integer, intent(in) :: num_inflaton

      write(this%ci, '(I2)'), num_inflaton
      this%ci = adjustl(this%ci)
      this%array_fmt = '(a25,'//trim(this%ci)//'es10.3)'

    end subroutine make_formatting


end module modpk_output

MODULE ode_path
  use modpkparams, only : dp
  implicit none

  INTEGER*4 :: nok,nbad,kount
  LOGICAL, SAVE :: save_steps=.false.
  LOGICAL :: ode_underflow
  LOGICAL :: ode_ps_output
  LOGICAL :: ode_infl_end
  LOGICAL :: infl_ended
  real(dp) :: dxsav
  real(dp), DIMENSION(:), POINTER :: xp
  real(dp), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path


MODULE internals
  use modpkparams, only : num_inflaton, dp
  IMPLICIT NONE
  real(dp), PARAMETER :: PI=3.141592653589793238462643383279502884197e0_dp
  real(dp) :: h_ik
  !MULTIFIELD
  integer :: index_ptb_y, index_ptb_vel_y, index_tensor_y, index_uzeta_y
  !END MULTIFIELD
  real(dp) :: k, a_ik


END MODULE internals


!Module that holds objects for observables, power spectra, etc
module modpk_observables
  use modpkparams, only : dp, num_inflaton
  implicit none

  integer*4 :: ik
  real(dp) :: eval_ps,k_start, useq_ps

  !Power spectrum type, used for one k
  !Simplifies all the various defns for isocurv ptbs
  type :: power_spectra
    !Mode
    real(dp) :: k
    !Ptb spectra
    complex(dp), dimension(:,:), allocatable :: matrix
    !Proj ptb spectra onto adiab direction
    real(dp) :: adiab
    !Tensor mode spectrum
    real(dp) :: tensor
    !Approximated zeta spectrum
    real(dp) :: powz
    !Proj ptb spectra onto directions perpend to adiab
    real(dp) :: isocurv
    !Cross-spectrum for isocurv + adiab
    real(dp) :: cross_ad_iso
    !Total non-adiab pressure ptb spectrum
    real(dp) :: pnad
    !Total entropic ptb spectrum, proportional to non-adiab pressure ptb
    real(dp) :: entropy
    !Expansion scalar for field space bundle width
    real(dp) :: bundle_exp_scalar

  end type

  !For temporary calc of spectra in odeint
  type(power_spectra) :: power_internal

  type, public :: KahanSum
    real(dp) :: summand = 0e0_dp
    real(dp) :: remainder = 0e0_dp
    contains
      procedure :: add => kahan_summation
      procedure :: clear => kahan_clear_mem
  end type

  !Type to save the ICs and observs. Add new observables here
  type :: observables
    real(dp), dimension(:), allocatable :: ic
    !Spectra amplitudes
    real(dp) :: As
    real(dp) :: A_iso
    real(dp) :: A_pnad
    real(dp) :: A_ent
    real(dp) :: A_cross_ad_iso
    !Bundle width from arXiv:1203.2635
    real(dp) :: A_bundle
    !Spectral indices
    real(dp) :: ns
    real(dp) :: nt
    real(dp) :: n_iso
    real(dp) :: n_pnad
    real(dp) :: n_ent
    !Tensor-to-scalar
    real(dp) :: r
    !Running, etc
    real(dp) :: alpha_s
    real(dp) :: runofrun
    !Non-Gaussianity
    real(dp) :: f_NL
    real(dp) :: tau_NL
    contains
      procedure :: printout => ic_print_observables
      procedure :: set_zero => set_observs_to_zero
  end type observables

  contains

    !Procedures for Kahan summation objects

    !Algorithm for Kahan summation
    subroutine kahan_summation(this,val)

      class(KahanSum) :: this
      real(dp), intent(in) :: val
      real(dp) :: ytemp, ttemp

      ytemp = val - this%remainder
      ttemp = this%summand + ytemp
      this%remainder = (ttemp - this%summand) - ytemp
      this%summand = ttemp

    end subroutine kahan_summation

    subroutine kahan_clear_mem(this)

      class(KahanSum) :: this

      this%remainder=0e0_dp
      this%summand=0e0_dp

    end subroutine kahan_clear_mem


    !Print the cosmo observables to file
    subroutine ic_print_observables(this, outunit)

      class(observables) :: this
      integer, intent(in) :: outunit

      if (num_inflaton*2 +12 > 120000) then
        print*, "Don't be silly."
        print*, "Too many fields to print out properly."
        print*, "Fix formatting."
        stop
      end if

      write(outunit, '(120000E18.10)') this%ic(:), this%As, this%ns,&
        this%r, this%nt, this%alpha_s, this%A_iso, this%A_pnad,&
        this%A_ent, this%A_bundle, this%n_iso, this%n_pnad, this%n_ent, &
        this%A_cross_ad_iso

    end subroutine ic_print_observables

    subroutine set_observs_to_zero(this)
      class(observables) :: this

        this%As = 0e0_dp
        this%A_iso = 0e0_dp
        this%A_pnad = 0e0_dp
        this%A_ent = 0e0_dp
        this%A_cross_ad_iso = 0e0_dp
        this%A_bundle = 0e0_dp
        this%ns = 0e0_dp
        this%nt = 0e0_dp
        this%n_iso = 0e0_dp
        this%n_pnad = 0e0_dp
        this%n_ent = 0e0_dp
        this%r = 0e0_dp
        this%alpha_s = 0e0_dp
        this%runofrun = 0e0_dp
        this%f_NL = 0e0_dp
        this%tau_NL = 0e0_dp

    end subroutine set_observs_to_zero

end module modpk_observables
