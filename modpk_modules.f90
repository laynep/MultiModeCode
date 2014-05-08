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

END MODULE modpkparams

module modpk_output
  implicit none

  type :: print_options
    logical :: modpkoutput
    logical :: save_traj
    logical :: output_badic
    logical :: fields_horiz
    logical :: fields_end_infl
    logical :: spectra

    !File unit numbers for output
    integer :: trajout
    integer :: spectraout
    integer :: fields_h_out
    integer :: fields_end_out

    contains
      procedure :: open_files => output_file_open

  end type

  type(print_options) :: out_opt

  contains

    subroutine output_file_open(this)
      class(print_options) :: this

      !Open output files
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

    end subroutine output_file_open

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


module powersp
  use modpkparams, only : dp
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


end module powersp
