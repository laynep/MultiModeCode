MODULE camb_interface
  !Simple module to tell an external program, such as CAMB or the
  !multimodecode_driver, that a given set of parameters does not give an
  !appropriate inflationary realization.
  IMPLICIT NONE
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
END MODULE camb_interface

MODULE modpkparams
  !Module defining many "global" variables that various cosmology portions of
  !the code will need access to.  A variable added here can be seen in most
  !places.

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
  INTEGER, PARAMETER :: nsteps=1e6
  real(dp), PARAMETER :: M_Pl=1.0e0_dp
  real(dp), PARAMETER :: Mpc2Mpl=2.6245e-57_dp
  real(dp) :: k_pivot, N_pivot, N_tot, H_pivot
  real(dp) :: a_end, a_pivot
  real(dp) :: a_init
  real(dp) :: h_init, rescale_factor

  real(dp), ALLOCATABLE :: phidot_sign(:)
  real(dp) :: Nefold_max=100000.e0_dp
  real(dp) :: t_max=1.e100_dp
  real(dp), dimension(nsteps*2) :: lna, hubarr, log_aharr, epsarr, dtheta_dN
  LOGICAL :: slowroll_infl_end
  LOGICAL :: slowroll_start=.false.

  !MULTIFIELD
  integer :: num_inflaton
  real(dp), dimension(:,:), allocatable :: vparams
  real(dp), dimension(:), allocatable :: auxparams
  integer :: numb_auxparams
  real(dp), allocatable :: phi_init0(:), phi_init(:)
  real(dp), allocatable :: dphi_init0(:), dphi_init(:), dphidt_init0(:)
  real(dp), allocatable:: phi_pivot(:), dphi_pivot(:), phi_infl_end(:), dphi_infl_end(:)

  real(dp), allocatable :: phiarr(:,:), dphiarr(:,:) !The first index is the multifield index
  real(dp), allocatable :: param_arr(:)
  real(dp) :: sigma_arr(nsteps)
  real(dp) :: delsigma = M_Pl      !specify total field distance travelled before inflation ends
  !END MULTIFIELD

  real(dp) :: findiffdphi

  real(dp) :: modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r

  !Flags for analytical calculations
  logical :: use_deltaN_SR
  logical :: evaluate_modes

  !Technical options
  type :: tech_options
    !Accuracy options
    integer :: accuracy_setting
    !Which integrator to use
    logical :: use_dvode_integrator
    logical :: use_integ_with_t
    !Employ scalar constraints on auxiliary parameters
    logical :: use_ode_constraints
    !Whether we have an analytic Jacobian (not yet ready)
    logical :: use_analytical_jacobian
    !Use option to find appropriate IC for single-field inflation
    logical :: automate_singlefield_ic

    real(dp) :: rk_accuracy_modes, rk_accuracy_back
    real(dp) :: dvode_rtol_modes, dvode_rtol_back
    real(dp), dimension(10000) :: dvode_atol_modes_real, dvode_atol_modes_imag, dvode_atol_back
    real(dp) :: dvode_dN_r, dvode_dN_c, dvode_dt

  end type tech_options

  type(tech_options) :: tech_opt

END MODULE modpkparams


MODULE ode_path
  !Module for control parameters for the numerical integration of either the
  !background equations or the mode equations.
  use modpkparams, only : dp
  implicit none

  INTEGER*4 :: nok,nbad,kount, kount_t
  LOGICAL, SAVE :: save_steps=.false.
  LOGICAL :: ode_underflow
  LOGICAL :: ode_ps_output
  LOGICAL :: ode_infl_end
  LOGICAL :: infl_ended
  real(dp) :: dxsav
  real(dp), DIMENSION(:), POINTER :: xp, xp_t
  real(dp), DIMENSION(:), POINTER :: param_p, param_p_t
  real(dp), DIMENSION(:,:), POINTER :: yp, yp_t
END MODULE ode_path


MODULE internals
  !Module that defines some variables for internal use in the code.
  use modpkparams, only : num_inflaton, dp
  IMPLICIT NONE
  real(dp), PARAMETER :: PI=3.141592653589793238462643383279502884197e0_dp
  real(dp) :: h_ik
  !MULTIFIELD
  !integer :: index_ptb_y, index_ptb_vel_y, index_tensor_y, index_uzeta_y
  !END MULTIFIELD
  real(dp) :: k, a_ik


END MODULE internals
