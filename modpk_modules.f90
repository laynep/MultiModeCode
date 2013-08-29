MODULE camb_interface
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
  LOGICAL :: modpkoutput=.false.
END MODULE camb_interface

MODULE modpkparams
  IMPLICIT NONE

  !Double precision.
  INTEGER, parameter :: DP = selected_real_kind(15, 307)

  LOGICAL :: use_modpk, vnderivs, instreheat

! increase max_vparams to use more potential parameters
  INTEGER*4, parameter :: max_vparams = 9

  INTEGER :: potential_choice

  INTEGER*4 :: nactual_bg, nactual_mode
  INTEGER, PARAMETER :: nsteps=100000
  real(dp), PARAMETER :: M_Pl=1.0d0
  real(dp), PARAMETER :: Mpc2Mpl=2.6245d-57
  real(dp) :: k_pivot, N_pivot, N_tot, H_pivot
  real(dp) :: a_end, a_pivot
  real(dp) :: a_init
  real(dp) :: h_init, rescale_factor

  real(dp), ALLOCATABLE :: phidot_sign(:)
  !DEBUG
  !real(dp) :: Nefold_max=10000.
  real(dp) :: Nefold_max=100000.
  real(dp) :: lna(nsteps)
  real(dp) :: hubarr(nsteps), aharr(nsteps), epsarr(nsteps), dtheta_dN(nsteps)
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

MODULE ode_path
  use modpkparams, only : dp
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
  real(dp), PARAMETER :: PI=3.141592653589793238462643383279502884197
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
    !Total non-adiab pressure ptb spectrum
    real(dp) :: pnad
    !Total entropic ptb spectrum, proportional to non-adiab pressure ptb
    real(dp) :: entropy
    !Expansion scalar for field space bundle width
    real(dp) :: bundle_exp_scalar

  end type

  !For temporary calc of spectra in odeint
  type(power_spectra) :: power_internal

end module powersp
