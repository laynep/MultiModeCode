MODULE camb_interface
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
  LOGICAL :: modpkoutput=.false.
END MODULE camb_interface

MODULE ode_path
  INTEGER*4 :: nok,nbad,kount
  LOGICAL, SAVE :: save_steps=.false.
  LOGICAL :: ode_underflow
  LOGICAL :: ode_ps_output
  LOGICAL :: ode_infl_end
  LOGICAL :: infl_ended
  DOUBLE PRECISION :: dxsav
  DOUBLE PRECISION, DIMENSION(:), POINTER :: xp
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

MODULE modpkparams
  IMPLICIT NONE

  LOGICAL :: use_modpk, vnderivs, instreheat

! increase max_vparams to use more potential parameters
  INTEGER*4, parameter :: max_vparams = 9

  INTEGER :: potential_choice

  INTEGER*4 :: nactual_bg, nactual_mode
  INTEGER*4, PARAMETER :: nsteps=10000
  DOUBLE PRECISION, PARAMETER :: M_Pl=1.0d0
  DOUBLE PRECISION, PARAMETER :: Mpc2Mpl=2.6245d-57
  DOUBLE PRECISION :: k_pivot, N_pivot, N_tot, H_pivot
  DOUBLE PRECISION :: a_end, a_pivot
  DOUBLE PRECISION :: a_init
  DOUBLE PRECISION :: h_init, rescale_factor

  DOUBLE PRECISION, ALLOCATABLE :: phidot_sign(:)
  DOUBLE PRECISION :: Nefold_max=10000.
  DOUBLE PRECISION :: lna(nsteps)
  DOUBLE PRECISION :: hubarr(nsteps), aharr(nsteps), epsarr(nsteps), dtheta_dN(nsteps)
  LOGICAL :: slowroll_infl_end
  LOGICAL :: slowroll_start=.false.

  !MULTIFIELD
  INTEGER, parameter :: DP = selected_real_kind(kind(1.0d0))
  INTEGER :: num_inflaton
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vparams
  DOUBLE PRECISION, ALLOCATABLE :: phi_init0(:), phi_init(:)
  DOUBLE PRECISION, ALLOCATABLE:: phi_pivot(:), dphi_pivot(:), phi_infl_end(:)
  
  DOUBLE PRECISION, ALLOCATABLE :: phiarr(:,:), dphiarr(:,:) !The first index is the multifield indice
  DOUBLE PRECISION :: sigma_arr(nsteps)
  DOUBLE PRECISION :: delsigma = M_Pl      !specify total field distance travelled before inflation ends
  !END MULTIFIELD

  DOUBLE PRECISION :: findiffdphi

  DOUBLE PRECISION :: modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r

END MODULE modpkparams


MODULE internals
  use modpkparams, only : num_inflaton
  IMPLICIT NONE
  REAL, PARAMETER :: PI=3.141592653589793238462643383279502884197
  DOUBLE PRECISION :: h_ik, powt_ik, powz_ik
  ![ LP: ] MULTIFIELD
  double precision :: pow_adiab_ik, pow_isocurv_ik
  double precision, dimension(:,:), allocatable :: pow_ptb_ij
  integer :: index_ptb_y, index_ptb_vel_y, index_tensor_y, index_uzeta_y
  !END MULTIFIELD
  DOUBLE PRECISION :: k, a_ik


END MODULE internals


MODULE powersp
  IMPLICIT NONE
  INTEGER*4 :: ik
  DOUBLE PRECISION :: eval_ps,k_start, useq_ps
END MODULE powersp
