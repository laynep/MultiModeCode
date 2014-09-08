!Module for handling exceptions and warnings.
!More info is sometimes printed to screen where these occur.
module modpk_errorhandling
  use modpk_io, only : out_opt
  implicit none

  private
  public :: raise, run_outcome

  !General class for errors
  type :: error
    contains
      procedure, public :: fatal_code => raise_exception_code
      procedure, public :: fatal_cosmo => raise_exception_cosmo
      procedure, public :: warning => raise_warning
  end type

  !Type that contains indices for classifying bad sets of parameters
  type :: run_outcome_type
    integer :: success = 0
    integer :: infl_didnt_start = 1
    integer :: pivot_didnt_leaveH = 2
    integer :: cant_set_modeIC = 3
    integer :: cant_init_scalefact = 4
    integer :: ref_efold_didnt_leaveH = 5
    integer :: bad_reheat = 6
    integer :: underflow = -1

    contains
      procedure, public :: is_fatal => compare_params_fatal
      procedure, public :: print_outcome => output_the_outcome

  end type

  type(error) :: raise
  type(run_outcome_type) :: run_outcome

  contains

    !Code-related warnings handled here.
    subroutine raise_warning(self, msg, fname, line)
      class(error) :: self
      character(*), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line

      if (out_opt%modpkoutput) then
        if (.not. out_opt%output_reduced) then

          if (present(fname)) &
            print*, "WARNING: Occured in file ", fname
          if (present(fname)) &
            print*, "WARNING: At line number ", line

        end if

        print*, "MODECODE: ", trim(msg)

      end if

    end subroutine raise_warning

    subroutine raise_exception_code(self, msg, fname, line)
      class(error) :: self
      character(*), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line


      !Print out even if asking for no output (modpkoutput=.false.)
      print*, "**********************************************"
      print*, "ERROR: Encountered a fatal exception."
      if (present(fname)) &
        print*, "ERROR: Occured in file ", fname
      if (present(line)) &
        print*, "ERROR: At line number ", line
      print*, "ERROR: The error is code-related."
      print*, "**********************************************"

      print*, "MODECODE: ", trim(msg)

      stop

    end subroutine raise_exception_code


    subroutine raise_exception_cosmo(self, msg, fname, line)
      class(error) :: self
      character(*), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line


      !Print out even if asking for no output (modpkoutput=.false.)
      print*, "**********************************************"
      print*, "ERROR: Encountered a fatal exception."
      if (present(fname)) &
        print*, "ERROR: Occured in file ", fname
      if (present(line)) &
        print*, "ERROR: At line number ", line
      print*, "ERROR: The error is cosmology-related."
      print*, "**********************************************"

      print*, "MODECODE: ", trim(msg)

      stop

    end subroutine raise_exception_cosmo

    subroutine compare_params_fatal(self, outcome)

      class(run_outcome_type) :: self
      integer, intent(in) :: outcome
      logical :: is_fatal

      if (outcome < 0) then
        !Fatal

        is_fatal = .true.

        call raise%fatal_cosmo(&
          "This set of parameters gives a fatal error &
          and it's not worth trying new parameters.",&
          __FILE__, __LINE__)

      else
        !Not fatal
        is_fatal = .false.
      end if

    end subroutine compare_params_fatal

    subroutine output_the_outcome(self, outcome)

      class(run_outcome_type) :: self
      integer, intent(in) :: outcome
      logical :: is_fatal

      if (outcome == self%success) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Successful set of parameters.", outcome

      else if (outcome == self%infl_didnt_start) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Inflation didn't start.", outcome

      else if (outcome == self%pivot_didnt_leaveH) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: The pivot scale didn't leave the horizon.", outcome

      else if (outcome == self%cant_set_modeIC) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: A modes' IC couldn't be consistently set.", outcome

      else if (outcome == self%cant_init_scalefact) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Too many e-folds; can't initialize the scale factor.", outcome

      else if (outcome == self%ref_efold_didnt_leaveH) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Trying to save some reference N, &
          but got less evolution than that.", outcome

      else if (outcome == self%bad_reheat) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Didn't satisfy reheating bounds.", outcome

      else if (outcome == self%underflow) then
        if (out_opt%modpkoutput) &
          print*, "MODECODE: Numerical underflow error in odeint.", outcome

      else

        print*, "MODECODE: pk_bad=", outcome
        call raise%fatal_code(&
          "The outcome flag pk_bad was set to a value that &
          wasn't caught by our run_outcome type.  You shouldn't &
          be seeing this.", &
          __FILE__, __LINE__)

      end if

      !Check if the whole run should be stopped because of the outcome
      call self%is_fatal(outcome)

      print*, "............... RESTARTING ..............."

    end subroutine output_the_outcome

end module modpk_errorhandling
