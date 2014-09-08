!Module for handling fatal exceptions.
!Also, code-related warnings are handled here, while cosmology-related warnings
!are printed to screen where they occur.
module modpk_errorhandling
  implicit none

  private
  public :: raise

  type :: error
    contains
      procedure, public :: fatal_code => raise_exception_code
      procedure, public :: fatal_cosmo => raise_exception_cosmo
      procedure, public :: warning => raise_warning
  end type

  type(error) :: raise

  character(30) :: msg_fmt
  integer :: ii

  contains

    !Code-related warnings handled here.
    subroutine raise_warning(self, msg, fname, line)
      class(error) :: self
      character(100), dimension(:), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line


      print*, "**********************************************"
      if (present(fname)) &
        print*, "WARNING: Occured in file ", fname
      if (present(fname)) &
        print*, "WARNING: At line number ", line

      do ii=1, size(msg)
        print*, "WARNING: ", trim(msg(ii))
      end do

    end subroutine raise_warning

    subroutine raise_exception_code(self, msg, fname, line)
      class(error) :: self
      character(*), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line


      print*, "**********************************************"
      print*, "MODECODE: Encountered a fatal exception."
      if (present(fname)) &
        print*, "MODECODE: Occured in file ", fname
      if (present(fname)) &
        print*, "MODECODE: At line number ", line
      print*, "MODECODE: Error is code-related."
      print*, "**********************************************"

      print*, "MODECODE: ", trim(msg)

      !Fatal
      stop

    end subroutine raise_exception_code


    subroutine raise_exception_cosmo(self, msg, fname, line)
      class(error) :: self
      character(*), intent(in) :: msg
      character(*), intent(in), optional :: fname
      integer, intent(in), optional :: line


      print*, "**********************************************"
      print*, "MODECODE: Encountered a fatal exception."
      if (present(fname)) &
        print*, "MODECODE: Occured in file ", fname
      if (present(fname)) &
        print*, "MODECODE: At line number ", line
      print*, "MODECODE: Error is cosmology-related."
      print*, "**********************************************"

      print*, "MODECODE: ", trim(msg)

      !Fatal
      stop

    end subroutine raise_exception_cosmo

end module modpk_errorhandling
