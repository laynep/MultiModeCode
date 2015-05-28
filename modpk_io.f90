module modpk_io
  !Behold the beauty that is Fortran IO.
  !Controls output options for the run.
  implicit none

  private
  public :: out_opt

  type :: detailed_print_options
    !Do we need more detailed options
    logical :: needmoreopts

    !Print out other relevant parameters
    logical :: write_Npiv
    logical :: write_num_inflaton
    logical :: write_vparams
    logical :: write_auxparams
  end type

  type :: print_options

    !Verbosity
    logical :: modpkoutput
    logical :: output_reduced

    !Things to print
    logical :: save_traj
    logical :: output_badic
    logical :: fields_horiz
    logical :: fields_end_infl
    logical :: spectra
    logical :: modes
    logical :: phi0

    !File unit numbers for output
    integer :: trajout
    integer :: spectraout
    integer :: fields_h_out
    integer :: fields_end_out
    integer :: outsamp
    integer :: outsamp_SR
    integer :: outsamp_N_iso
    integer :: outsamp_N_iso_SR
    integer, dimension(4) :: modeout

    !If first write, then make column headers
    logical :: first_trajout = .true.
    logical :: first_spectraout = .true.
    logical :: first_fields_h_out = .true.
    logical :: first_fields_end_out = .true.
    logical :: first_outsamp = .true.
    logical :: first_outsamp_SR = .true.
    logical :: first_outsamp_N_iso = .true.
    logical :: first_outsamp_N_iso_SR = .true.
    logical :: first_modeout = .true.

    !Writing fmts
    character(16) :: e_fmt = '(a25, 900es12.4)'
    character(36) :: e2_fmt = '(a25, es17.9, a3, es16.9, a1)'
    character(16) :: i_fmt = '(a25,I3)'
    character(16) :: array_fmt
    character(len=2) :: ci

    !Detailed output options
    type (detailed_print_options) :: detailed

    contains

      procedure, public :: open_files => output_file_open
      procedure, public :: close_files => output_file_close
      procedure, public :: formatting => make_formatting

  end type

  !Main global instance for printing options
  type(print_options) :: out_opt

  contains

    !Open output files
    subroutine output_file_open(self,ICs,SR)
      class(print_options) :: self
      logical, intent(in), optional :: ICs, SR

      if (self%save_traj) &
        open(newunit=self%trajout, &
          file="out_trajectory.csv")
      if (self%spectra) &
        open(newunit=self%spectraout, &
          file="out_powerspectra.csv")
      if (self%fields_horiz) &
        open(newunit=self%fields_h_out, &
          file="out_fields_horizon_cross.csv")
      if (self%fields_end_infl) &
        open(newunit=self%fields_end_out, &
          file="out_fields_infl_end.csv")
      if (self%modes) then
        open(newunit=self%modeout(1), &
          file="out_modes_1.csv")
        open(newunit=self%modeout(2), &
          file="out_modes_2.csv")
        open(newunit=self%modeout(3), &
          file="out_modes_3.csv")
        open(newunit=self%modeout(4), &
          file="out_modes_4.csv")
      end if

      if (present(ICs) .and. ICs) then
        open(newunit=self%outsamp,&
          file="out_ic_eqen.csv")
        open(newunit=self%outsamp_N_iso,&
          file="out_ic_isoN.csv")
      end if

      if (present(SR) .and. SR) then
        open(newunit=self%outsamp_SR,&
          file="out_ic_eqen_SR.csv")
        open(newunit=self%outsamp_N_iso_SR,&
          file="out_ic_isoN_SR.csv")
      end if

    end subroutine output_file_open

    !Close output files
    subroutine output_file_close(self,ICs,SR)
      class(print_options) :: self
      logical, intent(in), optional :: ICs, SR

      if (self%save_traj) &
        close(self%trajout)
      if (self%spectra) &
        close(self%spectraout)
      if (self%fields_horiz) &
        close(self%fields_h_out)
      if (self%fields_end_infl) &
        close(self%fields_end_out)
      if (self%modes) then
        close(self%modeout(1))
        close(self%modeout(2))
        close(self%modeout(3))
        close(self%modeout(4))
      end if

      if (present(ICs) .and. ICs) then
        close(self%outsamp)
        close(self%outsamp_N_iso)
      end if

      if (present(SR) .and. SR) then
        close(self%outsamp_SR)
        close(self%outsamp_N_iso_SR)
      end if

    end subroutine output_file_close

    subroutine make_formatting(self, num_inflaton)
      class(print_options) :: self
      integer, intent(in) :: num_inflaton

      write(self%ci, '(I2)'), num_inflaton
      self%ci = adjustl(self%ci)
      self%array_fmt = '(a25,'//trim(self%ci)//'es10.3)'

    end subroutine make_formatting


end module modpk_io
