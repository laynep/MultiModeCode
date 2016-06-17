module modpk_observables
  !Module that defines the various observables one could calculate around the
  !pivot scale.  Defines objects for observables, power spectra, etc.
  use modpkparams, only : dp, num_inflaton, &
    vparams, auxparams, numb_auxparams, N_pivot, num_inflaton
  use modpk_errorhandling, only : raise
  use modpk_io, only : out_opt
  use csv_file, only : csv_write
  implicit none

  integer*4 :: ik
  real(dp) :: eval_ps,k_start, useq_ps

  !Power spectrum type, used for one k
  !Simplifies all the various defns for isocurv ptbs
  type :: power_spectra
    !Mode
    real(dp) :: k
    !Ptb spectra
    complex(dp), dimension(:,:), allocatable :: phi_ij
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
    !Total pressure ptb spectrum
    real(dp) :: pressure
    !Adiab pressure ptb spectrum
    real(dp) :: press_ad
    !Total entropic ptb spectrum, proportional to non-adiab pressure ptb
    real(dp) :: entropy
    !Expansion scalar for field space bundle width
    real(dp) :: bundle_exp_scalar

    contains

      procedure, public :: init => spectra_initializer

  end type

  !For temporary calc of spectra in odeint
  type(power_spectra) :: power_internal

  type, public :: KahanSum
    real(dp) :: summand = 0e0_dp
    real(dp) :: remainder = 0e0_dp
    contains
      procedure, public :: add => kahan_summation
      procedure, public :: clear => kahan_clear_mem
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
      procedure, public :: printout => ic_print_observables
      procedure, public :: print_header => ic_print_headers
      procedure, public :: set_zero => set_observs_to_zero
      procedure, public :: set_finite_diff => calculate_observs_finitediff
  end type observables

  private :: ic_print_observables, set_observs_to_zero, &
    calculate_observs_finitediff,kahan_summation,kahan_clear_mem

  contains

    !Procedures for Kahan summation objects

    !Algorithm for Kahan summation
    subroutine kahan_summation(self,val)

      class(KahanSum) :: self
      real(dp), intent(in) :: val
      real(dp) :: ytemp, ttemp

      ytemp = val - self%remainder
      ttemp = self%summand + ytemp
      self%remainder = (ttemp - self%summand) - ytemp
      self%summand = ttemp

    end subroutine kahan_summation

    subroutine kahan_clear_mem(self)

      class(KahanSum) :: self

      self%remainder=0e0_dp
      self%summand=0e0_dp

    end subroutine kahan_clear_mem

    !Write the cosmo observables headers to file
    subroutine ic_print_headers(self, outunit)
      class(observables) :: self
      integer, intent(in) :: outunit

      character(1024) :: cname
      integer :: ii, jj

      !First write any desired auxiliary parameters
      if (out_opt%detailed%needmoreopts) then

        if (out_opt%detailed%write_Npiv) &
          call csv_write(outunit, 'N_*',advance=.false.)

        if (out_opt%detailed%write_num_inflaton) &
          call csv_write(outunit, 'N_fields',advance=.false.)

        if (out_opt%detailed%write_vparams) then
          do ii=1,size(vparams,2); do jj=1,size(vparams,1)
            write(cname, "(A5,I4.4,A1,I4.4,A1)") "vpar(", jj, ",", ii, ")"
            call csv_write(&
              outunit,&
              trim(cname), &
              advance=.false.)
          end do; end do

        end if

        if (out_opt%detailed%write_auxparams) then
          do ii=1,numb_auxparams
            write(cname, "(A8,I4.4)") "auxparam", ii
            call csv_write(&
              outunit,&
              trim(cname), &
              advance=.false.)
          end do
        end if

      end if

      !First columns for IC
      if (out_opt%phi0) then
        do ii=1,size(self%ic)
          write(cname, "(A3,I4.4)") "phi_piv", ii
          call csv_write(&
            outunit,&
            trim(cname), &
            advance=.false.)
        end do
      end if

      !Remaining columns
      call csv_write(outunit,&
        (/character(len=10) ::&
          'As', 'ns', 'r', 'nt', 'alpha_s', &
          'A_iso', 'A_pnad', 'A_ent', 'A_bundle', &
          'n_iso', 'n_pnad', 'n_ent', &
          'A_cross', 'f_NL', 'tau_NL'/),&
        advance=.true.)

    end subroutine ic_print_headers


    !Print the cosmo observables to file
    subroutine ic_print_observables(self, outunit)

      class(observables) :: self
      integer, intent(in) :: outunit

      !First write the desired auxiliary parameters
      if (out_opt%detailed%needmoreopts) then

        if (out_opt%detailed%write_Npiv) &
          call csv_write(outunit,N_pivot,advance=.false.)

        if (out_opt%detailed%write_num_inflaton) &
          call csv_write(outunit,num_inflaton,advance=.false.)

        if (out_opt%detailed%write_vparams) then
          call csv_write(outunit,&
            (/ vparams/), advance=.false.)
        end if

        if (out_opt%detailed%write_auxparams) then
          call csv_write(outunit,&
            (/ auxparams(1:numb_auxparams)/), advance=.false.)
        end if

      end if

      if (out_opt%phi0) &
        call csv_write(outunit,&
          (/ self%ic(:)/), advance=.false.)

      call csv_write(outunit,&
        (/self%As, &
        self%ns,&
        self%r, &
        self%nt, &
        self%alpha_s, &
        self%A_iso, &
        self%A_pnad,&
        self%A_ent, &
        self%A_bundle, &
        self%n_iso, &
        self%n_pnad, &
        self%n_ent, &
        self%A_cross_ad_iso, &
        self%f_NL, &
        self%tau_NL /), &
        advance=.true.)

    end subroutine ic_print_observables

    subroutine set_observs_to_zero(self)
      class(observables) :: self

        self%As = 0e0_dp
        self%A_iso = 0e0_dp
        self%A_pnad = 0e0_dp
        self%A_ent = 0e0_dp
        self%A_cross_ad_iso = 0e0_dp
        self%A_bundle = 0e0_dp
        self%ns = 0e0_dp
        self%nt = 0e0_dp
        self%n_iso = 0e0_dp
        self%n_pnad = 0e0_dp
        self%n_ent = 0e0_dp
        self%r = 0e0_dp
        self%alpha_s = 0e0_dp
        self%runofrun = 0e0_dp
        self%f_NL = 0e0_dp
        self%tau_NL = 0e0_dp

    end subroutine set_observs_to_zero

    subroutine calculate_observs_finitediff(self, dlnk, &
        pk0, pklow1, pkhigh1, &
        pklow2, pkhigh2, &
        bundle_width)
      class(observables) :: self
      type(power_spectra), intent(in) :: pk0, pklow1, pkhigh1
      type(power_spectra), intent(in), optional :: pklow2, pkhigh2
      real(dp), intent(in) :: dlnk
      real(dp), intent(in), optional :: bundle_width
      logical :: runofrun

      if (present(pklow2)) then
        runofrun = .true.
      else
        runofrun = .false.
      endif

      !Amplitudes
      self%As = pk0%adiab
      self%A_iso=pk0%isocurv
      self%A_pnad=pk0%pnad
      self%A_ent=pk0%entropy
      self%A_cross_ad_iso = pk0%cross_ad_iso

      !Bundle width
      self%A_bundle=bundle_width

      !Finite difference evaluation of spectral indices
      self%ns = 1.e0_dp+log(pkhigh1%adiab/pklow1%adiab)/dlnk/2.e0_dp
      self%nt = log(pkhigh1%tensor/pklow1%tensor)/dlnk/2.e0_dp
      self%n_iso=log(pkhigh1%isocurv/pklow1%isocurv)/dlnk/2.e0_dp
      self%n_pnad=log(pkhigh1%pnad/pklow1%pnad)/dlnk/2.e0_dp
      self%n_ent=log(pkhigh1%entropy/pklow1%entropy)/dlnk/2.e0_dp

      !Tensor-to-scalar
      self%r = pk0%tensor/pk0%adiab

      if (runofrun) then

        !alpha_s from 5-pt stencil
        self%alpha_s = (1.0e0_dp/12.0e0_dp/dlnk**2)*&
          (-log(pkhigh2%adiab) + 16.0e0_dp*log(pkhigh1%adiab) - &
          30.0e0_dp*log(pk0%adiab) + 16.0e0_dp*log(pklow1%adiab) - &
          log(pklow2%adiab))

        self%runofrun = (1.0e0_dp/2.0e0_dp/dlnk**3)*&
          (log(pkhigh2%adiab) -2.0e0_dp* log(pkhigh1%adiab) &
          + 2.0e0_dp*log(pklow1%adiab) -log(pklow2%adiab))

      else

        self%alpha_s = log(pkhigh1%adiab*pklow1%adiab/pk0%adiab**2)/dlnk**2

        !Default
        self%runofrun = 0e0_dp

      end if

    end subroutine calculate_observs_finitediff

    subroutine spectra_initializer(self)
      class(power_spectra) :: self

      self%k = 0.0e0_dp
      self%adiab = 0.0e0_dp
      self%tensor = 0.0e0_dp
      self%powz = 0.0e0_dp
      self%isocurv = 0.0e0_dp
      self%cross_ad_iso = 0.0e0_dp
      self%pnad = 0.0e0_dp
      self%pressure = 0.0e0_dp
      self%press_ad = 0.0e0_dp
      self%entropy = 0.0e0_dp
      self%bundle_exp_scalar = 0.0e0_dp

      if (allocated(self%phi_ij)) &
        deallocate(self%phi_ij)

    end subroutine spectra_initializer

end module modpk_observables
