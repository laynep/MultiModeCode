!A module that implements perturbative reheating after the end of
!inflation.

module modpk_reheat
  use modpkparams, only : dp, slowroll_infl_end, vparams
  use potential, only : getH, geteps, getkineticenergy, getw
  use modpk_errorhandling, only : raise

  implicit none

  logical :: use_reheat=.false.
  type :: reheat_model_flags
    integer :: perturbative = 1
  end type reheat_model_flags
  integer :: reheat_model

  contains

    subroutine reheat_ending_check(phi,dphidN)

      real(dp), intent(in), dimension(:) :: phi, dphidN
      integer :: ii

      !print*, "testing reheat_ending_check"

      call consistency_checks()

      select case(reheat_model)
      case default

        call raise%fatal_code("Please specify your reheating model.",__FILE__,__LINE__)

      case(1)

      end select


      !print*, "this is rho:", &
      print*,((1.0e0_dp/vparams(2,1))*vparams(1,ii)*abs(phi(ii))**(vparams(2,1)) + &
        0.5e0_dp*dphidN(ii)**2*getH(phi,dphidN)**2,ii=1,size(phi))

      !print*,sum((1.0e0_dp/vparams(2,1))*vparams(1,:)*abs(phi(:))**(vparams(2,1)) + &
      !  0.5e0_dp*dphidN(:)**2*getH(phi,dphidN)**2)

      !print*, "eps", geteps(phi,dphidN)
      !print*, "w", getw(phi,dphidN)
      !print*, "m/H",sqrt(vparams(1,:))/getH(phi,dphidN)


      !print*, geteps(phi,dphidN), getw(phi,dphidN)

      !print*, "phi0", phi
      !print*, "dphi0", dphidN

      !if (geteps(phi,dphidN)>1.0) stop




    end subroutine reheat_ending_check

    subroutine consistency_checks()

      if (use_reheat) then
        slowroll_infl_end = .false.
      end if

    end subroutine consistency_checks

end module modpk_reheat
