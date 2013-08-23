!Some routines for post-processing the IC data
module modpk_postprocessing
  use modpkparams, only : dp, num_inflaton
  use modpk_icsampling, only : ic_and_observables
  implicit none

  contains

    function calc_clumping_penalty(ics, sigma) result(E)

      type(ic_and_observables), dimension(:), intent(in) :: ics
      real(dp), intent(in) :: sigma
      !real(dp),dimension(size(ics),size(ics)) :: Eij
      real(dp),dimension(size(ics)) :: E
      real(dp),dimension(2*num_inflaton) :: xi,xj,rad

      integer :: i, j

      E=0e0_dp
      do i=1,size(ics)
        xi = ics(i)%ic(:)
        do j=1, size(ics)
          xj = ics(j)%ic(:)
          rad = xi-xj
          E(i) = E(i) + exp(-1e0_dp/2.0e0_dp/sigma**2/sum(rad*rad))
        end do
      end do

    end function calc_clumping_penalty

end module modpk_postprocessing
