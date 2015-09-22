!Module that lets us compare the mode evolution to the slow-roll expectation for
!sum-separable potentials, using the results of Battefeld-Easther astro-ph/0610296
module modpk_deltaN
  use modpkparams, only : dp, vparams, N_pivot, potential_choice
  use modpk_observables, only : power_spectra
  use internals, only : pi
  use potential, only : pot, dVdphi, d2Vdphi2, d3Vdphi3
  use modpk_errorhandling, only : raise, assert
  implicit none

  !The horizon crossing approximation
  logical :: HC_approx=.false.

  !Tolerance value for checking division by zero
  real(dp), private, parameter :: div_tol=1e-20

  !Zero checking for regularization
  interface check_div_zero
  	module procedure check_div_zero_real
  	module procedure check_div_zero_array
  end interface

  contains

    !Evalutes adiabatic power spectrum under SR approximation at field positions
    !phi_end, given that the mode of interest crossed the horizon at phi_pivot
    !Assumes a massless, uncorrelated mode subhorizon
    function PR_SR(phi_pivot,phi_end, spectrum) result(PR)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) ::PR
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp) :: H_piv, V_piv, P_dphi
      real(dp), dimension(size(phi_pivot)) :: eps_i, u_i
      type(power_spectra), intent(in), optional :: spectrum
      integer :: ii, jj

      V_piv = pot(phi_pivot)
      H_piv = sqrt(V_piv/3.0e0_dp)
      dN = dNdphi_SR(phi_pivot,phi_end)

      !eps_i = eps_SR(phi_pivot)
      !u_i = (V_i_sum_sep(phi_pivot)+Z_i_BE(phi_end))/V_piv

      P_dphi = (H_piv/2.0e0_dp/pi)**2

      PR = sum(dN*dN)*P_dphi

      if (present(spectrum)) then
        PR=0e0_dp
        do ii=1,size(dN); do jj=1,size(dN)
          PR = PR+ &
            dN(ii)*dN(jj)*spectrum%phi_ij(ii,jj)
        end do; end do
      end if

      call assert%check(.not. isnan(PR), __FILE__, __LINE__)

    end function PR_SR

    function r_SR(phi_pivot,phi_end, spectrum) result(r)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: r
      real(dp), dimension(size(phi_pivot)) :: dN
      type(power_spectra), intent(in), optional :: spectrum
      real(dp) :: P_tens, P_scal, H, V

      integer :: ii, jj

      dN = dNdphi_SR(phi_pivot,phi_end)

      r = 8.0e0_dp/sum(dN*dN)

      if (present(spectrum)) then
        V = pot(phi_pivot)
        H = sqrt(V/3.0e0_dp)
        P_tens = 8.0e0_dp*(H/2.0e0_dp/pi)**2

        P_scal=0e0_dp
        do ii=1,size(dN); do jj=1,size(dN)
          P_scal = P_scal + &
            dN(ii)*dN(jj)*spectrum%phi_ij(ii,jj)
        end do; end do

        r = P_tens/P_scal

      end if

      call assert%check(.not. isnan(r), __FILE__, __LINE__)


    end function r_SR

    function nt_SR(phi_pivot) result(nt)
      real(dp), dimension(:), intent(in) :: phi_pivot
      real(dp) :: nt, eps(size(phi_pivot))

      eps = eps_SR(phi_pivot)

      nt = -2.0e0_dp*sum(eps)
      !nt = -2.0e0_dp*sum(eps)/(1.0e0_dp - sum(eps))

      call assert%check(.not. isnan(nt), __FILE__, __LINE__)

    end function nt_SR

    !Assumes independent, GRFs at horizon exit
    function fNL_SR(phi_pivot,phi_end) result(fnl)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: fnl
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2N
      integer :: ii, jj

      dN = dNdphi_SR(phi_pivot,phi_end)
      d2N = d2Ndphi2_SR(phi_pivot,phi_end)

      fnl=0e0_dp
      do ii=1,size(dN); do jj=1,size(dN)
        fnl = fnl + dN(ii)*dN(jj)*d2N(ii,jj)
      end do; end do
      fnl = fnl*(-5.0e0_dp/6.0e0_dp)/(sum(dN*dN))**2

      call assert%check(.not. isnan(fnl), __FILE__, __LINE__)


    end function fNL_SR

    !Assumes independent, GRFs at horizon exit
    !From Eq 41 in astro-ph/0611075
    function tauNL_SR(phi_pivot,phi_end) result(taunl)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: taunl
      real(dp), dimension(size(phi_pivot)) :: dN
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2N
      integer :: aa, bb, cc

      dN = dNdphi_SR(phi_pivot,phi_end)
      d2N = d2Ndphi2_SR(phi_pivot,phi_end)

      taunl=0e0_dp
      do aa=1,size(dN); do bb=1,size(dN); do cc=1,size(dN)
        taunl = taunl + d2N(aa,bb)*d2N(aa,cc)*dN(bb)*dN(cc)
      end do; end do; end do
      taunl = taunl*(1.0e0_dp)/(sum(dN*dN))**3

      call assert%check(.not. isnan(taunl), __FILE__, __LINE__)

    end function tauNL_SR

    function ns_SR(phi_pivot,phi_end) result(ns)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: ns, eps_piv, V
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2V
      real(dp), dimension(size(phi_pivot)) :: dV, dN
      integer :: ii, jj
      real(dp), dimension(size(phi_pivot)) :: eps_i, eta_i, u_i

      !dV = dVdphi(phi_pivot)
      d2V = d2Vdphi2(phi_pivot)
      dN = dNdphi_SR(phi_pivot,phi_end)
      !eps_i = eps_SR(phi_pivot)
      !eta_i = eta_SR(phi_pivot)
      eps_piv = sum(eps_SR(phi_pivot))
      V = pot(phi_pivot)

      !u_i = (V_i_sum_sep(phi_pivot)+Z_i_BE(phi_end))/V

      !ns = 1.0e0_dp - 2.0e0_dp*eps_piv &
      !  - (4.0e0_dp/sum(u_i**2/eps_i))*&
      !  (1.0e0_dp - sum(eta_i*u_i**2/2.0e0_dp/eps_i))

      ns = 1.0e0_dp &
        - 2.0e0_dp*eps_piv &
        - (2.0e0_dp/sum(dN*dN))

      do ii=1,size(phi_pivot); do jj=1,size(phi_pivot)
        ns = ns +&
          2.0e0_dp*d2V(ii,jj)*dN(ii)*dN(jj)/V/sum(dN*dN)
      end do; end do

      call assert%check(.not. isnan(ns), __FILE__, __LINE__)

    end function ns_SR

    !Running of ns
    !Formula as in Eq 6.14 1203.3792
    function alpha_s_SR(phi_pivot,phi_end) result(alpha_s)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp) :: alpha_s


      real(dp), dimension(size(phi_pivot)) :: V_i, Z_i, u_i
      real(dp) :: V, eps_piv
      real(dp), dimension(size(phi_pivot)) :: eps_i_piv, eta_i_piv, &
        xi_i_piv
      real(dp) :: sum_ui_over_epsi
      real(dp) :: term1, term2, term3, term4, term5, term6

      V = pot(phi_pivot)
      eps_i_piv = eps_SR(phi_pivot)
      eps_piv = sum(eps_i_piv)
      xi_i_piv = xi_SR(phi_pivot)
      eta_i_piv = eta_SR(phi_pivot)

      V_i = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      u_i = (V_i + Z_i)/V

      !Check for division by zero
      call check_div_zero(eps_i_piv)

      sum_ui_over_epsi = sum(u_i**2/eps_i_piv)


      term1 = -8.0e0_dp * eps_piv**2
      term2 = 4.0e0_dp * sum(eps_i_piv*eta_i_piv)

      term3 = (-16.0e0_dp/sum_ui_over_epsi**2)*&
        (1.0e0_dp - sum(eta_i_piv*u_i**2/2.0e0_dp/eps_i_piv))**2

      term4 = (-8.0e0_dp/sum_ui_over_epsi)*&
        sum(eta_i_piv*u_i*(1.0e0_dp - &
          eta_i_piv*u_i**2/2.0e0_dp/eps_i_piv))

      term5 = (4.0e0_dp*eps_piv/sum_ui_over_epsi)*&
        sum(eta_i_piv*u_i**2/eps_i_piv)

      term6 = (-2.0e0_dp/sum_ui_over_epsi)*&
        sum(xi_i_piv*u_i**2/eps_i_piv)

      alpha_s = term1 + term2 + term3 + term4 + term5 + term6

      call assert%check(.not. isnan(alpha_s), __FILE__, __LINE__)



    end function alpha_s_SR

    !Deriv of N wrt phi on horiz cross surface
    !Battefeld-Easther eq 29
    function dNdphi_SR(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_pivot)) :: dNdphi_SR
      real(dp), dimension(size(phi_pivot)) :: eps, V_i, Z_i
      real(dp) :: V

      integer :: ii

      eps = eps_SR(phi_pivot)
      V_i = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      V = pot(phi_pivot)

      !Check for division by zero
      call check_div_zero(eps)

      dNdphi_SR = ((1.0e0_dp/sqrt(2.0e0_dp*eps))/V)*(V_i + Z_i)

      call assert%check(.not. any(isnan(dNdphi_SR)), __FILE__, __LINE__)

    end function dNdphi_SR

    !2nd Deriv of N wrt phi on horiz cross surface
    !Battefeld-Easther eq 32
    function d2Ndphi2_SR(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_pivot),size(phi_pivot)) :: d2Ndphi2_SR
      real(dp), dimension(size(phi_end),size(phi_end)) :: delta

      real(dp), dimension(size(phi_end)) :: eta_piv, eps_piv, V_i_piv, Z_i
      real(dp) :: V_piv
      real(dp), dimension(size(phi_end),size(phi_end)) :: dZ_ij

      integer :: ll, kk

      eta_piv = eta_SR(phi_pivot)
      eps_piv = eps_SR(phi_pivot)
      V_i_piv = V_i_sum_sep(phi_pivot)
      Z_i = Z_i_BE(phi_end)
      V_piv = pot(phi_pivot)
      dZ_ij = dZdphi_ij_BE(phi_pivot,phi_end)

      !Identity matrix
      delta=0e0_dp
      do ll=1, size(phi_end)
        delta(ll,ll) = 1.0e0_dp
      end do

      !Check for division by zero
      call check_div_zero(eps_piv)
      call check_div_zero(V_piv)

      do ll=1, size(phi_end); do kk=1,size(phi_end)
        d2Ndphi2_SR(ll,kk) = delta(kk,ll)*&
            (1.0e0_dp - &
              (eta_piv(ll)/2.0e0_dp/eps_piv(ll))*&
              (V_i_piv(ll)+Z_i(ll))/V_piv&
             ) +&
             (1.0e0_dp/sqrt(2.0e0_dp*eps_piv(ll))/V_piv)*dZ_ij(ll,kk)
      end do; end do

      call assert%check(.not. any(isnan(d2Ndphi2_SR)), __FILE__, __LINE__)

    end function d2Ndphi2_SR

    function eps_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: eps_SR
      real(dp), dimension(size(phi)) :: dV
      real(dp) :: V

      dV = dVdphi(phi)
      V = pot(phi)
      eps_SR = 0.5e0_dp*( dV**2/V**2)

      call assert%check(.not. any(isnan(eps_SR)), __FILE__, __LINE__)

    end function eps_SR

    function eta_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: eta_SR
      real(dp), dimension(size(phi),size(phi)) :: d2V
      real(dp) :: V
      integer :: i

      d2V = d2Vdphi2(phi)
      V = pot(phi)

      do i=1, size(eta_SR)
        eta_SR(i) = d2V(i,i)**2/V
      end do

      call assert%check(.not. any(isnan(eta_SR)), __FILE__, __LINE__)

    end function eta_SR

    function xi_SR(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: xi_SR
      real(dp), dimension(size(phi),size(phi),size(phi)) :: d3V
      real(dp) :: V
      real(dp), dimension(size(phi)) :: dV
      integer :: i

      d3V = d3Vdphi3(phi)
      V = pot(phi)
      dV = dVdphi(phi)

      do i=1, size(xi_SR)
        xi_SR(i) = sqrt(dV(i)*d3V(i,i,i)**2/V**2)
      end do

      call assert%check(.not. any(isnan(xi_SR)), __FILE__, __LINE__)

    end function xi_SR

    !The function Z_i from Battefeld-Easther that encodes all details from the
    !end of inflation surface.  Eq. 31
    function Z_i_BE(phi_end)
      real(dp), dimension(:), intent(in) :: phi_end
      real(dp), dimension(size(phi_end)) :: Z_i_BE

      real(dp), dimension(size(phi_end)) :: eps_end
      real(dp) :: eps, V

      if (HC_approx) then
        Z_i_BE=0e0_dp
        return
      end if

      eps_end = eps_SR(phi_end)
      eps = sum(eps_end)
      V = pot(phi_end)

      !Check for division by zero
      call check_div_zero(eps)

      Z_i_BE = V*eps_end/eps - V_i_sum_sep(phi_end)

      call assert%check(.not. any(isnan(Z_i_BE)), __FILE__, __LINE__)

    end function Z_i_BE

    !Deriv of Z_i wrt fields at horiz cross
    !Battefeld-Easther Eq. 33
    function dZdphi_ij_BE(phi_pivot,phi_end)
      real(dp), dimension(:), intent(in) :: phi_pivot, phi_end
      real(dp), dimension(size(phi_end),size(phi_end)) :: dZdphi_ij_BE

      real(dp), dimension(size(phi_end)) :: eps_end, eps_piv, eta_end
      real(dp) :: eps_t_end, V_end, V_piv
      real(dp), dimension(size(phi_end),size(phi_end)) :: delta

      integer :: ii, jj, ll, kk

      if (HC_approx) then
        dZdphi_ij_BE=0e0_dp
        return
      end if

      eps_end = eps_SR(phi_end)
      eps_piv = eps_SR(phi_pivot)
      eps_t_end = sum(eps_end)
      V_end = pot(phi_end)
      V_piv = pot(phi_pivot)
      eta_end = eta_SR(phi_end)

      !Check for division by zero
      call check_div_zero(eps_piv)
      call check_div_zero(eps_t_end)
      call check_div_zero(V_piv)

      !Identity matrix
      delta=0e0_dp
      do ll=1, size(phi_end)
        delta(ll,ll) = 1.0e0_dp
      end do

      !Summation over jj
      dZdphi_ij_BE = 0e0_dp
      do ll=1, size(phi_end); do kk=1,size(phi_end)
        do jj=1,size(phi_end)
          dZdphi_ij_BE(ll,kk) = dZdphi_ij_BE(ll,kk) + &
            (-V_end**2/V_piv)*sqrt(2.0e0_dp/eps_piv(kk))*&
            (&
              eps_end(jj)*&
              ((eps_end(ll)/eps_t_end) - delta(ll,jj))*&
              ((eps_end(kk)/eps_t_end) - delta(kk,jj))*&
              (1.0e0_dp - (eta_end(jj)/eps_t_end))&
            )
        end do
      end do; end do

      call assert%check(.not. any(isnan(dZdphi_ij_BE)), __FILE__, __LINE__)


    end function dZdphi_ij_BE

    !For a sum-separable potential V=\sum_i V_i.  This returns only the V_i part
    !For sum-separable potentials with different potential for each field, have
    !to put this in by hand.
    function V_i_sum_sep(phi)
      real(dp), dimension(:), intent(in) :: phi
      real(dp), dimension(size(phi)) :: V_i_sum_sep

      real(dp), dimension(:,:), allocatable :: vparams_temp
      integer :: vrows, jj
      real(dp) :: V

      !Sum-separable, but not same potential
      !Really, really ugly...
      if (potential_choice == 21) then
        allocate(vparams_temp(size(vparams,1), size(vparams,2)))
        vparams_temp = vparams

        !First field
        vparams(1,2) = 0.e0_dp
        V_i_sum_sep(1) = pot(phi)
        vparams = vparams_temp

        !Second field
        vparams(2,1) = 0.e0_dp
        V_i_sum_sep(2) = pot(phi)

        vparams = vparams_temp
        return
      end if


      !The idea: make temp copy of vparams; change vparams as if it had only the
      !one field; get V; restore vparams
      !NB: vparams(vrows,num_inflaton)

      vrows = size(vparams,1)

      allocate(vparams_temp(size(vparams,1), size(vparams,2)))
      vparams_temp = vparams

      do jj=1,size(phi)
        deallocate(vparams)
        allocate(vparams(vrows,1))
        vparams(:,1) = vparams_temp(:,jj)
        V_i_sum_sep(jj) = pot((/phi(jj)/))
      end do

      deallocate(vparams)
      allocate(vparams(size(vparams_temp,1), size(vparams_temp,2)))
      vparams = vparams_temp
      deallocate(vparams_temp)

    end function V_i_sum_sep

    !Subroutines for regularizing a 0/0 error in the summations above
    pure subroutine check_div_zero_array(array)
      real(dp), dimension(:), intent(inout) :: array
      integer :: ii

      if (any(abs(array)<div_tol)) then
        do ii=1, size(array)
          if (abs(array(ii))<div_tol) then
            array(ii) = div_tol
          end if
        end do
      end if

    end subroutine check_div_zero_array

    pure subroutine check_div_zero_real(real_)
      real(dp), intent(inout) :: real_

      if (abs(real_)<div_tol) then
        real_ = div_tol
      end if

    end subroutine check_div_zero_real



end module modpk_deltaN
