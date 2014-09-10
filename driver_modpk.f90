program driver_modpk
  use camb_interface
  use modpkparams
  use access_modpk, ONLY: potinit, evolve
  implicit none
  integer*4 :: i
  real(dp) :: kin, pow, powt, kmin, kmax, dlnk, pmin, pmax
  real(dp) :: pspivot, ptpivot
  real(dp) :: ps0, ps1, ps2, pt0, pt1, pt2

  !modpkoutput=.true.

  potential_choice = 2
  vparams_num = 2
  vnderivs=.false.

  slowroll_infl_end =.true.
  instreheat=.false.
  phi_infl_end=0.

  kmin = log(5.d-4)
  kmax = log(5.d0)

  vparams(1) = -2.2
  vparams(2) = 0.8

  phi_init0 = 20.
  k_pivot=0.05d0
  N_pivot=55.d0
  findiffdphi = epsilon(1.d0)

  call potinit()

  if(pk_bad==0) then

     call evolve(k_pivot,pspivot,ptpivot)

     do i=1,500
        kin=exp(kmin+real(i)*(kmax-kmin)/real(500-1))
        call evolve(kin,pow,powt)
        write(*,*) kin,pow,powt
     end do

     dlnk = 0.1d0
     call evolve(k_pivot,ps1,pt1)
     write(*,*)
     write(*,*) 'A_s =', ps1
     write(*,*) 'r =', pt1/ps1
     call evolve(k_pivot,ps0,pt0)
     call evolve(k_pivot*exp(-dlnk),ps1,pt1)
     call evolve(k_pivot*exp(dlnk),ps2,pt2)
     write(*,*) 'n_s =', 1.d0+log(ps2/ps1)/dlnk/2.d0
     write(*,*) 'n_t =', log(pt2/pt1)/dlnk/2.d0
     write(*,*) 'dn_s/dlnk =', log(ps1*ps2/ps0**2)/dlnk**2

  else
     write(*,*) 'pk_bad =', pk_bad
  endif

end program driver_modpk
