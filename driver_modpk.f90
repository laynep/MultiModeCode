PROGRAM driver_modpk
  USE camb_interface
  USE modpkparams
  USE access_modpk, ONLY: potinit, evolve
  IMPLICIT NONE
  INTEGER*4 :: i
  DOUBLE PRECISION :: kin, pow, powt, kmin, kmax, dlnk, pmin, pmax
  DOUBLE PRECISION :: pspivot, ptpivot
  DOUBLE PRECISION :: ps0, ps1, ps2, pt0, pt1, pt2

  !modpkoutput=.true.

  potential_choice = 2
  vparams_num = 2
  vnderivs=.false.

  slowroll_infl_end =.true.
  instreheat=.false.
  phi_infl_end=0.
'
  kmin = log(5.d-4)
  kmax = log(5.d0)

  vparams(1) = -2.2
  vparams(2) = 0.8

  phi_init0 = 20.
  k_pivot=0.05d0
  N_pivot=55.d0
  findiffdphi = epsilon(1.d0)

  CALL potinit()

  IF(pk_bad==0) THEN 

     call evolve(k_pivot,pspivot,ptpivot)
     
     DO i=1,500
        kin=exp(kmin+REAL(i)*(kmax-kmin)/REAL(500-1))
        CALL evolve(kin,pow,powt)
        WRITE(*,*) kin,pow,powt
     END DO

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

  ELSE
     write(*,*) 'pk_bad =', pk_bad
  ENDIF

  STOP
END PROGRAM driver_modpk
