



      External fcn,bsstep,qtrap,pdf
      include 'header'

      integer argnum

      Integer errflg,loop1,loop2
      REAL*8 eps,h1,hmin,hmax,x1,x2,start,finish,deltaM,delta,initial,H 
      REAL*8 v(NMAX)
      REAL*8 potential,kinetic,pdf,work,bisect,partial,test,spread


      beta = 0.5d0
      mass4 = 1d0
      nval = 500


      if(beta.ne.0d0) then ! beta = 0 all masses identical


        
                                ! set mass terms
         msqd(1) = (1d0-sqrt(beta))**2


         deltaM = ((1d0+sqrt(beta))**2 - (1d0-sqrt(beta))**2)/
     &        real(nval-1)
         
         delta = 1d0/real(nval-1)
         

         do loop1 = 2,nval-1
            start = msqd(loop1-1)
            finish = bisect(start,delta,start,start+deltaM)
            msqd(loop1) = finish
         enddo
         
         msqd(nval) = (1d0+sqrt(beta))**2
        
      
         do loop1 =1,nval
            msqd(loop1) = mass4**2 * msqd(loop1)
         enddo
 

      else         
         do loop1=1,nval
            msqd(loop1) = mass4**2
         enddo
      endif


      open(31,file = 'masses.dat')  ! write masses to a file.  
      
      do loop1=1,NVAL
         write(31,*) msqd(loop1)
      enddo
      close(31)

      end





c----------------------------------------------------------------------------------------------
c
c
c
c
      real*8 function partial(x,y)
      external qtrap,pdf

      include 'header'


      real*8 x,y,pdf,s


      call  qtrap(pdf,x,y,s)
      partial = s
      return
      end

c---------------------------------------------------------------------------------------------
c     marchenko-pastur distribution, with modification to exclude delta
c     -fn at origin. (Multiplied by beta)

      real*8 function pdf(x)

      include 'header'

      real*8 x,sb
      
      sb = sqrt(beta)


      if(x.lt.(1d0-sb)**2) then
         pdf = 0d0
         return
      endif

      if(x.gt.(1d0+sb)**2) then
         pdf = 0d0
         return
      endif
      pdf = sqrt((x-(1d0-sb)**2)*((1d0+sb)**2-x))/(2d0 *pi* x)/beta
      return
      end


cc----------------------------------------------------------------------------------------------
c
c     a quick and dirty bisection rountine.
c

      real*8 function bisect(base,target,lostart,histart)
      
      include 'header'
      real*8 hi,lo,base,lostart,histart,target

      real*8 midpt,midval,hival,loval,partial



      lo = lostart
      hi = histart

      hival = partial(base,hi)-target

      loval = partial(base,lo)-target


      do while(hival.lt.0d0)
         hi = hi + (hi-lo)
         hival = partial(base,hi) -target
      enddo


      do while( (hi-lo).gt.10d-10)
         midpt = lo+(hi-lo)/2d0
         midval = partial(base,midpt)-target
         if(loval*midval.ge.0d0) then ! loval and midval have same sign
            loval = midval
            lo =midpt
         else
            hival=midval
            hi = midpt
         endif
      enddo
      
      bisect = lo + (hi-lo)/2d0
      return
      end



      SUBROUTINE qtrap(func,a,b,s)
      INTEGER JMAX
      DOUBLE PRECISION a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-8, JMAX=200)
CU    USES trapzd
      INTEGER j
      DOUBLE PRECISION olds
      olds=-1.d30


      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        if (s.eq.0.d0.and.olds.eq.0.d0.and.j.gt.6) return
        olds=s
11    continue
      write(6,*) 'too many steps in qtrap'
      stop
      END


      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x



      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else

        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
