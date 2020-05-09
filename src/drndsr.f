      real*8 function drndsr()
C  YOKOYA's subroutine for
C  RANDOM NUMBER GENERATOR OF SYNCHROTRON RADIATION SPECTRUM.
C  GENERATES X(1),X(2).....X(N).
C  X=(PHOTON ENERGY)/(CRITICAL ENERGY)
C   PROBABILITY THAT X BE IN THE INTERVAL (Y,Y+DY) IS
C     PP(Y)*DY= DY*( INTEGRAL (3/(5*PI))*K(5/3,Z) FROM Y TO INFINITY )
C  METHOD: USE INVERSE function Y(P) OF P(Y)
C                 P(Y)=INTEGRAL PP(Y)*DY FROM 0 TO Y.
C          GENERATE UNIFORM RANDOM DEVIATE 'P' IN (0,1) BY THE
C          CONGRUENCE METHOD AND COMPUTE Y(P).
C  ACCURACY: RELATIVE ERROR OF Y(P) IS LESS THAN 2.0*E-4
C     EVERYWHERE IN THE RANGE 0<P<1-2**(-31).
C       <Y>=(8/(15*SQRT(3))*1.0000136
C       <Y**2>=(11/27)*1.0000153
C  CPU-TIME: 4.9*N+5.6 (MICRO SEC) ON M200H.OPT(3) (1984.4.24)
      implicit none
      real*8 p,tran,p2,q,t
      P=tran()
      IF(P.GE.0.658D0) GOTO 100
      P2=P**2
      drndsr=(((0.41841680D0*P2+0.14175015D0)*P2+0.30528148D0)*P2
     %       +0.53520052D0)*P2*P
      return
  100 Q=1D0-P
      IF(Q.LE.0.0297d0) GOTO 200
      drndsr=((-0.32809490D0*Q+0.20653610D0)*Q+0.011920005D0)/
     %    (((Q+0.88765651D0)*Q+0.19270658D0)*Q+0.0033141692D0)
      return
  200 T=-LOG(Q)
      drndsr=(674.99117D0*T+148.32871D0)/
     %  ((T-225.45670D0)*T-692.23986D0)+T
      RETURN
      END

      real*8 function tdusr(anp,an)
c a random energy loss of SR using drndsr
c anp: expected number of photons
c returns: a random energy loss due to SR in unit of uc
c cutoff: threshold to ingnore multiple photons beyond this probability 
c npmax: the maximum number of photons to try
      implicit none
      real*8 anp,p,p0,pn,tran,drndsr,an,x
      real*8 , parameter :: cutoff=0.9999d0
      integer*4 np
      integer*4, parameter:: npmax=10000
      p0=exp(-anp)
      pn=p0
      p=tran()
      tdusr=0.d0
      an=0.d0
      x=1.d0
      do np=1,npmax
        if(p .le. pn)then
          return
        endif
        tdusr=tdusr+drndsr()
        an=an+1.d0
        x=x*anp/an
        pn=pn+x*p0
        if(pn .gt. cutoff)then
          exit
        endif
      enddo
      return
      end

      subroutine tdusrn(anp,dph,r1,r2,an)
c a random energy loss of SR using drndsr
c anp: expected number of photons
c returns: dph = sum of photon energies/uc, 
c  r1 = sum(r_n u_n)/sum(u_n), r2 = sum(r_n^2 u_n)/sum(r_n u_n)/2
c  an = number of actually emitted photons
c cutoff: threshold to ingnore multiple photons beyond this probability 
c npmax: the maximum number of photons to try
      implicit none
      real*8 , intent(in)::anp
      real*8 ,intent(out)::dph,an,r1,r2
      real*8 p,p0,pn,tran,drndsr,x,r,dph1
      real*8 , parameter :: cutoff=0.9999d0
      integer*4 , parameter :: npmax=10000
      integer*4 k
      p0=exp(-anp)
      pn=p0
      p=tran()
      r1=0.d0
      r2=0.d0
      an=0.d0
      x=1.d0
      dph=0.d0
      do k=1,npmax
        if(p .le. pn)then
          exit
        endif
        an=an+1.d0
        x=x*anp/an
        pn=pn+x*p0
        dph1=drndsr()
        dph=dph+dph1
        r=tran()
        r1=r1+r*dph1
        r2=r2+r**2*dph1
        if(pn .gt. cutoff)then
          exit
        endif
      enddo
      if(r1 .ne. 0.d0)then
        r2=r2/r1*.5d0
        r1=r1/dph
      endif
      return
      end

      subroutine tdusrnpl(anp,dph,r1,r2,an,dpr,rph)
c
c Same as tdusrn, with photonlist generation
c      
      implicit none
      integer*4 , parameter :: npmax=10000
      real*8 , intent(in)::anp
      real*8 ,intent(out)::dph,an,r1,r2,dpr(npmax),rph(npmax)
      real*8 p,p0,pn,tran,drndsr,x,r,dph1
      real*8 , parameter :: cutoff=0.9999d0
      integer*4 k
      p0=exp(-anp)
      pn=p0
      p=tran()
      r1=0.d0
      r2=0.d0
      an=0.d0
      x=1.d0
      dph=0.d0
      do k=1,npmax
        if(p .le. pn)then
          exit
        endif
        an=an+1.d0
        x=x*anp/an
        pn=pn+x*p0
        dph1=drndsr()
        dph=dph+dph1
        r=tran()
        r1=r1+r*dph1
        r2=r2+r**2*dph1
        dpr(k)=dph1
        rph(k)=r
        if(pn .gt. cutoff)then
          exit
        endif
      enddo
      if(r1 .ne. 0.d0)then
        r2=r2/r1*.5d0
        r1=r1/dph
      endif
      return
      end
