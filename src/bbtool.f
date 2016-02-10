C   09/02/93 305021512  MEMBER NAME  GAUINV   *.FORT     M  E2FORT
      FUNCTION GAUINV(P0)
C  INVERSE OF (INTEGRATED) NORMAL DISTRIBUTION FUNCTION
C              1         X= Y
C     P(Y)=-----------* INTEGRAL EXP(-X**2/2) DX
C          SQRT(2*PI)    X= -INF
C     IF P(Y)=P0, THEN GAUINV(P0)=Y.
C        0 < P0 < 1 ,   -INF < Y < +INF
C  IF THIS ROUTINE IS USED TO CONVERT UNIFORM RANDOM NUMBERS TO
C  GAUSSIAN, MAXIMUM RELATIVE ERROR IN THE DISTRIBUTION FUNCTION
C  DP/DX=EXP(-X**2/2)/SQRT(2*PI) IS LESS THAN 0.640E-3 EVERYWHERE
C  IN THE RANGE  2**(-31) < P0 < 1-2**31.  (MINIMAX APPROXIMATION)
C
      IMPLICIT REAL*8(A-H,O-Z)
C------------------------
      DATA PP1/0.334624883253D0/, QQ2/0.090230446775D0/,
     1     QQ3/0.049905685242D0/, QQ4/0.027852994157D0/,
     2     QQ5/0.015645650215D0/
      DATA A3/ 4.5585614D+01/, A2/ 2.1635544D+00/, A1/ 2.7724523D+00/,
     1     A0/ 2.5050240D+00/,
     2     B4/ 4.0314354D+02/, B3/-2.7713713D+02/, B2/ 7.9731883D+01/,
     3     B1/-1.4946512D+01/, B0/ 2.2157257D+00/,
     4     C4/ 4.1394487D+03/, C3/-1.5585873D+03/, C2/ 2.4648581D+02/,
     5     C1/-2.4719139D+01/, C0/ 2.4335936D+00/,
     6     D4/ 4.0895693D+04/, D3/-8.5400893D+03/, D2/ 7.4942805D+02/,
     7     D1/-4.1028898D+01/, D0/ 2.6346872D+00/,
     8     E4/ 3.9399134D+05/, E3/-4.6004775D+04/, E2/ 2.2566998D+03/,
     9     E1/-6.8317697D+01/, E0/ 2.8224654D+00/,
     O     F0/-8.1807613D-02/, F1/-2.8358733D+00/, F2/ 1.4902469D+00/
C------------------------
      GAUINV=0.D0
      P=P0-0.5D0
      P1=ABS(P)
      IF(P1.GE.PP1) GOTO 120
      P2=P**2
      GAUINV=(((A3*P2+A2)*P2+A1)*P2+A0)*P
      RETURN
  120 Q=0.5D0-P1
      IF(Q.LE.QQ2) GOTO 140
      GAUINV=(((B4*Q+B3)*Q+B2)*Q+B1)*Q+B0
      GOTO 200
  140 IF(Q.LE.QQ3) GOTO 150
      GAUINV=(((C4*Q+C3)*Q+C2)*Q+C1)*Q+C0
      GOTO 200
  150 IF(Q.LE.QQ4) GOTO 160
      GAUINV=(((D4*Q+D3)*Q+D2)*Q+D1)*Q+D0
      GOTO 200
  160 IF(Q.LE.QQ5) GOTO 170
      GAUINV=(((E4*Q+E3)*Q+E2)*Q+E1)*Q+E0
      GOTO 200
  170 IF(Q.LE.0D0) GOTO 900
      T=SQRT(-2D0*LOG(Q))
      GAUINV=T+F0+F1/(F2+T)
  200 IF(P.LT.0D0) GAUINV=-GAUINV
      RETURN
  900 WRITE(6,910) P0
  910 FORMAT(' (FUNC.GAUINV) INVALID INPUT ARGUMENT ',1PD20.13)
      RETURN
      END
C
C CERRF *************************************************
      SUBROUTINE CERRF(XX, YY, WX, WY)
*----------------------------------------------------------------------*
* Purpose:                                                             *
*   Modification of WWERF, double precision complex error function,    *
*   written at CERN by K. Koelbig.                                     *
* Input:                                                               *
*   XX, YY    (real)    Argument to CERF.                              *
* Output:                                                              *
*   WX, WY    (real)    Function result.                               *
*----------------------------------------------------------------------*

*---- Single precision version.
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (I-N)
C     REAL RWX, RWY
      PARAMETER         (MWFLT = 1, MREAL = 3)
*     IMPORTANT: MWNAM must be an integer multiple of MWFLT.
      PARAMETER         (MCWRD = 4)
      PARAMETER         (MCNAM = 16, MWNAM = MCNAM / MCWRD)

      PARAMETER         (CC     = 1.12837 91670 9551D0)
      PARAMETER         (ONE    = 1.D0)
      PARAMETER         (TWO    = 2.D0)
      PARAMETER         (XLIM   = 5.33D0)
      PARAMETER         (YLIM   = 4.29D0)
      DIMENSION         RX(33), RY(33)

c     XX=REAL(ZI)
c     YY=AIMAG(ZI)

      X = ABS(XX)
      Y = ABS(YY)

      IF (Y .LT. YLIM  .AND.  X .LT. XLIM) THEN
         Q = (ONE - Y / YLIM) * SQRT(ONE - (X/XLIM)**2)
         H = ONE / (3. 2* Q)
         NC =  7+ INT(23.0*Q)
         XL = H**( 1- NC)
         XH = Y + 0.5/H
         YH = X
         NU =  10+ INT(21.0*Q)
         RX(NU+1) = 0.
         RY(NU+1) = 0.

         DO 10 N = NU, 1, -1
            TX = XH + N * RX(N+1)
            TY = YH - N * RY(N+1)
            TN = TX*TX + TY*TY
            RX(N) = 0. 5* TX / TN
            RY(N) = 0. 5* TY / TN
   10    CONTINUE

         SX = 0.
         SY = 0.

         DO 20 N = NC, 1, -1
            SAUX = SX + XL
            SX = RX(N) * SAUX - RY(N) * SY
            SY = RX(N) * SY + RY(N) * SAUX
            XL = H * XL
   20    CONTINUE

         WX = CC * SX
         WY = CC * SY
      ELSE
         XH = Y
         YH = X
         RX(1) = 0.
         RY(1) = 0.

         DO 30 N = 9, 1, -1
            TX = XH + N * RX(1)
            TY = YH - N * RY(1)
            TN = TX*TX + TY*TY
            RX(1) = 0. 5* TX / TN
            RY(1) = 0. 5* TY / TN
   30    CONTINUE

         WX = CC * RX(1)
         WY = CC * RY(1)
      ENDIF

      IF(Y .EQ. 0.) WX = EXP(-X**2)
      IF(YY .LT. 0.) THEN
         WX = TWO * EXP(Y*Y-X*X) * COS(TWO*X*Y) - WX
         WY = - TWO * EXP(Y*Y-X*X) * SIN(TWO*X*Y) - WY
         IF(XX .GT. 0.) WY = -WY
      ELSE
         IF(XX .LT. 0.) WY = -WY
      ENDIF
*      print *,'wx= ',wx,' wy= ',wy
C     RWX = WX
C     RWY = WY
*      print *,'rwx= ',rwx,' rwy= ',rwy

      RETURN
      END
!
!
!
      real*8 FUNCTION bessi0(x)
      implicit none
      REAL*8 x
!Returns the modified Bessel function I0(x) for any real x. 
      REAL*8 ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     & q8,q9,y !Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     &     1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     &     0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1, 
     &     0.2635537d-1,-0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75) then 
         y=(x/3.75)**2
         bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))) 
      else
         ax=abs(x)
         y=3.75/ax
         bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     &        +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      END



      real*8 FUNCTION bessk0(x)
      implicit none
      REAL*8 x
! USESbessi0
! Returns the modified Bessel function K0(x) for positive real x. 
      REAL*8 bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y 
! Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     &     0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/ 
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     &     -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0) then ! Polynomial fit. 
         y=x*x/4.0
         bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+
     &        y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         y=(2.0/x)
         bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     &        y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      END



      real*8 FUNCTION bessi1(x)
      implicit none
      REAL*8 x
! Returns the modified Bessel function I1(x) for any real x. 
      REAL*8 ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
! Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, 
     &     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/ 
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, 
     &     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     &     -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75) then !Polynomial fit. 
         y=(x/3.75)**2
         bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         ax=abs(x)
         y=3.75/ax
         bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     &        y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
         if(x.lt.0.)bessi1=-bessi1
      endif
      return
      END


      real*8 FUNCTION bessk1(x) 
      implicit none
      REAL*8 x
! USES bessi1
! Returns the modified Bessel function K1(x) for positive real x. 
      REAL*8 bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y 
! Accumulate polynomials in double precision.
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0, 
     &     -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     &     0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0) then !Polynomial fit. 
         y=x*x/4.0
         bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+
     &        y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else 
         y=2.0/x
         bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     &        y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      END

      real*8 FUNCTION bessk(n,x)
      implicit none
      INTEGER n
      REAL*8 x
! USESbessk0,bessk1
! Returns the modified Bessel function Kn(x) for positive x and n ≥ 2. 
      INTEGER j
      REAL*8 bk,bkm,bkp,tox,bessk0,bessk1
      if(n.lt.2) then
         write(*,*) 'bad argument n in bessk'
         stop
      endif
      tox=2.0/x
      bkm=bessk0(x)
      bk=bessk1(x)
      do j=1,n-1
         bkp=bkm+j*tox*bk
         bkm=bk
         bk=bkp
      enddo
      bessk=bk
      return
      END


      real*8 FUNCTION bessi(n,x)
      implicit none
      INTEGER n,IACC
      REAL*8 x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.0e10,BIGNI=1.0e-10)
!C USESbessi0Returns the modified Bessel function 
!In(x) for any real x and n ≥ 2. 
      INTEGER j,m
      REAL*8 bi,bim,bip,tox,bessi0
      if (n.lt.2) then
         write(*,*) 'bad argument n in bessi'
         stop
      endif
      if (x.eq.0.) then
         bessi=0.
      else
         tox=2.0/abs(x)
         bip=0.0
         bi=1.0
         bessi=0.
         m=2*((n+int(sqrt(float(IACC*n))))) !Downward recurrence from even m.
         ! Make IACC larger to increase accuracy. 
         do j=m,1,-1 
            bim=bip+float(j)*tox*bi
            bip=bi
            bi=bim
            if (abs(bi).gt.BIGNO) then
               bessi=bessi*BIGNI
               bi=bi*BIGNI
               bip=bip*BIGNI
            endif
            if (j.eq.n) bessi=bip 
         enddo
         bessi=bessi*bessi0(x)/bi
         if (x.lt.0..and.mod(n,2).eq.1) bessi=-bessi
      endif
      return
      END



      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      implicit none
      INTEGER MAXIT
      REAL*8 ri,rip,rk,rkp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     &     PI=3.141592653589793d0)
!

      INTEGER i,l,nl 
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,
     &     fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,
     &     qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,
     &     rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) then
         write(*,'(A,1P,2E11.3,A)') 
     &        'bad arguments in bessik   ',x,xnu,'<0.'
         stop
      endif
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN) h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do i=1,MAXIT
         b=b+xi2
         d=1.d0/(b+d)
         c=b+1.d0/c
         del=c*d
         h=del*h
         if(abs(del-1.d0).lt.EPS)goto 1
      enddo 

      write(*,*) '***** x too large in bessik; try asymp. expan. *****'

 1    continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi 
      do l=nl,1,-1
         ritemp=fact*ril+ripl
         fact=fact-xi
         ripl=fact*ritemp+ril
         ril=ritemp
      enddo
      f=ripl/ril
      if(x.lt.XMIN) then
         x2=.5d0*x
         pimu=PI*xmu
         if(abs(pimu).lt.EPS)then
            fact=1.d0
         else
            fact=pimu/sin(pimu)
         endif
         d=-log(x2)
         e=xmu*d
         if(abs(e).lt.EPS)then
            fact2=1.d0
         else
            fact2=sinh(e)/e
         endif
!         write(*,*) xmu
         call beschb(xmu,gam1,gam2,gampl,gammi)
         ff=fact*(gam1*cosh(e)+gam2*fact2*d)
         sum=ff
         e=exp(e)
         p=0.5d0*e/gampl
         q=0.5d0/(e*gammi)
         c=1.d0
         d=x2*x2
         sum1=p
         do i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2)
            c=c*d/i
            p=p/(i-xmu)
            q=q/(i+xmu)
            del=c*ff
            sum=sum+del
            del1=c*(p-i*ff)
            sum1=sum1+del1
            if(abs(del).lt.abs(sum)*EPS) goto 2
         enddo
         write(*,*) '****** bessk series failed to converge ******'
 2       continue
         rkmu=sum
         rk1=sum1*xi2
      else
         b=2.d0*(1.d0+x)
         d=1.d0/b
         delh=d
         h=delh
         q1=0.d0
         q2=1.d0
         a1=.25d0-xmu2
         c=a1
         q=c
         a=-a1
         s=1.d0+q*delh
         do i=2,MAXIT
            a=a-2*(i-1)
            c=-a*c/i
            qnew=(q1-b*q2)/a
            q1=q2
            q2=qnew
            q=q+c*qnew
            b=b+2.d0
            d=1.d0/(b+a*d)
            delh=(b*d-1.d0)*delh
            h=h+delh
            dels=q*delh
            s=s+dels
            if(abs(dels/s).lt.EPS) goto 3
         enddo
         write(*,*) '***** bessik: failure to converge in cf2 *****'
 3       continue
         h=a1*h
         rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
         rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do i=1,nl
         rktemp=(xmu+i)*xi2*rk1+rkmu
         rkmu=rk1
         rk1=rktemp
      enddo
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END



      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
! USESchebev
! Evaluates NC1 and NC2 by Chebyshev expansion for |x| ≤ 1/2. 
! Also returns 1/NC(1 + x) and1/NC(1 − x). If converting to 
! double precision, set NUSE1 = 7, NUSE2 = 8. 
      REAL*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371168d0,6.5165112670737d-3,3.087090173086d-4,
     &     -3.4706269649d-6,6.9437664d-9,3.67795d-11,-1.356d-13/ 
      DATA c2/1.843740587300905d0,-7.68528408447867d-2,
     &     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,
     &     2.423096d-10,-1.702d-13,-1.49d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END


      real*8 FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 a,b,x,c(m)
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) 
     &     write(*,*) '**** x not in range in chebev ****'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do j=m,2,-1
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
      enddo 
      chebev=y*d-dd+0.5*c(1)
      return
      END

