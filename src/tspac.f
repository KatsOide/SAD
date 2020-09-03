      subroutine tspac(np,x,px,y,py,z,g,dv,pz,al,al1)
      use tmacro
      implicit real*8 (a-h,o-z)
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),pz(np)
c      real*8 tr(4,4)
      complex*16 fxy
c      complex*16 fxygau
      complex*16 bbkick1
      sxx=0.d0
      sx=0.d0
      syy=0.d0
      sy=0.d0
      sxy=0.d0
      szz=0.d0
      sz=0.d0
      sp=0.d0
      alcs=al*.5d0-al1
      do i=1,np
        xc=x(i)+px(i)*alcs
        yc=y(i)+py(i)*alcs
        sx=sx+xc
        sxx=sxx+xc**2
        sy=sy+yc
        syy=syy+yc**2
        sxy=sxy+xc*yc
        sz=sz+z(i)
        szz=szz+z(i)**2
c        sp=sp+g(i)*(2.d0+g(i))
        sp=sp+g(i)
      enddo
      sx=sx/np
      sy=sy/np
      sz=sz/np
      sxx=sxx/np-sx**2
      syy=syy/np-sy**2
      sxy=sxy/np-sx*sy
      tilt=atan2(2.d0*sxy,sxx-syy)
      tx=cos(tilt)
      ty=sin(tilt)
      sigz=sqrt(szz/np-sz**2)
      sig1=.5d0*(sxx+syy)
      sig2=.5d0*hypot(sxx-syy,2.d0*sxy)
c      sig2=.5d0*sqrt((sxx-syy)**2+4.d0*sxy**2)
      sigu=sqrt(sig1+sig2)
      sigv=sqrt(abs(sig1-sig2))
      pb=p0*(1.d0+sp/np)
c      fs=2.d0*pbunch*rclassic/sigz/(sigu+sigv)*al/h0**3
      fs=-pbunch*rclassic/sigz*al/p0/(1.d0+pb**2)/sqrt(2.d0*pi)
      do i=1,np
c        p=(1.d0+g(i))**2
        p=1.d0+g(i)
        xc=x(i)+px(i)*alcs-sx
        yc=y(i)+py(i)*alcs-sy
        u= xc*tx+yc*ty
        v=-xc*ty+yc*tx
c        fxy=fxygau(u/sigu,v/sigu,1.d0,sigv/sigu)
        fxy=bbkick1(u,v,sigu,sigv)
c        call bbkick(dcmplx(sigu,sigv),
c     $       dcmplx(u,v),fxy,1,tr)
        fp=fs/p*exp(-.5d0*((z(i)-sz)/sigz)**2)
        dpx=fp*(dble(fxy)*tx-imag(fxy)*ty)
        dpy=fp*(dble(fxy)*ty+imag(fxy)*tx)
        x(i)=x(i)-dpx*alcs
        y(i)=y(i)-dpy*alcs
        px(i)=px(i)+dpx
        py(i)=py(i)+dpy
      enddo
c      write(*,*)u,v,dpx,dpy
      return
      end

      function tfgaussiancoulomb(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      complex*16 f
c      real*8 tr(4,4)
c      complex*16 fxygau
      complex*16 bbkick1
      if(isp .ne. isp1+4)then
        go to 9000
      endif
      do i=isp1+1,isp
        if(ktfnonrealq(ktastk(i)))then
          go to 9000
        endif
      enddo
c      f=fxygau(
c     $     rlist(isp1+1),
c     $     rlist(isp1+2),
c     $     rlist(isp-1),
c     $     rlist(isp))
c      call bbkick(dcmplx(vstk(ivstkoffset+isp-1),
c     $     vstk(ivstkoffset+isp)),
c     $     dcmplx(vstk(ivstkoffset+isp1+1),
c     $     vstk(ivstkoffset+isp1+2)),f,1,tr)
      f=bbkick1(rtastk(isp1+1),rtastk(isp1+2),
     $     rtastk(isp-1),rtastk(isp))
      kx=kxcalocv(-1,dble(f),imag(f))
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"x, y, sigx, sigy"')
      kx=dxnullo
      return
      end

      function tfgaussiancoulombu(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      real*8 twspu
      if(isp .ne. isp1+4)then
        go to 9000
      endif
      do i=isp1+1,isp
        if(ktfnonrealq(dtastk(i)))then
          go to 9000
        endif
      enddo
      kx%x(1)=twspu(rtastk(isp1+1),rtastk(isp1+2),
     $     rtastk(isp-1),rtastk(isp),1.d-10,1.d-13)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"x, y, sigx, sigy"')
      kx=dxnullo
      return
      end

      subroutine tfgaussiancoulombfitted(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: klr
      integer*4 isp1,irtc,i,itfmessage
      real*8 fx,fy,fu,fxx,fyy,fxy
      if(isp .ne. isp1+4)then
        go to 9000
      endif
      do i=isp1+1,isp
        if(ktfnonrealq(dtastk(i)))then
          go to 9000
        endif
      enddo
      call twspfu(rtastk(isp1+1),rtastk(isp1+2),
     $     rtastk(isp-1),rtastk(isp),fx,fy,fu,fxx,fyy,fxy)
      kx=kxavaloc(-1,6,klr)
      klr%rbody(1)=fx
      klr%rbody(2)=fy
      klr%rbody(3)=fu
      klr%rbody(4)=fxx
      klr%rbody(5)=fyy
      klr%rbody(6)=fxy
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"x, y, sigx, sigy"')
      return
      end

C   14/04/87 304210516  MEMBER NAME  FXYGAU   *.FORT     M  FORTRAN77
C******************** FXYGAU *******************************************
C   Yokoya's funtion to calculate force of gaussian beam.
c 
      COMPLEX*16 function FXYGAU(X,Y,SX,SY)
C   ELECTROMAGNETIC FORCE OF GAUSSIAN BEAM
C   NOT USED IN TRACKING. FOR COMPARISON ONLY.
      use macmath
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CC1,CC2
      PARAMETER (CC1 = m_sqrt2 * m_2_sqrtpi, CC2 = m_1_sqrt2)
C             CC1=SQRT(8/PI)=SQRT(2) * (2 / SQRT(PI))  CC2=SQRT(1/2)
      COMPLEX*16 Z,DFXY
      COMPLEX*16 CWERF
      IF(ABS(SX-SY).GE.1D-5*SX) GOTO 200
      EX=0.5d0*(X**2+Y**2)/SX**2
      IF(EX.LE.1D-5) GOTO 100
      Z=DCMPLX(X,-Y)
      FAC=CC1*SX
      IF(EX.LE.40.d0) FAC=FAC*(1.d0-EXP(-EX))
      FXYGAU=FAC/Z
      RETURN
  100 FAC=CC1/(2.d0*SX)*(1.d0-0.5d0*EX)
      FXYGAU=DCMPLX(FAC*X,FAC*Y)
      RETURN
  200 XA=ABS(X)
      YA=ABS(Y)
      D1=1.d0/SQRT(SX**2-SY**2)
      D2=CC2*D1
      Z=D2*DCMPLX(XA,YA)
      DFXY=CWERF(Z)
      EX=0.5d0*((X/SX)**2+(Y/SY)**2)
      IF(EX.GE.40.d0) GOTO 300
      R=SX/SY
      Z=DCMPLX(D2*XA/R,D2*YA*R)
      DFXY=DFXY-EXP(-EX)*CWERF(Z)
  300 DFXY=D1*(SX+SY)*DFXY
      IF(X.LT.0.d0) DFXY=CONJG(DFXY)
      FXYGAU=DCMPLX(IMAG(DFXY),DREAL(DFXY))
      IF(Y.LT.0.d0) FXYGAU=CONJG(FXYGAU)
      RETURN
      END
C   14/04/87 304210518  MEMBER NAME  CWERF    *.FORT     M  FORTRAN77
C******************** CWERF ****************************************
      COMPLEX*16 function CWERF(Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z
      LOGICAL B
      DATA CONST/1.12837 91670 9551D0/
      H2=0.d0
      XX=dble(Z)
      YY=imag(Z)
      X=ABS(XX)
      Y=ABS(YY)
      IF(Y .LT. 4.29D0 .AND. X .LT. 5.33D0) GO TO 1
      H=0.
      NC=0
      NU=8
      ALAMDA=0.
      B=.TRUE.
      GO TO 2
    1 S=(1D0-Y/4.29D0)*SQRT(1D0-X**2/28.41D0)
      H=1.6D0*S
      H2=2D0*H
      NC=6+INT(23D0*S)
      NU=9+INT(21D0*S)
      ALAMDA=H2**NC
      B= ALAMDA .EQ. 0D0
    2 R1=0D0
      R2=0D0
      S1=0D0
      S2=0D0
      N=NU+1
    3 N=N-1
      FN=N+1
      T1=Y+H+FN*R1
      T2=X-FN*R2
      C=0.5D0/(T1**2+T2**2)
      R1=C*T1
      R2=C*T2
      IF(H .LE. 0.0 .OR. N .GT. NC) GO TO 4
      T1=ALAMDA+S1
      S1=R1*T1-R2*S2
      S2=R2*T1+R1*S2
      ALAMDA=ALAMDA/H2
    4 IF(N .GT. 0) GO TO 3
      IF(B) GO TO 6
      RS1=S1
      RS2=S2
      GO TO 7
    6 RS1=R1
      RS2=R2
    7 RS1=CONST*RS1
      IF(Y .EQ. 0D0) RS1=EXP(-X**2)
      CWERF=DCMPLX(RS1,CONST*RS2)
      IF(YY .LT. 0D0) GO TO 8
      IF(XX .LT. 0D0) CWERF=CONJG(CWERF)
      RETURN
    8 CWERF=2D0*EXP(-DCMPLX(X,Y)**2)-CWERF
      IF(XX .GT. 0D0) CWERF=CONJG(CWERF)
      RETURN
      END
C******************** BBKICK *******************************************00016300
      SUBROUTINE BBKICK(S,D,P,K,T)
C  BEAM-BEAM KICK ANGLE AND TRANSFER MATRIX FOR 2-DIM GAUSSIAN BUNCH.
C  INPUT:
C    S(1),S(2)= HORIZONTAL AND VERTICAL RMS BEAM SIZEOF DRIVING BUNCH.
C    D(1),D(2)= HOR. AND VER. POSITION OF KICKED PARTICLE W.R.T.
C               THE DRIVING BUNCH. (SAME UNIT AS S'S)
C    K          IF K=1, KICK ANGLE ONLY. IF K=2, KICK ANGLE AND
C               TRANSFER MATRIX
C  OUTPUT:
C    P(1),P(2)= HOR. AND VER. KICK ANGLE. (MUST BE MULTIPLIED BY
C               A=(# OF PARTICLES IN THE DRIVING BUNCH)*
C                 (CLASSICAL RADIUS)/(BEAM ENERGY/REST ENERGY)
C               FOR SAME SIGN CHARGE A IS NEGATIVE.
C    T(4,4)     4*4 TRANSFER MATRIX NEAR (D(1),D(2)) IS GIVEN BY
C                M=I+A*T   (I=4*4 UNIT MATRIX. A=DEFINED ABOVE)
C               THE VECTOR IS (X,X',Y,Y').
C
      use macmath
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RTPI
      PARAMETER (RTPI = m_sqrtpi)
      COMPLEX*16 CWERF,ddf,ddz1,ddz2,ddf2
      DIMENSION S(2),D(2),P(2),T(4,4)
      DIMENSION FRI(2)
      EQUIVALENCE (DDF,FRI)
C
      IF(K.EQ.2) THEN
        T=0.d0
c        DO 100 J=1,4
c          DO 100 I=1,4
c 100        T(I,J)=0D0
      ENDIF
      EX=0.5*((D(1)/S(1))**2+(D(2)/S(2))**2)
      IF(EX.LE.1D-12) GOTO 300
      G=0D0
      IF(EX.LE.100D0) G=EXP(-EX)
      IF(ABS(S(1)-S(2))/S(1).LE.1D-5) GOTO 200
      I=1
      IF(S(1).LT.S(2)) I=2
      J=3-I
      R1=S(J)/S(I)
      R2=S(I)/S(J)
      C0=SQRT(2D0*(S(I)**2-S(J)**2))
      C1=2D0*RTPI/C0
      C2=2D0*C1/C0
      DDZ1=DCMPLX(ABS(D(I)),ABS(D(J)))/C0
      DDZ2=DCMPLX(R1*ABS(D(I)),R2*ABS(D(J)))/C0
      DDF=CWERF(DDZ1)-G*CWERF(DDZ2)
      IF(D(I).EQ.0D0) FRI(2)=0D0
      IF(D(J).EQ.0D0) FRI(1)=0D0
      P(I)=-C1*IMAG(DDF)
      P(J)=-C1*DREAL(DDF)
      IF(D(I).LT.0D0) P(I)=-P(I)
      IF(D(J).LT.0D0) P(J)=-P(J)
      IF(K.EQ.1) RETURN
      I2=I*2  
      I1=I2-1
      J2=J*2
      J1=J2-1
      T(I2,I1)=C2*(IMAG(DDZ1*DDF)-(1D0-R1*G)/RTPI)
      T(J2,J1)=-C2*(IMAG(DDZ1*DDF)-(1D0-R2*G)/RTPI)
      T(I2,J1)=C2*DREAL(DDZ1*DDF)
      IF(D(I)*D(J).LT.0D0) T(I2,J1)=-T(I2,J1)
      T(J2,I1)=T(I2,J1)
      RETURN
  200 SS=0.5D0*(S(1)+S(2))
      DDF=1D0/DCMPLX(D(1),-D(2))
      C1=2D0*(1D0-G)
      P(1)=-C1*DREAL(DDF)
      P(2)=-C1*IMAG(DDF)
      IF(K.EQ.1) RETURN
      DDF2=DDF**2
      C2=2D0*G/SS**2
      T(2,1)=+C1*DREAL(DDF2)-C2*D(1)*DREAL(DDF)
      T(4,3)=-C1*DREAL(DDF2)-C2*D(2)*IMAG(DDF)
      T(2,3)=C1*IMAG(DDF2)-C2*D(1)*D(2)/(D(1)**2+D(2)**2)
      T(4,1)=T(2,3)
      RETURN
  300 T21=-2D0/(S(1)*(S(1)+S(2)))
      T43=-2D0/(S(2)*(S(1)+S(2)))
      P(1)=T21*D(1)
      P(2)=T43*D(2)
      IF(K.EQ.1) RETURN
      T(2,1)=T21
      T(4,3)=T43
      RETURN
      END

C******************** BBKICK1 *******************************************00016300
      complex*16 function bbkick1(dx,dy,sx,sy)
c Modified from Yokoya's
C  BEAM-BEAM KICK ANGLE AND TRANSFER MATRIX FOR 2-DIM GAUSSIAN BUNCH.
C  INPUT:
C    S(1),S(2)= HORIZONTAL AND VERTICAL RMS BEAM SIZEOF DRIVING BUNCH.
C    D(1),D(2)= HOR. AND VER. POSITION OF KICKED PARTICLE W.R.T.
C               THE DRIVING BUNCH. (SAME UNIT AS S'S)
C    K          IF K=1, KICK ANGLE ONLY. IF K=2, KICK ANGLE AND
C               TRANSFER MATRIX
C  OUTPUT:
C    P(1),P(2)= HOR. AND VER. KICK ANGLE. (MUST BE MULTIPLIED BY
C               A=(# OF PARTICLES IN THE DRIVING BUNCH)*
C                 (CLASSICAL RADIUS)/(BEAM ENERGY/REST ENERGY)
C               FOR SAME SIGN CHARGE A IS NEGATIVE.
C    T(4,4)     4*4 TRANSFER MATRIX NEAR (D(1),D(2)) IS GIVEN BY
C                M=I+A*T   (I=4*4 UNIT MATRIX. A=DEFINED ABOVE)
C               THE VECTOR IS (X,X',Y,Y').
C
      use macmath
      implicit none
      REAL*8 RTPI
      PARAMETER (RTPI = m_sqrtpi)
      real*8 ex,g,r1,r2,c0,c1,ss,t21,t43
      real*8 dx,dy,sx,sy,px,py
      COMPLEX*16 CWERF,ddf,ddz1,ddz2,bbkick1rec
      real*8 FRI(2)
      EQUIVALENCE (DDF,FRI)
C
      if(sx .lt. sy)then
        bbkick1=dcmplx(0.d0,1.d0)*conjg(bbkick1rec(dy,dx,sy,sx))
        return
      endif
      EX=0.5*((dx/sx)**2+(dy/sy)**2)
      IF(EX.LE.1D-12) GOTO 300
      G=0D0
      IF(EX.LE.100D0) G=EXP(-EX)
      IF(ABS(sx-sy)/sx.LE.1D-5) GOTO 200
      R1=Sy/Sx
      R2=Sx/Sy
      C0=SQRT(2D0*(Sx**2-Sy**2))
      C1=2D0*RTPI/C0
      DDZ1=DCMPLX(ABS(Dx),ABS(Dy))/C0
      DDZ2=DCMPLX(R1*ABS(Dx),R2*ABS(Dy))/C0
      DDF=CWERF(DDZ1)-G*CWERF(DDZ2)
      IF(Dx.EQ.0D0) FRI(2)=0D0
      IF(Dy.EQ.0D0) FRI(1)=0D0
      IF(Dx.EQ.0D0) FRI(2)=0D0
      IF(Dy.EQ.0D0) FRI(1)=0D0
      Px=-C1*IMAG(DDF)
      Py=-C1*DREAL(DDF)
      IF(Dx.LT.0D0) Px=-Px
      IF(Dy.LT.0D0) Py=-Py
      bbkick1=dcmplx(px,py)
      RETURN
  200 SS=0.5D0*(Sx+Sy)
      DDF=1D0/DCMPLX(Dx,-Dy)
      C1=2D0*(1D0-G)
      Px=-C1*DREAL(DDF)
      Py=-C1*IMAG(DDF)
      bbkick1=dcmplx(px,py)
      RETURN
  300 T21=-2D0/(Sx*(Sx+Sy))
      T43=-2D0/(Sy*(Sx+Sy))
      Px=T21*Dx
      Py=T43*Dy
      bbkick1=dcmplx(px,py)
      RETURN
      END

      complex*16 function bbkick1rec(dx,dy,sx,sy)
      implicit none
      real*8 dx,dy,sx,sy
      COMPLEX*16 bbkick1
      bbkick1rec=bbkick1(dx,dy,sx,sy)
      return
      end
