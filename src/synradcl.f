      SUBROUTINE SYNRADCL(E0,B0,DS,FLAG,NG,EG,thx,THY,XI3,XI2,IRTN)
C
c  1/9/2016 Modified for different CHARGE.
c  1/6/2003 Received from K. Yokoya
c  1/7/2003 modified for f77, replace RANDCAIN with tran
c  suppress no radiation return
c  1/9/2003 restored thx for consistency with recoil
c
C  Generator of classical synchrotron radiation (electron/positron)
C   Needs random number routine RANDCAIN (uniform in (0,1)).
C
C Input parameter (kept intact)
C   E0:   Initial electron energy (eV)
C   B0:   Transverse magnetic field (Tesla)
C   DS:   Interaction length (m)
C   FLAG  0:  photon energy only
C         1:  photon energy and vertical angle
C         2:  photon energy, vertical angle and linear/circular polarization
C Output parameter
C   NG    0:  no radiation
C         1:  one photon emitted
C   EG    Photon energy (eV). 
C         Never exceeds E0 (classical formula slightly modified.)
C   THY   Vertical angle (radian)  (Not used if FLAG=0)
c   thx   horizontal angle (radian)  (Not used if FLAG=0)
C   XI3   Linear polarization (Not used if FLAG<=1)
C         (1+XI3)/2 : Probability to be polarized in the orbit plane.
C         (1-XI3)/2 : Probability to be perpendicular to the orbit plane.
C   XI2   Circular polarization  (Not used if FLAG<=1)
C         (1+XI2)/2 : Probability to be right-handed.
C         (1-XI2)/2 : Probability to be left-handed.
C         The coordinate axes are defined as
C           z:  direction of electron motion at the moment of radiation
C           x:  normal to the orbit in the orbit plane (>0 outward)
C           y:  perpendicular to the orbit plane ((x,y,z): right-handed)
C         XI3 and XI2 are the Stokes parameters based on this (x,y) axis.
C         XI1 is zero.
C   IRTN    0: normal
C           1: emission probability > 0.1  (warning only)
C           2: Upsilon > 0.05 (warning only)
C         101: emission probability > 1  (no calculation)
C         102: Upsilon > 1  (no calculation)
C 
C Caution
C  (1)  The emission probability is 
C          p=5/(2*sqrt(3))*alpha*gamma*DS/rho = 6 * B (Tesla) * DS(m)
C       This must be less than 1 (This routine generates at most one photon.)
C  (2)  Even when p<1, there is a statistical problem if p is not small enough.
C       The probability of emitting two photons during DS is p^2/2, which is
C       ignored here. Although the average number of photons and average
C       energy loss is correct with any p<1, the square average may differ.
C       It is safe to choose p<0.1.
C  (3)  The x-angle of photon is not computed. It is meaningless because
C       the orbit deflection by the field is (assumed to be) much larger
C       than 1/gamma. If one needs the x-angle, one should choose a point
C       in DS randomly where the photon is emitted and take electron orbit
C       angle at that point as the x-angle of the photon.
C  (4)  MASS and COMPTON have to be changed for proton radiation.
C    
      use tmacro
      IMPLICIT NONE
c     Including physical constant
c     CVELOC:	Speed of light
c     FINEST:	Fine-structure constant
c     PLANKR:	Dirac's constant(Planck's constant over 2Pi)
c     ELEMCH:	Elementary charge
c     ELMASS:	Electron mass energy equivalent in eV
c     Including math constant
c     M_SQRT3:	SQRT(3.D0)
      real*8 MASS,COMPTON,C1
c      parameter (MASS=elmass)
c      parameter (COMPTON=PLANKR*CVELOC/(ELEMCH*MASS))
c      parameter (C1=5.D0/(2.D0*M_SQRT3)*FINEST)
      INTEGER FLAG,NG,IRTN
      real*8 E0,B0,DS,EG,THY,XI3,XI2
      real*8 B,GAM,RHO,UPSILON,P00,P1,U,V,X,TH,PHI,
     %   THYG,THYG2,Z,R
      real*8 tran,thx,SRFUNCU,SRFUNCY,BK1323

      MASS=amass
      COMPTON=rclassic/rcratio
      C1=5.D0/(2.D0*M_SQRT3)*rcratio
      NG=0
      IRTN=0
      IF(E0.LE.0.OR.B0.EQ.0) RETURN
      B=ABS(B0)
      GAM=E0/MASS
c      RHO=E0/(CVELOC*B)
      RHO=brho/B
      UPSILON=COMPTON*GAM**2/RHO
      IF(UPSILON.GE.1D0) THEN
        IRTN=102
        RETURN
      ENDIF
      P00=C1*GAM/RHO*DS
      IF(P00.GE.1D0) THEN
        IRTN=101
        RETURN
      ENDIF
c     No-radiation return is suppressed, should be handled outside
c      P1=tran()
c      IF(P1.GE.P00) RETURN
c
      P1=tran()
      IF(P1.EQ.0) RETURN
      IF(FLAG.LE.0) THEN
        X=1.5D0*UPSILON*SRFUNCY(P1)
        EG=X/(1+X)*E0
      ELSE
        U=SRFUNCU(P1)
        IF(U.EQ.0) RETURN
        P1=tran()
        V=2*SQRT(2D0)*SIN(ASIN(5*P1/(4*SQRT(2D0)))/3)
        X=1.5D0*UPSILON*U*V**3
        EG=X/(1+X)*E0
        PHI=2*PI*tran()
        TH=SQRT(MAX(0D0,1-V**2))/(V*GAM)
        THY=TH*SIN(PHI)
C          TH*COS(PHI) can be x-angle
        thx=th*cos(phi)
        IF(FLAG.GE.2) THEN
          THYG=GAM*THY
          THYG2=SQRT(1+THYG**2)
          Z=X/(3*UPSILON)*THYG2**3
          R=-THYG/THYG2*BK1323(Z)
          XI2=2*R/(1+R**2)
          XI3=(1-R**2)/(1+R**2)
        ENDIF
      ENDIF
      NG=1
      IF(UPSILON.GE.0.05) IRTN=2
      IF(P00.GE.0.1D0) IRTN=1
         RETURN
      END

      FUNCTION SRFUNCU(P)
C   Inverse function of
C         p = integral  x*K(1/3,x) dx over 0<x<u   --->  u=SRFUNCU(p)
C ( 0<u<infty, 0<p<1)
C   Relative error of SRFUNCU is < 6*10^(-5).
C   Return 0 if P<=0 or P>=1
      IMPLICIT NONE
      real*8 SRFUNCU,P
      real*8 PC1,pc2
      data pc1 /0.3190D0/,PC2/0.7495D0/
      real*8 C00
      data c00 /0.84628437532D0/
C         C00=3/sqrt(4*pi)
      real*8 C1(6)
      data c1 /
     % 1.0202936D+00, 4.3043378D-01, 2.3886659D-01, 9.0502838D-01,
     %-1.4086504D+00, 1.7219123D+00/
      real*8 C2(6)
      data c2/
     %-9.9376219D-02, 4.6221963D+00,-1.1694960D+01, 2.6245046D+01,
     %-2.8102206D+01, 1.3195570D+01/
      real*8 C3(5)
      data c3 /
     % 2.8020354D-02, 9.2900218D-01,-1.4673954D+00, 1.6974039D+00,
     %-8.3557103D-01/
      INTEGER I
      real*8 P15,P25,F,CP,CPI

      F=0
      IF(P.LE.PC1) THEN
        IF(P.GT.0) THEN
          P15=P**0.2D0
          P25=P15**2
          F=0
          DO I=6,1,-1
            F=(F+C1(I))*P25
          ENDDO
          F=F*P15
        ENDIF
      ELSEIF(P.LE.PC2) THEN
        F=0
        DO I=6,1,-1
          F=F*P+C2(I)
        ENDDO
      ELSEIF(P.LE.0.999999999D0) THEN
        CP=-LOG(C00*(1-P))
        CPI=1/CP
        F=0
        DO I=5,1,-1
          F=F*CPI+C3(I)
        ENDDO
        F=F+CP+0.5D0*LOG(CP)
      ENDIF
      SRFUNCU=F
      RETURN
      END

      FUNCTION SRFUNCY(P)
C  Convert uniform random number (0,1) into photon energy.
C   Photon energy = SRFUNCY * critical_energy
C  y=SRFUNC(p) is the inverse of
C      p = \int_0^y dx \int_x^\infty [9\sqrt(3)/(8\pi)] *  K(5/3,z) dz
C  Return 0 if P<=0 or P>=1
      IMPLICIT NONE
      real*8 SRFUNCY,P
      real*8 PC1,pc2
      data pc1 /0.658D0/,PC2/0.9703/
      real*8 ya1,ya2,ya3,ya4,yb0,yb1,yb2,yc0,yc1,yc2,
     $     yd0,yd1,ye0,ye1
      data 
     %  YA1/0.53520052D0/, YA2/0.30528148D0/, YA3/0.14175015D0/,
     %  YA4/0.41841680D0/,
     %  YB0/0.011920005D0/,YB1/0.20653610D0/, YB2/-0.32809490D0/,
     %  YC0/0.33141692D-2/,YC1/0.19270658D0/, YC2/0.88765651D0/,
     %  YD0/148.32871D0/,  YD1/674.99117D0/,
     %  YE0/-692.23986D0/, YE1/-225.45670D0/
      real*8 Y,P2,Q,T
      
      Y=0
      IF(P.LE.PC1) THEN
        IF(P.GT.0) THEN
          P2=P**2
          Y=(((YA4*P2+YA3)*P2+YA2)*P2+YA1)*P2*P
        ENDIF
      ELSEIF(P.LE.PC2) THEN
        Q=1-P
        Y=((YB2*Q+YB1)*Q+YB0)/(((Q+YC2)*Q+YC1)*Q+YC0)
      ELSEIF(P.LT.1D0) THEN
        T=-LOG(1-P)
        Y=T+(YD1*T+YD0)/((T+YE1)*T+YE0)
      ENDIF
      SRFUNCY=Y
      RETURN
      END

      real*8 function BK1323(X)
C   BesselK(1/3,x) / BesselK(2/3,x) (X>=0)
C   Relative error does not exceed 3.83e-5
C   Return 0 if X<0.
      IMPLICIT NONE
      real*8 X
      real*8 XC1,xc2
      data xc1/0.3551D0/,XC2/0.8579D0/
      real*8 c11,c12,c13,c14,c15,c21,c22,c23,c24,c25,c31,c32,c33,c34
      data C11/1.5701681D0/, C12/-1.4937488D0/,C13/1.7410516D0/,
     % C14/-1.6504717D0/,C15/0.81410785D0/,
     $     C21/ 0.55852340D0/, C22/1.0262830D0/,C23/-1.4970125D0/,
     % C24/1.1671599D0/, C25/-0.37018753D0/,
     $     C31/-0.16539043D0/,C32/0.083425502D0/,C33/-0.042281332D0/,
     % C34/0.010874519D0/
      real*8 X1,X2

      IF(X.GE.XC2) THEN
        X1=1/X
        BK1323=(((C34*X1+C33)*X1+C32)*X1+C31)*X1+1D0
      ELSEIF(X.GE.XC1) THEN
        BK1323=(((C25*X+C24)*X+C23)*X+C22)*X+C21
      ELSEIF(X.GT.0D0) THEN
        X1=X**(0.333333333333333333D0)
        X2=X1**2
        BK1323=((((C15*X2+C14)*X2+C13)*X2+C12)*X2+C11)*X1
      ELSE
        BK1323=0
      ENDIF
      RETURN
      END

      subroutine radangle(GAM,RHO,U,THX,THY,XI3,XI2)
      use tmacro
      implicit none
      real*8 , intent(in) :: GAM,RHO,U
      real*8 , intent(out)::THX,THY,XI3,XI2
      real*8 P1,tran,V,X,COMPTON,UPSILON,
     $     PHI,TH,THYG,THYG2,Z,R,BK1323
      COMPTON=rclassic/rcratio
      UPSILON=COMPTON*GAM**2/RHO
      P1=tran()
      V=2*SQRT(2D0)*SIN(ASIN(5*P1/(4*SQRT(2D0)))/3)
      X=1.5D0*UPSILON*U*V**3
      PHI=2*PI*tran()
      TH=SQRT(MAX(0D0,1-V**2))/(V*GAM)
      THY=TH*SIN(PHI)
C          TH*COS(PHI) can be x-angle
      thx=th*cos(phi)
      THYG=GAM*THY
      THYG2=SQRT(1+THYG**2)
      Z=X/(3*UPSILON)*THYG2**3
      R=-THYG/THYG2*BK1323(Z)
      XI2=2*R/(1+R**2)
      XI3=(1-R**2)/(1+R**2)
      RETURN
      END
