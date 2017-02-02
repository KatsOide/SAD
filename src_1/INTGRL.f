      subroutine INTGRL(LATT,TWISS,IDP,DP,LFNO)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      PARAMETER (ML=30)
*     IMPLICIT REAL*8 'inc/A-H,O-Z.inc'   *** includeD IN TFMACRO ***
      DIMENSION LATT(2,NLAT),TWISS(NLAT,-NDIM:NDIM,ntwissfun)
*    TWP(1)=AX    (2)=BX    (3)=PHX
*    TWP(4)=AY    (5)=BY    (6)=PHY
*    TWP(7)=EX    (8)=EPX   (9)=EY   (10)=EPY
*    TWP(11)=R1   (12)=R2   (13)=R3  (14)=R4
*    TWP(15)=DX   (16)=DXP  (17)=DY  (18)=DYP
      COMMON /TWISSP/TWPIN(ntwissfun),TWPOT(ntwissfun),TWPMG(ntwissfun)
      COMMON /ELEM/ANG,E1,E2,TETA,RLE,RKL1,RKL2
      COMMON /DELEM/DL,DKL0,DKL1,DKL2,CKL0X,CKL0Y,CKL1X,CKL1Y
      DIMENSION RI(15),DI(15),DIEG(15)
      COMMON /FUNS/FUNI(ML,15)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      DATA COEFF/7.0394513D-6/
      DATA CO4/1.2113938D-3/,CO5/1.4674773D-6/
*   NRI  NUMBER OF RADIATION INTEGRALS  MUST < 15 (DIM OF FUNI,DI,RI)
      DATA NRI/12/,IPOUT/0/,IB,IQ/0,0/
c
      save
c
      DO 1 I=1,NRI
    1 RI(I)=0.
      RFV=0.
      RLT=0.
      TOTANG=0.
      ENERGY=0.511D-3*H0
      WRITE(LFNO,'(A/)')
     +'****************************************************************'
      WRITE(LFNO,'(A)')
     + '     RADIATION INTEGRALS OF A LATTICE WITH X-Y COUPLING'
      WRITE(LFNO,'(/A)')
     +'     CONTRIBUTIONS FROM SOLENOID MAGNETS ARE NOT includeD '
      WRITE(LFNO,'(A)') '              IN THIS VERSION.   1991.2.26 '
      WRITE(LFNO,'(A)')
     +'****************************************************************'
      WRITE(LFNO,'(/ A,I6)') 'NLAT =',NLAT
      WRITE(LFNO,'(// A,A)') '           I1U      I2       I3   ',
     +         '   I4U      I5U      I1V      I4V      I5V'
      DO 10 J=1,NLAT
         DO 11 I=1,ntwissfun
         TWPIN(I)=TWISS(J,0,I)
         TWPMG(I)=TWPIN(I)
         IF(J.LT.NLAT) THEN
           TWPOT(I)=TWISS(J+1,0,I)
         ELSE
           TWPOT(I)=TWPIN(I)
         END IF
   11    CONTINUE
         DO 12 I=1,15
            DI(I)=0.
   12       DIEG(I)=0.
         CKL0X=0.
         CKL0Y=0.
         CKL1X=0.
         CKL1Y=0.
         ID=IDTYPE(LATT(1,J))
         RLE=RLIST(IDVAL(LATT(1,J))+1)
         RLT=RLT+RLE
*           WRITE(LFNO,'(A,I5,2F13.5)') ' ELEMENT  ',ID,RLE,RLT
         IF(ID.EQ.2) THEN
            IB=IB+1
            ANG=RLIST(LATT(2,J)+2)
            E1=RLIST(LATT(2,J)+3)*ANG
            E2=RLIST(LATT(2,J)+4)*ANG
            TETA=RLIST(LATT(2,J)+5)
            IF(ANG.EQ.0) GOTO 10
            IR=DABS(ANG)/0.002
            IF(IR.EQ.0) IR=1
            IF(IR.GT.14) IR=14
            IRP=2*IR
            DL=RLE/IRP
            DKL0=ANG/IRP
            TOTANG=TOTANG+ANG
*      WRITE(LFNO,'(A,6F12.5)') ' BM ',ANG,E1,E2,RLE,TOTANG,RLT
            IF(IPOUT.EQ.1.AND.IB.LE.10) THEN
            WRITE(LFNO,'(// A,3F12.5)') PNAME(LATT(1,J)),RLE,RKL1,RLT
            WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))') ' TWPIN ',J,TWPIN
            WRITE(LFNO,'(A,I4/ 10F8.3/ 1P,(10D8.1))') ' TWPOT ',J,TWPOT
            WRITE(LFNO,'(/ A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
            END IF
            IF(E1.NE.0.) THEN
               CALL BEDGE(E1)
               CALL TEDGE(DIEG)
               IF(IPOUT.EQ.1.AND.IB.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
            ENDIF
            DO 100 I=1,IRP
               CALL BMAG(TETA)
               CALL TCAL(I,0)
               IF(IPOUT.EQ.1.AND.IB.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
100         CONTINUE
               CALL TCAL(IRP+1,1)
            IF(E2.NE.0.) THEN
               CALL BEDGE(E2)
               CALL TEDGE(DIEG)
               IF(IPOUT.EQ.1.AND.IB.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
            ENDIF
         ELSE IF(ID.EQ.4) THEN
            IQ=IQ+1
            RKL1=RLIST(LATT(2,J)+2)
            TETA=RLIST(LATT(2,J)+4)
            IF(IPOUT.EQ.1.AND.IQ.LE.10) THEN
            WRITE(LFNO,'(// A,3F12.5)') PNAME(LATT(1,J)),RLE,RKL1,RLT
            WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))') ' TWPIN ',J,TWPIN
            WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))') ' TWPOT ',J,TWPOT
            END IF
            IF(RKL1.EQ.0) GOTO 10
            IR=DABS(RKL1)/0.02
            IF(IR.EQ.0) IR=1
            IF(IR.GT.14) IR=14
            IRP=2*IR
            DL=RLE/IRP
            DKL1=RKL1/IRP
            DO 200 I=1,IRP
               IF(IPOUT.EQ.1.AND.IQ.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
               CALL QUADRU(TETA,TWPMG(15),TWPMG(17))
               CALL TCAL(I,0)
200         CONTINUE
               CALL TCAL(IRP+1,1)
               IF(IPOUT.EQ.1.AND.IQ.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
         ELSE IF(ID.EQ.6) THEN
            RKL2=RLIST(LATT(2,J)+2)
            TETA=RLIST(LATT(2,J)+4)
*           WRITE(LFNO,'(A,3F12.5)') ' S  ',RLE,RKL2,RLT
            CALL SEXTU(TETA,TWPMG(15),TWPMG(17))
            IR=0
            CALL TCAL(1,0)
            IF(RLE.EQ.0.) THEN
               DL=0.2
               WRITE(LFNO,'(A)') ' L(SEXTUPOLE) EQ 0. SET L=0.2'
            ENDIF
         ELSE IF(ID.EQ.31) THEN
            RFV=RFV+RLIST(LATT(2,J)+2)
            RHARM=RLIST(LATT(2,J)+3)
            FREQ=RLIST(LATT(2,J)+5)
*           WRITE(LFNO,*) RLIST(LATT(2,J)+3),RLIST(LATT(2,J)+5)
         END IF
         IF(ID.EQ.2.OR.ID.EQ.4.OR.ID.EQ.6) THEN
*      WRITE(LFNO,'( 7F10.3)') AX0,AX(J+1),BX0,BX(J+1),GX0,EX0,EX(J+1)
            DO 50 I=1,NRI
               DI(I)=SYMPS(FUNI(1,I),IR)+DIEG(I)
50             RI(I)=RI(I)+DI(I)
            WRITE(LFNO,'( A,8(1PD9.2))') PNAME(LATT(1,J)),
     +           (DI(I),I=1,8)
         END IF
   10 CONTINUE
      QX=TWPOT(3)/(2*PI)
      QY=TWPOT(6)/(2*PI)
      WRITE(LFNO,'(//A)') '    SYNCHROTRON RADIATION INTEGRALS'
      IF(RHARM.EQ.0) RHARM=FREQ*RLT/CVELOC
      IF(FREQ.EQ.0) FREQ=RHARM*CVELOC/RLT
*     WRITE(LFNO,'(A,E12.4,A,E12.4)') ' P0 = ',P0,'   H0 = ',H0
*     BFIELD=ENERGY*1.0E9/(CVELOC*RHO)
*     RAMBDAC=18.64/(BFIELD*ENERGY*ENERGY)
      WRITE(LFNO,'(/ A,F9.3)') ' ENERGY (GEV)               = ',ENERGY
      WRITE(LFNO,'(A,F9.3)') ' DP/P                       = ',DP0
      WRITE(LFNO,'(A,F9.3)') ' TOTAL LENGTH (M)           = ',RLT
      TOTANG=TOTANG*180./PI
      WRITE(LFNO,'(A,F9.3)') ' TOTAL BENDING ANGLE        = ',TOTANG
*     WRITE(LFNO,'(A,F9.4)') ' BFIELD (TESLA)             = ',BFIELD
*     WRITE(LFNO,'(A,F6.3)') ' LAMBDAC (A)               = ',RAMBDAC
      RFV=RFV/1000.
      WRITE(LFNO,'(A,F15.3)') ' RF VOLTS   (KV)            = ',RFV
      FREQ=FREQ*1.E-6
      WRITE(LFNO,'(A,F6.0)') ' HARMONIC NUMBER            = ',RHARM
      WRITE(LFNO,'(A,F10.4)') ' RF FREQUENCY (MHZ)        = ',FREQ
      WRITE(LFNO,'(/ A)') ' RADIATION INTEGRALS  '
      WRITE(LFNO,'(2(A,1PD12.3))') '   I1X = ',RI(1),'    I1Y = ',RI(6)
      WRITE(LFNO,'(2(A,1PD12.3))') '   I2  = ',RI(2)
      WRITE(LFNO,'(2(A,1PD12.3))') '   I3  = ',RI(3)
      WRITE(LFNO,'(2(A,1PD12.3))') '   I4U = ',RI(4),'    I4V = ',RI(7)
      WRITE(LFNO,'(2(A,1PD12.3))') '   I5U = ',RI(5),'    I5V = ',RI(8)
      WRITE(LFNO,'(A,1PD12.3)') '        TOTAL BENDING ANGL = ',RI(11)
      ELOSS=2.0D6*COEFF*RI(2)*ENERGY**4
      ALPHAX=RI(1)/RLT
      ALPHAY=RI(6)/RLT
      DAMPU=1.-RI(4)/RI(2)
      DAMPV=1.-RI(7)/RI(2)
      CO2=CVELOC*ENERGY**3*COEFF/RLT
      TAUU=1000.0/(CO2*(RI(2)-RI(4)))
      TAUV=1000.0/(CO2*(RI(2)-RI(7)))
      TAUE=1000.0/(CO2*(2*RI(2)+RI(4)+RI(7)))
      SIG=CO4*ENERGY*DSQRT(RI(3)/(2.*RI(2)+RI(4)+RI(7)))
      EMITU=CO5*ENERGY*ENERGY*RI(5)/(RI(2)-RI(4))
      EMITV=CO5*ENERGY*ENERGY*RI(8)/(RI(2)-RI(7))
      CHRX=-RI(9)/4./PI
      CHRY=-RI(10)/4./PI
      WRITE(LFNO,'(/ A,F8.5,A,F8.5)')
     +  ' NU-X-  = ',QX,'    NU-Y- = ',QY
      WRITE(LFNO,'(A,F12.6)') ' MOMENTUM COMPACTION FAC X  = ',ALPHAX
      WRITE(LFNO,'(A,F12.6)') ' MOMENTUM COMPACTION FAC Y  = ',ALPHAY
      WRITE(LFNO,'(A,F12.6)') ' ENERGY LOSS PER TURN (KV)  = ',ELOSS
      WRITE(LFNO,'(A,F12.8)') ' DU=1-I4U/I2                = ',DAMPU
      WRITE(LFNO,'(A,F12.8)') ' DV=1-I4V/I2                = ',DAMPV
      WRITE(LFNO,'(/ A)') ' DAMPING TIME CONSTANTS'
      WRITE(LFNO,'(A,F12.6)') '                TAUU (MSEC) = ',TAUU
      WRITE(LFNO,'(A,F12.6)') '                TAUV (MSEC) = ',TAUV
      WRITE(LFNO,'(A,F12.6)') '                TAUE (MSEC) = ',TAUE
      WRITE(LFNO,'(A,1PE12.4)') ' ENERGY SPREAD (SIG) DE/E   = ',SIG
      WRITE(LFNO,'(A,1PE12.4)') '            EMITTANCE U     = ',EMITU
      WRITE(LFNO,'(A,1PE12.4)') '            EMITTANCE V     = ',EMITV
      IF (ELOSS.LE.RFV) THEN
         SIPH=ELOSS/RFV
         COPH=SQRT(1.-SIPH*SIPH)
         PHIS=PI-DATAN(SIPH/COPH)
         ALPHA=ALPHAX+ALPHAY
         QS=0.00039894228D0*DSQRT(RHARM*ALPHA*RFV*COPH/ENERGY)
         SIGL=1.D6*SIG*RLT*DSQRT((ALPHA*ENERGY)/(2*PI*RHARM*RFV*COPH))
         OVERV=RFV/ELOSS
         OVERV=SQRT(OVERV*OVERV-1)
         RKQ=2.0*(OVERV-ATAN(OVERV))
         EPMAX=SQRT(1.D-6*ELOSS*RKQ/(PI*ALPHA*RHARM*ENERGY))
         XI=0.5*(EPMAX/SIG)**2
         WRITE(LFNO,'(/ A,F12.6)') ' NU-S                       = ',QS
         WRITE(LFNO,'(A,F12.6)') ' BUCKET HIGHT   DEMAX/E     = ',EPMAX
         WRITE(LFNO,'(A,F12.6)') ' NATURAL BUNCH LENGTH (MM)  = ',SIGL
         IF(XI.LT.100) THEN
         TAUQ=TAUE/1000.*DEXP(XI)/2./XI
         TAUQ=EXP(2.3026*TAUQ)
         WRITE(LFNO,'(A,1PD12.4)') ' QUANTUM LIFETIME(HOURS)    = ',TAUQ
         ELSE
         WRITE(LFNO,'(A,1PD12.4)') ' Q.LIFE IS TOO LONG      XI = ',XI
         END IF
      ELSE
         WRITE(LFNO,*) ' NOT ENOUGH RF VOLTS;  ELOSS,RFV = ',ELOSS,RFV
      END IF
      WRITE(LFNO,'(/ A,F12.6)') ' LINEAR CHROMATICITY KX*BU  = ',CHRX
      WRITE(LFNO,'(A,F12.6)') ' LINEAR CHROMATICITY KY*BV  = ',CHRY
      END
