      module radint
      use tffitcode,only:ntwissfun
      integer*4 , parameter ::NRI=12,IPOUT=0,ML=30
      real*8 TWPIN(ntwissfun),TWPOT(ntwissfun),TWPMG(ntwissfun)
      real*8 ANG,E1,E2,TETA,RLE,RKL1,RKL2
      real*8  DL,DKL0,DKL1,DKL2,CKL0X,CKL0Y,CKL1X,CKL1Y
      real*8 FUNI(ML,NRI)
      real*8 TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)

      contains
      subroutine INTGRL(LFNO)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer
      use kyparam
      use macphys
      implicit none
      integer*4 ,intent(in)::LFNO
      integer*4 IB,IQ
*    TWP(1)=AX    (2)=BX    (3)=PHX
*    TWP(4)=AY    (5)=BY    (6)=PHY
*    TWP(7)=EX    (8)=EPX   (9)=EY   (10)=EPY
*    TWP(11)=R1   (12)=R2   (13)=R3  (14)=R4
*    TWP(16)=DX   (17)=DXP  (18)=DY  (19)=DYP
      type (sad_comp) ,pointer ::cmp
      real*8  RI(NRI),DI(NRI),DIEG(NRI)
      real*8 TOTANG,ENERGY,ALPHA,ALPHAX,ALPHAY,CHRX,CHRY,CO2,COPH,
     $     DAMPU,DAMPV,ELOSS,EMITU,EMITV,EPMAX,FREQ,PHIS,QX,QY,QS,
     $     RFV,RHARM,RKQ,RLT,SIG,SIGL,SIPH,TAUE,TAUQ,TAUU,TAUV,
     $     XI,OVERV,SYMPS
      integer*4 J,idc,ID,I,IR,IRP
      real*8 ,parameter ::
     $     COEFF=(elradi/(elmass**3)/3)*1e27,
     $     CO4=sqrt(55.d0/32.d0/sqrt(3.d0)*elradi/
     $     finest)/elmass*1e9,CO5=CO4**2
c      DATA COEFF/7.0394513D-6/
c      DATA CO4/1.2113938D-3/,CO5/1.4674773D-6/
*   NRI  NUMBER OF RADIATION INTEGRALS  MUST < 15 (DIM OF FUNI,DI,RI)
c
      IB=0
      IQ=0
      RI=0.d0
      RFV=0.d0
      RLT=0.d0
      TOTANG=0.d0
      ENERGY=amass/1e9*H0
      WRITE(LFNO,'(A/)')
     +'****************************************************************'
      WRITE(LFNO,'(A)')
     + '     RADIATION INTEGRALS OF A LATTICE WITH X-Y COUPLING'
      WRITE(LFNO,'(/A)')
     +'     CONTRIBUTIONS FROM SOLENOID MAGNETS ARE NOT included '
      WRITE(LFNO,'(A)') '              IN THIS VERSION.   1991.2.26 '
      WRITE(LFNO,'(A)') '                  Refurbished:  2019.10.10 '
      WRITE(LFNO,'(A)')
     +'****************************************************************'
c      WRITE(LFNO,'(/ A,I6)') 'NLAT =',NLAT
c      WRITE(LFNO,'(// A,A)') '           I1U      I2       I3   ',
c     +         '   I4U      I5U      I1V      I4V      I5V'
      DO 10 J=1,NLAT-1
        call compelc(j,cmp)
        idc=cmp%id
        ID=idtype(idc)
        DO 11 I=1,ntwissfun
         TWPIN(I)=TWISS(J,0,I)
         TWPMG(I)=TWPIN(I)
         IF(J.LT.NLAT) THEN
           TWPOT(I)=TWISS(J+1,0,I)
         ELSE
           TWPOT(I)=TWPIN(I)
         END IF
   11    CONTINUE
         DI=0.d0
         DIEG=0.d0
         CKL0X=0.
         CKL0Y=0.
         CKL1X=0.
         CKL1Y=0.
c         ID=IDTYPE(LATT(1,J))
         RLE=cmp%value(kytbl(kwL,ID))
c         RLE=RLIST(IDVAL(LATT(1,J))+1)
         RLT=RLT+RLE
c         WRITE(LFNO,'(A,2I5,1p2g13.5)') ' ELEMENT  ',J,ID,RLE,RLT
         select case (ID)
         case(icBEND)
            IB=IB+1
            ANG=cmp%value(ky_ANGL_BEND)
            E1=cmp%value(ky_E1_BEND)*ANG
            E2=cmp%value(ky_E2_BEND)*ANG
            TETA=cmp%value(ky_ROT_BEND)
c            ANG=RLIST(LATT(2,J)+2)
c            E1=RLIST(LATT(2,J)+3)*ANG
c            E2=RLIST(LATT(2,J)+4)*ANG
c            TETA=RLIST(LATT(2,J)+5)
            IF(ANG.EQ.0) GOTO 10
            IR=min(14,max(1,int(ABS(ANG)/0.002d0)))
            IRP=2*IR
            DL=RLE/IRP
            DKL0=ANG/IRP
            TOTANG=TOTANG+ANG
*      WRITE(LFNO,'(A,6F12.5)') ' BM ',ANG,E1,E2,RLE,TOTANG,RLT
            IF(IPOUT.EQ.1.AND.IB.LE.10) THEN
              WRITE(LFNO,'(// A,3F12.5)') PNAME(idc),RLE,RKL1,RLT
              WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))')
     $             ' TWPIN ',J,TWPIN
              WRITE(LFNO,'(A,I4/ 10F8.3/ 1P,(10D8.1))')
     $             ' TWPOT ',J,TWPOT
              WRITE(LFNO,'(/ A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
            END IF
            IF(E1.NE.0.) THEN
               CALL BEDGE(E1)
               CALL TEDGE(DIEG)
c               write(*,*)'bend-e1 ',i,E1,tu(1,1),tu(2,2)
               IF(IPOUT.EQ.1.AND.IB.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
            ENDIF
            DO 100 I=1,IRP
               CALL BMAG(TETA)
c               write(*,*)'bend-i ',i,tu(1,1),tu(2,2)
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
          case (icQUAD)
            IQ=IQ+1
            RKL1=cmp%value(ky_K1_QUAD)
            TETA=cmp%value(ky_ROT_QUAD)
c            RKL1=RLIST(LATT(2,J)+2)
c            TETA=RLIST(LATT(2,J)+4)
            IF(IPOUT.EQ.1.AND.IQ.LE.10) THEN
              WRITE(LFNO,'(// A,3F12.5)') PNAME(idc),RLE,RKL1,RLT
              WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))')
     $             ' TWPIN ',J,TWPIN
              WRITE(LFNO,'(/A,I4/ 10F8.3/ 1P,(10D8.1))')
     $             ' TWPOT ',J,TWPOT
            END IF
            IF(RKL1.EQ.0) GOTO 10
            IR=int(ABS(RKL1)/0.02d0)
            IF(IR.EQ.0) IR=1
            IF(IR.GT.14) IR=14
            IRP=2*IR
            DL=RLE/IRP
            DKL1=RKL1/IRP
            DO 200 I=1,IRP
               IF(IPOUT.EQ.1.AND.IQ.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
               CALL QUADRU(TETA,TWPMG(mfitdx),TWPMG(mfitdy))
               CALL TCAL(I,0)
200         CONTINUE
               CALL TCAL(IRP+1,1)
               IF(IPOUT.EQ.1.AND.IQ.LE.10)
     +         WRITE(LFNO,'(A/ 10F8.3/ 1P,(10D8.1))') ' TWPMG ',TWPMG
         case (icSEXT)
            RKL2=cmp%value(ky_K_THIN)
            TETA=cmp%value(ky_ROT_THIN)
c            RKL2=RLIST(LATT(2,J)+2)
c            TETA=RLIST(LATT(2,J)+4)
*           WRITE(LFNO,'(A,3F12.5)') ' S  ',RLE,RKL2,RLT
            CALL SEXTU(TETA,TWPMG(mfitdx),TWPMG(mfitdy))
            IR=0
            CALL TCAL(1,0)
            IF(RLE.EQ.0.) THEN
               DL=0.2
               WRITE(LFNO,'(A)') ' L(SEXTUPOLE) EQ 0. SET L=0.2'
            ENDIF
         CASE (icCAVI)
            RFV=cmp%value(ky_VOLT_CAVI)
            RHARM=cmp%value(ky_HARM_CAVI)
            FREQ=cmp%value(ky_FREQ_CAVI)
c            RFV=RFV+RLIST(LATT(2,J)+2)
c            RHARM=RLIST(LATT(2,J)+3)
c            FREQ=RLIST(LATT(2,J)+5)
*           WRITE(LFNO,*) RLIST(LATT(2,J)+3),RLIST(LATT(2,J)+5)
         END SELECT
         IF(ID.EQ.icBEND.OR.ID.EQ.icQUAD.OR.ID.EQ.icSEXT) THEN
*      WRITE(LFNO,'( 7F10.3)') AX0,AX(J+1),BX0,BX(J+1),GX0,EX0,EX(J+1)
           DO I=1,NRI
             DI(I)=SYMPS(FUNI(1,I),IR)+DIEG(I)
             RI(I)=RI(I)+DI(I)
           enddo
c            WRITE(LFNO,'( A,8(1PD9.2))') PNAME(idc),
c     +           (DI(I),I=1,8)
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
      WRITE(LFNO,'(/ A,F9.3)') ' ENERGY (GeV)               = ',ENERGY
      WRITE(LFNO,'(A,F9.3)')   ' DP/P                       = ',DP0
      WRITE(LFNO,'(A,F9.3)')   ' TOTAL LENGTH (m)           = ',RLT
      TOTANG=TOTANG*180./PI
      WRITE(LFNO,'(A,F9.3)')   ' TOTAL BENDING ANGLE (DEG)  = ',TOTANG
      WRITE(LFNO,'(A,1pg15.7)')' TOTAL BENDING ANGLE (RAD)  = ',RI(11)
*     WRITE(LFNO,'(A,F9.4)') ' BFIELD (TESLA)             = ',BFIELD
*     WRITE(LFNO,'(A,F6.3)') ' LAMBDAC (A)               = ',RAMBDAC
      RFV=RFV/1000.
      WRITE(LFNO,'(A,1pg15.7)')' RF VOLTAGE (kV)            = ',RFV
      FREQ=FREQ*1.E-6
      WRITE(LFNO,'(A,F6.0)')   ' HARMONIC NUMBER            = ',RHARM
      WRITE(LFNO,'(A,F10.4)')  ' RF FREQUENCY (MHz)         = ',FREQ
      WRITE(LFNO,'(/ A)') ' RADIATION INTEGRALS  '
      WRITE(LFNO,'(2(A,1pg15.7))') '   I1X = ',RI(1),'    I1Y = ',RI(6)
      WRITE(LFNO,'(2(A,1pg15.7))') '   I2  = ',RI(2)
      WRITE(LFNO,'(2(A,1pg15.7))') '   I3  = ',RI(3)
      WRITE(LFNO,'(2(A,1pg15.7))') '   I4U = ',RI(4),'    I4V = ',RI(7)
      WRITE(LFNO,'(2(A,1pg15.7))') '   I5U = ',RI(5),'    I5V = ',RI(8)
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
      CHRX=-RI(9)/4.d0/PI
      CHRY=-RI(10)/4.d0/PI
      WRITE(LFNO,'(/ A,F8.5,A,F8.5)')
     +  ' NU-X-  = ',QX,'    NU-Y- = ',QY
      WRITE(LFNO,'(A,1pg12.4)') ' MOMENTUM COMPACTION FAC X  = ',ALPHAX
      WRITE(LFNO,'(A,1pg12.4)') ' MOMENTUM COMPACTION FAC Y  = ',ALPHAY
      WRITE(LFNO,'(A,1pG12.4)') ' ENERGY LOSS PER TURN (kV)  = ',ELOSS
      WRITE(LFNO,'(A,F12.8)') ' DU=1-I4U/I2                = ',DAMPU
      WRITE(LFNO,'(A,F12.8)') ' DV=1-I4V/I2                = ',DAMPV
      WRITE(LFNO,'(/ A)') ' DAMPING TIME CONSTANTS'
      WRITE(LFNO,'(A,1pg12.4)') '                TAUU (msec) = ',TAUU
      WRITE(LFNO,'(A,1pg12.4)') '                TAUV (msec) = ',TAUV
      WRITE(LFNO,'(A,1pg12.4)') '                TAUE (msec) = ',TAUE
      WRITE(LFNO,'(A,1Pg12.4)') ' ENERGY SPREAD (SIG) DE/E   = ',SIG
      WRITE(LFNO,'(A,1Pg12.4)') '            EMITTANCE U     = ',EMITU
      WRITE(LFNO,'(A,1Pg12.4)') '            EMITTANCE V     = ',EMITV
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
         WRITE(LFNO,'(/ A,F12.6)')' NU-S                       = ',QS
         WRITE(LFNO,'(A,F12.6)')  ' BUCKET HIGHT   DEMAX/E     = ',EPMAX
         WRITE(LFNO,'(A,F12.6)')  ' NATURAL BUNCH LENGTH (mm)  = ',SIGL
         IF(XI.LT.100) THEN
         TAUQ=TAUE/1000.*DEXP(XI)/2./XI
         TAUQ=EXP(2.3026*TAUQ)
         WRITE(LFNO,'(A,1PD12.4)') ' QUANTUM LIFETIME(hours)    = ',TAUQ
         ELSE
         WRITE(LFNO,'(A,1PD12.4)') ' Q.LIFE IS TOO LONG      XI = ',XI
         END IF
      ELSE
         WRITE(LFNO,*) ' NOT ENOUGH RF VOLTS;  ELOSS,RFV = ',ELOSS,RFV
      END IF
      WRITE(LFNO,'(/ A,F12.6)') ' LINEAR CHROMATICITY GX     = ',CHRX
      WRITE(LFNO,'(A,F12.6)')   ' LINEAR CHROMATICITY GY     = ',CHRY
      END subroutine

      subroutine TCAL(I,IFLG)
      use tffitcode
      implicit none
      real*8, save:: R0(4,4),R0I(4,4)
      integer*4 I,IFLG
      real*8 CURV,CURV2,CURV3,GU,GV,RK1
      associate (
     $     AU=>TWPMG(mfitax),
     $     BU=>TWPMG(mfitbx),
     $     PU=>TWPMG(mfitnx),
     $     AV=>TWPMG(mfitay),
     $     BV=>TWPMG(mfitby),
     $     PV=>TWPMG(mfitny),
     $     EU=>TWPMG(mfitex),
     $     EPU=>TWPMG(mfitepx),
     $     EV=>TWPMG(mfitey),
     $     EPV=>TWPMG(mfitepy),
     $     R1=>TWPMG(mfitr1),
     $     R2=>TWPMG(mfitr2),
     $     R3=>TWPMG(mfitr3),
     $     R4=>TWPMG(mfitr4),
     $     DETR=>TWPMG(mfitdetr),
     $     DU=>TWPMG(mfitdx),
     $     DPU=>TWPMG(mfitdpx),
     $     DV=>TWPMG(mfitdy),
     $     DPV=>TWPMG(mfitdpy))
      IF(DL.EQ.0.) WRITE(6,'(A)') '  DL=0. IN TCAL '
      RK1=CKL1X/DL
      CURV=DSQRT(CKL0X*CKL0X+CKL0Y*CKL0Y)/DL
      GU=(1+AU*AU)/BU
      GV=(1+AV*AV)/BV
      IF(I.EQ.1) CALL R0CAL(R0,R0I)
      CURV2=CURV*CURV
      CURV3=CURV2*ABS(CURV)
      FUNI(I,1)=CKL0X*
     +          (R0I(1,1)*EU+R0I(1,2)*EPU+R0I(1,3)*EV+R0I(1,4)*EPV)
      FUNI(I,2)=CURV2*DL
      FUNI(I,3)=CURV3*DL
      IF(CURV.GT.0) THEN
      FUNI(I,4)=(CKL0X+2*CKL1X/CURV)*CURV2*(R0I(1,1)*EU+R0I(1,2)*EPU)
     +         +(CKL0Y-2*CKL1Y/CURV)*CURV2*(R0I(3,1)*EU+R0I(3,2)*EPU)
      ELSE
      FUNI(I,4)=0.
      END IF
      FUNI(I,5)=CURV3*(GU*EU*EU+2.*AU*EPU*EU+BU*EPU*EPU)*DL
c      write(*,*)'tcal-I5 ',funi(i,5),CURV3,DL
* 6: I1Y    7:I4Y     8:I5Y
      FUNI(I,6)=CKL0Y*
     +          (R0I(3,1)*EU+R0I(3,2)*EPU+R0I(3,3)*EV+R0I(3,4)*EPV)
      IF(CURV.GT.0) THEN
      FUNI(I,7)=(CKL0X+2*CKL1X/CURV)*CURV2*(R0I(1,3)*EV+R0I(1,4)*EPV)
     +         +(CKL0Y-2*CKL1Y/CURV)*CURV2*(R0I(3,3)*EV+R0I(3,4)*EPV)
      ELSE
      FUNI(I,7)=0.
      END IF
      FUNI(I,8)=CURV3*(GV*EV*EV+2.*AV*EPV*EV+BV*EPV*EPV)*DL
* 9:  CHROMX   10: CHROMY
c      FUNI(I,9)=(RK1+CURV2*0.d0)*BU*DL
c      FUNI(I,10)=-RK1*BV*DL
      FUNI(I,9) = GU*DL
      FUNI(I,10)= GV*DL
c      write(*,'(a,1p6g15.7)')'chrom ',RK1,CURV2,BU,BV,DL
      FUNI(I,11)=CKL0X
      FUNI(I,12)=CKL0Y
      IF(IFLG.EQ.1) RETURN
      CALL TMATRS(R0,R0I)
      CALL TWTRANS()
      end associate
      END subroutine

      subroutine BMAG(TETA1)
      implicit none
      real*8 TETA1,SI,CO,CURV
      CURV=DKL0/DL
      CALL ZEROV4(VM)
      CALL UNITM4(TM,1.d0)
      SI=DSIN(DKL0)
      CO=DCOS(DKL0)
      TM(1,1)=CO
      IF(CURV.NE.0.) THEN
        TM(1,2)=SI/CURV
        TM(2,1)=-SI*CURV
        VM(1)=(1-CO)/CURV
        VM(2)= SI
      ELSE
        TM(1,2)=DL
      END IF
      TM(2,2)=TM(1,1)
      TM(3,4)=DL
      CKL0X=DCOS(TETA1)*DKL0
      CKL0Y=-DSIN(TETA1)*DKL0
      CKL1X=0.
      CKL1Y=0.
      CALL TROT(TM,VM,TETA1,TO,VO)
      RETURN
      END subroutine

      subroutine TMATRS(R0,R0I)
      use tffitcode
      implicit none
      real*8 A(2,2),B(2,2),C(2,2),D(2,2),AI(2,2),CAI(2,2),ACI(2,2)
      real*8 RA1(2,2),RB1(2,2),RC1(2,2),RD1(2,2)
      real*8 TR(4,4),R0(4,4),R0I(4,4)
      real*8 DETCAI,SMU1,RMU1
      integer*4 J,K
      TR=matmul(TO,R0I)
c      WRITE(6,'(A)') 'TMATRS: TO, R0I, TR, R0, TU'
c      WRITE(6,'(8F10.4)') TO
c      WRITE(6,'(8F10.4)') R0I
c      WRITE(6,'(8F10.4)') TR
      CALL DM42(TR,A,B,C,D)
      CALL MINV2(A,AI)
      CALL MMUL2(C,AI,CAI)
      DETCAI=CAI(1,1)*CAI(2,2)-CAI(1,2)*CAI(2,1)
      IF(DETCAI.LT. 1.d0) THEN
         SMU1=1./(1.+DETCAI)
         RMU1=DSQRT(SMU1)
         CALL FMUL2(-RMU1,CAI,RC1)
         CALL SYMTR2(RC1,RB1)
         CALL UNITM2(RA1,RMU1)
         CALL UNITM2(RD1,RMU1)
         CALL CM24(RA1,RB1,RC1,RD1,R0)
         CALL R1INV(R0,R0I)
         DO J=1,2
           DO K=1,2
             TWPMG(10+K+2*(J-1))=RC1(J,K)
           enddo
         enddo
      ELSE
         SMU1=1./(1.+1./DETCAI)
         RMU1=DSQRT(SMU1)
         CALL MINV2(CAI,ACI)
         CALL FMUL2(-RMU1,ACI,RD1)
         CALL SYMTR2(RD1,RA1)
         CALL UNITM2(RB1,RMU1)
         CALL UNITM2(RC1,RMU1)
         CALL CM24(RA1,RB1,RC1,RD1,R0)
         CALL R2INV(R0,R0I)
         DO J=1,2
           DO K=1,2
             TWPMG(10+K+2*(2-J))=-RD1(J,K)
           enddo
         enddo
      END IF
      TU=matmul(R0,TR)
c      WRITE(6,'(8F10.4)') R0
c      WRITE(6,'(8F10.4)') TU
      VU=matmul(R0,VO)
c      CALL VMUL4(R0,VO,VU)
      RETURN
      END subroutine

      subroutine TWTRANS()
      use tffitcode
      use tmacro
      implicit none
      real*8 TTX(3,3),TTY(3,3),DX1(4)
      real*8 GX,GY,BX1,AX1,BY1,AY1,EX1,EPX1,EY1,EPY1
      associate (
     $     AX=>TWPMG(mfitax),
     $     BX=>TWPMG(mfitbx),
     $     PX=>TWPMG(mfitnx),
     $     AY=>TWPMG(mfitay),
     $     BY=>TWPMG(mfitby),
     $     PY=>TWPMG(mfitny),
     $     EX=>TWPMG(mfitex),
     $     EPX=>TWPMG(mfitepx),
     $     EY=>TWPMG(mfitey),
     $     EPY=>TWPMG(mfitepy),
     $     R1=>TWPMG(mfitr1),
     $     R2=>TWPMG(mfitr2),
     $     R3=>TWPMG(mfitr3),
     $     R4=>TWPMG(mfitr4),
     $     DETR=>TWPMG(mfitdetr),
     $     DX=>TWPMG(mfitdx:mfitdpy))
      TTX(1,1)=TU(1,1)*TU(1,1)
      TTX(1,2)=-2.*TU(1,1)*TU(1,2)
      TTX(1,3)=TU(1,2)*TU(1,2)
      TTX(2,1)=-TU(1,1)*TU(2,1)
      TTX(2,2)=1.+2.*TU(1,2)*TU(2,1)
      TTX(2,3)=-TU(1,2)*TU(2,2)
      TTX(3,1)=TU(2,1)*TU(2,1)
      TTX(3,2)=-2.*TU(2,2)*TU(2,1)
      TTX(3,3)=TU(2,2)*TU(2,2)
      TTY(1,1)=TU(3,3)*TU(3,3)
      TTY(1,2)=-2.*TU(3,3)*TU(3,4)
      TTY(1,3)=TU(3,4)*TU(3,4)
      TTY(2,1)=-TU(3,3)*TU(4,3)
      TTY(2,2)=1.+2.*TU(3,4)*TU(4,3)
      TTY(2,3)=-TU(3,4)*TU(4,4)
      TTY(3,1)=TU(4,3)*TU(4,3)
      TTY(3,2)=-2.*TU(4,4)*TU(4,3)
      TTY(3,3)=TU(4,4)*TU(4,4)
      GX=(1+AX*AX)/BX
      GY=(1+AY*AY)/BY
      BX1=TTX(1,1)*BX+TTX(1,2)*AX+TTX(1,3)*GX
      AX1=TTX(2,1)*BX+TTX(2,2)*AX+TTX(2,3)*GX
*     GX1=TTX(3,1)*BX+TTX(3,2)*AX+TTX(3,3)*GX
*     GX1=(1+AX1*AX1)/BX1
      BY1=TTY(1,1)*BY+TTY(1,2)*AY+TTY(1,3)*GY
      AY1=TTY(2,1)*BY+TTY(2,2)*AY+TTY(2,3)*GY
*     GY1=TTY(3,1)*BY+TTY(3,2)*AY+TTY(3,3)*GY
*     GY1=(1+AY1*AY1)/BY1
      EX1=TU(1,1)*EX+TU(1,2)*EPX+VU(1)
      EPX1=TU(2,1)*EX+TU(2,2)*EPX+VU(2)
      EY1=TU(3,3)*EY+TU(3,4)*EPY+VU(3)
      EPY1=TU(4,3)*EY+TU(4,4)*EPY+VU(4)
      AX=AX1
      AY=AY1
      BX=BX1
      BY=BY1
      EX=EX1
      EPX=EPX1
      EY=EY1
      EPY=EPY1
      DX1=matmul(TO,DX)
c      CALL VMUL4(TO,DX,DX1)
      DX(1)=DX1(1)+VO(1)*DP0
*    +      DP0*(R0I(1,1)*EX+R0I(1,2)*EPX+R0I(1,3)*EY+R0I(1,4)*EPY)
      DX(2)=DX1(2)+VO(2)*DP0
*    +      DP0*(R0I(2,1)*EX+R0I(2,2)*EPX+R0I(2,3)*EY+R0I(2,4)*EPY)
      DX(3)=DX1(3)+VO(3)*DP0
*    +      DP0*(R0I(3,1)*EX+R0I(3,2)*EPX+R0I(3,3)*EY+R0I(3,4)*EPY)
      DX(4)=DX1(4)+VO(4)*DP0
*    +      DP0*(R0I(4,1)*EX+R0I(4,2)*EPX+R0I(4,3)*EY+R0I(4,4)*EPY)
      end associate
      END subroutine

      subroutine QUADRU(TETA1,DX,DY)
      implicit none
      real*8 ,intent(in)::TETA1,DX,DY
      real*8 RK,P,S1,S2,S3,S4,TET2,DR,SQRK
      CALL ZEROV4(VM)
      CALL UNITM4(TM,1.d0)
      RK=DKL1/DL
      IF(RK.GE.0.) SQRK=DSQRT(RK)
      IF(RK.LT.0.) SQRK=DSQRT(-RK)
      P=SQRK*DL
      S1=DSIN(P)
      S2=DCOS(P)
      S3=DSINH(P)
      S4=DCOSH(P)
      IF (RK.GE.0.) THEN
*     QF TYPE
        TM(1,1)=S2
        TM(1,2)=S1/SQRK
        TM(2,1)=-SQRK*S1
        TM(2,2)=S2
        TM(3,3)=S4
        TM(3,4)=S3/SQRK
        TM(4,3)=SQRK*S3
        TM(4,4)=S4
        ELSE
*     QD TYPE
        TM(1,1)=S4
        TM(1,2)=S3/SQRK
        TM(2,1)=SQRK*S3
        TM(2,2)=S4
        TM(3,3)=S2
        TM(3,4)=S1/SQRK
        TM(4,3)=-SQRK*S1
        TM(4,4)=S2
      END IF
      TET2=TETA1*2.
      CKL0Y=-DKL1*(DX*DSIN(TET2)+DY*DCOS(TET2))
      CKL0X=DKL1*(DX*DCOS(TET2)-DY*DSIN(TET2))
*     CKL1X=DKL1*DCOS(TETA1)
*     CKL1Y=-CKL1X
      DR=DSQRT(DX*DX+DY*DY)
      IF(DR.GT.0.) THEN
        CKL1X=DKL1*DX/DR
        CKL1Y=-DKL1*DY/DR
      ELSE
        CKL1X=DKL1
        CKL1Y=-DKL1
      END IF
      CALL TROT(TM,VM,TETA1,TO,VO)
      RETURN
      END subroutine

      subroutine R0CAL(R0,R0I)
      use tffitcode
      implicit none
      real*8 RA0(2,2),RB0(2,2),RC0(2,2),RD0(2,2)
      real*8 R0(4,4),R0I(4,4)
      real*8 DETR,RMU0
      integer*4 J,K
c
      DO J=1,2
        DO K=1,2
          RC0(J,K)=TWPMG(10+K+2*(J-1))
        enddo
      enddo
      DETR=RC0(1,1)*RC0(2,2)-RC0(1,2)*RC0(2,1)
      IF(DETR.LE.1.) THEN
         RMU0=DSQRT(1.-DETR)
         CALL SYMTR2(RC0,RB0)
         CALL UNITM2(RA0,RMU0)
         CALL UNITM2(RD0,RMU0)
         CALL CM24(RA0,RB0,RC0,RD0,R0)
         CALL R1INV(R0,R0I)
      ELSE
*    DET(R)>1 CASE
         RA0(1,1)=TWPMG(13)/DETR
         RA0(1,2)=TWPMG(14)/DETR
         RA0(2,1)=TWPMG(11)/DETR
         RA0(2,2)=TWPMG(12)/DETR
         CALL SYMTR2(RA0,RD0)
         RMU0=DSQRT(1.-1./DETR)
         CALL UNITM2(RB0,RMU0)
         CALL UNITM2(RC0,RMU0)
         CALL CM24(RA0,RB0,RC0,RD0,R0)
         CALL R2INV(R0,R0I)
      ENDIF
      RETURN
      END subroutine

      subroutine TEDGE(DIEG)
      use tffitcode
      implicit none
      real*8 DIEG(*)
      real*8 R0(4,4),R0I(4,4)
      real*8 CURVX,CURVY
      associate (
     $     AU=>TWPMG(mfitax),
     $     BU=>TWPMG(mfitbx),
     $     PU=>TWPMG(mfitnx),
     $     AV=>TWPMG(mfitay),
     $     BV=>TWPMG(mfitby),
     $     PV=>TWPMG(mfitny),
     $     EU=>TWPMG(mfitex),
     $     EPU=>TWPMG(mfitepx),
     $     EV=>TWPMG(mfitey),
     $     EPV=>TWPMG(mfitepy),
     $     R1=>TWPMG(mfitr1),
     $     R2=>TWPMG(mfitr2),
     $     R3=>TWPMG(mfitr3),
     $     R4=>TWPMG(mfitr4),
     $     DETR=>TWPMG(mfitdetr),
     $     DU=>TWPMG(mfitdx),
     $     DPU=>TWPMG(mfitdpx),
     $     DV=>TWPMG(mfitdy),
     $     DPV=>TWPMG(mfitdpy))
      CURVX=CKL0X/DL
      CURVY=CKL0Y/DL
      CALL R0CAL(R0,R0I)
      DIEG(4)=DIEG(4)-(TO(2,1)*CURVX*(R0I(1,1)*EU+R0I(1,2)*EPU)
     +          -TO(4,3)*CURVY*(R0I(3,1)*EU+R0I(3,2)*EPU))
      DIEG(7)=DIEG(7)-(TO(2,1)*CURVX*(R0I(1,3)*EV+R0I(1,4)*EPV)
     +          -TO(4,3)*CURVY*(R0I(3,3)*EV+R0I(3,4)*EPV))
      DIEG(9)=DIEG(9)-TM(2,1)*BU
      DIEG(10)=DIEG(10)-TM(4,3)*BV
      CALL TMATRS(R0,R0I)
      CALL TWTRANS()
      end associate
      END subroutine

      subroutine BEDGE(EANG)
      implicit none
      real*8 EANG,CURV
      CURV=DKL0/DL
      CALL ZEROV4(VM)
      CALL UNITM4(TM,1.d0)
      TM(2,1)=CURV*DTAN(EANG)
      TM(4,3)=-CURV*DTAN(EANG)
      CALL TROT(TM,VM,TETA,TO,VO)
      CKL0X=DCOS(TETA)*DKL0
      CKL0Y=-DSIN(TETA)*DKL0
      CKL1X=0.
      CKL1Y=0.
      RETURN
      END subroutine

      subroutine SEXTU(TETA1,DX,DY)
      implicit none
      real*8 TETA1,DX,DY,SINT,COST
      CALL UNITM4(TM,1.d0)
      CALL UNITM4(TO,1.d0)
      CALL ZEROV4(VM)
      CALL ZEROV4(VO)
      SINT=DSIN(3.*TETA1)
      COST=DCOS(3.*TETA1)
      CKL0Y=-DKL2*(DX*DY*COST           +(DX*DX-DY*DY)/2.*SINT)
      CKL0X=DKL2*((DX*DX-DY*DY)/2.*COST-DX*DY*SINT)
*     CKL1X=DKL2*(2.*(DX-DY)*COST-2.*DX*SINT)
*     CKL1Y=-CKL1X
        CKL1X=DKL2*DX
        CKL1Y=-DKL2*DY
      RETURN
      END subroutine

      subroutine TROT (TM,VM,TETA,TO,VO)
      implicit none
      real*8 TM(4,4),VM(4),TO(4,4),VO(4)
      real*8 ROT(4,4),RINV(4,4),TS(4,4)
      real*8 TETA,CTET,STET
      CTET=DCOS(TETA)
      STET=DSIN(TETA)
      CALL UNITM4(ROT,CTET)
      CALL UNITM4(RINV,CTET)
      ROT(1,3)=-STET
      ROT(2,4)=-STET
      ROT(3,1)=STET
      ROT(4,2)=STET
      RINV(1,3)=STET
      RINV(2,4)=STET
      RINV(3,1)=-STET
      RINV(4,2)=-STET
      TS=matmul(TM,ROT)
      TO=matmul(RINV,TS)
      VO=matmul(RINV,VM)
c      CALL VMUL4(RINV,VM,VO)
      RETURN
      END subroutine

      end module
