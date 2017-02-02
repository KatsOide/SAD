      subroutine TEDGE(DIEG)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer*4 ntwissfun
      parameter (ntwissfun=20)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      COMMON /TWISSP/TWPIN(ntwissfun),TWPOT(ntwissfun),
     +   AU,BU,PU,AV,BV,PV,EU,EPU,EV,EPV,R1,R2,R3,R4,DU,DUP,DV,DVP
      COMMON /DELEM/DL,DKL0,DKL1,DKL2,CKL0X,CKL0Y,CKL1X,CKL1Y
      DIMENSION DIEG(*)
      DIMENSION R0(4,4),R0I(4,4)
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
      CALL TWTRANS(R0,R0I)
      END
