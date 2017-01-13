      subroutine TMATRS(R0,R0I)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer*4 ntwissfun
      parameter (ntwissfun=20)
      COMMON /TWISSP/TWPIN(ntwissfun),TWPOT(ntwissfun),TWPMG(ntwissfun)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      DIMENSION A(2,2),B(2,2),C(2,2),D(2,2),AI(2,2),CAI(2,2),ACI(2,2)
      DIMENSION RA1(2,2),RB1(2,2),RC1(2,2),RD1(2,2)
      DIMENSION TR(4,4),R0(4,4),R0I(4,4)
      CALL MMUL4(TO,R0I,TR)
*     WRITE(6,'(8F10.4)') R0
*     WRITE(6,'(8F10.4)') R0I
*     WRITE(6,'(8F10.4)') TR
      CALL DM42(TR,A,B,C,D)
      CALL MINV2(A,AI)
      CALL MMUL2(C,AI,CAI)
      DETCAI=CAI(1,1)*CAI(2,2)-CAI(1,2)*CAI(2,1)
      IF(DETCAI.GT.-1.) THEN
         SMU1=1./(1.+DETCAI)
         RMU1=DSQRT(SMU1)
         CALL FMUL2(-RMU1,CAI,RC1)
         CALL SYMTR2(RC1,RB1)
         CALL UNITM2(RA1,RMU1)
         CALL UNITM2(RD1,RMU1)
         CALL CM24(RA1,RB1,RC1,RD1,R0)
         CALL R1INV(R0,R0I)
         DO 20 J=1,2
         DO 20 K=1,2
  20       TWPMG(10+K+2*(J-1))=RC1(J,K)
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
         DO 21 J=1,2
         DO 21 K=1,2
  21       TWPMG(10+K+2*(2-J))=-RD1(J,K)
      END IF
      CALL MMUL4(R0,TR,TU)
*     WRITE(6,'(A)') ' R1,TU'
*     WRITE(6,'(8F10.4)') R0
*     WRITE(6,'(8F10.4)') TU
      CALL VMUL4(R0,VO,VU)
      RETURN
      END
