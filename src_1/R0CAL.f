      subroutine R0CAL(R0,R0I)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer*4 ntwissfun
      parameter (ntwissfun=20)
      COMMON /TWISSP/TWPIN(ntwissfun),TWPOT(ntwissfun),TWPMG(ntwissfun)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      DIMENSION RA0(2,2),RB0(2,2),RC0(2,2),RD0(2,2)
      DIMENSION R0(4,4),R0I(4,4)
c
      save
c
      DO 10 J=1,2
      DO 10 K=1,2
  10  RC0(J,K)=TWPMG(10+K+2*(J-1))
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
      END
