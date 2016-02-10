      subroutine QUADRU(TETA,DX,DY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      COMMON /DELEM/DL,DKL0,DKL1,DKL2,CKL0X,CKL0Y,CKL1X,CKL1Y
c
      save
c
      CALL ZEROV4(VM)
      CALL UNITM4(TM,1.)
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
      TET2=TETA*2.
      CKL0Y=-DKL1*(DX*DSIN(TET2)+DY*DCOS(TET2))
      CKL0X=DKL1*(DX*DCOS(TET2)-DY*DSIN(TET2))
*     CKL1X=DKL1*DCOS(TETA)
*     CKL1Y=-CKL1X
      DR=DSQRT(DX*DX+DY*DY)
      IF(DR.GT.0.) THEN
        CKL1X=DKL1*DX/DR
        CKL1Y=-DKL1*DY/DR
      ELSE
        CKL1X=DKL1
        CKL1Y=-DKL1
      END IF
      CALL TROT(TM,VM,TETA,TO,VO)
      RETURN
      END
