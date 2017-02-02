      subroutine BMAG(TETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      COMMON /DELEM/DL,DKL0,DKL1,DKL2,CKL0X,CKL0Y,CKL1X,CKL1Y
c
      save
c
      CURV=DKL0/DL
      CALL ZEROV4(VM)
      CALL UNITM4(TM,1.)
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
         CKL0X=DCOS(TETA)*DKL0
         CKL0Y=-DSIN(TETA)*DKL0
         CKL1X=0.
         CKL1Y=0.
      CALL TROT(TM,VM,TETA,TO,VO)
      RETURN
      END
