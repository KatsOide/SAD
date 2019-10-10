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
      END
