      subroutine VMUL4(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4,4), B(4),C(4)
      DO 10 I=1,4
      C(I)=0.
      DO 10 K=1,4
      C(I)=C(I)+A(I,K)*B(K)
10    CONTINUE
      RETURN
      END
