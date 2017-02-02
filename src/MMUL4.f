      subroutine MMUL4(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4,4), B(4,4),C(4,4)
      DO 10 I=1,4
      DO 10 J=1,4
      C(I,J)=0.
      DO 10 K=1,4
      C(I,J)=C(I,J)+A(I,K)*B(K,J)
10    CONTINUE
      RETURN
      END
