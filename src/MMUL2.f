      subroutine MMUL2(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,2), B(2,2),C(2,2)
      DO 10 I=1,2
      DO 10 J=1,2
      C(I,J)=0.
      DO 10 K=1,2
      C(I,J)=C(I,J)+A(I,K)*B(K,J)
10    CONTINUE
      RETURN
      END
