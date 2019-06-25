      subroutine FMUL2(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(2,2),C(2,2)
c
      save
c
      DO 10 I=1,2
      DO 10 J=1,2
      C(I,J)=A*B(I,J)
10    CONTINUE
      RETURN
      END
