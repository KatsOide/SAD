      subroutine ZEROV4(A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4)
      DO 10 I=1,4
10    A(I)=0.D0
      RETURN
      END
