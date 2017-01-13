      subroutine ZEROM4(A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(16)
      DO 10 I=1,16
10    A(I)=0.D0
      RETURN
      END
