      subroutine UNITM2(A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4)
      CALL ZEROM2(A)
      DO 10 J=1,2
 10   A(1+3*(J-1))=X
      RETURN
      END
