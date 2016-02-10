      subroutine UNITM4(A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(16)
      CALL ZEROM4(A)
      DO 10 J=1,4
 10   A(1+5*(J-1))=X
      RETURN
      END
