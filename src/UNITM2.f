      subroutine UNITM2(A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4)
      A=0.d0
      A(1)=X
      A(4)=X
c$$$      CALL ZEROM2(A)
c$$$      DO 10 J=1,2
c$$$ 10   A(1+3*(J-1))=X
      RETURN
      END
