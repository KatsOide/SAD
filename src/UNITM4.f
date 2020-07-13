      subroutine UNITM4(A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(16)
      A=0.d0
      A(1)=X
      A(6)=X
      A(11)=X
      A(16)=X
c$$$      CALL ZEROM4(A)
c$$$      DO 10 J=1,4
c$$$ 10   A(1+5*(J-1))=X
      RETURN
      END
