      subroutine SYMTR2(A,B)
*  B=-SI*A*S
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,2),B(2,2)
      B(1,1)=-A(2,2)
      B(1,2)=A(1,2)
      B(2,1)=A(2,1)
      B(2,2)=-A(1,1)
      RETURN
      END
