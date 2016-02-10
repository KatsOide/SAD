      subroutine CM24(A,B,C,D,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(4,4), A(2,2),B(2,2),C(2,2),D(2,2)
c
      save
c     
      DO 10 I=1,2
      DO 10 J=1,2
      R(I,J)=A(I,J)
      R(I,J+2)=B(I,J)
      R(I+2,J)=C(I,J)
      R(I+2,J+2)=D(I,J)
10    CONTINUE
      RETURN
      END
