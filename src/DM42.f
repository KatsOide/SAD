      subroutine DM42(R,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(4,4), A(2,2),B(2,2),C(2,2),D(2,2)
c
      save
c
      DO 10 I=1,2
      DO 10 J=1,2
      A(I,J)=R(I,J)
      B(I,J)=R(I,J+2)
      C(I,J)=R(I+2,J)
      D(I,J)=R(I+2,J+2)
10    CONTINUE
      RETURN
      END
