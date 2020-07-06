      subroutine DM42(R,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(4,4), A(2,2),B(2,2),C(2,2),D(2,2)
c
      save
c
      A=R(1:2,1:2)
      B=R(1:2,3:4)
      C=R(3:4,1:2)
      D=R(3:4,3:4)
c$$$      DO 10 I=1,2
c$$$      DO 10 J=1,2
c$$$      A(I,J)=R(I,J)
c$$$      B(I,J)=R(I,J+2)
c$$$      C(I,J)=R(I+2,J)
c$$$      D(I,J)=R(I+2,J+2)
c$$$10    CONTINUE
      RETURN
      END
