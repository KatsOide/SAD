      subroutine R2INV(R0,R0I)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R0(4,4),R0I(4,4)
      R0I(1:2,1:2)=-R0(3:4,3:4)
      R0I(3:4,3:4)=-R0(1:2,1:2)
      R0I(1:2,3:4)= R0(1:2,3:4)
      R0I(3:4,1:2)= R0(3:4,1:2)
c$$$      DO 10 I=1,2
c$$$      DO 10 J=1,2
c$$$         R0I(I,J)=-R0(I+2,J+2)
c$$$         R0I(I+2,J+2)=-R0(I,J)
c$$$         R0I(I,J+2)=R0(I,J+2)
c$$$  10     R0I(I+2,J)=R0(I+2,J)
      RETURN
      END
