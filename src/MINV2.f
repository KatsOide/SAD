      subroutine MINV2(A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,2),B(2,2)
      DETA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DETA.EQ.0.) THEN
           WRITE(6,'(A)') 'DETERMINANT=0 MATRIX-INV NOT FOUND'
           WRITE(6,'(4F10.4)') A
           STOP
      END IF
      B(1,1)=A(2,2)/DETA
      B(1,2)=-A(1,2)/DETA
      B(2,1)=-A(2,1)/DETA
      B(2,2)=A(1,1)/DETA
      RETURN
      END
