      function SYMPS(Z,N)
      PARAMETER (ML=30)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(ML)
      IF(N.EQ.0) THEN
         SYMPS=Z(1)
         RETURN
      END IF
      S2=0.
      S3=0.
      DO 10 I=1,N-1
        S2=S2+Z(I*2+1)
10    CONTINUE
      DO 20 I=1,N
        S3=S3+Z(I*2)
20    CONTINUE
      SYMPS=(Z(1)+Z(2*N+1)+2.*S2+4.*S3)/3.
      END
