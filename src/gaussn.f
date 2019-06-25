      real*8 function gaussn(IG0)
C  GENERATE GAUSSIAN RANDUM DEVIATES WITH ZERO MEAN AND UNIT S.D.
C   INPUT: IG0=ODD INTEGER
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 xx
      parameter (xx=6.550464984D-19/3.063100571D-00)
      INTEGER*4 IRC
      data IRC /48828125/
C
 1000 continue
      IF (IG0 .EQ. 0) go to 9999
        IG0=IRc*IG0
        IGA=ABS(IG0)
        Y=DBLE(IGA)
        IF(IGA .LT. 1288490189) THEN
          YY=Y*Y
          gaussn=((xx*YY+1.d0)*3.06310057134D-29*YY
     %            +5.84190096639D-10)*Y
        ELSE IF(IGA .LT. 1997159793) THEN
          gaussn=((2.15247545000D9-Y)*Y+5.13646191373D17)
     %      /((0.448498722421D0*Y-3.55927532674D9)*Y+5.77424700846D18)
        ELSE IF(IGA .LT. 2126008812) THEN
          gaussn=((Y-5.53451026774D9)*Y+7.30305470154D18)
     %      /((1.68096862343D0*Y-7.78160540942D9)*Y+8.96788425914D18)
        ELSE
          T=DSQRT(-DLOG((0.5-2.3283064365387D-10*Y)**2))
          gaussn=(T-(T+7.0309866D0)
     %                /((0.18604881D0*T+3.0556532D0)*T+3.2217897D0))
        ENDIF
        IF (IG0 .LT. 0) gaussn=-gaussn
      RETURN
 9999 continue
        call errmsg('gaussn',
     &              'random seed becomes 0',0,0)
        ig0=irc**2+2
        go to 1000
      END
