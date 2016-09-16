      real*8 function atof(token,slen)
      implicit none
      character*(*) token
      integer*4 slen,m
      real*8 eval1
      integer*4 MAXLLEN,LLEN
      parameter(MAXLLEN=255)
      parameter(LLEN=MAXLLEN)
c
      atof=0.d0
      if (slen .gt. LLEN)  then
       call errmsg('atof',
     &             'Data is too long.',0,0)
       call errmsg('atof',token,0,16)
      else if (slen .le.0) then
       call errmsg('atof',
     &             'Null string appears as argument.',0,16)
      else       
c        write(fmt0,'( 2H(F ,I2.2, 3H.0)  )')slen
c        read(token(:slen),fmt0) atof
        atof=eval1(token,len(token),1,m)
      endif
      return
      end
