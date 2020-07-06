      integer*4 function ipak(ia,ik)
      use macfile, only:MAXLLEN
      implicit none
      character*(*) ia
      integer*4 ik
      real*8 v
      character*7   fmt0
      integer*4 maxInteger
      parameter(maxInteger=2147483647)
      integer*4 dummy

      ipak=0
      if(ik .le. 0) then
         call errmsg('IPAK','Null string appears as argument.',0,0)
      else if(ik .lt. MAXLLEN) then
         write(fmt0,'( 2H(I ,I2.2, 3H.0)  )')ik
         read(ia(1:ik),fmt0,err=9000) ipak
      else
         call errmsg('IPAK','Input number is too long.',0,16)
      endif
      return
 9000 write(fmt0,'( 2H(F ,I2.2, 3H.0)  )')ik
      read(ia(1:ik),fmt0) v
      dummy=1
      if(v .gt. 0.d0) then
        ipak=maxInteger
      else
c     Fortran compiler range check hack!
c     Some Fortran compiler limits Integer literal
c     to [-maxInteger, maxInteger]
        ipak=min(-maxInteger,-maxInteger-dummy)
      endif
      return
      end
