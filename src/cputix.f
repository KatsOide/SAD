      subroutine cputix
      implicit none
      include 'inc/MACFILE.inc'
      real*8 ctime
      real*8 ctime0,dt
      integer*4 mfalloc,irtn
      data ctime0/0.0d0/
      save ctime0
c
      call cputime(ctime,irtn)
      dt=(ctime-ctime0)*1.0d-6
      write(errfl,100)
     $     ctime*1.d-6,dt,
     &      mfalloc(-1)
  100 format(1h ,
     &         ' cpu time=',1pg11.4,'(sec) dt=',3pf11.3,'(msec)',
     &         ' free area::',i8)
      return
      end

      subroutine cputime(ctime,irtc)
      implicit none
      real*8 ctime,second
      integer irtc
      ctime=second()*1.d6
      irtc=0
      return
      end
