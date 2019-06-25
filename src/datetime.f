      subroutine datetime(dat,tim)
      implicit none
      character*(*) dat,tim
      character*8 date
      character*10 time
      character*5 zone
      integer value(8)
      call Date_and_Time(date,time,zone,value)
      dat=date
      tim=time(1:6)
      return
      end

      subroutine datetime1(value)
      implicit none
      character*8 date
      character*10 time
      character*5 zone
      integer value(8)
      call Date_and_Time(date,time,zone,value)
      return
      end
