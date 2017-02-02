      subroutine fdate1(dat)
      implicit none
      character*(*) dat
      character*8 d
      character*12 tim
      character*9 day
      integer lene

      call datetime(d,tim)
      if(d(7:7) .eq. ' ')then
        d(7:7)='0'
      endif

      call tday(day)
      dat=tim(1:2)//":"//tim(3:4)//":"//
     $     tim(5:6)//' '//day(1:lene(day))//' '//
     $     d(5:6)//"/"//d(7:8)//"/"//d(1:4)
      return
      end
