      module version
        character*19, parameter ::
c                     /'1234567890123456789'/
     $     versionid  ='1.1.9.0.7k64       ',
     $     versiondate='9/13/2020 22:00:00 '

        character*25 builtdate
        character*30 startdat
      end module

      program MAIN
      use version
      use maccbk
      use tfmem
      use tftok
      implicit none
c
      call fdate1(startdat)
      call buildinfo_get_string('Built:Date', builtdate)
c
      write(*,*)
     $     '*** Welcome to SAD Ver.',versionid(1:len_trim(versionid)),
     $     ' built at ',builtdate(1:len_trim(builtdate)),' ***'
      write(*,*)'*** Today: ',startdat(1:len_trim(startdat)),' ***'
c
      call talocinit
      call inifil
      call initbl
      call tftokinit
      call ktfinitshare
c
      call toplvl
      end

      subroutine fdate1(dat)
      implicit none
      character*(*) ,intent(out):: dat
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
