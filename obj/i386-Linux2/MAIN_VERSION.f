      implicit none
      include 'inc/MACFILE.inc'
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACCBK.inc'
c
      integer*4 lnblnk
      character*30 dat
c
      character*19 versionid,versiondate
      character*25 builtdate
      common /version/ versionid,versiondate
c                      /'         1111111111'/
c                      /'1234567890123456789'/
      data versionid   /'1.0.10.4.9a20      '/
      data versiondate /'12/28/2010 23:00:00'/
      data builtdate   /'2010-12-28 17:08:30 -0800'/
c
      call fdate1(dat)
      write(*,*)
     $     '*** Welcome to SAD Ver.',versionid(1:lnblnk(versionid)),
     $     ' built at ',builtdate(1:lnblnk(builtdate)),' ***'
      write(*,*)'*** Today: ',dat(1:lnblnk(dat)),' ***'
c
      call inimem
      call inifil
      call initbl
c
      call toplvl
      end
