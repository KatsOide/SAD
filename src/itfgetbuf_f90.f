      integer*4 function itfgetbuf(lfn,buf,limit,irtc)
      implicit none
      integer*4 lfn,limit,irtc
      character(len=*) :: buf
      integer :: bufmax,nc,status
      bufmax = min(len(buf),limit)
      read(lfn,fmt='(a)',advance='no',end=9000,
     $     size=nc,iostat=status)buf(:bufmax)
      if(status .gt. 0)then
         irtc=1
         itfgetbuf=-999
         return
      endif
      irtc=0
      itfgetbuf=nc
      return
 9000 continue
      irtc=-1
      itfgetbuf=-99
      return
      end

c     CAUTION:
c     This fgetc implementation by non-advancing input
c     can't handle '\r' character correctly!
c
c     for Examples:
c     On some Fortran compiler, '\r\n' 2-character sequence returns
c     1-character '\n', because CR+LF('\r\n') is assumed as EOL.
c
c     Some broken fortran compiler stops non-advancing input at '\r'
c     same as '\n' detection and drops next 1-character.
c     Ex). un-seek-able input stream with gfortran 4.2.4 20080220 (prerelease)
c
      integer*4 function fgetc(lfn,buf)
      implicit none
      use ISO_C_BINDING
      integer*4 lfn
      character(len=*) :: buf
      integer :: nc,status
      read(lfn,fmt='(a1)',advance='no',end=9000,
     $     size=nc,iostat=status)buf
      if(status .gt. 0)then
         fgetc=1
         return
      endif
      if(status .lt. 0)then
         buf=C_NEW_LINE
      endif
      fgetc=0
      return
 9000 continue
      fgetc=-1
      return
      end
