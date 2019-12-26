      integer*4 function itfgetbuf(lfn,buf,limit,irtc)
      implicit none
      integer*4 lfn,limit,i,fgetc,irtc,nc
      character*(*) buf
c      write(*,*)'itfgetbuf ',lfn
      read(lfn,'(q,a)',end=9000,err=9100)nc,buf(:min(nc,limit))
      if(nc .lt. limit)then
        nc=nc+1
        buf(nc)=char(10)
      endif
      irtc=0
      itfgetbuf=nc
      return
 9000 irtc=-1
      itfgetbuf=-99
      return
 9100 irtc=1
      itfgetbuf=-999
      return
      end
