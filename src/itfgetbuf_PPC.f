      integer*4 function itfgetbuf(lfn,buf,limit,irtc)
      implicit none
      integer*4 lfn,limit,i,fgetc,irtc
      character buf(limit)
      do i=1,limit
        irtc=fgetc(lfn,buf(i))
        if(irtc .ne. 0)then
          itfgetbuf=i-1
          if(irtc .lt. 0)then
            if(i .gt. 1)then
              irtc=0
            else
              itfgetbuf=-99
            endif
          else
            itfgetbuf=-999
          endif
          return
        endif
        if(buf(i) .eq. char(10))then
          itfgetbuf=i-1
          if(i .gt. 1)then
            if(buf(i-1) .eq. char(13))then
              itfgetbuf=i-2
            endif
          endif
          return
        endif
      enddo
      itfgetbuf=limit
      return
      end
