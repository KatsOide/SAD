      subroutine preabuf(outc,in,stat)
      logical eof,err,nextword,stat
      character*132 buf
      character*(*) outc
      data nextword/.false./
      save nextword,ipo
c
      eof=.false.
      err=.false.
c----- get next word or line
      if(stat) goto 1
      if(nextword) then
        goto 2
      else
        goto 1
      endif
c----- get one line
 1    call readstr(in,buf,irtc)
      if(irtc .ne. 0)then
        if(irtc .lt. 0)then
          go to 991
        else
          go to 990
        endif
      endif
c      read(in,'(a)',end=991,err=990) buf
c     ... neglect the part preceded by '!' ....
      ipo=ifany(buf,'!',1)
      if(ipo.ne.0) then
        if(ipo.gt.1) then
          buf=buf(1:ipo-1)
        else
          goto 1
        endif
      endif
c----- get word
      l=lene(buf)
      ipo=1
 2    ipo=notany(buf,' ',ipo)
      if(ipo.ne.0) then
        is=ipo
        ipo=ifany(buf,' ',ipo)
        if(ipo.ne.0) then
          if=ipo
          nextword=.true.
        else
          if=l+1
          nextword=.false.
        endif
      else
        goto 1
      endif
c
      outc=buf(is:if-1)
      stat=.false.
      return
 990  err=.true.
      stat=.true.
      return
 991  eof=.true.
      stat=.true.
      return
      end
