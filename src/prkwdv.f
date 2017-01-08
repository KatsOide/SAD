      subroutine prkwdv(head,kwd,val)
      use maccbk
      use macfile
      implicit none
      character head*(*),kwd*(*)
      real*8 val
      integer*4 chrcnt,lene
      integer*4 llen,idxl,slen,len,i
      character*80 line
      data chrcnt/0/
      save chrcnt
c
      if (val .eq. 0.d0) return
c
      slen=lene(kwd)
      if (chrcnt+slen .gt. 80-23) then
        write(outfl,*) ' ',head(1:lene(head)),line(1:chrcnt)
        chrcnt=0
      endif
      write(line(chrcnt+1:chrcnt+slen+23),
     &      '(A,''='',1PG21.14,'' '')')kwd(:slen),val
      chrcnt=chrcnt+slen+23
      return
c
      entry prkwnm(head,kwd)
c
      slen=lene(head)
      slen=slen+lene(kwd)
      slen=slen+3
      if (chrcnt+slen .gt. 80) then
          write(outfl,*) line(1:chrcnt)
          chrcnt=0
      endif
c      write ( line(chrcnt+1:),"(A,A,A,A)")
c     &        head(1:lene(head))," ",kwd(1:lene(kwd)),"=("
c      chrcnt=chrcnt+slen
      do i=1,lene(head)
         chrcnt=chrcnt+1
         line(chrcnt:chrcnt)=head(i:i)
      end do
      chrcnt=chrcnt+1
      line(chrcnt:)=" "
      do i=1,lene(kwd)
         chrcnt=chrcnt+1
         line(chrcnt:)=kwd(i:i)
      end do
      chrcnt=chrcnt+1
      line(chrcnt:)="="      
      chrcnt=chrcnt+1
      line(chrcnt:)="("      
c
      return
c
      entry prkw(kwd)
c      slen=lene(kwd)
      slen=len(kwd)
      if (chrcnt+slen .gt. 80) then
          write(outfl,*) line(:chrcnt)
          chrcnt=0
      endif
      write (line(chrcnt+1:),'(A)') kwd(:slen)
      chrcnt=chrcnt+slen
      return
c
      entry prkwot(kwd)
      if (chrcnt+len(kwd) .gt. 80) then
        write(outfl,'(A)') line(:chrcnt)
        write(outfl,'(A)') kwd
      else
        write(outfl,'(A,A)') line(:chrcnt),kwd
      endif
      chrcnt=0
      return
c
      Entry Lprnum(llen,idxl,head)
c
      do 1000  i=0,llen-1
        if (chrcnt+11 .gt. 80) then
           write(outfl,*) line(:chrcnt)
           chrcnt=0
        endif
        write(line(chrcnt+1:chrcnt+11),
     &        '(1PG10.4,'','')')rlist(i+idxl)
        chrcnt=chrcnt+11
 1000 continue
      if (chrcnt .gt. 1) then
        chrcnt=chrcnt-1
c       write(outfl,*) line(:chrcnt-1)
c       chrcnt=0
      endif
      return
c
      end
