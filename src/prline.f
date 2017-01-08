      subroutine prline(idxl,head)
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      implicit real*8 (a-h,o-z)
      integer idxl
      character*(*) head
c
      character    head1*6,cwork*80
      integer slen,llen
      integer*8 listp,iw
      character*17 fmt
      data fmt /'(A**,A4,A8,:,A20)'/
c
      slen=len(head)
      write(fmt(3:4),'(I2.2)') slen
      cwork(:slen)=head
      llen=ilist(1,idval(idxl))
      listp=idval(idxl)
      write(outfl,'(1H ,A8,''made of '',I3,'' elements=('')')
     &        pname(idxl),llen
      do 100 iw=listp+1,listp+llen
        if (ilist(1,iw) .ne. 1) then
          write(head1,'(I3,1H*)')
     &               ilist(1,iw)
        else
          head1=' '
        endif
        if (idtype(ilist(2,iw)) .eq. icLINE) then
            write(outfl,fmt) head,head1(:4)
     &                         ,pname(ilist(2,iw)),' '
        else if (idtype(ilist(2,iw)) .eq. icNULL) then
          write(outfl,fmt) head,head1(:4)
     &                  ,pname(ilist(2,iw)),':not defined yet'
        else if (idtype(ilist(2,iw)) .lt. icMXEL) then
          cwork(slen+1:slen+4)=head1
          call prelem(ilist(2,iw),cwork(:slen+4))
        endif
 100  continue
      write(outfl,'(T9,'')'' )')
      return
      end
