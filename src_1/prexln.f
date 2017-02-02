      subroutine prexln(idxl,head)
      use maccbk
      implicit real*8 (a-h,o-z)
      integer idxl
      character*(*) head
c
      include 'inc/MACCODE.inc'
      include 'inc/MACFILE.inc'
c
      character*12 head1,cwork*80
      integer slen
c
c for debug
      call ptrace('prexln   '//pname(idxl)//'!',1)
c end debug
      if (ilist(2,idval(idxl)) .le. 0) then
        call errmsg('prexln',
     &       pname(idxl)//' is expanded now!',0,0)
        call expnln(idxl)
      endif
c
      slen=len(head)
      cwork(:slen)=head
      ip=ilist(2,idval(idxl))
      llen=ilist(1,ip)
      write(outfl,'(1H ,A8,'' expanded(length='',I4,'')!=('')')
     &             pname(idxl),llen
      do 100 iw=ip+1,ip+llen
          head1='            '
        if (idtype(ilist(1,iw)) .eq. icLINE) then
            write(outfl,*) head,head1
     &                         ,pname(ilist(1,iw))
        else if (idtype(ilist(1,iw)) .eq. icNULL) then
          write(outfl,*) head,
     &                   pname(ilist(1,iw)),':not defined yet'
        else if (idtype(ilist(1,iw)) .lt. icMXEL) then
          cwork(slen+1:slen+12)=head1
c         call prelem(ilist(1,iw),cwork(:slen+12))
          call prelm0(ilist(1,iw),ilist(2,iw),cwork(:slen+12))
        else
c.........for debug
c         print *,'element ',pname(ilist(1,iw)),'!'
c.........end debug
          call errmsg('prexln',
     &         'Woom, something wrong occured?',0,0)
          call errmsg('prexln',
     &         pname(ilist(1,iw)),0,0)
        endif
 100  continue
      write(outfl,'(T9,'')'' )')
c for debug
      call ptrace('prexln '//pname(idxl)//'!',-1)
c end debug
      return
      end
