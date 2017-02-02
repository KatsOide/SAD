      subroutine doread(dummy)
      use maccbk
      implicit none
      include 'inc/MACTTYP.inc'
      include 'inc/MACFILE.inc'
      include 'inc/MACMISC.inc'
      character*(MAXSTR) token
c      character*(22) cfmsg
      integer slen,ival,ttype
      real*8 rval
      logical*4 skipch,ok
      integer*4 f,flmgr,maxstk
      parameter (maxstk=255)
      integer*4 fstk(maxstk),stkpt,nextfn,ipak,i,itemp,dummy
c
c      data cfmsg/'read data from file=**'/
c     
      ok=.true.
      f=0
      stkpt=0
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
         if (.not. (ttype .eq. ttypNM .or. ttype .eq. ttypST .or.
     $        ttype .eq.ttypID)) then
            call errmsg('doread',
     &           'syntax error:'//LCURL//' is expected.',0,16)
         endif
      endif
 2000 if( skipch(COMMA,token,slen,ttype,rval,ival) ) go to 2000
      if(f .gt. 0 .and. f .lt. 99)then
         stkpt=stkpt+1
         if (stkpt .gt. maxstk) then
            call errmsg('doread',
     &           'stack over flow',0,0)
            stkpt=maxstk
         else
            fstk(stkpt)=f
         end if
         f=0
      endif
      if((token(:slen) .eq. RCURL) .or.
     &     (token(:slen) .eq. SEMIC))    go to 3000
c                                              EXIT
      if (ttype .eq. ttypID) then
         if ((token(:3) .eq. 'STD') .or. (token(:3) .eq. 'std')) then
            if(ok) f=STDIN
         else
            if(ok) then
               f=nextfn(0)
               if (f .ne. 0)
     $              open(f,file=token(:slen),status='OLD',err=9000)
c     for debug
            print *,'input file is redirected to ',token(:slen)
c     end debug
            endif
         endif
      else if(ttype .eq. ttypST) then
         if(ok) then
            f=nextfn(0)
            if( f .ne. 0)
     $           open(f,file=token(:slen),status='OLD',err=9000)
c     for debug
         print *,'input file is redirected to ',token(:slen)
c     end debug
         end if
      else if(ttype .eq. ttypNM) then
         if(ok) then
            f=ipak(token(:slen),slen)
            if((f .le. 0) .or. (f .ge. 99))then
               call errmsg('doread',
     &              'invalid file number '//token(:slen),0,0)
            endif
         endif
      else
         call errmsg('doread',
     &        'invalid input',0,0)
      endif
      call gettok(token,slen,ttype,rval,ival)
      go to 2000
c     
 3000 continue
      f=infl
      do i=stkpt,1,-1
         itemp=flmgr(f)
         f=fstk(i)
      end do
      if(f .ne. infl)then
         call myfflush
         infl=f 
      endif
c.....for debug
c     call ptrace('doread',-1)
c.....end debug
      return
c
 9000 continue
      call errmsg('doread','Cannot open file '//token(:slen),0,0)
      ok=.false.
      call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
      end
c
      subroutine redirectInput(fn,fd)
      character*255 fn
      integer fd, f,flmgr
      include 'inc/MACFILE.inc'
c
      fd=nextfn(fd)
      if (fd .eq. 0) then
         call errmsg('redirectInput','No more file descriptor',0,0)
         return 
      endif
      open(fd,file=fn,status='OLD',err=9000)
      f=flmgr(infl)
      infl=fd
      return
c
 9000 continue
      call errmsg('redirect','Cannot open file '//fn,0,0)
      return
      end
c     
      integer*4 function nextfn(i)
      use tfrbuf
      implicit none
      integer*4 i
c
      integer*4 ios,f, is,ie
      logical*4 od
      character*255 msg
c
      f=0
      is=11
      ie=nbuf
      do f=is,ie
         inquire(f,IOSTAT=ios,err=9000,OPENED=od)
         if( .not. od) then
            nextfn=f
            return
         endif
      end do
c
 1000 continue
      nextfn=0
      return
c
 9000 continue
      call perror(msg)
      go to 1000
c
      end
