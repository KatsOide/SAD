      subroutine doread(dummy)
      use tfrbuf
      use maccbk
      use macttyp
      use macfile
      use macmisc
      use tfcsi, only:ipoint,lrecl
      use ffsfile
      implicit none
      character*(MAXSTR) token
c      character*(22) cfmsg
      integer slen,ival,ttype
      real*8 rval
      logical*4 skipch,ok
      integer*4 f,maxstk,ifd1
      integer*8 flmgr,itemp,kfile,ksize,mapallocfd
      parameter (maxstk=255)
      integer*4 fstk(maxstk),stkpt,ipak,i,dummy,irtc
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
      if(f .gt. 0 .and. f .lt. 99 .and. f .ne. infl)then
        if(lfnp .ge. maxlfn)then
          call errmsg('doread','lfn stack over flow',0,0)
        endif
      endif
      if((token(:slen) .eq. RCURL) .or.
     &     (token(:slen) .eq. SEMIC))    go to 3000
c                                              EXIT
      if (ttype .eq. ttypID) then
         if ((token(:3) .eq. 'STD') .or. (token(:3) .eq. 'std')) then
            if(ok) f=STDIN
         else
            if(ok) then
               f=nextfn(moderead)
               if (f .ne. 0)
     $              open(f,file=token(:slen),status='OLD',err=9000)
c     for debug
            print *,'input file is redirected to ',token(:slen)
c     end debug
            endif
         endif
      else if(ttype .eq. ttypST) then
        if(ok) then
          f=nextfn(moderead)
          open(f,FILE=token(:slen),STATUS='OLD',err=9000)
          ifd1=fnum(f)
          kfile=mapallocfd(ifd1,ksize,irtc)
          if(irtc .eq. 0)then
            call irbopen1(f,kfile/8,ksize+modemapped,ifd1)
            call trbassign(f)
            ipoint=1
            lrecl=0
          else
            write(*,*)'???-input file open error: ',
     $           token(:slen)
          endif
c     if( f .ne. 0)
c     $           open(f,file=token(:slen),status='OLD',err=9000)
c     for debug
          print *,'input file is redirected to #',f,token(:slen)
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
c      write(*,*)'doread-return ',f,infl
      if(f .ne. 0 .and. f .ne. infl)then
        lfnp=lfnp+1
        lfnstk(lfnp)=infl
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
      use tfrbuf, only:nextfn,moderead
      use macfile
      character*(MAXLLEN) fn
      integer fd
      integer*8 flmgr,f
c
      fd=nextfn(moderead)
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
 
