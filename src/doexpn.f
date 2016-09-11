      subroutine doexpn(dummy)
      use maccbk
      use mackw
      use macttyp
      use macmisc
      implicit none
c
      character*(MAXSTR) token
      integer*4 slen,ival,ttype,hsrchz,idx,dummy
      real*8 rval
      logical skipch,skipped
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
         if (ttype .ne. ttypID) then
            call errmsg('doexpn',
     &           'syntax error:'//LCURL//' is expected.',0,16)
         endif
      endif
c     strat line expand
 2000 continue
      if (ttype .eq. ttypID) then
        idx =hsrchz(token(:slen))
        if (idtype(idx) .eq. icNULL) then
          call errmsg('doexpn',
     &         token(:slen)//' is not defined yet',0,0)
        else if (idtype(idx) .lt. icMXEL) then
          call errmsg('doexpn',
     &        token(:slen)//' is '//pname(kytbl(0,idtype(idx)))
     &        ,0,0)
          call errmsg('doexpn',
     &        'unable to expand.'
     &        ,0,0)
        else if (idtype(idx) .eq. icLINE) then
          call expnln(idx)
        else
          call errmsg('doepxn',
     &      'syntax error: '//token(:slen)//
     $         ' is already used as an element name.'
     &        ,0,16)
        endif
      else   if((token(:slen) .eq. RCURL) .or.
     &   (token(:slen) .eq. SEMIC))    then
        return
      else
        call errmsg('doexpn',
     &       'syntax error:expand{<line name> ,<line name>!}',0,16)
      endif
      call gettok(token,slen,ttype,rval,ival)
      skipped=skipch(COMMA,token,slen,ttype,rval,ival)
      go to 2000
c
      end
