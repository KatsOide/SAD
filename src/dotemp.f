      subroutine dotemp(dummy)
      use maccbk
      implicit none
      integer*4 dummy

      include 'inc/CBKMAC.inc'
      character*(MAXSTR) token
      integer slen,ival,ttype,hsrchz
      real*8 rval
      logical skipch,skipped
      integer idx

c for debug
      call ptrace('dotemp',1)
c end debug
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .eq. ttypID) then
c         call errmsg('dotemp',
c    &              'missing '''//LCURL//''' is assumed',0,0)
        else
          call errmsg('dotemp',
     &      ' syntax error:temp{<line name>}',
     &                0,16)
        endif
      endif
c
 2000 continue
        skipped=skipch(COMMA,token,slen,ttype,rval,ival)
        if((token(:slen) .eq. RCURL) .or.
     &     (token(:slen) .eq. SEMIC)) go to 3000
        if (ttype .eq. ttypID) then
          idx=hsrchz(token(:slen))
          if (idtype(idx) .eq. icLINE) then
            if (ilist(2,idval(idx)) .le. 0) then
              call expnln(idx)
            endif
              call temp(ilist(2,idval(idx)))
          else
            call errmsg('dotemp',
     &                   token(:slen)//'is not line name.',0,0)
            call errmsg('dotemp',
     &                   token(:slen)//'is ignored.',0,0)
          endif
        else
          call errmsg('dotemp',
     &         'syntax error:illegal token as input',0,0)
          call errmsg('dotemp',
     &         'syntax:tempics{<linename>}',0,16)
        endif
        call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
 3000 continue
c for debug
      call ptrace('dotemp',-1)
c end debug
      return
      end
