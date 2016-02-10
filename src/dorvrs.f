      subroutine dorvrs(dummy)
      use maccbk
      implicit none
      integer*4 dummy

      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
      character*(MAXSTR) token
      integer slen,ival,ttype,hsrchz,lpname,idx
      real*8 rval
      logical skipch
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .ne. ttypID) then
          call errmsg('dorvrs',
     &                ' syntax error:print{<line or element name',
     &                0,16)
        endif
      endif
c
 2000 continue
      if((token(:slen) .eq. RCURL) .or.
     &   (token(:slen) .eq. SEMIC)) then
        return
      else if (ttype .eq. ttypID) then
        idx=hsrchz(token(:slen))
        call gettok(token,slen,ttype,rval,ival)
 2100   if (skipch(COMMA,token,slen,ttype,rval,ival) ) go to 2100
        if (idtype(idx) .eq. icNULL) then
          call errmsg('dorvrs',
     &         'element name  '//
     $          pname(idx)(:lpname(idx))//
     $          '! is not defined yet.',0,0)
c
        else if (idtype(idx) .eq. icLINE) then
            call rvrsln(idx)
        else
          call errmsg('dorvrs',
     &         'syntax error:line name expected.'//
     $          pname(idx)(:lpname(idx))
     $          ,0,0)
        endif
      else
        call errmsg('dorvrs',
     &       'syntax error:illegal token as input',0,0)
        call errmsg('dorvrs',
     &       'syntax: print <Line name>;',0,16)
      endif
      go to 2000
c
      end
