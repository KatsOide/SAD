      subroutine doprin(dummy)
      use maccbk
      implicit none
      real*8  dummy
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACFILE.inc'
      character*(MAXSTR) token
      integer*4 slen,ival,ttype,hsrchz,hsrch,lpname
      real*8 rval
      logical skipch
c     
      integer*4 idx,ioldf,i
c     
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
         if (ttype .eq. ttypID) then
c     call errmsg('doprint',
c     &              'missing '''//LCURL//''' is assumed',0,0)
         else
            call errmsg('doprint',
     &           ' syntax error:print{<line or element name',
     &           0,16)
         endif
      endif
c     
 2000 continue
      if((token(:slen) .eq. RCURL) .or.
     &     (token(:slen) .eq. SEMIC)) then
         return
      else if (ttype .eq. ttypID) then
         idx=hsrchz(token(:slen))
         call gettok(token,slen,ttype,rval,ival)
 2100    if (skipch(COMMA,token,slen,ttype,rval,ival) ) go to 2100
         if (idtype(idx) .eq. icNULL) then
            if (pname(idx)(:lpname(idx)) .eq. 'GLOBALS') then
               ioldf=outfl
               outfl=21
               do 2110 i=1,HTMAX
                  if((idtype(i) .ge. icGLI) .and.
     &                 (idtype(i) .le. icGLR))
     &                 call iprglb(i)
 2110          continue
               outfl=ioldf
            else if (pname(idx)(:lpname(idx)) .eq. 'FLAGS') then
               ioldf=outfl
               outfl=21
               do 2120 i=1,HTMAX
                  if(idtype(i) .eq. icFLAG)
     &                 call prnflg(pname(i))
 2120          continue
               outfl=ioldf
            else
               idx=hsrch(pname(idx)(:lpname(idx))//'$')
               if (idtype(idx) .eq. icFLAG) then
                  call prnflg(pname(idx))
               else
                  call errmsg('doprin',
     $                 'element name  '//
     $                 pname(idx)(:lpname(idx))//'! is not defined yet.'
     $                 ,0,4)
               endif
            endif
c     
         else if (idtype(idx) .eq. icLINE) then
            if (skipch('LONG',token,slen,ttype,rval,ival) ) then
               call prline(idx,' ')
            else if (skipch('EXPAND',token,slen,ttype,rval,ival) ) then
               call prexln(idx,' ')
            else if (skipch('SHORT',token,slen,ttype,rval,ival) ) then
               call sprexl(idx)
            else if (skipch('SAD',token,slen,ttype,rval,ival) ) then
               call prSad(idx)
            else if (skipch('MATRIX',token,slen,ttype,rval,ival) ) then
               call prntmt(idx)
            else if (skipch('COD',token,slen,ttype,rval,ival)
     &              .or. skipch('cod',token,slen,ttype,rval,ival) )then
               if(ilist(2,idval(idx)) .ne. 0) then
                  if (skipch('PLOT',token,slen,ttype,rval,ival) .or.
     $                 skipch('plot',token,slen,ttype,rval,ival) )then
                     call pltcod(ilist(2,idval(idx)))
                  else
                     call prncod(ilist(2,idval(idx)))
                  endif
               else
                  call errmsg('doprin',
     &                 pname(idx)(:lpname(idx))//
     $                 ' is not expanded yet.',0,4)
               endif
            else
               call sprlin(idx)
            endif
         else if (idtype(idx) .lt. icMXEL) then
            call prelem(idx,' ')
         else if ((idtype(idx) .ge. icGLI )
     &           .and.  (idtype(idx) .le. icGLR))      then
            call iprglb(idx)
         else
            idx=hsrch(pname(idx)(:lpname(idx))//'$')
            if (idtype(idx) .eq. icFLAG) then
               call prnflg(pname(idx))
            else
               call errmsg('doprin',
     &              'syntax error:non element name is given',0,4)
            endif
         endif
      else
         call errmsg('doprin',
     &        'syntax error:illegal token as input',0,0)
         call errmsg('doprin',
     &        'syntax: print{<name of line or element> , ...!}',0,16)
      endif
      go to 2000
c     
      end
