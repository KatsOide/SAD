      function Lrdnum(idxl,token,slen,ttype)
c     $Header: /SAD/cvsroot/oldsad/src/Lrdnum.f,v 1.2.14.2 2012/09/15 08:14:55 oide Exp $
      use maccbk
      implicit none
c
      include 'inc/MACCODE.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACTTYP.inc'
c
      integer*4 Lrdnum,scan, yyparse,yypushtoken,yypoptoken,yyreturn
      integer idxl,idxm
      character*(MAXSTR) token
      integer slen,ival,ttype
      real*8 rval,yyval
Cmacro cbkmac
c
      logical skipch,skiped
c
      Lrdnum=0
 1000 continue
      if (skipch(RPAR,token,slen,ttype,rval,ival)) return
c     
      skiped=skipch(COMMA,token,slen,ttype,rval,ival)
c for debug
c     call errmsg('Lrdnum',
c     &              'for debug :'//token(:slen),0,0)
c end debug
      if ( ttype .ne. ttypNM .and.
     $     ttype .ne. ttypID .and. ttype .ne. ttypDM) then
         call errmsg('Lrdnum',
     &        'math-expression is expected as input.('
     &        //token(:slen)//')'
     &        ,0,0)
 1100   continue
        if (skipch(RPAR,token,slen,ttype,rval,ival)) return
        go to 1100
      else
         Lrdnum=Lrdnum+1
         idxm=idxl+Lrdnum
         scan=yypushtoken(token(:slen),slen,ttype)
         scan=yyparse(yyreturn,yyval)
c     print *,' rdexpr: ', scan, yyreturn,yyval
         rlist(idxm)=yyval
         scan=yypoptoken(token,slen,ttype)
      endif
      go to 1000
c     
      end
