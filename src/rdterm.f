      Subroutine rdterm(token,slen,ttype,val,status)
c     $Header: /SAD/cvsroot/oldsad/src/rdterm.f,v 1.2.14.2 2012/09/15 08:14:57 oide Exp $
      use maccbk
      implicit none
      character*(*) token
      integer slen,ttype,status
      integer*4 scan, yyparse,yypushtoken,yypoptoken,yyreturn
      real*8 val,yyval

      include 'inc/MACCODE.inc'
      include 'inc/MACTTYP.inc'
c
c      print *,'rdterm:',slen,token(:slen)
      status=0
      scan=yypushtoken(token,slen,ttype)
      scan=yyparse(yyreturn,yyval)
      val=yyval
c yyparse reads one more token.
      scan=yypoptoken(token,slen,ttype)
      return
      end
