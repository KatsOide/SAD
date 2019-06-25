      Subroutine rdterm(token,slen,ttype,val,status)
c     $Header: /SAD/cvsroot/oldsad/src/rdterm.f,v 1.2.14.4 2016/09/09 04:38:31 oide Exp $
      use maccbk
      use mackw
      use macttyp
      implicit none
      character*(*) token
      integer slen,ttype,status
      integer*4 scan, yyparse,yypushtoken,yypoptoken,yyreturn
      real*8 val,yyval
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
