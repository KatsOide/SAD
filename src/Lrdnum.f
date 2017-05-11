      function Lrdnum(idxl,token,slen,ttype)
c     $Header: /SAD/work/cvsclone/oldsad/src/Lrdnum.f,v 1.2.14.5 2016-09-09 12:08:57 oide Exp $
      use maccbk
      use maccode
      use macttyp
      use macmisc
      implicit none
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
