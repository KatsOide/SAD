      Subroutine doACT(cmdidx)
      use maccbk
      implicit none
      integer*4 cmdidx
c.....cmdidxbrings pointer to command name
      include 'inc/MACMISC.inc'
      include 'inc/MACCODE.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACVAR.inc'
      character*(MAXSTR) token
      integer*4 slen,ival,ttype,hsrch,hsrchz,mfalloc,mcfallo
      integer*4 scan, yyparse,yypushtoken,yypoptoken,yyreturn
      real*8 rval,yyval
      logical*4 skipch,ldummy,inlist
      integer*4 restme,membas,memuse,idx
      integer*4 i,lidx,iwork,ic
      integer*4 argp,argdfp,offset
c     
c     for debug
c     call ptrace('doAct '//pname(cmdidx)(2:)//'!',1)
c     end debug
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
         if (ttype .eq. ttypID) then
c     call errmsg('doAct',
c     &              'missing '''//LCURL//''' is assumed',0,0)
         else
           if((token(:slen) .eq. RCURL) .or.
     &          (token(:slen) .eq. SEMIC))then
           else
             call errmsg('doAct',
     &'syntax error:<command>{use <line name> <key word>=<value>}',
     &            0,16)
           endif
         endif
      endif
c
      argdfp=idval(cmdidx)
      if(ilist(2,argdfp) .le. 0) then
         argp=mcfallo(ilist(1,argdfp))
         ilist(2,argdfp)=argp
         ilist(1,argp)=ilist(1,argdfp)-1
         do i=1,ilist(1,argp)
            ilist(1,argp+i)=ilist(1,argdfp+i+1)
            ilist(2,argp+i)=0
         enddo
      else
         argp=ilist(2,argdfp)
      endif
c
 2000 if (skipch(COMMA,token,slen,ttype,rval,ival)) go to 2000
      inlist=.false.
      membas=0
      memuse=0
      restme=0
c
      if((token(:slen) .eq. RCURL) .or.
     &     (token(:slen) .eq. SEMIC))    go to  3000
c
c      write(*,*)"doACT:", token(:slen),slen,ttype
      if (ttype .eq. ttypID) then
         call capita(token(:slen))
         idx=hsrch(token(:slen))
         call gettok(token,slen,ttype,rval,ival)
         ldummy=skipch('=',token,slen,ttype,rval,ival)
         if (idtype(idx) .eq. icVAR) then
            offset=0
            do i=1,ilist(1,argdfp)
               if(ilist(1,argp+i) .eq. idx) then
                  offset=i
                  go to 2200
               endif
            enddo
c..........read right value
 2200       continue
            if (skipch(COMMA,token,slen,ttype,rval,ival)) go to 2200
            if (offset .eq. 0) then
               call errmsg(pname(cmdidx),'undefined parameter',0,0)
            else if( (.not. inlist) .and.
     &              skipch(LPAR,token,slen,ttype,rval,ival)) then
               inlist=.true.
               restme=mfalloc(-1)
               membas=mfalloc(restme)
               memuse=1
               ilist(2,argp+offset)=membas
            else if(inlist  .and.
     &              skipch(RPAR,token,slen,ttype,rval,ival)) then
               inlist=.false.
               if(restme .lt. memuse) then
                  call errmsg('doact'//pname(cmdidx)(2:),
     &                 ' break memory area.',0,0)
c     stop 9999
                  call pr_mem_map
               endif
               ilist(1,membas)=memuse
               call freeme(membas+memuse,restme-memuse)
            else if ((ttype .eq. ttypNM) .or.
     $              (ttype .eq. ttypID) .or.
     $              (ttype .eq. ttypDM) ) then
c............fordebug
c     print *,token(:slen),offset,pname(idx)
c............enddebug
               if (mod(idval(idx),VarLst) .eq. VarInt) then
                  scan=yypushtoken(token(:slen),slen,ttype)
                  scan=yyparse(yyreturn,yyval)
                  scan=yypoptoken(token,slen,ttype)
                  if (inlist) then
                     ilist(1,membas+memuse)=INT(yyval)
                     memuse=memuse+1
                  else
                     ilist(2,argp+offset)=INT(yyval)
                  endif
c                  print  *,' doACT/int:', INT(yyval)
               else if (mod(idval(idx),VarLst) .eq. VarRl) then
                  scan=yypushtoken(token(:slen),slen,ttype)
                  scan=yyparse(yyreturn,yyval)
                  scan=yypoptoken(token,slen,ttype)
                  if (inlist) then
                     rlist(membas+memuse)=yyval
                     memuse=memuse+1
                  else
                     ilist(2,argp+offset)=mcfallo(1)
                     rlist(ilist(2,argp+offset))=yyval
                  endif
c                  print  *,' doACT/real:', yyval
               else if(mod(idval(idx),varLst) .eq. VarPt) then
                  lidx=hsrchz(token(:slen))
                  if (inlist) then
                     ilist(2,membas+memuse)=lidx
                     memuse=memuse+1
                  else
                     ilist(2,argp+offset)=lidx
                  endif
                  call gettok(token,slen,ttype,rval,ival)
               else if(mod(idval(idx),VarLst).eq. VarLog) then
                  if (skipch('T',token,slen,ttype,rval,ival)) then
                     if (inlist) then
                        ilist(1,membas+memuse)=1
                        memuse=memuse+1
                     else
                        ilist(2,argp+offset)=1
                     endif
                  else if (skipch('F',token,slen,ttype,rval,ival)) then
                     if (inlist) then
                        ilist(1,membas+memuse)=0
                        memuse=memuse+1
                     else
                        ilist(2,argp+offset)=0
                     endif
                  else
                     call errmsg(pname(cmdidx),
     &                    pname(ilist(1,argp+offset))
     $                    //' is logical variable.'
     &                    //'range of variable is T or F',
     &                    0,0)
                     call gettok(token,slen,ttype,rval,ival)
                  endif
               else if (mod(idval(idx),VarLst) .eq. VarStr) then
                  if (inlist) then
                     iwork=membas+memuse
                     ilist(1,membas+memuse)=slen
                     ilist(2,membas+memuse)=membas+memuse+1
                     memuse=memuse+(slen+7)/8+1
                  else
                     iwork =mcfallo((slen+7)/8+1)
                     ilist(2,argp+offset)=iwork
                     ilist(1,iwork)=slen
                     ilist(2,iwork)=iwork+1
                  endif
                  do ic=1,(slen+7)/8
                     read(token(ic*8-7:ic*8),'(A8)')rlist(iwork+ic)
                  enddo
                  call gettok(token,slen,ttype,rval,ival)
               else
                  call errmsg(pname(cmdidx),
     &            'argument type does not match with the definition of'
     &                 //pname(ilist(1,argp+offset))//'.',
     &                0,0)
                  call errmsg(pname(cmdidx),
     &                 'argument given is '//token(:slen),
     &                 0,0)
               endif
            else
               call errmsg(pname(cmdidx),
     &              'unsupported type.'//token(:slen),0,0)
               call gettok(token,slen,ttype,rval,ival)
            endif
            if (inlist) goto 2200
         else
           call errmsg(pname(cmdidx),
     &         'syntax error:undefined variable '//pname(idx)
     &                                          ,0,0)
           call gettok(token,slen,ttype,rval,ival)
         endif
       else
         call errmsg(pname(cmdidx),'undefined Keyword',0,0)
         call errmsg(pname(cmdidx),token(:slen),0,0)
         call gettok(token,slen,ttype,rval,ival)
       endif
       go to 2000
c
 3000 continue
        call funcall1(ilist(1,argdfp+1),argp)
c for debug
c        call ptrace('doAct '//pname(cmdidx)(2:)//'!',-1)
c end debug
       return
       end
