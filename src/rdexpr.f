      real*8 function rdexpr(token,slen,ttype)
c     $Header: /SAD/cvsroot/oldsad/src/rdexpr.f,v 1.7.2.8 2016/09/19 04:55:35 oide Exp $
      use maccbk
      use mackw
      use macttyp
      use cbkmac
      use tfmem, only:ktaloc,tfree
      implicit none
      character*(*) token
      integer*4 slen,ttype
c     
      real*8 yyval
      integer*4 scan,yyparse,yypushtoken,yypoptoken,yyreturn
c     for debug
c     call ptrace('rdexpr '//token(:slen)//'!',1)
c      print *,'rdexpr: ',slen,' ',token(:slen)
c end debug
c
c   Oide changed 3 lines 8/20/1999
c      write(*,*)'rdexpr ',token(:slen),ttype,ttypDM
      rdexpr=0
      if (ttype .eq. ttypNM .or. ttype .eq. ttypID
     $     .or. token(:slen) .eq. '-'
     $     .or. token(:slen) .eq. '+'
     $     .or. token(:slen) .eq. LPAR) then
         scan=yypushtoken(token(:slen),slen,ttype)
         scan=yyparse(yyreturn,yyval)
c         print *,' rdexpr: ', scan, yyreturn,yyval
         rdexpr=yyval
c     yyparse reads one more token.
         scan=yypoptoken(token,slen,ttype)
c
c 4/22/2017 we abandon the error asignment scheme in the MAIN level.
c
c         print *,' rdexpr:',' next token is ',token(:slen),' ', ttype
c$$$      else if(token(:slen) .eq. LPAR ) then
c$$$c read list of numbers for error definiton
c$$$c       (<number1>,<number2>,......)
c$$$c or    (<distribution type>,<list of parameters>)
c$$$c       <distribution> = Normal, uniform etc.
c$$$c at first find end of error list
c$$$        idxerl=idval(elmidx)
c$$$        do while(idxerl .ne. 0)
c$$$          write(*,*)'rdexpr ',elmidx,idxerl
c$$$          idxerl=klist(idxerl)
c$$$        enddo
c$$$c        allmem=mtaloc(-1)
c$$$        membas=ktaloc(allmem)
c$$$        memuse=0
c$$$c.........fordebug
c$$$c     print *,'fallocated ',allmem,' from ',membas
c$$$c.........enddebug
c$$$c     
c$$$        klist(idxerl)=membas
c$$$        idxerl=klist(idxerl)
c$$$        ilist(2,idxerl-1)=kwcode
c$$$        klist(idxerl)=0
c$$$        memuse=memuse+2
c$$$        call gettok(token,slen,ttype,rval,ival)
c$$$        if (ttype .eq. ttypNM) then
c$$$           ilist(2,idxerl+1)=0
c$$$           ilist(1,idxerl+1)=Lrdnum(idxerl+1,token,slen,ttype)
c$$$           memuse=memuse+ilist(1,idxerl+1)
c$$$           rdexpr=rlist(idval(elmidx)+kytbl(kwcode,idtype(elmidx)))
c$$$        else if (ttype .eq. ttypID) then
c$$$           idxran=hsrch(token(:slen))
c$$$           if (idtype(idxran) .ne. icRAND) then
c$$$             token1=token(:slen)
c$$$              call errmsg('rdexor',
c$$$     &             token1(:slen)//' is not supported as '//
c$$$     &             'distribution',0,0)
c$$$              call gettok(token,slen,ttype,rval,ival)
c$$$ 3100         continue
c$$$              if ( .not. skipch(RPAR,token,slen,ttype,rval,ival))
c$$$     &             go to 3100
c$$$           else
c$$$c     ilist(2,idxerl+1)=idval(idxran) changed by NY apr.9,88
c$$$              ilist(2,idxerl+1)=idxran
c$$$              call gettok(token,slen,ttype,rval,ival)
c$$$              ilist(1,idxerl+1)=Lrdnum(idxerl+1,token,slen,ttype)
c$$$              memuse=memuse+ilist(1,idxerl+1)
c$$$           endif
c$$$           rdexpr=rlist(idval(elmidx)+kytbl(kwcode,idtype(elmidx)))
c$$$        endif
c$$$        if(allmem .gt. memuse+3)then
c$$$          ilist(1,membas+memuse)=allmem-memuse
c$$$          call tfree(membas+memuse+1)
c$$$        endif
c$$$c        call freeme(membas+memuse,allmem-memuse)
      else
         call errmsg('rdexpr',
     &        'syntax error:<number>  <unit> ! or (<num1>,...)',0,16)
      endif
c.....fordebug
c     call ptrace('rdexpr '//token(:slen)//'!',-1)
c.....end debug
      return
      end
