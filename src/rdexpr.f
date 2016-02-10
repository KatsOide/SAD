      real*8 function rdexpr(elmidx,kwcode,token,slen,ttype)
c     $Header: /SAD/cvsroot/oldsad/src/rdexpr.f,v 1.7.2.2 2012/09/15 08:14:57 oide Exp $
      use maccbk
      implicit none
      integer*4 elmidx,kwcode
      character*(*) token
      character*64 token1
      integer*4 slen,ival,ttype,mfalloc
      real*8 rval
      include 'inc/CBKMAC.inc'
c     
      real*8 yyval
      integer*4 Lrdnum,scan,
     $     yyparse,yypushtoken,yypoptoken,yyreturn
     $     ,allmem,membas,memuse,hsrch,idxerl,idxran
      logical skipch
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
     $     .or. token(:slen) .eq. '+') then
         scan=yypushtoken(token(:slen),slen,ttype)
         scan=yyparse(yyreturn,yyval)
c         print *,' rdexpr: ', scan, yyreturn,yyval
         rdexpr=yyval
c     yyparse reads one more token.
         scan=yypoptoken(token,slen,ttype)
c         print *,' rdexpr:',' next token is ',token(:slen),' ', ttype
      else if(token(:slen) .eq. LPAR ) then
c read list of numbers for error definiton
c       (<number1>,<number2>,......)
c or    (<distribution type>,<list of parameters>)
c       <distribution> = Normal, uniform etc.
c at first find end of error list
        idxerl=idval(elmidx)
 2100   continue
        if (ilist(2,idxerl).ne. 0) then
           idxerl=ilist(2,idxerl)
           go to 2100
        endif
c
        allmem=mfalloc(-1)
        membas=mfalloc(allmem)
        memuse=0
c.........fordebug
c     print *,'fallocated ',allmem,' from ',membas
c.........enddebug
c     
        ilist(2,idxerl)=membas
        idxerl=ilist(2,idxerl)
        ilist(1,idxerl)=kwcode
        ilist(2,idxerl)=0
        memuse=memuse+2
        call gettok(token,slen,ttype,rval,ival)
        if (ttype .eq. ttypNM) then
           ilist(2,idxerl+1)=0
           ilist(1,idxerl+1)=Lrdnum(idxerl+1,token,slen,ttype)
           memuse=memuse+ilist(1,idxerl+1)
           rdexpr=rlist(idval(elmidx)+kytbl(kwcode,idtype(elmidx)))
        else if (ttype .eq. ttypID) then
           idxran=hsrch(token(:slen))
           if (idtype(idxran) .ne. icRAND) then
             token1=token(:slen)
              call errmsg('rdexor',
     &             token1(:slen)//' is not supported as '//
     &             'distribution',0,0)
              call gettok(token,slen,ttype,rval,ival)
 3100         continue
              if ( .not. skipch(RPAR,token,slen,ttype,rval,ival))
     &             go to 3100
           else
c     ilist(2,idxerl+1)=idval(idxran) changed by NY apr.9,88
              ilist(2,idxerl+1)=idxran
              call gettok(token,slen,ttype,rval,ival)
              ilist(1,idxerl+1)=Lrdnum(idxerl+1,token,slen,ttype)
              memuse=memuse+ilist(1,idxerl+1)
           endif
           rdexpr=rlist(idval(elmidx)+kytbl(kwcode,idtype(elmidx)))
        endif
        call freeme(membas+memuse,allmem-memuse)
      else
         call errmsg('rdexpr',
     &        'syntax error:<number>  <unit> ! or (<num1>,...)',0,16)
      endif
c.....fordebug
c     call ptrace('rdexpr '//token(:slen)//'!',-1)
c.....end debug
      return
      end
