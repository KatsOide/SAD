      subroutine dolist(dummy)
      use maccbk
      implicit real*8 (a-h,o-z)
      include 'inc/MACFILE.inc'
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
      character*(MAXSTR) token
      integer slen,ival,ttype,hsrch,lpname
      real*8 rval
      logical skipch
      logical Ldump,Lline,Lelem(icMXEL),Lcmd,Lmem
      integer idx,lcount,used,iky
c
          Lmem=.false.
          Lcmd=.false.
          Ldump=.false.
          Lline=.false.
          do 100 i=1,icMXEL
            Lelem(i)=.false.
  100     continue
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .eq. ttypID) then
c         call errmsg('dolist',
c    &              'missing '''//LCURL//''' is assumed',0,0)
        else
          call errmsg('dolist',
     &      ' syntax error:list{ ALL!  LINE!  DUMP! ......}',
     &                0,16)
        endif
      endif
c
 2000 if (skipch(COMMA,token,slen,ttype,rval,ival)) go to 2000
      call CAPITA(token(:slen))
      if (ttype .eq. ttypID) then
        if(token(:slen) .eq. 'DUMP') then
          Ldump=.true.
        else if(token(:slen) .eq. 'ALL') then
          Lline=.true.
c         Ldrift=.true.
          do 2100 i=1,icMXEL
            Lelem(i)=.true.
 2100     continue
        else if(token(:slen) .eq. 'LINE') then
          Lline=.true.
        else if((idtype(hsrch(token(:slen)))  .eq. icDEF) .and.
     &          idval(hsrch(token(:slen))) .lt. icMXEL) then
          Lelem(idval(hsrch(token(:slen))))=.true.
        else if(token(:slen) .eq. 'COMMAND') then
          Lcmd=.true.
        else if(token(:slen) .eq. 'MEMORY' .or.
     &          token(:slen) .eq. 'MEM'    .or.
     &          token(:slen) .eq. 'mem'    .or.
     &          token(:slen) .eq. 'memmory'
     &                                    ) then
          Lmem=.true.
        endif
      else if((token(:slen) .eq. RCURL) .or.
     &        (token(:slen) .eq. SEMIC)) then
        go to 3000
      else
        call errmsg('dolist',
     &       'syntax error:illegal token as input',0,0)
        call errmsg('dolist',
     &       'syntax: print{<name of line or element> , ...!}',0,16)
      endif
      call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
 3000 continue
      lcount=0
      used=0
      do 3900 idx = 1,HTMAX
        if (mod(lcount , 55) .eq. 0) then
            write(errfl,*)'index  ','Name  ','Type  ','Value  ','other'
            lcount=lcount+1
        endif
        if(idtype(idx) .ne. icNULL) then
        used=used+1
        if ((idtype(idx) .eq. icLINE) .and. (Lline .or. ldump)) then
          write(errfl,*)
     $          idx, pname(idx)(:lpname(idx)) ,'LINE    ', idval(idx)
     &          , rlist(idval(idx))
          lcount=lcount+1
        elseif ((idtype(idx) .lt. icMXEL) .and.
     &         (Lelem(idtype(idx)) .or. ldump)) then
           iky=kytbl(0,idtype(idx))
           write(errfl,*   )
     $          idx, pname(idx)(:lpname(idx)), pname(iky)(:lpname(iky))
           lcount=lcount+1
        elseif ((idtype(idx) .eq. icRSVD) .and. (Lcmd
     &          .or. ldump)) then
           write(errfl,*) idx,pname(idx)(:lpname(idx))
     $          ,'COMMAND ',idval(idx)
           lcount=lcount+1
        elseif ((idtype(idx) .eq. icDEF) .and. (Lcmd
     &          .or. ldump)) then
           write(errfl,*) idx,pname(idx)(:lpname(idx))
     $          ,'ELEMENT ',idval(idx)
           lcount=lcount+1
        elseif ((idtype(idx) .eq. icACT) .and. (Lcmd
     &          .or. ldump)) then
           write(errfl,*   ) idx,pname(idx)(:lpname(idx))
     $          ,'Action  ',idval(idx)
     &          ,ilist(1,idval(idx))
     &          ,ilist(2,idval(idx))
           lcount=lcount+1
        elseif ((idtype(idx) .eq. icGLI) .and. ldump) then
           write(errfl,*   ) idx,pname(idx)(:lpname(idx))
     $          ,'Integer ',idval(idx)
           lcount=lcount+1
        elseif ((idtype(idx) .eq. icGLR) .and. ldump) then
           write(errfl,*   ) idx,pname(idx)(:lpname(idx))
     $          ,'Real    ',idval(idx), rlist(idval(idx))
           lcount=lcount+1
        else if (ldump) then
           if(idval(idx) .ne. 0) then
              write(errfl,*) idx,pname(idx)(:lpname(idx))
     $             ,idtype(idx) ,idval(idx)
     &             ,rlist(idval(idx))
           endif
           lcount=lcount+1
        endif
      endif
 3900 continue
      write(errfl,
     $     '(1H ,''used Hash table='',I5,''/'',I5)')
     $     used, HTMAX
      if (Lmem) then
         call pr_mem_map
      endif
      return
      end
