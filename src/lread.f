      subroutine lread(idxl,token,slen,ttype,rval,ival)
      use maccbk
      use tfstk
      use cbkmac
      implicit none
      integer*4 idxl,idx
      character*(MAXSTR) token,wkstr
      integer*4 slen,ival,ttype
      real*8 rval
c     
      integer*8 ia
      integer*4 hsrchz
      integer*4 n,isp0,mult
      logical skipch
c     
c     first character must be '='
      if (.not. skipch(EQCHR,token,slen,ttype,rval,ival)) then
         call errmsg('lread',
     &        'missing '''//EQCHR//''' is assumed.',0,0)
      endif
c     next is to be a '('
      if (.not. skipch(LPAR,token,slen,ttype,rval,ival)) then
         call errmsg('lread',
     &        'missing '''//LPAR//''' is assumed.',0,0)
      endif
c     intialize parameters
      isp0=isp
c     now read elemnt list
      w2000: do while(.true.)
        if (skipch(COMMA,token,slen,ttype,rval,ival))then
          cycle
        endif
        mult=1
        w2100: do while(.true.)
          if(skipch(PLUS,token,slen,ttype,rval,ival))then
            cycle
          else if(skipch(MINUS,token,slen,ttype,rval,ival)) then
            mult=-mult
            cycle
          else if (ttype .eq. ttypID) then
            idx=hsrchz(token(:slen))
            if ( (idtype(idx) .gt. icMXEL) .and.
     &           (idtype(idx) .ne. icLINE)      ) then
              wkstr(:slen)=token(:slen)
              call errmsg('lread',
     &             'syntax error: '//wkstr(:slen)
     &             //' is not line nor element.',0,16)
            else
              isp=isp+1
              itastk(1,isp)=mult
              itastk(2,isp)=idx
            endif
            exit
          else if(ttype .eq. ttypNM) then
            if (ival .eq. 0) then
              call errmsg('lread',
     &             'invalid muliplicity',0,16)
            else
              mult=mult*ival
            endif
            call gettok(token,slen,ttype,rval,ival)
            if (.not. skipch(STAR,token,slen,ttype,rval,ival)) then
              if (ttype .eq. ttypID) then
                call errmsg('lread',
     &               'missing '''//STAR//''' is assumed',0,0)
              else
                call errmsg('lread',
     &               'syntax error:<line elemnt>:: <num>*!<ID>',0,16)
              endif
            endif
            cycle
          else if (ttype .eq. ttypDM) then
            if(.not. skipch(RPAR,token,slen,ttype,rval,ival))then
              call errmsg('lread',
     &             'ivalid delimiter',0,0)
            endif
            exit w2000
          else
            call errmsg('lread',
     &           'syntax error: <num>*!<element name> is expected',0,16)
          endif
        enddo w2100
        call gettok(token,slen,ttype,rval,ival)
      enddo w2000
      n=isp-isp0
      ia=ktaloc(n+1)
      ilist(1,ia)=n
      ilist(2,ia)=0
      klist(ia+1:ia+n)=ktastk(isp0+1:isp0+n)
      idval(idxl)=ia
c      write(*,*)'lread ',ia,n,pname(idxl),idtype(idxl)
      return
      end
