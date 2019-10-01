      subroutine rdkwdl(elmidx)
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      implicit none
      integer elmidx
c
      character*(MAXSTR) token
      integer*4 slen,ival,ttype,hsrch,offset,idx
      real*8 rval,rdexpr
      logical skipch
c
c for debug
c     print *,'>enter into rdkwdl'
c end debug
c first character must be '='
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(EQCHR,token,slen,ttype,rval,ival)) then
           call errmsg('rdkwdl',
     &          'missing '''//EQCHR//''' is assumed.',0,0)
      endif
c next is to be a '{' or '('
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival)) then
        if (.not. skipch(LPAR,token,slen,ttype,rval,ival)) then
           call errmsg('rdkwdl',
     &          'missing '''//LCURL//''' is assumed.',0,0)
        endif
      endif
c
      ilist(1,idval(elmidx))=kytbl(kwMAX,idtype(elmidx))
c
c  now read list of <key word> = <expression>
 2000 continue
        if (ttype .eq. ttypID) then
c for debug
c        print *,'keyword: ',token(:slen)
c end debug
c         token(:8)=token(:slen)//NULSTR !deleted jul,16,'87 by N.Y
          call CAPITA(token(:slen))
          idx= hsrch(token(:slen))
          if (idtype(idx) .eq. icKWRD) then
            call gettok(token,slen,ttype,rval,ival)
            if (.not. skipch(EQCHR,token,slen,ttype,rval,ival)) then
              call errmsg('rdkwdl',
     &             'missing '''//EQCHR//''' is assumed.',0,0)
            endif
            offset=kytbl(idval(idx),idtype(elmidx))
c           print *,idtype(elmidx),kwMAX
            if (offset .eq. 0) then
              call errmsg('rdkwdl',
     &             'invalid key word name: '//
     &              pname(kytbl(idval(idx),0)),0,16)
            else
c for debug
c             print *,'>>goto rdexpr '//token(:slen)//'!'
c end debug
              rlist(idval(elmidx)+offset)=rdexpr(token,slen,ttype)
c rdexpr  return next token in token(:slen)
            endif
          else
            write(*,*)'rdkwl ',idx,idtype(idx)
            call errmsg('rdkwdl',
     &           'syntax error: Undefined keyword: '//token(:slen)
     $           ,0,16)
          endif
        else if ((token(:slen) .eq. RCURL)
     &       .or.(token(:slen) .eq. SEMIC)
     &       .or.(token(:slen) .eq. RPAR)) then
c for debug
c         print *,'<leave from rdkwdl'
c end debug
          return
        else
          call errmsg('rdkwdl',
     &         'syntax error:<key word>=<expression>.',0,16)
        endif
        if (token(:slen) .eq. COMMA)then
          call gettok(token,slen,ttype,rval,ival)
        endif
      go to 2000
c
      end
