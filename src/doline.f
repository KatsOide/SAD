      subroutine doline(incode)
      use trackbypass
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      use tfmem, only:tfree
      implicit none
      character*(MAXSTR) token
      integer*8 incode,icode,idx1,idm
      integer*4 index,idx,membas
      integer slen,ival,ttype,hsrchz,slen2,ttype2
      real*8 rval
      logical skipch
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .eq. ttypID) then
c         call errmsg('doline',
c    &              'missing '''//LCURL//''' is assumed',0,0)
        else
          call errmsg('doline',
     &                ' syntax error:'//LCURL//' is expected.',0,16)
        endif
      endif
      icode=incode
      go to 2000
c
c strat line definition
      entry dolin2(index,slen2,ttype2)
      slen=slen2
      ttype=ttype2
      token=pname(index)(:slen)
      icode=idtype(index)
c
 2000 if (skipch(COMMA,token,slen,ttype,rval,ival)) go to 2000
c.......for debug
c       print *,'read line defininition: ',slen,'"',token(:slen),'"'
c.......end debug
c
      if((token(:slen) .eq. RCURL) .or.
     &   (token(:slen) .eq. SEMIC)) then
        return
c
      else if (ttype .eq. ttypID) then
        idx=hsrchz(token(:slen))
        call gettok(token,slen,ttype,rval,ival)
        if ((idtype(idx) .eq. icNULL)
     &       .or. (idtype(idx) .eq. icode))then
          if(idtype(idx) .eq. icode)then
            call errmsg('doline',pname(idx)(:lpname(idx))//
     &           ' is redefined',0,-1)
          endif
          idtype(idx)=int(icode)
          idx1=idval(idx)
          if(idx1 .ne. 0) then
            membas=ilist(2,idx1)
            if(membas .gt. 0) then
              idm=idval(membas)
              if(idtype(membas) .eq. icLine)then
                if(idm .ne. lattuse)then
                  call tclrline(idm)
                  idtype(membas)=0
                else
                  lattredef=lattuse
                endif
              elseif(idtype(membas) .gt. 0)then
c                write(*,*)'doline ',idm,idtype(membas)
                ilist(1,idm)=ilist(1,idm)+1
                call tfree(idm+1)
                idtype(membas)=0
              endif
              ilist(2,idx1)=0
            endif
            call tfree(idx1)
            idval(idx)=0
          endif
c          print *,'doline-before-lread'
          call lread(idx,token,slen,ttype,rval,ival)
c          print *,'doline-lread: ',slen,ttype,idx,idtype(idx),
c     $         icode,'"',token(:slen),'"'
        else if (idtype(idx) .ge. icRSVD) then
          call errmsg('doline',
     &        'syntax error:you can not use reaserved words as a name'
     &        ,0,16)
        else
          call errmsg('doline',
     &         'syntax error: '//pname(idx)//' is already used.'
     &         ,0,16)
        endif
      else
        call errmsg('doline',
     &       'syntax error:LINE{<line name>=(line definition)}',0,16)
      endif
      go to 2000
c
      end
