      subroutine doflag(fval,dummy)
      use maccbk
      use mackw
      use macttyp
      use macmisc
      implicit none
      integer fval
      real*8 dummy
c
      character*(MAXSTR) token
      integer slen,ival,ttype,hsrchz,idx
      real*8 rval
      logical skipch,skipped
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .eq. ttypID) then
c         call errmsg('doflag',
c    &              'missing '''//LCURL//''' is assumed',0,0)
        else
          call errmsg('doflag',
     &                ' syntax error:'//LCURL//' is expected.',0,0)
          return
        endif
      endif
C
 2000 continue
c     write(*,*)'@continue'
        skipped=skipch(COMMA,token,slen,ttype,rval,ival)
        if(token(:slen) .eq. RCURL) return
        if(token(:slen) .eq. SEMIC) return
        if (ttype .eq. ttypID) then
          idx =hsrchz('$'//token(:slen)//'$')
c     write(*,*)'@hsrchz',idx,slen
          if(idtype(idx) .eq. icFLAG) then
c...........   for debug
c           call errmsg('onflag','$'//token(:slen)//'$',0,0)
c...........   end debug
            call IsetGL('$'//token(:slen)//'$',fval,idx)
c     write(*,*)'@isetgl'
          else
            call defglb('$'//token(:slen)//'$',icFLAG,idx)
c     write(*,*)'@defglb'
            call IsetGL('$'//token(:slen)//'$',fval,idx)
c     write(*,*)'@isetgl'
            call errmsg('doflag','New flag  '//
     &        token(:slen)//' is defined.',0,-1 )
          endif
        else
          call errmsg('doflag',
     &      'syntax error:on!off{<FLAG name LIST>}',0,16)
        endif
c
        call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
      end
c
      subroutine doonfl(dummy)
      use maccbk
      implicit none
      real*8 dummy
c
       call doflag(FLAGON,dummy)
       return
       end
c
      subroutine dooffl(dummy)
      use maccbk
      implicit none
      real*8 dummy
c     
      call doflag(FLAGOF,dummy)
      return
      end
