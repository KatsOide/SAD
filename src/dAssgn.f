      Subroutine dAssgn(token,slen,status)
      use maccbk
      use maccode
      use macttyp
      use macvar
      use macmisc
      implicit none
c     
      character*(MAXSTR) token
      integer slen,status

      real*8 rval, val
      logical skipch,skiped
      integer*4 idx,ival,ttype
      integer*8 newblk,allmem,membas,memptr,memuse,m
      integer*8 ktaloc
c     macro functions
      logical issign
      character char
      issign(char)=(char .eq. '+') .or. (char .eq. '-')
C     implementation
      status=0
      call defglb(token(:slen),icNULL,idx)
      if ((idtype(idx) .ne. icGLR) .and.
     &     (idtype(idx) .ne. icGLI) .and.
     &     (idtype(idx) .ne. icGLL) .and.
     &     (idtype(idx) .ne. icNULL))then
         status=-1
         return
      endif
c     
      call gettok(token,slen,ttype,rval,ival)
      skiped = skipch('=',token,slen,ttype,rval,ival)
      if(ttype .eq. ttypNM .or. ttype .eq. ttypID
     $     .or. issign(token(:slen)) ) then
         call rdterm(token,slen,ttype,rval,status)
c         write(*,*)'dAssgn ',ttype,rval,ival,status,token(1:slen)
         if (status .ne. 0) return
         if (idtype(idx) .eq. icGLI) then
            idval(idx)=INT8(rval)
         else if ((idtype(idx) .eq. icGLR) .or.
     &           (idtype(idx) .eq. icNULL)) then
            if(idval(idx) .le. 0)then
              idval(idx)=ktaloc(3)
c              idval(idx)=mfalloc(1)
            endif
            idtype(idx)=icGLR
            rlist(idval(idx))=rval
         else if (idtype(idx) .eq. icGLL) then
            call tfreem(idval(idx),ilist(1,idval(idx)))
            idval(idx)=ktaloc(3)
c            call freeme(idval(idx),ilist(1,idval(idx)))
c            idval(idx)=mfalloc(2)
            ilist(1,idval(idx))=1
            if(ilist(2,idval(idx)) .eq. icNULL)
     &           ilist(2,idval(idx))=icGLR
            rlist(idval(idx)+1)=rval
         endif
      else if(skipch(LPAR,token,slen,ttype,rval,ival)) then
         if(idtype(idx) .ne. icGLL .and. idtype(idx) .ne.icNULL) then
            status=-1
            call errmsg('dasgn', 'type mismatch',0,4)
            return
         endif
         if(idtype(idx) .eq. icGLL)then
           call tfreem(idval(idx),ilist(1,idval(idx)))
c           call freeme(idval(idx),ilist(1,idval(idx)))
         endif
         call defglb(pname(idx),icGLL,idx)
         allmem=pagesz/4
         membas=ktaloc(allmem)
c         membas=mfalloc(allmem)
         if(membas .eq. 0) then
            call errmsg('dAssgn',' cannot allocate memory',0,0)
            stop
         end if
c     Note
c     ilist(1,membas): size
c     ilist(2,membas): data type
c     See ilist(*,ptr)@src/LgetGL.f
         memptr=membas+1
 2000    if(skipch(COMMA,token,slen,ttype,rval,ival)) go to 2000
         if(skipch(RPAR,token,slen,ttype,rval,ival)) go to 3000
         if(token(:slen) .eq. SEMIC) go to 3000
         call rdterm(token,slen,ttype,val,status)
         if(status .ne. 0) then
            call errmsg('dAssgn',
     &           'sytax error.',0,4)
            return
         endif
         rlist(memptr)=val
         memptr=memptr+1
         if(memptr-membas .ge. allmem) then
               newblk=ktaloc(2*allmem)
c               newblk=mfalloc(2*allmem)
               if (newblk .eq. 0) then
                  call errmsg('dAssgn',
     &                 ' cannot extend working area.',32,0)
                  stop
               end if
               do m=1,memptr-membas-1
                  klist(newblk+m)=klist(membas+m)
               end do
               memptr=newblk+(memptr-membas)
               call tfree(membas)
c               call freeme(membas,allmem)
               allmem=2*allmem
               membas=newblk
         endif
         go to 2000
c     End of List.
 3000    continue
         if (token(:slen) .ne. SEMIC) then
            call errmsg('dAssgn',' Illegal delimiter ',0,4)
         endif
         memuse=memptr-membas
         if(allmem .lt. memuse) then
            call errmsg('dassgn',
     &           ' broken memory area.',0,16)
            stop 9999
         endif
         if(allmem .gt. memuse+3)then
           ilist(1,memptr)=int(allmem-memuse)
           call tfree(memptr+1)
         endif
c         call freeme(memptr,allmem-memuse)
         call LsetGL(pname(idx),idx,membas,memuse-1,icGLR)
      endif
c     for debug
c     call ptrace('dAssgn',-1)
c     end debug
      return
      end
