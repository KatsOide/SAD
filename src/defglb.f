      Subroutine defglb(id,type,idx)
      use maccbk
      use maccode
      use tfmem, only:ktaloc
      implicit none
      character*(*) id,idw*(MAXSTR)
      integer*4 type,idx,lenw
      real*8 rval
      integer*8 ival8
      integer*4 ival,idxx,hsrchz,idx1
      logical isglobal
      isglobal(idx)=((idtype(idx) .eq. icGLR)
     &         .or.  (idtype(idx) .eq. icGLI)
     &         .or.  (idtype(idx) .eq. icFLAG)
     &         .or.  (idtype(idx) .eq. icGLL))
c
      idx=hsrchz(id)
      if (idtype(idx) .eq. icNULL) then
        idtype(idx)=type
      else if(isglobal(idx)) then
         if(type .ne. icNULL .and. type .ne. idtype(idx)) then
            idw=id
            call errmsg('defglb',
     $           idw(:lenw(idw))//' cannot be used as '//
     &           'global variable name',0,4)
         endif
      else
         idw=id
         call errmsg('defglb',
     $       idw(:lenw(idw))//' cannot be used as '//
     &        'global variable name',0,4)
      endif
      return
c
      entry RsetGL(id,rval,idx1)
      idxx=hsrchz(id)
      if(idtype(idxx) .eq. icGLR .or. idtype(idxx) .eq. icNULL) then
        if (idval(idxx) .le. 0)then
          idval(idxx)=ktaloc(3)
c          idval(idxx)=mfalloc(1)
        endif
        rlist(idval(idxx))=rval
      else
        idw=id
        call errmsg('setglobal','type mismatch . '//idw//'!',0,0)
      endif
      return
c
      entry IsetGL(id,ival,idx1)
c
      idxx=hsrchz(id)
      if((idtype(idxx) .eq. icGLI)
     &   .or. (idtype(idxx) .eq. icFLAG))then
        idval(idxx)=int8(ival)
      else
        idw=id
        call errmsg('setglobal','type mismatch . '//idw//'!',0,0)
      endif
      return

      entry IsetGL8(id,ival8,idx1)
c
      idxx=hsrchz(id)
      if((idtype(idxx) .eq. icGLI)
     &   .or. (idtype(idxx) .eq. icFLAG))then
        idval(idxx)=ival8
      else
        idw=id
        call errmsg('setglobal','type mismatch . '//idw//'!',0,0)
      endif
      return
c
      end
