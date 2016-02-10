      Subroutine defglb(id,type,idx)
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      character*(*) id,idw*(MAXSTR)
      integer*4 type,idx,idx1,lenw
      real*8 rval
      integer*4 ival,idxx,hsrchz,mfalloc
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
        if (idval(idxx) .le. 0) idval(idxx)=mfalloc(1)
        rlist(idval(idxx))=rval
      else
        idw=id
        call errmsg('setglobal','type miss match . '//idw//'!',0,0)
      endif
      return
c
      entry IsetGL(id,ival,idx1)
c
      idxx=hsrchz(id)
      if((idtype(idxx) .eq. icGLI)
     &   .or. (idtype(idxx) .eq. icFLAG))then
        idval(idxx)=ival
      else
        idw=id
        call errmsg('setglobal','type miss match . '//idw//'!',0,0)
      endif
      return
c
      end
