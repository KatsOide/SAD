      integer*4 function IgetGL(id)
      use maccbk
      use maccode
      implicit none
      character*(*) id
      integer*4 idx
      character*(MAXPNAME) idw
      integer*4 hsrchz
      idx=hsrchz(id)
      if    ((idtype(idx) .eq. icGLI)
     &        .or. (idtype(idx) .eq. icFLAG))then
        IgetGL=int(idval(idx))
      else if(idtype(idx) .eq. icGLL) then
        IgetGL=int(idval(idx))
      else
        idw=id
        call errmsg('getglobal','type miss match . '//idw//'!',0,0)
        IgetGL=0
      endif
      return
      end
c
      integer*8 function IgetGL8(id)
      use maccbk
      use maccode
      implicit none
      character*(*) id
      integer*4 idx
      character*(MAXPNAME) idw
      integer*4 hsrchz
      idx=hsrchz(id)
      if    ((idtype(idx) .eq. icGLI)
     &        .or. (idtype(idx) .eq. icFLAG))then
        IgetGL8=idval(idx)
      else if(idtype(idx) .eq. icGLL) then
        IgetGL8=idval(idx)
      else
        idw=id
        call errmsg('getglobal','type miss match . '//idw//'!',0,0)
        IgetGL8=0
      endif
      return
      end
c
      real*8 function RgetGL(id,idx)
      use maccbk
      use maccode
      implicit none
      character*(*) id
      integer*4 idx
      character*(MAXPNAME) idw
      integer*4 hsrchz
      idx=hsrchz(id)
      if(idtype(idx) .eq. icGLR) then
        if(idval(idx) .gt. 0) then
          RgetGL=rlist(idval(idx))
        else
          idw=id
          call errmsg('RgetGL',
     &                'undefined global variable '//idw//'!',0,0)
          RgetGL=0.0d0
        endif
      else
        idw=id
        call errmsg('RgetGL','type mismatch '//idw//'!',0,0)
        RgetGL=0.0d0
      endif
      return
      end
