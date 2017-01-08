      subroutine LgetGL(id,idxx,ptr,n,itype)
      use maccbk
      use maccode
      implicit none

      character*(*) id,idw*(MAXPNAME)
      integer*8 ptr
      integer*4 idxx,n,itype
      integer*4 hsrchz
c
      idxx=hsrchz(id)
      if(idtype(idxx) .eq. icGLL) then
        ptr=idval(idxx)
        n= ilist(1,ptr)-1
        itype=ilist(2,ptr)
      else
        idw=id
        call errmsg('getglobal','type miss match . '//idw//'!',0,0)
        ptr=0
        n= 0
        itype=0
      endif
      return
c
      entry LsetGL(id,idxx,ptr,n,itype)
c
      idxx=hsrchz(id)
      if(idtype(idxx) .eq. icGLL) then
        idval(idxx)=ptr
        ilist(1,ptr)=n+1
        ilist(2,ptr)=itype
      else
        idw=id
        call errmsg('setglobal','type miss match . '//idw//'!',0,0)
      endif
      return
      end
