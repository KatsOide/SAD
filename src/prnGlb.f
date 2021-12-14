      Subroutine prnGlb(gname)
      use maccbk
      use maccode
      use macfile
      implicit real*8 (a-h,o-z)
      character*(*) gname
      integer*4  idx,ptr,n,itype,hsrch,ival
      real*8 rval
      idx=hsrch(gname)
      go to 1000
      Entry iprglb(idxp)
      idx=idxp
 1000 continue
      if(idtype(idx) .eq. icGLI) then
        ival=IgetGL(pname(idx),idx)
        if(outfl .ne. 6) print *,pname(idx),'=',ival,';'
        if(pname(idx)(1:1) .ne. '$') then
          write(outfl,*) pname(idx),'=',ival,';'
        else
          write(outfl,*) '!',pname(idx),'=',ival,';'
        endif
      else if(idtype(idx) .eq. icGLR) then
        rval=RgetGL(pname(idx),idx)
        if(outfl .ne. 6) print *,pname(idx),'=',rval,';'
        if(pname(idx)(1:1) .ne. '$') then
          write(outfl,*) pname(idx),'=',rval,';'
        else
          write(outfl,*) '!',pname(idx),'=',rval,';'
        endif
      else if(idtype(idx) .eq. icGLL) then
        call LgetGL(pname(idx),idx,ptr,n,itype)
        if(itype .eq. icGLR) then
          if(outfl .ne. 6)
     &          print *,pname(idx),'=(',(rlist(ptr+i),i=1,n),');'
          if(pname(idx)(1:1) .ne. '$') then
c            write(errfl,*) pname(idx),'size:',n,'ptr:',ptr
            write(outfl,*) pname(idx),'=(',(rlist(ptr+i),i=1,n),');'
          else
            write(outfl,*) '!',pname(idx),'=(',(rlist(ptr+i),i=1,n),');'
          endif
        else
          write(errfl,*) pname(idx),'has invalid element type id',itype
        endif
      else
        call errmsg('dassgn',
     &       ' Invalid global type'//pname(idx),0,16)
      endif
      return
      end
