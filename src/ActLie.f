      Subroutine ActLie(argp)
      use maccbk
      use macvar
      implicit real*8 (a-h,o-z)
      integer*4 argp
C
C      idummy=sethtb('LIE     ',icACT,mfalloc(3))
C      ilist(1,idval(idummy))=2
C      call setfnp(ilist(1,idval(idummy)+1),Act)
C      ilist(1,idval(idummy)+2)=hsrch('USE')
      integer len,idx,hsrchz
      integer*8 rslvin
c
      save
c
      len=ilist(1,argp)
      idx=ilist(2,argp+1)
      nidx=hsrchz('+'//pname(idx))
      idtype(nidx)=idtype(idx)
      if(idval(nidx) .eq. 0) then
        idval(nidx)=rslvin(idx)
      endif
      if (ilist(2,idval(nidx)) .le. 0) then
        call errmsg('ActLie',
     &       'Expanding '//pname(nidx)//' now.',0,0)
        call expnln(nidx)
      endif
c      call pr_mem_map
c     call prexln(nidx,'            ')
      call AAALIE(nidx)
c      call pr_mem_map
      return
      end
