      integer function flmgr(ifn)
      use maccbk
c
c     flmgr(0) returns last saved file number
c     flmgr(>0) push ifn onto stack.
      implicit real*8 (a-h,o-z)
      integer*4 ifn
      include 'inc/MACFILE.inc'
c
      integer*4 itemp
c
      if(ifn .eq. 0) then
        flmgr=ilist(1,inflpt)
        if(ilist(2,inflpt) .ne. 0) then
          itemp=inflpt
          inflpt=ilist(2,itemp)
          call freeme(itemp,1)
        endif
c        write(*,*)'@flmgr pop ',flmgr,inflpt
      else
        flmgr=mfalloc(1)
        ilist(1,flmgr)=ifn
        ilist(2,flmgr)=inflpt
        inflpt=flmgr
c        write(*,*)'@flmgr push ',ifn,inflpt
      endif
      return
      end
