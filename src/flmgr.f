      integer*8 function flmgr(ifn)
      use maccbk
      use macfile
      use tfmem, only:ktaloc,tfree
c
c     flmgr(0) returns last saved file number
c     flmgr(>0) push ifn onto stack.
      implicit none
      integer*4 ifn
c
      integer*8 itemp
c
      if(ifn .eq. 0) then
        flmgr=klist(inflpt)
        if(klist(inflpt+1) .ne. 0) then
          itemp=inflpt
          inflpt=klist(itemp+1)
          call tfree(itemp)
c          call freeme(itemp,1)
        endif
c        write(*,*)'@flmgr pop ',flmgr,inflpt
      else
        flmgr=ktaloc(5)
        klist(flmgr)=int8(ifn)
        klist(flmgr+1)=inflpt
        inflpt=flmgr
c        write(*,*)'@flmgr push ',ifn,inflpt
      endif
      return
      end
