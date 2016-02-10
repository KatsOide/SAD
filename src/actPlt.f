       Subroutine actPlt(argp)
       use maccbk
       implicit real*8 (a-h,o-z)
       integer argp
       integer plot$,ptype$,turns$,npart$,span$,cent$,nparm
       parameter (ptype$=1,turns$=3,npart$=2,span$=4,
     &            cent$=5,nparm=6)
       integer last$,newpl$
c......for debug
c      call ptrace('ActPlt',1)
c......end debug
c......fallocate memory for new ploting definition
       last$=0
       newpl$=IgetGL('$PLOT$',plot$)
c......find last of list.
 1000 continue
        if (newpl$ .eq. 0) goto 1100
        last$=newpl$
        newpl$=ilist(2,last$)
        go to 1000
 1100 continue
       newpl$=mcfallo(nparm)
       ilist(1,newpl$)=nparm
       ilist(2,newpl$)=0
       if (last$ .eq. 0) then
         call IsetGL('$PLOT$',newpl$,plot$)
       else
         ilist(2,last$)=newpl$
       endif
       last$=newpl$
c......get values of parameters
       ilist(1,last$+ptype$)=ilist(2,argp+ptype$)
       ilist(1,last$+npart$)=ilist(2,argp+npart$)
       ilist(1,last$+turns$)=ilist(2,argp+turns$)
       rlist(last$+span$)=rlist(ilist(2,argp+span$))
       rlist(last$+cent$)=rlist(ilist(2,argp+cent$))
c......for debug
c      print *,
c    & ilist(1,last$+ptype$),
c    & ilist(1,last$+turns$),
c    & ilist(1,last$+npart$),
c    & rlist(last$+span$),
c    & rlist(last$+cent$)
c......end debug
c......for debug
c      call ptrace('ActPlt',-1)
c......end debug
       return
       end
