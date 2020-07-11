       Subroutine actPlt(argp)
       use maccbk
       implicit real*8 (a-h,o-z)
       integer*8 argp,ktcaloc,newpl$,last$,IgetGL8
       integer plot$,ptype$,turns$,npart$,span$,cent$,nparm
       parameter (ptype$=1,turns$=3,npart$=2,span$=4,
     &            cent$=5,nparm=6)
c......for debug
c      call ptrace('ActPlt',1)
c......end debug
c......fallocate memory for new ploting definition
       last$=0
       newpl$=IgetGL8('$PLOT$')
c......find last of list.
 1000 continue
        if (newpl$ .eq. 0) goto 1100
        last$=newpl$
        newpl$=klist(last$)
        go to 1000
 1100 continue
       newpl$=ktcaloc(nparm)
       klist(newpl$)=0
       if (last$ .eq. 0) then
         call IsetGL8('$PLOT$',newpl$,plot$)
       else
         klist(last$)=newpl$
       endif
       last$=newpl$
c......get values of parameters
       ilist(1,last$+ptype$)=ilist(1,argp+ptype$*2)
       ilist(1,last$+npart$)=ilist(1,argp+npart$*2)
       ilist(1,last$+turns$)=ilist(1,argp+turns$*2)
       rlist(last$+span$)=rlist(ilist(1,argp+span$*2))
       rlist(last$+cent$)=rlist(ilist(1,argp+cent$*2))
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
