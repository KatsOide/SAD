      subroutine filaux(idx)
      use maccbk
      use mackw
      use macphys
      use macfile
      use macmisc
      implicit none
      integer*4 idx

      integer NAUX
      parameter (NAUX=6)
      integer*8 pexln,iaux,ktcaloc
      integer*4 llen,status,iw,idummy
c      integer*4 i,IgetGL
c      real*8 dist
      real*8 u0,Vc,codin(6),codou(6)
      real*8 g,p0,pmass,phis,RgetGL
C.....calculate total length of line
c      write(*,*)'filaux-0 ',idx
      pexln=idval(ilist(2,idval(idx)))
      if (ilist(2,pexln) .eq. 0)then
        iaux=ktcaloc(NAUX)
        ilist(1,iaux)=NAUX
        klist(pexln)=iaux
        rlist(iaux+1)=0.d0
      else if(ilist(1,klist(pexln)) .ne. NAUX) then
        call errmsg('filaux','Inconsistent memory fallocation.',0,4)
        return
      else
         iaux=klist(pexln)
      endif
      llen=ilist(1,pexln-1)-2
c      dist=0
c   K. Oide 10/26/1996
c      do 1100 i=pexln+1,pexln+llen
c          if (kytbl(kwL,idtype(ilist(1,i))) .ne. 0) then
c            dist=dist+
c     &           rlist(ilist(2,i)
c     &                 +kytbl(kwL,idtype(ilist(1,i))))
c          endif
c 1100   continue
c      rlist(iaux+1)=dist
c   K. Oide 11/22/1995
      ilist(1,iaux+5)=idx
c
      status = 0
      call v_clear(codin,6)
      call v_clear(codou,6)
c
c.....for radiation loss.
      p0=RgetGL('MOMENTUM',iw)
      pmass=RgetGL('$MASS$',idummy)
      g=p0/pmass
      g=1.0/(1.0+g*g)
      u0=0.d0
c      oldRF=IgetGL('$RFSW$',idummy)
c      call IsetGL('$RFSW$',FLAGOF,idummy)
c
c      call m_unit(tm,6)
c      if(igetgl('$COD$',idummy) .ne. 0)then
c        do 2100 i=1,llen
c          call v_clear(codin,6)
c          de=radlos(pexln+i,codin,p0,pmass,du)
c          u0=u0+du
c          call optelm(pexln+i
c    &                ,codin,tm,codou,tmw,status)
c          call m_copy_6d(tmw,tm)
c2100    continue
c      endif
c      call IsetGL('$RFSW$',oldRF,idummy)
c       write(errfl,'('' Radiation loss per revolution:'',1PG11.4,'//
c    &              ''' eV'')') u0
c
c     rlist(iaux+2)=-tm(5,6)/dist
c      write (errfl,*) 'Design orbit length:',rlist(iaux+1)
c     write (errfl,*) 'momentum compaction :',rlist(iaux+2)
c
      rlist(iaux+3)=u0
      vc=0.d0
      if(abs(u0) .ge. abs(vc)) then
c       call errmsg('expand','Insufficent RF power. ',0,4)
        phis=0.d0
      else
c       if(rlist(iaux+2) .gt. g) then
c         phis=atan2(-u0,-sqrt(vc*vc-u0*u0))
c       else
c         phis=atan2(-u0, sqrt(vc*vc-u0*u0))
c       endif
        phis=0.d0
c       write(errfl,'(A5,G13.6,A2,T25,2(A5,G13.6,A4))')
c    &             'U0=',-u0,'eV' ,'Vc=',vc,'eV  '
c    &             ,'PHIS=',180.0d0*phis/pi,'deg'
      endif
      call RsetGL('PHIS',phis,idummy)
      return
      end
