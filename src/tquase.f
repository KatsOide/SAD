c Obsolete Aug. 10, 2020
c
      subroutine tquase(trans,cod,beam,srot,al,ak,bz,
     1     dx,dy,theta,radlvl,
     1     fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0)
      use tfstk
      use ffs_flag
      use tmacro
      use temw, only:tsetr0
      use sol,only:tsolrote
      use kradlib, only:tradke      
      implicit none
      real*8 , intent(in)::ak,al,bz,dx,dy,theta,
     $     f1in,f2in,f1out,f2out,eps0,radlvl
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      integer*4 ,intent(in):: mfring
      logical*4 ,intent(in):: fringe
      real*8 bxs,bys,bzs,theta1,ak1,aln,akn
      integer*4 itgetqraddiv,i,ndiv
      integer*4 , parameter :: ndivmax=1000
      logical*4 enarad,krad
      if(al .eq. 0.d0)then
        call tthine(trans,cod,beam,srot,4,al,ak,dx,dy,theta,.false.,l)
        return
      endif
      enarad=radlvl /= 1.d0
      krad=enarad .and. al /= 0.d0
      cod(2)=cod(2)+.5d0*bz*dy
      cod(4)=cod(4)-.5d0*bz*dx
      if(ak*al .lt. 0.d0)then
        theta1=theta-m_pi_2
        ak1=-ak
      else
        theta1=theta
        ak1=ak
      endif
      call tsolrote(trans,cod,beam,srot,0.d0,0.d0,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.true.)
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring /= 2)then
        call tqfrie(trans,cod,beam,ak1,al,bz)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,f1in,f2in,bz)
      endif
      if(krad .and. f1in /= 0.d0)then
        call tradke(trans,cod,beam,srot,f1in,0.d0,bz*.5d0)
      else
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(krad)then
        if(eps0 .eq. 0.d0)then
          ndiv=min(ndivmax,max(1,merge(itgetqraddiv(cod,ak1,al,bzs*.5d0)
        else
          ndiv=int(dble(itgetqraddiv(cod,ak1,al,bzs*.5d0))/eps0)
        endif
      else
        ndiv=1
      endif
      aln=al/ndiv
      akn=ak1/ndiv
      do i=1,ndiv
        call tsolque(trans,cod,beam,srot,aln,akn,
     $       bz,0.d0,0.d0,eps0,krad,irad)
      enddo
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,-f1out,f2out,bz)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring /= 1)then
        call tqfrie(trans,cod,beam,-ak1,al,bz)
      endif
      if(krad .and. f1out /= 0.d0)then
        call tradke(trans,cod,beam,srot,f1out,0.d0,bzs*.5d0)
      endif
      call tsolrote(trans,cod,beam,srot,0.d0,0.d0,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.false.)
      cod(2)=cod(2)-.5d0*bz*dy
      cod(4)=cod(4)+.5d0*bz*dx
      return
      end
