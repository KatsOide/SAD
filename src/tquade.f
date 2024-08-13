      subroutine tquade(trans,cod,beam,srot,al,ak,bz,
     1     dx,dy,theta,krad,
     1     fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0,kin,achro,l)
      use tfstk
      use ffs_flag
      use tmacro
      use drife
      use temw, only:tsetr0
      use sol,only:tsolrote
      use kradlib, only:tradke      
      use tparastat,only:setndivelm
      implicit none
      real*8 , intent(in)::ak,al,bz,dx,dy,theta,
     $     f1in,f2in,f1out,f2out,eps0
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      integer*4 ,intent(in):: mfring,l
      logical*4 ,intent(in):: fringe
      real*8 bxs,bys,bzs,theta1,ak1,aln,akn
      integer*4 itgetqraddiv,i,ndiv
      integer*4 , parameter :: ndivmax=1000
      logical*4 ,intent(in):: krad,kin,achro
      if(al == 0.d0)then
        call tthine(trans,cod,beam,srot,4,al,ak,dx,dy,theta,.false.,l)
        return
      elseif(ak == 0.d0)then
        call tdrife0(trans,cod,beam,srot,al,bz,0.d0,.true.,.false.,irad)
        return
      endif
      if(ak*al .lt. 0.d0)then
        theta1=theta-m_pi_2
        ak1=-ak
      else
        theta1=theta
        ak1=ak
      endif
      call tsolrote(trans,cod,beam,srot,0.d0,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.true.)
c cod has canonical momenta!
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring /= 2)then
        call tqfrie(trans,cod,beam,ak1,al,bzs)
      endif
      if(mfring == 1 .or. mfring == 3)then
        call tqlfre(trans,cod,beam,al,ak1,f1in,f2in,bzs)
      endif
      if(krad .and. f1in /= 0.d0)then
        call tradke(trans,cod,beam,srot,f1in,0.d0,bzs*.5d0)
      else
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(krad)then
        if(eps0 == 0.d0)then
          ndiv=min(ndivmax,max(1,itgetqraddiv(cod,ak1,al,bzs*.5d0)))
        else
          ndiv=min(ndivmax,max(1,int(dble(itgetqraddiv(cod,ak1,al,bzs*.5d0)/eps0))))
        endif
      else
        ndiv=1
      endif
      call setndivelm(l,0)
      aln=al/ndiv
      akn=ak1/ndiv
      do i=1,ndiv
        call tsolque(trans,cod,beam,srot,aln,akn,
     $       bzs,0.d0,0.d0,eps0,krad,irad,l)
      enddo
      if(mfring == 2 .or. mfring == 3)then
        call tqlfre(trans,cod,beam,al,ak1,-f1out,f2out,bzs)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring /= 1)then
        call tqfrie(trans,cod,beam,-ak1,al,bzs)
      endif
      if(krad .and. f1out /= 0.d0)then
        call tradke(trans,cod,beam,srot,f1out,0.d0,bzs*.5d0)
      endif
      call tsolrote(trans,cod,beam,srot,0.d0,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.false.)
      return
      end
