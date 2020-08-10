      subroutine tquade(trans,cod,beam,srot,al,ak,bz,
     1     dx,dy,theta,krad,
     1     fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0,
     $     kin,achro)
      use tfstk
      use ffs_flag
      use tmacro
      use temw, only:tsetr0
      use sol,only:tsolrote
      use kradlib, only:tradke      
      implicit none
      real*8 , intent(in)::ak,al,bz,dx,dy,theta,
     $     f1in,f2in,f1out,f2out,eps0
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      integer*4 ,intent(in):: mfring
      logical*4 ,intent(in):: fringe
      real*8 bxs,bys,bzs,theta1,ak1,aln,akn
      integer*4 itgetqraddiv,i,ndiv
      integer*4 , parameter :: ndivmax=1000
      logical*4 ,intent(in):: krad,kin,achro
      if(al .eq. 0.d0)then
        call tthine(trans,cod,beam,srot,4,
     $       al,ak,dx,dy,theta,.false.)
        return
      elseif(ak .eq. 0.d0)then
        call tdrife(trans,cod,beam,srot,al,bz,
     $       0.d0,0.d0,0.d0,.true.,.false.,irad)
        return
      endif
      if(ak*al .lt. 0.d0)then
        theta1=theta-m_pi_2
        ak1=-ak
      else
        theta1=theta
        ak1=ak
      endif
      call tsolrote(trans,cod,beam,srot,al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.true.)
c cod has canonical momenta!
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring .ne. 2)then
        call tqfrie(trans,cod,beam,ak1,al,bzs)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,f1in,f2in,bzs)
      endif
      if(krad .and. f1in .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1in,0.d0,bzs*.5d0)
      else
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(krad)then
        if(eps0 .eq. 0.d0)then
          ndiv=max(1,itgetqraddiv(cod,ak1,al,bzs*.5d0))
        else
          ndiv=max(1,
     $         int(dble(itgetqraddiv(cod,ak1,al,bzs*.5d0))/eps0))
        endif
        ndiv=min(ndivmax,ndiv)
      else
        ndiv=1
      endif
      aln=al/ndiv
      akn=ak1/ndiv
      do i=1,ndiv
        call tsolque(trans,cod,beam,srot,aln,akn,
     $       bzs,0.d0,0.d0,eps0,krad,irad)
      enddo
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,-f1out,f2out,bzs)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring .ne. 1)then
        call tqfrie(trans,cod,beam,-ak1,al,bzs)
      endif
      if(krad .and. f1out .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1out,0.d0,bzs*.5d0)
      endif
      call tsolrote(trans,cod,beam,srot,al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.false.)
      return
      end
