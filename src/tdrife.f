      subroutine tdrife(trans,cod,beam,srot,al,bz,ak0x,ak0y,alr,
     $     dvon,enarad,irad)
      use element_drift_common
      use tmacro, only:bradprev
      use temw, only:tsetr0,tmulbs
      use sol, only:tsoldz
      use tspin
      use kradlib, only:tradke
      use mathfun, only: sqrtl
      implicit none
      integer*4 ,intent(in):: irad
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in)::  al,bz,ak0x,ak0y
      real*8 trans1(6,6),pr,pxi,pyi,pzi,a,ale,alz,dv,dvdp,bzh,alr
      logical*4 ,intent(in):: dvon,enarad
      if(al .eq. 0.d0)then
        cod(2)=cod(2)-ak0x
        cod(4)=cod(4)+ak0y
        return
      endif
      call tinitr(trans1)
      if(dvon)then
        call tgetdv(cod(6),dv,dvdp)
      else
        dv=0.d0
        dvdp=0.d0
      endif
      pr=1.d0+cod(6)
      if(abs(bz) .eq. 0.d0 .and.
     $     abs(ak0x) .eq. 0.d0 .and. abs(ak0y) .eq. 0.d0)then
        pxi=cod(2)
        pyi=cod(4)
        a=pxi**2+pyi**2
        pzi=pr*sqrtl(1.d0-a/pr**2)
        ale=al/pzi
        alz=ale/pzi**2
        trans1(1,2)=ale+pxi**2*alz
        trans1(1,4)=pxi*pyi*alz
        trans1(1,6)=-pxi*pr*alz
        trans1(3,2)=trans1(1,4)
        trans1(3,4)=ale+pyi**2*alz
        trans1(3,6)=-pyi*pr*alz
        trans1(5,2)=trans1(1,6)
        trans1(5,4)=trans1(3,6)
        trans1(5,6)=dvdp*al+a*alz
        trans(1:5:2,1:irad)=trans(1:5:2,1:irad)
     $       +matmul(trans1(1:5:2,2:6:2),trans(2:6:2,1:irad))
c        trans(1,1:irad)=trans(1,1:irad)+trans1(1,2)*trans(2,1:irad)
c     $       +trans1(1,4)*trans(4,1:irad)+trans1(1,6)*trans(6,1:irad)
c        trans(3,1:irad)=trans(3,1:irad)+trans1(3,2)*trans(2,1:irad)
c     $       +trans1(3,4)*trans(4,1:irad)+trans1(3,6)*trans(6,1:irad)
c        trans(5,1:irad)=trans(5,1:irad)+trans1(5,2)*trans(2,1:irad)
c     $       +trans1(5,4)*trans(4,1:irad)+trans1(5,6)*trans(6,1:irad)
        cod(1)=cod(1)+pxi/pzi*al
        cod(3)=cod(3)+pyi/pzi*al
        cod(5)=cod(5)-(a/(pr+pzi)/pzi+dv)*al
        bzh=0.d0
      else
        bzh=bz*.5d0
c cod has canonical momenta!
        cod(2)=cod(2)+bzh*cod(3)
        cod(4)=cod(4)-bzh*cod(1)
        call tsoldz(trans1,cod,al,ak0y/al,ak0x/al,bz,.true.)
        trans1(1:5,3)=trans1(1:5,3)+bzh*trans1(1:5,2)
        trans1(1:5,1)=trans1(1:5,1)-bzh*trans1(1:5,4)
        cod(2)=cod(2)-bzh*cod(3)
        cod(4)=cod(4)+bzh*cod(1)
        trans1(2,1:6)=trans1(2,1:6)-bzh*trans1(3,1:6)
        trans1(4,1:6)=trans1(4,1:6)+bzh*trans1(1,1:6)
        call tmultr5(trans,trans1,irad)
      endif
      if(irad .gt. 6)then
        call tmulbs(beam ,trans1,.true.)
      endif
      if(enarad)then
        call tradke(trans,cod,beam,srot,alr,0.d0,bzh)
      endif
      bradprev=0.d0
      return
      end
