      module drife
      use sad_basics
      use mathfun, only: sqrtl

      contains
      subroutine tdrife(trans,cod,beam,srot,al,bz,ak0x,ak0y,alr,dvon,enarad,irad)
      use element_drift_common
      use tmacro, only:bradprev
      use temw, only:tsetr0,tmulbs
      use tspin
      use kradlib, only:tradke
      use sad_basics
      implicit none
      integer*4 ,intent(in):: irad
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in)::  al,bz,ak0x,ak0y
      real*8 trans1(6,6),pr,dv,dvdp,bzh,alr,beam0(42)
      logical*4 ,intent(in):: dvon,enarad
      if(al == 0.d0)then
        cod(2)=cod(2)-ak0x
        cod(4)=cod(4)+ak0y
        return
      endif
      if(dvon)then
        call tgetdv(cod(6),dv,dvdp)
      else
        dv=0.d0
        dvdp=0.d0
      endif
      pr=1.d0+cod(6)
      if(abs(bz) == 0.d0 .and.
     $     abs(ak0x) == 0.d0 .and. abs(ak0y) == 0.d0)then
        call tdrife1(trans1,cod,dv,dvdp,al)
        trans(1:5:2,1:irad)=trans(1:5:2,1:irad)
     $       +matmul(trans1(1:5:2,2:6:2),trans(2:6:2,1:irad))
        bzh=0.d0
      else
        call tinitr(trans1)
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
c        writen(*,'(a,1p10g12.4)')'tdrife ',al,ak0y/al,ak0x/al,bz,trans1(1,2)
        call tmultr5(trans,trans1,irad)
      endif
      beam0=beam
      if(irad > 6)then
        call tmulbs(beam ,trans1,.true.)
      endif
      if(enarad)then
        call tradke(trans,cod,beam,srot,alr,0.d0,bzh)
      endif
      bradprev=0.d0
      return
      end

      subroutine tdrife0(trans,cod,beam,srot,al,bz,alr,dvon,enarad,irad)
      use element_drift_common
      use tmacro, only:bradprev
      use temw, only:tsetr0,tmulbs
      use tspin
      use kradlib, only:tradke
      implicit none
      integer*4 ,intent(in):: irad
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: al,bz
      real*8 trans1(6,6),pr,dv,dvdp,bzh,alr,beam0(42)
      logical*4 ,intent(in):: dvon,enarad
      if(al == 0.d0)then
        return
      endif
      if(dvon)then
        call tgetdv(cod(6),dv,dvdp)
      else
        dv=0.d0
        dvdp=0.d0
      endif
      pr=1.d0+cod(6)
      if(abs(bz) == 0.d0)then
        call tdrife1(trans1,cod,dv,dvdp,al)
        trans(1:5:2,1:irad)=trans(1:5:2,1:irad)+matmul(trans1(1:5:2,2:6:2),trans(2:6:2,1:irad))
        bzh=0.d0
      else
        call tinitr(trans1)
        bzh=bz*.5d0
c cod has canonical momenta!
        cod(2)=cod(2)+bzh*cod(3)
        cod(4)=cod(4)-bzh*cod(1)
        call tsoldz(trans1,cod,al,0.d0,0.d0,bz,.true.)
        trans1(1:5,3)=trans1(1:5,3)+bzh*trans1(1:5,2)
        trans1(1:5,1)=trans1(1:5,1)-bzh*trans1(1:5,4)
        cod(2)=cod(2)-bzh*cod(3)
        cod(4)=cod(4)+bzh*cod(1)
        trans1(2,1:6)=trans1(2,1:6)-bzh*trans1(3,1:6)
        trans1(4,1:6)=trans1(4,1:6)+bzh*trans1(1,1:6)
c        writen(*,'(a,1p10g12.4)')'tdrife ',al,ak0y/al,ak0x/al,bz,trans1(1,2)
        call tmultr5(trans,trans1,irad)
      endif
      beam0=beam
      if(irad > 6)then
        call tmulbs(beam ,trans1,.true.)
      endif
      if(enarad)then
        call tradke(trans,cod,beam,srot,alr,0.d0,bzh)
      endif
      bradprev=0.d0
      return
      end

      end module drife
