      subroutine tdrife(trans,cod,beam,srot,al,bz,ak0x,ak0y,alr,
     $     dvon,enarad,irad)
      use tfstk, only: sqrtl
      use element_drift_common
      use tmacro, only:bradprev
      use temw, only:tsetr0,bsi
      implicit none
      integer*4 irad,i,itmax
      parameter (itmax=10)
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6),srot(3,3)
      real*8 al,bz,ak0x,ak0y,pr,pxi,pyi,pzi,a,ale,alz,
     $     dv,dvdp,bzh,alr
      logical*4 dvon,enarad
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
        a=min(ampmax,pxi**2+pyi**2)
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
        do i=1,irad
          trans(1,i)=trans(1,i)+trans1(1,2)*trans(2,i)
     $         +trans1(1,4)*trans(4,i)+trans1(1,6)*trans(6,i)
          trans(3,i)=trans(3,i)+trans1(3,2)*trans(2,i)
     $         +trans1(3,4)*trans(4,i)+trans1(3,6)*trans(6,i)
          trans(5,i)=trans(5,i)+trans1(5,2)*trans(2,i)
     $         +trans1(5,4)*trans(4,i)+trans1(5,6)*trans(6,i)
        enddo
        if(irad .gt. 6)then
          call tmulbs(beam ,trans1,.true.,.true.)
        endif
        cod(1)=cod(1)+pxi/pzi*al
        cod(3)=cod(3)+pyi/pzi*al
        cod(5)=cod(5)-(a/(pr+pzi)/pzi+dv)*al
      else
c        if(enarad)then
c          br=tbrhoz()
c          bx=ak0y/al*br
c          by=ak0x/al*br
cc   dl/dx is temporarily set to zero, of course it is wrong...
c          call trade(trans,beam,cod,bx,by,bz*br,bz,
c     $         0.d0,0.d0,0.d0,0.d0,
c     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
c        else
c          br=0.d0
c        endif
        bzh=bz*.5d0
        if(enarad)then
          call tsetr0(trans(:,1:6),cod(1:6),bzh,0.d0)
        endif
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
        if(irad .gt. 6)then
          call tmulbs(beam ,trans1,.false.,.true.)
        endif
        if(enarad)then
          call tradke(trans,cod,beam,srot,alr,0.d0,bzh)
        endif
c        if(enarad)then
c          call trade(trans,beam,cod,bx,by,bz*br,bz,
c     $         0.d0,0.d0,0.d0,0.d0,
c     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
c        endif
c        if(irad .gt. 6)then
c          call tmulbs(beam ,trans1,.true.,.true.)
c        endif
      endif
      bradprev=0.d0
      return
      end
