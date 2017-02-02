      subroutine tdrife(trans,cod,beam,al,bz,ak0x,ak0y,
     $     dvon,enarad,calpol,irad,ld)
      use element_drift_common
      use tmacro, only:bradprev
      implicit none
      integer*4 irad,ld,i,itmax
      parameter (itmax=10)
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6)
      real*8 al,bz,ak0x,ak0y,pr,pxi,pyi,pzi,a,ale,alz,
     $     dv,dvdp,bzh,bx,by,br,tbrhoz
      logical*4 dvon,calpol,enarad
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
      if(bz .eq. 0.d0 .and. ak0x .eq. 0.d0 .and. ak0y .eq. 0.d0)then
        pxi=cod(2)
        pyi=cod(4)
        a=min(ampmax,pxi**2+pyi**2)
        pzi=pr*sqrt(1.d0-a/pr**2)
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
        if(irad .gt. 6 .or. calpol)then
          call tmulbs(beam ,trans1,.true.,.true.)
        endif
        if(calpol)then
          call polpar(10,ld,al,0.d0,0.d0,0.d0,0.d0,cod)
        endif
        cod(1)=cod(1)+pxi/pzi*al
        cod(3)=cod(3)+pyi/pzi*al
        cod(5)=cod(5)-(a/(pr+pzi)/pzi+dv)*al
      else
        if(enarad)then
          br=tbrhoz()
          bx=ak0y/al*br
          by=ak0x/al*br
c   dl/dx is temporarily set to zero, of course it is wrong...
          call trade(trans,beam,cod,bx,by,bz*br,bz,
     $         0.d0,0.d0,0.d0,0.d0,
     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
        else
          br=0.d0
        endif
        bzh=bz*.5d0
        cod(2)=cod(2)+bzh*cod(3)
        cod(4)=cod(4)-bzh*cod(1)
        call tsoldz(trans1,cod,al,ak0y/al,ak0x/al,bz,.true.)
        trans1(1,3)=trans1(1,3)+bzh*trans1(1,2)
        trans1(1,1)=trans1(1,1)-bzh*trans1(1,4)
        trans1(2,3)=trans1(2,3)+bzh*trans1(2,2)
        trans1(2,1)=trans1(2,1)-bzh*trans1(2,4)
        trans1(3,3)=trans1(3,3)+bzh*trans1(3,2)
        trans1(3,1)=trans1(3,1)-bzh*trans1(3,4)
        trans1(4,3)=trans1(4,3)+bzh*trans1(4,2)
        trans1(4,1)=trans1(4,1)-bzh*trans1(4,4)
        trans1(5,3)=trans1(5,3)+bzh*trans1(5,2)
        trans1(5,1)=trans1(5,1)-bzh*trans1(5,4)
        cod(2)=cod(2)-bzh*cod(3)
        cod(4)=cod(4)+bzh*cod(1)
        trans1(2,1)=trans1(2,1)-bzh*trans1(3,1)
        trans1(4,1)=trans1(4,1)+bzh*trans1(1,1)
        trans1(2,2)=trans1(2,2)-bzh*trans1(3,2)
        trans1(4,2)=trans1(4,2)+bzh*trans1(1,2)
        trans1(2,3)=trans1(2,3)-bzh*trans1(3,3)
        trans1(4,3)=trans1(4,3)+bzh*trans1(1,3)
        trans1(2,4)=trans1(2,4)-bzh*trans1(3,4)
        trans1(4,4)=trans1(4,4)+bzh*trans1(1,4)
        trans1(2,6)=trans1(2,6)-bzh*trans1(3,6)
        trans1(4,6)=trans1(4,6)+bzh*trans1(1,6)
        call tmultr5(trans,trans1,irad)
        if(enarad)then
          call trade(trans,beam,cod,bx,by,bz*br,bz,
     $         0.d0,0.d0,0.d0,0.d0,
     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
        endif
        if(irad .gt. 6 .or. calpol)then
          call tmulbs(beam ,trans1,.true.,.true.)
        endif
        if(calpol)then
          call polpar(10,ld,al,0.d0,0.d0,0.d0,0.d0,cod)
        endif
      endif
      bradprev=0.d0
      return
      end
