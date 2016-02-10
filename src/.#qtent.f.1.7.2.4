      subroutine qtent(trans,cod,twiss,k1,idp,disp)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k1,idp
      real*8 trans(4,5),twiss(nlat,-ndim:ndim,ntwissfun),
     $     cod(6),r1,r2,r3,r4,detr,cc
      logical*4 disp,normal
      r1=twiss(k1,idp,mfitr1)
      r2=twiss(k1,idp,mfitr2)
      r3=twiss(k1,idp,mfitr3)
      r4=twiss(k1,idp,mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twiss(k1,idp,mfitdetr) .lt. 1.d0)then
        trans(1,1)=cc
        trans(1,2)=0.d0
        trans(1,3)= r4
        trans(1,4)=-r2
        trans(2,1)=0.d0
        trans(2,2)=cc
        trans(2,3)=-r3
        trans(2,4)= r1
        trans(3,1)=-r1
        trans(3,2)=-r2
        trans(3,3)=cc
        trans(3,4)=0.d0
        trans(4,1)=-r3
        trans(4,2)=-r4
        trans(4,3)=0.d0
        trans(4,4)=cc
      else
        trans(1,1)= r4
        trans(1,2)=-r2
        trans(1,3)= cc
        trans(1,4)= 0.d0
        trans(2,1)=-r3
        trans(2,2)= r1
        trans(2,3)=0.d0
        trans(2,4)= cc
        trans(3,1)= cc
        trans(3,2)=0.d0
        trans(3,3)=-r1
        trans(3,4)=-r2
        trans(4,1)=0.d0
        trans(4,2)= cc
        trans(4,3)=-r3
        trans(4,4)=-r4
      endif
      if(disp)then
        if(normal)then
          trans(1,5)=cc*twiss(k1,idp,mfitex)
     1              +r4*twiss(k1,idp,mfitey)-r2*twiss(k1,idp,mfitepy)
          trans(2,5)=cc*twiss(k1,idp,mfitepx)
     1              -r3*twiss(k1,idp,mfitey)+r1*twiss(k1,idp,mfitepy)
          trans(3,5)=cc*twiss(k1,idp,mfitey)
     1              -r1*twiss(k1,idp,mfitex)-r2*twiss(k1,idp,mfitepx)
          trans(4,5)=cc*twiss(k1,idp,mfitepy)
     1              -r3*twiss(k1,idp,mfitex)-r4*twiss(k1,idp,mfitepx)
        else
          trans(1,5)=r4*twiss(k1,idp,mfitex)-r2*twiss(k1,idp,mfitepx)
     1              +cc*twiss(k1,idp,mfitey)
          trans(2,5)=-r3*twiss(k1,idp,mfitex)+r1*twiss(k1,idp,mfitepx)
     1              +cc*twiss(k1,idp,mfitepy)
          trans(3,5)=cc*twiss(k1,idp,mfitex)
     1              -r1*twiss(k1,idp,mfitey)-r2*twiss(k1,idp,mfitepy)
          trans(4,5)=cc*twiss(k1,idp,mfitepx)
     1              -r3*twiss(k1,idp,mfitey)-r4*twiss(k1,idp,mfitepy)
        endif
      else
        trans(1,5)=0.d0
        trans(2,5)=0.d0
        trans(3,5)=0.d0
        trans(4,5)=0.d0
      endif
      cod(1)=twiss(k1,idp,mfitdx)
      cod(2)=twiss(k1,idp,mfitdpx)
      cod(3)=twiss(k1,idp,mfitdy)
      cod(4)=twiss(k1,idp,mfitdpy)
      cod(5)=twiss(k1,idp,mfitdz)
      cod(6)=twiss(k1,idp,mfitddp)
      return
      end

      subroutine qcanontoangle(trans,cod)
      implicit none
      real*8 trans(4,5),cod(6),pr
      pr=1.d0+cod(6)
      cod(2)=cod(2)/pr
      cod(4)=cod(4)/pr
      trans(2,1)=trans(2,1)/pr
      trans(2,2)=trans(2,2)/pr
      trans(2,3)=trans(2,3)/pr
      trans(2,4)=trans(2,4)/pr
      trans(2,5)=(trans(2,5)-cod(2))/pr
      trans(4,1)=trans(4,1)/pr
      trans(4,2)=trans(4,2)/pr
      trans(4,3)=trans(4,3)/pr
      trans(4,4)=trans(4,4)/pr
      trans(4,5)=(trans(4,5)-cod(4))/pr
      return
      end

      subroutine qangletocanon(trans,cod)
      implicit none
      real*8 trans(4,5),cod(6),pr
      pr=1.d0+cod(6)
      trans(2,1)=trans(2,1)*pr
      trans(2,2)=trans(2,2)*pr
      trans(2,3)=trans(2,3)*pr
      trans(2,4)=trans(2,4)*pr
      trans(2,5)=trans(2,5)*pr+cod(2)
      trans(4,1)=trans(4,1)*pr
      trans(4,2)=trans(4,2)*pr
      trans(4,3)=trans(4,3)*pr
      trans(4,4)=trans(4,4)*pr
      trans(4,5)=trans(4,5)*pr+cod(4)
      cod(2)=cod(2)*pr
      cod(4)=cod(4)*pr
      return
      end

      subroutine qtentu(trans,cod,utwiss,disp)
      use ffs
      use tffitcode
      implicit none
      real*8 trans(4,5),cod(6),utwiss(*),
     $     r1,r2,r3,r4,detr,cc
      logical*4 disp,normal
      r1=utwiss(mfitr1)
      r2=utwiss(mfitr2)
      r3=utwiss(mfitr3)
      r4=utwiss(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      normal=utwiss(mfitdetr) .lt. 1.d0
      if(normal)then
        trans(1,1)=cc
        trans(1,2)=0.d0
        trans(1,3)= r4
        trans(1,4)=-r2
        trans(2,1)=0.d0
        trans(2,2)=cc
        trans(2,3)=-r3
        trans(2,4)= r1
        trans(3,1)=-r1
        trans(3,2)=-r2
        trans(3,3)=cc
        trans(3,4)=0.d0
        trans(4,1)=-r3
        trans(4,2)=-r4
        trans(4,3)=0.d0
        trans(4,4)=cc
      else
        trans(1,1)= r4
        trans(1,2)=-r2
        trans(1,3)= cc
        trans(1,4)= 0.d0
        trans(2,1)=-r3
        trans(2,2)= r1
        trans(2,3)=0.d0
        trans(2,4)= cc
        trans(3,1)= cc
        trans(3,2)=0.d0
        trans(3,3)=-r1
        trans(3,4)=-r2
        trans(4,1)=0.d0
        trans(4,2)= cc
        trans(4,3)=-r3
        trans(4,4)=-r4
      endif
      if(disp)then
        if(normal)then
          trans(1,5)=cc*utwiss(7)
     1              +r4*utwiss(9)-r2*utwiss(10)
          trans(2,5)=cc*utwiss(8)
     1              -r3*utwiss(9)+r1*utwiss(10)
          trans(3,5)=cc*utwiss(9)
     1              -r1*utwiss(7)-r2*utwiss(8)
          trans(4,5)=cc*utwiss(10)
     1              -r3*utwiss(7)-r4*utwiss(8)
        else
          trans(1,5)=-r1*utwiss(7)-r2*utwiss(8)
     1              +cc*utwiss(9)
          trans(2,5)=-r3*utwiss(7)-r4*utwiss(8)
     1              +cc*utwiss(10)
          trans(3,5)=cc*utwiss(7)
     1              +r4*utwiss(9)-r2*utwiss(10)
          trans(4,5)=cc*utwiss(8)
     1              -r3*utwiss(9)+r1*utwiss(10)
c             write(*,'(1P5G15.7)')((trans(ii,jj),jj=1,5),ii=1,4)
        endif
      else
        trans(1,5)=0.d0
        trans(2,5)=0.d0
        trans(3,5)=0.d0
        trans(4,5)=0.d0
      endif
      cod(1)=utwiss(mfitdx)
      cod(2)=utwiss(mfitdpx)
      cod(3)=utwiss(mfitdy)
      cod(4)=utwiss(mfitdpy)
      cod(5)=utwiss(mfitdz)
      cod(6)=utwiss(mfitddp)
      return
      end

      subroutine qinicanon(trans,cod)
      implicit none
      real*8 trans(6,6),cod(6),pr
      pr=1.d0+cod(6)
      trans(1,1)=1.d0
      trans(1,2)=0.d0
      trans(1,3)=0.d0
      trans(1,4)=0.d0
      trans(1,5)=0.d0
      trans(1,6)=0.d0
      trans(2,1)=0.d0
      trans(2,2)=pr
      trans(2,3)=0.d0
      trans(2,4)=0.d0
      trans(2,5)=0.d0
      trans(2,6)=cod(2)
      trans(3,1)=0.d0
      trans(3,2)=0.d0
      trans(3,3)=1.d0
      trans(3,4)=0.d0
      trans(3,5)=0.d0
      trans(3,6)=0.d0
      trans(4,1)=0.d0
      trans(4,2)=0.d0
      trans(4,3)=0.d0
      trans(4,4)=pr
      trans(4,5)=0.d0
      trans(4,6)=cod(4)
      trans(5,1)=0.d0
      trans(5,2)=0.d0
      trans(5,3)=0.d0
      trans(5,4)=0.d0
      trans(5,5)=1.d0
      trans(5,6)=0.d0
      trans(6,1)=0.d0
      trans(6,2)=0.d0
      trans(6,3)=0.d0
      trans(6,4)=0.d0
      trans(6,5)=0.d0
      trans(6,6)=1.d0
      cod(2)=cod(2)*pr
      cod(4)=cod(4)*pr
      return
      end
