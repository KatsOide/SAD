      subroutine qgettr(trans1,k1,k2,idp,disp,fold)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in):: k1,k2,idp
      real*8 ,intent(inout):: trans1(4,5)
      logical*4 ,intent(in):: disp,fold
      real*8 utwiss1(ntwissfun),utwiss2(ntwissfun)
c      do i=1,mfitdetr
        utwiss1(1:mfitdetr)=twiss(k1,idp,1:mfitdetr)
        utwiss2(1:mfitdetr)=twiss(k2,idp,1:mfitdetr)
c      enddo
      call qgettru(utwiss1,utwiss2,
     $     twiss(nlat,idp,mfitnx),twiss(nlat,idp,mfitny),
     $     trans1,k1,k2,disp,fold,trpt)
      return
      end

      subroutine qgettru(utwiss1,utwiss2,amux,amuy,
     $     trans1,k1,k2,disp,fold,trpt1)
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 k1,k2
      real*8 ,intent(in):: utwiss1(*),utwiss2(*),amux,amuy
      real*8 ,intent(inout):: trans1(4,5)
      real*8 trans(4,5)
      logical*4 ,intent(in):: disp,fold,trpt1
      real*8 sqrtb1,sqrtb2,dpsi1,cospsi,sinpsi,gr,
     $     detr,cc,r1,r2,r3,r4,th
      sqrtb1=sqrt(utwiss1(mfitbx))
      sqrtb2=sqrt(utwiss2(mfitbx))
      dpsi1=utwiss2(mfitnx)-utwiss1(mfitnx)
c next line modified by Kikuchi. 7/28/93
c     if(k2 .lt. k1)then
      if(.not. trpt1 .and.
     $     (k2 .lt. k1 .or. k2 .eq. k1 .and. fold))then
         dpsi1=amux+dpsi1
      endif
      th=tan(.5d0*dpsi1)
      sinpsi=2.d0*th/(1.d0+th**2)
      cospsi=1.d0-th*sinpsi
c      sinpsi=sin(dpsi1)
c      cospsi=cos(dpsi1)
      trans=0.d0
      trans(1,1)=sqrtb2/sqrtb1*(cospsi+utwiss1(mfitax)*sinpsi)
      trans(1,2)=sqrtb2*sqrtb1*sinpsi
      trans(2,1)=-((utwiss2(mfitax)-utwiss1(mfitax))*cospsi
     1             +(utwiss1(mfitax)*utwiss2(mfitax)+1.d0)*sinpsi)
     1           /sqrtb1/sqrtb2
      trans(2,2)=sqrtb1/sqrtb2*(cospsi-utwiss2(mfitax)*sinpsi)
      if(disp)then
        gr=1
        trans(1,5)=utwiss2(mfitex)
        trans(2,5)=utwiss2(mfitepx)
      else
        gr=sqrt(gammab(k2)/gammab(k1))
        trans(1,5)=utwiss2(mfitex)-gr*(trans(1,1)*utwiss1(mfitex)
     1                                +trans(1,2)*utwiss1(mfitepx))
        trans(2,5)=utwiss2(mfitepx)-gr*(trans(2,1)*utwiss1(mfitex)
     1                                +trans(2,2)*utwiss1(mfitepx))
      endif
      sqrtb1=sqrt(utwiss1(mfitby))
      sqrtb2=sqrt(utwiss2(mfitby))
      dpsi1=utwiss2(mfitny)-utwiss1(mfitny)
c      if(k2 .lt. k1)then
      if(.not. trpt1 .and.
     $     (k2 .lt. k1 .or. k2 .eq. k1 .and. fold))then
        dpsi1=amuy+dpsi1
      endif
      th=tan(.5d0*dpsi1)
      sinpsi=2.d0*th/(1.d0+th**2)
      cospsi=1.d0-th*sinpsi
c      cospsi=cos(dpsi1)
c      sinpsi=sin(dpsi1)
      trans(3,3)=sqrtb2/sqrtb1*(cospsi+utwiss1(mfitay)*sinpsi)
      trans(3,4)=sqrtb2*sqrtb1*sinpsi
      trans(4,3)=-((utwiss2(mfitay)-utwiss1(mfitay))*cospsi
     1             +(utwiss1(mfitay)*utwiss2(mfitay)+1.d0)*sinpsi)
     1           /sqrtb1/sqrtb2
      trans(4,4)=sqrtb1/sqrtb2*(cospsi-utwiss2(mfitay)*sinpsi)
      if(disp)then
        trans(3,5)=utwiss2(mfitey)
        trans(4,5)=utwiss2(mfitepy)
      else
        trans(3,5)=utwiss2(mfitey)-gr*(trans(3,3)*utwiss1(mfitey)
     1                                +trans(3,4)*utwiss1(mfitepy))
        trans(4,5)=utwiss2(mfitepy)-gr*(trans(4,3)*utwiss1(mfitey)
     1                                 +trans(4,4)*utwiss1(mfitepy))
      endif
      r1=utwiss2(mfitr1)
      r2=utwiss2(mfitr2)
      r3=utwiss2(mfitr3)
      r4=utwiss2(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(utwiss2(mfitdetr) .lt. 1.d0)then
        trans1(1,1)=cc*trans(1,1)+r4*trans(3,1)-r2*trans(4,1)
        trans1(2,1)=cc*trans(2,1)-r3*trans(3,1)+r1*trans(4,1)
        trans1(3,1)=cc*trans(3,1)-r1*trans(1,1)-r2*trans(2,1)
        trans1(4,1)=cc*trans(4,1)-r3*trans(1,1)-r4*trans(2,1)
        trans1(1,2)=cc*trans(1,2)+r4*trans(3,2)-r2*trans(4,2)
        trans1(2,2)=cc*trans(2,2)-r3*trans(3,2)+r1*trans(4,2)
        trans1(3,2)=cc*trans(3,2)-r1*trans(1,2)-r2*trans(2,2)
        trans1(4,2)=cc*trans(4,2)-r3*trans(1,2)-r4*trans(2,2)
        trans1(1,3)=cc*trans(1,3)+r4*trans(3,3)-r2*trans(4,3)
        trans1(2,3)=cc*trans(2,3)-r3*trans(3,3)+r1*trans(4,3)
        trans1(3,3)=cc*trans(3,3)-r1*trans(1,3)-r2*trans(2,3)
        trans1(4,3)=cc*trans(4,3)-r3*trans(1,3)-r4*trans(2,3)
        trans1(1,4)=cc*trans(1,4)+r4*trans(3,4)-r2*trans(4,4)
        trans1(2,4)=cc*trans(2,4)-r3*trans(3,4)+r1*trans(4,4)
        trans1(3,4)=cc*trans(3,4)-r1*trans(1,4)-r2*trans(2,4)
        trans1(4,4)=cc*trans(4,4)-r3*trans(1,4)-r4*trans(2,4)
        trans1(1,5)=cc*trans(1,5)+r4*trans(3,5)-r2*trans(4,5)
        trans1(2,5)=cc*trans(2,5)-r3*trans(3,5)+r1*trans(4,5)
        trans1(3,5)=cc*trans(3,5)-r1*trans(1,5)-r2*trans(2,5)
        trans1(4,5)=cc*trans(4,5)-r3*trans(1,5)-r4*trans(2,5)
      else
        trans1(1,1)=cc*trans(3,1)-r1*trans(1,1)-r2*trans(2,1)
        trans1(2,1)=cc*trans(4,1)-r3*trans(1,1)-r4*trans(2,1)
        trans1(3,1)=cc*trans(1,1)+r4*trans(3,1)-r2*trans(4,1)
        trans1(4,1)=cc*trans(2,1)-r3*trans(3,1)+r1*trans(4,1)
        trans1(1,2)=cc*trans(3,2)-r1*trans(1,2)-r2*trans(2,2)
        trans1(2,2)=cc*trans(4,2)-r3*trans(1,2)-r4*trans(2,2)
        trans1(3,2)=cc*trans(1,2)+r4*trans(3,2)-r2*trans(4,2)
        trans1(4,2)=cc*trans(2,2)-r3*trans(3,2)+r1*trans(4,2)
        trans1(1,3)=cc*trans(3,3)-r1*trans(1,3)-r2*trans(2,3)
        trans1(2,3)=cc*trans(4,3)-r3*trans(1,3)-r4*trans(2,3)
        trans1(3,3)=cc*trans(1,3)+r4*trans(3,3)-r2*trans(4,3)
        trans1(4,3)=cc*trans(2,3)-r3*trans(3,3)+r1*trans(4,3)
        trans1(1,4)=cc*trans(3,4)-r1*trans(1,4)-r2*trans(2,4)
        trans1(2,4)=cc*trans(4,4)-r3*trans(1,4)-r4*trans(2,4)
        trans1(3,4)=cc*trans(1,4)+r4*trans(3,4)-r2*trans(4,4)
        trans1(4,4)=cc*trans(2,4)-r3*trans(3,4)+r1*trans(4,4)
        trans1(1,5)=cc*trans(3,5)-r1*trans(1,5)-r2*trans(2,5)
        trans1(2,5)=cc*trans(4,5)-r3*trans(1,5)-r4*trans(2,5)
        trans1(3,5)=cc*trans(1,5)+r4*trans(3,5)-r2*trans(4,5)
        trans1(4,5)=cc*trans(2,5)-r3*trans(3,5)+r1*trans(4,5)
      endif
      return
      end
