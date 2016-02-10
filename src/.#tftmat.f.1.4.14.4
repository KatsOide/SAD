      subroutine tftmat(twiss,gammab,trans,l1,l2,idp,fold)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 l1,l2,idp,i
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
      real*8 trans(4,5)
      real*8 utwiss1(ntwissfun),utwiss2(ntwissfun)
      logical*4 fold
      do i=1,mfitdetr
        utwiss1(i)=twiss(l1,idp,i)
        utwiss2(i)=twiss(l2,idp,i)
      enddo
      call tftmatu(utwiss1,utwiss2,
     $     twiss(nlat,idp,mfitnx),twiss(nlat,idp,mfitny),
     $     gammab,trans,l1,l2,fold,trpt)
      return
      end

      subroutine tftmatu(utwiss1,utwiss2,amux,amuy,
     $     gammab,trans,l1,l2,fold,trpt1)
      use ffs
      use tffitcode
      implicit none
      integer*4 l1,l2
      real*8 utwiss1(*),utwiss2(*),gammab(*)
      real*8 trans(4,5),amux,amuy
      logical*4 fold,trpt1
      integer*4 i
      real*8 detr,cc,x,px,y,py,r1,r2,r3,r4
      call qgettru(utwiss1,utwiss2,amux,amuy,
     $     gammab,trans,l1,l2,.false.,fold,trpt1)
      r1=utwiss1(mfitr1)
      r2=utwiss1(mfitr2)
      r3=utwiss1(mfitr3)
      r4=utwiss1(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(utwiss1(mfitdetr) .lt. 1.d0)then
        do 20 i=1,4
          x=trans(i,1)
          px=trans(i,2)
          y=trans(i,3)
          py=trans(i,4)
          trans(i,1)=cc*x +r1*y+r3*py
          trans(i,2)=cc*px+r2*y+r4*py
          trans(i,3)=cc*y -r4*x+r3*px
          trans(i,4)=cc*py+r2*x-r1*px
20      continue
      else
        do 30 i=1,4
          x=trans(i,1)
          px=trans(i,2)
          y=trans(i,3)
          py=trans(i,4)
          trans(i,1)=cc*y  +(-r4*x+r3*px)
          trans(i,2)=cc*py +( r2*x-r1*px)
          trans(i,3)=cc*x  +( r1*y+r3*trans(i,4))
          trans(i,4)=cc*px +( r2*y+r4*trans(i,4))
30      continue
      endif
      return
      end
