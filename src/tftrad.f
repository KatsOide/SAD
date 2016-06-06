      subroutine tftrad(trdtbl,trval,lfno)
      use tfstk
      implicit none
      real*8 trdtbl(3,6),trval,range(3,3),r1(3),phi(3)
      integer*4 lfno
      range(1:2,1)=trdtbl(1:2,1)
      range(1:2,2)=trdtbl(1:2,2)
      r1=trdtbl(1:3,3)
      phi(1:2)=0.d0
      phi(3)=-0.5d0*pi
      call tfdapert1(range,r1,3,trval,
     $     phi,0.d0,3,1,lfno)
      return
      end
