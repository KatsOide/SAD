      subroutine tinitr(trans)
      implicit none
      real*8 ,intent(out):: trans(36)
      trans=(/
     $     1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     $     0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     $     0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,
     $     0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,0.d0,1.d0/)
      return
      end

      subroutine tinitr12(trans)
      implicit none
      real*8 ,intent(out):: trans(72)
      trans=0.d0
      trans(1 )=1.d0
      trans(8 )=1.d0
      trans(15)=1.d0
      trans(22)=1.d0
      trans(29)=1.d0
      trans(36)=1.d0
      return
      end
