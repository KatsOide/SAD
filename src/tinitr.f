      subroutine tinitr(trans)
      implicit none
      real*8 trans(36)
      trans(1 )=1.d0
      trans(2 )=0.d0
      trans(3 )=0.d0
      trans(4 )=0.d0
      trans(5 )=0.d0
      trans(6 )=0.d0
      trans(7 )=0.d0
      trans(8 )=1.d0
      trans(9 )=0.d0
      trans(10)=0.d0
      trans(11)=0.d0
      trans(12)=0.d0
      trans(13)=0.d0
      trans(14)=0.d0
      trans(15)=1.d0
      trans(16)=0.d0
      trans(17)=0.d0
      trans(18)=0.d0
      trans(19)=0.d0
      trans(20)=0.d0
      trans(21)=0.d0
      trans(22)=1.d0
      trans(23)=0.d0
      trans(24)=0.d0
      trans(25)=0.d0
      trans(26)=0.d0
      trans(27)=0.d0
      trans(28)=0.d0
      trans(29)=1.d0
      trans(30)=0.d0
      trans(31)=0.d0
      trans(32)=0.d0
      trans(33)=0.d0
      trans(34)=0.d0
      trans(35)=0.d0
      trans(36)=1.d0
      return
      end

      subroutine tinitr12(trans)
      implicit none
      real*8 trans(72)
      trans=0.d0
      trans(1 )=1.d0
      trans(8 )=1.d0
      trans(15)=1.d0
      trans(22)=1.d0
      trans(29)=1.d0
      trans(36)=1.d0
      return
      end
