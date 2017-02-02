      subroutine tgrot(a,geo0,geo1)
      implicit none
      real*8 geo0(3,3),geo1(3,3),a(3),chi2,a13,a33,chi1,a21,a22,chi3
      chi2=asin(
     1    -geo0(1,2)*geo1(1,3)-geo0(2,2)*geo1(2,3)-geo0(3,2)*geo1(3,3))
      a13= geo0(1,1)*geo1(1,3)+geo0(2,1)*geo1(2,3)+geo0(3,1)*geo1(3,3)
      if(a13 .eq. 0.d0)then
        chi1=0.d0
      else
        a33= geo0(1,3)*geo1(1,3)+geo0(2,3)*geo1(2,3)+geo0(3,3)*geo1(3,3)
        chi1=-atan2(a13,a33)
      endif
      a21= geo0(1,2)*geo1(1,1)+geo0(2,2)*geo1(2,1)+geo0(3,2)*geo1(3,1)
      if(a21 .eq. 0.d0)then
        chi3=0.d0
      else
        a22= geo0(1,2)*geo1(1,2)+geo0(2,2)*geo1(2,2)+geo0(3,2)*geo1(3,2)
        chi3=atan2(-a21,a22)
      endif
      a(1)=chi1
      a(2)=chi2
      a(3)=chi3
      return
      end
