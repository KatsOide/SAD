      subroutine tinse(trans,cod,beam,trx)
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      implicit none
      real*8 trans(6,12),cod(6),beam(42),trx(6,7),trans1(6,7),
     $     x1,px1,y1,py1
      trans1(:,1:6)=trx(:,1:6)
      call tmultr(trans,trans1(:,1:6),irad)
      call tmulbs(beam ,trans1(:,1:6),.true.)
       x1=trans1(1,1)*cod(1)+trans1(1,2)*cod(2)+
     $    trans1(1,3)*cod(3)+trans1(1,4)*cod(4)+trans1(1,6)*cod(6)+
     $    trans1(1,7)
      px1=trans1(2,1)*cod(1)+trans1(2,2)*cod(2)+
     $    trans1(2,3)*cod(3)+trans1(2,4)*cod(4)+trans1(2,6)*cod(6)+
     $    trans1(2,7)
       y1=trans1(3,1)*cod(1)+trans1(3,2)*cod(2)+
     $    trans1(3,3)*cod(3)+trans1(3,4)*cod(4)+trans1(3,6)*cod(6)+
     $    trans1(3,7)
      py1=trans1(4,1)*cod(1)+trans1(4,2)*cod(2)+
     $    trans1(4,3)*cod(3)+trans1(4,4)*cod(4)+trans1(4,6)*cod(6)+
     $    trans1(4,7)
      cod(5)=trans1(5,1)*cod(1)+trans1(5,2)*cod(2)+
     $       trans1(5,3)*cod(3)+trans1(5,4)*cod(4)+cod(5)
      cod(1)=x1
      cod(2)=px1
      cod(3)=y1
      cod(4)=py1
      return
      end
