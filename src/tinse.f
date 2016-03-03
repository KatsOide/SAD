      subroutine tinse(trans,cod,beam,trx,ld)
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 ld      
      real*8 trans(6,12),cod(6),beam(42),trx(6,7),trans1(6,13),
     $     x1,px1,y1,py1
      call tclr(trans1(1,7),36)
      call tmov(trx,trans1,36)
      call tmultr(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      if(calpol)then
        call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
      endif
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
