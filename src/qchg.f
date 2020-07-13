      subroutine qchg(trans,cod,dx,dy,theta,enter)
      implicit none
      real*8 ,intent(inout):: trans(4,5),cod(6)
      real*8 ,intent(in):: dx,dy,theta
      real*8 cost,sint,x(5),y(5),xi,pxi
      logical*4 ,intent(in):: enter
      if(enter)then
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
      endif
      if(theta .ne. 0.d0)then
        cost=cos(theta)
        sint=sin(theta)
c        do 10 i=1,5
          x=trans(1,:)
          trans(1,:)= cost*x-sint*trans(3,:)
          trans(3,:)= sint*x+cost*trans(3,:)
          y=trans(2,:)
          trans(2,:)= cost*y-sint*trans(4,:)
          trans(4,:)= sint*y+cost*trans(4,:)
c10      continue
        xi=cod(1)
        cod(1)= cost*xi-sint*cod(3)
        cod(3)= sint*xi+cost*cod(3)
        pxi=cod(2)
        cod(2)= cost*pxi-sint*cod(4)
        cod(4)= sint*pxi+cost*cod(4)
      endif
      if( .not. enter)then
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
      endif
      return
      end
