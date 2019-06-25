      subroutine qchg(trans,cod,dx,dy,theta,enter)
      implicit none
      real*8 trans(4,5),cod(6),dx,dy,theta,
     $     cost,sint,x,y,xi,pxi
      integer*4 i
      logical*4 enter
      if(enter)then
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
      endif
      if(theta .ne. 0.d0)then
        cost=cos(theta)
        sint=sin(theta)
        do 10 i=1,5
          x=trans(1,i)
          trans(1,i)= cost*x-sint*trans(3,i)
          trans(3,i)= sint*x+cost*trans(3,i)
          y=trans(2,i)
          trans(2,i)= cost*y-sint*trans(4,i)
          trans(4,i)= sint*y+cost*trans(4,i)
10      continue
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
