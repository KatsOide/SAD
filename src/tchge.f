      subroutine tchge(trans,cod,beam,dx,dy,theta,enter,ld)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 ld
      real*8 trans(6,12),cod(6),beam(42),trans1(6,13),
     $     dx,dy,theta,cost,sint,xi,pxi
      logical enter
      if(enter)then
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
      endif
      if(theta .ne. 0.d0)then
        call tinitr(trans1)
        cost=cos(theta)
        sint=sin(theta)
        trans1(1,1)= cost
        trans1(1,3)=-sint
        trans1(3,1)= sint
        trans1(3,3)= cost
        trans1(2,2)= cost
        trans1(2,4)=-sint
        trans1(4,2)= sint
        trans1(4,4)= cost
        call tmultr(trans,trans1,irad)
        call tmulbs(beam,trans1,.true.,.true.)
        if(calpol)then
          call polpar(0,ld,theta,0.d0,0.d0,0.d0,0.d0,cod)
        endif
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
