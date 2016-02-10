      subroutine tmuld(trans,trans1)
      implicit none
      real*8 trans(6,12),trans1(6,6),v1,v2,v3,v4,v5,v6
      integer*4 i
      do 10 i=7,12
        v1=trans(1,i)+trans(1,i-6)
        v2=trans(2,i)+trans(2,i-6)
        v3=trans(3,i)+trans(3,i-6)
        v4=trans(4,i)+trans(4,i-6)
        v5=trans(5,i)+trans(5,i-6)
        v6=trans(6,i)+trans(6,i-6)
        trans(2,i)=trans(2,i)+trans1(2,1)*v1
        trans(2,i)=trans(2,i)+trans1(2,2)*v2
        trans(2,i)=trans(2,i)+trans1(2,3)*v3
        trans(2,i)=trans(2,i)+trans1(2,4)*v4
        trans(2,i)=trans(2,i)+trans1(2,5)*v5
        trans(2,i)=trans(2,i)+trans1(2,6)*v6
        trans(4,i)=trans(4,i)+trans1(4,1)*v1
        trans(4,i)=trans(4,i)+trans1(4,2)*v2
        trans(4,i)=trans(4,i)+trans1(4,3)*v3
        trans(4,i)=trans(4,i)+trans1(4,4)*v4
        trans(4,i)=trans(4,i)+trans1(4,5)*v5
        trans(4,i)=trans(4,i)+trans1(4,6)*v6
        trans(6,i)=trans(6,i)+trans1(6,1)*v1
        trans(6,i)=trans(6,i)+trans1(6,2)*v2
        trans(6,i)=trans(6,i)+trans1(6,3)*v3
        trans(6,i)=trans(6,i)+trans1(6,4)*v4
        trans(6,i)=trans(6,i)+trans1(6,5)*v5
        trans(6,i)=trans(6,i)+trans1(6,6)*v6
10    continue
      return
      end

      subroutine tmuld6(trans,trans1)
      implicit none
      real*8 trans(6,12),trans1(6,6),v1,v2,v3,v4,v5,v6
      integer*4 i
      do 10 i=7,12
        v1=trans(1,i)+trans(1,i-6)
        v2=trans(2,i)+trans(2,i-6)
        v3=trans(3,i)+trans(3,i-6)
        v4=trans(4,i)+trans(4,i-6)
        v5=trans(5,i)+trans(5,i-6)
        v6=trans(6,i)+trans(6,i-6)
        trans(1,i)=trans(1,i)+trans1(1,1)*v1
        trans(1,i)=trans(1,i)+trans1(1,2)*v2
        trans(1,i)=trans(1,i)+trans1(1,3)*v3
        trans(1,i)=trans(1,i)+trans1(1,4)*v4
        trans(1,i)=trans(1,i)+trans1(1,5)*v5
        trans(1,i)=trans(1,i)+trans1(1,6)*v6
        trans(2,i)=trans(2,i)+trans1(2,1)*v1
        trans(2,i)=trans(2,i)+trans1(2,2)*v2
        trans(2,i)=trans(2,i)+trans1(2,3)*v3
        trans(2,i)=trans(2,i)+trans1(2,4)*v4
        trans(2,i)=trans(2,i)+trans1(2,5)*v5
        trans(2,i)=trans(2,i)+trans1(2,6)*v6
        trans(3,i)=trans(3,i)+trans1(3,1)*v1
        trans(3,i)=trans(3,i)+trans1(3,2)*v2
        trans(3,i)=trans(3,i)+trans1(3,3)*v3
        trans(3,i)=trans(3,i)+trans1(3,4)*v4
        trans(3,i)=trans(3,i)+trans1(3,5)*v5
        trans(3,i)=trans(3,i)+trans1(3,6)*v6
        trans(4,i)=trans(4,i)+trans1(4,1)*v1
        trans(4,i)=trans(4,i)+trans1(4,2)*v2
        trans(4,i)=trans(4,i)+trans1(4,3)*v3
        trans(4,i)=trans(4,i)+trans1(4,4)*v4
        trans(4,i)=trans(4,i)+trans1(4,5)*v5
        trans(4,i)=trans(4,i)+trans1(4,6)*v6
        trans(5,i)=trans(5,i)+trans1(5,1)*v1
        trans(5,i)=trans(5,i)+trans1(5,2)*v2
        trans(5,i)=trans(5,i)+trans1(5,3)*v3
        trans(5,i)=trans(5,i)+trans1(5,4)*v4
        trans(5,i)=trans(5,i)+trans1(5,5)*v5
        trans(5,i)=trans(5,i)+trans1(5,6)*v6
        trans(6,i)=trans(6,i)+trans1(6,1)*v1
        trans(6,i)=trans(6,i)+trans1(6,2)*v2
        trans(6,i)=trans(6,i)+trans1(6,3)*v3
        trans(6,i)=trans(6,i)+trans1(6,4)*v4
        trans(6,i)=trans(6,i)+trans1(6,5)*v5
        trans(6,i)=trans(6,i)+trans1(6,6)*v6
10    continue
      return
      end
