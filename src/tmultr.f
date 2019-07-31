      subroutine tmultr(trans,trans1,n)
      implicit none
      integer*4 n
      real*8 trans(6,n),trans1(6,6),
     $     v1(n),v2(n),v3(n),v4(n),v5(n),v6(n)
c      do 10 i=1,n
        v1=trans(1,:)
        v2=trans(2,:)
        v3=trans(3,:)
        v4=trans(4,:)
        v5=trans(5,:)
        v6=trans(6,:)
        trans(1,:)=trans1(1,1)*v1+trans1(1,2)*v2+trans1(1,3)*v3
     1              +trans1(1,4)*v4+trans1(1,5)*v5+trans1(1,6)*v6
        trans(2,:)=trans1(2,1)*v1+trans1(2,2)*v2+trans1(2,3)*v3
     1              +trans1(2,4)*v4+trans1(2,5)*v5+trans1(2,6)*v6
        trans(3,:)=trans1(3,1)*v1+trans1(3,2)*v2+trans1(3,3)*v3
     1              +trans1(3,4)*v4+trans1(3,5)*v5+trans1(3,6)*v6
        trans(4,:)=trans1(4,1)*v1+trans1(4,2)*v2+trans1(4,3)*v3
     1              +trans1(4,4)*v4+trans1(4,5)*v5+trans1(4,6)*v6
        trans(5,:)=trans1(5,1)*v1+trans1(5,2)*v2+trans1(5,3)*v3
     1              +trans1(5,4)*v4+trans1(5,5)*v5+trans1(5,6)*v6
        trans(6,:)=trans1(6,1)*v1+trans1(6,2)*v2+trans1(6,3)*v3
     1              +trans1(6,4)*v4+trans1(6,5)*v5+trans1(6,6)*v6
c10    continue
      return
      end

      subroutine tmultr5(trans,trans1,n)
      implicit none
      integer*4 n
      real*8 trans(6,n),trans1(6,6),
     $     v1(n),v2(n),v3(n),v4(n),v6(n)
c      if(trans1(1,5) .ne. 0.d0
c     $     .or. trans1(2,5) .ne. 0.d0
c     $     .or. trans1(3,5) .ne. 0.d0
c     $     .or. trans1(4,5) .ne. 0.d0
c     $     .or. trans1(5,5) .ne. 1.d0)then
c       write(*,*)trans(6,1000000)
c      endif
          v1=trans(1,:)
          v2=trans(2,:)
          v3=trans(3,:)
          v4=trans(4,:)
          v6=trans(6,:)
          trans(1,:)=trans1(1,1)*v1+trans1(1,2)*v2+trans1(1,3)*v3
     1         +trans1(1,4)*v4+trans1(1,6)*v6
          trans(2,:)=trans1(2,1)*v1+trans1(2,2)*v2+trans1(2,3)*v3
     1         +trans1(2,4)*v4+trans1(2,6)*v6
          trans(3,:)=trans1(3,1)*v1+trans1(3,2)*v2+trans1(3,3)*v3
     1         +trans1(3,4)*v4+trans1(3,6)*v6
          trans(4,:)=trans1(4,1)*v1+trans1(4,2)*v2+trans1(4,3)*v3
     1         +trans1(4,4)*v4+trans1(4,6)*v6
          trans(5,:)=trans1(5,1)*v1+trans1(5,2)*v2+trans1(5,3)*v3
     1         +trans1(5,4)*v4+trans(5,:)+trans1(5,6)*v6
      return
      end
