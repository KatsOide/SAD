      subroutine trotg(geo,omega,cost,sint)
      implicit none
      integer*4 i
      real*8 geo(3,3),omega(3),cost,sint,
     $     o1,o2,o3,gomega,gt1,gt2,gt3
      o1=omega(1)
      o2=omega(2)
      o3=omega(3)
      do 10 i=1,3
        gomega=geo(1,i)*o1+geo(2,i)*o2+geo(3,i)*o3
        gt1=geo(3,i)*o2-geo(2,i)*o3
        gt2=geo(1,i)*o3-geo(3,i)*o1
        gt3=geo(2,i)*o1-geo(1,i)*o2
        geo(1,i)=cost*(geo(1,i)-gomega*o1)+sint*gt1+gomega*o1
        geo(2,i)=cost*(geo(2,i)-gomega*o2)+sint*gt2+gomega*o2
        geo(3,i)=cost*(geo(3,i)-gomega*o3)+sint*gt3+gomega*o3
10    continue
      return
      end
