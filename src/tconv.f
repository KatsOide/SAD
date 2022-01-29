      subroutine tconv(x1,x2,idir)
      use tmacro
      use mathfun
      implicit none
      integer*4 ,intent(in):: idir
      real*8 ,intent(in):: x1(8)
      real*8 ,intent(out):: x2(8)
      real*8 pr,h1,p1
      if(idir .ge. 0)then
c        x2(6)=x1(6)*(2.d0+x1(6))
        x2(6)=x1(6)
        pr=1.d0+x2(6)
        x2(1)=x1(1)
        x2(2)=x1(2)*pr
        x2(3)=x1(3)
        x2(4)=x1(4)*pr
        x2(5)=x1(5)
      else
        pr=1.d0+x1(6)
        x2(1)=x1(1)
        x2(2)=x1(2)/pr
        x2(3)=x1(3)
        x2(4)=x1(4)/pr
        x2(5)=x1(5)
        p1=p0*pr
        h1=p2h(p1)
c        h1=p1*sqrt(1.d0+1.d0/p1**2)
        x2(7)=-x1(6)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
        x2(6)=x1(6)
c        x2(6)=x1(6)/(1.d0+sqrt(pr))
c       a=min(.95d0,x2(2)**2+x2(4)**2)
c       x2(8)=-a/(1.d0+sqrt(1.d0-a))
      endif
      return
      end
