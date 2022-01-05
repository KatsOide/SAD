      subroutine gaus3(x,p,y,dydp,mp)
      implicit none
      real*8 x,p,y,dydp
      integer mp
      dimension p(mp),dydp(mp)

      dydp(1)=exp(-0.5*(x-p(3))**2/p(2)**2)
      y=p(1)*dydp(1)
      dydp(2)=p(1)*(x-p(3))**2/p(2)**3 *dydp(1)
      dydp(3)=p(1)*(x-p(3))/p(2)**2 *dydp(1)
      return
      end

