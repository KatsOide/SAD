      subroutine tmap(x1,x2,idir)
      use temw, only: rn=>r, ri
      implicit none
      real*8 ,intent(in):: x1(6)
      real*8 ,intent(out):: x2(6)
      integer*4 idir
c      real*8 x(6)

      if(idir >= 0)then
        x2=matmul(ri,x1)
      else
        x2=matmul(rn,x1)
      endif
      return
      end
