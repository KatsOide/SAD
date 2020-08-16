      subroutine tconvm(np,px,py,g,dv,idir)
      use tmacro
      use mathfun
      implicit none
      integer*4 ,intent(in):: np,idir
      real*8 ,intent(inout):: px(np),py(np),dv(np),g(np)
      real*8 pr(np),h1(np),p1(np)
      if(idir .ge. 0)then
        px=px+g*px
        py=py+g*py
      else
        pr=1.d0+g
        px=px/pr
        py=py/pr
        p1=p0*pr
        h1=p2h(p1)
        dv=-g*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
      endif
      return
      end
