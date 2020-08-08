      subroutine tconvm(np,px,py,g,dv,idir)
      use tmacro
      use mathfun
      implicit none
      integer*4 ,intent(in):: np,idir
      integer*4 i
      real*8 ,intent(inout):: px(np),py(np),dv(np),g(np)
      real*8 pr,h1,p1
      if(idir .ge. 0)then
        px=px+g*px
        py=py+g*py
c        do 10 i=1,np
cc          g(i)=g(i)*(2.d0+g(i))
c          pr=1.d0+g(i)
c          px(i)=px(i)*pr
c          py(i)=py(i)*pr
c10      continue
      else
        do concurrent (i=1:np)
          pr=1.d0+g(i)
          px(i)=px(i)/pr
          py(i)=py(i)/pr
          p1=p0*pr
          h1=p2h(p1)
c          h1=p1*sqrt(1.d0+1.d0/p1**2)
          dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
c          g(i)=g(i)/(1.d0+sqrt(pr))
c         a=min(.95d0,px(i)**2+py(i)**2)
c         pz(i)=-a/(1.d0+sqrt(1.d0-a))
        enddo
      endif
      return
      end
