      subroutine tconvm(np,px,py,g,dv,idir)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 np,idir,i
      real*8 px(np),py(np),dv(np),g(np),pr,h1,p1
      if(idir .ge. 0)then
        do 10 i=1,np
c          g(i)=g(i)*(2.d0+g(i))
          pr=1.d0+g(i)
          px(i)=px(i)*pr
          py(i)=py(i)*pr
10      continue
      else
        do 20 i=1,np
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
20      continue
      endif
      return
      end
