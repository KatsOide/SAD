      subroutine ttcav(np,x,px,y,py,z,g,dv,sx,sy,sz,
     1                 al,ak,harm,phi,freq,dx,dy,theta,krad)
      use tmacro
      use ffs_flag, only:calpol
      use mathfun
      use kradlib
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np)
      real*8 al,ak,harm,phi,freq,dx,dy,theta,cost,sint,w,phic,
     $     dp1r,p1r,p1,h1,t,ph,dh,h2,a,dpr,dp2r,p2r,xi,pxi,v
      logical*4 krad
      include 'inc/TENT.inc'
      if(al .ne. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,dv,al*.5d0)
        if(krad)then
          pxr0=px
          pyr0=py
          zr0=z
        endif
      endif
      if(harm .eq. 0.d0)then
        w=pi2*freq/c
      else
        w=omega0*harm/c
      endif
      phic=phi*sign(1.d0,charge)
      v=ak*p0*w
      do 10 i=1,np
c        dp1r=g(i)*(2.d0+g(i))
        dp1r=g(i)
        p1r=1.d0+dp1r
        p1=p0*p1r
        if(dv(i) .gt. 0.1d0)then
          h1=p2h(p1)
c          h1=sqrt(1.d0+p1**2)
        else
          h1=p1r*h0/(1.d0-dv(i)+dvfs)
        endif
        t=min(1.d4,max(-1.d4,-z(i)*h1/p1))
        ph=w*t+phic
        dh=v*x(i)*cos(ph)
        px(i)=px(i)-ak*sin(ph)/p1r
        h2=h1+dh
        a=max(-1.d0,dh*(h1+h2)/p1**2)
        dpr=a*(.5d0-a*(.125d0-a*.0625d0))
        dpr=(dpr**2+a)/(2.d0+2.d0*dpr)
        dpr=(dpr**2+a)/(2.d0+2.d0*dpr)
        dp2r=dp1r+p1r*dpr
        p2r=1.d0+dp2r
c        gg=dp2r*(.5d0-dp2r*(.125d0-dp2r*.0625d0))
c        gg=(gg**2+dp2r)/(2.d0+2.d0*gg)
c        g(i)=(gg**2+dp2r)/(2.d0+2.d0*gg)
        g(i)=dp2r
        dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs
        px(i)=px(i)*p1r/p2r
        py(i)=py(i)*p1r/p2r
        z(i)=-t*p2r*p0/h2
10    continue
      if(al .ne. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,dv,al*.5d0)
        if(krad)then
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al,0.d0)
        endif
      endif
      include 'inc/TEXIT.inc'
      return
      end
