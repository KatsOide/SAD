cvs $Hearder$
      subroutine tsdrad(np,x,px,y,py,z,g,dv,al,rho)
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 np,i,n,ndiv
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     al,rho,b,uave,gi,pr,bzx,pxi,pyi,s,brad,rhor,al1,
     $     alx,prob,ep,dpzi,pzi,rhox,phi,sinphi,a12,a14,a22,a24,dp,
     $     h1,pxf,tran,drndsr,gi0
      b=(brhoz/rho)**2
      uave=8.d0/15.d0/sqrt(3.d0)
      do i=1,np
        gi=g(i)
        pr=1.d0+gi
        bzx=1.d0/rho/pr
        pxi=px(i)+bzx*y(i)*.5d0
        pyi=py(i)-bzx*x(i)*.5d0
        s=min(.999d0,pxi**2+pyi**2)
        brad=b*s
        if(brad .eq. 0.d0)then
          z(i)=z(i)-dv(i)*al
          cycle
        endif
        rhor=brhoz/max(1.d-20,sqrt(brad))
        al1=.1d0*rhor/anrad/p0
        ndiv=1+int(al/al1)
        alx=al/ndiv
        prob=.1d0*alx/al1
        ep=urad/rhor*p0**3
        dpzi=-s/(1.d0+sqrt(1.d0-s))
        pzi=1.d0+dpzi
        gi0=0.d0
        do n=1,ndiv
          if(trpt)then
            if(tran() .le. 1.d0-prob)then
              dp=0.d0
            else
              dp=-ep*pr**2*drndsr()
            endif
          else
            if(rfluct)then
              dp=-ep*pr**2*drndsr()*prob
            else
              dp=-ep*pr**2*uave*prob
            endif
          endif
          gi=gi+dp*.5d0
          pr=1.d0+gi
          if(gi .ne. gi0)then
            gi0=gi
            rhox=rho*pr
            phi=alx/rhox/pzi
            sinphi=sin(phi)
            a12=sinphi*rhox
            a14= 2.d0*sin(phi*.5d0)**2*rhox
            a22=cos(phi)
            a24= sinphi
            x(i)=x(i)+a12*pxi+a14*pyi
            pxf =     a22*pxi+a24*pyi
            y(i)=y(i)-a14*pxi+a12*pyi
            pyi =    -a24*pxi+a22*pyi
            pxi=pxf
          endif
          if(dp .ne. 0.d0)then
            gi=gi+dp*.5d0
          endif
        enddo
        g(i)=gi
        bzx=1.d0/rho/(1.d0+gi)
        px(i)=pxi-bzx*y(i)*.5d0
        py(i)=pyi+bzx*x(i)*.5d0
        z(i)=z(i)+(dpzi/pzi-dv(i))*al
      enddo
      do i=1,np
        pr=1.d0+g(i)
        h1=sqrt(1.d0+(p0*pr)**2)
        dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
      enddo
      return
      end
