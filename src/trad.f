      subroutine trad(np,x,px,y,py,g,dv,by,bx,dbydx,dldx,dldxe,al,
     $     f1,f2,als,ala,dir)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 gmin,dbydx,dldx,dldxe,al,dir,brad,bx,by,
     $     pr,p,hh,dp,de,h1,er,alc,ald,dldxx,dldpx,
     $     als,ala,alr,alr1,f1,f2
      parameter (gmin=-0.99d0)
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np)
      call tradel(al,f1,f2,als,ala,alr,alr1)
c        write(*,*)al,f1,als,ala,alr,alr1
      alc=alr*crad
      ald=dir*al
      dldxx=dldx+dldxe
c      dldpx=dldx*ald
      dldpx=dldx*ald*.5d0
      if(rfluct)then
        er=c/amass*erad*alr1/alr
        call tran_array(dv(1),np)
        do i=1,np
          dv(i)=(dv(i)-.5d0)*3.46410161513775461d0
        enddo
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          pr=1.d0+g(i)
          p=p0*pr
          hh=1.d0+p**2
          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          de=er*sqrt(hh*brad)/p*dp*hh
          g(i)=max(gmin,g(i)+dp+sqrt(abs(de))*dv(i))
        enddo
      else
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          hh=1.d0+(p0*(1.d0+g(i)))**2
          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          g(i)=max(gmin,g(i)+dp)
        enddo
      endif
      do i=1,np
        pr=1.d0+g(i)
        h1=p2h(p0*pr)
c        h1=sqrt(1.d0+(p0*pr)**2)
        dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
c        g(i)=g(i)/(1.d0+sqrt(pr))
      enddo
      return
      end
