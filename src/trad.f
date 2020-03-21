c  Obsolete 22 Jan 2020
c
      subroutine trad(np,x,px,y,py,g,dv,by,bx,dbydx,dldx,dldxe,al,
     $     f1,f2,als,ala,dir)
      use ffs_flag
      use tmacro
      use mathfun
      implicit none
      integer*4 np,i
      real*8 dbydx,dldx,dldxe,al,dir,brad,bx,by,tdusr,dg,
     $     pr,p,hh,dp,h1,alc,ald,dldxx,dldpx,al1,alr1,
     $     als,ala,alr2,alr3,f1,f2,uc,rhoinv0,anp,an,h
      real*8, parameter:: gmin=-0.9999d0
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np)
      call tradel(al,f1,f2,als,ala,alr1,alr2,alr3)
c      write(*,'(a,1p7g14.6)')'trad ',al,f1,f2,als,ala,alr,alr1
      ald=dir*al
      dldxx=dldx+dldxe
c      dldpx=dldx*ald
      dldpx=dldx*ald*.5d0
      if(rfluct .and. al .gt. 0.d0)then
c        er=c/amass*erad*alr1/alr
c        call tran_array(dv(1),np)
c        do i=1,np
c          dv(i)=(dv(i)-.5d0)*3.46410161513775461d0
c        enddo
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          pr=1.d0+g(i)
          p=p0*pr
          h=p2h(p)
          al1=(1.d0+x(i)*dldxx+px(i)*dldpx)*alr1
     1         *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          rhoinv0=sqrt(brad)/brhoz
          anp=anrad*p0*rhoinv0*al1
          uc=cuc*h**3/p0*rhoinv0
c          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
c     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
c          de=er*sqrt(hh*brad)/p*dp*hh
c          g(i)=max(gmin,g(i)+dp+sqrt(abs(de))*dv(i))
c          if(i .eq. 1)then
c            write(*,'(a,1p8g15.7)')'trad ',anp,uc,1.d0/rhoinv,
c     $           al1*rhoinv,dg
c          endif
          dg=-uc*tdusr(anp,an)
          if(dg .ne. 0.d0)then
            g(i)=max(gmin,g(i)+dg)
            pr=1.d0+g(i)
            h1=p2h(p0*pr)
            dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
          endif
        enddo
      else
        alc=alr2*crad
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          hh=1.d0+(p0*(1.d0+g(i)))**2
          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          g(i)=max(gmin,g(i)+dp)
          pr=1.d0+g(i)
          h1=p2h(p0*pr)
          dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
        enddo
      endif
      return
      end

      subroutine tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phi0)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      implicit none
      integer*4 np,i
      real*8 ,intent(inout)::
     $     x(np),px(np),y(np),py(np),dv(np),z(np),g(np),
     $     sx(np),sy(np),sz(np)
      real*8 , intent(in)::al,phi0
      cphi0=cos(phi0)
      sphi0=sin(phi0)
      if(rfluct .and. al .ne. 0.d0)then
        do i=1,np
          call tradkf1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),
     $         sx(i),sy(i),sz(i),
     $         pxr0(i),pyr0(i),zr0(i),bsi(i),al,i)
        enddo
      else
        call tradkn(np,x,px,y,py,z,g,dv,
     $       sx,sy,sz,al)
c        do i=1,np
c          call tradk1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),
c     $         sx(i),sy(i),sz(i),
c     $         px0(i),py0(i),zr0(i),bsi(i),al)
c        enddo
      endif
      pxr0=px
      pyr0=py
      zr0=z
      bsi=0.d0
      return
      end
