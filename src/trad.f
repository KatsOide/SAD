      module tspin
      use macphys

      type spin
      sequence
      real*8 sx,sy,sz
      end type

      contains
        subroutine tradkf1(x,px,y,py,z,g,dv,sp,px0,py0,bsi,al,dldx)
        use tfstk, only:pxy2dpz,p2h
        use ffs_flag
        use tmacro
        implicit none
        real*8 x,px,y,py,z,g,dv,px0,py0,bsi,al,dldx,
     $       dpx,dpy,pz,pz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,h2,tdusr
        type (spin) sp
        real*8, parameter:: gmin=-0.9999d0,
     $       cave=8.d0/15.d0/sqrt(3.d0),thetamax=0.005d0
        dpx=px-px0
        dpy=py-py0
c     theta=abs(dcmplx(dpx,dpy))
c     if(theta .gt. thetamax)then
        pz=1.d0+pxy2dpz(px,py)
        pz0=1.d0+pxy2dpz(px0,py0)
        ppx=py0*pz-pz0*py
        ppy=pz*px-px0*pz
        ppz=px0*py-py0*px
        theta=sqrt(ppx**2+ppy**2+ppz**2)
c     endif
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        anp=anrad*p*theta
        dg=tdusr(anp)
        if(dg .ne. 0.d0)then
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          al1=al*(1.d0+(dldx*pxm**2+pym**2)*.5d0)
          uc=cuc*(1.d0+p**2)*theta/al1*pr
          dg=-dg*uc
          g=max(gmin,g+dg)
          ddpx=-.5d0*dpx*dg
          ddpy=-.5d0*dpy*dg
          x=x+.5d0*ddpx*al
          y=y+.5d0*ddpy*al
          px=px+ddpx
          py=py+ddpy
          pr=1.d0+g
          h2=p2h(p0*pr)
          dv=-g*(1.d0+pr)/h2/(h2+pr*h0)+dvfs
          z=z*p0*pr/h2*h1/p
          if(calpol)then
            call sprot(sp,pxm,pym,ppx,ppy,ppz,bsi,h2)
          endif
        elseif(calpol)then
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sp,px0,py0,ppx,ppy,ppz,bsi,h1)
        endif
        return
        end subroutine

        subroutine tradk1(x,px,y,py,z,g,dv,sp,px0,py0,bsi,al,dldx)
        use tfstk, only:pxy2dpz,p2h
        use ffs_flag
        use tmacro
        implicit none
        real*8 x,px,y,py,z,g,dv,px0,py0,bsi,al,dldx,
     $       pz,pz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,
     $       pxm,pym,al1,uc,ddpx,ddpy,h2,h1
        type (spin) sp
        real*8, parameter:: gmin=-0.9999d0,
     $       cave=8.d0/15.d0/sqrt(3.d0),thetamax=0.005d0
        dpx=px-px0
        dpy=py-py0
        pz=1.d0+pxy2dpz(px,py)
        pz0=1.d0+pxy2dpz(px0,py0)
        ppx=py0*pz-pz0*py
        ppy=pz*px-px0*pz
        ppz=px0*py-py0*px
        theta=sqrt(ppx**2+ppy**2+ppz**2)
        pxm=px0+dpx*.5d0
        pym=py0+dpy*.5d0
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        al1=al*(1.d0+(dldx*pxm**2+pym**2)*.5d0)
        anp=anrad*p*theta
        uc=cuc*(1.d0+p**2)*theta/al1*pr
        dg=-cave*anp*uc
        g=max(gmin,g+dg)
        ddpx=-.5d0*dpx*dg
        ddpy=-.5d0*dpy*dg
        x=x+.5d0*ddpx*al
        y=y+.5d0*ddpy*al
        px=px+ddpx
        py=py+ddpy
        pr=1.d0+g
        h2=p2h(p0*pr)
        dv=-g*(1.d0+pr)/h2/(h2+pr*h0)+dvfs
        z=z*p0*pr/h2*h1/p
        if(calpol)then
          call sprot(sp,pxm,pym,ppx,ppy,ppz,bsi,h2)
        endif
        return
        end subroutine

        subroutine sprot(sp,pxm,pym,ppx,ppy,ppz,bsi,h)
        use tfstk,only:pxy2dpz
        use tmacro
        implicit none
        type (spin) sp
        real*8 pxm,pym,ppx,ppy,ppz,bsi,h,pzm,
     $       bx,by,bz,bp,blx,bly,blz,btx,bty,btz,ct,
     $       gx,gy,gz,g,
     $       gnx,gny,gnz,
     $       ux,uy,uz,u,
     $       unx,uny,unz,
     $       vnx,vny,vnz,
     $       su,sv,sw,su1,sv1,cosu,sinu
        real*8 , parameter :: cl=1.d0+gspin
        pzm=1.d0+pxy2dpz(pxm,pym)
        bx=pym*ppz-pzm*ppy
        by=pzm*ppx-pxm*ppz
        bz=pxm*ppy-pym*ppx+bsi
        bp=bx*pxm+by*pym+bz*pzm
        blx=bp*pxm
        bly=bp*pym
        blz=bp*pzm
        btx=bx-blx
        bty=by-bly
        btz=bz-blz
        ct=1.d0+h*gspin
        gx=ct*btx+cl*blx
        gy=ct*bty+cl*bly
        gz=ct*btz+cl*blz
        g=abs(dcmplx(gx,abs(dcmplx(gy,gz))))
        if(g .eq. 0.d0)then
          return
        endif
        gnx=gx/g
        gny=gy/g
        gnz=gz/g
        ux=sp%sy*gnz-sp%sz*gny
        uy=sp%sz*gnx-sp%sx*gnz
        uz=sp%sx*gny-sp%sy*gnx
        u=abs(dcmplx(ux,abs(dcmplx(uy,uz))))
        if(u .eq. 0.d0)then
          return
        endif
        unx=ux/u
        uny=uy/u
        unz=uz/u
        vnx=gny*unz-gnz*uny
        vny=gnz*unx-gnx*unz
        vnz=gnx*uny-gny*unx
        su=sp%sx*unx+sp%sy*uny+sp%sz*unz
        sv=sp%sx*vnx+sp%sy*vny+sp%sz*vnz
        sw=sp%sx*gnx+sp%sy*gny+sp%sz*gnz
        sinu=u*g
        cosu=sqrt(1.d0-sinu**2)
        su1=cosu*su-sinu*sv
        sv1=sinu*su+cosu*sv
        sp%sx=su1*unx+sv1*vnx+sw*gnx
        sp%sy=su1*uny+sv1*vny+sw*gny
        sp%sz=su1*unz+sv1*vnz+sw*gnz
        return
        end subroutine

        subroutine texspin(np,sp)
        implicit none
        integer*4 np,i
        real*8 sx
        type (spin) sp(np)
        do i=1,np
          sx=sp(i)%sx
          sp(i)%sx=sp(i)%sy
          sp(i)%sy=sx
        enddo
        return
        end subroutine

      end module

      subroutine trad(np,x,px,y,py,g,dv,by,bx,dbydx,dldx,dldxe,al,
     $     f1,f2,als,ala,dir)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 dbydx,dldx,dldxe,al,dir,brad,bx,by,tdusr,dg,
     $     pr,p,hh,dp,h1,alc,ald,dldxx,dldpx,al1,alr1,
     $     als,ala,alr2,alr3,f1,f2,uc,rhoinv0,anp
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
          al1=(1.d0+x(i)*dldxx+px(i)*dldpx)*alr1
     1         *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          rhoinv0=sqrt(brad)/brhoz
          anp=anrad*p0*rhoinv0*al1
          uc=cuc*(1.d0+p**2)*rhoinv0
c          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
c     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
c          de=er*sqrt(hh*brad)/p*dp*hh
c          g(i)=max(gmin,g(i)+dp+sqrt(abs(de))*dv(i))
c          if(i .eq. 1)then
c            write(*,'(a,1p8g15.7)')'trad ',anp,uc,1.d0/rhoinv,
c     $           al1*rhoinv,dg
c          endif
          dg=-uc*tdusr(anp)
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

      subroutine tradk(np,x,px,y,py,z,g,dv,sp,px0,py0,bsi,al)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),dv(np),z(np),g(np),
     $     px0(np),py0(np),bsi(np),al
      type (spin) sp(np)
      if(rfluct .and. al .gt. 0.d0)then
        do i=1,np
          call tradkf1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),sp(i),
     $         px0(i),py0(i),bsi(i),al,1.d0)
        enddo
      else
        do i=1,np
          call tradk1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),sp(i),
     $         px0(i),py0(i),bsi(i),al,1.d0)
        enddo
      endif
      return
      end

      subroutine tradki(np,x,px,y,py,z,g,dv,sp,px0,py0,bsi,al)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),dv(np),z(np),g(np),
     $     px0(np),py0(np),bsi(np),al(np)
      type (spin) sp(np)
      if(rfluct)then
        do i=1,np
          call tradkf1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),sp(i),
     $         px0(i),py0(i),bsi(i),al(i),0.d0)
        enddo
      else
        do i=1,np
          call tradk1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),sp(i),
     $         px0(i),py0(i),bsi(i),al(i),0.d0)
        enddo
      endif
      return
      end
