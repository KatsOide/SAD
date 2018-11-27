      subroutine trad(np,x,px,y,py,g,dv,by,bx,dbydx,dldx,dldxe,al,
     $     f1,f2,als,ala,dir)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 dbydx,dldx,dldxe,al,dir,brad,bx,by,tdusr,dg,
     $     pr,p,hh,dp,h1,alc,ald,dldxx,dldpx,al1,alr1,
     $     als,ala,alr2,alr3,f1,f2,cuc,uc,rhoinv0,anp
      real*8, parameter:: gmin=-0.9999d0
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np)
      call tradel(al,f1,f2,als,ala,alr1,alr2,alr3)
c      write(*,'(a,1p7g14.6)')'trad ',al,f1,f2,als,ala,alr,alr1
      ald=dir*al
      dldxx=dldx+dldxe
c      dldpx=dldx*ald
      dldpx=dldx*ald*.5d0
      if(rfluct .and. al .gt. 0.d0)then
        cuc=1.5d0*rclassic/rcratio
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

      subroutine tradk(np,x,px,y,py,px0,py0,g,dv,al)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 al,tdusr,pr,p,al1,cuc,uc,anp,
     $     dpx,dpy,pxm,pym,theta,dg,h1,ddpx,ddpy,
     $     pz,pz0,ppx,ppy,ppz
      real*8, parameter:: gmin=-0.9999d0,
     $     cave=8.d0/15.d0/sqrt(3.d0),thetamax=0.005d0
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np),
     $     px0(np),py0(np)
      cuc=1.5d0*rclassic/rcratio
      if(rfluct .and. al .gt. 0.d0)then
        do i=1,np
          dpx=px(i)-px0(i)
          dpy=py(i)-py0(i)
          theta=abs(dcmplx(dpx,dpy))
          if(theta .gt. thetamax)then
            pz=1.d0+pxy2dpz(px(i),py(i))
            pz0=1.d0+pxy2dpz(px0(i),py0(i))
            ppx=py0(i)*pz-pz0*py(i)
            ppy=pz*px(i)-px0(i)*pz
            ppz=px0(i)*py(i)-py0(i)*px(i)
            theta=sqrt(ppx**2+ppy**2+ppz**2)
          endif
          pr=1.d0+g(i)
          p=p0*pr
          anp=anrad*p*theta
          dg=tdusr(anp)
          if(dg .ne. 0.d0)then
            pxm=px0(i)+dpx*.5d0
            pym=py0(i)+dpy*.5d0
            al1=al*(1.d0+(pxm**2+pym**2)*.5d0)
            uc=cuc*(1.d0+p**2)*theta/al1*pr
            dg=-dg*uc
c            write(*,*)'tradk ',i,dg,rhoinv,anp,uc
            g(i)=max(gmin,g(i)+dg)
            ddpx=-.5d0*dpx*dg
            ddpy=-.5d0*dpy*dg
            x(i)=x(i)+.5d0*ddpx*al
            y(i)=y(i)+.5d0*ddpy*al
            px(i)=px(i)+ddpx
            py(i)=py(i)+ddpy
            pr=1.d0+g(i)
            h1=p2h(p0*pr)
            dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
          endif
        enddo
      else
        do i=1,np
          dpx=px(i)-px0(i)
          dpy=py(i)-py0(i)
          theta=abs(dcmplx(dpx,dpy))
          if(theta .gt. thetamax)then
            pz=1.d0+pxy2dpz(px(i),py(i))
            pz0=1.d0+pxy2dpz(px0(i),py0(i))
            ppx=py0(i)*pz-pz0*py(i)
            ppy=pz*px(i)-px0(i)*pz
            ppz=px0(i)*py(i)-py0(i)*px(i)
            theta=sqrt(ppx**2+ppy**2+ppz**2)
          endif
          pxm=px0(i)+dpx*.5d0
          pym=py0(i)+dpy*.5d0
          pr=1.d0+g(i)
          p=p0*pr
          al1=al*(1.d0+(pxm**2+pym**2)*.5d0)
          anp=anrad*p*theta
          uc=cuc*(1.d0+p**2)*theta/al1*pr
          dg=-cave*anp*uc
          g(i)=max(gmin,g(i)+dg)
          ddpx=-.5d0*dpx*dg
          ddpy=-.5d0*dpy*dg
          x(i)=x(i)+.5d0*ddpx*al
          y(i)=y(i)+.5d0*ddpy*al
          px(i)=px(i)+ddpx
          py(i)=py(i)+ddpy
          pr=1.d0+g(i)
          h1=p2h(p0*pr)
          dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
        enddo
      endif
      return
      end

      subroutine tradki(np,x,px,y,py,px0,py0,g,dv,al)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 tdusr,pr,p,al1,cuc,uc,anp,
     $     dpx,dpy,pym,theta,dg,h1,ddpx,ddpy,
     $     pz,pz0,ppx,ppy,ppz
      real*8, parameter:: gmin=-0.9999d0,
     $     cave=8.d0/15.d0/sqrt(3.d0),thetamax=0.005d0
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np),
     $     px0(np),py0(np),al(np)
      cuc=1.5d0*rclassic/rcratio
      if(rfluct)then
        do i=1,np
          dpx=px(i)-px0(i)
          dpy=py(i)-py0(i)
          theta=abs(dcmplx(dpx,dpy))
          if(theta .gt. thetamax)then
            pz=1.d0+pxy2dpz(px(i),py(i))
            pz0=1.d0+pxy2dpz(px0(i),py0(i))
            ppx=py0(i)*pz-pz0*py(i)
            ppy=pz*px(i)-px0(i)*pz
            ppz=px0(i)*py(i)-py0(i)*px(i)
            theta=sqrt(ppx**2+ppy**2+ppz**2)
          endif
          pr=1.d0+g(i)
          p=p0*pr
          anp=anrad*p*theta
          dg=tdusr(anp)
          if(dg .ne. 0.d0)then
            pym=py0(i)+dpy*.5d0
            al1=al(i)*(1.d0+pym**2*.5d0)
            uc=cuc*(1.d0+p**2)*theta/al1*pr
            dg=-dg*uc
c            write(*,*)'tradk ',i,dg,rhoinv,anp,uc
            g(i)=max(gmin,g(i)+dg)
            ddpx=-.5d0*dpx*dg
            ddpy=-.5d0*dpy*dg
            x(i)=x(i)+.5d0*ddpx*al(i)
            y(i)=y(i)+.5d0*ddpy*al(i)
            px(i)=px(i)+ddpx
            py(i)=py(i)+ddpy
            pr=1.d0+g(i)
            h1=p2h(p0*pr)
            dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
          endif
        enddo
      else
        do i=1,np
          dpx=px(i)-px0(i)
          dpy=py(i)-py0(i)
          theta=abs(dcmplx(dpx,dpy))
          if(theta .gt. thetamax)then
            pz=1.d0+pxy2dpz(px(i),py(i))
            pz0=1.d0+pxy2dpz(px0(i),py0(i))
            ppx=py0(i)*pz-pz0*py(i)
            ppy=pz*px(i)-px0(i)*pz
            ppz=px0(i)*py(i)-py0(i)*px(i)
            theta=sqrt(ppx**2+ppy**2+ppz**2)
          endif
          pym=py0(i)+dpy*.5d0
          pr=1.d0+g(i)
          p=p0*pr
          al1=al(i)*(1.d0+pym**2*.5d0)
          anp=anrad*p*theta
          uc=cuc*(1.d0+p**2)*theta/al1*pr
          dg=-cave*anp*uc
          g(i)=max(gmin,g(i)+dg)
          ddpx=-.5d0*dpx*dg
          ddpy=-.5d0*dpy*dg
          x(i)=x(i)+.5d0*ddpx*al(i)
          y(i)=y(i)+.5d0*ddpy*al(i)
          px(i)=px(i)+ddpx
          py(i)=py(i)+ddpy
          pr=1.d0+g(i)
          h1=p2h(p0*pr)
          dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
        enddo
      endif
      return
      end
