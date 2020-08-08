      subroutine spkick(np,x,px,y,py,z,g,dv,al,
     $     radius,alx,kturn,kptbl)
      use tfstk
      use ffs
      use tffitcode
      use tmacro, only:l_track
      use iso_c_binding
      implicit none
      integer*4 nzmax,ntheta
      parameter (nzmax=1000,ntheta=12)
      integer*8 itt,irho,iphi,iq,iphis,iey,iez
      integer*4 ,intent(inout):: np,kptbl(np0,6)
      integer*4 kturn
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     $     sx(np0),sy(np0),sz(np0),zz(np0)
      real*8 ,pointer ::rho(:,:,:),phi(:,:,:),q(:,:),phis(:),
     $     ex(:,:,:),ey(:,:,:),ez(:,:,:)
      real*8 al,radius,alx
      integer*4 i,npz,nx,nz,nz1,nz2,nza,nrho,nphi,nq
      real*8 v0,gamma0,dx,dy,zc(nzmax),zf(nzmax)
      if(radius .eq. 0.d0)then
        return
      endif
c      call spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
      call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl,np,kturn,
     $     radius,radius,
     $     0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      if(radius .le. 0.d0)then
        return
      endif
      if(np .le. 1)then
        return
      endif
      call tftclupdate(7)
      call spzz(np,px,py,z,g,zz,v0,gamma0)
      itt=ktaloc((np+1)/2)
      do i=1,np
        ilist(i,itt)=i
      enddo
      call spsort(np,ilist(1:np,itt),zz)
      call spmesh(nx,nz,nz1,nz2,dx,dy,zc,zf,np,npz,x,y,zz,radius,
     $     ilist(1:np,itt),nzmax)
      nza=nz2-nz1+1
      nrho=(2*nx)**2*nza
      irho=ktaloc(nrho)
      call c_f_pointer(c_loc(rlist(irho)),rho,[2*nx,2*nx,nza])
      rho=0.d0
      call sprho(rho,nx,nz,nz1,nz2,dx,dy,zc,
     $     np,npz,x,y,zz,ilist(1:np,itt))
      nphi=(2*nx+2)**2*(nza+2)
      iphi=ktaloc(nphi)
      call c_f_pointer(c_loc(rlist(iphi)),phi,[2*nx+2,2*nx,nza+2])
      phi=0.d0
      nq=ntheta*nz
      iq=ktaloc(nq)
      call c_f_pointer(c_loc(rlist(iq)),q,[ntheta,nz])
      q=0.d0
      iphis=ktaloc(nq)
      call c_f_pointer(c_loc(rlist(iphis)),phis,[nq])
      phis=0.d0
      call spphi(rho,phi,q,phis,
     $     nx,nz,nz1,nz2,np,dx,dy,zc,zf,radius,ntheta)
      call c_f_pointer(c_loc(rlist(irho)),ex,[2*nx,2*nx,nza])
      iey=ktaloc(nrho)
      call c_f_pointer(c_loc(rlist(iey)),ey,[2*nx,2*nx,nza])
      iez=ktaloc(nrho)
      call c_f_pointer(c_loc(rlist(iez)),ez,[2*nx,2*nx,nza])
      call spfield(phi,ex,ey,ez,
     $     nx,nz,nz1,nz2,dx,dy,zc,gamma0)
      call spdeflect(np,x,px,y,py,z,g,dv,zz,
     $     ex,ey,ez,
     $     nx,nz,nz1,nz2,npz,ilist(1:np,itt),dx,dy,zc,
     $     al,alx,v0)
      call tfree(int8(iphis))
      call tfree(int8(iq))
      call tfree(int8(iey))
      call tfree(int8(iez))
      call tfree(int8(iphi))
      call tfree(int8(irho))
      call tfree(int8(itt))
      return
      end

      subroutine spdeflect(np,x,px,y,py,z,g,dv,zz,
     $     ex,ey,ez,
     $     nx,nz,nz1,nz2,npz,itab,dx,dy,zc,al,alx,v0)
      use ffs
      use tffitcode
      use mathfun
      implicit none
      real*8 eps,almin
      parameter (eps=1.d0,almin=0.01d0)
      integer*4 np,npz,nx,nz,i,j,nz1,nz2,itab(np),k,lx,ly,jj,j1
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),zz(np)
      real*8 zc(nz),
     $     ex(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     ey(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     ez(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     dx,dy,al,v0,alpha,v00,x0,y0,ax,ay,az,
     $     ax1,ay1,az1,exii,eyii,ezii,bxii,byii,v,fxii,fyii,dp,
     $     p,h,dh,v1,dz,p1,pr,pr0,alx,
     $     alpha1
c     real*8 spx,spxpx,spy,spypy,sp,spp
c      real*8 exm,eym,ezm,gamma0,sigpx,sigpy,sigp
      v00=p0/h0
      alpha1=charge**2*e/4.d0/pi/ep0/p0/amass*pbunch/np0
c      write(*,*)'spdeflect ',pbunch,np,np0
c      spx=0.d0
c      spxpx=0.d0
c      spy=0.d0
c      spypy=0.d0
c      sp=0.d0
c      spp=0.d0
c      do i=1,np
c        dp=g(i)*(2.d0+g(i))
c        g(i)=dp
c        spx=spx+px(i)
c        spxpx=spxpx+px(i)**2
c        spy=spy+py(i)
c        spypy=spypy+py(i)**2
c        sp=sp+dp
c        spp=spp+dp**2
c      enddo
c      sigpx=sqrt(spxpx/np-(spx/np)**2)
c      sigpy=sqrt(spypy/np-(spy/np)**2)
c      sigp=sqrt(spp/np-(sp/np)**2)
c      alx=min(al,max(almin,min(
c     $     eps*sigpx/(alpha1*exm/gamma0/v00),
c     $     eps*sigpy/(alpha1*eym/gamma0/v00),
c     $     eps*sigp/(alpha1*ezm/v00))))
      alx=al
      alpha=alpha1*alx
      x0=(.5d0-nx)*dx
      y0=(.5d0-nx)*dy

      do j=nz1,nz2-1
        dz=zc(j+1)-zc(j)
        j1=(j-nz1)*npz
        do k=1,min(npz,np-j1)
          jj=j1+k
          i=itab(jj)
          az=(zz(i)-zc(j))/dz
          ax=(x(i)-x0)/dx
          lx=floor(ax)
          ax=ax-lx
          lx=lx-nx
          ay=(y(i)-y0)/dy
          ly=floor(ay)
          ay=ay-ly
          ly=ly-nx
          ax1=1.d0-ax
          ay1=1.d0-ay
          az1=1.d0-az
          exii=ax1*(
     $         ay1*(az1*ex(lx,  ly,  j  )+az*ex(lx,  ly,  j+1))+
     $         ay *(az1*ex(lx,  ly+1,j  )+az*ex(lx,  ly+1,j+1)))+
     $         ax *(
     $         ay1*(az1*ex(lx+1,ly,  j  )+az*ex(lx+1,ly,  j+1))+
     $         ay *(az1*ex(lx+1,ly+1,j  )+az*ex(lx+1,ly+1,j+1)))
          eyii=ax1*(
     $         ay1*(az1*ey(lx,  ly,  j  )+az*ey(lx,  ly,  j+1))+
     $         ay *(az1*ey(lx,  ly+1,j  )+az*ey(lx,  ly+1,j+1)))+
     $         ax *(
     $         ay1*(az1*ey(lx+1,ly,  j  )+az*ey(lx+1,ly,  j+1))+
     $         ay *(az1*ey(lx+1,ly+1,j  )+az*ey(lx+1,ly+1,j+1)))
          ezii=ax1*(
     $         ay1*(az1*ez(lx,  ly,  j  )+az*ez(lx,  ly,  j+1))+
     $         ay *(az1*ez(lx,  ly+1,j  )+az*ez(lx,  ly+1,j+1)))+
     $         ax *(
     $         ay1*(az1*ez(lx+1,ly,  j  )+az*ez(lx+1,ly,  j+1))+
     $         ay *(az1*ez(lx+1,ly+1,j  )+az*ez(lx+1,ly+1,j+1)))
          bxii=-v0*eyii
          byii= v0*exii
          v=(1.d0-dv(i)+dvfs)*v00
          fxii=exii-v*byii
          fyii=eyii+v*bxii
c          dp=g(i)*(2.d0+g(i))
          dp=g(i)
          pr0=1.d0+dp
          px(i)=px(i)+alpha*fxii/v
          py(i)=py(i)+alpha*fyii/v
          p=pr0*p0
          h=p2h(p)
c          h=p*sqrt(1.d0+1.d0/p**2)
          dh=p0*(dp*(p+p0)/(h+h0)+alpha*ezii)
          h=h0+dh
          p1=h2p(h)
c          p1=h*sqrt(1.d0-1.d0/h**2)
c          p1=sqrt((h-1.d0)*(h+1.d0))
          dp=dh*(h+h0)/(p1+p0)/p0
          pr=1.d0+dp
c          g(i)=dp/(1.d0+sqrt(pr))
          g(i)=dp
          dv(i)=-dp*(1.d0+pr)/h/(h+pr*h0)+dvfs
          v1=(1.d0-dv(i)+dvfs)*v00
          z(i)=v1/v*z(i)
          px(i)=px(i)/pr
          py(i)=py(i)/pr
        enddo
      enddo
      return
      end

      subroutine spfield(phi,ex,ey,ez,nx,nz,nz1,nz2,
     $     dx,dy,zc,gamma0)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nx,nz,i,j,m,nz1,nz2
      real*8 phi(-nx-1:nx,-nx-1:nx,nz1-1:nz2+1),zc(nz),
     $     ex(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     ey(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     ez(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     dx,dy,gamma0,ax,ay,az0,az1,az2,dz0,dz1,dz
      ax=-.5d0/dx
      ay=-.5d0/dy
      do m=nz1,nz2
        dz0=zc(m)-zc(m-1)
        dz1=zc(m+1)-zc(m)
        dz=dz0+dz1
        az0= dz1/dz0/dz/gamma0
        az2=-dz0/dz1/dz/gamma0
        az1=-(1.d0/dz0-1.d0/dz1)/gamma0
        do i=-nx,nx-1
          do j=-nx,nx-1
            ex(i,j,m)=ax*(phi(i+1,j,m)-phi(i-1,j,m))
            ey(i,j,m)=ay*(phi(i,j+1,m)-phi(i,j-1,m))
            ez(i,j,m)=az0*phi(i,j,m-1)
     $           +az1*phi(i,j,m)+az2*phi(i,j,m+1)
          enddo
        enddo
      enddo
c      m=(nz1+nz2)/2
c      do i=-nx,nx-1
c        write(*,'(:1p4(g15.7,a))')
c     $       (i+.5d0)*dx,char(9),ex(i,0,m),char(9),
c     $       ey(0,i,m),char(9),ez(i,i,m)
c      enddo
c      do m=nz1,nz2
c        write(*,'(:1p4(g15.7,a))')
c     $       zc(m),char(9),ex(0,0,m),char(9),
c     $       ey(0,0,m),char(9),ez(0,0,m)
c      enddo
      return
      end

      subroutine spphi(rho,phi,q,phis,nx,nz,nz1,nz2,np,
     $     dx,dy,zc,zf,radius,ntheta)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 itmax
      parameter (itmax=20)
      integer*4 nx,nz,np,it,nz1,nz2,ntheta,i
      real*8 eps
      parameter (eps=0.01d0)
      real*8 dx,dy,zc(nz),zf(nz),radius,dtheta,
     $     rho(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     phi(-nx-1:nx,-nx-1:nx,nz1-1:nz2+1),
     $     q(ntheta,nz),phis(ntheta,nz),sx(ntheta),sy(ntheta),
     $     dqabs,drtheta(1-ntheta:ntheta-1)
      dtheta=pi/ntheta
      do i=1-ntheta,ntheta-1
        drtheta(i)=abs(2.d0*radius*sin((i+.5d0)*dtheta))
      enddo
      do i=1,ntheta
        sx(i)=radius*cos((i-.5d0)*2.d0*dtheta)
        sy(i)=radius*sin((i-.5d0)*2.d0*dtheta)
      enddo
      call spphis(rho,phis,nx,nz,nz1,nz2,dx,dy,zc,zf,sx,sy,ntheta)
      it=0
 1    call spphi1(q,phis,nz,nz1,nz2,
     $     zc,zf,drtheta,radius,dqabs,ntheta)
      if(dqabs .gt. eps*np)then
        it=it+1
        if(it .lt. itmax)then
          go to 1
        endif
      endif
      call spphi2(q,rho,phi,nx,nz,nz1,nz2,
     $     dx,dy,zc,radius,ntheta)
      return
      end

      subroutine spphi2(q,rho,phi,nx,nz,nz1,nz2,
     $     dx,dy,zc,radius,ntheta)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 itmax
      real*8 eps,acc
      parameter (eps=0.0001d0,itmax=100,acc=1.5d0)
      integer*4 nx,nz,nz1,nz2,ntheta,i,j,k,l,m,n,
     $     i1,j1,k1,nd,it,l1
      real*8 dx,dy,zc(nz),y1,y2,x1,x2,d1,d2,dxy1,dxy2,
     $     alpha1,alpha2,sphi,sdphi,phic,
     $     x,y,dxh,dyh,dtheta,radius,a,
     $     az0(nz1:nz2),az1(nz1:nz2),az2(nz1:nz2),
     $     ax0,ax1,ay0,ay1,v,dphic,
     $     rho(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     phi(-nx-1:nx,-nx-1:nx,nz1-1:nz2+1),
     $     q(ntheta,nz),dz,sx,sy,z1,z2,asinh,dmin
c      real*8 dxy1a,dxy1b,dxy1c,dxy1d,dxy2a,dxy2b,dxy2c,dxy2d
      dmin=((dx-dy)**3*dx*dy/(dx**2+dy**2)**2)**2
      dxh=dx*.5d0
      dyh=dy*.5d0
      nd=2*nx+1
      do i=0,nx-1
        x=i*dx+dxh
        i1=-1-i
        do j=0,nx-1
          y=j*dy+dyh
          j1=-1-j
          do k=-nx-1,nx,nd
            k1=-1-k
            x1=k*dx+dxh
            y2=k*dy+dyh
            do l=-nx,nx-1
              l1=-1-l
              y1=l*dy+dyh
              x2=l*dx+dxh
              dxy1=sqrt((x1-x)**2+(y1-y)**2+dmin)
c              dxy1a=sqrt((x1-x+dxh)**2+(y1-y+dyh)**2)
c              dxy1b=sqrt((x1-x-dxh)**2+(y1-y+dyh)**2)
c              dxy1c=sqrt((x1-x+dxh)**2+(y1-y-dyh)**2)
c              dxy1d=sqrt((x1-x-dxh)**2+(y1-y-dyh)**2)
              dxy2=sqrt((x2-x)**2+(y2-y)**2+dmin)
c              dxy2a=sqrt((x2-x+dxh)**2+(y2-y+dyh)**2)
c              dxy2b=sqrt((x2-x-dxh)**2+(y2-y+dyh)**2)
c              dxy2c=sqrt((x2-x+dxh)**2+(y2-y-dyh)**2)
c              dxy2d=sqrt((x2-x-dxh)**2+(y2-y-dyh)**2)
              do n=nz1-1,nz2+1
                do m=nz1,nz2
                  z1=(zc(m-1)+zc(m))*.5d0-zc(n)
                  z2=(zc(m+1)+zc(m))*.5d0-zc(n)
                  alpha1=(asinh(z2/dxy1)-asinh(z1/dxy1))/(z2-z1)
c     $                 asinh(z2/dxy1a)-asinh(z1/dxy1a)+
c     $                 asinh(z2/dxy1b)-asinh(z1/dxy1b)+
c     $                 asinh(z2/dxy1c)-asinh(z1/dxy1c)+
c     $                 asinh(z2/dxy1d)-asinh(z1/dxy1d)
                  alpha2=(asinh(z2/dxy2)-asinh(z1/dxy2))/(z2-z1)
c     $                 asinh(z2/dxy2a)-asinh(z1/dxy2a)+
c     $                 asinh(z2/dxy2b)-asinh(z1/dxy2b)+
c     $                 asinh(z2/dxy2c)-asinh(z1/dxy2c)+
c     $                 asinh(z2/dxy2d)-asinh(z1/dxy2d)
                  phi(k, l, n)=phi(k, l, n)+alpha1*rho(i, j, m)
                  phi(k1,l, n)=phi(k1,l, n)+alpha1*rho(i1,j, m)
                  phi(k, l1,n)=phi(k, l1,n)+alpha1*rho(i, j1,m)
                  phi(k1,l1,n)=phi(k1,l1,n)+alpha1*rho(i1,j1,m)
                  phi(l, k, n)=phi(l, k, n)+alpha2*rho(i, j, m)
                  phi(l1,k, n)=phi(l1,k, n)+alpha2*rho(i1,j, m)
                  phi(l, k1,n)=phi(l, k1,n)+alpha2*rho(i, j1,m)
                  phi(l1,k1,n)=phi(l1,k1,n)+alpha2*rho(i1,j1,m)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      dtheta=2.d0*pi/ntheta
      do i=1,ntheta
        sx=radius*cos((i-1)*dtheta)
        sy=radius*sin((i-1)*dtheta)
        do m=1,nz
          do n=nz1-1,nz2+1
            dz=(zc(n)-zc(m))**2+dmin
            do k=-nx-1,nx,nd
              x1=dx*k+dxh
              y2=dy*k+dyh
              d1=(x1-sx)**2+dz
              d2=(y2-sy)**2+dz
              do l=-nx,nx-1
                y1=dy*l+dyh
                x2=dx*l+dxh
                phi(k,l,n)=phi(k,l,n)+q(i,m)/sqrt((y1-sy)**2+d1)
                phi(l,k,n)=phi(l,k,n)+q(i,m)/sqrt((x2-sx)**2+d2)
              enddo
            enddo
          enddo
        enddo
      enddo
      it=0
      do m=nz1,nz2
        az0(m)=2.d0/(zc(m)-zc(m-1))/(zc(m+1)-zc(m-1))
        az1(m)=2.d0/(zc(m)-zc(m-1))/(zc(m+1)-zc(m))
        az2(m)=2.d0/(zc(m+1)-zc(m))/(zc(m+1)-zc(m-1))
      enddo
      ax0=1.d0/dx**2
      ax1=2.d0/dx**2
      ay0=1.d0/dy**2
      ay1=2.d0/dy**2
 1    sphi=0.d0
      sdphi=0.d0
      do m=nz1,nz2
        a=1.d0/(az1(m)+ax1+ay1)
        v=8.d0*pi/(dx*dy*(zc(m+1)-zc(m-1)))
        do i=-nx,nx-1
          do j=-nx,nx-1
            phic=(az0(m)*phi(i,j,m-1)+az2(m)*phi(i,j,m+1)
     $           +(phi(i-1,j,m)+phi(i+1,j,m))*ax0
     $           +(phi(i,j-1,m)+phi(i,j+1,m))*ay0
     $           +v*rho(i,j,m))*a
            dphic=phic-phi(i,j,m)
            sdphi=sdphi+abs(dphic)
            sphi=sphi+abs(phic)
            phi(i,j,m)=phi(i,j,m)+acc*dphic
          enddo
        enddo
      enddo
c      write(*,*)'spphi2 ',it,sdphi,sphi
      if(sdphi .gt. sphi*eps)then
        it=it+1
        if(it .lt. itmax)then
          go to 1
        endif
      endif
      return
      end

      subroutine spphi1(q,phis,nz,nz1,nz2,
     $     zc,zf,drtheta,radius,dqabs,ntheta)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nz,i,j,l,m,jj,ntheta,i1,nz1,nz2
      real*8 zc(nz),zf(nz),alpha,dqa,r,q1,dz,
     $     q(ntheta,nz),phis(ntheta,nz),fact,radius,asinh,
     $     dqabs,drtheta(1-ntheta:ntheta-1),
     $     dq(ntheta)
      dqa=0.d0
      fact=-.4d0*pi/asinh(.5d0*(zc(nz2)-zc(nz1))/radius)/ntheta
      do m=1,nz
c        alpha=-sqrt(drtheta(0)**2+(zc(m)-zf(m))**2)/64.d0
        alpha=fact*(zc(m)-zf(m))
        do i=1,ntheta
          if(i .eq. 1)then
            i1=ntheta
          else
            i1=i-1
          endif
          q1=alpha*(phis(i,m)+phis(i1,m))
          dqa=dqa+abs(q1)
          dq(i)=q1
          q(i,m)=q(i,m)+q1
        enddo
c        write(*,*)'spphi1 ',m,alpha,q1,
c     $       phis(ntheta,m),phis(ntheta-1,m)
        do l=1,nz
          dz=(zc(m)-zf(l))**2
          do i=0,ntheta-1
            r=sqrt(drtheta(i)**2+dz)
            do j=1,ntheta
              jj=i+j
              if(jj .gt. ntheta)then
                jj=jj-ntheta
              endif
              phis(jj,l)=phis(jj,l)+dq(j)/r
            enddo
          enddo
        enddo
      enddo
      dqabs=dqa
      return
      end

      subroutine spphis(rho,phis,nx,nz,nz1,nz2,
     $     dx,dy,zc,zf,sx,sy,ntheta)
      implicit none
      integer*4 nx,nz,nz1,nz2,ntheta,i,j,k,l,m
      real*8 dx,dy,zc(nz),zf(nz),
     $     rho(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     phis(ntheta,nz),sx(ntheta),sy(ntheta),
     $     dxh,dyh,x,y,d
c      logical*4 isnan
      dxh=dx*.5d0
      dyh=dy*.5d0
      do i=-nx,nx-1
        x=i*dx+dxh
        do j=-nx,nx-1
          y=j*dy+dyh
          do k=1,ntheta
            d=(x-sx(k))**2+(y-sy(k))**2
            do m=nz1,nz2
              do l=1,nz
                phis(k,l)=phis(k,l)+
     $               rho(i,j,m)/sqrt(d+(zc(m)-zf(l))**2)
              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine sprho(rho,nx,nz,nz1,nz2,dx,dy,zc,
     $     np,npz,x,y,zz,itab)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nx,nz1,nz2,np,npz,itab(np),
     $     ix1,ix2,iy1,iy2,iz1,iz2,
     $     i,j,jj,nz,j1
      real*8 dx,zc(nz),x(np),y(np),zz(np),
     $     rho(-nx:nx-1,-nx:nx-1,nz1:nz2),
     $     ax,ay,rx1,rx2,ry1,ry2,rz1,rz2,dz,x0,y0,dy,
     $     dmin
c      logical*4 isnan
      dmin=abs(((dx-dy)**3*dx*dy/(dx**2+dy**2)**2))
      x0=dx*(-nx+0.5d0)
      y0=dy*(-nx+0.5d0)
      do iz1=nz1,nz2-1
        iz2=iz1+1
        dz=max(zc(iz2)-zc(iz1),dmin)
        j1=(iz1-nz1)*npz
        do jj=1,min(npz,np-j1)
          j=j1+jj
          i=itab(j)
          ax=(x(i)-x0)/dx
          ix1=min(2*nx-2,max(0,int(ax)))
          rx2=ax-ix1
          rx1=1.d0-rx2
          ix1=ix1-nx
          ix2=ix1+1
          ay=(y(i)-y0)/dy
          iy1=min(2*nx-2,max(0,int(ay)))
          ry2=ay-iy1
          ry1=1.d0-ry2
          iy1=iy1-nx
          iy2=iy1+1
          rz2=(zz(i)-zc(iz1))/dz
          rz1=1.d0-rz2
          rho(ix1,iy1,iz1)=rho(ix1,iy1,iz1)+rx1*ry1*rz1
          rho(ix1,iy1,iz2)=rho(ix1,iy1,iz2)+rx1*ry1*rz2
          rho(ix1,iy2,iz1)=rho(ix1,iy2,iz1)+rx1*ry2*rz1
          rho(ix1,iy2,iz2)=rho(ix1,iy2,iz2)+rx1*ry2*rz2
          rho(ix2,iy1,iz1)=rho(ix2,iy1,iz1)+rx2*ry1*rz1
          rho(ix2,iy1,iz2)=rho(ix2,iy1,iz2)+rx2*ry1*rz2
          rho(ix2,iy2,iz1)=rho(ix2,iy2,iz1)+rx2*ry2*rz1
          rho(ix2,iy2,iz2)=rho(ix2,iy2,iz2)+rx2*ry2*rz2
        enddo
      enddo
      return
      end

      subroutine spmesh(nx,nz,nz1,nz2,dx,dy,zc,zf,
     $     np,npz,x,y,zz,radius,itab,nzmax)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nxmin,nzmax,nzmin
      parameter (nxmin=4,nzmin=10)
      integer*4 np,i,itab(np),nx,nz,npz,nza,na1,na2,nz1,nz2
      real*8 x(np),y(np),zz(np),dx,dy,zc(nzmax),xm,ym,
     $     radius,dz1,dz2,zf(nzmax)
      nza=min(nzmax-100,max(nint(5.d0*dble(np/1000)**0.67d0),nzmin))
      npz=(np+nza-1)/nza
c      write(*,*)'spmesh ',np,nzmax,npz,nza,itab(npz+1),itab(1),
c     $     itab(np),itab((nza-1)*npz+1)
      dz1=zz(itab(npz+1))-zz(itab(1))
      dz2=zz(itab(np))-zz(itab((nza-1)*npz+1))
      na1=max(1,floor(radius*3.d0/max(dz1,radius/10.d0)))
      nz1=na1+1
      nz2=nz1+nza
      na2=max(1,floor(radius*3.d0/max(dz2,radius/10.d0)))
      nz=nz2+na2
      if(nz .gt. nzmax)then
        write(*,*)'spmesh-too large mesh'
      endif
      do i=0,nza-1
        zc(nz1+i)=zz(itab(i*npz+1))
      enddo
      zc(nz2)=zz(itab(np))
      do i=1,na1
        zc(i)=(i-na1-1)*dz1+zc(nz1)
      enddo
      do i=1,na2
        zc(nz2+i)=i*dz2+zc(nz2)
      enddo
      zf(1)=zc(1)-dz1*.5d0
      do i=2,nz
        zf(i)=(zc(i)+zc(i-1))*.5d0
      enddo
      xm=0.d0
      ym=0.d0
      do i=1,np
        xm=max(xm,abs(x(i)))
        ym=max(ym,abs(y(i)))
      enddo
      xm=min(radius/sqrt(2.d0),xm+1.d-9)
      ym=min(radius/sqrt(2.d0),ym+1.d-9)
      nx=max(nxmin,nint(sqrt(sqrt(dble(npz)))/2.d0))
      dx=2.d0*xm/(2*nx-1)
      dy=2.d0*ym/(2*nx-1)
      if(max(dx,dy) .gt. 31.d0*min(dx,dy))then
        write(*,*)'spmesh-too flat beam: ',max(dx,dy)/min(dx,dy)
      endif
      return
      end

      subroutine spzz(np,px,py,z,g,zz,v0,gamma0)
      use ffs
      use tffitcode
      use mathfun
      implicit none
      integer*4 np,i
      real*8 px(np),py(np),z(np),g(np),zz(np),v0,sdp,gamma0,p,h
      sdp=0.d0
      do i=1,np
c        sdp=sdp+g(i)*(2.d0+g(i))
        sdp=sdp+g(i)
      enddo
      sdp=sdp/np
      p=p0*(1.d0+sdp)
      h=p2h(p)
c      h=p*sqrt(1.d0+1.d0/p**2)
c      h=sqrt(1.d0+p**2)
      v0=p/h
      gamma0=h
      do i=1,np
        zz(i)=z(i)*sqrt(
     $       max(0.d0,(1.d0-py(i))*(1.d0+py(i))-px(i)**2))*h
      enddo
      return
      end

      recursive subroutine spsort(n,itab,zz)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 n,itab(n),m,i1,i2,is,im,ip1,ip2
      real*8 zz(n)
      if(n .le. 1)then
        return
      endif
      i1=itab(1)
      i2=itab(n)
      if(zz(i1) .gt. zz(i2))then
        is=i1
        i1=i2
        i2=is
      endif
      if(n .eq. 2)then
        itab(1)=i1
        itab(2)=i2
        return
      endif
      m=(n+1)/2
      im=itab(m)
      if(zz(i1) .lt. zz(im))then
        if(zz(im) .gt. zz(i2))then
          is=im
          im=i2
          i2=is
        endif
      elseif(zz(i1) .ne. zz(im))then
        is=im
        im=i1
        i1=is
      endif
      itab(1)=i1
      itab(n)=i2
      itab(m)=im
      if(n .eq. 3)then
        return
      endif
      itab(m)=itab(2)
      itab(2)=im
      ip1=3
      ip2=n-1
      do while(ip1 .le. ip2)
        do while(zz(itab(ip1)) .le. zz(im) .and. ip1 .le. ip2)
          ip1=ip1+1
        enddo
        do while(zz(im) .le. zz(itab(ip2)) .and. ip1 .le. ip2)
          ip2=ip2-1
        enddo
        if(ip2 .gt. ip1)then
          is=itab(ip1)
          itab(ip1)=itab(ip2)
          itab(ip2)=is
          ip1=ip1+1
          ip2=ip2-1
        endif
      enddo
      ip1=ip1-1
      is=itab(ip1)
      itab(ip1)=im
      itab(2)=is
      call spsort(ip1-1,itab,zz)
      call spsort(n-ip1,itab(ip1+1),zz)
      return
      end

c     drift in the free space
      subroutine spdrift_free(np,x,px,y,py,z,g,dv,sx,sy,sz,al,
     $     radius,kturn,kptbl)
      use tfstk
      use ffs
      use tffitcode
      use tmacro,only:l_track
      implicit none
      integer*4 nzmax
      real*8 alstep
      parameter (nzmax=1000,alstep=0.05d0)
      integer*4 ,intent(inout):: np,kptbl(np0,6)
      integer*4 ,intent(in):: kturn
      real*8 x(np0),px(np0),y(np),py(np0),z(np0),g(np0),dv(np0),
     $     sx(np0),sy(np0),sz(np0)
      real*8 al,radius
      integer*4 ndiv,i
      real*8 aln,alx
      ndiv=max(1,nint(abs(al)/alstep))
      aln=al/ndiv

      call tdrift_free(np,x,px,y,py,z,dv,aln*.5d0)

      call spkick(np,x,px,y,py,z,g,dv,aln,
     $     radius,alx,kturn,kptbl)
      do i=2,ndiv
        call tdrift_free(np,x,px,y,py,z,dv,aln)
        call spkick(np,x,px,y,py,z,g,dv,aln,
     $       radius,alx,kturn,kptbl)
      enddo

      call tdrift_free(np,x,px,y,py,z,dv,aln*.5d0)
c      call spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
      if(radius .ne. 0.d0)then
         call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $       kptbl,np,kturn,
     $        radius,radius,
     $        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      return
      end

c     drift in the parallel solenoid
      subroutine spdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,radius,kturn,kptbl)
      use tfstk
      use ffs
      use tffitcode
      use tmacro,only:l_track
      implicit none
      integer*4 nzmax
      real*8 alstep
      parameter (nzmax=1000,alstep=0.05d0)
      integer*4 ,intent(inout):: np,kptbl(np0,6)
      integer*4 ,intent(in):: kturn
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     $     sx(np0),sy(np0),sz(np0)
      real*8 al,bz,radius
      integer*4 ndiv,i
      real*8 aln,alx
      ndiv=max(1,nint(abs(al)/alstep),nint(abs(bz*al)/1.5d0))
      aln=al/ndiv

      call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     aln*.5d0,bz,.false.)

      call spkick(np,x,px,y,py,z,g,dv,aln,
     $     radius,alx,kturn,kptbl)
      do i=2,ndiv
        call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       aln,bz,.false.)
        call spkick(np,x,px,y,py,z,g,dv,aln,
     $       radius,alx,kturn,kptbl)
      enddo

      call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     aln*.5d0,bz,.false.)
c      call spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
      if(radius .ne. 0.d0)then
         call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $       kptbl,np,kturn,
     $        radius,radius,
     $        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      return
      end

      subroutine spdrift(np,
     $     x,px,y,py,z,g,dv,sx,sy,sz,al,bz,ak0x,ak0y,radius,kptbl)
cProbably obsolete
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nzmax
      real*8 alstep
      parameter (nzmax=1000,alstep=0.05d0)
      integer*4 ,intent(inout):: np,kptbl(np0,6)
      integer*4 kturn
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     $     sx(np0),sy(np0),sz(np0)
      real*8 al,bz,ak0x,ak0y,radius
      integer*4 ndiv,i
      real*8 aln,akxn,akyn,alx
      ndiv=max(1,nint(abs(al)/alstep),nint(abs(bz*al)/1.5d0))
      aln=al/ndiv
      akxn=ak0x/ndiv
      akyn=ak0y/ndiv

      call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     aln*.5d0,bz,akxn*.5d0,akyn*.5d0,.false.)

      call spkick(np,x,px,y,py,z,g,dv,aln,
     $     radius,alx,kturn,kptbl)
      do i=2,ndiv
        call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       aln,bz,akxn,akyn,.false.)
        call spkick(np,x,px,y,py,z,g,dv,aln,
     $       radius,alx,kturn,kptbl)
      enddo

      call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     aln*.5d0,bz,akxn*.5d0,akyn*.5d0,.false.)
c      call spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
      if(radius .ne. 0.d0)then
         call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $       kptbl,np,kturn,
     $        radius,radius,
     $        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      return
      end

      subroutine spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
cObsolete
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 np,kptbl(np0,6),i,j,k
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     radius,rsq,x1,px1,y1,py1,z1,g1,dv1
      if(radius .eq. 0.d0)then
        return
      endif
      rsq=radius**2
      do i=1,np
 20     if(i .gt. np)then
          return
        endif
c
c     Don't flip aperture conditional!!
c     This conditional is designed to drop NaN particle.
c     In IEEE754 standard, any comparision operation
c     with NaN operand is defined as ``False''.
        if(.not. (x(i)**2+y(i)**2 .le. rsq
     $       .and. abs(z(i)) .le. zlost))then
c          write(*,*)'spapert ',i,x(i),y(i),
c     $         z(i),radius,zlost
          j=kptbl(np,3)
          k=kptbl(i,3)
          kptbl(k,1)=np
          kptbl(j,1)=i
          kptbl(np,3)=k
          kptbl(i,3)=j
            x1=x(i)
            px1=px(i)
            y1=y(i)
            py1=py(i)
            z1=z(i)
            g1=g(i)
            dv1=dv(i)
            x(i)=x(np)
            px(i)=px(np)
            y(i)=y(np)
            py(i)=py(np)
            z(i)=z(np)
            g(i)=g(np)
            dv(i)=dv(np)
            x(np)=x1
            px(np)=px1
            y(np)=y1
            py(np)=py1
            z(np)=z1
            g(np)=g1
            dv(np)=dv1
          np=np-1
          go to 20
        endif
      enddo
      return
      end
