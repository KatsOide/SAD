      subroutine tsteer(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phib,dx,dy,
     1     theta,cost,sint,
     1     cosp1,sinp1,cosp2,sinp2,
     $     fb1,fb2,fringe,eps,enarad)
      use ffs_flag
      use tmacro
      use tspin
      use bendib, only:tbendal
      use photontable, only:tsetphotongeo
      use mathfun, only:pxy2dpz
      implicit none
      real*8 , parameter :: a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1     a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1     a13=231.d0/13312.d0,a15=143.d0/10240.d0
      integer*4 , parameter :: ndivmax=1000
      integer*4 np,i,ndiv,n,n1,n2
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      real*8 al,phib,dx,dy,theta,cost,sint,rhob,eps,aln,
     $     dxfr1,dyfr1,dyfra1,dyi,alc,alx,alr,
     $     dxfr2,dyfr2,dyfra2,f1r,f2r,phic,
     $     dp,p,rhoe,pxi,s,dpv1,pv1,dpv2,pv2,fa,f,ff,
     $     dpz1,pz1,dpz2,pz2,phsq,u,w,dl,dpx,pyi,xi,pxf,d,
     $     cosp1,sinp1,cosp2,sinp2,tanp1,tanp2,fb1,fb2
      logical*4 fringe,enarad,krad
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,2,al,-phib,
     1             dx,dy,theta,cost,sint,1.d0,.false.)
        return
      elseif(phib .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,dv,al)
        return
      endif
      include 'inc/TENT.inc'
      tanp1=sinp1/cosp1
      tanp2=0.d0
      rhob=al/phib
      krad=enarad .and. rad
      f1r=0.d0
      f2r=0.d0
      dxfr1=0.d0
      dyfr1=0.d0
      dyfra1=0.d0
      dxfr2=0.d0
      dyfr2=0.d0
      dyfra2=0.d0
      n1=1
      n2=0
      if(fb1 .ne. 0.d0)then
        dxfr1=-fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        if(fringe)then
          dyfra1=4.d0*dyfr1/fb1**2
        endif
        if(krad)then
          f1r=0.5d0*fb1
          n1=-1
        endif
      endif
      if(fb2 .ne. 0.d0 .and. krad)then
        f2r=.5d0*fb2
        n2=2
      endif
      phic=phib*(1.d0-(f1r+f2r)/al)
      if(krad)then        
        px0=px
        py0=py
        zr0=z
        ndiv=min(ndivmax,max(1,ndivrad(phic,0.d0,0.d0,eps)))
      else
        ndiv=1
      endif
      n2=ndiv+n2
      alc=al*phic/phib
      aln=alc/ndiv
      do n=n1,n2
        call tbendal(n,ndiv,f1r,f2r,aln,alx,alr)
        if(n .eq. n1+1)then
          dxfr1=0.d0
          dyfr1=0.d0
          dyfra1=0.d0
          tanp1=0.d0
        endif
        if(n .eq. n2 .and. fb2 .ne. 0.d0)then
          dxfr2=-fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          if(fringe)then
            dyfra2=4.d0*dyfr2/fb2**2
          endif
          tanp2=sinp2/cosp2
        endif
        do i=1,np
          dp=g(i)
          p=1.d0+dp
          rhoe=rhob*p
          pxi=px(i)-tanp1*x(i)/rhoe
          x(i)=x(i)+dxfr1*dp/p
          py(i)=py(i)+((dyfr1-dyfra1*y(i)**2)/p**2+tanp1/rhoe)*y(i)
          z(i)=z(i)+(dxfr1*pxi+
     $         (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/p)/p
          dpv1=pxy2dpz(pxi,0.d0)
          pv1=1.d0+dpv1
          fa=y(i)/rhoe/pv1
          f=(1.d0-(y(i)/rhob)**2/6.d0)*fa
          ff=(.5d0-(y(i)/rhob)**2/24.d0)*y(i)*fa/pv1**2
          pyi=py(i)+pxi*f
          xi  =x(i)-ff
          z(i)=z(i)+ff*pxi
          dpz1=pxy2dpz(pxi,pyi)
          pz1=1.d0+dpz1
          dpx=alx/rhoe
          pxf=pxi+dpx
          dpz2=pxy2dpz(pxf,pyi)
          pz2=1.d0+dpz2
          d=pxf*pz1+pxi*pz2
          if(d .eq. 0.d0 .or. pxi*pxf .lt. 0.d0)then
            phsq=pxi**2+pz1**2
            u=(dpx+pxf*dpz1-pxi*dpz2)/phsq
            w=(pxf*dpz1-pxi*dpz2+dpx*pyi**2)/phsq
          else
            u=dpx*(pxf+pxi)/d
            w=-dpx*(pxf*dpz1+pxi*dpz2)/d
          endif
          s=u**2
          if(s .gt. 2.d-2)then
            dl=rhoe*((asin(u)/u-1.d0)*u+w)
          else
            if(s .gt. 2.d-4)then
              dl=rhoe*(s*(a3+s*(a5+s*(a7+s*(a9+s*(a11+s*(a13+s*a15))))))
     $             *u+w)
            else
              dl=rhoe*(s*(a3+s*(a5+s*(a7+s*a9)))*u+w)
            endif
          endif
          x(i)=xi+alx*(pxf+pxi)/(pz2+pz1)
          dyi=pyi*(alx+dl)
          bsi(i)=-dyi/rhob
          y(i)=y(i)+dyi
          dpv2=pxy2dpz(pxf,0.d0)
          pv2=1.d0+dpv2
          fa=y(i)/rhoe/pv2
          f=(1.d0-(y(i)/rhob)**2/6.d0)*fa
          ff=(.5d0-(y(i)/rhob)**2/24.d0)*y(i)*fa/pv2**2
          py(i)=pyi-pxf*f
          x(i)=x(i)+ff
          px(i)=pxf-x(i)*tanp2/rhoe
          z(i)=z(i)-(dl+dv(i)*alx+ff*pxf)
          x(i)=x(i)-dxfr2*dp/p
          py(i)=py(i)+((dyfr2-dyfra2*y(i)**2)/p**2+tanp2/rhoe)*y(i)
          z(i)=z(i)-(dxfr2*pxf-
     $         (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/p)/p
        enddo
        if(krad)then
          if(photons)then
            if(n .eq. 1)then
              call tsetphotongeo(alx,0.d0,theta,.true.)
            else
              call tsetphotongeo(alx,0.d0,theta,.false.)
            endif
          endif
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         px0,py0,zr0,bsi,alr,0.d0)
        endif
      enddo
      include 'inc/TEXIT.inc'
      return
      end
