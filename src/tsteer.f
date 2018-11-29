      subroutine tsteer(np,x,px,y,py,z,g,dv,pz,l,al,phib,dx,dy,
     1     theta,cost,sint,
     1     cosp1,sinp1,cosp2,sinp2,
     $     fb1,fb2,mfring,fringe,enarad,eps0)
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      implicit none
      real*8 a3,a5,a7,a9,a11,a13,a15
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      integer*4 np,mfring,i,l
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np),
     $     px0(np),py0(np),ds(np)
      real*8 al,phib,dx,dy,theta,cost,sint,eps0,rhob,
     $     dxfr1,dyfr1,dyfra1,f1r,f2r,
     $     dxfr2,dyfr2,dyfra2,
     $     dp,p,rhoe,pxi,s,dpv1,pv1,dpv2,pv2,fa,f,ff,
     $     dpz1,pz1,dpz2,pz2,phsq,u,w,dl,brad,dpx,pyi,xi,pxf,d,
     $     cosp1,sinp1,cosp2,sinp2,tanp1,tanp2,fb1,fb2
      logical*4 fringe,enarad,enrad
      enrad=enarad .and. rad
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,pz,2,l,al,-phib,
     1             dx,dy,theta,cost,sint,1.d0,.false.)
        return
      elseif(phib .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        return
      endif
      if(enrad .and. trpt)then
        call tstrad(np,x,px,y,py,z,g,dv,pz,l,al,phib,dx,dy,
     1       theta,cost,sint,
     1       cosp1,sinp1,cosp2,sinp2,
     $       fb1,fb2,mfring,eps0)
        return
      endif
      include 'inc/TENT.inc'
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      if(enrad)then
        px0=px
        py0=py
c        if(iprev(l) .eq. 0)then
c          f1r=fb1
c        else
c          f1r=0.d0
c        endif
c        if(inext(l) .eq. 0)then
c          f2r=fb2
c        else
c          f2r=0.d0
c        endif
c        brad=brhoz/rhob
c        call trad(np,x,px,y,py,g,dv,brad,0.d0,0.d0,
c     1             0.d0,-tanp1*2.d0/al,.5d0*al,
c     $       f1r,f2r,0.d0,al,1.d0)
      endif
      dxfr1=0.d0
      dyfr1=0.d0
      dyfra1=0.d0
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=-fb1**2/rhob/24.d0
          dyfr1=fb1/rhob**2/6.d0
          if(fringe)then
            dyfra1=4.d0*dyfr1/fb1**2
          endif
        endif
      endif
      dxfr2=0.d0
      dyfr2=0.d0
      dyfra2=0.d0
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=-fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          if(fringe)then
            dyfra2=4.d0*dyfr2/fb2**2
          endif
        endif
      endif
      do 100 i=1,np
c        dp=g(i)*(2.d0+g(i))
        dp=g(i)
        p=1.d0+dp
        rhoe=rhob*p
        pxi=px(i)-tanp1*x(i)/rhoe
        x(i)=x(i)+dxfr1*dp/p
        py(i)=py(i)+((dyfr1-dyfra1*y(i)**2)/p**2+tanp1/rhoe)*y(i)
        z(i)=z(i)+(dxfr1*pxi+
     $       (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/p)/p
        s=min(.9d0,pxi**2)
        dpv1=s*(-.5d0-s*(.125d0+s*.0625d0))
        dpv1=(dpv1**2-s)/(2.d0+2.d0*dpv1)
        dpv1=(dpv1**2-s)/(2.d0+2.d0*dpv1)
        pv1=1.d0+dpv1
        fa=y(i)/rhoe/pv1
        f=(1.d0-(y(i)/rhob)**2/6.d0)*fa
        ff=(.5d0-(y(i)/rhob)**2/24.d0)*y(i)*fa/pv1**2
        pyi=py(i)+pxi*f
        xi  =x(i)-ff
        z(i)=z(i)+ff*pxi
        s=min(.9d0,pxi**2+pyi**2)
        dpz1=s*(-.5d0-s*(.125d0+s*(.0625d0+s*.0390625d0)))
        dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
        dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
        pz1=1.d0+dpz1
        dpx=al/rhoe
        pxf=pxi+dpx
        s=min(.9d0,pxf**2+pyi**2)
        dpz2=s*(-.5d0-s*(.125d0+s*(.0625d0+s*.0390625d0)))
        dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
        dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
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
     $           *u+w)
          else
            dl=rhoe*(s*(a3+s*(a5+s*(a7+s*a9)))*u+w)
          endif
        endif
        x(i)=xi+al*(pxf+pxi)/(pz2+pz1)
        if(enrad)then
          ds(i)=al+dl
        endif
        y(i)=y(i)+pyi*(al+dl)
        s=min(.9d0,pxf**2)
        dpv2=s*(-.5d0-s*(.125d0+s*.0625d0))
        dpv2=(dpv2**2-s)/(2.d0+2.d0*dpv2)
        dpv2=(dpv2**2-s)/(2.d0+2.d0*dpv2)
        pv2=1.d0+dpv2
        fa=y(i)/rhoe/pv2
        f=(1.d0-(y(i)/rhob)**2/6.d0)*fa
        ff=(.5d0-(y(i)/rhob)**2/24.d0)*y(i)*fa/pv2**2
        py(i)=pyi-pxf*f
        x(i)=x(i)+ff
        px(i)=pxf-x(i)*tanp2/rhoe
        z(i)=z(i)-(dl+dv(i)*al+ff*pxf)
        x(i)=x(i)-dxfr2*dp/p
        py(i)=py(i)+((dyfr2-dyfra2*y(i)**2)/p**2+tanp2/rhoe)*y(i)
        z(i)=z(i)-(dxfr2*pxf-
     $       (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/p)/p
100   continue
      if(enrad)then
        call tradki(np,x,px,y,py,px0,py0,g,dv,ds)
c        call trad(np,x,px,y,py,g,dv,brad,0.d0,0.d0,
c     1            0.d0,-tanp2*2.d0/al,.5d0*al,
c     $       f1r,f2r,al,al,-1.d0)
      endif
      include 'inc/TEXIT.inc'
      return
      end
