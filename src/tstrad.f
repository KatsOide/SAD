c  Obsolete
c
      subroutine tstrad(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     l,al,phib,dx,dy,
     1     theta,cost,sint,
     1     cosp1,sinp1,cosp2,sinp2,
     $     fb1,fb2,mfring,eps0)
      use ffs_flag
      use tmacro
      implicit none
      integer*4 l,np,mfring,i,nx,ngamma,n
      real*8 al,phib,dx,dy,theta,cost,sint,
     1     cosp1,sinp1,cosp2,sinp2,fb1,fb2,
     $     eps0,rhob,dxfr1,dyfr1,dxfr2,dyfr2,tanp1,tanp2,
     $     dp,p,rhoe,f,ff,rhoba,eps,aln,phin,ur,an,sp,
     $     dpz1,pz1,dpx,pxf,s,dpz2,pz2,d,phsq,u,w,dl,
     $     dprad,dpradx,dprady,bya,h,alr,p1,prob,pyi,pxi,xi
      real*8, parameter :: a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np)
      real*8 tran
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,2,l,al,-phib,
     1             dx,dy,theta,1.d0,.false.)
        return
      elseif(phib .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,dv,al)
        return
      endif
      include 'inc/TENT.inc'
      rhob=al/phib
      if(mfring .gt. 0 .or. mfring .eq. -1)then
        dxfr1=-fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
      else
        dxfr1=0.d0
        dyfr1=0.d0
      endif
      if(mfring .gt. 0 .or. mfring .eq. -2)then
        dxfr2=-fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
      else
        dxfr2=0.d0
        dyfr2=0.d0
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      do 210 i=1,np
c        dp=g(i)*(2.d0+g(i))
        dp=g(i)
c        g(i)=dp
        p=1.d0+dp
        rhoe=rhob*p
        px(i)=px(i)-x(i)*tanp1/rhoe
        f=y(i)/rhoe
        py(i)=py(i)+(px(i)+tanp1)*f
        ff=y(i)*f*.5d0
        x(i)=x(i)-ff
        z(i)=z(i)+ff*px(i)
        x(i)=x(i)+dxfr1*dp/p
        py(i)=py(i)+dyfr1*y(i)/p**2
        z(i)=z(i)+(dxfr1*px(i)+.5d0*dyfr1*y(i)**2/p)/p
210   continue
      rhoba=abs(rhob)
      nx=int(al/rhoba*anrad*p0/.07d0)+1
      if(radlight)then
        if(eps0 .le. 0.d0)then
          eps=1.d-5
        else
          eps=1.d-5*eps0
        endif
        ngamma=int(h0*abs(phib)*1.d-4/eps)+1
        nx=max(nx,ngamma)
      endif
      aln=al/nx
      phin=phib/nx
      ur=urad*p0**3
      an=anrad*p0
      sp=0.d0
      do 10 n=1,nx
        do 100 i=1,np
          pxi=px(i)
          pyi=py(i)
          dp=g(i)
          p=1.d0+dp
          alr=aln*(1.d0+(pxi**2+pyi**2)*.5d0)
          prob=an*alr/rhoba
          if(tran() .gt. 1.d0-prob)then
            bya=brhoz/rhob
            call tsynchrad(p,alr,0.d0,bya,
     $                 dprad,dpradx,dprady,
     $                 i,l,aln*(n-1),0.d0,0.d0,theta,
     $                 x(i),y(i),pxi,pyi)
            dp=dp-dprad*.5d0
            pxi=pxi-dpradx*.5d0
            pyi=pyi-dprady*.5d0
            p=1.d0+dp
            sp=sp-dprad
          else
            dprad=0.d0
            dpradx=0.d0
            dprady=0.d0
          endif
          rhoe=rhob*p
          s=min(.9d0,pxi**2+pyi**2)
          dpz1=s*(-.5d0-s*(.125d0+s*.0625d0))
c         if(s .gt. 1.d-4)then
            dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
            dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
c         endif
          pz1=1.d0+dpz1
          dpx=aln/rhoe
          pxf=pxi+dpx
          s=min(.9d0,pxf**2+pyi**2)
          dpz2=s*(-.5d0-s*(.125d0+s*.0625d0))
c         if(s .gt. 1.d-4)then
            dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
            dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
c         endif
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
          x(i)=x(i)+aln*(pxf+pxi)/(pz2+pz1)
          px(i)=pxf-dpradx*.5d0
          y(i)=y(i)+pyi*(aln+dl)
          py(i)=pyi-dprady*.5d0
          z(i)=z(i)-(dl+dv(i)*aln)
          g(i)=dp-dprad*.5d0
          if(radlight)then
            call tlstore(np,x,y,z,dv,0.d0,aln,0.d0,0.d0,
     $           p0/h0*c,dvfs,.true.)
          endif
100     continue
10    continue
      if(radcod)then
        sp=0.d0
      else
        sp=sp/np
      endif
      do 220 i=1,np
        dp=g(i)-sp
        p=1.d0+dp
        x(i)=x(i)-dxfr2*dp/p
        py(i)=py(i)+dyfr2*y(i)/p**2
        z(i)=z(i)-(dxfr2*px(i)-.5d0*dyfr2*y(i)**2/p)/p
        rhoe=rhob*p
        f=y(i)/rhoe
        py(i)=py(i)+(tanp2-px(i))*f
        ff=y(i)*f*.5d0
        x(i)=x(i)+ff
        z(i)=z(i)-ff*px(i)
        px(i)=px(i)-x(i)*tanp2/rhoe
c        g(i)=dp/(1.d0+sqrt(p))
        g(i)=dp
        p1=p*p0
        h=sqrt(1.d0+p1**2)
        dv(i)=-dp*(1.d0+p)/h/(h+p*h0)+dvfs
220   continue
      include 'inc/TEXIT.inc'
      return
      end
