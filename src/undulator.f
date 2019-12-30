C   18/01/93 303061342  MEMBER NAME  UNDULATOR *.FORT     M  E2FORT
        subroutine undulator(np,x,px,y,py,z,g,dv,sx,sy,sz,ulist)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c  Undulator tracking subroutine.                                 c
c  K.Ohmi and E.Forest        KEK Report 92-14                    c
c                                                                 c
c                           written by K.Ohmi                     c
c                                                                 c
c                                        Dec  1994                c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use tmacro
      use ffs_flag, only:calpol
      implicit real*8 (a-h,o-z)
c      implicit none
      parameter (nsli=1,npole=2,len=3,i_F0=4,i_G0=5,i_dph=6,
     &   i_Kz=7,i_Kx=8,i_Ky=9,i_Qx=10,i_Qy=11,lambda=12,i_ds=13)
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np)
      real*8 ulist(*)
      real*8 sinkz,sinkzp
      real*8 kx1,kx2,kx3,kx4,ff,gg,fx,fy,gx,gy,dpds
      real*8 det1a,fgx,fgy,tmp
      integer i,nslice
c
      nslice=ulist(1)+0.01
c
c   Lorentz factor gamma is stored in h0.
c   Factor gammabeta is stored in p0.
c     gamma=h0
c     gammabet=p0
c
c
c     Initialize dx/dy/theta(cost/sint) for inc/T(ENT|EXIT).inc
      dx=0.d0
      dy=0.d0
      theta=0.d0
c
      include 'inc/TENT.inc'
      do 10 i=1,np
c
c   Get Canonical variables
c       g(i)=pn-1.d0   Not good for precision
c       g(i)=g(i)*(g(i)+2.d0)
c   pn=|p|/p0=1+delta
       pn=1.d0+g(i)
       px(i)=px(i)*pn
       py(i)=py(i)*pn
 10    continue

      do 11 i=1,np
       do 21 j=1,nslice
         kx1=ulist(i_Kx)*x(i)
         kx2=ulist(i_Ky)*y(i)
         kx3=ulist(i_Qx)*x(i)
         kx4=ulist(i_Qy)*y(i)
         sinkz=sin(ulist(i_Kz)*j*ulist(i_ds))
         sinkzp=sin(ulist(i_Kz)*j*ulist(i_ds)+ulist(i_dph))
         ff=ulist(i_F0)/ulist(i_Kz)*cos(kx1)*cosh(kx2)*sinkz
     &       -ulist(i_Qy)*ulist(i_G0)/ulist(i_Qx)/ulist(i_Kz)*
     &        sinh(kx3)*sin(kx4)*sinkzp
         gg=-ulist(i_G0)/ulist(i_Kz)*cos(kx4)*cosh(kx3)*sinkzp
     &       +ulist(i_Kx)*ulist(i_F0)/ulist(i_Ky)/ulist(i_Kz)*
     &        sinh(kx2)*sin(kx1)*sinkz
c
         fx=-ulist(i_Kx)*ulist(i_F0)/
     &        ulist(i_Kz)*sin(kx1)*cosh(kx2)*sinkz
     &       -ulist(i_Qy)*ulist(i_G0)/
     &        ulist(i_Kz)*cosh(kx3)*sin(kx4)*sinkzp
         fy=ulist(i_Ky)*ulist(i_F0)/
     &        ulist(i_Kz)*cos(kx1)*sinh(kx2)*sinkz
     &       -ulist(i_Qy)*ulist(i_Qy)*ulist(i_G0)/ulist(i_Qx)/
     &        ulist(i_Kz)*sinh(kx3)*cos(kx4)*sinkzp
c
         gx=-ulist(i_G0)*ulist(i_Qx)/
     &        ulist(i_Kz)*cos(kx4)*sinh(kx3)*sinkzp
     &       +ulist(i_Kx)*ulist(i_Kx)*ulist(i_F0)/ulist(i_Ky)/
     &        ulist(i_Kz)*sinh(kx2)*cos(kx1)*sinkz
         gy=ulist(i_G0)*ulist(i_Qy)/
     &        ulist(i_Kz)*sin(kx4)*cosh(kx3)*sinkzp
     &       +ulist(i_Kx)*ulist(i_F0)/
     &        ulist(i_Kz)*cosh(kx2)*sin(kx1)*sinkz
         dpds=ulist(i_ds)/(1.+g(i))
c
         fx=dpds*fx
         gx=dpds*gx
         fy=dpds*fy
         gy=dpds*gy
         det1a=1.-fx-gy+fx*gy-fy*gx
c
         fgx=dpds*(ff*fx+gg*gx)
         fgy=dpds*(ff*fy+gg*gy)
c
         tmp=fx
         fx=(1.-gy)/det1a
         gx=gx/det1a
         fy=fy/det1a
         gy=(1.-tmp)/det1a
c
         tmp=px(i)-fgx
         px(i)=fx*tmp+gx*(py(i)-fgy)
         py(i)=fy*tmp+gy*(py(i)-fgy)
c
         fgx=(px(i)-ff)
         fgy=(py(i)-gg)
         x(i)=x(i)+dpds*fgx
         y(i)=y(i)+dpds*fgy
         z(i)=z(i)-0.5*(fgx*fgx+fgy*fgy)*(dpds*dpds)/ulist(i_ds)
 21    continue
 11    continue
c
c   Return to SAD variables
c
c   pn=|p|/p0=1+delta
      do 19 i=1,np
      pn=1.d0+g(i)
      px(i)=px(i)/pn
      py(i)=py(i)/pn
c   When you change momentum, change dv(i) as
      h1=sqrt(1.d0+(pn*p0)**2)
      dv(i)=-g(i)*(1.d0+pn)/h1/(h1+pn*h0)+dvfs
c   g(i)=sqrt(pn)-1.d0  Not good for precision
c      g(i)=g(i)/(1.d0+sqrt(pn))
 19   continue

      include 'inc/TEXIT.inc'
      return
      end

      subroutine undinit(p_in,ulist)
      use tmacro
      implicit real*8 (a-h,o-z)
      parameter (nsli=1,npole=2,len=3,i_F0=4,i_G0=5,i_dph=6,
     &   i_Kz=7,i_Kx=8,i_Ky=9,i_Qx=10,i_Qy=11,lambda=12,i_ds=13)
      parameter (nulist=20)
      real*8 p_in(70),ulist(nulist)
c
c  p_in(1) length , (2) Bx     , (3) By     , (4) Kx , (5) Qy ,
c      (6) dphase , (7) slice  , (8) Pole
c
      ulist(nsli)=p_in(7)
      ulist(npole)=p_in(8)
      ulist(len)=p_in(1)
      ulist(i_F0)=p_in(2)/pgev*c
      ulist(i_G0)=p_in(3)/pgev*c
      ulist(i_dph)=p_in(6)*pi/180.
c
      ulist(i_Kz)=pi2*ulist(npole)/ulist(len)
      ulist(i_Kx)=p_in(4)
      ulist(i_Ky)=sqrt(ulist(i_Kz)*ulist(i_Kz)
     &            +ulist(i_Kx)*ulist(i_Kx))
      ulist(i_Qy)=p_in(5)
      ulist(i_Qx)=sqrt(ulist(i_Kz)*ulist(i_Kz)
     &            +ulist(i_Qy)*ulist(i_Qy))
      ulist(lambda)=ulist(len)/ulist(npole)
c
      ulist(i_ds)=ulist(len)/ulist(nsli)
c
      return
      end




