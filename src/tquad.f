      subroutine tquad(np,x,px,y,py,z,g,dv,pz,l,al,ak,
     1                 dx,dy,theta,cost,sint,radlvl,chro,
     1                 fringe,f1in,f2in,f1out,f2out,mfring,eps0,kin)
      use ffs_flag
      use tmacro
c      use ffs_pointer, only:inext,iprev
      use tfstk, only:pxy2dpz,sqrt1
      implicit none
      logical*4 enarad,chro,fringe,kin
      integer*4 np,l,i,mfring
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np),
     $     px0(np),py0(np),
     $     al,ak,dx,dy,theta,cost,sint,radlvl,eps0,alr,
     $     f1in,f1out,f2in,f2out,p,a,ea,b,pxi,pxf,pyf,xi
      real*8, parameter :: ampmax=0.9999d0
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,pz,4,l,0.d0,ak,
     $             dx,dy,theta,cost,sint, 1.d0,.false.)
        return
      elseif(ak .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        return
      endif
      enarad=rad .and. radlvl .ne. 1.d0
      if(trpt .and. enarad)then
        call tqrad(np,x,px,y,py,z,g,dv,pz,l,al,ak,dx,dy,theta,
     1             cost,sint,radlvl,f1in,f2in,f1out,f2out,mfring)
        return
      endif
      include 'inc/TENT.inc'
      if(enarad)then
        px0=px
        py0=py
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 2)then
        call ttfrin(np,x,px,y,py,z,g,4,ak,al,0.d0)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        do 2110 i=1,np
c          p=(1.d0+g(i))**2
          p=1.d0+g(i)
          a=f1in/p
          ea=exp(a)
          b=f2in/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
2110    continue
      endif
c$$$      if(enarad)then
c$$$        if(iprev(l) .eq. 0)then
c$$$          f1r=sqrt(abs(24.d0*f1in/ak*al))
c$$$        else
c$$$          f1r=0.d0
c$$$        endif
c$$$        if(inext(l) .eq. 0)then
c$$$          f2r=sqrt(abs(24.d0*f1out/ak*al))
c$$$        else
c$$$          f2r=0.d0
c$$$        endif
c$$$        b1=brho*ak/al
c$$$        call trad(np,x,px,y,py,g,dv,0.d0,0.d0,
c$$$     1       b1,0.d0,0.d0,.5d0*al,
c$$$     $       f1r,f2r,0.d0,al,1.d0)
c$$$      endif
      if(enarad)then
        call tsolqur(np,x,px,y,py,z,g,dv,pz,al,ak,
     $       0.d0,0.d0,0.d0,eps0,px0,py0,alr)
      else
        call tsolqu(np,x,px,y,py,z,g,dv,pz,al,ak,0.d0,0.d0,0.d0,eps0)
      endif
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        do 2120 i=1,np
c          p=(1.d0+g(i))**2
          p=1.d0+g(i)
          a=-f1out/p
          ea=exp(a)
          b= f2out/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
2120    continue
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 1)then
        call ttfrin(np,x,px,y,py,z,g,4,-ak,al,0.d0)
      endif
      if(enarad)then
        call tradk(np,x,px,y,py,px0,py0,g,dv,alr)
      endif
      include 'inc/TEXIT.inc'
      return
      end
c
      subroutine tthin(np,x,px,y,py,z,g,dv,pz,nord,l,al,ak,
     1                 dx,dy,theta,cost,sint,radlvl,fringe)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
c     alpha=1/sqrt(12),beta=1/6-alpha/2,gamma=1/40-1/24/sqrt(3)
      integer*4 nmult
      parameter (nmult=21)
      real*8 alpha,beta,gamma,alpha1
      parameter (alpha=2.88675134594812882d-1,
     1           beta =2.23290993692602255d-2,
     1           gamma=9.43738783765593145d-4,
     1           alpha1=.5d0-alpha)
      integer*4 l,nord,np,kord,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),pz(np),
     $     px0(np),py0(np)
      real*8 fact(0:nmult)
      real*8 theta,sint,cost,dx,dy,al,ak,
     $     ala,alb,aki,akf,dpz,al1,radlvl,
     $     f1,f2,f3,f4,f5,xi,yi,zi,pr,r,rcx1,pxi,rk1,rk
      complex*16 cx
      logical enarad,fringe
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2432902008176640000.d0,51090942171709440000.d0/
c     begin initialize for preventing compiler warning
      ala=0.d0
      alb=0.d0
      aki=0.d0
c     end   initialize for preventing compiler warning
      if(ak .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        return
      endif
      enarad=rad .and. radlvl .eq. 0.d0 .and. al .ne. 0.d0
      if(enarad .and. trpt .and. rfluct)then
        call tthinrad(np,x,px,y,py,z,g,dv,pz,nord,l,al,ak,
     1                 dx,dy,theta,cost,sint,fringe)
        return
      endif
      include 'inc/TENT.inc'
      if(enarad)then
        px0=px
        py0=py
      endif
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,ak,al,0.d0)
      endif
      kord=nord/2-1
      akf=ak/fact(kord)
      if(al .ne. 0.d0)then
        ala=al*alpha1
        alb=al*alpha
        do 10 i=1,np
c          a=px(i)**2+py(i)**2
          dpz=pxy2dpz(px(i),py(i))
c          dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
          al1=ala/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*ala
10      continue
        akf=akf*.5d0
      endif
      if(kord .eq. 2)then
        do 1210 i=1,np
c     aki=akf/(1.d0+g(i))**2
          aki=akf/(1.d0+g(i))
          px(i)=px(i)-aki*(x(i)-y(i))*(x(i)+y(i))
          py(i)=py(i)+2.d0*aki*x(i)*y(i)
 1210   continue
      else
        do 1020 i=1,np
c     pr=(1.d0+g(i))**2
          pr=1.d0+g(i)
          aki=akf/pr
          cx=dcmplx(x(i),-y(i))**kord
          px(i)=px(i)-aki*dble(cx)
          py(i)=py(i)-aki*imag(cx)
 1020   continue
      endif
      if(al .ne. 0.d0)then
        if(kord .le. 0)then
          do 1030 i=1,np
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)*2.d0
            x(i)=x(i)+px(i)*al1
            y(i)=y(i)+py(i)*al1
            z(i)=z(i)+dpz  *al1-dv(i)*alb*2.d0
1030      continue
        elseif(kord .eq. 2)then
          f1=beta*2.d0*al
          f2=4.d0*gamma*al**2
          f3=gamma*6.d0*al**2
          f4=.5d0*beta*al
          f5=f2
          do 1031 i=1,np
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
c            pr=(1.d0+g(i))**2
            pr=1.d0+g(i)
            aki=akf/pr*2.d0
            r=xi  **2+yi  **2
            cx=dcmplx(xi  ,-yi  )**2
            rcx1=xi  *dble(cx)+yi  *imag(cx)
            px(i)=px(i)+aki**2*(f1*r*xi
     1        -aki*(f2*xi  *rcx1
     1             +f3*r*dble(cx)) )
            py(i)=py(i)+aki**2*(f1*r*yi
     1        -aki*(f2*yi  *rcx1
     1             +f3*r*imag(cx)) )
            zi  =zi  +aki**2*r*(f4*r-f5*aki*rcx1)
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
1031      continue
        elseif(kord .eq. 1)then
          f1=beta*al
          f3=gamma*2.d0*al**2
          f4=.5d0*f1
          f5=f3
          do i=1,np
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
            r=xi**2+yi**2
            rcx1=(xi-yi)*(xi+yi)
            px(i)=px(i)+aki**2*(f1-aki*f3)*xi
            py(i)=py(i)+aki**2*(f1+aki*f3)*yi
            zi  =zi  +aki**2*(f4*r-f5*aki*rcx1)
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
          enddo
        else
          f1=beta*kord*al
          f2=2.d0*gamma*kord*(kord-1)*al**2
          f3=gamma*kord*(kord+1)*al**2
          f4=.5d0*beta*al
          f5=2.d0*gamma*kord*al**2
          do 1032 i=1,np
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
c            pr=(1.d0+g(i))**2
            pr=1.d0+g(i)
            aki=akf/pr*2.d0
            r=xi  **2+yi  **2
            rk1=r**(kord-2)
            rk=r*rk1
            cx=dcmplx(xi  ,-yi  )**kord
            rcx1=xi  *dble(cx)+yi  *imag(cx)
            px(i)=px(i)+aki**2*(f1*rk*xi
     1        -aki*(f2*rk1*xi  *rcx1
     1             +f3*rk*dble(cx)) )
            py(i)=py(i)+aki**2*(f1*rk*yi
     1        -aki*(f2*rk1*yi  *rcx1
     1             +f3*rk*imag(cx)) )
            zi  =zi  +aki**2*rk*(f4*r-f5*aki*rcx1)
c            a=px(i)**2+py(i)**2
            dpz=pxy2dpz(px(i),py(i))
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
1032      continue
        endif
        if(kord .eq. 2)then
          do 1220 i=1,np
c            aki=akf/(1.d0+g(i))**2
            aki=akf/(1.d0+g(i))
            px(i)=px(i)-aki*(x(i)-y(i))*(x(i)+y(i))
            py(i)=py(i)+2.d0*aki*x(i)*y(i)
1220      continue
        else
          do 1040 i=1,np
c            aki=akf/(1.d0+g(i))**2
            aki=akf/(1.d0+g(i))
            cx=dcmplx(x(i),-y(i))**kord
            px(i)=px(i)-aki*dble(cx)
            py(i)=py(i)-aki*imag(cx)
1040      continue
        endif
        do 1050 i=1,np
c          a=px(i)**2+py(i)**2
          dpz=pxy2dpz(px(i),py(i))
c          dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
          al1=ala/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*ala
1050    continue
      endif
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,-ak,al,0.d0)
      endif
      if(enarad)then
        call tradk(np,x,px,y,py,px0,py0,g,dv,al)
      endif
      include 'inc/TEXIT.inc'
      return
      end
c
      subroutine tthinrad(np,x,px,y,py,z,g,dv,pz,nord,l,al,ak,
     1                 dx,dy,theta,cost,sint,fringe)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 nmult,n,ndivmax
      parameter (nmult=21)
      real*8 ampmax,eps00
      parameter (ampmax=0.05d0,eps00=0.005d0,ndivmax=1000)
      integer*4 l,np,i,kord,ndiv,nord
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),pz(np)
      real*8 fact(0:nmult)
      real*8 dx,dy,theta,cost,sint,an,sp,akfl,akf,alsum,aln,
     $     al0,pr,thr,alx,al1,aki,dpx,dpy,alr,prob,bxa,bya,
     $     delp,xi,pxi,ak,al,h1,dpz
      real*8 dprad,dpradx,dprady
      real*8 tran
      complex*16 cx
      logical fringe
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2432902008176640000.d0,51090942171709440000.d0/
      include 'inc/TENT.inc'
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,ak,al,0.d0)
      endif
      n=nord/2
      kord=n-1
      akf=ak/fact(kord)
      ndiv=min(ndivmax,max(2,int(
     $     sqrt(ampmax**(n-1)/6.d0/eps00*abs(akf)*al))+1))
      an=anrad*p0
      sp=0.d0
      akfl=akf/al
      do i=1,np
        alsum=0.d0
        aln=al/ndiv
        al0=0.d0
c        delp=g(i)*(2.d0+g(i))
        delp=g(i)
        pr=1.d0+delp
 100    cx=dcmplx(x(i),-y(i))**kord
        thr=abs(akfl)*abs(cx)
        if(thr .eq. 0.d0)then
          alx=aln
        else
          alx=min(aln,0.08d0/an/thr)
        endif
        al1=al0+alx*.5d0
        aki=akfl*al1/pr
        dpx=-aki*dble(cx)
        dpy=-aki*imag(cx)
        px(i)=px(i)+dpx
        py(i)=py(i)+dpy
        alr=al1*(1.d0+(px(i)**2+py(i)**2)*.5d0)
        thr=thr*alr
        if(thr .ne. 0.d0)then
          prob=an*thr
          if((tran()-.5d0)*3.46410161513775461d0 .gt. .5d0-prob)then
            bya=-dpx*brhoz/al1
            bxa= dpy*brhoz/al1
            call tsynchrad(pr,alr,bxa,bya,
     $           dprad,dpradx,dprady,
     $           i,l,alsum,-al0,0.d0,theta,
     $           x(i),y(i),px(i),py(i))
            px(i)=px(i)-dpradx
            py(i)=py(i)-dprady
          else
            dprad=0.d0
          endif
        else
          dprad=0.d0
        endif
        delp=delp-dprad
        pr=1.d0+delp
        sp=sp-dprad
        if(alx .ne. 0.d0)then
c          a=px(i)**2+py(i)**2
          dpz=pxy2dpz(px(i),py(i))
c          dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c          dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
          al1=alx/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*alx
        endif
        al0=alx*.5d0
        alsum=alsum+alx
        if(aln .gt. 0)then
          aln=min(al/ndiv,max(al-alsum,0.d0))
          go to 100
        endif
c        h1=sqrt(1.d0+(p0*pr)**2)
        h1=p2h(p0*pr)
c        g(i)=delp/(1.d0+sqrt(max(.01d0,pr)))
        g(i)=delp
        dv(i)=-delp*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
      enddo
      if(.not. radcod  )then
        sp=sp/np
        do 1110 i=1,np
c          dp=(2.d0+g(i))*g(i)-sp
c          g(i)=dp/(1.d0+sqrt(1.d0+dp))
          g(i)=g(i)-sp
c     Change of dv(i) is ignored here.
 1110   continue
      endif
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,-ak,al,0.d0)
      endif
      include 'inc/TEXIT.inc'
      return
      end
c
      subroutine ttfrin(np,x,px,y,py,z,g,nord,ak,al,bz)
      implicit none
      real*8 xlimit,plimit
      parameter (xlimit=10.d0,plimit=0.99d0)
      integer*4 np,nord,i,kord
      real*8 x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ak,al,akk,aki,a,b,ab,t,dx1,dy1,d,xx,yy,
     $     h,f,px1,py1,xi,yi,bz,bzph,pr,pxa,pya
      complex*16 cx,cp,cz,cz1,cx1,cp1,ca
      integer*4 nmult
      parameter (nmult=22)
      real*8 fact(0:nmult),an(nmult+1)
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2.432902008176640000d18,5.1090942171709440000d19,
     $     1.124000727777607680000d21/
      data an /1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,
     $     11.d0,12.d0,13.d0,14.d0,15.d0,16.d0,17.d0,18.d0,19.d0,
     $     20.d0,21.d0,22.d0,23.d0/
      if(al .eq. 0.d0 .or. ak .eq. 0.d0)then
        return
      endif
      if(nord .eq. 4)then
        akk=ak/al/4.d0
        if(bz .eq. 0.d0)then
          do i=1,np
c            aki=akk/(1.d0+g(i))**2
            aki=akk/(1.d0+g(i))
            xi=x(i)
            yi=y(i)
            a=aki*xi**2
            b=aki*yi**2
            ab=aki*(xi-yi)*(xi+yi)
            t=ab**2/6.d0
            dx1= xi*(a/3.d0+b+t)
            dy1=-yi*(a+b/3.d0-t)
            x(i)=xi+dx1
            y(i)=yi+dy1
            d=a+b
            xx=1.d0+d+ab*(5.d0*a-b)/6.d0
            yy=1.d0-d+ab*(a-5.d0*b)/6.d0
            h=2.d0*aki*xi*yi*(1.d0-ab/3.d0)
            f=xx*yy+h**2
            px1=(px(i)*yy+py(i)*h)/f
            py1=(py(i)*xx-px(i)*h)/f
            z(i)=z(i)-px1*(dx1+t*xi)-py1*(dy1+t*yi)
            px(i)=px1
            py(i)=py1
          enddo
        else
          do i=1,np
c            pr=(1.d0+g(i))**2
            pr=(1.d0+g(i))
            aki=akk/pr
            bzph=.5d0*bz/pr
            xi=x(i)
            yi=y(i)
            a=aki*xi**2
            b=aki*yi**2
            ab=aki*(xi-yi)*(xi+yi)
            t=ab**2/6.d0
            dx1= xi*(a/3.d0+b+t)
            dy1=-yi*(a+b/3.d0-t)
            x(i)=xi+dx1
            y(i)=yi+dy1
            d=a+b
            xx=1.d0+d+ab*(5.d0*a-b)/6.d0
            yy=1.d0-d+ab*(a-5.d0*b)/6.d0
            h=2.d0*aki*xi*yi*(1.d0-ab/3.d0)
            f=xx*yy+h**2
            pxa=px(i)+bzph*y(i)
            pya=py(i)-bzph*x(i)
            px(i)=(pxa*yy+pya*h)/f-bzph*yi
            py(i)=(pya*xx-pxa*h)/f+bzph*xi
            z(i)=z(i)-px(i)*(dx1+t*xi)-py(i)*(dy1+t*yi)
          enddo
        endif
      elseif(nord .eq. 2)then
        akk=ak/al
        do i=1,np
          aki=akk/(1.d0+g(i))
          x(i)=x(i)+.5d0*aki*y(i)**2
          py(i)=py(i)-aki*px(i)*y(i)
          z(i)=z(i)-.5d0*aki*px(i)*y(i)**2
        enddo
      elseif(nord .eq. 6)then
c        akk=ak/al/12.d0
        akk=ak/al/24.d0
        do i=1,np
c          aki=akk/(1.d0+g(i))**2
          aki=akk/(1.d0+g(i))
          cx=dcmplx(x(i),y(i))
          cp=dcmplx(px(i),-py(i))
          cz1=cx*cx
          cz=cz1*cx
          a=aki*imag(cz)*2.d0
          ca=-aki*cx*(cz/4.d0-conjg(cz))
          cx1=cx+ca
c          d=1.d0+a**2-(aki*4.d0)**2*
          d=1.d0+a**2-(aki*3.d0)**2*
     1         (dble(cz)**2+imag(cz)**2)
c          cp1=(dcmplx(1.d0,a)*cp-(aki*4.d0)*cz1*conjg(cx*cp))/d
          cp1=(dcmplx(1.d0,a)*cp-(aki*3.d0)*cz1*conjg(cx*cp))/d
          z(i)=z(i)-dble(ca*cp1)
          x(i)=dble(cx1)
          y(i)=imag(cx1)
          px(i)=dble(cp1)
          py(i)=-imag(cp1)
        enddo
      else
        kord=nord/2
c        akk=ak/al/fact(kord)/2.d0
        akk=ak/al/fact(kord)/4.d0
        do i=1,np
c          aki=akk/(1.d0+g(i))**2
          aki=akk/(1.d0+g(i))
          cx=dcmplx(x(i),y(i))
          cp=dcmplx(px(i),-py(i))
          cz1=cx**(kord-1)
          cz=cz1*cx
          a=aki*imag(cz)*2.d0
          ca=-aki*cx*(cz/an(kord+1)-conjg(cz))
c          write(*,*)'ttfrin ',kord,ca,an(kord+1),cz
          cx1=cx+ca
          d=1.d0+a**2-(aki*an(kord))**2*
     1         (dble(cz)**2+imag(cz)**2)
          cp1=(dcmplx(1.d0,a)*cp-(aki*an(kord))*cz1*conjg(cx*cp))/d
          z(i)=z(i)-dble(ca*cp1)
          x(i)=dble(cx1)
          y(i)=imag(cx1)
          px(i)=dble(cp1)
          py(i)=-imag(cp1)
        enddo
      endif
      return
      end

      subroutine ttfrins(np,x,px,y,py,z,g,nord,ak,al,bz)
      use macmath
      implicit none
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 bz,aka
      integer*4 nord,i
      real*8 ak(2),al,x0,px0,theta,cost,sint
      if((ak(1) .ne. 0.d0 .or. ak(2) .ne. 0.d0) .and. al .ne. 0.d0)then
c        theta=pi/nord
        if(ak(1) .eq. 0.d0)then
          theta=pi/nord
          aka=ak(2)
        elseif(ak(2) .eq. 0.d0)then
          theta=0.d0
          aka=ak(1)
        else
          theta=atan2(ak(2),ak(1))*2.d0/nord
          aka=sqrt(ak(1)**2+ak(2)**2)
        endif
        if(theta .ne. 0.d0)then
          cost=cos(theta)
          sint=sin(theta)
          do i=1,np
            x0   =x(i)
            x(i) = cost*x0 -sint*y(i)
            y(i) = sint*x0 +cost*y(i)
            px0  =px(i)
            px(i)= cost*px0-sint*py(i)
            py(i)= sint*px0+cost*py(i)
          enddo
          call ttfrin(np,x,px,y,py,z,g,nord,aka,al,bz)
          do i=1,np
            x0   =x(i)
            x(i) = cost*x0 +sint*y(i)
            y(i) =-sint*x0 +cost*y(i)
            px0  =px(i)
            px(i)= cost*px0+sint*py(i)
            py(i)=-sint*px0+cost*py(i)
          enddo
        else
          call ttfrin(np,x,px,y,py,z,g,nord,aka,al,bz)
        endif
      endif
      return
      end

      subroutine tqlfri(np,x,px,y,py,z,g,f1,f2,bz)
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 bz,f1,f2,xf,yf,pxf,pyf,a,b,bb,bzph,ea,f,p
      do i=1,np
c        p=(1.d0+g(i))**2
        p=(1.d0+g(i))
        bzph=.5d0*bz/p
        a=f1/p
        ea=exp(a)
        b=f2/p
        bb=b*bzph/p
        f=1.d0/((1.d0-bb)*(1.d0+bb))
        xf=(ea*x(i)+b*(px(i)+bzph*y(i)/ea-bb*py(i)))*f
        yf=(y(i)/ea-b*(py(i)-bzph*x(i)*ea-bb*px(i)))*f
        pxf=(px(i)+bzph*yf)/ea
        pyf=(py(i)-bzph*xf)*ea
        z(i)=z(i)-(a*ea*x(i)+b*ea*(1.d0+.5d0*a)*pxf)*pxf
     $       +(a*y(i)/ea+b/ea*(1.d0-.5d0*a)*pyf)*pyf
        px(i)=pxf-bzph*y(i)
        py(i)=pyf+bzph*x(i)
        x(i)=xf
        y(i)=yf
      enddo
      return
      end
