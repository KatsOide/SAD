      module bendib
      use tfstk
      implicit none
      real*8 xi,pxi,yi,pyi,zi,p,dp,rhoe,rho0,rhob,drhob,drhop,
     $     rhosq,
     $     akk,akxsq,akysq,akx,aky,dcx,aksx,dcy,aksy,phix,phiy,
     $     spx,spy,sxkx,syky,dcxkx,xsxkx
      real*8 a3,a5,a7,a9,a11,a13,a15,psqmax,epsbend
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0,
     $     psqmax=0.9999d0,epsbend=1.d-3)

      contains

      subroutine tbendiinit(ak1,al)
      implicit none
      real*8, intent (in):: ak1,al
      real*8 xsin,xsinh
      rhosq=rho0*rhoe
      drhop=(rhoe-rho0)/rhosq
      akk=ak1/al
      akxsq=-1.d0/rhosq-akk/p
      akysq=akk/p
      if(akxsq .gt. 0.d0)then
        akx=sqrt(akxsq)
        phix=akx*al
        dcx=2.d0*sinh(.5d0*phix)**2
        spx=sinh(phix)
        aksx=akx*spx
        dcxkx=dcx/akxsq
        sxkx=spx/akx
        xsxkx=xsinh(phix)/akx/akxsq
      elseif(akxsq .lt. 0.d0)then
        akx=sqrt(-akxsq)
        phix=akx*al
        dcx=-2.d0*sin(.5d0*phix)**2
        spx=sin(phix)
        aksx=-akx*spx
        dcxkx=dcx/akxsq
        sxkx=spx/akx
        xsxkx=-xsin(phix)/akx/akxsq
      else
        akx=0.d0
        phix=0.d0
        dcx=0.d0
        spx=0.d0
        aksx=0.d0
        dcxkx=0.5d0*al**2
        sxkx=al
        xsxkx=-1.d0/6.d0*al**3
      endif
      if(akysq .gt. 0.d0)then
        aky=sqrt(akysq)
        phiy=aky*al
        dcy=2.d0*sinh(.5d0*phiy)**2
        spy=sinh(phiy)
        aksy=aky*spy
        syky=spy/aky
      elseif(akysq .lt. 0.d0)then
        aky=sqrt(-akysq)
        phiy=aky*al
        dcy=-2.d0*sin(.5d0*phiy)**2
        spy=sin(phiy)
        aksy=-aky*spy
        syky=spy/aky
      else
        aky=0.d0
        phiy=0.d0
        dcy=0.d0
        spy=0.d0
        aksy=0.d0
        syky=al
      endif
      return
      end subroutine

      subroutine tbendicorr(ak1,al,phi)
      implicit none
      real*8 , intent(in) :: ak1,al,phi
      real*8 pzi,dpzi,xr,s
      s=min(psqmax,pxi**2+pyi**2)
      dpzi=sqrt1(-s)
      pzi=1.d0+dpzi
      xi =(-al*dpzi*pxi+pzi*xi)/(pzi-phi*pxi)
      pxi=pxi+phi*dpzi
      yi =yi+(-al*dpzi+xi*phi)*pyi/pzi
      zi =zi-phi*dpzi*xi/pzi+al*dpzi*(.5d0*pzi+.5d0-1.d0/pzi)
      xr=xi/rho0
      pxi=pxi-ak1*xi*xr*(0.5d0-xr*(2.d0-xr)/12.d0)/p
      return
      end subroutine

      subroutine tbendibody(al)
      implicit none
      real*8, intent(in):: al
      real*8 dxf,dpxf,dyf,dpyf,hi,xiksq
c      dxf =(drhopak+xi)*dcx+aksx/akxsq*pxi
c      dpxf= aksx*(drhopak+xi)+pxi*dcx
c      dyf = dcy*yi+aksy/akysq*pyi
c      dxf = drhop*dcxkx+xi*dcx+sxkx*pxi
      xiksq = drhop    +akxsq*xi
      dxf = dcxkx*xiksq+sxkx*pxi
      dpxf= sxkx*xiksq +dcx*pxi
      dyf = dcy*yi     +syky*pyi
      dpyf= aksy*yi    +dcy*pyi
      hi=.5d0*(pxi**2+pyi**2-akxsq*xi**2-akysq*yi**2)
     $     -drhop*xi
      zi =zi-.25d0*(dxf*pxi+xi*dpxf+dxf*dpxf
     $     +dyf*pyi+yi*dpyf+dyf*dpyf
     $     -(5.d0/rho0-1.d0/rhoe)*(drhop*xsxkx-sxkx*xi-dcxkx*pxi))
     $     -hi*al*.5d0
      xi =xi +dxf
      pxi=pxi+dpxf
      yi =yi +dyf
      pyi=pyi+dpyf
      return
      end subroutine

      end module

      subroutine tbendi(np,x,px,y,py,z,g,dv,sx,sy,sz,l,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb1,fb2,mfring,enarad,fringe,eps0)
      use tbendcom, only:tbrot
      use bendib
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      use tspin
      implicit none
      integer*4 np,mfring,i,ndiv,n,l
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np),px0(np),py0(np),zr0(np),bsi(np),
     $     al,phib,phi0,cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,cost,sint,fb1,fb2,eps0,
     $     tanp1,tanp2,aind,b,dxfr1,dyfr1,dyfra1,pr,eps,
     $     af,f,fpx,ff,akn,aln,phin,f1r,f2r,
     $     dxfr2,dyfr2,dyfra2,dtheta,cphin,sphin,sx0
      logical*4 enarad,fringe,krad
      include 'inc/TENT.inc'
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      rho0=al/phi0
      aind=rho0/phi0*ak
      krad=rad .and. enarad .and. al .ne. 0.d0
      if(krad)then
        px0=px
        py0=py
        zr0=z
        bsi=0.d0
        if(iprev(l) .eq. 0)then
          f1r=fb1
        else
          f1r=0.d0
        endif
        if(inext(l) .eq. 0)then
          f2r=fb2
        else
          f2r=0.d0
        endif
        b=brhoz/rhob
c        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
c     1             1.d0/rho0,-tanp1*2.d0/al,.5d0*al,
c     $       f1r,f2r,0.d0,al,1.d0)
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 2)then
        call ttfrin(np,x,px,y,py,z,g,4,ak,al,0.d0)
      endif
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fb1**2/rhob/24.d0
          dyfr1=fb1/rhob**2/6.d0
          dyfra1=4.d0*dyfr1/fb1**2
          do i=1,np
c            dp=g(i)*(2.d0+g(i))
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)+dxfr1*dp/pr
            py(i)=py(i)+(dyfr1-dyfra1*y(i)**2)*y(i)/pr**2
            z(i)=z(i)+(dxfr1*px(i)+
     $           (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
      if(eps0 .le. 0.d0)then
        eps=epsbend
      else
        eps=epsbend*eps0
      endif
      ndiv=1+int(abs(phi0/eps))
      if(fringe)then
        af=1.d0
      else
        af=0.d0
      endif
      akn=ak/ndiv
      aln=al/ndiv
      phin=phi0/ndiv
      cphin=cos(phin)
      sphin=sin(phin)
      do i=1,np
        dp=g(i)
        p=1.d0+dp
        rhoe=rhob*p
        yi=y(i)
        f=yi/rhoe
        fpx=af*px(i)
        pyi=py(i)-(tanp1+fpx)*f
        ff  =af*yi*f*.5d0
        xi=x(i)+ff
        zi=z(i)-ff*fpx
        pxi=px(i)+tanp1*xi/rhoe
        call tbendiinit(akn,aln)
        call tbendicorr(akn*.5d0,aln*.5d0,phin*.5d0)
        call tbendibody(aln)
        if(krad)then
          bsi(i)=bsi(i)+akn/aln*xi*yi
          if(rfluct)then
            call tradkf1(xi,pxi,yi,pyi,zi,dp,dv(i),sx(i),sy(i),sz(i),
     $           px0(i),py0(i),zr0(i),cphin,sphin,bsi(i),aln)
          else
            call tradk1(xi,pxi,yi,pyi,zi,dp,dv(i),sx(i),sy(i),sz(i),
     $           px0(i),py0(i),zr0(i),cphin,sphin,bsi(i),aln)
          endif
          px0(i)=pxi
          py0(i)=pyi
          zr0(i)=zi
          bsi(i)=0.d0
        endif
        do n=2,ndiv
          call tbendicorr(akn,aln,phin)
          call tbendibody(aln)
          if(krad .and. n .ne. ndiv)then
            if(rfluct)then
              call tradkf1(xi,pxi,yi,pyi,zi,dp,dv(i),sx(i),sy(i),sz(i),
     $             px0(i),py0(i),zr0(i),cphin,sphin,bsi(i),aln)
            else
              call tradk1(xi,pxi,yi,pyi,zi,dp,dv(i),sx(i),sy(i),sz(i),
     $             px0(i),py0(i),zr0(i),cphin,sphin,bsi(i),aln)
            endif
            px0(i)=pxi
            py0(i)=pyi
            zr0(i)=zi
          endif
        enddo
        call tbendicorr(akn*.5d0,aln*.5d0,phin*.5d0)
        zi=zi-dv(i)*al
        px(i)=pxi+tanp2*xi/rhoe
        y(i)=yi
        f=yi/rhoe
        fpx=af*px(i)
        py(i)=pyi-(tanp2-fpx)*f
        ff  =af*yi*f*.5d0
        x(i)=xi-ff
        z(i)=zi+ff*fpx
      enddo
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          dyfra2=4.d0*dyfr2/fb2**2
          do i=1,np
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)-dxfr2*dp/pr
            py(i)=py(i)+(dyfr2-dyfra2*y(i)**2)*y(i)/pr**2
            z(i)=z(i)-(dxfr2*px(i)-
     $           (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
c        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
c     1             1.d0/rho0,-tanp2*2.d0/al,.5d0*al,
c     $       f1r,f2r,al,al,-1.d0)
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 1)then
        call ttfrin(np,x,px,y,py,z,g,4,-ak,al,0.d0)
      endif
      if(krad)then
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       px0,py0,zr0,cphin,sphin,bsi,aln)
      endif
c      if(dphiy .ne. 0.d0)then
c        do i=1,np
c          pr=(1.d0+g(i))**2
c          pr=1.d0+g(i)
c          px(i)=px(i)+dphix/pr
c          py(i)=py(i)+dphiy/pr
c        enddo
c      endif
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      include 'inc/TEXIT.inc'
      return
      end
