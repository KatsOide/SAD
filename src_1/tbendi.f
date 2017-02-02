      module bendib
      use tfstk
      implicit none
      real*8 xi,pxi,yi,pyi,zi,p,dp,rhoe,rho0,rhob,drhob,drhop,
     $     drhopak,rhosq,
     $     akk,akxsq,akysq,akx,aky,dcx,aksx,dcy,aksy,phix,phiy
      real*8 a3,a5,a7,a9,a11,a13,a15,psqmax,epsbend
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0,
     $     psqmax=0.9999d0,epsbend=1.d-3)

      contains

      subroutine tbendiinit(ak1,al)
      implicit none
      real*8, intent (in):: ak1,al
      rhosq=rho0*rhoe
      drhop=(rhoe-rho0)/rhosq
      akk=ak1/al
      akxsq=(-p/rhosq-akk)/p
      drhopak=drhop/akxsq
      akysq=akk/p
      if(akxsq .ge. 0.d0)then
        akx=sqrt(akxsq)
        phix=akx*al
        dcx=2.d0*sinh(.5d0*phix)**2
        aksx=akx*sinh(phix)
      else
        akx=sqrt(-akxsq)
        phix=akx*al
        dcx=-2.d0*sin(.5d0*phix)**2
        aksx=-akx*sin(phix)
      endif
      if(akysq .ge. 0.d0)then
        aky=sqrt(akysq)
        phiy=aky*al
        dcy=2.d0*sinh(.5d0*phiy)**2
        aksy=aky*sinh(phiy)
      else
        aky=sqrt(-akysq)
        phiy=aky*al
        dcy=-2.d0*sin(.5d0*phiy)**2
        aksy=-aky*sin(phiy)
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
      real*8 dxf,dpxf,dyf,dpyf,hi
      dxf =(drhopak+xi)*dcx+aksx/akxsq*pxi
      dpxf= aksx*(drhopak+xi)+pxi*dcx
      dyf = dcy*yi+aksy/akysq*pyi
      dpyf=aksy*yi+dcy*pyi
      hi=.5d0*(pxi**2+pyi**2-akxsq*xi**2-akysq*yi**2)
     $     -drhop*xi
      zi =zi-.25d0*(dxf*pxi+xi*dpxf+dxf*dpxf
     $     +dyf*pyi+yi*dpyf+dyf*dpyf
     $     -(5.d0/rho0-1.d0/rhoe)*(drhopak*al-dpxf/akxsq))
     $     -hi*al*.5d0
      xi =xi +dxf
      pxi=pxi+dpxf
      yi =yi +dyf
      pyi=pyi+dpyf
      return
      end subroutine

      end module

      subroutine tbendi(np,x,px,y,py,z,g,dv,pz,l,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1     fb1,fb2,mfring,enarad,fringe,eps0)
      use bendib
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 np,mfring,i,ndiv,n,l
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np),
     $     al,phib,phi0,cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dphix,dphiy,cost,sint,fb1,fb2,eps0,
     $     tanp1,tanp2,aind,b,dxfr1,dyfr1,dyfra1,pr,eps,
     $     af,f,fpx,ff,akn,aln,phin,f1r,f2r,
     $     dxfr2,dyfr2,dyfra2
      logical*4 enarad,fringe
      include 'inc/TENT.inc'
      if(dphiy .ne. 0.d0)then
        do i=1,np
c          pr=(1.d0+g(i))**2
          pr=1.d0+g(i)
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
        enddo
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      rho0=al/phi0
      aind=rho0/phi0*ak
      if(rad .and. enarad)then
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
        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
     1             1.d0/rho0,-tanp1*2.d0/al,.5d0*al,
     $       f1r,f2r,0.d0,al,1.d0)
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
        do n=2,ndiv
          call tbendicorr(akn,aln,phin)
          call tbendibody(aln)
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
      if(rad .and. enarad)then
        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
     1             1.d0/rho0,-tanp2*2.d0/al,.5d0*al,
     $       f1r,f2r,al,al,-1.d0)
      endif
      if(dphiy .ne. 0.d0)then
        do i=1,np
c          pr=(1.d0+g(i))**2
          pr=1.d0+g(i)
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
        enddo
      endif
      include 'inc/TEXIT.inc'
      return
      end
