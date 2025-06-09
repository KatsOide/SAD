      module bendeb
      use tfstk
      use bendib
      use temw,only:bsir0,tmulbs
      use tspin
      implicit none
      logical*4 ,save :: tbinit
      real*8 ,save:: al,b1,b,aind,trans1(6,6),cx,cy,
     $     drhopp,akxsqp,dcxp,sxkxp,dcyp,sykyp,phixp,phiyp,
     $     akysqp,dcxkxp,aksxp,aksyp,xsxkxp

      contains
      subroutine tbendeinit(ak1,al)
      implicit none
      real*8 , intent(in) :: ak1,al
      p=1.d0+dp
      rhoe=rhob*p
      call tbendiinit(ak1,al,.true.)
      cx=1.d0+dcx
      cy=1.d0+dcy
      drhopp=1.d0/rhoe/p
      akxsqp=-akxsq/p
      phixp=-0.5d0*phix/p
      if(akxsq > 0.d0)then
        dcxp=spx*phixp
        dcxkxp=dcxp/akxsq+dcxkx/p
        sxkxp= 0.5d0*sxkx/p+cx*phixp/akx
        xsxkxp=(-dcx*phixp/akxsq-1.5d0*xsxkx/p)/akx
        aksxp=-0.5d0*aksx/p+akx*cx*phixp
      elseif(akxsq < 0.d0)then
        dcxp=-spx*phixp
        dcxkxp=dcxp/akxsq+dcxkx/p
        sxkxp= 0.5d0*sxkx/p+cx*phixp/akx
        xsxkxp=(-dcx*phixp/akxsq-1.5d0*xsxkx/p)/akx
        aksxp=-0.5d0*aksx/p-akx*cx*phixp
      else
        dcxp=0.d0
        dcxkxp=0.d0
        sxkxp=0.d0
        xsxkxp=0.d0
        aksxp=0.d0
      endif
      akysqp=-akysq/p
      phiyp=-0.5d0*phiy/p
      if(akysq > 0.d0)then
        dcyp=phiyp*spy
        sykyp= 0.5d0*syky/p+cy*phiyp/aky
        aksyp=-0.5d0*aksy/p+aky*cy*phiyp
      elseif(akysq < 0.d0)then
        dcyp=-phiyp*spy
        sykyp= 0.5d0*syky/p+cy*phiyp/aky
        aksyp=-0.5d0*aksy/p-aky*cy*phiyp
      else
        dcyp=0.d0
        sykyp=0.d0
      endif
      tbinit=.false.
      return
      end subroutine

      subroutine tbendecorr(trans,cod,beam,ak1,al)
      use tmacro
      use sad_basics
      use mathfun
      use sad_basics
      implicit none
      real*8 , intent(in) :: ak1,al
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42)
      real*8 pzi,dpzi,xr,s,a,trans2(6,6),dk,pxf,xf,yf,zf,phi0
      phi0=al/rho0
      p  =1.d0+cod(6)
      xi =cod(1)
      pxi=min(p,max(-p,cod(2)))/p
      yi =cod(3)
      pyi=min(p,max(-p,cod(4)))/p
      s=min(psqmax,pxi**2+pyi**2)
      dpzi=sqrt1(-s)
      pzi=1.d0+dpzi
      a=1.d0/(pzi-phi0*pxi)
      xf =(-al*dpzi*pxi+pzi*xi)*a
      pxf=pxi+phi0*dpzi
      yf =yi+(-al*dpzi+xf*phi0)*pyi/pzi
c      zf =cod(5)-phi0*dpzi*xf/pzi+al*dpzi*(.5d0*pzi+.5d0-1.d0/pzi)
      zf=cod(5)+(.5d0*s+(1.d0+xi/rho0)*dpzi/pzi)*al
      xr=xf/rho0
      pxf=pxf-ak1*xf*xr*(0.5d0-xr*(2.d0-xr)/12.d0)/p
      dk=-ak1*xr*(1.d0-xr*(1.d0-2.d0*xr))
      call tinitr(trans2)
      trans2(1,1)=pzi*a
      trans2(1,2)=((al+phi0*xi)*(pzi**2+pxi**2)
     $     -al*(phi0*pxi**3+pzi**3))*a**2/pzi/p
      trans2(1,4)=(al*(1.d0-phi0*pxi)+phi0*xi)*pxi*pyi*a**2/pzi/p
      trans2(1,6)=pxi/pzi*a*(al*(-a*(1.d0-phi0*pxi)+pzi**2)-a*phi0*xi)
      trans2(2,1)=dk*trans2(1,1)
      trans2(2,2)=1.d0-phi0*pxi/pzi+dk*trans2(1,2)
      trans2(2,4)=-phi0*pyi/pzi+dk*trans2(1,4)
      trans2(2,6)=-phi0/pzi*dpzi+dk*trans2(1,6)
      trans2(3,1)=phi0*pyi*a
      trans2(3,2)=(-al*phi0*(pzi**2+pxi**2)
     $     +(al+phi0*xi)*(pxi+phi0*pzi))*a**2*pyi/pzi/p
      trans2(3,4)=(-al*(pzi**2*(pzi-phi0*pxi)+phi0*pxi*pyi**2)
     $     +(pyi**2+pzi*(pzi-phi0*pxi))*(al+phi0*xi))*a**2/p/pzi
      trans2(3,6)=-(al*(dpzi*(1.d0+pzi)*phi0*pxi-pzi**3)+(al+phi0*xi))
     $     *a**2*pyi/pzi/p
      trans2(5,1)=-trans2(1,1)*trans2(2,6)
      trans2(5,2)=-trans2(1,2)*trans2(2,6)+trans2(2,2)*trans2(1,6)
      trans2(5,4)=trans2(3,6)
c      zf =cod(5)-phi0*dpzi*xf/pzi+al*dpzi*(.5d0*pzi+.5d0-1.d0/pzi)
c      trans2(5,6)=(-al*pzi*(a*phi0*pxi+pzi)+
c     $     a*(al*(-phi0*pxi*dpzi*(phi0*pxi-1.d0)*a/pzi**2+1.d0)+
c     $     phi0*xi*(pxi*phi0*dpzi*a/pzi**2+1.d0)))/p
c      zf=cod(5)+(.5d0*s+(1.d0+xi/rho0)*dpzi/pzi)*al
      trans2(5,6)=(-s/p+(1.d0+xi/rho0)*s/pzi**3)*al
      call tmultr5(trans,trans2,irad)
      call tmulbs(beam ,trans2,.true.)
      cod(1)=xf
      cod(2)=min(p,max(-p,pxf*p))
      cod(3)=yf
      cod(5)=zf
      return
      end subroutine

      subroutine tbendebody(trans,cod,beam,srot,al,phin,
     $     ak1,alr,bsi1,bsi2,enarad)
      use tmacro
      use kradlib, only:tradke
      use temw,only:tmulbs
      use sad_basics
      use mathfun, only:pxy2dpz
      implicit none
      real*8, intent(in):: ak1,al,bsi1,bsi2,alr
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 pxf,xf,phin,dxpx,
     $     dxf,dpxf,dyf,dpyf,hi,hip,pyf,yf,zf,xiksq
      logical*4 ,intent(in):: enarad
      dp=cod(6)
      p=1.d0+dp
      xi=cod(1)
      pxi=min(p,max(-p,cod(2)))/p
      yi=cod(3)
      pyi=min(p,max(-p,cod(4)))/p
      zi=cod(5)
      if(tbinit)then
        call tbendeinit(ak1,al)
      endif
      xiksq=drhop+akxsq*xi
      dxf = dcxkx*xiksq+sxkx*pxi
      dpxf= sxkx*xiksq +dcx*pxi
      dyf = dcy*yi     +syky*pyi
      dpyf= aksy*yi    +dcy*pyi
      hi=.5d0*(pxi**2+pyi**2-akxsq*xi**2-akysq*yi**2)
     $     -drhop*xi
      hip=.5d0*(-2.d0*(pxi**2+pyi**2)/p-akxsqp*xi**2-akysqp*yi**2)
     $     -drhopp*xi
      dxpx=dxf*pxi+xi*dpxf+dxf*dpxf+dyf*pyi+yi*dpyf+dyf*dpyf
      zf =zi-.25d0*(dxpx
c     $     -(5.d0/rho0-1.d0/rhoe)*(drhopak*al-dpxf/akxsq))
     $     -(5.d0/rho0-1.d0/rhoe)*(drhop*xsxkx-sxkx*xi-dcxkx*pxi))
     $     -hi*al*.5d0
      xf =xi +dxf
      pxf=pxi+dpxf
      yf =yi +dyf
      pyf=pyi+dpyf
      trans1(1,1)=cx
      trans1(1,2)=sxkx/p
c       dcxkx*(drhop+akxsq*xi)+sxkx*pxi
      trans1(1,6)=drhopp*dcxkx+drhop*dcxkxp
     $     +xi*dcxp+(sxkxp-sxkx/p)*pxi
c      trans1(1,6)=drhopp*dcxkx+(drhop+akxsqp*xi)*dcxkxp
c     $     +(sxkxp-0.5d0*sxkx/p)*pxi
      trans1(2,1)=aksx*p
      trans1(2,2)=cx
c      trans1(2,6)=aksxp*(drhopak+xi)*p
c     $     +aksx*(drhopakp*p+drhopak+xi)+pxi*dcxp
      trans1(2,6)=p*(drhopp*sxkx+drhop*sxkxp)+drhop*sxkx
     $     +(p*aksxp+aksx)*xi+dcxp*pxi
      trans1(3,3)=cy
      trans1(3,4)=syky/p
c      trans1(3,6)=yi*dcyp+(sykyp-0.5d0*syky/p)*pyi
      trans1(3,6)=yi*dcyp+(sykyp-syky/p)*pyi
      trans1(4,3)=aksy*p
      trans1(4,4)=cy
      trans1(4,6)=(aksyp*p+aksy)*yi+pyi*dcyp
      trans1(5,1)=-trans1(1,1)*trans1(2,6)+trans1(2,1)*trans1(1,6)
      trans1(5,2)=-trans1(1,2)*trans1(2,6)+trans1(2,2)*trans1(1,6)
      trans1(5,3)=-trans1(3,3)*trans1(4,6)+trans1(4,3)*trans1(3,6)
      trans1(5,4)=-trans1(3,4)*trans1(4,6)+trans1(4,4)*trans1(3,6)
c      dpzi=p*pxy2dpz(pxi,pyi)
c      pzi=p+dpzi
c      az=1.d0/(pzi-p*phin*pxi)
c      trans1(5,6)=(-al*pzi*(az*p**2*phin*pxi + pzi)/p**3 +
c     $     az*(al*(-phin*p*pxi*dpzi*(phin*p*pxi - p)*az/pzi**2 + 1.d0) +
c     $     phin*xi*(p**2*pxi*phin*dpzi*az/pzi**2 + 1.d0)))
c     $     +hip*al+h0/h1emit**3*al
c      zf =zi-.25d0*(dxpx
c     $     -(5.d0/rho0-1.d0/rhoe)*(drhop*xsxkx-sxkx*xi-dcxkx*pxi))
c     $     -hi*al*.5d0
c      dxpx=dxf*pxi+xi*dpxf+dxf*dpxf+dyf*pyi+yi*dpyf+dyf*dpyf
      trans1(5,6)=-0.25d0*(-dxpx/p
     $     +trans1(1,6)*pxf+trans1(2,6)/p*xf
     $     +trans1(3,6)*pyf+trans1(4,6)/p*yf
     $     -1.d0/rhoe/p*(drhop*xsxkx-sxkx*xi-dcxkx*pxi)
     $     -(5.d0/rho0-1.d0/rhoe)
c     $     *(drhopp*xsxkx+drhop*xsxkxp-sxkxp*xi-dcxkxp*pxi))
     $     *(drhopp*xsxkx+drhop*xsxkxp-sxkxp*xi-dcxkxp*pxi+dcxkx*pxi/p))
     $     -hip*al*.5d0+h0/h1emit**3*al
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.)
      cod(1)=xf
      cod(2)=min(p,max(-p,pxf*p))
      cod(3)=yf
      cod(4)=min(p,max(-p,pyf*p))
      cod(5)=zf-dvemit*al
      if(enarad)then
        bsir0=bsir0+bsi1*(yi/rhob+ak1*xi*yi)
     $       -bsi2*(yf/rhob+ak1*xf*yf)
        call tradke(trans,cod,beam,srot,alr,phin,0.d0)
        tbinit=.true.
      endif
      return
      end subroutine

      subroutine tbendebody0(trans,cod,beam,srot,aln,
     $     phi0n,sn,xsn,cn,dcn,alr,
     $     bsi1,bsi2,enarad)
      use tmacro
      use kradlib, only:tradke
      use temw,only:tmulbs
      use sad_basics
      use mathfun
      implicit none
      real*8, intent(in):: phi0n,sn,xsn,cn,dcn,aln,bsi1,bsi2,alr
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 xi,pxi,yi,pyi,dp,pxf,pza,pzf,pzi,s1,u,xf,yf,zf,
     $     dpxf,dpzf,dpzi,omega,phsq,da,dm22,
     $     dodp,dodpx,dodpy,dodx,
     $     dpzadp,dpzadpx,dpzadpy,dpzadx,dpza,ddodp,ddpzadp
      logical*4 , intent(in) ::enarad
      xi=cod(1)
      dp=cod(6)
      p=1.d0+dp
      pxi=min(p,max(-p,cod(2)))
      yi=cod(3)
      pyi=min(p,max(-p,cod(4)))
      dpzi=p*pxy2dpz(pxi/p,pyi/p)
      pzi=p+dpzi
      dpza=dpzi+(drhob-xi)/rhob
      pza=dpza+dp
c         pzi-p+p-1+(rhob-rho-xi)/rhob
      dpxf=dcn*pxi+sn*pza
      pxf=pxi+dpxf
      dpzf=p*pxy2dpz(pxf/p,pyi/p)
      pzf=p+dpzf
      s1=pxf/pzf
c      xf=rhob*(sn*pxi-cn*pza+p+dpzf-rho/rhob)
c      xf=rhob*(sn*pxi-cn*pza+p+dpzf-1+drhob/rhob)
c      xf=rhob*(sn*pxi-cn*pza+dp+dpzf+drhob/rhob)
c      xf=rhob*(sn*pxi-cn*dpza-dcn*dp+dpzf+drhob/rhob)
c      xf=rhob*(sn*pxi-cn*dpzi-cn*(drhob-xi)/rhob-dcn*dp+dpzf+drhob/rhob)
c      xf=rhob*(sn*pxi-cn*dpzi+cn*xi/rhob-dcn*dp+dpzf-dcn*drhob/rhob)
      xf=rhob*(sn*pxi-cn*dpzi-dcn*dp+dpzf)+cn*xi-dcn*drhob
c      if(abs(rhob) > 1.d10)then
c        write(*,'(a,1p10g12.4)')'tbendb0 ',xf,sn*pxi*rhob,-cn*dpza*rhob,-rhob*dcn*dp,rhob*dpzf,drhob
c      endif
      phsq=pzf**2+pxf**2
      da=asin((-p*dpxf-pxf*dpzi+pxi*dpzf)/phsq)
      omega=phi0n+da
      yf=yi+rhob*omega*pyi
      zf=cod(5)-da*p*rhob-phi0n*(drhob+dp*rhob)
      dpzadx =-1.d0/rhob
      dpzadpx=-pxi/pzi
      dpzadpy=-pyi/pzi
      ddpzadp=dpzi/pzi
c      dpzadp = 1.d0+ddpzadp
      dpzadp = p/pzi
      u=(s1+dpzadpx)/phsq
      trans1(2,1)=sn*dpzadx
      dm22=dcn+sn*dpzadpx
      trans1(2,2)=1.d0+dm22
      trans1(2,4)=sn*dpzadpy
      trans1(2,6)=sn*dpzadp
      dodx =-trans1(2,1)/pzf
c      dodpx=1.d0/pzi-trans1(2,2)/pzf
      dodpx=((dpzf-dpzi)/pzi-dm22)/pzf
      dodpy=-u*pyi-trans1(2,4)/pzf
      ddodp = u*p -sn*ddpzadp/pzf
      dodp  = ddodp-sn/pzf
c      dodp=u*p-trans1(2,6)/pzf
      trans1(1,1)=cn+sn*s1
      trans1(1,2)=rhob*(sn-cn*dpzadpx-s1*trans1(2,2))
      trans1(1,4)=rhob*(-cn*dpzadpy-s1*trans1(2,4)-pyi/pzf)
      trans1(1,6)=rhob*(-cn*dpzadp-s1*trans1(2,6)+p/pzf)
      trans1(3,1)=rhob*pyi*dodx
      trans1(3,2)=rhob*pyi*dodpx
      trans1(3,4)=rhob*(omega+pyi*dodpy)
      trans1(3,6)=rhob*pyi*dodp
      trans1(5,1)=-trans1(1,1)*trans1(2,6)+trans1(2,1)*trans1(1,6)
      trans1(5,2)=-trans1(1,2)*trans1(2,6)+trans1(2,2)*trans1(1,6)
      trans1(5,4)=-trans1(1,4)*trans1(2,6)+trans1(2,4)*trans1(1,6)+trans1(3,6)
c      trans1(5,6)=-rhob*(omega+dodp*p)+h0/h1emit**3*aln
      trans1(5,6)=-rhob*(da+xsn+sn*dpzf/pzf+ddodp*p)+h0/h1emit**3*aln
c      if(abs(rhob) > 1.e10)then
c        write(*,*)'tbb0 ',-rhob*(omega+dodp*p),-rhob*(da+xsn+sn*dpzf/pzf+ddodp*p)
c      endif
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.)
      cod(1)=xf
      cod(3)=yf
      cod(2)=min(p,max(-p,pxf))
      cod(4)=min(p,max(-p,pyi))
      cod(5)=zf-dvemit*aln
      if(enarad)then
        bsir0=bsir0+bsi1*yi/rhob-bsi2*cod(3)/rhob
        call tradke(trans,cod,beam,srot,alr,phi0n,0.d0)
      endif
      return
      end subroutine

      subroutine tbendef1(trans,cod,beam,srot,al0,phi0,fb,rb0,enarad)
      use mathfun, only:xsincos
      implicit none
      real*8 , intent(inout)::trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 , intent(in)::al0,phi0,rb0,fb
      logical*4 , intent(in)::enarad
      real*8 phib1,sn,xsn,cn,dcn,bsi1,bsi2,rb,fl
      integer*4 i
      logical*4 en
      rb=rb0
      bsi2=0.d0
      bsi1=merge(1.d0,0.d0,rb < .5d0)
      en=enarad
      do i=1,2
        fl=0.5d0*fb*rb
        phib1=phi0*fl/al0
        call xsincos(phib1,sn,xsn,cn,dcn)
        call tbendebody0(trans,cod,beam,srot,fl,
     $     phib1,sn,xsn,cn,dcn,.5d0*fb,bsi1,bsi2,
     $     en)
        rb=1.d0/6.d0/rb
        bsi1=0.d0
        if(rb < .5d0)then
          bsi2=1.d0
          en=.false.
        endif
      enddo
      return
      end subroutine

      end module

      subroutine tbende(trans,cod,beam,srot,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,dz,theta,dtheta,dchi2,alg,phig,
     $     fb1,fb2,mfring,fringe,eps0,enarad,alcorr,next,l)
      use bendeb
      use tfstk
      use ffs_flag
      use tmacro
      use multa, only:nmult
      use kradlib, only:tradke,dphipol
      use temw,only:tsetr0
      use chg,only:tchge
      use tparastat,only:setndivelm
      use sad_basics
      use mathfun, only:xsincos
      implicit none
      integer*4 ,intent(in)::  mfring,l
      integer*4 ndiv,nrad,n,n1,n2
      real*8 ,intent(in):: al0,phib,phi0,apsi1,apsi2,ak,dx,dy,dz,
     $     theta,dtheta,fb1,fb2,eps0,dchi2,alg,phig
      real*8 phibl,bsi1,bsi2,alx,alr,
     $     dxfr1,dyfr1,dxfr2,dyfr2,phi1,
     $     eps,akn,tanp1,tanp2,f,
     $     dyfra1,dyfra2,psi1,psi2,cod11,
     $     cn,sn,xsn,dcn,phin,aln,alx0,akx0,epsr1,
     $     rbc,akc,alc,phic,f1r,f2r
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      complex*16 akm(0:nmult)
      logical*4 ,intent(in):: enarad,alcorr,fringe,next
      logical*4 prev
      if(alcorr .and.
     $     mfring /= 0 .and. al0 /= 0.d0
     $     .and. phi0 /= 0.d0)then
        al=al0-((phi0*fb1)**2+(phi0*fb2)**2)/48.d0/al0
     1       *sin(.5d0*(phi0*(1.d0-psi1-psi2)-apsi1-apsi2))
     $       /sin(.5d0*phi0)
      else
        al=al0
      endif
      if(phi0 == 0.d0)then
        if(ak == 0.d0)then
          call tsteee(trans,cod,beam,srot,al0,-phib,dx,dy,theta,enarad,
     $         apsi1,apsi2,
     $         fb1,fb2,mfring,fringe,eps0,next,l)
        elseif(phib == phi0)then
          call tquade(trans,cod,beam,srot,al0,ak,0.d0,
     1     dx,dy,theta,enarad,fringe,0.d0,0.d0,0.d0,0.d0,0,eps0,
     $     .true.,.false.,l)
        else
          akm=(0.d0,0.d0)
          akm(0)=phib-phi0
          akm(1)=ak
          call tmulte(trans,cod,beam,srot,l,al,akm,0.d0,
     $         0.d0,psi1,psi2,apsi1,apsi2,
     1         dx,dy,0.d0,0.d0,0.d0,theta,dtheta,
     $         eps0,enarad,fringe,
     $         0.d0,0.d0,0.d0,0.d0,mfring,fb1,fb2,.true.,
     $         0.d0,0.d0,0.d0,0.d0,0.d0,
     $         .false.,.false.)
        endif
        return
      elseif(phib == 0.d0)then
        call tchge(trans,cod,beam,srot,
     $       dx,dy,dz,theta,dtheta,dchi2,alg,phig,.true.)
        call tbdrifte(trans,cod,beam,srot,al,phi0,h0,h1emit,dvemit,
     $       irad)
        call tchge(trans,cod,beam,srot,
     $       dx,dy,dz,theta,dtheta,dchi2,alg-al0,phig-phi0,.false.)
        return
      elseif(al == 0.d0)then
        call tbthie(trans,cod,beam,srot,phib,phi0,dx,dy,dz,theta,dtheta,dchi2)
        return
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,dz,theta,dtheta,dchi2,alg,phig,.true.)
c      write(*,'(a,1p10g12.4)')'tbende-2 ',cod(1:5)
      if(enarad)then
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      phibl=phib/al
      rhob=1.d0/phibl
      rho0=al/phi0
      prev=bradprev /= 0.d0
      f1r=0.d0
      f2r=0.d0
      if(fb1 /= 0.d0)then
        if(mfring > 0 .or. mfring == -1)then
          dxfr1=fb1**2*phibl/24.d0
          dyfr1=fb1*phibl**2/6.d0
          dyfra1=merge(4.d0*dyfr1/fb1**2,0.d0,fringe)
          call tblfre(trans,cod,beam,dxfr1,dyfr1,dyfra1)
          f1r=0.5d0*fb1
        endif
      endif
      if(fb2 /= 0.d0 .and. mfring > 0 .or. mfring == -2)then
        f2r=0.5d0*fb2
      endif
      rbc=1.d0-(f1r+f2r)/al0
      phic=phi0*rbc
      eps=merge(epsbend,epsbend*eps0,eps0 == 0.d0)
      drhob=rhob-rho0
      aind=rho0/phi0*ak
      b=brhoz/rhob
      b1=b*aind/rhob
      ndiv=merge(1,1+int(abs(phic/eps)),ak == 0.d0)
      if(enarad)then
        epsr1=merge(epsrad,epsrad*eps0,eps0 == 0.d0)
        nrad=1+int(abs(al0*rbc/epsr1*crad*(h0*b)**2))
        ndiv=max(ndiv,int(nrad*emidiv*emidib),
     $       1+int(abs(phib*h0*anrad)/epsr1/1.d6*emidiv*emidib))
c        write(*,*)'tbende ',ndiv,nrad,phib,b,epsrad,crad
      endif
      if(calpol)then
        ndiv=max(ndiv,1+int(max(abs(phic),1.d-6)*h0*gspin/dphipol))
      endif
      call setndivelm(l,ndiv)
      call tinitr(trans1)
      tanp1=tan(psi1*phi0+apsi1)
      tanp2=tan(psi2*phi0+apsi2)
      f=1.d0/rho0
      call tbedge(trans,cod,beam,al,phib,psi1*phi0+apsi1,.true.)
      cod11=cod(1)
      akc=ak*rbc
      alc=al*rbc
      akn=akc/ndiv
      aln=alc/ndiv
      phin=phic/ndiv
      if(ak == 0.d0)then
        call xsincos(phin,sn,xsn,cn,dcn)
        bsi1=1.d0
        bsi2=0.d0
        n1=1
        n2=ndiv
        if(f1r /= 0.d0)then
          n1=0
        endif
        if(f2r /= 0.d0)then
          n2=ndiv+1
        endif
        do n=n1,n2
          if(n == 0)then
            call tbendef1(trans,cod,beam,srot,al0,phi0,fb1,rbl,enarad)
          elseif(n == ndiv+1)then
            call tbendef1(trans,cod,beam,srot,al0,phi0,fb2,rbh,enarad)
            alr=f2r
            phi1=alr*rbl/al0*phi0
          else
            if(n == n2)then
              bsi2=1.d0
            endif
            call tbendebody0(trans,cod,beam,srot,aln,
     $           phin,sn,xsn,cn,dcn,aln,bsi1,bsi2,
     $           enarad .and. n /= n2)
            alr=aln
            phi1=phin
          endif
          bsi1=0.d0
        enddo
      else
        tbinit=.true.
        n1=1
        n2=ndiv
        if(f1r /= 0.d0)then
          n1=-1
        endif
        if(f2r /= 0.d0)then
          n2=ndiv+2
        endif
        alx0=0.d0
        akx0=0.d0
        bsi1=1.d0
        bsi2=0.d0
        do n=n1,n2
          call tbendal(n,ndiv,f1r,f2r,aln,alx,alr)
          akx=ak*alx/al0
          phi1=phi0*alx/al0
          call tbendecorr(trans,cod,beam,
     $         (akx+akx0)*.5d0,(alx+alx0)*.5d0)
          akx0=akx
          alx0=alx
          if(n == n2)then
            bsi2=1.d0
          endif
          call tbendebody(trans,cod,beam,srot,alx,phi1,
     $         akx,alr,bsi1,bsi2,
     $         enarad .and. n /= n2)
          if(n <= 0 .or. n >= ndiv)then
            tbinit=.true.
          endif
          bsi1=0.d0
        enddo
        call tbendecorr(trans,cod,beam,akx0*.5d0,alx0*.5d0)
      endif
      if(.not. next)then
        bradprev=0.d0
      endif
      call tbedge(trans,cod,beam,al,phib,psi2*phi0+apsi2,.false.)
      if(f2r /= 0.d0)then
        dxfr2=-fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        dyfra2=merge(4.d0*dyfr2/fb2**2,0.d0,fringe)
        call tblfre(trans,cod,beam,dxfr2,dyfr2,dyfra2)
      endif
      if(enarad)then
        call tradke(trans,cod,beam,srot,alr,phi1,0.d0)
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,dz,theta,dtheta,dchi2,alg-al0,phig-phi0,.false.)
      return
      end

      subroutine tbdrifte(trans,cod,beam,srot,al,phi0,
     $     h0,h1emit,dvemit,irad)
      use drife
      use temw,only:tmulbs
      use sad_basics
      use mathfun
      implicit none
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: phi0,al,h0,h1emit,dvemit
      real*8 cp,sp,pr,pxi,pzf,
     $     trans1(6,6),xi,pyi,pzi,pxf,xf,dpzipxi,dpzipyi,dpzip,
     $     dpzfpxi,dpzfpyi,dpzfp,rho0,dl,dcp,xs,a
      integer*4 ,intent(in):: irad
      real*8 ,parameter ::psqmax=0.9999d0
      call xsincos(phi0,sp,xs,cp,dcp)
      rho0=al/phi0
      call tdrife0(trans,cod,beam,srot,rho0*sp,0.d0,0.d0,.true.,.false.,irad)
      xi=cod(1)-rho0*dcp
      pr=1.d0+cod(6)
      pxi=cod(2)
      pyi=cod(4)
      a=pxi**2+pyi**2
      pzi=pr*sqrtl(1.d0-a/pr**2)
c      pzi=sqrt(max(1.d-4,(pr-pxi)*(pr+pxi)-pyi**2))
      pzf=pzi*cp-pxi*sp
      pxf=pzi*sp+pxi*cp
      xf=xi*pzi/pzf
      dpzipxi=-pxi/pzi
      dpzfpxi= cp*dpzipxi-sp
      dpzipyi=-pyi/pzi
      dpzfpyi= cp*dpzipyi
      dpzip  = pr/pzi
      dpzfp  = cp*dpzip
      trans1(1,1)=pzi/pzf
      trans1(1,2)=(dpzipxi-dpzfpxi*pzi/pzf)/pzf*xi
      trans1(1,3)=0.d0
      trans1(1,4)=(dpzipyi-dpzfpyi*pzi/pzf)/pzf*xi
      trans1(1,6)=(dpzip  -dpzfp  *pzi/pzf)/pzf*xi
      trans1(2,1)=0.d0
      trans1(2,2)=cp+sp*dpzipxi
      trans1(2,3)=0.d0
      trans1(2,4)=   sp*dpzipyi
      trans1(2,6)=   sp*dpzip
      trans1(3,1)=sp*pyi/pzf
      trans1(3,2)=-xi*sp*pyi/pzf**2*dpzfpxi
      trans1(3,3)=1.d0
      trans1(3,4)= xi*sp*(1.d0-pyi*dpzfpyi/pzf)/pzf
      trans1(3,6)=-xi*sp*pyi/pzf**2*dpzfp
      trans1(5,1)=-pr/pzf*sp
      trans1(5,2)= pr/pzf**2*xi*sp*dpzfpxi
      trans1(5,3)=0.d0
      trans1(5,4)= pr/pzf**2*xi*sp*dpzfpyi
      dl=rho0*xsin(phi0)
      trans1(5,6)=-xi*sp*(1.d0-pr*dpzfp  /pzf)/pzf
     $     +h0/h1emit**3*dl
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.)
      cod(1)=xf
      cod(2)=pxf
      cod(3)=cod(3)+xi*sp*pyi/pzf
      cod(5)=cod(5)-pr/pzf*xi*sp-dl*dvemit
      return
      end

      subroutine qbend(trans,cod,
     1     al0,phib,phi0,psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,dz,theta,dtheta,dchi2,alg,phig,fb1,fb2,mfring,fringe,eps0,coup,l)
      use sad_basics
      implicit none
      integer*4 ,intent(in):: mfring,l
      real*8 ,intent(inout):: trans(4,5),cod(6)
      real*8 transe(6,12),beam(42),srot(3,9)
      real*8 ,intent(in):: dx,dy,dz,theta,fb1,fb2,al0,phib,phi0,
     $     psi1,psi2,ak,dtheta,eps0,apsi1,apsi2,dchi2,alg,phig
      logical*4 ,intent(out):: coup
      logical*4 ,intent(in):: fringe
      call tinitr(transe)
      call tbende(transe,cod,beam,srot,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,dz,theta,dtheta,dchi2,alg,phig,
     $     fb1,fb2,mfring,fringe,eps0,.false.,.true.,.false.,l)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) /= 0.d0 .or. trans(1,4) /= 0.d0 .or.
     $     trans(2,3) /= 0.d0 .or. trans(2,4) /= 0.d0
      return
      end

      subroutine tblfre(trans,cod,beam,dxfr,dyfr,dyfra)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      use sad_basics
      implicit none
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42)
      real*8 ,intent(in):: dxfr,dyfr,dyfra
      real*8 trans1(6,6),pr,ysq,dyfraysq
      pr=1.d0+cod(6)
      ysq=cod(3)**2
      dyfraysq=ysq*dyfra
      call tinitr(trans1)
      trans1(1,6)=dxfr/pr**2
      trans1(4,3)=(dyfr-3.d0*dyfraysq)/pr
      trans1(4,6)=-dyfr*cod(3)/pr**2
      trans1(5,2)=trans1(1,6)
      trans1(5,3)=-trans1(4,6)
      trans1(5,6)=-(2.d0*dxfr*cod(2)+
     $     (dyfr-.5d0*dyfraysq)*ysq)/pr**3
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam,trans1,.true.)
      cod(1)=cod(1)+dxfr*cod(6)/pr
      cod(4)=cod(4)+(dyfr-dyfraysq)*cod(3)/pr
      cod(5)=cod(5)+(dxfr*cod(2)+
     $     (.5d0*dyfr-.25d0*dyfraysq)*ysq)/pr**2
      return
      end

      subroutine tbthie(trans,cod,beam,srot,phib,phi0,
     1                 dx,dy,dz,theta,dtheta,dchi2)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      use chg,only:tchge
      use sad_basics
      implicit none
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: phib,phi0,dx,dy,dz,theta,dtheta,dchi2
      real*8 trans1(6,6)
      call tchge(trans,cod,beam,srot,
     $     dx,dy,dz,theta,dtheta,dchi2,0.d0,.5d0*phi0,.true.)
      call tinitr(trans1)
      trans1(2,6)=phi0
      trans1(5,1)=-phi0
      trans(:,1:irad)=matmul(trans1,trans(:,1:irad))
      call tmulbs(beam ,trans1,.true.)
      cod(2)=cod(2)+(phi0-phib)+cod(6)*phi0
      cod(5)=cod(5)-phi0*cod(1)
      call tchge(trans,cod,beam,srot,
     $     dx,dy,dz,theta,dtheta,dchi2,0.d0,-.5d0*phi0,.false.)
      return
      end
