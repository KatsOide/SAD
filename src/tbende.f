      module bendeb
      use tfstk
      use bendib
      use temw
      use tspin
      implicit none
      logical*4 tbinit
      real*8 al,b1,b,aind,trans1(6,13),cx,cy,
     $     drhopp,akxsqp,dcxp,sxkxp,dcyp,sykyp,phixp,phiyp,
     $     akysqp,dcxkxp,aksxp,aksyp,xsxkxp

      contains
      subroutine tbendeinit(ak1,al)
      implicit none
      real*8 , intent(in) :: ak1,al
      rhoe=rhob*p
      call tbendiinit(ak1,al)
      cx=1.d0+dcx
      cy=1.d0+dcy
      drhopp=1.d0/rhoe/p
      akxsqp=-akxsq/p
      phixp=-0.5d0*phix/p
      if(akxsq .gt. 0.d0)then
        dcxp=spx*phixp
        dcxkxp=dcxp/akxsq+dcxkx/p
        sxkxp= 0.5d0*sxkx/p+cx*phixp/akx
        xsxkxp=(-dcx*phixp/akxsq-1.5d0*xsxkx/p)/akx
        aksxp=-0.5d0*aksx/p+akx*cx*phixp
      elseif(akxsq .lt. 0.d0)then
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
      if(akysq .gt. 0.d0)then
        dcyp=phiyp*spy
        sykyp= 0.5d0*syky/p+cy*phiyp/aky
        aksyp=-0.5d0*aksy/p+aky*cy*phiyp
      elseif(akysq .lt. 0.d0)then
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
      use mathfun
      implicit none
      real*8 , intent(in) :: ak1,al
      real*8 trans(6,12),cod(6),beam(42)
      real*8 pzi,dpzi,xr,s,a,trans2(6,6),dk,pxf,xf,yf,zf,phi
      phi=al/rho0
      p  =1.d0+cod(6)
      xi =cod(1)
      pxi=min(p,max(-p,cod(2)))/p
      yi =cod(3)
      pyi=min(p,max(-p,cod(4)))/p
      s=min(psqmax,pxi**2+pyi**2)
      dpzi=sqrt1(-s)
      pzi=1.d0+dpzi
      a=1.d0/(pzi-phi*pxi)
      xf =(-al*dpzi*pxi+pzi*xi)*a
      pxf=pxi+phi*dpzi
      yf =yi+(-al*dpzi+xf*phi)*pyi/pzi
      zf =cod(5)-phi*dpzi*xf/pzi+al*dpzi*(.5d0*pzi+.5d0-1.d0/pzi)
      xr=xf/rho0
      pxf=pxf-ak1*xf*xr*(0.5d0-xr*(2.d0-xr)/12.d0)/p
      dk=-ak1*xr*(1.d0-xr*(1.d0-2.d0*xr))
      trans2(1,1)=pzi*a
      trans2(1,2)=((al+phi*xi)*(pzi**2+pxi**2)
     $     -al*(phi*pxi**3+pzi**3))*a**2/pzi/p
      trans2(1,3)=0
      trans2(1,4)=(al*(1.d0-phi*pxi)+phi*xi)*pxi*pyi*a**2/pzi/p
      trans2(1,5)=0.d0
      trans2(1,6)=pxi/pzi*a*(al*(-a*(1.d0-phi*pxi)+pzi**2)-a*phi*xi)
      trans2(2,1)=dk*trans2(1,1)
      trans2(2,2)=1.d0-phi*pxi/pzi+dk*trans2(1,2)
      trans2(2,3)=0.d0
      trans2(2,4)=-phi*pyi/pzi+dk*trans2(1,4)
      trans2(2,5)=0.d0
      trans2(2,6)=-phi/pzi*dpzi+dk*trans2(1,6)
      trans2(3,1)=phi*pyi*a
      trans2(3,2)=(-al*phi*(pzi**2+pxi**2)
     $     +(al+phi*xi)*(pxi+phi*pzi))*a**2*pyi/pzi/p
      trans2(3,3)=1.d0
      trans2(3,4)=(-al*(pzi**2*(pzi-phi*pxi)+phi*pxi*pyi**2)
     $     +(pyi**2+pzi*(pzi-phi*pxi))*(al+phi*xi))*a**2/p/pzi
      trans2(3,5)=0.d0
      trans2(3,6)=-(al*(dpzi*(1.d0+pzi)*phi*pxi-pzi**3)+(al+phi*xi))
     $     *a**2*pyi/pzi/p
      trans2(4,1)=0.d0
      trans2(4,2)=0.d0
      trans2(4,3)=0.d0
      trans2(4,4)=1.d0
      trans2(4,5)=0.d0
      trans2(4,6)=0.d0
      trans2(5,1)=-trans2(1,1)*trans2(2,6)
      trans2(5,2)=-trans2(1,2)*trans2(2,6)+trans2(2,2)*trans2(1,6)
      trans2(5,3)=0.d0
      trans2(5,4)=trans2(3,6)
      trans2(5,5)=1.d0
      trans2(5,6)=(-al*pzi*(a*phi*pxi+pzi)+
     $     a*(al*(-phi*pxi*dpzi*(phi*pxi-1.d0)*a/pzi**2+1.d0)+
     $     phi*xi*(pxi*phi*dpzi*a/pzi**2+1.d0)))/p
      trans2(6,1)=0.d0
      trans2(6,2)=0.d0
      trans2(6,3)=0.d0
      trans2(6,4)=0.d0
      trans2(6,5)=0.d0
      trans2(6,6)=1.d0
      call tmultr5(trans,trans2,irad)
      call tmulbs(beam ,trans2,.true.,.true.)
      cod(1)=xf
      cod(2)=min(p,max(-p,pxf*p))
      cod(3)=yf
      cod(5)=zf
      return
      end subroutine

      subroutine tbendebody(trans,cod,beam,srot,al,phin,
     $     ak1,alr,bsi1,bsi2,enarad)
      use tmacro
      implicit none
      real*8, intent(in):: ak1,al,bsi1,bsi2,alr
      real*8 trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 pxf,xf,phin,dxpx,
     $     dxf,dpxf,dyf,dpyf,hi,hip,pyf,yf,zf,
     $     xiksq
      logical*4 enarad
c      if(enarad)then
c        dpini=cod(6)
c        xr=cod(1)/rho0
c        xe=xi+xi*xr*(.5d0-xr*(2.d0-xr)/12.d0)
c        dxe=1.d0+xr*(1.d0-xr*(.5d0-xr/3.d0))
c        bx=b1*cod(3)
c        by=b+b1*xe
c        bxy=b1*dxe
c        call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $       0.d0,bxy,0.d0,fal,alnr,als,ala,fb1,fb2,prev,next)
c        tbinit=tbinit .or. cod(6) .ne. dpini
c      endif
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
c      dxf =(drhopak+xi)*dcx+aksx/akxsq*pxi
c      dpxf= aksx*(drhopak+xi)+pxi*dcx
c      dyf = dcy*yi+aksy/akysq*pyi
c      dpyf=aksy*yi+dcy*pyi
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
      trans1(1,6)=drhopp*dcxkx+drhop*dcxkxp
     $     +xi*dcxp+(sxkxp-0.5d0*sxkx/p)*pxi
      trans1(2,1)=aksx*p
      trans1(2,2)=cx
c      trans1(2,6)=aksxp*(drhopak+xi)*p
c     $     +aksx*(drhopakp*p+drhopak+xi)+pxi*dcxp
      trans1(2,6)=p*(drhopp*sxkx+drhop*sxkxp)+drhop*sxkx
     $     +(p*aksxp+aksx)*xi+dcxp*pxi
      trans1(3,3)=cy
      trans1(3,4)=syky/p
      trans1(3,6)=yi*dcyp+(sykyp-0.5d0*syky/p)*pyi
      trans1(4,3)=aksy*p
      trans1(4,4)=cy
      trans1(4,6)=(aksyp*p+aksy)*yi+pyi*dcyp
      trans1(5,1)=-trans1(1,1)*trans1(2,6)+trans1(2,1)*trans1(1,6)
      trans1(5,2)=-trans1(1,2)*trans1(2,6)+trans1(2,2)*trans1(1,6)
      trans1(5,3)=-trans1(3,3)*trans1(4,6)+trans1(4,3)*trans1(3,6)
      trans1(5,4)=-trans1(3,4)*trans1(4,6)+trans1(4,4)*trans1(3,6)
      trans1(5,6)=-0.25d0*(-dxpx/p
     $     +trans1(1,6)*pxf+trans1(2,6)/p*xf
     $     +trans1(3,6)*pyf+trans1(4,6)/p*yf
     $     -1.d0/rhoe/p*(drhop*xsxkx-sxkx*xi-dcxkx*pxi)
     $     -(5.d0/rho0-1.d0/rhoe)
     $     *(drhopp*xsxkx+drhop*xsxkxp-sxkxp*xi-dcxkxp*pxi))
     $     -hip*al*.5d0+h0/h1emit**3*al
c      write(*,'(a,1p8g15.7)')'tmulte ',drhopp,drhop,xsxkx,xsxkxp,
c     $     sxkxp,dcxkxp
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      cod(1)=xf
      cod(2)=min(p,max(-p,pxf*p))
      cod(3)=yf
      cod(4)=min(p,max(-p,pyf*p))
      cod(5)=zf-dvemit*al
      if(enarad)then
        bsi=bsi+bsi1*(yi/rhob+ak1*xi*yi)
     $       -bsi2*(yf/rhob+ak1*xf*yf)
        call tradke(trans,cod,beam,srot,alr,phin,0.d0)
        tbinit=.true.
      endif
      return
      end subroutine

      subroutine tbendebody0(trans,cod,beam,srot,aln,
     $     phi0n,snphi0,sinsq0,csphi0,alr,
     $     bsi1,bsi2,enarad)
      use tmacro
      use temw
      use mathfun
      implicit none
      real*8, intent(in):: phi0n,snphi0,sinsq0,csphi0,aln,bsi1,bsi2,alr
      real*8 trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 xi,pxi,yi,pyi,dp,pr,rhoe,
     $     dpzinv,phsq,dtn,s,dpz1,pz1,drho,dpx,pxf,
     $     pz2,d,sinda,da,dpz2,xf
      logical*4 , intent(in) ::enarad
      xi=cod(1)
      dp=cod(6)
      pr=1.d0+dp
      pxi=min(pr,max(-pr,cod(2)))/pr
      yi=cod(3)
      pyi=min(pr,max(-pr,cod(4)))/pr
      rhoe=rhob*pr
      s=min(psqmax,pxi**2+pyi**2)
      dpz1=sqrt1(-s)
      pz1=1.d0+dpz1
      drho=drhob+rhoe*dpz1+rhob*dp
      dpx=-(xi-drho)/rhoe*snphi0-sinsq0*pxi
      pxf=pxi+dpx
      s=min(psqmax,pxf**2+pyi**2)
      dpz2=sqrt1(-s)
      pz2=1.d0+dpz2
      d=pxf*pz1+pxi*pz2
      if(d .eq. 0.d0)then
        sinda=min(1.d0,max(-1.d0,2.d0*pxf*pz2/(pxf**2+pz2**2)))
      else
        sinda=min(1.d0,max(-1.d0,dpx*(pxf+pxi)/d))
      endif
      s=sinda**2
      if(s .gt. 2.d-4)then
        da=sinda
     1       *(1.d0+s*(a3+s*(a5+s*(a7+s*(a9+s*(a11+s*(a13+s*a15)))))))
      else
        da=sinda*(1.d0+s*(a3+s*(a5+a7*s)))
      endif
      trans1(2,1)=-snphi0/rhob
      trans1(2,2)=csphi0-snphi0*pxi/pz1
      trans1(2,6)=snphi0/pz1
      trans1(2,4)=-pyi*trans1(2,6)
      dpzinv=dpx/pz1/pz2*(pxi+pxf)/(pz1+pz2)
      phsq=(1.d0-pyi)*(1.d0+pyi)
      if(d .eq. 0.d0)then
        dtn=-2.d0*pxi/pz1
      else
        dtn= phsq*sinda/pz1/pz2
      endif
      trans1(1,1)=csphi0+snphi0*pxf/pz2
      trans1(1,2)=rhob*(snphi0-dtn*csphi0+pxi*pxf/pz1/pz2*snphi0)
      trans1(1,6)=rhob*(dpzinv+sinsq0/pz1-pxf/pz2*trans1(2,6))
      trans1(1,4)=-pyi*trans1(1,6)
      trans1(5,1)=rhob*trans1(2,1)/pz2
      trans1(5,2)=rhob*(dpzinv-sinsq0/pz2-snphi0*pxi/pz1/pz2)
      trans1(5,6)=rhob*(trans1(2,6)/pz2-dtn/phsq)
      trans1(5,4)=-pyi*trans1(5,6)
      trans1(3,1)=-pyi*trans1(5,1)
      trans1(3,2)=-pyi*trans1(5,2)
      trans1(3,4)=rhob*(phi0n-da)-pyi*trans1(5,4)
      trans1(3,6)=trans1(5,4)
      trans1(5,6)=trans1(5,6)-(phi0n-da)*rhob+h0/h1emit**3*aln
      xf=xi*csphi0+rhoe*(snphi0*pxi-dpx*(pxi+pxf)/(pz1+pz2))
     1     +drho*sinsq0
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      cod(1)=xf
      cod(3)=yi+pyi*rhoe*(phi0n-da)
      cod(2)=min(pr,max(-pr,pxf*pr))
      cod(4)=min(pr,max(-pr,pyi*pr))
      cod(5)=cod(5)-phi0n*(dp*rhob+drhob)+da*rhoe-dvemit*aln
      if(enarad)then
        bsi=bsi+bsi1*yi/rhob-bsi2*cod(3)/rhob
        call tradke(trans,cod,beam,srot,alr,phi0n,0.d0)
      endif
      return
      end subroutine

      subroutine tbrote(trans1,cod,srot,phi0,dtheta)
      implicit none
      real*8 trans1(6,6),cod(6),phi0,dtheta,chi1,chi2,chi3,sphi0,
     $     coschi2,sdt,cphi0,rr(3,3),srot(3,9)
      integer*4 i
      cphi0=cos(phi0*.5d0)
      sphi0=sin(phi0*.5d0)
      sdt=sin(dtheta)
      chi2=asin(sdt*sphi0)
      coschi2=cos(chi2)
      chi1=asin(sin(dtheta*.5d0)**2*2.d0*sphi0*cphi0/coschi2)
      chi3=asin(sdt*cphi0/coschi2)
c      write(*,*)'tbrote ',chi1,chi2,chi3
      call tsrote(trans1,cod,rr,chi1,chi2,chi3)
      do i=1,9
        srot(1,i)=dot_product(rr(1,:),srot(:,i))
        srot(2,i)=dot_product(rr(2,:),srot(:,i))
        srot(3,i)=dot_product(rr(3,:),srot(:,i))
      enddo
      return
      end

      end module

      subroutine tbende(trans,cod,beam,srot,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,theta,dtheta,
     $     fb1,fb2,mfring,fringe,eps0,enarad,alcorr,next,l)
      use bendeb
      use tfstk
      use ffs_flag
      use tmacro
      use multa, only:nmult
      use temw
      implicit none
      integer*4 mfring,ndiv,nrad,n,l,n1
      real*8 al0,phib,phi0,psi1,psi2,ak,dx,dy,theta,dtheta,
     $     fb1,fb2,eps0,phibl,bsi2,phib1,phib2,phir2,alx,alr,
     $     dxfr1,dyfr1,dxfr2,dyfr2,
     $     eps,akn,tanp1,tanp2,f,
     $     dyfra1,dyfra2,apsi1,apsi2,cod11,
     $     csphin,snphin,sinsqn,phin,aln,alx1,akx1,
     $     rbc,akc,alc,csphib1,csphib2,phic,rb1,rb2,
     $     sinsq1,sinsq2,snphib1,snphib2,alr2
      real*8 trans(6,12),cod(6),beam(42),srot(3,9)
      complex*16 akm(0:nmult)
      logical*4 enarad,alcorr,fringe,next,prev
      if(alcorr .and.
     $     mfring .ne. 0 .and. al0 .ne. 0.d0
     $     .and. phi0 .ne. 0.d0)then
        al=al0-((phi0*fb1)**2+(phi0*fb2)**2)/48.d0/al0
     1       *sin(.5d0*(phi0*(1.d0-psi1-psi2)-apsi1-apsi2))
     $       /sin(.5d0*phi0)
      else
        al=al0
      endif
      if(phi0 .eq. 0.d0)then
        if(ak .eq. 0.d0)then
          call tsteee(trans,cod,beam,srot,al0,-phib,dx,dy,theta,enarad,
     $         apsi1,apsi2,
     $         fb1,fb2,mfring,fringe,next)
        elseif(phib .eq. phi0)then
          call tquade(trans,cod,beam,srot,al0,ak,
     1     dx,dy,theta,enarad,fringe,0.d0,0.d0,0.d0,0.d0,0,eps0,
     $     .true.,.false.,next)
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
      elseif(phib .eq. 0.d0)then
        call tchge(trans,cod,beam,srot,
     $       dx,dy,theta,dtheta,phi0,.true.)
        call tbdrifte(trans,cod,beam,srot,al,phi0,h0,h1emit,dvemit,
     $       irad)
        call tchge(trans,cod,beam,srot,
     $       -dx,-dy,-theta,-dtheta,-phi0,.false.)
        return
      elseif(al .eq. 0.d0)then
        call tbthie(trans,cod,beam,srot,phib,phi0,dx,dy,theta,dtheta)
        return
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,dtheta,phi0,.true.)
      if(enarad)then
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      phibl=phib/al
      rhob=1.d0/phibl
      rho0=al/phi0
      prev=bradprev .ne. 0.d0
      rb1=0.d0
      rb2=0.d0
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fb1**2*phibl/24.d0
          dyfr1=fb1*phibl**2/6.d0
          if(fringe)then
            dyfra1=4.d0*dyfr1/fb1**2
          else
            dyfra1=0.d0
          endif
          call tblfre(trans,cod,beam,dxfr1,dyfr1,dyfra1)
          rb1=0.5d0*fb1/al0
        endif
      endif
      if(fb2 .ne. 0.d0 .and.
     $       mfring .gt. 0 .or. mfring .eq. -2)then
        rb2=0.5d0*fb2/al0
      endif
      rbc=1.d0-rb1-rb2
      phic=phi0*rbc
      if(eps0 .eq. 0.d0)then
        eps=epsbend
      else
        eps=epsbend*eps0
      endif
      drhob=rhob-rho0
      aind=rho0/phi0*ak
      b=brhoz/rhob
      b1=b*aind/rhob
      if(ak .eq. 0.d0)then
        ndiv=1
      else
        ndiv=1+int(abs(phic/eps))
      endif
      if(enarad)then
        nrad=int(abs(al0*rbc/epsrad*crad*(h0*b)**2))
        ndiv=max(ndiv,int(nrad*emidiv*emidib),
     1       int(abs(phib*h0*anrad)/epsrad/1.d6*emidiv*emidib))
        if(ndiv .eq. 1 .and. rb2 .eq. 0.d0)then
          bsi2=1.d0
        else
          bsi2=0.d0
        endif
      endif
      trans1(1,3)=0.d0
      trans1(1,4)=0.d0
      trans1(1,5)=0.d0
      trans1(2,3)=0.d0
      trans1(2,4)=0.d0
      trans1(2,5)=0.d0
      trans1(3,1)=0.d0
      trans1(3,2)=0.d0
      trans1(3,5)=0.d0
      trans1(4,1)=0.d0
      trans1(4,2)=0.d0
      trans1(4,5)=0.d0
      trans1(5,5)=1.d0
      trans1(6,1)=0.d0
      trans1(6,2)=0.d0
      trans1(6,3)=0.d0
      trans1(6,4)=0.d0
      trans1(6,5)=0.d0
      trans1(6,6)=1.d0
      tanp1=tan(psi1*phi0+apsi1)
      tanp2=tan(psi2*phi0+apsi2)
      f=1.d0/rho0
      call tbedge(trans,cod,beam,al,phib,psi1*phi0+apsi1,.true.)
c      if(enarad .and. fb1 .ne. 0.d0)then
c        call tradke(trans,cod,beam,srot,fb1,0.d0,0.d0)
c      else
c        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
c      endif
      cod11=cod(1)
      akc=ak*rbc
      alc=al*rbc
      akn=akc/ndiv
      aln=alc/ndiv
      phin=phic/ndiv
c      write(*,'(a,i5,1p8g15.7)')'tbende ',ndiv,rb1,rb2,al0,alc,aln,phin
      if(ak .eq. 0.d0)then
        trans1(3,3)=1.d0
        trans1(4,3)=0.d0
        trans1(4,4)=1.d0
        trans1(4,6)=0.d0
        trans1(5,3)=0.d0        
        csphin=cos(phin)
        snphin=sin(phin)
        if(csphin .ge. 0.d0)then
          sinsqn=snphin**2/(1.d0+csphin)
        else
          sinsqn=1.d0-csphin
        endif
        if(rb1 .eq. 0.d0)then
          call tbendebody0(trans,cod,beam,srot,aln,
     $         phin,snphin,sinsqn,csphin,aln,1.d0,bsi2,
     $         enarad .and. ndiv .ne. 1)
          n1=2
        else
          phib1=phi0*rb1
          csphib1=cos(phib1)
          snphib1=sin(phib1)
          if(csphib1 .ge. 0.d0)then
            sinsq1=snphib1**2/(1.d0+csphib1)
          else
            sinsq1=1.d0-csphib1
          endif
          call tbendebody0(trans,cod,beam,srot,al0*rb1,
     $         phib1,snphib1,sinsq1,csphib1,al0*rb1*1.5d0,1.d0,bsi2,
     $         enarad)
          n1=1
        endif
        bsi2=0.d0
        do n=n1,ndiv
          if(n .eq. ndiv .and. rb2 .ne. 0.d0)then
            bsi2=1.d0
          endif
          call tbendebody0(trans,cod,beam,srot,aln,
     $         phin,snphin,sinsqn,csphin,aln,0.d0,bsi2,
     $         enarad .and. (n .ne. ndiv .or. rb2 .ne. 0.d0))
        enddo
        if(rb2 .ne. 0.d0)then
          phib2=phi0*rb2
          csphib2=cos(phib2)
          snphib2=sin(phib2)
          if(csphib2 .ge. 0.d0)then
            sinsq2=snphib2**2/(1.d0+csphib2)
          else
            sinsq2=1.d0-csphib2
          endif
          alr2=al0*rb2*1.5d0
          phir2=phib2
          call tbendebody0(trans,cod,beam,srot,al0*rb2,
     $         phib2,snphib2,sinsq1,csphib2,alr2,0.d0,1.d0,
     $         .false.)
        else
          alr2=aln
          phir2=phin
        endif
      else
        tbinit=.true.
        if(rb1 .eq. 0.d0)then
          alx=aln
          akx=akn
          phix=phin
          alr=aln
          n1=2
        else
          alx=al0*rb1
          akx=ak*rb1
          phix=phi0*rb1
          alr=alx*1.5d0
          n1=1
        endif
        call tbendecorr(trans,cod,beam,akx*.5d0,alx*.5d0)
        call tbendebody(trans,cod,beam,srot,alx,phix,
     $       akx,alr,1.d0,bsi2,
     $       enarad .and. (ndiv .ge. n1 .or. rb2 .ne. 0.d0))
        bsi2=0.d0
        do n=n1,ndiv
          akx=(akn+akx)*.5d0
          alx=(aln+alx)*.5d0
          if(n .eq. ndiv .and. rb2 .eq. 0.d0)then
            bsi2=1.d0
          endif
          call tbendecorr(trans,cod,beam,akx,alx)
          akx=akn
          alx=aln
          call tbendebody(trans,cod,beam,srot,aln,phin,
     $         akn,aln,0.d0,bsi2,
     $         enarad .and. (n .ne. ndiv .or. rb2 .ne. 0.d0))
        enddo
        if(rb2 .ne. 0.d0)then
          alx1=al0*rb2
          akx1=ak*rb2
          phir2=phi0*rb2
          alr2=alx1*1.5d0
          call tbendecorr(trans,cod,beam,(akx+akx1)*.5d0,
     $         (alx+alx1)*.5d0)
          call tbendebody(trans,cod,beam,srot,alx1,phir2,
     $         akx1,alr,0.d0,1.d0,.false.)
          alx=alx1
          akx=akx1
        else
          phir2=phin
          alr2=aln
        endif
        call tbendecorr(trans,cod,beam,akx*5d0,alx*.5d0)
      endif
      if(.not. next)then
        bradprev=0.d0
      endif
      call tbedge(trans,cod,beam,al,phib,psi2*phi0+apsi2,.false.)
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=-fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          if(fringe)then
            dyfra2=4.d0*dyfr2/fb2**2
          else
            dyfra2=0.d0
          endif
          call tblfre(trans,cod,beam,dxfr2,dyfr2,dyfra2)
        endif
      endif
      if(enarad)then
        call tradke(trans,cod,beam,srot,alr2,phir2,0.d0)
      endif
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,-dtheta,-phi0,.false.)
      return
      end

      subroutine tbdrifte(trans,cod,beam,srot,al,phi0,
     $     h0,h1emit,dvemit,irad)
      use mathfun
      implicit none
      real*8 trans(6,12),cod(6),beam(42),srot(3,9),
     $     phi0,al,cp,sp,pr,pxi,pzf,
     $     trans1(6,6),xi,pyi,pzi,pxf,xf,dpzipxi,dpzipyi,dpzip,
     $     dpzfpxi,dpzfpyi,dpzfp,rho0,h0,dl,dcp,
     $     h1emit,dvemit,a,psqmax
      integer*4 irad
      parameter (psqmax=0.9999d0)
      cp=cos(phi0)
      sp=sin(phi0)
      if(cp .ge. 0.d0)then
        dcp=sp**2/(1.d0+cp)
      else
        dcp=1.d0-cp
      endif
      rho0=al/phi0
      call tdrife(trans,cod,beam,srot,rho0*sp,0.d0,0.d0,0.d0,0.d0,
     $     .true.,.false.,irad)
      xi=cod(1)+rho0*dcp
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
      call tmulbs(beam ,trans1,.true.,.true.)
      cod(1)=xf
      cod(2)=pxf
      cod(3)=cod(3)+xi*sp*pyi/pzf
      cod(5)=cod(5)-pr/pzf*xi*sp-dl*dvemit
      return
      end

      subroutine qbend(trans,cod,
     1     al0,phib,phi0,psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,theta,dtheta,fb1,fb2,mfring,fringe,eps0,coup)
      implicit none
      integer*4 mfring
      real*8 trans(4,5),cod(6),transe(6,12),beam(42),srot(3,9),
     $     dx,dy,theta,fb1,fb2,al0,phib,phi0,
     $     psi1,psi2,ak,dtheta,eps0,apsi1,apsi2
      logical*4 coup,fringe
      call tinitr(transe)
      call tbende(transe,cod,beam,srot,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,theta,dtheta,
     $     fb1,fb2,mfring,fringe,eps0,.false.,.true.,.false.,1)
      call qcopymat(trans,transe,.false.)
c      write(*,*)'qbend ',cod,al0,phi0,phib
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end

      subroutine tblfre(trans,cod,beam,dxfr,dyfr,dyfra)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      real*8 trans(6,12),cod(6),beam(42),trans1(6,13),
     $     dxfr,dyfr,pr,dyfra,ysq,dyfraysq
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
      call tmulbs(beam,trans1,.true.,.true.)
      cod(1)=cod(1)+dxfr*cod(6)/pr
      cod(4)=cod(4)+(dyfr-dyfraysq)*cod(3)/pr
      cod(5)=cod(5)+(dxfr*cod(2)+
     $     (.5d0*dyfr-.25d0*dyfraysq)*ysq)/pr**2
      return
      end

      subroutine tbthie(trans,cod,beam,srot,phib,phi0,
     1                 dx,dy,theta,dtheta)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      real*8 trans(6,12),cod(6),beam(42),phib,phi0,dx,dy,theta,
     $     trans1(6,13),dtheta,srot(3,9)
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,dtheta,phi0,.true.)
      call tinitr(trans1)
      trans1(2,6)=phi0
      trans1(5,1)=-phi0
      call tmultr(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      cod(2)=cod(2)+(phi0-phib)+cod(6)*phi0
      cod(5)=cod(5)-phi0*cod(1)
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,-dtheta,-phi0,.false.)
      return
      end
