      subroutine tbende(trans,cod,beam,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,theta,dtheta,
     $     fb1,fb2,mfring,fringe,eps0,enarad,alcorr,ld)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      real*8 a3,a5,a7,a9,ptmax
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,ptmax=0.9999d0)
      integer*4 ld,mfring,ndiv,nrad,k,n
      real*8 al0,phib,phi0,psi1,psi2,ak,dx,dy,theta,dtheta,
     $     fb1,fb2,eps0,
     $     al,dphix,rho0,
     $     dxfr1,dyfr1,dxfr2,dyfr2,
     $     eps,drhob,aind,f1r,
     $     b,b1,aln,akn,phibn,phi0n,csphi0,snphi0,sqphi0,tanp1,tanp2,
     $     f,xi,xr,xe,ak1,bx,by,bxy,dp,pr,pxi,yi,pyi,rhoe,s,
     $     dpz2,pz2,d,sinda,da,dpzinv,phsq,xf,xfr,dxe,t21,t43,
     $     dphiy,rhob,dpz1,pz1,drho,dpx,pxf,dtn,dyfra1,dyfra2,
     $     apsi1,apsi2,als
      real*8 trans(6,12),cod(6),beam(42),trans1(6,13)
      logical*4 enarad,alcorr,fringe
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
        call tsteee(trans,cod,beam,al0,-phib,dx,dy,theta,enarad,
     $       apsi1,apsi2,
     $       fb1,fb2,mfring,fringe,ld)
        return
      elseif(phib .eq. 0.d0)then
        call tchge(trans,cod,beam,-dx,-dy,theta,.true.,ld)
        call tbdrifte(trans,cod,beam,al,phi0,h0,h1emit,dvemit,
     $       irad,calpol,ld)
        call tchge(trans,cod,beam,dx,dy,-theta,.false.,ld)
        return
      elseif(al .le. 0.d0)then
        call tbthie(trans,cod,beam,phib,phi0,dx,dy,theta,ld)
        return
      endif
      call tchge(trans,cod,beam,-dx,-dy,theta,.true.,ld)
      if(dtheta .ne. 0.d0)then
        dphix=      phi0*sin(.5d0*dtheta)**2
        dphiy= .5d0*phi0*sin(dtheta)
        cod(2)=cod(2)+dphix
        cod(4)=cod(4)+dphiy
      else
        dphix=0.d0
        dphiy=0.d0
      endif
      rhob=al/phib
      rho0=al/phi0
      f1r=0.d0
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fb1**2/rhob/24.d0
          dyfr1=fb1/rhob**2/6.d0
          if(fringe)then
            dyfra1=4.d0*dyfr1/fb1**2
          else
            dyfra1=0.d0
          endif
          call tblfre(trans,cod,beam,dxfr1,dyfr1,dyfra1,ld)
          f1r=fb1
        endif
      endif
      if(eps0 .eq. 0.d0)then
        eps=1.d-6
      else
        eps=1.d-6*eps0
      endif
      drhob=rhob-rho0
      aind=rho0/phi0*ak
      b=brhoz/rhob
      b1=b*aind/rhob
      ndiv=1+int(sqrt(abs(ak*al)/eps/12.d0))
      if(enarad)then
        nrad=int(al/epsrad*crad*(h0*b)**2)
        ndiv=max(ndiv,int(nrad*emidiv*emidib),
     1       int(min(pi2,abs(phib))/epsrad/1.d3*emidiv*emidib))
      endif
      aln=al/ndiv
      akn=ak/ndiv
      phibn=phib/ndiv
      phi0n=phi0/ndiv
      csphi0=cos(phi0n)
      snphi0=sin(phi0n)
      if(csphi0 .ge. 0.d0)then
        sqphi0=snphi0**2/(1.d0+csphi0)
      else
        sqphi0=1.d0-csphi0
      endif
      tanp1=tan(psi1*phi0+apsi1)
      tanp2=tan(psi2*phi0+apsi2)
      f=1.d0/rho0
      call tbedge(trans,cod,beam,al,phib,psi1*phi0+apsi1,.true.,ld)
      trans1(1,5)=0.d0
      trans1(2,5)=0.d0
      trans1(3,5)=0.d0
      trans1(4,1)=0.d0
      trans1(4,2)=0.d0
      trans1(4,4)=1.d0
      trans1(4,5)=0.d0
      trans1(4,6)=0.d0
      trans1(5,5)=1.d0
      trans1(6,1)=0.d0
      trans1(6,2)=0.d0
      trans1(6,3)=0.d0
      trans1(6,4)=0.d0
      trans1(6,5)=0.d0
      trans1(6,6)=1.d0
      als=0.d0
      do 100 n=1,ndiv
        xi=cod(1)
        xr=xi/rho0
        xe=xi+xi*xr*(.5d0-xr*(2.d0-xr)/12.d0)
        dxe=1.d0+xr*(1.d0-xr*(.5d0-xr/3.d0))
        if(n .eq. 1)then
          ak1=akn*.5d0
          if(enarad)then
            bx=b1*cod(3)
            by=b+b1*xe
            bxy=b1*dxe
            call trade(trans,beam,cod,bx,by,0.d0,0.d0,
     $           0.d0,bxy,0.d0,
     1           f*aln*.5d0-tanp1,fb1,fb2,0.d0,al,aln*.5d0)
          endif
        else
          ak1=akn
          if(enarad)then
            bx=b1*cod(3)
            by=b+b1*xe
            bxy=b1*dxe
            call trade(trans,beam,cod,bx,by,0.d0,0.d0,
     $           0.d0,bxy,0.d0,f*aln,
     $           fb1,fb2,als,al,aln)
          endif
        endif
        als=als+aln
        dp=cod(6)
        pr=1.d0+dp
        pxi=min(ptmax,max(-ptmax,(cod(2)-ak1*xe)/pr))
        yi=cod(3)
        pyi=min(ptmax,max(-ptmax,(cod(4)+ak1*yi)/pr))
        rhoe=rhob*pr
        s=min(ptmax,pxi**2+pyi**2)
        dpz1=-s/(1.d0+sqrt(1.d0-s))
        pz1=1.d0+dpz1
        drho=drhob+rhoe*dpz1+rhob*dp
        dpx=-(xi-drho)/rhoe*snphi0-sqphi0*pxi
        pxf=pxi+dpx
        s=min(ptmax,pxf**2+pyi**2)
        dpz2=-s/(1.d0+sqrt(1.d0-s))
        pz2=1.d0+dpz2
        d=pxf*pz1+pxi*pz2
        if(d .eq. 0.d0)then
          sinda=min(1.d0,max(-1.d0,2.d0*pxf*pz2/(pxf**2+pz2**2)))
        else
          sinda=min(1.d0,max(-1.d0,dpx*(pxf+pxi)/d))
        endif
        da=asin(sinda)
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
        trans1(1,6)=rhob*(dpzinv+sqphi0/pz1-pxf/pz2*trans1(2,6))
        trans1(1,4)=-pyi*trans1(1,6)
        trans1(5,1)=rhob*trans1(2,1)/pz2
        trans1(5,2)=rhob*(dpzinv-sqphi0/pz2-snphi0*pxi/pz1/pz2)
        trans1(5,6)=rhob*(trans1(2,6)/pz2-dtn/phsq)
        trans1(5,4)=-pyi*trans1(5,6)
        trans1(3,1)=-pyi*trans1(5,1)
        trans1(3,2)=-pyi*trans1(5,2)
        trans1(3,4)=rhob*(phi0n-da)-pyi*trans1(5,4)
        trans1(3,6)=trans1(5,4)
        trans1(5,6)=trans1(5,6)-(phi0n-da)*rhob+h0/h1emit**3*aln
        if(ak1 .eq. 0.d0)then
          trans1(1,3)=0.d0
          trans1(2,3)=0.d0
          trans1(3,3)=1.d0
          trans1(4,3)=0.d0
          trans1(5,3)=0.d0
        else
          t21=-ak1*dxe
          t43= ak1
          trans1(1,1)=trans1(1,1)+trans1(1,2)*t21
          trans1(1,3)=            trans1(1,4)*t43
          trans1(2,1)=trans1(2,1)+trans1(2,2)*t21
          trans1(2,3)=            trans1(2,4)*t43
          trans1(3,1)=trans1(3,1)+trans1(3,2)*t21
          trans1(3,3)=1.d0       +trans1(3,4)*t43
          trans1(4,3)=            trans1(4,4)*t43
          trans1(5,1)=trans1(5,1)+trans1(5,2)*t21
          trans1(5,3)=            trans1(5,4)*t43
        endif
        xf=xi*csphi0+rhoe*(snphi0*pxi-dpx*(pxi+pxf)/(pz1+pz2))
     1       +drho*sqphi0
        if(n .eq. ndiv .and. akn .ne. 0.d0)then
          xfr=xf/rho0
          dxe=1.d0+xfr*(1.d0-xfr*(.5d0-xfr/3.d0))
          t21=-akn*.5d0*dxe
          t43= akn*.5d0
          do 120 k=1,6
            trans1(2,k)=trans1(2,k)+trans1(1,k)*t21
            trans1(4,k)=trans1(4,k)+trans1(3,k)*t43
120       continue
        else
          xfr=0.d0
        endif
        call tmultr5(trans,trans1,irad)
        call tmulbs(beam ,trans1,.true.,.true.)
        if(calpol)then
          if(n .eq. 1)then
            call polpar(21,ld,aln,phi0n,phibn,ak,0.d0,cod)
          else
            call polpar(22,ld,aln,phi0n,phibn,ak,0.d0,cod)
          endif
        endif
        cod(1)=xf
        cod(3)=yi+pyi*rhoe*(phi0n-da)
        if(n .eq. ndiv .and. akn .ne. 0.d0)then
          xe=xf+xf*xfr*(.5d0-xfr*(2.d0-xfr)/12.d0)
          cod(2)=pxf*pr-akn*.5d0*xe
          cod(4)=pyi*pr+akn*cod(3)*.5d0
        else
          cod(2)=pxf*pr
          cod(4)=pyi*pr
        endif
        cod(5)=cod(5)-phi0n*(dp*rhob+drhob)+da*rhoe-dvemit*aln
        if(n .eq. ndiv)then
          if(calpol)then
            if(irad .eq. 6)then
              npelm=npelm+1
            else
              ipelm=ipelm+1
              call tinitr(trans1)
              call tmov(trans1,rlist(ipoltr+(ipelm-1)*36),36)
              call polpar(23,ld,aln,phi0n,phibn,ak,0.d0,cod)
            endif
          endif
          if(enarad)then
            if(fb2 .ne. 0.d0 .and.
     $           mfring .gt. 0 .or. mfring .eq. -2)then
              f1r=fb2
            else
              f1r=0.d0
            endif
            bx=b1*cod(3)
            by=b+b1*xe
            bxy=b1*dxe
            call trade(trans,beam,cod,bx,by,0.d0,0.d0,
     $           0.d0,bxy,0.d0,
     $           f*aln*.5d0-tanp2,fb1,fb2,al,al,aln*.5d0)
          endif
        endif
100   continue
      call tbedge(trans,cod,beam,al,phib,psi2*phi0+apsi2,.false.,ld)
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=-fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          if(fringe)then
            dyfra2=4.d0*dyfr2/fb2**2
          else
            dyfra2=0.d0
          endif
          call tblfre(trans,cod,beam,dxfr2,dyfr2,dyfra2,ld)
        endif
      endif
      if(dtheta .ne. 0.d0)then
        cod(2)=cod(2)+dphix
        cod(4)=cod(4)+dphiy
      endif
      call tchge(trans,cod,beam,dx,dy,-theta,.false.,ld)
      return
      end

      subroutine tbdrifte(trans,cod,beam,al,phi0,
     $     h0,h1emit,dvemit,irad,calpol,ld)
      implicit none
      real*8 trans(6,12),cod(6),beam(42),phi0,al,cp,sp,pr,pxi,pzf,
     $     trans1(6,6),xi,pyi,pzi,pxf,xf,dpzipxi,dpzipyi,dpzip,
     $     dpzfpxi,dpzfpyi,dpzfp,rho0,h0,dl,xsin,dcp,
     $     h1emit,dvemit,a,ptmax
      integer*4 ld,irad
      logical*4 calpol
      parameter (ptmax=0.999d0)
      cp=cos(phi0)
      sp=sin(phi0)
      if(cp .ge. 0.d0)then
        dcp=sp**2/(1.d0+cp)
      else
        dcp=1.d0-cp
      endif
      rho0=al/phi0
      call tdrife(trans,cod,beam,rho0*sp,0.d0,0.d0,0.d0,
     $     .true.,.false.,calpol,irad,ld)
      xi=cod(1)+rho0*dcp
      pr=1.d0+cod(6)
      pxi=cod(2)
      pyi=cod(4)
      a=min(ptmax,pxi**2+pyi**2)
      pzi=pr*sqrt(1.d0-a/pr**2)
c      pzi=sqrt(max(1.d-4,(pr-pxi)*(pr+pxi)-pyi**2))
      pzf=pzi*cp-pxi*sp
      pxf=pzi*sp+pxi*cp
      xf=xi*pzi/pzf
      call tinitr(trans1)
      dpzipxi=-pxi/pzi
      dpzfpxi= cp*dpzipxi-sp
      dpzipyi=-pyi/pzi
      dpzfpyi= cp*dpzipyi
      dpzip  = pr/pzi
      dpzfp  = cp*dpzip
      trans1(1,1)=pzi/pzf
      trans1(1,2)=(dpzipxi-dpzfpxi*pzi/pzf)/pzf*xi
      trans1(1,4)=(dpzipyi-dpzfpyi*pzi/pzf)/pzf*xi
      trans1(1,6)=(dpzip  -dpzfp  *pzi/pzf)/pzf*xi
      trans1(2,2)=cp+sp*dpzipxi
      trans1(2,4)=   sp*dpzipyi
      trans1(2,6)=   sp*dpzip
      trans1(3,1)=sp*pyi/pzf
      trans1(3,2)=-xi*sp*pyi/pzf**2*dpzfpxi
      trans1(3,4)= xi*sp*(1.d0-pyi*dpzfpyi/pzf)/pzf
      trans1(3,6)=-xi*sp*pyi/pzf**2*dpzfp
      trans1(5,1)=-pr/pzf*sp
      trans1(5,2)= pr/pzf**2*xi*sp*dpzfpxi
      trans1(5,4)= pr/pzf**2*xi*sp*dpzfpyi
      dl=rho0*xsin(phi0)
      trans1(5,6)=-xi*sp*(1.d0-pr*dpzfp  /pzf)/pzf
     $     +h0/h1emit**3*dl
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      if(calpol)then
        call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
      endif
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
      real*8 trans(4,5),cod(6),transe(6,12),beam(42),
     $     dx,dy,theta,fb1,fb2,al0,phib,phi0,
     $     psi1,psi2,ak,dtheta,eps0,apsi1,apsi2
      logical*4 coup,fringe
      call tinitr(transe)
      call tbende(transe,cod,beam,al0,phib,phi0,
     $     psi1,psi2,apsi1,apsi2,ak,
     1     dx,dy,theta,dtheta,
     $     fb1,fb2,mfring,fringe,eps0,.false.,.true.,1)
      call qcopymat(trans,transe,.false.)
c      write(*,*)'qbend ',cod,al0,phi0,phib
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end

      subroutine tblfre(trans,cod,beam,dxfr,dyfr,dyfra,ld)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 ld
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
      if(calpol)then
        call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
      endif
      cod(1)=cod(1)+dxfr*cod(6)/pr
      cod(4)=cod(4)+(dyfr-dyfraysq)*cod(3)/pr
      cod(5)=cod(5)+(dxfr*cod(2)+
     $     (.5d0*dyfr-.25d0*dyfraysq)*ysq)/pr**2
      return
      end

      subroutine tbthie(trans,cod,beam,phib,phi0,
     1                 dx,dy,theta,ld)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 ld
      real*8 trans(6,12),cod(6),beam(42),phib,phi0,dx,dy,theta,
     $     trans1
      common /tem/ trans1(6,13)
      call tchge(trans,cod,beam,-dx,-dy,theta,.true.,ld)
      call tinitr(trans1)
      trans1(2,6)=phi0
      trans1(5,1)=-phi0
      call tmultr(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.,.true.)
      if(calpol)then
        call polpar(20,ld,0.d0,phi0,phib,0.d0,0.d0,cod)
      endif
      cod(2)=cod(2)+(phi0-phib)+cod(6)*phi0
      cod(5)=cod(5)-phi0*cod(1)
      call tchge(trans,cod,beam,dx,dy,-theta,.false.,ld)
      return
      end
