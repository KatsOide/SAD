      subroutine tmulta(
     $     np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak0,phi,
     $     psi1,psi2,bz,
     1     dx,dy,theta,dtheta,
     $     eps0,enarad,fb1,fb2,mfring,fringe)
      use ffs_flag, only:rad,ndivrad,calpol
      use tmacro
      use multa
      use tbendcom, only:tbrot,tbshift
      use tspin, only:tradke,pxr0,pyr0,zr0,bsi
      use photontable,only:tsetphotongeo,pp
      use mathfun
      implicit none
      integer*4 ndivmax
      real*8 ampmax,eps00
      parameter (ampmax=0.05d0,eps00=0.005d0,ndivmax=2000)
      integer*4 np,mfring,i,n,mfr,ndiv,nmmax,m,m1,k,nmmin
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     al,phi,psi1,psi2,bz,dx,dy,theta,eps0,fb1,fb2,
     $     dtheta,pr,cost,sint,rho0,rhob,
     $     sinp1,sinp2,cosp1,cosp2,phin,aln,cosw,sinw,sqwh,sinwp1,
     $     eps,w,r,rk(0:nmult),als,ak0r,ak1r,ak1n,
     $     phib,phibn,an(0:nmult+1)
      real*8 sx(np),sy(np),sz(np)
      complex*16 ak0(0:nmult),ak(0:nmult),akn(0:nmult),
     $     cx1,csl,csr,cl,cr,cg,cx
      logical*4 enarad,fringe,krad
      real*8 fact(0:nmult+1)
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2432902008176640000.d0,51090942171709440000.d0,
     $     1124000727777607680000.d0/
      data an/1.d0,1.d0,
     $0.5d0,
     $0.33333333333333333333d0,
     $0.25d0,
     $0.2d0,
     $0.166666666666666666667d0,
     $0.142857142857142857143d0,
     $0.125d0,
     $0.111111111111111111111d0,
     $0.1d0,
     $0.090909090909090909091d0,
     $0.083333333333333333333d0,
     $0.076923076923076923077d0,
     $0.071428571428571428571d0,
     $0.066666666666666666667d0,
     $0.0625d0,
     $0.058823529411764705882d0,
     $0.055555555555555555556d0,
     $0.052631578947368421053d0,
     $0.05d0,
     $0.047619047619047619048d0,
     $0.045454545454545454545d0/
      if(bz .ne. 0.d0)then
        write(*,*)
     $       'MULT with nonzero ANGLE and BZ is not yet supported.'
        call abort
      endif
      if(gknini)then
        call gkninit
      endif
      cost=cos(theta)
      sint=sin(theta)
      call tbshift(np,x,px,y,py,z,dx,dy,phi,cost,sint,.true.)
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi,dtheta)
      endif
      if(eps0 .eq. 0.d0)then
        eps=eps00
      else
        eps=eps00*eps0
      endif
      nmmax=0
      ak=ak0
      ak0r=dble(ak(0))
      ak1r=dble(ak(1))
      ak(0)=dcmplx(0.d0,imag(ak(0)))
      ak(1)=dcmplx(0.d0,imag(ak(1)))
      do n=nmult,0,-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmax=n
          exit
        endif
      enddo
      nmmin=nmmax
      do n=0,nmmax-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmin=n
          exit
        endif
      enddo
      ndiv=1
      do n=nmmin,nmmax
        ndiv=max(ndiv,
     $int(sqrt(ampmax**(n-1)/6.d0/fact(n-1)/eps*abs(ak(n)*al)))+1)
      enddo
      phib=phi+ak0r
      krad=rad .and. enarad
      if(krad)then
        ndiv=max(ndiv,ndivrad(phib,ak1r,0.d0,eps0))
      endif
      ndiv=min(ndivmax,ndiv)
      aln=al/ndiv
      if(fb1 .ne. 0.d0)then
        aln=aln-(phi*fb1)**2/al/48.d0
     1       *sin(.5d0*(phi-psi1-psi2))/sin(.5d0*phi)/ndiv
      endif
      if(fb2 .ne. 0.d0)then
        aln=aln-(phi*fb2)**2/al/48.d0
     1       *sin(.5d0*(phi-psi1-psi2))/sin(.5d0*phi)/ndiv
      endif
      rho0=al/phi
      rhob=al/phib
      phin=phi/ndiv
      phibn=phib/ndiv
      ak1n=ak1r/ndiv
      sinp1=sin(psi1)
      cosp1=cos(psi1)
      sinp2=sin(psi2)
      cosp2=cos(psi2)
      do m=nmmin,nmmax
        akn(m)=ak(m)/(fact(m+1)*ndiv)
      enddo
      als=aln*.5d0
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
        if(calpol)then
          bsi=0.d0
        endif
        pp%theta=theta
      endif
      do n=1,ndiv
        if(n .eq. 1)then
          w=phin*.5d0-psi1
          cosw=cos(w)
          sinw=sin(w)
          if(cosw .ge. 0.d0)then
            sqwh=sinw**2/(1.d0+cosw)
          else
            sqwh=1.d0-cosw
          endif
          sinwp1=sin(phin*.5d0)
          mfr=0
          if(mfring .eq. 2)then
            mfr=0
          elseif(mfring .ne. 0)then
            mfr=-1
          endif
          if(krad)then
            do i=1,np
              cx1=dcmplx(x(i),y(i))
              cx=0.d0
              do k=nmmax,2,-1
                cx=(cx+ak(k))*cx1*an(k+1)
              enddo
              cx=.5d0*cx*cx1**2
              bsi(i)=imag(cx)/al
            enddo
          endif
          call tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         als,phibn*.5d0,phin*.5d0,psi1,0.d0,
     1         cosp1,sinp1,1.d0,0.d0,
     1         ak1n,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $         fb1,fb2,mfr,fringe,cosw,sinw,sqwh,sinwp1,
     1         enarad,eps0,.false.,1)
          w=phin
          cosw=cos(w)
          sinw=sin(w)
          if(cosw .ge. 0.d0)then
            sqwh=sinw**2/(1.d0+cosw)
          else
            sqwh=1.d0-cosw
          endif
          sinwp1=sinw
        else
          call tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         aln,phibn,phin,0.d0,0.d0,
     1         1.d0,0.d0,1.d0,0.d0,
     1         ak1n,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $         0.d0,0.d0,0,.false.,cosw,sinw,sqwh,sinwp1,
     1         enarad,eps0,.false.,2)
          als=als+aln
        endif
        do i=1,np
          pr=(1.d0+g(i))
          cx1=dcmplx(x(i),y(i))
          csl=(0.d0,0.d0)
          csr=(0.d0,0.d0)
          r=rho0+x(i)
          rk(0)=1.d0+sqrt1(x(i)/rho0)
          do k=1,nmult-nmmin
            rk(k)=rk(k-1)/r
          enddo
          do m=nmult,0,-1
            cl=(0.d0,0.d0)
            cr=(0.d0,0.d0)
            do k=max(0,m-nmmax),m-nmmin
              m1=m-k
              cg=gkn(m1,k)*rk(k)*akn(m1)
              cl=cl+(m+1)*cg
              cr=cr+(.5d0-k)*cg
            enddo
            csl=csl*cx1+cl
            csr=csr*cx1+cr
          enddo
          px(i)=px(i)-(dble(csr*cx1)/r+dble(csl))/pr
          py(i)=py(i)+imag(csl)/pr
        enddo
        if(krad)then
          pxr0=px
          pyr0=py
          zr0=z
        endif
      enddo
      w=phin*.5d0-psi2
      cosw=cos(w)
      sinw=sin(w)
      if(cosw .ge. 0.d0)then
        sqwh=sinw**2/(1.d0+cosw)
      else
        sqwh=1.d0-cosw
      endif
      sinwp1=sinw
      mfr=0
      if(mfring .eq. 1)then
        mfr=0
      elseif(mfring .ne. 0)then
        mfr=-2
      endif
      if(krad)then
        do i=1,np
          cx1=dcmplx(x(i),y(i))
          cx=0.d0
          do k=nmmax,2,-1
            cx=(cx+ak(k))*cx1*an(k+1)
          enddo
          cx=.5d0*cx*cx1**2
          bsi(i)=-imag(cx)/al
        enddo
      endif
      call tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     aln*.5d0,phibn*.5d0,phin*.5d0,0.d0,psi2,
     1     1.d0,0.d0,cosp2,sinp2,
     1     ak1n,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $     fb1,fb2,mfr,fringe,cosw,sinw,sqwh,sinwp1,
     1     enarad,eps0,.false.,2)
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi,-dtheta)
      endif
      call tbshift(np,x,px,y,py,z,-dx,-dy,-phi,cost,-sint,.false.)
      return
      end

      subroutine tmultae(trans,cod,beam,srot,al,ak0,
     $     phi,psi1,psi2,apsi1,apsi2,bz,
     1     dx,dy,theta,dtheta,
     $     eps0,enarad,fringe,fb1,fb2,mfring,l)
      use tfstk
      use tmacro
      use multa
      use tspin, only:tradke      
      implicit none
      integer*4 ndivmax
      real*8 ampmax,eps00
      parameter (ampmax=0.05d0,eps00=0.005d0,ndivmax=2000)
      integer*4 mfring,n,mfr,ndiv,nmmax,m,m1,k,nmmin,l
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6),srot(3,9),
     $     al,phi,psi1,psi2,bz,dx,dy,theta,eps0,fb1,fb2,
     $     dphix,dphiy,dtheta,rho0,
     $     phin,aln,eps,r,rk(0:nmult),apsi1,apsi2,
     $     psi1n,psi2n,ak0r,ak1r,ak1n,phib,phibn
      complex*16 ak(0:nmult),ak0(0:nmult),akn(0:nmult),
     $     cx1,csl,csr,cl,cr,cg,
     $     csxx,csxy,csyy,cxx,cxy,cyy
      logical*4 enarad,fringe
      real*8 fact(0:nmult+1)
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2432902008176640000.d0,51090942171709440000.d0,
     $     1124000727777607680000.d0/
      if(bz .ne. 0.d0)then
        write(*,*)
     $       'MULT with nonzero ANGLE and BZ is not yet supported.'
        call abort
      endif
      if(gknini)then
        call gkninit
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,dtheta,phi,.true.)
      if(dtheta .ne. 0.d0)then
        dphix=      phi*sin(.5d0*dtheta)**2
        dphiy= .5d0*phi*sin(dtheta)
        cod(2)=cod(2)+dphix
        cod(4)=cod(4)+dphiy
      else
        dphix=0.d0
        dphiy=0.d0
      endif
      if(eps0 .eq. 0.d0)then
        eps=eps00
      else
        eps=eps00*eps0
      endif
      ak=ak0
      ak0r=dble(ak(0))
      ak1r=dble(ak(1))
      ak(0)=dcmplx(0.d0,imag(ak(0)))
      ak(1)=dcmplx(0.d0,imag(ak(1)))
      nmmax=0
      do n=nmult,0,-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmax=n
          exit
        endif
      enddo
      nmmin=nmmax
      do n=0,nmmax-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmin=n
          exit
        endif
      enddo
      ndiv=1
      do n=nmmin,nmmax
        ndiv=max(ndiv,
     $int(sqrt(ampmax**(n-1)/6.d0/fact(n-1)/eps*abs(ak(n)*al)))+1)
      enddo
      ndiv=min(ndivmax,ndiv)
      psi1n=2.d0*psi1*ndiv
      psi2n=2.d0*psi2*ndiv
      aln=al/ndiv
      if(fb1 .ne. 0.d0 .and. (mfring .eq. 1 .or. mfring .eq. 3))then
        aln=aln-(phi*fb1)**2/al/48.d0
     1       *sin(.5d0*(phi*(1.d0-psi1-psi2)-apsi1-apsi2))
     $       /sin(.5d0*phi)/ndiv
      endif
      if(fb2 .ne. 0.d0 .and. (mfring .eq. 2 .or. mfring .eq. 3))then
        aln=aln-(phi*fb2)**2/al/48.d0
     1       *sin(.5d0*(phi*(1.d0-psi1-psi2)-apsi1-apsi2))
     $       /sin(.5d0*phi)/ndiv
      endif
      phib=phi+ak0r
      phin=phi/ndiv
      phibn=phib/ndiv
      ak1n=ak1r/ndiv
      rho0=aln/phin
      do m=nmmin,nmmax
        akn(m)=ak(m)/(fact(m+1)*ndiv)
      enddo
      call tinitr(trans1)
      do n=1,ndiv
        if(n .eq. 1)then
          mfr=0
          if(mfring .eq. 2)then
            mfr=0
          elseif(mfring .ne. 0)then
            mfr=-1
          endif
          call tbende(trans,cod,beam,srot,aln*.5d0,phibn*.5d0,phin*.5d0,
     $         psi1n,0.d0,apsi1,0.d0,ak1n*.5d0,
     $         0.d0,0.d0,0.d0,0.d0,
     1         fb1,fb2,mfr,fringe,eps0,enarad,.false.,.false.,l)
        else
          call tbende(trans,cod,beam,srot,aln,phibn,phin,
     $         0.d0,0.d0,0.d0,0.d0,ak1n,
     $         0.d0,0.d0,0.d0,0.d0,
     1         0.d0,0.d0,0,.false.,eps0,enarad,.false.,.false.,l)
        endif
        cx1=dcmplx(cod(1),cod(3))
        csl=(0.d0,0.d0)
        csr=(0.d0,0.d0)
        csxx=(0.d0,0.d0)
        csxy=(0.d0,0.d0)
        csyy=(0.d0,0.d0)
        r=rho0+cod(1)
        rk(0)=sqrt(r/rho0)
        do k=1,nmult-nmmin
          rk(k)=rk(k-1)/r
        enddo
        do m=nmult,0,-1
          cl=(0.d0,0.d0)
          cr=(0.d0,0.d0)
          cxx=(0.d0,0.d0)
          cxy=(0.d0,0.d0)
          cyy=(0.d0,0.d0)
          do k=max(0,m-nmmax),m-nmmin
            m1=m-k
            cg=gkn(m1,k)*rk(k)*akn(m1)
            cl=cl+(m+1)*cg
            cr=cr+(.5-k)*cg
            cxx=cxx+(.25d0*(4*k**2-1)*cx1-(2*k-1)*(m+1)*r)*cg
            if(m .eq. 0)then
              cxy=cxy+(1-2*k)*cg
            else
              cxy=cxy+(m+1)*((1-2*k)*cx1+2*m*r)*cg
              cyy=cyy+m*(m+1)*cg
            endif
          enddo
          csl=csl*cx1+cl
          csr=csr*cx1+cr
          csxx=csxx*cx1+cxx
          if(m .ne. 0)then
            csxy=csxy*cx1+cxy
            csyy=csyy*cx1+cyy
          else
            csxy=csxy+cxy
          endif
          if(enarad)then
            call tradke(trans,cod,beam,srot,aln,0.d0,0.d0)
          endif
        enddo
c        write(*,*)'tmultae ',dble(csr*cx1)/r,dble(csl),nmmin
        cod(2)=cod(2)-(dble(csr*cx1)/r+dble(csl))
        cod(4)=cod(4)+imag(csl)
        trans1(2,1)=-(dble(csxx)/r**2+dble(csyy))
        trans1(2,3)=imag(csxy)*.5d0/r
        trans1(4,1)=trans1(2,3)
        trans1(4,3)=dble(csyy)
        call tmultr5(trans,trans1,irad)
        call tmulbs(beam,trans1,.true.)
      enddo
      mfr=0
      if(mfring .eq. 1)then
        mfr=0
      elseif(mfring .ne. 0)then
        mfr=-2
      endif
      call tbende(trans,cod,beam,srot,aln*.5d0,phibn*.5d0,phin*.5d0,
     $     0.d0,psi2n,0.d0,apsi2,ak1n*.5d0,
     $     0.d0,0.d0,0.d0,0.d0,
     1     fb1,fb2,mfr,fringe,eps0,enarad,.false.,.false.,l)
      if(dtheta .ne. 0.d0)then
        cod(2)=cod(2)+dphix
        cod(4)=cod(4)+dphiy
      endif
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,-dtheta,-phi,.false.)
      return
      end
