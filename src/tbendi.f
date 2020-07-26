      subroutine tbendi(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb1,fb2,mfring,enarad,fringe,eps0)
      use tbendcom, only:tbrot,tbshift
      use bendib
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      use kradlib
      use tspin, only:cphi0,sphi0
      use photontable
      implicit none
      integer*4 , parameter :: ndivmax=1000
      integer*4 np,mfring,i,ndiv,n,n1,n2
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np),
     $     alx(-1:ndivmax+2),alr(-1:ndivmax+2),
     $     akxn(-1:ndivmax+2),phixn(-1:ndivmax+2),
     $     cphixn(-1:ndivmax+2),sphixn(-1:ndivmax+2),
     $     al,phib,phi0,cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,cost,sint,fb1,fb2,eps0,
     $     tanp1,tanp2,aind,b,dxfr1,dyfr1,dyfra1,pr,eps,
     $     af,f,fpx,ff,aln,f1r,f2r,akx0,alx0,
     $     phix0,alc,
     $     dxfr2,dyfr2,dyfra2,dtheta
      logical*4 enarad,fringe,krad
      call tbshift(np,x,px,y,py,z,dx,dy,phi0,cost,sint,.true.)
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      rho0=al/phi0
      aind=rho0/phi0*ak
      if(eps0 .le. 0.d0)then
        eps=epsbend
      else
        eps=epsbend*eps0
      endif
      ndiv=1+int(abs(phi0/eps))
      krad=rad .and. enarad .and. al .ne. 0.d0
      f1r=0.d0
      f2r=0.d0
      n1=1
      n2=0
      alc=al
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
        if(calpol)then
          bsi=0.d0
        endif
        if(iprev(l_track) .eq. 0)then
          f1r=.5d0*fb1
          n1=-1
        endif
        if(inext(l_track) .eq. 0)then
          f2r=.5d0*fb2
          n2=2
        endif
        b=brhoz/rhob
        alc=al-f1r-f2r
        ndiv=max(ndiv,ndivrad(phib*alc/al,ak*alc/al,0.d0,eps0))
      endif
      ndiv=min(ndivmax,ndiv)
      n2=ndiv+n2
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 2)then
        call ttfrin(np,x,px,y,py,z,g,4,ak,al,0.d0)
      endif
      if(fb1 .ne. 0.d0)then
        dxfr1=fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        dyfra1=4.d0*dyfr1/fb1**2
        do i=1,np
c     dp=g(i)*(2.d0+g(i))
          dp=g(i)
          pr=1.d0+dp
          x(i)=x(i)+dxfr1*dp/pr
          py(i)=py(i)+(dyfr1-dyfra1*y(i)**2)*y(i)/pr**2
          z(i)=z(i)+(dxfr1*px(i)+
     $         (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/pr)/pr
        enddo
      endif
      if(fringe)then
        af=1.d0
      else
        af=0.d0
      endif
      aln=alc/ndiv
      do n=n1,n2
        call tbendal(n,ndiv,f1r,f2r,aln,alx(n),alr(n))
        akxn(n)=ak*alx(n)/al
        phixn(n)=phi0*alx(n)/al
        cphixn(n)=cos(phixn(n))
        sphixn(n)=sin(phixn(n))
      enddo
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
        akx0=0.d0
        alx0=0.d0
        phix0=0.d0
        if(photons .and. krad)then
          call tsetpcvt(l_track,dx,dy,theta,dtheta,phi0,al)
        endif
        do n=n1,n2
          call tbendiinit(akxn(n),alx(n),n .eq. n1)
          call tbendicorr((akxn(n)+akx0)*.5d0,(alx(n)+alx0)*.5d0,
     $         (phixn(n)+phix0)*.5d0)
          akx0=akxn(n)
          alx0=alx(n)
          phix0=phixn(n)
          call tbendibody(alx(n))
          if(krad)then
            bsi(i)=bsi(i)+akxn(n)/alx(n)*xi*yi
            if(n .ne. n2)then
              cphi0=cphixn(n)
              sphi0=sphixn(n)
              if(rfluct)then
                call tradkf1(xi,pxi,yi,pyi,zi,dp,dv(i),
     $               sx(i),sy(i),sz(i),
     $               pxr0(i),pyr0(i),zr0(i),bsi(i),alr(n),i)
              else
                call tradk1(xi,pxi,yi,pyi,zi,dp,dv(i),
     $               sx(i),sy(i),sz(i),
     $               pxr0(i),pyr0(i),zr0(i),bsi(i),alr(n))
              endif
              pxr0(i)=pxi
              pyr0(i)=pyi
              zr0(i)=zi
              bsi(i)=0.d0
              if(photons)then
                pcvt%fr0=pcvt%fr0+alx(n)/al
c                call tsetphotongeo(alx(n),phixn(n),theta,.false.)
              endif
            endif
          endif
        enddo
        call tbendicorr(akx0*.5d0,alx0*.5d0,phix0*.5d0)
        zi=zi-dv(i)*al
        px(i)=pxi+tanp2*xi/rhoe
        y(i)=yi
        f=yi/rhoe
        fpx=af*px(i)
        py(i)=pyi-(tanp2-fpx)*f
        ff  =af*yi*f*.5d0
        x(i)=xi-ff
        z(i)=zi+ff*fpx
        g(i)=dp
c        write(*,'(a,1p6g15.7)')'tbendi-2 ',zi,z(i),xi,x(i),ff,fpx
      enddo
      if(fb2 .ne. 0.d0)then
        dxfr2=fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        dyfra2=4.d0*dyfr2/fb2**2
        do i=1,np
          dp=g(i)
          pr=1.d0+dp
          x(i)=x(i)-dxfr2*dp/pr
          py(i)=py(i)+(dyfr2-dyfra2*y(i)**2)*y(i)/pr**2
          z(i)=z(i)-(dxfr2*px(i)-
     $         (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/pr)/pr
        enddo
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 1)then
        call ttfrin(np,x,px,y,py,z,g,4,-ak,al,0.d0)
      endif
      if(krad)then
        if(photons)then
          call tsetpcvt(l_track,dx,dy,theta,dtheta,phi0,al)
          pcvt%fr0=1.d0-alr(n2)/al
        endif
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,alr(n2),phixn(n2))
      endif
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      call tbshift(np,x,px,y,py,z,-dx,-dy,-phi0,cost,-sint,.false.)
      return
      end
