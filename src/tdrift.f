c     common module for drift element implementation
      module element_drift_common
      implicit none

c     Maximum amplitude of (px/p0)^2 + (py/p0)^2
      real(8), public, parameter :: ampmax = 0.9999d0

      contains
c     drift in the free space
      subroutine tdrift_free(np,x,px,y,py,z,dv,al)
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 ,intent(in):: np
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np)
      real*8 ,intent(in):: al
      real*8 al1,dpz
      integer*4 i
      do concurrent (i=1:np)
        dpz=pxy2dpz(px(i),py(i))
        al1=al/(1.d0+dpz)
        x(i)=x(i)+px(i)*al1
        y(i)=y(i)+py(i)*al1
        z(i)=z(i)+dpz*al1-dv(i)*al
      enddo
c      real*8 al1(np),dpz(np)
c      dpz=pxy2dpz(px,py)
c      al1=al/(1.d0+dpz)
c      x=x+px*al1
c      y=y+py*al1
c      z=z+dpz  *al1-dv*al
      return
      end

c     drift in the parallel solenoid
      subroutine tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,enarad)
      use ffs_flag, only:rfluct,photons,calpol
      use photontable
      use kradlib
      use tspin, only:cphi0,sphi0
      use mathfun, only: sqrtl,pxy2dpz,xsincos
      use tmacro, only:l_track
      implicit none
      integer*4 ,intent(in):: np
      real*8 ,intent(in):: al,bz
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     g(np),dv(np),sx(np),sy(np),sz(np)
      logical*4 ,intent(in):: enarad
      real*8 px1,py1
      integer*4 i
      real*8 pr,bzp,pxi,pyi
      real*8 dpzi,pzi,al1,zi
      real*8 phi,a22,a24,a12,a14
      if(enarad)then
        if(calpol)then
          bsi=0.d0
        endif
        cphi0=1.d0
        sphi0=0.d0
        if(rfluct .and. photons)then
          call tsetpcvt(l_track,0.d0,0.d0,0.d0,0.d0,0.d0,al)
        endif
      endif
      do i=1,np
         pr=(1.d0+g(i))
         bzp=bz/pr
         pxi=px(i)+bzp*y(i)*.5d0
         pyi=py(i)-bzp*x(i)*.5d0
         zi=z(i)

         dpzi=pxy2dpz(pxi,pyi)
         pzi=1.d0+dpzi
         al1=al/pzi

c         Definition of a* coefficient
c         phi := al1 * bzp
c         a22 := cos(phi)
c         a24 := sin(phi)
c         a12 := (sin(phi)       / phi) * al1
c         a14 := ((1 - cos(phi)) / phi) * al1
c              = (sin(phi)       / phi) * al1 * (sin(phi) / (1 + cos(phi)))
c
c         sin(phi) / phi = \sum_{i=0}^{N-1} (-1)^i \frac{x^{2i}}{(2i-1)!}
c                          + R_2N(x) / x
c         |R_2N(x) / x| =< \frac{x^{2N}}{(2N)!}
c
         phi=al1*bzp
         call xsincos(phi,a24,a12,a22,a14)
c         a22=cos(phi)
c         a24=sin(phi)
c         if(abs(phi) .gt. 1.d-4)then
c            a12=a24/bzp
c         else
c            a12=(1.d0-(1.d0/6.d0)*phi*phi)*al1
c         endif
c         if(a22 .ge. 0.d0)then
c            a14=a12*a24/(1.d0+a22)
c         else
c            a14=(1.d0-a22)/bzp
c         endif
         x(i) =x(i)+(a24*pxi-a14*pyi)/bzp
         y(i) =y(i)+(a14*pxi+a24*pyi)/bzp
         z(i) =z(i)+dpzi *al1-dv(i)*al
         px1=     a22*pxi+a24*pyi
         py1=    -a24*pxi+a22*pyi
         if(enarad)then
           bsi(i)=bz
           if(rfluct)then
             call tradkf1(x(i),px1,y(i),py1,z(i),g(i),dv(i),
     $            sx(i),sy(i),sz(i),
     $            pxi,pyi,zi,bsi(i),al,i)
           else
             call tradk1(x(i),px1,y(i),py1,z(i),g(i),dv(i),
     $            sx(i),sy(i),sz(i),
     $            pxi,pyi,zi,bsi(i),al)
           endif
         endif
         px(i)=px1-bzp*y(i)*.5d0
         py(i)=py1+bzp*x(i)*.5d0
      enddo
c      write(*,'(a,106g15.7)')'td_sol ',x(1),px(1),y(1),py(1),z(1),g(1)
      return
      end

      subroutine tsolconv(pz0,plz,pbz,bpr,
     $     phi,sinphi,xsinphi,cosphi,dcosphi,ndiag)
      use mathfun
      implicit none
      integer*4 ,parameter :: itmax=20
      real*8 ,parameter :: conv=1.d-15
      real*8 ,intent(in):: pz0,plz,pbz,bpr
      real*8 ,intent(out):: phi,sinphi,xsinphi,cosphi,dcosphi
      real*8 dphi,s,u,phi0
      integer*4 ,intent(inout):: ndiag
      integer*4 j
      phi=asinz((bpr-pbz)/hypot(pz0,pbz))+atan(pbz,pz0)
      if(plz .eq. 0.d0)then
        call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
      else
        do j=1,itmax
          call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
          s=plz*xsinphi+pz0*sinphi-pbz*dcosphi
          u=-plz*dcosphi+pz0*cosphi+pbz*sinphi
          dphi=merge((bpr-s)/u,(bpr-s)/pz0,u .ne. 0.d0)
          phi0=phi
          phi=phi+dphi
          if(phi0 .eq. phi .or. abs(dphi) .le. conv*abs(phi))then
            return
          endif
        enddo
        if(ndiag .ge. 0)then
          ndiag=ndiag-1
          write(*,'(a,1p8g13.5)')'tsolconv convergence error: ',
     $         phi,dphi,bpr,s,u,plz,pz0,pbz
          if(ndiag .eq. -1)then
            write(*,*)
     $           'Further tsolconv messages will be suppressed.'
          endif
        endif
      endif
      return
      end

      subroutine tsoldz(trans,cod,al,bxs0,bys0,bzs0,drift)
      use mathfun
      implicit none
      integer*4 ,save::ndiag=15
      integer*4 ,parameter ::itmax=15
      real*8 ,intent(inout):: trans(6,6),cod(6)
      real*8 ,intent(in):: al,bxs0,bys0,bzs0
      real*8 bxs,bys,bzs,pxi,pyi,pz0,
     $     dpz0dpx,dpz0dpy,dpz0dp,phi,
     $     dphidz,dphidpx,dphidpy,dphidp,a24,a12,a22,a14,
     $     da12,da14,pr,dpz0,dv,dvdp,phix,phiy,phiz,babs,
     $     alb,pbx,pby,pbz,pl,dpl,dphizsq,a,r,
     $     dpldpx,dpldpy,dpldp,dplz,plx,ply,plz,ptx,pty,ptz,
     $     cosphi,sinphi,dcosphi,dphi,
     $     xsinphi,ax,ay,az,cx,cy,albabs,ap,bpr,db
      real*8 ,parameter ::conv=1.d-15,bzthre=1.d-20,ptmax=0.9999d0
      logical*4 ,intent(in):: drift
      bxs=bxs0
      bys=bys0
      bzs=bzs0
      babs=norm2([bzs,bxs,bys])
      if(abs(babs) .lt. bzthre)then
        bxs=0.d0
        bys=0.d0
        bzs=0.d0
        babs=0.d0
      endif
      call tinitr(trans)
      pr=1.d0+cod(6)
c cod does NOT have canonical momenta!
      pxi=cod(2)
      pyi=cod(4)
      a=pxi**2+pyi**2
      dpz0=-a/pr/(1.d0+sqrtl(1.d0-a/pr**2))
      pz0=pr+dpz0
      r=al/pz0
      dpz0dpx= -pxi/pz0
      dpz0dpy= -pyi/pz0
      dpz0dp =   pr/pz0
      if(bxs .eq. 0.d0 .and. bys .eq. 0.d0)then
        phi=bzs*r
        dphidz  = 1.d0/pz0
        dphidpx = -r/pz0*dpz0dpx
        dphidpy = -r/pz0*dpz0dpy
        dphidp  = -r/pz0*dpz0dp
        if(bzs .eq. 0.d0)then
          a24=0.d0
          a12=r
          a22=1.d0
          a14=0.d0
          da12=1.d0
          da14=0.d0
        else
          a24=sin(phi)
          a12=a24/bzs
          a22=cos(phi)
          a14=merge(a24**2/(1.d0+a22),1.d0-a22,a22 .eq. 0.d0)/bzs
          da12=a22
          da14=a24
        endif
        cod(1)=cod(1)+a12*pxi+a14*pyi
        cod(3)=cod(3)-a14*pxi+a12*pyi
        cod(2)=       a22*pxi+a24*pyi
        cod(4)=      -a24*pxi+a22*pyi
        trans(1,2)=( da12*pxi+da14*pyi)*dphidpx+a12
        trans(1,4)=( da12*pxi+da14*pyi)*dphidpy+a14
        trans(1,6)=( da12*pxi+da14*pyi)*dphidp
        trans(3,2)=(-da14*pxi+da12*pyi)*dphidpx-a14
        trans(3,4)=(-da14*pxi+da12*pyi)*dphidpy+a12
        trans(3,6)=(-da14*pxi+da12*pyi)*dphidp
        trans(2,2)=( -a24*pxi+ a22*pyi)*bzs*dphidpx+a22
        trans(2,4)=( -a24*pxi+ a22*pyi)*bzs*dphidpy+a24
        trans(2,6)=( -a24*pxi+ a22*pyi)*bzs*dphidp
        trans(4,2)=( -a22*pxi- a24*pyi)*bzs*dphidpx-a24
        trans(4,4)=( -a22*pxi- a24*pyi)*bzs*dphidpy+a22
        trans(4,6)=( -a22*pxi- a24*pyi)*bzs*dphidp
        trans(5,2)=r*pr/pz0*dpz0dpx
        trans(5,4)=r*pr/pz0*dpz0dpy
        if(drift)then
          call tgetdv(cod(6),dv,dvdp)
          cod(5)=cod(5)+al*(dpz0/pz0-dv)
          trans(5,6)=al*(a/pz0**3+dvdp)
        else
          trans(1,5)=( da12*pxi+da14*pyi)*dphidz
          trans(3,5)=(-da14*pxi+da12*pyi)*dphidz
          trans(2,5)=( -a24*pxi+ a22*pyi)*bzs*dphidz
          trans(4,5)=( -a22*pxi- a24*pyi)*bzs*dphidz
          cod(5)=cod(5)-r*pr
          trans(5,5)=-pr/pz0
          trans(5,6)= r*a/pz0**2
        endif
      else
        phix=bxs/babs
        phiy=bys/babs
        phiz=bzs/babs
        alb=1.d0/babs
        albabs=al*babs
        dphizsq=phix**2+phiy**2
        dpl=pxi*phix+pyi*phiy+dpz0*phiz
        pl=pr*phiz+dpl
        dpldpx=phix+phiz*dpz0dpx
        dpldpy=phiy+phiz*dpz0dpy
        dpldp =     phiz*dpz0dp
        plx=pl*phix
        ply=pl*phiy
        plz=pl*phiz
        ptx=pxi-plx
        pty=pyi-ply
        ptz=dpz0 -dpl*phiz+pr*dphizsq
        pbx=pty*phiz-ptz*phiy
        pby=ptz*phix-ptx*phiz
        pbz=ptx*phiy-pty*phix
        bpr=albabs/pr
        db=bpr-pbz
        ap=hypot(pz0,pbz)
        dphi=-atan(pbz,pz0)
        phi=asinz(db/ap)-dphi
        if(al .ne. 0.d0)then
          if(plz .eq. 0.d0)then
            call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
          else
            call tsolconv(pz0,plz,pbz,bpr,
     $           phi,sinphi,xsinphi,cosphi,dcosphi,ndiag)
          endif
        else
          phi=0.d0
          xsinphi=0.d0
          sinphi=0.d0
          dcosphi=0.d0
          cosphi=1.d0
        endif
        dplz=-pr*dphizsq+dpl*phiz
        ax=pxi+ptx*dcosphi+pbx*sinphi
        ay=pyi+pty*dcosphi+pby*sinphi
        az=pz0+ptz*dcosphi+pbz*sinphi
        dphidpx=-(dpldpx*phiz*xsinphi+dpz0dpx*sinphi
     $       -phiy*dcosphi)/az
        dphidpy=-(dpldpy*phiz*xsinphi+dpz0dpy*sinphi
     $       +phix*dcosphi)/az
        dphidp =-(dpldp *phiz*xsinphi+dpz0dp*sinphi)/az
        dphidz =babs/az
        cod(1)=cod(1)+(plx*phi+ptx*sinphi-pbx*dcosphi)*alb
        cod(3)=cod(3)+(ply*phi+pty*sinphi-pby*dcosphi)*alb
        cod(2)=pxi+ptx*dcosphi+pbx*sinphi
        cod(4)=pyi+pty*dcosphi+pby*sinphi
c cod does NOT have canonical momenta!
        trans(1,2)=alb*(dpldpx*phix*xsinphi+sinphi
     $       +dpz0dpx*phiy       *dcosphi+ax*dphidpx)
        trans(1,4)=alb*(dpldpy*phix*xsinphi
     $       -(phiz-dpz0dpy*phiy)*dcosphi+ax*dphidpy)
        trans(1,6)=alb*(dpldp *phix*xsinphi
     $       +dpz0dp *phiy       *dcosphi+ax*dphidp )
        trans(3,2)=alb*(dpldpx*phiy*xsinphi
     $       -(dpz0dpx*phix-phiz)*dcosphi+ay*dphidpx)
        trans(3,4)=alb*(dpldpy*phiy*xsinphi+sinphi
     $       -dpz0dpy*phix       *dcosphi+ay*dphidpy)
        trans(3,6)=alb*(dpldp *phiy*xsinphi
     $       -dpz0dp *phix       *dcosphi+ay*dphidp )
        cx=-ptx*sinphi+pbx*cosphi
        cy=-pty*sinphi+pby*cosphi
        trans(2,2)=cosphi-dpldpx*phix*dcosphi
     $       -dpz0dpx*phiy*sinphi+cx*dphidpx
        trans(2,4)=      -dpldpy*phix*dcosphi
     $       +(phiz-dpz0dpy*phiy)*sinphi+cx*dphidpy
        trans(2,6)=      -dpldp *phix*dcosphi
     $       -dpz0dp *phiy*sinphi+cx*dphidp
        trans(4,2)=      -dpldpx*phiy*dcosphi
     $       +(dpz0dpx*phix-phiz)*sinphi+cy*dphidpx
        trans(4,4)=cosphi-dpldpy*phiy*dcosphi
     $       +dpz0dpy*phix*sinphi+cy*dphidpy
        trans(4,6)=      -dpldp *phiy*dcosphi
     $       +dpz0dp *phix*sinphi+cy*dphidp
        trans(5,2)= -pr*alb*dphidpx
        trans(5,4)= -pr*alb*dphidpy
        if(drift)then
          call tgetdv(cod(6),dv,dvdp)
          cod(5)=cod(5)+((dpl*phiz-dphizsq*pr)*xsinphi
     $         +dpz0*sinphi-pbz*dcosphi)*alb-dv*al
          trans(5,6)= -alb*(phi+pr*dphidp)+al*dvdp
        else
          trans(1,5)= alb*ax*dphidz
          trans(3,5)= alb*ay*dphidz
          trans(2,5)= cx*dphidz
          trans(4,5)= cy*dphidz
          cod(5)=cod(5)-pr*phi*alb
          trans(5,6)= -alb*(phi+pr*dphidp)
          trans(5,5)= -pr*alb*dphidz
        endif
      endif
      return
      end subroutine

      end module element_drift_common

      subroutine tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,ak0x,ak0y,enarad)
      use element_drift_common
      use kradlib
      use tspin, only:cphi0,sphi0
      use ffs_flag, only:rfluct,photons,calpol
      use photontable
      use mathfun
      use tmacro, only:l_track
      implicit none
      integer*4 np,i,itmax,ndiag
      real*8 conv
      parameter (itmax=15,conv=1.d-15)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np)
      real*8 al,bz,pr,bzp,phi,px1,py1,px0,py0,z0,
     $     sinphi,cosphi,ak0x,ak0y,b,phix,phiy,phiz,
     $     dphizsq,dpz0,pz0,plx,ply,plz,ptx,pty,ptz,
     $     pbx,pby,pbz,dphi,dcosphi,pl,dpl,alb,
     $     xsinphi,bpr,bsi0,ap,db
      logical*4 enarad
      data ndiag/15/
      if(ak0x .eq. 0.d0 .and. ak0y .eq. 0.d0)then
        if(bz .eq. 0.d0)then
          call tdrift_free(np,x,px,y,py,z,dv,al)
          return
        else
          call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al,bz,enarad)
          return
        endif
      else
        if(enarad)then
          cphi0=1.d0
          sphi0=0.d0
          if(calpol)then
            bsi=0.d0
          endif
          if(rfluct .and. photons)then
            call tsetpcvt(l_track,0.d0,0.d0,0.d0,0.d0,0.d0,al)
          endif
        endif
c        b=hypot(hypot(ak0x,ak0y),bz*al)
        b=norm2([ak0x,ak0y,bz*al])
        phix=ak0y/b
        phiy=ak0x/b
        phiz=bz*al/b
        dphizsq=phix**2+phiy**2
        do i=1,np
          bsi0=bz+ak0x*y(i)+ak0y*x(i)
          pr=(1.d0+g(i))
          alb=al*pr/b
          bzp=bz/pr
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          px0=px(i)
          py0=py(i)
          z0=z(i)
c          s=min(ampmax,px(i)**2+py(i)**2)
c          dpz0=-s/(1.d0+sqrtl(1.d0-s))
          dpz0=pxy2dpz(px(i),py(i))
          pz0=1.d0+dpz0
          dpl=px(i)*phix+py(i)*phiy+dpz0*phiz
          pl=phiz+dpl
          plx=pl*phix
          ply=pl*phiy
          plz=pl*phiz
          ptx=px(i)-plx
          pty=py(i)-ply
          ptz=dpz0 -dpl*phiz+dphizsq
          pbx=pty*phiz-ptz*phiy
          pby=ptz*phix-ptx*phiz
          pbz=ptx*phiy-pty*phix
          bpr=b/pr
          db=bpr-pbz
          ap=hypot(pz0,pbz)
          dphi=-atan(pbz,pz0)
          phi=asin(min(1.d0,max(-1.d0,db/ap)))-dphi
          call tsolconv(pz0,plz,pbz,bpr,
     $     phi,sinphi,xsinphi,cosphi,dcosphi,ndiag)
          x(i)=x(i)+(plx*phi+ptx*sinphi-pbx*dcosphi)*alb
          y(i)=y(i)+(ply*phi+pty*sinphi-pby*dcosphi)*alb
          z(i)=z(i)+((dpl*phiz-dphizsq)*xsinphi
     $         +dpz0*sinphi-pbz*dcosphi)*alb-dv(i)*al
          px1=px0+ptx*dcosphi+pbx*sinphi
          py1=py0+pty*dcosphi+pby*sinphi
          if(enarad)then
            bsi(i)=bsi(i)+bsi0-ak0x*y(i)-ak0y*x(i)
            if(rfluct)then
              call tradkf1(x(i),px1,y(i),py1,z(i),g(i),dv(i),
     $         sx(i),sy(i),sz(i),
     $         px0,py0,z0,bsi(i),al,i)
            else
              call tradk1(x(i),px1,y(i),py1,z(i),g(i),dv(i),
     $         sx(i),sy(i),sz(i),
     $         px0,py0,z0,bsi(i),al)
            endif
          endif
          px(i)=px1-bzp*y(i)*.5d0
          py(i)=py1+bzp*x(i)*.5d0
        enddo
      endif
      return
      end
