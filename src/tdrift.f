c     common module for drift element implementation
      module element_drift_common
      implicit none
      private

c     Maximum amplitude of (px/p0)^2 + (py/p0)^2
      real(8), public, parameter :: ampmax = 0.9999d0

      end module element_drift_common

c     drift in the free space
      subroutine tdrift_free(np,x,px,y,py,z,dv,al)
      use element_drift_common
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np)
      real*8 al,al1,dpz
      do i=1,np
        dpz=pxy2dpz(px(i),py(i))
        al1=al/(1.d0+dpz)
        x(i)=x(i)+px(i)*al1
        y(i)=y(i)+py(i)*al1
        z(i)=z(i)+dpz  *al1-dv(i)*al
      enddo
      return
      end

c     drift in the parallel solenoid
      subroutine tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,enarad)
      use element_drift_common
      use ffs_flag, only:rfluct,photons,calpol
      use photontable
      use kradlib
      use tspin, only:cphi0,sphi0
      use mathfun, only: sqrtl,pxy2dpz,xsincos
      use tmacro, only:l_track
      implicit none
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np),px1,py1
      real*8 al,bz
      integer*4 i
      real*8 pr,bzp,pxi,pyi
      real*8 dpzi,pzi,al1,zi
      real*8 phi,a22,a24,a12,a14
      logical*4 enarad
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

      subroutine tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,ak0x,ak0y,enarad)
      use element_drift_common
      use kradlib
      use tspin, only:cphi0,sphi0
      use sol, only:tsolconv
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
        b=hypot3(ak0x,ak0y,bz*al)
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
