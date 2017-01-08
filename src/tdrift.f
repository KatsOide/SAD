c     common module for drift element implementation
      module element_drift_common
      implicit none
      private

c     Maximum amplitude of (px/p0)^2 + (py/p0)^2
      real(8), public, parameter :: ampmax = 0.9999d0

      end module element_drift_common

c     drift in the free space
      subroutine tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
      use element_drift_common
      use tfstk
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),pz(np)
      real*8 al,al1,dpz
      do i=1,np
c        s=px(i)**2+py(i)**2
        dpz=pxy2dpz(px(i),py(i))
c        dpz=s*(-.5d0-s*(.125d0+s*.0625d0))
c        dpz=(dpz**2-s)/(2.d0+2.d0*dpz)
c        dpz=(dpz**2-s)/(2.d0+2.d0*dpz)
        al1=al/(1.d0+dpz)
        x(i)=x(i)+px(i)*al1
        y(i)=y(i)+py(i)*al1
        z(i)=z(i)+dpz  *al1-dv(i)*al
      enddo
      return
      end

c     drift in the parallel solenoid
      subroutine tdrift_solenoid(np,x,px,y,py,z,g,dv,pz,al,bz)
      use tfstk, only: sqrtl
      use element_drift_common
      implicit none
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),pz(np)
      real*8 al,bz
      integer*4 i
      real*8 pr,bzp,pxi,pyi
      real*8 s,dpzi,pzi,al1
      real*8 phi,a22,a24,a12,a14
      do i=1,np
c         pr=(1.d0+g(i))**2
         pr=(1.d0+g(i))
         bzp=bz/pr
         pxi=px(i)+bzp*y(i)*.5d0
         pyi=py(i)-bzp*x(i)*.5d0

         s=min(ampmax,pxi**2+pyi**2)
         dpzi=-s/(1.d0+sqrtl(1.d0-s))
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
         a22=cos(phi)
         a24=sin(phi)
         if(abs(phi) .gt. 1.d-4)then
            a12=a24/bzp
         else
            a12=(1.d0-(1.d0/6.d0)*phi*phi)*al1
         endif
         if(a22 .ge. 0.d0)then
            a14=a12*a24/(1.d0+a22)
         else
            a14=(1.d0-a22)/bzp
         endif

         x(i) =x(i)+a12*pxi+a14*pyi
         y(i) =y(i)-a14*pxi+a12*pyi
         px(i)=     a22*pxi+a24*pyi-bzp*y(i)*.5d0
         py(i)=    -a24*pxi+a22*pyi+bzp*x(i)*.5d0
         z(i) =z(i)+dpzi *al1-dv(i)*al
      enddo
      return
      end

      subroutine tdrift(np,x,px,y,py,z,g,dv,pz,
     $     al,bz,ak0x,ak0y)
      use element_drift_common
      use tfstk, only:sqrtl
      implicit none
      integer*4 np,i,j,itmax,ndiag
      real*8 conv
      parameter (itmax=15,conv=1.d-15)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),pz(np),g(np)
      real*8 al,bz,pr,bzp,s,phi,
     $     sinphi,ak0x,ak0y,b,phix,phiy,phiz,
     $     dphizsq,dpz0,pz0,plx,ply,plz,ptx,pty,ptz,
     $     pbx,pby,pbz,phi0,dphi,dcosphi,pl,dpl,alb,
     $     xsinphi,xsin,r,bpr
      data ndiag/15/
      if(ak0x .eq. 0.d0 .and. ak0y .eq. 0.d0)then
        if(bz .eq. 0.d0)then
          call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
          return
        else
          call tdrift_solenoid(np,x,px,y,py,z,g,dv,pz,al,bz)
          return
        endif
      else
c        b=hypot(hypot(ak0x,ak0y),bz*al)
        b=abs(dcmplx(abs(dcmplx(ak0x,ak0y)),bz*al))
        phix=ak0y/b
        phiy=ak0x/b
        phiz=bz*al/b
        dphizsq=phix**2+phiy**2
        do i=1,np
c          pr=(1.d0+g(i))**2
          pr=(1.d0+g(i))
          alb=al*pr/b
          bzp=bz/pr
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          s=min(ampmax,px(i)**2+py(i)**2)
          dpz0=-s/(1.d0+sqrtl(1.d0-s))
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
          phi=asin(min(1.d0,max(-1.d0,bpr/pz0)))
          dphi=0.d0
          do j=1,itmax
            sinphi=sin(phi)
            dcosphi=2.d0*sin(.5d0*phi)**2
            xsinphi=xsin(phi)
            s=plz*xsinphi+pz0*sinphi+pbz*dcosphi
            r=pz0-ptz*dcosphi+pbz*sinphi
            if(r .ne. 0.d0)then
              dphi=(bpr-s)/r
            endif
            phi0=phi
            phi=phi+dphi
            if(phi0 .eq. phi .or. abs(dphi) .le. conv*abs(phi))then
              go to 100
            endif
          enddo
          if(ndiag .ge. 0)then
            ndiag=ndiag-1
            write(*,*)'tdrift convergence error',phi,dphi
            if(ndiag .eq. -1)then
              write(*,*)
     $             'Further tdrift messages will be suppressed.'
            endif
          endif
 100      x(i)=x(i)+(plx*phi+ptx*sinphi+pbx*dcosphi)*alb
          y(i)=y(i)+(ply*phi+pty*sinphi+pby*dcosphi)*alb
          z(i)=z(i)+((dpl*phiz-dphizsq)*xsinphi
     $         +dpz0*sinphi+pbz*dcosphi)*alb-dv(i)*al
          px(i)=px(i)-ptx*dcosphi+pbx*sinphi-bzp*y(i)*.5d0
          py(i)=py(i)-pty*dcosphi+pby*sinphi+bzp*x(i)*.5d0
        enddo
      endif
      return
      end
