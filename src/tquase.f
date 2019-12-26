      subroutine tquase(trans,cod,beam,srot,al,ak,bz,
     1     dx,dy,theta,radlvl,
     1     fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0)
      use tfstk
      use ffs_flag
      use tmacro
      use temw, only:tsetr0
      use tspin, only:tradke      
      implicit none
      real*8 , intent(in)::ak,al,bz,dx,dy,theta,
     $     f1in,f2in,f1out,f2out,eps0
      real*8 trans(6,12),cod(6),beam(42),srot(3,9),
     $    radlvl,aln,akn
      integer*4 mfring
      logical*4 fringe
      real*8 bxs,bys,bzs,theta1,ak1
      integer*4 itgetqraddiv
      integer*4 i,ndiv
      integer*4 , parameter :: ndivmax=1000
      logical*4 enarad,krad
      if(al .eq. 0.d0)then
        call tthine(trans,cod,beam,srot,4,
     $       al,ak,dx,dy,theta,.false.)
        return
      endif
      enarad=radlvl .ne. 1.d0
      krad=enarad .and. al .ne. 0.d0
      cod(2)=cod(2)+.5d0*bz*dy
      cod(4)=cod(4)-.5d0*bz*dx
      if(ak .lt. 0.d0)then
        theta1=theta-m_pi_2
        ak1=-ak
      else
        theta1=theta
        ak1=ak
      endif
      call tsolrot(trans,cod,beam,al,0.d0,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.true.)
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring .ne. 2)then
        call tqfrie(trans,cod,beam,ak1,al,bz)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,f1in,f2in,bz)
      endif
      if(krad .and. f1in .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1in,0.d0,bz*.5d0)
      else
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
      endif
c      if(ifv .eq. 0)then
        if(krad)then
          if(eps0 .eq. 0.d0)then
            ndiv=max(1,itgetqraddiv(cod,ak1,al,bzs*.5d0))
          else
            ndiv=max(1,
     $           int(dble(itgetqraddiv(cod,ak1,al,bzs*.5d0))/eps0))
          endif
          ndiv=min(ndivmax,ndiv)
        else
          ndiv=1
        endif
c        write(*,*)'tquase ',ndiv,krad,irad,radlvl
        aln=al/ndiv
        akn=ak1/ndiv
        do i=1,ndiv
          call tsolque(trans,cod,beam,srot,aln,akn,
     $         bz,0.d0,0.d0,eps0,krad,irad)
        enddo
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak1,-f1out,f2out,bz)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring .ne. 1)then
        call tqfrie(trans,cod,beam,-ak1,al,bz)
      endif
      if(krad .and. f1out .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1out,0.d0,bzs*.5d0)
      endif
      call tsolrot(trans,cod,beam,al,0.d0,dx,dy,0.d0,
     $     0.d0,0.d0,theta1,bxs,bys,bzs,.false.)
      cod(2)=cod(2)-.5d0*bz*dy
      cod(4)=cod(4)+.5d0*bz*dx
      return
      end

      subroutine tsolrot(trans,cod,beam,al,bz,dx,dy,dz,
     $     chi1,chi2,theta,bxs,bys,bzs,ent)
      use mathfun, only: sqrtl
      implicit none
      integer*4 i,itgetirad
      real*8 trans(6,12),cod(6),beam(21),trans1(6,6),
     $     trans2(6,6)
      real*8 bz,dx,dy,theta,cost,sint,x0,px0,bzh,dz,chi1,chi2,
     $     al,s0,cchi1,schi1,cchi2,schi2,dcchi1,dcchi2,
     $     ds1,pr,pz0,dpz0dpx,dpz0dpy,dpz0dp,
     $     pz1,ds2,bxs,bys,bzs,bxs0,bzh1,pzmin,a,ptmax
      logical*4 ent
      parameter (pzmin=1.d-10,ptmax=0.9999d0)
      call tinitr(trans1)
      if(ent)then
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or.
     $       chi2 .ne. 0.d0)then
          bzh1=bz*.5d0
          cod(2)=cod(2)+bzh1*cod(3)
          cod(4)=cod(4)-bzh1*cod(1)
          s0=-.5d0*al
          cchi1=cos(chi1)
          schi1=sin(chi1)
          if(cchi1 .ge. 0.d0)then
            dcchi1=schi1**2/(1.d0+cchi1)
          else
            dcchi1=(1.d0-cchi1)
          endif
          cchi2=cos(chi2)
          schi2=sin(chi2)
          if(cchi2 .ge. 0.d0)then
            dcchi2=schi2**2/(1.d0+cchi2)
          else
            dcchi2=(1.d0-cchi2)
          endif
          bzs=cchi2*cchi1*bz
          bxs=-schi1*bz
          bys=-schi2*cchi1*bz
          cod(1)=cod(1)-dx
          cod(3)=cod(3)-dy
          pr=1.d0+cod(6)
          a=cod(2)**2+cod(4)**2
          pz0=pr*sqrtl(1.d0-a/pr**2)
c          pz0=sqrt(max(pzmin,(pr-cod(2))*(pr+cod(2))-cod(4)**2))
          dpz0dpx=-cod(2)/pz0
          dpz0dpy=-cod(4)/pz0
          dpz0dp = pr/pz0
          pz1   = schi1*cod(2)+cchi1*pz0
          cod(2)= cchi1*cod(2)-schi1*pz0
          cod(4)= cchi2*cod(4)-schi2*pz1
          trans1(2,2)= cchi1      -schi1*dpz0dpx
          trans1(2,4)=            -schi1*dpz0dpy
          trans1(2,6)=            -schi1*dpz0dp
          trans1(4,2)=-schi2*schi1-schi2*cchi1*dpz0dpx
          trans1(4,4)= cchi2      -schi2*cchi1*dpz0dpy
          trans1(4,6)=            -schi2*cchi1*dpz0dp
          ds1   = schi1*cod(1)-dcchi1*s0-cchi1*dz
          cod(1)= cchi1*cod(1)-schi1*(s0-dz)
          ds2   = schi2*cod(3)+cchi2*ds1-dcchi2*s0
          cod(3)= cchi2*cod(3)-schi2*(ds1+s0)
          trans1(1,1)= cchi1
          trans1(3,1)=-schi2*schi1
          trans1(3,3)= cchi2
          trans1(5,1)=-cchi2*schi1
          trans1(5,3)=-schi2
          trans1(5,5)=0.d0
          call tsoldz(trans2,cod,-ds2,bxs,bys,bzs,.false.)
          call tmultr(trans1,trans2,6)
          trans1(5,5)=1.d0
          bzh=bzs*.5d0
          cod(2)=cod(2)-bzh*cod(3)
          cod(4)=cod(4)+bzh*cod(1)
          do i=1,6
            trans1(i,1)=trans1(i,1)-bzh1*trans1(i,4)
            trans1(i,3)=trans1(i,3)+bzh1*trans1(i,2)
          enddo
          do i=1,6
            trans1(2,i)=trans1(2,i)-bzh*trans1(3,i)
            trans1(4,i)=trans1(4,i)+bzh*trans1(1,i)
          enddo
        else
          bzh=bz*.5d0
          cod(1)=cod(1)-dx
          cod(3)=cod(3)-dy
          cod(2)=cod(2)+bzh*dy
          cod(4)=cod(4)-bzh*dx
          bxs=0.d0
          bys=0.d0
          bzs=bz
        endif
      endif
      if(theta .ne. 0.d0)then
        cost=cos(theta)
        if(ent)then
          sint=sin(theta)
        else
          sint=-sin(theta)
        endif
        x0=cod(1)
        cod(1)=cost*x0-sint*cod(3)
        cod(3)=sint*x0+cost*cod(3)
        px0=cod(2)
        cod(2)=cost*px0-sint*cod(4)
        cod(4)=sint*px0+cost*cod(4)
        do i=1,6
          x0=trans1(1,i)
          trans1(1,i)=cost*x0-sint*trans1(3,i)
          trans1(3,i)=sint*x0+cost*trans1(3,i)
          px0=trans1(2,i)
          trans1(2,i)=cost*px0-sint*trans1(4,i)
          trans1(4,i)=sint*px0+cost*trans1(4,i)
        enddo
        bxs0=bxs
        bxs=cost*bxs0-sint*bys
        bys=sint*bxs0+cost*bys
      endif
      if(.not. ent)then
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or.
     $       chi2 .ne. 0.d0)then
          s0=-.5d0*al
          cchi1=cos(chi1)
          schi1=sin(chi1)
          if(cchi1 .ge. 0.d0)then
            dcchi1=schi1**2/(1.d0+cchi1)
          else
            dcchi1=(1.d0-cchi1)
          endif
          cchi2=cos(chi2)
          schi2=sin(chi2)
          if(cchi2 .ge. 0.d0)then
            dcchi2=schi2**2/(1.d0+cchi2)
          else
            dcchi2=(1.d0-cchi2)
          endif
          call tinitr(trans2)
          pr=1.d0+cod(6)
          bzh=bzs*.5d0
          cod(2)=cod(2)+bzh*cod(3)
          cod(4)=cod(4)-bzh*cod(1)
          a=cod(2)**2+cod(4)**2
          pz0=pr*sqrtl(1.d0-a/pr**2)
c          pz0=sqrt(max(pzmin,(pr-cod(2))*(pr+cod(2))-cod(4)**2))
          dpz0dpx=-cod(2)/pz0
          dpz0dpy=-cod(4)/pz0
          dpz0dp = pr/pz0
          ds1   =-schi2*cod(3)+dcchi2*s0
          cod(3)= cchi2*cod(3)-schi2*s0
          ds2   =-schi1*cod(1)+cchi1*ds1+dcchi1*s0+dz
          cod(1)= cchi1*cod(1)+schi1*(ds1-s0)
          pz1   =-schi2*cod(4)+cchi2*pz0
          cod(2)= cchi1*cod(2)+schi1*pz1
          cod(4)= cchi2*cod(4)+schi2*pz0
          trans2(1,1)= cchi1
          trans2(1,3)=-schi1*schi2
          trans2(3,3)= cchi2
          trans2(5,1)= schi1
          trans2(5,3)= cchi1*schi2
          trans2(2,2)= cchi1         +schi1*cchi2*dpz0dpx
          trans2(2,4)=-schi1*schi2   +schi1*cchi2*dpz0dpy
          trans2(2,3)= trans2(2,2)*bzh
          trans2(2,1)=-trans2(2,4)*bzh
          trans2(2,6)=                schi1*cchi2*dpz0dp
          trans2(4,2)=                schi2*dpz0dpx
          trans2(4,4)= cchi2         +schi2*dpz0dpy
          trans2(4,6)=                schi2*dpz0dp
          trans2(4,3)= trans2(4,2)*bzh
          trans2(4,1)=-trans2(4,4)*bzh
          trans1(5,5)=0.d0
          call tmultr5(trans1,trans2,6)
          call tsoldz(trans2,cod,-ds2,0.d0,0.d0,bz,.false.)
          call tmultr(trans1,trans2,6)
          trans1(5,5)=1.d0
          cod(1)=cod(1)+dx
          cod(3)=cod(3)+dy
          bzh1=bz*.5d0
          cod(2)=cod(2)-bzh1*cod(3)
          cod(4)=cod(4)+bzh1*cod(1)
          do i=1,6
            trans1(2,i)=trans1(2,i)-bzh1*trans1(3,i)
            trans1(4,i)=trans1(4,i)+bzh1*trans1(1,i)
          enddo
        else
          bzh=bz*.5d0
          cod(2)=cod(2)-bzh*dy
          cod(4)=cod(4)+bzh*dx
          cod(1)=cod(1)+dx
          cod(3)=cod(3)+dy
        endif
      endif
c      write(*,'(a/,6(1p,6g11.4/))')
c     $     'tsolrot ',((trans1(i,j),j=1,6),i=1,6)
      call tmultr5(trans,trans1,itgetirad())
      call tmulbs(beam,trans1,.false.,.true.)
      return
      end

      subroutine tsoldz(trans,cod,al,bxs0,bys0,bzs0,drift)
      use mathfun
      implicit none
      integer*4 j,itmax,ndiag
      parameter (itmax=15)
      real*8 trans(6,6),cod(6),al,bxs,bys,bzs,pxi,pyi,pz0,
     $     dpz0dpx,dpz0dpy,dpz0dp,phi,
     $     dphidz,dphidpx,dphidpy,dphidp,a24,a12,a22,a14,
     $     da12,da14,pr,dpz0,dv,dvdp,s,r,phix,phiy,phiz,babs,
     $     alb,pbx,pby,pbz,pl,dpl,dphizsq,a,
     $     dpldpx,dpldpy,dpldp,dplz,plx,ply,plz,ptx,pty,ptz,
     $     cosphi,sinphi,dcosphi,phi0,dphi,
     $     xsinphi,ax,ay,az,cx,cy,conv,albabs,u,
     $     bxs0,bys0,bzs0,bzthre,ptmax
      parameter (conv=1.d-15,bzthre=1.d-20,ptmax=0.9999d0)
      logical*4 drift
      data ndiag/15/
      bxs=bxs0
      bys=bys0
      bzs=bzs0
      babs=sqrt(bzs**2+bxs**2+bys**2)
      if(abs(babs) .lt. bzthre)then
        bxs=0.d0
        bys=0.d0
        bzs=0.d0
        babs=0.d0
      endif
      call tinitr(trans)
      pr=1.d0+cod(6)
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
          if(a22 .ge. 0.d0)then
            a14=a24**2/(1.d0+a22)/bzs
          else
            a14=(1.d0-a22)/bzs
          endif
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
        if(al .ne. 0.d0)then
          phi=asin(min(1.d0,max(-1.d0,albabs/pz0)))
          dphi=0.d0
          do j=1,itmax
            sinphi=sin(phi)
            dcosphi=2.d0*sin(.5d0*phi)**2
            xsinphi=xsin(phi)
            s=plz*xsinphi+pz0*sinphi+pbz*dcosphi
            u=pz0-ptz*dcosphi+pbz*sinphi
            if(u .ne. 0.d0)then
              dphi=(albabs-s)/u
            else
              dphi=0.d0
            endif
            phi0=phi
            phi=phi+dphi
            if(phi0 .eq. phi .or. abs(dphi) .le. conv*abs(phi))then
              go to 100
            endif
          enddo
          if(ndiag .ge. 0)then
            ndiag=ndiag-1
            write(*,*)'tsoldz convergence error',phi,dphi,u,babs
            if(ndiag .eq. -1)then
              write(*,*)
     $             'Further tsoldz message will be suppressed.'
            endif
          endif
        else
          phi=0.d0
          xsinphi=0.d0
          sinphi=0.d0
          dcosphi=0.d0
        endif
 100    dplz=-pr*dphizsq+dpl*phiz
        cosphi=cos(phi)
        ax=pxi-ptx*dcosphi+pbx*sinphi
        ay=pyi-pty*dcosphi+pby*sinphi
        az=pz0-ptz*dcosphi+pbz*sinphi
        dphidpx=-(dpldpx*phiz*xsinphi+dpz0dpx*sinphi
     $       +phiy*dcosphi)/az
        dphidpy=-(dpldpy*phiz*xsinphi+dpz0dpy*sinphi
     $       -phix*dcosphi)/az
        dphidp =-(dpldp *phiz*xsinphi+dpz0dp*sinphi)/az
        dphidz =babs/az
        cod(1)=cod(1)+(plx*phi+ptx*sinphi+pbx*dcosphi)*alb
        cod(3)=cod(3)+(ply*phi+pty*sinphi+pby*dcosphi)*alb
        cod(2)=pxi-ptx*dcosphi+pbx*sinphi
        cod(4)=pyi-pty*dcosphi+pby*sinphi
        trans(1,2)=alb*(dpldpx*phix*xsinphi+sinphi
     $       -dpz0dpx*phiy       *dcosphi+ax*dphidpx)
        trans(1,4)=alb*(dpldpy*phix*xsinphi
     $       +(phiz-dpz0dpy*phiy)*dcosphi+ax*dphidpy)
        trans(1,6)=alb*(dpldp *phix*xsinphi
     $       -dpz0dp *phiy       *dcosphi+ax*dphidp )
        trans(3,2)=alb*(dpldpx*phiy*xsinphi
     $       +(dpz0dpx*phix-phiz)*dcosphi+ay*dphidpx)
        trans(3,4)=alb*(dpldpy*phiy*xsinphi+sinphi
     $       +dpz0dpy*phix       *dcosphi+ay*dphidpy)
        trans(3,6)=alb*(dpldp *phiy*xsinphi
     $       +dpz0dp *phix       *dcosphi+ay*dphidp )
        cx=-ptx*sinphi+pbx*cosphi
        cy=-pty*sinphi+pby*cosphi
        trans(2,2)=cosphi+dpldpx*phix*dcosphi
     $       -dpz0dpx*phiy*sinphi+cx*dphidpx
        trans(2,4)=       dpldpy*phix*dcosphi
     $       +(phiz-dpz0dpy*phiy)*sinphi+cx*dphidpy
        trans(2,6)=       dpldp *phix*dcosphi
     $       -dpz0dp *phiy*sinphi+cx*dphidp
        trans(4,2)=       dpldpx*phiy*dcosphi
     $       +(dpz0dpx*phix-phiz)*sinphi+cy*dphidpx
        trans(4,4)=cosphi+dpldpy*phiy*dcosphi
     $       +dpz0dpy*phix*sinphi+cy*dphidpy
        trans(4,6)=       dpldp *phiy*dcosphi
     $       +dpz0dp *phix*sinphi+cy*dphidp
        trans(5,2)= -pr*alb*dphidpx
        trans(5,4)= -pr*alb*dphidpy
        if(drift)then
          call tgetdv(cod(6),dv,dvdp)
          cod(5)=cod(5)+((dpl*phiz-dphizsq*pr)*xsinphi
     $         +dpz0*sinphi+pbz*dcosphi)*alb-dv*al
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
      end
