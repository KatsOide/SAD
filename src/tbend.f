      recursive subroutine tbend(np,x,px,y,py,z,g,dv,pz,l,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1     fb1,fb2,mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     enarad,alb,ale,ala,eps)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      use multa, only:nmult
      implicit none
      integer*4 np,mfring,i,l,ndiv,mfr1,mfr2,ndivmax
      parameter (ndivmax=1024)
      real*8 al,phib,phi0,cosp1,sinp1,cosp2,sinp2,ak,dx,dy,theta,
     $     dphix,dphiy,cost,sint,cosw,sinw,sqwh,sinwp1,eps,
     $     tanp1,tanp2,rhob,rho0,b,drhob,dp,p,
     $     pinv,rhoe,pxi,pyi,s,dpzi,pzi,sp1,x1,dz1,y1,z1,px1,
     $     py1,pv1sqi,f,ff,x2,py2,z2,dph2,ph2,dpx2,pz2,drho,
     $     t2,dpx3,px3,dpz3,pz3,t3,x3,da,y3,z3,pv2sqi,x4,py4,z4,dpz4,
     $     dz4,pr,xi,dxfr1,dyfr1,dzfr1,dxfr2,dyfr2,dzfr2,dpz32,
     $     dyfra1,dyfra2,fa,t4,dpx3a,dpz2,t2t3,dcosp,px1px3,
     $     psi1,psi2,wn1,wn2,wnc,aln,phibn,phi0n,alb,ale,ala,als,
     $     fb1,fb2,eps1,f1r,f2r,
     $     coswn1,sinwn1,sqwhn1,sinwp1n1,
     $     coswnc,sinwnc,sqwhnc,sinwp1nc,
     $     coswn2,sinwn2,sqwhn2,sinwp1n2
      real*8 a3,a5,a7,a9,a11,a13,a15
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      real*8 smax,smin,rphidiv
      parameter (smax=0.99d0,smin=0.01d0,rphidiv=3e-3)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np)
      complex*16 akm(0:nmult)
      logical*4 enarad,fringe
      if(phi0 .eq. 0.d0)then
        if(ak .eq. 0.d0)then
          call tsteer(np,x,px,y,py,z,g,dv,pz,l,al,-phib,
     1         dx,dy,theta,cost,sint,
     1         cosp1,sinp1,cosp2,sinp2,
     $         fb1,fb2,mfring,fringe,eps)
        elseif(phib .eq. phi0)then
          call tquad(np,x,px,y,py,z,g,dv,pz,l,al,ak,
     1         dx,dy,theta,cost,sint,0.d0,.true.,
     1         fringe,0.d0,0.d0,0,eps,.true.)
        else
          akm=0.d0
          akm(0)=phib-phi0
          akm(1)=ak
          call tmulti(np,x,px,y,py,z,g,dv,pz,
     $         al,ak,0.d0,0.d0,
     $         psi1,psi2,
     $         dx,dy,0.d0,0.d0,0.d0,theta,0.d0,
     $         eps,enarad,fringe,
     $         0.d0,0.d0,0.d0,0.d0,
     $         mfring,fb1,fb2,
     $         0.d0,0.d0,0.d0,0.d0,0.d0,
     $         .false.,.false.,
     $         int8(0),int8(0),int8(0),int8(0))
        endif
        return
      elseif(phib .eq. 0.d0)then
        call tbdrift(np,x,px,y,py,z,dv,pz,al,phi0)
        return
      elseif(al .eq. 0.d0)then
        call tbthin(np,x,px,y,py,z,g,pz,phib,phi0,dx,dy,
     1              theta,cost,sint)
        return
      elseif(rad .and. enarad .and. trpt)then
        call tbrad(np,x,px,y,py,z,g,dv,pz,l,al,phib,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1       fb1,fb2,mfring,
     1       fringe,eps)
        return
      elseif(ak .ne. 0.d0)then
        call tbendi(np,x,px,y,py,z,g,dv,pz,l,al,phib,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1       fb1,fb2,mfring,enarad,fringe,eps)
        return
      endif
      include 'inc/TENT.inc'
      if(dphiy .ne. 0.d0)then
        do i=1,np
          pr=1.d0+g(i)
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
        enddo
      endif
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
        eps1=eps
        if(eps1 .eq. 0.d0)then
          eps1=1.d0
        endif
        ndiv=min(ndivmax,int(1+(abs(phi0)*
     $       max(h0/100.d0,1.d0/rphidiv)/eps1)))
        if(ndiv .gt. 1)then
          aln=al/ndiv
          phibn=phib/ndiv
          phi0n=phi0/ndiv
          psi1=atan2(sinp1,cosp1)
          psi2=atan2(sinp2,cosp2)
          wn1=phi0n-psi1
          coswn1=cos(wn1)
          sinwn1=sin(wn1)
          if(coswn1 .gt. 0.d0)then
            sqwhn1=sinwn1**2/(1.d0+coswn1)
          else
            sqwhn1=1.d0-coswn1
          endif
          sinwp1n1=sin(phi0n)
          wn2=phi0n-psi2
          coswn2=cos(wn2)
          sinwn2=sin(wn2)
          if(coswn2 .gt. 0.d0)then
            sqwhn2=sinwn2**2/(1.d0+coswn2)
          else
            sqwhn2=1.d0-coswn2
          endif
          sinwp1n2=sinwn2
          wnc=phi0n
          coswnc=cos(wnc)
          sinwnc=sin(wnc)
          if(coswnc .gt. 0.d0)then
            sqwhnc=sinwnc**2/(1.d0+coswnc)
          else
            sqwhnc=1.d0-coswnc
          endif
          sinwp1nc=sinwnc
          if(mfring .gt. 0 .or. mfring .eq. -1)then
            mfr1=-1
          else
            mfr1=0
          endif
          if(mfring .gt. 0 .or. mfring .eq. -2)then
            mfr2=-2
          else
            mfr2=0
          endif
          als=alb+aln
          call tbend(np,x,px,y,py,z,g,dv,pz,l,aln,phibn,phi0n,
     1         cosp1,sinp1,1.d0,0.d0,
     1         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     1         fb1,fb2,mfr1,fringe,coswn1,sinwn1,sqwhn1,sinwp1n1,
     1         .true.,alb,als,ala,eps)
          do i=2,ndiv-1
            call tbend(np,x,px,y,py,z,g,dv,pz,l,aln,phibn,phi0n,
     1           1.d0,0.d0,1.d0,0.d0,
     1           0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     1           fb1,fb2,0,.false.,coswnc,sinwnc,sqwhnc,sinwp1nc,
     1           .true.,als,als+aln,ala,eps)
            als=als+aln
          enddo
          call tbend(np,x,px,y,py,z,g,dv,pz,l,aln,phibn,phi0n,
     1         1.d0,0.d0,cosp2,sinp2,
     1         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     1         fb1,fb2,mfr2,fringe,coswn2,sinwn2,sqwhn2,sinwp1n2,
     1         .true.,als,ale,ala,eps)
          if(dphiy .ne. 0.d0)then
            do i=1,np
              pr=1.d0+g(i)
              px(i)=px(i)+dphix/pr
              py(i)=py(i)+dphiy/pr
            enddo
          endif
          include 'inc/TEXIT.inc'
          return
        endif
      endif
      rhob=al/phib
      rho0=al/phi0
      if(rad .and. enarad)then
        tanp1=sinp1/cosp1
        b=brhoz/rhob
        if(alb .ne. 0.d0)then
          als=alb+al*.25d0
        else
          als=0.d0
        endif
        call trad(np,x,px,y,py,g,dv,b,0.d0,0.d0,
     1       1.d0/rho0,-tanp1*2.d0/al,.5d0*al,
     $       f1r,f2r,als,ala,1.d0)
      endif
      dxfr1=0.d0
      dyfr1=0.d0
      dyfra1=0.d0
      dzfr1=0.d0
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fb1**2/rhob/24.d0
          dyfr1=fb1/rhob**2/6.d0
          if(fringe)then
            dyfra1=4.d0*dyfr1/fb1**2
          endif
          dzfr1=dxfr1*sinp1
        endif
      endif
      dxfr2=0.d0
      dyfr2=0.d0
      dyfra2=0.d0
      dzfr2=0.d0
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          if(fringe)then
            dyfra2=4.d0*dyfr2/fb2**2
          endif
          dzfr2=dxfr2*sinp2
        endif
      endif
      if(cosp1*cosp2 .gt. 0.d0)then
        dcosp=(sinp2-sinp1)*(sinp2+sinp1)/(cosp1+cosp2)
      else
        dcosp=cosp1-cosp2
      endif
      drhob=rhob-rho0
      do 100 i=1,np
        dp=g(i)
        p=1.d0+dp
        pinv=1.d0/p
        rhoe=rhob*p
        pxi=px(i)
        pyi=py(i)
        s=min(smax,pxi**2+pyi**2)
c        dpzi=-s/(1.d0+sqrt(1.d0-s))
        dpzi=sqrt1(-s)
        pzi=1.d0+dpzi
        sp1=sinp1/pzi
        x1=x(i)/(cosp1-pxi*sp1)
        dz1=x1*sp1
        y1=y(i)+pyi*dz1
        z1=z(i)-dz1
        px1= pxi*cosp1+pzi*sinp1
        x1=x1+dxfr1*dp*pinv
        py1=pyi+(dyfr1-dyfra1*y1**2)*y1*pinv**2
        z1=z1+(dxfr1*px1+
     $       (.5d0*dyfr1-.25d0*dyfra1*y1**2)*y1**2*pinv)*pinv-dzfr1
        pv1sqi=1.d0/max(smin,1.d0-px1**2)
        fa=y1/rhoe*sqrt(pv1sqi)
        f=(1.d0-(y1/rhob)**2/6.d0)*fa
        ff=(.5d0-(y1/rhob)**2/24.d0)*y1*fa*pv1sqi
        x2=x1+ff
        py2=py1-px1*f
        z2=z1-px1*ff
        s=min(smax,py2**2)
c        dph2=-s/(1.d0+sqrt(1.d0-s))
        dph2=sqrt1(-s)
        ph2=1.d0+dph2
        dpx2=pxi*cosp1+(dpzi-dph2)*sinp1
        s=min(smax,px1**2+py2**2)
c        dpz2=-s/(1.d0+sqrt(1.d0-s))
        dpz2=sqrt1(-s)
        pz2=1.d0+dpz2
c        pz2=sqrt(max(smin,1.d0-py2**2-px1**2))
        drho=drhob+rhoe*dph2+rhob*dp
        t2=(px1+ph2*sinp1)/(pz2+ph2*cosp1)
        dpx3a=(x2*sinw-drho*(sinp2+sinwp1))/rhoe
        dpx3=dpx3a-dpx2*(cosw-sinw*t2)
        px3=ph2*sinp2+dpx3
        px1px3=ph2*(sinp2+sinp1)+dpx3a+dpx2*(sqwh+sinw*t2)
        s=min(smax,px3**2+py2**2)
c        dpz3=-s/(1.d0+sqrt(1.d0-s))
        dpz3=sqrt1(-s)
        pz3=1.d0+dpz3
        dpz32=px1px3*(px1-px3)/(pz2+pz3)
        t3=(px3+ph2*sinp2)/(pz3+ph2*cosp2)
c        t2t3=ph2*(sinp2/(pz3+ph2*cosp2)+sinp1/(pz2+ph2*cosp1))
c     $       +((ph2*sinp2+dpx3a-dpx2*(cosw-sinw*t2))*(pz2+ph2*cosp1)
c     $       +(dpx2+ph2*sinp1)*(pz3+ph2*cosp2))/
c     $       (pz3+ph2*cosp2)/(pz2+ph2*cosp1)
        t2t3=(ph2*sinp2+px1px3)/(pz3+ph2*cosp2)
     $       +ph2*sinp1/(pz2+ph2*cosp1)
     $       +px1*(dpz32-ph2*dcosp)
     $       /(pz3+ph2*cosp2)/(pz2+ph2*cosp1)
c        write(*,*)t2,t3,t2+t3,t2t3,px1+px3,px1px3
        t4=(cosp2+t3*sinp2)*(pz2*cosp1+px1*sinp1)
c        x3=x2*cosw+
c     1     (rho0*((cosw*t2+sinw)*dpx2-dpx3*t3)-
c     1      drho*((px3-px1)*(px3+px1)/(pz3+pz2)-sqwh*pz2-sinw*px1))/ph2
        x3=x2*(cosw-rho0/rhoe*t3*sinw/ph2)
     1       +(rho0*(cosw*t2t3+sinw*(1.d0-t2*t3))*dpx2-
     1       drho*(-(sinp2+sinwp1)*rho0/rhoe*t3
     $       -dpz32-sqwh*pz2-sinw*px1))/ph2
c        da=asin(min(1.d0,max(-1.d0,
c     $       (dpx2*((cosp1+t2*sinp1)*(pz3*cosp2+px3*sinp2)
c     $       -t4*(cosw-sinw*t2))
c     1       +dpx3a*t4)/ph2**2)))
        da=asin(min(1.d0,max(-1.d0,
     $       (dpx2*(
     $       sinp1*(t2*(pz3*cosp2+px3*sinp2)-px1*cosp2)
     $       -sinp2*(t3*(pz2*cosp1+px1*sinp1)-px3*cosp1)
     $       +cosp1*cosp2*dpz32
     $       +t4*(sqwh+sinw*t2))
     1       +dpx3a*t4)/ph2**2)))
        y3=y1+py2*rhoe*(phi0+da)
        z3=z2-phi0*(dp*rhob+drhob)-da*rhoe-dv(i)*al
        pv2sqi=1.d0/max(smin,1.d0-px3**2)
        fa=y3/rhoe*sqrt(pv2sqi)
        f=(1.d0-(y3/rhob)**2/6.d0)*fa
        ff=(.5d0-(y3/rhob)**2/24.d0)*y3*fa*pv2sqi
        x4=x3-ff
        py4=py2-px3*f
        z4=z3-px3*ff
        x4=x4-dxfr2*dp*pinv
        py4=py4+(dyfr2-dyfra2*y3**2)*y3*pinv**2
        z4=z4+(dxfr2*px3+
     $       (.5d0*dyfr2-.25d0*dyfra2*y3**2)*y3**2*pinv)*pinv-dzfr2
        s=min(smax,px3**2+py4**2)
c        dpz4=-s/(1.d0+sqrt(1.d0-s))
        dpz4=sqrt1(-s)
        px(i)=-cosp2*dpx3+sinp2*(dpz4-dpz3-dpx3*t3)
        dz4=x4*sinp2/(cosp2*(1.d0+dpz4)+sinp2*px3)
        x(i)=x4*cosp2+px(i)*dz4
        py(i)=py4
        y(i)=y3+py4*dz4
        z(i)=z4-dz4
100   continue
      if(rad .and. enarad)then
        tanp2=sinp2/cosp2
        if(ale .eq. ala)then
          als=ale
        else
          als=ale-al*.25d0
        endif
        call trad(np,x,px,y,py,g,dv,b,0.d0,0.d0,
     1       1.d0/rho0,-tanp2*2.d0/al,.5d0*al,
     $       f1r,f2r,als,ala,-1.d0)
      endif
      if(dphiy .ne. 0.d0)then
        do i=1,np
          pr=1.d0+g(i)
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
        enddo
      endif
      include 'inc/TEXIT.inc'
      return
      end

      subroutine tbdrift(np,x,px,y,py,z,dv,pz,al,phi0)
      use tfstk
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),pz(np),
     $     al,phi0,cp,sp,rho0,dx,xi,pzi,pzf,dl,xsin,dcp
      cp=cos(phi0)
      sp=sin(phi0)
      if(cp .ge. 0.d0)then
        dcp=sp**2/(1.d0+cp)
      else
        dcp=1.d0-cp
      endif
      rho0=al/phi0
      call tdrift_free(np,x,px,y,py,z,0.d0,dv,pz,rho0*sp)
      dx=rho0*dcp
      dl=rho0*xsin(phi0)
      do i=1,np
        xi=x(i)+dx
        pzi=1.d0+sqrt1(-px(i)**2-py(i)**2)
c        pzi=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
        pzf=pzi*cp-px(i)*sp
        x(i)=xi*pzi/pzf
        y(i)=y(i)+xi*sp*py(i)/pzf
        z(i)=z(i)-xi*sp/pzf+(1.d0-dv(i))*dl
        px(i)=px(i)*cp+pzi*sp
      enddo
      return
      end

      subroutine tbthin(np,x,px,y,py,z,g,pz,phib,phi0,dx,dy,
     1                 theta,cost,sint)
      implicit none
      integer*4 np,i
      real*8 phib,phi0,dx,dy,theta,cost,sint,xi,pxi
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),pz(np)
      include 'inc/TENT.inc'
      do 10 i=1,np
c        px(i)=px(i)+phi0-phib/(1.d0+g(i))**2
        px(i)=px(i)+phi0-phib/(1.d0+g(i))
        z(i)=z(i)-x(i)*phi0
10    continue
      include 'inc/TEXIT.inc'
      return
      end
