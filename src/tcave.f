      subroutine tcave(trans,cod,beam,srot,l,
     $     al,vc,harm,phi,freq,dx,dy,theta,v10,v20,v11,v02,
     $     fringe,mfring,autophi)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only: gammab
      use temw,only:tmulbs
      use mathfun
      implicit none
      real*8 eps,oneev
      parameter (eps=1.d-2)
      parameter (oneev=1.d0+3.83d-12)
      integer*4 l,ndiv,n,mfring
      real*8 trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 trans1(6,6)
      real*8 al,vc,harm,phi,freq,w,v,vn,vcn,p1,h1,
     $     aln,phis,phic,dhg,v1,t,phii,dh,a,
     $     h2,p2,pf,v2,dgb,rg2,
     $     dx,dy,theta,v10,v20,v11,v02,
     $     v10a,v11a,v20a,v02a,ddhdx,ddhdy,ddhdz,ddhdp,
     $     wi,offset1,va,sp,cp,av,dpxa,dpya,dav,davdz,davdp,
     $     dpx,dpy,dv,s0
      logical*4 fringe,autophi
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,0.d0,0.d0,.true.)
      if(harm .eq. 0.d0)then
        w=pi2*freq/c
      else
        w=omega0*harm/c
      endif
      if(w .eq. 0.d0)then
        wi=0.d0
      else
        wi=1.d0/w
      endif
      p0=gammab(l)
      h0=p2h(p0)
c      h0=p0*sqrt(1.d0+1.d0/p0**2)
      p1=gammab(l+1)
      h1=p2h(p1)
c      h1=p1*sqrt(1.d0+1.d0/p1**2)
c      h1=sqrt(p1**2+1.d0)
      v=vc/amass*abs(charge)
      ndiv=1+int(min(abs(w*al),sqrt((v*(1.d0/h0+1.d0/h1))**2/3.d0/eps)))
      vn=v/ndiv
      vcn=vc/ndiv
      aln=al/ndiv
      phis=trf0*w
      phic=phi*sign(1.d0,charge)
c      write(*,'(a,1p5g15.7)')'tcave ',phi,phis,phic,trf0
      v10a=v10/ndiv/amass*abs(charge)
      v11a=v11/ndiv/amass*abs(charge)
      v20a=v20/ndiv/amass*abs(charge)+vn*(w*(.5d0/p0+.5d0/p1))**2/4.d0
      v02a=v02/ndiv/amass*abs(charge)+vn*(w*(.5d0/p0+.5d0/p1))**2/4.d0
      if(p1 .ne. p0)then
        dhg=(p1-p0)*(p1+p0)/(h1+h0)/ndiv
      else
        dhg=0.d0
      endif
      vc0=vc0+vc
      if(omega0 .ne. 0.d0)then
        hvc0=hvc0+(c*w)/omega0*vc
      endif
      if(rfsw)then
        if(trpt .or. radcod .or. autophi)then
          offset1=0.d0
          s0=0.d0
        else
          s0=sin(phis)
          offset1=sin(phis)
        endif
        if(al .ne. 0.d0 .and. fringe
     $       .and. mfring .ge. 0 .and. mfring .ne. 2)then
          call tcavfrie(trans,cod,beam,al,v,w,phic,phis-phic,s0,p0,
     $         irad,irad .gt. 6 .or. calpol,autophi)
        endif
        call tinitr(trans1)
        dgb=0.d0
        do n=1,ndiv
         if(al .ne. 0.d0)then
            if(n .eq. 1)then
              call tdrife(trans,cod,beam,srot,aln*.5d0,
     $             0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
            else
              call tdrife(trans,cod,beam,srot,aln,
     $             0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
              call tgetdvh(dgb,dv)
              cod(5)=cod(5)+dv*aln
            endif
          endif
          p1=p0*(1.d0+cod(6))
          h1=p2h(p1)
c          h1=sqrt(1.d0+p1**2)
          h1=p1+1.d0/(h1+p1)
          v1=p1/h1
          t=-cod(5)/v1
          va=vn+(v10a+v20a*cod(1)+v11a*cod(3))*cod(1)
     $           +v02a*cod(3)**2
          if(autophi)then
            phii=phic
            sp=sin(phii)
            cp=cos(phii)
          else
            phii=w*t+phic-phis
            sp=sin(phii)
            cp=cos(phii)
            dvcacc=dvcacc+vcn*cp*w
            ddvcacc=ddvcacc+vcn*sp*w**2
          endif
          vcacc=vcacc-vcn*sp
          dh=max(oneev-h1,-va*(sp+offset1))
c          write(*,'(a,1p6g14.6)')'tcave ',vcacc,vcn,sp
          h2=h1+dh
          p2=h2p(h2)
c          p2=h2*sqrt(1.d0-1.d0/h2**2)
          pf    =(h2+h1)/(p2+p1)*dh
          p2=p1+pf
          v2=p2/h2
          a=-va*w*cp
          av=-(cp*wi-offset1*t)/p0
          dpxa=v10a+2.d0*v20a*cod(1)+v11a*cod(3)
          dpya=v11a*cod(1)+2.d0*v02a*cod(3)
          dpx=dpxa*av
          dpy=dpya*av
          dav=-(sp+offset1)/p0
          davdz=dav/v1
          davdp=dav/p1/h1**2*t*p0
          trans1(2,1)=2.d0*v20a*av
          trans1(2,3)=v11a*av
          trans1(2,5)=dpxa*davdz
          trans1(2,6)=dpxa*davdp
          trans1(4,1)=v11a*av
          trans1(4,3)=2.d0*v02a*av
          trans1(4,5)=dpya*davdz
          trans1(4,6)=dpya*davdp
          ddhdx=-(2.d0*v20a*cod(1)+v11a*cod(3))*(sp+offset1)
          ddhdy=-(2.d0*v02a*cod(3)+v11a*cod(1))*(sp+offset1)
          ddhdz=-a/v1
          ddhdp=-a*t/p1/h1**2
          trans1(6,1)=ddhdx/v2/p0
          trans1(6,3)=ddhdy/v2/p0
          trans1(6,5)=ddhdz/v2/p0
          trans1(6,6)=(v1+ddhdp)/v2
          trans1(5,1)=-ddhdx*t/p2/h2**2
          trans1(5,3)=-ddhdy*t/p2/h2**2
          trans1(5,5)=v2/v1-ddhdz*t/p2/h2**2
c          write(*,*)'tcave ',trans1(5,5),v1,v2
          trans1(5,6)=-t*p0*(ddhdp-dh*(h2*(h2+h1)+p1**2)/p1/h1**2)
     $         /p2/h2**2
c          trans1(5,5)=(p2+a*t/p2/h2)/h2/v1
c          trans1(5,6)=t*p0/h1**2/h2**2/p1/p2*(dh*(h2*(h2+h1)+p1**2)+a*t)
c          trans1(6,5)=-a/v1/v2/p0
c          trans1(6,6)=(p1-a*t/p1/h1)/h1/v2
          cod(2)=cod(2)+dpx
          cod(4)=cod(4)+dpy
          cod(6)=cod(6)+pf/p0
          cod(5)=-t*v2
          trans(:,1:irad)=matmul(trans1,trans(:,1:irad))
c          call tmultr(trans,trans1,irad)
          if(irad .gt. 6)then
            call tmulbs(beam,trans1,.true.)
          endif
          dgb=dgb+dhg
        enddo
        if(al .ne. 0.d0)then
          call tdrife(trans,cod,beam,srot,aln*.5d0,
     $         0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
          call tgetdvh(dgb,dv)
          cod(5)=cod(5)+dv*aln*.5d0
          if(fringe .and. mfring .ge. 0 .and. mfring .ne. 1)then
            call tcavfrie(trans,cod,beam,al,-v,w,phic,phis-phic,s0,p0,
     $           irad,irad .gt. 6 .or. calpol,autophi)
          endif
        endif
      else
        phii=phic
        sp=sin(phii)
        cp=cos(phii)
        dvcacc=dvcacc+vc*cp*w
        ddvcacc=ddvcacc+vc*sp*w**2
        vcacc=vcacc-vc*sp
        if(al .ne. 0.d0)then
          call tdrife(trans,cod,beam,srot,al,
     $         0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
        endif
      endif
      if(dhg .ne. 0.d0)then
        rg2=p0/gammab(l+1)
c        rg=sqrt(rg2)
        call tinitr(trans1)
        trans1(2,2)=rg2
        trans1(4,4)=rg2
        trans1(6,6)=rg2
        trans(2,1:irad)=trans(2,1:irad)*rg2
        trans(4,1:irad)=trans(4,1:irad)*rg2
        trans(6,1:irad)=trans(6,1:irad)*rg2
        if(irad .gt. 6)then
          call tmulbs(beam,trans1,.true.)
        endif
        cod(2)=cod(2)*rg2
        cod(4)=cod(4)*rg2
        cod(6)=(cod(6)+1.d0)*rg2-1.d0
        pgev=gammab(l+1)*amass
        call tphyzp
        call tesetdv(cod(6))
      endif
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,0.d0,0.d0,.false.)
c      write(*,'(a,i5,1p6g15.7)')'tcave ',l+1,dhg,rg2,
c     $     trans(5,5),trans(5,6),trans(6,5),trans(6,6)
      return
      end

      subroutine tcavfrie(trans,cod,beam,al,v,w,phic,dphis,s0,p0,
     $     irad,calb,autophi)
      use temw,only:tmulbs
      use mathfun, only:p2h
      implicit none
      real*8 trans(6,12),cod(6),trans1(6,6),beam(42),
     $     v,al,p0,vf,dp1r,p1r,p1,h1,v1,t,phic,
     $     ph,w,dphis,sph,cph,dpt,s0,wc,dha,dh,h2,a,dpr,
     $     dp2r,p2r,v2,at,az
      integer*4 irad
      logical*4 calb,autophi
      vf=v/al/p0*.5d0
      dp1r=cod(6)
      p1r=1.d0+dp1r
      p1=p0*p1r
      h1=p2h(p1)
c      h1=p1*sqrt(1.d0+1.d0/p1**2)
      v1=p1/h1
      t=-cod(5)/v1
      if(autophi)then
        ph=phic
      else
        ph=w*t-dphis
      endif
      sph=sin(ph)
      cph=cos(ph)
      dpt=vf*(sph+s0)
c      write(*,'(a,1p6g15.7)')'tcavfrie ',dpt,ph,phic,w,t,dphis
      cod(2)=cod(2)+cod(1)*dpt
      cod(4)=cod(4)+cod(3)*dpt
      wc=-w*vf
      dha=wc*(cod(1)**2+cod(3)**2)*.5d0
      dh=dha*cph*p0
      h2=h1+dh
      a=dh*(h1+h2)/p1**2
      dpr=a/(1.d0+sqrt(1.d0+a))
      dp2r=dp1r+p1r*dpr
      p2r=1.d0+dp2r
      v2=p2r*p0/h2
      cod(6)=dp2r
      call tinitr(trans1)
c      write(*,*)'tcavfrie ',dpt
      at=cod(5)/p1r/h1**2
      trans1(2,1)=dpt
      trans1(2,5)=-cod(1)*vf*w/v1*cph
      trans1(2,6)=-trans1(2,5)*at
      trans1(4,3)=dpt
      trans1(4,5)=-cod(3)*vf*w/v1*cph
      trans1(4,6)=-trans1(4,5)*at
      trans1(6,1)=wc*cod(1)/v2*cph
      trans1(6,3)=wc*cod(3)/v2*cph
      trans1(6,5)=w*dha*sph/v1/v2
      trans1(6,6)=v1/v2-trans1(6,5)*at
      az=cod(5)/v1/h2**3*p0
      trans1(5,1)=az*trans1(6,1)
      trans1(5,3)=az*trans1(6,3)
      trans1(5,5)=v2/v1+az*trans1(6,5)
      trans1(5,6)=az*trans1(6,6)-at*v2/v1
      cod(5)=cod(5)*v2/v1
      trans(:,1:irad)=matmul(trans1,trans(:,1:irad))
c      call tmultr(trans,trans1,irad)
      if(calb)then
        call tmulbs(beam,trans1,.true.)
      endif
      return
      end
