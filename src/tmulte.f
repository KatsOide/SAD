      subroutine tmulte(trans,cod,beam,srot,l,al,ak,bz,
     $     phia,psi1,psi2,apsi1,apsi2,
     1     dx,dy,dz,chi1,chi2,theta,dtheta,
     $     eps0,enarad,fringe,
     $     f1in,f2in,f1out,f2out,mfring,
     $     fb1,fb2,bfrm,vc,harm,phi,freq,wakew1,
     $     rtaper,autophi)
      use ffs_flag, only:calpol,radcod,rfsw,trpt
      use ffs_pointer , only:gammab
      use tmacro, only:amass,c,charge,ddvcacc,dvcacc,e,h0,hvc0,
     $     irad,omega0,p0,pbunch,pgev,trf0,vc0,vcacc
      use multa, only:nmult
      use temw,only:bsir0,tsetr0,tmulbs,code
      use sol, only:tsolrote
      use kradlib, only:tradke
      use mathfun
      use multa, only:fact,an
      use macmath
      implicit none
      integer*4 ,parameter ::ndivmax=300
      real*8 ,parameter::ampmax=0.05d0,pmax=0.9999d0,oneev=1.d0+3.83d-12
      integer*4 ,intent(in):: mfring,l
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: f1in,f2in,f1out,f2out,
     $     al,vc,harm,phi,freq,bz,dx,dy,dz,chi1,chi2,theta,dtheta,
     $     fb1,fb2,rtaper,wakew1,eps0,phia,psi1,psi2,
     $     apsi1,apsi2
      complex*16 ,intent(in):: ak(0:nmult)
      logical*4 ,intent(in):: enarad,fringe,autophi,bfrm
      real*8 bxs,bys,bzs,al1,p1,h1,v1,t,phii,a,dh,
     $     h2,p2,pf,v2,eps,v,w,aln,vn,phis,phic,ak1,vcn,veff,
     $     dhg,rg2,dgb,w1n,theta2,v10a,v11a,v02a,v20a,offset1,va,sp,cp,
     $     av,dpxa,dpya,dpx,dpy,dav,davdz,davdp,ddhdx,ddhdy,ddhdp,
     $     ddhdz,wi,dv,s0
      integer*4 n,ndiv,m,kord,nmmax,nmmin,itgetqraddiv
      real*8 trans1(6,6)
      complex*16 cx,cx0,cx2,cr,cr1
      complex*16 akn(0:nmult),ak0n,akn0
      logical*4 acc,krad
      if(phia .ne. 0.d0)then
        call tmultae(trans,cod,beam,srot,al,ak,
     $       phia,psi1,psi2,apsi1,apsi2,bz,
     1       dx,dy,theta,dtheta,
     $       eps0,enarad,fringe,fb1,fb2,mfring,l)
        return
      endif
      code=cod
      theta2=theta+dtheta+akang(ak(1),al,cr1)
      call tsolrote(trans,cod,beam,srot,al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.true.)
      akn0=(ak(0)*cr1+dcmplx(bys,bxs)*al)*rtaper
      krad=enarad .and. al .ne. 0.d0
      if(krad)then
        call tsetr0(trans(:,1:6),cod,bzs*.5d0,0.d0)
      endif
      do n=nmult,1,-1
        if(ak(n) .ne. 0.d0)then
          nmmax=n
          go to 1
        endif
      enddo
      if(vc .ne. 0.d0 .or. gammab(l+1) .ne. gammab(l))then
        nmmax=0
      else
        call tdrife(trans,cod,beam,srot,
     $       al,bzs,dble(akn0),imag(akn0),al,
     $       .true.,krad,irad)
        dhg=0.d0
        go to 1000
      endif
 1    if(eps0 .eq. 0.d0)then
        eps=5.d-3
      else
        eps=5.d-3*eps0
      endif
      ndiv=1
      do n=2,nmmax
        ndiv=max(ndiv,
     $       int(sqrt(ampmax**(n-1)
     $       /6.d0/fact(n-1)/eps*abs(ak(n)*al)))+1)
      enddo
      if(krad)then
        ndiv=min(ndivmax,max(ndiv,
     $       itgetqraddiv(cod,dble(ak(1)),al,bzs*.5d0)))
      else
        ndiv=min(ndiv,ndivmax)
      endif
      aln=al/ndiv
      acc=vc .ne. 0.d0 .and. rfsw
      p0=gammab(l)
      h0=p2h(p0)
      p1=gammab(l+1)
      h1=p2h(p1)
      offset1=0.d0
      vn=0.d0
      vcn=0.d0
      phis=0.d0
      v20a=0.d0
      phic=0.d0
      w=0.d0
      wi=0.d0
      v02a=0.d0
      if(acc)then
        if(harm .eq. 0.d0)then
          w=m_2pi*freq/c
        else
          w=omega0*harm/c
        endif
        if(w .eq. 0.d0)then
          wi=0.d0
        else
          wi=1.d0/w
        endif
        vc0=vc0+vc
        if(omega0 .ne. 0.d0)then
          hvc0=hvc0+(c*w)/omega0*vc
        endif
        if(rfsw)then
          v=vc/amass*abs(charge)
          ndiv=min(ndivmax,max(ndiv,1+int(min(abs(w*al),
     $         sqrt((v*(1.d0/h0+1.d0/h1))**2/3.d0/eps)))))
          aln=al/ndiv
          vn=v/ndiv
          vcn=vc/ndiv
          phis=trf0*w
          phic=phi*sign(1.d0,charge)
          v10a=0.d0
          v11a=0.d0
          v20a=vn*(w*(.5d0/p0+.5d0/p1))**2/4.d0
          v02a=vn*(w*(.5d0/p0+.5d0/p1))**2/4.d0
          if(trpt .or. radcod .or. autophi)then
            s0=0.d0
            offset1=0.d0
          else
            s0=sin(phis)
            offset1=sin(phis)
          endif
        else
          phii=phic
          sp=sin(phii)
          cp=cos(phii)
          dvcacc=dvcacc+vc*cp*w
          ddvcacc=ddvcacc+vc*sp*w**2
          vcacc=vcacc-vc*sp
        endif
      endif
      if(p1 .ne. p0)then
        dhg=(p1-p0)*(p1+p0)/(h1+h0)/ndiv
      else
        dhg=0.d0
      endif
      cr=cr1*rtaper
      akn(0)=akn0/ndiv
      do n=1,max(nmmax,1)
        cr=cr*cr1
        akn(n)=(ak(n)*cr)/ndiv
      enddo
      akn(1)=dble(akn(1))
      ak1=dble(akn(1))*.5d0
      al1=aln*.5d0
      ak0n=akn(0)*.5d0
      if(al .ne. 0.d0)then
        if(fringe .and. mfring .ne. 2)then
          if(acc)then
            call tcavfrie(trans,cod,beam,al,v,w,phic,phis-phic,s0,p0,
     $           irad,irad .gt. 6 .or. calpol,autophi)
          endif
          if(ak1 .ne. 0.d0)then
            call tqfrie(trans,cod,beam,ak1,al1,bzs)
          endif
        endif
        if(bfrm .and. ak0n .ne. (0.d0,0.d0))then
          if(mfring .eq. 1 .or. mfring .eq. 3)then
            call tbfrme(trans,cod,beam,ak0n/al1,fb1,.true.)
          elseif(mfring .ne. 2)then
            call tbfrme(trans,cod,beam,ak0n/al1,0.d0,.true.)
          endif
        endif
        if(mfring .eq. 1 .or. mfring .eq. 3)then
          call tqlfre(trans,cod,beam,al1,ak1,f1in,f2in,bzs)
        endif
        nmmin=2
        if(krad)then
          if(f1in .ne. 0.d0)then
            call tradke(trans,cod,beam,srot,f1in,0.d0,bzs*.5d0)
          else
            call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
          endif
        endif
      else
        call tsetr0(trans(:,1:6),cod(1:6),bzs*.5d0,0.d0)
        nmmin=1
      endif
      call tinitr(trans1)
      w1n=pbunch*abs(charge)*e*wakew1/amass/p0/ndiv
      dgb=0.d0
      do m=1,ndiv
        if(nmmin .eq. 2)then
          call tsolque(trans,cod,beam,srot,al1,ak1,
     $         bzs,dble(ak0n),imag(ak0n),
     $         eps0,krad,irad)
          call tgetdvh(dgb,dv)
          cod(5)=cod(5)+dv*al1
        endif
        ak1=dble(akn(1))
        al1=aln
        ak0n=akn(0)
        cx0=dcmplx(cod(1),cod(3))
        cx=(0.d0,0.d0)
        cx2=(0.d0,0.d0)
        do kord=nmmax,nmmin,-1
          cx=(cx+akn(kord))*cx0*an(kord)
          cx2=cx2*cx0*an(kord)+akn(kord)
        enddo
        if(nmmin .eq. 2)then
          cx=cx*cx0
          cx2=cx2*cx0
        else
          cx=cx+akn(0)
        endif
        cx=dcmplx(min(pmax,max(-pmax,dble(cx))),
     $       min(pmax,max(-pmax,imag(cx))))
        trans1(2,1)=-dble(cx2)+w1n
        trans1(2,3)=imag(cx2)
        trans1(4,1)=imag(cx2)
        trans1(4,3)= dble(cx2)+w1n
        cod(2)=cod(2)-dble(cx)+w1n*cod(1)
        cod(4)=cod(4)+imag(cx)+w1n*cod(3)
        if(m .eq. 1)then
          bsir0=bsir0+imag(cx)/al1
        else
          bsir0=0.d0
        endif
        if(m .eq. ndiv)then
          bsir0=bsir0-imag(cx)/al1
        endif
        if(acc)then
          p1=p0*(1.d0+cod(6))
          h1=p2h(p1)
          h1=p1+1.d0/(h1+p1)
          v1=p1/h1
          t=-cod(5)/v1
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
          va=vn+(v10a+v20a*cod(1)+v11a*cod(3))*cod(1)
     $           +v02a*cod(3)**2
          dh=max(oneev-h1,-va*(sp+offset1))
          veff=vcn
          vc0=vc0+veff
          if(omega0 .ne. 0.d0)then
            hvc0=hvc0+(c*w)/omega0*veff
          endif
          h2=h1+dh
          p2=h2p(h2)
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
          trans1(1,1)=1.d0
          trans1(2,1)=trans1(2,1)+2.d0*v20a*av
          trans1(2,2)=1.d0
          trans1(2,3)=trans1(2,3)+v11a*av
          trans1(2,5)=dpxa*davdz
          trans1(2,6)=dpxa*davdp
          trans1(3,3)=1.d0
          trans1(4,1)=trans1(4,1)+v11a*av
          trans1(4,3)=trans1(4,3)+2.d0*v02a*av
          trans1(4,4)=1.d0
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
          trans1(5,6)=-t*p0*(ddhdp-dh*(h2*(h2+h1)+p1**2)/p1/h1**2)
     $         /p2/h2**2
          cod(2)=cod(2)+dpx
          cod(4)=cod(4)+dpy
          cod(6)=cod(6)+pf/p0
          cod(5)=-t*v2
          trans(:,1:irad)=matmul(trans1,trans(:,1:irad))
          dgb=dgb+dhg
c$$$        elseif(irad .eq. 6)then
c$$$          trans(2,1)=trans(2,1)+trans1(2,1)*trans(1,1)
c$$$     $         +trans1(2,3)*trans(3,1)
c$$$          trans(4,1)=trans(4,1)+trans1(4,1)*trans(1,1)
c$$$     $         +trans1(4,3)*trans(3,1)
c$$$          trans(2,2)=trans(2,2)+trans1(2,1)*trans(1,2)
c$$$     $         +trans1(2,3)*trans(3,2)
c$$$          trans(4,2)=trans(4,2)+trans1(4,1)*trans(1,2)
c$$$     $         +trans1(4,3)*trans(3,2)
c$$$          trans(2,3)=trans(2,3)+trans1(2,1)*trans(1,3)
c$$$     $         +trans1(2,3)*trans(3,3)
c$$$          trans(4,3)=trans(4,3)+trans1(4,1)*trans(1,3)
c$$$     $         +trans1(4,3)*trans(3,3)
c$$$          trans(2,4)=trans(2,4)+trans1(2,1)*trans(1,4)
c$$$     $         +trans1(2,3)*trans(3,4)
c$$$          trans(4,4)=trans(4,4)+trans1(4,1)*trans(1,4)
c$$$     $         +trans1(4,3)*trans(3,4)
c$$$          trans(2,5)=trans(2,5)+trans1(2,1)*trans(1,5)
c$$$     $         +trans1(2,3)*trans(3,5)
c$$$          trans(4,5)=trans(4,5)+trans1(4,1)*trans(1,5)
c$$$     $         +trans1(4,3)*trans(3,5)
c$$$          trans(2,6)=trans(2,6)+trans1(2,1)*trans(1,6)
c$$$     $         +trans1(2,3)*trans(3,6)
c$$$          trans(4,6)=trans(4,6)+trans1(4,1)*trans(1,6)
c$$$     $         +trans1(4,3)*trans(3,6)
        else
          trans(2,1:irad)=trans(2,1:irad)+trans1(2,1)*trans(1,1:irad)
     $         +trans1(2,3)*trans(3,1:irad)
          trans(4,1:irad)=trans(4,1:irad)+trans1(4,1)*trans(1,1:irad)
     $         +trans1(4,3)*trans(3,1:irad)
        endif
        if(irad .gt. 6)then
          call tmulbs(beam ,trans1,.true.)
        endif
      enddo
      if(nmmin .eq. 2)then
        call tsolque(trans,cod,beam,srot,al1*.5d0,ak1*.5d0,
     $       bzs,dble(ak0n)*.5d0,imag(ak0n)*.5d0,
     $       eps0,krad,irad)
        call tgetdvh(dgb,dv)
        cod(5)=cod(5)+dv*al1*.5d0
      endif
      if(al .ne. 0.d0)then
        if(mfring .eq. 2 .or. mfring .eq. 3)then
          call tqlfre(trans,cod,beam,al1,ak1,-f1out,f2out,bzs)
        endif
        if(bfrm .and. ak0n .ne. (0.d0,0.d0))then
          if(mfring .eq. 2 .or. mfring .eq. 3)then
            call tbfrme(trans,cod,beam,-ak0n/al1,fb2,.false.)
          elseif(mfring .ne. 1)then
            call tbfrme(trans,cod,beam,-ak0n/al1,0.d0,.false.)
          endif
        endif
        if(fringe .and. mfring .ne. 1)then
          if(ak1 .ne. 0.d0)then
            call tqfrie(trans,cod,beam,-ak1,al1,bzs)
          endif
          if(acc)then
            call tcavfrie(trans,cod,beam,al,-v,w,phic,phis-phic,s0,p0,
     $           irad,irad .gt. 6,autophi)
          endif
        endif
      endif
      if(krad .and. f1out .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1out,0.d0,bzs*.5d0)
      endif
 1000 continue
c      if(ktfenanq(cod(5)) .or. ktfenanq(trans(5,6)))then
c        write(*,'(a,2i5)')'tmulte-1000 ',l,ndiv
c        write(*,'(1p6g15.7)')(trans(i,1:6),i=1,6),cod
c      endif
      call tsolrote(trans,cod,beam,srot,al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.false.)
      if(dhg .ne. 0.d0)then
        rg2=p0/gammab(l+1)
        call tinitr(trans1)
        trans1(2,2)=rg2
        trans1(4,4)=rg2
        trans1(6,6)=rg2
        trans(:,1:irad)=matmul(trans1,trans(:,1:irad))
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
c      if(ktfenanq(cod(5)) .or. ktfenanq(trans(5,6)))then
c        write(*,'(a,i5)')'tmulte-end ',l
c        write(*,'(1p6g15.7)')(trans(i,1:6),i=1,6),cod
c      endif
      return
      end
