      subroutine tmulte(trans,cod,beam,l,al,ak,bz0,
     $     phia,psi1,psi2,apsi1,apsi2,
     1     dx,dy,dz,chi1,chi2,theta,dtheta,
     $     eps0,enarad,fringe,
     $     f1in,f2in,f1out,f2out,mfring,
     $     fb1,fb2,bfrm,vc,harm,phi,freq,wakew1,
     $     rtaper,autophi,ld)
      use tfstk
      use ffs_flag
      use ffs_pointer , only:gammab
      use tmacro
      use multa, only:nmult
      implicit none
      integer*4 ndivmax
      real*8 ampmax,oneev,pmax
      parameter (ampmax=0.05d0,ndivmax=300,pmax=0.9999d0)
      parameter (oneev=1.d0+3.83d-12)
      integer*4 mfring,ld,l,n,ndiv,m,kord,i,nmmax,nmmin
      real*8 f1in,f2in,f1out,f2out,
     $     al,vc,harm,phi,freq,bz,dx,dy,dz,chi1,chi2,theta,
     $     eps0,bxs,bys,bzs,al1,p1,h1,v1,t,phii,a,dh,dtheta,
     $     h2,p2,pf,v2,eps,v,w,aln,vn,phis,phic,ak1,vcn,veff,
     $     dhg,rg2,dgb,wakew1,w1n,theta1,phia,psi1,psi2,
     $     apsi1,apsi2,bz0,v10a,v11a,v20a,v02a,offset1,va,sp,cp,
     $     av,dpxa,dpya,dpx,dpy,dav,davdz,davdp,ddhdx,ddhdy,ddhdp,
     $     ddhdz,wi,dv,s0,fb1,fb2,rtaper,cod60,cod10,cod30,
     $     trans10(6,6)
      real*8 trans(6,12),trans1(6,6),cod(6),beam(42)
      complex*16 cx,cx0,cx2,cr,cr1
      real*8 fact(0:nmult),an(0:nmult)
      complex*16 ak(0:nmult),akn(0:nmult),ak0n
      logical*4 enarad,fringe,acc,bfrm,autophi
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1     720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $     39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $     1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $     6402373705728000.d0,121645100408832000.d0,
     $     2432902008176640000.d0,51090942171709440000.d0/
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
     $0.047619047619047619048d0/
      if(phia .ne. 0.d0)then
        call tmultae(trans,cod,beam,al,ak,
     $       phia,psi1,psi2,apsi1,apsi2,bz0,
     1       dx,dy,theta,dtheta,
     $       eps0,enarad,fringe,fb1,fb2,mfring,l,ld)
        return
      endif
      if(imag(ak(1)) .eq. 0.d0)then
        theta1=0.d0
      else
        theta1=atan2(imag(ak(1)),dble(ak(1)))*.5d0
      endif
      call tsolrot(trans,cod,beam,al,bz0,dx,dy,dz,
     $     chi1,chi2,theta+theta1,bxs,bys,bzs,.true.,ld)
      cr1=dcmplx(cos(theta1),-sin(theta1))
      akn(0)=(ak(0)*cr1+dcmplx(bys,bxs)*al)*rtaper
      bz=bz0
      do n=nmult,0,-1
        if(ak(n) .ne. 0.d0)then
          nmmax=n
          go to 1
        endif
      enddo
      if(vc .ne. 0.d0 .or. gammab(l+1) .ne. gammab(l))then
        nmmax=0
      else
        call tdrife(trans,cod,beam,al,bzs,dble(akn(0)),imag(akn(0)),
     $       .true.,enarad,calpol,irad,ld)
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
      ndiv=min(ndiv,ndivmax)
c      write(*,*)'tmulte-ndiv ',ndiv
      acc=vc .ne. 0.d0 .and. rfsw
      p0=gammab(l)
      h0=p2h(p0)
      p1=gammab(l+1)
      h1=p2h(p1)
      if(vc .ne. 0.d0)then
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
        vc0=vc0+vc
        if(omega0 .ne. 0.d0)then
          hvc0=hvc0+(c*w)/omega0*vc
        endif
        if(rfsw)then
          v=vc/amass*abs(charge)
          ndiv=max(ndiv,1+int(min(abs(w*al),
     $         sqrt((v*(1.d0/h0+1.d0/h1))**2/3.d0/eps))))
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
      else
        aln=al/ndiv
c     begin initialize for preventing compiler warning
        phic=0.d0
        phis=0.d0
        offset1=0.d0
        w=0.d0
        wi=0.d0
        vn=0.d0
        vcn=0.d0
        v20a=0.d0
        v02a=0.d0
c     end   initialize for preventing compiler warning
      endif
      if(p1 .ne. p0)then
        dhg=(p1-p0)*(p1+p0)/(h1+h0)/ndiv
      else
        dhg=0.d0
      endif
      cr=cr1*rtaper
      akn(0)=akn(0)/ndiv
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
            call tqfrie(trans,cod,beam,ak1,al1,ld,bzs)
          endif
        endif
        if(bfrm .and. ak0n .ne. (0.d0,0.d0))then
          if(mfring .eq. 1 .or. mfring .eq. 3)then
            call tbfrme(trans,cod,beam,ak0n/al1,fb1,.true.,ld)
          elseif(mfring .ne. 2)then
            call tbfrme(trans,cod,beam,ak0n/al1,0.d0,.true.,ld)
          endif
        endif
        if(mfring .eq. 1 .or. mfring .eq. 3)then
          call tqlfre(trans,cod,beam,al1,ak1,f1in,f2in,bzs,ld)
        endif
        nmmin=2
      else
        nmmin=1
      endif
      call tinitr(trans1)
      w1n=pbunch*abs(charge)*e*wakew1/amass/p0/ndiv
      dgb=0.d0
      do m=1,ndiv
        if(nmmin .eq. 2)then
          cod10=cod(1)
          cod30=cod(3)
          cod60=cod(6)
          trans10=trans(:,1:6)
          call tsolque(trans,cod,beam,al1,ak1,
     $         bzs,dble(ak0n),imag(ak0n),
     $         eps0,enarad,radcod,calpol,irad,ld)
          call tgetdvh(dgb,dv)
          cod(5)=cod(5)+dv*al1
c          if(abs(trans(2,3)).gt. 2.d0 .or.
c     $         abs(trans(1,3)) .gt. 2.d0)then
c            write(*,'(a,i5,1p8g14.6)')'tmulte-02 ',m,
c     $           cod(1),cod(3),cod(6),dy,ak1,ak(2)
c            write(*,'(1p6g14.6)')trans(:,1:6)
c            write(*,'(1p6g14.6)')trans10
c          endif
        endif
        ak1=dble(akn(1))
        al1=aln
        ak0n=akn(0)
        if(calpol)then
          call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
        endif
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
        if(acc)then
          p1=p0*(1.d0+cod(6))
          h1=p2h(p1)
c          h1=sqrt(1.d0+p1**2)
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
c          write(*,*)'tmulte-dh ',dh*amass,phii,phic,phis
          veff=vcn
          vc0=vc0+veff
          if(omega0 .ne. 0.d0)then
            hvc0=hvc0+(c*w)/omega0*veff
          endif
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
          call tmultr(trans,trans1,irad)
          dgb=dgb+dhg
        elseif(irad .eq. 6)then
            trans(2,1)=trans(2,1)+trans1(2,1)*trans(1,1)
     $           +trans1(2,3)*trans(3,1)
            trans(4,1)=trans(4,1)+trans1(4,1)*trans(1,1)
     $           +trans1(4,3)*trans(3,1)
            trans(2,2)=trans(2,2)+trans1(2,1)*trans(1,2)
     $           +trans1(2,3)*trans(3,2)
            trans(4,2)=trans(4,2)+trans1(4,1)*trans(1,2)
     $           +trans1(4,3)*trans(3,2)
            trans(2,3)=trans(2,3)+trans1(2,1)*trans(1,3)
     $           +trans1(2,3)*trans(3,3)
            trans(4,3)=trans(4,3)+trans1(4,1)*trans(1,3)
     $           +trans1(4,3)*trans(3,3)
            trans(2,4)=trans(2,4)+trans1(2,1)*trans(1,4)
     $           +trans1(2,3)*trans(3,4)
            trans(4,4)=trans(4,4)+trans1(4,1)*trans(1,4)
     $           +trans1(4,3)*trans(3,4)
            trans(2,5)=trans(2,5)+trans1(2,1)*trans(1,5)
     $           +trans1(2,3)*trans(3,5)
            trans(4,5)=trans(4,5)+trans1(4,1)*trans(1,5)
     $           +trans1(4,3)*trans(3,5)
            trans(2,6)=trans(2,6)+trans1(2,1)*trans(1,6)
     $           +trans1(2,3)*trans(3,6)
            trans(4,6)=trans(4,6)+trans1(4,1)*trans(1,6)
     $           +trans1(4,3)*trans(3,6)
        else
          do i=1,irad
            trans(2,i)=trans(2,i)+trans1(2,1)*trans(1,i)
     $           +trans1(2,3)*trans(3,i)
            trans(4,i)=trans(4,i)+trans1(4,1)*trans(1,i)
     $           +trans1(4,3)*trans(3,i)
          enddo
        endif
        if(irad .gt. 6 .or. calpol)then
          call tmulbs(beam ,trans1,.true.,.true.)
        endif
        if(calpol)then
          call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
        endif
      enddo
      if(nmmin .eq. 2)then
        call tsolque(trans,cod,beam,al1*.5d0,ak1*.5d0,
     $       bzs,dble(ak0n)*.5d0,imag(ak0n)*.5d0,
     $       eps0,enarad,radcod,calpol,irad,ld)
        call tgetdvh(dgb,dv)
        cod(5)=cod(5)+dv*al1*.5d0
      endif
      if(al .ne. 0.d0)then
        if(mfring .eq. 2 .or. mfring .eq. 3)then
          call tqlfre(trans,cod,beam,al1,ak1,-f1out,f2out,bzs,ld)
        endif
        if(bfrm .and. ak0n .ne. (0.d0,0.d0))then
          if(mfring .eq. 2 .or. mfring .eq. 3)then
            call tbfrme(trans,cod,beam,-ak0n/al1,fb2,.false.,ld)
          elseif(mfring .ne. 1)then
            call tbfrme(trans,cod,beam,-ak0n/al1,0.d0,.false.,ld)
          endif
        endif
        if(fringe .and. mfring .ne. 1)then
          if(ak1 .ne. 0.d0)then
            call tqfrie(trans,cod,beam,-ak1,al1,ld,bzs)
          endif
          if(acc)then
            call tcavfrie(trans,cod,beam,al,-v,w,phic,phis-phic,s0,p0,
     $           irad,irad .gt. 6 .or. calpol,autophi)
          endif
        endif
      endif
 1000 continue
      call tsolrot(trans,cod,beam,al,bz,dx,dy,dz,
     $     chi1,chi2,theta+theta1,bxs,bys,bzs,.false.,ld)
      if(dhg .ne. 0.d0)then
        rg2=p0/gammab(l+1)
c        rg=sqrt(rg2)
        call tinitr(trans1)
        trans1(2,2)=rg2
        trans1(4,4)=rg2
        trans1(6,6)=rg2
        call tmultr(trans,trans1,irad)
        if(irad .gt. 6 .or. calpol)then
          call tmulbs(beam,trans1,.true.,.true.)
        endif
        cod(2)=cod(2)*rg2
        cod(4)=cod(4)*rg2
        cod(6)=(cod(6)+1.d0)*rg2-1.d0
        pgev=gammab(l+1)*amass
        call tphyzp
        call tesetdv(cod(6))
      endif
      return
      end
