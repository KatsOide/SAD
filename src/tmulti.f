      subroutine tmulti(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,ak,akr0,bz,phia,psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,theta2,
     $     eps0,krad,fringe,f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,dofr,
     $     radius,rtaper,ndiv0,nmmax,
     $     kturn,kptbl)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      use kradlib
      use sol,only:tsolrot
      use photontable,only:tsetpcvt,pcvt
      use mathfun
      use multa, only:fact,aninv
      use kyparam, only:nmult
c      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 , parameter ::itmax=10,ndivmax=1000
      real*8 ,parameter :: conv=3.d-16,oneev=1.d0+1.d-6,
     $     alstep=0.05d0,pmin=1.d-10,arad=0.01d0
      integer*4 ,intent(inout):: np
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np),
     $     dv(np),sx(np),sy(np),sz(np),bz
      real*8 ,intent(in):: al,phia,psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,theta2,
     $     eps0,f1in,f2in,f1out,f2out,fb1,fb2,radius,rtaper
      complex*16 ,intent(in):: ak(0:nmult),akr0(0:nmult)
      integer*4 ,intent(in):: mfring,kturn,nmmax,ndiv0
      integer*4 ,intent(inout):: kptbl(np0,6)
      logical*4 ,intent(in):: fringe,krad
      logical*1 ,intent(in):: dofr(0:nmult)
      logical*4 spac1,nzleng
      integer*4 i,m,n,ndiv,ibsi
      real*8 bxs,bys,bzs,b,a,epsr,wi,
     $     akr1,ak1,al1,p,ea,pxf,pyf,sv,wm,alx
      complex*16 akr(0:nmult),akrm(0:nmult),cx,cx1,ak01,b0
      if(phia .ne. 0.d0)then
        call tmulta(
     $       np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,ak,phia,
     $       psi1,psi2,bz,
     1       dx,dy,theta,dtheta,
     $       eps0,krad,fb1,fb2,mfring,fringe)
        return
      endif
      b0=0.d0
      nzleng=al .ne. 0.d0
      spac1=.false.
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.true.)
      akr(0)=(akr0(0)+dcmplx(bys*al,bxs*al))*rtaper
      if(nmmax .eq. 0 .and. .not. spac)then
        call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,bzs,dble(akr(0)),imag(akr(0)),krad)
      else
c     Zero-clear akr(1) for case: nmmax .eq. 0
        akr(1)=0.d0
        akr(1:nmmax)=akr0(1:nmmax)*rtaper
c     Im[ark(1)] should be 0, because akr(1) := ak(1)*cr1*cr1,
c     cr1 := Exp[-theta1], ak(1) = Abs[ak(1)] * Exp[2 theta1]
        akr1=dble(akr(1))
        akr(1)=akr1
        ndiv=ndiv0
        epsr=merge(1.d0,eps0,eps0 .eq. 0.d0)
        if(nzleng)then
          if(spac)then
            spac1 = radius .ne. 0.d0
            ndiv=max(ndiv,nint(abs(al)/(alstep*epsr)),
     $           nint(abs(bzs*al)/1.5d0/epsr))
          endif
          if(krad)then
            pxr0=px
            pyr0=py
            zr0=z
            ndiv=max(ndiv,ndivrad(abs(akr(0)),akr1,bz,eps0))
            if(photons)then
              call tsetpcvt(l_track,dx,dy,theta2,0.d0,0.d0,al)
            endif
          endif
          ndiv=min(ndivmax,ndiv)
          if(fringe)then
            if(mfring .ne. 2)then
              do n=0,nmmax
                if(dofr(n))then
                  call ttfrins(np,x,px,y,py,z,g,n*2+2,akr(n),
     $                 al,bzs)
                endif
              enddo
            endif
          endif
          if(mfring .eq. 1 .or. mfring .eq. 3)then
            if(akr(0) .ne. (0.d0,0.d0) .and. fb1 .ne. 0.d0)then
              call tblfri(np,x,px,y,py,z,g,al,akr(0),fb1)
            endif
            if(f1in .ne. 0.d0 .or. f2in .ne. 0.d0)then
              do concurrent (i=1:np)
                p=(1.d0+g(i))
                a=f1in/p
                ea=exp(a)
                b=f2in/p
                pxf=px(i)/ea
                pyf=py(i)*ea
                z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $               +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
                x(i)=ea*x(i)+b*px(i)
                y(i)=y(i)/ea-b*py(i)
                px(i)=pxf
                py(i)=pyf
              enddo
            endif
          endif
        endif
        wi=1.d0/ndiv
        ak1=akr1*wi*.5d0
        al1=al*wi*.5d0
        ak01=akr(0)*wi*.5d0
        sv=0.d0
        ibsi=1
        do m=1,ndiv
          akrm(0:nmmax)=akr(0:nmmax)*wi
          if(nzleng)then
            call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           al1,ak1,bzs,dble(ak01),imag(ak01),ibsi,eps0)
            if(krad)then
              if(m .eq. 1 .and. calpol)then
                do concurrent (i=1:np)
                  cx1=dcmplx(x(i),y(i))
                  cx=0.d0
                  do n=nmmax,2,-1
                    cx=(cx+akr(n))*cx1*aninv(n+1)
                  enddo
                  bsi(i)=bsi(i)+imag(.5d0*cx*cx1**2)/al
                enddo
              endif
              call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
              pcvt%fr0=pcvt%fr0+al1/al
            endif
            ibsi=0
            wm=wi
            if(m .eq. ndiv)then
              wm=wi*.5d0
            endif
            al1=al*wm
            ak1=akr1*wm
            ak01=akr(0)*wm
            if(nmmax .ge. 2)then
              do concurrent (i=1:np)
                cx1=dcmplx(x(i),y(i))
                cx=akrm(nmmax)*cx1*aninv(nmmax)
                do n=nmmax-1,2,-1
                  cx=(cx+akrm(n))*cx1*aninv(n)
                enddo
                cx=cx*cx1/(1.d0+g(i))
                px(i)=px(i)-dble(cx)
                py(i)=py(i)+imag(cx)
              enddo
            endif
            if(spac1)then
              call spkick(np,x,px,y,py,z,g,dv,sx,sy,sz,al*wm,radius,
     $             alx,kturn,kptbl)
            endif
          elseif(nmmax .ge. 1)then
            do concurrent (i=1:np)
              cx1=dcmplx(x(i),y(i))
              cx=akrm(nmmax)*cx1*aninv(nmmax)
              do n=nmmax-1,1,-1
                cx=(cx+akrm(n))*cx1*aninv(n)
              enddo
              cx=(cx+akrm(0))/(1.d0+g(i))
              px(i)=px(i)-dble(cx)
              py(i)=py(i)+imag(cx)
            enddo
          endif
        enddo
        if(nzleng)then
          call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al1,ak1,bzs,dble(ak01),imag(ak01),2,eps0)
          if(mfring .eq. 2 .or. mfring .eq. 3)then
            if(f1out .ne. 0.d0 .or. f2out .ne. 0.d0)then
              do concurrent (i=1:np)
                p=(1.d0+g(i))
                a=-f1out/p
                ea=exp(a)
                b= f2out/p
                pxf=px(i)/ea
                pyf=py(i)*ea
                z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $               +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
                x(i)=ea*x(i)+b*px(i)
                y(i)=y(i)/ea-b*py(i)
                px(i)=pxf
                py(i)=pyf
              enddo
            endif
            if(fb2 .ne. 0.d0 .and. akr(0) .ne. (0.d0,0.d0))then
              call tblfri(np,x,px,y,py,z,g,al,-akr(0),fb2)
            endif
          endif
          if(fringe .and. mfring .ne. 1)then
            do n=nmmax,0,-1
              if(dofr(n))then
                call ttfrins(np,x,px,y,py,z,g,n*2+2,-akr(n),
     $               al,bzs)
              endif
            enddo
          endif
          if(krad)then
            if(calpol .and. nmmax .ge. 2)then
              do concurrent (i=1:np)
                cx1=dcmplx(x(i),y(i))
                cx=akr(nmmax)*cx1*aninv(nmmax+1)
                do n=nmmax-1,2,-1
                  cx=(cx+akr(n))*cx1*aninv(n+1)
                enddo
                bsi(i)=bsi(i)-imag(.5d0*cx*cx1**2)/al
              enddo
            endif
            call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
            pcvt%fr0=pcvt%fr0+al1/al
          endif
        endif
        if(spac1)then
          call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $         kptbl,np,kturn,
     $         radius,radius,
     $         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
        endif
      endif
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.false.)
      return
      end

      subroutine tmultiacc(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,ak,akr0,bz,phia,psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,theta2,
     $     eps0,krad,fringe,f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,dofr,
     $     vc,w,phirf,dphirf,vnominal,
     $     radius,rtaper,autophi,ndiv0,nmmax,
     $     kturn,kptbl)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      use kradlib
      use sol,only:tsolrot
      use photontable,only:tsetpcvt,pcvt
      use mathfun
      use multa, only:fact,aninv
      use kyparam, only:nmult
c      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 , parameter ::itmax=10,ndivmax=1000
      real*8 ,parameter :: conv=3.d-16,oneev=1.d0+1.d-6,
     $     alstep=0.05d0,pmin=1.d-10,arad=0.01d0
      integer*4 ,intent(inout):: np
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np),
     $     dv(np),sx(np),sy(np),sz(np),bz
      real*8 ,intent(in):: al,phia,psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,theta2,
     $     eps0,f1in,f2in,f1out,f2out,fb1,fb2,w,
     $     vc,phirf,dphirf,vnominal,radius,rtaper
      complex*16 ,intent(in):: ak(0:nmult),akr0(0:nmult)
      integer*4 ,intent(in):: mfring,kturn,nmmax,ndiv0
      integer*4 ,intent(inout):: kptbl(np0,6)
      logical*4 ,intent(in):: fringe,autophi,krad
      logical*1 ,intent(in):: dofr(0:nmult)
      logical*4 spac1,nzleng
      integer*4 i,m,n,ndiv,ibsi
      real*8 bxs,bys,bzs,b,a,epsr,wi,v,phis,r,wl,
     $     r1,we,phic,dphis,offset,offset1,
     $     tlim,akr1,ak1,al1,p,ea,pxf,pyf,sv,asinh,ws1,wm,
     $     h2,p2,dp2,pr2,dvn,dzn,dp1r,p1r,p1,h1,t,ph,dh,dpr,dp2r,p2r,
     $     alx,dp2p2,dp,dp1,pr1,
     $     he,vcorr,v20a,v02a,v1a,v11a,av,dpx,dpy,pe,ah
      real*8 ws(ndivmax+1)
      complex*16 akr(0:nmult),akrm(0:nmult),cx,cx1,ak01,b0
      if(phia .ne. 0.d0)then
        call tmulta(
     $       np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,ak,phia,
     $       psi1,psi2,bz,
     1       dx,dy,theta,dtheta,
     $       eps0,krad,fb1,fb2,mfring,fringe)
        return
      endif
      dphis=0.d0
      b0=0.d0
      nzleng=al .ne. 0.d0
      spac1=.false.
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.true.)
      akr(0)=(akr0(0)+dcmplx(bys*al,bxs*al))*rtaper
      if(nmmax .eq. 0 .and. vc .eq. 0.d0 .and. .not. spac)then
        call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,bzs,dble(akr(0)),imag(akr(0)),krad)
        go to 1000
      endif
c     Zero-clear akr(1) for case: nmmax .eq. 0
      akr(1)=0.d0
      akr(1:nmmax)=akr0(1:nmmax)*rtaper
c     Im[ark(1)] should be 0, because akr(1) := ak(1)*cr1*cr1,
c     cr1 := Exp[-theta1], ak(1) = Abs[ak(1)] * Exp[2 theta1]
      akr1=dble(akr(1))
      akr(1)=akr1
      ndiv=ndiv0
      epsr=merge(1.d0,eps0,eps0 .eq. 0.d0)
      if(nzleng)then
        if(spac)then
          spac1 = radius .ne. 0.d0
          ndiv=max(ndiv,nint(abs(al)/(alstep*epsr)),
     $         nint(abs(bzs*al)/1.5d0/epsr))
        endif
        if(krad)then
          pxr0=px
          pyr0=py
          zr0=z
          ndiv=max(ndiv,ndivrad(abs(akr(0)),akr1,bz,eps0))
          if(photons)then
            call tsetpcvt(l_track,dx,dy,theta2,0.d0,0.d0,al)
          endif
        endif
        ndiv=min(ndivmax,ndiv)
      endif
      wi=merge(0.d0,1.d0/w,w .eq. 0.d0)
      v=vc*abs(charge)/amass
      he=h0+vnominal
      pe=h2p(he)
      vcorr=v*(w*(.5d0/p0+.5d0/pe))**2/4.d0
      v20a=vcorr
      v02a=vcorr
      v1a=0.d0
      v11a=0.d0
      r=abs(vnominal*(1.d0/h0+1.d0/he))
      ndiv=min(ndivmax,max(ndiv,
     $     1+int(min(abs(w*al),
     $     sqrt(r**2/3.d0/(epsr*1.d-3))))))
      if(abs(r) .gt. 1.d-6)then
        wl=asinh(r*(2.d0+r)/2.d0/(1.d0+r))
        ndiv=min(ndivmax,max(ndiv,nint(wl*10.d0/epsr)))
        r1=wl/ndiv/2.d0
        ws1=2.d0*exp(r1)*sinh(r1)
        we=1.d0+ws1
        ws(1)=ws1/r
        do i=2,ndiv
          ws(i)=ws(i-1)*we
        enddo
      else
        ws(1:ndiv)=1.d0/ndiv
      endif
      phic=(phirf+dphirf)*charge
      phis=merge(0.d0,merge(phic,w*trf0,autophi),trpt)
      dphis=phis-phic
      if(rad .or. trpt .or. autophi)then
        offset=sin(dphis)
        offset1=0.d0
      else
        offset=sin(dphis)-sin(phis)
        offset1=sin(phis)
      endif
      tlim=1.d4
      ws(ndiv+1)=0.d0
      ak1=akr1*ws(1)*.5d0
      al1=al*ws(1)*.5d0
      ak01=akr(0)*ws(1)*.5d0
      if(nzleng)then
        if(fringe)then
          if(mfring .ne. 2)then
            call tcavfrin(np,x,px,y,py,z,g,dv,al,v,w,p0,h0,
     $           dphis,dvfs,offset)
            do n=0,nmmax
              if(dofr(n))then
                call ttfrins(np,x,px,y,py,z,g,n*2+2,akr(n),
     $               al,bzs)
              endif
            enddo
          endif
        endif
        if(mfring .eq. 1 .or. mfring .eq. 3)then
          if(akr(0) .ne. (0.d0,0.d0) .and. fb1 .ne. 0.d0)then
            call tblfri(np,x,px,y,py,z,g,al,akr(0),fb1)
          endif
          if(f1in .ne. 0.d0 .or. f2in .ne. 0.d0)then
            do concurrent (i=1:np)
              p=(1.d0+g(i))
              a=f1in/p
              ea=exp(a)
              b=f2in/p
              pxf=px(i)/ea
              pyf=py(i)*ea
              z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
              x(i)=ea*x(i)+b*px(i)
              y(i)=y(i)/ea-b*py(i)
              px(i)=pxf
              py(i)=pyf
            enddo
          endif
        endif
      endif
      sv=0.d0
      ibsi=1
      do m=1,ndiv
        akrm(0:nmmax)=akr(0:nmmax)*ws(m)
        if(nzleng)then
          call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al1,ak1,bzs,dble(ak01),imag(ak01),ibsi,eps0)
          if(krad)then
            if(m .eq. 1 .and. calpol)then
              do concurrent (i=1:np)
                cx1=dcmplx(x(i),y(i))
                cx=0.d0
                do n=nmmax,2,-1
                  cx=(cx+akr(n))*cx1*aninv(n+1)
                enddo
                bsi(i)=bsi(i)+imag(.5d0*cx*cx1**2)/al
              enddo
            endif
            call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
            pcvt%fr0=pcvt%fr0+al1/al
          endif
          ibsi=0
          wm=.5d0*(ws(m)+ws(m+1))
          al1=al*wm
          ak1=akr1*wm
          ak01=akr(0)*wm
          if(nmmax .ge. 2)then
            do concurrent (i=1:np)
              cx1=dcmplx(x(i),y(i))
              cx=akrm(nmmax)*cx1*aninv(nmmax)
              do n=nmmax-1,2,-1
                cx=(cx+akrm(n))*cx1*aninv(n)
              enddo
              cx=cx*cx1/(1.d0+g(i))
              px(i)=px(i)-dble(cx)
              py(i)=py(i)+imag(cx)
            enddo
          endif
          if(spac1)then
            call spkick(np,x,px,y,py,z,g,dv,sx,sy,sz,al*ws(m),radius,
     $           alx,kturn,kptbl)
          endif
        elseif(nmmax .ge. 1)then
          do concurrent (i=1:np)
            cx1=dcmplx(x(i),y(i))
            cx=akrm(nmmax)*cx1*aninv(nmmax)
            do n=nmmax-1,1,-1
              cx=(cx+akrm(n))*cx1*aninv(n)
            enddo
            cx=(cx+akrm(0))/(1.d0+g(i))
            px(i)=px(i)-dble(cx)
            py(i)=py(i)+imag(cx)
          enddo
        endif
        if(vnominal .ne. 0.d0)then
          sv=sv+vnominal*ws(m)
          h2=h0+sv
          p2=h2p(h2)
          dp2=sv*(h2+h0)/(p2+p0)/p0
          pr2=1.d0+dp2
          dvn=-dp2*(1.d0+pr2)/h2/(h2+pr2*h0)
          dzn=-dvn*al*ws(m)*.5d0
        else
          dzn=0.d0
        endif
        do concurrent (i=1:np)
          dp1r=g(i)
          p1r=1.d0+dp1r
          p1=p0*p1r
          h1=merge(p2h(p1),p1r*h0/(1.d0-dv(i)+dvfs),
     $         dv(i) .gt. 0.1d0)
          z(i)=z(i)-dzn
          t=min(tlim,max(-tlim,-z(i)*h1/p1))
          ph=.5d0*w*t
          dh=max(oneev-h1,
     $         (v+(v1a+v20a*x(i)+v11a*y(i))*x(i)+v02a*y(i)**2)
     $         *(-2.d0*sin(ph)*cos(ph-dphis)+offset)*ws(m))
          h2=h1+dh
          ah=max(dh*(h1+h2)/p1**2,-1.d0+pmin)
          dpr=sqrt1(ah)
          dp2r=max(dp1r+p1r*dpr,pmin-1.d0)
          p2r=1.d0+dp2r
          g(i)=dp2r
          dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs
          z(i)=-t*p2r*p0/h2-dzn
          av=-(cos(2.d0*ph-dphis)*wi-offset1*t)*ws(m)/p0
          dpx=(v1a+2.d0*v20a*x(i)+v11a*y(i))*av
          dpy=(v11a*x(i)+2.d0*v02a*y(i))*av
          px(i)=(px(i)*p1r+dpx)/p2r
          py(i)=(py(i)*p1r+dpy)/p2r
        enddo
      enddo
      if(nzleng)then
        call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al1,ak1,bzs,dble(ak01),imag(ak01),2,eps0)
        if(mfring .eq. 2 .or. mfring .eq. 3)then
          if(f1out .ne. 0.d0 .or. f2out .ne. 0.d0)then
            do concurrent (i=1:np)
              p=(1.d0+g(i))
              a=-f1out/p
              ea=exp(a)
              b= f2out/p
              pxf=px(i)/ea
              pyf=py(i)*ea
              z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
              x(i)=ea*x(i)+b*px(i)
              y(i)=y(i)/ea-b*py(i)
              px(i)=pxf
              py(i)=pyf
            enddo
          endif
          if(fb2 .ne. 0.d0 .and. akr(0) .ne. (0.d0,0.d0))then
            call tblfri(np,x,px,y,py,z,g,al,-akr(0),fb2)
          endif
        endif
        if(fringe .and. mfring .ne. 1)then
          do n=nmmax,0,-1
            if(dofr(n))then
              call ttfrins(np,x,px,y,py,z,g,n*2+2,-akr(n),
     $             al,bzs)
            endif
          enddo
          call tcavfrin(np,x,px,y,py,z,g,dv,al,-v,w,p0,h0,
     $         dphis,dvfs,offset)
        endif
        if(krad)then
          if(calpol .and. nmmax .ge. 2)then
            do concurrent (i=1:np)
              cx1=dcmplx(x(i),y(i))
              cx=akr(nmmax)*cx1*aninv(nmmax+1)
              do n=nmmax-1,2,-1
                cx=(cx+akr(n))*cx1*aninv(n+1)
              enddo
              bsi(i)=bsi(i)-imag(.5d0*cx*cx1**2)/al
            enddo
          endif
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
          pcvt%fr0=pcvt%fr0+al1/al
        endif
      endif
      if(spac1)then
        call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $       kptbl,np,kturn,
     $       radius,radius,
     $       0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
 1000 continue
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.false.)
      if(vnominal .ne. 0.d0)then
        h2=h0+vnominal
        p2=h2p(h2)
        dp2p2=vnominal*(h2+h0)/(p2+p0)/p2
        do concurrent (i=1:np)
          dp=max(g(i),pmin-1.d0)
          dp1=dp*p0/p2-dp2p2
          pr1=1.d0+dp1
          g(i)=dp1
          h1=p2h(p2*pr1)
          dv(i)=-dp1*(1.d0+pr1)/h1/(h1+pr1*h2)+dvfs
        enddo
        bz=bz*p0/p2
        pgev=p2*amass
        call tphyzp
      endif
      return
      end

      subroutine tblfri(np,x,px,y,py,z,g,al,ck0,fb1)
      implicit none
      integer*4 ,intent(in):: np
      integer*4 i
      complex*16 ,intent(in):: ck0
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ,intent(in):: al,fb1
      real*8 dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray,p,rhob
      if(dble(ck0) .ne. 0.d0)then
        rhob=al/dble(ck0)
        dxfrx=fb1**2/rhob/24.d0
        dyfrx=fb1/rhob**2/6.d0
        dyfrax=4.d0*dyfrx/fb1**2
      else
        dxfrx=0.d0
        dyfrx=0.d0
        dyfrax=0.d0
      endif
      if(imag(ck0) .ne. 0.d0)then
        rhob=al/imag(ck0)
        dyfry=fb1**2/rhob/24.d0
        dxfry=fb1/rhob**2/6.d0
        dxfray=4.d0*dxfry/fb1**2
      else
        dyfry=0.d0
        dxfry=0.d0
        dxfray=0.d0
      endif
      do concurrent (i=1:np)
        p=1.d0+g(i)
        x(i)=x(i)+dxfrx*g(i)/p
        y(i)=y(i)-dyfry*g(i)/p
        px(i)=px(i)+(dxfry-dxfray*x(i)**2)/p**2*x(i)
        py(i)=py(i)+(dyfrx-dyfrax*y(i)**2)/p**2*y(i)
        z(i)=z(i)+(dxfrx*px(i)-dyfry*py(i)+
     $       ((.5d0*dyfrx-.25d0*dyfrax*y(i)**2)*y(i)**2
     $       +(.5d0*dxfry-.25d0*dxfray*x(i)**2)*x(i)**2)/p)/p
      enddo
      return
      end
