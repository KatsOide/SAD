      subroutine tmulti(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,ak,bz,phia,psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,theta2,cr1,
     $     eps0,krad,fringe,f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,
     $     vc,w,phirf,dphirf,vnominal,
     $     radius,rtaper,autophi,
     $     kturn,latt,kptbl)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      use sol,only:tsolrot
      use photontable,only:tsetphotongeo
      use mathfun
c      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 , parameter ::nmult=21,itmax=10,ndivmax=1000
      real*8 ,parameter :: conv=3.d-16,oneev=1.d0+1.d-6,
     $     ampmax=0.05d0,alstep=0.05d0,eps00=0.005d0,pmin=1.d-10,
     $     arad=0.01d0
c      parameter (oneev=1.d0+3.83d-12)
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np)
      real*8 sx(np),sy(np),sz(np)
      real*8 al,f1in,f2in,f1out,f2out
      complex*16 ak(0:nmult)
      real*8 bz,phia,psi1,psi2,dx,dy,dz,chi1,chi2,theta,dtheta,eps0
      logical*4 fringe,autophi
      integer*4 mfring
      real*8 fb1,fb2,vc,phirf,dphirf,radius
      integer*8 latt(nlat)
      integer*4 kturn,kptbl(np0,6)
      logical*4 acc,spac1,dofr(0:nmult),krad
      integer*4 i,m,n,ndiv,nmmin,nmmax,ibsi
      real*8 , intent(in)::theta2
      real*8 bxs,bys,bzs,vnominal,b,a,eps,w,wi,v,phis,r,wl,
     $     r1,we,wsn,phic,dphis,offset,offset1,
     $     tlim,akr1,ak1,al1,p,ea,pxf,pyf,sv,wsm,asinh,ws1,wm,
     $     h2,p2,dp2,pr2,dvn,dzn,dp1r,p1r,p1,h1,t,ph,dh,dpr,dp2r,p2r,
     $     alx,dp2p2,dp,dp1,pr1,rtaper,
     $     he,vcorr,v20a,v02a,v1a,v11a,av,dpx,dpy,pe,ah
      real*8 ws(ndivmax)
      complex*16 , intent(in)::cr1
      complex*16 akr(0:nmult),cr,cx,cx1,ak01,b0
      real*8 fact(0:nmult),an(0:nmult+1)
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
     $0.047619047619047619048d0,
     $0.045454545454545454545d0/
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
c      theta2=theta+dtheta+akang(ak(1),al,cr1)
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.true.)
      akr(0)=(ak(0)*cr1+dcmplx(bys*al,bxs*al))*rtaper
      do n=nmult,0,-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmax=n
          go to 1
        endif
      enddo
      if(vc .ne. 0.d0 .or. spac)then
        nmmax=0
      else
        call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,bzs,dble(akr(0)),imag(akr(0)),krad)
        go to 1000
      endif
 1    cr=cr1*rtaper
c     Zero-clear akr(1) for case: nmmax .eq. 0
      akr(1)=0.d0
      do n=1,nmmax
        cr=cr*cr1
        akr(n)=ak(n)*cr
      enddo
c     Im[ark(1)] MIGHT be 0, because akr(1) := ak(1)*cr1*cr1,
c     cr1 := Exp[-theta1], ak(1) = Abs[ak(1)] * Exp[2 theta1]
      akr1=dble(akr(1))
      akr(1)=akr1
      if(eps0 .eq. 0.d0)then
        eps=eps00
      else
        eps=eps00*eps0
      endif
      ndiv=1
      do n=2,nmmax
        ndiv=max(ndiv,
     $int(sqrt(ampmax**(n-1)/6.d0/fact(n-1)/eps*abs(akr(n)*al)))+1)
      enddo
      if(spac)then
        ndiv=max(ndiv,nint(abs(al)/(alstep*eps/eps00)),
     $       nint(eps00/eps*abs(bzs*al)/1.5d0))
      endif
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
        ndiv=max(ndiv,ndivrad(abs(akr(0)),akr1,bz,eps0))
      endif
      ndiv=min(ndivmax,ndiv)
      acc=(trpt .or. rfsw) .and. vc .ne. 0.d0
      if(acc)then
        if(w .eq. 0.d0)then
          wi=0.d0
        else
          wi=1.d0/w
        endif
        v=vc*abs(charge)/amass
        he=h0+vnominal
        pe=h2p(he)
c        pe=sqrt((he-1.d0)*(he+1.d0))
        vcorr=v*(w*(.5d0/p0+.5d0/pe))**2/4.d0
        v20a=vcorr
        v02a=vcorr
        v1a=0.d0
        v11a=0.d0
        r=abs(vnominal*(1.d0/h0+1.d0/he))
        ndiv=min(ndivmax,max(ndiv,
     $       1+int(min(abs(w*al),
     $       sqrt(r**2/3.d0/(eps/eps00*1.d-3))))))
        if(abs(r) .gt. 1.d-6)then
          wl=asinh(r*(2.d0+r)/2.d0/(1.d0+r))
          ndiv=min(ndivmax,max(ndiv,
     $         nint(wl*10.d0*eps00/eps)))
          r1=wl/ndiv/2.d0
          ws1=2.d0*exp(r1)*sinh(r1)
          we=1.d0+ws1
          ws(1)=ws1/r
          do i=2,ndiv
            ws(i)=ws(i-1)*we
          enddo
        else
          wsn=1.d0/ndiv
          do i=1,ndiv
            ws(i)=wsn
          enddo
        endif            
        phic=(phirf+dphirf)*charge
        if(trpt)then
c          vnominal=v*sin(-phirf*charge)
          phis=0.d0
        else
c          vnominal=0.d0
          if(autophi)then
            phis=phic
          else
            phis=w*trf0
          endif
        endif
        dphis=phis-phic
        if(rad .or. trpt .or. autophi)then
          offset=sin(dphis)
          offset1=0.d0
        else
          offset=sin(dphis)-sin(phis)
          offset1=sin(phis)
        endif
        tlim=1.d4
      else
c     begin initialize for preventing compiler warning
        w=0.d0
        wi=0.d0
        v=0.d0
        v20a=0.d0
        v02a=0.d0
        dphis=0.d0
        offset=0.d0
        offset1=0.d0
c     end   initialize for preventing compiler warning
c        vnominal=0.d0
        wsn=1.d0/ndiv
        do i=1,ndiv
          ws(i)=wsn
        enddo
      endif
      ak1=akr1*ws(1)*.5d0
      al1=al*ws(1)*.5d0
      ak01=akr(0)*ws(1)*.5d0
      if(al .ne. 0.d0)then
        nmmin=2
      else
        nmmin=1
      endif
      if(al .ne. 0.d0)then
        if(fringe)then
          if(mfring .ne. 2 .and. acc)then
            call tcavfrin(np,x,px,y,py,z,g,dv,al,v,w,p0,h0,
     $           dphis,dvfs,offset)
          endif
          do n=0,nmmax
            if(akr(n) .ne. (0.d0,0.d0))then
              if(n .le. 2)then
                dofr(n)=.true.
              else
                dofr(n)=abs(akr(n))/al*ampmax**(n+1)*an(n+1) .gt. eps
              endif
              if(dofr(n) .and. mfring .ne. 2)then
                call ttfrins(np,x,px,y,py,z,g,n*2+2,akr(n),
     $               al,bzs)
              endif
            else
              dofr(n)=.false.
            endif
          enddo
        endif
        if(mfring .eq. 1 .or. mfring .eq. 3)then
          if(akr(0) .ne. (0.d0,0.d0) .and. fb1 .ne. 0.d0)then
            call tblfri(np,x,px,y,py,z,g,al,akr(0),fb1)
          endif
          if(f1in .ne. 0.d0 .or. f2in .ne. 0.d0)then
            do i=1,np
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
      spac1 = spac .and. radius .ne. 0.d0
      sv=0.d0
      ibsi=1
c      if(l_track .gt. 7300 .and. l_track .lt. 7310)then
c        write(*,'(a,2i5,1p4g15.7)')'tmulti-2 ',
c     $       l_track,ndiv,al1,ak1,x(1),px(1)
c      endif
      do m=1,ndiv
        if(nmmin .eq. 2)then
         call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al1,ak1,bzs,dble(ak01),imag(ak01),ibsi,eps0)
          if(krad)then
            if(m .eq. 1)then
              do i=1,np
                cx1=dcmplx(x(i),y(i))
                cx=0.d0
                do n=nmmax,2,-1
                  cx=(cx+akr(n))*cx1*an(n+1)
                enddo
                cx=.5d0*cx*cx1**2
                bsi(i)=bsi(i)+imag(cx)/al
              enddo
            endif
            if(photons)then
              if(m .eq. 1)then
                call tsetphotongeo(al1,0.d0,theta2,.true.)
              else
                call tsetphotongeo(al1,0.d0,0.d0,.false.)
              endif
            endif
            call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
          endif
          ibsi=0
        endif
        wsm=ws(m)
        if(m .eq. ndiv)then
          wm=wsm*.5d0
        else
          wm=(wsm+ws(m+1))*.5d0
        endif
        al1=al*wm
        ak1=akr1*wm
        ak01=akr(0)*wm
        if(nmmin .eq. 2)then
          do i=1,np
            cx1=dcmplx(x(i),y(i))
            cx=(0.d0,0.d0)
            do n=nmmax,2,-1
              cx=(cx+(akr(n)*wsm))*cx1*an(n)
            enddo
            cx=cx*cx1/(1.d0+g(i))
            px(i)=px(i)-dble(cx)
            py(i)=py(i)+imag(cx)
          enddo
        else
          do i=1,np
            cx1=dcmplx(x(i),y(i))
            cx=(0.d0,0.d0)
            do n=nmmax,1,-1
              cx=(cx+(akr(n)*wsm))*cx1*an(n)
            enddo
            cx=(cx+akr(0)*wsm)/(1.d0+g(i))
            px(i)=px(i)-dble(cx)
            py(i)=py(i)+imag(cx)
          enddo
        endif
        if(acc)then
          if(vnominal .ne. 0.d0)then
            sv=sv+vnominal*wsm
            h2=h0+sv
            p2=h2p(h2)
c            p2=sqrt((h2-1.d0)*(h2+1.d0))
            dp2=sv*(h2+h0)/(p2+p0)/p0
            pr2=1.d0+dp2
            dvn=-dp2*(1.d0+pr2)/h2/(h2+pr2*h0)
            dzn=-dvn*al*wsm*.5d0
          else
            dzn=0.d0
          endif
          do i=1,np
            dp1r=g(i)
            p1r=1.d0+dp1r
            p1=p0*p1r
            if(dv(i) .gt. 0.1d0)then
              h1=p2h(p1)
c              h1=sqrt(p1**2+1.d0)
            else
              h1=p1r*h0/(1.d0-dv(i)+dvfs)
            endif
            z(i)=z(i)-dzn
            t=min(tlim,max(-tlim,-z(i)*h1/p1))
            ph=.5d0*w*t
            dh=max(oneev-h1,
     $           (v+(v1a+v20a*x(i)+v11a*y(i))*x(i)+v02a*y(i)**2)
     $           *(-2.d0*sin(ph)*cos(ph-dphis)+offset)*wsm)
            h2=h1+dh
            ah=max(dh*(h1+h2)/p1**2,-1.d0+pmin)
            dpr=sqrt1(ah)
c            dpr=ah/(1.d0+sqrt(1.d0+ah))
            dp2r=max(dp1r+p1r*dpr,pmin-1.d0)
            p2r=1.d0+dp2r
            g(i)=dp2r
            dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs
            z(i)=-t*p2r*p0/h2-dzn
            av=-(cos(2.d0*ph-dphis)*wi-offset1*t)*wsm/p0
            dpx=(v1a+2.d0*v20a*x(i)+v11a*y(i))*av
            dpy=(v11a*x(i)+2.d0*v02a*y(i))*av
            px(i)=(px(i)*p1r+dpx)/p2r
            py(i)=(py(i)*p1r+dpy)/p2r
          enddo
        endif
        if(spac1)then
          call spkick(np,x,px,y,py,z,g,dv,sx,sy,sz,al*wsm,radius,alx,
     $          kturn,l_track,latt,kptbl)
        endif
      enddo
      if(nmmin .eq. 2)then
        call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al1,ak1,bzs,dble(ak01),imag(ak01),2,eps0)
      endif
c      if(l_track .gt. 7300 .and. l_track .lt. 7310)then
c        write(*,'(a,2i5,1p4g15.7)')'tmulti-2 ',
c     $       l_track,ndiv,al1,ak1,x(1),px(1)
c      endif
      if(spac1)then
        call tapert(l_track,latt,x,px,y,py,z,g,dv,sx,sy,sz,
     $       kptbl,np,kturn,
     $       radius,radius,
     $       0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      if(al .ne. 0.d0)then
        if(mfring .eq. 2 .or. mfring .eq. 3)then
          if(f1out .ne. 0.d0 .or. f2out .ne. 0.d0)then
            do i=1,np
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
          if(acc)then
            call tcavfrin(np,x,px,y,py,z,g,dv,al,-v,w,p0,h0,
     $         dphis,dvfs,offset)
          endif
        endif
        if(krad)then
          do i=1,np
            cx1=dcmplx(x(i),y(i))
            cx=0.d0
            do n=nmmax,2,-1
              cx=(cx+akr(n))*cx1*an(n+1)
            enddo
            cx=.5d0*cx*cx1**2
            bsi(i)=bsi(i)-imag(cx)/al
          enddo
          if(photons)then
            call tsetphotongeo(al1,0.d0,0.d0,.false.)
          endif
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al1,0.d0)
        endif
      endif
 1000 continue
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,theta2,bxs,bys,bzs,.false.)
      if(vnominal .ne. 0.d0)then
        h2=h0+vnominal
        p2=h2p(h2)
c        p2=sqrt((h2-1.d0)*(h2+1.d0))
        dp2p2=vnominal*(h2+h0)/(p2+p0)/p2
        do i=1,np
          dp=max(g(i),pmin-1.d0)
          dp1=dp*p0/p2-dp2p2
          pr1=1.d0+dp1
          g(i)=dp1
          h1=p2h(p2*pr1)
c          h1=sqrt(1.d0+(p2*pr1)**2)
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
      integer*4 np,i
      complex*16 ck0
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),
     $     dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray,
     $     p,al,fb1,rhob
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
      do i=1,np
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
