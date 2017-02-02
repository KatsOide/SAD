      subroutine tmulti(np,x,px,y,py,z,g,dv,pz,
     $     al,ak,bz,phia,
     $     psi1,psi2,
     $     dx,dy,dz,chi1,chi2,theta,dtheta,
     $     eps0,enarad,fringe,f1,f2,mfring,fb1,fb2,
     $     vc,w,phirf,dphirf,radius,rtaper,autophi,
     $     kturn,l,latt,kptbl)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 nmult,itmax,ndivmax
      real*8 conv,ampmax,alstep,eps00,oneev,pmin
      parameter (nmult=21,itmax=10,ndivmax=1000,conv=3.d-16,
     $     ampmax=0.05d0,alstep=0.05d0,eps00=0.005d0,pmin=1.d-10)
c      parameter (oneev=1.d0+3.83d-12)
      parameter (oneev=1.d0+1.d-6)
      integer*4 np
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 al,f1r,f2r
      complex*16 ak(0:nmult)
      real*8 bz,phia,psi1,psi2,dx,dy,dz,chi1,chi2,theta,dtheta,eps0
      logical*4 enarad,fringe,autophi
      real*8 f1,f2
      integer*4 mfring
      real*8 fb1,fb2,vc,phirf,dphirf,radius
      integer*8 latt(nlat)
      integer*4 kturn,l,kptbl(np0,6)
      logical acc,spac1,dofr(0:nmult)
      integer*4 i,j,m,n,ndiv,nmmin,nmmax
      real*8 pz0,s0,bxs,bys,bzs,
     $     vnominal,theta1,theta2,
     $     cchi1,schi1,b,b1,
     $     phix,phiy,phiz,dphizsq,pr,ds1,ds2,pz1,
     $     dcchi1,cchi2,schi2,bzp,alb,s,dpz0,
     $     dpl,pl,plx,ply,plz,ptx,pty,ptz,pbx,pby,pbz,phi,
     $     sinphi,dcosphi,xsinphi,xsin,dphi,phi0,pz2,a,
     $     fx,fy,cost,sint,x0,px0,bxs0,eps,w,wi,v,phis,r,wl,
     $     dcchi2,radlvl,r1,we,wsn,phic,dphis,offset,offset1,
     $     tlim,akr1,ak1,al1,p,ea,pxf,pyf,sv,wsm,asinh,ws1,wm,
     $     h2,p2,dp2,pr2,dvn,dzn,dp1r,p1r,p1,h1,t,ph,dh,dpr,dp2r,p2r,
     $     alx,pyi,pxi,y1,x1,a24,a12,a22,a14,dp2p2,dp,dp1,pr1,
     $     he,vcorr,v20a,v02a,v1a,v11a,av,dpx,dpy,pe,ah,z00,
     $     rtaper
      real*8 ws(ndivmax)
      complex*16 akr(0:nmult),cr,cr1,cx,cx1,ak01,b0
      real*8 fact(0:nmult),an(0:nmult)
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
        call tmulta(
     $       np,x,px,y,py,z,g,dv,pz,
     $       l,al,ak,phia,
     $       psi1,psi2,bz,
     1       dx,dy,theta,dtheta,
     $       eps0,enarad,fb1,fb2,mfring,fringe)
        return
      endif
      radlvl=1.d0
      b0=0.d0
      z00=z(1)
      if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or. chi2 .ne. 0.d0)then
c     begin initialize for preventing compiler warning
        fx=0.d0
        fy=0.d0
c     end   initialize for preventing compiler warning
        s0=-al*.5d0
        cchi1=cos(chi1)
        schi1=sin(chi1)
        if(cchi1 .ge. 0.d0)then
          dcchi1=schi1**2/(1.d0+cchi1)
        else
          dcchi1=(1.d0-cchi1)
        endif
        cchi2=cos(chi2)
        schi2=sin(chi2)
        if(cchi1 .ge. 0.d0)then
          dcchi2=schi2**2/(1.d0+cchi2)
        else
          dcchi2=(1.d0-cchi2)
        endif
        if(bz .ne. 0.d0)then
          b=abs(bz)
          phix=(-schi1)         *sign(1.d0,bz)
          phiy=( cchi1)*(-schi2)*sign(1.d0,bz)
          phiz=( cchi1)*( cchi2)*sign(1.d0,bz)
          dphizsq=phix**2+phiy**2
          bxs=phix*b
          bys=phiy*b
          bzs=phiz*b
          do i=1,np
            pr=(1.d0+g(i))
            px(i)=px(i)+bz*y(i)/pr*.5d0
            py(i)=py(i)-bz*x(i)/pr*.5d0
            x(i)=x(i)-dx
            y(i)=y(i)-dy
            ds1  = schi1*x(i)-dcchi1*s0-cchi1*dz
            x(i) = cchi1*x(i)-schi1*(s0-dz)
            ds2  =-(schi2*y(i)+cchi2*ds1-dcchi2*s0)
            y(i) = cchi2*y(i)-schi2*(ds1+s0)
            pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c            pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
            pz1  = schi1*px(i)+cchi1*pz0
            px(i)= cchi1*px(i)-schi1*pz0
            py(i)= cchi2*py(i)-schi2*pz1
            bzp=bzs/pr
            alb=pr/b
            s=px(i)**2+py(i)**2
            dpz0=sqrt1(-s)
c            dpz0=-s/(1.d0+sqrt(1.d0-s))
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
            if(ds2 .ne. 0.d0)then
              phi=ds2/alb/pz0
              do j=1,itmax
                sinphi=sin(phi)
                dcosphi=2.d0*sin(.5d0*phi)**2
                xsinphi=xsin(phi)
                s=(plz*xsinphi+pz0*sinphi+pbz*dcosphi)*alb
                dphi=(ds2-s)/alb/(pz0-ptz*dcosphi+pbz*sinphi)
                phi0=phi
                phi=phi+dphi
                if(phi0 .eq. phi .or. abs(dphi/phi) .lt. conv)then
                  go to 100
                endif
              enddo
c              write(*,*)'tmulti convergence error',phi,dphi
            else
              phi=0.d0
              sinphi=0.d0
              dcosphi=0.d0
            endif
 100        x(i)=x(i)+(plx*phi+ptx*sinphi+pbx*dcosphi)*alb
            y(i)=y(i)+(ply*phi+pty*sinphi+pby*dcosphi)*alb
            z(i)=z(i)-phi*alb
            px(i)=px(i)-ptx*dcosphi+pbx*sinphi-bzp*y(i)*.5d0
            py(i)=py(i)-pty*dcosphi+pby*sinphi+bzp*x(i)*.5d0
          enddo        
        else
          bxs=0.d0
          bys=0.d0
          bzs=0.d0
          do i=1,np
            pr=(1.d0+g(i))
            x(i)=x(i)-dx
            y(i)=y(i)-dy
            pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c            pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
            pz1  = schi1*px(i)+cchi1*pz0
            px(i)= cchi1*px(i)-schi1*pz0
            py(i)= cchi2*py(i)-schi2*pz1
            ds1  = schi1*x(i)-dcchi1*s0-cchi1*dz
            x(i) = cchi1*x(i)-schi1*(s0-dz)
            ds2  = schi2*y(i)+cchi2*ds1-dcchi2*s0
            y(i) = cchi2*y(i)-schi2*(ds1+s0)
            pz2=1.d0+sqrt1(-px(i)**2-py(i)**2)
c            pz2=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
            a=ds2/pz2
            x(i) =x(i)-a*px(i)
            y(i) =y(i)-a*py(i)
            z(i) =z(i)+a
          enddo
        endif
      else
c     begin initialize for preventing compiler warning: is it necessary?
        cchi1=0.d0
        cchi2=0.d0
        schi1=0.d0
        schi2=0.d0
        dcchi1=0.d0
        dcchi2=0.d0
        s0=0.d0
c     end   initialize for preventing compiler warning
        bzs=bz
        bxs=0.d0
        bys=0.d0
        fx= bzs*dy*.5d0
        fy=-bzs*dx*.5d0
        do i=1,np
          pr=(1.d0+g(i))
          x(i)=x(i)-dx
          y(i)=y(i)-dy
          px(i)=px(i)+fx/pr
          py(i)=py(i)+fy/pr
        enddo
      endif
      if(imag(ak(1)) .eq. 0.d0)then
        theta1=0.d0
      else
        theta1=atan2(imag(ak(1)),dble(ak(1)))*.5d0
      endif
      theta2=theta+theta1
      if(theta2 .ne. 0.d0)then
        cost=cos(theta2)
        sint=sin(theta2)
        do i=1,np
          x0=x(i)
          x(i)=cost*x0-sint*y(i)
          y(i)=sint*x0+cost*y(i)
          px0=px(i)
          px(i)=cost*px0-sint*py(i)
          py(i)=sint*px0+cost*py(i)
        enddo
        bxs0=bxs
        bxs=cost*bxs0-sint*bys
        bys=sint*bxs0+cost*bys
      else
c     begin initialize for preventing compiler warning
        cost=1.d0
        sint=0.d0
c     end   initialize for preventing compiler warning
      endif
c      write(*,'(2a,1p7g12.5)')'tmulti ',pname(latt(1,l))(1:8),rtaper,
c     $     x(np),px(np),y(np),py(np),z(np),g(np)
      cr1=dcmplx(cos(theta1),-sin(theta1))
      akr(0)=(ak(0)*cr1+dcmplx(bys*al,bxs*al))*rtaper
      do n=nmult,0,-1
        if(ak(n) .ne. (0.d0,0.d0))then
          nmmax=n
          go to 1
        endif
      enddo
      if(vc .ne. 0.d0 .or. spac)then
        nmmax=0
        go to 1
      else
        call tdrift(np,x,px,y,py,z,g,dv,pz,al,
     $       bzs,dble(akr(0)),imag(akr(0)))
      endif
      vnominal=0.d0
      go to 1000
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
     $int(sqrt(ampmax**(n-1)/6.d0/fact(n-1)/eps*abs(akr(n))*al))+1)
      enddo
      ndiv=min(ndivmax,ndiv)
      if(spac)then
        ndiv=max(ndiv,nint(al/(alstep*eps/eps00)),
     $       nint(eps00/eps*abs(bzs*al)/1.5d0))
      endif
      acc=(trpt .or. rfsw) .and. vc .ne. 0.d0
      if(acc)then
        if(w .eq. 0.d0)then
          wi=0.d0
        else
          wi=1.d0/w
        endif
        v=vc*abs(charge)/amass
        if(trpt)then
          vnominal=v*sin(-phirf*charge)
          phis=0.d0
        else
          vnominal=0.d0
          if(autophi)then
            phis=phic
          else
            phis=w*trf0
          endif
        endif
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
        phic=(phirf+dphirf)*charge-vcphic
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
        vnominal=0.d0
        wsn=1.d0/ndiv
        do i=1,ndiv
          ws(i)=wsn
        enddo
      endif
      ak1=akr1*ws(1)*.5d0
      al1=al*ws(1)*.5d0
      ak01=akr(0)*ws(1)*.5d0
      if(al .gt. 0.d0)then
        nmmin=2
      else
        nmmin=1
      endif
      if(al .gt. 0.d0)then
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
          if(f1 .ne. 0.d0 .or. f2 .ne. 0.d0)then
            do i=1,np
              p=(1.d0+g(i))
              a=f1/p
              ea=exp(a)
              b=f2/p
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
        if(rad .and. enarad .and.
     $       (akr1 .ne. 0.d0 .or. akr(0) .ne. (0.d0,0.d0)))then
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
          radlvl=0.d0
          b1=brhoz*akr1/al
          b0=brhoz*akr(0)/al
          call trad(np,x,px,y,py,g,dv,dble(b0),-imag(b0),
     1         b1,0.d0,0.d0,.5d0*al,
     $         f1r,f2r,0.d0,al,1.d0)
        endif
      endif
      spac1 = spac .and. radius .ne. 0.d0
      sv=0.d0
      do m=1,ndiv
        if(nmmin .eq. 2)then
          call tsolqu(np,x,px,y,py,z,g,dv,pz,al1,ak1,
     $         bzs,dble(ak01),imag(ak01),eps0)
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
          call spkick(np,x,px,y,py,z,g,dv,pz,al*wsm,radius,alx,
     $          kturn,l,latt,kptbl)
        endif
      enddo
      if(nmmin .eq. 2)then
        call tsolqu(np,x,px,y,py,z,g,dv,pz,al1,
     $       ak1,bzs,dble(ak01),imag(ak01),eps0)
      endif
      if(spac1)then
c        call spapert(np,x,px,y,py,z,g,dv,radius,kptbl)
        call tapert(l,latt,x,px,y,py,z,g,dv,pz,kptbl,np,kturn,
     $       radius,radius,
     $       0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      if(al .gt. 0.d0)then
        if(radlvl .eq. 0.d0)then
          call trad(np,x,px,y,py,g,dv,dble(b0),-imag(b0),
     1         b1,0.d0,0.d0,.5d0*al,
     $         f1r,f1r,al,al,-1.d0)
        endif
        if(mfring .eq. 2 .or. mfring .eq. 3)then
          if(f1 .ne. 0.d0 .or. f2 .ne. 0.d0)then
            do i=1,np
              p=(1.d0+g(i))
              a=-f1/p
              ea=exp(a)
              b= f2/p
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
      endif
 1000 if(theta2 .ne. 0.d0)then
        do i=1,np
          x0=x(i)
          x(i)= cost*x0+sint*y(i)
          y(i)=-sint*x0+cost*y(i)
          px0=px(i)
          px(i)= cost*px0+sint*py(i)
          py(i)=-sint*px0+cost*py(i)
        enddo
      endif
      if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or. chi2 .ne. 0.d0)then
        if(bz .ne. 0.d0)then
          do i=1,np
            pr=(1.d0+g(i))
            px(i)=px(i)+bzs/pr*y(i)*.5d0
            py(i)=py(i)-bzs/pr*x(i)*.5d0
            pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c            pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
            pz1  =-schi2*py(i)+cchi2*pz0
            pyi  = cchi2*py(i)+schi2*pz0
            pz2  =-schi1*px(i)+cchi1*pz1
            pxi  = cchi1*px(i)+schi1*pz1
            y1   = cchi2*y(i)-schi2*s0
            ds1  =-schi2*y(i)+dcchi2*s0
            ds2  =-schi1*x(i)+cchi1*ds1+dcchi1*s0+dz
            x1   = cchi1*x(i)+schi1*(ds1-s0)
            bzp=bz/pr
            phi=-bzp*ds2/pz2
            a24=sin(phi)
            a12=a24/bzp
            a22=cos(phi)
            if(a22 .ge. 0.d0)then
              a14=a24**2/(1.d0+a22)/bzp
            else
              a14=(1.d0-a22)/bzp
            endif
            x(i) =x1 +a12*pxi+a14*pyi+dx
            y(i) =y1 -a14*pxi+a12*pyi+dy
            px(i)=    a22*pxi+a24*pyi-bzp*y(i)*.5d0
            py(i)=   -a24*pxi+a22*pyi+bzp*x(i)*.5d0
            z(i) =z(i)+ds2/pz2
          enddo
        else
          do i=1,np
c            pr=(1.d0+g(i))**2
            pr=(1.d0+g(i))
            pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c            pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
            pz1  =-schi2*py(i)+cchi2*pz0
            pyi  = cchi2*py(i)+schi2*pz0
            pz2  =-schi1*px(i)+cchi1*pz1
            pxi  = cchi1*px(i)+schi1*pz1
            y1   = cchi2*y(i)-schi2*s0
            ds1  =-schi2*y(i)+dcchi2*s0
            ds2  =-schi1*x(i)+cchi1*ds1+dcchi1*s0+dz
            x1   = cchi1*x(i)+schi1*(ds1-s0)
            a=ds2/pz2
            px(i)=pxi
            py(i)=pyi
            x(i) =x1-a*pxi+dx
            y(i) =y1-a*pyi+dy
            z(i) =z(i)+a
          enddo
        endif
      else
        do i=1,np
          pr=(1.d0+g(i))
          px(i)=px(i)-fx/pr
          py(i)=py(i)-fy/pr
          x(i)=x(i)+dx
          y(i)=y(i)+dy
        enddo
      endif
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
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),
     $     dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray,
     $     p,ck0(2),al,fb1,rhob
      if(ck0(1) .ne. 0.d0)then
        rhob=al/ck0(1)
        dxfrx=fb1**2/rhob/24.d0
        dyfrx=fb1/rhob**2/6.d0
        dyfrax=4.d0*dyfrx/fb1**2
      else
        dxfrx=0.d0
        dyfrx=0.d0
        dyfrax=0.d0
      endif
      if(ck0(2) .ne. 0.d0)then
        rhob=al/ck0(2)
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
