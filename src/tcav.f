      subroutine tcav(np,x,px,y,py,z,g,dv,sx,sy,sz,al,vc,
     $     w,phi,dphi,vnominal,
     $     lwl,wakel,lwt,waket,
     $     dx,dy,theta,v1,v20,v11,v02,
     $     fringe,mfring,autophi)
      use ffs_flag
      use tmacro
      use mathfun
      implicit none
      integer*4, parameter ::ndivmax=1000
      real*8, parameter::eps=1.d-3,oneev=1.d0+3.83d-12
      integer*4 ,intent(in):: np,mfring,lwl,lwt
      integer*4 ndiv,i,n,itab(np),izs(np)
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     g(np),dv(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: wakel(2,lwl),waket(2,lwt),
     $     al,vc,w,phi,dphi,vnominal,dx,dy,theta,v1,v20,v11,v02
      real*8 he,wi,cost,sint,v,v1a,ws(ndivmax),
     $     v20a,v11a,v02a,phis,r,wl,r1,ws1,we,wsn,phic,
     $     dphis,offset,offset1,tlim,sv,a,dpz,al1,h2,p2,
     $     dp2,pr2,dvn,dzn,dp1r,p1r,p1,h1,dp2r,p2r,av,dpx,fw,
     $     alx,dpepe,xi,pxi,fw0,
     $     asinh,t,ph,dh,dpr,dpy,dp,dp1,pr1,pe,vcorr
      logical*4 fringe,autophi
      if(w .eq. 0.d0)then
        wi=0.d0
      else
        wi=1.d0/w
      endif
      include 'inc/TENT.inc'
      phic=(phi+dphi)*charge
      v=vc/amass*abs(charge)
      v1a=v1/amass*abs(charge)
      v20a=.5d0*v20/amass*abs(charge)
      v11a=v11/amass*abs(charge)
      v02a=.5d0*v02/amass*abs(charge)
      phis=merge(0.d0,merge(phic,w*trf0,autophi),trpt)
      he=h0+vnominal
      pe=h2p(he)
c      pe=he*sqrt(1.d0-1.d0/he**2)
c      pe=sqrt((he-1.d0)*(he+1.d0))
      dpepe=vnominal*(he+h0)/(pe+p0)/pe
      vcorr=v*(w*(.5d0/p0+.5d0/pe))**2/4.d0
      v20a=v20a+vcorr
      v02a=v02a+vcorr
      r=abs(vnominal*(1.d0/h0+1.d0/he))
      ndiv=1+int(min(abs(w*al),
     $     sqrt(r**2/3.d0/eps)))
      if(r .gt. 1.d-6)then
        wl=asinh(r*(2.d0+r)/2.d0/(1.d0+r))
        ndiv=min(ndivmax,max(ndiv,
     $       nint(wl*10.d0)))
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
      dphis=phis-phic
      if(rad .or. trpt .or. autophi)then
        offset=sin(dphis)
        offset1=0.d0
      else
        offset=sin(dphis)-sin(phis)
        offset1=sin(phis)
      endif
      if(rfsw)then
        tlim=1.d4
      else
        tlim=0.d0
      endif
      if(al .ne. 0.d0 .and. fringe .and.
     $     mfring .ge. 0 .and. mfring .ne. 2)then
        call tcavfrin(np,x,px,y,py,z,g,dv,al,v,w,p0,h0,
     $     dphis,dvfs,offset)
      endif
      fw0=(abs(charge)*e*pbunch*anbunch/amass)/np0
      sv=0.d0
      do 110 n=1,ndiv
        wsn=ws(n)
        if(al .ne. 0.d0)then
          if(n .eq. 1)then
            alx=al*wsn*.5d0
          else
            alx=al*(wsn+ws(n-1))*.5d0
          endif
          do 10 i=1,np
            a=px(i)**2+py(i)**2
            dpz=sqrt1(-a)
c            dpz=-a/(1.d0+sqrt(1.d0-a))
            al1=alx/(1.d0+dpz)
            x(i)=x(i)+px(i)*al1
            y(i)=y(i)+py(i)*al1
            z(i)=z(i)+dpz  *al1-dv(i)*alx
10        continue
          if(vnominal .ne. 0.d0)then
            sv=sv+vnominal*wsn
            h2=h0+sv
            p2=h2p(h2)
c            p2=h2*sqrt(1.d0-1.d0/h2**2)
c            p2=sqrt((h2-1.d0)*(h2+1.d0))
            dp2=sv*(h2+h0)/(p2+p0)/p0
            pr2=1.d0+dp2
            dvn=-dp2*(1.d0+pr2)/h2/(h2+pr2*h0)
            dzn=-dvn*al*wsn*.5d0
          else
            dzn=0.d0
          endif
        else
          dzn=0.d0
        endif
        do 20 i=1,np
          dp1r=g(i)
          p1r=1.d0+dp1r
          p1=p0*p1r
          if(dv(i) .gt. 0.1d0)then
            h1=p2h(p1)
c            h1=p1*sqrt(1.d0+1.d0/p1**2)
c            h1=sqrt(1.d0+p1**2)
          else
            h1=p1r*h0/(1.d0-dv(i)+dvfs)
          endif
          z(i)=z(i)-dzn
          t=min(tlim,max(-tlim,-z(i)*h1/p1))
          ph=.5d0*w*t
          dh=max(oneev-h1,
     $         (v+(v1a+v20a*x(i)+v11a*y(i))*x(i)+v02a*y(i)**2)
     $         *(-2.d0*sin(ph)*cos(ph-dphis)+offset)*wsn)
c          if(autophi)then
c            write(*,'(a,1p4g15.7)')'tcav ',dh,ph,dphis,offset
c          endif
          h2=h1+dh
          a=max(dh*(h1+h2)/p1**2,-1.d0)
          dpr=sqrt1(a)
c          dpr=a/(1.d0+sqrt(1.d0+a))
          dp2r=dp1r+p1r*dpr
          p2r=1.d0+dp2r
          g(i)=dp2r
          dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs
          av=-(cos(2.d0*ph-dphis)*wi-offset1*t)*wsn/p0
          dpx=(v1a+2.d0*v20a*x(i)+v11a*y(i))*av
          dpy=(v11a*x(i)+2.d0*v02a*y(i))*av
          px(i)=(px(i)*p1r+dpx)/p2r
          py(i)=(py(i)*p1r+dpy)/p2r
          z(i)=-t*p2r*p0/h2-dzn
20      continue
        if(lwake .or. twake)then
          fw=fw0*wsn
          call txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         0.d0,0.d0,0.d0,
     $         int(anbunch),
     $         fw,lwl,wakel,lwt,waket,p0,h0,itab,izs,n .eq. 1)
        endif
110   continue
c      write(*,'(a,1p5g15.7)')'tcav-2 ',y(1),py(1),z(1),g(1),dv(1)
      if(al .ne. 0.d0)then
        alx=al*ws(ndiv)*.5d0
        do 30 i=1,np
          a=px(i)**2+py(i)**2
          dpz=sqrt1(-a)
c          dpz=-a/(1.d0+sqrt(1.d0-a))
          al1=alx/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*alx
30      continue
        if(fringe .and. mfring .ge. 0 .and. mfring .ne. 1)then
          call tcavfrin(np,x,px,y,py,z,g,dv,al,-v,w,p0,h0,
     $         dphis,dvfs,offset)
        endif
      endif
      if(vnominal .ne. 0.d0)then
c        write(*,*)'tcav ',vnominal,p0,pe
        do 210 i=1,np
          dp=g(i)
          dp1=dp*p0/pe-dpepe
          pr1=1.d0+dp1
          g(i)=dp1
          h1=p2h(pe*pr1)
c          h1=pe*pr1*sqrt(1.d0+1.d0/(pe*pr1)**2)
c          h1=sqrt(1.d0+(pe*pr1)**2)
          dv(i)=-dp1*(1.d0+pr1)/h1/(h1+pr1*he)+dvfs
210     continue
        pgev=pe*amass
        call tphyzp
      endif
      include 'inc/TEXIT.inc'
      return
      end

      subroutine tcavfrin(np,x,px,y,py,z,g,dv,al,v,w,p0,h0,
     $     dphis,dvfs,offset)
      use mathfun
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np)
      real*8 al,v,p0,dphis,offset,h0,w,s0,dvfs,
     $     vf,dp1r,p1r,p1,h1,t,ph,dpt,dh,h2,a,dpr,dp2r,p2r,pmin
      parameter (pmin=1.e-10)
      vf=v/al/p0*.5d0
      s0=sin(dphis)-offset
      do i=1,np
        dp1r=g(i)
        p1r=1.d0+dp1r
        p1=p0*p1r
        if(dv(i) .gt. 0.1d0)then
          h1=p2h(p1)
c          h1=p1*sqrt(1.d0+1.d0/p1**2)
c          h1=sqrt(1.d0+p1**2)
        else
          h1=p1r*h0/(1.d0-dv(i)+dvfs)
        endif
        t=-z(i)*h1/p1
        ph=w*t
        dpt=vf*(sin(ph-dphis)+s0)/p1r
        px(i)=px(i)+x(i)*dpt
        py(i)=py(i)+y(i)*dpt
        dh=-w*p0*vf*(x(i)**2+y(i)**2)*.5d0*cos(ph-dphis)
        h2=h1+dh
        a=max(dh*(h1+h2)/p1**2,-1.d0)
        dpr=sqrt1(a)
c        dpr=a/(1.d0+sqrt(1.d0+a))
        dp2r=max(dp1r+p1r*dpr,pmin-1.d0)
        p2r=1.d0+dp2r
        g(i)=dp2r
        z(i)=z(i)*p2r/p1r*h1/h2
        dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs
      enddo
      return
      end
