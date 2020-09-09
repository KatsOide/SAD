      module photontable
      implicit none
      type photonp
        sequence
        real*8 al,phi,theta,rho,cost,sint,chi,geo1(3,4)        
        integer*4 l
      end type
      type photoncvt
        sequence
        real*8 dx,dy,theta,dtheta,phi0,al,cost,sint,fr0
        integer*4 l
      end type
      integer*4 ntable,ltable
      parameter (ntable=256,ltable=100000)
      integer*8 kphtable(ntable)
      integer*4 nt,itp,ilp,lt
      type (photonp) pp
      type (photoncvt) pcvt

      contains

      subroutine tphotoninit()
      use tfstk
      implicit none
      nt=ntable
      lt=ltable
      itp=0
      ilp=0
      kphtable(1)=0
      return
      end subroutine

      subroutine tsetpgfrac(l,al)
      use maccbk, only:idtype
      use sad_main, only:sad_comp
      use ffs
      use ffs_pointer,only:compelc
      use geolib
      implicit none
      integer*4 ,intent(in):: l
      real*8 ,intent(in):: al
      real*8 ,parameter :: almin=1.d-12
      type (sad_comp) ,pointer ::cmp
      integer*4 k,irtc,it
      real*8 fr,al0
      fr=0.d0
      call compelc(l,cmp)
      it=idtype(cmp%id)
      k=kytbl(kwL,it)
      if(k .ne. 0)then
        al0=cmp%value(k)
        if(al0 .ne. 0.d0)then
          fr=max(al,almin)/al0
        endif
      endif
      pp%geo1=tfgeofrac(l,fr,irtc)
      return
      end

      subroutine tsetphotongeo(al0,phi0,theta0,ini)
      use tmacro, only:l_track
      use ffs_pointer, only:geo
      implicit none
      real*8, intent(in):: al0,phi0,theta0
      logical*4 , intent(in)::ini
      real*8 gv(3,4),cp0,sp0,cost,sint,r1,r2
      associate(l=>pp%l,al=>pp%al,phi=>pp%phi,theta=>pp%theta,
     $     x1=>pp%geo1(1,1),x2=>pp%geo1(2,1),x3=>pp%geo1(3,1),
     $     y1=>pp%geo1(1,2),y2=>pp%geo1(2,2),y3=>pp%geo1(3,2),
     $     z1=>pp%geo1(1,3),z2=>pp%geo1(2,3),z3=>pp%geo1(3,3),
     $     gx0=>pp%geo1(1,4),gy0=>pp%geo1(2,4),gz0=>pp%geo1(3,4),
     $     rho=>pp%rho,chi=>pp%chi,geo1=>pp%geo1)
      l=l_track
      gv=merge(geo(:,:,l),geo1,ini)
      al=al0
      theta=theta0
      phi=phi0
      cost=cos(theta)
      sint=sin(theta)
      x1= cost*gv(1,1)-sint*gv(1,2)
      x2= cost*gv(2,1)-sint*gv(2,2)
      x3= cost*gv(3,1)-sint*gv(3,2)
      y1= sint*gv(1,1)+cost*gv(1,2)
      y2= sint*gv(2,1)+cost*gv(2,2)
      y3= sint*gv(3,1)+cost*gv(3,2)
      z1=gv(1,3)
      z2=gv(2,3)
      z3=gv(3,3)
      if(phi .eq. 0.d0)then
        rho=0.d0
        gx0=gv(1,4)+z1*al
        gy0=gv(2,4)+z2*al
        gz0=gv(3,4)+z3*al
      else
        rho=al/phi
        sp0=sin(phi)
        cp0=cos(phi)
        r1=rho*sp0
        r2=rho*merge(sp0**2/(1.d0+cp0),1.d0-cp0,cp0 .ge. 0.d0)
        gx0=gv(1,4)+(r1*z1-r2*x1)
        gy0=gv(2,4)+(r1*z2-r2*x2)
        gz0=gv(3,4)+(r1*z3-r2*x3)
        z1=-sp0*x1+cp0*gv(1,3)
        x1= cp0*x1+sp0*gv(1,3)
        z2=-sp0*x2+cp0*gv(2,3)
        x2= cp0*x2+sp0*gv(2,3)
        z3=-sp0*x3+cp0*gv(3,3)
        x3= cp0*x3+sp0*gv(3,3)
c        write(*,'(a,1p12g10.2)')'tsetphgv ',gx0,gy0,gz0,
c     $       x1,x2,x3,y1,y2,y3,z1,z2,z3
      endif
      chi=merge(0.d0,2.d0*atan2(x3,-y3),x3 .eq. 0.d0)
      return
      end associate
      end subroutine

      subroutine tphrec(xi,pxi,yi,pyi,dp,dpr,ppx,ppy,p1,h1,rho,ds,k)
      use tfstk
      use tmacro, only:p0
      use mathfun, only:pxy2dpz
      use geolib
      implicit none
      integer*4 , intent(in)::k
      real*8 ,intent(in):: xi,pxi,yi,pyi,dp,dpr,ppx,ppy,p1,h1,ds,rho
      integer*4 irtc
      integer*8 kp
      real*8 cod(6),dpa(3),pxir,pzi,pyir,gv(3,4),gv1(3,4),fr,pp,
     $     thu,thv,xi30,xi1,xi2,xi3,c1,c2,s1,s2,ppx1,ppy1
      real*8 ,parameter :: frmin=1.d-12
      associate(l=>pcvt%l,cost=>pcvt%cost,sint=>pcvt%sint,al=>pcvt%al,
     $     fr0=>pcvt%fr0)
      fr=merge(fr0,fr0+ds/al,al .eq. 0.d0)
      gv=tfgeofrac(l,fr,irtc)
      if(irtc .ne. 0)then
        return
      endif
      call radangle(h1,rho*p1/p0,dpr,thu,thv,xi30,xi2)
      pp=hypot(ppx,ppy)
      ppx1=ppx/pp
      ppy1=ppy/pp
      pxir=pxi+( thv*ppx1+thu*ppy1)
      pyir=pyi+(-thv*ppy1+thu*ppx1)
      cod=tphchge([xi,pxir,yi,pyir,0.d0,0.d0])
c      cod=[xi,pxir,yi,pyir,0.d0,0.d0]
      gv1=tforbitgeo(gv,cod)
      pzi=1.d0+pxy2dpz(cod(2),cod(4))
      dpa=dp*matmul(gv(:,1:3),[cod(2),cod(4),pzi])
      c1=( cost*ppx1+sint*ppy1)/pp
      s1=(-sint*ppx1+cost*ppy1)/pp
      c2=(c1-s1)*(c1+s1)
      s2=2.d0*c1*s1
      xi3=s2*xi30
      xi1=c2*xi30
      if(ilp .eq. 0)then
        itp=itp+1
        kphtable(itp)=ktaloc(10*lt)
        ilp=1
      endif
      kp=kphtable(itp)+(ilp-1)*10
      ilist(1,kp)=k
      ilist(2,kp)=l
      rlist(kp+1:kp+3)=gv1(:,4)
      rlist(kp+4:kp+6)=dpa
      rlist(kp+7)=xi1
      rlist(kp+8)=xi2
      rlist(kp+9)=xi3
      ilp=ilp+1
      if(ilp .gt. lt)then
        ilp=0
      endif
      return
      end associate
      end subroutine

      subroutine tphotonconv(xi,pxi,yi,pyi,dp,dpr,p1,h1,ds,k)
      use tfstk
      use tmacro, only:p0
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 , intent(in)::k
      integer*8 kp
      real*8 xi,pxi,yi,pyi,dp,p1,h1,xi3a,gx,gy,gz,dpr,
     $     dpgx,dpgy,dpgz,ds,xir,pxir,zir,pzi,pyir,
     $     phi1,cp,sp,thu,thv,xi30,xi1,xi2,xi3,pxia
      associate(l=>pp%l,al=>pp%al,phi=>pp%phi,theta=>pp%theta,
     $     x1=>pp%geo1(1,1),x2=>pp%geo1(2,1),x3=>pp%geo1(3,1),
     $     y1=>pp%geo1(1,2),y2=>pp%geo1(2,2),y3=>pp%geo1(3,2),
     $     z1=>pp%geo1(1,3),z2=>pp%geo1(2,3),z3=>pp%geo1(3,3),
     $     gx0=>pp%geo1(1,4),gy0=>pp%geo1(2,4),gz0=>pp%geo1(3,4),
     $     rho=>pp%rho,chi=>pp%chi,geo1=>pp%geo1)
      call radangle(h1,rho*p1/p0,dpr,thu,thv,xi30,xi2)
      pxir=pxi+thu
      pyir=pyi+thv
      pzi=1.d0+pxy2dpz(pxir,pyir)
      if(phi .ne. 0.d0)then
        phi1=ds/rho
        cp=cos(phi1)
        sp=sin(phi1)
        xir=xi*cp
        zir=rho*sp
        xi3=(cp-sp)*(cp+sp)*xi30
        xi1=-2.d0*cp*sp*xi30
        pxia=pxir
        pxir=pxia*cp-pzi*sp
        pzi=pxia*sp+pzi*cp
      else
        xir=xi
        zir=ds
        xi3=xi30
        xi1=0.d0
      endif
      gx=gx0+xir*x1+yi*y1+zir*z1
      gy=gy0+xir*x2+yi*y2+zir*z2
      gz=gz0+xir*x3+yi*y3+zir*z3
      dpgx=dp*(pzi*z1+pxir*x1+pyir*y1)
      dpgy=dp*(pzi*z2+pxir*x2+pyir*y2)
      dpgz=dp*(pzi*z3+pxir*x3+pyir*y3)
c      write(*,'(a,2i5,1p5g12.4)')'phconv ',k,l,
c     $     pxia,pxir*x2,pyir*y2,pzi*z2,dpgy/dp
      xi3a=cos(chi)*xi3+sin(chi)*xi1
      xi1=-sin(chi)*xi3+cos(chi)*xi1
      xi3=xi3a
      if(ilp .eq. 0)then
        itp=itp+1
        kphtable(itp)=ktaloc(10*lt)
        ilp=1
      endif
c      write(*,*)'phconv ',itp,ilp
      kp=kphtable(itp)+(ilp-1)*10
      ilist(1,kp)=k
      ilist(2,kp)=l
      rlist(kp+1)=gx
      rlist(kp+2)=gy
      rlist(kp+3)=gz
      rlist(kp+4)=dpgx
      rlist(kp+5)=dpgy
      rlist(kp+6)=dpgz
      rlist(kp+7)=xi1
      rlist(kp+8)=xi2
      rlist(kp+9)=xi3
      ilp=ilp+1
      if(ilp .gt. lt)then
        ilp=0
      endif
      end associate
      end subroutine

      subroutine tgswap(l)
      use ffs_pointer, only:geo
      implicit none
      integer*4 , intent(in)::l
      real*8 v(3)
      v=geo(:,1,l)
      geo(:,1,l)=geo(:,2,l)
      geo(:,2,l)=v
      return
      end subroutine

      subroutine tphotonlist()
      use tfstk
      use ffs_pointer,only:gammab
      use tmacro
      implicit none
      type (sad_dlist), pointer ::klx
      type (sad_rlist), pointer ::klri
      integer*4 nitem
      parameter (nitem=12)
      integer*8 kax, kp,kt
      integer*4 nph,i
      real*8 dp
      integer*8 ,save ::kphlist=0
      if(kphlist .eq. 0)then
        kphlist=ktfsymbolz('`PhotonList',11)-4
      endif
      call tflocal(klist(kphlist))
      if(itp .le. 0)then
        dlist(kphlist)=dxnulll
      else
        nph=(itp-1)*lt+max(ilp-1,0)
c        write(*,*)'phlist ',itp,nph
        kax=ktadaloc(-1,nph,klx)
        klx%attr=ior(klx%attr,lconstlist)
        itp=1
        ilp=0
        kt=kphtable(1)
        do i=1,nph
          ilp=ilp+1
          if(ilp .gt. lt)then
            ilp=1
            itp=itp+1
            kt=kphtable(itp)
          endif
          kp=kt+(ilp-1)*10
          klx%dbody(i)%k=ktflist+ktavaloc(0,nitem,klri)
          klri%attr=lconstlist
          dp=norm2([rlist(kp+4),rlist(kp+5),rlist(kp+6)])
c          dp=hypot(rlist(kp+4),hypot(rlist(kp+5),rlist(kp+6)))
c          dp=sqrt(rlist(kp+4)**2+rlist(kp+5)**2
c     $         +rlist(kp+6)**2)
          klri%rbody(1)=dp*amass*gammab(ilist(2,kp))
          klri%rbody(2)=rlist(kp+1)
          klri%rbody(3)=rlist(kp+2)
          klri%rbody(4)=rlist(kp+3)
          klri%rbody(5)=rlist(kp+4)/dp
          klri%rbody(6)=rlist(kp+5)/dp
          klri%rbody(7)=rlist(kp+6)/dp
          klri%rbody(8)=rlist(kp+7)
          klri%rbody(9)=rlist(kp+8)
          klri%rbody(10)=rlist(kp+9)
          klri%rbody(11)=ilist(1,kp)
          klri%rbody(12)=ilist(2,kp)
        enddo
        do i=1,itp
          if(kphtable(i) .ne. 0)then
            call tfree(kphtable(i))
          endif
        enddo
        klist(kphlist)=ktflist+ktfcopy1(kax)
      endif
c      call tfdebugprint(klist(kphlist),'phlist',1)
c      write(*,*)'with ',itp,ilp
      return
      end subroutine

      real*8 function tphchge(cod) result(codx)
      use tmacro, only:irad
      implicit none
      dimension codx(6)
      real*8 ,intent(in):: cod(6)
      real*8 trans(6,12),beam(42),srot(3,9)
      integer*4 ir0
      ir0=irad
      irad=0
      codx=cod
      call tchge(trans,codx,beam,srot,
     $     -pcvt%dx,-pcvt%dy,-pcvt%theta,-pcvt%dtheta,-pcvt%phi0,
     $     .false.)
      irad=ir0
      return
      end function 

      subroutine tsetpcvt(l,dx,dy,theta,dtheta,phi0,al)
      implicit none
      integer*4 ,intent(in):: l
      real*8 ,intent(in):: dx,dy,theta,dtheta,phi0,al
      pcvt%l=l
      pcvt%dx=dx
      pcvt%dy=dy
      pcvt%theta=theta
      pcvt%dtheta=dtheta
      pcvt%phi0=phi0
      pcvt%al=al
      if(phi0 .eq. 0.d0)then
        pcvt%cost=cos(theta+dtheta)
        pcvt%sint=sin(theta+dtheta)
      else
        pcvt%cost=cos(theta)
        pcvt%sint=sin(theta)
      endif
      pcvt%fr0=0.d0
      return
      end subroutine
      
      end module

      module kradlib
      real*8 , allocatable :: pxr0(:),pyr0(:),zr0(:),bsi(:)

      contains
        subroutine tradkf1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     px00,py0,zr00,bsi,al,k)
        use ffs_flag
        use tmacro
        use photontable, only:tphrec
        use mathfun, only:pxy2dpz,p2h
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0,sprot
        implicit none
        integer*4 ,parameter :: npmax=10000
        integer*4 , intent(in)::k
        integer*4 i
        real*8 , intent(inout)::x,px,y,py,z,g,dv
        real*8 , intent(in)::px00,py0,zr00,bsi,al
        real*8 dpx,dpy,dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,p2,h2,sx,sy,sz,
     $       ppa,an,dph,r1,r2,px0,xr,yr,rho
        real*8 dpr(npmax),rph(npmax)
        dpz0=pxy2dpz(px00,py0)
        px0= cphi0*px00+sphi0*(1.d0+dpz0)
        dpz0=cphi0*dpz0-sphi0*px00
        dpx=px-px0
        dpy=py-py0
        dpz=pxy2dpz(px,py)
        dpz0=pxy2dpz(px0,py0)
        ppx=py*dpz0-dpz*py0+dpy
        ppy=dpz*px0-px*dpz0-dpx
        ppz=px*py0-py*px0
c        ppa=hypot(ppx,hypot(ppy,ppz))
        ppa=norm2([ppx,ppy,ppz])
        theta=asin(min(1.d0,max(-1.d0,ppa)))
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        anp=anrad*h1*theta
        al1=al-z+zr00
        rho=al1/theta
        if(photons)then
          call tdusrnpl(anp,dph,r1,r2,an,dpr,rph)
        else
          call tdusrn(anp,dph,r1,r2,an)
        endif
        if(an .ne. 0.d0)then
          uc=cuc*h1**3/p0/rho
          if(photons)then
            do i=1,int(an)
              dg=dpr(i)*uc
              xr=x-rph(i)*(px-.5d0*dpx*rph(i))*al
              yr=y-rph(i)*(py-.5d0*dpy*rph(i))*al
              call tphrec(xr,px,yr,py,dg,dpr(i),
     $             ppx,ppy,p,h1,rho,al-rph(i)*al,k)
c              call tphotonconv(xr,px,yr,py,dg,
c     $             dpr(i),p,h1,-rph(i)*al,k)
            enddo
          endif
          dg=-dph*uc
          dg=dg/(1.d0-2.d0*dg)
          g=max(gmin,g+dg)
          ddpx=-r1*dpx*dg
          ddpy=-r1*dpy*dg
          x=x+r2*ddpx*al1
          y=y+r2*ddpy*al1
          px=px+ddpx
          py=py+ddpy
          pr=1.d0+g
          p2=p0*pr
          h2=p2h(p2)
          dv=-g*(1.d0+pr)/h2/(h2+p2)+dvfs
          z=z*p2/h2*h1/p
          if(calpol)then
            pxm=px0+dpx*.5d0
            pym=py0+dpy*.5d0
            call sprot(sx,sy,sz,pxm,pym,
     $           ppx,ppy,ppz,bsi,merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0),
     $           h1,p2*h2/al1,an)
          endif
        elseif(calpol)then
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi,
     $         merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0),h1,
     $         p*h1/al1,-1.d0)
        endif
        return
        end subroutine

        subroutine tradkfn(np,xn,pxn,yn,pyn,zn,gn,dvn,sxn,syn,szn,al)
        use ffs_flag
        use tmacro
        use photontable, only:tphrec
        use mathfun, only:pxy2dpz,p2h
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0,sprot
        implicit none
        integer*4 ,parameter :: npmax=10000
        integer*4 , intent(in)::np
        integer*4 i,k
        real*8 , intent(inout)::
     $       xn(np),pxn(np),yn(np),pyn(np),zn(np),gn(np),dvn(np),
     $       sxn(np),syn(np),szn(np)
        real*8 , intent(in)::al
        real*8 dpx,dpy,dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,p2,h2,
     $       ppa,an,dph,r1,r2,px0,xr,yr,rho
        real*8 dpr(npmax),rph(npmax)
        do k=1,np
          dpz0=pxy2dpz(pxr0(k),pyr0(k))
          px0= cphi0*pxr0(k)+sphi0*(1.d0+dpz0)
          dpz0=cphi0*dpz0-sphi0*pxr0(k)
          dpx=pxn(k)-px0
          dpy=pyn(k)-pyr0(k)
          dpz=pxy2dpz(pxn(k),pyn(k))
          dpz0=pxy2dpz(px0,pyr0(k))
          ppx=pyn(k)*dpz0-dpz*pyr0(k)+dpy
          ppy=dpz*px0-pxn(k)*dpz0-dpx
          ppz=pxn(k)*pyr0(k)-pyn(k)*px0
          ppa=norm2([ppx,ppy,ppz])
          theta=asin(min(1.d0,max(-1.d0,ppa)))
          pr=1.d0+gn(k)
          p=p0*pr
          h1=p2h(p)
          anp=anrad*h1*theta
          al1=al-zn(k)+zr0(k)
          rho=al1/theta
          if(photons)then
            call tdusrnpl(anp,dph,r1,r2,an,dpr,rph)
          else
            call tdusrn(anp,dph,r1,r2,an)
          endif
          if(an .ne. 0.d0)then
            uc=cuc*h1**3/p0/rho
            if(photons)then
              do i=1,int(an)
                dg=dpr(i)*uc
                xr=xn(k)-rph(i)*(pxn(k)-.5d0*dpx*rph(i))*al
                yr=yn(k)-rph(i)*(pyn(k)-.5d0*dpy*rph(i))*al
                call tphrec(xr,pxn(k),yr,pyn(k),dg,dpr(i),
     $               ppx,ppy,p,h1,rho,al-rph(i)*al,k)
c                call tphotonconv(xr,pxn(k),yr,pyn(k),dg,
c     $               dpr(i),p,h1,-rph(i)*al,k)
              enddo
            endif
            dg=-dph*uc
            dg=dg/(1.d0-2.d0*dg)
            gn(k)=max(gmin,gn(k)+dg)
            ddpx=-r1*dpx*dg
            ddpy=-r1*dpy*dg
            xn(k)=xn(k)+r2*ddpx*al1
            yn(k)=yn(k)+r2*ddpy*al1
            pxn(k)=pxn(k)+ddpx
            pyn(k)=pyn(k)+ddpy
            pr=1.d0+gn(k)
            p2=p0*pr
            h2=p2h(p2)
            dvn(k)=-gn(k)*(1.d0+pr)/h2/(h2+p2)+dvfs
            zn(k)=zn(k)*p2/h2*h1/p
            if(calpol)then
              pxm=px0    +dpx*.5d0
              pym=pyr0(k)+dpy*.5d0
              call sprot(sxn(k),syn(k),szn(k),pxm,pym,
     $             ppx,ppy,ppz,bsi(k),
     $             merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0),h1,
     $             p2*h2/al1,an)
            endif
          elseif(calpol)then
            pxm=px0    +dpx*.5d0
            pym=pyr0(k)+dpy*.5d0
            call sprot(sxn(k),syn(k),szn(k),pxm,pym,ppx,ppy,ppz,
     $           bsi(k),merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0),
     $           h1,p*h1/al1,-1.d0)
          endif
        enddo
        return
        end subroutine

        subroutine tradk1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     px00,py0,zr00,bsi0,al)
        use ffs_flag
        use tmacro
        use mathfun, only:pxy2dpz,p2h
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0,sprot
        implicit none
        real*8 , intent(inout)::x,px,y,py,z,g,dv
        real*8 , intent(in)::px00,py0,zr00,bsi0,al
        real*8 dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,
     $       px0,pxm,pym,al1,uc,ddpx,ddpy,h2,h1,sx,sy,sz,ppa,p2
        dpz0=pxy2dpz(px00,py0)
        px0= cphi0*px00+sphi0*(1.d0+dpz0)
        dpz0=cphi0*dpz0-sphi0*px00
        dpx=px-px0
        dpy=py-py0
        dpz=pxy2dpz(px,py)
        ppx=py*dpz0-dpz*py0+dpy
        ppy=dpz*px0-px*dpz0-dpx
        ppz=px*py0-py*px0
        ppa=norm2([ppx,ppy,ppz])
        theta=asin(min(1.d0,max(-1.d0,ppa)))
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        al1=al-z+zr00
        anp=anrad*h1*theta
        uc=cuc*h1**3/p0*theta/al1
        dg=-cave*anp*uc
        dg=dg/(1.d0-2.d0*dg)
        g=max(gmin,g+dg)
        ddpx=-.5d0*dpx*dg
        ddpy=-.5d0*dpy*dg
        x=x+ddpx*al1/3.d0
        y=y+ddpy*al1/3.d0
        px=px+ddpx
        py=py+ddpy
        pr=1.d0+g
        p2=p0*pr
        h2=p2h(p2)
        dv=-g*(1.d0+pr)/h2/(h2+p2)+dvfs
        z=z*p2/h2*h1/p
        if(calpol)then
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi0,
     $         merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0),
     $         h2,p2*h2/al1,anp)
        endif
        return
        end subroutine

        subroutine tradkn(np,xn,pxn,yn,pyn,zn,gn,dvn,sxn,syn,szn,al)
        use ffs_flag
        use tmacro
        use mathfun, only:pxy2dpz,p2h
        use tspin
        implicit none
        integer*4 , intent(in)::np
        real*8 , intent(inout)::
     $       xn(np),pxn(np),yn(np),pyn(np),zn(np),gn(np),dvn(np),
     $       sxn(np),syn(np),szn(np)
        real*8 , intent(in)::al
        real*8 dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,
     $       px0,pxm,pym,al1,uc,ddpx,ddpy,h2,h1,ppa,p2
        integer*4 i
        do i=1,np
          dpz0=pxy2dpz(pxr0(i),pyr0(i))
          px0= cphi0*pxr0(i)+sphi0*(1.d0+dpz0)
          dpz0=cphi0*dpz0-sphi0*pxr0(i)
          dpx=pxn(i)-px0
          dpy=pyn(i)-pyr0(i)
          dpz=pxy2dpz(pxn(i),pyn(i))
          ppx=pyn(i)*dpz0-dpz*pyr0(i)+dpy
          ppy=dpz*px0-pxn(i)*dpz0-dpx
          ppz=pxn(i)*pyr0(i)-pyn(i)*px0
          ppa=norm2([ppx,ppy,ppz])
          theta=asin(min(1.d0,max(-1.d0,ppa)))
          pr=1.d0+gn(i)
          p=p0*pr
          h1=p2h(p)
          al1=al-zn(i)+zr0(i)
          anp=anrad*h1*theta
          uc=cuc*h1**3/p0*theta/al1
          dg=-cave*anp*uc
          dg=dg/(1.d0-2.d0*dg)
          gn(i)=max(gmin,gn(i)+dg)
          ddpx=-.5d0*dpx*dg
          ddpy=-.5d0*dpy*dg
          xn(i)=xn(i)+ddpx*al1/3.d0
          yn(i)=yn(i)+ddpy*al1/3.d0
          pxn(i)=pxn(i)+ddpx
          pyn(i)=pyn(i)+ddpy
          pr=1.d0+gn(i)
          p2=p0*pr
          h2=p2h(p2)
          dvn(i)=-gn(i)*(1.d0+pr)/h2/(h2+p2)+dvfs
          zn(i)=zn(i)*p2/h2*h1/p
          if(calpol)then
            pxm=px0    +dpx*.5d0
            pym=pyr0(i)+dpy*.5d0
            call sprot(sxn(i),syn(i),szn(i),pxm,pym,ppx,ppy,ppz,
     $           bsi(i),merge(theta/ppa*pr,0.d0,ppa .ne. 0.d0)
     $           ,h2,p2*h2/al1,anp)
          endif
        enddo
        return
        end subroutine

        subroutine tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phi0)
        use tfstk
        use ffs_flag, only:calpol,rfluct
        use tspin
        implicit none
        integer*4 , intent(in)::np
        real*8 ,intent(inout)::
     $       x(np),px(np),y(np),py(np),dv(np),z(np),g(np),
     $       sx(np),sy(np),sz(np)
        real*8 , intent(in)::al,phi0
        if(al .ne. 0.d0)then
          cphi0=cos(phi0)
          sphi0=sin(phi0)
          if(rfluct)then
            call tradkfn(np,x,px,y,py,z,g,dv,sx,sy,sz,al)
          else
            call tradkn(np,x,px,y,py,z,g,dv,sx,sy,sz,al)
          endif
        endif
        pxr0=px
        pyr0=py
        zr0=z
        if(calpol)then
          bsi=0.d0
        endif
        return
        end subroutine

        subroutine tradke(trans,cod,beam,srot,al,phir0,bzh)
        use tmacro
        use temw,only:codr0,bzhr0,bsir0,calint,tinv6,gintd,transr,
     $       tmulbs
        use ffs_flag,only:radcod,calpol
        use mathfun, only:pxy2dpz,p2h
        use tspin, only:cave,cl,cuu,gmin,sflc
        implicit none
        real*8 , intent(inout)::trans(6,12),cod(6),beam(42),
     $       srot(3,9)
        real*8 , intent(in)::al,bzh,phir0
        real*8 transi(6,6),tr1(6,6),dxpa(6),tr2(6,6),
     $       ddpz(6),dal(6),duc(6),dddpx(6),dddpy(6),ddg(6),
     $       dtheta(6),danp(6),dbeam(21),dpxi(6),dpyi(6),
     $       c1,dpx,dpy,ddpx,ddpy,pxr0,ct,pz00,das,bt,
     $       pr,px,py,pz,pz0,xpx,xpy,xpz,xpa,theta,th,
     $       p,h1,al1,anp,uc,dg,g,pr1,pxi,pyi,
     $       p2,h2,de,cp,sp,b,pxm,pym,gi,dh1r,
     $       pxh,pyh,pzh,xpzb,btx,bty,btz,dct,sinu,cosu,dcosu,
     $       gx,gy,gz,blx,bly,blz,
     $       sx(9),sy(9),sz(9),sux(9),suy(9),suz(9),sw(9),
     $       dpxh(6),dpyh(6),dpzh(6),bp,dbp(6),dpxr0(6),dpz0(6),
     $       dxpx(6),dxpy(6),dxpz(6),dxpzb(6),dblx(6),dbly(6),dblz(6),
     $       dbtx(6),dbty(6),dbtz(6),dgx(6),dgy(6),dgz(6),dpz00(6)
        gi=codr0(6)
        pr=1.d0+gi
        th=tan(.5d0*phir0)
        sp=2.d0*th/(1.d0+th**2)
        cp=1.d0-th*sp
c     codr0 has canonical momenta!
        pxi=codr0(2)+bzhr0*codr0(3)
        pyi=codr0(4)-bzhr0*codr0(1)
        pz00=pr*(1.d0+pxy2dpz(pxi/pr,pyi/pr))
        pxr0= cp*pxi+sp*pz00
        pz0 =-sp*pxi+cp*pz00
        px=cod(2)+bzh*cod(3)
        py=cod(4)-bzh*cod(1)
        pz=pr*(1.d0+pxy2dpz(px/pr,py/pr))
        dpx=px-pxr0
        dpy=py-pyi
        xpx=(py*pz0 -pz*pyi)
        xpy=(pz*pxr0-px*pz0)
        xpz=(px*pyi-py*pxr0)
        xpa=abs(dcmplx(xpx,abs(dcmplx(xpy,xpz))))/pr**2
        theta=asin(min(1.d0,xpa))
        p=p0*pr
        h1=p2h(p)
        al1=al-cod(5)+codr0(5)
        anp=anrad*h1*theta
        uc=cuc*h1**3/p0*theta/al1
        dg=-cave*anp*uc
        u0=u0-dg
        g=max(gmin,gi+dg)
        pr1=1.d0+g
        ddpx=.5d0*dpx*dg
        ddpy=.5d0*dpy*dg
        c1=al1/pr/3.d0
        if(radcod)then
          cod(1)=cod(1)+ddpx*c1
          cod(3)=cod(3)+ddpy*c1
          cod(2)=px*pr1/pr+ddpx-bzh*cod(3)
          cod(4)=py*pr1/pr+ddpy+bzh*cod(1)
          cod(6)=g
          p2=p0*pr1
          h2=p2h(p2)
          cod(5)=cod(5)*p2/h2*h1/p
          call tesetdv(g)
        else
          p2=p
          h2=h1
        endif
        if(irad .gt. 6)then
          transi=tinv6(transr)
c          call tinv6(transr,transi)
          transi=matmul(trans(:,1:6),transi)
c          call tmultr(transi,trans(:,1:6),6)
          tr2=transi
          if(bzh .ne. 0.d0)then
            tr2(2,:)=tr2(2,:)+bzh*tr2(3,:)
            tr2(4,:)=tr2(4,:)-bzh*tr2(1,:)
          endif
          ddpz=(tr2(6,:)*pr-tr2(2,:)*px-tr2(4,:)*py)/pz
          dpxi=(/0.d0,1.d0,bzhr0,0.d0,0.d0,0.d0/)
          dpyi=(/-bzhr0,0.d0,0.d0,1.d0,0.d0,0.d0/)
          dpz00=(-pxi*dpxi-pyi*dpyi)/pz00
          dpz00(6)=dpz00(6)+pr/pz00
          dpxr0=cp*dpxi+sp*dpz00
          dpz0=(-pxr0*dpxr0-pyi*dpyi)/pz0
          dpz0(6)=dpz0(6)+pr/pz0
          dxpx=tr2(4,:)*pz0+py*dpz0-ddpz*pyi-pz*dpyi
          dxpy=ddpz*pxr0+pz*dpxr0-tr2(2,:)*pz0-px*dpz0
          dxpz=tr2(2,:)*pyi+px*dpyi-tr2(4,:)*pxr0-py*dpxr0
          dh1r=p*p0/h1**2
          if(xpa .ne. 0.d0)then
            dxpa=(xpx*dxpx+xpy*dxpy+xpz*dxpz)/xpa/pr**2
            dxpa(6)=dxpa(6)-2.d0*xpa/pr
            dal=-tr2(5,:)
            dal(5)=dal(5)+1.d0
            das=1.d0/sqrt(1.d0-xpa**2)
            dtheta=dxpa*das
            danp=anrad*h1*dtheta
            danp(6)=danp(6)+anp*dh1r
            duc=uc*(dtheta/theta-dal/al1)
            duc(6)=duc(6)+3.d0*uc*dh1r
            ddg=-cave*(danp*uc+anp*duc)
            dddpx=.5d0*((tr2(2,:)-dpxr0)*dg+ddpx*ddg)
            dddpy=.5d0*((tr2(4,:)-dpyi )*dg+ddpy*ddg)
            tr1(1,:)=c1*dddpx
            tr1(1,6)=tr1(1,6)-ddpx/pr
            tr1(3,:)=c1*dddpy
            tr1(3,6)=tr1(3,6)-ddpy/pr
            tr1(2,:)=(tr2(2,:)*dg+px*ddg)/pr+dddpx
            tr1(2,6)=tr1(2,6)-px*dg/pr
            tr1(4,:)=(tr2(4,:)*dg+py*ddg)/pr+dddpy
            tr1(4,6)=tr1(4,6)-py*dg/pr
c     derivative of dz has been ignored.
            tr1(5,:)=0.d0
            tr1(6,:)=ddg
c     write(*,'(a,1p8g15.7)')'tradke  ',tr2(2,:)
c     write(*,'(a,1p8g15.7)')' ddg    ',ddg,dg
c     write(*,'(a,1p8g15.7)')' danp   ',danp,anp
c     write(*,'(a,1p8g15.7)')' duc    ',duc,uc
c     write(*,'(a,1p8g15.7)')' dtheta ',dtheta,theta
c     write(*,'(a,1p8g15.7)')' dxpy   ',dxpy,xpy
c     write(*,'(a,1p8g15.7)')' dpxr0  ',dpxr0,pxr0
c     do i=1,6
c     write(*,'(1p6g15.7)')tr1(i,:)
c     enddo
            if(bzh .ne. 0.d0)then
              tr1(2,:)=tr1(2,:)-bzh  *tr1(3,:)
              tr1(4,:)=tr1(4,:)+bzh  *tr1(1,:)
            endif
            call tmuld6(trans,tr1)
            tr1(1,1)=tr1(1,1)+1.d0
            tr1(2,2)=tr1(2,2)+1.d0
            tr1(3,3)=tr1(3,3)+1.d0
            tr1(4,4)=tr1(4,4)+1.d0
            tr1(5,5)=tr1(5,5)+1.d0
            tr1(6,6)=tr1(6,6)+1.d0
            call tmulbs(beam,tr1,calint)
            de=anp*uc**2*cuu
            pxm=pxi+px
            pym=pyi+py
            b=bzh*.5d0
            dbeam=0.d0
            dbeam(3)=(beam(3)+b*(2.d0*beam(5)+b*beam(6))
     $           +(pxm**2+pxi**2+px**2)/6.d0)*de
            dbeam(8) =(beam(8)-b*(beam(2)-beam(10)+b*beam(4))
     $           +(pxm*pym+pxi*pyi+px*py)/6.d0)*de
            dbeam(10)=(beam(10)+b*(-2.d0*beam(7)+b*beam(1))
     $           +(pym**2+pyi**2+py**2)/6.d0)*de
            dbeam(17)=pxm*de*.5d0
            dbeam(19)=pym*de*.5d0
            dbeam(21)=de
            beam(1:21)=beam(1:21)+dbeam
            if(calint)then
              beam(22:42)=beam(22:42)+dbeam
            endif
          endif
          if(calpol)then
            xpzb=xpz+(bsir0+bzh*2.d0*al)*pr**2
            dxpzb=dxpz
            dxpzb(6)=dxpzb(6)+2.d0*(bsir0+bzh*2.d0*al)*pr
            pxh=(pxr0+px)/pr*.5d0
            pyh=(pyi+py)/pr*.5d0
            dpxh=(tr2(2,:)+dpxr0)/pr*.5d0
            dpxh(6)=dpxh(6)-pxh/pr
            dpyh=(tr2(4,:)+dpyi)/pr*.5d0
            dpyh(6)=dpyh(6)-pyh/pr
            pzh=1.d0+pxy2dpz(pxh,pyh)
            dpzh=-(pxh*dpxh+pyh*dpyh)/pzh
            bp=(xpx*pxh+xpy*pyh+xpzb*pzh)/pr
            dbp=(dxpx*pxh+xpx*dpxh+dxpy*pyh
     $           +xpy*dpyh+dxpzb*pzh+xpzb*dpzh)/pr
            dbp(6)=dbp(6)-bp/pr
            blx=bp*pxh
            bly=bp*pyh
            blz=bp*pzh
            btx=xpx/pr-blx
            bty=xpy/pr-bly
            btz=xpzb/pr-blz
            dblx=dbp*pxh+bp*dpxh
            dbly=dbp*pyh+bp*dpyh
            dblz=dbp*pzh+bp*dpzh
            dbtx=dxpx/pr-dblx
            dbtx(6)=dbtx(6)-xpx/pr**2
            dbty=dxpy/pr-dbly
            dbty(6)=dbty(6)-xpy/pr**2
            dbtz=dxpzb/pr-dblz
            dbtz(6)=dbtz(6)-xpzb/pr**2
            ct=h1*gspin
            dct=ct*dh1r
            ct=ct+1.d0
            gx=ct*btx+cl*blx
            gy=ct*bty+cl*bly
            gz=ct*btz+cl*blz
            dgx=ct*dbtx+cl*dblx
            dgx(6)=dgx(6)+dct*btx
            dgy=ct*dbty+cl*dbly
            dgy(6)=dgy(6)+dct*bty
            dgz=ct*dbtz+cl*dblz
            dgz(6)=dgz(6)+dct*btz
            srot(1,4:9)=srot(1,4:9)
     $           +dgx(1)*transr(1,:)+dgx(2)*transr(2,:)
     $           +dgx(3)*transr(3,:)+dgx(4)*transr(4,:)
     $           +dgx(5)*transr(5,:)+dgx(6)*transr(6,:)
            srot(2,4:9)=srot(2,4:9)
     $           +dgy(1)*transr(1,:)+dgy(2)*transr(2,:)
     $           +dgy(3)*transr(3,:)+dgy(4)*transr(4,:)
     $           +dgy(5)*transr(5,:)+dgy(6)*transr(6,:)
            srot(3,4:9)=srot(3,4:9)
     $           +dgz(1)*transr(1,:)+dgz(2)*transr(2,:)
     $           +dgz(3)*transr(3,:)+dgz(4)*transr(4,:)
     $           +dgz(5)*transr(5,:)+dgz(6)*transr(6,:)
            g=norm2([gx,gy,gz])
            if(g .ne. 0.d0)then
              bt=abs(dcmplx(btx,abs(dcmplx(bty,btz))))
              th=tan(.5d0*g)
              sinu=2.d0*th/(1.d0+th**2)
              dcosu=th*sinu
c              sinu=sin(g)
c              dcosu=2.d0*sin(g*.5d0)**2
              cosu=1.d0-dcosu
              sinu=sinu/g
              sx=srot(1,:)
              sy=srot(2,:)
              sz=srot(3,:)
              sw=(sx*gx+sy*gy+sz*gz)/g
              gintd=gintd+sw(1:3)*sflc*anp*(bt*h1*p/al1)**2
              sw=sw*dcosu/g
              sux=sy*gz-sz*gy
              suy=sz*gx-sx*gz
              suz=sx*gy-sy*gx
              sx       =cosu*sx+sinu*sux+sw*gx
              srot(2,:)=cosu*sy+sinu*suy+sw*gy
              sz       =cosu*sz+sinu*suz+sw*gz
              srot(1,:)= cp*sx+sp*sz
              srot(3,:)=-sp*sx+cp*sz
            else
              srot(1,:)=  cp*srot(1,:)+sp*srot(3,:)
              srot(3,:)=(-sp*srot(1,:)+srot(3,:))/cp
            endif
          endif
        endif
        codr0(1:6)=cod(1:6)
        transr=trans(:,1:6)
        bzhr0=bzh
        bsir0=0.d0
        return
        end subroutine

        subroutine tallocvar(var,np)
        implicit none
        integer*4 ,intent(in):: np
        real*8 ,allocatable :: var(:)
        if(.not. allocated(var))then
          allocate(var(np))
        elseif(sizeof(var) .lt. np*8)then
          deallocate(var)
          allocate(var(np))
        endif
        return
        end

        subroutine tallocrad(np)
        use ffs_flag, only:rad
        implicit none
        integer*4 ,intent(in):: np
        call tallocvar(bsi,np)
        if(rad)then
          call tallocvar(pxr0,np)
          call tallocvar(pyr0,np)
          call tallocvar(zr0,np)
        endif
       return
       end

      end module
