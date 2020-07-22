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
      use tmacro
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
      if(ini)then
        gv=geo(:,:,l)
c        call tggeol(l,gv)
      else
        gv=geo1
      endif
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
        if(cp0 .ge. 0.d0)then
          r2=rho*sp0**2/(1.d0+cp0)
        else
          r2=rho*(1.d0-cp0)
        endif
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
      if(x3 .eq. 0.d0)then
        chi=0.d0
      else
        chi=2.d0*atan2(x3,-y3)
      endif
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
      if(al .eq. 0.d0)then
        fr=fr0
      else
        fr=fr0+ds/al
      endif
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
      use mathfun, only:hypot3
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
          dp=hypot3(rlist(kp+4),rlist(kp+5),rlist(kp+6))
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
