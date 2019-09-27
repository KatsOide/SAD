      module multa
        integer*4, parameter :: nmult=21
        logical*4 :: gknini=.true.
        real*8 :: gkn(0:nmult,0:nmult)=0.d0

        contains
        subroutine gkninit
        implicit none
        integer*4 n,k
        do n=0,nmult
          gkn(n,0)=1.d0
          do k=1,nmult-n
            gkn(n,k)=gkn(n,k-1)
     $           *dble((2*k-3)*(2*k+1))/8.d0/dble(k*(k+n+1))
          enddo
        enddo
        gknini=.false.
        return
        end subroutine
      end module

      module tbendcom
        real*8 rho0,rhob,f1r,f2r,fb1,fb2

        contains
        subroutine tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
        use tfstk
        use ffs_flag, only:calpol
        use mathfun
        implicit none
        integer*4 np,i
        real*8 phi0,dtheta,
     $       r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $       pxi,pyi,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf,
     $       cphi0,sphi0,cdt,sdt,sdth2,sxf,syf
        real*8 x(np),px(np),y(np),py(np),z(np),sx(np),sy(np),sz(np)
        cphi0=cos(phi0*.5d0)
        sphi0=sin(phi0*.5d0)
        sdth2=sin(dtheta*.5d0)**2
        cdt=1.d0-2.d0*sdth2
        sdt=sin(dtheta)
        r11=cdt*cphi0**2+sphi0**2
        r12=-cphi0*sdt
        r13=-2.d0*sdth2*sphi0*cphi0
        r21=-r12
        r22=cdt
        r23=sphi0*sdt
        r31=r13
        r32=-r23
        r33=cdt*sphi0**2+cphi0**2
        if(calpol)then
          do i=1,np
            xi=x(i)
            yi=y(i)
            pxi=px(i)
            pyi=py(i)
            pzi=1.d0+pxy2dpz(pxi,pyi)
            xf =r11*xi +r12*yi
            yf =r21*xi +r22*yi
            zf =r31*xi +r32*yi
            pxf=r11*pxi+r12*pyi+r13*pzi
            pyf=r21*pxi+r22*pyi+r23*pzi
            pzf=r31*pxi+r32*pyi+r33*pzi
            px(i)=pxf
            py(i)=pyf
            x(i)=xf-pxf/pzf*zf
            y(i)=yf-pyf/pzf*zf
            z(i)=z(i)+zf/pzf
            sxf  =r11*sx(i)+r12*sy(i)+r13*sz(i)
            syf  =r21*sx(i)+r22*sy(i)+r23*sz(i)
            sz(i)=r31*sx(i)+r32*sy(i)+r33*sz(i)
            sx(i)=sxf
            sy(i)=syf
          enddo
        else
          do i=1,np
            xi=x(i)
            yi=y(i)
            pxi=px(i)
            pyi=py(i)
            pzi=1.d0+pxy2dpz(pxi,pyi)
            xf =r11*xi +r12*yi
            yf =r21*xi +r22*yi
            zf =r31*xi +r32*yi
            pxf=r11*pxi+r12*pyi+r13*pzi
            pyf=r21*pxi+r22*pyi+r23*pzi
            pzf=r31*pxi+r32*pyi+r33*pzi
            px(i)=pxf
            py(i)=pyf
            x(i)=xf-pxf/pzf*zf
            y(i)=yf-pyf/pzf*zf
            z(i)=z(i)+zf/pzf
          enddo
        endif
        return
        end subroutine

      end module


      module photontable
      implicit none
      type photonp
        sequence
        real*8 al,phi,theta,rho,sp0,cp0,
     $       r1,r2,cost,sint,chi,geo1(3,4)        
        integer*4 l
      end type
      integer*4 ntable,ltable
      parameter (ntable=256,ltable=100000)
      integer*8 kphtable(ntable)
      integer*4 nt,itp,ilp,lt
      type (photonp) pp

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

      subroutine tsetphotongeo(al0,phi0,theta0,ini)
      use tmacro, only:l_track
      use ffs_pointer,only:ibzl,geo
      implicit none
      real*8, intent(in):: al0,phi0,theta0
      logical*4 , intent(in)::ini
      integer*4 ls
      real*8 gv(3,4)
      associate(l=>pp%l,al=>pp%al,phi=>pp%phi,theta=>pp%theta,
     $     x1=>pp%geo1(1,1),x2=>pp%geo1(2,1),x3=>pp%geo1(3,1),
     $     y1=>pp%geo1(1,2),y2=>pp%geo1(2,2),y3=>pp%geo1(3,2),
     $     z1=>pp%geo1(1,3),z2=>pp%geo1(2,3),z3=>pp%geo1(3,3),
     $     gx0=>pp%geo1(1,4),gy0=>pp%geo1(2,4),gz0=>pp%geo1(3,4),
     $     rho=>pp%rho,sp0=>pp%sp0,cp0=>pp%cp0,
     $     r1=>pp%r1,r2=>pp%r2,cost=>pp%cost,sint=>pp%sint,
     $     chi=>pp%chi,geo1=>pp%geo1)
      l=l_track
      if(ini)then
        call tggeol(l,gv)
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
        write(*,'(a,1p12g10.2)')'tsetphgv ',gx0,gy0,gz0,
     $       x1,x2,x3,y1,y2,y3,z1,z2,z3
      endif
      if(x3 .eq. 0.d0)then
        chi=0.d0
      else
        chi=2.d0*atan2(x3,-y3)
      endif
      return
      end associate
      end subroutine

      subroutine tphotonconv(xi,pxi,yi,pyi,dp,dpr,p1,h1,ds,k)
      use tfstk
      use tmacro, only:brho,p0
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 , intent(in)::k
      integer*8 kp
      real*8 xi,pxi,yi,pyi,dp,p1,h1,xi3a,gx,gy,gz,dpr,
     $     dpz,dpgx,dpgy,dpgz,ds,xir,pxir,zir,pzi,pyir,
     $     al1,phi1,cp,sp,thu,thv,xi30,xi1,xi2,xi3,pxia
      associate(l=>pp%l,al=>pp%al,phi=>pp%phi,theta=>pp%theta,
     $     x1=>pp%geo1(1,1),x2=>pp%geo1(2,1),x3=>pp%geo1(3,1),
     $     y1=>pp%geo1(1,2),y2=>pp%geo1(2,2),y3=>pp%geo1(3,2),
     $     z1=>pp%geo1(1,3),z2=>pp%geo1(2,3),z3=>pp%geo1(3,3),
     $     gx0=>pp%geo1(1,4),gy0=>pp%geo1(2,4),gz0=>pp%geo1(3,4),
     $     rho=>pp%rho,sp0=>pp%sp0,cp0=>pp%cp0,
     $     r1=>pp%r1,r2=>pp%r2,cost=>pp%cost,sint=>pp%sint,
     $     chi=>pp%chi,geo1=>pp%geo1)
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

      end module

      subroutine tphotonlist()
      use photontable
      use tfstk
      use tmacro
      implicit none
      type (sad_dlist), pointer ::klx
      type (sad_rlist), pointer ::klri
      integer*4 nitem
      parameter (nitem=12)
      integer*8 kax, kp,kt
      integer*4 nph,i
      real*8 dp
      integer*8 kphlist
      data kphlist/0/
      if(kphlist .eq. 0)then
        kphlist=ktfsymbolz('`PhotonList',11)-4
      endif
      call tflocal(klist(kphlist))
      if(itp .le. 0)then
        kax=kxnulll
      else
        nph=(itp-1)*lt+max(ilp-1,0)
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
          dp=sqrt(rlist(kp+4)**2+rlist(kp+5)**2
     $         +rlist(kp+6)**2)
          klri%rbody(1)=dp*amass
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
      endif
      klist(kphlist)=ktflist+ktfcopy1(kax)
      return
      end

      subroutine tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb10,fb20,mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     enarad,eps)
      use ffs_flag,only:rad
      use tspin
      implicit none
      integer*4 np,mfring,l
      real*8 al,phib,phi0,cosp1,sinp1,cosp2,sinp2,ak,dx,dy,theta,
     $     cost,sint,cosw,sinw,sqwh,sinwp1,eps,
     $     fb10,fb20,dtheta
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      logical*4 enarad,fringe
      if(rad .and. enarad)then
        bsi=0.d0
      endif
      call tbend0(np,x,px,y,py,z,g,dv,sx,sy,sz,px0,py0,zr0,bsi,
     $     al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb10,fb20,mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     enarad,eps,.true.,0)
      return
      end

      subroutine tbend0(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,
     $     al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb10,fb20,mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     enarad,eps,ini,iniph)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev,geo
      use multa, only:nmult
      use tbendcom
      use tspin
      use photontable
      implicit none
      integer*4 np,mfring,i,ndiv,ndivmax,iniph
      parameter (ndivmax=1024)
      real*8 al,phib,phi0,cosp1,sinp1,cosp2,sinp2,ak,dx,dy,theta,
     $     cost,sint,cosw,sinw,sqwh,sinwp1,eps,
     $     xi,pxi,psi1,psi2,fb10,fb20,dtheta
      real*8 a3,a5,a7,a9,a11,a13,a15
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      real*8 smax,smin,rphidiv
      parameter (smax=0.99d0,smin=0.01d0,rphidiv=3e-3)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      complex*16 akm(0:nmult)
      logical*4 enarad,fringe,ini
      if(phi0 .eq. 0.d0)then
        if(ak .eq. 0.d0)then
          call tsteer(np,x,px,y,py,z,g,dv,sx,sy,sz,al,-phib,
     1         dx,dy,theta+dtheta,cos(theta+dtheta),sin(theta+dtheta),
     1         cosp1,sinp1,cosp2,sinp2,
     $         fb10,fb20,mfring,fringe,eps,enarad)
        elseif(phib .eq. phi0)then
          call tquad(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,
     1         dx,dy,theta+dtheta,
     $         cos(theta+dtheta),sin(theta+dtheta),0.d0,.true.,
     1         fringe,0.d0,0.d0,0,eps,.true.)
        else
          akm=0.d0
          akm(0)=phib-phi0
          akm(1)=ak
          call tmulti(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al,ak,0.d0,0.d0,
     $         psi1,psi2,
     $         dx,dy,0.d0,0.d0,0.d0,theta+dtheta,0.d0,
     $         eps,enarad,fringe,
     $         0.d0,0.d0,0.d0,0.d0,
     $         mfring,fb10,fb20,
     $         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     $         .false.,.false.,
     $         i00,i00,i00,i00)
        endif
        return
c      elseif(rad .and. enarad .and. trpt)then
c        call tbrad(np,x,px,y,py,z,g,dv,sx,sy,sz,l,al,phib,phi0,
c     1       cosp1,sinp1,cosp2,sinp2,
c     1       ak,dx,dy,theta,dtheta,cost,sint,
c     1       fb10,fb20,mfring,
c     1       fringe,eps)
c        return
      elseif(ak .ne. 0.d0)then
        call tbendi(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phib,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       ak,dx,dy,theta,dtheta,cost,sint,
     1       fb10,fb20,mfring,enarad,fringe,eps)
        return
      endif
      include 'inc/TENT.inc'
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      if(phib .eq. 0.d0)then
        call tbdrift(np,x,px,y,py,z,dv,sx,sz,al,phi0)
        go to 9000
      elseif(al .eq. 0.d0)then
        call tbthin(np,x,px,y,py,z,g,sx,sy,sz,phib,phi0,dx,dy,
     1              theta,dtheta,cost,sint)
        go to 9000
      endif
      rhob=al/phib
      rho0=al/phi0
      fb1=fb10
      fb2=fb20
      if(rad .and. enarad)then
        if(ini)then
          px0=px
          py0=py
          zr0=z
        endif
        if(iprev(l_track) .eq. 0)then
          f1r=fb1
        else
          f1r=0.d0
        endif
        if(inext(l_track) .eq. 0)then
          f2r=fb2
        else
          f2r=0.d0
        endif
        ndiv=min(ndivmax,ndivrad(phib,0.d0,0.d0,eps))
        if(photons)then
          select case (iniph)
            case (0)
              call tsetphotongeo(al/ndiv,phi0/ndiv,theta,.true.)
            case (1)
              call tsetphotongeo(al/ndiv,phi0/ndiv,pp%theta,.true.)
            case default
              call tsetphotongeo(al/ndiv,phi0/ndiv,0.d0,.false.)
          end select
        endif
      else
        ndiv=1
      endif
c      write(*,*)'tbend0 ',ndiv
      if(ndiv .gt. 1)then
        call tbendr(np,x,px,y,py,z,g,dv,sx,sy,sz,px0,py0,zr0,bsi,
     $       al,phib,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       mfring,fringe,ndiv)
      else
        call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       px0,py0,zr0,bsi,
     $       al,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       mfring,fringe,
     $       cosw,sinw,sqwh,sinwp1,
     1       enarad,1.d0,1.d0)
      endif
 9000 if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      include 'inc/TEXIT.inc'
      return
      end

      subroutine tbdrift(np,x,px,y,py,z,dv,sx,sz,al,phi0)
      use tfstk
      use ffs_flag, only:calpol
      use mathfun
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),
     $     sx(np),sz(np),
     $     al,phi0,cp,sp,rho0,dx,xi,pzi,pzf,dl,dcp,th
      th=tan(.5d0*phi0)
      sp=2.d0*th/(1.d0+th**2)
      dcp=th*sp
      cp=1.d0-dcp
c      cp=cos(phi0)
c      sp=sin(phi0)
c      if(cp .ge. 0.d0)then
c        dcp=sp**2/(1.d0+cp)
c      else
c        dcp=1.d0-cp
c      endif
      rho0=al/phi0
      call tdrift_free(np,x,px,y,py,z,dv,rho0*sp)
      dx=rho0*dcp
      dl=rho0*xsin(phi0)
      if(calpol)then
        do i=1,np
          xi=x(i)+dx
          pzi=1.d0+pxy2dpz(px(i),py(i))
          pzf=pzi*cp-px(i)*sp
          x(i)=xi*pzi/pzf
          y(i)=y(i)+xi*sp*py(i)/pzf
          z(i)=z(i)-xi*sp/pzf+(1.d0-dv(i))*dl
          px(i)=px(i)*cp+pzi*sp
          sx(i)= cp*sx(i)+sp*sz(i)
          sz(i)=(sz(i)-sp*sx(i))/cp
        enddo
      else
        do i=1,np
          xi=x(i)+dx
          pzi=1.d0+pxy2dpz(px(i),py(i))
          pzf=pzi*cp-px(i)*sp
          x(i)=xi*pzi/pzf
          y(i)=y(i)+xi*sp*py(i)/pzf
          z(i)=z(i)-xi*sp/pzf+(1.d0-dv(i))*dl
          px(i)=px(i)*cp+pzi*sp
        enddo
      endif
      return
      end

      subroutine tbthin(np,x,px,y,py,z,g,sx,sy,sz,phib,phi0,dx,dy,
     1                 theta,dtheta,cost,sint)
      use tbendcom, only:tbrot
      implicit none
      integer*4 np,i
      real*8 phib,phi0,dx,dy,theta,cost,sint,xi,pxi,dtheta
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),sx(np),sy(np),sz(np)
      include 'inc/TENT.inc'
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      do 10 i=1,np
c        px(i)=px(i)+phi0-phib/(1.d0+g(i))**2
        px(i)=px(i)+phi0-phib/(1.d0+g(i))
        z(i)=z(i)-x(i)*phi0
10    continue
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      include 'inc/TEXIT.inc'
      return
      end

      subroutine tbendr(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,
     $     al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     mfring,fringe,ndiv)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      use photontable
      implicit none
      integer*4 np,mfring,i,ndiv,mfr1,mfr2,l
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      real*8 al,phib,phi0,cosp1,sinp1,cosp2,sinp2,rho0,
     $     psi1,psi2,wn1,wn2,wnc,aln,phibn,phi0n,
     $     coswn1,sinwn1,sqwhn1,sinwp1n1,
     $     coswnc,sinwnc,sqwhnc,sinwp1nc,
     $     coswn2,sinwn2,sqwhn2,sinwp1n2
      logical*4 fringe,ph
      rho0=al/phi0
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
      call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,
     $     aln,phi0n,
     1     cosp1,sinp1,1.d0,0.d0,
     1     mfr1,fringe,
     $     coswn1,sinwn1,sqwhn1,sinwp1n1,
     1     .true.,1.d0,0.d0)
      if(photons)then
        call tsetphotongeo(aln,phi0n,0.d0,.false.)
      endif
      do i=2,ndiv-1
        call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       px0,py0,zr0,bsi,
     $       aln,phi0n,
     1       1.d0,0.d0,1.d0,0.d0,
     1       0,.false.,
     $       coswnc,sinwnc,sqwhnc,sinwp1nc,
     1       .true.,0.d0,0.d0)
        if(photons)then
          call tsetphotongeo(aln,phi0n,0.d0,.false.)
        endif
      enddo
      call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,
     $     aln,phi0n,
     1     1.d0,0.d0,cosp2,sinp2,
     1     mfr2,fringe,
     $     coswn2,sinwn2,sqwhn2,sinwp1n2,
     1     .true.,0.d0,1.d0)
      return
      end

      subroutine tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,px0,py0,zr0,bsi,
     $     al,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     enarad,bsi1,bsi2)
      use tfstk
      use ffs_flag
      use tmacro
      use multa, only:nmult
      use tbendcom
      use tspin
      use photontable,only:pp
      use mathfun
      implicit none
      integer*4 np,mfring,i
      real*8 al,phi0,cosp1,sinp1,cosp2,sinp2,
     $     cosw,sinw,sqwh,sinwp1,drhob,dp,p,
     $     pinv,rhoe,pxi,pyi,dpzi,pzi,sp1,x1,dz1,y1,z1,px1,
     $     py1,pv1sqi,f,ff,x2,py2,z2,dph2,ph2,dpx2,pz2,drho,
     $     t2,dpx3,px3,dpz3,pz3,t3,x3,da,y3,z3,pv2sqi,x4,py4,z4,dpz4,
     $     dz4,dxfr1,dyfr1,dzfr1,dxfr2,dyfr2,dzfr2,dpz32,
     $     dyfra1,dyfra2,fa,t4,dpx3a,t2t3,dcosp,px1px3,
     $     phi0a,bsi1,bsi2
      real*8, parameter :: smin=1.d-4
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      logical*4 enarad,fringe
c      write(*,*)'tbendcore-1 ',z(1:3)
      if((mfring .gt. 0 .or. mfring .eq. -1) .and. fb1 .ne. 0.d0)then
        dxfr1=fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        dzfr1=dxfr1*sinp1
        if(fringe)then
          dyfra1=4.d0*dyfr1/fb1**2
        else
          dyfra1=0.d0
        endif
      else
        dxfr1=0.d0
        dyfr1=0.d0
        dyfra1=0.d0
        dzfr1=0.d0
      endif
      if((mfring .gt. 0 .or. mfring .eq. -2) .and. fb2 .ne. 0.d0)then
        dxfr2=fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        dzfr2=dxfr2*sinp2
        if(fringe)then
          dyfra2=4.d0*dyfr2/fb2**2
        else
          dyfra2=0.d0
        endif
      else
        dxfr2=0.d0
        dyfr2=0.d0
        dyfra2=0.d0
        dzfr2=0.d0
      endif
      if(cosp1*cosp2 .gt. 0.d0)then
        dcosp=(sinp2-sinp1)*(sinp2+sinp1)/(cosp1+cosp2)
      else
        dcosp=cosp1-cosp2
      endif
      drhob=rhob-rho0
      do 100 i=1,np
        bsi(i)=bsi(i)+bsi1*y(i)/rhob
        dp=g(i)
        p=1.d0+dp
        pinv=1.d0/p
        rhoe=rhob*p
        pxi=px(i)
        pyi=py(i)
        dpzi=pxy2dpz(pxi,pyi)
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
        ff=.25d0*(f+fa)*y1*pv1sqi
        x2=x1+ff
        py2=py1-px1*f
        z2=z1-px1*ff
        dph2=pxy2dpz(0.d0,py2)
        ph2=1.d0+dph2
        dpx2=pxi*cosp1+(dpzi-dph2)*sinp1
        pz2=1.d0+pxy2dpz(px1,py2)
        drho=drhob+rhoe*dph2+rhob*dp
        t2=(px1+ph2*sinp1)/(pz2+ph2*cosp1)
        dpx3a=(x2*sinw-drho*(sinp2+sinwp1))/rhoe
        dpx3=dpx3a-dpx2*(cosw-sinw*t2)
        px3=ph2*sinp2+dpx3
        px1px3=ph2*(sinp2+sinp1)+dpx3a+dpx2*(sqwh+sinw*t2)
        dpz3=pxy2dpz(px3,py2)
        pz3=1.d0+dpz3
        dpz32=px1px3*(px1-px3)/(pz2+pz3)
        t3=(px3+ph2*sinp2)/(pz3+ph2*cosp2)
        t2t3=(ph2*sinp2+px1px3)/(pz3+ph2*cosp2)
     $       +ph2*sinp1/(pz2+ph2*cosp1)
     $       +px1*(dpz32-ph2*dcosp)
     $       /(pz3+ph2*cosp2)/(pz2+ph2*cosp1)
c        write(*,*)t2,t3,t2+t3,t2t3,px1+px3,px1px3
        t4=(cosp2+t3*sinp2)*(pz2*cosp1+px1*sinp1)
        x3=x2*(cosw-rho0/rhoe*t3*sinw/ph2)
     1       +(rho0*(cosw*t2t3+sinw*(1.d0-t2*t3))*dpx2-
     1       drho*(-(sinp2+sinwp1)*rho0/rhoe*t3
     $       -dpz32-sqwh*pz2-sinw*px1))/ph2
        da=asin(min(1.d0,max(-1.d0,
     $       (dpx2*(
     $       sinp1*(t2*(pz3*cosp2+px3*sinp2)-px1*cosp2)
     $       -sinp2*(t3*(pz2*cosp1+px1*sinp1)-px3*cosp1)
     $       +cosp1*cosp2*dpz32
     $       +t4*(sqwh+sinw*t2))
     1       +dpx3a*t4)/ph2**2)))
        phi0a=phi0+da
        y3=y1+py2*rhoe*phi0a
        z3=z2-phi0*(dp*rhob+drhob)-da*rhoe-dv(i)*al
        pv2sqi=1.d0/max(smin,1.d0-px3**2)
        fa=y3/rhoe*sqrt(pv2sqi)
        f=(1.d0-(y3/rhob)**2/6.d0)*fa
        ff=.25d0*(f+fa)*y3*pv2sqi
        py4=py2-px3*f
        z4=z3-px3*ff
        x4=x3-ff-dxfr2*dp*pinv
        py4=py4+(dyfr2-dyfra2*y3**2)*y3*pinv**2
        z4=z4+(dxfr2*px3+
     $       (.5d0*dyfr2-.25d0*dyfra2*y3**2)*y3**2*pinv)*pinv-dzfr2
        dpz4=pxy2dpz(px3,py4)
        px(i)=-cosp2*dpx3+sinp2*(dpz4-dpz3-dpx3*t3)
        dz4=x4*sinp2/(cosp2*(1.d0+dpz4)+sinp2*px3)
        x(i)=x4*cosp2+px(i)*dz4
        py(i)=py4
        y(i)=y3+py4*dz4
        z(i)=z4-dz4
        bsi(i)=bsi(i)-bsi2*y(i)/rhob
100   continue
      if(rad .and. enarad)then
        write(*,*)'tbendcore-tradk ',pp%l,pp%al,pp%phi
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       px0,py0,zr0,cos(phi0),sin(phi0),bsi,al)
        px0=px
        py0=py
        zr0=z
        bsi=0.d0
      endif
      return
      end
