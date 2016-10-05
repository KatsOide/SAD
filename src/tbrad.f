      subroutine tbrad(np,x,px,y,py,z,g,dv,pz,
     $     l,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1     fs1,fs2,mfring,fringe,eps0)
      use tfstk, only: sqrtl
      use ffs_flag
      use tmacro
      use bendeb, only:epsbend
      implicit none
c      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
c     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
c     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      integer*4 np,l,i,mfring,ndiv,nx,ngamma,n
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np)
      real*8 al,phib,phi0,tanp1,tanp2,ak,dx,dy,theta,dphix,dphiy,
     $     cost,sint,fs1,fs2,eps0,bxa,bya,alr,dpradx,dprady,
     $     alsum,brad0,brad1,brad2,rhob,rho0,af,aind,dxfr1,
     $     dyfr1,dyfra1,p,f,eps,fpx,ff,ur,an,phin,cs00,sn00,
     $     akk0,by0,byx,drhob,sp,ak1,dp,xi,yi,pxi,pyi,s,dpz1,xr,
     $     brad,rho,alx,prob,pz1,dpx,pxf,dpz2,pz2,d,
     $     sinda,da,al0,cosp1,sinp1,cosp2,sinp2,akn,aln,
     $     bx00,by00,csphi0,dprad,drho,dxfr2,p1,phix,pr,
     $     rhoe,sinsq0,snphi0,sq00,tran,dyfr2,dyfra2,h
      logical*4 fringe
c     begin initialize for preventing compiler warning
      brad0=0.d0
      brad1=0.d0
      brad2=0.d0
c     end  initialize for preventing compiler warning
      include 'inc/TENT.inc'
      if(dphiy .ne. 0.d0)then
        do 3510 i=1,np
c          pr=(1.d0+g(i))**2
          pr=(1.d0+g(i))
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
3510    continue
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      rho0=al/phi0
      aind=rho0/phi0*ak
      if(fs1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fs1**2/rhob/24.d0
          dyfr1=fs1/rhob**2/6.d0
          dyfra1=4.d0*dyfr1/fs1**2
          do i=1,np
c            dp=g(i)*(2.d0+g(i))
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)+dxfr1*dp/pr
            py(i)=py(i)+(dyfr1-dyfra1*y(i)**2)*y(i)/pr**2
            z(i)=z(i)+(dxfr1*px(i)+
     $           (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
      if(fringe)then
        af=1.d0
      else
        af=0.d0
      endif
      do 1100 i=1,np
c        g(i)=g(i)*(2.d0+g(i))
        p=1.d0+g(i)
        rhoe=rhob*p
        f=y(i)/rhoe
        fpx=af*px(i)
        py(i)=py(i)-(tanp1+fpx)*f
        ff  =af*y(i)*f*.5d0
        x(i)=x(i)+ff
        z(i)=z(i)-ff*fpx
        px(i)=px(i)+tanp1*x(i)/rhoe
1100  continue
      if(eps0 .le. 0.d0)then
        eps=epsbend
      else
        eps=epsbend*eps0
      endif
      ur=urad*p0**3
      an=anrad*p0
      ndiv=1+int(sqrt(abs(ak*al)/eps/12.d0))
      nx=int(al/abs(rhob)*an/.18d0)+1
      ndiv=max(ndiv,nx)
      if(radlight)then
        ngamma=int(h0*abs(phib)*1.d-2*anrad/eps)+1
        ndiv=max(ndiv,ngamma)
      endif
      phin=phi0/ndiv
      cs00=cos(phin)
      sn00=sin(phin)
      if(cs00 .ge. 0.d0)then
        sq00=sn00**2/(1.d0+cs00)
      else
        sq00=1.d0-cs00
      endif
c      sq00=2.d0*sin(phin*.5d0)**2
      akn=ak/ndiv
      aln=al/ndiv
      akk0=ak/al
      by0=brhoz/rhob
      byx=brhoz*akk0
      drhob=rhob-rho0
      sp=0.d0
      do 1120 n=1,ndiv
        if(n .eq. 1)then
          ak1=akn*.5d0
        else
          ak1=akn
        endif
        do 1130 i=1,np
          dp=g(i)
          p=1.d0+dp
          xi=x(i)
          yi=y(i)
          if(byx .eq. 0.d0)then
            pxi=px(i)
            pyi=py(i)
            if(n .eq. 1)then
              s=pxi**2+pyi**2
              dpz1=-s/(1.d0+sqrtl(1.d0-s))
            else
              dpz1=pz(i)
            endif
          else
            xr=xi/rho0
            pxi=px(i)-ak1*(xi+xi*xr*(.5d0-xr*(2.d0-xr)/12.d0))/p
            pyi=py(i)+ak1*yi/p
            s=min(.95d0,pxi**2+pyi**2)
            dpz1=-s/(1.d0+sqrtl(1.d0-s))
          endif
          al0=aln
          alsum=aln*(n-1)
110       if(byx .eq. 0.d0)then
            brad=by0**2
          else
            bx00=byx*yi
            by00=by0+byx*xi
            brad0=by00**2+bx00**2
            brad1=byx*(by00*pxi+bx00*pyi)*.5d0
            brad2=byx*(byx*(pxi**2+pyi**2)
     1                -akk0*(by00*xi-bx00*yi))/12.d0
            brad=brad0+al0*(brad1+al0*brad2)
          endif
          if(brad .ne. 0.d0)then
            if(byx .eq. 0.d0)then
              rho=abs(rhob)
              alx=min(al0,.19d0*rho/an)
 111          alr=alx*(1.d0+xi/rho0+(pxi**2+pyi**2)*.5d0)
              prob=.5d0*alr*an/rho
c              write(*,*)'tbrad ',n,i,prob,alx,byx,an
              if(prob .gt. .1d0)then
                alx=alx*(.095d0/prob)
                go to 111
              endif
            else
              rho=brho/sqrt(brad)
              alx=min(al0,.19d0*rho/an)
112           if(alx .ne. al0)then
                brad=brad0+alx*(brad1+alx*brad2)
                rho=brho/sqrt(brad)
              endif
              alr=alx*(1.d0+xi/rho0+(pxi**2+pyi**2)*.5d0)*.5d0
              prob=alr*an/rho
              if(prob .gt. .1d0)then
                alx=alx*(.095d0/prob)
                go to 112
              endif
            endif
            if(tran() .gt. 1.d0-prob)then
              bxa=byx*(yi*(1.d0+akk0*alx**2/6.d0)+pyi*alx*.5d0)
              bya=by0+byx*(xi*(1.d0-akk0*alx**2/6.d0)+pxi*alx*.5d0)
              call tsynchrad(p,alr,bxa,bya,
     $             dprad,dpradx,dprady,
     $             i,l,alsum,0.d0,phi0*alsum/al,theta,
     $             xi,yi,pxi,pyi)
              dp=max(-.999d0,dp-dprad)
              p=1.d0+dp
              sp=sp-dprad
              pxi=pxi-dpradx
              pyi=pyi-dprady
            endif
          else
            alx=al0
          endif
          alsum=alsum+alx
          if(alx .eq. aln)then
            phix=phin
            csphi0=cs00
            snphi0=sn00
            sinsq0=sq00
          else
c           write(*,'(1p2g12.4,2i5)')aln,alx,ndiv,n
            phix=phi0*alx/al
            csphi0=cos(phix)
            snphi0=sin(phix)
            if(csphi0 .ge. 0.d0)then
              sinsq0=snphi0**2/(1.d0+csphi0)
            else
              sinsq0=1.d0-csphi0
            endif
c            sinsq0=2.d0*sin(phix*.5d0)**2
          endif
          rhoe=rhob*p
          pz1=1.d0+dpz1
          drho=drhob+rhoe*dpz1+rhob*dp
          dpx=-(xi-drho)/rhoe*snphi0-sinsq0*pxi
          pxf=pxi+dpx
          s=min(.95d0,pxf**2+pyi**2)
          dpz2=-s/(1.d0+sqrtl(1.d0-s))
          pz2=1.d0+dpz2
          d=pxf*pz1+pxi*pz2
          if(d .eq. 0.d0)then
            sinda=2.d0*pxf*pz2/(pxf**2+pz2**2)
          else
            sinda=dpx*(pxf+pxi)/d
          endif
          da=asin(sinda)
          xi=xi*csphi0+rhoe*(snphi0*pxi-dpx*(pxi+pxf)/(pz1+pz2))
     1         +drho*sinsq0
          pxi=pxf
          yi=yi+pyi*rhoe*(phix-da)
          z(i)=z(i)-phix*(dp*rhob+drhob)+da*rhoe-dv(i)*alx
          dpz1=dpz2
          if(byx .eq. 0.d0)then
            brad=by0**2
          else
            bx00=byx*yi
            by00=by0+byx*xi
            brad0=by00**2+bx00**2
            brad1=byx*(by00*pxi+bx00*pyi)*.5d0
            brad2=byx*(byx*(pxi**2+pyi**2)
     1                -akk0*(by00*xi-bx00*yi))/12.d0
            brad=brad0+alx*(brad1+alx*brad2)
          endif
          if(brad .ne. 0.d0)then
            if(byx .eq. 0.d0)then
              rho=abs(rhob)
            else
              rho=brho/sqrt(brad)
            endif
            alr=alx*(1.d0+xi/rho0+(pxi**2+pyi**2)*.5d0)*.5d0
            prob=alr*an/rho
            if(tran() .gt. 1.d0-prob)then
              bxa=byx*(yi*(1.d0+akk0*alx**2/6.d0)+pyi*alx*.5d0)
              bya=by0+byx*(xi*(1.d0-akk0*alx**2/6.d0)+pxi*alx*.5d0)
              call tsynchrad(p,alr,bxa,bya,
     $             dprad,dpradx,dprady,
     $             i,l,alsum,-alr,phi0*alsum/al,theta,
     $             xi,yi,pxi,pyi)
              dp=max(-.999d0,dp-dprad)
              p=1.d0+dp
              sp=sp-dprad
              pxi=pxi-dpradx
              pyi=pyi-dprady
            endif
          endif
          al0=al0-alx
          if(al0 .gt. 0.d0)then
            go to 110
          endif
          x(i)=xi
          px(i)=pxi
          y(i)=yi
          py(i)=pyi
          pz(i)=dpz1
          g(i)=dp
1130    continue
        if(radlight)then
          call tlstore(np,x,y,z,dv,theta,0.d0,phin*n,rho0,
     $         p0/h0*c,dvfs,n .eq. ndiv)
        endif
1120  continue
      if(radcod)then
        sp=0.d0
      else
        sp=sp/np
      endif
      do 1140 i=1,np
        if(akn .ne. 0.d0)then
          pr=1.d0+g(i)
          xr=x(i)/rho0
          px(i)=px(i)-
     1          akn*(x(i)+x(i)*xr*(.5d0-xr*(2.d0-xr)/12.d0))/pr*.5d0
          py(i)=py(i)+akn*y(i)/pr*.5d0
        endif
        dp=max(-.99d0,g(i)-sp)
        pr=1.d0+dp
c        g(i)=dp/(1.d0+sqrt(pr))
        g(i)=dp
        p1=pr*p0
        h=p1*sqrt(1.d0+1.d0/p1**2)
        dv(i)=-dp*(1.d0+pr)/h/(h+pr*h0)+dvfs
1140  continue
      do 1150 i=1,np
c        p=(1.d0+g(i))**2
        p=(1.d0+g(i))
        rhoe=rhob*p
        px(i)=px(i)+tanp2*x(i)/rhoe
        f=y(i)/rhoe
        fpx=af*px(i)
        py(i)=py(i)-(tanp2-fpx)*f
        ff  =af*y(i)*f*.5d0
        x(i)=x(i)-ff
        z(i)=z(i)+ff*fpx
1150  continue
      if(fs2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=fs2**2/rhob/24.d0
          dyfr2=fs2/rhob**2/6.d0
          dyfra2=4.d0*dyfr2/fs2**2
          do i=1,np
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)-dxfr2*dp/pr
            py(i)=py(i)+(dyfr2-dyfra2*y(i)**2)*y(i)/pr**2
            z(i)=z(i)-(dxfr2*px(i)-
     $           (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
      if(dphiy .ne. 0.d0)then
        do 3520 i=1,np
c          pr=(1.d0+g(i))**2
          pr=(1.d0+g(i))
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
3520    continue
      endif
      include 'inc/TEXIT.inc'
      return
      end

      subroutine tlinit(np,h0,geo)
      use tfmem, only:ktaloc
      implicit none
      integer*4 np
      real*8 h0,geo(3,4)
      integer*8 katbl
      integer*4 np0,ltbl,lpoint
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
      np0=np
      ltbl=int(h0*10.d0)+1
      katbl=ktaloc(np*4*ltbl)
      lpoint=0
      gx=geo(1,4)
      gy=geo(2,4)
      gz=geo(3,4)
      gt=0.d0
      geo0(:,3)=geo(:,3)
      geo0(1,4)=gx
      geo0(2,4)=gy
      geo0(3,4)=gz
      t0=gt
      return
      end

      subroutine tlstore(np,x,y,z,dv,theta,dl,dphi,rho0,v0,dvfs,keep)
      use tfstk
      implicit none
      integer*8 katbl,new,kai
      integer*4 np,np0,ltbl,lpoint,m,i
      real*8 x(np),y(np),z(np),dv(np),theta,dphi,rho0,v0,
     $     x1,x2,x3,y1,y2,y3,sp0,cp0,r1,r2,
     $     geoi(3,3),cost,sint,dl
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt,dvfs
      logical*4 keep
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
c     begin initialize for preventing compiler warning
      cost=1.d0
      sint=0.d0
      y1=0.d0
      y2=0.d0
      y3=0.d0
c     end   initialize for preventing compiler warning
c      write(*,*)'tlstore ',lpoint,dphi,dl,theta,rho0,v0
      if(dphi .eq. 0.d0)then
        gx=geo0(1,4)+dl*geo0(1,3)
        gy=geo0(2,4)+dl*geo0(2,3)
        gz=geo0(3,4)+dl*geo0(3,3)
        geoi=geo0(:,1:3)
        gt=t0+dl/v0
      else
        if(theta .ne. 0.d0)then
          cost=cos(theta)
          sint=sin(theta)
          x1= cost*geo0(1,1)-sint*geo0(1,2)
          x2= cost*geo0(2,1)-sint*geo0(2,2)
          x3= cost*geo0(3,1)-sint*geo0(3,2)
          y1= sint*geo0(1,1)+cost*geo0(1,2)
          y2= sint*geo0(2,1)+cost*geo0(2,2)
          y3= sint*geo0(3,1)+cost*geo0(3,2)
        else
          x1= geo0(1,1)
          x2= geo0(2,1)
          x3= geo0(3,1)
        endif
        sp0=sin(dphi)
        cp0=cos(dphi)
        r1=rho0*sp0
        if(cp0 .ge. 0.d0)then
          r2=rho0*sp0**2/(1.d0+cp0)
        else
          r2=rho0*(1.d0-cp0)
        endif
        gx=geo0(1,4)+(r1*geo0(1,3)-r2*x1)
        gy=geo0(2,4)+(r1*geo0(2,3)-r2*x2)
        gz=geo0(3,4)+(r1*geo0(3,3)-r2*x3)
        geoi(1,1)= cp0*x1+sp0*geo0(1,3)
        geoi(1,3)=-sp0*x1+cp0*geo0(1,3)
        geoi(2,1)= cp0*x2+sp0*geo0(2,3)
        geoi(2,3)=-sp0*x2+cp0*geo0(2,3)
        geoi(3,1)= cp0*x3+sp0*geo0(3,3)
        geoi(3,3)=-sp0*x3+cp0*geo0(3,3)
        geoi(1,2)=geo0(1,2)
        geoi(2,2)=geo0(2,2)
        geoi(3,2)=geo0(3,2)
        gt=t0+abs(rho0*dphi)/v0
      endif
      lpoint=lpoint+1
      if(lpoint .gt. ltbl)then
c        write(*,*)'tlspect ',ltbl,np0,ltbl*8*np0
        new=ktaloc(ltbl*8*np0)
        if(new .gt. 0)then
          call tmov(rlist(katbl),rlist(new),ltbl*np0)
          call tmov(rlist(katbl+ltbl),rlist(new+ltbl*2),ltbl*np0)
          call tmov(rlist(katbl+ltbl*2),rlist(new+ltbl*4),ltbl*np0)
          call tmov(rlist(katbl+ltbl*3),rlist(new+ltbl*6),ltbl*np0)
          ltbl=ltbl*2
          call tfree(katbl)
          katbl=new
        else
          write(*,*)'Too long ltbl',ltbl
          call abort
        endif
      endif
      m=ltbl*4
      do i=1,np
        kai=(i-1)*4*ltbl+katbl+lpoint-1
        rlist(kai       )=gt-z(i)/(v0*(1.d0-dv(i)+dvfs))
        rlist(kai+ltbl  )=gx+x(i)*geoi(1,1)+y(i)*geoi(1,2)
        rlist(kai+ltbl*2)=gy+x(i)*geoi(2,1)+y(i)*geoi(2,2)
        rlist(kai+ltbl*3)=gz+x(i)*geoi(3,1)+y(i)*geoi(3,2)
      enddo
      if(keep)then
        geo0(:,3)=geoi(:,3)
        if(dphi .ne. 0 .and. theta .ne. 0.d0)then
          geo0(1,2)=-sint*geo0(1,1)+cost*y1
          geo0(1,1)= cost*geo0(1,1)+sint*y1
          geo0(2,2)=-sint*geo0(2,1)+cost*y2
          geo0(2,1)= cost*geo0(2,1)+sint*y2
          geo0(3,2)=-sint*geo0(3,1)+cost*y3
          geo0(3,1)= cost*geo0(3,1)+sint*y3
        endif
        geo0(1,4)=gx
        geo0(2,4)=gy
        geo0(3,4)=gz
        t0=gt
      endif
      return
      end

      subroutine tlresult(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: klx,klij,kli
      integer*8 katbl,kati
      integer*4 np0,ltbl,lpoint,i,j
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
c      write(*,*)ltbl,lpoint,katbl,ilist(1,katbl-1)
      kx=kxadaloc(-1,np0,klx)
      do i=1,np0
        kati=katbl+ltbl*4*(i-1)
        klx%dbody(i)=kxadaloc(0,4,kli)
        do j=1,4
          kli%dbody(j)=kxavaloc(0,lpoint,klij)
          call tmov(rlist(kati+(j-1)*ltbl),klij%rbody(1),lpoint)
        enddo
      enddo
      call tfree(katbl)
      return
      end

      subroutine tlfield(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_list), pointer :: klt
      integer*8 kx, ke
      integer*4 isp1,irtc,n,m,itfmessage
      real*8 tx,ty,tz
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      call tfmsize(ktastk(isp1+1),n,m,irtc)
      if(irtc .ne. 0)then
        return
      elseif(n .ne. 4)then
        irtc=itfmessage(9,'General::wrongleng','"#1","4"')
        return
      endif
      if(.not. tfreallistq(ktastk(isp),klt))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of Reals for #2"')
        return
      endif
      if(klt%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongleng','"#2","3"')
        return
      endif
      tx=klt%rbody(1)
      ty=klt%rbody(2)
      tz=klt%rbody(3)
      call tlcalc(ktastk(isp1+1),ke,m,tx,ty,tz)
      kx=ktflist+ke
      return
      end

      subroutine tlcalc(ka,ke,lpoint,tx,ty,tz)
      use tfstk
      use tmacro
      implicit none
      integer*8 ka,ke
      integer*4 lpoint,i,j
      real*8 xtbl(lpoint,4),f(lpoint,10),
     $     tx,ty,tz,dxi,dyi,dzi,ri,r0,dri,coeff,dt,f1,z0
      call tfl2m(klist(ktfaddr(ka)-3),xtbl,lpoint,4,.true.)
      z0=4.d-7*2.d0*asin(1.d0)*c
      r0=sqrt(tx**2+ty**2+tz**2)
      do i=1,lpoint
        dxi=tx-xtbl(i,2)
        dyi=ty-xtbl(i,3)
        dzi=tz-xtbl(i,4)
        ri=sqrt(dxi**2+dyi**2+dzi**2)
        dri=-(xtbl(i,2)*(tx+dxi)+xtbl(i,3)*(ty+dyi)
     $       +xtbl(i,4)*(tz+dzi))/(r0+ri)
        xtbl(i,1)=xtbl(i,1)+dri/c
        xtbl(i,2)=dxi/ri
        xtbl(i,3)=dyi/ri
        xtbl(i,4)=dzi/ri
        f(i,8)=ri/c
      enddo
      do i=2,4
        call spline(lpoint,xtbl(1,1),xtbl(1,i),f(1,i),f(1,1),0,0.d0)
      enddo
      do i=5,7
        do j=1,lpoint
          f(j,10)=xtbl(j,i-3)/f(j,8)**2
        enddo
        call spline(lpoint,xtbl(1,1),f(1,10),f(1,i),f(1,1),0,0.d0)
        f1=f(lpoint-1,i)
        do j=1,lpoint-1
          dt=xtbl(j+1,1)-xtbl(j,1)
          f(j,i)=f(j,10)+f(j,8)*
     $         ((f(j+1,10)-f(j,10))/dt-dt*(2.d0*f(j,i)+f(j+1,i)))
        enddo
        dt=xtbl(lpoint,1)-xtbl(lpoint-1,1)
        f(lpoint,i)=f(lpoint,10)+f(lpoint,8)*
     $       ((f(lpoint,10)-f(lpoint-1,10))/dt
     $       +dt*(f1+2.d0*f(lpoint,i)))
      enddo
      coeff=1.d-7*e*charge
      do i=1,lpoint
        f(i,1)=xtbl(i,1)
        f(i,2)=coeff*(6.d0*f(i,2)+f(i,5))
        f(i,3)=coeff*(6.d0*f(i,3)+f(i,6))
        f(i,4)=coeff*(6.d0*f(i,4)+f(i,7))
        f(i,5) =(xtbl(i,3)*f(i,4)-xtbl(i,4)*f(i,3))/z0
        f(i,6) =(xtbl(i,4)*f(i,2)-xtbl(i,2)*f(i,4))/z0
        f(i,7) =(xtbl(i,2)*f(i,3)-xtbl(i,3)*f(i,2))/z0
        f(i,8) =f(i,3)*f(i,7)-f(i,4)*f(i,6)
        f(i,9) =f(i,4)*f(i,5)-f(i,2)*f(i,7)
        f(i,10)=f(i,2)*f(i,6)-f(i,3)*f(i,5)
      enddo
      ke=ktfaddr(kxm2l(f,lpoint,10,lpoint,.true.))
      return
      end

      subroutine tlspect(isp1,kx,nparallel,irtc)
      use tfstk
      use macphys
      implicit none
      type (sad_list), pointer :: klp,kl1,kl2,klt,klx
      real*8 speedoflight
      parameter (speedoflight=cveloc)
      integer*8 kx,k1,k2,kac,kas,kal,kack,kad
      integer*4 isp1,irtc,m,na,i,k,
     $     na1,nastep,nparallel,ipr,fork_worker,getpid,j,
     $     wait,irw,naa,isw,ipid,itfmessage
      real*8 al0,al1,dal,w,c,s,a0,a1,x,ak,ak0,ak1,dak
      integer*4 ichpid(100)
      character*64 fn
c     begin initialize for preventing compiler warning
      ipid=0
      naa=0
      kal=0
      kas=0
c     end   initialize for preventing compiler warning
      irtc=0
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. tflistqk(ktastk(isp1+1),klp))then
        irtc=itfmessage(9,'General::wrongtype','"List for #1"')
        return
      endif
      if(klp%nl .ne. 2)then
        irtc=itfmessage(9,'General::wrongleng','"#1","2"')
        return
      endif
      k1=klp%body(1)
      k2=klp%body(2)
      if(.not. tfreallistq(k1,kl1) .or. .not. tfreallistq(k2,kl2))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of List of Reals for #1"')
        return
      endif
      m=kl1%nl
      if(m .ne. kl2%nl)then
        irtc=itfmessage(9,'General::equalleng','"#1[[1]] and #[[2]]"')
        return
      endif
      if(.not. tfreallistq(ktastk(isp),klt))then
        irtc=itfmessage(9,'General::wrongtype','"List of Reals for #2"')
        return
      endif
      if(klt%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongleng','"#","3"')
        return
      endif
      al0=klt%rbody(1)
      al1=klt%rbody(2)
      dal=klt%rbody(3)
      na=int((al1-al0)/dal+2.5d0)
      ak0=2.d0*pi/al1
      ak1=2.d0*pi/al0
      dak=(ak1-ak0)/(na-1)
      na1=1
      if(nparallel .gt. 1 .and. na .ge. nparallel)then
        ipid=getpid()
        nastep=nparallel
        ipr=1
        do while(na1 .lt. nparallel .and. ipr .gt. 0)
          ipr=fork_worker()
          if(ipr .gt. 0)then
            ichpid(na1)=ipr
            na1=na1+1
          endif
        enddo
      else
        ipr=1
        nastep=1
      endif
      if(ipr .eq. 0)then
        naa=(na-na1)/nastep+1
        kac=ktaloc(naa*3)
      else
        kac=ktavaloc(0,na)
        kas=ktavaloc(0,na)
        kal=ktavaloc(0,na)
      endif
      do k=na1,na,nastep
        ak=ak0+dak*(k-1)
        w=ak*speedoflight
        c=0.d0
        s=0.d0
        do i=2,m-1
          x=w*kl1%rbody(i)
          a0=(kl2%rbody(i)-kl2%rbody(i-1))/
     $         (kl1%rbody(i)-kl1%rbody(i-1))
          a1=(kl2%rbody(i+1)-kl2%rbody(i))/
     $         (kl1%rbody(i+1)-kl1%rbody(i))
          c=c+cos(x)*(a0-a1)
          s=s+sin(x)*(a0-a1)
        enddo
        a0=(kl2%rbody(m)-kl2%rbody(m-1))/
     $       (kl1%rbody(m)-kl1%rbody(m-1))
        x=w*kl1%rbody(m)
        c=c/w+cos(x)*a0/w+kl2%rbody(m)*sin(x)
        s=s/w+sin(x)*a0/w-kl2%rbody(m)*cos(x)
        a1=(kl2%rbody(2)-kl2%rbody(1))/
     $       (kl1%rbody(2)-kl1%rbody(1))
        x=w*kl1%rbody(1)
        c=(c-cos(x)*a1/w-kl2%rbody(1)*sin(x))/w
        s=(s-sin(x)*a1/w+kl2%rbody(1)*cos(x))/w
        if(ipr .eq. 0)then
          kack=kac+(k-na1)/nastep*3
          rlist(kack  )=ak
          rlist(kack+1)=c
          rlist(kack+2)=s
        else
          rlist(kal+k)=ak
          rlist(kac+k)=c
          rlist(kas+k)=s
        endif
      enddo
      if(nparallel .gt. 1 .and. na .ge. nparallel)then
        if(ipr .eq. 0)then
          write(fn,'(a,i12.12)')'/tmp/',(ipid*nparallel)+na1
          open(999,file=fn,form='unformatted',status='unknown')
          write(999)(rlist(kac+i),i=0,naa*3-1)
          stop
        else
          kad=ktaloc(((na-1)/nastep+1)*3)
          do j=1,nparallel-1
 3001       irw=wait(isw)
            do k=1,nparallel-1
              if(irw .eq. ichpid(k))then
                na1=k
                ichpid(k)=0
                go to 3010
              endif
            enddo
            write(*,*)'RadiationSpectrum Strange child: ',irw,isw
            go to 3001
 3010       write(fn,'(a,i12.12)')'/tmp/',(ipid*nparallel)+na1
            open(999,file=fn,form='unformatted',status='old')
            naa=(na-na1)/nastep+1
            read(999)(rlist(kad+i),i=0,naa*3-1)
            close(999,status='DELETE')
            do k=na1,na,nastep
              kack=kad+(k-na1)/nastep*3
              rlist(kal+k)=rlist(kack  )
              rlist(kac+k)=rlist(kack+1)
              rlist(kas+k)=rlist(kack+2)
            enddo
          enddo
          call tfree(kad)
        endif
      endif
      kx=ktflist+ktadaloc(-1,3,klx)
      klx%body(1)=ktflist+kal
      klx%body(2)=ktflist+kac
      klx%body(3)=ktflist+kas
      return
      end

      subroutine tsynchrad(p,alr,bxa,bya,
     $     dprad,dpradx,dprady,
     $     mp,l,al,s0,phi,theta,xi,yi,pxi,pyi)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 irtn,ng,l,mp
      real*8 p,alr,bxa,bya,dprad,dpradx,dprady,tran,
     $     b0,e0,eg,thx,thy,xi30,thu,thv,
     $     c1,c2,s1,s2,xi3,xi2,xi1,s0,phi1,
     $     thx1,thy1,ds,thu1,xi,yi,pxi,pyi,al,al1,
     $     phi,theta,gx,gy,gz,dpgx,dpgy,dpgz
      e0=amass*sqrt((p*p0)**2+1.d0)
      b0=sqrt(bxa**2+bya**2)
      call SYNRADCL(E0,B0,alr,2,NG,EG,thu,thv,XI30,XI2,IRTN)
      if(irtn .gt. 100 .or. ng .eq. 0)then
        dprad=0.d0
        dpradx=0.d0
        dprady=0.d0
      else
        dprad=eg/amass/p0
        c1=bya/b0
        s1=bxa/b0
        thx=thu*c1+thv*s1
        thy=-thu*s1+thv*c1
        dpradx=thx*dprad
        dprady=thy*dprad
        if(photons)then
          c2=(c1-s1)*(c1+s1)
          s2=2.d0*c1*s1
          xi3=c2*xi30
          xi1=-s2*xi30
          ds=tran()*alr
          thu1=thu-(s0+ds)*b0/brho
          thx1=thu1*c1+thv*s1
          thy1=-thu1*s1+thv*c1
          al1=al+ds+s0
          if(al .eq. 0.d0)then
            phi1=0.d0
          else
            phi1=phi*al1/al
          endif
          call tphotonconv(al1,phi1,theta,
     $         geo(:,:,l),xi,yi,
     $         dprad,pxi-thx1,pyi-thy1,
     $         gx,gy,gz,dpgx,dpgy,dpgz,xi1,xi3)
          call tphotonstore(mp,l,gx,gy,gz,
     $         dpgx*p0,dpgy*p0,dpgz*p0,
     $         xi1,xi2,xi3)
        endif
      endif
      return
      end

      module photontable
      implicit none
      integer*4 ntable,ltable
      parameter (ntable=256,ltable=100000)
      integer*8 kphtable(ntable)
      integer*4 nt,itp,ilp,lt
      end module

      subroutine tphotoninit()
      use photontable
      use tfstk
      use tmacro
      implicit none
      nt=ntable
      lt=ltable
      itp=0
      ilp=0
      kphtable(1)=0
      return
      end

      subroutine tphotonconv(al,phi,theta,geo1,xi,yi,dp,dpx,dpy,
     $     gx,gy,gz,dpgx,dpgy,dpgz,xi1,xi3)
      use tfstk, only: sqrtl
      implicit none
      real*8 al,phi,theta,geo1(3,4),xi,yi,dp,dpx,dpy,gx,gy,gz,
     $     dpz,x1,x2,x3,y1,y2,y3,z1,z2,z3,rho0,sp0,cp0,r1,r2,
     $     dpgx,dpgy,dpgz,cost,sint,xi1,xi3,chi,xi3a
      cost=cos(theta)
      sint=sin(theta)
      x1= cost*geo1(1,1)-sint*geo1(1,2)
      x2= cost*geo1(2,1)-sint*geo1(2,2)
      x3= cost*geo1(3,1)-sint*geo1(3,2)
      y1= sint*geo1(1,1)+cost*geo1(1,2)
      y2= sint*geo1(2,1)+cost*geo1(2,2)
      y3= sint*geo1(3,1)+cost*geo1(3,2)
      if(phi .eq. 0.d0)then
        gx=geo1(1,4)+geo1(1,3)*al
        gy=geo1(2,4)+geo1(2,3)*al
        gz=geo1(3,4)+geo1(3,3)*al
        z1=geo1(1,3)
        z2=geo1(2,3)
        z3=geo1(3,3)
      else
        rho0=al/phi
        sp0=sin(phi)
        cp0=cos(phi)
        r1=rho0*sp0
        if(cp0 .ge. 0.d0)then
          r2=rho0*sp0**2/(1.d0+cp0)
        else
          r2=rho0*(1.d0-cp0)
        endif
        gx=geo1(1,4)+(r1*geo1(1,3)-r2*x1)
        gy=geo1(2,4)+(r1*geo1(2,3)-r2*x2)
        gz=geo1(3,4)+(r1*geo1(3,3)-r2*x3)
        z1=-sp0*x1+cp0*geo1(1,3)
        x1= cp0*x1+sp0*geo1(1,3)
        z2=-sp0*x2+cp0*geo1(2,3)
        x2= cp0*x2+sp0*geo1(2,3)
        z3=-sp0*x3+cp0*geo1(3,3)
        x3= cp0*x3+sp0*geo1(3,3)
      endif
      gx=gx+xi*x1+yi*y1
      gy=gy+xi*x2+yi*y2
      gz=gz+xi*x3+yi*y3
      dpz=dp*sqrtl(1.d0-dpx**2-dpy**2)
      dpgx=dpz*z1+dp*(dpx*x1+dpy*y1)
      dpgy=dpz*z2+dp*(dpx*x2+dpy*y2)
      dpgz=dpz*z3+dp*(dpx*x3+dpy*y3)
      if(x3 .eq. 0.d0)then
        chi=0.d0
      else
        chi=2.d0*atan2(x3,-y3)
      endif
      xi3a=cos(chi)*xi3+sin(chi)*xi1
      xi1=-sin(chi)*xi3+cos(chi)*xi1
      xi3=xi3a
      return
      end

      subroutine tphotonstore(mp,l,gx,gy,gz,dpgx,dpgy,dpgz,
     $     xi1,xi2,xi3)
      use photontable
      use tfstk
      use tmacro
      implicit none
      integer*8 kp
      integer*4 mp,l
      real*8 gx,gy,gz,dpgx,dpgy,dpgz,xi1,xi2,xi3
      if(ilp .eq. 0)then
        itp=itp+1
        kphtable(itp)=ktaloc(10*lt)
        ilp=1
      endif
      kp=kphtable(itp)+(ilp-1)*10
      ilist(1,kp)=mp
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
      return
      end

      subroutine tphotonlist()
      use photontable
      use tfstk
      use tmacro
      implicit none
      type (sad_list), pointer ::kli,klx
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
          klx%body(i)=ktflist+ktavaloc(0,nitem,kli)
          kli%attr=lconstlist
          dp=sqrt(rlist(kp+4)**2+rlist(kp+5)**2
     $         +rlist(kp+6)**2)
          kli%rbody(1)=dp*amass
          kli%rbody(2)=rlist(kp+1)
          kli%rbody(3)=rlist(kp+2)
          kli%rbody(4)=rlist(kp+3)
          kli%rbody(5)=rlist(kp+4)/dp
          kli%rbody(6)=rlist(kp+5)/dp
          kli%rbody(7)=rlist(kp+6)/dp
          kli%rbody(8)=rlist(kp+7)
          kli%rbody(9)=rlist(kp+8)
          kli%rbody(10)=rlist(kp+9)
          kli%rbody(11)=ilist(1,kp)
          kli%rbody(12)=ilist(2,kp)
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
