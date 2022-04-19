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
      if(k /= 0)then
        al0=cmp%value(k)
        if(al0 /= 0.d0)then
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
      if(irtc /= 0)then
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
      if(ilp > lt)then
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
      real*8 ,intent(in):: xi,pxi,yi,pyi,dp,p1,h1,ds
      real*8 xi3a,gx,gy,gz,dpr,
     $     dpgx,dpgy,dpgz,xir,pxir,zir,pzi,pyir,
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
      if(phi /= 0.d0)then
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
      if(ilp > lt)then
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
        kax=ktadaloc(-1,nph,klx)
        klx%attr=ior(klx%attr,lconstlist)
        itp=1
        ilp=0
        kt=kphtable(1)
        do i=1,nph
          ilp=ilp+1
          if(ilp > lt)then
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
          if(kphtable(i) /= 0)then
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
     $     -pcvt%dx,-pcvt%dy,pcvt%theta,pcvt%dtheta,0.d0,0.d0,
     $     pcvt%phi0,.false.)
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


      module tspin
      use macphys

      real*8, parameter :: pst=8.d0*sqrt(3.d0)/15.d0,
     $     sflc=.75d0*(elradi/finest)**2,cfpd=5.d0*sqrt(3.d0)/8.d0
      real*8, parameter:: gmin=-0.9999d0,cave=8.d0/15.d0/sqrt(3.d0)
      real*8, parameter:: cuu=11.d0/27.d0,cl=1.d0+gspin

      integer*4 ,parameter :: mord=6,lind=13
      integer*4 , parameter ::
     $     mlen(mord) =(/18,105, 392,1134, 2772, 6006/),
     $     mleni(mord)=(/24,234,1456,6825,26208,86632/)

      real*8 cphi0,sphi0

      type scmat
        complex*16 , pointer :: cmat(:,:,:)
        integer*4 , pointer :: ind(:,:)
        integer*4 , pointer :: ias(:)
        integer*4 nind,iord,id,maxi
      end type

      contains
        subroutine spinitrm(rm,nord,id,l)
        implicit none
        type (scmat) , intent(inout):: rm
        integer*4 , intent(in)::nord,id,l
        allocate(rm%cmat(3,3,l),rm%ind(lind,l),rm%ias(l))
        rm%nind=0
        rm%iord=nord
        rm%id=id
        rm%maxi=l
        return
        end

        integer*4 function iaind(rm,ind) result(ia)
        type (scmat) , intent(inout):: rm
        integer*4 , intent(in):: ind(lind)
        integer*4 i1,i2,im,k
        logical*4 found
        i1=1
        i2=rm%nind
        im=0
        do while(i2 .ge. i1)
          im=(i1+i2)/2
          ia=rm%ias(im)
          found=.true.
          do k=1,lind
            if(rm%ind(k,ia) > ind(k))then
              i2=im-1
              found=.false.
              exit
            elseif(rm%ind(k,ia) < ind(k))then
              i1=im+1
              found=.false.
              exit
            endif
          enddo
          if(found)then
            return
          endif
        enddo
        ia=rm%nind+1
        if(i1 == im+1)then
          im=im+1
        endif
        if(im < ia)then
          rm%ias(im+1:ia)=rm%ias(im:rm%nind)
        endif
        rm%ias(im)=ia
        rm%nind=ia
        if(rm%nind > rm%maxi)then
          write(*,*)'Insufficient matrix table ',rm%id,rm%nind,
     $         rm%iord,ind
          stop
        endif
        rm%ind(:,ia)=ind
        rm%cmat(:,:,ia)=(0.d0,0.d0)
        return
        end function

        integer*4 function indn(i,i1,i2,idx,imx,is)
        implicit none
        dimension indn(lind)
        integer*4 , intent(in)::i,i1,i2,idx,imx,is
        indn=0
        indn(i*2-1:i*2)=(/i1,i2/)
        indn(i+6)=idx
        indn(i+9)=imx
        indn(lind)=is
        return
        end function

        integer*4 function ind2(i,j,i1,i2,j1,j2)
        implicit none
        dimension ind2(lind)
        integer*4 , intent(in)::i,j,i1,i2,j1,j2
        ind2=0
        ind2(i*2-1:i*2)=(/i1,i2/)
        ind2(j*2-1:j*2)=(/j1,j2/)
        return
        end function

        integer*4 function ind3(i,i1,i2)
        implicit none
        dimension ind3(lind)
        integer*4 , intent(in)::i,i1,i2
        ind3(1:6)=1
        ind3(i*2-1:i*2)=(/i1,i2/)
        ind3(7:13)=0
        return
        end function

        subroutine spsetrm(am,ind,rm)
        implicit none
        type (scmat) , intent(inout)::rm
        integer*4 , intent(in) :: ind(lind)
        complex*16 , intent(in) ::am(3,3)
        rm%cmat(:,:,iaind(rm,ind))=am
        return
        end subroutine

        subroutine spdotrm(rma,rmb,rmc)
        implicit none
        type (scmat) , intent(in):: rma,rmb
        type (scmat) , intent(inout) :: rmc
        integer*4 i,j,ia
        do i=1,rma%nind
          do j=1,rmb%nind
            ia=iaind(rmc,rma%ind(:,i)+rmb%ind(:,j))
            rmc%cmat(:,:,ia)=rmc%cmat(:,:,ia)
     $           +matmul(rma%cmat(:,:,i),rmb%cmat(:,:,j))
          enddo
        enddo
        return
        end subroutine

        subroutine spaddrm(rma,rmb,rmc)
        implicit none
        type (scmat) ,intent(in):: rma,rmb
        type (scmat) ,intent(out):: rmc
        integer*4 ind(lind),i,ic
        do i=1,rma%nind
          ind=rma%ind(:,i)
          ic=iaind(rmc,ind)
          rmc%cmat(:,:,ic)=rmc%cmat(:,:,ic)+rma%cmat(:,:,i)
        enddo
        do i=1,rmb%nind
          ind=rmb%ind(:,i)
          ic=iaind(rmc,ind)
          rmc%cmat(:,:,ic)=rmc%cmat(:,:,ic)+rmb%cmat(:,:,i)
        enddo
        return
        end subroutine

        subroutine spcopyrm(rma,rmb)
        implicit none
        type (scmat) rma,rmb
        integer*4 n
        n=rma%nind
        rmb%cmat(:,:,1:n)=rma%cmat(:,:,1:n)
        rmb%ind(:,1:n)=rma%ind(:,1:n)
        rmb%ias(1:n)=rma%ias(1:n)
        rmb%nind=n
        return
        end subroutine

        subroutine spintrm(rm,dx,dy,dz,amx,amy,amz,ams,ssprd)
        use gammaf
        use mathfun
        use macmath
        implicit none
        type (scmat) , intent(inout) :: rm
        real*8 , intent(in) :: dx,dy,dz,amx,amy,amz,ams,ssprd
        real*8 a
        complex*16 cm(3,3),cei,cphi,cb
        integer*4 i,ia,ind(lind)
        do i=1,rm%nind
          ind=rm%ind(:,i)
          cphi=dcmplx(-ind(7)*dx-ind(8)*dy-ind(9)*dz,
     $         ind(10)*amx+ind(11)*amy+ind(12)*amz+ind(lind)*ams)
          cei=exp(-cphi)
          a=m_sqrt2*ssprd
c          cb=(1.d0-cei)/a
          cb=-cexp1(-cphi)/a
          cm=-rm%cmat(:,:,i)*m_2pi*m_sqrtpi
     $         *exp(cb**2-cphi)*cerfc(cb)/a
          rm%cmat(:,:,i)=-cm
          ind(7:lind)=0
          ia=iaind(rm,ind)
          rm%cmat(:,:,ia)=rm%cmat(:,:,ia)+cm
        enddo
        return
        end subroutine

        subroutine spmulrm(rm,cv,indv,rd)
        implicit none
        type (scmat) , intent(in) :: rm
        type (scmat) , intent(inout) :: rd
        complex*16 , intent(in):: cv
        integer*4 , intent(in):: indv(lind)
        integer*4 i,ia
        do i=1,rm%nind
          ia=iaind(rd,indv+rm%ind(:,i))
          rd%cmat(:,:,ia)=rd%cmat(:,:,ia)+cv*rm%cmat(:,:,i)
        enddo
        return
        end subroutine

        subroutine sprmulrm(rm,r)
        implicit none
        type (scmat) , intent(inout) :: rm
        real*8 , intent(in):: r
        rm%cmat(:,:,1:rm%nind)=r*rm%cmat(:,:,1:rm%nind)
        return
        end subroutine

        subroutine spcalcres(rm,rmi,m,dx,amx,ams,ssprd)
        implicit none
        type (scmat), intent(inout):: rm(mord),rmi(mord)
        integer*4 , intent(in)::m
        real*8 , intent(in)::dx(3),amx(3),ams,ssprd
        integer*4 k
        call spdotrm(rm(1),rm(m-1),rm(m))
        call sprmulrm(rm(m),1.d0/dble(m))
        call spcopyrm(rm(m),rmi(m))
        do k=1,m-1
          call spdotrm(rm(k),rmi(m-k),rmi(m))
        enddo
        call spintrm(rmi(m),dx(1),dx(2),dx(3),amx(1),amx(2),amx(3),ams,ssprd)
        return
        end

        subroutine spdepol(gxr,gxi,gyr,gyi,gzr,gzi,e1,e2,em,dx,ssprd,amx,ams,rmd)
        implicit none
        type (scmat) rm(mord),rmi(mord)
        real*8 , intent(in)::gxr(3),gxi(3),gyr(3),gyi(3),gzr(3),gzi(3),
     $       dx(3),amx(3),ams,e1(3),e2(3),em(3),ssprd
        real*8 , intent(out) :: rmd(3,3,3)
        integer*4 i,ia1,ia2,ia3,ia4,ia5,ia6,
     $       ia20,ia21,ia40,ia41,ia42,
c     $       ia60,ia61,
     $       ia62,ia63,i1,i2,k,
     $       ib40,ib41,ic40,ic41,ib60,ib61,ic60,ic61,id60,id61,
     $       ie61,ie62,if61,if62
        complex*16 gx,gy,gz,gxc,gyc,gzc
        complex*16 ,parameter :: cI=(0.d0,1.d0),c0=(0.d0,0.d0),c1=(1.d0,0.d0)
        real*8 de12,se12,sesq,e1e2
        do i=1,mord
          call spinitrm(rm(i),i,i,mlen(i))
          call spinitrm(rmi(i),i,i+mord,mleni(i))
        enddo
        do i=1,3
          gx=dcmplx(gxr(i),gxi(i))
          gxc=conjg(gx)
          gy=dcmplx(gyr(i)+gzi(i),gzr(i)-gyi(i))
          gyc=conjg(gy)
          gz=dcmplx(gyr(i)-gzi(i),gzr(i)+gyi(i))
          gzc=conjg(gz)
          ia1=iaind(rm(1),indn(i,0,1,1,-1,-1))
          rm(1)%cmat(:,:,ia1)=RESHAPE(0.25d0*gzc*(/
     $         c0,-cI,c1,
     $         cI,c0,c0,
     $         -c1,c0,c0/),(/3,3/))
          ia2=iaind(rm(1),indn(i,0,1,1,-1,0))
          rm(1)%cmat(:,:,ia2)=RESHAPE(0.5d0*gxc*(/
     $         c0,c0, c0,
     $         c0,c0,-c1,
     $         c0,c1, c0/),(/3,3/))
          ia3=iaind(rm(1),indn(i,0,1,1,-1,1))
          rm(1)%cmat(:,:,ia3)=RESHAPE(0.25d0*gyc*(/
     $         c0,cI,c1,
     $         -cI,c0,c0,
     $         -c1,c0,c0/),(/3,3/))
          ia4=iaind(rm(1),indn(i,1,0,1,1,-1))
          rm(1)%cmat(:,:,ia4)=conjg(rm(1)%cmat(:,:,ia3))
          ia5=iaind(rm(1),indn(i,1,0,1,1,0))
          rm(1)%cmat(:,:,ia5)=conjg(rm(1)%cmat(:,:,ia2))
          ia6=iaind(rm(1),indn(i,1,0,1,1,1))
          rm(1)%cmat(:,:,ia6)=conjg(rm(1)%cmat(:,:,ia1))
        enddo

        call spcopyrm(rm(1),rmi(1))
        call spintrm(rmi(1),dx(1),dx(2),dx(3),amx(1),amx(2),amx(3),ams,ssprd)
c        do k=2,mord
        do k=2,2
          call spcalcres(rm,rmi,k,dx,amx,ams,ssprd)
        enddo

        rmd=0.d0
        do i=1,3
          ia20=iaind(rmi(2),indn(i,0,2,0,0,0))
          ia21=iaind(rmi(2),indn(i,1,1,0,0,0))

          ia40=iaind(rmi(4),indn(i,0,4,0,0,0))
          ia41=iaind(rmi(4),indn(i,1,3,0,0,0))
          ia42=iaind(rmi(4),indn(i,2,2,0,0,0))
          i1=mod(i,3)+1
          ib40=iaind(rmi(4),ind2(i,i1,0,2,1,1))
          ib41=iaind(rmi(4),ind2(i,i1,1,1,1,1))
          i2=mod(i1,3)+1
          ic40=iaind(rmi(4),ind2(i,i2,0,2,1,1))
          ic41=iaind(rmi(4),ind2(i,i2,1,1,1,1))

c          ia60=iaind(rmi(6),indn(i,0,6,0,0,0))
c          ia61=iaind(rmi(6),indn(i,1,5,0,0,0))
          ia62=iaind(rmi(6),indn(i,2,4,0,0,0))
          ia63=iaind(rmi(6),indn(i,3,3,0,0,0))
          ib60=iaind(rmi(6),ind2(i,i1,0,2,2,2))
          ib61=iaind(rmi(6),ind2(i,i1,1,1,2,2))
          ic60=iaind(rmi(6),ind2(i,i2,0,2,2,2))
          ic61=iaind(rmi(6),ind2(i,i2,1,1,2,2))
          id60=iaind(rmi(6),ind3(i,0,2))
          id61=iaind(rmi(6),ind3(i,1,1))
          ie61=iaind(rmi(6),ind2(i,i1,1,3,1,1))
          ie62=iaind(rmi(6),ind2(i,i1,2,2,1,1))
          if61=iaind(rmi(6),ind2(i,i2,1,3,1,1))
          if62=iaind(rmi(6),ind2(i,i2,2,2,1,1))

          se12=(e1(i)+e2(i))*.5d0
          de12=(e1(i)-e2(i))*.5d0
c always de12 > 0
          sesq=(e1(i)**2+e2(i)**2)*.25d0
          e1e2=e1(i)*e2(i)*.25d0
          rmd(:,:,1)=rmd(:,:,1)
c     $         +e1(i)*dble(rmi(2)%cmat(:,:,ia21)+2.d0*rmi(2)%cmat(:,:,ia20))
c     $         +e2(i)*imag(rmi(2)%cmat(:,:,ia21)+2.d0*rmi(2)%cmat(:,:,ia20))
     $         +se12*dble(rmi(2)%cmat(:,:,ia21))
     $         +2.d0*de12*dble(rmi(2)%cmat(:,:,ia20))

          rmd(:,:,2)=rmd(:,:,2)
     $     +2.d0*(
     $         em(i)*(
     $         se12*dble(rmi(4)%cmat(:,:,ia42))
     $         +2.d0*de12*dble(rmi(4)%cmat(:,:,ia41)))
     $        +em(i1)*(
     $         se12*dble(rmi(4)%cmat(:,:,ib41))
     $         +2.d0*de12*dble(rmi(4)%cmat(:,:,ib40)))
     $        +em(i2)*(
     $         se12*dble(rmi(4)%cmat(:,:,ic41))
     $         +2.d0*de12*dble(rmi(4)%cmat(:,:,ic40))))
c     $     +de12*(6.d0*de12*dble(rmi(4)%cmat(:,:,ia40))
c     $           +(6.d0*se12+4.d0*em(i))
c     $              *dble(rmi(4)%cmat(:,:,ia41)))
c     $     +(2.d0*(em(i)*se12+e1e2)+3.d0*sesq)
c     $        *dble(rmi(4)%cmat(:,:,ia42))
 
          rmd(:,:,3)=rmd(:,:,3)
     $     +8.d0*(
     $         em(i )**2*(se12*dble(rmi(6)%cmat(:,:,ia63))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,ia62)))
     $        +em(i1)**2*(se12*dble(rmi(6)%cmat(:,:,ib61))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,ib60)))
     $        +em(i2)**2*(se12*dble(rmi(6)%cmat(:,:,ic61))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,ic60))))
     $     +4.d0*(
     $         em(i1)*em(i2)*(se12*dble(rmi(6)%cmat(:,:,id61))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,id60)))
     $        +em(i )*em(i1)*(se12*dble(rmi(6)%cmat(:,:,ie62))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,ie61)))
     $        +em(i )*em(i2)*(se12*dble(rmi(6)%cmat(:,:,if62))
     $         +2.d0*de12*dble(rmi(6)%cmat(:,:,if61))))
c     $     +de12*(3.d0*de12*(
c     $       10.d0*de12*dble(rmi(6)%cmat(:,:,ia60))
c     $       +(4.d0*em(i)+10.d0*se12)*dble(rmi(6)%cmat(:,:,ia61)))
c     $       +(em(i)*(16.d0*em(i)+12.d0*se12)+30.d0*sesq+36.d0*e1e2)
c     $           *dble(rmi(6)%cmat(:,:,ia62)))
c     $     +(em(i)*(8.d0*em(i)*se12+6.d0*sesq+4.d0*e1e2)
c     $       +3.d0*se12*(5.d0*sesq-2.d0*e1e2))
c     $      *dble(rmi(6)%cmat(:,:,ia63))
c     $
        enddo

        do i=1,mord
          deallocate(rm(i)%cmat,rmi(i)%cmat,rm(i)%ind,rmi(i)%ind,rm(i)%ias,rmi(i)%ias)
        enddo
        return
        end subroutine

        subroutine spnorm(srot,sps,smu,sdamp)
        use temw, only:gintd
        use mathfun, only:outer
        use macmath, only:m_2pi
        implicit none
        real*8 , intent(inout) :: srot(3,9)
        real*8 , intent(out) :: sps(3,3),smu,sdamp
        real*8 s,a(3,3),w(3,3),eig(2,3),dr(3),dsps(3),
     $       cm,sm,spsa1(3)
        real*8 , parameter :: smin=1.d-4
        integer*4 i
c        s=abs(dcmplx(srot(1,2),abs(dcmplx(srot(2,2),srot(3,2)))))
        s=norm2(srot(:,2))
        srot(:,2)=srot(:,2)/s
        s=srot(1,1)*srot(1,2)+srot(2,1)*srot(2,2)+srot(3,1)*srot(3,2)
        srot(:,1)=srot(:,1)-s*srot(:,2)
        s=abs(dcmplx(srot(1,1),abs(dcmplx(srot(2,1),srot(3,1)))))
        srot(:,1)=srot(:,1)/s
        srot(:,3)=outer(srot(:,1),srot(:,2))
        sps(1,1)=srot(2,3)-srot(3,2)
        sps(2,1)=srot(3,1)-srot(1,3)
        sps(3,1)=srot(1,2)-srot(2,1)
        s=abs(dcmplx(sps(1,1),abs(dcmplx(sps(2,1),sps(3,1)))))
        if(s < smin)then
          a=srot(:,1:3)
          call teigen(a,w,eig,3,3)
          do i=1,3
            if(eig(2,i) == 0.d0)then
              sps(:,1)=a(:,i)
              s=abs(dcmplx(sps(1,1),abs(dcmplx(sps(2,1),sps(3,1)))))
              exit
            endif
          enddo
        endif
        sps(:,1)=sps(:,1)/s
        dr=sps(:,1)-srot(:,1)*sps(1,1)-srot(:,2)*sps(2,1)
     $       -srot(:,3)*sps(3,1)
        a=srot(:,1:3)
        a(1,1)=a(1,1)-1.d0
        a(2,2)=a(2,2)-1.d0
        a(3,3)=a(3,3)-1.d0
        call tsolvg(a,dr,dsps,3,3,3)
        sps(:,1)=sps(:,1)+dsps
c        s=abs(dcmplx(sps(1,1),abs(dcmplx(sps(2,1),sps(3,1)))))
        s=norm2(sps(:,1))
        sps(:,1)=sps(:,1)/s
        if(abs(min(sps(1,1),sps(2,1),sps(3,1)))
     $       > abs(max(sps(1,1),sps(2,1),sps(3,1))))then
          sps(:,1)=-sps(:,1)
        endif
        dr=sps(:,1)-srot(:,1)*sps(1,1)-srot(:,2)*sps(2,1)
     $       -srot(:,3)*sps(3,1)
        sps(:,2)=0.d0
        if(abs(sps(1,1)) > abs(sps(2,1)))then
          sps(2,2)=1.d0
          s=sps(2,1)
        else
          sps(1,2)=1.d0
          s=sps(1,1)
        endif
        sps(:,2)=sps(:,2)-s*sps(:,1)
c        s=abs(dcmplx(sps(1,2),abs(dcmplx(sps(2,2),sps(3,2)))))
        s=norm2(sps(:,2))
        sps(:,2)=sps(:,2)/s
        sps(:,3)=outer(sps(:,1),sps(:,2))
        spsa1=srot(:,1)*sps(1,2)+srot(:,2)*sps(2,2)+srot(:,3)*sps(3,2)
        cm=dot_product(spsa1,sps(:,2))
        sm=dot_product(spsa1,sps(:,3))
        smu=atan(-sm,cm)
        sdamp=dot_product(sps(:,1),gintd)
        return
        end subroutine

        subroutine serot(xx,xy,yy,c1,c2,d1,d2,e1,e2,em1,em2)
        implicit none
        real*8 ,intent(in):: xx,xy,yy
        real*8 ,intent(inout):: c1,c2,d1,d2,e1,e2
        real*8 ,intent(out):: em1,em2
        real*8 tx,c,s
        tx=.5d0*atan(2.d0*xy,xx-yy)
        c=cos(tx)
        s=sin(tx)
        em1=c**2*xx+s**2*yy+2.d0*c*s*xy
        em2=c**2*yy+s**2*xx-2.d0*c*s*xy
        c1=  c*c1+s*c2
        c2=(-s*c1+  c2)/c
        d1=  c*d1+s*d2
        d2=(-s*d1+  d2)/c
        e1=  c*e1+s*e2
        e2=(-s*e1+  e2)/c
        return
        end

        subroutine srequpol(srot,sps,params,demit,sdamp,rm1,equpol)
        use temw,only:ipdampx,nparams,ipdampz,ipemx,ipemz,ipnup,
     $       ipnx,ipnz,r,dsg
        use macmath
        implicit none
        real*8 , intent(in)::srot(3,9),demit(21),sps(3,3),
     $       params(nparams),sdamp
        real*8 , intent(out)::equpol(3),rm1(3,3)
        real*8 drot(3,6),emit1(3),damp(3),ssprd,
     $       d1,d2,d3,d4,d5,d6,e1,e2,e3,e4,e5,e6,
     $       c1,c2,c3,c4,c5,c6,smu,
     $       dex1,dex2,dey1,dey2,dez1,dez2,
     $       rm(3,3),epol(3,3),b(3),rmd(3,3,3)
        integer*4 i
        smu=params(ipnup)*m_2pi
        drot=matmul(transpose(sps),matmul(srot(:,4:9),r))
        c1=drot(1,1)
        c2=drot(1,2)
        c3=drot(1,3)
        c4=drot(1,4)
        c5=drot(1,5)
        c6=drot(1,6)
        d1=drot(2,1)
        d2=drot(2,2)
        d3=drot(2,3)
        d4=drot(2,4)
        d5=drot(2,5)
        d6=drot(2,6)
        e1=drot(3,1)
        e2=drot(3,2)
        e3=drot(3,3)
        e4=drot(3,4)
        e5=drot(3,5)
        e6=drot(3,6)
        call serot(demit(1), demit(2), demit(3), c1,c2,d1,d2,e1,e2,dex1,dex2)
        call serot(demit(6), demit(9), demit(10),c3,c4,d3,d4,e3,e4,dey1,dey2)
        call serot(demit(15),demit(20),demit(21),c5,c6,d5,d6,e5,e6,dez1,dez2)
        emit1=params(ipemx:ipemz)
        damp=abs(params(ipdampx:ipdampz))
        ssprd=dot_product(drot(1,:),sqrt([emit1(1),emit1(1),emit1(2),emit1(2),emit1(3),emit1(3)]))
        call spdepol(
     $       (/c1,c3,c5/),(/c2,c4,c6/),
     $       (/d1,d3,d5/),(/d2,d4,d6/),
     $       (/e1,e3,e5/),(/e2,e4,e6/),
     $       (/dex1,dey1,dez1/),(/dex2,dey2,dez2/),
     $       emit1,damp,ssprd,
     $       params(ipnx:ipnz)*m_2pi,smu,rmd)
c        rm1=reshape([dsg(1),dsg(2),dsg(4),dsg(2),dsg(3),dsg(5),dsg(4),dsg(5),dsg(6)],[3,3])
        rm1=0.d0
        do i=1,1
c        do i=1,3
          rm1=rm1+rmd(:,:,i)
          rm=rm1
          rm(1,1)=rm(1,1)-sdamp
          b=(/-sdamp*pst,0.d0,0.d0/)
          call tsolva(rm,b,epol(:,i),3,3,3,min(1.d-8,sdamp/100.d0))
        enddo
        rm1(1,1)=rm1(1,1)-sdamp
        equpol=epol(1,:)
c        write(*,'(a,1p10g12.4)')'epol ',epol
c        write(*,'(a,1p10g12.4)')'rm1  ',rm1
c        write(*,'(a,1p10g12.4)')'pols ',matmul(rm1,epol)
        return
        end subroutine

        subroutine srotinit(srot)
        use temw,only:gintd
        implicit none
        real*8 ,intent(out):: srot(3,9)
        srot=0.d0
        srot(1,1)=1.d0
        srot(2,2)=1.d0
        srot(3,3)=1.d0
        gintd=0.d0
        return
        end subroutine

      end module

      module kradlib
      real*8 ,parameter :: dphipol=1.d-2
      real*8 , allocatable :: pxr0(:),pyr0(:),zr0(:),bsi(:)

      contains
        subroutine tradkf1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     px00,py0,zr00,bsi,al,k)
        use ffs_flag
        use tmacro
        use photontable, only:tphrec
        use mathfun, only:pxy2dpz,p2h,asinz
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0
        implicit none
        integer*4 ,parameter :: npmax=10000
        integer*4 , intent(in)::k
        integer*4 i
        real*8 , intent(inout)::x,px,y,py,z,g,dv
        real*8 , intent(in)::px00,py0,zr00,bsi,al
        real*8 dpx,dpy,dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,p2,h2,sx,sy,sz,
     $       ppa,an,dph,r1,r2,px0,xr,yr,rho,de
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
        theta=asinz(ppa)
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
        if(an /= 0.d0)then
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
          de=dph*uc
          dg=-de/(1.d0+2.d0*de)
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
            call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi,h1,p,de*uc)
          endif
        elseif(calpol)then
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi,h1,p,0.d0)
        endif
        return
        end subroutine

        subroutine tradkfn(np,xn,pxn,yn,pyn,zn,gn,dvn,sxn,syn,szn,al)
        use ffs_flag
        use tmacro
        use photontable, only:tphrec
        use mathfun, only:pxy2dpz,p2h,asinz
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0
        implicit none
        integer*4 ,parameter :: npmax=10000
        integer*4 , intent(in)::np
        real*8 , intent(inout)::
     $       xn(np),pxn(np),yn(np),pyn(np),zn(np),gn(np),dvn(np),
     $       sxn(np),syn(np),szn(np)
        real*8 , intent(in)::al
        integer*4 i,k
        real*8 dpx,dpy,dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,p2,h2,de,
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
          theta=asinz(ppa)
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
          if(an /= 0.d0)then
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
            de=dph*uc
            dg=-de/(1.d0+2.d0*de)
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
              call sprot(sxn(k),syn(k),szn(k),pxm,pym,ppx,ppy,ppz,bsi(k),h1,p,uc*de)
            endif
          elseif(calpol)then
            pxm=px0    +dpx*.5d0
            pym=pyr0(k)+dpy*.5d0
            call sprot(sxn(k),syn(k),szn(k),pxm,pym,ppx,ppy,ppz,bsi(k),h1,p,0.d0)
          endif
        enddo
        return
        end subroutine

        subroutine tradk1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     px00,py0,zr00,bsi0,al)
        use ffs_flag
        use tmacro
        use mathfun, only:pxy2dpz,p2h,asinz
        use tspin, only:cave,cl,cuu,gmin,sflc,cphi0,sphi0
        implicit none
        real*8 , intent(inout)::x,px,y,py,z,g,dv
        real*8 , intent(in)::px00,py0,zr00,bsi0,al
        real*8 dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,de,
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
        theta=asinz(ppa)
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        al1=al-z+zr00
        anp=anrad*h1*theta
        uc=cuc*h1**3/p0*theta/al1
        de=cave*anp*uc
        dg=-de/(1.d0+2.d0*de)
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
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi0,h1,p,uc*de)
        endif
        return
        end subroutine

        subroutine tradkn(np,xn,pxn,yn,pyn,zn,gn,dvn,sxn,syn,szn,al)
        use ffs_flag
        use tmacro
        use mathfun, only:pxy2dpz,p2h,asinz
        use tspin
        implicit none
        integer*4 , intent(in)::np
        real*8 , intent(inout)::
     $       xn(np),pxn(np),yn(np),pyn(np),zn(np),gn(np),dvn(np),
     $       sxn(np),syn(np),szn(np)
        real*8 , intent(in)::al
        real*8 dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,de,
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
          theta=asinz(ppa)
          pr=1.d0+gn(i)
          p=p0*pr
          h1=p2h(p)
          al1=al-zn(i)+zr0(i)
          anp=anrad*h1*theta
          uc=cuc*h1**3/p0*theta/al1
          de=cave*anp*uc
          dg=-de/(1.d0+2.d0*de)
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
            call sprot(sxn(i),syn(i),szn(i),pxm,pym,ppx,ppy,ppz,bsi(i),h1,p,uc*de)
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
        if(al /= 0.d0)then
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

        subroutine sprot(sx,sy,sz,pxm,pym,bx0,by0,bz0,bsi,h,p,ucde)
        use tmacro
        use ffs_flag, only:radpol
        use mathfun,only:pxy2dpz,sqrt1,xsincos
        use tspin
        implicit none
        real*8 ,intent(in):: pxm,pym,bsi,bx0,by0,bz0,h,p,ucde
        real*8 ,intent(inout):: sx,sy,sz
        real*8 bx,by,bz,bp,blx,bly,blz,btx,bty,btz,ct,pzm,
     $       gx,gy,gz,g,sux,suy,suz,bt,st,dst,
     $       sw,cosu,sinu,dcosu,xsinu
        pzm=1.d0+pxy2dpz(pxm,pym)
        bx=bx0/pzm
        by=by0/pzm
        bz=bz0/pzm+bsi
        bp=bx*pxm+by*pym+bz*pzm
        blx=bp*pxm
        bly=bp*pym
        blz=bp*pzm
        btx=bx-blx
        bty=by-bly
        btz=bz-blz
        ct=1.d0+h*gspin
        gx=ct*btx+cl*blx
        gy=ct*bty+cl*bly
        gz=ct*btz+cl*blz
        if(ucde > 0.d0 .and. radpol)then
          bt=norm2([btx,bty,btz])
          if(bt /= 0.d0)then
            st=(sx*btx+sy*bty+sz*btz)/bt
            dst=(pst-st)*cfpd*(p*p0/h**2)**2*ucde/bt
            sx=sx+dst*btx
            sy=sy+dst*bty
            sz=sz+dst*btz
c$$$            if(st /= 1.d0)then
c$$$              dr=dst/(1.d0-st**2)/bt
c$$$              dsx=dr*sx
c$$$              dsy=dr*sy
c$$$              dsz=dr*sz
c$$$              gx=gx+dsy*btz-dsz*bty
c$$$              gy=gy+dsz*btx-dsx*btz
c$$$              gz=gz+dsx*bty-dsy*btx
c$$$            else
c$$$              st1=st-dst
c$$$              sl1=sqrt(1-st1**2)
c$$$              sx=sl1*pxm+st1*btx/bt
c$$$              sy=sl1*pym+st1*bty/bt
c$$$              sz=sl1*pzm+st1*btz/bt
c$$$            endif
          endif
        endif
        g=norm2([gx,gy,gz])
        if(g /= 0.d0)then
          call xsincos(g,sinu,xsinu,cosu,dcosu)
          sw=-(sx*gx+sy*gy+sz*gz)*dcosu/g**2
          sinu=sinu/g
          sux=sy*gz-sz*gy
          suy=sz*gx-sx*gz
          suz=sx*gy-sy*gx
          sx=cosu*sx+sinu*sux+sw*gx
          sy=cosu*sy+sinu*suy+sw*gy
          sz=cosu*sz+sinu*suz+sw*gz
        endif
        sx= sx*cphi0+sz*sphi0
        sz=(sz-sx*sphi0)/cphi0
        return
        end subroutine

        subroutine tradke(trans,cod,beam,srot,al,phir0,bzh)
        use tmacro
        use temw,only:codr0,bzhr0,bsir0,calint,tinv6,gintd,transr,
     $       tmulbs,dsg
        use ffs_flag,only:radcod,calpol
        use mathfun, only:pxy2dpz,p2h,asinz,xsincos
        use tspin, only:cave,cl,cuu,gmin,sflc
        implicit none
        real*8 , intent(inout)::trans(6,12),cod(6),beam(42),srot(3,9)
        real*8 , intent(in)::al,bzh,phir0
        real*8 transi(6,6),tr1(6,6),dxpa(6),tr2(6,6),
     $       ddpz(6),dal(6),duc(6),dddpx(6),dddpy(6),ddg(6),
     $       dtheta(6),danp(6),dbeam(21),dpxi(6),dpyi(6),
     $       c1,dpx,dpy,ddpx,ddpy,pxr0,ct,pz00,das,bt,
     $       pr,px,py,pz,pz0,xpx,xpy,xpz,xpa,theta,th,
     $       p,h1,al1,anp,uc,dg,g,pr1,pxi,pyi,
     $       p2,h2,dee,cp,sp,b,pxm,pym,gi,dh1r,
     $       pxh,pyh,pzh,xpzb,btx,bty,btz,dct,sinu,cosu,dcosu,xsinu,
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
        xpa=norm2([xpx,xpy,xpz])/pr**2
        theta=asinz(xpa)
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
        if(irad > 6)then
          transi=matmul(trans(:,1:6),tinv6(transr))
c          call tmultr(transi,trans(:,1:6),6)
          tr2=transi
          if(bzh /= 0.d0)then
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
          if(xpa /= 0.d0)then
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
            if(bzh /= 0.d0)then
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
            dee=anp*uc**2*cuu
            pxm=pxi+px
            pym=pyi+py
            b=bzh*.5d0
            dbeam=0.d0
            dbeam(3)=(beam(3)+b*(2.d0*beam(5)+b*beam(6))
     $           +(pxm**2+pxi**2+px**2)/6.d0)*dee
            dbeam(8) =(beam(8)-b*(beam(2)-beam(10)+b*beam(4))
     $           +(pxm*pym+pxi*pyi+px*py)/6.d0)*dee
            dbeam(10)=(beam(10)+b*(-2.d0*beam(7)+b*beam(1))
     $           +(pym**2+pyi**2+py**2)/6.d0)*dee
            dbeam(17)=pxm*dee*.5d0
            dbeam(19)=pym*dee*.5d0
            dbeam(21)=dee
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
            pyh=(pyi +py)/pr*.5d0
            dpxh=(tr2(2,:)+dpxr0)/pr*.5d0
            dpxh(6)=dpxh(6)-pxh/pr
            dpyh=(tr2(4,:)+dpyi)/pr*.5d0
            dpyh(6)=dpyh(6)-pyh/pr
            pzh=1.d0+pxy2dpz(pxh,pyh)
            dpzh=-(pxh*dpxh+pyh*dpyh)/pzh
            bp =(xpx*pxh+xpy*pyh+xpzb*pzh)/pr**2
            dbp=(dxpx*pxh+xpx*dpxh+dxpy*pyh
     $           +xpy*dpyh+dxpzb*pzh+xpzb*dpzh)/pr**2
            dbp(6)=dbp(6)-2.d0*bp/pr
            blx=bp*pxh
            bly=bp*pyh
            blz=bp*pzh
            btx=xpx/pr**2-blx
            bty=xpy/pr**2-bly
            btz=xpzb/pr**2-blz
            dblx=dbp*pxh+bp*dpxh
            dbly=dbp*pyh+bp*dpyh
            dblz=dbp*pzh+bp*dpzh
            dbtx=dxpx/pr**2-dblx
            dbtx(6)=dbtx(6)-2.d0*xpx/pr**3
            dbty=dxpy/pr**2-dbly
            dbty(6)=dbty(6)-2.d0*xpy/pr**3
            dbtz=dxpzb/pr-dblz
            dbtz(6)=dbtz(6)-2.d0*xpzb/pr**3
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
c          dxpy=ddpz*pxr0+pz*dpxr0-tr2(2,:)*pz0-px*dpz0
c            if(abs(theta) < 1e-4)then
c              write(*,'(a,1p11g11.3)')'tradke ',al,theta,ct,dct,ct*dbty(6),dct*bty,srot(2,4),srot(2,9)
c              write(*,'(a,1p11g11.3)')'$dxpy       ',dxpy
c              write(*,'(a,1p11g11.3)')'$tr2(2      ',tr2(2,:)
c            endif
c            dsr2=dgy(1)*transr(1,:)+dgy(2)*transr(2,:)
c     $           +dgy(3)*transr(3,:)+dgy(4)*transr(4,:)
c     $           +dgy(5)*transr(5,:)+dgy(6)*transr(6,:)
c            if(abs(theta) < 1e-4)then
c              write(*,'(a,1p11g11.3)')'$dsr2      ',dgy(1),transr(1,6),dgy(1)*transr(1,6),dgy(6),dgy(6)*transr(6,6),dsr2(6)
c              write(*,'(a,1p11g11.3)')'$dgy       ',dgy
c            endif
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
            call tadddsg(dsg,srot,dgx,dgy,dgx,dbeam)
            g=norm2([gx,gy,gz])
            if(g /= 0.d0)then
              bt=norm2([btx,bty,btz])
              call xsincos(g,sinu,xsinu,cosu,dcosu)
c              th=tan(.5d0*g)
c              sinu=2.d0*th/(1.d0+th**2)
c              dcosu=th*sinu
cc              sinu=sin(g)
cc              dcosu=2.d0*sin(g*.5d0)**2
c              cosu=1.d0-dcosu
              sinu=sinu/g
              sx=srot(1,:)
              sy=srot(2,:)
              sz=srot(3,:)
              sw=(sx*gx+sy*gy+sz*gz)/g
              gintd=gintd+sw(1:3)*sflc*anp*(bt*h1*p/al1)**2
              sw=-sw*dcosu/g
              sux=sy*gz-sz*gy
              suy=sz*gx-sx*gz
              suz=sx*gy-sy*gx
              srot(1,:)=cosu*sx+sinu*sux+sw*gx
              srot(2,:)=cosu*sy+sinu*suy+sw*gy
              srot(3,:)=cosu*sz+sinu*suz+sw*gz
            endif
            srot(1,:)=  cp*srot(1,:)+sp*srot(3,:)
            srot(3,:)=(-sp*srot(1,:)+srot(3,:))/cp
          endif
        endif
        codr0(1:6)=cod(1:6)
        transr=trans(:,1:6)
        bzhr0=bzh
        bsir0=0.d0
        return
        end subroutine

        subroutine tadddsg(dsg,srot,dgx0,dgy0,dgz0,dbeam)
        implicit none
        real*8 ,intent(inout):: dsg(6)
        real*8 ,intent(in):: dgx0(6),dgy0(6),dgz0(6),dbeam(21),srot(3,9)
        real*8 dg0(3,6),sri(3,3),dgx(6),dgy(6),dgz(6)
        sri=transpose(srot(:,1:3))
        dg0(1,:)=dgx0
        dg0(2,:)=dgy0
        dg0(3,:)=dgz0
        dg0=matmul(sri,dg0)
        dgx=dg0(1,:)
        dgy=dg0(2,:)
        dgz=dg0(3,:)
        dsg(1)=dsg(1)
     $       +dgx(2)*(dgx(2)*dbeam(3 )+2.d0*dgx(4)*dbeam(8)+2.d0*dgx(6)*dbeam(17))
     $       +dgx(4)*(dgx(4)*dbeam(10)+2.d0*dgx(6)*dbeam(19))+dgx(6)*dgx(6)*dbeam(21)
        dsg(2)=dsg(2)
     $       +dgx(2)*(dgy(2)*dbeam(3 )+dgy(4)*dbeam(8 )+dgy(6)*dbeam(17))
     $       +dgx(4)*(dgy(2)*dbeam(8 )+dgy(4)*dbeam(10)+dgy(6)*dbeam(19))
     $       +dgx(6)*(dgy(2)*dbeam(17)+dgy(4)*dbeam(19)+dgy(6)*dbeam(21))
        dsg(3)=dsg(3)
     $       +dgy(2)*(dgy(2)*dbeam(3 )+2.d0*dgy(4)*dbeam(8)+2.d0*dgy(6)*dbeam(17))
     $       +dgy(4)*(dgy(4)*dbeam(10)+2.d0*dgy(6)*dbeam(19))+dgy(6)*dgy(6)*dbeam(21)
        dsg(4)=dsg(4)
     $       +dgx(2)*(dgz(2)*dbeam(3 )+dgz(4)*dbeam(8 )+dgz(6)*dbeam(17))
     $       +dgx(4)*(dgz(2)*dbeam(8 )+dgz(4)*dbeam(10)+dgz(6)*dbeam(19))
     $       +dgx(6)*(dgz(2)*dbeam(17)+dgz(4)*dbeam(19)+dgz(6)*dbeam(21))
        dsg(5)=dsg(5)
     $       +dgy(2)*(dgz(2)*dbeam(3 )+dgz(4)*dbeam(8 )+dgz(6)*dbeam(17))
     $       +dgy(4)*(dgz(2)*dbeam(8 )+dgz(4)*dbeam(10)+dgz(6)*dbeam(19))
     $       +dgy(6)*(dgz(2)*dbeam(17)+dgz(4)*dbeam(19)+dgz(6)*dbeam(21))
        dsg(6)=dsg(6)
     $       +dgz(2)*(dgz(2)*dbeam(3 )+2.d0*dgz(4)*dbeam(8)+2.d0*dgz(6)*dbeam(17))
     $       +dgz(4)*(dgz(4)*dbeam(10)+2.d0*dgz(6)*dbeam(19))+dgz(6)*dgz(6)*dbeam(21)
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
