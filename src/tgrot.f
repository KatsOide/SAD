      subroutine tgrot(a,geo0,geo1)
      implicit none
      real*8 geo0(3,3),geo1(3,3),a(3),chi2,a13,a33,chi1,a21,a22,chi3
      chi2=asin(
     1    -geo0(1,2)*geo1(1,3)-geo0(2,2)*geo1(2,3)-geo0(3,2)*geo1(3,3))
      a13= geo0(1,1)*geo1(1,3)+geo0(2,1)*geo1(2,3)+geo0(3,1)*geo1(3,3)
      if(a13 .eq. 0.d0)then
        chi1=0.d0
      else
        a33= geo0(1,3)*geo1(1,3)+geo0(2,3)*geo1(2,3)+geo0(3,3)*geo1(3,3)
        chi1=-atan2(a13,a33)
      endif
      a21= geo0(1,2)*geo1(1,1)+geo0(2,2)*geo1(2,1)+geo0(3,2)*geo1(3,1)
      if(a21 .eq. 0.d0)then
        chi3=0.d0
      else
        a22= geo0(1,2)*geo1(1,2)+geo0(2,2)*geo1(2,2)+geo0(3,2)*geo1(3,2)
        chi3=atan2(-a21,a22)
      endif
      a(1)=chi1
      a(2)=chi2
      a(3)=chi3
      return
      end

      subroutine tggeo1(k,geo1)
      use tfstk
      use ffs_pointer
      use kyparam
      implicit none
      integer*4 k,i
      integer*8 lp
      real*8 geo1(3,3),chi3,cschi3,snchi3,g1,
     $     chi1,chi2,cschi1,snchi1,cschi2,snchi2,
     $     tfchi
      lp=latt(k)
      chi1=rlist(lp+ky_DPX_SOL)
      chi2=rlist(lp+ky_DPY_SOL)
      cschi1=cos(chi1)
      snchi1=sin(chi1)
      cschi2=cos(chi2)
      snchi2=sin(chi2)
      chi3=tfchi(geo(1,1,k),3)
      cschi3=cos(chi3)
      snchi3=sin(chi3)
      do 110 i=1,3
        geo1(i,1)= geo(i,1,k)*cschi3+geo(i,2,k)*snchi3
        geo1(i,2)=-geo(i,1,k)*snchi3+geo(i,2,k)*cschi3
        g1       = geo1(i,2)*snchi2+geo(i,3,k)*cschi2
        geo1(i,2)= geo1(i,2)*cschi2-geo(i,3,k)*snchi2
        geo1(i,3)= g1*cschi1+geo1(i,1)*snchi1
        geo1(i,1)=-g1*snchi1+geo1(i,1)*cschi1
 110  continue
      return
      end

      subroutine tfgeofrac(lxp,gv)
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 lxp,lxs,ibz
      real*8 tfbzs,bzs,cod(6),f,xi,pxi,yi,pyi,ds,gi,geo2(3,3),
     $     pos0,gam0,gv(12),geo1(12)
      call tmov(geo(1,1,lxp+1),geo1,12)
      pos0=pos(lxp+1)
      gam0=gammab(lxp+1)
      lxs=ibzl(2,lxp)
c      write(*,*)'tfgeofrac ',lxs
      if(lxs .eq. 0)then
        call tfgeo1(lxp,lxp+1,.true.,.false.)
      else
        call tggeo1(lxs,geo2)
        cod(1:6)=twiss(lxp,0,mfitdx:mfitddp)
        bzs=tfbzs(lxp,ibz)
        f=0.5d0*bzs
        xi=cod(1)
        yi=cod(3)
        pxi=cod(2)+f*yi
        pyi=cod(4)-f*xi
        ds=0.d0
        gi=cod(6)
        call tsgeo1(lxp,xi,pxi,yi,pyi,ds,gi,bzs,geo2,1,1.d0)
      endif
      call tmov(geo(1,1,lxp+1),gv,12)
      call tmov(geo1,geo(1,1,lxp+1),12)
      gammab(lxp+1)=gam0
      pos(lxp+1)=pos0
      return
      end

