      subroutine qdtwis(dtwiss,ctrans,iclast,
     $     k0,l,idp,iv,nfam,nut,disp,dzfit)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use sad_main
      use tffitcode
      use temw, only:tmultr45
      implicit none
      type (sad_comp), pointer::cmp
      integer*4 nfam,nut,
     $     iclast(-nfam:nfam),k0,k,ke,l,idp,iv,
     $     iutk,iutl,i,k1
      real*8 dtwiss(mfittry),dcod(6),dcod2(6),dtrans(4,5),
     $     trans(4,5),trans1(4,5),dcod1(6),dtrans1(4,5),
     $     ctrans(4,7,-nfam:nfam),g1,dir,psi1,psi2,dt1,dt2,
     $     r,gr,detp,sqrdet,ddetp,dsqr,
     $     x11,x12,x21,x22,dx11,dx12,dx21,dx22,
     $     ax0,bx0,gx0,dax,dbx,
     $     y11,y12,y21,y22,dy11,dy12,dy21,dy22,
     $     ay0,by0,gy0,day,dby,
     $     bxx,cosmx,sinmx,byy,cosmy,sinmy,
     $     r1,r2,r3,r4,dr1,dr2,dr3,dr4
      logical*4 nzcod,disp,dzfit,normal
      k=k0
      g1=gammab(1)
      ke=idtypec(k)
      call compelc(k,cmp)
      dir=direlc(k)
      iutk=itwissp(k)
      iutl=itwissp(l)
      dcod(5)=0.d0
c      go to (1110,1120,1900,1140,1900,1160,1900,1160,1900,1160,
c     $       1900,1160,1900,1900,1900,1900,1900,1900,1900,1320,
c     $       1900,1320,1900,1900,1900,1900,1900,1900,1900,1900,
c     $       1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,
c     $       4100),ke
      select case (ke)
      case (icDRFT)
        call qddrif(dtrans,dcod,utwiss(1,idp,iutk))
      case (icBEND)
        if(iv .eq. ky_ANGL_BEND .or.
     $       iv .eq. ky_K1_BEND)then
          if(dir .gt. 0.d0)then
            psi1=cmp%value(ky_E1_BEND)
            psi2=cmp%value(ky_E2_BEND)
          else
            psi2=cmp%value(ky_E1_BEND)
            psi1=cmp%value(ky_E2_BEND)
          endif
          call qdbend(dtrans,dcod,cmp%value(ky_L_BEND),
     1         cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND),
     $         cmp%value(ky_ANGL_BEND),
     1         psi1,psi2,cmp%value(ky_K1_BEND),
     1         utwiss(1,idp,iutk),
     1         cmp%value(ky_DX_BEND),cmp%value(ky_DY_BEND),
     1         cmp%value(ky_ROT_BEND),iv)
          dt1=dtrans(1,5)
          dt2=dtrans(2,5)
        else
          call qdtrans(ke,iutk,k,k+1,
     $         iv,dtrans,dcod,idp)
        endif
      case (icQUAD)
        if(iv .eq. ky_K1_QUAD .or.
     $       iv .eq. ky_ROT_QUAD)then
          call qdquad(dtrans,dcod,
     $         cmp%value(ky_L_QUAD),cmp%value(ky_K1_QUAD),
     $         k,idp,cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     $         cmp%value(ky_ROT_QUAD),iv,nfam,nut)
          if(geocal .and. cmp%value(ky_DX_QUAD) .ne. 0.d0
     $         .or. cmp%value(ky_DY_QUAD) .ne. 0.d0)then
            call qdquad(dtrans1,dcod1,
     $           cmp%value(ky_L_QUAD),cmp%value(ky_K1_QUAD),
     $           k,0,cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     $           cmp%value(ky_ROT_QUAD),iv,nfam,nut)
            dcod(1:4)=dcod(1:4)-dcod1(1:4)
          endif
        else
          call qdtrans(ke,iutk,k,k+1,
     $         iv,dtrans,dcod,idp)
        endif
      case (icSEXT,icOCTU,icDECA,icDODECA)
        if(iv .eq. 2 .or. iv .eq. 4)then
          call qdthin(dtrans,dcod,ke,
     $         cmp%value(ky_L_THIN),cmp%value(ky_K_THIN),
     1         k,idp,cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     $         cmp%value(ky_ROT_THIN),iv,nfam,nut)
        else
          call qdtrans(ke,iutk,k,k+1,iv,dtrans,dcod,idp)
        endif
      case (icMULT,icSOL)
        call qdtrans(ke,iutk,k,k+1,iv,dtrans,dcod,idp)
      case (icMARK)
        if(k .eq. 1)then
          if(iv .eq. mfitddp)then
            call tftmatu(utwiss(1,idp,1),utwiss(1,idp,iutl),
     $           0.d0,0.d0,dtrans,1,l,.false.,trpt)
            r=sqrt(gammab(k)/gammab(l))
            do i=mfitdx,mfitdpy
              dtwiss(i)=dtrans(i-mfitdx+1,5)*r
            enddo
            go to 9000
          elseif(.not. cell)then
            k=min(nlat-1,max(2,1+int(
     $           rlist(latt(1)+ky_OFFSET_MARK))))
            gammab(1)=gammab(k)
            iutk=itwissp(k)
            call qdini(utwiss(1:ntwissfun,idp,1),
     $           utwiss(1:ntwissfun,idp,iutk),k,dtrans,dcod,iv)
            k1=k
            go to 2002
          else
            dtwiss=0.d0
            go to 9000
          endif
        else
          dtwiss=0.d0
          go to 9000
        endif
      case default
        dtwiss=0.d0
        go to 9000
      end select
      k1=k+1
 2002 gr=sqrt(gammab(k1)/gammab(l))
      dcod2(1:4)=dcod(1:4)*gr
      dcod2(5)=dcod(5)
      dcod2(6)=dcod(6)*gr**2
      dtrans(1:4,5)=dtrans(1:4,5)/gr
      call tftmatu(utwiss(1,idp,itwissp(k1)),
     $     utwiss(1,idp,iutl),
     $     utwiss(mfitnx,idp,nut),utwiss(mfitny,idp,nut),
     $     trans,k1,l,.false.,trpt)
      nzcod=.false.
      trans1=matmul(trans(:,1:4),dtrans)
      dcod1(1:4)=matmul(trans(:,1:4),dcod2(1:4))
      nzcod=nzcod .or. maxval(abs(dcod1(1:4))) .ne. 0.d0
      dcod1(5)=dcod2(5)
      if(dzfit .or. (nzcod .and. disp .and. .not. cell))then
        call qddtwiss(k,k1,l,
     $       trans,trans1,dcod,idp,
     $       ctrans(1,1,idp),iclast(idp),trpt)
        dcod1(5)=dcod(5)
      endif
      call qgettru(utwiss(1:ntwissfun,idp,iutk),
     $     utwiss(1:ntwissfun,idp,iutl),
     $     utwiss(3,idp,nut),utwiss(6,idp,nut),
     $     trans,k,l,.true.,.true.,trpt)
      r1=utwiss(mfitr1,idp,iutl)
      r2=utwiss(mfitr2,idp,iutl)
      r3=utwiss(mfitr3,idp,iutl)
      r4=utwiss(mfitr4,idp,iutl)
      detp=r1*r4-r2*r3
      sqrdet=sqrt(1.d0-detp)
      ddetp=trans1(1,1)*trans(2,2)+trans(1,1)*trans1(2,2)
     1     -trans1(1,2)*trans(2,1)-trans(1,2)*trans1(2,1)
      normal=utwiss(mfitdetr,idp,iutl) .lt. 1.d0
      if(normal)then
        dsqr=.5d0*ddetp/sqrdet
        x11=trans(1,1)/sqrdet
        x12=trans(1,2)/sqrdet
        x21=trans(2,1)/sqrdet
        x22=trans(2,2)/sqrdet
        dx11=(trans1(1,1)-x11*dsqr)/sqrdet
        dx12=(trans1(1,2)-x12*dsqr)/sqrdet
        dx21=(trans1(2,1)-x21*dsqr)/sqrdet
        dx22=(trans1(2,2)-x22*dsqr)/sqrdet
        dr1=-trans1(3,1)*x22-trans(3,1)*dx22
     1       +trans1(3,2)*x21+trans(3,2)*dx21
        dr2= trans1(3,1)*x12+trans(3,1)*dx12
     1       -trans1(3,2)*x11-trans(3,2)*dx11
        dr3=-trans1(4,1)*x22-trans(4,1)*dx22
     1       +trans1(4,2)*x21+trans(4,2)*dx21
        dr4= trans1(4,1)*x12+trans(4,1)*dx12
     1       -trans1(4,2)*x11-trans(4,2)*dx11
        dtwiss(mfitex)=trans(1,5)*dsqr+trans1(1,5)*sqrdet
     1           -dr4*trans(3,5)-r4*trans1(3,5)
     1           +dr2*trans(4,5)+r2*trans1(4,5)
        dtwiss(mfitepx)=trans(2,5)*dsqr+trans1(2,5)*sqrdet
     1           +dr3*trans(3,5)+r3*trans1(3,5)
     1           -dr1*trans(4,5)-r1*trans1(4,5)
        dtwiss(mfitdetr)=dr1*r4+r1*dr4-dr2*r3-r2*dr3
      else
        dsqr=-.5d0*ddetp/sqrdet
        x11=trans(3,1)/sqrdet
        x12=trans(3,2)/sqrdet
        x21=trans(4,1)/sqrdet
        x22=trans(4,2)/sqrdet
        dx11=(trans1(3,1)-x11*dsqr)/sqrdet
        dx12=(trans1(3,2)-x12*dsqr)/sqrdet
        dx21=(trans1(4,1)-x21*dsqr)/sqrdet
        dx22=(trans1(4,2)-x22*dsqr)/sqrdet
        dr1=(-trans1(1,1)*x22-trans(1,1)*dx22
     1       +trans1(1,2)*x21+trans(1,2)*dx21)
        dr2=( trans1(1,1)*x12+trans(1,1)*dx12
     1       -trans1(1,2)*x11-trans(1,2)*dx11)
        dr3=(-trans1(2,1)*x22-trans(2,1)*dx22
     1       +trans1(2,2)*x21+trans(2,2)*dx21)
        dr4=( trans1(2,1)*x12+trans(2,1)*dx12
     1       -trans1(2,2)*x11-trans(2,2)*dx11)
        dtwiss(mfitex)=trans(3,5)*dsqr+trans1(3,5)*sqrdet
     1       +(-dr4*trans(1,5)-r4*trans1(1,5)
     1       +dr2*trans(2,5)+r2*trans1(2,5))
        dtwiss(mfitepx)=trans(4,5)*dsqr+trans1(4,5)*sqrdet
     1       +( dr3*trans(1,5)+r3*trans1(1,5)
     1       -dr1*trans(2,5)-r1*trans1(2,5))
        dtwiss(mfitdetr)=-dr1*r4-r1*dr4+dr2*r3+r2*dr3
      endif
      dtwiss(mfitr1)=dr1
      dtwiss(mfitr2)=dr2
      dtwiss(mfitr3)=dr3
      dtwiss(mfitr4)=dr4
      dtwiss(mfitpex) =trans1(1,5)
      dtwiss(mfitpepx)=trans1(2,5)
      ax0=utwiss(mfitax,idp,iutk)
      bx0=utwiss(mfitbx,idp,iutk)
      gx0=(1.d0+ax0**2)/bx0
      dax      =(dx11*x22+x11*dx22+dx12*x21+x12*dx21)*ax0
     1          -(dx11*x21+x11*dx21)*bx0-(dx12*x22+x12*dx22)*gx0
      dtwiss(mfitax)=dax/(1.d0+utwiss(mfitax,idp,iutl)**2)
      dbx      =2.d0*(-(dx11*x12+x11*dx12)*ax0
     1               +dx11*x11*bx0+dx12*x12*gx0)
      dtwiss(mfitbx)=dbx/utwiss(mfitbx,idp,iutl)
      dtwiss(mfitnx)=(dx12*(bx0*x11-ax0*x12)-x12*(bx0*dx11-ax0*dx12))
     1          /(x12**2+(bx0*x11-ax0*x12)**2)
      dtwiss(mfitdx)=dcod1(1)
      dtwiss(mfitdpx)=dcod1(2)
      if(normal)then
        y11=trans(3,3)/sqrdet
        y12=trans(3,4)/sqrdet
        y21=trans(4,3)/sqrdet
        y22=trans(4,4)/sqrdet
        dy11=(trans1(3,3)-y11*dsqr)/sqrdet
        dy12=(trans1(3,4)-y12*dsqr)/sqrdet
        dy21=(trans1(4,3)-y21*dsqr)/sqrdet
        dy22=(trans1(4,4)-y22*dsqr)/sqrdet
        dtwiss(mfitey)=trans(3,5)*dsqr+trans1(3,5)*sqrdet
     1            +dr1*trans(1,5)+r1*trans1(1,5)
     1            +dr2*trans(2,5)+r2*trans1(2,5)
        dtwiss(mfitepy)=trans(4,5)*dsqr+trans1(4,5)*sqrdet
     1            +dr3*trans(1,5)+r3*trans1(1,5)
     1            +dr4*trans(2,5)+r4*trans1(2,5)
      else
        y11=trans(1,3)/sqrdet
        y12=trans(1,4)/sqrdet
        y21=trans(2,3)/sqrdet
        y22=trans(2,4)/sqrdet
        dy11=(trans1(1,3)-y11*dsqr)/sqrdet
        dy12=(trans1(1,4)-y12*dsqr)/sqrdet
        dy21=(trans1(2,3)-y21*dsqr)/sqrdet
        dy22=(trans1(2,4)-y22*dsqr)/sqrdet
        dtwiss(mfitey)=trans(1,5)*dsqr+trans1(1,5)*sqrdet
     1       +( dr1*trans(3,5)+r1*trans1(3,5)
     1       +dr2*trans(4,5)+r2*trans1(4,5))
        dtwiss(mfitepy)=trans(2,5)*dsqr+trans1(2,5)*sqrdet
     1       +( dr3*trans(3,5)+r3*trans1(3,5)
     1       +dr4*trans(4,5)+r4*trans1(4,5))
      endif
      dtwiss(mfitpey) =trans1(3,5)
      dtwiss(mfitpepy)=trans1(4,5)
      ay0=utwiss(mfitay,idp,iutk)
      by0=utwiss(mfitby,idp,iutk)
      gy0=(1.d0+ay0**2)/by0
      day      =(dy11*y22+y11*dy22+dy12*y21+y12*dy21)*ay0
     1          -(dy11*y21+y11*dy21)*by0-(dy12*y22+y12*dy22)*gy0
      dtwiss(mfitay)=day/(1.d0+utwiss(mfitay,idp,iutl)**2)
      dby      =2.d0*(-(dy11*y12+y11*dy12)*ay0
     1               +dy11*y11*by0+dy12*y12*gy0)
      dtwiss(mfitby)=dby/utwiss(mfitby,idp,iutl)
      dtwiss(mfitny)=(dy12*(by0*y11-ay0*y12)-y12*(by0*dy11-ay0*dy12))
     1          /(y12**2+(by0*y11-ay0*y12)**2)
      dtwiss(mfitdy)=dcod1(3)
      dtwiss(mfitdpy)=dcod1(4)
      dtwiss(mfitdz)=dcod1(5)
      dtwiss(mfitgmx)=2.d0*dtwiss(mfitax)*ax0/(1.d0+ax0**2)-
     $     dtwiss(mfitbx)
      dtwiss(mfitgmy)=2.d0*dtwiss(mfitay)*ay0/(1.d0+ay0**2)-
     $     dtwiss(mfitby)
c      dtwiss(mfitgmz)=2.d0*dtwiss(mfitaz)*az0/(1.d0+az0**2)-
c     $     dtwiss(mfitbz)
      if(l .eq. nlat)then
        bxx=sqrt(utwiss(mfitbx,idp,iutl)/utwiss(mfitbx,idp,1))
        cosmx=cos(utwiss(mfitnx,idp,iutl))
        sinmx=sin(utwiss(mfitnx,idp,iutl))
        dtwiss(mfittrx)=
     $       bxx*(.5d0*(cosmx+utwiss(mfitax,idp,1)*sinmx)*
     $         dtwiss(mfitbx)+
     1       (-sinmx+utwiss(mfitax,idp,1)*cosmx)*dtwiss(mfitnx))+
     1       (-.5d0*(cosmx-utwiss(mfitax,idp,iutl)*sinmx)*
     $         dtwiss(mfitbx)+
     1       (-sinmx-utwiss(mfitax,idp,iutl)*cosmx)*dtwiss(mfitnx)-
     1         sinmx*dax)/bxx
        byy=sqrt(utwiss(mfitby,idp,iutl)/utwiss(mfitby,idp,1))
        cosmy=cos(utwiss(mfitny,idp,iutl))
        sinmy=sin(utwiss(mfitny,idp,iutl))
        dtwiss(mfittry)=
     $       byy*(.5d0*(cosmy+utwiss(mfitay,idp,1)*sinmy)*
     $         dtwiss(mfitby)+
     1       (-sinmy+utwiss(mfitay,idp,1)*cosmy)*dtwiss(mfitny))+
     1       (-.5d0*(cosmy-utwiss(mfitay,idp,iutl)*sinmy)*
     $         dtwiss(mfitby)+
     1       (-sinmy-utwiss(mfitay,idp,iutl)*cosmy)*dtwiss(mfitny)-
     1       sinmy*day)/byy
      endif
 9000 gammab(1)=g1
      call resetnan(dtwiss)
      return
      end

      subroutine qdini(utwiss1,utwiss2,k2,dtrans,dcod,iv)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 iv,k2 
      real*8 utwiss1(ntwissfun),utwiss2(ntwissfun),
     $     trans(4,5),dtrans(4,5),dcod(6)
      real*8 b,detr,damu,dir
      dtrans=0.d0
      dcod=0.d0
      dir=direlc(1)
      if(iv .ge. mfitr1 .and. iv .le. mfitr4)then
        detr=utwiss1(mfitdetr)
      else
c     begin initialize for preventing compiler warning
        detr=0
c     end   initialize for preventing compiler warning
        call qgettru(utwiss1,utwiss2,0.d0,0.d0,
     $       trans,1,k2,.true.,.false.,.true.)
      endif
c      go to (100,200,300,400,500,600,700,800,900,1000,
c     $     1100,1200,1300,1400,1500,1600,1700,1800,1900),iv
c      return
      select case (iv)
      case (mfitax)
        b=-dir/utwiss1(mfitbx)
        dtrans(1:4,1)= trans(1:4,2)*b
      case (mfitbx)
        b=.5d0/utwiss1(mfitbx)
        dtrans(1:4,1)= trans(1:4,1)*b
        dtrans(1:4,2)=-trans(1:4,2)*b
      case (mfitay)
        b=-dir/utwiss1(mfitby)
        dtrans(1:4,3)= trans(1:4,4)*b
      case (mfitby)
        b=.5d0/utwiss1(mfitby)
        dtrans(1:4,3)= trans(1:4,3)*b
        dtrans(1:4,4)=-trans(1:4,4)*b
      case (mfitex)
        dtrans(1:4,5)=trans(1:4,1)
      case (mfitepx)
        dtrans(1:4,5)=dir*trans(1:4,2)
      case (mfitey)
        dtrans(1:4,5)=trans(1:4,3)
      case (mfitepy)
        dtrans(1:4,5)=dir*trans(1:4,4)
      case (mfitr1)
        if(detr .lt. 1.d0)then
          damu=-.5d0*utwiss1(mfitr4)/sqrt(1.d0-detr)
          dtrans(1,1)=damu
          dtrans(3,1)=-1.d0
          dtrans(2,2)=damu
          dtrans(3,3)=damu
          dtrans(2,4)=1.d0
          dtrans(4,4)=damu
        else
        endif
      case (mfitr2)
        if(detr .lt. 1.d0)then
          damu= .5d0*dir*utwiss1(mfitr3)/sqrt(1.d0-detr)
          dtrans(1,1)=damu
          dtrans(3,2)=-dir
          dtrans(2,2)=damu
          dtrans(3,3)=damu
          dtrans(1,4)=-dir
          dtrans(4,4)=damu
        else
        endif
      case (mfitr3)
        if(detr .lt. 1.d0)then
          damu= .5d0*dir*utwiss1(mfitr2)/sqrt(1.d0-detr)
          dtrans(1,1)=damu
          dtrans(4,1)=-dir
          dtrans(2,2)=damu
          dtrans(3,3)=damu
          dtrans(2,3)=-dir
          dtrans(4,4)=damu
        else
        endif
      case (mfitr4)
        if(detr .lt. 1.d0)then
          damu=-.5d0*utwiss1(mfitr1)/sqrt(1.d0-detr)
          dtrans(1,1)=damu
          dtrans(4,2)=-1.d0
          dtrans(2,2)=damu
          dtrans(3,3)=damu
          dtrans(1,3)=1.d0
          dtrans(4,4)=damu
        else
        endif
      case (mfitdx)
        dcod(1)=1.d0
      case (mfitdpx)
        dcod(2)=dir
      case (mfitdy)
        dcod(3)=1.d0
      case (mfitdpy)
        dcod(4)=dir
      case (mfitdz)
        dcod(5)=1.d0
      case (mfitnx,mfitny,mfitddp)
      case default
        return
      end select
      dcod(6)=0.d0
      return
      end

      subroutine nancheck(a,str)
      use tfstk, only: ktfenanq
      implicit none
      real*8 a(4,5)
      character*(*) str
      integer*4 i,j
      do i=1,5
        do j=1,4
          if(ktfenanq(a(j,i)))then
            write(*,*)str,' ',j,i
            return
          endif
        enddo
      enddo
      return
      end

      subroutine qddtwiss(k,k1,l,trans,dtrans,dcod,idp,
     $     ctrans,iclast,trpt)
      use mackw
      use tfstk
      use ffs_pointer
      use ffs_fit, only:nut
      use tffitcode
      use ffs, only:ffs_bound
      use temw, only:tmultr45
      implicit none
      type (ffs_bound) fbound,fbound1
      real*8 eps
      parameter (eps=1.d-4)
      integer*4 k,l,idp,itwk,itwl,
     $     iclast,k1,itwk1,ibg,ibe,itwe,ibe1,itwbe
      real*8 dtrans(4,5),dcod(6),trans(4,5),
     $     trans2(4,5),trans3(4,5),trans1(4,5),transe(4,5),
     $     transe2(4,5),cod2(6),code(6),dcode(6),
     $     w,ctrans(4,7)
      logical*4 over,trpt
c      equivalence (trans2,trans2s)
      itwk=itwissp(k)
      itwk1=itwissp(k1)
      itwl=itwissp(l)
      w=eps/(abs(dcod(1))+abs(dcod(2))+abs(dcod(3))+abs(dcod(4))
     $     +abs(dcod(5)))
      cod2(6)=utwiss(mfitddp,idp,itwk1)
      call tfbndsol(k,ibg,ibe)
      call tffsbound1(k1,l,fbound)
      if(iclast .gt. 0 .and. iclast .le. fbound%le .and.
     $     (iclast .ne. fbound%le .or. ctrans(4,7) .le. fbound%fe))then
        cod2(1:4)=ctrans(:,6)
        cod2(5)=ctrans(1,7)
        fbound1=fbound
        fbound1%lb=iclast
        fbound1%fb=ctrans(4,7)
        call qcod(1,fbound1,trans3,cod2,.true.,over)
        trans2=tmultr45(ctrans(:,1:5),trans3)
      else
        cod2(1:5)=utwiss(mfitdx:mfitdz, idp,itwk1)+w*dcod(1:5)
        if(ibg .eq. 0 .or. l .le. max(ibg,ibe))then
          call qcod(1,fbound,trans2,cod2,.true.,over)
        else
          if(ibg .lt. ibe)then
            ibe1=ibe+1
            itwe=itwissp(ibe1)
            code(2)=utwiss(mfitdpx,idp,itwe)
            code(4)=utwiss(mfitdpy,idp,itwe)
            fbound1=fbound
            fbound1%le=ibe1
            fbound1%fe=0.d0
            call qcod(1,fbound1,transe,cod2,.true.,over)
            transe(2,5)=transe(2,5)-cod2(2)+code(2)
            transe(4,5)=transe(4,5)-cod2(4)+code(4)
            call tftmatu(utwiss(1,idp,itwissp(ibe1)),
     $           utwiss(1,idp,itwl),
     $           0.d0,0.d0,
     $           transe2,ibe1,l,.false.,trpt)
            cod2(1:5)=utwiss(mfitdx:mfitdz, idp,itwk1)
            call qcod(1,fbound1,trans3,cod2,.true.,over)
            trans3(2,5)=trans3(2,5)-cod2(2)+code(2)
            trans3(4,5)=trans3(4,5)-cod2(4)+code(4)
            trans=tmultr45(trans3,transe2)
          else
            itwbe=itwissp(ibe)
            call tftmatu(utwiss(1,idp,itwk1),
     $           utwiss(1,idp,itwbe),
     $           0.d0,0.d0,
     $           transe,k1,ibe,.false.,trpt)
            dcode(1)=
     $            transe(1,1)*dcod(1)+transe(1,2)*dcod(2)
     $           +transe(1,3)*dcod(3)+transe(1,4)*dcod(4)
            dcode(2)=
     $            transe(2,1)*dcod(1)+transe(2,2)*dcod(2)
     $           +transe(2,3)*dcod(3)+transe(2,4)*dcod(4)
            dcode(3)=
     $            transe(3,1)*dcod(1)+transe(3,2)*dcod(2)
     $           +transe(3,3)*dcod(3)+transe(3,4)*dcod(4)
            dcode(4)=
     $            transe(4,1)*dcod(1)+transe(4,2)*dcod(2)
     $           +transe(4,3)*dcod(3)+transe(4,4)*dcod(4)
            cod2(1:4)=utwiss(mfitdx:mfitdpy, idp,itwbe)-w*dcode(1:4)
            cod2(5)=utwiss(mfitdz, idp,itwbe)
            fbound1=fbound
            fbound1%lb=ibe
            fbound1%fb=0.d0
            fbound1%le=k1
            fbound1%fe=0.d0
            call qcod(1,fbound1,trans2,cod2,.true.,over)
            transe(2,5)=transe(2,5)-w*dcode(2)
            transe(4,5)=transe(4,5)-w*dcode(4)
            transe=tmultr45(transe,trans2)
            call tftmatu(utwiss(1,idp,itwk1),
     $           utwiss(1,idp,itwl),
     $           0.d0,0.d0,
     $           transe2,k1,l,.false.,trpt)
          endif
          cod2(1:5)=utwiss(mfitdx:mfitdz, idp,itwl)
          trans2=tmultr45(transe,transe2)
        endif
      endif
c      call tmov(trans,transe2,20)
c      write(*,*)'qddtws-2 '
c      write(*,'(1p5g13.5)')((transe2(i,j),j=1,5),i=1,4)
c      write(*,'(1p5g13.5)')((trans2(i,j),j=1,5),i=1,4)
      iclast=fbound%le
      ctrans(4,7)=fbound%fe
      ctrans(:,1:5)=trans2
      trans2=(trans2-trans)/w
c      trans2s(1:20)=(trans2s(1:20)-trans(1:20))/w
      dcod(5)=(cod2(5)-utwiss(mfitdz,idp,itwl))/w
      ctrans(:,6)=cod2(1:4)
      ctrans(1:2,7)=cod2(5:6)
      call qgettru(utwiss(1:ntwissfun,idp,itwk),
     $     utwiss(1:ntwissfun,idp,itwk1),
     $     utwiss(3,idp,nut),utwiss(6,idp,nut),
     $     trans1,k,k1,.true.,.true.,trpt)
      dtrans=tmultr45(trans1,trans2)+dtrans
c      transe=tmultr45(trans1,trans2)
c      call tadd(transe,dtrans,dtrans,20)
      return
      end

      subroutine qdtrans(ke,kk1,j,je,
     $     iv,dtrans,dcod,idp)
      use kyparam
      use tfstk
      use ffs_pointer
      use sad_main
      use tffitcode
      use mackw
      use ffs_seg
      use temw, only:tmultr45
      implicit none
      type (sad_comp), pointer :: cmp
      real*8 eps,vmin
      parameter (eps=1.d-6,vmin=1.d-6)
      integer*4 ke,iv,idp,j,kk1,je
      real*8 dtrans(4,5),dcod(6),v0,wv,dv,trans2(4,5),
     $     cod2(20),trans1(4,5),cod1(20),trans(4,5)
      logical*4 over
      call compelc(j,cmp)
      v0=tfvalvar(j,iv)
      wv=1.d0
c      go to (
c     $     4900, 200,4900, 400,4900, 600,4900, 600,4900, 600,
c     $      600,4900,4900,4900,4900,4900,4900,4900,4900,6000,
c     $     4900,2200,4900,4900,4900,4900,4900,4900,4900,4900,
c     $     4900,4900,4900,4900,4900,4900,4900,4900,4900,4900,
c     $     4900),ke
      select case (ke)
      case (icBEND)
        wv=1.d0
      case (icQUAD)
        wv=1.d0
      case (icSEXT,icOCTU,icDECA,icDODECA)
        wv=1.d0
      case (icMULT)
        if(iv .ge. ky_K1_MULT)then
          wv=10.d0**((iv-ky_K1_MULT)/2)
        else
          wv=1.d0
        endif
      case default
        dtrans=0.d0
        dcod=0.d0
        return  
      end select        
      dv=max(abs(eps*v0),abs(vmin*wv))
c      call tfsetcmp(v0+dv,cmp,iv)
      cmp%value(iv)=v0+dv
      cod2(1:6)=utwiss(mfitdx:mfitddp,idp,kk1)
      call qtwiss1(0.d0,idp,j,je,trans2,cod2,.true.,over)
c      write(*,'(a,i5,1p12g13.5)')'qdtrans ',iv,v0+dv,trans2(1:4)
      cmp%value(iv)=v0-dv
      cod1(1:6)=utwiss(mfitdx:mfitddp,idp,kk1)
      call qtwiss1(0.d0,idp,j,je,trans1,cod1,.true.,over)
c      write(*,'(a,i5,1p12g13.5)')'qdtrans ',iv,trans1(1:4)
      trans2=(trans2-trans1)/(2.d0*dv)
      dcod(1:5)=(cod2(1:5)-cod1(1:5))/(2.d0*dv)
      cmp%value(iv)=v0
c      write(*,'(a,i5,1p8g15.7)')'qdtrans ',iv,dv,trans2(1:4)
      call qtentu(trans,cod1,utwiss(1,idp,kk1),.true.)
      dtrans=tmultr45(trans,trans2)
      return
      end
