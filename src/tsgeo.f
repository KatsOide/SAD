      subroutine tsgeo(k,ke,ke1,sol)
      use kyparam
      use tfstk
      use ffs
      use sad_main
      use ffs_pointer
      use tffitcode
      use ffs_seg
      use geolib
      implicit none
      real*8, parameter:: conv=3.d-16
      integer*4 ,intent(in)::k
      integer*4 ,intent(inout)::ke,ke1
      integer*4 i,kg,k1,k2,idir,kbz
      real*8 geos(3,4),bzs,
     $     chi1,chi2,chi3,
     $     xi,yi,pxi,pyi,
     $     dir,tfbzs,pos0,ds,a(3),dg4(3),
     $     gi,geot(3,3)
      integer*8 led
      logical*4 ,intent(inout):: sol
      type (sad_comp),pointer :: cmp,cmp1,cmp2,cmpp,cmpe
      sol=.false.
      call compelc(k,cmp1)
      kg=0
      if(cmp1%value(ky_BND_SOL) .ne. 0.d0)then
        if(cmp1%value(ky_GEO_SOL) .ne. 0.d0)then
          kg=k
        endif
        do i=k+1,nlat-1
          if(idtypec(i) .eq. icSOL)then
            call compelc(i,cmp)
            if(cmp%value(ky_BND_SOL) .ne. 0.d0)then
              if(cmp%value(ky_GEO_SOL) .ne. 0.d0)then
                if(kg .ne. 0)then
                  write(*,*)' Duplicated GEO of Solenoid ',
     $                 pname(idelc(k))(1:lpnamec(k))
                  return
                endif
                kg=i
              endif
              ke=i
              sol=.true.
              exit
            endif
          endif
        enddo
      endif
      if(.not. sol)then
        write(*,*)' Missing BOUND of Solenoid ',
     $       pname(idelc(k))(1:lpnamec(k))
        return
      endif
      if(kg .eq. 0)then
        write(*,*)' Missing GEO of Solenoid ',
     $     pname(idelc(k))(1:lpnamec(k))
        sol=.false.
        return
      endif
      call compelc(ke,cmp2)
      ke1=ke+1
      if(kg .eq. k)then
        k1=k
        k2=ke
        cmpp=>cmp1
        cmpe=>cmp2
        dir=1.d0
        idir=1
      else
        k1=ke
        k2=k
        cmpp=>cmp2
        cmpe=>cmp1
        dir=-1.d0
        idir=-1
      endif
      cmpp%value(ky_DZ_SOL)=0.d0
      chi1=cmpp%value(ky_DPX_SOL)
      chi2=cmpp%value(ky_DPY_SOL)
      if(idir .gt. 0)then
        pos0=0
        bzs=tfbzs(k1,kbz)
        chi3=tfchi(geo(:,:,k),3)
        geo(:,1:3,k+1)=tfderotgeo(geo(:,:,k),(/chi1,chi2,chi3/))
        geo(:,4,k+1)=geo(:,4,k)
        cmp1%value(ky_CHI1_SOL:ky_CHI3_SOL)=
     $       tgrot(geo(:,:,k),geo(:,:,k+1))
        pos(k+1)=pos(k)
      else
        geos=geo(:,:,k)
        pos0=pos(k)
        bzs=tfbzs(k1-1,kbz)
        geo(:,:,ke)=geoini
        geo(:,:,ke1)=tfchitogeo((/chi1,chi2,0.d0/))
        cmp2%value(ky_CHI1_SOL:ky_CHI3_SOL)=
     $       tgrot(geo(:,1:3,ke),geo(:,:,ke1))
        pos(ke:ke1)=0.d0
        geo(:,4,ke:ke1)=0.d0
      endif
      xi=cmpp%value(ky_DX_SOL)
      yi=cmpp%value(ky_DY_SOL)
      pxi=-sin(chi1)*cos(chi2)
      pyi=-sin(chi2)
      led=idvalc(k2)
      ds=0.d0
      gi=0.d0
      if(cmpp%value(ky_FRIN_SOL) .eq. 0.d0)then
        call tsfrin(1,xi,pxi,yi,pyi,ds,gi,bzs)
        ds=-ds*dir
      endif
      pxi=pxi+yi*bzs*.5d0
      pyi=pyi-xi*bzs*.5d0
      do i=k1+idir,k2,idir
        call tsgeo1(i,xi,pxi,yi,pyi,ds,gi,bzs,idir,dir)
      enddo
      rlist(led+ky_DX_SOL:led+ky_DZ_SOL)=(/xi,yi,ds/)
      cmpe%value(ky_DX_SOL:ky_DZ_SOL)=(/xi,yi,ds/)
      if(idir .gt. 0)then
        chi3=tfchi(geo(:,1:3,ke1),3)
        geo(:,1:3,ke1)=tfderotgeo(geo(:,1:3,ke1),(/0.d0,0.d0,chi3/))
        cmp2%value(ky_CHI1_SOL:ky_CHI3_SOL)=
     $       tgrot(geo(:,1:3,ke),geo(:,:,ke1))
      else
        a=tgrot(geo(:,:,k),geo(:,:,ke1))
        geot=matmul(tfderotgeo(geos,(/0.d0,0.d0,a(3)/)),
     $       transpose(geo(:,1:3,k)))
        pos(k:ke1)=pos(k:ke1)+pos0-pos(k)
        dg4=geos(:,4)-matmul(geot,geo(:,4,k))
        do i=k+1,ke1
          geo(:,:,i)=matmul(geot,geo(:,:,i))
          geo(:,4,i)=geo(:,4,i)+dg4
        enddo
        geo(:,:,k)=geos
        cmp1%value(ky_CHI1_SOL:ky_CHI3_SOL)=
     $       tgrot(geos,geo(:,:,k+1))
      endif
      return
      end

      subroutine tsgeo1(i,xi,pxi,yi,pyi,ds,gi,bzs,idir,dir)
      use kyparam
      use tfstk
      use ffs
      use sad_main
      use ffs_pointer,only:pos,geo,idtypec,direlc,compelc,setdirelc,
     $     tsetfringep
      use tffitcode
      use ffs_seg
      use mathfun, only:akang
      use geolib
      implicit none
      type (sad_comp), pointer :: cmp
      type (sad_dlist), pointer :: lsegp
      integer*4 , intent(in):: i,idir
      integer*4 i0,i1,lt,irtc,mfr,kbz
      real*8 , intent(in)::dir
      real*8 , intent(inout):: xi,yi,pxi,pyi,ds,gi,bzs
      real*8 phi,ak,sinphi,a14,a12,a22,al,pzi,
     $     a24,dx,pxf,dy,pyf,dl,theta,phix,phiy,f,xf,
     $     zf,gf,dvf,ak1,ftable(4),bzs0,tfbzs,db,yf,
     $     sxf,syf,szf,theta2,
     $     trans(6,12),cod(6),beam(42),srot(3,9)
      logical*4 seg,dirf
      complex*16 cr1
      real*8 ,save::dummy(256)=0.d0
      i0=i+(1-idir)/2
      i1=2*i+1-i0
      lt=idtypec(i)
      call compelc(i,cmp)
      seg=tcheckseg(cmp,lt,al,lsegp,irtc)
      if(irtc .ne. 0)then
        call tffserrorhandle(i,irtc)
        return
      endif
      select case (lt)
      case (icDRFT)
        pzi=1.d0
     1       -(pxi**2+pyi**2)/(1.d0+sqrt((1.d0-pyi)*(1.d0+pyi)-pxi**2))
        phi=al*bzs/pzi*dir
        ak=bzs
        sinphi=sin(phi)
        a22= cos(phi)
        if(ak .eq. 0.d0)then
          a12=al/pzi*dir
          a14=0.d0
        else
          a12= sinphi/ak
          a14=merge(a12*sinphi/(1.d0+a22),(1.d0-a22)/ak,
     $         a22 .ge. 0.d0)
c     a14= 2.d0*sin(phi*.5d0)**2/ak
        endif
        a24= sinphi
        dx = a12*pxi+a14*pyi
        pxf= a22*pxi+a24*pyi
        dy =-a14*pxi+a12*pyi
        pyf=-a24*pxi+a22*pyi
        dl=(pxi**2+pyi**2)/pzi/(1.d0+pzi)*al
        xi=xi+dx
        yi=yi+dy
      case (icBEND)
        phi=cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND)
        theta=cmp%value(ky_ROT_BEND)
        phiy=phi*cos(theta)
        phix=phi*sin(theta)
        f=.5d0*bzs
        xf=xi
        pxf=(pxi-f*yi)*dir
        yf=yi
        pyf=(pyi+f*xi)*dir
        zf=0.d0
        gf=0.d0
        dvf=0.d0
        call tdrift(1,xf,pxf,yf,pyf,zf,gf,dvf,sxf,syf,szf,
     $       al,bzs*dir,phiy,phix,.false.)
        pxf=pxf*dir+f*yf
        pyf=pyf*dir-f*xf
        dl=-zf
        dx=xf-xi
        dy=yf-yi
        xi=xf
        yi=yf
      case (icQUAD)
        xf=xi
        f=.5d0*bzs
        pxf=(pxi-f*yi)*dir
        yf=yi
        pyf=(pyi+f*xi)*dir
        zf=0.d0
        theta=cmp%value(ky_ROT_QUAD)
        ak1=cmp%value(ky_K1_QUAD)
        theta2=theta+akang(dcmplx(ak1,0.d0),al,cr1)
        call setdirelc(i,direlc(i)*dir)
        dirf=direlc(i) .gt. 0.d0
        if(dirf)then
          mfr=nint(cmp%value(ky_FRMD_QUAD))
        else
          mfr=nint(cmp%value(ky_FRMD_QUAD))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        gf=0.d0
        dvf=0.d0
        if(al .ne. 0.d0)then
          call tsetfringep(cmp,icQUAD,ak1/al,ftable)
        else
          ftable=0.d0
        endif
        call tquad(1,xf,pxf,yf,pyf,zf,gf,dvf,
     $       sxf,syf,szf,
     $       al,ak1,bzs*dir,
     $       cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),theta,
     $       theta2,
     1       .false.,.true.,cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,cmp%value(ky_EPS_QUAD),.true.)
        call setdirelc(i,direlc(i)*dir)
        pxf=pxf*dir+f*yf
        pyf=pyf*dir-f*xf
        dl=-zf
        dx=xf-xi
        dy=yf-yi
        xi=xf
        yi=yf
      case (icMULT)
        f=bzs*.5d0
        cod=[xi,(pxi-f*yi)*dir,yi,(pyi+f*xi)*dir,0.d0,0.d0]
        call setdirelc(i,direlc(i)*dir)
        dirf=direlc(i) .gt. 0.d0
        call tinitr12(trans)
        if(seg)then
          call tmulteseg(trans,cod,beam,srot,i,cmp,bzs*dir,
     $         lsegp,.false.,1.d0)
        else
          call tmulte1(trans,cod,beam,srot,i,cmp,bzs*dir,.false.,1.d0)
        endif
        call setdirelc(i,direlc(i)*dir)
        xf=cod(1)
        yf=cod(3)
        pxf=cod(2)*dir+f*yf
        pyf=cod(4)*dir-f*xf
        dl=-cod(5)
        dx=xf-xi
        dy=yf-yi
        xi=xf
        yi=yf
      case (icCAVI)
        f=bzs*.5d0
        cod(1)=xi
        cod(2)=(pxi-f*yi)*dir
        cod(3)=yi
        cod(4)=(pyi+f*xi)*dir
        cod(5:6)=0.d0
        call setdirelc(i,direlc(i)*dir)
        dirf=direlc(i) .gt. 0.d0
        call tinitr12(trans)
        call tmulte(trans,cod,beam,srot,i,al,
     $       dummy,
     $       bzs*dir,
     $       0.d0,0.d0,0.d0,0.d0,
     1       cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $       0.d0,0.d0,0.d0,
     $       cmp%value(ky_ROT_CAVI),
     $       0.d0,0.d0,.false.,
     $       cmp%value(ky_FRIN_CAVI) .eq. 0.d0,
     $       0.d0,0.d0,0.d0,0.d0,
     $       cmp%value(ky_FRMD_CAVI),0.d0,0.d0,
     $       .true.,
     $       (cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI))*dir,
     $       cmp%value(ky_HARM_CAVI),
     $       cmp%value(ky_PHI_CAVI),cmp%value(ky_FREQ_CAVI),
     $       0.d0,1.d0,
     $       cmp%value(ky_APHI_CAVI) .ne. 0.d0,
     $       i)
        call setdirelc(i,direlc(i)*dir)
        xf=cod(1)
        yf=cod(3)
        pxf=cod(2)*dir+f*yf
        pyf=cod(4)*dir-f*xf
        dl=-cod(5)
        dx=xf-xi
        dy=yf-yi
        xi=xf
        yi=yf
      case (icSOL)
        al=0.d0
        dx=0.d0
        dy=0.d0
        dl=0.d0
        if(cmp%value(ky_BND_SOL) .eq. 0.d0)then
          bzs0=bzs
          bzs=tfbzs(merge(i-1,i,idir .lt. 0),kbz)
          db=bzs-bzs0
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            pxf=pxi-yi*bzs0*.5d0
            pyf=pyi+xi*bzs0*.5d0
            gi=0.d0
            xf=xi
            yf=yi
            call tsfrin(1,xf,pxf,yf,pyf,dl,gi,db)
            dx=xf-xi
            dy=yf-yi
            xi=xf
            yi=yf
            dl=-dl
            pxf=pxf+yi*bzs*.5d0
            pyf=pyf-xi*bzs*.5d0
          else
            pxf=pxi+yi*db*.5d0
            pyf=pyi-xi*db*.5d0
          endif
          cmp%value(ky_DX_SOL:ky_GEO_SOL)=0.d0
        else
          pxf=pxi-yi*bzs*.5d0
          pyf=pyi+xi*bzs*.5d0
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            gi=0.d0
            xf=xi
            yf=yi
            call tsfrin(1,xf,pxf,yf,pyf,dl,gi,-bzs)
            dx=xf-xi
            dy=yf-yi
            xi=xf
            yi=yf
            dl=-dl
          endif
          geo(:,:,i1)=tforbitgeo(geo(:,:,i0),(/xi,pxf,yi,pyf/))
          pxi=pxf
          pyi=pyf
          pos(i1)=pos(i0)+dl*dir
          ds=ds+dl
          return
        endif
      case default
        al=0.d0
        dx=0.d0
        pxf=pxi
        dy=0.d0
        pyf=pyi
        dl=0.d0
      end select
      geo(:,4,i1)=geo(:,4,i0)+geo(:,3,i0)*al*dir
      geo(:,1:3,i1)=geo(:,1:3,i0)
      pxi=pxf
      pyi=pyf
      pos(i1)=pos(i0)+(al+dl)*dir
      ds=ds+dl
      return
      end
