      module tracklim
      real*8, parameter ::plimit=0.9999d0,zlimit=1.d10,vmax=.9d0,
     $     ampmax=0.9999d0
      real*8 xlimit
      end module

      subroutine tturn(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n)
      use tfstk
      use tmacro
      use tspin
      implicit none
      integer*4 np,n,kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      real*8 sx(np),sy(np),sz(np)
      logical*4 normal
      call tturn0(np,1,nlat,x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl,n,normal)
      return
      end

      subroutine tturn0(np,lb,le,x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl,n,normal)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs, only:ffs_bound
      use sad_main
      use ffs_pointer, only:compelc
      use tspin
      use tracklim
      implicit none
      type (ffs_bound) fbound
      type (sad_comp), pointer ::cmp
      type (sad_descriptor) :: dsave(kwMAX)
      integer*4 np,n,la,ls,nvar,lb,le
c      integer*4 isb,itwb,itwb1,itwb2,itwb3,itwb4,ntw
      real*8 pgev00
      integer*4 kptbl(np0,6),lv,itfdownlevel,irtc
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      real*8 sx(np0),sy(np0),sz(np0)
      logical*4 sol,chg,tfinsol,normal
      pgev00=pgev
      sol=tfinsol(lb)
      novfl=0
      la=0
      if(radlight)then
        call tlstore(np,x,y,z,dv,0.d0,0.d0,0.d0,0.d0,
     $       p0/h0*c,dvfs,.true.)
      endif
      call tffsbound1(lb,le,fbound)
      normal=fbound%lb .lt. fbound%le .or.
     $     fbound%lb .eq. fbound%le .and. fbound%fb .le. fbound%fe
      if(.not. normal)then
        return
      endif
      levele=levele+1
      if(fbound%lb .eq. fbound%le)then
        call compelc(fbound%lb,cmp)
        call qfracsave(fbound%lb,dsave,nvar,.true.)
        call qfracseg(cmp,cmp,fbound%fb,fbound%fe,chg,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(fbound%lb,irtc)
        else
          call tturn1(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n,
     $         sol,la,fbound%lb,fbound%lb)
        endif
        if(chg)then
          call qfracsave(fbound%lb,dsave,nvar,.false.)
        endif
      else
        call compelc(fbound%lb,cmp)
        if(fbound%fb .ne. 0.d0)then
          call qfracsave(fbound%lb,dsave,nvar,.true.)
          call qfracseg(cmp,cmp,fbound%fb,1.d0,chg,irtc)
          if(irtc .ne. 0)then
            call tffserrorhandle(fbound%lb,irtc)
          else
            call tturn1(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n,
     $           sol,la,fbound%lb,fbound%lb)
          endif
          if(chg)then
            call qfracsave(fbound%lb,dsave,nvar,.false.)
          endif
          ls=fbound%lb+1
        else
          ls=fbound%lb
        endif
        if(fbound%le .gt. ls+1 .or.
     $       fbound%le .eq. ls+1 .and. fbound%fe .eq. 0.d0)then
         call tturn1(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n,
     $         sol,la,ls,fbound%le-1)
        endif
        if(fbound%fe .ne. 0.d0)then
          call qfracsave(fbound%le,dsave,nvar,.true.)
          call compelc(fbound%le,cmp)
          call qfracseg(cmp,cmp,0.d0,fbound%fe,chg,irtc)
c          call qfraccomp(fbound%le,0.d0,fbound%fe,ideal,chg)
          call tturn1(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n,
     $         sol,la,fbound%le,fbound%le)
          if(chg)then
            call qfracsave(fbound%le,dsave,nvar,.false.)
          endif
        endif
      endif
      lv=itfdownlevel()
      if(trpt .and. codplt)then
c        call ttstat(np,x,px,y,py,z,g,dv,0.d0,
c     1       ' ',sa,ss,0.d0,
c     1       .false.,.false.,0)
c        call ttstat(np,x,px,y,py,z,g,dv,rlist(ilist(2,iwakepold+4)),
c     1       ' ',sa,ss,0.d0,
c     1       .false.,.false.,0)
c        itwb=ilist(1,iwakepold+6)
c     1       +ilist(1,iwakepold+5)*ilist(2,iwakepold+5)
c        ntw=(2*ilist(1,iwakepold+5)+1)*ilist(2,iwakepold+5)
c        itwb1=itwb+ntw*14-1
c        itwb2=itwb+ntw*15-1
c        itwb3=itwb+ntw*16-1
c        itwb4=itwb+ntw*17-1
c        isb=ilist(2,iwakepold+6)
c        rlist(itwb1+nlat)=sa(1)
c        rlist(itwb2+nlat)=sa(2)
c        rlist(itwb3+nlat)=sa(3)
c        rlist(itwb4+nlat)=sa(4)
c        call tt6621(ss,rlist(isb+21*(nlat-1)))
      endif
      if(rfsw)then
      else
        z(1:np)=0.d0
      endif
      if(pgev .ne. pgev00)then
        pgev=pgev00
        call tphyzp
      endif
      return
      end

      subroutine tturn1(np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl,n,
     $     sol,la,lbegin,lend)
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use ffs_wake
      use sad_main
      use ffs_pointer, only: direlc,compelc,twiss,pos,beamsize
      use tparastat
      use tfcsi, only:icslfno
      use ffs_seg
      use tspin
      use kradlib
      use mathfun
      use tracklim
      implicit none
      integer*4,parameter :: la1=15
      type (sad_comp), pointer:: cmp
      type (sad_dlist) , pointer ::lsegp
      integer*4 ,intent(inout):: np,la
      integer*4 ,intent(in):: n,lbegin,lend
      integer*4 kdx,kdy,krot
      integer*4 kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),
     $     g(np0),dv(np0),dpz,al1
      real*8 sx(np0),sy(np0),sz(np0)
      real*8 bz,al,ak0,ak1,tgauss,ph,harmf,
     $     sspac0,sspac,fw,dx,dy,rot,sspac1,sspac2,ak,rtaper,
     $     cod(6)
      integer*4 l,lele,i,ke,lwl,lwt,lwlc,lwtc,irtc,
     $     nextwake,nwak,itab(np),izs(np)
      integer*8 iwpl,iwpt,iwplc,iwptc
      logical*4 sol,out,autophi,seg,krad,wspaccheck
      if(np .le. 0)then
        return
      endif
      if(wspaccheck())then
        return
      endif
      out=.true.
      nextwake=0
      iwpt=0
      iwpl=0
      iwplc=0
      iwptc=0
      sspac1=0.d0
      if(wake)then
        do i=1,nwakep
          if(iwakeelm(i) .ge. lbegin)then
            nwak=i
            nextwake=iwakeelm(i)
            exit
          endif
        enddo
      endif
c      itwb=ilist(1,iwakepold+6)
c     1     +ilist(1,iwakepold+5)*ilist(2,iwakepold+5)
c      ntw=(2*ilist(1,iwakepold+5)+1)*ilist(2,iwakepold+5)
c      itwb1=itwb+ntw*14-1
c      itwb2=itwb+ntw*15-1
c      itwb3=itwb+ntw*16-1
c      itwb4=itwb+ntw*17-1
c      isb=ilist(2,iwakepold+6)
      xlimit=alost*3.d0
      sspac0=pos(lbegin-1)
      call tsetdvfs
      if(rad)then
        allocate(pxr0(np0))
        allocate(pyr0(np0))
        allocate(zr0(np0))
      endif
      allocate(bsi(np0))
      bsi=0.d0
      do l=lbegin,lend
        l_track=l
c        if(abs(x(1))+abs(y(1)) .gt. 0.01d0)then
c        write(*,*)'tturn1-l ',l,np,x(np),y(np),g(np)
c          write(*,*)'tturn1 ',l,np,kptbl(1,1:6)
c        endif
c        call tfmemcheckprint('tturn',l,.false.,irtc)
c        if(trpt .and. codplt)then
c          call ttstat(np,x,px,y,py,z,g,dv,0.d0,
c     1         ' ',sa,ss,0.d0,
c     1         .false.,.false.,0)
c        endif
        if(la .le. 0)then
          call limitnan(x(1:np),-xlimit,xlimit,xlimit)
          call limitnan(px(1:np),-plimit,plimit,plimit)
          call limitnan(y(1:np),-xlimit,xlimit,xlimit)
          call limitnan(py(1:np),-plimit,plimit,plimit)
          call limitnan(z(1:np),-zlimit,zlimit,zlimit)
          call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     1         kptbl,np,n,
     $         0.d0,0.d0,0.d0,0.d0,
     $         -alost,-alost,alost,alost,0.d0,0.d0,0.d0,0.d0)
          if(np .le. 0)then
            go to 9000
          endif
          la=la1
        endif
        call compelc(l,cmp)
        lele=idtype(cmp%id)
        if(sol)then
          if(l .eq. lbegin)then
            call tsol(np,x,px,y,py,z,g,dv,sx,sy,sz,l,lend,
     $           ke,sol,kptbl,la,n,nwak,nextwake,out)
          endif
          if(np .le. 0)then
            go to 9000
          endif
          sol=l .lt. ke
          go to 1020
        endif
        if(l .eq. nextwake)then
          iwpl=abs(kwaketbl(1,nwak))
          if(iwpl .ne. 0)then
            lwl=(ilist(1,iwpl-1)-2)/2
          else
            lwl=0
          endif
          iwpt=abs(kwaketbl(2,nwak))
          if(iwpt .ne. 0)then
            lwt=(ilist(1,iwpt-1)-2)/2
          else
            lwt=0
          endif
          if(lele .ne. icCAVI)then
            fw=(abs(charge)*e*pbunch*anbunch/amass)/np0*.5d0
            kdx=kytbl(kwDX,lele)
            if(kdx .ne. 0)then
              dx=cmp%value(kdx)
            else
              dx=0.d0
            endif
            kdy=kytbl(kwDY,lele)
            if(kdy .ne. 0)then
              dy=cmp%value(kdy)
            else
              dy=0.d0
            endif
            krot=kytbl(kwROT,lele)
            if(krot .ne. 0)then
              rot=cmp%value(krot)
            else
              rot=0.d0
            endif
            call txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           dx,dy,rot,int(anbunch),
     $           fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $           p0,h0,itab,izs,.true.)
          endif
        endif
        if(wspac .or. pspac)then
          sspac=(pos(l)+pos(l+1))*.5d0
          if(sspac .ne. sspac0 .and. wspac)then
c            if(abs(px(1))+abs(py(1)) .gt. 0.01d0)then
c              write(*,*)'tturn1-wspac ',l,np,px(1),py(1),g(1)
c            endif
            cod=twiss(l,0,mfitdx:mfitddp)
            call twspac(np,x,px,y,py,z,g,dv,sx,sy,sz,sspac-sspac0,
     $           cod,beamsize(:,l),l)
c            if(abs(px(1))+abs(py(1)) .gt. 0.01d0)then
c              write(*,*)'tturn1-wspac ',l,np,px(1),py(1),g(1)
c            endif
c            write(*,*)'twspac-end'
          endif
          sspac0=sspac
       endif
       seg=tcheckseg(cmp,lele,al,lsegp,irtc)
       if(irtc .ne. 0)then
         call tffserrorhandle(l,irtc)
         go to 1010
       endif
       select case (lele)
       case (icDRFT)
         if(spac)then
           call spdrift_free(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $          cmp%value(ky_L_DRFT),
     $          cmp%value(ky_RADI_DRFT),n,kptbl)
         else
           if(cmp%value(ky_KIN_DRFT) .eq. 0.d0)then
             do concurrent (i=1:np)
               dpz=pxy2dpz(px(i),py(i))
               al1=al/(1.d0+dpz)
               x(i)=x(i)+px(i)*al1
               y(i)=y(i)+py(i)*al1
               z(i)=z(i)+dpz  *al1-dv(i)*al
             enddo
           else
             x(1:np)=x(1:np)+px(1:np)*al
             y(1:np)=y(1:np)+py(1:np)*al
             z(1:np)=z(1:np)-((px(1:np)**2+py(1:np)**2)*.5d0
     $            +dv(1:np))*al
           endif
           go to 1011
         endif

       case (icBEND)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         if(cmp%value(ky_RANK_BEND) .ne. 0.d0)then
           ak0=cmp%value(ky_ANGL_BEND)
     $          +cmp%value(ky_K0_BEND)
     $          +cmp%value(ky_RANK_BEND)*tgauss()
         else
           ak0=cmp%value(ky_ANGL_BEND)
     $          +cmp%value(ky_K0_BEND)
         endif
         ak1=cmp%value(ky_K1_BEND)
         krad=rad .and. cmp%value(ky_RAD_BEND) .eq. 0.d0 .and.
     $        cmp%value(p_L_BEND) .ne. 0.d0
         if(rad)then
           if(radcod .and. radtaper)then
             rtaper=1.d0-dp0
     $            +(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
             ak0=ak0*rtaper
             ak1=ak1*rtaper
           endif
           if(krad .and. calpol)then
             bsi=0.d0
           endif
         endif

         call tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $        cmp%value(p_L_BEND),ak0,
     $        cmp%value(ky_ANGL_BEND),
     $        cmp%value(p_PSI1_BEND),cmp%value(p_PSI2_BEND),
     1        cmp%value(p_COSPSI1_BEND),cmp%value(p_SINPSI1_BEND),
     1        cmp%value(p_COSPSI2_BEND),cmp%value(p_SINPSI2_BEND),
     1        ak1,
     1        cmp%value(ky_DX_BEND),cmp%value(ky_DY_BEND),
     1        cmp%value(p_THETA_BEND),cmp%value(ky_DROT_BEND),
c     $       cmp%value(p_DPHIX_BEND),cmp%value(p_DPHIY_BEND),
     1        cmp%value(p_COSTHETA_BEND),cmp%value(p_SINTHETA_BEND),
     $        cmp%value(p_FB1_BEND),cmp%value(p_FB2_BEND),
     $        int(cmp%value(ky_FRMD_BEND)),
     $        cmp%value(ky_FRIN_BEND) .eq. 0.d0,
     1        cmp%value(p_COSW_BEND),cmp%value(p_SINW_BEND),
     $        cmp%value(p_SQWH_BEND),cmp%value(p_SINWP1_BEND),
     1        krad,cmp%value(ky_EPS_BEND),.true.,0)

       case (icQUAD)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         rtaper=1.d0
         if(rad .and. radcod .and. radtaper)then
           rtaper=(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
     $          -dp0
         endif
         call tquad(np,x,px,y,py,z,g,dv,sx,sy,sz,al,
     1        cmp%value(ky_K1_QUAD)*rtaper,0.d0,
     $        cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1        cmp%value(ky_ROT_QUAD),
     1        cmp%value(p_THETA2_QUAD),
     1        rad .and. cmp%value(ky_RAD_QUAD) .eq. 0.d0 .and.
     $        al .ne. 0.d0,
     $        cmp%value(ky_CHRO_QUAD) .eq. 0.d0,
     1        cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $        cmp%value(p_AKF1F_QUAD)*rtaper,
     $        cmp%value(p_AKF2F_QUAD)*rtaper,
     $        cmp%value(p_AKF1B_QUAD)*rtaper,
     $        cmp%value(p_AKF2B_QUAD)*rtaper,
     1        int(cmp%value(p_FRMD_QUAD)),cmp%value(ky_EPS_QUAD),
     $        cmp%value(ky_KIN_QUAD) .eq. 0.d0)

       case (icSEXT,icOCTU,icDECA,icDODECA)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         ak1=cmp%value(ky_K_THIN)
         if(rad .and. radcod .and. radtaper)then
           ak1=ak1*((
     $          2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))
     $          *.5d0-dp0)
         endif
         call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $        lele,al,ak1,
     1        cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     1        cmp%value(p_THETA_THIN),
     $        rad .and. cmp%value(ky_RAD_THIN) .eq. 0.d0 .and.
     $        al .ne. 0.d0,
     1        cmp%value(ky_FRIN_THIN) .eq. 0.d0)

       case (icUND)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         call undulator(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $        cmp%value(p_PARAM_UND))
         
       case (icWIG)
         call twig(np,x,px,y,py,z,g,dv,al,cmp%value(ky_BMAX_WIG),
     1        int(cmp%value(ky_PRD_WIG)),
     $        cmp%value(ky_DX_WIG),cmp%value(ky_DY_WIG),
     1        cmp%value(ky_ROT_WIG),cmp%value(p_PARAM_WIG))

       case (icSOL)
         call tsol(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $        l,lend,
     $        ke,sol,kptbl,la,n,nwak,nextwake,out)
         if(np .le. 0)then
           go to 9000
         endif

       case (icST)
         write(*,*)'Use BEND with ANGLE=0 for STEER.'
         call abort
         
       case (icMULT)
         rtaper=1.d0
         if(rad .and. radcod .and. radtaper)then
           rtaper=1.d0-dp0
     $          +(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
         endif
         bz=0.d0
         if(seg)then
           call tmultiseg(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $          cmp,lsegp,bz,rtaper,n,kptbl)
         else
           call tmulti1(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $          cmp,bz,rtaper,n,kptbl)
         endif

       case (icCAVI)
         if(tparacheck(icCAVI,cmp))then
           call tpara(cmp)
         endif
         ak=cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI)
         if(cmp%value(ky_RANV_CAVI) .ne. 0.d0)then
           ak=ak+cmp%value(ky_RANV_CAVI)*tgauss()
         endif
         if(cmp%value(ky_RANP_CAVI) .eq. 0.d0)then
           ph=cmp%value(ky_DPHI_CAVI)
         else
           ph=cmp%value(ky_DPHI_CAVI)+
     $          cmp%value(ky_RANP_CAVI)*tgauss()
         endif
         autophi=cmp%value(ky_APHI_CAVI) .ne. 0.d0
         if(autophi)then
           ph=ph+gettwiss(mfitdz,l)*cmp%value(p_W_CAVI)
         endif
         if(twake .or. lwake)then
           if(l .eq. nextwake)then
             iwplc=max(iwpl-1,0)
             iwptc=max(iwpt-1,0)
             lwlc=lwl
             lwtc=lwt
           else
c     iwplc=abs(ilist(1,lp+ky_LWAK_CAVI))
c     if(iwplc .eq. 0 .or. .not. lwake)then
             lwlc=0
c     else
c     lwlc=(ilist(1,iwplc-1)-2)/2
c     endif
c     iwptc=abs(ilist(1,lp+ky_TWAK_CAVI))
c     if(iwptc .eq. 0 .or. .not. twake)then
             lwtc=0
c     else
c     lwtc=(ilist(1,iwptc-1)-2)/2
c     endif
           endif
           call tcav(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,
     1          cmp%value(p_W_CAVI),cmp%value(ky_PHI_CAVI),ph,
     $          cmp%value(p_VNOMINAL_CAVI),
     $          lwlc,rlist(iwplc+1),lwtc,rlist(iwptc+1),
     1          cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $          cmp%value(ky_ROT_CAVI),
     $          cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $          cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $          cmp%value(ky_FRIN_CAVI) .eq. 0.d0,
     $          int(cmp%value(p_FRMD_CAVI)),autophi)
         else
c           write(*,*)'tturn-tcav-0 ',cmp%value(p_W_CAVI),
c     $          cmp%value(p_VNOMINAL_CAVI)
           call tcav(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,
     1          cmp%value(p_W_CAVI),cmp%value(ky_PHI_CAVI),ph,
     $          cmp%value(p_VNOMINAL_CAVI),
     $          0,0.d0,0,0.d0,
     1          cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $          cmp%value(ky_ROT_CAVI),
     $          cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $          cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $          cmp%value(ky_FRIN_CAVI) .eq. 0.d0,
     $          int(cmp%value(p_FRMD_CAVI)),autophi)
         endif

       case (icTCAV)
         if(rfsw)then
           if(cmp%value(ky_RANK_TCAV) .eq. 0.d0)then
             ak=cmp%value(ky_K0_TCAV)
           else
             ak=cmp%value(ky_K0_TCAV)+cmp%value(ky_RANK_TCAV)*tgauss()
           endif
           if(cmp%value(ky_RANP_TCAV) .eq. 0.d0)then
             ph=cmp%value(ky_PHI_TCAV)
           else
             ph=cmp%value(ky_PHI_TCAV)+cmp%value(ky_RANP_TCAV)*tgauss()
           endif
           harmf=cmp%value(ky_HARM_TCAV)-int(cmp%value(ky_HARM_TCAV))
           ph=ph+harmf*(n-1)*pi2
           call ttcav(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,
     $          cmp%value(ky_HARM_TCAV),ph,cmp%value(ky_FREQ_TCAV),
     1          cmp%value(ky_DX_TCAV),cmp%value(ky_DY_TCAV),
     $          cmp%value(ky_ROT_TCAV),rad .and.
     $          cmp%value(ky_RAD_TCAV) .eq. 0.d0)
         else
           call tdrift_free(np,x,px,y,py,z,dv,al)
         endif

       case (icMAP)
         call temap(np,np0,x,px,y,py,z,g,dv,l,n,kptbl)
         go to 1010

       case (icINS)
         call tins(np,x,px,y,py,z,g,cmp%value(ky_DIR_INS+1))
         go to 1010

       case (icCOORD)
         call tcoord(np,x,px,y,py,z,
     1        cmp%value(ky_DX_COORD),cmp%value(ky_DY_COORD),
     $        cmp%value(ky_DZ_COORD),cmp%value(ky_CHI1_COORD),
     $        cmp%value(ky_CHI2_COORD),cmp%value(ky_CHI3_COORD),
     1        cmp%value(ky_DIR_COORD) .eq. 0.d0)
         go to 1010

       case (icBEAM)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         call beambeam(np,x,px,y,py,z,g,dv,sx,sy,sz,cmp%value(1),
     $        cmp%value(p_PARAM_BEAM),n)
         go to 1020

       case (icProt)
         if(.not. cmp%update)then
           call tpara(cmp)
         endif
         call phsrot(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $        cmp%value(p_PARAM_Prot))
         go to 1010

       case (icSPCH)
         if(pspac) then
           sspac2=(rlist(ifpos+l-1)+rlist(ifpos+l))*.5d0
c     print *,'tturn l sspac2',l,sspac2
           call tpspac(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $          pbunch, amass, p0, h0, sspac2-sspac1,
     $          pspac_nx,pspac_ny,pspac_nz,
     $          pspac_dx,pspac_dy,pspac_dz,
     $          l,lend,pspac_nturn,pspac_nturncalc)
           sspac1=sspac2
         endif
         go to 1020

       case (icMARK)
         go to 1010

       case (icAPRT)
         call tapert1(x,px,y,py,z,g,dv,sx,sy,sz,kptbl,np,n)
         if(np .le. 0)then
           go to 9000
         endif
         la=la1
         go to 1010

       case default
         go to 1010
       end select
 1020   la=la-1
 1011   if(radlight)then
          if(lele .eq. icBEND)then
          elseif(lele .eq. icMULT .and.
     $           cmp%value(ky_ANGL_MULT) .ne. 0.d0)then
          elseif(lele .lt. 32 .and. al .ne. 0.d0)then
            call tlstore(np,x,y,z,dv,0.d0,al,0.d0,0.d0,
     $           p0/h0*c,dvfs,.true.)
          endif
        endif
 1010   continue
        if(l .eq. nextwake)then
          if(lele .ne. icCAVI)then
            call txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           dx,dy,rot,int(anbunch),
     $           fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $           p0,h0,itab,izs,.false.)
          endif
          nwak=nwak+1
          if(nwak .gt. nwakep)then
            nextwake=0
          else
            nextwake=iwakeelm(nwak)
          endif
        endif
      enddo
      call limitnan(x(1:np),-xlimit,xlimit,xlimit)
      call limitnan(px(1:np),-plimit,plimit,plimit)
      call limitnan(y(1:np),-xlimit,xlimit,xlimit)
      call limitnan(py(1:np),-plimit,plimit,plimit)
      call limitnan(z(1:np),-zlimit,zlimit,zlimit)
      call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     1     kptbl,np,n,
     $     0.d0,0.d0,0.d0,0.d0,
     $     -alost,-alost,alost,alost,0.d0,0.d0,0.d0,0.d0)
      if(np .le. 0)then
        go to 9000
      endif
      la=la1
c      call tfmemcheckprint('tturn',1,.false.,irtc)
      if(wspac.or.pspac)then
        sspac=pos(lend)
        if(sspac .ne. sspac0)then
           if(wspac)then
             cod=twiss(lend,0,mfitdx:mfitddp)
             call twspac(np,x,px,y,py,z,g,dv,sx,sy,sz,sspac-sspac0,
     $            cod,beamsize(:,lend),lend)
           endif
           if(pspac) then
              call tpspac(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             pbunch, amass, p0, h0, sspac-sspac0,
     $             pspac_nx,pspac_ny,pspac_nz,
     $             pspac_dx,pspac_dy,pspac_dz,
     $            lend,lend,pspac_nturn,pspac_nturncalc)
           endif
        endif
      endif
 9000 deallocate(bsi)
      if(rad)then
        deallocate(zr0)
        deallocate(pyr0)
        deallocate(pxr0)
      endif
      return
      end

      subroutine tt6621(s,b)
      implicit none
      real*8 s(6,6),b(21)
      b( 1)=s(1,1)
      b( 2)=s(2,1)
      b( 3)=s(2,2)
      b( 4)=s(3,1)
      b( 5)=s(3,2)
      b( 6)=s(3,3)
      b( 7)=s(4,1)
      b( 8)=s(4,2)
      b( 9)=s(4,3)
      b(10)=s(4,4)
      b(11)=s(5,1)
      b(12)=s(5,2)
      b(13)=s(5,3)
      b(14)=s(5,4)
      b(15)=s(5,5)
      b(16)=s(6,1)
      b(17)=s(6,2)
      b(18)=s(6,3)
      b(19)=s(6,4)
      b(20)=s(6,5)
      b(21)=s(6,6)
      return
      end

      subroutine tmultiseg(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     cmp,lsegp,bz,rtaper,n,kptbl)
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_seg
      use tspin
      implicit none
      type (sad_comp) :: cmp
      integer*4 np,n
      integer*4 kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),
     $     g(np0),dv(np0)
      real*8 sx(np0),sy(np0),sz(np0)
      real*8 bz,rtaper
      type (sad_dlist) :: lsegp
      type (sad_dlist), pointer :: lal,lk
      type (sad_rlist), pointer :: lak,lkv
      real*8 :: rsave(cmp%ncomp2)
      integer*4 i,nseg,i1,i2,istep,k,k1,k2,nk
      integer*8 kk
      integer*4 , parameter :: nc=ky_PROF_MULT-1
      rsave(1:nc)=cmp%value(1:nc)
      nk=lsegp%nl
      call descr_sad(lsegp%dbody(1),lal)
      call descr_sad(lal%dbody(2),lak)
      nseg=lak%nl
      if(cmp%ori)then
        i1=1
        i2=nseg
        istep=1
      else
        i1=nseg
        i2=1
        istep=-1
      endif
      do i=i1,i2,istep
        do concurrent (k=1:nc)
          if(integv(k,icMULT))then
            cmp%value(k)=rsave(k)*lak%rbody(i)
          endif
        enddo
        do k=1,nk
          call descr_sad(lsegp%dbody(k),lk)
          call descr_sad(lk%dbody(2),lkv)
          kk=ktfaddr(lsegp%dbody(k))
          k1=ilist(1,kk+1)
          k2=ilist(2,kk+1)
          if(k1 .eq. k2)then
            cmp%value(k1)=0.d0
          endif
          cmp%value(k1)=cmp%value(k1)+rsave(k2)*lkv%rbody(i)
        enddo
        cmp%update=.false.
        call tmulti1(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       cmp,bz,rtaper,n,kptbl)
      enddo
      cmp%value(1:nc)=rsave(1:nc)
      return
      end

      subroutine tmulti1(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     cmp,bz,rtaper,n,kptbl)
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use tparastat
      use tspin
      implicit none
      type (sad_comp) :: cmp
      integer*4 np,n
      integer*4 kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      real*8 sx(np0),sy(np0),sz(np0)
      real*8 ph,bz,rtaper
      logical*4 autophi
      if(tparacheck(icMULT,cmp))then
        call tpara(cmp)
      endif
      autophi=cmp%value(ky_APHI_MULT) .ne. 0.d0
      ph=cmp%value(ky_DPHI_MULT)
      if(autophi)then
        ph=ph+gettwiss(mfitdz,l_track)*cmp%value(p_W_MULT)
      endif
      call tmulti(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     cmp%value(p_L_MULT),cmp%value(ky_K0_MULT),
     $     bz,cmp%value(p_ANGL_MULT),
     $     cmp%value(p_PSI1_MULT),cmp%value(p_PSI2_MULT),
     1     cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $     cmp%value(ky_DZ_MULT),
     $     cmp%value(p_CHI1_MULT),cmp%value(p_CHI2_MULT),
     $     cmp%value(ky_ROT_MULT),
     $     cmp%value(ky_DROT_MULT),
     $     cmp%value(p_THETA2_MULT),
     $     cmp%value(p_CR1_MULT),
     $     cmp%value(ky_EPS_MULT),
     $     rad .and. cmp%value(ky_RAD_MULT) .eq. 0.d0 .and.
     $     cmp%value(p_L_MULT) .ne. 0.d0,
     $     cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $     cmp%value(p_AKF1F_MULT)*rtaper,
     $     cmp%value(p_AKF2F_MULT)*rtaper,
     $     cmp%value(p_AKF1B_MULT)*rtaper,
     $     cmp%value(p_AKF2B_MULT)*rtaper,
     $     int(cmp%value(p_FRMD_MULT)),
     $     cmp%value(p_FB1_MULT),cmp%value(p_FB2_MULT),
     $     cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT),
     $     cmp%value(p_W_MULT),
     $     cmp%value(ky_PHI_MULT),ph,cmp%value(p_VNOMINAL_MULT),
     $     cmp%value(ky_RADI_MULT),rtaper,autophi,
     $     n,kptbl)
      return
      end
