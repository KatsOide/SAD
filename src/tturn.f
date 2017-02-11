      subroutine tturn(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n)
      use tfstk
      use tmacro
      implicit none
      real*8 plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,kptbl(np0,6)
      integer*8 latt(nlat)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      logical*4 normal
      call tturn0(np,latt,1,nlat,x,px,y,py,z,g,dv,pz,kptbl,n,normal)
      return
      end

      subroutine tturn0(np,latt,lb,le,x,px,y,py,z,g,dv,pz,
     $     kptbl,n,normal)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs, only:ffs_bound
      implicit none
      type (ffs_bound) fbound
      real*8 plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,la,ls,nvar,lb,le
c      integer*4 isb,itwb,itwb1,itwb2,itwb3,itwb4,ntw
      real*8 pgev00
      integer*4 kptbl(np0,6)
      integer*8 latt(nlat)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sa(6),ss(6,6),vsave(100)
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
      write(*,*)'tturn0 ',fbound%lb,fbound%fb,fbound%le,fbound%fe
      if(fbound%lb .eq. fbound%le)then
        call qfracsave(fbound%lb,vsave,nvar,.true.)
        call qfraccomp(fbound%lb,fbound%fb,fbound%fe,ideal,chg)
        call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $       sol,la,fbound%lb,fbound%lb)
        if(chg)then
          call qfracsave(fbound%lb,vsave,nvar,.false.)
        endif
      else
        if(fbound%fb .ne. 0.d0)then
          call qfracsave(fbound%lb,vsave,nvar,.true.)
          call qfraccomp(fbound%lb,fbound%fb,1.d0,ideal,chg)
          call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $         sol,la,fbound%lb,fbound%lb)
          if(chg)then
            call qfracsave(fbound%lb,vsave,nvar,.false.)
          endif
          ls=fbound%lb+1
        else
          ls=fbound%lb
        endif
        if(fbound%le .gt. ls+1)then
          call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $         sol,la,ls,fbound%le-1)
        endif
        if(fbound%fe .ne. 0.d0)then
          call qfracsave(fbound%le,vsave,nvar,.true.)
          call qfraccomp(fbound%le,0.d0,fbound%fe,ideal,chg)
          call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $         sol,la,fbound%le,fbound%le)
          if(chg)then
            call qfracsave(fbound%le,vsave,nvar,.false.)
          endif
        endif
      endif
      if(trpt .and. codplt)then
        call ttstat(np,x,px,y,py,z,g,dv,0.d0,
     1       ' ',sa,ss,0.d0,
     1       .false.,.false.,0)
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

      subroutine tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     sol,la,lbegin,lend)
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use ffs_wake
      use sad_main
      use ffs_pointer, only: direlc,compelc
      implicit none
      integer*4 la1
      parameter (la1=15)
      real*8 xlimit,plimit,zlimit,vmax,ampmax
      parameter (plimit=0.7d0,zlimit=1.d10,ampmax=0.9999d0)
      parameter (vmax=.9d0)
      type (sad_comp), pointer:: cmp
      integer*4 np,n,la,lbegin,lend,kdx,kdy,krot,kyl
      integer*4 kptbl(np0,6)
      integer*8 latt(nlat)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),
     $     g(np0),dv(np0),pz(np0),dpz,al1
      real*8 sa(6),ss(6,6),phi,bz,harm,w,
     $     al,ak0,ak1,psi1,psi2,tgauss,ph,harmf,
     $     sspac0,sspac,fw,dx,dy,rot,sspac1,sspac2,
     $     fb1,fb2,chi1,chi2,ak,rtaper,vnominal
      integer*4 l,lele,i,mfr,ke,lwl,lwt,lwlc,lwtc,
     $     nextwake,nwak,itab(np),izs(np)
      integer*8 iwpl,iwpt,iwplc,iwptc,itp
      logical*4 sol,out,autophi
      if(np .le. 0)then
        return
      endif
      call wspaccheck
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
      sspac0=rlist(ifpos+lbegin-1)
      call tsetdvfs
      do l=lbegin,lend
c        if(mod(l,100) .eq. 1)then
c          write(*,*)'ttrun1 ',l,z(1),dvfs
c        endif
        if(trpt .and. codplt)then
          call ttstat(np,x,px,y,py,z,g,dv,0.d0,
     1         ' ',sa,ss,0.d0,
     1         .false.,.false.,0)
        endif
        if(la .le. 0)then
          call tapert(l,latt,x,px,y,py,z,g,dv,pz,
     1         kptbl,np,n,
     $         0.d0,0.d0,0.d0,0.d0,
     $         -alost,-alost,alost,alost,0.d0,0.d0,0.d0,0.d0)
          if(np .le. 0)then
            return
          endif
          la=la1
        endif
        call compelc(l,cmp)
        lele=idtype(cmp%id)
        if(sol)then
          if(l .eq. lbegin)then
            call tsol(np,x,px,y,py,z,g,dv,pz,latt,l,lend,
     $           ke,sol,kptbl,la,n,nwak,nextwake,out)
          endif
          if(np .le. 0)then
            return
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
            call txwake(np,x,px,y,py,z,g,dv,
     $           dx,dy,rot,int(anbunch),
     $           fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $           p0,h0,itab,izs,.true.)
          endif
        endif
        if(wspac .or. pspac)then
          sspac=(rlist(ifpos+l-1)+rlist(ifpos+l))*.5d0
          if(sspac .ne. sspac0)then
             if(wspac) then
               call twspac(np,x,px,y,py,z,g,dv,pz,sspac-sspac0,
     $              gettwiss(mfitdx,l),
c     $              rlist(iftwis+((mfitdx-1)*(2*ndim+1)+ndim)*nlat
c     $              +l-1),
     $              rlist(ifsize+(l-1)*21))
             endif
          endif
          sspac0=sspac
       endif
       itp=cmp%param
       kyl=kytbl(kwL,lele)
       if(kyl .ne. 0.d0)then
         al=cmp%value(kyl)
       else
         al=0.d0
       endif
c        if(l .eq. lbegin .or. l .eq. lend)then
c          write(*,'(1X,I5,I8,I5,I6,1X,A,1P7G11.3)')
c     1       l,itp,lele,lp,pname(ilist(2,latt(l)))(1:8),
c     $       x(4),px(4),y(4),py(4),
c     1       z(4),g(4)*(2.d0+g(4))
c        endif
       go to (1100,1200,1010,1400,1010,1600,1010,1600,1010,1600,
     1      1010,1600,1010,1010,1010,1010,1010,1800,1900,2000,
     1      2100,2200,1010,1010,1010,1010,1010,1010,1010,1010,
     1      3100,3200,3300,3400,3500,3600,3700,3800,1010,1010,
     1      4100,1010,4300),lele
       go to 1010
 1100   if(spac)then
          call spdrift_free(np,x,px,y,py,z,g,dv,pz,
     $         cmp%value(ky_L_DRFT),
     $         cmp%value(ky_RADI_DRFT),n,l,latt,kptbl)
          go to 1020
        else
          if(cmp%value(kytbl(kwKIN,lele)) .eq. 0.d0)then
c$$$            do i=1,np
c$$$              s=min(ampmax,px(i)**2+py(i)**2)
c$$$              dpz=sqrt1(-s)
c$$$              pzi=1.d0+dpz
c$$$              al1=al/pzi
c$$$              x(i)=x(i)+px(i)*al1
c$$$              y(i)=y(i)+py(i)*al1
c$$$              z(i)=z(i)+dpz*al1-dv(i)*al
c$$$            enddo
c            a=px(1:np)**2+py(1:np)**2
c            dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
c            al1=al/(1.d0+dpz)
c            x(1:np)=x(1:np)+px(1:np)*al1
c            y(1:np)=y(1:np)+py(1:np)*al1
c            z(1:np)=z(1:np)+dpz  *al1-dv(1:np)*al
            do i=1,np
              dpz=pxy2dpz(px(i),py(i))
c              dpz=s*(-.5d0-s*(.125d0+s*.0625d0))
c              dpz=(dpz**2-s)/(2.d0+2.d0*dpz)
c              dpz=(dpz**2-s)/(2.d0+2.d0*dpz)
              al1=al/(1.d0+dpz)
              x(i)=x(i)+px(i)*al1
              y(i)=y(i)+py(i)*al1
              z(i)=z(i)+dpz  *al1-dv(i)*al
            enddo
          else
c            do i=1,np
              x(1:np)=x(1:np)+px(1:np)*al
              y(1:np)=y(1:np)+py(1:np)*al
              z(1:np)=z(1:np)-((px(1:np)**2+py(1:np)**2)*.5d0
     $             +dv(1:np))*al
c            enddo
          endif
          go to 1011
        endif
 1200   continue
        al=rlist(itp)
        if(cmp%value(ky_RANK_BEND) .ne. 0.d0)then
          ak0=cmp%value(ky_ANGL_BEND)
     $         +cmp%value(ky_K0_BEND)
     $         +cmp%value(ky_RANK_BEND)*tgauss()
        else
          ak0=cmp%value(ky_ANGL_BEND)
     $         +cmp%value(ky_K0_BEND)
        endif
        if(direlc(l) .gt. 0.d0)then
          fb1=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB1_BEND)
          fb2=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB2_BEND)
        else
          fb2=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB1_BEND)
          fb1=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB2_BEND)
        endif
        ak1=cmp%value(ky_K1_BEND)
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          ak0=ak0*rtaper
          ak1=ak1*rtaper
        endif
        call tbend(np,x,px,y,py,z,g,dv,pz,l,al,ak0,
     $       cmp%value(ky_ANGL_BEND),
     1       rlist(itp+1),rlist(itp+2),rlist(itp+3),rlist(itp+4),
     1       ak1,
     1       cmp%value(ky_DX_BEND),
     $       cmp%value(ky_DY_BEND),
     1       rlist(itp+13),rlist(itp+11),rlist(itp+12),
     1       rlist(itp+5),rlist(itp+6),
     $       fb1,fb2,
     $       int(cmp%value(ky_FRMD_BEND)),
     $       cmp%value(ky_FRIN_BEND) .eq. 0.d0,
     1       rlist(itp+7),rlist(itp+8),rlist(itp+9),rlist(itp+10),
     1       cmp%value(ky_RAD_BEND) .eq. 0.d0,
     $       0.d0,al,al,
     1       cmp%value(ky_EPS_BEND))
        go to 1020
 1400   continue
        if(direlc(l) .gt. 0.d0)then
          mfr=int(cmp%value(ky_FRMD_QUAD))
        else
          mfr=int(cmp%value(ky_FRMD_QUAD))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tquad(np,x,px,y,py,z,g,dv,pz,l,al,
     1       cmp%value(ky_K1_QUAD)*rtaper,
     $       cmp%value(ky_DX_QUAD),
     $       cmp%value(ky_DY_QUAD),
     1       rlist(itp+4),rlist(itp+2),rlist(itp+3),
     1       cmp%value(ky_RAD_QUAD),
     $       cmp%value(ky_CHRO_QUAD) .eq. 0.d0,
     1       cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $       rlist(itp+6)*rtaper,rlist(itp+7)*rtaper,
     $       rlist(itp+8)*rtaper,rlist(itp+9)*rtaper,
     1       mfr,cmp%value(ky_EPS_QUAD),
     $       cmp%value(ky_KIN_QUAD) .eq. 0.d0)
        go to 1020
 1600   ak1=cmp%value(ky_K_THIN)
        if(rad .and. radcod .and. radtaper)then
          ak1=ak1*(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tthin(np,x,px,y,py,z,g,dv,pz,lele,l,al,ak1,
     1       cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     $       rlist(itp+4),rlist(itp+2),rlist(itp+3),
     $       cmp%value(ky_RAD_THIN),
     1       cmp%value(ky_FRIN_THIN) .eq. 0.d0)
        go to 1020
 1800   call undulator(np,x,px,y,py,z,g,dv,pz,rlist(itp))
        go to 1020
 1900   call twig(np,x,px,y,py,z,g,dv,al,cmp%value(ky_BMAX_WIG),
     1       int(cmp%value(ky_PRD_WIG)),
     $       cmp%value(ky_DX_WIG),cmp%value(ky_DY_WIG),
     1       cmp%value(ky_ROT_WIG),itp)
        go to 1020
 2000   call tsol(np,x,px,y,py,z,g,dv,pz,latt,l,lend,
     $       ke,sol,kptbl,la,n,nwak,nextwake,out)
        if(np .le. 0)then
          return
        endif
        go to 1020
 2100   write(*,*)'Use BEND with ANGLE=0 for STEER.'
        call abort
 2200   phi=cmp%value(ky_ANGL_MULT)
        mfr=nint(cmp%value(ky_FRMD_MULT))
        if(direlc(l) .gt. 0.d0)then
          psi1=phi*cmp%value(ky_E1_MULT)
     $         +cmp%value(ky_AE1_MULT)
          psi2=phi*cmp%value(ky_E2_MULT)
     $         +cmp%value(ky_AE2_MULT)
          fb1=cmp%value(ky_FB1_MULT)
          fb2=cmp%value(ky_FB2_MULT)
          chi1=cmp%value(ky_CHI1_MULT)
          chi2=cmp%value(ky_CHI2_MULT)
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=phi*cmp%value(ky_E2_MULT)
     $         +cmp%value(ky_AE2_MULT)
          psi2=phi*cmp%value(ky_E1_MULT)
     $         +cmp%value(ky_AE1_MULT)
          fb2=cmp%value(ky_FB1_MULT)
          fb1=cmp%value(ky_FB2_MULT)
          chi1=-cmp%value(ky_CHI1_MULT)
          chi2=-cmp%value(ky_CHI2_MULT)
        endif
        bz=0.d0
        if(trpt)then
          vnominal=cmp%value(ky_VOLT_MULT)/amass*abs(charge)
     $         *sin(-cmp%value(ky_PHI_MULT)*sign(1.d0,charge))
        else
          vnominal=0.d0
        endif
        harm=cmp%value(ky_HARM_MULT)
        if(harm .eq. 0.d0)then
          w=pi2*cmp%value(ky_FREQ_MULT)/c
        else
          w=omega0*harm/c
        endif
        autophi=cmp%value(ky_APHI_MULT) .ne. 0.d0
        ph=cmp%value(ky_DPHI_MULT)
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tmulti(np,x,px,y,py,z,g,dv,pz,al,
     $       cmp%value(ky_K0_MULT),
     $       bz,
     $       phi,psi1,psi2,
     1       cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $       cmp%value(ky_DZ_MULT),
     $       chi1,chi2,cmp%value(ky_ROT_MULT),
     $       cmp%value(ky_DROT_MULT),
     $       cmp%value(ky_EPS_MULT),
     $       cmp%value(ky_RAD_MULT) .eq. 0.d0,
     $       cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $       rlist(itp+1)*rtaper,rlist(itp+2)*rtaper,
     $       rlist(itp+3)*rtaper,rlist(itp+4)*rtaper,
     $       mfr,fb1,fb2,
     $       cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT),
     $       w,
     $       cmp%value(ky_PHI_MULT),ph,vnominal,
     $       cmp%value(ky_RADI_MULT),rtaper,autophi,
     $       n,l,latt,kptbl)
        go to 1020
 3100   harm=cmp%value(ky_HARM_CAVI)
        if(harm .eq. 0.d0)then
          w=pi2*cmp%value(ky_FREQ_CAVI)/c
        else
          w=omega0*harm/c
        endif
        if(trpt)then
          vnominal=cmp%value(ky_VOLT_CAVI)/amass*abs(charge)
     $         *sin(-cmp%value(ky_PHI_CAVI)*sign(1.d0,charge))
        else
          vnominal=0.d0
        endif
        ak=cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI)
        if(cmp%value(ky_RANV_CAVI) .ne. 0.d0)then
          ak=ak+cmp%value(ky_RANV_CAVI)*tgauss()
        endif
        if(cmp%value(ky_RANP_CAVI) .eq. 0.d0)then
          ph=cmp%value(ky_DPHI_CAVI)
        else
          ph=cmp%value(ky_DPHI_CAVI)+
     $         cmp%value(ky_RANP_CAVI)*tgauss()
        endif
        autophi=cmp%value(ky_APHI_CAVI) .ne. 0.d0
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
        endif
        mfr=nint(cmp%value(ky_FRMD_CAVI))
        if(direlc(l) .gt. 0.d0)then
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        if(twake .or. lwake)then
          if(l .eq. nextwake)then
            iwplc=max(iwpl-1,0)
            iwptc=max(iwpt-1,0)
            lwlc=lwl
            lwtc=lwt
          else
c            iwplc=abs(ilist(1,lp+ky_LWAK_CAVI))
c            if(iwplc .eq. 0 .or. .not. lwake)then
              lwlc=0
c            else
c              lwlc=(ilist(1,iwplc-1)-2)/2
c            endif
c            iwptc=abs(ilist(1,lp+ky_TWAK_CAVI))
c            if(iwptc .eq. 0 .or. .not. twake)then
              lwtc=0
c            else
c              lwtc=(ilist(1,iwptc-1)-2)/2
c            endif
          endif
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,cmp%value(ky_PHI_CAVI),ph,vnominal,
     $         lwlc,rlist(iwplc+1),lwtc,rlist(iwptc+1),
     1         cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $         cmp%value(ky_ROT_CAVI),
     $         cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $         cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $         cmp%value(ky_FRIN_CAVI) .eq. 0.d0,mfr,autophi)
        else
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,cmp%value(ky_PHI_CAVI),ph,vnominal,
     $         0,0.d0,0,0.d0,
     1         cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $         cmp%value(ky_ROT_CAVI),
     $         cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $         cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $         cmp%value(ky_FRIN_CAVI) .eq. 0.d0,mfr,autophi)
        endif
        go to 1020
 3200   if(rfsw)then
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
          call ttcav(np,x,px,y,py,z,g,dv,pz,al,ak,
     $         cmp%value(ky_HARM_TCAV),ph,cmp%value(ky_FREQ_TCAV),
     1         cmp%value(ky_DX_TCAV),cmp%value(ky_DY_TCAV),
     $         cmp%value(ky_ROT_TCAV))
        else
          call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        endif
        go to 1020
 3300   continue
c        write(*,*)'ttrun1-temaxp ',np,np0,i,n
        call temap(np,np0,x,px,y,py,z,g,dv,l,n,kptbl)
        go to 1010
 3400   call tins(np,x,px,y,py,z,g,cmp%value(ky_DIR_INS+1))
        go to 1010
 3500   call tcoord(np,x,px,y,py,z,
     1       cmp%value(ky_DX_COORD),cmp%value(ky_DY_COORD),
     $       cmp%value(ky_DZ_COORD),cmp%value(ky_CHI1_COORD),
     $       cmp%value(ky_CHI2_COORD),cmp%value(ky_CHI3_COORD),
     1       cmp%value(ky_DIR_COORD) .eq. 0.d0)
        go to 1010
 3600   call beambeam(np,x,px,y,py,z,g,dv,pz,cmp%value(1),rlist(itp),n)
c        write(*,*)'beambeam-end ',cmp%param
        go to 1010
 3700   call phsrot(np,x,px,y,py,z,g,dv,pz,rlist(itp))
        go to 1010
 3800   if(pspac) then
           sspac2=(rlist(ifpos+l-1)+rlist(ifpos+l))*.5d0
c           print *,'tturn l sspac2',l,sspac2
           call tpspac(np,x,px,y,py,z,g,dv,pz,
     $          pbunch, amass, p0, h0, sspac2-sspac1,
     $          pspac_nx,pspac_ny,pspac_nz,
     $          pspac_dx,pspac_dy,pspac_dz,l,lend,
     $          pspac_nturn,pspac_nturncalc)
           sspac1=sspac2
        endif
        go to 1010
 4100   go to 1010
 4300   call tapert1(l,latt,x,px,y,py,z,g,dv,pz,
     1       kptbl,np,n)
        if(np .le. 0)then
          return
        endif
        la=la1
        go to 1010
 1020   la=la-1
        do i=1,np
          x(i)=min(xlimit,max(-xlimit,x(i)))
          y(i)=min(xlimit,max(-xlimit,y(i)))
          px(i)=min(plimit,max(-plimit,px(i)))
          py(i)=min(plimit,max(-plimit,py(i)))
          z(i)=min(zlimit,max(-zlimit,z(i)))
        enddo
 1011   if(radlight)then
          if(lele .eq. icBEND)then
          elseif(lele .eq. icMULT .and.
     $           cmp%value(ky_ANGL_MULT) .ne. 0.d0)then
          elseif(lele .lt. 32 .and. al .gt. 0.d0)then
            call tlstore(np,x,y,z,dv,0.d0,al,0.d0,0.d0,
     $           p0/h0*c,dvfs,.true.)
          endif
        endif
 1010   continue
        if(l .eq. nextwake)then
          if(lele .ne. icCAVI)then
            call txwake(np,x,px,y,py,z,g,dv,
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
      if(wspac.or.pspac)then
        sspac=rlist(ifpos+lend-1)
        if(sspac .ne. sspac0)then
           if(wspac)then
              call twspac(np,x,px,y,py,z,g,dv,pz,sspac-sspac0,
     $            gettwiss(mfitdx,lend),
c     $             rlist(iftwis+((mfitdx-1)*(2*ndim+1)+ndim)*nlat
c     $             +lend-1),
     $             rlist(ifsize+(lend-1)*21))
           endif
           if(pspac) then
              call tpspac(np,x,px,y,py,z,g,dv,pz,
     $             pbunch, amass, p0, h0, sspac-sspac0,
     $             pspac_nx,pspac_ny,pspac_nz,
     $             pspac_dx,pspac_dy,pspac_dz)
           endif
        endif
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
