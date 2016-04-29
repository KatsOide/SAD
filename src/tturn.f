      subroutine tturn(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      real*8 plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,latt(2,nlat),kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      call tturn0(np,latt,1,nlat,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     .false.,0,0,0)
      return
      end

      subroutine tturn0(np,latt,lb,le,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     wake,iwakeelm,kwaketbl,nwakep)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      real*8 plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,la,ls,lbegin,lend,nvar,lb,le
      integer*4 isb,itwb,itwb1,itwb2,itwb3,itwb4,ntw
      real*8 frbegin,frend,pgev00
      integer*4 latt(2,nlat),kptbl(np0,6),nwakep,iwakeelm(nwakep)
      integer*8 kwaketbl(2,nwakep)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sa(6),ss(6,6),vsave(100)
      logical*4 sol,chg,wake,tfinsol
      pgev00=pgev
      sol=tfinsol(lb)
      novfl=0
      la=0
      if(radlight)then
        call tlstore(np,x,y,z,dv,0.d0,0.d0,0.d0,0.d0,
     $       p0/h0*c,dvfs,.true.)
      endif
      call tffsbound1(nlat,latt,lb,le,lbegin,frbegin,lend,frend)
      if(frbegin .ne. 0.d0)then
        call qfracsave(latt(1,lbegin),vsave,nvar,ideal,.true.)
        call qfraccomp(latt(1,lbegin),frbegin,1.d0,ideal,chg)
        call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $       sol,la,lbegin,lbegin,
     $       wake,iwakeelm,kwaketbl,nwakep)
        if(chg)then
          call qfracsave(latt(1,lbegin),vsave,nvar,ideal,.false.)
        endif
        ls=lbegin+1
      else
        ls=lbegin
      endif
      call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     sol,la,ls,lend-1,
     $     wake,iwakeelm,kwaketbl,nwakep)
      if(frend .ne. 0.d0)then
        call qfracsave(latt(1,lend),vsave,nvar,ideal,.true.)
        call qfraccomp(latt(1,lend),0.d0,frend,ideal,chg)
        call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $       sol,la,lend,lend,
     $     wake,iwakeelm,kwaketbl,nwakep)
        if(chg)then
          call qfracsave(latt(1,lend),vsave,nvar,ideal,.false.)
        endif
      endif
      if(trpt .and. codplt)then
        call ttstat(np,x,px,y,py,z,g,dv,rlist(ilist(2,iwakepold+4)),
     1       ' ',
     1       sa,ss,0.d0,
     1       .false.,.false.,0)
        itwb=ilist(1,iwakepold+6)
     1       +ilist(1,iwakepold+5)*ilist(2,iwakepold+5)
        ntw=(2*ilist(1,iwakepold+5)+1)*ilist(2,iwakepold+5)
        itwb1=itwb+ntw*14-1
        itwb2=itwb+ntw*15-1
        itwb3=itwb+ntw*16-1
        itwb4=itwb+ntw*17-1
        isb=ilist(2,iwakepold+6)
        rlist(itwb1+nlat)=sa(1)
        rlist(itwb2+nlat)=sa(2)
        rlist(itwb3+nlat)=sa(3)
        rlist(itwb4+nlat)=sa(4)
        call tt6621(ss,rlist(isb+21*(nlat-1)))
      endif
      if(rfsw)then
      else
        call tclr(z,np)
      endif
      if(pgev .ne. pgev00)then
        pgev=pgev00
        call tphyzp
      endif
      return
      end

      subroutine tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     sol,la,lbegin,lend,
     $     wake,iwakeelm,kwaketbl,nwakep)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 la1
      parameter (la1=15)
      real*8 xlimit,plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,la,lbegin,lend,kdx,kdy,krot
      integer*4 latt(2,nlat),kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sa(6),ss(6,6),phi,bz,harm,w,
     $     al,a,dpz,al1,ak0,ak1,psi1,psi2,tgauss,ph,harmf,
     $     sspac0,sspac,fw,dx,dy,rot,sspac1,sspac2,
     $     fb1,fb2,chi1,chi2,ak,rtaper
      integer*4 itwb,ntw,itwb1,itwb2,itwb3,itwb4,
     $     isb,l,lele,lp,itp,i,mfr,ke,lwl,lwt,lwlc,lwtc,
     $     nwakep,iwakeelm(nwakep),
     $     nextwake,nwak,itab(np),izs(np)
      integer*8 kwaketbl(2,nwakep),iwpl,iwpt,iwplc,iwptc
      logical*4 sol,out,wake,autophi
      if(np .le. 0)then
        return
      endif
      call wspaccheck
      out=.true.
      nextwake=0
      if(wake)then
        do i=1,nwakep
          if(iwakeelm(i) .ge. lbegin)then
            nwak=i
            nextwake=iwakeelm(i)
            exit
          endif
        enddo
      endif
      itwb=ilist(1,iwakepold+6)
     1     +ilist(1,iwakepold+5)*ilist(2,iwakepold+5)
      ntw=(2*ilist(1,iwakepold+5)+1)*ilist(2,iwakepold+5)
      itwb1=itwb+ntw*14-1
      itwb2=itwb+ntw*15-1
      itwb3=itwb+ntw*16-1
      itwb4=itwb+ntw*17-1
      isb=ilist(2,iwakepold+6)
      xlimit=alost*3.d0
      sspac0=rlist(ifpos+lbegin-1)
      call tsetdvfs
      do l=lbegin,lend
c        if(mod(l,100) .eq. 1)then
c          write(*,*)'ttrun1 ',l,z(1),dvfs
c        endif
        if(trpt .and. codplt)then
          call ttstat(np,x,px,y,py,z,g,dv,rlist(ilist(2,iwakepold+4)),
     1         ' ',
     1         sa,ss,0.d0,
     1         .false.,.false.,0)
          rlist(itwb1+l)=sa(1)
          rlist(itwb2+l)=sa(2)
          rlist(itwb3+l)=sa(3)
          rlist(itwb4+l)=sa(4)
          call tt6621(ss,rlist(isb+21*(l-1)))
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
        lele=idtype(latt(1,l))
        if(sol)then
          if(l .eq. lbegin)then
            call tsol(np,x,px,y,py,z,g,dv,pz,latt,l,lend,
     $           ke,sol,kptbl,la,n,
     $           nwakep,iwakeelm,kwaketbl,nwak,nextwake,out)
          endif
          if(np .le. 0)then
            return
          endif
          sol=l .lt. ke
          go to 1020
        endif
        lp=latt(2,l)
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
              dx=rlist(lp+kdx)
            else
              dx=0.d0
            endif
            kdy=kytbl(kwDY,lele)
            if(kdy .ne. 0)then
              dy=rlist(lp+kdy)
            else
              dy=0.d0
            endif
            krot=kytbl(kwROT,lele)
            if(krot .ne. 0)then
              rot=rlist(lp+krot)
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
!             if(pspac) then
!                call tpspac(np,x,px,y,py,z,g,dv,pz,
!     $           pbunch, amass, p0, h0, sspac-sspac0,
!     $           pspac_nx,pspac_ny,pspac_nz,
!     $           pspac_dx,pspac_dy,pspac_dz)
!             endif
          endif
          sspac0=sspac
       endif
       itp=ilist(2,lp)
       al=rlist(lp+1)
c        if(l .eq. lbegin .or. l .eq. lend)then
c          write(*,'(1X,I5,I8,I5,I6,1X,A,1P7G11.3)')
c     1       l,itp,lele,lp,pname(latt(1,l))(1:8),
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
          call spdrift_free(np,x,px,y,py,z,g,dv,pz,rlist(lp+1),
     $         rlist(lp+kytbl(kwRADI,icDRFT)),n,l,latt,kptbl)
          go to 1020
        else
          if(rlist(lp+2) .eq. 0.d0)then
            do i=1,np
              a=px(i)**2+py(i)**2
              dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
              dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
              dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
              al1=al/(1.d0+dpz)
              x(i)=min(xlimit,max(-xlimit,x(i)+px(i)*al1))
              y(i)=min(xlimit,max(-xlimit,y(i)+py(i)*al1))
              z(i)=z(i)+dpz  *al1-dv(i)*al
            enddo
          else
            do i=1,np
              x(i)=x(i)+px(i)*al
              y(i)=y(i)+py(i)*al
              z(i)=z(i)-((px(i)**2+py(i)**2)*.5d0+dv(i))*al
            enddo
          endif
          go to 1011
        endif
 1200   continue
        al=rlist(itp)
        if(rlist(lp+14) .ne. 0.d0)then
          ak0=rlist(lp+2)+rlist(lp+11)+rlist(lp+14)*tgauss()
        else
          ak0=rlist(lp+2)+rlist(lp+11)
        endif
        if(rlist(lp+ilist(1,lp)) .gt. 0.d0)then
          fb1=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB1,icBEND))
          fb2=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB2,icBEND))
        else
          fb2=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB1,icBEND))
          fb1=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB2,icBEND))
        endif
        ak1=rlist(lp+kytbl(kwK1,icBEND))
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          ak0=ak0*rtaper
          ak1=ak1*rtaper
        endif
        call tbend(np,x,px,y,py,z,g,dv,pz,l,al,ak0,
     $       rlist(lp+kytbl(kwANGL,icBEND)),
     1       rlist(itp+1),rlist(itp+2),rlist(itp+3),rlist(itp+4),
     1       ak1,
     1       rlist(lp+kytbl(kwDX,icBEND)),
     $       rlist(lp+kytbl(kwDY,icBEND)),
     1       rlist(itp+13),rlist(itp+11),rlist(itp+12),
     1       rlist(itp+5),rlist(itp+6),
     $       fb1,fb2,
     $       nint(rlist(lp+kytbl(kwFRMD,icBEND))),
     $       rlist(lp+kytbl(kwFRIN,icBEND)) .eq. 0.d0,
     1       rlist(itp+7),rlist(itp+8),rlist(itp+9),rlist(itp+10),
     1       rlist(lp+kytbl(kwRAD,icBEND)) .eq. 0.d0,
     $       0.d0,al,al,
     1       rlist(lp+13))
        go to 1020
 1400   continue
        al=rlist(lp+1)
        if(rlist(lp+ilist(1,lp)) .gt. 0.d0)then
          mfr=nint(rlist(lp+12))
        else
          mfr=nint(rlist(lp+12))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tquad(np,x,px,y,py,z,g,dv,pz,l,al,
     1       rlist(lp+kytbl(kwK1,icQUAD))*rtaper,
     $       rlist(lp+5),rlist(lp+6),
     1       rlist(itp+7),rlist(itp+2),rlist(itp+3),
     1       rlist(lp+kytbl(kwRAD,icQUAD)),rlist(lp+8) .eq. 0.d0,
     1       rlist(lp+9) .eq. 0.d0,
     $       rlist(itp+5)*rtaper,rlist(itp+6)*rtaper,
     1       mfr,rlist(lp+kytbl(kwF1,icQuad)),rlist(lp+13),
     $       rlist(lp+14) .eq. 0.d0)
        go to 1020
 1600   ak1=rlist(lp+2)
        if(rad .and. radcod .and. radtaper)then
          ak1=ak1*(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tthin(np,x,px,y,py,z,g,dv,pz,lele,l,al,ak1,
     1       rlist(lp+5),rlist(lp+6),rlist(itp+4),
     1       rlist(itp+2),rlist(itp+3),rlist(lp+7),
     1       rlist(lp+8) .eq. 0.d0)
        go to 1020
 1800   call undulator(np,x,px,y,py,z,g,dv,pz,rlist(itp))
        go to 1020
 1900   call twig(np,x,px,y,py,z,g,dv,al,rlist(lp+2),
     1       int(rlist(lp+3)),rlist(lp+5),rlist(lp+6),
     1       rlist(lp+4),itp)
        go to 1020
 2000   call tsol(np,x,px,y,py,z,g,dv,pz,latt,l,lend,
     $       ke,sol,kptbl,la,n,
     $       nwakep,iwakeelm,kwaketbl,nwak,nextwake,out)
        if(np .le. 0)then
          return
        endif
        go to 1020
 2100   write(*,*)'Use BEND with ANGLE=0 for STEER.'
        call forcesf()
 2200   phi=rlist(lp+kytbl(kwANGL,icMULT))
        mfr=nint(rlist(lp+14))
        if(rlist(lp+ilist(1,lp)) .gt. 0.d0)then
          psi1=phi*rlist(lp+kytbl(kwE1,icMULT))
     $         +rlist(lp+kytbl(kwAE1,icMULT))
          psi2=phi*rlist(lp+kytbl(kwE2,icMULT))
     $         +rlist(lp+kytbl(kwAE2,icMULT))
          fb1=rlist(lp+kytbl(kwFB1,icMULT))
          fb2=rlist(lp+kytbl(kwFB2,icMULT))
          chi1=rlist(lp+kytbl(kwCHI1,icMULT))
          chi2=rlist(lp+kytbl(kwCHI2,icMULT))
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=phi*rlist(lp+kytbl(kwE2,icMULT))
     $         +rlist(lp+kytbl(kwAE2,icMULT))
          psi2=phi*rlist(lp+kytbl(kwE1,icMULT))
     $         +rlist(lp+kytbl(kwAE1,icMULT))
          fb2=rlist(lp+kytbl(kwFB1,icMULT))
          fb1=rlist(lp+kytbl(kwFB2,icMULT))
          chi1=-rlist(lp+kytbl(kwCHI1,icMULT))
          chi2=-rlist(lp+kytbl(kwCHI2,icMULT))
        endif
        bz=0.d0
        harm=rlist(lp+kytbl(kwHARM,icMULT))
        if(harm .eq. 0.d0)then
          w=pi2*rlist(lp+kytbl(kwFREQ,icMULT))/c
        else
          w=omega0*harm/c
        endif
        autophi=rlist(lp+kytbl(kwAPHI,icMULT)) .ne. 0.d0
        ph=rlist(lp+kytbl(kwDPHI,icMULT))
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tmulti(np,x,px,y,py,z,g,dv,pz,al,
     $       rlist(lp+kytbl(kwK0,icMULT)),
     $       bz,
     $       phi,psi1,psi2,
     1       rlist(lp+3),rlist(lp+4),rlist(lp+5),
     $       chi1,chi2,rlist(lp+8),
     $       rlist(lp+kytbl(kwDROT,icMULT)),
     $       rlist(lp+9),rlist(lp+10) .eq. 0.d0,
     $       rlist(lp+11) .eq. 0.d0,
     $       rlist(itp+1)*rtaper,rlist(itp+2)*rtaper,mfr,fb1,fb2,
     $       rlist(lp+15),w,rlist(lp+17),ph,
     $       rlist(lp+kytbl(kwRADI,icMULT)),rtaper,autophi,
     $       n,l,latt,kptbl)
        go to 1020
 3100   harm=rlist(lp+kytbl(kwHARM,icCAVI))
        if(harm .eq. 0.d0)then
          w=pi2*rlist(lp+kytbl(kwFREQ,icCAVI))/c
        else
          w=omega0*harm/c
        endif
        if(rlist(lp+9) .eq. 0.d0)then
          ak=rlist(lp+2)
        else
          ak=rlist(lp+2)+rlist(lp+9)*tgauss()
        endif
        if(rlist(lp+10) .eq. 0.d0)then
          ph=rlist(lp+kytbl(kwDPHI,icCAVI))
        else
          ph=rlist(lp+kytbl(kwDPHI,icCAVI))+
     $         rlist(lp+kytbl(kwRANP,icCAVI))*tgauss()
        endif
        autophi=rlist(lp+kytbl(kwAPHI,icCAVI)) .ne. 0.d0
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
c          write(*,*)'tturn ',l,gettwiss(mfitdz,l),ph
        endif
        mfr=nint(rlist(lp+kytbl(kwFRMD,icCAVI)))
        if(rlist(lp+ilist(1,lp)) .gt. 0.d0)then
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
            iwplc=abs(ilist(1,lp+kytbl(kwLWAK,icCAVI)))
            if(iwplc .eq. 0 .or. .not. lwake)then
              lwlc=0
            else
              lwlc=(ilist(1,iwplc-1)-2)/2
            endif
            iwptc=abs(ilist(1,lp+kytbl(kwTWAK,icCAVI)))
            if(iwptc .eq. 0 .or. .not. twake)then
              lwtc=0
            else
              lwtc=(ilist(1,iwptc-1)-2)/2
            endif
          endif
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,rlist(lp+kytbl(kwPHI,icCAVI)),ph,
     $         lwlc,rlist(iwplc+1),lwtc,rlist(iwptc+1),
     1         rlist(lp+13),rlist(lp+14),rlist(lp+15),
     $         rlist(lp+16),rlist(lp+17),rlist(lp+18),rlist(lp+19),
     $         rlist(lp+kytbl(kwFRIN,icCAVI)) .eq. 0.d0,mfr,autophi)
        else
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,rlist(lp+kytbl(kwPHI,icCAVI)),ph,
     $         0,0.d0,0,0.d0,
     1         rlist(lp+13),rlist(lp+14),rlist(lp+15),
     $         rlist(lp+16),rlist(lp+17),rlist(lp+18),rlist(lp+19),
     $         rlist(lp+kytbl(kwFRIN,icCAVI)) .eq. 0.d0,mfr,autophi)
        endif
        go to 1020
 3200   if(rfsw)then
          if(rlist(lp+9) .eq. 0.d0)then
            ak=rlist(lp+2)
          else
            ak=rlist(lp+2)+rlist(lp+9)*tgauss()
          endif
          if(rlist(lp+10) .eq. 0.d0)then
            ph=rlist(lp+4)
          else
            ph=rlist(lp+4)+rlist(lp+10)*tgauss()
          endif
          harmf=rlist(lp+3)-int(rlist(lp+3))
          ph=ph+harmf*(n-1)*pi2
          call ttcav(np,x,px,y,py,z,g,dv,pz,al,ak,rlist(lp+3),
     1         ph,rlist(lp+5),
     1         rlist(lp+6),rlist(lp+7),rlist(lp+8))
        else
          call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        endif
        go to 1020
 3300   continue
        call temap(np,np0,x,px,y,py,z,g,dv,l,n,kptbl)
        go to 1010
 3400   call tins(np,x,px,y,py,z,g,rlist(lp+20))
        go to 1010
 3500   call tcoord(np,x,px,y,py,z,
     1       rlist(lp+1),rlist(lp+2),rlist(lp+3),
     1       rlist(lp+4),rlist(lp+5),rlist(lp+6),
     1       rlist(lp+7) .eq. 0.d0)
        go to 1010
 3600   call beambeam(np,x,px,y,py,z,g,dv,pz,rlist(lp+1),rlist(itp),n)
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
     $           rlist(lp+kytbl(kwANGL,icMULT)) .ne. 0.d0)then
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
