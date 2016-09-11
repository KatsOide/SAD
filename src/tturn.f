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
      call tturn0(np,latt,1,nlat,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     .false.,0,0,0)
      return
      end

      subroutine tturn0(np,latt,lb,le,x,px,y,py,z,g,dv,pz,kptbl,n)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      real*8 plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      integer*4 np,n,la,ls,lbegin,lend,nvar,lb,le
c      integer*4 isb,itwb,itwb1,itwb2,itwb3,itwb4,ntw
      real*8 frbegin,frend,pgev00
      integer*4 kptbl(np0,6)
      integer*8 latt(nlat)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sa(6),ss(6,6),vsave(100)
      logical*4 sol,chg,tfinsol
      pgev00=pgev
      sol=tfinsol(lb)
      novfl=0
      la=0
      if(radlight)then
        call tlstore(np,x,y,z,dv,0.d0,0.d0,0.d0,0.d0,
     $       p0/h0*c,dvfs,.true.)
      endif
      call tffsbound1(lb,le,lbegin,frbegin,lend,frend)
      if(frbegin .ne. 0.d0)then
        call qfracsave(lbegin,vsave,nvar,.true.)
        call qfraccomp(lbegin,frbegin,1.d0,ideal,chg)
        call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $       sol,la,lbegin,lbegin)
        if(chg)then
          call qfracsave(lbegin,vsave,nvar,.false.)
        endif
        ls=lbegin+1
      else
        ls=lbegin
      endif
      call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $     sol,la,ls,lend-1)
      if(frend .ne. 0.d0)then
        call qfracsave(lend,vsave,nvar,.true.)
        call qfraccomp(lend,0.d0,frend,ideal,chg)
        call tturn1(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n,
     $       sol,la,lend,lend)
        if(chg)then
          call qfracsave(lend,vsave,nvar,.false.)
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
      use tfstk
      use ffs
      use tffitcode
      use ffs_wake
      use sad_main
      use ffs_pointer, only: direlc,compelc
      implicit none
      integer*4 la1
      parameter (la1=15)
      real*8 xlimit,plimit,zlimit,vmax
      parameter (plimit=0.7d0,zlimit=1.d10)
      parameter (vmax=.9d0)
      type (sad_comp), pointer:: cmp
      integer*4 np,n,la,lbegin,lend,kdx,kdy,krot
      integer*4 kptbl(np0,6)
      integer*8 latt(nlat)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sa(6),ss(6,6),phi,bz,harm,w,
     $     al,a,dpz,al1,ak0,ak1,psi1,psi2,tgauss,ph,harmf,
     $     sspac0,sspac,fw,dx,dy,rot,sspac1,sspac2,
     $     fb1,fb2,chi1,chi2,ak,rtaper
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
       al=cmp%value(kytbl(kwL,lele))
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
          call spdrift_free(np,x,px,y,py,z,g,dv,pz,cmp%value(1),
     $         cmp%value(kytbl(kwRADI,icDRFT)),n,l,latt,kptbl)
          go to 1020
        else
          if(cmp%value(kytbl(kwKIN,lele)) .eq. 0.d0)then
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
        if(cmp%value(14) .ne. 0.d0)then
          ak0=cmp%value(2)+cmp%value(11)+cmp%value(14)*tgauss()
        else
          ak0=cmp%value(2)+cmp%value(11)
        endif
        if(direlc(l) .gt. 0.d0)then
          fb1=cmp%value(kytbl(kwF1,icBEND))
     $         +cmp%value(kytbl(kwFB1,icBEND))
          fb2=cmp%value(kytbl(kwF1,icBEND))
     $         +cmp%value(kytbl(kwFB2,icBEND))
        else
          fb2=cmp%value(kytbl(kwF1,icBEND))
     $         +cmp%value(kytbl(kwFB1,icBEND))
          fb1=cmp%value(kytbl(kwF1,icBEND))
     $         +cmp%value(kytbl(kwFB2,icBEND))
        endif
        ak1=cmp%value(kytbl(kwK1,icBEND))
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          ak0=ak0*rtaper
          ak1=ak1*rtaper
        endif
        call tbend(np,x,px,y,py,z,g,dv,pz,l,al,ak0,
     $       cmp%value(kytbl(kwANGL,icBEND)),
     1       rlist(itp+1),rlist(itp+2),rlist(itp+3),rlist(itp+4),
     1       ak1,
     1       cmp%value(kytbl(kwDX,icBEND)),
     $       cmp%value(kytbl(kwDY,icBEND)),
     1       rlist(itp+13),rlist(itp+11),rlist(itp+12),
     1       rlist(itp+5),rlist(itp+6),
     $       fb1,fb2,
     $       nint(cmp%value(kytbl(kwFRMD,icBEND))),
     $       cmp%value(kytbl(kwFRIN,icBEND)) .eq. 0.d0,
     1       rlist(itp+7),rlist(itp+8),rlist(itp+9),rlist(itp+10),
     1       cmp%value(kytbl(kwRAD,icBEND)) .eq. 0.d0,
     $       0.d0,al,al,
     1       cmp%value(13))
        go to 1020
 1400   continue
        if(direlc(l) .gt. 0.d0)then
          mfr=nint(cmp%value(12))
        else
          mfr=nint(cmp%value(12))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tquad(np,x,px,y,py,z,g,dv,pz,l,al,
     1       cmp%value(kytbl(kwK1,icQUAD))*rtaper,
     $       cmp%value(5),cmp%value(6),
     1       rlist(itp+4),rlist(itp+2),rlist(itp+3),
     1       cmp%value(kytbl(kwRAD,icQUAD)),cmp%value(8) .eq. 0.d0,
     1       cmp%value(9) .eq. 0.d0,
     $       rlist(itp+6)*rtaper,rlist(itp+7)*rtaper,
     $       rlist(itp+8)*rtaper,rlist(itp+9)*rtaper,
     1       mfr,cmp%value(13),
     $       cmp%value(14) .eq. 0.d0)
        go to 1020
 1600   ak1=cmp%value(2)
        if(rad .and. radcod .and. radtaper)then
          ak1=ak1*(2.d0+gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tthin(np,x,px,y,py,z,g,dv,pz,lele,l,al,ak1,
     1       cmp%value(5),cmp%value(6),rlist(itp+4),
     1       rlist(itp+2),rlist(itp+3),cmp%value(7),
     1       cmp%value(8) .eq. 0.d0)
        go to 1020
 1800   call undulator(np,x,px,y,py,z,g,dv,pz,rlist(itp))
        go to 1020
 1900   call twig(np,x,px,y,py,z,g,dv,al,cmp%value(2),
     1       int(cmp%value(3)),cmp%value(5),cmp%value(6),
     1       cmp%value(4),itp)
        go to 1020
 2000   call tsol(np,x,px,y,py,z,g,dv,pz,latt,l,lend,
     $       ke,sol,kptbl,la,n,nwak,nextwake,out)
        if(np .le. 0)then
          return
        endif
        go to 1020
 2100   write(*,*)'Use BEND with ANGLE=0 for STEER.'
        call forcesf()
 2200   phi=cmp%value(kytbl(kwANGL,icMULT))
        mfr=nint(cmp%value(14))
        if(direlc(l) .gt. 0.d0)then
          psi1=phi*cmp%value(kytbl(kwE1,icMULT))
     $         +cmp%value(kytbl(kwAE1,icMULT))
          psi2=phi*cmp%value(kytbl(kwE2,icMULT))
     $         +cmp%value(kytbl(kwAE2,icMULT))
          fb1=cmp%value(kytbl(kwFB1,icMULT))
          fb2=cmp%value(kytbl(kwFB2,icMULT))
          chi1=cmp%value(kytbl(kwCHI1,icMULT))
          chi2=cmp%value(kytbl(kwCHI2,icMULT))
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=phi*cmp%value(kytbl(kwE2,icMULT))
     $         +cmp%value(kytbl(kwAE2,icMULT))
          psi2=phi*cmp%value(kytbl(kwE1,icMULT))
     $         +cmp%value(kytbl(kwAE1,icMULT))
          fb2=cmp%value(kytbl(kwFB1,icMULT))
          fb1=cmp%value(kytbl(kwFB2,icMULT))
          chi1=-cmp%value(kytbl(kwCHI1,icMULT))
          chi2=-cmp%value(kytbl(kwCHI2,icMULT))
        endif
        bz=0.d0
        harm=cmp%value(kytbl(kwHARM,icMULT))
        if(harm .eq. 0.d0)then
          w=pi2*cmp%value(kytbl(kwFREQ,icMULT))/c
        else
          w=omega0*harm/c
        endif
        autophi=cmp%value(kytbl(kwAPHI,icMULT)) .ne. 0.d0
        ph=cmp%value(kytbl(kwDPHI,icMULT))
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
        endif
        rtaper=1.d0
        if(rad .and. radcod .and. radtaper)then
          rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
        endif
        call tmulti(np,x,px,y,py,z,g,dv,pz,al,
     $       cmp%value(kytbl(kwK0,icMULT)),
     $       bz,
     $       phi,psi1,psi2,
     1       cmp%value(3),cmp%value(4),cmp%value(5),
     $       chi1,chi2,cmp%value(8),
     $       cmp%value(kytbl(kwDROT,icMULT)),
     $       cmp%value(9),cmp%value(10) .eq. 0.d0,
     $       cmp%value(11) .eq. 0.d0,
     $       rlist(itp+1)*rtaper,rlist(itp+2)*rtaper,
     $       rlist(itp+3)*rtaper,rlist(itp+4)*rtaper,
     $       mfr,fb1,fb2,
     $       cmp%value(15),w,cmp%value(17),ph,
     $       cmp%value(kytbl(kwRADI,icMULT)),rtaper,autophi,
     $       n,l,latt,kptbl)
        go to 1020
 3100   harm=cmp%value(kytbl(kwHARM,icCAVI))
        if(harm .eq. 0.d0)then
          w=pi2*cmp%value(kytbl(kwFREQ,icCAVI))/c
        else
          w=omega0*harm/c
        endif
        if(cmp%value(9) .eq. 0.d0)then
          ak=cmp%value(2)
        else
          ak=cmp%value(2)+cmp%value(9)*tgauss()
        endif
        if(cmp%value(10) .eq. 0.d0)then
          ph=cmp%value(kytbl(kwDPHI,icCAVI))
        else
          ph=cmp%value(kytbl(kwDPHI,icCAVI))+
     $         cmp%value(kytbl(kwRANP,icCAVI))*tgauss()
        endif
        autophi=cmp%value(kytbl(kwAPHI,icCAVI)) .ne. 0.d0
        if(autophi)then
          ph=ph+gettwiss(mfitdz,l)*w
c          write(*,*)'tturn ',l,gettwiss(mfitdz,l),ph
        endif
        mfr=nint(cmp%value(kytbl(kwFRMD,icCAVI)))
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
c            iwplc=abs(ilist(1,lp+kytbl(kwLWAK,icCAVI)))
c            if(iwplc .eq. 0 .or. .not. lwake)then
              lwlc=0
c            else
c              lwlc=(ilist(1,iwplc-1)-2)/2
c            endif
c            iwptc=abs(ilist(1,lp+kytbl(kwTWAK,icCAVI)))
c            if(iwptc .eq. 0 .or. .not. twake)then
              lwtc=0
c            else
c              lwtc=(ilist(1,iwptc-1)-2)/2
c            endif
          endif
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,cmp%value(kytbl(kwPHI,icCAVI)),ph,
     $         lwlc,rlist(iwplc+1),lwtc,rlist(iwptc+1),
     1         cmp%value(13),cmp%value(14),cmp%value(15),
     $         cmp%value(16),cmp%value(17),cmp%value(18),cmp%value(19),
     $         cmp%value(kytbl(kwFRIN,icCAVI)) .eq. 0.d0,mfr,autophi)
        else
          call tcav(np,x,px,y,py,z,g,dv,al,ak,
     1         w,cmp%value(kytbl(kwPHI,icCAVI)),ph,
     $         0,0.d0,0,0.d0,
     1         cmp%value(13),cmp%value(14),cmp%value(15),
     $         cmp%value(16),cmp%value(17),cmp%value(18),cmp%value(19),
     $         cmp%value(kytbl(kwFRIN,icCAVI)) .eq. 0.d0,mfr,autophi)
        endif
        go to 1020
 3200   if(rfsw)then
          if(cmp%value(9) .eq. 0.d0)then
            ak=cmp%value(2)
          else
            ak=cmp%value(2)+cmp%value(9)*tgauss()
          endif
          if(cmp%value(10) .eq. 0.d0)then
            ph=cmp%value(4)
          else
            ph=cmp%value(4)+cmp%value(10)*tgauss()
          endif
          harmf=cmp%value(3)-int(cmp%value(3))
          ph=ph+harmf*(n-1)*pi2
          call ttcav(np,x,px,y,py,z,g,dv,pz,al,ak,cmp%value(3),
     1         ph,cmp%value(5),
     1         cmp%value(6),cmp%value(7),cmp%value(8))
        else
          call tdrift_free(np,x,px,y,py,z,g,dv,pz,al)
        endif
        go to 1020
 3300   continue
c        write(*,*)'ttrun1-temaxp ',np,np0,i,n
        call temap(np,np0,x,px,y,py,z,g,dv,l,n,kptbl)
        go to 1010
 3400   call tins(np,x,px,y,py,z,g,cmp%value(20))
        go to 1010
 3500   call tcoord(np,x,px,y,py,z,
     1       cmp%value(1),cmp%value(2),cmp%value(3),
     1       cmp%value(4),cmp%value(5),cmp%value(6),
     1       cmp%value(7) .eq. 0.d0)
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
     $           cmp%value(kytbl(kwANGL,icMULT)) .ne. 0.d0)then
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
