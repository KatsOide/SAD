      subroutine tturne(trans,cod,beam,
     $     iatr,iacod,iabmi,plot,update,rt)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,ffs_bound
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      use temw, only:normali
      implicit none
      type (ffs_bound) fbound
      real*8 codmax,demax
      parameter (codmax=1.d4,demax=.5d0)
      integer*8 iatr,iacod,iabmi
      real*8 trans(6,12),cod(6),beam(42)
      real*8 z0,pgev00,alambdarf,dzmax,phis
      logical*4 plot,update,rt
      pgev00=pgev
      vc0=0.d0
      u0=0.d0
      hvc0=0.d0
      vcacc=0.d0
      dvcacc=0.d0
      ddvcacc=0.d0
      z0=cod(5)
      if(calint)then
        touckl(:) = 0.d0
        touckm(:,:,:) = 0.d0
        if(.not. allocated(toucke))then
          allocate(toucke(ntouckl,nlat))
        endif
        if(size(toucke, 2) .lt. nlat)then
          deallocate(toucke)
          allocate(toucke(ntouckl,nlat))
        endif
        toucke(:,:) = 0.d0
      endif
      if(irad .eq. 6)then
        npelm=0
      endif
      ipelm=0
      normali=.true.
      call tffsbound(fbound)
      call tturne0(trans,cod,beam,fbound,
     $     iatr,iacod,iabmi,0,plot,rt)
      if(update)then
        if(vcacc .ne. 0.d0)then
          wrfeff=sqrt(abs(ddvcacc/vcacc))
        elseif(vc0 .ne. 0.d0)then
          wrfeff=abs(dvcacc/vc0)
        else
          wrfeff=0.d0
        endif
        if(wrfeff .eq. 0.d0 .and. vc0 .ne. 0.d0)then
          wrfeff=hvc0/vc0*omega0/c
        endif
        if(wrfeff .ne. 0.d0)then
          alambdarf=pi2/wrfeff
          vceff=abs(dcmplx(vcacc,dvcacc/wrfeff))
        else
          alambdarf=circ
          vceff=0.d0
        endif
        if(vceff .eq. 0.d0)then
          vceff=vc0
        endif
        if(vcacc .eq. 0.d0)then
          vcacc=vceff*sin(trf0*wrfeff)
        endif
        if(trpt)then
          trf0=0.d0
          vcalpha=1.d0
        else
          if(vc0 .ne. 0.d0)then
            vcalpha=vceff/vc0
          else
            vcalpha=0.d0
          endif
          if(vceff .ne. 0.d0)then
            dzmax=alambdarf*.24d0
            phis=asin(abs(vcacc/vceff))
            if(vceff .gt. u0*pgev)then
              if(trans(5,6) .lt. 0.d0)then
                trf0=trf0+(asin(u0*pgev/vceff)-phis)/wrfeff
              else
                trf0=trf0+
     $                 (pi-asin(u0*pgev/vceff)-phis)/wrfeff
              endif
            else
              trf0=trf0+(.5*pi-phis)/wrfeff
            endif
          endif
          if(trf0 .lt. 0.d0)then
            trf0=-mod(-trf0+0.5d0*alambdarf,alambdarf)+alambdarf*0.5d0
          else
            trf0= mod(trf0-0.5d0*alambdarf,alambdarf)+alambdarf*0.5d0
          endif
        endif
        call RsetGL1('DTSYNCH',trf0)
        call RsetGL1('PHICAV',trf0*wrfeff)
        call RsetGL1('EFFVCRATIO',vcalpha)
        call RsetGL1('EFFVC',vceff)
        call RsetGL1('EFFRFFREQ',wrfeff*c/pi2)
      endif
      if(pgev00 .ne. pgev)then
        pgev=pgev00
        call tphyzp
      endif
      return
      end

      subroutine tturne0(trans,cod,beam,fbound,
     $     iatr,iacod,iabmi,idp,plot,rt)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,ffs_bound
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      implicit none
      type (ffs_bound) fbound
      type (sad_list), pointer :: kli
      integer*8 iatr,iacod,iabmi
      integer*4 ls,l,nvar,lx,idp,le1
      real*8 trans(6,12),cod(6),beam(42)
      real*8 trans1(6,12),cod1(6),beam1(42)
      real*8 vsave(kwMAX)
      real*8 r,xp,xb,xe,fr,fra,frb,tffselmoffset
      logical*4 sol,plot,chg,sol1,cp0,int0,rt
      sol=.false.
      if(fbound%fb .ne. 0.d0)then
        call qfracsave(fbound%lb,vsave,nvar,.true.)
        call qfraccomp(fbound%lb,fbound%fb,1.d0,ideal,chg)
        call tturne1(trans,cod,beam,
     $       iatr,iacod,iabmi,idp,plot,sol,rt,fbound%lb,fbound%lb)
        if(chg)then
          call qfracsave(fbound%lb,vsave,nvar,.false.)
        endif
        ls=fbound%lb+1
      else
        ls=fbound%lb
      endif
      if(fbound%fe .eq. 0.d0)then
        le1=min(nlat-1,fbound%le)
        call tturne1(trans,cod,beam,
     $       iatr,iacod,iabmi,idp,plot,sol,rt,ls,le1)
        if(plot)then
          call tfsetplot(trans,cod,beam,0,
     $         le1+1,iatr,iacod,.false.,idp)
        endif
      else
        call tturne1(trans,cod,beam,
     $       iatr,iacod,iabmi,idp,plot,sol,rt,ls,fbound%le-1)
        call qfracsave(fbound%le,vsave,nvar,.true.)
        call qfraccomp(fbound%le,0.d0,fbound%fe,ideal,chg)
        call tturne1(trans,cod,beam,
     $       iatr,iacod,iabmi,idp,plot,sol,rt,fbound%le,nlat-1)
        if(chg)then
          call qfracsave(fbound%le,vsave,nvar,.false.)
        endif
        if(plot)then
          call tfsetplot(trans,cod,beam,0,
     $         nlat,iatr,iacod,.false.,idp)
        endif
      endif
      if(plot)then
        if(codplt)then
          call tfadjustn(idp,mfitnx)
          call tfadjustn(idp,mfitny)
        endif
        xb=fbound%lb+fbound%fb
        xe=fbound%le+fbound%fe
        do l=1,nlat-1
          xp=min(xe,max(xb,tffselmoffset(l)))
          if(xp .ne. dble(l))then
            lx=int(xp)
            fr=xp-lx
 8101       if(fr .eq. 0.d0)then
              if(iatr .ne. 0)then
                if(iatr .gt. 0)then
                  if(l .ge. fbound%lb .and. l .le. fbound%le)then
                    call tflocal(klist(iatr+l))
                  endif
                  klist(iatr+l)=ktfcopy1(klist(iatr+lx))
                endif
                if(iacod .gt. 0)then
                  if(l .ge. fbound%lb .and. l .le. fbound%le)then
                    call tflocal(klist(iacod+l))
                  endif
                  klist(iacod+l)=ktfcopy(klist(iacod+lx))
                endif
              endif
              if(codplt)then
                twiss(l,idp,mfitdx:mfitddp)=twiss(lx,idp,mfitdx:mfitddp)
                if(irad .ge. 12)then
                  beamsize(:,l)=beamsize(:,lx)
                endif
              elseif(rt)then
                twiss(l,idp,mfitddp)=twiss(lx,idp,mfitddp)
              endif
            else
              if(lx .eq. fbound%lb)then
                fra=fbound%fb
                frb=max(fr,fra)
              elseif(lx .eq. fbound%le)then
                fra=0.d0
                frb=min(fbound%fe,fr)
              else
                fra=0.d0
                frb=fr
              endif
              if(fra .eq. frb)then
                fr=0.d0
                go to 8101
              endif
c     below is incorrect for fra <> 0
              call qfracsave(lx,vsave,nvar,.true.)
              call qfraccomp(lx,fra,frb,ideal,chg)
              if(.not. chg)then
                fr=0.d0
                go to 8101
              endif
              cod1=0.d0
              if(iatr .ne. 0)then
                if(iatr .gt. 0 .and. ktflistqd(dlist(iatr+lx),kli))then
                  call tfl2m(kli,trans1,6,6,.false.)
                else
                  call tinitr(trans1)
                endif
                if(iacod .gt. 0 .and.
     $               ktflistqd(dlist(iacod+lx),kli))then
                  cod1=kli%rbody(1:6)
                endif
              else
                call tinitr(trans1)
              endif
              trans1(:,7:12)=0.d0
              if(codplt)then
                cod1=twiss(lx,idp,mfitdx:mfitddp)
                if(irad .ge. 12)then
                  beam1(1:21)=beamsize(:,lx)
                  if(calint)then
                    beam1(22:42)=beamsize(:,lx)
                  endif
                endif
              endif
              sol1=.false.
              cp0=codplt
              codplt=.false.
              int0=calint
              calint=.false.
              call tturne1(trans1,cod1,beam1,
     $             int8(0),int8(0),int8(0),idp,.false.,sol1,rt,
     $             lx,lx)
              codplt=cp0
              calint=int0
              call qfracsave(lx,vsave,nvar,.false.)
c              write(*,*)'tturne0 ',lx,l,fbound%lb,fbound%le,idp,
c     $             gammab(lx)/(gammab(lx)*(1.d0-frb)+gammab(lx+1)*frb)
              call tfsetplot(trans1,cod1,beam1,lx,
     $             l,iatr,iacod,
     $             l .ge. fbound%lb .and. l .le. fbound%le,idp)
            endif
          endif
        enddo
      elseif(radtaper .and. radcod)then
        if(fbound%le .eq. 1)then
          r=1.d0
        else
          r=gammab(fbound%le-1)/gammab(fbound%le)
        endif
        twiss(fbound%le,idp,mfitddp)=cod(6)*r
      endif
      return
      end

      subroutine tturne1(trans,cod,beam,
     $     iatr,iacod,iabmi,idp,plot,sol,rt,ibegin,iend)
      use kyparam
      use tfstk
      use tffitcode
      use ffs, only: gettwiss
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      implicit none
      real*8 codmax,demax
      parameter (codmax=1.d4,demax=.5d0)
      type (sad_comp), pointer :: cmp
      integer*8 iatr,iacod,iabmi,kbmz,kbmzi,lp
      integer*4 idp,i,l1
      real*8 trans(6,12),cod(6),beam(42),bmir(6,6),
     $     bmi(21),bmh(21)
      real*8 psi1,psi2,apsi1,apsi2,alid,
     $     r,dir,al,alib,dtheta,theta0,ftable(4),
     $     phi,fb1,fb2,chi1,chi2,ak0,ak1,rtaper
      integer*4 l,ld,lele,kl,mfr,ibegin,iend,ke
      logical*4 sol,plot,bmaccum,plotib,isnan,rt,next
      save kbmz
      data kbmz /0/
      if(kbmz .eq. 0)then
        kbmz=ktadaloc(0,6)
        kbmzi=ktraaloc(-1,6)
        do i=1,6
          klist(kbmz+i)=ktflist+ktfcopy1(kbmzi)
        enddo
      endif
      bmaccum=.false.
      plotib=plot .and. iabmi .ne. 0
      alid=0.d0
      call tsetdvfs
      call tesetdv(cod(6))
      bradprev=0.d0
      do l=ibegin,iend
        next=inext(l) .ne. 0
        if(isnan(cod(1)) .or. isnan(cod(3)))then
          if(isnan(cod(1)))then
            cod(1)=0.d0
          endif
          if(isnan(cod(3)))then
            cod(3)=0.d0
          endif
          call tinitr(trans)
          if(.not. plot)then
            return
          endif
        endif
        if(sol)then
          sol=l .lt. ke
          alid=0.d0
          bmaccum=.false.
          go to 1010
        endif
        if(plot)then
          if(iatr .ne. 0)then
            if(iatr .gt. 0)then
              dlist(iatr+l)=
     $             dtfcopy1(kxm2l(trans,6,6,6,.false.))
            endif
            if(iacod .gt. 0)then
              dlist(iacod+l)=
     $             dtfcopy1(kxm2l(cod,0,6,1,.false.))
            endif
          endif
          if(codplt)then
c            if(l .eq. 1)then
c              r=1.d0
c            else
c              r=gammab(l-1)/gammab(l)
c            endif
            call tsetetwiss(trans,cod,beam,0,l,idp)
c            write(*,'(a,i5,1p6g15.7)')'tturne1 ',l,twiss(l,idp,1:6)
c            et=twiss(l,0,1:mfitzpy)
c            call checketwiss(trans,et)
          endif
        elseif(radtaper .and. radcod)then
          if(l .eq. 1)then
            r=1.d0
          else
            r=gammab(l-1)/gammab(l)
          endif
          twiss(l,idp,mfitddp)=cod(6)*r
        endif
        ld=idelc(l)
        lele=idtype(ld)
        if(ideal)then
          lp=idval(ld)
        else
          lp=elatt%comp(l)
        endif
        call loc_comp(lp,cmp)
        dir=direlc(l)
        kl=kytbl(kwL,lele)
        if(kl .ne. 0)then
          al=cmp%value(kl)
        else
          al=0.d0
        endif
        if(calint)then
          if(al .ne. 0.d0)then
            if(lele .eq. icDRFT)then
              alib=al*.25d0+alid
              alid=al*.25d0
            else
              alib=al*.5d0+alid
              alid=al*.5d0
            endif
            call tintrb(trans,cod,beam,bmi,alib,alid,l)
            if(plotib)then
              if(lele .eq. icDRFT)then
                if(bmaccum)then
                  bmh=bmh+bmi
                else
                  bmh=bmi
                endif
                bmaccum=.true.
              else
                if(bmaccum)then
                  bmi=bmi+bmh
                endif
                call tconvbm(bmi,bmir)
                dlist(iabmi+l)=
     $               dtfcopy1(kxm2l(bmir,6,6,6,.false.))
                bmaccum=.false.
              endif
            endif
          elseif(plotib)then
            klist(iabmi+l)=ktflist+ktfcopy1(kbmz)
          endif
        endif
c        WRITE(*,*)lele,' ',PNAME(ILIST(2,LATT(L)))(1:16)
c        if(l .lt. 5)then
c          write(*,*)'tturne1-l ',l,beam(6)
c        endif
        go to (1100,1200,1010,1400,1010,1600,1010,1600,1010,1600,
     $       1010,1600,1010,1010,1010,1010,1010,1010,1010,3000,
     $       3100,3200,1010,1010,1010,1010,1010,1010,1010,1010,
     $       4100,4200,4300,4400,4500),lele
        go to 5000
1100    continue
        if(calint)then
          if(al .ne. 0.d0)then
            al=al*.5d0
            call tdrife(trans,cod,beam,al,0.d0,0.d0,0.d0,
     $           .true.,.false.,calpol,irad,ld)
            call tintrb(trans,cod,beam,bmi,al,al*.5d0,l)
            if(plotib)then
              bmi=bmi*0.5d0
              bmh=bmh+bmi
              call tconvbm(bmh,bmir)
              bmh=bmi
              dlist(iabmi+l)=
     $             dtfcopy1(kxm2l(bmir,6,6,6,.false.))
            endif
          endif
        endif
        call tdrife(trans,cod,beam,al,0.d0,0.d0,0.d0,
     $       .true.,.false.,calpol,irad,ld)
        go to 1010
1200    continue
        if(dir .gt. 0.d0)then
          psi1=cmp%value(ky_E1_BEND)
          psi2=cmp%value(ky_E2_BEND)
          apsi1=cmp%value(ky_AE1_BEND)
          apsi2=cmp%value(ky_AE2_BEND)
          fb1=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB1_BEND)
          fb2=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB2_BEND)
        else
          psi1=cmp%value(ky_E2_BEND)
          psi2=cmp%value(ky_E1_BEND)
          apsi1=cmp%value(ky_AE2_BEND)
          apsi2=cmp%value(ky_AE1_BEND)
          fb2=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB1_BEND)
          fb1=cmp%value(ky_F1_BEND)
     $         +cmp%value(ky_FB2_BEND)
        endif
        dtheta=cmp%value(ky_DROT_BEND)
        theta0=cmp%value(ky_ROT_BEND)+dtheta
        ak0=cmp%value(ky_K0_BEND)
     $       +cmp%value(ky_ANGL_BEND)
        ak1=cmp%value(ky_K1_BEND)
        if(radcod .and. radtaper)then
          if(rt)then
            l1=nextl(l)
            ak0=ak0*(4.d0+3.d0*cod(6)+gettwiss(mfitddp,l1))*.25d0
            ak1=ak1*(4.d0+3.d0*cod(6)+gettwiss(mfitddp,l1))*.25d0
          else
            ak0=ak0*(1.d0+cod(6))
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        call tbende(trans,cod,beam,al,
     1       min(pi2,max(-pi2,ak0)),
     $       cmp%value(ky_ANGL_BEND),
     $       psi1,psi2,apsi1,apsi2,ak1,
     1       cmp%value(ky_DX_BEND),
     $       cmp%value(ky_DY_BEND),theta0,dtheta,
     $       fb1,fb2,
     $       nint(cmp%value(ky_FRMD_BEND)),
     $       cmp%value(ky_FRIN_BEND) .eq. 0.d0,
     $       cmp%value(ky_EPS_BEND),
     1       cmp%value(ky_RAD_BEND) .eq. 0.d0,.true.,
     $       next,l,ld)
        go to 1010
1400    continue
        if(dir .gt. 0.d0)then
          mfr=nint(cmp%value(ky_FRMD_QUAD))
        else
          mfr=nint(cmp%value(ky_FRMD_QUAD))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        ak1=cmp%value(ky_K1_QUAD)
        if(radcod .and. radtaper)then
          if(rt)then
            l1=nextl(l)
            ak1=ak1*(2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0
          else
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        call tsetfringepe(cmp,icQUAD,dir,ftable)
        call tquade(trans,cod,beam,al,ak1,
     $       cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1       cmp%value(ky_ROT_QUAD),
     $       cmp%value(ky_RAD_QUAD) .eq. 0.d0,
     1       cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,cmp%value(ky_EPS_QUAD),
     $       cmp%value(ky_KIN_QUAD) .eq. 0.d0,next,ld)
        go to 1010
 1600   ak1=cmp%value(ky_K_THIN)
        if(radcod .and. radtaper)then
          if(rt)then
            l1=nextl(l)
            ak1=ak1*(2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0
          else
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        call tthine(trans,cod,beam,lele,al,ak1,
     1       cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     $       cmp%value(ky_ROT_THIN),.false.,ld)
        go to 1010
 3000   call tsole(trans,cod,beam,l,ke,sol,
     1       iatr,iacod,iabmi,idp,plot,rt)
        alid=0.d0
        go to 1010
 3100   write(*,*)'Use BEND with ANGLE=0 for ST.'
        call forcesf()
 3200   phi=cmp%value(ky_ANGL_MULT)
        mfr=nint(cmp%value(ky_FRMD_MULT))
        if(dir .gt. 0.d0)then
          psi1=cmp%value(ky_E1_MULT)
          psi2=cmp%value(ky_E2_MULT)
          apsi1=cmp%value(ky_AE1_MULT)
          apsi2=cmp%value(ky_AE2_MULT)
          fb1=cmp%value(ky_FB1_MULT)
          fb2=cmp%value(ky_FB2_MULT)
          chi1=cmp%value(ky_CHI1_MULT)
          chi2=cmp%value(ky_CHI2_MULT)
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=cmp%value(ky_E2_MULT)
          psi2=cmp%value(ky_E1_MULT)
          apsi1=cmp%value(ky_AE2_MULT)
          apsi2=cmp%value(ky_AE1_MULT)
          fb2=cmp%value(ky_FB1_MULT)
          fb1=cmp%value(ky_FB2_MULT)
          chi1=-cmp%value(ky_CHI1_MULT)
          chi2=-cmp%value(ky_CHI2_MULT)
        endif
        rtaper=1.d0
        if(radcod .and. radtaper)then
          if(rt)then
            l1=nextl(l)
            rtaper=(2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0
          else
            rtaper=1.d0+cod(6)
          endif
        endif
        call tsetfringepe(cmp,icMULT,dir,ftable)
        call tmulte(trans,cod,beam,l,al,
     $       cmp%value(ky_K0_MULT),
     $       0.d0,
     $       phi,psi1,psi2,apsi1,apsi2,
     1       cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $       cmp%value(ky_DZ_MULT),
     $       chi1,chi2,cmp%value(ky_ROT_MULT),
     $       cmp%value(ky_DROT_MULT),
     $       cmp%value(ky_EPS_MULT),
     $       cmp%value(ky_RAD_MULT) .eq. 0.d0,
     $       cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,fb1,fb2,
     $       cmp%value(ky_K0FR_MULT) .eq. 0.d0,
     $       cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT),
     $       cmp%value(ky_HARM_MULT),
     $       cmp%value(ky_PHI_MULT),cmp%value(ky_FREQ_MULT),
     $       cmp%value(ky_W1_MULT),rtaper,
     $       cmp%value(ky_APHI_MULT) .ne. 0.d0,
     $       ld)
        go to 1010
 4100   mfr=nint(cmp%value(ky_FRMD_CAVI))
        if(dir .gt. 0.d0)then
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
c        write(*,*)'tturne-tcave',cod
        call tcave(trans,cod,beam,l,al,
     1       cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI),
     $       cmp%value(ky_HARM_CAVI),
     1       cmp%value(ky_PHI_CAVI),cmp%value(ky_FREQ_CAVI),
     1       cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $       cmp%value(ky_ROT_CAVI),
     $       cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $       cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $       cmp%value(ky_FRIN_CAVI) .eq. 0.d0,mfr,
     $       cmp%value(ky_APHI_CAVI) .ne. 0.d0,
     $       ld)
c        write(*,*)'tturne-tcave-1',cod
        go to 1010
 4200   call ttcave(trans,cod,beam,al,
     1       cmp%value(ky_K0_TCAV),cmp%value(ky_HARM_TCAV),
     1       cmp%value(ky_PHI_TCAV),cmp%value(ky_FREQ_TCAV),
     1       cmp%value(ky_DX_TCAV),cmp%value(ky_DY_TCAV),
     $       cmp%value(ky_ROT_TCAV),ld)
        go to 1010
 4300   call temape(trans,cod,beam,l)
        go to 1010
 4400   call tinse(trans,cod,beam,cmp%value(ky_DIR_INS+1),ld)
        go to 1010
 4500   call tcoorde(trans,cod,beam,
     1       cmp%value(ky_DX_COORD),cmp%value(ky_DY_COORD),
     $       cmp%value(ky_DZ_COORD),cmp%value(ky_CHI1_COORD),
     $       cmp%value(ky_CHI2_COORD),cmp%value(ky_CHI3_COORD),
     1       cmp%value(ky_DIR_COORD) .eq. 0.d0,ld)
        go to 1010
 5000   go to 1010
 1010   continue
      enddo
      if(calint)then
        if(alid .ne. 0.d0)then
          call tintrb(trans,cod,beam,bmi,alid,alid,l)
          alid=0.d0
        else
          bmi=0.d0
        endif
        if(plotib)then
          if(bmaccum)then
            call tadd(bmh,bmi,bmi,21)
          endif
          call tconvbm(bmi,bmir)
          dlist(iabmi+iend+1)=
     $         dtfcopy1(kxm2l(bmir,6,6,6,.false.))
        endif
      endif
      return
      end

      subroutine tconvbm(b,br)
      implicit none
      real*8 b(21),br(6,6)
      integer*4 i,j,n
      n=0
      do i=1,6
        do j=1,i
          n=n+1
          br(i,j)=b(n)
          br(j,i)=b(n)
        enddo
      enddo
      return
      end

      subroutine tfadjustn(idp,m)
      use tffitcode
      use ffs_pointer, only:twiss
      use tmacro
      implicit none
      integer*4 idp,m,l
      real*8 phi0
      real*8,parameter :: toln=2.d-9
      phi0=0.d0
      do l=2,nlat
        twiss(l,idp,m)=twiss(l,idp,m)+phi0
        if(twiss(l,idp,m)+toln .lt. twiss(l-1,idp,m))then
          phi0=phi0+pi2
          twiss(l,idp,m)=twiss(l,idp,m)+pi2
        endif
      enddo
      return
      end

      subroutine tfsetplot(trans,cod,beam,lorg,
     $     l,iatr,iacod,local,idp)
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tffitcode
      use tmacro
      implicit none
      integer*8 iatr,iacod
      integer*4 l,idp,lorg
      real*8 trans(6,12),cod(6),beam(21)
      logical*4 local
      if(iatr .ne. 0)then
        if(iatr .gt. 0)then
          if(local)then
            call tflocal(klist(iatr+l))
          endif
          dlist(iatr+l)=
     $         dtfcopy1(kxm2l(trans,6,6,6,.false.))
        endif
        if(iacod .gt. 0)then
          if(local)then
            call tflocal(klist(iacod+l))
          endif
          dlist(iacod+l)=
     $         dtfcopy1(kxm2l(cod,0,6,1,.false.))
        endif
      endif
      if(codplt)then
        call tsetetwiss(trans,cod,beam,lorg,l,idp)
c        write(*,'(a,2i5,1p6g15.7)')'tsetplot  ',lorg,l,
c     $       twiss(l,idp,mfitzx:mfitzpy)
      elseif(radcod .and. radtaper)then
        twiss(l,idp,mfitddp)=cod(6)
      endif
      return
      end

      subroutine tsetetwiss(trans,cod,beam,lorg,l,idp)
      use ffs
      use ffs_pointer
      use tffitcode
      use tmacro
      use temw
      implicit none
      integer*4 l,idp,lorg,l0
      real*8 trans(6,6),ti(6,6),twi(ntwissfun),cod(6),
     $     beam(21),ril(6,6),gr,tr0(6,6)
      logical*4 norm
      real*8,parameter :: toln=-2.d-9
      if(trpt)then
        gr=gammab(l)/gammab(max(1,lorg-1))
        tr0=trans*sqrt(gr)
        call tinv6(tr0,ti)
      else
        call tinv6(trans,ti)
      endif
      if(lorg .eq. 0)then
        call tmultr(ti,ri,6)
        norm=normali
        l0=1
      else
        l0=lorg
        twi=twiss(lorg,idp,1:ntwissfun)
        call etwiss2ri(twi,ril,norm)
        call tmultr(ti,ril,6)
      endif
      call tfetwiss(ti,cod,twi,norm)
      if(l .eq. 1)then
        twi(mfitnx)=0.d0
        twi(mfitny)=0.d0
        twi(mfitnz)=0.d0
      else
        if(twi(mfitnx) .lt. toln)then
          twi(mfitnx)=twiss(l0,idp,mfitnx)+twi(mfitnx)+pi2
        else
          twi(mfitnx)=twiss(l0,idp,mfitnx)+twi(mfitnx)
        endif
        if(twi(mfitny) .lt. toln)then
          twi(mfitny)=twiss(l0,idp,mfitny)+twi(mfitny)+pi2
        else
          twi(mfitny)=twiss(l0,idp,mfitny)+twi(mfitny)
        endif
c        write(*,*)'setetwiss ',l,l0,idp,
c     $       twi(mfitny),twiss(l0,idp,mfitny)
        twi(mfitnz)=twiss(l0,idp,mfitnz)+twi(mfitnz)
      endif
c      twi(mfitdpx)=twi(mfitdpx)*rgb
c      twi(mfitdpy)=twi(mfitdpy)*rgb
c      twi(mfitddp)=twi(mfitddp)*rgb
c      write(*,*)'setetwiss ',twi(mfitddp),rgb
      twiss(l,idp,1:ntwissfun)=twi
      if(irad .ge. 12)then
        beamsize(:,l)=beam
      endif
      return
      end

      subroutine checketwiss(trans,tw1)
      use ffs_pointer
      use tffitcode
      use tmacro
      use temw
      implicit none
      integer*4 i
      real*8 tw1(ntwissfun),ra(6,6),trans(6,6),ti(6,6)
      logical*4 normal
      call etwiss2ri(tw1,ra,normal)
      ti=r
      call tmultr(ti,trans,6)
      call tmultr(ti,ra,6)
      write(*,*)'checketwiss ',tw1(mfitdetr)
      do i=1,6
        write(*,'(1p6g15.7)')ti(i,:)
      enddo
      return
      end

      subroutine tesetdv(dp)
      use tfstk
      use tmacro
      implicit none
      real*8 dp,p1
      p1=p0*(1.d0+dp)
      h1emit=p2h(p1)
c      h1emit=p1*sqrt(1.d0+1.d0/p1**2)
c      h1emit=p1+1.d0/(sqrt(1.d0+p1**2)+p1)
      dvemit=-(p1+p0)/h1emit/(p0*h1emit+p1*h0)*dp+dvfs
      return
      end

      subroutine tsetdvfs
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      real*8 rgetgl1
      if(rfsw)then
        dvfs=rgetgl1('FSHIFT')
      else
        dvfs=0.d0
      endif
      return
      end

      subroutine tgetdv(dp,dv,dvdp)
      use tfstk
      use tmacro
      implicit none
      real*8 dp,dv,dvdp,pr,p1,h1
      pr=1.d0+dp
      p1=p0*pr
      h1=p2h(p1)
c      h1=p1*sqrt(1.d0+1.d0/p1**2)
c      h1=p1+1.d0/(sqrt(1.d0+p1**2)+p1)
      dv=-(pr+1.d0)/h1/(h1+pr*h0)*dp+dvfs
      dvdp=h0/h1**3
      return
      end

      subroutine tgetdvh(dh,dv)
      use tfstk
      use tmacro
      implicit none
      real*8 dh,dv,h1,p1
      if(dh .ne. 0.d0)then
        h1=h0+dh
        p1=h2p(h1)
c        p1=h1*sqrt(1.d0-1.d0/h1**2)
c        p1=h1-1.d0/(sqrt(h1**2-1.d0)+h1)
        dv=-(p0+p1)/(p1*h0+p0*h1)/p1/p0*dh
      else
        dv=0.d0
      endif
      return
      end
