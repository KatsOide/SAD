      subroutine tturne(trans,cod,beam,srot,iae,plot,update,rt,optics)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,ffs_bound
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      use temw, only:normali,calint,iaemit
      use tfcsi, only:icslfno
      implicit none
      type (iaemit),intent(in):: iae
      type (ffs_bound) fbound
      real*8 codmax,demax
      parameter (codmax=1.d4,demax=.5d0)
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 z0,pgev00,alambdarf,dzmax,phis,vcacc1
      logical*4 ,intent(in):: plot,update,rt,optics
      pgev00=pgev
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
      normali=.true.
      call tffsbound(fbound)
      call tturneg(trans,cod,beam,srot,fbound,iae,plot,rt,optics)
      if(update)then
        if(wrfeff .ne. 0.d0)then
          alambdarf=pi2/wrfeff
          vceff=abs(dcmplx(vcacc,dvcacc/wrfeff))
        else
          alambdarf=circ
          vceff=0.d0
        endif
        if(vceff .ne. 0.d0)then
          dzmax=alambdarf*.24d0
          phis=asin(abs(vcacc1/vceff))
          if(radcod)then
c     trf0=-(cod(5)+z0)*0.5d0
          else
            if(vceff .gt. u0*pgev)then
              if(trans(5,6) .lt. 0.d0)then
                trf0=(asin(u0*pgev/vceff))/wrfeff
              else
                trf0=(pi-asin(u0*pgev/vceff))/wrfeff
              endif
c         write(*,'(a,1p6g15.7)')'ttrune ',u0*pgev,vceff,wrfeff,
c     $             trans(5,6),trf0
c         write(*,'(1p6g15.7)')(trans(i,1:6),i=1,6),cod
            else
              trf0=(.5*pi)/wrfeff
            endif
            if(trf0 .lt. 0.d0)then
              trf0=-mod(-trf0+0.5d0*alambdarf,alambdarf)
     $             +alambdarf*0.5d0
            else
              trf0= mod(trf0-0.5d0*alambdarf,alambdarf)
     $             +alambdarf*0.5d0
            endif
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

      subroutine tturneg(trans,cod,beam,srot,fbound,
     $     iae,plot,rt,optics)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,ffs_bound
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      use temw, only:normali,iaemit
      use tfcsi, only:icslfno
      implicit none
      type (iaemit) ,intent(in):: iae
      type (ffs_bound) ,intent(in):: fbound
      real*8 ,parameter::codmax=1.d4,demax=.5d0
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 z0,vcacc1
      logical*4 ,intent(in):: plot,rt,optics
      vc0=0.d0
      z0=cod(5)
      u0=0.d0
      hvc0=0.d0
      vcacc=0.d0
      dvcacc=0.d0
      ddvcacc=0.d0
      normali=.true.
      call tturne0(trans,cod,beam,srot,fbound,iae,0,plot,rt,optics)
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
        vceff=abs(dcmplx(vcacc,dvcacc/wrfeff))
      else
        vceff=0.d0
      endif
      if(vceff .eq. 0.d0)then
        vceff=vc0
      endif
      vcacc1=vcacc
      if(vcacc1 .eq. 0.d0)then
        vcacc1=vceff*sin(trf0*wrfeff)
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
      endif
      return
      end

      subroutine tturne0(trans,cod,beam,srot,fbound,
     $     iae,idp,plot,rt,optics)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,ffs_bound
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      use temw, only:calint,iaemit,iaez,beamplt
      implicit none
      type (iaemit) ,intent(in):: iae
      type (ffs_bound) ,intent(in):: fbound
      type (sad_dlist), pointer :: kli
      type (sad_rlist), pointer :: klir
      type (sad_comp) , pointer :: cmp
      integer*4 ,intent(in):: idp
      integer*4 ls,l,nvar,lx,le1,lv,itfdownlevel,irtc
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 trans1(6,12),cod1(6),beam1(42)
      type (sad_descriptor) dsave(kwMAX)
      real*8 r,xp,xb,xe,fr,fra,frb,tffselmoffset
      logical*4 ,intent(in):: plot,rt,optics
      logical*4 sol,chg,sol1,cp0,int0
      sol=.false.
      levele=levele+1
c      write(*,*)'tturne1-1 ',fbound%lb,fbound%fb,fbound%le,fbound%fe
      if(fbound%fb .ne. 0.d0)then
        call compelc(fbound%lb,cmp)
        call qfracsave(fbound%lb,dsave,nvar,.true.)
        call qfracseg(cmp,cmp,fbound%fb,1.d0,chg,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(l,irtc)
        else
c        call qfraccomp(fbound%lb,fbound%fb,1.d0,ideal,chg)
          call tturne1(trans,cod,beam,srot,
     $         iae,idp,plot,sol,rt,optics,
     $         fbound%lb,fbound%lb)
        endif
        if(chg)then
          call qfracsave(fbound%lb,dsave,nvar,.false.)
        endif
        ls=fbound%lb+1
      else
        ls=fbound%lb
      endif
      if(fbound%fe .eq. 0.d0)then
        le1=min(nlat-1,fbound%le)
        call tturne1(trans,cod,beam,srot,
     $       iae,idp,plot,sol,rt,optics,
     $       ls,le1)
        if(plot)then
          call tfsetplot(trans,cod,beam,fbound%lb,
     $         le1+1,iae%iatr,iae%iacod,.false.,idp)
c          write(*,'(a,l2,i5,1p6g15.7)')'tturne0-le1 ',codplt,le1+1,
c     $       twiss(le1+1,idp,mfitnx)/m_2pi
        endif
      else
        if(ls .lt. fbound%le)then
          call tturne1(trans,cod,beam,srot,
     $         iae,idp,plot,sol,rt,optics,
     $         ls,fbound%le-1)
        endif
        call compelc(fbound%le,cmp)
        call qfracsave(fbound%le,dsave,nvar,.true.)
        call qfracseg(cmp,cmp,0.d0,fbound%fe,chg,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(fbound%le,irtc)
        else
          call tturne1(trans,cod,beam,srot,
     $         iae,idp,plot,sol,rt,optics,
     $         fbound%le,min(nlat-1,fbound%le))
        endif
        if(chg)then
          call qfracsave(fbound%le,dsave,nvar,.false.)
        endif
        if(plot)then
          call tfsetplot(trans,cod,beam,0,
     $         nlat,iae%iatr,iae%iacod,.false.,idp)
        endif
      endif
      if(plot)then
        call tfadjustn(idp,mfitnx)
        call tfadjustn(idp,mfitny)
        xb=fbound%lb+fbound%fb
        xe=fbound%le+fbound%fe
        do l=1,nlat-1
          xp=min(xe,max(xb,tffselmoffset(l)))
          if(xp .ne. dble(l))then
            lx=int(xp)
            fr=xp-lx
 8101       if(fr .eq. 0.d0)then
              if(iae%iatr .ne. 0)then
                if(iae%iatr .gt. 0)then
                  if(l .ge. fbound%lb .and. l .le. fbound%le)then
                    call tflocal(klist(iae%iatr+l))
                  endif
                  klist(iae%iatr+l)=ktfcopy1(klist(iae%iatr+lx))
                endif
                if(iae%iacod .gt. 0)then
                  if(l .ge. fbound%lb .and. l .le. fbound%le)then
                    call tflocal(klist(iae%iacod+l))
                  endif
                  klist(iae%iacod+l)=ktfcopy(klist(iae%iacod+lx))
                endif
              endif
              if(plot)then
                twiss(l,idp,mfitdx:mfitddp)=twiss(lx,idp,mfitdx:mfitddp)
                if(irad .ge. 12 .and. beamplt)then
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
              call compelc(lx,cmp)
              call qfracsave(lx,dsave,nvar,.true.)
              call qfracseg(cmp,cmp,fra,frb,chg,irtc)
              if(irtc .ne. 0)then
                call tffserrorhandle(l,irtc)
              else
                if(.not. chg)then
                  fr=0.d0
                  go to 8101
                endif
                cod1=0.d0
                if(iae%iatr .ne. 0)then
                  if(iae%iatr .gt. 0 .and.
     $                 ktflistq(dlist(iae%iatr+lx),kli))then
                    call tfl2m(kli,trans1,6,6,.false.)
                  else
                    call tinitr(trans1)
                  endif
                  if(iae%iacod .gt. 0 .and.
     $                 tfreallistq(dlist(iae%iacod+lx),klir))then
                    cod1=klir%rbody(1:6)
                  endif
                else
                  call tinitr(trans1)
                endif
                trans1(:,7:12)=0.d0
                if(plot)then
                  cod1=twiss(lx,idp,mfitdx:mfitddp)
                  if(irad .ge. 12 .and. beamplt)then
                    beam1(1:21)=beamsize(:,lx)
                    beam1(22:42)=0.d0
                  endif
                endif
                sol1=.false.
                cp0=codplt
                codplt=.false.
                int0=calint
                calint=.false.
                call tturne1(trans1,cod1,beam1,srot,
     $               iaez,idp,.false.,sol1,rt,
     $               optics,
     $               lx,lx)
                codplt=cp0
                calint=int0
              endif
              call qfracsave(lx,dsave,nvar,.false.)
              call tfsetplot(trans1,cod1,beam1,lx,
     $             l,iae%iatr,iae%iacod,
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
      lv=itfdownlevel()
      return
      end

      subroutine tturne1(trans,cod,beam,srot,
     $     iae,idp,plot,sol,rt,optics,
     $     ibegin,iend)
      use kyparam
      use tfstk
      use tffitcode
      use ffs, only: gettwiss,limitcod
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      use tfcsi, only:icslfno
      use ffs_seg
      use temw,only:tsetr0,calint,iaemit
      implicit none
      real*8 , parameter:: demax=.5d0,tapmax=0.3d0
      type (iaemit) ,intent(in):: iae
      type (sad_comp), pointer :: cmp
      type (sad_dlist), pointer :: lsegp
      integer*8 kbmzi,lp
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 bmir(6,6),bmi(21),bmh(21),trans1(6,6)
      real*8 psi1,psi2,apsi1,apsi2,alid,alid1,
     $     r,dir,al,alib,dtheta,theta0,ftable(4),
     $     fb1,fb2,ak0,ak1,rtaper,als
      integer*4 ,intent(in):: idp,ibegin,iend
      integer*4 l,ld,lele,mfr,ke,irtc,i,l1
      logical*4 ,intent(in):: plot,rt,optics
      logical*4 ,intent(inout):: sol
      logical*4 bmaccum,plotib,next,seg,coup,err,inin
      integer*8 ,save:: kbmz=0
      if(kbmz .eq. 0)then
        kbmz=ktadaloc(0,6)
        kbmzi=ktraaloc(-1,6)
        do i=1,6
          klist(kbmz+i)=ktflist+ktfcopy1(kbmzi)
        enddo
      endif
      bmaccum=.false.
      plotib=plot .and. iae%iabmi .ne. 0
      alid=0.d0
      call tsetdvfs
      call tesetdv(cod(6))
      bradprev=0.d0
      do l=ibegin,iend
        call tsetr0(trans,cod,0.d0,0.d0)
        next=inext(l) .ne. 0
        inin=.false.
        do i=1,6
          if(ktfenanq(cod(i)))then
            inin=.true.
            cod(i)=0.d0
          endif
        enddo
        if(inin)then
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
          if(iae%iatr .ne. 0)then
            if(iae%iatr .gt. 0)then
              dlist(iae%iatr+l)=
     $             dtfcopy1(kxm2l(trans,6,6,6,.false.))
            endif
            if(iae%iacod .gt. 0)then
              dlist(iae%iacod+l)=
     $             dtfcopy1(kxm2l(cod,0,6,1,.false.))
            endif
          endif
c            if(l .eq. 1)then
c              r=1.d0
c            else
c              r=gammab(l-1)/gammab(l)
c            endif
c          write(*,*)'tturne1-setetwiss ',l,beamplt,idp,cod(5)
          call tsetetwiss(trans,cod,beam,ibegin,l,idp)
c          write(*,*)'tturne1-setetwiss-end ',l,plot,cod(5)
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
        seg=tcheckseg(cmp,lele,al,lsegp,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(l,irtc)
          go to 1010
        endif
        dir=direlc(l)
        if(calint)then
          als=al
          if(als .ne. 0.d0)then
            if(lele .eq. icDRFT)then
              alib=als*.25d0+alid
              alid1=als*.25d0
            else
              alib=als*.5d0+alid
              alid1=als*.5d0
            endif
            call tintrb(trans,cod,beam,bmi,alib,alid,optics,l)
            alid=alid1
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
                dlist(iae%iabmi+l)=
     $               dtfcopy1(kxm2l(bmir,6,6,6,.false.))
                bmaccum=.false.
              endif
            endif
          elseif(plotib)then
            klist(iae%iabmi+l)=ktflist+ktfcopy1(kbmz)
          endif
        endif
c        WRITE(*,*)lele,' ',PNAME(ILIST(2,LATT(L)))(1:16)
c        if(mod(l,100) .eq. 0)then
c          write(*,*)'tturne1-l ',l,lele
c        endif
c        go to (1100,1200,1010,1400,1010,1600,1010,1600,1010,1600,
c     $       1010,1600,1010,1010,1010,1010,1010,1010,1010,3000,
c     $       3100,3200,1010,1010,1010,1010,1010,1010,1010,1010,
c     $       4100,4200,4300,4400,4500),lele
c        go to 5000

        select case (lele)
        case (icDRFT)
          if(calint)then
            if(al .ne. 0.d0)then
              al=al*.5d0
              call tdrife(trans,cod,beam,srot,al,0.d0,0.d0,0.d0,0.d0,
     $             .true.,.false.,irad)
              call tintrb(trans,cod,beam,bmi,al,al*.5d0,optics,l)
              if(plotib)then
                bmi=bmi*0.5d0
                bmh=bmh+bmi
                call tconvbm(bmh,bmir)
                bmh=bmi
                dlist(iae%iabmi+l)=
     $               dtfcopy1(kxm2l(bmir,6,6,6,.false.))
              endif
            endif
          endif
          call tdrife(trans,cod,beam,srot,al,0.d0,0.d0,0.d0,0.d0,
     $         .true.,.false.,irad)

        case (icBEND)
          if(dir .gt. 0.d0)then
            psi1=cmp%value(ky_E1_BEND)
            psi2=cmp%value(ky_E2_BEND)
            apsi1=cmp%value(ky_AE1_BEND)
            apsi2=cmp%value(ky_AE2_BEND)
            fb1=cmp%value(ky_F1_BEND)
     $           +cmp%value(ky_FB1_BEND)
            fb2=cmp%value(ky_F1_BEND)
     $           +cmp%value(ky_FB2_BEND)
          else
            psi1=cmp%value(ky_E2_BEND)
            psi2=cmp%value(ky_E1_BEND)
            apsi1=cmp%value(ky_AE2_BEND)
            apsi2=cmp%value(ky_AE1_BEND)
            fb2=cmp%value(ky_F1_BEND)
     $           +cmp%value(ky_FB1_BEND)
            fb1=cmp%value(ky_F1_BEND)
     $           +cmp%value(ky_FB2_BEND)
          endif
          dtheta=cmp%value(ky_DROT_BEND)
          theta0=cmp%value(ky_ROT_BEND)
          ak0=cmp%value(ky_K0_BEND)
     $         +cmp%value(ky_ANGL_BEND)
          ak1=cmp%value(ky_K1_BEND)
c          if(l .eq. 4551)then
c            write(*,*)'ttrune-icBend ',ak0,cod(6),
c     $           gettwiss(mfitddp,nextl(l))
c          endif
          if(radcod .and. radtaper)then
            if(rt)then
              l1=nextl(l)
              rtaper=((4.d0+3.d0*cod(6)+gettwiss(mfitddp,l1))*.25d0-dp0)
            else
              rtaper=(1.d0-dp0+cod(6))
            endif
            rtaper=min(1.d0+tapmax,max(1.d0-tapmax,rtaper))
            ak0=ak0*rtaper
            ak1=ak1*rtaper
          endif
          call tbende(trans,cod,beam,srot,al,
     1         min(pi2,max(-pi2,ak0)),
     $         cmp%value(ky_ANGL_BEND),
     $         psi1,psi2,apsi1,apsi2,ak1,
     1         cmp%value(ky_DX_BEND),
     $         cmp%value(ky_DY_BEND),theta0,dtheta,
     $         fb1,fb2,
     $         nint(cmp%value(ky_FRMD_BEND)),
     $         cmp%value(ky_FRIN_BEND) .eq. 0.d0,
     $         cmp%value(ky_EPS_BEND),
     1         cmp%value(ky_RAD_BEND) .eq. 0.d0,.true.,
     $         next,l)

        case (icQUAD)
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
              rtaper=((2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0-dp0)
            else
              rtaper=1.d0-dp0+cod(6)
            endif
            ak1=ak1*min(1.d0+tapmax,max(1.d0-tapmax,rtaper))
          endif
          call tsetfringepe(cmp,icQUAD,ftable)
          call tquade(trans,cod,beam,srot,al,ak1,0.d0,
     $         cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1         cmp%value(ky_ROT_QUAD),
     $         cmp%value(ky_RAD_QUAD) .eq. 0.d0 .and. al .ne. 0.d0,
     1         cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $         ftable(1),ftable(2),ftable(3),ftable(4),
     $         mfr,cmp%value(ky_EPS_QUAD),
     $         cmp%value(ky_KIN_QUAD) .eq. 0.d0,
     $         cmp%value(ky_CHRO_QUAD) .ne. 0.d0)

        case (icSEXT,icOCTU,icDECA,icDODECA)
          ak1=cmp%value(ky_K_THIN)
          if(radcod .and. radtaper)then
            if(rt)then
              l1=nextl(l)
              rtaper=((2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0-dp0)
            else
              rtaper=(1.d0-dp0+cod(6))
            endif
            ak1=ak1*min(1.d0+tapmax,max(1.d0-tapmax,rtaper))
          endif
          call tthine(trans,cod,beam,srot,lele,al,ak1,
     1         cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     $         cmp%value(ky_ROT_THIN),.false.)

        case(icSOL)
          call tsole(trans,cod,beam,srot,l,ke,sol,
     1         iae,idp,plot,rt)
          alid=0.d0

        case (icST)
          write(*,*)'Use BEND with ANGLE=0 for ST.'
          call forcesf()

        case(icMULT)
          rtaper=1.d0
          if(radcod .and. radtaper)then
            if(rt)then
              l1=nextl(l)
              rtaper=(2.d0+cod(6)+gettwiss(mfitddp,l1))*.5d0-dp0
            else
              rtaper=1.d0-dp0+cod(6)
            endif
          endif
          rtaper=min(1.d0+tapmax,max(1.d0-tapmax,rtaper))
          if(seg)then
c            call tfevals('Print["PROF-TTE-0: ",LINE["PROFILE","Q1"]]',
c     $       kxx,irtc)
            call tmulteseg(trans,cod,beam,srot,
     $           l,cmp,0.d0,lsegp,.true.,rtaper)
c            call tfevals('Print["PROF-TTE-1: ",LINE["PROFILE","Q1"]]',
c     $       kxx,irtc)
          else
            call tmulte1(trans,cod,beam,srot,l,cmp,0.d0,.true.,rtaper)
          endif

        case (icCAVI)
          mfr=nint(cmp%value(ky_FRMD_CAVI))
          if(dir .gt. 0.d0)then
          else
            mfr=mfr*(11+mfr*(2*mfr-9))/2
          endif
c     write(*,*)'tturne-tcave',cod
          call tcave(trans,cod,beam,srot,l,al,
     1         cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI),
     $         cmp%value(ky_HARM_CAVI),
     1         cmp%value(ky_PHI_CAVI),cmp%value(ky_FREQ_CAVI),
     1         cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $         cmp%value(ky_ROT_CAVI),
     $         cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $         cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $         cmp%value(ky_FRIN_CAVI) .eq. 0.d0,mfr,
     $         cmp%value(ky_APHI_CAVI) .ne. 0.d0,
     $         ld)
c     write(*,*)'tturne-tcave-1',cod

        case (icTCAV)
          call ttcave(trans,cod,beam,srot,al,
     1         cmp%value(ky_K0_TCAV),cmp%value(ky_HARM_TCAV),
     1         cmp%value(ky_PHI_TCAV),cmp%value(ky_FREQ_TCAV),
     1         cmp%value(ky_DX_TCAV),cmp%value(ky_DY_TCAV),
     $         cmp%value(ky_ROT_TCAV),
     $         cmp%value(ky_RAD_TCAV) .eq. 0.d0)

        case (icMAP)
          if(optics)then
            call qemap(trans1,cod,l,coup,err)
            trans(:,1:6)=matmul(trans1,trans(:,1:6))
c            call tmultr(trans,trans1,6)
          else
            call temape(trans,cod,beam,l)
          endif

        case(icINS)
          call tinse(trans,cod,beam,cmp%value(ky_DIR_INS+1))

        case (icCOORD)
          call tcoorde(trans,cod,beam,
     1         cmp%value(ky_DX_COORD),cmp%value(ky_DY_COORD),
     $         cmp%value(ky_DZ_COORD),cmp%value(ky_CHI1_COORD),
     $         cmp%value(ky_CHI2_COORD),cmp%value(ky_CHI3_COORD),
     1         cmp%value(ky_DIR_COORD) .eq. 0.d0)

        case default
        end select
 1010   continue
        call limitcod(cod)
      enddo
      call limitnan(cod,-1.d10,1.d10)
c      call tfmemcheckprint('tturne-end0',0,.true.,irtc)
      if(calint)then
        if(alid .ne. 0.d0)then
          call tintrb(trans,cod,beam,bmi,alid,alid,optics,iend+1)
          alid=0.d0
        else
          bmi=0.d0
        endif
        if(plotib)then
          if(bmaccum)then
            bmi=bmi+bmh
c            call tadd(bmh,bmi,bmi,21)
          endif
          call tconvbm(bmi,bmir)
          dlist(iae%iabmi+iend+1)=
     $         dtfcopy1(kxm2l(bmir,6,6,6,.false.))
        endif
      endif
c      call tfmemcheckprint('tturne-end1',0,.true.,irtc)
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
      use temw, only:toln
      implicit none
      integer*4 idp,m,l
      real*8 phi0
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
      real*8 trans(6,12),cod(6),beam(42)
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
      use temw, only:tfetwiss,etwiss2ri,normali,ri,toln,tinv6,
     $     beamplt
      use macmath
      implicit none
      integer*4 l,idp,lorg,l0
      real*8 trans(6,6),ti(6,6),twi(ntwissfun),cod(6),beam(42),gr
      logical*4 norm
      if(trpt)then
        gr=gammab(l)/gammab(max(1,lorg-1))
        ti=tinv6(trans*sqrt(gr))
      else
        ti=tinv6(trans)
      endif
      if(lorg .le. 1)then
        ti=matmul(ri,ti)
        norm=normali
        l0=1
      else
        l0=lorg
        twi=twiss(lorg,idp,1:ntwissfun)
        ti=matmul(etwiss2ri(twi,norm),ti)
      endif
      twi=tfetwiss(ti,cod,norm)
      if(l .eq. 1)then
        twi(mfitnx)=0.d0
        twi(mfitny)=0.d0
        twi(mfitnz)=0.d0
      else
        if(twi(mfitnx) .lt. -toln)then
          twi(mfitnx)=twiss(l0,idp,mfitnx)+twi(mfitnx)+pi2
        else
          twi(mfitnx)=twiss(l0,idp,mfitnx)+twi(mfitnx)
        endif
        if(twi(mfitny) .lt. -toln)then
          twi(mfitny)=twiss(l0,idp,mfitny)+twi(mfitny)+pi2
        else
          twi(mfitny)=twiss(l0,idp,mfitny)+twi(mfitny)
        endif
        twi(mfitnz)=twiss(l0,idp,mfitnz)+twi(mfitnz)
      endif
      twiss(l,idp,1:ntwissfun)=twi
c      write(*,'(a,2i5,1p6g15.7)')'setetwiss-end ',l,l0,
c     $     twi(mfitnx)/m_2pi,toln
      if(irad .ge. 12 .and. beamplt)then
        beamsize(:,l)=beam(1:21)+beam(22:42)
      endif
      return
      end

      subroutine checketwiss(trans,tw1)
      use ffs_pointer
      use tffitcode
      use tmacro
      use temw,only:r,etwiss2ri
      implicit none
      integer*4 i
      real*8 tw1(ntwissfun),trans(6,6),ti(6,6)
      logical*4 normal
      ti=matmul(etwiss2ri(tw1,normal),matmul(trans,r))
c      ti=r
c      call tmultr(ti,trans,6)
c      call tmultr(ti,ra,6)
      write(*,*)'checketwiss ',tw1(mfitdetr)
      do i=1,6
        write(*,'(1p6g15.7)')ti(i,:)
      enddo
      return
      end

      subroutine tesetdv(dp)
      use tmacro
      use mathfun
      implicit none
      real*8 ,intent(in):: dp
      real*8 p1
      p1=p0*(1.d0+dp)
      h1emit=p2h(p1)
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
      use tmacro
      use mathfun
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
      use tmacro
      use mathfun
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

      subroutine tmulteseg(trans,cod,beam,srot,
     $     l,cmp,bzs,lsegp,enarad,rtaper)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer, only:tsetfringepe
      use tffitcode
      use sad_main
      implicit none
      type (sad_comp) :: cmp
      type (sad_dlist) :: lsegp
      type (sad_dlist), pointer :: lal,lk
      type (sad_rlist), pointer :: lak,lkv
      real*8 :: rsave(cmp%ncomp2)
      real*8 trans(6,12),cod(6),beam(42),srot(3,9),rtaper,bzs
      integer*4 i,nseg,i1,i2,istep,k,l,k1,k2,nk
      integer*8 kk
      integer*4 , parameter :: nc=ky_PROF_MULT-1
      logical*4 , intent(in)::enarad
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
        do k=1,nc
          if(integv(k,icMULT))then
            cmp%value(k)=rsave(k)*lak%rbody(i)
          else
            cmp%value(k)=rsave(k)
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
        call tmulte1(trans,cod,beam,srot,l,cmp,bzs,enarad,rtaper)
      enddo
      cmp%value(1:nc)=rsave(1:nc)
      return
      end

      subroutine tmulte1(trans,cod,beam,srot,l,cmp,bzs,enarad,rtaper)
      use kyparam
      use tfstk
      use tffitcode
      use ffs, only: gettwiss
      use ffs_pointer,only:tsetfringepe
      use ffs_flag
      use tmacro
      use sad_main
      implicit none
      type (sad_comp) ,intent(in):: cmp
      integer*4 ,intent(in):: l
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      integer*4 mfr
      real*8 ,intent(in):: bzs,rtaper
      real*8 phi,al,ftable(4),
     $     psi1,psi2,apsi1,apsi2,fb1,fb2,chi1,chi2
      logical*4 , intent(in):: enarad
      al=cmp%value(ky_L_MULT)
      phi=cmp%value(ky_ANGL_MULT)
      mfr=nint(cmp%value(ky_FRMD_MULT))
      if(cmp%ori)then
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
      call tsetfringepe(cmp,icMULT,ftable)
      call tmulte(trans,cod,beam,srot,l,al,
     $     cmp%value(ky_K0_MULT),
     $     bzs,
     $     phi,psi1,psi2,apsi1,apsi2,
     1     cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $     cmp%value(ky_DZ_MULT),
     $     chi1,chi2,cmp%value(ky_ROT_MULT),
     $     cmp%value(ky_DROT_MULT),
     $     cmp%value(ky_EPS_MULT),
     $     enarad .and. cmp%value(ky_RAD_MULT) .eq. 0.d0,
     $     cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $     ftable(1),ftable(2),ftable(3),ftable(4),
     $     mfr,fb1,fb2,
     $     cmp%value(ky_K0FR_MULT) .eq. 0.d0,
     $     cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT),
     $     cmp%value(ky_HARM_MULT),
     $     cmp%value(ky_PHI_MULT),cmp%value(ky_FREQ_MULT),
     $     cmp%value(ky_W1_MULT),rtaper,
     $     cmp%value(ky_APHI_MULT) .ne. 0.d0)
      return
      end

      logical*4 function nanm(a)
      use tfstk, only:ktfenanzeroq
      real*8 , intent(in):: a(6,6)
      integer*4 i,j
      do i=1,6
        do j=1,6
          if(ktfenanzeroq(a(i,j)))then
            nanm=.true.
            return
          endif
        enddo
      enddo
      nanm=.false.
      return
      end
