      subroutine tturne(latt,trans,cod,beam,
     $     twiss,bsize,gammab,
     $     iatr,iacod,iabmi,ndim,plot,update,rt)
      use touschek_table
      use tfstk
      use tffitcode
      use ffs, only: gettwiss
      implicit none
      include 'inc/TMACRO1.inc'
      type (sad_list), pointer :: kli
      real*8 codmax,demax
      parameter (codmax=1.d4,demax=.5d0)
      integer*8 iatr,iacod,iabmi
      integer*4 ndim,lbegin,lend,ls,l,nvar,lx
      integer*4 latt(2,nlat)
      real*8 trans(6,12),cod(6),beam(42)
      real*8 trans1(6,12),cod1(6),beam1(42)
      real*8 twiss(nlat,-ndim:ndim,*),bsize(21,nlat),
     $     gammab(nlat),vsave(100)
      real*8 pgev00,frbegin,frend,phis,dvcphic,r,
     $     xp,xb,xe,fr,fra,frb,tffselmoffset,z0
      logical*4 sol,plot,update,chg,sol1,cp0,int0,rt
      pgev00=pgev
      vc0=0.d0
      u0=0.d0
      hvc0=0.d0
      vccos=0.d0
      vcsin=0.d0
      z0=cod(5)
      sol=.false.
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
      call tffsbound(nlat,latt,lbegin,frbegin,lend,frend)
      if(frbegin .ne. 0.d0)then
        call qfracsave(latt(1,lbegin),vsave,nvar,ideal,.true.)
        call qfraccomp(latt(1,lbegin),frbegin,1.d0,ideal,chg)
        call tturne1(latt,trans,cod,beam,twiss,bsize,gammab,
     $       iatr,iacod,iabmi,ndim,plot,sol,rt,lbegin,lbegin)
        if(chg)then
          call qfracsave(latt(1,lbegin),vsave,nvar,ideal,.false.)
        endif
        ls=lbegin+1
      else
        ls=lbegin
      endif
      call tturne1(latt,trans,cod,beam,twiss,bsize,gammab,
     $     iatr,iacod,iabmi,ndim,plot,sol,rt,ls,lend-1)
      if(frend .ne. 0.d0)then
        call qfracsave(latt(1,lend),vsave,nvar,ideal,.true.)
        call qfraccomp(latt(1,lend),0.d0,frend,ideal,chg)
        call tturne1(latt,trans,cod,beam,twiss,bsize,gammab,
     $       iatr,iacod,iabmi,ndim,plot,sol,rt,lend,lend)
        if(chg)then
          call qfracsave(latt(1,lend),vsave,nvar,ideal,.false.)
        endif
      endif
      if(plot)then
        call tfsetplot(trans,cod,beam,twiss,bsize,ndim,
     $       lend,iatr,iacod,.false.,
     $       gammab(lend-1)/gammab(lend))
        xb=lbegin+frbegin
        xe=lend+frend
        do l=1,nlat
          xp=min(xe,max(xb,tffselmoffset(latt,l,nlat)))
          if(xp .ne. dble(l))then
            lx=int(xp)
            fr=xp-lx
 8101       if(fr .eq. 0.d0)then
              if(iatr .ne. 0)then
                if(iatr .gt. 0)then
                  if(l .ge. lbegin .and. l .le. lend)then
                    call tflocal(klist(iatr+l))
                  endif
                  klist(iatr+l)=ktfcopy1(klist(iatr+lx))
                endif
                if(iacod .gt. 0)then
                  if(l .ge. lbegin .and. l .le. lend)then
                    call tflocal(klist(iacod+l))
                  endif
                  klist(iacod+l)=ktfcopy(klist(iacod+lx))
                endif
              endif
              if(codplt)then
                twiss(l,0,mfitdx:mfitddp)=twiss(lx,0,mfitdx:mfitddp)
                bsize(:,l)=bsize(:,lx)
              endif
            else
              if(lx .eq. lbegin)then
                fra=frbegin
                frb=max(fr,fra)
              elseif(lx .eq. lend)then
                fra=0.d0
                frb=min(frend,fr)
              else
                fra=0.d0
                frb=fr
              endif
              if(fra .eq. frb)then
                fr=0.d0
                go to 8101
              endif
c     below is incorrect for fra <> 0
              call qfracsave(latt(1,lx),vsave,nvar,ideal,.true.)
              call qfraccomp(latt(1,lx),fra,frb,ideal,chg)
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
              call tclr(trans1(1,7),36)
              if(codplt)then
                cod1=twiss(lx,0,mfitdx:mfitddp)
                beam1(1:21)=bsize(:,lx)
                if(calint)then
                  beam1(22:42)=bsize(:,lx)
                endif
              endif
              sol1=.false.
              cp0=codplt
              codplt=.false.
              int0=calint
              calint=.false.
              call tturne1(latt,trans1,cod1,beam1,twiss,bsize,gammab,
     $             int8(0),int8(0),int8(0),ndim,.false.,sol1,rt,
     $             lx,lx)
              codplt=cp0
              calint=int0
              call qfracsave(latt(1,lx),vsave,nvar,ideal,.false.)
              call tfsetplot(trans1,cod1,beam1,twiss,bsize,ndim,
     $             l,iatr,iacod,
     $             l .ge. lbegin .and. l .le. lend,
     $             gammab(lx)/(gammab(lx)*(1.d0-frb)+gammab(lx+1)*frb))
            endif
          endif
        enddo
      elseif(radtaper .and. codplt)then
        if(lend .eq. 1)then
          r=1.d0
        else
          r=gammab(lend-1)/gammab(lend)
        endif
        twiss(lend,0,mfitddp)=cod(6)*r
      endif
      if(vc0 .ne. 0.d0 .and. update)then
        vceff=sign(sqrt(vccos**2+vcsin**2),vc0)
        if(trpt)then
          trf0=0.d0
          vcphic=0.d0
          vcalpha=1.d0
        else
          dvcphic=atan2(vcsin,vccos)
          if(abs(dvcphic) .lt. pi*.5d0)then
            vcphic=vcphic+dvcphic
          endif
          vcalpha=vceff/vc0
          phis=asin(min(1.d0,max(-1.d0,u0*p0*amass/sign(vceff,vccos))))
          trf0=phis*c*p0/h0/omega0/hvc0*vceff
        endif
        call RsetGL1('DTSYNCH',trf0)
        call RsetGL1('PHICAV',vcphic)
        call RsetGL1('EFFVCRATIO',vcalpha)
      endif
      if(pgev00 .ne. pgev)then
        pgev=pgev00
        call tphyzp
      endif
      return
      end

      subroutine tturne1(latt,trans,cod,beam,twiss,bsize,gammab,
     $     iatr,iacod,iabmi,ndim,plot,sol,rt,ibegin,iend)
      use tfstk
      use tffitcode
      use ffs, only: gettwiss
      implicit none
      include 'inc/TMACRO1.inc'
      real*8 codmax,demax
      parameter (codmax=1.d4,demax=.5d0)
      integer*8 iatr,iacod,iabmi,kbmz,kbmzi
      integer*4 ndim,latt(2,nlat),i
      real*8 trans(6,12),cod(6),beam(42),bmir(6,6),
     $     bmi(21),bmh(21)
      real*8 twiss(nlat,-ndim:ndim,*),bsize(21,nlat),gammab(nlat),
     $     psi1,psi2,apsi1,apsi2,alid,cod60,
     $     r,dir,al,alib,dtheta,theta0,
     $     phi,fb1,fb2,chi1,chi2,ak0,ak1,rtaper
      integer*4 l,ld,lele,lp,kl,
     $     mfr,ibegin,iend,ke
      logical*4 sol,plot,bmaccum,plotib,isnan,rt
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
      do l=ibegin,iend
        if(isnan(cod(1)) .or. isnan(cod(3)))then
          cod=0.d0
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
            if(l .eq. 1)then
              r=1.d0
            else
              r=gammab(l-1)/gammab(l)
            endif
            twiss(l,0,mfitdx )=cod(1)
            twiss(l,0,mfitdpx)=cod(2)*r
            twiss(l,0,mfitdy )=cod(3)
            twiss(l,0,mfitdpy)=cod(4)*r
            twiss(l,0,mfitdz )=cod(5)
            twiss(l,0,mfitddp)=cod(6)*r
            bsize(:,l)=beam(1:21)
          endif
        elseif(radtaper .and. codplt)then
          if(l .eq. 1)then
            r=1.d0
          else
            r=gammab(l-1)/gammab(l)
          endif
          twiss(l,0,mfitddp)=cod(6)*r
        endif
        ld=latt(1,l)
        lele=idtype(ld)
        if(ideal)then
          lp=idval(ld)
        else
          lp=latt(2,l)
        endif
        dir=rlist(latt(2,l)+ilist(1,latt(2,l)))
        kl=kytbl(kwL,lele)
        if(kl .ne. 0)then
          al=rlist(lp+kl)
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
                  call tadd(bmh,bmi,bmh,21)
                else
                  call tmov(bmi,bmh,21)
                endif
                bmaccum=.true.
              else
                if(bmaccum)then
                  call tadd(bmh,bmi,bmi,21)
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
c        WRITE(*,*)lele,' ',PNAME(LATT(1,L))(1:16)
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
              call tmul(bmi,.5d0,bmi,21)
              call tadd(bmh,bmi,bmh,21)
              call tconvbm(bmh,bmir)
              call tmov(bmi,bmh,21)
              dlist(iabmi+l)=
     $             dtfcopy1(kxm2l(bmir,6,6,6,.false.))
            endif
          endif
        endif
        cod60=cod(6)
        call tdrife(trans,cod,beam,al,0.d0,0.d0,0.d0,
     $       .true.,.false.,calpol,irad,ld)
        go to 1010
1200    continue
        if(dir .gt. 0.d0)then
          psi1=rlist(lp+kytbl(kwE1,icBEND))
          psi2=rlist(lp+kytbl(kwE2,icBEND))
          apsi1=rlist(lp+kytbl(kwAE1,icBEND))
          apsi2=rlist(lp+kytbl(kwAE2,icBEND))
          fb1=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB1,icBEND))
          fb2=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB2,icBEND))
        else
          psi1=rlist(lp+kytbl(kwE2,icBEND))
          psi2=rlist(lp+kytbl(kwE1,icBEND))
          apsi1=rlist(lp+kytbl(kwAE2,icBEND))
          apsi2=rlist(lp+kytbl(kwAE1,icBEND))
          fb2=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB1,icBEND))
          fb1=rlist(lp+kytbl(kwF1,icBEND))
     $         +rlist(lp+kytbl(kwFB2,icBEND))
        endif
        dtheta=rlist(lp+kytbl(kwDROT,icBEND))
        theta0=rlist(lp+kytbl(kwROT,icBEND))+dtheta
        ak0=rlist(lp+kytbl(kwK0,icBEND))+rlist(lp+kytbl(kwANGL,icBEND))
        ak1=rlist(lp+kytbl(kwK1,icBEND))
        if(radcod .and. radtaper)then
          if(rt)then
            ak0=ak0*(4.d0+3.d0*cod(6)+gettwiss(mfitddp,l+1))*.25d0
            ak1=ak1*(4.d0+3.d0*cod(6)+gettwiss(mfitddp,l+1))*.25d0
          else
            ak0=ak0*(1.d0+cod(6))
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        cod60=cod(6)
        call tbende(trans,cod,beam,al,
     1       min(pi2,max(-pi2,ak0)),
     $       rlist(lp+2),
     $       psi1,psi2,apsi1,apsi2,ak1,
     1       rlist(lp+kytbl(kwDX,icBEND)),
     $       rlist(lp+kytbl(kwDY,icBEND)),theta0,dtheta,
     $       fb1,fb2,
     $       nint(rlist(lp+kytbl(kwFRMD,icBEND))),
     $       rlist(lp+kytbl(kwFRIN,icBEND)) .eq. 0.d0,
     $       rlist(lp+kytbl(kwEPS,icBEND)),
     1       rlist(lp+kytbl(kwRAD,icBEND)) .eq. 0.d0,.true.,ld)
        go to 1010
1400    continue
        if(dir .gt. 0.d0)then
          mfr=nint(rlist(lp+12))
        else
          mfr=nint(rlist(lp+12))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        ak1=rlist(lp+kytbl(kwK1,icQUAD))
        if(radcod .and. radtaper)then
          if(rt)then
            ak1=ak1*(2.d0+cod(6)+gettwiss(mfitddp,l+1))*.5d0
          else
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        call tquade(trans,cod,beam,al,
     1       ak1,rlist(lp+5),rlist(lp+6),
     1       rlist(lp+4),rlist(lp+kytbl(kwRAD,icQUAD)) .eq. 0.d0,
     1       rlist(lp+9) .eq. 0.d0,
     1       rlist(lp+10),rlist(lp+11),mfr,
     $       rlist(lp+13),rlist(lp+14) .eq. 0.d0,ld)
        go to 1010
 1600   ak1=rlist(lp+2)
        if(radcod .and. radtaper)then
          if(rt)then
            ak1=ak1*(2.d0+cod(6)+gettwiss(mfitddp,l+1))*.5d0
          else
            ak1=ak1*(1.d0+cod(6))
          endif
        endif
        call tthine(trans,cod,beam,lele,al,ak1,
     1             rlist(lp+5),rlist(lp+6),rlist(lp+4),.false.,ld)
        go to 1010
 3000   call tsole(trans,cod,beam,latt,l,ke,sol,
     1       twiss,bsize,gammab,iatr,iacod,iabmi,ndim,plot)
        alid=0.d0
        go to 1010
 3100   write(*,*)'Use BEND with ANGLE=0 for ST.'
        stop
 3200   phi=rlist(lp+kytbl(kwANGL,icMULT))
        mfr=nint(rlist(lp+kytbl(kwFRMD,icMULT)))
        if(dir .gt. 0.d0)then
          psi1=rlist(lp+kytbl(kwE1,icMULT))
          psi2=rlist(lp+kytbl(kwE2,icMULT))
          apsi1=rlist(lp+kytbl(kwAE1,icMULT))
          apsi2=rlist(lp+kytbl(kwAE2,icMULT))
          fb1=rlist(lp+kytbl(kwFB1,icMULT))
          fb2=rlist(lp+kytbl(kwFB2,icMULT))
          chi1=rlist(lp+kytbl(kwCHI1,icMULT))
          chi2=rlist(lp+kytbl(kwCHI2,icMULT))
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=rlist(lp+kytbl(kwE2,icMULT))
          psi2=rlist(lp+kytbl(kwE1,icMULT))
          apsi1=rlist(lp+kytbl(kwAE2,icMULT))
          apsi2=rlist(lp+kytbl(kwAE1,icMULT))
          fb2=rlist(lp+kytbl(kwFB1,icMULT))
          fb1=rlist(lp+kytbl(kwFB2,icMULT))
          chi1=-rlist(lp+kytbl(kwCHI1,icMULT))
          chi2=-rlist(lp+kytbl(kwCHI2,icMULT))
        endif
        rtaper=1.d0
        if(radcod .and. radtaper)then
          if(rt)then
            rtaper=(2.d0+cod(6)+gettwiss(mfitddp,l+1))*.5d0
          else
            rtaper=1.d0+cod(6)
          endif
        endif
        call tmulte(trans,cod,beam,gammab,l,al,
     $       rlist(lp+kytbl(kwK0,icMULT)),
     $       0.d0,
     $       phi,psi1,psi2,apsi1,apsi2,
     1       rlist(lp+3),rlist(lp+4),rlist(lp+5),
     $       chi1,chi2,rlist(lp+8),
     $       rlist(lp+kytbl(kwDROT,icMULT)),
     $       rlist(lp+9),rlist(lp+10) .eq. 0.d0,
     $       rlist(lp+11) .eq. 0.d0,
     $       rlist(lp+12),rlist(lp+13),mfr,fb1,fb2,
     $       rlist(lp+kytbl(kwK0FR,icMULT)) .eq. 0.d0,
     $       rlist(lp+15),rlist(lp+16),rlist(lp+17),rlist(lp+18),
     $       rlist(lp+kytbl(kwW1,icMULT)),rtaper,
     $       ld)
        go to 1010
 4100   mfr=nint(rlist(lp+kytbl(kwFRMD,icCAVI)))
        if(dir .gt. 0.d0)then
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        call tcave(trans,cod,beam,gammab,l,al,
     1       rlist(lp+2),rlist(lp+3),
     1       rlist(lp+4),rlist(lp+5),
     1       rlist(lp+13),rlist(lp+14),rlist(lp+15),
     $       rlist(lp+16),rlist(lp+17),rlist(lp+18),rlist(lp+19),
     $       rlist(lp+kytbl(kwFRIN,icCAVI)) .eq. 0.d0,mfr,
     $       ld)
        go to 1010
 4200   call ttcave(trans,cod,beam,al,
     $       rlist(lp+2),rlist(lp+3),
     1       rlist(lp+4),rlist(lp+5),
     1       rlist(lp+6),rlist(lp+7),rlist(lp+8),ld)
        go to 1010
 4300   call temape(trans,cod,beam,l)
        go to 1010
 4400   call tinse(trans,cod,beam,rlist(lp+20),ld)
        go to 1010
 4500   call tcoorde(trans,cod,beam,
     1       rlist(lp+1),rlist(lp+2),rlist(lp+3),
     1       rlist(lp+4),rlist(lp+5),rlist(lp+6),
     1       rlist(lp+7) .eq. 0.d0,ld)
        go to 1010
 5000   go to 1010
 1010   continue
      enddo
      if(calint)then
        if(alid .ne. 0.d0)then
          call tintrb(trans,cod,beam,bmi,alid,alid,l)
          alid=0.d0
        else
          call tclr(bmi,21)
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

      subroutine tfsetplot(trans,cod,beam,twiss,bsize,ndim,
     $     l,iatr,iacod,local,r)
      use tfstk
      use tffitcode
      implicit none
      include 'inc/TMACRO1.inc'
      integer*8 iatr,iacod
      integer*4 l,ndim
      real*8 trans(6,12),cod(6),beam(21),bsize(21,*),r,
     $     twiss(nlat,-ndim:ndim,*)
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
        twiss(l,0,mfitdx )=cod(1)
        twiss(l,0,mfitdpx)=cod(2)*r
        twiss(l,0,mfitdy )=cod(3)
        twiss(l,0,mfitdpy)=cod(4)*r
        twiss(l,0,mfitdz )=cod(5)
        twiss(l,0,mfitddp)=cod(6)*r
        call tmov(beam,bsize(1,l),21)
      endif
      return
      end

      subroutine tesetdv(dp)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
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
      implicit none
      include 'inc/TMACRO1.inc'
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
      implicit none
      include 'inc/TMACRO1.inc'
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
      implicit none
      include 'inc/TMACRO1.inc'
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
