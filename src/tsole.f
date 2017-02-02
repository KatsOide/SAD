      subroutine tsole(trans,cod,beam,k,ke,sol,
     1     iatr,iacod,iabmi,idp,plot,rt)
      use kyparam
      use tfstk
      use tffitcode
      use ffs, only:gettwiss
      use ffs_pointer
      use sad_main
      use ffs_flag
      use tmacro
      implicit none
      real*8 conv
      parameter (conv=3.d-16)
      integer*8 iatr,iacod,iabmi,iatrl,iacodl,iabmilz
      integer*4 k,ke,i,l,idp
      real*8 trans(6,12),cod(6),beam(42),bmir(6,6),rtaper
      real*8 r
      logical*4 sol,plot,rt
      save iabmilz
      data iabmilz /0/
      sol=.true.
      do i=k+1,nlat-1
        if(idtypec(i) .eq. icSOL)then
          if(rlist(idvalc(i)+8) .ne. 0.d0)then
            ke=i
            go to 20
          endif
        endif
      enddo
      sol=.false.
      return
 20   do 130 l=k,ke
        rtaper=1.d0
        if(radtaper)then
          if(rt)then
            rtaper=(2.d0+cod(6)+gettwiss(mfitddp,nextl(l)))*.5d0
          else
            rtaper=(1.d0+cod(6))
          endif
        endif
        call tsole1(trans,cod,beam,l,rtaper,.true.,.false.)
        if(plot)then
          if(iatr .ne. 0)then
            if(iatr .gt. 0)then
              call tflocal(klist(iatr+l+1))
              iatrl=ktfaddr(kxm2l(trans,6,6,6,.false.))
              klist(iatr+l+1)=ktflist+ktfcopy1(iatrl)
            endif
            if(iacod .gt. 0)then
              call tflocal(klist(iacod+l+1))
              iacodl=ktfaddr(kxm2l(cod,0,6,1,.false.))
              klist(iacod+l+1)=ktflist+ktfcopy1(iacodl)
            endif
          endif
          if(codplt)then
            if(l .eq. 1)then
              r=1.d0
            else
              r=gammab(l)/gammab(l+1)
            endif
            call tsetetwiss(trans,cod,beam,0,l+1,idp,r)
            if(irad .gt. 6)then
              beamsize(:,l+1)=beam
            endif
          endif
          if(calint .and. iabmi .ne. 0)then
            if(iabmilz .eq. 0)then
              bmir=0.d0
              iabmilz=ktfaddr(kxm2l(bmir,6,6,6,.false.))
            endif
            call tflocal(klist(iabmi+l))
            klist(iabmi+l)=ktflist+ktfcopy1(iabmilz)
          endif
        elseif(radtaper .and. radcod)then
          if(l .eq. 1)then
            r=1.d0
          else
            r=gammab(l)/gammab(l+1)
          endif
          twiss(l+1,idp,mfitddp)=cod(6)*r
        endif
130   continue
      return
      end

      subroutine tsole1(trans,cod,beam,l,rtaper,enarad,qsol)
      use kyparam
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      implicit none
      integer*4 l,ld,lt,mfr,kb,kyl
      integer*8 lp
      type (sad_comp), pointer ::cmp
      real*8 trans(6,12),cod(6),beam(42),al,theta,
     $     phi,phix,phiy,bzs,cod1(6),trans1(6,6),trans2(6,6),
     $     tfbzs,radlvl,bzs0,fb1,fb2,chi1,chi2,
     $     psi1,psi2,apsi1,apsi2,f1,rtaper,ftable(4),ak1
      logical*4 enarad,dir,ent,qsol,coup,err,enarad1
      ld=idelc(l)
      lt=idtype(ld)
      lp=elatt%comp(l)
      call loc_comp(lp,cmp)
      kyl=kytbl(kwL,lt)
      if(kyl .eq. 0)then
        al=0.d0
      else
        al=cmp%value(kyl)
      endif
      bzs=tfbzs(l,kb)
      if(lt .eq. icDRFT)then
        call tdrife(trans,cod,beam,al,
     $       bzs,0.d0,0.d0,.true.,
     $       enarad .and. cmp%value(ky_RAD_DRFT) .eq. 0.d0,
     $       calpol,irad,ld)
      elseif(lt .eq. icBEND)then
        theta=cmp%value(ky_ROT_BEND)
     $       +cmp%value(ky_DROT_BEND)
        phi=cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND)
        phiy= phi*cos(theta)
        phix= phi*sin(theta)
        call tdrife(trans,cod,beam,al,
     $       bzs,phiy,phix,.true.,
     $       enarad .and. cmp%value(ky_RAD_BEND) .eq. 0.d0,
     $       calpol,irad,ld)
      elseif(lt .eq. icQUAD)then
        dir=direlc(l) .gt. 0.d0
        if(dir)then
          mfr=nint(cmp%value(ky_FRMD_QUAD))
        else
          mfr=nint(cmp%value(ky_FRMD_QUAD))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        if(enarad)then
          radlvl=cmp%value(ky_RAD_QUAD)
        else
          radlvl=1.d0
        endif
        ak1=cmp%value(ky_K1_QUAD)
        call tsetfringepe(cmp,icQUAD,direlc(l),ftable)
        call tquase(trans,cod,beam,
     $       al,ak1,bzs,
     $       cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     $       cmp%value(ky_ROT_QUAD),
     1       radlvl,cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,cmp%value(ky_EPS_QUAD),l,dir,ld)
      elseif(lt .eq. icMULT)then
        dir=direlc(l).gt. 0.d0
        phi=cmp%value(ky_ANGL_MULT)
        mfr=nint(cmp%value(ky_FRMD_MULT))
        if(dir)then
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
        call tsetfringepe(cmp,icMULT,direlc(l),ftable)
        call tmulte(trans,cod,beam,l,al,
     $       cmp%value(ky_K0_MULT),bzs,
     $       phi,psi1,psi2,apsi1,apsi2,
     1       cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $       cmp%value(ky_DZ_MULT),
     $       chi1,chi2,cmp%value(ky_ROT_MULT),
     $       cmp%value(ky_DROT_MULT),
     $       cmp%value(ky_EPS_MULT),
     $       cmp%value(ky_RAD_MULT) .eq. 0.d0 .and. enarad,
     $       cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,fb1,fb2,
     $       cmp%value(ky_K0FR_MULT) .eq. 0.d0,
     $       cmp%value(ky_VOLT_MULT),cmp%value(ky_HARM_MULT),
     $       cmp%value(ky_PHI_MULT),cmp%value(ky_FREQ_MULT),
     $       cmp%value(ky_W1_MULT),rtaper,
     $       cmp%value(ky_APHI_MULT) .ne. 0.d0,
     $       ld)
      elseif(lt .eq. icSOL)then
        enarad1=enarad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0
        if(rlist(idval(ld)+ky_BND_SOL) .ne. 0.d0)then
          ent=direlc(l) .gt. 0.d0 .and. l .eq. kb
     $         .or. direlc(l) .lt. 0.d0 .and. l .ne. kb
          if(calpol)then
            cod1=cod
          endif
          bzs0=tfbzs(l-1,kb)
          if(enarad1 .and. .not. ent)then
            f1=cmp%value(ky_F1_SOL)
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,-bzs0,0.d0,f1,brhoz)
            endif
          endif
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            if(ent)then
              call tsconv(trans1,cod,lp,.true.)
              call tsfrie(trans2,cod,bzs)
            else
              call tsfrie(trans1,cod,-bzs0)
              call tsconv(trans2,cod,lp,.false.)
            endif
            call tmultr5(trans1,trans2,6)
          else
            call tsconv(trans1,cod,lp,ent)
          endif            
          call tmultr5(trans,trans1,irad)
          call tmulbs(beam ,trans1,.true.,.true.)
          if(calpol)then
            call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod1)
          endif
          if(enarad1 .and. ent)then
            f1=cmp%value(ky_F1_SOL)
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,0.d0,bzs,f1,brhoz)
            endif
          endif
        else
          bzs0=tfbzs(l-1,kb)
          if(enarad1)then
            f1=cmp%value(ky_F1_SOL)
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,bzs0,bzs,f1,brhoz)
            endif
          endif
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            if(calpol)then
              cod1=cod
            endif
            call tsfrie(trans1,cod,bzs-bzs0)
            call tmultr5(trans,trans1,irad)
            call tmulbs(beam ,trans1,.true.,.true.)
            if(calpol)then
              call polpar(0,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod1)
            endif
          endif
        endif
      elseif(lt .eq. icMAP)then
        if(qsol)then
          call qemap(trans,cod,l,coup,err)
        else
          call temape(trans,cod,beam,l)
        endif
      endif
      return
      end

      subroutine qsol(trans,cod,k,coup)
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 k
      real*8 trans(4,5),cod(6),transe(6,12),beam(42)
      logical*4 coup,radtaper0
      radtaper0=radtaper
      radtaper=.false.
      call tinitr(transe)
      call tsole1(transe,cod,beam,k,1.d0,.false.,.true.)
      radtaper=radtaper0
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
