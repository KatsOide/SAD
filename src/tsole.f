      subroutine tsole(trans,cod,beam,k,ke,sol,
     1     iatr,iacod,iabmi,plot,rt)
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
      integer*4 k,ke,i,l
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
            rtaper=(2.d0+cod(6)+gettwiss(mfitddp,l+1))*.5d0
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
            r=gammab(l)/gammab(l+1)
            twiss(l+1,0,mfitdx )=cod(1)
            twiss(l+1,0,mfitdpx)=cod(2)*r
            twiss(l+1,0,mfitdy )=cod(3)
            twiss(l+1,0,mfitdpy)=cod(4)*r
            twiss(l+1,0,mfitdz )=cod(5)
            twiss(l+1,0,mfitddp)=cod(6)*r
            beamsize(:,l+1)=beam
          endif
          if(calint .and. iabmi .ne. 0)then
            if(iabmilz .eq. 0)then
              bmir=0.d0
              iabmilz=ktfaddr(kxm2l(bmir,6,6,6,.false.))
            endif
            call tflocal(klist(iabmi+l))
            klist(iabmi+l)=ktflist+ktfcopy1(iabmilz)
          endif
        endif
130   continue
      return
      end

      subroutine tsole1(trans,cod,beam,l,rtaper,enarad,qsol)
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tmacro
      use sad_main
      implicit none
      integer*4 l,ld,lt,mfr,kb
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
      al=cmp%value(1)
      bzs=tfbzs(l,kb)
      if(lt .eq. icDRFT)then
        call tdrife(trans,cod,beam,al,
     $       bzs,0.d0,0.d0,.true.,
     $       enarad .and. cmp%value(kytbl(kwRAD,icDRFT)) .eq. 0.d0,
     $       calpol,irad,ld)
      elseif(lt .eq. icBEND)then
        theta=cmp%value(kytbl(kwROT,icBEND))
     $       +cmp%value(kytbl(kwDROT,icBEND))
        phi=cmp%value(2)+cmp%value(kytbl(kwK0,icBEND))
        phiy= phi*cos(theta)
        phix= phi*sin(theta)
        call tdrife(trans,cod,beam,al,
     $       bzs,phiy,phix,.true.,
     $       enarad .and. cmp%value(kytbl(kwRAD,icBEND)) .eq. 0.d0,
     $       calpol,irad,ld)
      elseif(lt .eq. icQUAD)then
        dir=direlc(l) .gt. 0.d0
        if(dir)then
          mfr=nint(cmp%value(12))
        else
          mfr=nint(cmp%value(12))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        if(enarad)then
          radlvl=cmp%value(7)
        else
          radlvl=1.d0
        endif
        ak1=cmp%value(kytbl(kwK1,icQUAD))
        call tsetfringepe(cmp,icQUAD,direlc(l),ftable)
        call tquase(trans,cod,beam,
     $       al,ak1,bzs,
     $       cmp%value(5),cmp%value(6),cmp%value(4),
     1       radlvl,cmp%value(9) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,cmp%value(13),l,dir,ld)
      elseif(lt .eq. icMULT)then
        dir=direlc(l).gt. 0.d0
        phi=cmp%value(kytbl(kwANGL,icMULT))
        mfr=nint(cmp%value(kytbl(kwFRMD,icMULT)))
        if(dir)then
          psi1=cmp%value(kytbl(kwE1,icMULT))
          psi2=cmp%value(kytbl(kwE2,icMULT))
          apsi1=cmp%value(kytbl(kwAE1,icMULT))
          apsi2=cmp%value(kytbl(kwAE2,icMULT))
          fb1=cmp%value(kytbl(kwFB1,icMULT))
          fb2=cmp%value(kytbl(kwFB2,icMULT))
          chi1=cmp%value(kytbl(kwCHI1,icMULT))
          chi2=cmp%value(kytbl(kwCHI2,icMULT))
        else
          mfr=mfr*(11+mfr*(2*mfr-9))/2
          psi1=cmp%value(kytbl(kwE2,icMULT))
          psi2=cmp%value(kytbl(kwE1,icMULT))
          apsi1=cmp%value(kytbl(kwAE2,icMULT))
          apsi2=cmp%value(kytbl(kwAE1,icMULT))
          fb2=cmp%value(kytbl(kwFB1,icMULT))
          fb1=cmp%value(kytbl(kwFB2,icMULT))
          chi1=-cmp%value(kytbl(kwCHI1,icMULT))
          chi2=-cmp%value(kytbl(kwCHI2,icMULT))
        endif
        call tsetfringepe(cmp,icMULT,direlc(l),ftable)
        call tmulte(trans,cod,beam,l,al,
     $       cmp%value(kytbl(kwK0,icMULT)),bzs,
     $       phi,psi1,psi2,apsi1,apsi2,
     1       cmp%value(3),cmp%value(4),cmp%value(5),
     $       chi1,chi2,cmp%value(8),
     $       cmp%value(kytbl(kwDROT,icMULT)),
     $       cmp%value(9),
     $       cmp%value(kytbl(kwRAD,icMULT)) .eq. 0.d0 .and. enarad,
     $       cmp%value(11) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,fb1,fb2,
     $       cmp%value(kytbl(kwK0FR,icMULT)) .eq. 0.d0,
     $       cmp%value(15),cmp%value(16),cmp%value(17),cmp%value(18),
     $       cmp%value(kytbl(kwW1,icMULT)),rtaper,
     $       cmp%value(kytbl(kwAPHI,icMULT)) .ne. 0.d0,
     $       ld)
      elseif(lt .eq. icSOL)then
        enarad1=enarad .and. cmp%value(kytbl(kwRAD,icSOL)) .eq. 0.d0
        if(rlist(idval(ld)+kytbl(kwBND,icSOL)) .ne. 0.d0)then
          ent=direlc(l) .gt. 0.d0 .and. l .eq. kb
     $         .or. direlc(l) .lt. 0.d0 .and. l .ne. kb
          if(calpol)then
            cod1=cod
          endif
          bzs0=tfbzs(l-1,kb)
          if(enarad1 .and. .not. ent)then
            f1=cmp%value(kytbl(kwF1,icSOL))
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,-bzs0,0.d0,f1,brhoz)
            endif
          endif
          if(cmp%value(kytbl(kwFRIN,icSOL)) .eq. 0.d0)then
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
            f1=cmp%value(kytbl(kwF1,icSOL))
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,0.d0,bzs,f1,brhoz)
            endif
          endif
        else
          bzs0=tfbzs(l-1,kb)
          if(enarad1)then
            f1=cmp%value(kytbl(kwF1,icSOL))
            if(f1 .ne. 0.d0)then
              call trades(trans,beam,cod,bzs0,bzs,f1,brhoz)
            endif
          endif
          if(cmp%value(kytbl(kwFRIN,icSOL)) .eq. 0.d0)then
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
