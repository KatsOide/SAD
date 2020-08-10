      subroutine tsole(trans,cod,beam,srot,k,ke,sol,iae,idp,plot,rt)
      use kyparam
      use tfstk
      use tffitcode
      use ffs, only:gettwiss,mfitddp
      use ffs_pointer
      use sad_main
      use ffs_flag
      use tmacro
      use ffs_seg
      use temw, only:calint,tmulbs,iaemit,beamplt
      use kradlib, only:tradke
      implicit none
      real*8 conv
      parameter (conv=3.d-16)
      type (iaemit) ,intent(in):: iae
      integer*8 iatrl,iacodl,iabmilz
      integer*4 ,intent(in):: k,idp
      integer*4 ,intent(out):: ke
      integer*4 i,l
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 rtaper, bmir(6,6),r
      logical*4 ,intent(out):: sol
      logical*4 ,intent(in):: plot,rt
      save iabmilz
      data iabmilz /0/
      sol=.true.
      do i=k+1,nlat-1
        if(idtypec(i) .eq. icSOL)then
          if(rlist(idvalc(i)+ky_BND_SOL) .ne. 0.d0)then
            ke=i
            go to 20
          endif
        endif
      enddo
      sol=.false.
      return
 20   do l=k,ke
        rtaper=1.d0
        if(radtaper)then
          if(rt)then
            rtaper=(2.d0+cod(6)+gettwiss(mfitddp,nextl(l)))*.5d0-dp0
          else
            rtaper=(1.d0-dp0+cod(6))
          endif
        endif
        call tsole1(trans,cod,beam,srot,l,rtaper,.true.,.false.)
c          write(*,*)'tsole-l ',l,sol,rtaper
c          write(*,'(1p6g15.7)')(trans(i,1:6),i=1,6),cod
        if(plot)then
          if(iae%iatr .ne. 0)then
            if(iae%iatr .gt. 0)then
              call tflocal(klist(iae%iatr+l+1))
              iatrl=ktfaddr(kxm2l(trans,6,6,6,.false.))
              klist(iae%iatr+l+1)=ktflist+ktfcopy1(iatrl)
            endif
            if(iae%iacod .gt. 0)then
              call tflocal(klist(iae%iacod+l+1))
              iacodl=ktfaddr(kxm2l(cod,0,6,1,.false.))
              klist(iae%iacod+l+1)=ktflist+ktfcopy1(iacodl)
            endif
          endif
          if(codplt)then
            call tsetetwiss(trans,cod,beam,0,l+1,idp)
            if(irad .gt. 6 .and. beamplt)then
              beamsize(:,l+1)=beam
            endif
          endif
          if(calint .and. iae%iabmi .ne. 0)then
            if(iabmilz .eq. 0)then
              bmir=0.d0
              iabmilz=ktfaddr(kxm2l(bmir,6,6,6,.false.))
            endif
            call tflocal(klist(iae%iabmi+l))
            klist(iae%iabmi+l)=ktflist+ktfcopy1(iabmilz)
          endif
        elseif(radtaper .and. radcod)then
          if(l .eq. 1)then
            r=1.d0
          else
            r=gammab(l)/gammab(l+1)
          endif
          twiss(l+1,idp,mfitddp)=cod(6)*r
        endif
      enddo
      return
      end

      subroutine tsole1(trans,cod,beam,srot,l,rtaper,enarad,qsol)
      use kyparam
      use tfstk
      use ffs_pointer
      use ffs_flag
      use ffs, only:mfitddp,gettwiss
      use tmacro
      use sad_main
      use ffs_seg
      use temw, only:tsetr0,tmulbs
      use kradlib, only:tradke
      use tspin
      implicit none
      integer*4 ,intent(in):: l
      integer*4 ld,lt,mfr,kb,irtc
      integer*8 lp
      type (sad_comp), pointer ::cmp
      type (sad_dlist), pointer :: lsegp
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: rtaper
      real*8 rr(3,3),al,theta,
     $     phi,phix,phiy,bzs,trans1(6,6),trans2(6,6),
     $     tfbzs,radlvl,bzs0,f1,ftable(4),ak1
      logical*4 ,intent(in):: enarad,qsol
      logical*4 dir,ent,coup,err,krad,seg
      real*8,save::dummy(256)=0.d0
      ld=idelc(l)
      lt=idtype(ld)
      lp=elatt%comp(l)
      call loc_comp(lp,cmp)
      seg=tcheckseg(cmp,lt,al,lsegp,irtc)
      if(irtc .ne. 0)then
        call tffserrorhandle(l,irtc)
        return
      endif
      bzs=tfbzs(l,kb)
      select case(lt)
      case (icDRFT)
        call tsetr0(trans,cod,bzs,0.d0)
        call tdrife(trans,cod,beam,srot,al,
     $       bzs,0.d0,0.d0,al,.true.,
     $       enarad .and. cmp%value(ky_RAD_DRFT) .eq. 0.d0,
     $       irad)
      case (icBEND)
        call tsetr0(trans,cod,bzs,0.d0)
        theta=cmp%value(ky_ROT_BEND)
     $       +cmp%value(ky_DROT_BEND)
        phi=cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND)
        phiy= phi*cos(theta)
        phix= phi*sin(theta)
        call tdrife(trans,cod,beam,srot,al,
     $       bzs,phiy,phix,al,.true.,
     $       enarad .and. cmp%value(ky_RAD_BEND) .eq. 0.d0,
     $       irad)
      case(icQUAD)
        dir=direlc(l) .gt. 0.d0
        if(dir)then
          mfr=nint(cmp%value(ky_FRMD_QUAD))
        else
          mfr=nint(cmp%value(ky_FRMD_QUAD))
          mfr=mfr*(11+mfr*(2*mfr-9))/2
        endif
        krad=enarad .and. cmp%value(ky_RAD_QUAD) .eq. 0.d0
     $       .and. al .ne. 0.d0
        ak1=cmp%value(ky_K1_QUAD)
        call tsetfringepe(cmp,icQUAD,ftable)
        call tquade(trans,cod,beam,srot,
     $       al,ak1,bzs,
     $       cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     $       cmp%value(ky_ROT_QUAD),
     1       krad,cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $       ftable(1),ftable(2),ftable(3),ftable(4),
     $       mfr,cmp%value(ky_EPS_QUAD),
     $       cmp%value(ky_KIN_QUAD) .eq. 0.d0,
     $       cmp%value(ky_CHRO_QUAD) .ne. 0.d0)
      case(icMULT)
        if(seg)then
          call tmulteseg(trans,cod,beam,srot,l,cmp,bzs,lsegp,
     $         enarad,rtaper)
        else
          call tmulte1(trans,cod,beam,srot,l,cmp,bzs,enarad,rtaper)
        endif
      case(icCAVI)
        call tmulte(trans,cod,beam,srot,l,al,
     $       dummy,
     $       bzs,
     $       0.d0,0.d0,0.d0,0.d0,0.d0,
     1       cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $       0.d0,0.d0,0.d0,
     $       cmp%value(ky_ROT_CAVI),
     $       0.d0,0.d0,enarad,
     $       cmp%value(ky_FRIN_CAVI) .eq. 0.d0,
     $       0.d0,0.d0,0.d0,0.d0,
     $       cmp%value(ky_FRMD_CAVI),0.d0,0.d0,
     $       .true.,
     $       cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI),
     $       cmp%value(ky_HARM_CAVI),
     $       cmp%value(ky_PHI_CAVI),cmp%value(ky_FREQ_CAVI),
     $       0.d0,1.d0,
     $       cmp%value(ky_APHI_CAVI) .ne. 0.d0,
     $       ld)
      case(icSOL)
        f1=cmp%value(ky_F1_SOL)
        krad=enarad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0
     $       .and. f1 .ne. 0.d0
        if(cmp%value(ky_BND_SOL) .ne. 0.d0)then
          ent=direlc(l) .gt. 0.d0 .and. l .eq. kb
     $         .or. direlc(l) .lt. 0.d0 .and. l .ne. kb
          bzs0=tfbzs(l-1,kb)
          if(krad)then
            if(ent)then
              call tsconv(trans1,cod,rr,lp,.true.)
              call tmultr5(trans,trans1,irad)
              call tmulbs(beam,trans1,.true.)
              if(calpol)then
                srot=matmul(rr,srot)
              endif
              call tsetr0(trans(:,1:6),cod,0.d0,0.d0)
              if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
                call tsfrie(trans1,cod,bzs)
                call tmultr5(trans,trans1,irad)
                call tmulbs(beam,trans1,.true.)
              endif
              call tradke(trans,cod,beam,srot,f1,0.d0,bzs*.5d0)
            else
              call tsetr0(trans(:,1:6),cod(1:6),bzs0*.5d0,0.d0)
              if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
                call tsfrie(trans1,cod,-bzs0)
                call tmultr5(trans,trans1,irad)
                call tmulbs(beam,trans1,.true.)
              endif
              call tradke(trans,cod,beam,srot,f1,0.d0,0.d0)
              call tsconv(trans1,cod,rr,lp,.false.)
              call tmultr5(trans,trans1,irad)
              call tmulbs(beam,trans1,.true.)
              if(calpol)then
                srot=matmul(rr,srot)
              endif
            endif
          else
            if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
              if(ent)then
                call tsconv(trans1,cod,rr,lp,.true.)
                call tsfrie(trans2,cod,bzs)
              else
                call tsfrie(trans1,cod,-bzs0)
                call tsconv(trans2,cod,rr,lp,.false.)
              endif
              call tmultr5(trans1,trans2,6)
            else
              call tsconv(trans1,cod,rr,lp,ent)
            endif
            call tmultr5(trans,trans1,irad)
            call tmulbs(beam ,trans1,.true.)
            if(calpol)then
              srot=matmul(rr,srot)
            endif
          endif
        else
          bzs0=tfbzs(l-1,kb)
          if(krad)then
            call tsetr0(trans(:,1:6),cod(1:6),bzs0*.5d0,0.d0)
c              call trades(trans,beam,cod,bzs0,bzs,f1,brhoz)
          endif
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            call tsfrie(trans1,cod,bzs-bzs0)
            call tmultr5(trans,trans1,irad)
            call tmulbs(beam ,trans1,.true.)
          endif
          if(krad)then
            call tradke(trans,cod,beam,srot,f1,0.d0,bzs*.5d0)
          endif
        endif
      case(icMAP)
        if(qsol)then
          call qemap(trans1,cod,l,coup,err)
          trans(:,1:6)=matmul(trans1,trans(:,1:6))
c          call tmultr(trans,trans1,6)
        else
          call temape(trans,cod,beam,l)
        endif
      end select
c      write(*,'(a,i5,1p6g15.7/16x,1p6g15.7)')'tsole1-end ',l,code,cod
      return
      end

      subroutine qsol(trans,cod,k,coup)
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in):: k
      real*8 ,intent(inout):: trans(4,5),cod(6),beam(42),srot(3,9)
      real*8 transe(6,12)
      logical*4 ,intent(out):: coup
      logical*4 radtaper0,calpol0
      radtaper0=radtaper
      radtaper=.false.
      calpol0=calpol
      calpol=.false.
      call tinitr(transe)
      call tsole1(transe,cod,beam,srot,k,1.d0,.false.,.true.)
      radtaper=radtaper0
      calpol=calpol0
      call qcopymatg(trans,transe,k)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
