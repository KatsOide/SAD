      subroutine tfemit(isp1,kx,irtc)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use iso_c_binding
      use tfcsi, only:icslfno
      use temw, only:tfinitemip,nparams,tfetwiss,iptwiss,tfinibeam,
     $     iaemit,iaez
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: kl1,kl2,klr
      type (ffs_bound) fbound
      type (iaemit) iae
      integer*8 kaparam
      integer*4 isp1,irtc,narg,mode,itgetfpe,itfmessage,lno,
     $     is,ie,nel
      real*8 param(nparams),trans(6,12),cod(6),beam(42),btr(441),sx,
     $     srot(3,9)
      real*8 dummy(6,6),ris(6,6)
      logical*4 stab,rt,bi
      call tfinitemip
      narg=isp-isp1
      codin=0.d0
      beamin=0.d0
      if(narg .gt. 4 .or. narg .lt. 2)then
        irtc=itfmessage(9,'General::narg','"2, 3, or 4"')
        return
      endif
      is=1
      ie=nlat
      if(narg .eq. 4)then
        if(tfnonreallistq(dtastk(isp),klr))then
          irtc=itfmessage(9,'General::wongval','"{begin, end} for #4"')
          return
        endif
        is=int(klr%rbody(1))
        ie=merge(nlat-int(abs(klr%rbody(2)))+1,int(klr%rbody(2)),
     $       klr%rbody(2) .lt. 0.d0)
c        write(*,'(a,1p2g15.7,3i5)')'tfemit ',klr%rbody(1:2),is,ie,nlat
      endif
      nel=ie-is+1
      bi=.false.
      if(narg .ge. 3)then
        if(tflistq(ktastk(isp1+3),kl))then
          if(tfreallistq(kl%dbody(1),kl1))then
            if(kl1%nl .ne. 6)then
              irtc=itfmessage(9,'General::wrongval',
     $      '"{x, px, y, py, z, dp} for InitialOrbit"')
              return
            endif
            codin=kl1%rbody(1:6)
          elseif(nel .ne. nlat)then
            codin=twiss(is,0,mfitdx:mfitddp)
          endif
          if(tfreallistq(kl%dbody(2),kl2))then
            beamin=kl2%rbody(1:21)
            bi=.true.
          endif
        endif
      endif
      if(.not. bi .and. trpt)then
        beamin=tfinibeam(is)
      endif
      if(ktfnonrealq(ktastk(isp1+1),mode))then
        irtc=itfmessage(9,'General::wrongtype','"Real for #1"')
        return
      endif
      if(mode .lt. -1 .or. mode .gt. 3)then
        irtc=itfmessage(9,'General::wrongnum','"-1, 0, 1, 2, or, 3"')
        return
      endif
      lno=int(rtastk(isp1+2))
      if(lno .eq. -1)then
        lno=icslfno()
      endif
      irtc=0
      iae=iaez
      if(mode .ge. 1)then
        iae%iamat=ktadalocnull(-1,6)
        if(mode .ge. 2)then
          iae%iacod=ktadalocnull(-1,nel)
          if(mode .eq. 3)then
            iae%iatr=ktadalocnull(-1,nel)
            if(intra)then
              iae%iabmi=ktadalocnull(-1,nel)
            endif
          endif
        endif
      endif
      call tfsetparam
      if(ifsize .eq. 0 .and. codplt)then
        ifsize=ktaloc(nel*21)
        call c_f_pointer(c_loc(rlist(ifsize)),beamsize,[21,nel])
        modesize=0
      endif
c      write(*,*)'tfemit-4 ',codplt,ifsize,nel
      if(nel .eq. nlat)then
        cod=codin
        call temit(trans,cod,beam,btr,
     $       .not. trpt .and. mode .ge. 0,iae,codplt,param,stab,lno)
      else
        call tffsbound1(is,ie-1,fbound)
        cod=codin
        rt=radcod .and. radtaper
        param=dnotanumber
        irad=12
        beam(1:21)=beamin
        beam(22:42)=0.d0
        call tinitr12(trans)
        call tturneg(trans,cod,beam,srot,fbound,
     $     iae,.true.,rt,.false.)
        call setparams(param,cod)
        call rotri(is,ris)
        dummy=dnotanumber
        call setiamat(iae%iamat,ris,codin,beam,dummy,trans)
        param(iptwiss:iptwiss+ntwissfun-1)=tfetwiss(ris,codin,
     $       twiss(is,0,mfitdetr) .lt. 1.d0)
      endif
      kx=kxadaloc(-1,merge(6,2+max(0,mode),mode .eq. 3 .and. intra),
     $     klx)
      kaparam=ktfaddr(kxm2l(param,0,nparams,1,.false.))
      if(itgetfpe() .gt. 0)then
        stab=.false.
        call tclrfpe
      endif
      sx=merge(1.d0,0.d0,stab)
c      write(*,*)mode,iax,iabmi,iamat,iaparam,nparams
      klx%rbody(1)=sx
      klx%dbody(2)%k=ktflist+ktfcopy1(kaparam)
      if(mode .ge. 1)then
        klx%dbody(3)%k=ktflist+ktfcopy1(iae%iamat)
        if(mode .ge. 2)then
          klx%dbody(4)%k=ktflist+ktfcopy1(iae%iacod)
          if(mode .eq. 3)then
            klx%dbody(5)%k=ktflist+ktfcopy1(iae%iatr)
            if(intra)then
              klx%dbody(6)%k=ktflist+ktfcopy1(iae%iabmi)
            endif
          endif
        endif
      endif
      return
      end

      subroutine tfemits(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      use temw, only:nparams,tfinitemip
      implicit none
      type (sad_descriptor) kx
      integer*8 kparams
      integer*4 isp1,irtc,mphi2,i,itfmessage
      real*8 arg(4),emxe,emye,rese
      call tfinitemip
      if(isp .ne. isp1+4)then
        go to 9001
      endif
      do i=1,4
        if(ktfnonrealq(ktastk(isp1+i)))then
          go to 9001
        endif
        arg(i)=rtastk(isp1+i)
      enddo
      mphi2=int(max(1.d0,min(32.d0,arg(4))))
      call tfgeo(.true.)
      kparams=ktaloc(nparams)
      codin=0.d0
      if(ifsize .eq. 0 .and. codplt)then
        ifsize=ktaloc(nlat*21)
c        ilist(2,iwakepold+6)=ifsize
      endif
      kx%k=ktfoper+mtfnull
      call temits(
     $     mphi2,
     $     arg(1)*pi2,arg(2)*pi2,arg(3)*pi2,
     $     emxe,emye,rese,rlist(kparams),
     $     0,kx,irtc)
      call tfree(kparams)
      if(.not. codplt)then
        call tfgeo(.true.)
      endif
      return
 9001 irtc=itfmessage(9,'General::wrongtype',
     $     '"[{nusstart,nusstop,nusstep},mphi]"')
      return
      end

      subroutine rotri(is,ris)
      use ffs_pointer
      use tffitcode
      use temw,only:r,ri,tinv6
      implicit none
      integer*4 ,intent(in)::is
      real*8 , intent(out)::ris(6,6)
      real*8 mu(3),c(3),s(3),trans(6,6)
      if(is .ne. 1)then
        call tftmat6(trans,1,is)
        ris=tinv6(matmul(trans,r))
c        call tmultr(rs,trans,6)
c        call tinv6(rs,ris)
        mu=twiss(is,0,(/mfitnx,mfitny,mfitnz/))
        c=cos(mu)
        s=sin(mu)
        trans=0.d0
        trans(1,1)=c(1)
        trans(1,2)=s(1)
        trans(2,1)=-s(1)
        trans(2,2)=c(1)
        trans(3,3)=c(2)
        trans(3,4)=s(2)
        trans(4,3)=-s(2)
        trans(4,4)=c(2)
        trans(5,5)=c(3)
        trans(5,6)=s(3)
        trans(6,5)=-s(3)
        trans(6,6)=c(3)
        ris=matmul(trans,ris)
c        call tmultr(ris,trans,6)
      else
        ris=ri
      endif
      return
      end
