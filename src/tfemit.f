      subroutine tfemit(isp1,kx,irtc)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use iso_c_binding
      use tfcsi, only:icslfno
      use temw, only:tfinitemip,nparams
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: kl1,kl2
      integer*8 iatr,iacod,iamat,kaparam,iabmi
      integer*4 isp1,irtc,narg,mode,itgetfpe,itfmessage,lno
      real*8 param(nparams),trans(6,12),cod(6),beam(42),btr(441),sx
      logical*4 stab
      call tfinitemip
      narg=isp-isp1
      codin=0.d0
      beamin=0.d0
      if(narg .eq. 3)then
        if(tflistq(ktastk(isp),kl))then
          if(tfreallistq(kl%dbody(1),kl1))then
            if(kl1%nl .ne. 6)then
              irtc=itfmessage(9,'General::wrongval',
     $      '"{x, px, y, py, z, dp} for InitialOrbit"')
              return
            endif
            codin=kl1%rbody(1:6)
          endif
          if(tfreallistq(kl%dbody(2),kl2))then
            beamin=kl2%rbody(1:21)
          endif
        endif
      elseif(narg .ne. 2)then
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
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
      iatr=-1
      iacod=-1
      iamat=-1
      iabmi=0
      if(mode .ge. 1)then
        iamat=ktadalocnull(-1,6)
        if(mode .ge. 2)then
          iacod=ktadalocnull(-1,nlat)
          if(mode .eq. 3)then
            iatr=ktadalocnull(-1,nlat)
            if(intra)then
              iabmi=ktadalocnull(-1,nlat)
            endif
          endif
        endif
      endif
      call tfsetparam
      if(ifsize .eq. 0 .and. codplt)then
        ifsize=ktaloc(nlat*21)
        call c_f_pointer(c_loc(rlist(ifsize)),beamsize,[21,nlat])
        updatesize=.false.
c        ilist(2,iwakepold+6)=ifsize
      endif
      call temit(trans,cod,beam,btr,
     $     mode .ge. 0,iatr,iacod,iabmi,iamat,
     $     .true.,param,stab,lno)
      if(mode .eq. 3 .and. intra)then
        kx=kxadaloc(-1,6,klx)
      else
        kx=kxadaloc(-1,2+max(0,mode),klx)
      endif
      kaparam=ktfaddr(kxm2l(param,0,nparams,1,.false.))
      if(itgetfpe() .gt. 0)then
        stab=.false.
        call tclrfpe
      endif
      if(stab)then
        sx=1.d0
      else
        sx=0.d0
      endif
c      write(*,*)mode,iax,iabmi,iamat,iaparam,nparams
      klx%rbody(1)=sx
      klx%dbody(2)%k=ktflist+ktfcopy1(kaparam)
      if(mode .ge. 1)then
        klx%dbody(3)%k=ktflist+ktfcopy1(iamat)
        if(mode .ge. 2)then
          klx%dbody(4)%k=ktflist+ktfcopy1(iacod)
          if(mode .eq. 3)then
            klx%dbody(5)%k=ktflist+ktfcopy1(iatr)
            if(intra)then
              klx%dbody(6)%k=ktflist+ktfcopy1(iabmi)
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
