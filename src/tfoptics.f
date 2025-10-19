      subroutine tfoptics(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      use tfcsi, only:icslfno
      use ffs_fit, only:ffs_stat
      use cellm, only:qcell1
      implicit none
      type (ffs_bound) fbound
      type (ffs_stat) optstat
      type (sad_dlist), pointer :: klx1,klx
      type (sad_rlist), pointer :: klx2,klx3,klxi
      integer*8 kx,kax,kax1,kax2,kax3,kaxi,kaini,itoff
      integer*4 isp1,narg,irtc,idim,k,i,itfloc,itfmessage,lout
      logical*4 cell0,geo,calc6d0
      narg=isp-isp1
      if(narg /= 5)then
        if(narg /= 6)then
          irtc=itfmessage(9,'General::narg','"5 or 6"')
          return
        endif
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real(GEOCAL) for #6"')
          return
        endif
        geo=ktastk(isp) /= 0
      else
        geo=.false.
      endif
      if(ktfnonrealq(ktastk(isp1+5)))then
        irtc=itfmessage(9,'General::wrongtype','"Real(IDIM) for #5"')
        return
      endif
      idim=int(rtastk(isp1+5))
      if(idim /= 2 .and. idim /= 3)then
        irtc=itfmessage(9,'General::wrongnum','"IDIM must be 2 or 3"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+4)))then
        irtc=itfmessage(9,'General::wrongtype','"Real(CELL) for #4"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+3)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List(inicond) for #3"')
        return
      endif
      kaini=ktfaddr(ktastk(isp1+3))
      if(ilist(2,kaini-1) /= ntwissfun)then
        irtc=itfmessage(9,'General::wrongleng','"#3","28"')
        return
      endif
      if(ktfnonreallistq(kaini))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of Reals for #3"')
        return
      endif
      fbound%lb=itfloc(ktastk(isp1+1),irtc)
      if(irtc /= 0)then
        return
      endif
      fbound%le=itfloc(ktastk(isp1+2),irtc)
      if(irtc /= 0)then
        return
      endif
      cell0=cell
      cell=rtastk(isp1+4) /= 0.d0
      do k=1,ntwissfun
        itoff=(2*ndim+1)*nlat*(k-1)+ndim*nlat+iftwis+fbound%lb-1
        rlist(itoff)=rlist(kaini+k)
      enddo
      if(geo)then
        call tfgeo(.true.)
      endif
      lout=icslfno()
      fbound%fb=0.d0
      fbound%fe=0.d0
      calc6d0=calc6d
      if(idim .eq. 2)then
        calc6d=.false.
      else
        calc6d=.true.
      endif
      call qcell1(fbound,0,optstat,.false.,lout)
      calc6d=calc6d0
      cell=cell0
      kax=ktadaloc(-1,3,klx)
      kax2=ktraaloc(0,3,klx2)
      klx%dbody(2)%k=ktflist+kax2
      if(optstat%stabx)then
        klx2%rbody(1)=1.d0
      endif
      if(optstat%staby)then
        klx2%rbody(2)=1.d0
      endif
      if(optstat%stabz)then
        klx2%rbody(3)=1.d0
      endif
      kax3=ktraaloc(0,3,klx3)
      klx%dbody(3)%k=ktflist+kax3
      klx3%rbody(1)=optstat%tracex
      klx3%rbody(2)=optstat%tracey
      klx3%rbody(3)=optstat%tracez
      kax1=ktadaloc(0,fbound%le-fbound%lb+1,klx1)
      klx%dbody(1)%k=ktflist+kax1
      do i=fbound%lb,fbound%le
        kaxi=ktraaloc(0,ntwissfun,klxi)
        klx1%dbody(i-fbound%lb+1)%k=ktflist+kaxi
        do k=1,ntwissfun
          itoff=(2*ndim+1)*nlat*(k-1)+ndim*nlat+iftwis+i-1
          klxi%rbody(k)=rlist(itoff)
        enddo
      enddo
      kx=ktflist+kax
      return
      end
