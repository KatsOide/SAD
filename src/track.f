      subroutine track(latt,iparam)
      use tfstk
      use ffs
      use touschek_table, only: initialize_tampl
      use trackbypass, only: bypasstrack, lattuse
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      use tparastat
      implicit none
      type (sad_descriptor) kx
      integer*8 ikptbl,ig,ipz,ix,ixx,iy,iyy,iz,izz,ifz,imt,kzx,
     $     latt,iparam,lscal
      integer*4 irtc,l,isp1,
     $     nts,itfdownlevel,naz,ltpara,igetgl1
      character*20 title
      logical*4, save :: trackinit=.false.
      real*8 ol,trval,dt1,df,rgetgl1,dt0,phi(3)
c      write(*,*)'track-0'
      if(bypasstrack)then
        write(*,*)
     $       '??? FFS, EMIT, TRACK in GetMAIN are bypassed. ???'
        return
      endif
      open(0,file='/dev/null',status='UNKNOWN',err=10)
 10   call tsetupsig
      if(.not. trackinit)then
        trackinit=.true.
c        write(*,*)'track-0.0 ',klist( 1621700052)
        call tffsvinit
        iffssave=0
c        write(*,*)'track-0.1 ',klist( 1621700052)
        call ffs_init_flag
c        write(*,*)'track-0.2 ',klist( 1621700052)
        call csinit(0,1,'!',.false.)
c        write(*,*)'track-0.3 ',klist( 1621700052)
        call tfinitn
c        write(*,*)'track-0.4'
        call tfevals('CONVERGENCE=1E-9;ExponentOfResidual=2;'//
     $       'OffMomentumWeight=1;MatchingResidual=0;'//
     $       'NetResidual=0;StabilityLevel=0;'//
     $       'FFS$NumericalDerivative=False;'//
     $       'DP=0.01;DPM=0;XIX=0;XIY=0;TITLE="";CASE="";'//
     $       'NFAMP=3;'//
     $       'DP0:=LINE["DDP",1];(DP0=v_)^:=(LINE["DDP",1]=v);'//
     $       'Protect[DP0];'//
     $       'System$Names=Select[Names["*"],'//
     $       'ToUpperCase[#[1]]==#[1]&]',
     $       kx,irtc)
        initmessage=0
        ifibzl=0
        ifgamm=0
        tparaed=.false.
      else
        levele=1
      endif
c      isp1=isp-1
c      call tfmemcheck(isp1,kx,irtc)
c      l=itfdownlevel()
      levele=levele+1
      novfl=0
      call cputime(dt0,irtc)
      pltfl=8
      ilattp=latt
      call loc_el(ilattp,elatt)
      lattuse=latt
      lscal =klist(iparam+1)
      np0   =ilist(1,lscal+1)
      nturn =ilist(1,lscal+2)
      ltpara=ilist(2,elatt%aux)
      charge=rlist(lscal+4)
      call rsetgl1('NP',dble(np0))
      call rsetgl1('TURNS',dble(nturn))
      call rsetgl1('CHARGE',charge)
      amass =rgetgl1('MASS')
      pgev  =rgetgl1('MOMENTUM')
      pbunch=rgetgl1('PBUNCH')
      anbunch=rgetgl1('NBUNCH')
      coumin=rgetgl1('MINCOUP')
      emidiv=rgetgl1('EMITDIV')
      emidib=rgetgl1('EMITDIVB')
      emidiq=rgetgl1('EMITDIVQ')
      emidis=rgetgl1('EMITDIVS')
      fridiv=rgetgl1('FRINGDIV')
      alost =rgetgl1('LOSSAMPL')
      zlost =rgetgl1('LOSSDZ')
      trf0  =rgetgl1('DTSYNCH')
      vcalpha=rgetgl1('EFFVCRATIO')
      nlat  =elatt%nlat1-1
c      write(*,*)'track (np0,nturn,nlat) =',np0,nturn,nlat
      df    =rgetgl1('FSHIFT')
      isynch=igetgl1('$RFSW$'  )
      intra =igetgl1('$INTRA$' ) .ne. 0
      calpol=igetgl1('$POL$'   ) .ne. 0
      rad   =igetgl1('$RAD$'   ) .ne. 0
      calcod=igetgl1('$COD$'   ) .ne. 0
      trpt  =igetgl1('$TRPT$'  ) .ne. 0
      radcod=igetgl1('$RADCOD$') .ne. 0
      emiout=igetgl1('$EMIOUT$') .ne. 0
      dapert=igetgl1('$DAPERT$') .ne. 0
      rfluct=igetgl1('$FLUC$'  ) .ne. 0
      k64   =igetgl1('$K64$'   ) .ne. 0
      fourie=igetgl1('$FOURIE$') .ne. 0
      smearp=igetgl1('$SMEAR$' ) .ne. 0
      geocal=igetgl1('$GEOCAL$' ) .ne. 0
      calc6d=igetgl1('$CALC6D$') .ne. 0
      intres=igetgl1('$INTRES$') .ne. 0
      halfres=igetgl1('$HALFRES$') .ne. 0
      sumres=igetgl1('$SUMRES$') .ne. 0
      diffres=igetgl1('$DIFFRES$') .ne. 0
      photons=igetgl1('$PHOTONS$' ) .ne. 0
      nparallel=max(1,int(rgetgl1('NPARA')))
      calc6d=.false.
      radlight=.false.
      ffsprmpt=.false.
      rfsw  =isynch .ne. 0
      call tsetgcut
      call tphyzp
      call tsetdvfs
      ol=rlist(elatt%aux+1)
      if(nturn .lt. 0 .and. np0 .eq. 0)then
        if(ol .le. 0.d0)then
          omega0=0.d0
        else
          write(*,*)'Design orbit length =',ol
        endif
      else
        call tfinitgeo
      endif
      jitter=.true.
      trgauss=.true.
      gauss=.false.
      bunchsta=.false.
      spac=.false.
      wspac=.false.
      selfcod=.false.
      orbitcal=.true.
      dp0   =0.d0
      call initialize_tampl()
      call tclrpara(elatt,nlat-1)
      call tclrfpe
      write(*,'(a)')
     1' RFSW RADCOD RAD  FLUC  INTRA'//
     1' POL   COD  DAPER EMIOU K64   FOURI SMEAR'
      write(*,9001)rfsw,radcod,rad,rfluct,
     1             intra,calpol,calcod,dapert,emiout,k64,fourie,
     1             smearp
9001  format(1x,12(L3,3X))
      if(nturn .eq. 0)then
        write(*,*)'Use EMIT_TANCE or Emittance[] within FFS.'
        go to 8001
      endif
      if(nturn .lt. 0)then
        if(np0 .eq. 0)then
          trpt=.false.
          rfsw=.true.
          call tffs
          title='FFS'
          go to 8001
        elseif(np0 .eq. -1)then
          title='UNDEFINED'
          go to 8001
        elseif(nturn .eq. -1)then
          nturn=1
          trpt=.true.
          rfsw=.true.
        else
        endif
      else
        trpt=.false.
      endif
      nturn=abs(nturn)
      ikptbl=ktaloc(np0*3)
      call tkptblini(ilist(1,ikptbl))
      ig=ktaloc(np0)
      ipz=ktaloc(np0)
      if(.not. dapert .or. trpt)then
        if(trpt)then
          ix=1
          ixx=1
          iy=1
          iyy=1
          iz=1
          izz=1
          ifz=0
          nts=1
        else
          ix=ktaloc(np0)
          ixx=ktaloc(np0)
          iy=ktaloc(np0)
          iyy=ktaloc(np0)
          iz=ktaloc(np0)
          izz=ktaloc(np0)
          nts=2**int(log(dble(nturn+1))/log(2.d0))
          if(nts .lt. nturn+1)then
            nts=nts*2
          endif
          if(nts .ge. 2*(nturn+1))then
            nts=nts/2
          endif
          if(smearp)then
            naz=np0*nts*6
            ifz=ktaloc(naz)
            if(ifz .gt. 0)then
              call tclr(rlist(ifz),naz)
            else
              ifz=0
            endif
          else
            ifz=0
          endif
        endif
        call tracka(rlist(ilattp+1),rlist(ikptbl),
     1        rlist(ilist(2,iparam+16)),rlist(ilist(2,iparam+17)),
     1        rlist(ilist(2,iparam+18)),rlist(ilist(2,iparam+19)),
     1        rlist(ilist(2,iparam+20)),rlist(ilist(2,iparam+21)),
     1        rlist(ig),rlist(ipz),
     1        rlist(ix),rlist(ixx),rlist(iy),rlist(iyy),
     1        rlist(iz),rlist(izz),
     1        rlist(ifz),max(1,nts),ifz .gt. 0)
        if(.not. trpt)then
          if(ifz .gt. 0)then
            call tfree(ifz)
          endif
          call tfree(izz)
          call tfree(iz)
          call tfree(iyy)
          call tfree(iy)
          call tfree(ixx)
          call tfree(ix)
        endif
      else
        imt=ktaloc(np0/2+1)
        kzx=ktaloc(np0)
        trval=0.d0
        rlist(ipz)=0.d0
        phi=0.d0
        call trackd(ilist(1,ilattp+1),ilist(1,ikptbl),
     1        rlist(ilist(2,iparam+16)),rlist(ilist(2,iparam+17)),
     1        rlist(ilist(2,iparam+18)),rlist(ilist(2,iparam+19)),
     1        rlist(ilist(2,iparam+20)),rlist(ilist(2,iparam+21)),
     1        rlist(ig),rlist(ipz),
     1        ilist(1,kzx),ilist(1,imt),trval,phi,0.d0,0.d0,3,1,outfl)
        call tfree(kzx)
        call tfree(imt)
      endif
      call tsptrm
      call tfree(ig)
      call tfree(ikptbl)
      title='Tracking'
 8001 call rsetgl1('LOSSAMPL',alost )
      call rsetgl1('LOSSDZ',zlost )
      call isetgl1('$RFSW$',isynch)
      call isetgl1('$INTRA$',intra )
      call isetgl1('$POL$',calpol)
      call isetgl1('$RAD$',rad   )
      call isetgl1('$COD$',calcod)
      call isetgl1('$RADCOD$',radcod)
      call isetgl1('$EMIOUT$',emiout)
      call isetgl1('$DAPERT$',dapert)
      call isetgl1('$FLUC$',rfluct)
      call isetgl1('$K64$',k64)
      call isetgl1('$FOURIE$',fourie)
      call isetgl1('$SMEAR$',smearp)
      call isetgl1('$GEOCAL$',geocal)
      call isetgl1('$INTRES$',intres)
      call isetgl1('$HALFRES$',halfres)
      call isetgl1('$SUMRES$',sumres)
      call isetgl1('$DIFFRES$',diffres)
      call isetgl1('$PHOTONS$',photons)
      nlat  =elatt%nlat1-1
      call tclrpara(elatt,nlat-1)
      call cputime(dt1,irtc)
      write(*,'(1X,2A,F10.3,A)')
     1     title,' end:  CPU time =',(dt1-dt0)*1.d-6,' sec'
      l=itfdownlevel()
      levele=levele+1
      isp1=isp-1
      l=itfdownlevel()
      close(0)
      return
      end

      integer*8 function itfilattp()
      use tfstk
      use tmacro
      implicit none
      itfilattp=ilattp
      return
      end

      real*8 function rgetgl1(vname)
      implicit none
      integer*4 ia
      character*(*) vname
      real*8 rgetgl
      ia=0
      rgetgl1=rgetgl(vname,ia)
      return
      end

      integer*4 function igetgl1(vname)
      implicit none
      integer*4 ia,igetgl
      character*(*) vname
      ia=0
      igetgl1=igetgl(vname,ia)
      return
      end

      subroutine rsetgl1(vname,val)
      implicit none
      integer*4 ia
      character*(*) vname
      real*8 val
      ia=0
      call rsetgl(vname,val,ia)
      return
      end

      subroutine isetgl1(vname,ival)
      implicit none
      integer*4 ia,ival
      character*(*) vname
      ia=0
      call isetgl(vname,ival,ia)
      return
      end

      subroutine twigp
      return
      end

      subroutine twig
      return
      end

      subroutine liemap
      return
      end

      subroutine tdimad
      return
      end

      subroutine aaalie
      return
      end
