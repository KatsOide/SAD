      subroutine track(latt,iparam)
      use tfstk
      use ffs
      use touschek_table, only: initialize_tampl
      use trackbypass, only: bypasstrack, lattuse
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      use tparastat
c      use tfcsi, only:ipoint,lrecl,lfni
      implicit none
      type (sad_descriptor) kx
      integer*8 ikptbl,ig,ipz,ix,ixx,iy,iyy,iz,izz,ifz,
     $     latt,iparam,lscal
      integer*4 irtc,l,isp1,
     $     nts,itfdownlevel,naz,ltpara,IgetGL
      character*20 title
      logical*4, save :: trackinit=.false.
      real*8 ol,dt1,df,rgetgl1,dt0
      if(bypasstrack)then
        write(*,*)
     $       '??? FFS, EMIT, TRACK in GetMAIN are bypassed. ???'
        return
      endif
      open(0,file='/dev/null',status='UNKNOWN',err=10)
 10   call tsetupsig
      if(.not. trackinit)then
        trackinit=.true.
        call tffsvinit
        iffssave=0
        call ffs_init_flag
        convcase=.true.
        call tfinitn
        call tfinittws
        call tfevals('CONVERGENCE=1E-9;ExponentOfResidual=2;'//
     $       'OffMomentumWeight=1;MatchingResidual='//
     $       'NetResidual=StabilityLevel=0;'//
     $       'FFS$NumericalDerivative=False;'//
     $       'DP=0.01;DPM=XIX=XIY=0;TITLE=CASE="";NFAMP=4;'//
     $       '(DP0=v_)^:=(LINE["DDP",1]=v);'//
     $       'System$Names=Select[Names["*"],'//
     $       'ToUpperCase[#[1]]==#[1]&];Protect[DP0];',
     $       kx,irtc)
        initmessage=0
        ifibzl=0
        ifgamm=0
        tparaed=.false.
      else
        levele=1
      endif
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
      nlat  =elatt%nlat0+1
      df    =rgetgl1('FSHIFT')
      isynch=IgetGL('$RFSW$'  )
      intra =IgetGL('$INTRA$' ) .ne. 0
      calpol=IgetGL('$POL$'   ) .ne. 0
      rad   =IgetGL('$RAD$'   ) .ne. 0
      calcod=IgetGL('$COD$'   ) .ne. 0
      trpt  =IgetGL('$TRPT$'  ) .ne. 0
      radcod=IgetGL('$RADCOD$') .ne. 0
      radpol=IgetGL('$RADPOL$') .ne. 0
      emiout=IgetGL('$EMIOUT$') .ne. 0
      dapert=IgetGL('$DAPERT$') .ne. 0
      rfluct=IgetGL('$FLUC$'  ) .ne. 0
      k64   =IgetGL('$K64$'   ) .ne. 0
      fourie=IgetGL('$FOURIE$') .ne. 0
      smearp=IgetGL('$SMEAR$' ) .ne. 0
      geocal=IgetGL('$GEOCAL$' ) .ne. 0
      calc6d=IgetGL('$CALC6D$') .ne. 0
      intres=IgetGL('$INTRES$') .ne. 0
      halfres=IgetGL('$HALFRES$') .ne. 0
      sumres=IgetGL('$SUMRES$') .ne. 0
      diffres=IgetGL('$DIFFRES$') .ne. 0
      photons=IgetGL('$PHOTONS$' ) .ne. 0
      nparallel=max(1,int(rgetgl1('NPARA')))
      keepexp=.true.
      calexp=.true.
      calc6d=.false.
      radlight=.false.
      ffsprmpt=.false.
      rfsw  =isynch .ne. 0
      suspend=.true.
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
      spac=.false.
      wspac=.false.
      selfcod=.false.
      orbitcal=.true.
      calopt=.true.
      dp0   =0.d0
      call initialize_tampl()
      call tclrpara
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
        write(*,*)'Obsolete -- use DynamicApertureSurvey in FFS.'
c$$$        imt=ktaloc(np0/2+1)
c$$$        kzx=ktaloc(np0)
c$$$        trval=0.d0
c$$$        rlist(ipz)=0.d0
c$$$        phi=0.d0
c$$$        call trackd(ilist(1,ikptbl),
c$$$     1        rlist(ilist(2,iparam+16)),rlist(ilist(2,iparam+17)),
c$$$     1        rlist(ilist(2,iparam+18)),rlist(ilist(2,iparam+19)),
c$$$     1        rlist(ilist(2,iparam+20)),rlist(ilist(2,iparam+21)),
c$$$     1        rlist(ig),rlist(ipz),
c$$$     1        ilist(1,kzx),ilist(1,imt),trval,phi,0.d0,0.d0,3,1,outfl)
c$$$        call tfree(kzx)
c$$$        call tfree(imt)
      endif
      call tsptrm
      call tfree(ig)
      call tfree(ikptbl)
      title='Tracking'
 8001 call rsetgl1('LOSSAMPL',alost )
      call rsetgl1('LOSSDZ',zlost )
      call isetgll('$RFSW$',rfsw)
      call isetgll('$INTRA$',intra )
      call isetgll('$POL$',calpol)
      call isetgll('$RAD$',rad   )
      call isetgll('$COD$',calcod)
      call isetgll('$RADCOD$',radcod)
      call isetgll('$RADPOL$',radpol)
      call isetgll('$EMIOUT$',emiout)
      call isetgll('$DAPERT$',dapert)
      call isetgll('$FLUC$',rfluct)
      call isetgll('$K64$',k64)
      call isetgll('$FOURIE$',fourie)
      call isetgll('$SMEAR$',smearp)
      call isetgll('$GEOCAL$',geocal)
      call isetgll('$INTRES$',intres)
      call isetgll('$HALFRES$',halfres)
      call isetgll('$SUMRES$',sumres)
      call isetgll('$DIFFRES$',diffres)
      call isetgll('$PHOTONS$',photons)
      nlat  =elatt%nlat0+1
      call tclrpara
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

      subroutine isetgll(vname,lval)
      implicit none
      logical*4 , intent(in)::lval
      integer*4 ia,ival
      character*(*) vname
      ia=0
      ival=lval
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
