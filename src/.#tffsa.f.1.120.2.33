      subroutine tffsa(lfnb,kffs,irtcffs)
      use tfstk
      use ffs
      use trackbypass
      use ffslocal
      use tffitcode
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ndimmax,nflaga,nflagn,maxrpt,maxlfn,hsrchz,nwakep
      integer*8 kffs,k,kx,ktaloc,kwakep,kwakeelm,itwisso,
     $     ifvalvar2,iparams,kax,iutwiss
      integer*4 kk,i,lfnb,ia,iflevel,j,nqcol1,ielm,igelm,k1,
     $     ii,irtc0,it,itemon,itmon,itestr,itstr,itt,
     $     itt1,itt2,itt3,itt4,itt5,iuse,l,itfuplevel,icslrecl,
     $     levelr,lfnl0,lpw,meas0,mfpnta,igetgl1,lenw,
     $     mphi2,newcor,next,nextt,nfam,nfam1,nfcol,nfp,nfr,nmon,
     $     nqcol,nster,nrpt1,ns,nut,icsmrk,itfpeeko,
     $     itfgetrecl,nl
      real*8 rmax,amus0,amus1,amusstep,apert,axi,ayi,ctime1,
     $     dpm2,dpxi,dpyi,em,emxe,emye,epxi,epyi,pspan,r,r2i,r3i,
     $     trval,rese,v,wa,wd,wl,xa,ya,xxa,xya,yya,getva,rgetgl1,
     $     wp,getvad,tgetgcut
      parameter (rmax=1.d35,ndimmax=500)
      parameter (nflaga=44,nflagn=nflag-nflaga)
      parameter (maxrpt=32,maxlfn=128)
      character*255 word,wordp,title,case,tfgetstrv,tfgetstrs,tfgetstr
      character*8 nlist(mfit1)
      character*(MAXPNAME) ename
      character*(MAXPNAME+8) name
      character*16 autofg
      character*20 str
      integer*4 irtcffs,irtc,nc,nfc0
      integer*4 lfnstk(maxlfn),lfret(maxlfn),
     $     lfrecl(maxlfn),lflinep(maxlfn),nrpt(maxrpt),
     $     irptp(maxrpt),iuid(-ndimmax:ndimmax),
     $     jfam(-ndimmax:ndimmax),kfam(-ndimmax:ndimmax)
      real*8 scale(mfit1),dp(-ndimmax:ndimmax),
     $     tracex(-ndimmax:ndimmax),tracey(-ndimmax:ndimmax),
     $     dfam(4,-ndimmax:ndimmax),residual(-ndimmax:ndimmax),
     $     df(maxcond),chi0(3),trdtbl(3,6),
     $     uini(28,-ndimmax:ndimmax)
      logical*4 hstab(-ndimmax:ndimmax),vstab(-ndimmax:ndimmax),
     $     geomet,err,new,cmd,open98,abbrev,ftest,
     $     fitflg,frefix,exist,init,trpt0,expnd,chguse,visit,
     $     byeall,expndc,tfvcomp,tffsinitialcond,inicond,wake,
     $     geocal0,busy
      logical*4 lfopen(maxlfn)
      save lfnstk,lfret,lfrecl,lfopen,open98,exist,init,scale,
     $     nmon,nster
      save busy
      data busy /.false./
      common /tt/ itt1,itt2,itt3,itt4,itt5
      data nlist /
     $     'AX      ','BX      ','NX      ','AY      ',
     1     'BY      ','NY      ','EX      ','EPX     ',
     1     'EY      ','EPY     ','R1      ','R2      ',
     1     'R3      ','R4      ','DETR    ',
     $     'DX      ','DPX     ',
     1     'DY      ','DPY     ','DZ      ','DDP     ',
     1     'PEX     ','PEPX    ','PEY     ','PEPY    ',
     $     'TRX     ','TRY     ','LENG    ','GX      ',
     $     'GY      ','GZ      ','CHI1    ','CHI2    ',
     $     'CHI3    ','DEX     ','DEPX    ','DEY     ',
     $     'DEPY    ','DDX     ','DDPX    ','DDY     ',
     $     'DDPY    ','PDEX    ','PDEPX   ','PDEY    ',
     $     'PDEPY   '/
      itwisso(kk,i,j)=iftwis+kk+nlat*(i+ndim+(j-1)*ndima)-1
      flv%mcommon=int((sizeof(flv)+7)/8)
c      write(*,*)'tffsa ',flv%mcommon
      kffs=ktfoper+mtfnull
      irtcffs=0
      l=itfuplevel()
      chguse=.false.
c     begin initialize for preventing compiler warning
      levelr=0
c     end   initialize for preventing compiler warning
 101  if(lfnb .le. 1 .or. chguse)then
        call tffsalloc(ilist(1,ilattp+1),ilist(1,ilattp))
        if(.not. chguse)then
          call cputime(flv%ctime0,irtc0)
          flv%ctime2=flv%ctime0
        endif
        flv%iut=0
        flv%nvar=0
        flv%ntouch=0
        call tclr(trdtbl,18)
        itt1=0
        flv%itmax=40
c        bzero=1.0d0
        id1=1
        id2=nlat
        iorgr=1
        geo0=0.d0
        geo0(1,1)=1.d0
        geo0(2,2)=1.d0
        geo0(3,3)=1.d0
        chi0(1)=0.d0
        chi0(2)=0.d0
        chi0(3)=0.d0
        if(geocal .or. chguse)then
          geocal0=geocal
          geocal=.true.
          call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $         rlist(ifgamm),.true.)
          geocal=geocal0
        endif
        if(lfnb .le. 0)then
          go to 8900
        endif
        call tffsinitparam(ilist(1,ilattp+1))
c     kikuchi ... next 1 line added     (11/13/'91)
        call corinit(newcor,nster,nmon,itstr,itestr,itmon,itemon)
c     
        flv%measp=nlat
        mfpnt=nlat
        mfpnt1=nlat
        flv%nfc=0
        call tfinitcalc
        call tmast(ilist(1,ilattp+1),ilist(1,ifmast),
     $       ilist(1,ifele1),rlist(ifpos))
        call twmov(ilist(2,ilattp+1),rlist(iftwis),nlat,ndim,.true.)
        if(.not. chguse)then
          do i=1,mfit1
            scale(i)=1.d0
          enddo
          scale(mfitnx)=pi2
          scale(mfitny)=pi2
          scale(mfitchi1)=pi/180.d0
          scale(mfitchi2)=pi/180.d0
          scale(mfitchi3)=pi/180.d0
          xixf=0.d0
          xiyf=0.d0
          flv%rsconv=1.d-9
          convgo=.false.
          trsize=.false.
          cellstab=.true.
          canon=.true.
          simulate=.true.
          absweit=.true.
          jitter=.true.
          lwake=.false.
          twake=.false.
          bipol=.true.
          cell=.false.
          gauss=.false.
          open98=.false.
          geomet=.false.
          fseed=.true.
          ideal=.false.
          codplt=.false.
          np0=1000
          ia=0
          call rsetgl('NP',dble(np0),ia)
          nturn=1
          trval=0.d0
          dp0=0.d0
          call tfsetsymbolr('XIX',3,xixf)
          call tfsetsymbolr('XIY',3,xiyf)
          call tfsetsymbolr('DPM',3,0.d0)
          call tfsetsymbolr('DP',2,dpmax)
          call tfsetsymbolr('ExponentOfResidual',18,2.d0)
          lfnp=1
          lfnstk(1)=5
          lfopen(1)=.false.
          lfret(1)=0
          lfrecl(1)=1
          lflinep(1)=1
          if(infl .ne. 5)then
            lfnp=2
            lfnstk(2)=infl
            lfret(2)=0
            lfrecl(2)=1
            lflinep(2)=1
          endif
        else
        endif
      else
        lfnp=lfnb
        lfnstk(lfnp)=0
        lfrecl(lfnp)=icslrecl()
        lflinep(lfnp)=icslrecl()
      endif
      iffserr=0
      if(chguse)then
        ios=0
        chguse=.false.
        go to 10
      endif
      levelr=0
      iflevel=0
      lfret(lfnp)=0
      lfopen(lfnp)=.false.
      lfni=lfnstk(lfnp)
      lfno=outfl
 2    if(lfnb .eq. 1)then
        if(lfni .ne. 5)then
          lfn1=lfno
        else
          if(igetgl1('$LOG$') .eq. 0)then
            lfn1=0
          else
            lfn1=lfno
          endif
        endif
      else
        lfn1=0
      endif
      call csrst(lfn1)
 10   continue
      if(iffserr .ne. 0)then
        if(lfnb .gt. 1)then
          ios=1
        else
          iffserr=0
        endif
      endif
      if(ios .gt. 0)then
        ios=0
        call tfclose(lfnp,lfnp,lfnstk,lfopen,lfret,lfrecl,
     $       lflinep,maxlfn,lfni,lfnb)
        if(lfnp .lt. lfnb)then
          go to 9000
        endif
        if(lfni .ne. 5)then
          lfn1=lfno
        else
          if(igetgl1('$LOG$') .eq. 0)then
            lfn1=0
          else
            lfn1=lfno
          endif
        endif
      elseif(ios .lt. 0)then
        ios=0
      endif
      call getwrd(word)
      if(ios .ne. 0)then
        go to 10
      endif
 12   if(word .eq. ' ')then
        go to 10
      elseif(word(1:1) .eq. '!')then
        go to 2
      endif
      call cssets(0)
      call tfprint(word,lfno,.false.,itt,nextt,exist)
c      write(*,*)'tffsa-1 ',word(1:lenw(word)),exist,itt,ios
      if(exist .or. ios .ne. 0)then
        go to 10
      endif
      if(iffserr .ne. 0 .and. lfnb .gt. 1)then
        go to 10
      endif
      if(word(1:1) .eq. '!')then
        go to 2
      endif
      if(word .eq. 'UNTIL')then
        if(levelr .le. 0)then
          call termes(lfno,'?UNTIL without REPEAT.',' ')
          go to 2
        endif
        ftest=getva(exist) .ne. 0.d0
        if(exist)then
          if(ftest)then
            levelr=levelr-1
            go to 4010
          endif
        endif
        nrpt(levelr)=nrpt(levelr)-1
        if(nrpt(levelr) .gt. 0)then
          call cssetp(irptp(levelr))
        else
          levelr=levelr-1
        endif
        ios=0
 4010   if(levelr .eq. 0)then
c          call cssetlinep(icslrecl())
          call cssetrec(.false.)
        endif
        go to 10
      elseif(abbrev(word,'REP_EAT','_'))then
        levelr=levelr+1
        nrpt1=int(getva(exist))
        if(.not. exist)then
          nrpt1=65535
        endif
        nrpt(levelr)=max(1,nrpt1)
        irptp(levelr)=icsmrk()
        call cssetrec(.true.)
        go to 10
      endif
      call tfif(word,iflevel,lfno,exist)
      if(exist)then
        go to 10
      endif
      call tflag(word,exist)
      if(exist)then
        go to 10
      endif
      call tgetfv
     1     (word,nlist,flv%nfc,lfno,ilist(1,ilattp+1),flv%icalc,
     $     flv%ncalc,
     1     flv%kfit,flv%fitval,flv%mfitp,flv%ifitp,flv%ifitp1,scale,
     $     maxcond,exist,err)
      if(err)then
        go to 2
      endif
      if(exist)then
        go to 10
      endif
      call terror(word,ilist(1,ilattp),rlist(ifpos),
     $     ilist(1,ifival),
     $     rlist(ifcoup),ilist(1,ifele),ilist(1,ifele1),
     $     ilist(1,ifklp),
     1     rlist(iferrk),ilist(1,ifmult),ilist(1,ifmast),
     $     new,lfno,exist,err)
      if(err)then
        go to 2
      endif
      if(exist)then
        go to 10
      elseif(new)then
        go to 12
      endif
      call tffile(word,lfnstk,lfopen,lfret,lfrecl,lflinep,
     $     maxlfn,lfni,lfno,lfnb,init,exist)
      if(lfnp .lt. lfnb)then
        go to 9000
      endif
      if(init)then
        call cssets(0)
        go to 2
      endif
      if(exist)then
        go to 10
      endif
      call tftrak(word,trdtbl,trval,lfno,exist)
      if(exist)then
        go to 10
      endif
      call tfgetlineps(word,lenw(word),nl,kax,1,irtc)
      if(irtc .eq. 0)then
        if(nl .gt. 0)then
          go to 7000
        endif
      else
        if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
          call tfreseterror
        endif
        irtc=0
      endif

      cmd=.false.
      if(word .eq. 'STOP' .or. word .eq. 'EXIT')then
        call tffsadjustvar
        call tfsave(word,.false.,ilist(1,ilattp+1),
     $       ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,
     1       ilist(1,ifklp),ilist(1,ifival),rlist(iftwis),
     $       rlist(iferrk))
        call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),
     $       rlist(ifpos),rlist(ifgamm),.false.)
        call corfree(newcor,nster,nmon,itstr,itestr,itmon,itemon)
        go to 8900
      elseif(word .eq. 'QUIT')then
        call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $       rlist(ifgamm),.false.)
        call corfree(newcor,nster,nmon,itstr,itestr,itmon,itemon)
        go to 8900
      elseif(word .eq. 'ABORT')then
        call tfresetsharedmap()
        stop
      elseif(word .eq. 'USE' .or. word .eq. 'VISIT')then
        visit=word .eq. 'VISIT'
        call peekwd(word,next)
        if(abbrev(word,'NOEXP_AND','_'))then
          expnd=.false.
          call cssetp(next)
        elseif(abbrev(word,'EXP_AND','_'))then
          expnd=.true.
          call cssetp(next)
        else
          expnd=.true.
        endif
        it=itfpeeko(k,next)
        call tfbeamline(k,iuse,ename,irtc)
        if(iuse .eq. 0)then
          if(ktfstringq(k) .or. ktfsymbolq(k))then
            word=tfgetstrs(k,nc)
            iuse=hsrchz(word(1:nc))
          else
            call peekwd(word,next)
            iuse=hsrchz(word(1:lenw(word)))
          endif
        endif
        if(idtype(iuse) .ne. icLINE)then
          if(irtc .ne. 0)then
            call tfaddmessage(ename,lenw(ename),6)
          endif
          if(visit)then
            call termes(lfno,
     $           'Missing beamline in VISIT.',' ')
          else
            call termes(lfno,
     $           'Missing beamline in USE.',' ')
          endif
          go to 2
        else
          if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
            call tfreseterror
          endif
          call tfclearlinep()
          call cssetp(next)
          call tffssaveparams(2,ilattp,0,err)
          expnd=expnd .and. .not. err
          call tclrpara(ilist(1,ilattp+1),nlat-1)
          if(visit)then
            call tffssaveparams(3,0,0,err)
            if(err)then
              call termes(lfno,'Too deep VISIT.',' ')
              go to 2
            endif
            call tffssaveparams(0,ilattp,flv%mcommon,err)
            call tfblocksym('`FitFunction',12)
            call tfblocksym('`FitValue',9)
            call tfblocksym('`ElementValues',14)
          else
            call tfunblocksym('`FitFunction',12,.false.)
            call tfunblocksym('`FitValue',9,.false.)
            call tfunblocksym('`ElementValues',14,.false.)
            call tffsfree
            call tffsfreebzl
            if(lattuse .eq. lattredef)then
              call tclrline(lattredef)
              lattredef=0
            endif
          endif
          if(ilist(2,idval(iuse)) .le. 0 .or. expnd)then
            call expnln(iuse)
          endif
          call filaux(idval(iuse))
          ilattp=ilist(2,idval(iuse))
          lattuse=ilattp
          call tclrpara(ilist(1,ilattp+1),ilist(1,ilattp))
          dleng =rlist(ilist(2,ilattp)+1)*rgetgl1('FSHIFT')
          ename=pname(iuse)
          call tmovb(ename,flv%blname,MAXPNAME)
          call tfsetbeamlinename(ename)
          chguse=.true.
          go to 101
        endif
      elseif(word .eq. 'BYE')then
        call peekwd(word,next)
        if(word .eq. 'ALL')then
          call cssetp(next)
          byeall=.true.
        else
          byeall=.false.
        endif
        call tffssaveparams(4,0,0,err)
        if(err)then
          if(byeall)then
            go to 10
          else
            call termes(lfno,'BYE without VISIT.',' ')
            go to 2
          endif
        endif
        call tclrpara(ilist(1,ilattp+1),ilist(1,ilattp))
        call tffsfree
        call tffsfreebzl
        if(byeall)then
          call tffssaveparams(-1,0,2,err)
        endif
        call tffssaveparams(1,ilattp,flv%mcommon,err)
        lattuse=ilattp
        nlat=ilist(1,ilattp)+1
        dleng =rlist(ilist(2,ilattp)+1)*rgetgl1('FSHIFT')
        call tmovb(flv%blname,ename,MAXPNAME)
        call tfsetbeamlinename(ename)
        call tfclearlinep()
        call tfunblocksym('`FitFunction',12,.true.)
        call tfunblocksym('`FitValue',9,.true.)
        call tfunblocksym('`ElementValues',14,.true.)
      elseif(word .eq. 'SPLIT')then
        call termes(lfno,
     $       'SPLIT is obsolete.   Use OFFSET of a marker.',' ')
      elseif(abbrev(word,'MAXI_TERATION','_'))then
        call tfgeti(flv%itmax,1.d0,word,lfno,exist)
        go to 31
      elseif(abbrev(word,'ATT_RIBUTE','_'))then
        call tfattr(word,ilist(1,ilattp+1),ilist(1,ifele),
     $       ilist(1,ifele1),ilist(1,ifival),
     $       rlist(ifcoup),rlist(ifaux),rlist(iferrk),
     $       ilist(1,ifklp),ilist(1,ifmult),lfno,exist,
     $       kffs,irtcffs,lfnb .gt. 1)
        go to 30
      elseif(word .eq. 'SAVE')then
        call tffsadjustvar
        call tfsave(word,.true.,ilist(1,ilattp+1),
     $       ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,
     1       ilist(1,ifklp),ilist(1,ifival),rlist(iftwis),
     $       rlist(iferrk))
      elseif(abbrev(word,'EXPAND','_'))then
        call tffsadjustvar
        call tfsave(word,.true.,ilist(1,ilattp+1),
     $       ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,
     1       ilist(1,ifklp),ilist(1,ifival),rlist(iftwis),
     $       rlist(iferrk))
      elseif(word .eq. 'RESET')then
        call tfrst(word,.true.,ilist(1,ilattp+1),
     $       ilist(1,ifvarele),ilist(1,ifvvar),
     $       ilist(1,ifivcomp),rlist(ifvalvar),flv%nvar,
     $       ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,
     1       ilist(1,ifklp),ilist(1,ifival),ilist(1,ifele),
     $       ilist(1,ifele1),rlist(iferrk),
     $       rlist(ifcoup))
      elseif(word .eq. 'FREE' .or. word .eq. 'FIX')then
        frefix=word .eq. 'FREE'
        call peekwd(word,next)
        if(word .eq. ' ')then
ckikuchi ... next 5 lines added     (8/17/'90)
          call mcmess(lfno)
          if(frefix) then
            call mcfre(newcor)
          else
            call mcfix(ilist(1,ilattp+1),rlist(iftwis),
     $           rlist(ifgamm),newcor)
          endif
          go to 10
        endif
        call tffsfreefix(ilist(1,ilattp+1),frefix,
     $       ilist(1,ifvarele),ilist(1,ifvvar),
     $       ilist(1,ifivcomp),rlist(ifvalvar),
     $       flv%nvar,nve,
     $       ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,
     $       ilist(1,ifival),ilist(1,ifklp),nele,
     $       rlist(iferrk),ilist(1,ifmult),ilist(1,ifele),
     $       ilist(1,ifele1),
     $       nlat,nlist,lfno)
      elseif(word .eq. 'FIT')then
        call getwdl2(word,wordp)
        if(word .eq. 'ALL')then
          mfpnt=1
          mfpnt1=nlat
        else
          mfpnt=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
          if(.not. exist)then
            mfpnt1=mfpnt
            go to 12
          endif
          call getwdl2(word,wordp)
          mfpnta=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
          if(exist)then
            mfpnt1=max(mfpnt,mfpnta)
            mfpnt=min(mfpnt,mfpnta)
          else
            mfpnt1=mfpnt
            go to 12
          endif
        endif
      elseif(abbrev(word,'FITP_OINTS','_'))then
        call tfgeti(nfp,1.d0,word,lfno,exist)
        if(.not. exist)then
          go to 2
        endif
        do kk=1,flv%nfc
          if(flv%kfit(kk) .le. mfit)then
            if(abs(flv%mfitp(kk)) .gt. 2)then
              flv%mfitp(kk)=sign(abs(nfp)+1,flv%mfitp(kk))
            endif
          endif
        enddo
      elseif(abbrev(word,'CAL_CULATE','_'))then
        call tfcalc(word,nlist,flv%icalc,flv%ncalc,mfpnt,mfpnt1,exist)
        if(exist)then
          expndc=.not. tfvcomp(ilist(1,ifivcomp),
     $         ilist(1,ifvarele),ilist(1,ifvvar),
     $         ilist(1,ifival),nele,flv%nvar)
          if(abbrev(word,'NOEXP_AND','_'))then
            word=' '
            expndc=.false.
          elseif(abbrev(word,'EXP_AND','_'))then
            word=' '
            expndc=.true.
          endif
          fitflg=.false.
          go to 1000
        else
          go to 12
        endif
      elseif(word .eq. 'GO')then
        expndc=.not. tfvcomp(ilist(1,ifivcomp),
     $       ilist(1,ifvarele),ilist(1,ifvvar),
     $       ilist(1,ifival),nele,flv%nvar)
        call peekwd(word,next)
        if(abbrev(word,'NOEXP_AND','_'))then
          call cssetp(next)
          expndc=.false.
        elseif(abbrev(word,'EXP_AND','_'))then
          call cssetp(next)
          expndc=.true.
        endif
        word=' '
        fitflg=.true.
        convgo=.false.
        go to 1000
      elseif(abbrev(word,'REC_OVER','_'))then
        ifvalvar2=ifvalvar+nve
        do i=1,flv%nvar
          v=rlist(ifvalvar+i-1)
          rlist(ifvalvar+i-1)=rlist(ifvalvar2+i-1)
          rlist(ifvalvar2+i-1)=v
        enddo
        call tfsetv(ilist(1,ilattp+1),ilist(1,ifvarele),
     $       ilist(1,ifvvar),ilist(1,ifivcomp),
     $       rlist(ifvalvar),flv%nvar,nele,
     $       ilist(1,ifklp),ilist(1,ifival),rlist(ifcoup),
     $       rlist(iferrk),ilist(1,ifele),
     $       ilist(1,ifele1),nlat)
      elseif(abbrev(word,'T_YPE','_'))then
        call tfsetparam
        call tftype(ilist(1,ilattp+1),ilist(1,ifival),
     $       ilist(1,ifele1),rlist(iferrk),ilist(1,ifklp),lfno,word,
     $       rlist(iftwis),emx,emy,dpmax,nlat,nele,ndim)
        go to 12
      elseif(abbrev(word,'DISP_LAY','_'))then
        call tfsetparam
        call tfdisp(word,id1,id2,
     1       ilist(1,ilattp+1),rlist(iftwis),rlist(ifpos),rlist(ifgeo),
     $       rlist(ifgamm),dp0,rlist(ifsize),
     $       ilist(1,ifele),ilist(1,ifele1),ilist(1,ifival),
     $       scale,ilist(1,ifmult),lfno,exist)
        go to 30
      elseif(word .eq. 'SCALE')then
        call tscale(nlist,scale,lfno)
      elseif(word.eq.'MAPANA') then
        call gosadpls(ilist(1,ilattp+1),ilist(1,ifklp),
     $       ilist(1,ifele),lfno)
      elseif(abbrev(word,'AP_ERTURE','_'))then
        call tfsetparam
        call tfaprt(ilist(1,ilattp+1),
     $       ilist(1,ifmult),lfno,word,rlist(iftwis),rlist(ifgamm))
        go to 12
      elseif(abbrev(word,'MAT_RIX','_'))then
        call tfmat(ilist(1,ilattp+1),ilist(1,ifmult),rlist(iftwis),
     $       rlist(ifgamm),lfno,word,wordp,exist)
        go to 30
      elseif(abbrev(word,'CHRO_MATICITY','_'))then
        call tfsetparam
        i=0
        call tfchro(ilist(1,ilattp+1),ilist(1,ifmult),
     1       rlist(itwisso(1,i,1)),rlist(itwisso(1,i,2)),
     $       rlist(itwisso(1,i,3)),rlist(itwisso(1,i,7)),
     1       rlist(itwisso(1,i,4)),rlist(itwisso(1,i,5)),
     $       rlist(itwisso(1,i,6)),lfno)
      elseif(word .eq. 'RADINT')then
        call tfsetparam
        call intgrl( ilist(1,ilattp+1),rlist(iftwis),0,1.d0,lfno)
      elseif(word .eq. 'DIMAD')then
        call tfsetparam
        call tdimad(ilist(1,ilattp+1),ilist(1,ifmult),lfno)
      elseif(word .eq. 'ZAP')then
        call tfsetparam
        apert=getva(exist)
        if(.not. exist)then
          apert=1.d0
        endif
        i=0
        call tfzap(apert,rlist(ifpos),
     1       rlist(itwisso(1,i,1)),rlist(itwisso(1,i,2)),
     $       rlist(itwisso(1,i,7)),rlist(itwisso(1,i,8)),
     1       rlist(itwisso(1,i,4)),rlist(itwisso(1,i,5)),
     $       lfno)
c     Next several lines are added by N. Yamamoto Apr. 25,'93
C     23/06/92 212101550  MEMBER NAME  TRCOD    *.FORT     M  E2FORT
      else if(abbrev(word,'TRC_OD','_'))then
         call getwdl2(word,wordp)
c     write(6,*)'trcod scan ', word ,abbrev(word,'FIXCOD','_')
         if(abbrev(word,'FIXCOD','_')) then
            call getwdl2(word,wordp)
            meas0=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
            if(exist)then
               flv%measp=meas0
            endif
            title=tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodf(ilist(1,ilattp+1),flv%measp,name,
     $           ilist(1,ifmult),
     &           .false.,word,title,case,exist,
     &           rlist(iftwis),lfno)
            if(.not. exist)then
               go to 12
            endif
         else if(abbrev(word,'MVCOD','_')) then
            call getwdl2(word,wordp)
            meas0=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
            if(.NOT. exist)then
               meas0=flv%measp
            endif
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodM(ilist(1,ilattp+1),meas0,name,ilist(1,ifmult),
     &           .false.,word,title,case,exist,
     &           rlist(iftwis),lfno)
            if(.not. exist)then
               go to 12
            endif
         else if(abbrev(word,'STN_COD','_')) then
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodS(ilist(1,ilattp+1),meas0,name,ilist(1,ifmult),
     &           .false.,word,title,case,exist,
     &           rlist(iftwis),lfno)
            go to 10
         else if(abbrev(word,'ALIGN','_')) then
            call getwdl(word)
            call nalign(ilist(1,ilattp+1),ilist(1,ifmult),
     $           ilist(1,ifmast),word,lfno,exist)
            go to 10
         else
            meas0=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
            if(exist)then
               flv%measp=meas0
            endif
            call elname(ilist(1,ilattp+1),flv%measp,
     $           ilist(1,ifmult),name)
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcoda(ilist(1,ilattp+1),flv%measp,name,
     $           ilist(1,ifmult),.false.,word,title,case,exist,
     &           rlist(iftwis),xa,ya,xxa,xya,yya,lfno)
         end if
         if(.not. exist)then
            go to 12
         endif
      else if(abbrev(word,'PUSHS_EED ','_'))then
         call tfevalb('NISTACK$OBJ@Push[]',18,kx,irtc)
         go to 10
      else if(abbrev(word,'POPS_EED ','_'))then
         call tfevalb('NISTACK$OBJ@Pop[]',17,kx,irtc)
         go to 10
      else if(abbrev(word,'EXCGS_EED ','_'))then
         call tfevalb('NISTACK$OBJ@Exchange[]',22,kx,irtc)
         go to 10
      else if(abbrev(word,'PEEKS_EED ','_') .OR.
     &        abbrev(word,'PKS_EED','_')        )then
         call tfgeti(i,1.d0,word,lfno,exist)
         if( .not. exist) i=0
         write(word,'(A,I12,A)')'NISTACK$OBJ@Peek[',i,']'
         call tfevalb(word,len_trim(word),kx,irtc)
         go to 31
      else if(abbrev(word,'DELW_AVE ','_')) then
c    syntax DELWave wave_length, amplitude, initial_phase,
c                   directio_vector(angle)
        wl=getva(exist)
        if(.not. exist) then
          call termes(lfno,'? syntax:DELW <WAVE LENGTH> <amplitude>'//
     &     '<initial_phase><direction>',' ')
          go to 10
        end if
        wa=getva(exist)
        if(.not. exist) then
          call termes(lfno,'? syntax:DELW <wave length> <AMPLITUDE>'//
     &     '<initial_phase><direction>',' ')
          go to 10
        end if
        wp=getva(exist)
        if(.not. exist) then
          wp=0.d0
        end if
        wd=getva(exist)
        if(.not. exist) then
          wd=0.d0
        end if
        call ndelw(wl,wa,wp,wd,
     &             ilist(1,ilattp+1),ilist(1,ifele),rlist(ifcoup),
     $       rlist(ifgeo),rlist(ifpos),ilist(1,ifmast))
        go to 10
c End of the lines added by N. Yamamoto Apr. 25, '93
      else
        cmd=.true.
      endif
      if(.not. cmd)then
        go to 10
      endif
      cmd=.false.
      if(word .eq. 'VARY')then
        call tfchgv(ilist(1,ilattp+1),ilist(1,ifival),
     $       ilist(1,ifklp),lfno)
        call tfinitvar(ilist(1,ifvarele),ilist(1,ifvvar),
     $       ilist(1,ifivcomp),
     $       rlist(ifvalvar),flv%nvar,ilist(1,ilattp+1),ilist(1,ifklp),
     $       ilist(1,ifival),rlist(iferrk),nlat,nele)
        go to 10
      elseif(abbrev(word,'REJ_ECT','_'))then
        call tfrej
     1   (word,nlist,flv%nfc,mfpnt,mfpnt1,
     1    flv%kfit,flv%mfitp,flv%ifitp,flv%ifitp1,exist)
        go to 30
      elseif(abbrev(word,'COUP_LE','_'))then
        call tfcoup(ilist(1,ilattp+1),rlist(ifcoup),
     $       ilist(1,ifele),ilist(1,ifele1),
     $       ilist(1,ifklp),ilist(1,ifival),rlist(iferrk),
     $       ilist(1,ifmult),lfno,exist)
        call tffsadjustvar
        if(.not. exist)then
          go to 2
        endif
        go to 10
      elseif(abbrev(word,'DECOUP_LE','_'))then
        call tfdecoup(ilist(1,ilattp+1),rlist(ifcoup),
     $       ilist(1,ifele),ilist(1,ifele1),
     $       ilist(1,ifklp),ilist(1,ifmult),lfno)
        call tffsadjustvar
        go to 10
      elseif(abbrev(word,'INDEP_ENDENT','_'))then
        call tfindep(ilist(1,ilattp+1),ilist(1,ifele),ilist(1,ifmult))
        call tffsadjustvar
        go to 10
      elseif(word .eq. 'SBUNCH')then
        iwakepold=ifwakep
        call tfgetr(rlist(iwakepold+1),1.d0,word,lfno,exist)
        go to 31
      elseif(word .eq. 'DBUNCH')then
        iwakepold=ifwakep
        call tfdbun(word,rlist(ilist(2,iwakepold)),ilist(1,iwakepold),
     1              lfno,err)
        go to 32
      elseif(word .eq. 'SLICE')then
        iwakepold=ifwakep
        ns=ilist(1,iwakepold+2)
        call tfgeti(ns,1.d0,word,lfno,exist)
        if(ns .le. 0)then
          call termes(lfno,'?Parameter out of range in SLICE.',' ')
          go to 2
        endif
        ilist(1,iwakepold+2)=ns
        go to 31
      elseif(abbrev(word,'ORI_GIN','_') .or. word .eq. 'ORG')then
        iorgr=igelm(ilist(1,ilattp+1),word,ilist(1,ifmult),exist)
        if(.not. exist)then
          iorgr=1
          go to 12
        endif
        do i=1,3
          geo0(i,4)=getvad(exist,word,geo0(i,4))
          if( .not. exist)then
            go to 12
          endif
        enddo
        do i=1,3
          chi0(i)=getvad(exist,word,chi0(i))
          if( .not. exist)then
            exit
          endif
        enddo
        call tsetg(geo0,chi0)
        if(.not. exist)then
          go to 12
        endif
        go to 10
      elseif(word .eq. 'SHOW')then
        call tshow(flv,ilist(1,ilattp+1),
     $       scale,nlist,ilist(1,ifmult),
     $       kffs,irtcffs,lfnb .gt. 1,lfno)
        go to 10
      elseif(abbrev(word,'STAT_US','_'))then
        call cputime(ctime1,irtc0)
        call tfsetparam
        lpw=min(itfgetrecl()-1,255)
        call twbuf(word,lfno,1,0,0,0)
        call twbuf('CTIME='//autofg((ctime1-flv%ctime0)*1.d-6,'8.3'),
     1          lfno,1,lpw,8,1)
        call twbuf('sec  DT='//autofg((ctime1-flv%ctime2)*1.d-6,'8.3'),
     1          lfno,1,lpw,8,1)
        call twbuf('sec',lfno,1,lpw,8,1)
        flv%ctime2=ctime1
        call twbuf(word,lfno,1,0,0,-1)
        do i=1,nflaga
          if(fname(i) .ne. ' ')then
            if(flags(i))then
              call twbuf(fname(i),lfno,1,lpw,8,1)
            else
              if(sino(i) .ne. ' ')then
                call twbuf(sino(i),lfno,1,lpw,8,1)
              else
                call twbuf('NO'//fname(i),lfno,1,lpw,8,1)
              endif
            endif
          endif
        enddo
        call twbuf(word,lfno,1,0,0,-1)
        call twbuf('CHARGE='//autofg(charge,'S7.4'),lfno,1,lpw,8,1)
        flv%rsconv=rfromd(kxsymbolv('CONVERGENCE',10))
c        flv%rsconv=rlist(ktlookup('CONVERGENCE'))
        call twbuf('CONVERGENCE='//autofg(flv%rsconv,'S8.5'),
     $       lfno,1,lpw,8,1)
        call twbuf('DP='//autofg(dpmax,'S8.5'),lfno,1,lpw,8,1)
        call twbuf('DP0='//autofg(dp0,'S8.5'),lfno,1,lpw,8,1)
        call twbuf('EMITDIV='//autofg(emidiv,'S7.4'),lfno,1,lpw,8,1)
        call twbuf('EMITDIVB='//autofg(emidib,'S7.4'),lfno,1,lpw,8,1)
        call twbuf('EMITDIVQ='//autofg(emidiq,'S7.4'),lfno,1,lpw,8,1)
        call twbuf('EMITDIVS='//autofg(emidis,'S7.4'),lfno,1,lpw,8,1)
        call twbuf('EMITX='//autofg(emx,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('EMITY='//autofg(emy,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('GCUT='//autofg(tgetgcut(),'S7.4'),lfno,1,lpw,8,1)
        call twbuf('MINCOUP='//autofg(coumin,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('MOMENTUM='//autofg(pgev,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('PBUNCH='//autofg(pbunch,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('NBUNCH='//autofg(anbunch,'S9.6'),lfno,1,lpw,8,1)
        call twbuf('XIX='//autofg(xixf/pi2,'S7.4'),lfno,1,lpw,8,1)
        call twbuf('XIY='//autofg(xiyf/pi2,'S7.4'),lfno,1,lpw,8,1)
        write(word,'(A,I6)')'NP=',np0
        call twbuf(word,lfno,1,lpw,8,1)
        write(word,'(A,I4)')'MAXITERATION=',flv%itmax
        call twbuf(word,lfno,1,lpw,8,1)
        call twbuf(word,lfno,1,lpw,8,-1)
        call twelm(lfno,ilist(1,ilattp+1),ilist(1,ifmult),
     $       mfpnt,mfpnt1,'FIT',lpw,8)
        call twelm(lfno,ilist(1,ilattp+1),ilist(1,ifmult),
     $       flv%measp,0,'MEA_SURE',lpw,8)
        call twelm(lfno,ilist(1,ilattp+1),ilist(1,ifmult),
     $       iorgr,0,'ORG',lpw,8)
        call twelm(lfno,ilist(1,ilattp+1),ilist(1,ifmult),
     $       id1,id2,'DISP_LAY',lpw,8)
        call twbuf(word,lfno,1,lpw,8,-1)
        go to 10
      elseif(abbrev(word,'VAR_IABLES','_') .or. word .eq. 'VARS')then
        call tfinitvar(ilist(1,ifvarele),ilist(1,ifvvar),
     $       ilist(1,ifivcomp),
     $       rlist(ifvalvar),flv%nvar,ilist(1,ilattp+1),ilist(1,ifklp),
     $       ilist(1,ifival),rlist(iferrk),nlat,nele)
        call tfvars(ilist(1,ifvarele),ilist(1,ifvvar),
     $       rlist(ifvalvar),ilist(1,ifivcomp),
     $       flv%nvar,ilist(1,ilattp+1),ilist(1,ifklp),ilist(1,ifmult),
     $       ilist(1,ifival),rlist(ifaux),
     $       ilist(1,ifele),ilist(1,ifele1),
     $       rlist(ifcoup),
     $       kffs,irtcffs,lfnb .gt. 1,lfno)
        go to 10
      elseif(abbrev(word,'RENUM_BER','_'))then
        call tffsrenumber(ilist(1,ilattp+1),ilist(1,ifmult),
     $       ilist(1,ifele),ilist(1,ifele1),ilist(1,ifele2),
     $       ilist(1,ifklp),lfno)
        go to 10
      elseif(abbrev(word,'REV_ERSE','_'))then
        axi=-axi
        ayi=-ayi
        epxi=-epxi
        epyi=-epyi
        dpxi=-dpxi
        dpyi=-dpyi
        r2i=-r2i
        r3i=-r3i
        go to 10
      elseif(word .eq. 'DRAW')then
        call tfsetparam
        call tfevalb('CANVASDRAW[]',12,kx,irtc)
        if(irtc .ne. 0 .or. ktfnonstringq(kx))then
          title=Tfgetstrv('TITLE')
          case=Tfgetstrv('CASE')
          call twsdrw(ilist(1,ilattp+1),rlist(ifpos),ilist(1,ifele),
     $         ilist(1,ifmult),word,wordp,lfno,
     1         rlist(iftwis),rlist(ifgamm),
     1         0,rlist(itmon),rlist(itemon),nmon,
     1         title,case,exist)
        else
          word=tfgetstr(kx,nc)
          exist=nc .eq. 0
        endif
        go to 30
      elseif(word .eq. 'GEO')then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call geodrw(ilist(1,ilattp+1),rlist(ifgeo),
     $       word,lfno,title,case)
        go to 10
      elseif(word .eq. 'TDR')then
        call ttdr(lfno,err)
        go to 32
      elseif(abbrev(word,'PRI_NT','_'))then
        word=' '
        call tfprint(word,lfno,.true.,itt,nextt,exist)
        go to 30
      else
        cmd=.true.
      endif
      if(.not. cmd)then
        go to 10
      endif

      cmd=.false.
      if(abbrev(word,'MEA_SURE','_'))then
        call tfsetparam
        call getwdl2(word,wordp)
        meas0=ielm(ilist(1,ilattp+1),wordp,1,ilist(1,ifmult),exist)
        if(exist)then
          flv%measp=meas0
        endif
        call elname(ilist(1,ilattp+1),flv%measp,ilist(1,ifmult),name)
ckikuchi ... 1 line modified
        if(twake .or. lwake)then
          iwakepold=ifwakep
          ns=ilist(1,iwakepold+2)
          np0=max(1,np0/ilist(1,iwakepold)/ns/8)*ns*
     $         ilist(1,iwakepold)*8
          if(ifsize .eq. 0)then
            ifsize=ktaloc(nlat*21)
            ilist(2,iwakepold+6)=int(ifsize)
          endif
        endif
        trpt=.true.
        trf0=0.d0
        vcphic=0.d0
        vcalpha=1.d0
        title=tfgetstrv('TITLE')
        case=tfgetstrv('CASE')
        call trackb(ilist(1,ilattp+1),flv%measp,name,.true.,
     $       word,title,case,exist,
     $       kffs,irtcffs,lfnb .gt. 1,
     1       xa,ya,xxa,xya,yya,lfno)
        if(.not. exist)then
          go to 12
        endif
        go to 10
      elseif(abbrev(word,'ORBITJ_ITTER','_'))then
        call tfojit(rlist(iftwis),lfno)
        go to 10
      elseif(abbrev(word,'PUTM_EASURE','_'))then
         call  nfputm(lfno,xa,ya,xxa,xya,yya)
         goto 10
      elseif(abbrev(word,'GET_MEASURE','_'))then
         call tfgetm(flv,flv%mfitp,xa,ya)
         go to 10
      elseif(abbrev(word,'LOGM_EASURE','_'))then
         write(lfno,9771)xa,ya,xxa,xya,yya
 9771    format(1p,5g15.7)
         go to 10
C   19/12/90 302031849  MEMBER NAME  CORRECT  *.FORT     M  E2FORT
      elseif(abbrev(word,'STE_ERING','_'))then
ckiku   call tfstr(word,ilist(1,ilattp+1),ist,nstr)
        call mcmess(lfno)
        call mcstr(word,ilist(1,ilattp+1),
     $       ilist(1,ifmult),ilist(1,ifmast),itstr,itestr,nster,
     $       .false.,lfno)
        go to 12
      elseif(word.eq.'COUPLESTE') then
        call mcoupste(word,ilist(1,ilattp+1),ilist(1,ifmult),
     $       rlist(itstr),nster,lfno)
        go to 12
      elseif(word.eq.'READSTE')then
        call mcstr(' ',ilist(1,ilattp+1),ilist(1,ifmult),
     $       ilist(1,ifmast),itstr,itestr,nster,.true.,lfno)
        call preadstr(word,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),rlist(itstr),nster,lfno)
        go to 12
      elseif(word.eq.'WRITESTE')then
        call pwrtstr(word,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),rlist(itstr),nster,lfno)
        go to 12
      elseif(abbrev(word,'MON_ITOR','_'))then
        call mcmess(lfno)
        call mcmon(word,ilist(1,ilattp+1),ilist(1,ifmult),rlist(ifpos),
     $       ilist(1,ifmast),itmon,itemon,nmon,.false.,
     &             lfno)
        go to 12
      elseif(word.eq.'READMON')then
        call mcmon(' ',ilist(1,ilattp+1),ilist(1,ifmult),rlist(ifpos),
     $       ilist(1,ifmast),itmon,itemon,nmon,.true.,
     &             lfno)
        call preadmon(word,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),rlist(itmon),nmon,lfno)
        go to 12
      elseif(word.eq.'WRITEMON')then
        call pwrtmon(word,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),rlist(itmon),nmon,lfno)
        go to 12
      elseif(word.eq.'KILL')then
        call pkill(word,wordp,ilist(1,ilattp+1),ilist(1,ifmult),
     $       rlist(itmon),nmon,rlist(itstr),
     1       nster,lfno)
        go to 12
      elseif(abbrev(word,'COR_RECT','_')) then
        call mcmess(lfno)
        call mccor(word,wordp,
     $       ilist(1,ilattp+1),rlist(ifpos),rlist(iftwis),
     $       rlist(ifgamm),ilist(1,ifmult),ilist(1,ifmast),flv%kfit,
     $       flv%ifitp,
     $       flv%mfitp,flv%fitval,flv%nfc,rlist(itstr),rlist(itestr),
     $       nster,
     $       rlist(itmon),rlist(itemon),nmon,newcor,lfno)
        go to 12
      elseif(abbrev(word,'KI_CK','_'))then
        call mckick(word,wordp,ilist(1,ilattp+1),ilist(1,ifmult))
        go to 12
      elseif(abbrev(word,'CLE_AR','_')) then
        call mclear(word,ilist(1,ilattp+1),rlist(itstr),nster)
        goto 12
      elseif(abbrev(word,'UNDO_CORRECT','_')) then
        call pundo(ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(ifgamm),0,1.d0+dp0)
        go to 10
      elseif(word.eq.'SUM+'.or.word.eq.'SUM-'.or.word.eq.'SUM'.or.
     1       abbrev(word,'SUMC_LR','_')) then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call pstati(word,ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(itmon),rlist(itemon),nmon,
     1       rlist(itstr),nster,title,case,lfno)
        go to 12
      elseif(abbrev(word,'STO_RE','_')) then
        call mstore(word,ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(itstr),nster,rlist(itmon),
     1       rlist(itemon),nmon)
        go to 12
      elseif(abbrev(word,'RECA_LL','_')) then
        call mrecal(word,ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(itstr),nster,rlist(itmon),rlist(itemon),nmon,lfno)
        go to 12
      elseif(abbrev(word,'DIR_ECTORY','_'))then
        call mfdir(word,lfno)
        goto 12
      elseif(abbrev(word,'DEL_ETE','_'))then
        call mfdel(word)
        goto 12
      elseif(word.eq.'ADD'.or.word.eq.'SUB'.or.word.eq.'SWAP'
     z       .or.abbrev(word,'FAC_TOR','_')
     z       .or.abbrev(word,'DIV_IDE','_')
     1       .or.word.eq.'PUSH'.or.word.eq.'DROP'            ) then
        call mstack(word,ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(itstr),nster,rlist(itmon),
     1       rlist(itemon),nmon)
        goto 12
      elseif(word.eq.'TRIM') then
        call ptrim(word,ilist(1,ilattp+1),rlist(iftwis),rlist(itmon),
     $       rlist(itemon),nmon,lfno)
        goto 12
      elseif(word.eq.'TCOD') then
        call mcrcod(ilist(1,ilattp+1),rlist(iftwis),ilist(1,ifmult),
     $       rlist(itmon),rlist(itemon),nmon,
     $       .false.,.false.,.true.,lfno)
        go to 10
      elseif(word.eq.'MCOD') then
        call mcrcod(ilist(1,ilattp+1),rlist(iftwis),ilist(1,ifmult),
     $       rlist(itmon),rlist(itemon),
     1       nmon,.true.,.false.,.true.,lfno)
        go to 10
      elseif(abbrev(word,'EPS_ILON','_')) then
        call mcepst(word,lfno)
        go to 12
      elseif(abbrev(word,'DPS_HIFT','_')) then
        call mdpmax(word,lfno)
        go to 12
      elseif(abbrev(word,'ECOR_R','_')) then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call pecorr(word,ilist(1,ilattp+1),rlist(ifpos),
     $       rlist(iferrk),title,case,lfno)
        goto 12
      elseif(word .eq. 'BUMP')then
        call getwdl(word)
        call pcbak(ilist(1,ilattp+1),rlist(iftwis))
        if(abbrev(word,'I_TERATION','_')) then
          call getwdl(word)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.true.)
          call pbump(ilist(1,ilattp+1),rlist(iftwis),
     $         rlist(ifgamm),0,ilist(1,ifmult),rlist(itstr),
     1               rlist(itestr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%fitval,flv%nfc,
     1               .true.,.false.,.false.,lfno)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.false.)
        else
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.true.)
          call pbump(ilist(1,ilattp+1),rlist(iftwis),
     $         rlist(ifgamm),0,ilist(1,ifmult),rlist(itstr),
     1               rlist(itestr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%fitval,flv%nfc,
     1               .false.,.false.,.false.,lfno)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.false.)
        endif
        go to 12
      elseif(word.eq.'LBUMP') then
        call getwdl(word)
        call pcbak(ilist(1,ilattp+1),rlist(iftwis))
        call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $       flv%nfc,.true.)
        call pbump(ilist(1,ilattp+1),rlist(iftwis),rlist(ifgamm),
     $       0,ilist(1,ifmult),rlist(itstr),rlist(itestr),
     1             nster,flv%kfit,flv%ifitp,flv%mfitp,flv%fitval,
     $       flv%nfc,.false.,.true.,
     1             .false.,lfno)
        call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $       flv%nfc,.false.)
        go to 12
      elseif(abbrev(word,'TOR_','_')) then
        call ptol(word,ilist(1,ilattp+1),ilist(1,ifmult),
     $       ilist(1,ifmast),rlist(iftwis),
     $       rlist(ifgamm),0,1.d0+dp0,
     1            rlist(itstr),nster,lfno)
        goto 12
      elseif(word.eq.'VBUMP') then
        call pvbump(word,wordp,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),ilist(1,ifmast),rlist(itstr),nster,
     $       nlist,flv%kfit,flv%ifitp,flv%mfitp,flv%fitval,flv%nfc,lfno)
        goto 12
      elseif(word.eq.'BTUNE') then
        if(ifsize .eq. 0)then
          ifsize=ktaloc(nlat*21)
          ilist(2,iwakepold+6)=int(ifsize)
        endif
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call petune(word,ilist(1,ilattp+1),rlist(iftwis),
     $       ilist(1,ifmult),rlist(ifgamm),rlist(ifsize),nlist,
     $       rlist(itstr),rlist(itestr),nster,title,case,lfno)
        goto 12
      elseif(word.eq.'BALS') then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call pasex(word,wordp,
     $       ilist(1,ilattp+1),rlist(ifpos),rlist(iftwis),
     $       ilist(1,ifmult),ilist(1,ifmast),rlist(ifgamm),title,case,
     $       rlist(itstr),rlist(itestr),nster,rlist(itmon),
     $       rlist(itemon),nmon,lfno)
        go to 12
      elseif(word.eq.'PETIN'.or.word.eq.'PETOUT') then
        call petcod(word,rlist(iftwis),rlist(itmon),nmon)
        goto 12
      elseif(abbrev(word,'WR_ITE','_')) then
        call getwdl(word)
        if(word.eq.'RMATQ'.or.word.eq.'RMAT') then
          call pwrite(word,ilist(1,ilattp+1),rlist(iftwis),
     $         rlist(ifgamm),rlist(ifpos),rlist(itmon),nmon)
        elseif(abbrev(word,'LAT_TICE','_').or.
     1         abbrev(word,'PS_NAME','_')) then
          call pwrlat(word,wordp,ilist(1,ilattp+1),ilist(1,ifmult),lfno)
        endif
        goto 12
      else
        cmd=.true.
      endif
      if(.not. cmd)then
        go to 10
      endif
      cmd=.false.
      if(word .eq. 'WAKE')then
        lfnl0=lfn1
        lfn1=lfno
        call tfwake(word,ilist(1,ilattp+1),ilist(1,ifmult),lfno,err)
        lfn1=lfnl0
        if(err)then
          go to 2
        endif
        go to 10
      elseif(abbrev(word,'DELC_OR','_'))then
        call tcorr(word,ilist(1,ilattp+1),rlist(ifpos),ilist(1,ifmult),
     $       ilist(1,ifmast),lfno)
        go to 10
      elseif(abbrev(word,'LTR_ACK','_'))then
        call tfltra(ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(ifgamm),word,lfno)
        go to 10
      elseif(abbrev(word,'EMIT_TANCE','_'))then
        pspan=getva(exist)
        if(exist)then
          ia=0
          call rsetgl('PSPAN',pspan,ia)
        endif
        trpt0=trpt
        trpt=.false.
        call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $       rlist(ifgamm),.true.)
        if(ifsize .eq. 0 .and. codplt)then
          ifsize=ktaloc(nlat*21)
          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call temitf(ilist(1,ilattp+1),rlist(iftwis),
     $       rlist(ifsize),rlist(ifgamm),ndim,codplt,lfni,lfno)
        trpt=trpt0
        if(.not. codplt)then
          call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $         rlist(ifgamm),.true.)
        endif
        emx=max(4.d-13/p0,rgetgl1('EMITX'))
        emy=max(4.d-13/p0,rgetgl1('EMITY'))
c        dpmax=rgetgl1('SIGE')
c        rlist(itlookup('DP',ivtype))=dpmax
        gauss=.true.
        go to 10
      elseif(abbrev(word,'SYNCHROB_ETA','_'))then
        amus0=pi2*getva(exist)
        if(.not. exist)then
          go to 7310
        endif
        amus1=pi2*getva(exist)
        if(.not. exist)then
          go to 7310
        endif
        amusstep=pi2*getva(exist)
        if(.not. exist)then
          go to 7310
        endif
        call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $       rlist(ifgamm),.true.)
        mphi2=9
        iparams=ktaloc(59)
        call tclr(codin,6)
        if(ifsize .eq. 0 .and. codplt)then
          ifsize=ktaloc(nlat*21)
          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call temits(ilist(1,ilattp+1),
     $     rlist(iftwis),rlist(ifsize),rlist(ifgamm),
     $     ndim,ntwissfun,
     $     mphi2,
     $     amus0,amus1,amusstep,
     $     emxe,emye,rese,rlist(iparams),
     $     lfni,lfno,int8(0),irtc)
        call tfree(iparams)
        if(.not. codplt)then
          call tfgeo(ilist(1,ilattp+1),rlist(ifgeo),rlist(ifpos),
     $         rlist(ifgamm),.true.)
        endif
        emx=max(4.d-13/p0,rgetgl1('EMITX'))
        emy=max(4.d-13/p0,rgetgl1('EMITY'))
c        dpmax=rgetgl1('SIGE')
c        rlist(itlookup('DP',ivtype))=dpmax
        gauss=.true.
        go to 10
 7310   call termes(lfno,
     $       'Usage: SYNCHROB_ETA nus_start nus_stop nus_step.',' ')
        go to 2
      elseif(abbrev(word,'BEAM_SIZE','_'))then
        if(ifsize .eq. 0)then
          ifsize=ktaloc(nlat*21)
          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call tfsetparam
        call tfsize(rlist(iftwis),rlist(ifgamm),rlist(ifsize))
        go to 10
      elseif(abbrev(word,'ALI_GN','_'))then
        call talign(ilist(1,ilattp+1),word,wordp,
     $       rlist(ifpos),ilist(1,ifmult),lfno,exist)
        go to 30
      endif
 7000 call tfgetv(word,ilist(1,ilattp+1),
     $     ilist(1,ifival),ilist(1,ifklp),ilist(1,ifele),
     $     ilist(1,ifele1),rlist(ifcoup),
     $     rlist(iferrk),rlist(ifaux),
     $     ilist(1,iftouchele),ilist(1,iftouchv),flv%ntouch,lfno,exist)
      if(.not. exist)then
        if(itt .ge. 0)then
          call cssetp(nextt)
          call tfprintout(lfno,itfgetrecl(),16,irtc)
          go to 10
        else
          call termes(lfno,
     1         '?Undefined command or element: ',word)
        endif
        if(ios .ne. 0)then
          go to 10
        endif
        go to 2
      endif
      go to 10
30    if(.not. exist)then
        go to 12
      else
        go to 10
      endif
31    if(.not. exist)then
        go to 2
      else
        go to 10
      endif
 32   if(err)then
        go to 2
      else
        go to 10
      endif
 1000 continue
      if(busy)then
        call termes(lfno,'?Recursive CAL or GO',' ')
        go to 2
c        irtcffs=itfmessage(9,'FFS::busy','""')
c        go to 8900
      endif
      busy=.true.
      call tfsetparam
      nfc0=flv%nfc
      if(cell)then
        call tfsetupcell(flv,nlat,maxcond)
c        write(*,*)'tffsa ',flv%nfc,nfc0
      endif
      geomet=.false.
      do i=1,flv%ncalc
        if(flv%icalc(3,i) .gt. mfitleng .and.
     $       flv%icalc(3,i) .le. mfitchi3)then
          geomet=.true.
        endif
      enddo
      nfcol=0
      do k1=1,flv%nfc
        if(flv%kfit(k1) .le. mfit .and. flv%mfitp(k1) .ne. 0)then
          if(flv%kfit(k1) .gt. mfitleng .and.
     $         flv%kfit(k1) .le. mfitchi3)then
            geomet=.true.
          endif
          nfcol=nfcol+1
          if(nfcol .gt. maxcond)then
            call termes(lfno,'?Too many fit conditions.',' ')
            go to 8810
          endif
          flv%kfitp(nfcol)=k1
        endif
      enddo
      if(.not. geomet)then
        do i=1,nele
          ii=(i-1)/2
          if(idtype(ilist(1,ilattp+ilist(i-ii*2,ifklp+ii)))
     $         .eq. 20)then
            geomet=.true.
            exit
          endif
        enddo
      endif
      call tffssetupcouple(lfno)
      if(expndc)then
        call tffsadjustvar
      else
        if(.not. expndc)then
          call termes(lfno,
     $         'Info-Element values are not expanded.',' ')
        endif
      endif
      if(fitflg)then
        call tmov(rlist(ifvalvar),rlist(ifvalvar+nve),flv%nvar)
      endif
      if(tffsinitialcond(iuid,uini,nfr,ndimmax,lfno,err))then
        inicond=.true.
        nfam=nfr
        if(iuid(-nfam) .lt. 0)then
          nfam1=1-nfam
        else
          nfam1=-nfam
        endif
        uini(mfitddp,0)=0.d0
        do i=nfam1,nfr
          kfam(i)=0
          jfam(i)=9999
          dp(i)=uini(mfitddp,i)
        enddo
        jfam(0)=0
      elseif(err)then
        call termes(lfno,'Error in InitialOrbits',' ')
        go to 8810
      else
        inicond=.false.
        nfr=-1
        do k=1,flv%nfc
          if(flv%kfit(k) .le. mfit .and. flv%mfitp(k) .ne. 0)then
            nfr=max((abs(flv%mfitp(k))-1)/2,nfr)
          endif
        enddo
        if(nfr .gt. 0 .and. dpmax .eq. 0.d0)then
          call termes(lfno,'?No momentum bandwidth.',' ')
          go to 8810
        endif
        if(nfr .lt. 0)then
          nfr=0
        endif
        dpm2=rfromd(kxsymbolv('DPM',3))
c        dpm2=rlist(ktlookup('DPM'))
        do i=-nfr,nfr
          dp(i)=dpmax*dble(i)/max(nfr,1)
          if(dpm2 .gt. 0.d0)then
            if(i .eq. -1)then
              dp(-1)=-dpm2
            elseif(i .eq. 1)then
              dp( 1)=+dpm2
            endif
          endif
        enddo
        em=abs(emx)+abs(emy)
        call tffamsetup(dp,dfam,jfam,kfam,ilist(2,ilattp+1),
     $       nfam,nfr,ndimmax,em)
        if(nfam .gt. nfr .and. kfam(-nfam) .le. 0)then
          nfam1=1-nfam
        else
          nfam1=-nfam
        endif
      endif
      wake=(twake .or. lwake) .and. trpt
      kwakep=0
      kwakeelm=0
      nwakep=0
      if(wake)then
        call tffssetupwake(ilist(1,ilattp+1),ilist(1,ifmult),
     $       kwakeelm,kwakep,nwakep,lfno,irtc)
        if(irtc .ne. 0)then
          call termes(lfno,'?Error in WakeFunction.',' ')
          go to 8810
        endif
        if(nwakep .eq. 0)then
          wake=.false.
        endif
      endif
      flv%iut=iutwiss(flv,ilist(1,iftwissp),ilist(1,ilattp+1),
     $     ilist(1,kwakeelm),nwakep,ilist(1,ifvarele),
     $     ilist(1,ifele),int(ifele1),int(ifele2),
     $     nlat,flv%nvar,nfcol,nfam,nut,.not. cell)
      if(flv%iut .le. 1)then
        call termes(lfno,
     $       '?Too many off-momentum or fit points.',' ')
        go to 8810
      endif
c      write(*,*)'tffsa-1 ',flv%iut,%LOC(rlist(flv%iut)),busy
      call tffsmatch(flv,rlist(iftwis),rlist(ifpos),rlist(ifgeo),
     $     rlist(ifgamm),rlist(flv%iut),
     $     ilist(1,ilattp+1),ilist(1,ifmult),ilist(1,ifele),
     $     ilist(1,ifele1),ilist(1,ifele2),rlist(ifcoup),
     $     ilist(1,iftwissp),
     $     ilist(1,ifival),ilist(1,ifklp),
     $     ilist(1,ifvarele),ilist(1,ifvvar),
     $     ilist(1,ifivcomp),rlist(ifvalvar),
     $     rlist(iferrk),rlist(ifaux),nlat,nele,ndim,nut,flv%nvar,
     $     dp(-nfam),tracex(-nfam),tracey(-nfam),
     $     hstab(-nfam),vstab(-nfam),df,nfr,nfam,nfam1,
     $     dfam(1,-nfam),jfam(-nfam),kfam(-nfam),
     $     inicond,iuid(-nfam),uini(1,-nfam),
     $     wake,ilist(1,kwakeelm),klist(kwakep),nwakep,
     $     nqcol,nqcol1,nfcol,nfc0,
     $     maxcond,nlist,brho,
     $     emx,emy,dpmax,dp0,coumin,r,residual(-nfam),absweit,
     $     cell,fitflg,geomet,cellstab,convgo,nparallel,orbitcal,
     $     lfno,irtc)
      call tclrfpe
      if(wake)then
        call tffsclearwake(kwakeelm,kwakep,nwakep)
      endif
      if(irtc .ne. 0)then
        call tmunmapp(flv%iut)
        go to 8810
      endif
      flv%nfc=nfc0
c      write(*,*)'tffsa-2 ',flv%iut,%LOC(rlist(flv%iut)),busy
      call tfshow(flv,ilist(1,ilattp+1),ilist(1,ifmult),
     $     nfr,nfam,nfam1,kfam(-nfam),
     1     scale,nlist,rlist(flv%iut),ilist(1,iftwissp),nut,
     $     rlist(ifgeo),rlist(ifpos),tracex(-nfam),tracey(-nfam),
     $     cellstab,residual(-nfam),hstab(-nfam),vstab(-nfam),
     $     dp(-nfam),df,inicond,iuid(-nfam),
     1     mfpnt,mfpnt1,nqcol,nqcol1,
     $     kffs,irtcffs,lfnb .gt. 1,lfno)
      call tmunmapp(flv%iut)
      call tffsclearcouple(ilist(1,ifele2))
      if(cell)then
        str=' '
        do i=nfam1,nfam
          if(.not. hstab(i))then
            str='Horizontal'
            exit
          endif
        enddo
        do i=nfam1,nfam
          if(.not. vstab(i))then
            if(str(1:1) .ne. ' ')then
              str='Horizontal/Vertical'
            else
              str='Vertical'
            endif
            exit
          endif
        enddo
        if(str(1:1) .ne. ' ')then
          write(lfno,*)str(1:len_trim(str)),' Unstable.'
        endif
      endif
      busy=.false.
      go to 12
 8810 busy=.false.
      flv%nfc=nfc0
      go to 10
 8900 call tffsfree
 9000 call tfconnectk(kffs,irtcffs)
      return
      end

      subroutine mcmess(lfno)
      implicit none
      integer*4 lfno
      write(lfno,*)
     $     'Orbit Correction commands will be removed soon.',
     $     'Use orbit correction functions instead.',
     $     'Please send to Katsunobu.Oide@kek.jp, ',
     $     'if you still need those commands.'
      return
      end

      integer*8 function iutwiss(flv,
     $     itwissp,latt,iwakeelm,nwakep,
     $     ivarele,iele,iele1,iele2,
     $     nlat,nvar,nfcol,nfam,nut,nonl)
      use tfstk
      use ffslocal, only:ffslocalv
      use tffitcode
      implicit none
      type (ffslocalv) flv
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nlat,nvar,nfcol,nfam,nut,nwakep,i2
      integer*4 itwissp(nlat),latt(2,nlat),ivarele(nvar),iele(nlat),
     $     iele1,iele2,iwakeelm(nwakep)
      integer*4 i,j,id,k
      integer*8 itmmapp
      real*8 tffsmarkoffset
      logical*4 nonl
      do i=3,nlat-1
        itwissp(i)=0
      enddo
      if(idtype(latt(1,nlat-1)) .eq. icMARK)then
        itwissp(nlat-1)=1
      endif
      if(nonl)then
        do i=2,nlat-1
          id=idtype(latt(1,i))
          if(id .ge. icSEXT .and. id .le. icDODECA
     $         .or. id .eq. icMULT .or. id .eq. icSOL)then
            itwissp(i)=1
            itwissp(i+1)=1
          endif
        enddo
      endif
      LOOP_I: do i=2,nlat-1
        do j=1,nvar
          if(ivarele(j) .eq. ilist(iele(i),iele1))then
            itwissp(i)=1
            itwissp(i+1)=1
            cycle LOOP_I
          endif
        enddo
        if(ilist(i,iele2) .ne. 0)then
          itwissp(i)=1
          itwissp(i+1)=1
        endif
      enddo LOOP_I
      do i=1,nwakep
        j=iwakeelm(i)
        itwissp(j)=1
        itwissp(j+1)=1
      enddo
      do k=1,nfcol
        j=flv%kfitp(k)
        if(flv%mfitp(j) .ne. 0 .and. flv%ifitp(j) .gt. 0)then
          itwissp(flv%ifitp(j))=1
          if(flv%ifitp1(j) .ne. flv%ifitp(j))then
            itwissp(flv%ifitp1(j))=1
            if(flv%mfitp(j) .lt. 0 .and. flv%kfit(j) .le. mfitdz
     $           .and. flv%kfit(j) .ne. mfitnx
     $           .and. flv%kfit(j) .ne. mfitny)then
              do i=min(flv%ifitp(j),flv%ifitp1(j)),
     $             max(flv%ifitp(j),flv%ifitp1(j))
                itwissp(i)=1
              enddo
            endif
          endif
        endif
      enddo
      do j=1,flv%ncalc
        if(flv%icalc(1,j) .gt. 0)then
          itwissp(flv%icalc(1,j))=1
          if(flv%icalc(2,j) .gt. 0)then
            itwissp(flv%icalc(2,j))=1
          endif
        endif
      enddo
      itwissp(1)=1
      i2=max(1,min(nlat,1+int(tffsmarkoffset(latt(1,1)))))
      itwissp(i2)=1
      itwissp(min(nlat,i2+1))=1
      itwissp(nlat)=1
      nut=1
      do i=2,nlat
        if(itwissp(i) .ne. 0)then
          nut=nut+1
          itwissp(i)=nut
        endif
      enddo
      iutwiss=itmmapp(nut*(2*nfam+1)*ntwissfun)
      return
      end

      subroutine tffssaveparams(icmd,ilattp,n,err)
      use tfstk
      use ffslocal
      use ffs, local_ilattp=>ilattp
      implicit none
      integer*4 nstk
      parameter (nstk=64)
c nlocal = mcommon in TFFSLOCAL.inc
      integer*8 ktaloc
      integer*4 icmd,ilattp,n,nr,istkp,i,nxh
      integer*8 iffsstk(nstk)
      logical*4 err
      save istkp,iffsstk
      data istkp /0/
      nxh=int((sizeof(ffv)+7)/8)
      nr=n
      err=.false.
      if(icmd .eq. 0)then
        istkp=istkp+1
        if(istkp .gt. nstk)then
          istkp=nstk
          err=.true.
          return
        endif
        iffsstk(istkp)=ktaloc(nr+nxh+1)
        ilist(1,iffsstk(istkp))=ilattp
        call tmov(ffv,rlist(iffsstk(istkp)+1),nxh)
        call tmov(flv,rlist(iffsstk(istkp)+nxh+1),nr)
        call tffsresetbzl
      elseif(icmd .eq. 1)then
        if(istkp .gt. 0)then
          ilattp=ilist(1,iffsstk(istkp))
          call tmov(rlist(iffsstk(istkp)+1),ffv,nxh)
          call tmov(rlist(iffsstk(istkp)+nxh+1),flv,nr)
          call tfree(iffsstk(istkp))
          istkp=istkp-1
        else
          err=.true.
        endif
      elseif(icmd .eq. 2)then
        err=.false.
        do i=istkp,1,-1
          if(ilist(1,iffsstk(i)) .eq. ilattp)then
            err=.true.
            return
          endif
        enddo
      elseif(icmd .eq. 3)then
        err=istkp .ge. nstk
        return
      elseif(icmd .eq. 4)then
        err=istkp .le. 0
        return
      else
        do i=istkp,n,-1
          call tmov(rlist(iffsstk(i)+1),ffv,nxh)
          call tffsfree
          call tffsfreebzl
          call tfree(iffsstk(i))
        enddo
        istkp=min(n-1,istkp)
      endif
      return
      end

      logical*4 function tfvcomp(ivcomp,
     $     ivarele,ivvar,ival,nele,nvar)
      implicit none
      integer*4 nvar,i,nele,
     $     ivcomp(nvar),ivarele(nvar),ivvar(nvar),ival(nele)
      tfvcomp=.false.
      do i=1,nvar
        if(ivcomp(i) .ne. 0 .and.
     $       ivvar(i) .ne. ival(ivarele(i)))then
          tfvcomp=.true.
          return
        endif
      enddo
      return
      end

      subroutine tffsadjustvar
      use tfstk
      use ffs
      use ffslocal
      use tffitcode
      implicit none
      call tffsadjust(ilist(1,iftouchele),ilist(1,iftouchv),
     $     ilist(1,ilattp+1),rlist(iferrk),rlist(ifcoup),
     $     ilist(1,ifele),ilist(1,ifele1),
     $     ilist(1,ifklp),ilist(1,ifival),flv%ntouch)
      call tfinitvar(ilist(1,ifvarele),ilist(1,ifvvar),
     $     ilist(1,ifivcomp),
     $     rlist(ifvalvar),flv%nvar,ilist(1,ilattp+1),ilist(1,ifklp),
     $     ilist(1,ifival),rlist(iferrk),nlat,nele)
      return
      end

      logical*4 function tffsinitialcond(
     $     iuid,uini,nfr,ndimmax,lfno,err)
      use tfstk
      use tffitcode
      implicit none
      type (sad_list), pointer :: klx,klj
      integer*8 kx
      integer*4 nfr,irtc,lfno,n,i,j,nfr1,
     $     ndimmax,iuid(-ndimmax:ndimmax)
      real*8 uini(28,-ndimmax:ndimmax)
      logical*4 err
      integer*8 iaini
      data iaini/0/
      tffsinitialcond=.false.
      err=.true.
      if(iaini .eq. 0)then
        iaini=ktfsymbolz('InitialOrbits',13)
      endif
      call tfsyeval(iaini,kx,irtc)
      if(irtc .ne. 0)then
        call tfemes(irtc,'InitialOrbits',1,lfno)
        return
      endif
      if(tflistqk(kx,klx))then
        n=klx%nl
        if(n .le. 0)then
          return
        endif
        nfr=(n+1)/2
        if(nfr .gt. ndimmax)then
          call termes(lfno,'Too many initial conditions.',' ')
          return
        endif
        nfr1=-(n/2)
        iuid(-nfr)=-1
        j=0
        do i=nfr1,nfr
          if(i .eq. 0)then
            call tclr(uini(1,0),28)
            iuid(0)=0
          else
            j=j+1
            if(tfreallistq(klx%body(j),klj))then
              if(klj%nl .eq. 6)then
                uini(mfitdx:mfitdx+5,i)=klj%rbody(1:6)
                uini(1,i)=-1.d0
                iuid(i)=j
              elseif(klj%nl .eq. 28)then
                uini(mfitax:mfitax+27,i)=klj%rbody(1:28)
                iuid(i)=j
              else
                return
              endif
            else
              return
            endif
          endif
        enddo
        tffsinitialcond=.true.
        err=.false.
      elseif(kx .eq. ktfoper+mtfnull .or. ktfsymbolq(kx))then
        err=.false.
        return
      endif
      return
      end

      subroutine tffsclearcouple(iele2)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 iele2(nlat),i
      do i=1,nlat-1
        if(iele2(i) .ne. 0)then
          call tfree(int8(iele2(i)))
        endif
      enddo
      return
      end

      subroutine tffssetupcouple(lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_list), pointer :: klx,kli,kle
      integer*8 ktaloc
      integer*4 lfno,i,j,k,nk,m,me,nc, ie,iet,ik,irtc
      character*(MAXPNAME) key,tfgetstrs
      type (sad_descriptor) itfelv,itfcoupk,kx,ki,kk,ke
      data itfelv%k,itfcoupk%k /0,0/
      ilist(1:nlat,ifele2)=0
      if(itfelv%k .eq. 0)then
        itfelv=kxsymbolz('`ElementValues',14)
        itfcoupk=kxsymbolz('`CoupledKeys',12)
      endif
      levele=levele+1
      call tfsyeval(itfcoupk,kx,irtc)
      call tfconnect(kx,irtc)
      if(irtc .ne. 0)then
        go to 9010
      endif
c      call tfdebugprint(kx,'setupcoup',3)
      if(.not. tflistqd(kx,klx))then
        go to 9000
      endif
      m=klx%nl
      do i=1,m
        ki=klx%dbody(i)
        if(.not. tflistqd(ki,kli))then
          go to 9000
        endif
        kk=kli%dbody(1)
        if(.not. ktfstringqd(kk))then
          go to 9000
        endif
        key=tfgetstrs(kk,nc)
        ke=kli%dbody(2)
        if(.not. tflistqd(ke,kle))then
          go to 9000
        endif
        me=kle%nl
        if(ktfreallistqo(kle))then
          do j=1,me
            ie=int(kle%rbody(j))
            call tfkeya(ie,key,ik)
            if(ik .lt. 0)then
              do k=1,nlat-1
                if(ilist(k,ifele1) .eq. -ie)then
                  iet=ilist(k,ifele2)
                  if(iet .eq. 0)then
                    iet=int(ktaloc(m+1))
                    ilist(1,iet)=0
                    ilist(k,ifele2)=iet
                  endif
                  nk=ilist(1,iet)+1
                  ilist(1,iet+nk)=i
                  ilist(2,iet+nk)=-ik
                  ilist(1,iet)=nk
                  ilist(nlat,ifele2)=1
                endif
              enddo
            elseif(ik .gt. 0)then
              iet=ilist(ie,ifele2)
              if(iet .eq. 0)then
                iet=int(ktaloc(m+1))
                ilist(1,iet)=0
                ilist(ie,ifele2)=iet
              endif
              nk=ilist(1,iet)+1
              ilist(1,iet+nk)=i
              ilist(2,iet+nk)=ik
              ilist(1,iet)=nk
              ilist(nlat,ifele2)=1
            endif
          enddo
        endif
      enddo
      return
 9010 if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
 9000 call termes(lfno,
     $     'ElementValues := { key[elem] :> f[ key1[elem1] ], ...}',
     $     ' ')
      return
      end

      subroutine tfkeya(i,keyword,ia)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,ia,it,kl,l
      character*(*) keyword
      character*8 tfkwrd,kw
      if(i .gt. 0)then
        it=idtype(ilist(1,ilattp+i))
      else
        kl=ilist(-i,ifklp)
        it=idtype(ilist(1,ilattp+kl))
      endif
      kw='-'
      l=0
      do while(kw .ne. ' ' .and. kw .ne. keyword)
        l=l+1
        kw=tfkwrd(it,l)
      enddo
      if(kw .eq. ' ')then
        ia=0
      else
        ia=sign(l,i)
      endif
      return
      end

      subroutine tfsetupcell(flv,nlat,maxcond)
      use ffslocal, only: ffslocalv
      use tffitcode
      implicit none
      type (ffslocalv) flv
      integer*4 maxcond,nlat,mfc
      integer*4 k
      if(flv%nfc .gt. maxcond-4)then
        return
      endif
      mfc=1
      do k=1,flv%nfc
        if(flv%kfit(k) .le. mfit)then
          mfc=max(abs(flv%mfitp(k)),mfc)
        endif
      enddo
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitax
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=0.d0
      flv%mfitp(flv%nfc)=mfc
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitay
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=0.d0
      flv%mfitp(flv%nfc)=mfc
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitbx
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=1.d0
      flv%mfitp(flv%nfc)=mfc
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitby
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=1.d0
      flv%mfitp(flv%nfc)=mfc
      return
      end

      subroutine tfsize(twiss,gammab,size)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),
     $     size(21,nlat)
      do 10 i=1,nlat
        call tfbeam(twiss,gammab,i,0.d0,size(1,i))
10    continue
      return
      end

      subroutine tscale(nlist,scale,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 lfno,i
      character*(*) nlist(mfit)
      real*8 scale(mfit)
      write(lfno,9001)(nlist(i),scale(i),i=1,mfit)
9001  format(:3(a,1pg15.7,1x))
      return
      end

      subroutine tclrline(line)
      use tfstk
      implicit none
      integer*4 line,i,ip,n
      do i=1,ilist(1,line)
        ip=ilist(2,line+i)
        n=ilist(1,ip)+1
        call tfreem(ip,n)
      enddo
      n=ilist(1,line)+1
      call tfreem(line,n)
      return
      end

      subroutine tfblocksym(str,nch)
      use tfstk
      implicit none
      character*(*) str
      integer*4 nch
      type (sad_descriptor) ks,kx
      type (sad_symbol), pointer :: sym
      ks=kxsymbolf(str,nch,.false.)
      call descr_sym(ks,sym)
      kx=kxnaloc1(max(0,sym%gen),sym%loc)
      return
      end

      subroutine tfunblocksym(str,nch,del)
      use tfstk
      implicit none
      character*(*) str
      integer*4 nch
      logical*4 del
      type (sad_descriptor) ks
      type (sad_symdef), pointer :: symd
      ks=kxsymbolf(str,nch,.false.)
      call descr_symdef(ks,symd)
      call tfdelete(symd,del,.true.)
      return
      end
