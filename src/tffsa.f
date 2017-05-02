      module track_tt
        integer*8 itt1,itt2,itt3,itt4,itt5,itt6
      end module

      subroutine tffsa(lfnb,kffs,irtcffs)
      use tfstk
      use ffs
      use ffs_pointer
      use trackbypass
      use tffitcode
      use ffs_fit
      use ffs_wake
      use sad_main
      use tfcsi
      use track_tt
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 maxrpt,maxlfn,hsrchz
      integer*8 kffs,k,kx,itwisso,
     $     ifvalvar2,iparams,kax,iutwiss
      integer*4 kk,i,lfnb,ia,iflevel,j,ielm,igelm,k1,
     $     ii,irtc0,it,itemon,itmon,itestr,itstr,itt,
     $     iuse,l,itfuplevel,
     $     levelr,lfnl0,lpw,meas0,mfpnta,igetgl1,lenw,
     $     mphi2,newcor,next,nextt,nfp,nmon,
     $     nster,nrpt1,itfpeeko,itfgetrecl,nl
      real*8 rmax,amus0,amus1,amusstep,apert,axi,ayi,ctime1,
     $     dpm2,dpxi,dpyi,em,emxe,emye,epxi,epyi,pspan,r,r2i,r3i,
     $     trval,rese,v,wa,wd,wl,xa,ya,xxa,xya,yya,getva,rgetgl1,
     $     wp,getvad,tgetgcut
      parameter (rmax=1.d35)
      parameter (maxrpt=32,maxlfn=128)
      character*256 word,wordp,title,case,tfgetstrv,tfgetstrs,tfgetstr
      character*(MAXPNAME) ename
      character*(MAXPNAME+8) name
      character*16 autofg
      character*20 str
      integer*4 irtcffs,irtc,nc
      integer*4 lfnstk(maxlfn),lfret(maxlfn),
     $     lfrecl(maxlfn),lflinep(maxlfn),nrpt(maxrpt),
     $     irptp(maxrpt),df(maxcond)
      real*8 chi0(3),trdtbl(3,6),rfromk
      logical*4 err,new,cmd,open98,abbrev,ftest,
     $     frefix,exist,init,trpt0,expnd,chguse,visit,
     $     byeall,expndc,tfvcomp,tffsinitialcond,
     $     geocal0,busy
      logical*4 lfopen(maxlfn)
      save lfnstk,lfret,lfrecl,lfopen,open98,exist,init,
     $     nmon,nster
      save busy
      data busy /.false./
      itwisso(kk,i,j)=iftwis+kk+nlat*(i+ndim+(j-1)*ndima)-1
      flv%mcommon=int((sizeof(flv)+7)/8)
c      write(*,*)'tffsa ',flv%mcommon
      kffs=ktfoper+mtfnull
      irtcffs=0
c      write(*,*)'tffsa-0 ',lfnb,ilattp,mstk
      l=itfuplevel()
c      write(*,*)'tffsa-1 ',l
      chguse=.false.
c     begin initialize for preventing compiler warning
      levelr=0
c     end   initialize for preventing compiler warning
 101  if(lfnb .le. 1 .or. chguse)then
        call tffsalloc()
        if(.not. chguse)then
          call cputime(flv%ctime0,irtc0)
          flv%ctime2=flv%ctime0
        endif
        flv%iut=0
        flv%nvar=0
        flv%ntouch=0
        trdtbl=0.d0
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
          call tfgeo(.true.)
          geocal=geocal0
        endif
        if(lfnb .le. 0)then
          go to 8900
        endif
        call tffsinitparam
c     kikuchi ... next 1 line added     (11/13/'91)
        call corinit(newcor,nster,nmon,itstr,itestr,itmon,itemon)
c     
        flv%measp=nlat
        mfpnt=nlat
        mfpnt1=nlat
        flv%nfc=0
        call tfinitcalc
        call tmast
        call twmov(1,twiss,nlat,ndim,.true.)
        if(.not. chguse)then
          scale=1.d0
          scale(mfitnx)=pi2
          scale(mfitny)=pi2
          scale(mfitnz)=pi2
          scale(mfitchi1)=pi/180.d0
          scale(mfitchi2)=pi/180.d0
          scale(mfitchi3)=pi/180.d0
          xixf=0.d0
          xiyf=0.d0
          flv%rsconv=1.d-9
          convgo=.false.
          trsize=.false.
          cellstab=.true.
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
     1     (word,flv%nfc,lfno,flv%icalc,flv%ncalc,
     1     flv%kfit,flv%fitval,flv%mfitp,flv%ifitp,flv%ifitp1,
     $     exist,err)
      if(err)then
        go to 2
      endif
      if(exist)then
        go to 10
      endif
      call terror(word,new,lfno,exist,err)
      if(err)then
        go to 2
      endif
      if(exist)then
        go to 10
      elseif(new)then
        go to 12
      endif
      call tffile(word,lfnstk,lfopen,lfret,lfnb,lfrecl,
     $     maxlfn,init,exist)
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
        call tfsave(word,.false.,flv%ntouch)
        call tfgeo(.false.)
        call corfree(newcor,nster,nmon,itstr,itestr,itmon,itemon)
        go to 8900
      elseif(word .eq. 'QUIT')then
        call tfgeo(.false.)
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
          call tffssaveparams(2,ilattp,err)
          expnd=expnd .and. .not. err
          call tclrpara(elatt,nlat-1)
          if(visit)then
            call tffssaveparams(0,ilattp,err)
            call tfblocksym('`FitFunction',12)
            call tfblocksym('`FitValue',9)
            call tfblocksym('`ElementValues',14)
          else
            call tfunblocksym('`FitFunction',12,.false.)
            call tfunblocksym('`FitValue',9,.false.)
            call tfunblocksym('`ElementValues',14,.false.)
            call tffsfree
            if(lattuse .eq. lattredef)then
              call tclrline(lattredef)
              lattredef=0
            endif
          endif
          if(ilist(2,idval(iuse)) .le. 0 .or. expnd)then
            call expnln(iuse)
          endif
          call filaux(iuse)
          ilattp=idval(ilist(2,idval(iuse)))
          call loc_el(ilattp,elatt)
          lattuse=ilattp
          call tclrpara(elatt,elatt%nlat1-2)
          dleng =rlist(elatt%aux+1)*rgetgl1('FSHIFT')
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
        call tffssaveparams(4,int8(0),err)
        if(err)then
          if(byeall)then
            go to 10
          else
            call termes(lfno,'BYE without VISIT.',' ')
            go to 2
          endif
        endif
        call tclrpara(elatt,elatt%nlat1-2)
        call tffsfree
        if(byeall)then
          call tffssaveparams(-2,int8(0),err)
        endif
        call tffssaveparams(1,ilattp,err)
        call loc_el(ilattp,elatt)
        lattuse=ilattp
        nlat=elatt%nlat1-1
        dleng =rlist(elatt%aux+1)*rgetgl1('FSHIFT')
        call ffs_init_pointer
        call ffs_twiss_pointer
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
        call tfattr(word,lfno,exist,kffs,irtcffs,lfnb .gt. 1)
        go to 30
      elseif(word .eq. 'SAVE')then
        call tffsadjustvar
        call tfsave(word,.true.,flv%ntouch)
      elseif(abbrev(word,'EXPAND','_'))then
        call tffsadjustvar
        call tfsave(word,.true.,flv%ntouch)
      elseif(word .eq. 'RESET')then
        call tfrst(word,.true.,flv%nvar,flv%ntouch)
      elseif(word .eq. 'FREE' .or. word .eq. 'FIX')then
        frefix=word .eq. 'FREE'
        call peekwd(word,next)
        if(word .eq. ' ')then
ckikuchi ... next 5 lines added     (8/17/'90)
          call mcmess(lfno)
          if(frefix) then
            call mcfre(newcor)
          else
            call mcfix(latt,twiss,
     $           gammab,newcor)
          endif
          go to 10
        endif
        call tffsfreefix(frefix,flv%nvar,flv%ntouch,lfno)
      elseif(word .eq. 'FIT')then
        call getwdl2(word,wordp)
        if(word .eq. 'ALL')then
          mfpnt=1
          mfpnt1=nlat
        else
          mfpnt=ielm(wordp,exist)
          if(.not. exist)then
            mfpnt1=mfpnt
            go to 12
          endif
          call getwdl2(word,wordp)
          mfpnta=ielm(wordp,exist)
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
          expndc=.not. tfvcomp()
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
        expndc=.not. tfvcomp()
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
        call tfsetv(flv%nvar)
      elseif(abbrev(word,'T_YPE','_'))then
        call tfsetparam
        call tftype(lfno,word)
        go to 12
      elseif(abbrev(word,'DISP_LAY','_'))then
        call tfsetparam
        call tfdisp(word,id1,id2,dp0,lfno,exist)
        go to 30
      elseif(word .eq. 'SCALE')then
        call tscale(nlist,scale,lfno)
      elseif(word.eq.'MAPANA') then
        call gosadpls(latt,ilist(1,ifklp),
     $       ilist(1,ifele),lfno)
      elseif(abbrev(word,'AP_ERTURE','_'))then
        call tfsetparam
        call tfaprt(lfno,word)
        go to 12
      elseif(abbrev(word,'MAT_RIX','_'))then
        call tfmat(lfno,word,wordp,exist)
        go to 30
      elseif(abbrev(word,'CHRO_MATICITY','_'))then
        call tfsetparam
        i=0
        call tfchro(latt,
     1       rlist(itwisso(1,i,1)),rlist(itwisso(1,i,2)),
     $       rlist(itwisso(1,i,3)),rlist(itwisso(1,i,7)),
     1       rlist(itwisso(1,i,4)),rlist(itwisso(1,i,5)),
     $       rlist(itwisso(1,i,6)),lfno)
      elseif(word .eq. 'RADINT')then
        call tfsetparam
        call intgrl( latt,twiss,0,1.d0,lfno)
      elseif(word .eq. 'DIMAD')then
        call tfsetparam
        call tdimad(latt,mult,lfno)
      elseif(word .eq. 'ZAP')then
        call tfsetparam
        apert=getva(exist)
        if(.not. exist)then
          apert=1.d0
        endif
        i=0
        call tfzap(apert,pos,
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
            meas0=ielm(wordp,exist)
            if(exist)then
               flv%measp=meas0
            endif
            title=tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodf(latt,flv%measp,name,
     $           mult,
     &           .false.,word,title,case,exist,
     &           twiss,lfno)
            if(.not. exist)then
               go to 12
            endif
         else if(abbrev(word,'MVCOD','_')) then
            call getwdl2(word,wordp)
            meas0=ielm(wordp,exist)
            if(.NOT. exist)then
               meas0=flv%measp
            endif
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodM(latt,meas0,name,mult,
     &           .false.,word,title,case,exist,
     &           twiss,lfno)
            if(.not. exist)then
               go to 12
            endif
         else if(abbrev(word,'STN_COD','_')) then
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcodS(latt,meas0,name,mult,
     &           .false.,word,title,case,exist,
     &           twiss,lfno)
            go to 10
         else if(abbrev(word,'ALIGN','_')) then
            call getwdl(word)
            call nalign(latt,mult,
     $           ilist(1,ifmast),word,lfno,exist)
            go to 10
         else
            meas0=ielm(wordp,exist)
            if(exist)then
               flv%measp=meas0
            endif
            call elname(flv%measp,name)
            title=Tfgetstrv('TITLE')
            case=Tfgetstrv('CASE')
            call trcoda(latt,flv%measp,name,
     $           mult,.false.,word,title,case,exist,
     &           twiss,xa,ya,xxa,xya,yya,lfno)
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
     &             latt,ilist(1,ifele),rlist(ifcoup),
     $       rlist(ifgeo),pos,ilist(1,ifmast))
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
        call tfchgv(lfno)
        call tfinitvar(flv%nvar)
        go to 10
      elseif(abbrev(word,'REJ_ECT','_'))then
        call tfrej
     1   (word,nlist,flv%nfc,mfpnt,mfpnt1,
     1    flv%kfit,flv%mfitp,flv%ifitp,flv%ifitp1,exist)
        go to 30
      elseif(abbrev(word,'COUP_LE','_'))then
        call tfcoup(lfno,exist)
        call tffsadjustvar
        if(.not. exist)then
          go to 2
        endif
        go to 10
      elseif(abbrev(word,'DECOUP_LE','_'))then
        call tfdecoup(lfno)
        call tffsadjustvar
        go to 10
      elseif(abbrev(word,'INDEP_ENDENT','_'))then
        call tfindep
        call tffsadjustvar
        go to 10
      elseif(word .eq. 'SBUNCH')then
        call termes(lfno,'?SBUNCH is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        call tfgetr(rlist(iwakepold+1),1.d0,word,lfno,exist)
c        go to 31
      elseif(word .eq. 'DBUNCH')then
        call termes(lfno,'?DBUNCH is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        call tfdbun(word,rlist(ilist(2,iwakepold)),ilist(1,iwakepold),
c     1              lfno,err)
c        go to 32
      elseif(word .eq. 'SLICE')then
        call termes(lfno,'?SLICE is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        ns=ilist(1,iwakepold+2)
c        call tfgeti(ns,1.d0,word,lfno,exist)
c        if(ns .le. 0)then
c          call termes(lfno,'?Parameter out of range in SLICE.',' ')
c          go to 2
c        endif
c        ilist(1,iwakepold+2)=ns
c        go to 31
      elseif(abbrev(word,'ORI_GIN','_') .or. word .eq. 'ORG')then
        iorgr=igelm(word,exist)
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
        call tshow(kffs,irtcffs,lfnb .gt. 1,lfno)
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
        do i=1,nflag
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
        call twelm(lfno,mfpnt,mfpnt1,'FIT',lpw,8)
        call twelm(lfno,flv%measp,0,'MEA_SURE',lpw,8)
        call twelm(lfno,iorgr,0,'ORG',lpw,8)
        call twelm(lfno,id1,id2,'DISP_LAY',lpw,8)
        call twbuf(word,lfno,1,lpw,8,-1)
        go to 10
      elseif(abbrev(word,'VAR_IABLES','_') .or. word .eq. 'VARS')then
        call tfinitvar(flv%nvar)
        call tfvars(flv%nvar,kffs,irtcffs,lfnb .gt. 1,lfno)
        go to 10
      elseif(abbrev(word,'RENUM_BER','_'))then
        call tffsrenumber(lfno)
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
          call twsdrw(latt,pos,ilist(1,ifele),
     $         word,wordp,lfno,
     1         twiss,0,rlist(itmon),rlist(itemon),nmon,
     1         title,case,exist)
        else
          word=tfgetstr(kx,nc)
          exist=nc .eq. 0
c          if(exist)then
c            call cssetp(icsmrk()+1)
c          endif
c          write(*,*)'DRAW ',exist,nc,ios,' ',word(1:nc)
        endif
        go to 30
      elseif(word .eq. 'GEO')then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call geodrw(rlist(ifgeo),word,lfno,title,case)
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
        meas0=ielm(wordp,exist)
        if(exist)then
          flv%measp=meas0
        endif
        call elname(flv%measp,name)
ckikuchi ... 1 line modified
        if(twake .or. lwake)then
c          iwakepold=ifwakep
c          ns=ilist(1,iwakepold+2)
c          np0=max(1,np0/ilist(1,iwakepold)/ns/8)*ns*
c     $         ilist(1,iwakepold)*8
          call ffs_init_sizep
        endif
        trpt=.true.
        trf0=0.d0
        vcphic=0.d0
        vcalpha=1.d0
        title=tfgetstrv('TITLE')
        case=tfgetstrv('CASE')
        call trackb(latt,flv%measp,name,.true.,
     $       word,title,case,exist,
     $       kffs,irtcffs,lfnb .gt. 1,
     1       xa,ya,xxa,xya,yya,lfno)
        if(.not. exist)then
          go to 12
        endif
        go to 10
      elseif(abbrev(word,'ORBITJ_ITTER','_'))then
        call tfojit(twiss,lfno)
        go to 10
      elseif(abbrev(word,'PUTM_EASURE','_'))then
         call  nfputm(lfno,xa,ya,xxa,xya,yya)
         goto 10
      elseif(abbrev(word,'GET_MEASURE','_'))then
         call tfgetm(flv%mfitp,xa,ya)
         go to 10
      elseif(abbrev(word,'LOGM_EASURE','_'))then
         write(lfno,9771)xa,ya,xxa,xya,yya
 9771    format(1p,5g15.7)
         go to 10
C   19/12/90 302031849  MEMBER NAME  CORRECT  *.FORT     M  E2FORT
      elseif(abbrev(word,'STE_ERING','_'))then
ckiku   call tfstr(word,latt,ist,nstr)
        call mcmess(lfno)
        call mcstr(word,latt,
     $       mult,ilist(1,ifmast),itstr,itestr,nster,
     $       .false.,lfno)
        go to 12
      elseif(word.eq.'COUPLESTE') then
        call mcoupste(word,latt,mult,
     $       rlist(itstr),nster,lfno)
        go to 12
      elseif(word.eq.'READSTE')then
        call mcstr(' ',latt,mult,
     $       ilist(1,ifmast),itstr,itestr,nster,.true.,lfno)
        call preadstr(word,latt,twiss,
     $       mult,rlist(itstr),nster,lfno)
        go to 12
      elseif(word.eq.'WRITESTE')then
        call pwrtstr(word,latt,twiss,
     $       mult,rlist(itstr),nster,lfno)
        go to 12
      elseif(abbrev(word,'MON_ITOR','_'))then
        call mcmess(lfno)
        call mcmon(word,latt,mult,pos,
     $       ilist(1,ifmast),itmon,itemon,nmon,.false.,
     &             lfno)
        go to 12
      elseif(word.eq.'READMON')then
        call mcmon(' ',latt,mult,pos,
     $       ilist(1,ifmast),itmon,itemon,nmon,.true.,
     &             lfno)
        call preadmon(word,latt,twiss,
     $       mult,rlist(itmon),nmon,lfno)
        go to 12
      elseif(word.eq.'WRITEMON')then
        call pwrtmon(word,latt,twiss,
     $       mult,rlist(itmon),nmon,lfno)
        go to 12
      elseif(word.eq.'KILL')then
        call pkill(word,wordp,latt,mult,
     $       rlist(itmon),nmon,rlist(itstr),
     1       nster,lfno)
        go to 12
      elseif(abbrev(word,'COR_RECT','_')) then
        call mcmess(lfno)
        call mccor(word,wordp,
     $       latt,pos,twiss,
     $       gammab,mult,ilist(1,ifmast),flv%kfit,
     $       flv%ifitp,
     $       flv%mfitp,flv%fitval,flv%nfc,rlist(itstr),rlist(itestr),
     $       nster,
     $       rlist(itmon),rlist(itemon),nmon,newcor,lfno)
        go to 12
      elseif(abbrev(word,'KI_CK','_'))then
        call mckick(word,wordp,latt,mult)
        go to 12
      elseif(abbrev(word,'CLE_AR','_')) then
        call mclear(word,latt,rlist(itstr),nster)
        goto 12
      elseif(abbrev(word,'UNDO_CORRECT','_')) then
        call pundo(latt,twiss,
     $       gammab,0,1.d0+dp0)
        go to 10
      elseif(word.eq.'SUM+'.or.word.eq.'SUM-'.or.word.eq.'SUM'.or.
     1       abbrev(word,'SUMC_LR','_')) then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call pstati(word,latt,twiss,
     $       rlist(itmon),rlist(itemon),nmon,
     1       rlist(itstr),nster,title,case,lfno)
        go to 12
      elseif(abbrev(word,'STO_RE','_')) then
        call mstore(word,latt,twiss,
     $       rlist(itstr),nster,rlist(itmon),
     1       rlist(itemon),nmon)
        go to 12
      elseif(abbrev(word,'RECA_LL','_')) then
        call mrecal(word,latt,twiss,
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
        call mstack(word,latt,twiss,
     $       rlist(itstr),nster,rlist(itmon),
     1       rlist(itemon),nmon)
        goto 12
      elseif(word.eq.'TRIM') then
        call ptrim(word,latt,twiss,rlist(itmon),
     $       rlist(itemon),nmon,lfno)
        goto 12
      elseif(word.eq.'TCOD') then
        call mcrcod(latt,twiss,mult,
     $       rlist(itmon),rlist(itemon),nmon,
     $       .false.,.false.,.true.,lfno)
        go to 10
      elseif(word.eq.'MCOD') then
        call mcrcod(latt,twiss,mult,
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
        call pecorr(word,latt,pos,
     $       rlist(iferrk),title,case,lfno)
        goto 12
      elseif(word .eq. 'BUMP')then
        call getwdl(word)
        call pcbak(latt,twiss)
        if(abbrev(word,'I_TERATION','_')) then
          call getwdl(word)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.true.)
          call pbump(latt,twiss,
     $         gammab,0,mult,rlist(itstr),
     1               rlist(itestr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%fitval,flv%nfc,
     1               .true.,.false.,.false.,lfno)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.false.)
        else
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.true.)
          call pbump(latt,twiss,
     $         gammab,0,mult,rlist(itstr),
     1               rlist(itestr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%fitval,flv%nfc,
     1               .false.,.false.,.false.,lfno)
          call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $         flv%nfc,.false.)
        endif
        go to 12
      elseif(word.eq.'LBUMP') then
        call getwdl(word)
        call pcbak(latt,twiss)
        call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $       flv%nfc,.true.)
        call pbump(latt,twiss,gammab,
     $       0,mult,rlist(itstr),rlist(itestr),
     1             nster,flv%kfit,flv%ifitp,flv%mfitp,flv%fitval,
     $       flv%nfc,.false.,.true.,
     1             .false.,lfno)
        call pbumps(rlist(itstr),nster,flv%kfit,flv%ifitp,flv%mfitp,
     $       flv%nfc,.false.)
        go to 12
      elseif(abbrev(word,'TOR_','_')) then
        call ptol(word,latt,mult,
     $       ilist(1,ifmast),twiss,
     $       gammab,0,1.d0+dp0,
     1            rlist(itstr),nster,lfno)
        goto 12
      elseif(word.eq.'VBUMP') then
        call pvbump(word,wordp,latt,twiss,
     $       mult,ilist(1,ifmast),rlist(itstr),nster,
     $       nlist,flv%kfit,flv%ifitp,flv%mfitp,flv%fitval,flv%nfc,lfno)
        goto 12
      elseif(word.eq.'BTUNE') then
        call ffs_init_sizep
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call petune(word,latt,twiss,
     $       mult,gammab,rlist(ifsize),nlist,
     $       rlist(itstr),rlist(itestr),nster,title,case,lfno)
        goto 12
      elseif(word.eq.'BALS') then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call pasex(word,wordp,
     $       latt,pos,twiss,
     $       mult,ilist(1,ifmast),gammab,title,case,
     $       rlist(itstr),rlist(itestr),nster,rlist(itmon),
     $       rlist(itemon),nmon,lfno)
        go to 12
      elseif(word.eq.'PETIN'.or.word.eq.'PETOUT') then
        call petcod(word,twiss,rlist(itmon),nmon)
        goto 12
      elseif(abbrev(word,'WR_ITE','_')) then
        call getwdl(word)
        if(word.eq.'RMATQ'.or.word.eq.'RMAT') then
          call pwrite(word,latt,twiss,
     $         gammab,pos,rlist(itmon),nmon)
        elseif(abbrev(word,'LAT_TICE','_').or.
     1         abbrev(word,'PS_NAME','_')) then
          call pwrlat(word,wordp,latt,mult,lfno)
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
        call tfevalb('WAKECOMMAND[]',13,kx,irtc)
        lfn1=lfnl0
        if(irtc .ne. 0 .or. kx .ne. ktftrue)then
          go to 2
        endif
        go to 10
      elseif(abbrev(word,'DELC_OR','_'))then
        call tcorr(word,latt,pos,ilist(1,ifmast),lfno)
        go to 10
      elseif(abbrev(word,'LTR_ACK','_'))then
        call tfltra(word,lfno)
        go to 10
      elseif(abbrev(word,'EMIT_TANCE','_'))then
        pspan=getva(exist)
        if(exist)then
          ia=0
          call rsetgl('PSPAN',pspan,ia)
        endif
        trpt0=trpt
        trpt=.false.
        call tfgeo(.true.)
        if(codplt)then
          call ffs_init_sizep
c          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call temitf(codplt,lfni,lfno)
        trpt=trpt0
        if(codplt)then
          updatesize=.true.
        else
          call tfgeo(.true.)
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
        call tfgeo(.true.)
        mphi2=9
        iparams=ktaloc(59)
        codin=0.d0
        if(codplt)then
          call ffs_init_sizep
c          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call temits(ndim,ntwissfun,mphi2,
     $     amus0,amus1,amusstep,
     $     emxe,emye,rese,rlist(iparams),
     $     lfni,lfno,int8(0),irtc)
        call tfree(iparams)
        if(codplt)then
          updatesize=.true.
        else
          call tfgeo(.true.)
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
        call tfsetparam
        call tfsize
        go to 10
      elseif(abbrev(word,'ALI_GN','_'))then
        call talign(latt,word,wordp,pos,lfno,exist)
        go to 30
      endif
 7000 call tfgetv(word,flv%ntouch,lfno,nextt,exist)
c      write(*,*)'tffsa ',nextt,exist,itt
      if(.not. exist)then
        if(itt .ge. 0)then
          call cssetp(nextt)
          call tfprintout(16.d0,irtc)
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
        call tfsetupcell(nlat,maxcond)
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
          call compelc(ilist(i-ii*2,ifklp+ii),cmp)
          if(idtype(cmp%id) .eq. 20)then
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
      if(tffsinitialcond(lfno,err))then
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
        call tffamsetup(1,em)
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
        call tffssetupwake(lfno,irtc)
        if(irtc .ne. 0)then
          call termes(lfno,'?Error in WakeFunction.',' ')
          go to 8810
        endif
        if(nwakep .eq. 0)then
          wake=.false.
        endif
      endif
      flv%iut=iutwiss(nlat,flv%nvar,nfcol,nfam,nut,.not. cell)
      if(flv%iut .le. 1)then
        call termes(lfno,
     $       '?Too many off-momentum or fit points.',' ')
        go to 8810
      endif
c      call tfevalb('Setup$FF[];Print["setupff ",FF$Orig]',36,kx,irtc)
      call tfevalb('Setup$FF[]',10,kx,irtc)
      call tffsmatch(df,dp0,r,nparallel,lfno,irtc)
      updatesize=.false.
      call tclrfpe
      if(wake)then
        call tffsclearwake
      endif
      if(irtc .ne. 0)then
        call tmunmapp(flv%iut)
        go to 8810
      endif
      call tfevalb('Reset$FF[]',10,kx,irtc)
      nqcol=nqcol-int(rfromk(kx))
      flv%nfc=nfc0
      call tfshow(cellstab,df,mfpnt,mfpnt1,
     $     kffs,irtcffs,lfnb .gt. 1,lfno)
      call tmunmapp(flv%iut)
      call tffsclearcouple(iele2)
      if(cell)then
        str=' '
        do i=nfam1,nfam
          if(.not. optstat(i)%stabx)then
            str='Horizontal'
            exit
          endif
        enddo
        do i=nfam1,nfam
          if(.not. optstat(i)%staby)then
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

      integer*8 function iutwiss(nlat,nvar,nfcol,nfam,nut,nonl)
      use tfstk
      use ffs, only:flv
      use ffs_pointer
      use tffitcode
      use ffs_wake
      use iso_c_binding
      use sad_main
      use mackw
      implicit none
      integer*4 nlat,nvar,nfcol,nfam,nut,i2
      integer*4 i,j,id,k
      integer*8 itmmapp
      real*8 tffsmarkoffset
      logical*4 nonl
      itwissp(3:nlat-1)=0
      if(idtypec(nlat-1) .eq. icMARK)then
        itwissp(nlat-1)=1
      endif
      if(nonl)then
        do i=2,nlat-1
          id=idtypec(i)
          if(id .ge. icSEXT .and. id .le. icDODECA
     $         .or. id .eq. icMULT .or. id .eq. icSOL)then
            itwissp(i)=1
            itwissp(i+1)=1
          endif
        enddo
      else
        do i=2,nlat-1
          id=idtypec(i)
          if(id .eq. icSOL)then
            itwissp(i)=1
            itwissp(i+1)=1
          endif
        enddo
      endif
      LOOP_I: do i=2,nlat-1
        do j=1,nvar
          if(ivarele(j) .eq. iele1(iele(i)))then
            itwissp(i)=1
            itwissp(i+1)=1
            cycle LOOP_I
          endif
        enddo
        if(iele2(i) .ne. 0)then
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
      i2=max(1,min(nlat,1+int(tffsmarkoffset(1))))
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
      call c_f_pointer(c_loc(rlist(iutwiss)),utwiss,
     $     [ntwissfun,2*nfam+1,nut])
      utwiss(1:ntwissfun,-nfam:nfam,1:nut)=>utwiss
      return
      end

      subroutine tffssaveparams(icmd,ilattp,err)
      use tfstk
      use ffs, local_ilattp=>ilattp
      implicit none
c nlocal = mcommon in TFFSLOCAL.inc
      integer*8 isave,ilattp
      integer*4 icmd,nxh
      logical*4 err
      nxh=int((sizeof(ffv)+7)/8)
      err=.false.
      if(icmd .eq. 0)then
        isave=ktaloc(nxh+2)
        klist(isave)=ilattp
        klist(isave+1)=iffssave
        call tmov(ffv,rlist(isave+2),nxh)
        iffssave=isave
      elseif(icmd .eq. 1)then
        if(iffssave .gt. 0)then
          isave=iffssave
          ilattp=klist(isave)
          call tmov(rlist(isave+2),ffv,nxh)
          iffssave=klist(isave+1)
          call tfree(isave)
        else
          err=.true.
        endif
      elseif(icmd .eq. 2)then
        err=.false.
        isave=iffssave
        do while(isave .gt. 0)
          if(klist(isave) .eq. ilattp)then
            err=.true.
            return
          endif
          isave=klist(isave+1)
        enddo
      elseif(icmd .eq. 3)then
        err=.false.
      elseif(icmd .eq. 4)then
        err=iffssave .eq. 0
      elseif(icmd .eq. -1)then
        do while(iffssave .gt. 0)
          isave=iffssave
          ilattp=klist(isave)
          call tmov(rlist(isave+2),ffv,nxh)
          iffssave=klist(isave+1)
          call tfree(isave)
        enddo
        err=.false.
      elseif(icmd .eq. -2)then
        isave=klist(iffssave+1)
        do while(isave .gt. 0)
          ilattp=klist(iffssave)
c          write(*,*)'tffssave -2: ',isave,ilattp
          call tmov(rlist(iffssave+2),ffv,nxh)
          call tfree(iffssave)
          iffssave=isave
          isave=klist(iffssave+1)
        enddo
      endif
      return
      end

      logical*4 function tfvcomp()
      use ffs_pointer
      use ffs, only:flv
      implicit none
      integer*4 i
      tfvcomp=.false.
      do i=1,flv%nvar
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
      use tffitcode
      implicit none
      call tffsadjust(flv%ntouch)
      call tfinitvar(flv%nvar)
      return
      end

      logical*4 function tffsinitialcond(lfno,err)
      use tfstk
      use ffs_fit
      use tffitcode
      implicit none
      type (sad_list), pointer :: klx,klj
      integer*8 kx
      integer*4 irtc,lfno,n,i,j,nfr1
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
            uini(:,0)=0.d0
c            call tclr(uini(1,0),28)
            iuid(0)=0
          else
            j=j+1
            if(tfreallistq(klx%body(j),klj))then
              if(klj%nl .eq. 6)then
                uini(mfitdx:mfitdx+5,i)=klj%rbody(1:6)
                uini(1,i)=-1.d0
                iuid(i)=j
              elseif(klj%nl .eq. ntwissfun)then
                uini(mfitax:mfitax+ntwissfun-1,i)=klj%rbody(1:ntwissfun)
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
      integer*8 iele2(nlat)
      integer*4 i
      do i=1,nlat-1
        if(iele2(i) .ne. 0)then
          call tfree(iele2(i))
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
      type (sad_symdef), pointer :: symd
      integer*8 iet
      integer*4 lfno,i,j,k,nk,m,me,nc, ie,ik,irtc
      character*(MAXPNAME) key,tfgetstrs
      type (sad_descriptor) itfelv,itfcoupk,kx,ki,kk,ke
      data itfelv%k,itfcoupk%k /0,0/
      klist(ifele2:ifele2+nlat-1)=0
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
        if(ktfsymbolqdef(kx%k,symd) .and. symd%value%k .eq. kx%k)then
          return
        endif
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
c            write(*,*)'setupcouple ',j,ik,key(1:10)
            if(ik .lt. 0)then
              do k=1,nlat-1
                if(ilist(k,ifele1) .eq. -ie)then
                  iet=klist(ifele2+k-1)
                  if(iet .eq. 0)then
                    iet=ktaloc(m+1)
                    ilist(1,iet)=0
                    klist(ifele2+k-1)=iet
                  endif
                  nk=ilist(1,iet)+1
c                  write(*,*)'setupcouple ',k,iet,ik,nk
                  ilist(1,iet+nk)=i
                  ilist(2,iet+nk)=-ik
                  ilist(1,iet)=nk
                  klist(ifele2+nlat-1)=1
                endif
              enddo
            elseif(ik .gt. 0)then
              iet=klist(ifele2+ie-1)
              if(iet .eq. 0)then
                iet=ktaloc(m+1)
                ilist(1,iet)=0
                klist(ifele2+ie-1)=iet
              endif
              nk=ilist(1,iet)+1
              ilist(1,iet+nk)=i
              ilist(2,iet+nk)=ik
              ilist(1,iet)=nk
              klist(ifele2+nlat-1)=1
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
      use ffs_pointer, only:idelc,idtypec
      use sad_main
      implicit none
      integer*4 i,ia,it,kl,l
      character*(*) keyword
      character*8 tfkwrd,kw
      if(i .gt. 0)then
        it=idtypec(i)
      else
        kl=ilist(-i,ifklp)
        it=idtypec(kl)
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

      subroutine tfsetupcell(nlat,maxcond)
      use ffs, only:flv
      use tffitcode
      implicit none
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

      recursive subroutine tclrline(line)
      use tfstk
      use maccode
      implicit none
      integer*4 i,n,idx
      integer*8 line
      do i=1,ilist(1,line)
        idx=ilist(2,line+i)
        if(idx .gt. 0 .and. idx .le. HTMAX)then
          if(idtype(idx) .eq. icLINE)then
            call tclrline(idval(idx))
          endif
        else
          write(*,*)'tclrline ',line,ilist(1,line),i,idx
        endif
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
