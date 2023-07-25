      module track_tt
        integer*8 itt1,itt2,itt3,itt4,itt5,itt6
      end module

      module ffsa
      use tfstk
      use tfcsi
      integer*8 ::idum=-1

      contains
      subroutine tffs
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      logical*4 err
      call tffsa(1,lfni,kx,irtc)
      if(irtc /= 0 .and. ierrorprint /= 0)then
        call tfreseterror
      endif
      call tffssaveparams(-1,idum,err)
      return
      end

      subroutine tffsa(lfnb,lfn,kffs,irtcffs)
      use ffs
      use ffs_pointer
      use trackbypass
      use tffitcode
      use ffs_fit
      use ffs_wake
      use sad_main
      use match
      use track_tt
      use tparastat
      use temw, only:nparams
      use tfrbuf
      use calc,only:twmov
      use tfshare, only:tfresetsharedmap,tmunmapp
      use ffsfile
      use radint
      use geolib
      use modul,only:tfunblocksym
      use iso_c_binding
      implicit none
      integer*4 maxrpt,hsrchz
      type (sad_descriptor) ,intent(out):: kffs
      type (sad_descriptor) kx,tfvars
      real*8 ,allocatable,dimension(:)::vparams
      integer*8 itwisso,kax,ildummy
      integer*4 kk,i,lfnb,ia,iflevel,j,ielm,ielme,igelme,k1,k,
     $     irtc0,it,itt,lfn,
     $     iuse,l,levelr,lfnl0,lpw,meas0,mfpnta,igetgl,lenw,
     $     mphi2,next,nextt,nfp,
     $     nrpt1,itfpeeko,itfgetrecl,nl
      real*8 rmax,amus0,amus1,amusstep,apert,axi,ayi,ctime1,
     $     dpm2,dpxi,dpyi,em,emxe,emye,epxi,epyi,pspan,r2i,r3i,
     $     trval,rese,v,wa,wd,wl,xa,ya,xxa,xya,yya,getva,rgetgl1,
     $     wp,getvad,tgetgcut
      type (ffs_res) r
      type (sad_string),pointer::strc
      parameter (rmax=1.d35)
      parameter (maxrpt=32)
      character*256 word,wordp,title,case,tfgetstrv,tfgetstrs,tfgetstr
      character*(MAXPNAME) ename
      character*(MAXPNAME) , pointer :: pename
      character*(MAXPNAME+8) name
      character*16 autofg
      character*20 str
      integer*4 irtcffs,irtc,nc,nrpt(maxrpt),irptp(maxrpt)
      real*8 chi0(3),trdtbl(3,6),df(maxcond)
      logical*4 err,new,cmd,open98,abbrev,ftest,
     $     frefix,exist,init,expnd,chguse,visit,byeall,geocal0,busy
      save open98,exist,init
      save busy
      data busy /.false./
      itwisso(kk,i,j)=iftwis+kk+nlat*(i+ndim+(j-1)*ndima)-1
      flv%mcommon=int((sizeof(flv)+7)/8)
      kffs=dxnullo
      irtcffs=0
      l=itfuplevel()
      chguse=.false.
c     begin initialize for preventing compiler warning
      levelr=0
c     end   initialize for preventing compiler warning
 101  if(lfnb <= 1 .or. chguse)then
        call tffsalloc()
        if(.not. chguse)then
          call cputime(flv%ctime0,irtc0)
          flv%ctime2=flv%ctime0
        endif
        flv%iut=0
        flv%nvar=0
        flv%ntouch=0
        flv%setref=.false.
        trdtbl=0.d0
        itt1=0
        flv%itmax=40
        id1=1
        id2=nlat
        iorgr=1
        geo0(:,:)=geoini
        chi0=0.d0
        if(geocal .or. chguse)then
          geocal0=geocal
          geocal=.true.
          call tfgeo(.true.)
          geocal=geocal0
        endif
        if(lfnb <= 0)then
          go to 8900
        endif
        call tffsinitparam
c     kikuchi ... next 1 line added     (11/13/'91)
c        call corinit(newcor,nster,nmon,itstr,itestr,itmon,itemon)
c     
        flv%measp=nlat
        mfpnt=nlat
        mfpnt1=nlat
        flv%nfc=0
        call tfinitcalc
        call tmast
        call twmov(1,twiss,nlat,ndim,.true.)
        if(.not. chguse)then
          call tfevals(
     $         '{EMITX,EMITY,EMITZ,SIGZ,SIGE}='//
     $         'LINE[{"EMITX","EMITY","EMITZ","SIGMAZ","SIGE"},1];'//
     $         'DP0=LINE["DDP",1];',kx,irtc)
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
          wakeopt=.false.
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
          itbuf(5)=moderead
          call trbreset(5)
          if(infl /= 5)then
            lfnp=2
            call trbnextl(infl)
            lfnstk(2)=infl
          endif
        else
        endif
      else
        lfnp=lfnb
        lfnstk(lfnp)=lfn
      endif
      iffserr=0
      if(chguse)then
        ios=0
        chguse=.false.
        go to 10
      endif
      levelr=0
      iflevel=0
      if(lfni /= lfnstk(lfnp))then
        lfni=lfnstk(lfnp)
        call trbassign(lfni)
      endif
      lfno=outfl
 2    lfn1=merge(merge(lfno,merge(0,lfno,igetgl('$LOG$') == 0),
     $     lfni /= 5),0,lfnb==1)
      call csrst(lfn1)
      kffs=dxnullo
 10   continue
      if(iffserr /= 0)then
        if(lfnb > 1)then
          ios=1
        else
          iffserr=0
        endif
      endif
      if(ios > 0)then
        ios=0
        call tfclose(lfnp,lfnb)
        if(lfnp .lt. lfnb)then
          go to 9000
        endif
        lfn1=merge(lfno,merge(0,lfno,igetgl('$LOG$') == 0),
     $       lfni /= 5)
      elseif(ios .lt. 0)then
        ios=0
      endif
      call getwrd(word)
      if(ios /= 0)then
        go to 10
      endif
 12   if(word == ' ')then
        go to 10
      elseif(word(1:1) == '!')then
        go to 2
      endif
      ios=0
      call tfprint(word,lfno,.false.,itt,nextt,exist)
      if(exist .or. ios /= 0)then
        go to 10
      endif
      if(iffserr /= 0 .and. lfnb > 1)then
        go to 10
      endif
      if(word(1:1) == '!')then
        go to 2
      endif
      if(word == 'UNTIL')then
        if(levelr <= 0)then
          call termes(lfno,'?UNTIL without REPEAT.',' ')
          go to 2
        endif
        ftest=getva(exist) /= 0.d0
        if(exist)then
          if(ftest)then
            levelr=levelr-1
            go to 4010
          endif
        endif
        nrpt(levelr)=nrpt(levelr)-1
        if(nrpt(levelr) > 0)then
          call cssetp(irptp(levelr))
        else
          levelr=levelr-1
        endif
        ios=0
 4010   if(levelr == 0)then
          rep=.false.
        endif
        go to 10
      elseif(abbrev(word,'REP_EAT','_'))then
        levelr=levelr+1
        nrpt1=int(getva(exist))
        if(.not. exist)then
          nrpt1=65535
        endif
        nrpt(levelr)=max(1,nrpt1)
        irptp(levelr)=ipoint
        rep=.true.
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
      call tffile(word,lfnb,init,exist)
      if(lfnp .lt. lfnb)then
        go to 9000
      endif
      if(init)then
        ios=0
        go to 2
      endif
      if(exist)then
        go to 10
      endif
c      call tftrak(word,trdtbl,trval,lfno,exist)
c      if(exist)then
c        go to 10
c      endif
      call tfgetlineps(word,lenw(word),nl,kax,1,irtc)
      if(irtc == 0)then
        if(nl > 0)then
          go to 7000
        endif
      else
        if(irtc > 0 .and. ierrorprint /= 0)then
          call tfreseterror
        endif
        irtc=0
      endif

      cmd=.false.
      if(word == 'STOP' .or. word == 'EXIT')then
        call tffsadjustvar
        call tfsave(word,.false.,flv%ntouch)
        call tfgeo(.false.)
c        call corfree(newcor,nster,nmon,itstr,itestr,itmon,itemon)
        go to 8900
      elseif(word == 'QUIT')then
        call tfgeo(.false.)
c        call corfree(newcor,nster,nmon,itstr,itestr,itmon,itemon)
        go to 8900
      elseif(word == 'ABORT')then
        call tfresetsharedmap()
        stop
      elseif(word == 'USE' .or. word == 'VISIT')then
        visit=word == 'VISIT'
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
        it=itfpeeko(kx,next)
c        call tfdebugprint(kx,'USE ',1)
        call tfbeamline(kx,iuse,ename,irtc)
        if(iuse == 0)then
          if(ktfstringq(kx) .or. ktfsymbolq(kx))then
            word=tfgetstrs(kx,nc)
            iuse=hsrchz(word(1:nc))
          else
            call peekwd(word,next)
            iuse=hsrchz(word(1:lenw(word)))
          endif
        endif
        if(idtype(iuse) /= icLINE)then
          if(irtc /= 0)then
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
          if(irtc > 0 .and. ierrorprint /= 0)then
            call tfreseterror
          endif
          call tfclearlinep()
          call cssetp(next)
          call tffssaveparams(2,ilattp,err)
          calexp=expnd
          expnd=expnd .and. .not. err
          call tclrpara
          if(visit)then
            call tffssaveparams(0,ilattp,err)
            call tfblocksym('`FitFunction',12)
            call tfblocksym('`FitValue',9)
            call tfblocksym('`ElementValues',14)
            call tfblocksym('`OpticsProlog',13)
            call tfblocksym('`OpticsEpilog',13)
          else
            call tfunblocksym('`FitFunction',12,.false.)
            call tfunblocksym('`FitValue',9,.false.)
            call tfunblocksym('`ElementValues',14,.false.)
            call tfunblocksym('`OpticsProlog',13,.false.)
            call tfunblocksym('`OpticsEpilog',13,.false.)
            call tffsfree
            if(lattuse == lattredef)then
              call tclrline(lattredef)
              lattredef=0
            endif
          endif
          if(ilist(2,idval(iuse)) <= 0 .or. expnd)then
            call expnln(iuse)
          endif
          ilattp=idval(ilist(2,idval(iuse)))
          call loc_el(ilattp,elatt)
          nlat=elatt%nlat0+1
          call filaux(iuse)
          lattuse=ilattp
          call tclrpara
          dleng =rlist(elatt%aux+1)*rgetgl1('FSHIFT')
          ename=pname(iuse)
          call c_f_pointer(c_loc(flv%blname),pename)
          pename=ename
c          call tmovb(ename,flv%blname,MAXPNAME)
          call tfsetbeamlinename(ename)
          chguse=.true.
          go to 101
        endif
      elseif(word == 'BYE')then
        call peekwd(word,next)
        if(word == 'ALL')then
          call cssetp(next)
          byeall=.true.
        else
          byeall=.false.
        endif
        call tffssaveparams(4,ildummy,err)
        if(err)then
          if(byeall)then
            go to 10
          else
            call termes(lfno,'BYE without VISIT.',' ')
            go to 2
          endif
        endif
        call tclrpara
        call tffsfree
        if(byeall)then
          call tffssaveparams(-2,ildummy,err)
        endif
        call tffssaveparams(1,ilattp,err)
        call loc_el(ilattp,elatt)
        lattuse=ilattp
        nlat=elatt%nlat0+1
        call tfresetparam
        dleng =rlist(elatt%aux+1)*rgetgl1('FSHIFT')
        call ffs_init_pointer
        call ffs_twiss_pointer
        call c_f_pointer(c_loc(flv%blname),pename)
        ename=pename
c        call tmovb(flv%blname,ename,MAXPNAME)
        call tfsetbeamlinename(ename)
        call tfclearlinep()
        call tfunblocksym('`FitFunction',12,.true.)
        call tfunblocksym('`FitValue',9,.true.)
        call tfunblocksym('`ElementValues',14,.true.)
        call tfunblocksym('`OpticsProlog',13,.true.)
        call tfunblocksym('`OpticsEpilog',13,.true.)
      elseif(word == 'SPLIT')then
        call termes(lfno,
     $       'SPLIT is obsolete.   Use OFFSET of a marker.',' ')
      elseif(abbrev(word,'MAXI_TERATION','_'))then
        call tfgeti(flv%itmax,1.d0,word,lfno,exist)
        go to 31
      elseif(abbrev(word,'ATT_RIBUTE','_'))then
        call tfattr(word,lfno,exist,kffs,irtcffs,lfnb > 1)
        go to 30
      elseif(word == 'SAVE')then
        call tffsadjustvar
        call tfsave(word,.true.,flv%ntouch)
      elseif(abbrev(word,'EXPAND','_'))then
        call tffsadjustvar
        call tfsave(word,.true.,flv%ntouch)
      elseif(word == 'RESET')then
        call tfrst(word,.true.)
      elseif(word == 'FREE' .or. word == 'FIX')then
        frefix=word == 'FREE'
        call peekwd(word,next)
        if(word == ' ')then
ckikuchi ... next 5 lines added     (8/17/'90)
c$$$          call mcmess(lfno)
c$$$          if(frefix) then
c$$$            call mcfre(newcor)
c$$$          else
c$$$            call mcfix(latt,twiss,
c$$$     $           gammab,newcor)
c$$$          endif
          go to 10
        endif
        call tffsfreefix(frefix,flv%nvar,lfno)
      elseif(word == 'FIT')then
        call getwdl2(word,wordp)
        if(word == 'ALL')then
          mfpnt=1
          mfpnt1=nlat
        else
          mfpnt=ielme(wordp,exist,lfno)
          if(.not. exist)then
            mfpnt1=mfpnt
            go to 12
          endif
          call getwdl2(word,wordp)
          mfpnta=ielme(wordp,exist,lfno)
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
          if(flv%kfit(kk) <= mfit)then
            if(abs(flv%mfitp(kk)) > 2)then
              flv%mfitp(kk)=sign(abs(nfp)+1,flv%mfitp(kk))
            endif
          endif
        enddo
      elseif(abbrev(word,'CAL_CULATE','_'))then
        call tfcalc(word,nlist,flv%icalc,flv%ncalc,mfpnt,mfpnt1,exist)
        if(exist)then
          if(abbrev(word,'NOEXP_AND','_'))then
            word=' '
            calexp=.false.
          elseif(abbrev(word,'EXP_AND','_'))then
            word=' '
            calexp=.true.
          elseif(.not. keepexp)then
            calexp=.not. tfvcomp()
          endif
          fitflg=.false.
          go to 1000
        else
          go to 12
        endif
      elseif(word == 'GO')then
        call peekwd(word,next)
        if(abbrev(word,'NOEXP_AND','_'))then
          call cssetp(next)
          calexp=.false.
        elseif(abbrev(word,'EXP_AND','_'))then
          call cssetp(next)
          calexp=.true.
        elseif(.not. keepexp)then
          calexp=.not. tfvcomp()
        endif
        word=' '
        fitflg=.true.
        convgo=.false.
        go to 1000
      elseif(abbrev(word,'REC_OVER','_'))then
        do i=1,flv%nvar
          v=nvevx(i)%valvar
          nvevx(i)%valvar=nvevx(i)%valvar2
          nvevx(i)%valvar2=v
        enddo
c        ifvalvar2=ifvalvar+nve
c        do i=1,flv%nvar
c          v=rlist(ifvalvar+i-1)
c          rlist(ifvalvar+i-1)=rlist(ifvalvar2+i-1)
c          rlist(ifvalvar2+i-1)=v
c        enddo
        call tfsetv(flv%nvar)
      elseif(abbrev(word,'T_YPE','_'))then
        call tfsetparam
        call tftype(lfno,word)
        go to 12
      elseif(abbrev(word,'DISP_LAY','_'))then
        call tfsetparam
        call tfdisp(word,id1,id2,dp0*gammab(1),lfno,exist)
        go to 30
      elseif(word == 'SCALE')then
        call tscale(nlist,scale,lfno)
      elseif(word=='MAPANA') then
c        call gosadpls(latt,ilist(1,ifklp),
c     $       ilist(1,ifele),lfno)
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
      elseif(word == 'RADINT')then
        call tfsetparam
        call intgrl(lfno)
      elseif(abbrev(word,'REF_ERENCE','_'))then
        call tfsetref
      elseif(word == 'DIMAD')then
        call tfsetparam
        call tdimad(latt,mult,lfno)
      elseif(word == 'ZAP')then
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
            call nalign(latt,mult,ilist(1,ifmast),word,exist)
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
         call tfevalb('NISTACK$OBJ@Push[]',kx,irtc)
         go to 10
      else if(abbrev(word,'POPS_EED ','_'))then
         call tfevalb('NISTACK$OBJ@Pop[]',kx,irtc)
         go to 10
      else if(abbrev(word,'EXCGS_EED ','_'))then
         call tfevalb('NISTACK$OBJ@Exchange[]',kx,irtc)
         go to 10
      else if(abbrev(word,'PEEKS_EED ','_') .OR.
     &        abbrev(word,'PKS_EED','_')        )then
         call tfgeti(i,1.d0,word,lfno,exist)
         if( .not. exist) i=0
         write(word,'(A,I12,A)')'NISTACK$OBJ@Peek[',i,']'
         call tfevalb(word,kx,irtc)
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
     &             latt,icomp,rlist(ifcoup),
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
      if(word == 'VARY')then
        call tfchgv(lfno)
        call tfinitvar
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
      elseif(word == 'SBUNCH')then
        call termes(lfno,'?SBUNCH is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        call tfgetr(rlist(iwakepold+1),1.d0,word,lfno,exist)
c        go to 31
      elseif(word == 'DBUNCH')then
        call termes(lfno,'?DBUNCH is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        call tfdbun(word,rlist(ilist(2,iwakepold)),ilist(1,iwakepold),
c     1              lfno,err)
c        go to 32
      elseif(word == 'SLICE')then
        call termes(lfno,'?SLICE is obsolete.',' ')
        go to 10
c        iwakepold=ifwakep
c        ns=ilist(1,iwakepold+2)
c        call tfgeti(ns,1.d0,word,lfno,exist)
c        if(ns <= 0)then
c          call termes(lfno,'?Parameter out of range in SLICE.',' ')
c          go to 2
c        endif
c        ilist(1,iwakepold+2)=ns
c        go to 31
      elseif(abbrev(word,'ORI_GIN','_') .or. word == 'ORG')then
        iorgr=igelme(word,exist,lfno)
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
        geo0(:,1:3)=tfrotgeo(geo0(:,1:3),-chi0*scale(mfitchi1:mfitchi3))
c        geo0(:,1:3)=tfchitogeo(-chi0*scale(mfitchi1:mfitchi3))
        if(.not. exist)then
          go to 12
        endif
        go to 10
      elseif(word == 'SHOW')then
        call tshow(kffs,irtcffs,lfnb > 1,lfno)
        go to 10
      elseif(abbrev(word,'STAT_US','_'))then
        call cputime(ctime1,irtc0)
        call tfsetparam
        lpw=min(itfgetrecl()-1,255)
        call twbuf(word,'',lfno,1,0,0,0)
        call twbuf('CTIME='//autofg((ctime1-flv%ctime0)*1.d-6,'S8.3'),
     $       ' sec',lfno,1,lpw,8,1)
        call twbuf('DT='//autofg((ctime1-flv%ctime2)*1.d-6,'S8.3'),
     $       ' sec',lfno,1,lpw,8,1)
        flv%ctime2=ctime1
        call twbuf(word,'',lfno,1,0,0,-1)
        do i=1,nflag
          if(fname(i) /= ' ')then
            if(flags(i))then
              call twbuf(fname(i),'',lfno,1,lpw,8,1)
            else
              if(sino(i) /= ' ')then
                call twbuf(sino(i),'',lfno,1,lpw,8,1)
              else
                call twbuf('NO'//fname(i),'',lfno,1,lpw,8,1)
              endif
            endif
          endif
        enddo
        call twbuf(word,'',lfno,1,0,0,-1)
        call twbuf('CHARGE='//autofg(charge,'S7.4'),'',lfno,1,lpw,8,1)
        call twbuf('MASS='//autofg(amass/1e6,'S10.7'),' MeV',
     $       lfno,1,lpw,8,1)
        call twbuf('MOMENTUM='//autofg(pgev/1e9,'S10.6'),' GeV',
     $       lfno,1,lpw,8,1)
        call twbuf('PBUNCH='//autofg(pbunch,'S9.6'),'',lfno,1,lpw,8,1)
        flv%rsconv=rfromd(kxsymbolv('CONVERGENCE',11))
        call twbuf('DP='//autofg(dpmax,'S8.5'),'',lfno,1,lpw,8,1)
        call twbuf('DP0='//autofg(dp0,'S8.5'),'',lfno,1,lpw,8,1)
        call twbuf('FSHIFT='//autofg(rgetgl1('FSHIFT'),'S9.6'),'',
     $       lfno,1,lpw,8,1)
        call twbuf('MINCOUP='//autofg(coumin,'S9.6'),'',lfno,1,lpw,8,1)
        call twbuf('EMITX='//autofg(emx,'S9.6'),' m',lfno,1,lpw,8,1)
        call twbuf('EMITY='//autofg(emy,'S9.6'),' m',lfno,1,lpw,8,1)
        call twbuf('EMITZ='//autofg(emz,'S9.6'),' m',lfno,1,lpw,8,1)
        call twbuf('SIGZ='//autofg(sigzs,'S9.6'),' m',lfno,1,lpw,8,1)
        call twbuf('SIGE='//autofg(sizedp,'S9.6'),'',lfno,1,lpw,8,1)
        call twbuf('XIX='//autofg(xixf/pi2,'S7.4'),'',lfno,1,lpw,8,1)
        call twbuf('XIY='//autofg(xiyf/pi2,'S7.4'),'',lfno,1,lpw,8,1)
        call twbuf('NBUNCH='//autofg(anbunch,''),'',lfno,1,lpw,8,1)
        call twbuf('NP='//autofg(dble(np0),'S8.1'),'',lfno,1,lpw,8,1)
        call twbuf('GCUT='//autofg(tgetgcut(),''),'',lfno,1,lpw,8,1)
        call tfevals('StandardForm[$FORM="S8.3";"CODCONV="//CODCONV]',kx,irtc)
        if(ktfstringq(kx,strc))then
          call twbuf(strc%str(1:strc%nch),'',lfno,1,lpw,8,1)
        endif
        call twbuf('EMITDIV='//autofg(emidiv,''),'',lfno,1,lpw,8,1)
        call twbuf('EMITDIVB='//autofg(emidib,''),'',lfno,1,lpw,8,1)
        call twbuf('EMITDIVQ='//autofg(emidiq,''),'',lfno,1,lpw,8,1)
        call twbuf('EMITDIVS='//autofg(emidis,''),'',lfno,1,lpw,8,1)
        call twbuf('CONVERGENCE='//autofg(flv%rsconv,''),'',
     $       lfno,1,lpw,8,1)
        call twbuf('MAXITERATION='//autofg(dble(flv%itmax),''),'',
     $       lfno,1,lpw,8,1)
        call twbuf('NPARA='//autofg(dble(nparallel),''),'',
     $       lfno,1,lpw,8,1)
        call twbuf(' ','',lfno,1,lpw,8,-1)
        call twelm(lfno,mfpnt,mfpnt1,'FIT',lpw,8)
        call twelm(lfno,id1,id2,'DISP_LAY',lpw,8)
        call twelm(lfno,iorgr,0,'ORG',lpw,8)
        call twelm(lfno,flv%measp,0,'MEA_SURE',lpw,8)
        call twbuf(' ','',lfno,1,lpw,8,-1)
        go to 10
      elseif(abbrev(word,'VAR_IABLES','_') .or. word == 'VARS')then
        call tfinitvar
        kffs=tfvars(flv%nvar,irtcffs,lfnb > 1,lfno)
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
      elseif(word == 'DRAW')then
        call tfsetparam
        call tfevalb('System`CANVASDRAW[]',kx,irtc)
c        call tfdebugprint(kx,'CANVASDRAW[]',1)
c        write(*,'(i8)')irtc
        if(irtc /= 0 .or. ktfnonstringq(kx%k))then
c          title=Tfgetstrv('TITLE')
c          case=Tfgetstrv('CASE')
c          call twsdrw(latt,pos,ilist(1,ifele),
c     $         word,wordp,lfno,
c     1         twiss,0,rlist(itmon),rlist(itemon),nmon,
c     1         title,case,exist)
        else
          word=tfgetstr(kx,nc)
          exist=nc == 0
        endif
        go to 30
      elseif(word == 'GEO')then
        title=Tfgetstrv('TITLE')
        case=Tfgetstrv('CASE')
        call geodrw(rlist(ifgeo),word,lfno,title,case)
        go to 10
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
        vcalpha=1.d0
        title=tfgetstrv('TITLE')
        case=tfgetstrv('CASE')
        call trackb(latt,flv%measp,name,
     $       kffs,irtcffs,lfnb > 1,
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
      else
        cmd=.true.
      endif
      if(.not. cmd)then
        go to 10
      endif
      cmd=.false.
      if(word == 'WAKE')then
        lfnl0=lfn1
        lfn1=lfno
        call tfevalb('WAKECOMMAND[]',kx,irtc)
        lfn1=lfnl0
        if(irtc /= 0 .or. kx%k /= ktftrue)then
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
        call tfgeo(.true.)
        if(codplt)then
          call ffs_init_sizep
c          ilist(2,iwakepold+6)=int(ifsize)
        endif
        codin=twiss(1,0,mfitdx:mfitddp)
        call temitf(codplt,lfno)
        if(codplt)then
          modesize=6
        else
          call tfgeo(.true.)
        endif
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
        allocate(vparams(nparams))
c        iparams=ktaloc(nparams)
        codin=0.d0
        if(codplt)then
          call ffs_init_sizep
c          ilist(2,iwakepold+6)=int(ifsize)
        endif
        call temits(mphi2,amus0,amus1,amusstep,
     $     emxe,emye,rese,vparams,
     $     lfno,i00,irtc)
        deallocate(vparams)
c        call tfree(iparams)
        if(codplt)then
          modesize=6
        else
          call tfgeo(.true.)
        endif
c        dpmax=rgetgl1('SIGE')
c        rlist(itlookup('DP',ivtype))=dpmax
        gauss=.true.
        go to 10
 7310   call termes(lfno,
     $       'Usage: SYNCHROB_ETA nus_start nus_stop nus_step.',' ')
        go to 2
      elseif(abbrev(word,'BEAM_SIZE','_'))then
        call tfsetparam
        call tfsize(.true.)
        go to 10
      elseif(abbrev(word,'ALI_GN','_'))then
        call talign(latt,word,wordp,pos,lfno,exist)
        go to 30
      endif
 7000 call tfgetv(word,lfno,nextt,exist)
      if(.not. exist)then
        if(itt .ge. 0)then
          call cssetp(nextt)
          call tfprintout(16.d0,irtc)
          go to 10
        else
          call termes(lfno,
     1         '?Undefined command or element: ',word)
        endif
        if(ios /= 0)then
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
      if(err)then
        go to 2
      else
        go to 10
      endif
 1000 continue
      if(busy)then
        call termes(lfno,'?Recursive CAL or GO',' ')
        go to 2
      endif
      busy=.true.
      call tfsetparam
      nfc0=flv%nfc
      if(cell)then
        call tfsetupcell(nlat,maxcond)
      endif
      geomet=.false.
      do i=1,flv%ncalc
        if(flv%icalc(3,i) > mfitleng .and.
     $       flv%icalc(3,i) <= mfitchi3)then
          geomet=.true.
        endif
      enddo
      nfcol=0
      do k1=1,flv%nfc
        if(flv%kfit(k1) <= mfit .and. flv%mfitp(k1) /= 0)then
          if(flv%kfit(k1) > mfitleng .and.
     $         flv%kfit(k1) <= mfitchi3)then
            geomet=.true.
          endif
          nfcol=nfcol+1
          if(nfcol > maxcond)then
            call termes(lfno,'?Too many fit conditions.',' ')
            go to 8810
          endif
          flv%kfitp(nfcol)=k1
        endif
      enddo
      if(.not. geomet)then
        do i=1,nele
          if(idtypec(nelvx(i)%klp) == icSOL)then
            geomet=.true.
            exit
          endif
        enddo
      endif
      call tffssetupcouple(lfno)
      if(calexp)then
        call tffsadjustvar
      else
        call termes(lfno,
     $         'Info-Element values are not expanded.',' ')
      endif
      if(fitflg)then
        nvevx(1:flv%nvar)%valvar2=nvevx(1:flv%nvar)%valvar
      endif
      if(tffsinitialcond(lfno,err))then
        inicond=.true.
        nfam=nfr
        nfam1=merge(1-nfam,-nfam,iuid(-nfam) .lt. 0)
        uini(mfitddp,0)=0.d0
c        do i=nfam1,nfr
          kfam(nfam1:nfr)=0
          jfam(nfam1:nfr)=9999
          dp(nfam1:nfr)=uini(mfitddp,nfam1:nfr)
c        enddo
        jfam(0)=0
      elseif(err)then
        call termes(lfno,'Error in InitialOrbits',' ')
        go to 8810
      else
        inicond=.false.
        nfr=-1
        do k=1,flv%nfc
          if(flv%kfit(k) <= mfit .and. flv%mfitp(k) /= 0)then
            nfr=max((abs(flv%mfitp(k))-1)/2,nfr)
          endif
        enddo
        if(nfr > 0 .and. dpmax == 0.d0)then
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
          if(dpm2 > 0.d0)then
            if(i == -1)then
              dp(-1)=-dpm2
            elseif(i == 1)then
              dp( 1)=+dpm2
            endif
          endif
        enddo
        em=abs(emx)+abs(emy)
        call tffamsetup(1,em)
        nfam1=merge(1-nfam,-nfam,
     $       nfam > nfr .and. kfam(-nfam) == 0)
      endif
      wake=(twake .or. lwake) .and. trpt .and. wakeopt
      kwakep=0
      kwakeelm=0
      nwakep=0
      if(wake)then
        call tffssetupwake(lfno,irtc)
c        write(*,*)'tffsa-setupwake-done ',nwakep
        if(irtc /= 0)then
          call termes(lfno,'?Error in WakeFunction.',' ')
          go to 8810
        endif
        if(nwakep == 0)then
          wake=.false.
        endif
      endif
      flv%iut=iutwiss(flv%nvar,nfcol,nfam,nut,.not. cell)
      if(flv%iut <= 1)then
        call termes(lfno,
     $       '?Too many off-momentum or fit points.',' ')
        go to 8810
      endif
      call tfevalb('Setup$FF[]',kx,irtc)
      call tffsmatch(df,dp0,r,lfno,irtc)
      if(.not. setref)then
        call tfsetref
      endif
      call tclrfpe
      if(wake)then
        call tffsclearwake
      endif
      if(irtc /= 0)then
        call tmunmapp(flv%iut)
        go to 8810
      endif
      call tfevalb('Reset$FF[]',kx,irtc)
      nqcol=nqcol-int(kx%x(1))
      flv%nfc=nfc0
      call tfshow(cellstab,df,mfpnt,mfpnt1,kffs,irtcffs,lfnb > 1,lfno)
      call tmunmapp(flv%iut)
      call tffsclearcouple
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
            str=merge('Horizontal/Vertical',
     $                'Vertical           ',
     $           str(1:1) /= ' ')    
            exit
          endif
        enddo
        if(str(1:1) /= ' ')then
          write(lfno,*)str(1:len_trim(str)),' Unstable.'
        endif
      endif
      busy=.false.
      go to 12
 8810 busy=.false.
      flv%nfc=nfc0
      go to 10
 8900 call tffsfree
 9000 call tfconnect(kffs,irtcffs)
      return
      end

      subroutine mcmess(lfno)
      implicit none
      integer*4 ,intent(in):: lfno
      write(lfno,*)
     $     'Orbit Correction commands will be removed soon.',
     $     'Use orbit correction functions instead.',
     $     'Please send to Katsunobu.Oide@kek.jp, ',
     $     'if you still need those commands.'
      return
      end

      integer*8 function iutwiss(nvar,nfcol,nfam,nut,nonl)
      use ffs, only:flv,nvevx
      use ffs_pointer
      use tffitcode
      use ffs_wake
      use iso_c_binding
      use sad_main
      use tfshare,only:itmmapp
      use mackw
      implicit none
      integer*4 ,intent(in):: nvar,nfcol,nfam
      integer*4 ,intent(out):: nut
      integer*4 i2,i,j,id,k
      logical*4 ,intent(in):: nonl
      itwissp(3:nlat-1)=0
      if(idtypec(nlat-1) == icMARK)then
        itwissp(nlat-1)=1
      endif
      if(nonl)then
        do i=2,nlat-1
          id=idtypec(i)
          if(id .ge. icSEXT .and. id <= icDODECA
     $         .or. id == icMULT .or. id == icSOL)then
            itwissp(i)=1
            itwissp(i+1)=1
          endif
        enddo
      else
        do i=2,nlat-1
          id=idtypec(i)
          if(id == icSOL)then
            itwissp(i)=1
            itwissp(i+1)=1
          endif
        enddo
      endif
      LOOP_I: do i=2,nlat-1
        do j=1,nvar
          if(nvevx(j)%ivarele == iele1(icomp(i)))then
            itwissp(i)=1
            itwissp(i+1)=1
            cycle LOOP_I
          endif
        enddo
        if(kele2(i) /= 0)then
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
        if(flv%mfitp(j) /= 0 .and. flv%ifitp(j) > 0)then
          itwissp(flv%ifitp(j))=1
          if(flv%ifitp1(j) /= flv%ifitp(j))then
            itwissp(flv%ifitp1(j))=1
            if(flv%mfitp(j) .lt. 0 .and. flv%kfit(j) <= mfit
     $           .and. flv%kfit(j) /= mfitnx
     $           .and. flv%kfit(j) /= mfitny)then
              do i=min(flv%ifitp(j),flv%ifitp1(j)),
     $             max(flv%ifitp(j),flv%ifitp1(j))
                itwissp(i)=1
              enddo
            endif
          endif
        endif
      enddo
      do j=1,flv%ncalc
        if(flv%icalc(1,j) > 0)then
          itwissp(flv%icalc(1,j))=1
          if(flv%icalc(2,j) > 0)then
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
        if(itwissp(i) /= 0)then
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
      use ffs, local_ilattp=>ilattp
      use iso_c_binding
      implicit none
c nlocal = mcommon in TFFSLOCAL.inc
      integer*8 ,intent(inout):: ilattp
      integer*8 isave
      integer*8 , pointer ::kffv(:)
      integer*4 ,intent(in):: icmd
      integer*4 nxh
      logical*4 ,intent(out):: err
      nxh=int((sizeof(ffv)+7)/8)
      call c_f_pointer(c_loc(ffv),kffv,[nxh])
      err=.false.
      select case (icmd)
      case(0)
        isave=ktaloc(nxh+2)
        klist(isave)=ilattp
        klist(isave+1)=iffssave
        klist(isave+2:isave+nxh+1)=kffv
c        call tmov(ffv,rlist(isave+2),nxh)
        iffssave=isave
      case(1)
        if(iffssave > 0)then
          isave=iffssave
          ilattp=klist(isave)
          kffv=klist(isave+2:isave+nxh+1)
c          call tmov(rlist(isave+2),ffv,nxh)
          iffssave=klist(isave+1)
          call tfree(isave)
        else
          err=.true.
        endif
      case(2)
        err=.false.
        isave=iffssave
        do while(isave > 0)
          if(klist(isave) == ilattp)then
            err=.true.
            return
          endif
          isave=klist(isave+1)
        enddo
      case(3)
        err=.false.
      case(4)
        err=iffssave == 0
      case(-1)
        do while(iffssave > 0)
          isave=iffssave
          ilattp=klist(isave)
          kffv=klist(isave+2:isave+nxh+1)
c          call tmov(rlist(isave+2),ffv,nxh)
          iffssave=klist(isave+1)
          call tfree(isave)
        enddo
        err=.false.
      case(-2)
        isave=klist(iffssave+1)
        do while(isave > 0)
          ilattp=klist(iffssave)
          kffv=klist(iffssave+2:iffssave+nxh+1)
c          call tmov(rlist(iffssave+2),ffv,nxh)
          call tfree(iffssave)
          iffssave=isave
          isave=klist(iffssave+1)
        enddo
      end select
      return
      end

      logical*4 function tfvcomp()
      use ffs_pointer
      use ffs, only:flv,nvevx,nelvx
      implicit none
      integer*4 i
      tfvcomp=.false.
      do i=1,flv%nvar
        if(nvevx(i)%ivcomp /= 0 .and.
     $       nvevx(i)%ivvar /= nelvx(nvevx(i)%ivarele)%ival)then
          tfvcomp=.true.
          return
        endif
      enddo
      return
      end

      logical*4 function tffsinitialcond(lfno,err)
      use ffs_fit
      use tffitcode
      use eeval
      implicit none
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klj
      type (sad_descriptor) kx,iaini
      integer*4 ,intent(in):: lfno
      integer*4 irtc,n,i,j,nfr1
      logical*4 ,intent(out):: err
      data iaini%k/0/
      tffsinitialcond=.false.
      err=.true.
      if(iaini%k == 0)then
        iaini%k=ktfsymbolz('InitialOrbits',13)
      endif
      kx=tfsyeval(iaini,irtc)
      if(irtc /= 0)then
        call tfemes(irtc,'InitialOrbits',1,lfno)
        return
      endif
      if(tflistq(kx,klx))then
        n=klx%nl
        if(n <= 0)then
          return
        endif
        nfr=(n+1)/2
        if(nfr > ndimmax)then
          call termes(lfno,'Too many initial conditions.',' ')
          return
        endif
        nfr1=-(n/2)
        iuid(-nfr)=-1
        j=0
        uini(mfitbx,nfr1:nfr)=0.d0
        do i=nfr1,nfr
          if(i == 0)then
            uini(:,0)=0.d0
c            call tclr(uini(1,0),28)
            iuid(0)=0
          else
            j=j+1
            if(tfreallistq(klx%dbody(j),klj))then
              if(klj%nl == 6)then
                uini(mfitdx:mfitddp,i)=klj%rbody(1:6)
                uini(mfitbx,i)=0.d0
                iuid(i)=j
              elseif(klj%nl == ntwissfun)then
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
      elseif(kx%k == ktfoper+mtfnull .or. ktfsymbolq(kx))then
        err=.false.
        return
      endif
      return
      end

      subroutine tffsclearcouple
      use ffs
      use tffitcode
      use ffs_pointer,only:kele2
      implicit none
      integer*4 i
      do i=1,nlat-1
        if(kele2(i) /= 0)then
          call tfree(kele2(i))
          kele2(i)=0
        endif
      enddo
      kele2(nlat)=0
      return
      end

      subroutine tffssetupcouple(lfno)
      use ffs
      use ffs_pointer, only: kele2
      use tffitcode
      use eeval
      implicit none
      type (sad_dlist), pointer :: klx,kli
      type (sad_rlist), pointer :: kle
      type (sad_symdef), pointer :: symd
      integer*8 iet
      integer*4 ,intent(in):: lfno
      integer*4 i,j,k,nk,m,me,nc, ie,ik,irtc
      character*(MAXPNAME) key,tfgetstrs
      type (sad_descriptor) kx,ki,kk,ke
      type (sad_descriptor) , save :: itfelv,itfcoupk
      data itfelv%k,itfcoupk%k /0,0/
      kele2(1:nlat)=0
      if(itfelv%k == 0)then
        itfelv=kxsymbolz('`ElementValues',14)
        itfcoupk=kxsymbolz('`CoupledKeys',12)
      endif
      levele=levele+1
c      if(lfni > 100)then
c        write(*,*)'setupcoup-syeval-0 ',lfni,ipoint,lrecl,ios
c      endif
      kx=tfsyeval(itfcoupk,irtc)
c      if(lfni > 100)then
c        write(*,*)'setupcoup-syeval ',lfni,ipoint,lrecl,ios
c      endif
      call tfconnect(kx,irtc)
      if(irtc /= 0)then
        go to 9010
      endif
c      call tfdebugprint(kx,'setupcoup',3)
      if(.not. tflistq(kx,klx))then
        if(ktfsymbolqdef(kx%k,symd) .and. symd%value%k == kx%k)then
          return
        endif
        go to 9000
      endif
      m=klx%nl
      do i=1,m
        ki=klx%dbody(i)
        if(.not. tflistq(ki,kli))then
          go to 9000
        endif
        kk=kli%dbody(1)
        if(.not. ktfstringq(kk))then
          go to 9000
        endif
        key=tfgetstrs(kk,nc)
        ke=kli%dbody(2)
        if(tfreallistq(ke,kle))then
          me=kle%nl
          do j=1,me
            ie=int(kle%rbody(j))
            call tfkeya(ie,key,ik)
c            write(*,*)'setupcouple ',j,ik,key(1:10)
            if(ik .lt. 0)then
              do k=1,nlat-1
                if(ilist(k,ifele1) == -ie)then
                  iet=kele2(k)
                  if(iet == 0)then
                    iet=ktaloc(m+1)
                    ilist(1,iet)=0
                    kele2(k)=iet
                  endif
                  nk=ilist(1,iet)+1
c                  write(*,*)'setupcouple ',k,iet,ik,nk
                  ilist(1,iet+nk)=i
                  ilist(2,iet+nk)=-ik
                  ilist(1,iet)=nk
                  kele2(nlat)=1
                endif
              enddo
            elseif(ik > 0)then
              iet=kele2(ie)
              if(iet == 0)then
                iet=ktaloc(m+1)
                ilist(1,iet)=0
                kele2(ie)=iet
              endif
              nk=ilist(1,iet)+1
              ilist(1,iet+nk)=i
              ilist(2,iet+nk)=ik
              ilist(1,iet)=nk
              kele2(nlat)=1
            endif
          enddo
        endif
      enddo
      return
 9010 if(irtc > 0 .and. ierrorprint /= 0)then
        call tfreseterror
      endif
 9000 call termes(lfno,
     $     'ElementValues := { key[elem] :> f[ key1[elem1] ], ...}',
     $     ' ')
      return
      end

      subroutine tfkeya(i,keyword,ia)
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
      use sad_main
      implicit none
      integer*4 ,intent(in):: i
      integer*4 ,intent(out):: ia
      integer*4 it,kl,l
      character*(*) ,intent(in):: keyword
      character*8 tfkwrd,kw
      if(i > 0)then
        it=idtypec(i)
      else
        kl=nelvx(-i)%klp
        it=idtypec(kl)
      endif
      kw='-'
      l=0
      do while(kw /= ' ' .and. kw /= keyword)
        l=l+1
        kw=tfkwrd(it,l)
      enddo
      ia=merge(0,sign(l,i),kw == ' ')
      return
      end

      subroutine tfsetupcell(nlat,maxcond)
      use ffs, only:flv
      use tffitcode
      implicit none
      integer*4 ,intent(in):: maxcond,nlat
      integer*4 mfc
      integer*4 k
      if(flv%nfc > maxcond-2)then
        return
      endif
      mfc=1
      do k=1,flv%nfc
        if(flv%kfit(k) <= mfit)then
          mfc=max(abs(flv%mfitp(k)),mfc)
        endif
      enddo
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitbmagx
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=1.d0
      flv%mfitp(flv%nfc)=mfc
      flv%nfc=flv%nfc+1
      flv%kfit(flv%nfc)=mfitbmagy
      flv%ifitp(flv%nfc)=1
      flv%ifitp1(flv%nfc)=nlat
      flv%fitval(flv%nfc)=1.d0
      flv%mfitp(flv%nfc)=mfc
c      write(*,*)'setupcell ',flv%nfc
c$$$      flv%nfc=flv%nfc+1
c$$$      flv%kfit(flv%nfc)=mfitax
c$$$      flv%ifitp(flv%nfc)=1
c$$$      flv%ifitp1(flv%nfc)=nlat
c$$$      flv%fitval(flv%nfc)=0.d0
c$$$      flv%mfitp(flv%nfc)=mfc
c$$$      flv%nfc=flv%nfc+1
c$$$      flv%kfit(flv%nfc)=mfitay
c$$$      flv%ifitp(flv%nfc)=1
c$$$      flv%ifitp1(flv%nfc)=nlat
c$$$      flv%fitval(flv%nfc)=0.d0
c$$$      flv%mfitp(flv%nfc)=mfc
c$$$      flv%nfc=flv%nfc+1
c$$$      flv%kfit(flv%nfc)=mfitbx
c$$$      flv%ifitp(flv%nfc)=1
c$$$      flv%ifitp1(flv%nfc)=nlat
c$$$      flv%fitval(flv%nfc)=1.d0
c$$$      flv%mfitp(flv%nfc)=mfc
c$$$      flv%nfc=flv%nfc+1
c$$$      flv%kfit(flv%nfc)=mfitby
c$$$      flv%ifitp(flv%nfc)=1
c$$$      flv%ifitp1(flv%nfc)=nlat
c$$$      flv%fitval(flv%nfc)=1.d0
c$$$      flv%mfitp(flv%nfc)=mfc
      return
      end

      subroutine tscale(nlist,scale,lfno)
      use ffs
      use tffitcode
      implicit none
      integer*4 ,intent(in):: lfno
      integer*4 i
      character*(*) ,intent(in):: nlist(mfit)
      real*8 scale(mfit)
      write(lfno,9001)(nlist(i),scale(i),i=1,mfit)
9001  format(:3(a,1pg15.7,1x))
      return
      end

      recursive subroutine tclrline(line)
      use maccode
      implicit none
      integer*4 i,n,idx
      integer*8 ,intent(in):: line
      do i=1,ilist(1,line)
        idx=ilist(2,line+i)
        if(idx > 0 .and. idx <= HTMAX)then
          if(idtype(idx) == icLINE)then
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
      implicit none
      character*(*) ,intent(in):: str
      integer*4 ,intent(in):: nch
      type (sad_descriptor) ks,kx
      type (sad_symbol), pointer :: sym
      ks=kxsymbolf(str,nch,.false.)
      call descr_sym(ks,sym)
      kx=kxnaloc1(max(0,sym%gen),sym%loc)
      return
      end

      subroutine tfsetref
      use ffs
      use ffs_pointer
      implicit none
      twiss(:,-1,:)=twiss(:,0,:)
      setref=.true.
      return
      end

      recursive subroutine tfffs(isp1,kx,irtc)
      use ffs
      use tffitcode
      use tfcsi
      use readbuf, only:trbopen,trbopenmap
      use tfrbuf, only:modestring,trbassign,trbclose
      use iso_c_binding
      use ffsfile, only:lfnp
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      type (csiparam) sav
      integer*4 outfl1,irtc,narg,lfn,isp1,itfmessage,itfmessagestr
      character*10 strfromis
      narg=isp-isp1
      if(narg > 2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      if(.not. ktfstringq(dtastk(isp1+1),str))then
        irtc=itfmessage(9,'General::wrongtype','"String for #1"')
        return
      endif
      call tftclupdate(7)
      outfl1=outfl
      if(narg == 2)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"File number for #2"')
          return
        endif
        outfl=int(rtastk(isp))
        if(outfl == -1)then
          outfl=icslfno()
        endif
      else
        outfl=0
      endif
      levele=levele+1
      sav=savep
      call trbopen(lfn,ktfaddr(ktastk(isp1+1)),
     $     int8(modestring),str%nch)
      if(lfn <= 0)then
        irtc=itfmessagestr(9,'FFS::lfn',str%str(1:str%nch))
        kx%k=ktfoper+mtfnull
      else
        call trbassign(lfn)
        ipoint=1
        lrecl=0
        call tffsa(lfnp+1,lfn,kx,irtc)
        call trbclose(lfn)
        call tclrfpe
        savep=sav
        call trbassign(lfni)
        outfl=outfl1
        if(irtc == 0 .and. iffserr /= 0)then
          irtc=itfmessagestr(9,'FFS::error',strfromis(int(iffserr)))
        endif
      endif
      call tfconnect(kx,irtc)
      return
      end

      end module
