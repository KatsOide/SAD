      subroutine tfdeffun3ep
      implicit none
      integer*4 i,itfunaloc,map(32),ieval(32)
c     Initialize map/ieval array
      do i=1,32
         map(i)=0
         ieval(i)=0
      enddo
c     Channel Access stuff
      i=itfunaloc('EPICS$CaOpen',3101,1,map,ieval,0)
      i=itfunaloc('EPICS$CaPendIO',3102,1,map,ieval,0)
      i=itfunaloc('EPICS$CaPendEvent',3103,1,map,ieval,0)
      i=itfunaloc('EPICS$CaPut',3104,2,map,ieval,0)
      i=itfunaloc('EPICS$CaInit',3105,0,map,ieval,0)
      i=itfunaloc('EPICS$CaFlushIO',3106,0,map,ieval,0)
      i=itfunaloc('EPICS$CaFieldType',3107,1,map,ieval,0)
      i=itfunaloc('EPICS$CaElementCount',3108,1,map,ieval,0)
      i=itfunaloc('EPICS$CaAddEvent',3109,3,map,ieval,0)
      i=itfunaloc('EPICS$CaDebugPrint',3110,1,map,ieval,0)
      i=itfunaloc('EPICS$CaClearChannel',3111,1,map,ieval,0)
      i=itfunaloc('EPICS$CaClearEvent',3112,1,map,ieval,0)
      i=itfunaloc('EPICS$CaHostName',3113,1,map,ieval,0)
      i=itfunaloc('EPICS$CaPutCB',3114,2,map,ieval,0)
c     Static Database Access stuff
      i=itfunaloc('EPICS$DbReadDatabase',3201,3,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteRecord',3202,1,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteMenu',3203,1,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteRecordType',3204,1,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteDevice',3205,1,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteDriver',3206,1,map,ieval,0)
      i=itfunaloc('EPICS$DbWriteBreaktable',3207,1,map,ieval,0)
      i=itfunaloc('EPICS$DbPath',3208,1,map,ieval,0)
      i=itfunaloc('EPICS$DbAddPath',3209,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetNRecordTypes',3210,1,map,ieval,0)
      i=itfunaloc('EPICS$DbFindRecordType',3211,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFirstRecordType',3212,0,map,ieval,0)
      i=itfunaloc('EPICS$DbNextRecordType',3213,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetRecordTypeName',3214,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetNFields',3215,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFirstField',3216,1,map,ieval,0)
      i=itfunaloc('EPICS$DbNextField',3217,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetFieldType',3218,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetFieldName',3219,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetDefault',3220,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetPrompt',3221,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetPromptGroup',3222,0,map,ieval,0)
      i=itfunaloc('EPICS$DbPutRecordAttribute',3223,0,map,ieval,0)
      i=itfunaloc('EPICS$DbdbGetRecordAttribute',3224,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetNRecords',3225,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFindRecord',3226,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFirstRecord',3227,0,map,ieval,0)
      i=itfunaloc('EPICS$DbNextRecord',3228,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetRecordName',3229,0,map,ieval,0)
      i=itfunaloc('EPICS$DbCreateRecord',3230,2,map,ieval,0)
      i=itfunaloc('EPICS$DbDeleteRecord',3231,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFreeRecord',3232,0,map,ieval,0)
      i=itfunaloc('EPICS$DbCopyRecord',3233,0,map,ieval,0)
      i=itfunaloc('EPICS$DbRenameRecord',3234,0,map,ieval,0)
      i=itfunaloc('EPICS$DbVisibleRecord',3235,0,map,ieval,0)
      i=itfunaloc('EPICS$DbInvisibleRecord',3236,0,map,ieval,0)
      i=itfunaloc('EPICS$DbIsVisibleRecord',3237,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFindField',3238,2,map,ieval,0)
      i=itfunaloc('EPICS$DbFoundField',3239,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetString',3240,1,map,ieval,0)
      i=itfunaloc('EPICS$DbPutString',3241,2,map,ieval,0)
      i=itfunaloc('EPICS$DbVerify',3242,2,map,ieval,0)
      i=itfunaloc('EPICS$DbGetRange',3243,0,map,ieval,0)
      i=itfunaloc('EPICS$DbIsDefaultValue',3244,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetNMenuChoices',3245,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetMenuChoices',3246,1,map,ieval,0)
      i=itfunaloc('EPICS$DbGetMenuIndex',3247,0,map,ieval,0)
      i=itfunaloc('EPICS$DbPutMenuIndex',3248,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetMenuStringFromIndex',3249,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetMenuIndexFromString',3250,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFindMenu',3251,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetNLinks',3252,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetLinkField',3253,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetLinkType',3254,0,map,ieval,0)
      i=itfunaloc('EPICS$DbCvtLinkToConstant',3255,0,map,ieval,0)
      i=itfunaloc('EPICS$DbCvtLinkToPvlink',3256,0,map,ieval,0)
      i=itfunaloc('EPICS$DbAllocForm',3257,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFreeForm',3258,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetFormPrompt',3259,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetFormValue',3260,0,map,ieval,0)
      i=itfunaloc('EPICS$DbPutForm',3261,0,map,ieval,0)
      i=itfunaloc('EPICS$DbVerifyForm',3262,0,map,ieval,0)
      i=itfunaloc('EPICS$DbGetRelatedField',3263,0,map,ieval,0)
      i=itfunaloc('EPICS$DbFindBrkTable',3264,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpPath',3265,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpRecord',3266,1,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpMenu',3267,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpRecordType',3268,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpFldDes',3269,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpDevice',3270,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpDrive',3271,0,map,ieval,0)
      i=itfunaloc('EPICS$DbDumpBreaktable',3272,0,map,ieval,0)
      i=itfunaloc('EPICS$DbPvdDump',3273,0,map,ieval,0)
      i=itfunaloc('EPICS$DbReportDeviceConfig',3274,0,map,ieval,0)
      return
      end

      subroutine tfefun3ep(isp1,id,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc, id
      
      go to(
     $   1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
     $   1110,1120,1130,1140
     $      ),id-3100
c        Srch,PnIO,PnEv, Put,Init,FlIO,FlTy,ElCn,AdEv,DbgP
c        ClCh,ClEv,HstN,PtCB
      go to(
     $   2901,2900,2990,2990,2990,2990,2990,2900,2900,2930,
     $   2990,2910,2930,2920,2940,2940,2940,2930,2920,2920,
     $   2990,2990,2990,2990,2990,2990,2990,2990,2990,2960,
     $   2990,2990,2990,2990,2990,2990,2990,2960,2990,2920,
     $   2970,2980,2990,2990,2930,2950,2990,2990,2990,2990,
     $   2990,2990,2990,2990,2990,2990,2990,2990,2990,2990,
     $   2990,2990,2990,2990,2990,2900,2990,2990,2990,2990,
     $   2990,2990,2990,2990,2990,2990,2990,2990,2990,2990,
     $   2990,2990,2990,2990,2990,2990,2990,2990,2990,2990
     $     )id-3200
c        RdDb,WtRc,    ,    ,    ,    ,    ,Path,AdPt,GNRT,
c            ,FsRt,NxRt,GtRN,GtNF,FsFd,NxFd,GtFT,GtFn,GtDf,
c            ,    ,    ,    ,    ,    ,    ,    ,    ,CrtR,
c            ,    ,    ,    ,    ,    ,    ,FnFl,    ,GtSt,
c        PtSt,Vrfy,    ,    ,GtNM,GtMC,    ,    ,    ,    ,
c            ,    ,    ,    ,    ,    ,    ,    ,    ,    ,
c            ,    ,    ,    ,    ,DpRc,    ,    ,    ,    ,
c            ,    ,    ,    ,    ,    ,    ,    ,    ,    ,
c            ,    ,    ,    ,    ,    ,    ,    ,    ,    ,
      goto 2990

 1010 call tfcasearch(isp1,kx,irtc)
      return
 1020 call tfcapendio(isp1,kx,irtc)
      return
 1030 call tfcapendevent(isp1,kx,irtc)
      return
 1040 call tfcaput(isp1,kx,irtc)
      return
 1050 call tfcainit(isp1,kx,irtc)
      return
 1060 call tfcaflushio(isp1,kx,irtc)
      return
 1070 call tfcafieldtype(isp1,kx,irtc)
      return
 1080 call tfcaelementcount(isp1,kx,irtc)
      return
 1090 call tfcaaddevent(isp1,kx,irtc)
      return
 1100 call tfcadebugprint(isp1,kx,irtc)
      return
 1110 call tfcaclearchannel(isp1,kx,irtc)
      return
 1120 call tfcaclearevent(isp1,kx,irtc)
      return
 1130 call tfcahostname(isp1,kx,irtc)
      return
 1140 call tfcaputcb(isp1,kx,irtc)
      return

 2900 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbstr(isp1,kx,irtc)
      isp=isp-1
      return
 2901 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbstrx3(isp1,kx,irtc)
      isp=isp-1
      return
 2910 isp=isp+1
      rtastk(isp)=id-3200
      call tfdb2int(isp1,kx,irtc)
      isp=isp-1
      return
 2920 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbint2str(isp1,kx,irtc)
      isp=isp-1
      return
 2930 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbint2int(isp1,kx,irtc)
      isp=isp-1
      return
 2940 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbintint2int(isp1,kx,irtc)
      isp=isp-1
      return
 2950 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbint2strs(isp1,kx,irtc)
      isp=isp-1
      return
 2960 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbintstr2int(isp1,kx,irtc)
      isp=isp-1
      return
 2970 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbintstr(isp1,kx,irtc)
      isp=isp-1
      return
 2980 isp=isp+1
      rtastk(isp)=id-3200
      call tfdbintstr2str(isp1,kx,irtc)
      isp=isp-1
      return
 2990 write(*,*) 'Unimplimented Db function is called.'
      kx=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfcasearch(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonstringq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      ka=ktfaddr(ktastk(isp))
      jlist(ilist(1,ka)+1,ka+1)=0
      call ecasearchandconnect(ilist(1,ka+1),kx,ilist(1,1))
      irtc=0
      return
      end

      subroutine tfcaaddevent(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,evtype,vtype,itfmessage
      real*8 chid
      if(isp .ne. isp1+3) then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+3)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      chid = rtastk(isp1+1)
      vtype = rtastk(isp1+2)
      evtype = rtastk(isp1+3)
      if (evtype.ne.0) call ecaaddarrayevent(chid,kx,vtype,evtype)
      irtc=0
      return
      end

      subroutine tfcapendio(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1, irtc,itfmessage
      real*4 t
      
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real Number"')
        return
      endif
      t = rtastk(isp)
      call ecapendio(t)
      kx=0
      irtc=0
      return
      end


      subroutine tfcaflushio(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc
      
c      integer*4 itfmessage
c      if(isp .ne. isp1)then
c        irtc=itfmessage(9,'General::narg','"0"')
c        return
c      endif
      call ecaflushio()
      kx=0
      irtc=0
      return
      end

      subroutine tfcapendevent(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1, irtc,itfmessage
      real*4 t
      
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real Number"')
        return
      endif
      t = rtastk(isp)
      call ecapendevent(t)
      kx=0
      irtc=0
      return
      end

      subroutine tfcaput(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,nret,itfmessage
      real*8 chid
      if (isp .ne. isp1+2) then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      chid=rtastk(isp1+1)
      nret=-1
      if(ktfrealq(ktastk(isp)))then
         call ecaput(chid,.false.,1,0,rtastk(isp),nret)
      elseif(ktfstringq(ktastk(isp)))then
         call ecaput(chid,.false.,0,0,ktfaddr(ktastk(isp)),nret)
      elseif(ktflistq(ktastk(isp)))then
         ka=ktfaddr(ktastk(isp))
         if(ilist(2,ka-1) .eq. 0)then
            nret=0
         else
           if(ktfreallistq(ka))then
             call ecaput(chid,.false.,1,ilist(2,ka-1),ka,nret)
           elseif(ktfrealq(klist(ka+1)))then
             call ecaput(chid,.false.,1,ilist(2,ka-1),ka,nret)
           else
             call ecaput(chid,.false.,0,ilist(2,ka-1),ka,nret)
           endif
         endif
      endif

      if(nret .eq. 0)then
        kx=ktfoper+mtfnull
        irtc=0
      elseif(nret .gt. 0)then
         irtc=itfmessage(9,'CA::Channel','"put error"')
      else
         irtc=itfmessage(9,'General::wrongtype',
     $        '"Real, String, List of Reals, List of Strings"')
      endif
      return
      end

      subroutine tfcaputcb(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,nret,itfmessage
      real*8 chid
      if (isp .ne. isp1+2) then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      chid = rtastk(isp1+1)
      nret=-1
      if(ktfrealq(ktastk(isp)))then
         call ecaput(chid,.true.,1,0,rtastk(isp),nret)
      elseif(ktfstringq(ktastk(isp)))then
         call ecaput(chid,.true.,0,0,ktfaddr(ktastk(isp)),nret)
      elseif(ktflistq(ktastk(isp)))then
         ka=ktfaddr(ktastk(isp))
         if(ilist(2,ka-1) .eq. 0)then
            nret=0
         else
            if(ktfreallistq(ka))then
               call ecaput(chid,.true.,1,ilist(2,ka-1),ka,nret)
            elseif(ktfrealq(klist(ka+1)))then
               call ecaput(chid,.true.,1,ilist(2,ka-1),ka,nret)
            else
               call ecaput(chid,.true.,0,ilist(2,ka-1),ka,nret)
            endif
         endif
      endif

      if(nret .eq. 0)then
        kx=ktfoper+mtfnull
        irtc=0
      else
         irtc=itfmessage(9,'General::wrongtype',
     $        '"Real, String, List of Reals, List of Strings"')
      endif
      return
      end

      subroutine tfcavaluecb(chid,stat,sev,t,type,nc,karray)
      use tfstk
      use efun
      use eeval
      implicit none
      type (sad_descriptor) kx
      integer*8 karray(nc)
      real*8 chid
      integer*4 stat,sev,type,nc
      real*8 t
      integer*4 isp0,isp2,irtc,i,l,itfdownlevel
      
      type (sad_descriptor) iaepicsvaluecb
      data iaepicsvaluecb%k /0/
      isp=isp+1
      isp0=isp
      if(iaepicsvaluecb%k .eq. 0)then
        iaepicsvaluecb%k=ktfsymbolz('EPICS$ValueCB',13)
      endif
      levele=levele+1
      dtastk(isp)=tfsyeval(iaepicsvaluecb,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      isp=isp+1
      rtastk(isp)=chid
      isp=isp+1
      if(type .eq. 0)then
         if(nc .eq. 1)then
            ktastk(isp)=karray(1)
         else
            isp2=isp
            do i=1,nc
               isp=isp+1
               ktastk(isp)=karray(i)
            enddo
            dtastk(isp2)=kxmakelist(isp2)
            isp=isp2
         endif
      else
         if(nc .eq. 1)then
           ktastk(isp)=karray(1)
         else
           dtastk(isp)=kxm2l(rlist(ksad_loc(karray(1)):),
     $          0,nc,nc,.false.)
         endif
      endif
      isp=isp+1
      isp2=isp
      isp=isp+1
      rtastk(isp)=stat
      isp=isp+1
      rtastk(isp)=sev
      isp=isp+1
      rtastk(isp)=t
      ktastk(isp2)=ktflist+ktfmakelist(isp2)
      isp=isp2
      kx=tfefunref(isp0,.true.,irtc)
 9000 if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
      l=itfdownlevel()
      isp=isp0-1
      return
      end

      subroutine tfcainit(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc
      call ecainit()
      kx=ktfoper + mtfnull
      irtc=0
      return
      end

      subroutine tfcafieldtype(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr
      integer*4 isp1,irtc,type,itfmessage
      real*8 chid
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      chid = rtastk(isp)
      call ecafieldtype(chid, type)
      kx=kfromr(dble(type))
      irtc=0
      return
      end

      subroutine tfcaelementcount(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr
      integer*4 isp1,irtc,count,itfmessage
      real*8 chid
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      chid = rtastk(isp)
      call ecaelementcount(chid, count)
      irtc=0
      kx = kfromr(dble(count))
      return
      end

      subroutine tfcadebugprint(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,bdbg,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      bdbg = rtastk(isp)
      call ecadebugprint(bdbg)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end


      subroutine tfcaclearchannel(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,itfmessage
      real*8 chid
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      chid = rtastk(isp)
      call ecaclearchannel(chid)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfcaclearevent(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,itfmessage
      real*8 evid
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      evid = rtastk(isp)
      call ecaclearevent(evid)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end


      subroutine tfcahostname(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,itfmessage
      real*8 chid
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      chid=rtastk(isp)
      call ecahostname(kx,chid)
      irtc=0
      return
      end

      subroutine tfdbstr(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,nfunc,iret,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      ka=ktfaddr(ktastk(isp1+1))
      jlist(ilist(1,ka)+1,ka+1)=0
      nfunc = rtastk(isp)
      call cdbstr(nfunc,ilist(1,ka+1),iret)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfdbstrx3(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka1,ka2,ka3
      integer*4 isp1,irtc,nfunc,iret,itfmessage
      if(isp .ne. isp1+4)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+3)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      ka1=ktfaddr(ktastk(isp1+1))
      jlist(ilist(1,ka1)+1,ka1+1)=0
      ka2=ktfaddr(ktastk(isp1+2))
      jlist(ilist(1,ka2)+1,ka2+1)=0
      ka3=ktfaddr(ktastk(isp1+3))
      jlist(ilist(1,ka3)+1,ka3+1)=0
      nfunc = rtastk(isp)
      call cdbstrx3(nfunc,ilist(1,ka1+1),ilist(1,ka2+1),
     $     ilist(1,ka3+1),iret)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfdb2int(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr
      integer*4 isp1,irtc,nfunc,iret,nret,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"0"')
        return
      endif
      nfunc = rtastk(isp)
      call cdb2int(nfunc,nret,iret)
      kx = kfromr(dble(nret))
      irtc=0
      return
      end

      subroutine tfdbint2str(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,nfunc,itfmessage,iarg
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      nfunc=rtastk(isp)
      iarg=rtastk(isp1+1)
      call cdbint2str(kx,nfunc,iarg)
      irtc=0
      return
      end

      subroutine tfdbint2int(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr
      integer*4 isp1,irtc,nfunc,iret,iarg,nret,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      iarg = rtastk(isp1+1)
      nfunc = rtastk(isp)
      call cdbint2int(nfunc,iarg,nret,iret)
      kx = kfromr(dble(nret))
      irtc=0
      return
      end

      subroutine tfdbintint2int(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr
      integer*4 isp1,irtc,nfunc,iret,iarg1,iarg2,nret,itfmessage
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      iarg1 = rtastk(isp1+1)
      iarg2 = rtastk(isp1+2)
      nfunc = rtastk(isp)
      call cdbintint2int(nfunc,iarg1,iarg2,nret,iret)
      kx = kfromr(dble(nret))
      irtc=0
      return
      end

      subroutine tfdbint2strs(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,isp0,isp2,nfunc,iarg,np,i,
     $     itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      nfunc=rtastk(isp)
      iarg=rtastk(isp1+1)
      isp0=isp
      call cdbint2strs(nfunc,iarg,itastk(1,isp0+1),np)
      isp=isp+np+1
      isp2=isp
      do i=1,np
         isp=isp+1
         ktastk(isp)=ktastk(isp0+i)
      enddo
      if(np .lt. 1)then
         isp=isp+1
         ktastk(isp)=ktfoper+mtfnull
      endif         
      kx=ktflist+ktfmakelist(isp2)
      isp=isp0
      irtc=0
      return
      end

      subroutine tfdbintstr2int(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kfromr,ka
      integer*4 isp1,irtc,nfunc,iret,iarg,nret,itfmessage
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      iarg = rtastk(isp1+1)
      ka=ktfaddr(ktastk(isp1+2))
      jlist(ilist(1,ka)+1,ka+1)=0
      nfunc = rtastk(isp)
      call cdbintstr2int(nfunc,iarg,ilist(1,ka+1),nret,iret)
      kx=kfromr(dble(nret))
      irtc=0
      return
      end

      subroutine tfdbintstr(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,nfunc,iret,iarg,itfmessage
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      iarg = rtastk(isp1+1)
      ka=ktfaddr(ktastk(isp1+2))
      jlist(ilist(1,ka)+1,ka+1)=0
      nfunc = rtastk(isp)
      call cdbintstr(nfunc,iarg,ilist(1,ka+1),iret)
      kx=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfdbintstr2str(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ka
      integer*4 isp1,irtc,nfunc,iarg, itfmessage
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"String"')
        return
      endif
      nfunc = rtastk(isp)
      iarg = rtastk(isp1+1)
      ka=ktfaddr(ktastk(isp1+2))
      call cdbintstr2str(kx,nfunc,iarg,ilist(1,ka+1))
      irtc=0
      return
      end
