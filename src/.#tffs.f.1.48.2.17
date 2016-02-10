      subroutine tffs
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      
      logical*4 err
      call tffsa(1,kx,irtc)
      if(irtc .ne. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
      call tffssaveparams(-1,0,1,err)
      return
      end

      subroutine tffsinitparam(latt)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kdp
      type (sad_symdef), pointer :: symddp
      integer*4 latt(2,1),j
      real*8 rgetgl1,sigz0
      j=idval(latt(1,1))
      emx=rlist(j+kytbl(kwEMIX,icMARK))
      if(emx .le. 0.d0)then
        emx=rgetgl1('EMITX')
      else
        call rsetgl1('EMITX',emx)
      endif
      emy=rlist(j+kytbl(kwEMIY,icMARK))
      if(emy .le. 0.d0)then
        emy=rgetgl1('EMITY')
      else
        call rsetgl1('EMITY',emy)
      endif
      dpmax=max(0.d0,rlist(j+kytbl(kwDP,icMARk)))
      kdp=kxsymbolz('DP',2,symddp)
      if(dpmax .le. 1.d-30)then
        dpmax=rfromd(symddp%value)
      endif
      if(dpmax .le. 1.d-30)then
        dpmax=0.01d0
      endif
      symddp%value=dfromr(dpmax)
      if(rlist(latt(2,1)+kytbl(kwBX,icMARK)) .le. 0.d0)then
        rlist(latt(2,1)+kytbl(kwBX,icMARK))=1.d0
      endif
      if(rlist(latt(2,1)+kytbl(kwBY,icMARK)) .le. 0.d0)then
        rlist(latt(2,1)+kytbl(kwBY,icMARK))=1.d0
      endif
      sigz0=max(0.d0,rlist(j+kytbl(kwSIGZ,icMARk)))
      call rsetgl1('SIGZ',sigz0)
      return
      end

      subroutine tffsalloc(latt,nl)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 ktaloc
      integer*4 nl,latt(2,nl),j,l,ntwis,k,i,itehash
      if(idtype(latt(1,1)) .ne. icMARK)then
        write(*,*)'The first element must be a MARK element.'
        stop
      endif
      marki=1
      nlat=ilist(1,ilattp)+1
      call tfhashelement(latt)
      ifmult =ktaloc(nlat/2+1)
      ifele1=ktaloc(nlat/2+1)
      nele=0
      LOOP_L: do l=1,nlat-1
        j=ielmhash+itehash(pname(ilist(1,ilattp+l)),MAXPNAME)
        LOOP_K: do k=1,ilist(1,j)
          i=ilist(1,k+ilist(2,j)-1)
          if(i .ge. l)then
            exit LOOP_K
          endif
          if(ilist(1,ilattp+i) .eq. ilist(1,ilattp+l))then
            ilist(l,ifele1)=ilist(i,ifele1)
            cycle LOOP_L
          endif
        enddo LOOP_K
        nele=nele+1
        ilist(nele,ifmult)=l
        ilist(l,ifele1)=nele
      enddo LOOP_L
      ifcoup=ktaloc(nlat)
      iferrk=ktaloc(nlat*2)
      ifele =ktaloc(nlat/2+1)
      ifele2=ktaloc(nlat/2+1)
      ifival=ktaloc(nele/2+1)
      ifklp =ktaloc(nele/2+1)
      ifaux =ktaloc(nele*2)
      ifmast =ktaloc(nlat/2+1)
      call tffsfreebzl
      ifibzl =ktaloc(nlat*3/2+2)
      call tfinit(latt,
     1     rlist(ifival),rlist(ifklp),
     1     rlist(ifaux),rlist(ifele),rlist(ifele1),rlist(ifele2),
     $     rlist(ifibzl),
     1     rlist(iferrk),rlist(ifcoup),rlist(ifmult))
      ndim=1
      ndima=ndim*2+1
      ntwis =nlat*ndima
      nve=(nele+nlat)/2+10
      iftwis=ktaloc(ntwis*ntwissfun)
      ifpos =ktaloc(nlat)
      ifgeo =ktaloc(nlat*12)
      ifgamm=ktaloc(nlat)
      iftwissp=ktaloc(nlat/2+1)
      ifvarele=ktaloc(nve)
      ifvvar=ktaloc(nve)
      ifvalvar=ktaloc(nve*4)
      ifivcomp=ktaloc(nve)
      iftouchele=ktaloc(nve)
      iftouchv=ktaloc(nve)
c      ifsize=ktaloc(nlat*21)
      ifsize=0
      iwakepold=ktaloc(16)
      ifwakep=iwakepold
      call tclr(rlist(iwakepold),16)
      ilist(1,iwakepold)=1
      ilist(2,iwakepold)=int(ktaloc(1))
      rlist(ilist(2,iwakepold))=0.d0
      rlist(iwakepold+1)=0.d0
      ilist(1,iwakepold+2)=1
      ilist(2,iwakepold+2)=int(ktaloc(1))
      ilist(1,iwakepold+5)=ndim
      ilist(2,iwakepold+5)=nlat
      ilist(1,iwakepold+6)=int(iftwis)
      ilist(2,iwakepold+6)=int(ifsize)
      return
      end

      subroutine tffsfree
      use tfstk
      use ffs
      use tffitcode
      implicit none
      call tfresethash
      call tfree(ilist(2,ifwakep))
      call tfree(ilist(2,ifwakep+2))
      call tfree(ifwakep)
      call tfree(iftouchv)
      call tfree(iftouchele)
      call tfree(ifivcomp)
      call tfree(ifvalvar)
      call tfree(ifvvar)
      call tfree(ifvarele)
      call tfree(iftwissp)
      if(ifsize .gt. 0)then
        call tfree(ifsize)
      endif
c      call tfree(ifgamm)
      call tfree(ifgeo)
      call tfree(ifpos)
      call tfree(iftwis)
c      call tfree(ifibzl)
      call tfree(ifmast)
      call tfree(ifmult)
      call tfree(ifcoup)
      call tfree(ifaux)
      call tfree(ifival)
      call tfree(ifklp)
      call tfree(ifele2)
      call tfree(ifele1)
      call tfree(ifele)
      call tfree(iferrk)
      return
      end

      subroutine tffsfreebzl
      use tfstk
      use ffs
      use tffitcode
      implicit none
      if(ifibzl .ne. 0)then
        call tfree(ifibzl)
        ifibzl=0
      endif
      if(ifgamm .ne. 0)then
        call tfree(ifgamm)
        ifgamm=0
      endif
      return
      end

      subroutine tffsresetbzl
      use tfstk
      use ffs
      use tffitcode
      implicit none
      ifibzl=0
      ifgamm=0
      return
      end

      subroutine tfffs(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 outfl1,irtc,narg,
     $     lfno1,lfni1,lfn11,lfret,lfrecl,icslfno,icslfni,
     $     icslfn1,icsmrk,icslrecl,isp1,itfmessage
      character*10 strfromis
      narg=isp-isp1
      if(narg .gt. 2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      if(.not. ktfstringqd(dtastk(isp1+1),str))then
        irtc=itfmessage(9,'General::wrongtype','"String for #1"')
        return
      endif
      call tftclupdate(7)
      outfl1=outfl
      if(narg .eq. 2)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"File number for #2"')
          return
        endif
        outfl=int(rtastk(isp))
        if(outfl .eq. -1)then
          outfl=icslfno()
        endif
      else
        outfl=0
      endif
      levele=levele+1
      lfno1=icslfno()
      lfni1=icslfni()
      lfn11=icslfn1()
      lfret=icsmrk()
      lfrecl=icslrecl()
      call cssetp(lfrecl)
      call setbuf(str%str,str%nch)
      call cssetp(lfrecl)
      call tffsa(lfnp+1,kx,irtc)
      call tclrfpe
      call cssetp(lfret)
      call cssetl(lfrecl)
      call cssetlfno(lfno1)
      call cssetlfni(lfni1)
      call cssetlfn1(lfn11)
      outfl=outfl1
      if(irtc .eq. 0 .and. iffserr .ne. 0)then
        irtc=itfmessage(9,'FFS::error',strfromis(iffserr))
      endif
      call tfconnect(kx,irtc)
      return
      end

      block data tfblockdata
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      data fname /
     $     'RAD     ','RFSW    ','RADCOD  ','COD     ',
     1     'INTRA   ','TRPT    ','EMIOUT  ','GAUSS   ',
     1     'BIPOL   ','CELL    ','FFSPRMPT','DAPERT  ',
     1     'FIXSEED ','IDEAL   ','CODPLOT ','CANON   ',
     1     'POL     ','FLUC    ','CMPLOT  ','FOURIER ',
     1     'TRACKSIZ','SIMULATE','ABSW    ','JITTER  ',
     1     'TRGAUSS ','LWAKE   ','TWAKE   ','BARYCOD ',
     1     'BUNCHSTA','CONV    ','STABLE  ','SPAC    ',
     $     'RADLIGHT','GEOCAL  ','PHOTONS ','WSPAC   ',
     $     'SELFCOD ','PSPAC   ','CONVCASE','PRSVCASE',
     $     'LOSSMAP ','ORBITCAL','RADTAPER','SORG     ',
     $  20*'        '/
      data sino  /
     $     '        ','        ','        ','        ',
     1     '        ','RING    ','        ','UNIFORM ',
     1     'UNIPOL  ','INS     ','        ','        ',
     1     'MOVESEED','REAL    ','        ','        ',
     1     '        ','DAMPONLY','        ','        ',
     1     '        ','OPERATE ','RELW    ','QUIET   ',
     1     'TRUNI   ','        ','        ','        ',
     1     'BATCHSTA','        ','UNSTABLE','        ',
     $     '        ','GEOFIX  ','        ','        ',
     $     '        ','        ','        ','        ',
     $     '        ','        ','        ','        ',
     $  20*'        '/
      end
