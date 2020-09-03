      subroutine tfinitn
      use tfstk
      use tfrbuf
      use tfcsi
      use macphys
      use version
      implicit none
      type (sad_symdef), pointer :: contd
      type (sad_dlist), pointer :: klx
      type (sad_string), pointer :: str
      character(len=1024) :: pkg
      character(len=64) :: env
      integer :: lpkg, lenv
      integer*8 ktfsymbolc,ktrvaloc,ktcontaloc,
     $     iaxsys,loc,ktcvaloc,kax,k1,k2,i
      integer*4 lpw,lenw,ifromstr
      call tfinfinit
      kinfinity=transfer(dinfinity,i00)
      kminfinity=transfer(-dinfinity,i00)
      knotanumber=transfer(dnotanumber,i00)
      call tfsinglechar
      levele=1
      itflocal=ktaloc(maxlevele+1)
      do i=itflocal,itflocal+maxlevele
        klist(i)=i
      enddo
      levelp=maxlevele
      lgeneration=0
      ltrace=0
      iordless=0
      modethrow=-1
      ncprolog=0
      initmessage=1
      itfcontroot=ktcontaloc(i00)
      call tfassigncont(itfcontroot,'`')
      call tfassigncont(itfcontroot,'Global`')
      iaxsys=ktfsymbolc('`System`',8,itfcontroot)
      call loc_symdef(iaxsys,contd)
      contd%sym%attr=iattrprotected
      itfcontext=ktfaddr(contd%value%k)
      itfcontextpath=ktaloc(2)
      ilist(2,itfcontextpath-1)=2
      klist(itfcontextpath)=itfcontext
      klist(itfcontextpath+1)=itfcontroot

      call tfdefop
      call tfdeffun
      dxnulll=kxaaloc(0,0)
      kxnulll=dxnulll%k
      dxnull=kxaaloc(0,0,klx)
      klx%head=dxnullo
      kxnull=dxnull%k
      dxnulls=kxsalocb(0,'  ',2,str)
      str%str(1:2)=char(0)//char(0)
      str%nch=0
      kxnulls=dxnulls%k
      iaxslotnull=ktaloc(nslots*2)
      k1=ktaaloc(0,0)
      klist(k1)=ktfoper+mtfslot
      k2=ktaaloc(0,0)
      klist(k2)=ktfoper+mtfslotseq
      klist(iaxslotnull  )=ktflist+k1
      klist(iaxslotnull+1)=ktflist+k2
      do i=2,nslots
        k1=ktavaloc(0,1)
        klist(k1)=ktfoper+mtfslot
        k2=ktavaloc(0,1)
        klist(k2)=ktfoper+mtfslotseq
        rlist(k1+1)=i
        rlist(k2+1)=rlist(k1+1)
        klist(iaxslotnull+(i-1)*2  )=ktflist+k1
        klist(iaxslotnull+(i-1)*2+1)=ktflist+k2
      enddo
      iaxslotpart=ktaloc(nslots+3)
      do i=-2,nslots
        dlist(iaxslotpart+i+2)=kxadaloc(0,2,klx)
        klx%head%k=ktfoper+mtfpart
        klx%dbody(1)=dtfcopy1(dlist(iaxslotnull))
        klx%rbody(2)=dble(i)
      enddo
      call tsvaloc('$Version',versionid)
      call tsvaloc('$CreationDate',versiondate)
      loc=-1
      dxliteral=kxsymbolz('Literal',7)
      dxeof=kxsymbolz('EndOfFile',9)
      iaxout=ktfsymbolz('Out',3)
      dxfailed=kxsymbolz('$Failed',7)
      kxmatrix=kxsymbolz('$Matrix',7)
      dxvect=kxsymbolz('Vector',6)
      dxvect1=kxsymbolf('Vector$',7,.true.)
      kxliteral=dxliteral%k
      kxeof=dxeof%k
      kxfailed=dxfailed%k
      kxvect=dxvect%k
      kxvect1=dxvect1%k
      call descr_symdef(kxsymbolz('System`$ReduceMath',18),
     $     redmath)
      
      kax=ktfsymbolz('SEED',4)

      iaxline=ktrvaloc('$Line',0.d0)
      iaximmediate=ktrvaloc('$Immediate',1.d0)
      iaxpriority=ktrvaloc('$Priority',0.d0)
      ierrorgen=ktrvaloc('$ErrorCount',0.d0)
      levelcompile=ktrvaloc('$CompileCount',1.d0)
      kax=ktrvaloc('$ReduceMath',0.d0)

c     Logical symbol
      kax=ktrvaloc('True',1.d0)
      kax=ktrvaloc('False',0.d0)

c     Math symbol
      kax=ktcvaloc('I',0.d0,1.d0)
      kax=ktrvaloc('NaN',dnotanumber)
      kax=ktrvaloc('INF',dinfinity)
      kax=ktrvaloc('Infinity',dinfinity)

c     Math constant
      kax=ktrvaloc('Pi',pi)
      kax=ktrvaloc('E',napier)
      kax=ktrvaloc('EulerGamma',euler)

c     Physical constant
      kax=ktrvaloc('SpeedOfLight',cveloc)
      kax=ktrvaloc('MKSAMu0',mu0)
      kax=ktrvaloc('MKSAEpsilon0',ep0)
      kax=ktrvaloc('SIMu0',mu0)
      kax=ktrvaloc('SIEpsilon0',ep0)
      kax=ktrvaloc('ElectronCharge',elemch)
      kax=ktrvaloc('FineStructureConstant',finest)
      kax=ktrvaloc('ElectronMass',elmass)
      kax=ktrvaloc('ElectronRadius',elradi)
      kax=ktrvaloc('ProtonMass',prmass)
      kax=ktrvaloc('ProtonRadius',prradi)
      kax=ktrvaloc('BoltzmannConstant',kboltzman)
      kax=ktrvaloc('ElectronGminus2over2',gspin)
      kax=ktrvaloc('PlanckConstant',plank)
      kax=ktrvaloc('PlanckHbar',plankr)

      ierrorth=0
      ierrorexp=0

c      write(*,*)'unlink ',shm_unlink('/bar'//char(0))
c      kshm=shm_map('/bar'//char(0),iunit,irtc)
c      write(*,*)'shm: ',iunit,irtc,kshm,rlist(kshm/8)
c      rlist(kshm/8)=m_pi/2
c      write(*,*)'shm: ',rlist(kshm/8)

      call get_environment_variable('COLUMNS',env)
      if(env .eq. ' ')then
        env='132'
      endif
      lpw=min(255,max(79,ifromstr(env(1:lenw(env)))-1))
      iavpw=ktrvaloc('System`PageWidth',dble(lpw))

      call get_environment_variable('SAD_PACKAGES',pkg)
      if(pkg .eq. ' ')then
        call get_environment_variable('SAD$PACKAGES',pkg)
        if(pkg .eq. ' ')then
          call buildinfo_get_string('Target:SAD_PKG_ROOT', pkg)
          if(pkg .eq. ' ')then
            pkg='/SAD/share/Packages/'
          endif
        endif
      endif
      call texpfn(pkg)
      lpkg=len_trim(pkg)
      if(pkg(lpkg:lpkg) .ne. '/')then
        lpkg=lpkg+1
        pkg(lpkg:lpkg)='/'
      endif

      call get_environment_variable('SAD_ENV',env)
      if(.not. (len_trim(env) .gt. 0 .and. index(env,'/') .ne. 1))then
        call buildinfo_get_string('Target:SAD_ENV', env)
        if(.not. (len_trim(env) .gt. 0 .and. index(env,'/') .ne. 1))then
          env='local'
        endif
      endif
c     index(env,'/') MUST be `0' or grater than `1'
      lenv=index(env,'/')-1
      if(lenv .le. 0)then
        lenv=len_trim(env)
      endif

      write(*,*) '*** SADScript Initialization: '//
     $     pkg(1:lpkg)//'init.n ***'
      lfno=0
      call tfgetf(pkg(1:lpkg)//'init.n')
c      write(*,*)'tfinitn 1 '
      klist(itfcontextpath)=itfcontroot
      klist(itfcontextpath+1)=itfcontext
      itfcontext=itfcontroot
      write(*,*) '*** Run time Environment:     '//
     $     pkg(1:lpkg)//'init.'//env(1:lenv)//'.n ***'
      call tfgetf(pkg(1:lpkg)//'init.'//env(1:lenv)//'.n')
c      write(*,*)'tfinitn-9 ',itfcontroot
      return
      end

      subroutine tfassigncont(kp,name)
      use tfstk
      implicit none
      type (sad_symdef), pointer :: contd
      integer*8 ,intent(in):: kp
      integer*8 ka,ktsydefc
      character*(*) ,intent(in):: name
      ka=ktsydefc(name,len(name),itfcontroot,.true.)
      call loc_symdef(ka,contd)
      call tflocald(contd%value)
      contd%sym%gen=-3
      contd%value%k=ktflist+ktfcopy1(kp)
      contd%sym%attr=iattrprotected
      if(klist(kp) .eq. 0)then
        klist(kp)=ktfsymbol+ktfcopy1(ka)
      endif
      return
      end

      subroutine tfsinglechar
      use tfstk
      implicit none
      type (sad_string), pointer :: str
      integer*8 j
      integer*4 i
      iaxschar=ktaloc(1280)
      do i=0,255
        j=iaxschar+5*i+3
        call loc_string(j,str)
        str%len=5
        str%override=-1
        str%alloc%k=ktfstring
        str%ref=1
        str%gen=0
        str%nch=1
        str%nc=0
        str%kstr(1)=0
        str%str(1:1)=char(i)
      enddo
      return
      end

      subroutine tfcontext(isp1,kx,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx,k
      type (sad_symbol),pointer :: sym
      type (sad_namtbl),pointer :: loc
      integer*8 ka
      integer*4 isp1,irtc,itfmessage
      if(isp1+1 .ne. isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      k=dtastk(isp)
      if(ktfoperq(k,ka))then
        call loc_sad(klist(ifunbase+ka),sym)
      elseif(.not. ktfsymbolq(k,sym))then
        irtc=itfmessage(9,'General::wrongtype','"Symbol"')
        return
      endif
      call sym_namtbl(sym,loc)
      kx=dlist(loc%cont)
      irtc=0
      return
      end

      subroutine tftocontext(isp1,kx,irtc)
      use tfstk
      use efun
      use funs
      use eeval
      implicit none
      type (sad_descriptor) kx,k,kc,ki,ks,ic
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      type (sad_dlist), pointer :: kl
      integer*8 kai,ktsydefc,ktfsymbolc,ktcontaloc,kaopt(1)
      integer*4 isp1,irtc,itfmessage,i,isp0,nc,ispopt,isp2
      character*4 optname(1)
      save kaopt,optname
      data kaopt /0/
      data optname /'Wrap'/
      isp0=isp
      call tfgetoptionstk(isp1,kaopt,optname,1,ispopt,irtc)
      ic%k=itfcontroot
      ks=dlist(itfcontroot)
      isp2=ispopt-1
      LOOP_I: do i=isp1+1,isp2
        ki=dtastk(i)
        if(ktfsymbolq(ki,sym))then
          call sym_symstr(sym,str)
        elseif(ktfstringq(ki,str))then
        elseif(ktfoperq(ki,kai))then
          if(kai .eq. mtfnull)then
            if(i .eq. isp1+1)then
              ic%k=0
            endif
            cycle LOOP_I
          endif
          call loc_symstr(klist(klist(ifunbase+kai)),str)
        elseif(ktflistq(ki,kl))then
          if(kl%head%k .eq. ktfoper+mtfslot)then
            k=tfslot(mtfslot,kl,.false.,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(ktfsymbolq(k,sym))then
              call sym_symstr(sym,str)
            elseif(ktfoperq(ki,kai))then
              if(kai .eq. mtfnull)then
                if(i .eq. isp1+1)then
                  ic%k=0
                endif
                cycle LOOP_I
              endif
              call loc_symstr(klist(klist(ifunbase+kai)),str)
            elseif(.not. ktfstringq(k,str))then
              go to 9000
            endif
          else
            go to 9000
          endif
        else
          go to 9000
        endif
        nc=str%nch
c        write(*,*)'tftocontext ',str%str(1:nc)
        if(i .ne. isp2 .and. ic%k .ne. 0)then
          if(str%str(nc:nc) .ne. '`')then
            str%str(nc+1:nc+1)='`'
            ks%k=ktfsymbol+ktsydefc(str%str,nc+1,ic%k,.true.)
            str%str(nc+1:nc+1)=char(0)
          else
            ks%k=ktfsymbol+ktsydefc(str%str,nc,ic%k,.true.)
          endif
          call descr_sad(ks,symd)
          symd%sym%gen=-3
          kc=symd%value
          if(ktfnonlistq(kc))then
            call tflocald(kc)
            ic%k=ktcontaloc(ktfaddrd(ks))
          else
            ic%k=ktfaddrd(kc)
          endif
        else
          ks%k=ktfsymbol+ktfsymbolc(str%str,nc,ic%k)
        endif
      enddo LOOP_I
      if(ispopt .gt. isp0)then
        isp=isp0
        kx=tfsyeval(ks,irtc)
      else
        isp=isp+1
        dtastk(isp)=ks
        kx=tfefunref(isp0+1,.true.,irtc)
        isp=isp0
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Symbol or Character-string"')
      return
      end

      subroutine tfsetcontext(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_symdef), pointer :: symd
      integer*4 isp1,irtc,itfmessage
      logical*4 tfcontextqk
      if(isp1+1 .ne. isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      kx=dtastk(isp)
      if(.not. tfcontextqk(kx%k))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Context (Symbol ending `)"')
        return
      endif
      call loc_sad(ktfaddr(kx),symd)
      itfcontext=ktfaddrd(symd%value)
      irtc=0
      return
      end

      subroutine tfsetcontextpath(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ki
      type (sad_dlist), pointer :: kl
      type (sad_symdef), pointer :: symd
      integer*8 ka1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,m,i,isp0
      logical*4 tfcontextqk
      if(isp1+1 .ne. isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      kx=dtastk(isp)
c      call tfdebugprint(kx,'setcpath',1)
      if(.not. tflistq(kx,kl) .or. ktfreallistq(kl))then
        go to 9000
      endif
      m=kl%nl
      isp0=isp
      do i=1,m
        ki=kl%dbody(i)
        if(.not. tfcontextqk(ki%k))then
          isp=isp0
          go to 9000
        endif
        call descr_sad(ki,symd)
        isp=isp+1
        ktastk(isp)=ktfaddrd(symd%value)
      enddo
      ka1=ktaloc(m)
      ilist(2,ka1-1)=m
      klist(ka1:ka1+m-1)=ktastk(isp0+1:isp0+m)
c      call tmov(ktastk(isp0+1),klist(ka1),m)
      isp=isp0
      call tfree(itfcontextpath)
      itfcontextpath=ka1
c      call tfdebugprint(kx,'setcontextpath',1)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"List of Contexts"')
      return
      end

      subroutine tfexit(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage,istatus
      if(isp .eq. isp1)then
        call tfresetsharedmap()
        call exit(0)
      elseif(isp .eq. isp1+1)then
        if(ktfrealq(ktastk(isp)))then
          call tfresetsharedmap()
          istatus=int(rtastk(isp1+1))
          call exit(istatus)
        elseif(ktastk(isp) .eq. ktfoper+mtfnull)then
          call tfresetsharedmap()
          call exit(0)
        endif
      endif
      kx%k=ktfoper+mtfnull
      irtc=itfmessage(9,'General::wrongtype','"Null"')
      return
      end

      subroutine tfhead(k,kx)
      use tfstk
      implicit none
      type (sad_descriptor) kxreal,kxsymbol,kxstring,kxpat,k,kx
      type (sad_symdef), pointer :: symd
      save kxreal,kxsymbol,kxstring,kxpat
      data kxreal%k /0/
      if(kxreal%k .eq. 0)then
        kxreal=kxsymbolz('Real',4)
        kxstring=kxsymbolz('String',6)
        kxpat%k=ktfoper+mtfcolon
        kxsymbol=kxsymbolz('Symbol',6,symd)
        kxsymbol=symd%value
      endif
      if(ktfrealq(k))then
        kx=kxreal
      else
        select case (ktftype(k%k))
        case (ktflist)
          kx=dlist(ktfaddr(k))
        case (ktfstring)
          kx=kxstring
        case (ktfsymbol)
          kx=kxsymbol
        case (ktfpat)
          kx=kxpat
        case (ktfoper)
          kx=kxsymbol
        end select
      endif
      return
      end

      function tfgetcommandline(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 narg,iargc,i,l,isp0,itfmessage
      character*256 arg
      kx=dxnullo
      if(isp .gt. isp1+1)then
        go to 9000
      elseif(isp .eq. isp1+1)then
        if(ktastk(isp) .ne. ktfoper+mtfnull)then
          go to 9000
        endif
      endif
      narg=iargc()
      irtc=0
      if(narg .le. 0)then
        kx=dxnulll
        return
      endif
      isp0=isp
      do i=1,narg
        call getarg(i,arg)
        l=len_trim(arg)
        isp=isp+1
        ktastk(isp)=ktfstring+ktsalocb(-1,arg,l)
      enddo
      kx=kxmakelist(isp0)
      isp=isp0
      return
 9000 irtc=itfmessage(9,'general::narg','"0"')
      return
      end

      integer*4 function itfunaloc(name,id,narg,map,ieval,immed)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_funtbl), pointer :: fun
      integer*8 k,ktalocr,loc,ktlookupc
      integer*4 ,intent(in):: id,narg,immed
      integer*4 maxnfun,i,nfun
      parameter (maxnfun=2048)
      character*(*) ,intent(in):: name
      logical*4 ev,nev
      integer*4 ,intent(in):: map(narg+1),ieval(narg+1)
      data nfun/0/
      if(nfun .eq. 0)then
        ifunbase=ktalocr(maxnfun)
      endif
      loc=ktlookupc(name,len(name),itfcontext,.true.)
      k=ktfaddr(kxnaloc(0,loc,max(narg,0)+2))
      call c_f_pointer(c_loc(klist(k-9)),fun)
      call tflocald(fun%def%value)
      fun%def%value%k=ktfoper+nfun
      fun%def%sym%attr=iattrconstant
      if(immed .gt. 0)then
        fun%def%sym%attr=fun%def%sym%attr+iattrimmediate
        if(immed .gt. 1)then
          fun%def%sym%attr=fun%def%sym%attr+iattrnumeric
        endif
      endif
      fun%def%sym%override=-2
      fun%narg=narg
      fun%id=id
      nev=.true.
      ev=.true.
      do i=1,narg+1
        fun%mapeval(1,i)=map(i)
        fun%mapeval(2,i)=ieval(i)
        nev=nev .and. map(i) .eq. 0 .and. ieval(i) .eq. 1
        ev=ev .and. ieval(i) .eq. 0
      enddo
      if(nev)then
        fun%mapeval(2,1)=-1
      elseif(ev)then
        fun%mapeval(2,1)=-2
      endif
      klist(ifunbase+nfun)=ktfcopy1(k)
      itfunaloc=nfun
      nfun=nfun+1
      if(nfun .gt. maxnfun)then
        write(*,*)'itfunaloc-too many functions.'
        call abort
      endif
      return
      end

      subroutine tfdefop
      use tfstk
      implicit none
      integer*4 i,itfunaloc,map(32),ieval(32)
c     Initialize map/ieval array
      map=0
      ieval=0
      i=itfunaloc('Null',-mtfnull,1,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     -iattrconstant
      i=itfunaloc('$f1$',-1,-1,map,ieval,0)
      i=itfunaloc('$f2$',-2,-1,map,ieval,0)
      i=itfunaloc('Plus',-mtfplus,1,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     +iattrorderless
      i=itfunaloc('$f4$',-4,-1,map,ieval,0)
      i=itfunaloc('Times',-mtftimes,1,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     +iattrorderless
      i=itfunaloc('$f6$',-6,-1,map,ieval,0)
      i=itfunaloc('InversePower',-mtfrevpower,1,map,ieval,2)
      i=itfunaloc('Power',-mtfpower,1,map,ieval,2)
      i=itfunaloc('Greater',-mtfgreater,1,map,ieval,2)

      i=itfunaloc('GreaterEqual',-mtfgeq,1,map,ieval,2)
      i=itfunaloc('LessEqual',-mtfleq,1,map,ieval,2)
      i=itfunaloc('Less',-mtfless,1,map,ieval,2)
      i=itfunaloc('Equal',-mtfequal,1,map,ieval,2)
      i=itfunaloc('Unequal',-mtfunequal,1,map,ieval,2)
      i=itfunaloc('And',-mtfand,-1,map,ieval,2)
      i=itfunaloc('Or',-mtfor,-1,map,ieval,2)
      i=itfunaloc('Not',-mtfnot,1,map,ieval,2)
      i=itfunaloc('SameQ',-mtfsame,1,map,ieval,2)
      i=itfunaloc('UnsameQ',-mtfunsame,1,map,ieval,2)

      i=itfunaloc('StringJoin',-mtfconcat,1,map,ieval,2)
      i=itfunaloc('$f21$',-mtfleftbra,-1,map,ieval,0)
      i=itfunaloc('$f22$',-mtfrightbra,-1,map,ieval,0)
      i=itfunaloc('List',-mtflist,-1,map,ieval,0)
      i=itfunaloc('$f24$',-mtfrightbrace,-1,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('SetDelayed',-mtfsetdelayed,2,map,ieval,0)
      i=itfunaloc('Set',-mtfset,-1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Complex',-mtfcomplex,1,map,ieval,2)
      i=itfunaloc('$f28$',-mtfleftparen,-1,map,ieval,0)
      i=itfunaloc('$f29$',-mtfrightparen,-1,map,ieval,0)

      i=itfunaloc('$f30$',-mtfcomma,-1,map,ieval,0)
      i=itfunaloc('CompoundExpression',-mtfcomp,-1,map,ieval,0)
      i=itfunaloc('Function',-mtffun,-1,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Pattern',-mtfcolon,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Rule',-mtfrule,-1,map,ieval,0)
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('RuleDelayed',-mtfruledelayed,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('ReplaceAll',-mtfreplace,2,map,ieval,0)
      i=itfunaloc('ReplaceRepeated',-mtfreplacerepeated,2,map,ieval,0)
      ieval(1)=1
      i=itfunaloc('UpSet',-mtfupset,2,map,ieval,0)
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('UpSetDelayed',-mtfupsetdelayed,2,map,ieval,0)

      i=itfunaloc('Unset',-mtfupsetdelayed,1,map,ieval,0)
      i=itfunaloc('PatternTest',-mtfpattest,-1,map,ieval,0)
      i=itfunaloc('FFSFlag',-mtfflag,1,map,ieval,0)
      i=itfunaloc('Slot',-mtfslot,-1,map,ieval,0)
      i=itfunaloc('SlotSequence',-mtfslotseq,-1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Dot',-mtfdot,2,map,ieval,2)
      i=itfunaloc('Alternatives',-mtfalt,1,map,ieval,0)
      i=itfunaloc('Map',-mtfmap,1,map,ieval,0)
      i=itfunaloc('MapAll',-mtfmapall,1,map,ieval,0)
      i=itfunaloc('Apply',-mtfapply,1,map,ieval,0)

      i=itfunaloc('Repeated',-mtfrepeated,-1,map,ieval,0)
      i=itfunaloc('RepeatedNull',-mtfrepeatednull,-1,map,ieval,0)
      i=itfunaloc('Inequality',-mtfinequality,1,map,ieval,2)
      ieval(1)=1
      i=itfunaloc('AddTo',-mtfaddto,2,map,ieval,0)
      i=itfunaloc('SubtractFrom',-mtfsubtractfrom,2,map,ieval,0)
      i=itfunaloc('TimesBy',-mtftimesby,2,map,ieval,0)
      i=itfunaloc('DivideBy',-mtfdivideby,2,map,ieval,0)
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Increment',-mtfincrement,2,map,ieval,0)
      i=itfunaloc('Decrement',-mtfdecrement,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Part',-mtfpart,-1,map,ieval,0)

      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Member',-mtfatt,2,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('MessageName',-mtfmessagename,1,map,ieval,0)
      ilist(1,klist(ifunbase+i)-3)=iattrholdall
      i=itfunaloc('TagSet',-mtftagset,-1,map,ieval,0)
      i=itfunaloc('$f63$',-mtfleftcomment,-1,map,ieval,0)
      i=itfunaloc('$f64$',-mtfrightcomment,-1,map,ieval,0)
      i=itfunaloc('Hold',-mtfhold,-1,map,ieval,0)
      i=itfunaloc('$f66$',-mtfend,-1,map,ieval,0)

      return
      end

      subroutine tfdeffun
      implicit none
c
c     id <= 1000
      call tfdeffun0
c     1000< id <=2000
      call tfdeffun1
c     2000< id <=3000 for EPICS. defined in the file tfefun2.f 
      call tfdeffun2
c     3000< id <=4000 for EPICS. defined in the file tfefun3ep.f
      call tfdeffun3ep
c     Call sadDefFunc_* collection
      call tfdeffuncs()
c
      return
      end
c
      subroutine tfdeffun0
      use tfstk
      implicit none
      integer*4 i,itfunaloc,map(32),ieval(32)
c     Initialize map/ieval array
      do i=1,32
         map(i)=0
         ieval(i)=0
      enddo
      map(1)=1
      i=itfunaloc('Sin',1,1,map,ieval,2)
      i=itfunaloc('Cos',2,1,map,ieval,2)
      i=itfunaloc('Tan',3,1,map,ieval,2)
      i=itfunaloc('Sinh',4,1,map,ieval,2)
      i=itfunaloc('Cosh',5,1,map,ieval,2)
      i=itfunaloc('Tanh',6,1,map,ieval,2)
      i=itfunaloc('Exp',7,1,map,ieval,2)
      i=itfunaloc('Log',8,1,map,ieval,2)
      map(1)=0
      i=itfunaloc('ArcTan',9,1,map,ieval,2)
      i=itfunaloc('Det',10,1,map,ieval,2)
      map(1)=1
      i=itfunaloc('Sqrt',11,1,map,ieval,2)
      i=itfunaloc('Floor',12,1,map,ieval,2)
      i=itfunaloc('Ceiling',13,1,map,ieval,2)
      map(1)=0
      i=itfunaloc('Min',14,0,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     +iattrorderless
      i=itfunaloc('Max',15,0,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     +iattrorderless
      i=itfunaloc('Mod',16,2,map,ieval,2)
      i=itfunaloc('StringLength',17,1,map,ieval,2)
      i=itfunaloc('Arg',18,1,map,ieval,2)
      map(1)=1
      i=itfunaloc('Sign',19,1,map,ieval,2)
      map(1)=0
      i=itfunaloc('Length',nfunlength,1,map,ieval,2)
      i=itfunaloc('Dimensions',21,1,map,ieval,2)
      i=itfunaloc('ReplacePart',nfunreppart,3,map,ieval,0)
      map(1)=1
      i=itfunaloc('ArcSin',23,1,map,ieval,2)
      i=itfunaloc('ArcCos',24,1,map,ieval,2)
      i=itfunaloc('ArcSinh',25,1,map,ieval,2)
      i=itfunaloc('ArcCosh',26,1,map,ieval,2)
      i=itfunaloc('ArcTanh',27,1,map,ieval,2)
      map(1)=0
      ieval(1)=1
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Table',28,2,map,ieval,2)
      i=itfunaloc('Do',nfundo,2,map,ieval,1)
      i=itfunaloc('Attributes',30,1,map,ieval,0)
      ieval(1)=0
      i=itfunaloc('Peek',31,1,map,ieval,0)
      map(1)=1
      i=itfunaloc('Abs',32,1,map,ieval,2)
      map(1)=0
      ieval(2)=0
      i=itfunaloc('Reverse',33,1,map,ieval,2)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('Module',nfunmodule,2,map,ieval,0)
      i=itfunaloc('Block',nfunblock,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('StringReplace',36,2,map,ieval,1)
      i=itfunaloc('SwitchCases',nfunswicases,2,map,ieval,0)
      i=itfunaloc('Flatten',38,2,map,ieval,1)
      ieval(2)=1
      ieval(3)=1
      ieval(4)=1
      ieval(5)=1
      i=itfunaloc('If',nfunif,4,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('Take',40,2,map,ieval,1)
      ieval(3)=0
      i=itfunaloc('Select',nfunselect,3,map,ieval,0)
      ieval(2)=1
      ieval(1)=1
      i=itfunaloc('While',42,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Join',43,2,map,ieval,1)
      i=itfunaloc('Append$',nfunappend,2,map,ieval,1)
      i=itfunaloc('Prepend$',nfunprepend,2,map,ieval,1)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('Clear',46,1,map,ieval,0)
      i=itfunaloc('Protect',47,1,map,ieval,0)
      i=itfunaloc('Unprotect',48,1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('Drop',49,2,map,ieval,1)
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('MapAt',50,3,map,ieval,0)
      i=itfunaloc('Inner',51,4,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Transpose',52,1,map,ieval,1)
      i=itfunaloc('SingularValues1',53,3,map,ieval,1)
      i=itfunaloc('DiagonalMatrix',54,1,map,ieval,1)
      i=itfunaloc('LinearSolveReal',55,3,map,ieval,1)
      i=itfunaloc('IdentityMatrix',56,1,map,ieval,1)
      i=itfunaloc('Eigensystem',57,1,map,ieval,1)
      i=itfunaloc('Operation',58,1,map,ieval,1)
      i=itfunaloc('Position',nfunposition,4,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Sum',60,2,map,ieval,1)
      i=itfunaloc('Product',61,2,map,ieval,1)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('Range$',62,3,map,ieval,2)
      map(1)=1
      i=itfunaloc('Re',63,1,map,ieval,2)
      i=itfunaloc('Im',64,1,map,ieval,2)
      i=itfunaloc('Conjugate',65,1,map,ieval,2)
      map(1)=0
      ieval(1)=0
      i=itfunaloc('ToString',66,2,map,ieval,0)
      i=itfunaloc('Depth',67,1,map,ieval,1)
      i=itfunaloc('Level',68,2,map,ieval,1)
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Write',69,2,map,ieval,0)
      i=itfunaloc('Get',70,1,map,ieval,0)
      i=itfunaloc('OpenWrite',71,1,map,ieval,0)
      i=itfunaloc('OpenAppend',72,1,map,ieval,0)
      i=itfunaloc('Close$',73,1,map,ieval,0)
      i=itfunaloc('Flush',74,1,map,ieval,0)
      ieval(1)=1
      i=itfunaloc('Print',75,1,map,ieval,0)
      ieval(1)=0
      i=itfunaloc('WriteString',76,2,map,ieval,0)
      ieval(2)=0
      ieval(3)=0
      i=itfunaloc('Return',77,1,map,ieval,0)
      i=itfunaloc('Head',78,1,map,ieval,1)
      i=itfunaloc('RealListQ',79,1,map,ieval,0)
      i=itfunaloc('Partition',80,3,map,ieval,1)
      i=itfunaloc('Throw',81,1,map,ieval,0)
      ieval(1)=1
      i=itfunaloc('Catch',82,1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('Thread',nfunthread,2,map,ieval,1)
      ieval(1)=1
      i=itfunaloc('SetAttributes',84,2,map,ieval,0)
      ieval(1)=0
      i=itfunaloc('MapIndexed',85,3,map,ieval,0)
      i=itfunaloc('FromCharacterCode',86,1,map,ieval,1)
      i=itfunaloc('ToCharacterCode',87,1,map,ieval,1)
      i=itfunaloc('ComplexQ',88,1,map,ieval,2)
      ieval(2)=0
      i=itfunaloc('Tr',89,1,map,ieval,2)
c      i=itfunaloc('SeedRandom',89,1,map,ieval,0)
      i=itfunaloc('SaveSharedMap',90,0,map,ieval,0)
c      i=itfunaloc('GaussRandom',90,1,map,ieval,0)
      ieval(2)=1
      i=itfunaloc('Switch',nfunswitch,1,map,ieval,1)
      ieval(2)=0
      i=itfunaloc('Sort',92,2,map,ieval,1)
      i=itfunaloc('Union1',93,2,map,ieval,1)
      ieval(2)=0
      i=itfunaloc('Order',94,2,map,ieval,1)
      ieval(1)=0
      i=itfunaloc('MemoryCheck',95,0,map,ieval,0)
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('Scan',nfunscan,3,map,ieval,0)
      i=itfunaloc('Identity',97,1,map,ieval,2)
      i=itfunaloc('TimeUsed',98,0,map,ieval,0)
      i=itfunaloc('NumberQ',99,0,map,ieval,1)
      i=itfunaloc('VectorQ',100,2,map,ieval,1)
      i=itfunaloc('AtomQ',101,1,map,ieval,1)
      i=itfunaloc('Outer',102,3,map,ieval,0)
      i=itfunaloc('MatchQ',103,2,map,ieval,0)
      ieval(1)=1
      ieval(2)=0
      i=itfunaloc('TracePrint',104,1,map,ieval,0)
      ieval(2)=1
      i=itfunaloc('Definition',105,1,map,ieval,0)
c-----Noboru addition     -----
      ieval(1)=0
      ieval(2)=0
      map(1)=0
      map(2)=0
      i=itfunaloc('READ',106,1,map,ieval,0)
c-----Noboru addition end -----
      i=itfunaloc('Intersection',107,1,map,ieval,1)
      i=itfunaloc('Complement',108,1,map,ieval,1)
      map(1)=1
      i=itfunaloc('Round',109,1,map,ieval,2)
      map(1)=0
      i=itfunaloc('InverseErf',110,1,map,ieval,2)
c      i=itfunaloc('SemCtrl',111,3,map,ieval,0)
c      i=itfunaloc('FromDate',111,1,map,ieval,1)
c      i=itfunaloc('ToDate',112,1,map,ieval,1)
      i=itfunaloc('ToInputString',113,1,map,ieval,0)
      ieval(3)=0
      i=itfunaloc('ReadString',114,3,map,ieval,0)
      i=itfunaloc('OpenRead',115,1,map,ieval,0)
      i=itfunaloc('ToExpression',116,1,map,ieval,0)
      i=itfunaloc('StringMatchQ',117,2,map,ieval,2)
      i=itfunaloc('StringPosition',118,3,map,ieval,2)
      i=itfunaloc('ToUpperCase',119,1,map,ieval,2)
      i=itfunaloc('Break',120,0,map,ieval,0)
      i=itfunaloc('Continue',121,0,map,ieval,0)
      i=itfunaloc('Goto',122,1,map,ieval,0)
      i=itfunaloc('Fourier',123,1,map,ieval,1)
      i=itfunaloc('InverseFourier',124,1,map,ieval,1)
      ieval(1)=1
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('Check',125,2,map,ieval,0)
      i=itfunaloc('Which',nfunwhich,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('MapFile',127,2,map,ieval,0)
      i=itfunaloc('UnmapFile',128,1,map,ieval,0)
c      i=itfunaloc('GetUID',129,1,map,ieval,0)
c      i=itfunaloc('GetGID',130,1,map,ieval,0)
      i=itfunaloc('ToLowerCase',131,1,map,ieval,2)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('Unevaluated$',nfununeval,1,map,ieval,1)      
      ieval(1)=0
      ieval(2)=0
      ieval(3)=0
      ieval(4)=0
      ieval(5)=0
      i=itfunaloc('Cases',nfuncases,2,map,ieval,0)      
      i=itfunaloc('DeleteCases',nfundelcases,2,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('Vectorize',135,1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
c      i=itfunaloc('BDPipe',136,1,map,ieval,0)
c      i=itfunaloc('BidirectionalPipe',136,1,map,ieval,0)
      i=itfunaloc('Names$',137,1,map,ieval,0)      
      i=itfunaloc('GarbageCollect',138,0,map,ieval,0)      
      i=itfunaloc('LogGamma1',139,1,map,ieval,2)
      i=itfunaloc('Factorial',140,1,map,ieval,2)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('With',nfunwith,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('SelectCases',nfunselcases,2,map,ieval,0)
      i=itfunaloc('Override',143,1,map,ieval,1)
      ieval(1)=1
      i=itfunaloc('AppendTo',144,2,map,ieval,0)
      i=itfunaloc('PrependTo',145,2,map,ieval,0)
      ieval(2)=1
      i=itfunaloc('FindRoot$',146,2,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('GammaRegularized',147,2,map,ieval,2)
      i=itfunaloc('GammaRegularizedQ',147,2,map,ieval,2)
      i=itfunaloc('GammaRegularizedP',148,2,map,ieval,2)
      map(1)=1
      i=itfunaloc('Erf',149,1,map,ieval,2)
      i=itfunaloc('Erfc',150,1,map,ieval,2)
      map(1)=0
      ieval(2)=1
      ieval(3)=1
      ieval(4)=1
      i=itfunaloc('Fit$',151,4,map,ieval,0)
      ieval(2)=0
      ieval(3)=0
      ieval(4)=0
      i=itfunaloc('Symbol',152,1,map,ieval,1)
      ieval(1)=1
      i=itfunaloc('SymbolName',153,1,map,ieval,1)
      ieval(1)=0
      i=itfunaloc('Extract',nfunextract,2,map,ieval,1)
      i=itfunaloc('Read',155,3,map,ieval,0)
      i=itfunaloc('Skip',156,3,map,ieval,0)
      i=itfunaloc('TemporaryName',157,0,map,ieval,0)
      i=itfunaloc('Exit',158,0,map,ieval,0)
      i=itfunaloc('StringFill',159,3,map,ieval,2)
      i=itfunaloc('Restrict',160,3,map,ieval,2)
      i=itfunaloc('MinMax',161,1,map,ieval,2)
      ilist(1,klist(ifunbase+i)-3)=ilist(1,klist(ifunbase+i)-3)
     $     +iattrorderless
      i=itfunaloc('Short',162,2,map,ieval,0)
      i=itfunaloc('$SetOMPNumThreads',163,1,map,ieval,0)
c      i=itfunaloc('Directory',164,0,map,ieval,0)
c      i=itfunaloc('SetDirectory',165,1,map,ieval,0)
c      i=itfunaloc('Wait',166,0,map,ieval,0)
      i=itfunaloc('BesselJ',167,2,map,ieval,2)
      i=itfunaloc('BesselY',168,2,map,ieval,2)
      i=itfunaloc('BesselI',169,2,map,ieval,2)
      i=itfunaloc('BesselK',170,2,map,ieval,2)
      i=itfunaloc('BaseForm',171,2,map,ieval,0)
      i=itfunaloc('StringTrim',172,1,map,ieval,2)
      i=itfunaloc('StringToStream',173,1,map,ieval,0)
      i=itfunaloc('EvenQ',174,1,map,ieval,2)
      i=itfunaloc('OddQ',175,1,map,ieval,2)
c      i=itfunaloc('DateString$',176,1,map,ieval,0)
      i=itfunaloc('Insert',nfuninsert,3,map,ieval,1)
      i=itfunaloc('Delete',nfundelete,2,map,ieval,1)
      i=itfunaloc('FlattenAt',179,3,map,ieval,1)
      i=itfunaloc('Replace',180,3,map,ieval,0)
c      i=itfunaloc('SetEnv',181,2,map,ieval,0)
      i=itfunaloc('Spline$',182,2,map,ieval,0)
      i=itfunaloc('FindIndex',183,2,map,ieval,0)
      i=itfunaloc('$SetContext',184,1,map,ieval,0)
      i=itfunaloc('$SetContextPath',185,1,map,ieval,0)
      ieval(1)=1
      ieval(2)=1
      i=itfunaloc('ToContext',186,1,map,ieval,0)
      i=itfunaloc('Context',187,1,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      i=itfunaloc('BitwiseAnd',188,2,map,ieval,2)
      i=itfunaloc('BitwiseOr',189,2,map,ieval,2)
      i=itfunaloc('BitwiseXor',190,2,map,ieval,2)
      ieval(2)=1
      ieval(3)=1
      i=itfunaloc('ReplaceMember',191,2,map,ieval,0)
      ieval(1)=1
      ieval(3)=0
      ieval(4)=1
      i=itfunaloc('MemberScan',192,4,map,ieval,0)
      i=itfunaloc('StandardForm',193,1,map,ieval,0)
      i=itfunaloc('Abort',194,0,map,ieval,0)
      i=itfunaloc('CheckAbort',195,0,map,ieval,0)
      ieval(1)=0
      ieval(2)=0
      ieval(4)=0
      ieval(5)=0
      i=itfunaloc('ReleaseHold',196,1,map,ieval,0)
      i=itfunaloc('NaNQ',197,1,map,ieval,0)
      i=itfunaloc('MapThread',198,2,map,ieval,0)
      i=itfunaloc('ScanThread',199,2,map,ieval,0)
      i=itfunaloc('Last',200,1,map,ieval,2)
      i=itfunaloc('First',201,1,map,ieval,2)
      i=itfunaloc('Second',202,1,map,ieval,2)
      i=itfunaloc('Third',203,1,map,ieval,2)
      i=itfunaloc('ObjectSymbol',204,1,map,ieval,0)
      i=itfunaloc('ProductLog',205,1,map,ieval,2)
      i=itfunaloc('GaussianCoulomb',206,4,map,ieval,2)
      i=itfunaloc('ClearMemberObject',207,1,map,ieval,0)
      i=itfunaloc('MemoryAllocate',208,1,map,ieval,0)
      i=itfunaloc('Duplicate',209,1,map,ieval,0)
      i=itfunaloc('$GetCommandLine',210,0,map,ieval,0)
      i=itfunaloc('SeekFile',211,3,map,ieval,0)
      i=itfunaloc('DigitQ',212,1,map,ieval,0)
      i=itfunaloc('LetterQ',213,1,map,ieval,0)
      i=itfunaloc('RealQ',214,1,map,ieval,0)
      i=itfunaloc('$NearlySameQ',215,4,map,ieval,0)
      i=itfunaloc('OpenShared',216,1,map,ieval,0)
      i=itfunaloc('ReadShared',217,1,map,ieval,0)
      i=itfunaloc('WriteShared',218,1,map,ieval,0)
      i=itfunaloc('SharedSize',219,1,map,ieval,0)
      ieval(1)=1
      i=itfunaloc('FBoundQ$',220,1,map,ieval,0)
      ieval(1)=0
      i=itfunaloc('GaussianCoulombU',221,4,map,ieval,2)
      i=itfunaloc('GaussianCoulombFitted',222,4,map,ieval,2)
      i=itfunaloc('Rest',223,1,map,ieval,1)
      i=itfunaloc('RotateRight1',224,1,map,ieval,1)
      i=itfunaloc('Difference',225,1,map,ieval,1)
      map(1)=1
      i=itfunaloc('Gamma0',226,1,map,ieval,2)
      map(1)=0
c      i=itfunaloc('CloseShared',227,1,map,ieval,0)
      i=itfunaloc('XSin',228,1,map,ieval,2)
      return
      end
c
c id > 1000
c
      subroutine tfdeffun1
      implicit none
      integer*4 i,itfunaloc,map(32),ieval(32)
c     Initialize map/ieval array
      do i=1,32
         map(i)=0
         ieval(i)=0
      enddo
      map(1)=1
      map(2)=1
      map(3)=1
      i=itfunaloc('Element$',1001,2,map,ieval,0)
      i=itfunaloc('Twiss',1002,2,map,ieval,0)
      i=itfunaloc('LINE',1003,2,map,ieval,0)
      map(1)=0
      i=itfunaloc('CalculateEmittance',1004,2,map,ieval,0)
      map(2)=0
      i=itfunaloc('TrackParticles',1005,3,map,ieval,0)
      map(3)=0
      map(4)=0
      map(5)=0
      map(6)=0
      ieval(4)=0
      ieval(5)=0
      ieval(6)=0
      i=itfunaloc('CalculateOptics',1006,5,map,ieval,0)
      i=itfunaloc('DynamicAperture',1007,5,map,ieval,0)
c-----Kikuchi added-----     
      map(1)=1
      map(2)=0
      map(3)=0
      map(4)=0
      map(5)=0
      map(6)=0
      map(7)=0
      map(8)=0
      ieval(7)=0
      ieval(8)=0
      i=itfunaloc('ResponseMatrix',1008,7,map,ieval,0)
      map(1)=0
      i=itfunaloc('Master',1009,0,map,ieval,0)
      map(1)=1
      i=itfunaloc('FLAG',1010,2,map,ieval,0)
      map(1)=0
      i=itfunaloc('LinearSolveConditioned',1011,7,map,ieval,0)
c-----Kikuchi addition end-----     
      map(1)=0
      i=itfunaloc('ExecDifferentialAlgebraPackage',1012,1,map,ieval,0)
      i=itfunaloc('InitializeDifferentialAlgebra',1013,1,map,ieval,0)
      i=itfunaloc('GetDifferentialAlgebraMap',1014,1,map,ieval,0)
      map(2)=0
      map(3)=0
      i=itfunaloc('FFS',1015,2,map,ieval,0)
      i=itfunaloc('RadiationField',1016,2,map,ieval,0)
      i=itfunaloc('RadiationSpectrum',1017,2,map,ieval,0)
      i=itfunaloc('FFSFlags',1018,1,map,ieval,0)
      i=itfunaloc('ExtractBeamLine',1019,1,map,ieval,0)
      i=itfunaloc('BeamLineName',1020,0,map,ieval,0)
      i=itfunaloc('SetElement$',1021,1,map,ieval,0)
      i=itfunaloc('Type$Key$',1022,1,map,ieval,0)
      i=itfunaloc('NormalCoordinate6',1023,1,map,ieval,0)
      i=itfunaloc('InitEmit',1024,0,map,ieval,0)
      i=itfunaloc('FFS$SHOW',1025,0,map,ieval,0)
c      ieval(1)=1
c      ieval(2)=1
      i=itfunaloc('ExpandBeamLine',1026,1,map,ieval,0)
c      ieval(1)=0
c      ieval(2)=0
      i=itfunaloc('GetMAIN',1027,1,map,ieval,0)
c      i=itfunaloc('Canvas3DClipTriangle',1028,1,map,ieval,0)
c      i=itfunaloc('Canvas3DLightTriangle',1029,2,map,ieval,0)
      i=itfunaloc('RGBColor',1030,3,map,ieval,1)
c      i=itfunaloc('Canvas3DProjection',1031,4,map,ieval,0)
c      i=itfunaloc('TclArg1',1032,1,map,ieval,0)
c      i=itfunaloc('CanvasSymbol',1033,1,map,ieval,0)
c      i=itfunaloc('TkCanvasCreateItemDirect',1034,1,map,ieval,0)
c      i=itfunaloc('CanvasSymbolDirect',1035,1,map,ieval,0)
c      i=itfunaloc('TclSetResult',1036,1,map,ieval,0)
      i=itfunaloc('SAD$LifeTrack',1037,1,map,ieval,0)
      i=itfunaloc('SynchroBetaEmittance1',1038,4,map,ieval,0)
      i=itfunaloc('CSRInit',1039,1,map,ieval,0)
      i=itfunaloc('CSRMatrix',1040,1,map,ieval,0)
      i=itfunaloc('CSRConvert',1041,1,map,ieval,0)
      i=itfunaloc('CSRTrack',1042,1,map,ieval,0)
      i=itfunaloc('CSRHaissin',1043,1,map,ieval,0)
      i=itfunaloc('CSRSetupOY',1044,1,map,ieval,0)
      i=itfunaloc('CSROYMatrix',1045,1,map,ieval,0)
      i=itfunaloc('SurvivedParticles',1046,1,map,ieval,0)
      i=itfunaloc('BBBrem1',1047,2,map,ieval,0)
      i=itfunaloc('MapParticles',1048,3,map,ieval,0)
      i=itfunaloc('UnmapParticles',1049,2,map,ieval,0)
      return
      end
