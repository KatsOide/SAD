
      module objsym
      use tfstk
      implicit none
      integer*4 nrecycle
      parameter (nrecycle=2048)
      integer*4 :: irec=0,irectbl(nrecycle)

      contains
      subroutine tfobjectsymbol(isp1,kx,irtc)
      use tfcode
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      integer*4 itfmessage,ls,l,n,n1,ifrac,id,n0
      character*1024 name
      character*32 buf
      character ch
      data id/0/
      if(isp1+1 /= isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonsymbolq(ktastk(isp),sym))then
        irtc=itfmessage(9,'General::wrongtype','"Symbol"')
        return
      endif
      call sym_namtbl(sym,loc)
      ls=min(1000,loc%str%nch)
      name(1:ls)=loc%str%str(1:ls)
      ls=ls+1
      name(ls:ls)="o"
      if(irec /= 0)then
        n0=irectbl(irec)
        irec=irec-1
      else
        id=id+1
        n0=id
      endif
      n=n0
      l=32
      do while(n /= 0)
        n1=n/62
        ifrac=n-n1*62
        if(ifrac .lt. 10)then
          ch=char(ichar('0')+ifrac)
        elseif(ifrac .lt. 36)then
          ch=char(ichar('a')+ifrac-10)
        else
          ch=char(ichar('A')+ifrac-36)
        endif
        buf(l:l)=ch
        l=l-1
        n=n1
      enddo
      name(ls+1:ls+32-l)=buf(l+1:32)
      ls=ls+33-l
      name(ls:ls)="$"
      kx=kxsymbolz(name,ls,symd)
      call sym_namtbl(symd%sym,loc)
      loc%kind=1
      loc%str%nc=n0
      irtc=0
      return
      end

      end module objsym

      module tfcx
      use tfstk
      integer*8, save :: icx=0,ithis,ithisloc,isyscont
      type (sad_descriptor) kxundefined,kxmemberobject,kxmemberobject2
      data kxundefined%k,kxmemberobject%k,kxmemberobject2%k /0,0,0/

      contains
      subroutine tfcxinit
      implicit none
      type (sad_dlist), pointer :: kl
      icx=ktadaloc(0,1,kl)
      kl%head=dtfcopy1(kxsymbolz('Class`cx$',9))
      kl%dbody(1)=dtfcopy1(kxsymbolz('Class`x$',8))
      ithis=ktfsymbolz('System`This',11)
      ithisloc=klist(ithis)
      isyscont=ktfaddr(klist(ktfsymbolz('`System`',8)-4))
      kxmemberobject=kxsymbolz('Class`MemberObject',18)
      kxmemberobject2=kxsymbolz('Class`MemberObject2',19)
      kxundefined   =kxsymbolz('System`Undefined',16)
      return
      end subroutine

      recursive subroutine tfclassmember(k10,k20,kx,eval,irtc)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(in):: k10,k20
      type (sad_descriptor) k1,k2,ks
      type (sad_symdef), pointer :: symd
      type (sad_symbol), pointer :: symh
      type (sad_dlist), pointer :: kl1,klx
      integer*8 ka2,k11,kv
      integer*4 isp1,l
      logical*4 ev
      logical*4 ,intent(in):: eval
      k1=k10
      k2=k20
      if(ktfoperq(k2,ka2))then
        k2%k=ktfsymbol+klist(ifunbase+ka2)
      endif
      if(ktfsymbolq(k2))then
        if(ktflistq(k1,kl1))then
          k11=kl1%head%k
          if(ktfsymbolqdef(k11,symd))then
            if(symd%sym%gen == -3)then
              go to 10
            elseif(symd%sym%override /= 0 .and.
     $             symd%downval /= 0)then
              call tfgetdefargp(kl1,ktfaddr(k11),kv,ev,irtc)
              if(irtc > 0)then
                return
              endif
              if(ev .and. irtc == 0)then
                call tfclassmember(dlist(kv),k2,kx,eval,irtc)
                return
              endif
            endif
          endif
        elseif(ktfsymbolqdef(k1%k,symd) .and.
     $         symd%sym%override /= 0)then
          if(ktflistq(symd%value,kl1) .and.
     $         ktfsymbolq(kl1%head%k,symh))then
            if(symh%gen == -3)then
              go to 10
            endif
          endif
        endif
      endif
      irtc=-1
      return
 10   if(kxmemberobject2%k == 0)then
        call tfcxinit
      endif
      if(kl1%nl /= 1)then
        irtc=-1
        return
      endif
      isp=isp+1
      isp1=isp
      levele=levele+1
      dtastk(isp)=kxmemberobject2
      isp=isp+1
      dtastk(isp)=kl1%dbody(1)
      call tfdeval(isp1,kxmemberobject2,ks,1,.false.,ev,irtc)
      if(irtc > 0)then
        isp=isp1-1
        l=itfdownlevel()
        return
      elseif(ev .and. irtc == 0)then
        dtastk(isp1)=ks
      else
        dtastk(isp1)=kxcompose(isp1)
        ks=kxmemberobject2
      endif
      isp=isp1+1
      dtastk(isp)=kl1%head
      isp=isp+1
      dtastk(isp)=k2
      call tfdeval(isp1,ks,kx,1,.false.,ev,irtc)
c     call tfdebugprint(kx,'classmember',1)
c     write(*,*)'with: ',irtc,ev,eval
      isp=isp1-1
      if(ev)then
        if(irtc /= 0)then
          l=itfdownlevel()
          return
        endif
        if(ktflistq(kx,klx))then
          if(klx%head%k == ktfoper+mtfhold)then
            kx=klx%dbody(1)
          endif
        elseif(kx%k == kxundefined%k)then
          irtc=-1
          l=itfdownlevel()
          return
        endif
        if(eval)then
          if(ktfsymbolq(kx))then
            kx=tfsyeval(kx,irtc)
            if(irtc /= 0)then
              l=itfdownlevel()
              return
            endif
          endif
        endif
        call tfconnect(kx,irtc)
      else
        l=itfdownlevel()
        irtc=-1
      endif
      return
      end

      recursive function tfsolvemember(list,reps,irtc) result(kx)
      use part,only:tfreplist
      use funs
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) k1
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: listx
      integer*4 ,intent(out):: irtc
      logical*4 eval
      logical*4 ,intent(out):: reps
      irtc=-1
      kx=sad_descr(list)
      reps=.false.
      if(list%head%k == ktfoper+mtfatt)then
        if(list%nl == 2 .and. ktfnonreallistqo(list))then
          call tfclassmember(list%dbody(1),list%dbody(2),kx,
     $         .false.,irtc)
        endif
        return
      elseif(list%head%k == ktfoper+mtfslot)then
        kx=tfslot(mtfslot,list,.false.,irtc)
        if(irtc > 0 .and. ierrorprint /= 0)then
          call tfreseterror
          irtc=-1
        endif
        reps=.true.
        return
      endif
      if(ktflistq(list%head,listx))then
        k1=tfsolvemember(listx,reps,irtc)
        if(irtc /= 0)then
          return
        endif
        if(reps)then
          listx=>tfduplist(list)
        else
          listx=>tfclonelist(list)
        endif
        call tfreplist(listx,0,k1,eval)
c        call tfloadlstk(list,listx)
c        listx%head=k1
c        call tfstk2l(listx,listx)
        kx=sad_descr(listx)
      endif
      return
      end

      subroutine tfatt(isp1,kx,eval,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k1
      integer*4 isp0,i
      logical*4 ,intent(in):: eval
      if(isp == isp1+2)then
        call tfclassmember(dtastk(isp1+1),dtastk(isp),kx,eval,irtc)
        return
      endif
      k1=dtastk(isp1+1)
      do i=isp1+2,isp
        call tfclassmember(k1,dtastk(i),k1,eval,irtc)
c        call tfdebugprint(k1,'tfatt',1)
c        write(*,*)'with: ',irtc,i
        if(irtc /= 0)then
          go to 10
        endif
      enddo
      kx=k1
      return
 10   if(irtc > 0)then
        return
      endif
      isp=isp+1
      isp0=isp
      ktastk(isp)=ktfoper+mtfatt
      isp=isp+1
      dtastk(isp)=k1
      dtastk(isp+1:isp+isp0-i)=dtastk(i:isp0-1)
      isp=isp+isp0-1
      kx=kxcompose(isp0)
      isp=isp0-1
      irtc=0
      return
      end

      recursive function tfrecompilearg(k,rep,irtc) result(kx)
      use tfcode
      implicit none
      type (sad_descriptor) kx,k,k1,k2,kd
      type (sad_dlist), pointer :: list,klx
      type (sad_rlist), pointer :: klr
      type (sad_pat), pointer :: pat
      type (sad_symbol), pointer :: sym2
      integer*8 ka1
      integer*4 irtc,i,m,isp1
      logical*4 rep,rep1,rep2
      irtc=0
      rep=.false.
      kx=k
      if(ktflistq(k,list))then
        if(iand(list%attr,lmemberlist) == 0)then
          return
        endif
        k1=list%head
        if(k1%k == ktfoper+mtfhold)then
          return
        endif
        k1=tfrecompilearg(k1,rep,irtc)
        if(ktfreallistq(list))then
          if(rep)then
            m=list%nl
            kx=kxavaloc(-1,m,klr)
            klr%rbody(1:m)=list%rbody(1:m)
c            call tmov(rlist(ka+1),rlist(kax+1),m)
            klr%attr=ior(larglist,list%attr)
            klr%head=dtfcopy(k1)
          endif
          return
        endif
        isp1=isp
        isp=isp+1
        dtastk(isp)=k1
        rep2=.false.
        do i=1,list%nl
          isp=isp+1
          dtastk(isp)=tfrecompilearg(list%dbody(i),rep1,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
          rep2=rep2 .or. rep1
        enddo
        if(list%head%k == ktfoper+mtfatt)then
          if(isp == isp1+3)then
            k2=dtastk(isp1+2)
            if(ktfsymbolq(k2,sym2))then
              if(sym2%override /= 0)then
                if(iand(sym2%attr,iattrdynamic) /= 0)then
                  go to 120
                endif
              endif
c              call tfdebugprint(ktastk(isp1+1),'rcmparg',3)
              call tfatt(isp1+1,kx,.false.,irtc)
              if(irtc > 0)then
                isp=isp1
                return
              elseif(irtc == 0)then
c                call tfdebugprint(kx,'==>',3)
                isp=isp1
                rep=.true.
                return
              endif
              irtc=0
            endif
          endif
        endif
 120    if(rep .or. rep2)then
          kx%k=ktflist+ktfcompose(isp1+1,klx)
          klx%attr=ior(larglist,list%attr)
          rep=.true.
        endif
        isp=isp1
      elseif(ktfpatq(k,pat))then
        k1=pat%expr
        if(ktfrefq(k1,ka1) .and. ka1 > 3)then
          k1=tfrecompilearg(k1,rep,irtc)
          if(irtc /= 0)then
            return
          endif
        endif
        kd=pat%default
        kd=tfrecompilearg(kd,rep1,irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        if(rep)then
          kx=kxpcopyss(k1,pat%head,pat%sym%alloc,kd)
        endif
      endif
      return
      end function

      subroutine tfclearmemberobject(isp1,kx,irtc)
      use objsym
      use tfcode
      use modul,only:tfdelete
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      integer*4 itfmessage
      if(kxmemberobject%k == 0)then
        call tfcxinit
      endif
      if(isp1+1 /= isp)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonsymbolq(ktastk(isp1+1),sym))then
        irtc=itfmessage(9,'General::wrongtype','"symbol for #1"')
        return
      endif
      irtc=0
      if(sym%override == 0)then
        sym=>tfsydef(sym)
      endif
      call sym_namtbl(sym,loc)
      call sym_symdef(sym,symd)
      call tfdelete(symd,.true.,.false.)
      irec=min(irec+1,nrecycle)
      irectbl(irec)=loc%str%nc
      kx%k=ktfoper+mtfnull
      return
      end

      subroutine tfmemberscan(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kc
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: listpc
      integer*8 kcv,kam
      integer*4 l,itfmessage,n,isp0,i
      logical*4 ev
      if(kxmemberobject%k == 0)then
        call tfcxinit
      endif
      if(isp /= isp1+4)then
        irtc=itfmessage(9,'General::narg','"4"')
        return
      endif
      call loc_sad(ktfaddr(ktastk(isp-1)),listpc)
      n=listpc%nl
      irtc=0
      if(n <= 0)then
        kx=kxundefined
        return
      endif
      kcv=ktastk(isp1+2)
      isp0=isp
      levele=levele+1
      if(ktastk(isp0) /= 0)then
        do i=1,n
          isp=isp0+1
          dtastk(isp)=kxmemberobject
          isp=isp+1
          dtastk(isp)=listpc%dbody(i)
          isp=isp+1
          ktastk(isp)=kcv
          call tftocontext(isp-2,kc,irtc)
          if(irtc /= 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          isp=isp0+2
          dtastk(isp)=kc
          kam=ktfcompose(isp-1)
          ktastk(isp-1)=ktflist+kam
          ktastk(isp)=ktastk(isp1+1)
          isp=isp+1
          dtastk(isp)=listpc%dbody(i)
          call tfdeval(isp-2,kxmemberobject,kx,1,.false.,ev,irtc)
          if(irtc /= 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          if(kx%k /= kxundefined%k)then
            isp=isp0
            call tfthrow(irtcret,kx,irtc)
            l=itfdownlevel()
            return
          endif
        enddo
      else
        isp=isp0+1
        dtastk(isp)=kxmemberobject
        isp=isp+1
        ktastk(isp)=kcv
        kam=ktfcompose(isp-1)
        do i=1,n
          isp=isp0+1
          ktastk(isp)=ktflist+kam
          isp=isp+1
          ktastk(isp)=ktastk(isp1+1)
          isp=isp+1
          dtastk(isp)=listpc%dbody(i)
          call tfdeval(isp-2,kxmemberobject,kx,1,.false.,ev,irtc)
          if(irtc /= 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          if(kx%k /= kxundefined%k)then
            isp=isp0
            call tfthrow(irtcret,kx,irtc)
            l=itfdownlevel()
            return
          endif
        enddo
      endif
      isp=isp0
      l=itfdownlevel()
      kx=kxundefined
      return
      end

      logical*4 function tfrepsymstk(sym,isp1,isp2,nrule1,nrule2,kx)
      use repl,only: tfmatchsymstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1,isp2,nrule1,nrule2
      type (sad_symbol) sym
      type (sad_dlist), pointer :: klx
      type (sad_namtbl), pointer :: loc
      integer*8 kaj,ka
      integer*4 j
      call loc_namtbl(sym%loc,loc)
      ka=sad_loc(sym%loc)
      if(tfmatchsymstk(ka,isp1,nrule1,j) .or.
     $     loc%cont /= itfcontroot .and. loc%cont /= isyscont
     $     .and. tfmatchsymstk(ka,isp2,nrule2,j))then
        kaj=ktfaddr(ktastk(ivstk2(1,j)))
        kx=kxadaloc(-1,2,klx)
        klx%head%k=ktfoper+mtfatt
        klx%dbody(1)%k=ktflist+ktfcopy1(icx)
        klx%dbody(2)%k=ktfsymbol+ktfcopy1(kaj)
        tfrepsymstk=.true.
      else
        tfrepsymstk=.false.
      endif
      return
      end

      recursive function tfreplacememberstk(isp1,isp2,isp3,
     $     nrule1,nrule2,k,m0,rep) result(kx)
      use tfcode
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,isp2,isp3,nrule1,nrule2,m0
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) k1,ki,kd
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: kl,klx,kli
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 i,m,isp0,n,m01
      logical*4 ,intent(out):: rep
      logical*4 rep1,rule
      if(icx == 0)then
        call tfcxinit
      endif
      rep=.false.
      kx=k
      if(ktflistq(k,kl))then
        k1=kl%head
        n=kl%nl
        rule=.false.
        if(k1%k == ktfoper+mtfatt)then
          m=1
          isp=isp+1
          ktastk(isp)=ktfoper+mtfatt
        else
          m=n
          if(tfclassq(k1))then
            isp=isp+1
            dtastk(isp)=k1
            rule=.true.
          else
            if(k1%k == kxliteral .and. n == 1)then
              kx=kl%dbody(1)
              rep=.true.
              return
            endif
            isp=isp+1
            if(.not. ktfstringq(k1))then
              dtastk(isp)=tfreplacememberstk(isp1,isp2,isp3,
     $             nrule1,nrule2,k1,0,rep)
            else
              dtastk(isp)=kl%head
            endif
          endif
        endif
        isp0=isp
        if(ktfreallistq(kl))then
          if(rep)then
            call loc_sad(ktavaloc(-1,n),klx)
            klx%head=dtfcopy(dtastk(isp0))
            klx%dbody(1:n)=kl%dbody(1:n)
            kx=sad_descr(klx)
          endif
          isp=isp0-1
          return
        endif
        LOOP_I: do i=1,n
          isp=isp+1
          dtastk(isp)=kl%dbody(i)
          if(i <= m .and. i /= m0)then
            ki=kl%dbody(i)
            m01=0
            if(rule .and. ktflistq(ki,kli))then
              if(kli%head%k == ktfoper+mtfrule .or.
     $             kli%head%k == ktfoper+mtfruledelayed)then
                m01=1
              endif
            endif
            if(.not. ktfstringq(ki) .and. ktfnonrealq(ki))then
              dtastk(isp)=tfreplacememberstk(isp1,isp2,isp3,
     $             nrule1,nrule2,ki,m01,rep1)
              rep=rep .or. rep1
            endif
          endif
        enddo LOOP_I
        if(rep)then
          kx=kxcompose(isp0)
        endif
        isp=isp0-1
        return
      elseif(ktfpatq(k,pat))then
        rep=.false.
        k1=pat%expr
        if(ktfobjq(k1))then
          k1=tfreplacememberstk(isp1,isp2,isp3,
     $         nrule1,nrule2,k1,0,rep)
        endif
        kd=pat%default
        if(ktfobjq(kd))then
          kd=tfreplacememberstk(isp1,isp2,isp3,
     $         nrule1,nrule2,kd,0,rep1)
          rep=rep .or. rep1
        endif
        if(rep)then
          kx=kxpcopyss(k1,pat%head,pat%sym%alloc,kd)
        endif
      else
        if(ktfoperq(k,ka))then
          call loc_sad(klist(ifunbase+ka),sym)
        elseif(ktfnonsymbolq(k,sym))then
          return
        endif
        if(sym%loc == ithisloc .and. sym%gen <= 0)then
          kx%k=ktflist+icx
          rep=.true.
        elseif(tfrepsymstk(sym,isp1,isp2,nrule1,nrule2,kx))then
          rep=.true.
        endif
      endif
      return
      end

      subroutine tfreplacemember(isp1,kx,irtc)
      use repl,only:tfsortsymbolstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ka,ki
      type (sad_dlist), pointer :: klm
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*4 itfmessage,ispa,ispb,i,ispc,nrule1,nrule2
      logical*4 rep
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"List for #1"')
        return
      endif
      if(ktfnonlistq(ktastk(isp),klm))then
        irtc=itfmessage(9,'General::wrongtype','"List for #2"')
        return
      endif
      ispa=isp
      call tfreplacedefarg(dtastk(isp1+1),ka,irtc)
      if(irtc /= 0)then
        return
      endif
      dtastk(isp+1:isp+2*klm%nl-1:2)=klm%dbody(1:klm%nl)
      isp=isp+2*klm%nl
c      do i=1,klm%nl
c        isp=isp+2
c        dtastk(isp-1)=klm%dbody(i)
c      enddo
      ispb=isp
      do i=ispa+1,ispb,2
        isp=isp+2
        ki=dtastk(i)
        if(.not. ktfsymbolq(ki,sym))then
          go to 9000
        endif
        call sym_symstr(sym,str)
        ktastk(isp-1)=ktfaddrd(ki)
        dtastk(i)=kxsymbolf(str%str,str%nch,.false.)
      enddo
      ispc=isp
c      write(*,*)'repmember ',ispa,ispc,ispb,klm%nl
      call tfsortsymbolstk(ispa,(ispb-ispa)/2,nrule1)
      call tfsortsymbolstk(ispb,(ispc-ispb)/2,nrule2)
      kx=tfreplacememberstk(ispa,ispb,ispc,nrule1,nrule2,ka,0,rep)
c      call tfdebugprint(ka,'repmember-9',3)
c      call tfdebugprint(kx,'==>',3)
      irtc=0
      isp=ispa
      return
 9000 irtc=itfmessage(9,'General::rongtype',
     $     '"List of Symbols for #2"')
      isp=ispa
      return
      end

      recursive logical*4 function tfclassq(k) result(lx)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) k1
      type (sad_dlist), pointer ::list
      type (sad_symdef), pointer ::symd
      lx=.false.
      if(ktflistq(k,list))then
        k1=list%head
        if(ktfsymbolqdef(k1%k,symd) .and. list%nl == 1)then
          if(symd%sym%gen == -3)then
            if(ktfnonreallistqo(list))then
              if(list%dbody(1)%k == k1%k)then
                lx=.true.
              endif
            endif
          endif
        endif
      elseif(ktfsymbolqdef(k%k,symd))then
        if(symd%sym%override /= 0 .and.
     $       ktflistq(symd%value))then
          lx=tfclassq(symd%value)
        endif
      endif
      return
      end

      subroutine tfreplacedefarg(k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kx1
      type (sad_dlist), pointer :: kl,klx
      integer*4 itfmessage
      integer*4 ,intent(out):: irtc
      if(ktfnonlistq(k,kl) .or. kl%head%k /= ktfoper+mtfhold
     $     .or. kl%nl /= 1)then
        go to 9000
      endif
      call tfreplacedefarg1(kl%dbody(1),kx1,irtc)
      if(irtc /= 0)then
        return
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head%k=ktfoper+mtfhold
      klx%dbody(1)=dtfcopy(kx1)
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Hold[definition]"')
      return
      end

      recursive subroutine tfreplacedefarg1(k1,kx,irtc)
      use dset,only:tfreplacedef
      implicit none
      type (sad_descriptor) ,intent(in):: k1
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kx1,kargr,kr
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kl1,klx1
      integer*8 ka
      integer*4 itfmessage,id
      if(ktfnonlistq(k1,kl1) .or. ktfnonoperq(kl1%head,ka))then
        go to 9100
      endif
      select case (int(ka))
      case (mtfset,mtfsetdelayed,mtfupset,mtfupsetdelayed)
        if(kl1%nl /= 2)then
          go to 9100
        endif
        call tfreplacedef(kl1%dbody(2),kl1%dbody(1),
     $       kr,kargr,irtc)
        if(irtc /= 0)then
          return
        endif
        if(kl1%dbody(2)%k == kr%k .and.
     $       kl1%dbody(1)%k == kargr%k)then
          kx=k1
          return
        endif
        kx=kxadaloc(-1,2,klx1)
        klx1%head%k=ktfoper+ka
        klx1%dbody(1)=dtfcopy(kargr)
        klx1%dbody(2)=dtfcopy(kr)
      case default
        id=iget_fun_id(ka)
        if(id == nfunif .or. id == nfunwith)then
          if(kl1%nl /= 2)then
            go to 9100
          endif
          call tfreplacedefarg1(kl1%dbody(2),kx1,irtc)
          if(irtc /= 0)then
            return
          endif
          kx=kxadaloc(-1,2,klx1)
          klx1%head%k=ktfoper+ka
          klx1%dbody(1)=dtfcopy(kl1%dbody(1))
          klx1%dbody(2)=dtfcopy(kx1)
        else
          go to 9100
        endif
      end select
      return
 9100 irtc=itfmessage(9,'General::wrongtype',
     $     '"Hold[f = g, f := g. f ^= g, f ^:= g, With[ .. ],'//
     $     ' or If[ .. ] ]"')
      return
      end

      end module tfcx
