      module tfcx
      use tfstk
      integer*8, save :: icx=0,ithis,ithisloc,isyscont
      type (sad_descriptor) kxundefined,kxmemberobject,kxmemberobject2
      data kxundefined%k,kxmemberobject%k,kxmemberobject2%k /0,0,0/

      contains
        subroutine tfcxinit
        use tfstk
        implicit none
        type (sad_list), pointer :: kl
        icx=ktadaloc(0,1,kl)
        kl%dbody(0)=dtfcopy1(kxsymbolz('Class`cx$',9))
        kl%dbody(1)=dtfcopy1(kxsymbolz('Class`x$',8))
        ithis=ktfsymbolz('System`This',11)
        ithisloc=klist(ithis)
        isyscont=ktfaddr(klist(ktfsymbolz('`System`',8)-4))
        kxmemberobject=kxsymbolz('Class`MemberObject',18)
        kxmemberobject2=kxsymbolz('Class`MemberObject2',19)
        kxundefined   =kxsymbolz('System`Undefined',16)
        return
        end subroutine

      end module

      recursive subroutine tfsolvemember(list,kx,reps,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      type (sad_list) list
      type (sad_list), pointer :: listx
      integer*4 irtc
      logical*4 eval,reps
      irtc=-1
      kx%k=ktflist+ksad_loc(list%head)
      reps=.false.
      if(list%head .eq. ktfoper+mtfatt)then
        if(list%nl .eq. 2 .and. ktfnonreallistqo(list))then
          call tfclassmember(list%dbody(1),list%dbody(2),kx,
     $         .false.,irtc)
        endif
        return
      elseif(list%head .eq. ktfoper+mtfslot)then
        call tfslot(int8(mtfslot),list,kx,.false.,irtc)
        if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
          call tfreseterror
          irtc=-1
        endif
        reps=.true.
        return
      endif
      if(ktflistq(list%head,listx))then
        call tfsolvemember(listx,k1,reps,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(reps)then
          call tfduplist(list,listx)
        else
          call tfclonelist(list,listx)
        endif
        call tfreplist(listx,0,k1,eval)
c        call tfloadlstk(list,listx)
c        listx%head=k1
c        call tfstk2l(listx,listx)
        kx%k=ktflist+ksad_loc(listx%head)
      endif
      return
      end

      subroutine tfatt(isp1,kx,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      integer*8 isp1,isp0,i
      integer*4 irtc
      logical*4 eval
      if(isp .eq. isp1+2)then
        call tfclassmember(dtastk(isp1+1),dtastk(isp),kx,eval,irtc)
        return
      endif
      k1=dtastk(isp1+1)
      do i=isp1+2,isp
        call tfclassmember(k1,dtastk(i),k1,eval,irtc)
c        call tfdebugprint(k1,'tfatt',1)
c        write(*,*)'with: ',irtc,i
        if(irtc .ne. 0)then
          go to 10
        endif
      enddo
      kx=k1
      return
 10   if(irtc .gt. 0)then
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

      recursive subroutine tfclassmember(k10,k20,kx,eval,irtc)
      use tfstk
      use tfcx
      implicit none
      type (sad_descriptor) k10,k20,kx,k1,k2,ks
      type (sad_symdef), pointer :: symd
      type (sad_symbol), pointer :: symh
      type (sad_list), pointer :: kl1,klx
      integer*8 ka2,k11,kv,isp1
      integer*4 irtc,l,itfdownlevel
      logical*4 ev,eval
      k1=k10
      k2=k20
      if(ktfoperqd(k2,ka2))then
        k2%k=ktfsymbol+klist(ifunbase+ka2)
      endif
      if(ktfsymbolqd(k2))then
        if(ktflistqd(k1,kl1))then
          k11=kl1%head
          if(ktfsymbolqdef(k11,symd))then
            if(symd%sym%gen .eq. -3)then
              go to 10
            elseif(symd%sym%override .ne. 0 .and.
     $             symd%downval .ne. 0)then
              call tfgetdefargp(kl1,ktfaddr(k11),kv,ev,irtc)
              if(irtc .gt. 0)then
                return
              endif
              if(ev .and. irtc .eq. 0)then
                call tfclassmember(dlist(kv),k2,kx,eval,irtc)
                return
              endif
            endif
          endif
        elseif(ktfsymbolqdef(k1%k,symd) .and.
     $         symd%sym%override .ne. 0)then
          if(ktflistqd(symd%value,kl1) .and.
     $         ktfsymbolq(kl1%head,symh))then
            if(symh%gen .eq. -3)then
              go to 10
            endif
          endif
        endif
      endif
      irtc=-1
      return
 10   if(kxmemberobject2%k .eq. 0)then
        call tfcxinit
      endif
      if(kl1%nl .ne. 1)then
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
      if(irtc .gt. 0)then
        isp=isp1-1
        l=itfdownlevel()
        return
      elseif(ev .and. irtc .eq. 0)then
        dtastk(isp1)=ks
      else
        dtastk(isp1)=kxcompose(isp1)
        ks=kxmemberobject2
      endif
      isp=isp1+1
      dtastk(isp)=kl1%dbody(0)
      isp=isp+1
      dtastk(isp)=k2
      call tfdeval(isp1,ks,kx,1,.false.,ev,irtc)
c      call tfdebugprint(kx,'classmember',1)
c      write(*,*)'with: ',irtc,ev,eval
      isp=isp1-1
      if(ev)then
        if(irtc .ne. 0)then
          l=itfdownlevel()
          return
        endif
        if(ktflistqd(kx,klx))then
          if(klx%head .eq. ktfoper+mtfhold)then
            kx=klx%dbody(1)
          endif
        elseif(kx%k .eq. kxundefined%k)then
          irtc=-1
          l=itfdownlevel()
          return
        endif
        if(eval)then
          if(ktfsymbolqd(kx))then
            call tfsyeval(kx,kx,irtc)
            if(irtc .ne. 0)then
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

      subroutine tfreplacemember(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ka,ki
      type (sad_list), pointer :: klm
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*8 isp1,ispa,ispb,i,ispc
      integer*4 irtc,itfmessage,nrule1,nrule2
      logical*4 rep
      if(isp .ne. isp1+2)then
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
      if(irtc .ne. 0)then
        return
      endif
c      call tfdebugprint(ka,'repmember',3)
      do i=1,klm%nl
        isp=isp+2
        ktastk(isp-1)=klm%body(i)
      enddo
      ispb=isp
      do i=ispa+1,ispb,2
        isp=isp+2
        ki=dtastk(i)
        if(.not. ktfsymbolqd(ki,sym))then
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
      call tfreplacememberstk(ispa,ispb,ispc,nrule1,nrule2,ka,kx,0,rep)
c      call tfdebugprint(ka,'repmember-9',3)
c      call tfdebugprint(kx,'==>',3)
      irtc=0
      isp=ispa
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"List of Symbols for #2"')
      isp=ispa
      return
      end

      recursive subroutine tfreplacememberstk(isp1,isp2,isp3,
     $     nrule1,nrule2,k,kx,m0,rep)
      use tfstk
      use tfcode
      use tfcx
      implicit none
      type (sad_descriptor) k,kx,k1,ki,kd
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: kl,klx,kli
      type (sad_symbol), pointer :: sym
      integer*8 ka,isp1,isp2,isp3,isp0
      integer*4 i,m,n,m0,m01,nrule1,nrule2
      logical*4 rep,rep1,tfclassq,rule,tfrepsymstk
      if(icx .eq. 0)then
        call tfcxinit
      endif
      rep=.false.
      kx=k
      if(ktflistqd(k,kl))then
        k1=kl%dbody(0)
        n=kl%nl
        rule=.false.
        if(k1%k .eq. ktfoper+mtfatt)then
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
            if(k1%k .eq. kxliteral .and. n .eq. 1)then
              kx=kl%dbody(1)
              rep=.true.
              return
            endif
            isp=isp+1
            if(.not. ktfstringqd(k1))then
              call tfreplacememberstk(isp1,isp2,isp3,
     $             nrule1,nrule2,k1,dtastk(isp),0,rep)
            else
              ktastk(isp)=kl%head
            endif
          endif
        endif
        isp0=isp
        if(ktfreallistqo(kl))then
          if(rep)then
            call loc_list(ktavaloc(-1,n),klx)
            klx%head=ktfcopy(ktastk(isp0))
            klx%body(1:n)=kl%body(1:n)
            kx%k=ktflist+ksad_loc(klx%head)
          endif
          isp=isp0-1
          return
        endif
        LOOP_I: do i=1,n
          isp=isp+1
          dtastk(isp)=kl%dbody(i)
          if(i .le. m .and. i .ne. m0)then
            ki=kl%dbody(i)
            m01=0
            if(rule .and. ktflistqd(ki,kli))then
              if(kli%head .eq. ktfoper+mtfrule .or.
     $             kli%head .eq. ktfoper+mtfruledelayed)then
                m01=1
              endif
            endif
            if(.not. ktfstringqd(ki) .and. ktfnonrealqd(ki))then
              call tfreplacememberstk(isp1,isp2,isp3,nrule1,nrule2,
     $             ki,dtastk(isp),m01,rep1)
              rep=rep .or. rep1
            endif
          endif
        enddo LOOP_I
        if(rep)then
          kx=kxcompose(isp0)
        endif
        isp=isp0-1
        return
      elseif(ktfpatqd(k,pat))then
        rep=.false.
        k1=pat%expr
        if(ktfobjqd(k1))then
          call tfreplacememberstk(isp1,isp2,isp3,
     $         nrule1,nrule2,k1,k1,0,rep)
        endif
        kd=pat%default
        if(ktfobjqd(kd))then
          call tfreplacememberstk(isp1,isp2,isp3,
     $         nrule1,nrule2,kd,kd,0,rep1)
          rep=rep .or. rep1
        endif
        if(rep)then
          kx=kxpcopyss(k1,pat%head,pat%sym%alloc,kd)
        endif
      else
        if(ktfoperqd(k,ka))then
          call loc_sad(klist(ifunbase+ka),sym)
        elseif(ktfnonsymbolqd(k,sym))then
          return
        endif
        if(sym%loc .eq. ithisloc .and. sym%gen .le. 0)then
          kx%k=ktflist+icx
          rep=.true.
        elseif(tfrepsymstk(sym,isp1,isp2,nrule1,nrule2,kx))then
          rep=.true.
        endif
      endif
      return
      end

      logical*4 function tfrepsymstk(sym,isp1,isp2,nrule1,nrule2,kx)
      use tfstk
      use tfcx
      implicit none
      type (sad_descriptor) kx
      type (sad_symbol) sym
      type (sad_list), pointer :: klx
      type (sad_namtbl), pointer :: loc
      integer*8 kaj,ka,isp1,isp2
      integer*4 j,nrule1,nrule2
      logical*4 tfmatchsymstk
      call loc_namtbl(sym%loc,loc)
      ka=sad_loc(sym%loc)
      if(tfmatchsymstk(ka,isp1,nrule1,j) .or.
     $     loc%cont .ne. itfcontroot .and. loc%cont .ne. isyscont
     $     .and. tfmatchsymstk(ka,isp2,nrule2,j))then
        kaj=ktfaddr(ktastk(ivstk2(1,j)))
        kx=kxadaloc(-1,2,klx)
        klx%head=ktfoper+mtfatt
        klx%body(1)=ktflist+ktfcopy1(icx)
        klx%body(2)=ktfsymbol+ktfcopy1(kaj)
        tfrepsymstk=.true.
      else
        tfrepsymstk=.false.
      endif
      return
      end

      subroutine tfreplacedefarg(k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx,kx1
      type (sad_list), pointer :: kl,klx
      integer*4 itfmessage,irtc
      if(ktfnonlistqd(k,kl) .or. kl%head .ne. ktfoper+mtfhold
     $     .or. kl%nl .ne. 1)then
        go to 9000
      endif
      call tfreplacedefarg1(kl%dbody(1),kx1,irtc)
      if(irtc .ne. 0)then
        return
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head=ktfoper+mtfhold
      klx%dbody(1)=dtfcopy(kx1)
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Hold[definition]"')
      return
      end

      recursive subroutine tfreplacedefarg1(k1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,kx,kx1,kargr,kr
      type (sad_list), pointer :: kl1,klx1
      integer*8 ka
      integer*4 itfmessage,irtc,id
      if(ktfnonlistqd(k1,kl1) .or. ktfnonoperqd(kl1%dbody(0),ka))then
        go to 9100
      endif
      select case (int(ka))
      case (mtfset,mtfsetdelayed,mtfupset,mtfupsetdelayed)
        if(kl1%nl .ne. 2)then
          go to 9100
        endif
        call tfreplacedef(kl1%dbody(2),kl1%dbody(1),
     $       kr,kargr,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(kl1%dbody(2)%k .eq. kr%k .and.
     $       kl1%dbody(1)%k .eq. kargr%k)then
          kx=k1
          return
        endif
        kx=kxadaloc(-1,2,klx1)
        klx1%head=ktfoper+ka
        klx1%dbody(1)=dtfcopy(kargr)
        klx1%dbody(2)=dtfcopy(kr)
      case default
        id=iget_fun_id(ka)
        if(id .eq. nfunif .or. id .eq. nfunwith)then
          if(kl1%nl .ne. 2)then
            go to 9100
          endif
          call tfreplacedefarg1(kl1%dbody(2),kx1,irtc)
          if(irtc .ne. 0)then
            return
          endif
          kx=kxadaloc(-1,2,klx1)
          klx1%head=ktfoper+ka
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

      recursive logical*4 function tfclassq(k) result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) k,k1
      type (sad_list), pointer ::list
      type (sad_symdef), pointer ::symd
      lx=.false.
      if(ktflistqd(k,list))then
        k1=list%dbody(0)
        if(ktfsymbolqdef(k1%k,symd) .and. list%nl .eq. 1)then
          if(symd%sym%gen .eq. -3)then
            if(ktfnonreallistqo(list))then
              if(list%body(1) .eq. k1%k)then
                lx=.true.
              endif
            endif
          endif
        endif
      elseif(ktfsymbolqdef(k%k,symd))then
        if(symd%sym%override .ne. 0 .and.
     $       ktflistqd(symd%value))then
          lx=tfclassq(symd%value)
        endif
      endif
      return
      end

      subroutine tfmemberscan(isp1,kx,irtc)
      use tfstk
      use tfcx
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: listpc
      integer*8 kcv,kam,kc,isp1,isp0
      integer*4 irtc,l,itfmessage,n,i,itfdownlevel
      logical*4 ev
      if(kxmemberobject%k .eq. 0)then
        call tfcxinit
      endif
      if(isp .ne. isp1+4)then
        irtc=itfmessage(9,'General::narg','"4"')
        return
      endif
      call loc_list(ktfaddr(ktastk(isp-1)),listpc)
      n=listpc%nl
      irtc=0
      if(n .le. 0)then
        kx=kxundefined
        return
      endif
      kcv=ktastk(isp1+2)
      isp0=isp
      levele=levele+1
      if(ktastk(isp0) .ne. 0)then
        do i=1,n
          isp=isp0+1
          dtastk(isp)=kxmemberobject
          isp=isp+1
          ktastk(isp)=listpc%body(i)
          isp=isp+1
          ktastk(isp)=kcv
          call tftocontext(isp-2,kc,irtc)
          if(irtc .ne. 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          isp=isp0+2
          ktastk(isp)=kc
          kam=ktfcompose(isp-1)
          ktastk(isp-1)=ktflist+kam
          ktastk(isp)=ktastk(isp1+1)
          isp=isp+1
          ktastk(isp)=listpc%body(i)
          call tfdeval(isp-2,kxmemberobject,kx,1,.false.,ev,irtc)
          if(irtc .ne. 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          if(kx%k .ne. kxundefined%k)then
            isp=isp0
            call tfthrow(0,kx,irtc)
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
          ktastk(isp)=listpc%body(i)
          call tfdeval(isp-2,kxmemberobject,kx,1,.false.,ev,irtc)
          if(irtc .ne. 0)then
            l=itfdownlevel()
            isp=isp0
            return
          endif
          if(kx%k .ne. kxundefined%k)then
            isp=isp0
            call tfthrow(0,kx,irtc)
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

      module objsym
      implicit none
      integer*4 nrecycle
      parameter (nrecycle=2048)
      integer*4 :: irec=0,irectbl(nrecycle)
      end module

      subroutine tfobjectsymbol(isp1,kx,irtc)
      use objsym
      use tfstk
      use tfcode
      implicit none
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      type (sad_descriptor) kx
      integer*8 isp1
      integer*4 irtc,itfmessage,ls,l,n,n1,ifrac,id,n0
      character*1024 name
      character*32 buf
      character ch
      data id/0/
      if(isp1+1 .ne. isp)then
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
      if(irec .ne. 0)then
        n0=irectbl(irec)
        irec=irec-1
      else
        id=id+1
        n0=id
      endif
      n=n0
      l=32
      do while(n .ne. 0)
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
      loc%str%nc=n0
      irtc=0
      return
      end

      subroutine tfclearmemberobject(isp1,kx,irtc)
      use objsym
      use tfstk
      use tfcode
      use tfcx
      implicit none
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage
      if(kxmemberobject%k .eq. 0)then
        call tfcxinit
      endif
      if(isp1+1 .ne. isp)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonsymbolq(ktastk(isp1+1),sym))then
        irtc=itfmessage(9,'General::wrongtype','"symbol for #1"')
        return
      endif
      irtc=0
      if(sym%override .eq. 0)then
        call tfsydef(sym,sym)
      endif
      call sym_namtbl(sym,loc)
      call sym_symdef(sym,symd)
      call tfdelete(symd,.true.,.false.)
      irec=min(irec+1,nrecycle)
      irectbl(irec)=loc%str%nc
      kx%k=ktfoper+mtfnull
      return
      end
