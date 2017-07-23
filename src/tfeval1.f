      recursive subroutine tfeval1(k1,k2,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      integer*4 irtc,iopc1,isp1
      select case (iopc1)
      case (mtfrevpower)
        call tfeval1(k2,k1,kx,mtfpower,irtc)
      case (mtfcomp)
        irtc=0
        kx=k2
      case (mtffun:mtfruledelayed,mtfpattest,mtfslot,mtfslotseq,
     $       mtfalt,mtfrepeated,mtfrepeatednull)
        call tfeexpr(k1,k2,kx,iopc1)
        irtc=0
      case (mtfnull:mtfdiv,mtfpower,mtfgreater:mtfleq,mtfnot:mtfor,
     $       mtfleftbra:mtfrightbrace,mtfcomplex:mtfcomma,
     $       mtfdot,mtfend)
        if(tfnumberqd(k2) .and.
     $     (tfnumberqd(k1) .or. iopc1 .eq. mtfnot))then
          call tfcmplx(k1,k2,kx,iopc1,irtc)
          return
        endif
        if(tflistq(k2) .or. tflistq(k1))then
          call tfearray(k1,k2,kx,iopc1,irtc)
        else
          call tfeexpr(k1,k2,kx,iopc1)
          irtc=0
        endif
      case default
        isp=isp+1
        isp1=isp
        ktastk(isp)=ktfoper+iopc1
        isp=isp+1
        dtastk(isp)=k1
        isp=isp+1
        dtastk(isp)=k2
        call tfefunref(isp1,kx,.true.,irtc)
        isp=isp1-1
      end select
      return
      end

      subroutine tfflagordef(isp1,kx,irtc)
      use tfstk
      type (sad_descriptor) kx
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,nc
      real*8 vx,fflogi
      logical*4 exist
      character*8 name
      if(isp .ne. isp1+1)then
        call tfdefinition(isp1,kx,irtc)
        return
      endif
      if(ktfsymbolqd(dtastk(isp),sym))then
        call sym_symstr(sym,str)
        nc=min(8,str%nch)
        name(1:nc)=str%str(1:nc)
        call capita(name(1:nc))
        vx=fflogi(name(1:nc),exist)
        if(exist)then
          kx=dfromr(vx)
          irtc=0
        else
          call tfdefinition(isp1,kx,irtc)
        endif
      else
        call tfdefinition(isp1,kx,irtc)
      endif
      return
      end

      subroutine tfeval1to(k1,k2,kx,iopc,old,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,kv,kr,ku,ks
      type (sad_list), pointer :: kl1
      integer*4 iopc,irtc
      logical*4 old
      if(ktflistq(k1,kl1))then
        call tfleval(kl1,kv,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      elseif(ktfsymbolqd(k1))then
        call tfsyeval(k1,kv,irtc)
        if(irtc .ne. 0)then
          return
        endif
      elseif(ktfpatqd(k1))then
        call tfpateval(k1,kv,irtc)
        if(irtc .ne. 0)then
          return
        endif
      else
        kv=k1
      endif
      if(iopc .eq. mtfaddto)then
        call tfeval1(kv,k2,kr,mtfplus,irtc)
      elseif(iopc .eq. mtftimesby)then
        call tfeval1(kv,k2,kr,mtfmult,irtc)
      elseif(iopc .eq. mtfsubtractfrom)then
        call tfeval1(-1.d0,k2,ku,mtfmult,irtc)
        if(irtc .ne. 0)then
          return
        endif
        call tfeval1(kv,ku,kr,mtfplus,irtc)
      else
        call tfeval1(k2,-1.d0,ku,mtfpower,irtc)
        if(irtc .ne. 0)then
          return
        endif
        call tfeval1(kv,ku,kr,mtfmult,irtc)
      endif
      if(irtc .ne. 0)then
        return
      endif
      call tfeevaldef(k1,ks,irtc)
      if(irtc .ne. 0)then
        return
      endif
      call tfset1(ks,kr,kx,mtfset,irtc)
      if(old)then
        kx=kv
      endif
      return
      end

      subroutine tfappendto(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k10,kr,k1
      type (sad_list), pointer :: kl
      integer*2, parameter :: nextra = 8
      integer*8 kp
      integer*4 isp1,irtc,itfmessage,isp0,mode
      logical*4 def,tfgetstoredp,ev
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      k10=dtastk(isp1+1)
      if(tfgetstoredp(k10,kp,def,irtc))then
        call tfappendto1(kp,dtastk(isp),kr,mode,ev,irtc)
        if(irtc .ne. 0)then
          return
        elseif(.not. ev)then
          kx=dlist(kp)
          return
        endif
      else
        if(irtc .gt. 0)then
          return
        endif
        if(ktflistq(k10,kl))then
          if(kl%head%k .eq. ktfoper+mtfpart)then
            call tfleval(kl,kx,.false.,irtc)
            if(irtc .ne. 0)then
              return
            endif
            call loc_list(ktfaddrd(kx),kl)
            call tfsetpart(kl,ktastk(isp),kx,mode,irtc)
            return
          endif
        endif
        call tfeevalref(k10,k1,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfnonlistqd(k1))then
          irtc=itfmessage(9,'General::wrongtype','"List"')
          return
        endif
        if(mode .eq. 1)then
          call tfappend(k1,dtastk(isp),kr,.true.,0,irtc)
        else
          call tfappend(k1,dtastk(isp),kr,.true.,1,irtc)
        endif
        if(irtc .ne. 0)then
          return
        endif
      endif
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfoper+mtfset
      isp=isp+1
      ktastk(isp)=ktastk(isp1+1)
      isp=isp+1
      dtastk(isp)=kr
      call tfefunref(isp0+1,kx,.true.,irtc)
      isp=isp0
      return
      end

      subroutine tfappendto1(kp,k2,kr,mode,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kr
      type (sad_list), pointer :: kl,klr
      integer*2, parameter :: nextra = 8
      integer*8 kp
      integer*4 irtc,itfmessage, i,n,mode
      logical*4 eval,ov
      eval=.true.
      k1=dlist(kp)
      if(ktfnonlistqd(k1,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      n=kl%nl
      ov=ktfovrwrtq(kl) .and. kl%ref .eq. 1
      if(ov .and. mode .eq. 1 .and. kl%lena .gt. 0)then
        klr=>kl
        klr%nl=n+1
        klr%lena=klr%lena-1
        klr%body(n+1)=0
        call tfreplist(klr,n+1,k2,eval)
        if(.not. eval)then
          return
        endif
      elseif(ov .and. mode .eq. 2 .and. kl%lenp .gt. 0)then
        call loc_list(ktfaddr(k1)-1,klr)
        klr%body(-3:0)=kl%body(-3:0)
        klr%lenp=klr%lenp-1
        klr%nl=n+1
        klr%body(1)=0
        call tfreplist(klr,1,k2,eval)
        klist(kp)=ktflist+ktfaddr(k1)-1
        if(.not. eval)then
          return
        endif
      elseif(mode .eq. 1)then
        kr%k=ktaalocsp(n+1,kl%lenp,nextra,klr)
        if(ktfreallistq(kl))then
          klr%head=dtfcopy(kl%head)
          klr%body(1:n)=kl%body(1:n)
        else
          klr%attr=lnonreallist
          do i=0,n
            klr%body(i)=ktfcopy(kl%body(i))
          enddo
        endif
        klr%body(n+1)=0
        call tfreplist(klr,n+1,k2,eval)
      else
        kr%k=ktaalocsp(n+1,nextra,kl%lena,klr)
        klr%head=dtfcopy(kl%head)
        if(ktfreallistq(kl))then
          klr%body(2:n+1)=kl%body(1:n)
        else
          klr%attr=lnonreallist
          do i=1,n
            klr%body(i+1)=ktfcopy(kl%body(i))
          enddo
        endif
        klr%body(1)=0
        call tfreplist(klr,1,k2,eval)
      endif
      if(eval)then
        call tfleval(klr,kr,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      else
        kr%k=ktflist+ksad_loc(klr%head%k)
      endif
      eval=.true.
      return
      end

      logical*4 function tfgetstoredp(ks0,kp,def,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ks0,ks
      type (sad_list), pointer :: lists
      type (sad_symdef), pointer :: symd
      integer*8 kp
      integer*4 irtc,itfmessageexp,itfmessage
      logical*4 def,ev
      def=.false.
      irtc=0
      tfgetstoredp=.false.
      ks=ks0
      do while(ktflistq(ks,lists))
        def=.true.
        ks=lists%head
      enddo
      if(ktfsymbolqdef(ks%k,symd))then
        if(symd%sym%override .eq. 0)then
          irtc=itfmessageexp(999,'General::invset',ks0)
          return
        endif
        if(ktfprotectedqo(symd%sym))then
          irtc=itfmessage(999,'General::protect','""')
          return
        endif
        if(def)then
          call loc_list(ktfaddrd(ks0),lists)
          call tfgetdefargp(lists,ktfaddr(ks),kp,ev,irtc)
          if(irtc .ne. 0)then
            return
          endif
        else
          kp=ksad_loc(symd%value%k)
        endif
      else
        return
      endif
      tfgetstoredp=.true.
      return
      end function 

      subroutine tfset1(k10,k20,kx,mopc,irtc)
      use tfstk
      use mackw
      implicit none
      type (sad_descriptor) k1,k2,kx,k10,k20,ks,ka
      type (sad_list),pointer :: list,kls1,kla,kls
      type (sad_symbol), pointer ::sym
      type (sad_symdef),pointer :: symd
      integer*8 ka1,kaa,kas,kar
      integer*4 irtc,mopc,itfmessage,itfmessageexp,i
      irtc=0
      k1=k10
      k2=k20
      kx=k2
      if(ktflistq(k1,list))then
        ka=list%head
        if(ktfoperqd(ka,kaa))then
          select case (kaa)
          case (mtflist)
            call tfearray(k1,k2,kx,mtfset,irtc)
            return
          case (mtfpart)
            call tfsetpart(list,k2,kx,0,irtc)
            return
          case(mtfmessagename)
            call loc_symdef(klist(ifunbase+mtfmessagename),symd)
            go to 12
          case (mtftagset)
            call tftagset(list,k2,kx,mopc,irtc)
            return
          end select
          go to 9900
        elseif(ktfsymbolqdef(ka%k,symd))then
          if(symd%sym%override .eq. -2)then
            if(ktflistq(symd%value,kls1) .and.
     $           kls1%head%k .eq. ktfoper+mtflist)then
              kas=ktadaloc(-1,list%nl+1,kls)
              kls%head%k=ktfoper+mtfpart
              kls%dbody(1)=dtfcopy1(ka)
              do i=1,list%nl
                kls%body(i+1)=ktfcopy(list%body(i))
              enddo
              call tfsetpart(kls,k2,kx,0,irtc)
              return
            endif
          else
            go to 9900
          endif
        elseif(ktflistq(ka,kla))then
          do while(ktflistq(ka,kla))
            ka=kla%head
          enddo
          if(.not. ktfsymbolqdef(ka%k,symd))then
            go to 9900
          endif
        else
          go to 9900
        endif
 12     if(symd%sym%override .eq. -2)then
          if(ktfprotectedqo(symd%sym) .and.
     $         symd%sym%gen .ne. -3)then
            irtc=itfmessage(999,'General::protect','""')
            return
          endif
          call tfdset(k2,symd%downval,kx,k1)
        else
          go to 9900
        endif
      elseif(ktfsymbolqd(k1,sym))then
        if(sym%override .eq. 0)then
          call tfsydef(sym,sym)
        endif
        if(ktfprotectedqo(sym))then
          irtc=itfmessage(999,'General::protect','""')
          return
        endif
        call sym_symdef(sym,symd)
        if(ktfrefqd(symd%value,kar))then
          if(ktfrealq(k2))then
            dlist(kar)=k2
            kx=k2
            return
          endif
          go to 9900
        endif
        if(k2%k .eq. ktfref)then
          ks%k=ktfcopy1(ktfsymbol+sad_loc(symd%sym%loc))
          kx%k=ktfoper+mtfnull
        else
          ks=dtfcopy(k2)
          kx=ks
        endif
        call tflocald(symd%value)
        symd%value=ks
      elseif(ktfrefqd(k1,ka1))then
c        if(ka1 .gt. 0 .and. ktfrealq(k2))then
        if(ka1 .gt. 0)then
          call tflocald(dlist(ka1))
          dlist(ka1)=dtfcopy(k2)
        else
          irtc=itfmessage(999,'General::invset',
     $         '"Illegal Momory Location"')
          return
        endif
      else
        go to 9900
      endif
      return
 9900 irtc=itfmessageexp(999,'General::invset',k10)
      return
      end

      subroutine tftagset(list,k2,kx,mopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ks,k2,kx,kl,kh
      type (sad_list) list
      type (sad_symbol), pointer :: syms
      type (sad_symdef),pointer :: symd
      type (sad_list), pointer :: list1,kls,kll,klh
      integer*4 irtc,mopc,itfmessage,itfmessageexp
      call tfeevaldef(list%body(1),ks,irtc)
      if(irtc .ne. 0)then
        return
      endif
      do while(ktflistq(ks,kls))
        ks=kls%head
      enddo
      if(.not. ktfsymbolqd(ks,syms))then
        irtc=itfmessage(9,'General::wrongtype','"Symbol"')
        return
      endif
      call tfeevaldef(list%dbody(2),kl,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ktfsymbolqd(kl))then
        if(ktfprotectedqo(syms))then
          irtc=itfmessage(9,'General::protect','""')
          return
        endif
        if(ks%k .eq. kl%k)then
          call tfset1(ks,k2,kx,mopc,irtc)
          return
        else
          irtc=itfmessage(9,'General::samesymbol',' ')
          return
        endif
      elseif(ktflistq(kl,kll))then
        if(ktfprotectedqo(syms) .and. syms%gen .ne. -3)then
          irtc=itfmessage(9,'General::protect','""')
          return
        endif
        kh=kll%head
        do while(ktflistq(kh,klh))
          kh=klh%head
        enddo
        call sym_symdef(syms,symd)
        if(kh%k .eq. ks%k)then
          call tfdset(k2,symd%downval,kx,kl)
        endif
        call descr_list(list%dbody(2),list1)
        call tfdsethead(list1,syms,kh)
        if(kh%k .ne. ktfref)then
          call tfdset(k2,symd%upval,kx,list%dbody(2))
        endif
      else
        irtc=itfmessageexp(9,'General::invset',kl)
        return
      endif
      kx=k2
      return
      end

      subroutine tfdsethead(list,sym,kh)
      use tfstk
      implicit none
      type (sad_descriptor) kh,ki,ki0
      type (sad_list) list
      type (sad_list), pointer :: kli
      type (sad_symbol) sym
      type (sad_symbol), pointer :: symi
      integer*4 i
      kh%k=ktfref
      do i=1,list%nl
        ki=list%dbody(i)
        ki0=ki
        if(ktflistq(ki,kli))then
          ki=kli%head
          if(ktflistq(ki))then
            ki0=ki
            do while(ktflistq(ki,kli))
              ki=kli%head
            enddo
          endif
        endif
        if(ktfsymbolqd(ki,symi))then
c          write(*,*)'dsethead ',symi%loc,symi%gen
          if(symi%loc .eq. sym%loc
     $         .and. max(0,symi%gen) .eq. max(0,sym%gen))then
            kh=ki0
            return
          endif
        endif
      enddo
      return
      end
