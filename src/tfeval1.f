      recursive function tfeval1(k1,k2,iopc1,irtc) result(kx)
      use tfstk
      use eexpr
      use efun
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      integer*4 isp1
      select case (iopc1)
      case (mtfrevpower)
        kx=tfeval1(k2,k1,mtfpower,irtc)
      case (mtfcomp)
        irtc=0
        kx=k2
      case (mtffun:mtfruledelayed,mtfpattest,mtfslot,mtfslotseq,
     $       mtfalt,mtfrepeated,mtfrepeatednull)
        kx=tfeexpr(k1,k2,iopc1)
        irtc=0
      case (mtfnull:mtfdiv,mtfpower,mtfgreater:mtfless,mtfand:mtfnot,
     $       mtfleftbra:mtfrightbrace,
     $       mtfcomplex:mtfcomma,mtfdot,mtfend)
        if(tfnumberq(k2) .and.
     $     (tfnumberq(k1) .or. iopc1 .eq. mtfnot))then
          kx=tfcmplx(k1,k2,iopc1,irtc)
          return
        endif
        if(tflistq(k2) .or. tflistq(k1))then
          kx=tfearray(k1,k2,iopc1,irtc)
        else
          kx=tfeexpr(k1,k2,iopc1)
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
        kx=tfefunref(isp1,.true.,irtc)
        isp=isp1-1
      end select
      return
      end

      subroutine tfflagordef(isp1,kx,irtc)
      use tfstk
      type (sad_descriptor) ,intent(out):: kx
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: isp1
      integer*4 nc
      real*8 vx,fflogi
      logical*4 exist
      character*8 name
      if(isp .ne. isp1+1)then
        call tfdefinition(isp1,kx,irtc)
        return
      endif
      if(ktfsymbolq(dtastk(isp),sym))then
        call sym_symstr(sym,str)
        nc=min(8,str%nch)
        name(1:nc)=str%str(1:nc)
        call capita(name(1:nc))
        vx=fflogi(name(1:nc),exist)
c        write(*,*)'flagordef-vx ',exist,vx,'"'//name(1:nc)//'"'
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

      function tfeval1to(k1,k2,iopc,old,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,kv,kr,ku,ks,tfeval1,tfset1
      type (sad_dlist), pointer :: kl1
      integer*4 ,intent(in):: iopc
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: old
      if(ktflistq(k1,kl1))then
        kv=tfleval(kl1,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      elseif(ktfsymbolq(k1))then
        kv=tfsyeval(k1,irtc)
        if(irtc .ne. 0)then
          return
        endif
      elseif(ktfpatq(k1))then
        kv=tfpateval(k1,irtc)
        if(irtc .ne. 0)then
          return
        endif
      else
        kv=k1
      endif
      if(iopc .eq. mtfaddto)then
        kr=tfeval1(kv,k2,mtfplus,irtc)
      elseif(iopc .eq. mtftimesby)then
        kr=tfeval1(kv,k2,mtfmult,irtc)
      elseif(iopc .eq. mtfsubtractfrom)then
        ku=tfeval1(sad_descr(-1.d0),k2,mtfmult,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kr=tfeval1(kv,ku,mtfplus,irtc)
      else
        ku=tfeval1(k2,sad_descr(-1.d0),mtfpower,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kr=tfeval1(kv,ku,mtfmult,irtc)
      endif
      if(irtc .ne. 0)then
        return
      endif
      ks=tfeevaldef(k1,irtc)
      if(irtc .ne. 0)then
        return
      endif
      kx=tfset1(ks,kr,mtfset,irtc)
      if(old)then
        kx=kv
      endif
      return
      end

      subroutine tfappendto(isp1,kx,mode,irtc)
      use tfstk
      use efun
      use eexpr
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) k10,kr,k1
      type (sad_dlist), pointer :: kl
      integer*2, parameter :: nextra = int2(8)
      integer*8 kp
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,isp0
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
            kx=tfleval(kl,.false.,irtc)
            if(irtc .ne. 0)then
              return
            endif
            call loc_sad(ktfaddrd(kx),kl)
            call tfsetpart(kl,dtastk(isp),kx,mode,irtc)
            return
          endif
        endif
        k1=tfeevalref(k10,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfnonlistq(k1))then
          irtc=itfmessage(9,'General::wrongtype','"List"')
          return
        endif
        kr=tfappend(k1,dtastk(isp),.true.,merge(0,1,mode .eq. 1),irtc)
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
      kx=tfefunref(isp0+1,.true.,irtc)
      isp=isp0
      return
      end

      subroutine tfappendto1(kp,k2,kr,mode,eval,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k2
      type (sad_descriptor) ,intent(out):: kr
      type (sad_descriptor) k1
      type (sad_dlist), pointer :: kl,klr
      integer*2, parameter :: nextra = int2(8)
      integer*8 ,intent(in):: kp
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i,n,mode
      logical*4 ,intent(out):: eval
      logical*4 ov
      eval=.true.
      k1=dlist(kp)
      if(ktfnonlistq(k1,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      n=kl%nl
      ov=ktfovrwrtq(kl) .and. kl%ref .eq. 1
      if(ov .and. mode .eq. 1 .and. kl%lena .gt. 0)then
        klr=>kl
        klr%nl=n+1
        klr%lena=klr%lena-int2(1)
        klr%dbody(n+1)%k=0
        call tfreplist(klr,n+1,k2,eval)
        if(.not. eval)then
          return
        endif
      elseif(ov .and. mode .eq. 2 .and. kl%lenp .gt. 0)then
        call loc_sad(ktfaddr(k1)-1,klr)
        klr%dbody(-3:0)=kl%dbody(-3:0)
        klr%lenp=klr%lenp-1
        klr%nl=n+1
        klr%dbody(1)%k=0
        call tfreplist(klr,1,k2,eval)
        klist(kp)=ktflist+ktfaddr(k1)-1
        if(.not. eval)then
          return
        endif
      elseif(mode .eq. 1)then
        kr%k=ktaalocsp(n+1,kl%lenp,nextra,klr)
        if(ktfreallistq(kl))then
          klr%head=dtfcopy(kl%head)
          klr%dbody(1:n)=kl%dbody(1:n)
        else
          klr%attr=lnonreallist
          do i=0,n
            klr%dbody(i)=dtfcopy(kl%dbody(i))
          enddo
        endif
        klr%dbody(n+1)%k=0
        call tfreplist(klr,n+1,k2,eval)
      else
        kr%k=ktaalocsp(n+1,nextra,kl%lena,klr)
        klr%head=dtfcopy(kl%head)
        if(ktfreallistq(kl))then
          klr%dbody(2:n+1)=kl%dbody(1:n)
        else
          klr%attr=lnonreallist
          do i=1,n
            klr%dbody(i+1)=dtfcopy(kl%dbody(i))
          enddo
        endif
        klr%dbody(1)%k=0
        call tfreplist(klr,1,k2,eval)
      endif
      if(eval)then
        kr=tfleval(klr,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      else
        kr=sad_descr(klr)
      endif
      eval=.true.
      return
      end

      logical*4 function tfgetstoredp(ks0,kp,def,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: ks0
      type (sad_descriptor) ks
      type (sad_dlist), pointer :: lists
      type (sad_symdef), pointer :: symd
      integer*8 kp
      integer*4 ,intent(out):: irtc
      integer*4 itfmessageexp,itfmessage
      logical*4 ,intent(out):: def
      logical*4 ev
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
          call loc_sad(ktfaddrd(ks0),lists)
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

      function  tfset1(k10,k20,mopc,irtc) result(kx)
      use tfstk
      use eexpr
      use mackw
      implicit none
      type (sad_descriptor) k1,k2,kx,ks,ka
      type (sad_descriptor) ,intent(in):: k10,k20
      type (sad_dlist),pointer :: list,kls1,kla,kls
      type (sad_symbol), pointer ::sym
      type (sad_symdef),pointer :: symd
      integer*8 ka1,kaa,kas,kar
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mopc
      integer*4 itfmessage,itfmessageexp,i
      irtc=0
      k1=k10
      k2=k20
      kx=k2
      if(ktflistq(k1,list))then
        ka=list%head
        if(ktfoperq(ka,kaa))then
          select case (kaa)
          case (mtflist)
            kx=tfearray(k1,k2,mtfset,irtc)
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
                kls%dbody(i+1)=dtfcopy(list%dbody(i))
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
      elseif(ktfsymbolq(k1,sym))then
        if(sym%override .eq. 0)then
          sym=>tfsydef(sym)
        endif
        if(ktfprotectedqo(sym))then
          irtc=itfmessage(999,'General::protect','""')
          return
        endif
        call sym_symdef(sym,symd)
        if(ktfrefq(symd%value,kar))then
          if(ktfrealq(k2))then
            dlist(kar)=k2
            kx=k2
            return
          endif
          go to 9900
        endif
        if(k2%k .eq. ktfref)then
          ks=dtfcopy1(sad_descr(symd%sym))
          kx%k=ktfoper+mtfnull
        else
          ks=dtfcopy(k2)
          kx=ks
        endif
        call tflocald(symd%value)
        symd%value=ks
      elseif(ktfrefq(k1,ka1))then
c        if(ka1 .gt. 0 .and. ktfrealq(k2))then
        if(ka1 .gt. 0)then
          call tflocald(dlist(ka1))
          dlist(ka1)=dtfcopy(k2)
        else
          irtc=itfmessage(999,'General::invset',
     $         '"Illegal Memory Location"')
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
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k2
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ks,kl,kh,tfset1
      type (sad_dlist) ,intent(in):: list
      type (sad_symbol), pointer :: syms
      type (sad_symdef),pointer :: symd
      type (sad_dlist), pointer :: list1,kls,kll,klh
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mopc
      integer*4 itfmessage,itfmessageexp
      ks=tfeevaldef(list%dbody(1),irtc)
      if(irtc .ne. 0)then
        return
      endif
      do while(ktflistq(ks,kls))
        ks=kls%head
      enddo
      if(.not. ktfsymbolq(ks,syms))then
        irtc=itfmessage(9,'General::wrongtype','"Symbol"')
        return
      endif
      kl=tfeevaldef(list%dbody(2),irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ktfsymbolq(kl))then
        if(ktfprotectedqo(syms))then
          irtc=itfmessage(9,'General::protect','""')
          return
        endif
        if(ks%k .eq. kl%k)then
          kx=tfset1(ks,k2,mopc,irtc)
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
        call descr_sad(list%dbody(2),list1)
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
      type (sad_descriptor) ,intent(out):: kh
      type (sad_descriptor) ki,ki0
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: kli
      type (sad_symbol) ,intent(in):: sym
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
        if(ktfsymbolq(ki,symi))then
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
