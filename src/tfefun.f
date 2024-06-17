      module eexpr
      use tfstk
      use complex
      type (sad_descriptor),external::tfset1

      contains
      function tfinsertsort(kl,ki) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: ki
      type (sad_dlist) ,intent(inout):: kl
      integer*4 isp0,isp1,isp2,i,ispm,itfcanonicalorder,isp3
      isp0=isp
      call tfgetllstkall(kl)
      isp1=isp0+1
      isp2=isp+1
      do while(isp2 > isp1)
        ispm=isp1+(isp2-isp1)/2
        i=itfcanonicalorder(ki,dtastk(ispm))
        if(i > 0)then
          isp1=ispm+1
        elseif(i == 0)then
          isp1=ispm
          isp2=ispm
        else
          isp2=ispm
        endif
      enddo
      isp2=isp
      isp3=isp+isp1-isp0-1
      ktastk(isp+1:isp3)=ktastk(isp0+1:isp1-1)
      isp=isp3+1
      dtastk(isp)=ki
      ktastk(isp+1:isp+isp2-isp1+1)=ktastk(isp1:isp2)
      isp=isp+isp2-isp1+1
      kx=kxcrelistm(isp-isp2,ktastk(isp2+1:isp),kl%head)
      isp=isp0
      return
      end function

      function tfappend(kl,k,eval,mode,irtc) result(kx)
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: kl,k
      type (sad_dlist), pointer :: list,listx
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mode
      integer*4 m,itfmessage,i
      logical*4 ,intent(in):: eval
      logical*4 ev
      if(.not. ktflistq(kl,list))then
        kx=dxnull
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      ev=eval .and. list%head%k /= ktfoper+mtflist .and.
     $     list%head%k /= ktfoper+mtfalt .and.
     $     list%head%k /= ktfoper+mtfnull
      m=list%nl
      call loc_sad(ktaaloc(-1,m+1),listx)
      listx%attr=list%attr
      if(ktfreallistq(list))then
        if(mode == 0)then
          listx%dbody(1:m)=list%dbody(1:m)
          if(ktfrealq(k))then
            listx%dbody(m+1)=k
          else
            listx%dbody(m+1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        else
          listx%dbody(2:m+1)=list%dbody(1:m)
          if(ktfrealq(k))then
            listx%dbody(1)=k
          else
            listx%dbody(1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        endif
      else
        if(mode == 0)then
          do i=1,m
            listx%dbody(i)=dtfcopy(list%dbody(i))
          enddo
          listx%dbody(m+1)=dtfcopy(k)
        else
          do i=1,m
            listx%dbody(i+1)=dtfcopy(list%dbody(i))
          enddo
          listx%dbody(1)=dtfcopy(k)
        endif
      endif
      listx%head=dtfcopy(list%head)
      if(iand(list%attr,kconstarg) /= 0)then
        if(.not. tfconstq(k%k))then
          listx%attr=ior(listx%attr-kconstarg,knoconstarg+lnoconstlist)
        endif
      endif
      if(ev)then
        kx=tfleval(listx,.true.,irtc)
      else
        kx=sad_descr(listx)
        irtc=0
      endif
      return
      end function

      subroutine tfappendto(isp1,kx,mode,irtc)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) k10,kr,k1,tfefunrefu
      type (sad_dlist), pointer :: kl
      integer*2, parameter :: nextra = int2(8)
      integer*8 kp
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,isp0
      logical*4 def,ev,tfgetstoredp
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      k10=dtastk(isp1+1)
      if(tfgetstoredp(k10,kp,def,irtc))then
        call tfappendto1(kp,dtastk(isp),kr,mode,ev,irtc)
        if(irtc /= 0)then
          return
        elseif(.not. ev)then
          kx=dlist(kp)
          return
        endif
      else
        if(irtc > 0)then
          return
        endif
        if(ktflistq(k10,kl))then
          if(kl%head%k == ktfoper+mtfpart)then
            kx=tfleval(kl,.false.,irtc)
            if(irtc /= 0)then
              return
            endif
            call loc_sad(ktfaddrd(kx),kl)
            call tfsetpart(kl,dtastk(isp),kx,mode,irtc)
            return
          endif
        endif
        k1=tfeevalref(k10,irtc)
        if(irtc /= 0)then
          return
        endif
        if(ktfnonlistq(k1))then
          irtc=itfmessage(9,'General::wrongtype','"List"')
          return
        endif
        kr=tfappend(k1,dtastk(isp),.true.,merge(0,1,mode == 1),irtc)
        if(irtc /= 0)then
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
      kx=tfefunrefu(isp0+1,irtc)
      isp=isp0
      return
      end

      recursive function tfset(isp1,upvalue,irtc) result(kx)
      use eeval,only:tfeevaldef
      implicit none
      type (sad_descriptor) kx,k1,k2,k110,k11
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_dlist), pointer :: kl11
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage,isp11,isp0
      logical*4 ,intent(in):: upvalue
      logical*4 euv
      if(isp1+1 .ge. isp)then
        kx=dxnull
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      k2=dtastk(isp)
      if(isp1+2 == isp)then
        k1=dtastk(isp1+1)
        if(ktfsymbolq(k1,sym))then
          if(sym%override == 0)then
            sym=>tfsydef(sym)
            k1=sad_descr(sym)
          endif
        elseif(ktflistq(k1))then
          k1=tfeevaldef(k1,irtc)
          if(irtc /= 0)then
            return
          endif
        endif
        if(k1%k /= ktastk(isp1+1) .and. upvalue)then
          k11=k1
          k110=k1
          if(ktflistq(k11,kl11))then
            k11=kl11%head
            if(ktflistq(k11))then
              k110=k11
              do while(ktflistq(k11,kl11))
                k11=kl11%head
              enddo
            endif
          endif
          if(ktfsymbolq(k11,sym))then
            if(sym%override == 0)then
              sym=>tfsydef(sym)
            endif
            call sym_symdef(sym,symd)
            if(symd%upval /= 0)then
              isp=isp+1
              isp11=isp
              dtastk(isp11)=dtastk(isp1)
              isp=isp+1
              dtastk(isp)=k1
              isp=isp+1
              dtastk(isp)=k2
              call tfdeval(isp11,dfromk(ksad_loc(sym%loc)),kx,0,
     $             .false.,euv,irtc)
              isp=isp11-1
              if(euv)then
                return
              endif
            endif
          endif
        endif
        kx=tfset1(k1,k2,int(ktfaddr(ktastk(isp1))),irtc)
      else
        isp0=isp
        kx=dtastk(isp)
        do i=isp-1,isp1+1,-1
          isp=isp0+3
          dtastk(isp-2)=dtastk(isp1)
          dtastk(isp-1)=dtastk(i)
          dtastk(isp  )=kx
          kx=tfset(isp0+1,upvalue,irtc)
          if(irtc /= 0)then
            return
          endif
        enddo
        isp=isp0
      endif
      return
      end

      subroutine tftagset(list,k2,kx,mopc,irtc)
      use eeval,only:tfeevaldef
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
      if(irtc /= 0)then
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
      if(irtc /= 0)then
        return
      endif
      if(ktfsymbolq(kl))then
        if(ktfprotectedqo(syms))then
          irtc=itfmessage(9,'General::protect','""')
          return
        endif
        if(ks%k == kl%k)then
          kx=tfset1(ks,k2,mopc,irtc)
          return
        else
          irtc=itfmessage(9,'General::samesymbol',' ')
          return
        endif
      elseif(ktflistq(kl,kll))then
        if(ktfprotectedqo(syms) .and. syms%gen /= -3)then
          irtc=itfmessage(9,'General::protect','""')
          return
        endif
        kh=kll%head
        do while(ktflistq(kh,klh))
          kh=klh%head
        enddo
        call sym_symdef(syms,symd)
        if(kh%k == ks%k)then
          call tfdset(k2,symd%downval,kx,kl)
        endif
        call descr_sad(list%dbody(2),list1)
        call tfdsethead(list1,syms,kh)
        if(kh%k /= ktfref)then
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
          if(symi%loc == sym%loc
     $         .and. max(0,symi%gen) == max(0,sym%gen))then
            kh=ki0
            return
          endif
        endif
      enddo
      return
      end

      subroutine tfsetpart(kln,k2,kx,mode,irtc)
      use part,only:tfpartrstk,tfreplist
      use eeval
      use funs,only:tfgetllstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ,intent(in):: k2
      integer*4 ,intent(in):: mode
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kr,k1
      type (sad_dlist) ,intent(in):: kln
      type (sad_dlist), pointer :: listl,list2,kli
      integer*8 kap
      integer*4 i,isp1,isp2,n,itfmessageexp
      logical*4 list,def,eval,evala,tfgetstoredp
      if(tfgetstoredp(kln%dbody(1),kap,def,irtc))then
        k1=dlist(kap)
        if(ktfnonlistq(k1,listl))then
          irtc=itfmessageexp(999,'General::invset',k1)
          return
        endif
      else
        if(irtc == 0)then
          irtc=itfmessageexp(999,'General::invset',sad_descr(kln))
        endif
        return
      endif
      isp1=isp
      call tfgetllstk(kln,2,-1)
      isp2=isp
      listl=>tfclonelist(listl)
      call tfpartrstk(listl,isp1,isp2,list,.false.,.true.,eval,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      if(mode == 0)then
        if(list)then
          if(tflistq(k2,list2))then
            n=list2%nl
            if(n == isp-isp2)then
              do i=1,n
                call loc_sad(ktastk(isp2+i),kli)
                call tfreplist(kli,itastk2(1,isp2+i),list2%dbody(i),evala)
                eval=eval .or. evala
              enddo
              go to 100
            endif
          endif
          do i=isp2+1,isp
            call loc_sad(ktastk(i),kli)
            call tfreplist(kli,itastk2(1,i),k2,evala)
            eval=eval .or. evala
          enddo
        elseif(isp == isp2+1)then
          call loc_sad(ktastk(isp),kli)
          call tfreplist(kli,itastk2(1,isp),k2,evala)
          eval=eval .or. evala
        endif
      else
        if(list)then
          if(tflistq(k2,list2))then
            n=list2%nl
            if(n == isp-isp2)then
              do i=1,n
                call tfappendto1(ktastk(isp2+i)+itastk2(1,isp2+i),
     $               list2%dbody(i),kr,mode,evala,irtc)
                if(irtc /= 0)then
                  go to 1000
                endif
                if(evala)then
                  call loc_sad(ktastk(isp2+i),kli)
                  call tfreplist(kli,itastk2(1,isp2+i),kr,evala)
                  eval=eval .or. evala
                endif
              enddo
              go to 100
            endif
          endif
          do i=1,isp-isp2
            call tfappendto1(ktastk(isp2+i)+itastk2(1,isp2+i),
     $           k2,kr,mode,evala,irtc)
            if(irtc /= 0)then
              go to 1000
            endif
            if(evala)then
              call loc_sad(ktastk(isp2+i),kli)
              call tfreplist(kli,itastk2(1,isp2+i),kr,evala)
              eval=eval .or. evala
            endif
          enddo
        elseif(isp == isp2+1)then
          call tfappendto1(ktastk(isp)+itastk2(1,isp),
     $         k2,kr,mode,evala,irtc)
          if(irtc /= 0)then
            go to 1000
          endif
          if(evala)then
            call loc_sad(ktastk(isp),kli)
            call tfreplist(kli,itastk2(1,isp),kr,evala)
            eval=eval .or. evala
          endif
        endif
      endif
 100  if(eval)then
        kx=tfleval(listl,.true.,irtc)
        if(irtc /= 0)then
          go to 1000
        endif
      else
        kx=sad_descr(listl)
      endif
      call tflocald(dlist(kap))
      dlist(kap)=dtfcopy(kx)
 1000 isp=isp1
      kx=k2
      return
      end

      subroutine tfappendto1(kp,k2,kr,mode,eval,irtc)
      use part,only:tfreplist
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
      ov=ktfovrwrtq(kl) .and. kl%ref == 1
      if(ov .and. mode == 1 .and. kl%lena > 0)then
        klr=>kl
        klr%nl=n+1
        klr%lena=klr%lena-int2(1)
        klr%dbody(n+1)%k=0
        call tfreplist(klr,n+1,k2,eval)
        if(.not. eval)then
          return
        endif
      elseif(ov .and. mode == 2 .and. kl%lenp > 0)then
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
      elseif(mode == 1)then
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
        if(irtc /= 0)then
          return
        endif
      else
        kr=sad_descr(klr)
      endif
      eval=.true.
      return
      end

      function tfnot(isp1,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,tfearray
      integer*4 isp1,irtc,iopc,itfmessage
      if(isp /= isp1+1)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(tfnumberq(dtastk(isp)))then
        kx=tfcmplx(dfromr(0.d0),dtastk(isp),iopc,irtc)
      elseif(tflistq(dtastk(isp)))then
        kx=tfearray(dfromr(0.d0),dtastk(isp),iopc,irtc)
      else
        kx=tfeexpr(dfromr(0.d0),dtastk(isp),iopc)
        irtc=0
      endif
      return
      end

      function tfplus(isp1,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,k1,k,ki
      type (sad_dlist), pointer :: klx
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 i,narg
      real*8 v,v1,vx,vi
      kx=dxnullo
      narg=isp-isp1
      if(narg == 2)then
        k1=dtastk(isp1+1)
        k =dtastk(isp)
        if(ktfrealq(k1,v1) .and. ktfrealq(k,v))then
          kx=merge(dfromr(v1+v),dfromr(v1*v),iopc == mtfplus)    
          irtc=0
        else
          if(tfnumberq(k) .and. tfnumberq(k1))then
            kx=tfcmplx(k1,k,iopc,irtc)
          else
            if(tflistq(k1) .or. tflistq(k))then
              kx=tfecmplxl(k1,k,iopc)
            else
              kx=tfeexpr(k1,k,iopc)
            endif
            irtc=0
          endif
        endif
        return
      elseif(narg == 0)then
        kx%k=merge(ktffalse,ktftrue,iopc == mtfplus)
      elseif(narg == 1)then
        if(ktfsymbolq(dtastk(isp)) .or.
     $       ktfpatq(dtastk(isp)))then
          kx=kxmakelist(isp1,klx)
          klx%head%k=ktfoper+iopc
        else
          kx=dtastk(isp1+1)
        endif          
        irtc=0
      else
        kx=dtastk(isp1+1)
        irtc=0
        do i=isp1+2,isp
          ki=dtastk(i)
          if(ktfrealq(ki,vi) .and. ktfrealq(kx,vx))then
            kx=merge(dfromr(vx+vi),dfromr(vx*vi),iopc == mtfplus)    
          else
            k1=kx
            if(tfnumberq(k1) .and. tfnumberq(ki))then
              kx=tfcmplx(k1,ki,iopc,irtc)
              if(irtc /= 0)then
                return
              endif
            elseif(tflistq(k1) .or. tflistq(ki))then
              kx=tfecmplxl(k1,ki,iopc)
              if(irtc /= 0)then
                return
              endif
            else
              kx=tfeexpr(k1,ki,iopc)
            endif
          endif
        enddo
      endif
      return
      end

      function tfpower(isp1,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,k1,ki
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      if(isp1 == isp)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      kx=dtastk(isp)
      irtc=0
      if(isp == isp1+1)then
        return
      endif
      do i=isp-1,isp1+1,-1
        ki=dtastk(i)
        k1=kx
        if(tfnumberq(k1) .and. tfnumberq(ki))then
          kx=tfcmplx(ki,k1,mtfpower,irtc)
          if(irtc /= 0)then
            return
          endif
        elseif(tflistq(k1) .or. tflistq(ki))then
          kx=tfecmplxl(ki,k1,mtfpower)
          irtc=0
        else
          kx=tfeexpr(ki,k1,mtfpower)
        endif
      enddo
      return
      end

      function tfrevpower(isp1,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,k1,ki
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      if(isp1 == isp)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      kx=dtastk(isp1+1)
      irtc=0
      if(isp == isp1+1)then
        return
      endif
      do i=isp1+2,isp
        ki=dtastk(i)
        k1=kx
        if(tfnumberq(k1) .and. tfnumberq(ki))then
          kx=tfcmplx(ki,k1,mtfpower,irtc)
          if(irtc /= 0)then
            return
          endif
        elseif(tflistq(k1) .or. tflistq(ki))then
          kx=tfecmplxl(ki,k1,mtfpower)
          irtc=0
        else
          kx=tfeexpr(ki,k1,mtfpower)
        endif
      enddo
      return
      end

      function tfequal(isp1,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp .lt. isp1+2)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(isp == isp1+2 .and.
     $     ktfstringq(dtastk(isp)) .and. ktfstringq(dtastk(isp1+1)))then
        kx%k=merge(ktftrue,ktffalse,
     $       tfsamestringq(dtastk(isp),dtastk(isp1+1)))
        if(iopc == mtfunequal)then
          kx%k=ktftrue-kx%k
        endif
        irtc=0
      else
        kx=tfrelation(isp1,iopc,irtc)
      endif
      return
      end

      recursive function tfrelation(isp1,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,tfearray
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,isp0,k
      kx=dxnullo
      if(isp .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(isp == isp1+2)then
        if(tfnumberq(dtastk(isp1+1)) .and.
     $       tfnumberq(dtastk(isp)))then
          kx=tfcmplx(dtastk(isp1+1),dtastk(isp),iopc,irtc)
        elseif(tflistq(dtastk(isp1+1))
     $         .or. tflistq(dtastk(isp)))then
          kx=tfearray(dtastk(isp1+1),dtastk(isp),iopc,irtc)
        else
          kx=tfeexpr(dtastk(isp1+1),dtastk(isp),iopc)
          irtc=0
        endif
      else
        isp0=isp
        do k=1,isp0-isp1-1
          ktastk(isp0+1)=ktastk(isp1+k)
          ktastk(isp0+2)=ktastk(isp1+k+1)
          isp=isp0+2
          kx=tfrelation(isp0,iopc,irtc)
          if(irtc /= 0)then
            isp=isp0
            return
          elseif(ktfnonrealq(kx))then
            irtc=-1
            return
          endif
          if(kx%k == 0)then
            isp=isp0
            return
          endif
        enddo
        isp=isp0
      endif
      return
      end

      function tfjoine(k1,k2,irtc) result(kx)
      use funs,only:tfgetllstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) k10,k20,ky1
      type (sad_dlist), pointer ::kl1,kl2
      integer*4 ,intent(out):: irtc
      integer*4 ma1,ma2,m,iopc,isp1
      if(ktfnonlistq(k1,kl1) .or. ktfnonlistq(k2,kl2)
     $     .or. .not. tfsameheadq(k1,k2))then
        kx=dxnullo
        irtc=-1
        return
      endif
      irtc=0
      iopc=int(ktfaddr(kl1%head%k))
      if(iopc == mtfplus .or. iopc == mtftimes)then
        if(.not. tfnumberq(kl1%dbody(1)))then
          kx=tfjoin2(k2,k1,.false.,irtc)
          go to 1000
        endif
        if(.not. tfnumberq(kl2%dbody(1)))then
          kx=tfjoin2(k1,k2,.false.,irtc)
          go to 1000
        endif
        ma1=kl1%nl
        ma2=kl2%nl
        m=ma1+ma2-1
        k10=kl1%dbody(1)
        k20=kl2%dbody(1)
        ky1=tfcmplx(k10,k20,iopc,irtc)
        if(irtc /= 0)then
          kx=dxnullo
          return
        endif
        if(ktfrealq(ky1))then
          if(iopc == mtfplus)then
            if(ky1%k == 0 .or. ky1%k == ktfmzero)then
              m=m-1
            endif
          else
            if(ky1%k == 0 .or. ky1%k == ktfmzero)then
              kx%k=0
              return
            elseif(ky1%k == ktftrue)then
              m=m-1
            endif
          endif
        endif
        if(m == 1)then
          kx=kl2%dbody(2)
          return
        endif
        isp1=isp
        if(m == ma1+ma2-2)then
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        else
          isp=isp+1
          dtastk(isp)=ky1
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        endif
        kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),k_descr(ktfoper+iopc))
        isp=isp1
      else
        kx=tfjoin2(k1,k2,.false.,irtc)
        return
      endif
 1000 isp=isp+1
      dtastk(isp)=kx
      call tfsort(isp-1,kx,0,irtc)
      isp=isp-1
      return
      end

      function tfjoin2(k1,k2,eval,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      integer*4 ,intent(out):: irtc
      integer*4 isp1
      logical*4 ,intent(in):: eval
      isp1=isp
      isp=isp+1
      dtastk(isp)=k1
      isp=isp+1
      dtastk(isp)=k2
      kx=tfjoin(isp1,eval,irtc)
      isp=isp1
      return
      end function

      function tfjoin(isp1,eval,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,tfefunrefu
      type (sad_descriptor) kf
      type (sad_dlist), pointer :: kl1,kli
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i,narg,isp0
      logical*4 ,intent(in):: eval
      logical*4 ev
      narg=isp-isp1
      if(ktfnonlistq(ktastk(isp1+1),kl1))then
        go to 9010
      endif
      if(narg .le. 1)then
        if(narg == 1)then
          kx=dtastk(isp1+1)
          irtc=0
        else
          kx=dxnullo
          irtc=itfmessage(9,'General::narg','"1 or more"')
        endif
        return
      endif
      kf=kl1%head
      isp=isp+1
      isp0=isp
      call tfgetllstkall(kl1)
      if(isp .ge. mstk)then
        kx=dxnullo
        irtc=itfmessage(9999,'General::stack','"Join"')
        return
      endif
      do i=isp1+2,isp0-1
        if(ktfnonlistq(ktastk(i),kli))then
          isp=isp0-1
          go to 9000
        endif
        if(.not. tfsameq(kli%head,kf))then
          go to 9100
        endif
        call tfgetllstkall(kli)
        if(isp .ge. mstk)then
          kx=dxnullo
          irtc=itfmessage(9999,'General::stack','"Join"')
          return
        endif
      enddo
      dtastk(isp0)=kf
      ev=eval
      if(ev .and. (kf%k == ktfoper+mtflist .or.
     $     kf%k == ktfoper+mtfalt .or.
     $     kf%k == ktfoper+mtfnull))then
          ev=.false.
      endif
      if(ev)then
        kx=tfefunrefu(isp0,irtc)
      else
        kx=kxcompose(isp0)
        irtc=0
      endif
      isp=isp0-1
      return
 9000 isp=isp0-1
 9010 irtc=itfmessage(9,'General::wrongtype',
     $     '"List or composition for all args"')
      kx=dxnullo
      return
 9100 irtc=itfmessage(9,'General::samehead',' ')
      kx=dxnullo
      isp=isp0-1
      return
      end

      subroutine tfclear(isp1,kx,irtc)
      use dset,only:tfcleardaloc
      use modul,only:tfdelete
      implicit none
      type (sad_descriptor) kx,ki
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_dlist), pointer :: kl
      integer*4 isp1,irtc,i,itfmessage
      LOOP_I: do i=isp1+1,isp
        ki=dtastk(i)
        if(ktfsymbolq(ki,sym))then
          if(sym%override == 0)then
            sym=>tfsydef(sym)
          endif
          if(ktfprotectedqo(sym))then
            irtc=itfmessage(999,'General::protect','""')
            return
          endif
          call sym_symdef(sym,symd)
          call tfdelete(symd,.false.,.false.)
        elseif(ktflistq(ki,kl))then
          ki=dtfcopy1(ki)
          call tfcleardef(kl,irtc)
          call tflocal1(ki%k)
          if(irtc /= 0)then
            return
          endif
        endif
      enddo LOOP_I
      kx%k=ktfoper+mtfnull
      irtc=0
      return
      end

      subroutine tfcleardef(kl,irtc)
      use tfcode
      use funs
      use eeval
      use dset,only:tfcleardaloc
      use iso_c_binding
      implicit none
      type (sad_descriptor) kh,kx
      type (sad_dlist) ,intent(inout):: kl
      type (sad_dlist), pointer :: klh
      type (sad_symbol), pointer :: symh
      type (sad_symdef), pointer :: def
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash
c      include 'DEBUG.inc'
      integer*8 kad,kad1,kadi,kadi1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,ik,i
      logical*4 discard,discard1
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfoper+mtfset
      isp=isp+1
      dtastk(isp)=sad_descr(kl)
      isp=isp+1
      ktastk(isp)=ktfref
      kx=tfset(isp0+1,.false.,irtc)
      isp=isp0
      if(irtc /= 0)then
        return
      endif
      kh=kl%head
      do while(ktflistq(kh,klh))
        kh=klh%head
      enddo
      if(ktfsymbolqdef(kh%k,def))then
        if(def%sym%override == 0)then
          symh=>tfsydef(def%sym)
          call sym_symdef(symh,def)
        endif
        kad=def%upval
        discard=.true.
        do ik=1,2
          do while(kad /= 0)
            call loc_defhash(kad,dhash)
            kad1=dhash%next
            discard1=.true.
            if(dhash%gen == maxgeneration)then
              do i=0,dhash%nhash
                kadi=dhash%dhash(i)%k
                do while(kadi /= 0)
                  call loc_deftbl(kadi,dtbl)
                  kadi1=dtbl%next
                  if(tfdefheadq(dtbl%arg,kl))then
                    call tfcleardaloc(kadi)
                    if(kadi1 /= 0)then
                      klist(kadi1+1)=dtbl%prev
                    endif
                    klist(dtbl%prev)=kadi1
                    call tfree(kadi)
                  else
                    discard1=.false.
                  endif
                  kadi=kadi1
                enddo
              enddo
            else
              call loc_deftbl(kad,dtbl)
              if(tfdefheadq(dtbl%arg,kl))then
                call tfcleardaloc(kad)
              else
                discard1=.false.
              endif
            endif
            if(discard1)then
              if(kad1 /= 0)then
                klist(kad1+1)=dhash%prev
              endif
              klist(dhash%prev)=kad1
              call tfree(kad)
            else
              discard=.false.
            endif
            kad=kad1
          enddo
          kad=def%downval
          if(discard)then
            if(ik == 1)then
              def%upval=0
            else
              def%downval=0
            endif
          endif
        enddo
      endif
      irtc=0
      return
      end

      logical function tfdefheadq(k,listp)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,intent(out):: listp
      type (sad_dlist), pointer :: kli,kla1,kla
      integer*4 i
      tfdefheadq=.true.
      call descr_sad(k,kla)
      if(tfsamelistqo(kla,listp))then
        return
      endif
      kla1=>kla
      do while(ktflistq(kla1%head,kla1))
        if(tfsamelistqo(kla1,listp))then
          return
        endif
      enddo
      do i=1,kla%nl
        if(ktflistq(kla%dbody(i),kli))then
          if(tfsamelistqo(kli,listp))then
            return
          endif
          do while(ktflistq(kli%head,kli))
            if(tfsamelistqo(kli,listp))then
              return
            endif
          enddo
        endif
      enddo
      tfdefheadq=.false.
      return
      end

      end module eexpr

      recursive function tfeexpr(k1,k,iopc1) result(ke)
      use tfstk
      use eexpr,only:tfjoine,tfcmplx,tfinsertsort,tfjoin2,tfappend
      use part,only:tfreplist
      use take,only:tftake
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k
      type (sad_descriptor) ke,kx,ky,ky1,k2,kx2,tfcompose
      type (sad_dlist), pointer :: listy,list1,listi,klx,kl1
      type (sad_rlist), pointer :: kle
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      type (sad_pat), pointer :: kp1
      type (sad_complex), pointer :: cx
      integer*8 ks
      integer*4 ,intent(in):: iopc1
      integer*4 m,irtc,i,isp0,iopc,itfcanonicalorder
      real*8 vx1,vy,v2,vx,x
      logical*4 eval
      iopc=iopc1
      ky=k
      select case (iopc)
      case(mtfplus,mtftimes)
        if(.not. tfnumberq(k1) .and. tfnumberq(ky))then
          kx=ky
          ky=k1
        else
          kx=k1
        endif
        if(ktflistq(ky,listy))then
          if(listy%head%k == ktfoper+iopc)then
            ke=tfjoine(kx,ky,irtc)
            if(irtc /= -1)then
              return
            endif
            irtc=0
            m=listy%nl
            if(tfnumberq(listy%dbody(1)))then
              if(tfnumberq(kx))then
                ky1=tfcmplx(listy%dbody(1),kx,iopc,irtc)
                if(irtc /= 0)then
                  ke%k=ktfoper+mtfnull
                  return
                endif
                if(ktfrealq(ky1))then
                  if(iopc == mtfplus)then
                    if(ky1%k == 0 .or. ky1%k == ktfmzero)then
                      if(m == 2)then
                        ke=listy%dbody(2)
                      else
                        ke=tftake(ky,dfromr(dble(-m+1)),.true.,.false.,irtc)
                      endif
                      return
                    endif
                  else
                    if(ky1%k == ktftrue)then
                      if(m == 2)then
                        ke=listy%dbody(2)
                      else
                        ke=tftake(ky,dfromr(dble(-m+1)),.true.,.false.,irtc)
                      endif
                      return
                    elseif((ky1%k == 0 .or. ky1%k == ktfmzero) .and.
     $                     redmath%value%k /= 0)then
                      ke%k=0
                      return
                    endif
                  endif
                endif
                listy=>tfclonelist(listy)
                call tfreplist(listy,1,ky1,eval)
                ke=sad_descr(listy)
              else
                ke=tfinsertsort(listy,kx)
              endif
              return
            endif
            if(ktfrealq(kx))then
              if(iopc == mtfplus .and. (kx%k == 0 .or. kx%k == ktfmzero))then
                ke=ky
                return
              elseif(iopc == mtftimes)then
                if(kx%k == ktftrue)then
                  ke=ky
                  return
                elseif((kx%k == 0 .or. kx%k == ktfmzero) .and.redmath%value%k /= 0)then
                  ke%k=0
                  return
                endif
              endif
            endif
            ke=tfinsertsort(listy,kx)
            return
          endif
        endif
        if(ktflistq(kx,klx))then
          if(klx%head%k == ktfoper+iopc)then
            ke=tfinsertsort(klx,ky)
            return
          endif
        elseif(ktfrealq(kx))then
          if(iopc == mtfplus .and. (kx%k == 0 .or. kx%k == ktfmzero))then
            ke=ky
            return
          elseif(iopc == mtftimes)then
            if(kx%k == ktftrue)then
              ke=ky
              return
            elseif((kx%k == 0 .or. kx%k == ktfmzero) .and. redmath%value%k /= 0)then
              ke%k=0
              return
            endif
          endif
        endif
        isp=isp+3
        ktastk(isp-2)=ktfoper+iopc
        if(itfcanonicalorder(ky,kx) .ge. 0)then
          dtastk(isp-1)=kx
          dtastk(isp  )=ky
        else
          dtastk(isp-1)=ky
          dtastk(isp  )=kx
        endif
        ke=kxcompose(isp-2)
        isp=isp-3
        return
      case (mtfnot)
        if(ktflistq(ky,listy))then
          if(listy%head%k == ktfoper+mtfnot)then
            ke=listy%dbody(1)
            return
          endif
        endif
        go to 5000
      case (mtfslot,mtfslotseq)
        if(ktfrealq(ky,vy))then
          ks=int8(vy)
          if(dble(ks) /= vy) then
            go to 5000
          endif
          if(ks .le. 0.or. ks > nslots)then
            go to 5000
          endif
        elseif(ky%k /= ktfoper+mtfnull)then
          go to 5000
        else
          ks=1
        endif
        if(iopc == mtfslot)then
          ke=dlist(iaxslotnull+(ks-1)*2)
        else
          ke=dlist(iaxslotnull+(ks-1)*2+1)
        endif
        return
      case (mtfflag)
        go to 5000
      case (mtfcomp,mtfconcat,mtfand,mtfor,mtfalt,mtfmessagename)
        if(ktflistq(k1,list1))then
          if(list1%head%k == ktfoper+iopc)then
            if(ktflistq(ky,listy))then
              if(listy%head%k == ktfoper+iopc)then
                ke=tfjoin2(k1,ky,.false.,irtc)
                return
              endif
            endif
            ke=tfappend(k1,ky,.false.,0,irtc)
            return
          endif
        endif
        if(ktflistq(ky,listy))then
          if(listy%head%k == ktfoper+iopc)then
            ke=tfappend(ky,k1,.false.,1,irtc)
            return
          endif
        endif
      case (mtfset,mtfpower)
        if(ktflistq(ky,listy))then
          if(listy%head%k == ktfoper+iopc)then
            ke=tfappend(ky,k1,.false.,1,irtc)
            return
          endif
        endif
        if(iopc == mtfpower)then
          if(ktfrealq(ky))then
            if(ky%k == ktftrue)then
              ke=k1
              return
            elseif((ky%k == 0 .or. ky%k == ktfmzero).and. redmath%value%k /= 0)then
              ke%k=ktftrue
              return
            endif
          endif
          if(ktflistq(k1,list1))then
            if(list1%head%k == ktfoper+mtfpower)then
              m=list1%nl
              if(m == 1)then
                kx=ky
              else
                if(m == 2)then
                  k2=list1%dbody(2)
                else
                  k2=tftake(k1,dfromr(dble(-m+1)),.true.,.false.,irtc)
                endif
                if(ktfrealq(k2,v2) .and. ktfrealq(ky,vy))then
                  kx=dfromr(v2*vy)
                else
                  kx=tfeexpr(k2,ky,mtftimes)
                endif
              endif
              ky=list1%dbody(1)
              ke=tfeexpr(ky,kx,mtfpower)
              return
            elseif(list1%head%k == ktfoper+mtftimes
     $             .and. ktfrealq(ky))then
              vx1=1.d0
              kx2%k=ktftrue
              m=list1%nl
              isp0=isp
              isp=isp+1
              ktastk(isp)=ktfoper+mtftimes
              do i=1,m
                isp=isp+1
                dtastk(isp)=list1%dbody(i)
                if(ktfrealq(ktastk(isp)))then
                  vx1=vx1*rtastk(isp)
                  isp=isp-1
                elseif(ktflistq(ktastk(isp),listi))then
                  if(listi%head%k == ktfoper+mtfcomplex)then
                    kx=tfeexpr(dtastk(isp),ky,mtfpower)
                    if(ktfrealq(kx,vx))then
                      vx1=vx1*vx
                    else
                      kx2=tfeexpr(kx2,kx,mtftimes)
                    endif
                    isp=isp-1
                  endif
                endif
              enddo
              if(isp /= isp0+m+1)then
                if(isp > isp0+2)then
                  kx=tfcompose(isp0+1,ktfoper+mtftimes,irtc)
                elseif(isp == isp0+2)then
                  kx=dtastk(isp)
                else
                  kx%k=ktftrue
                endif
                isp=isp0
                ke=tfeexpr(kx,ky,mtfpower)
                if(ktfrealq(kx2,v2))then
                  vx1=vx1*v2
                else
                  ke=tfeexpr(kx2,ke,mtftimes)
                endif
                if(vx1 /= 1.d0)then
                  kx=tfcmplx(sad_descr(vx1),ky,mtfpower,irtc)
                  ke=tfeexpr(kx,ke,mtftimes)
                endif
                return
              else
                isp=isp0
              endif
            endif
          elseif(ktfrealq(k1))then
            if(k1%k == ktftrue .and. redmath%value%k /= 0)then
              ke%k=ktftrue
              return
            endif
          endif
        endif
      case (mtfreplace,mtfreplacerepeated)
        if(ktflistq(k1,kl1))then
          if(kl1%head%k == ktfoper+iopc)then
            ke=tfappend(k1,ky,.false.,0,irtc)
            return
          endif
        endif
      case (mtfcolon)
        if(ktfsymbolq(k1,sym))then
          call sym_symstr(sym,str)
          ke=kxpalocb(str%str,str%nch,ky,transfer(ktfref,k))
          return
        elseif(ktfpatq(k1,kp1))then
          kp1%default=dtfcopy(ky)
          ke%k=ktfpat+ktfaddrd(k1)
          return
        endif
      case (mtfrevpower)
        ke=tfeexpr(k,k1,mtfpower)
        return
      case (mtfatt)
        if(ktfnonrealq(k1))then
          if(ktfrealq(k,x))then
            ke=kxavaloc(-1,1,kle)
            kle%rbody(1)=x
            kle%head=dtfcopy(k1)
            return
          elseif((ktfsymbolq(ky) .or. ktfoperq(ky)) .and.
     $           (ktfsymbolq(k1) .or. ktflistq(k1)) .or.
     $           ktfpatq(ky))then
            go to 4900
          elseif(ktflistq(k1,kl1) .and.
     $           kl1%head%k == ktfoper+mtfatt)then
            isp=isp+1
            isp0=isp
            ktastk(isp0)=ktfoper+mtfatt
            call tfgetllstkall(kl1)
            isp=isp+1
            dtastk(isp)=ky
            ke=kxcompose(isp0)
            isp=isp0-1
            return
          else
            isp=isp+2
            dtastk(isp-1)=k1
            dtastk(isp  )=ky
            ke=kxcompose(isp-1)
            isp=isp-2
            return
          endif
        endif
      case (mtfcomplex)
        if(ktfrealq(k))then
          if(k%k == 0 .or. k%k == ktfmzero)then
            ke=k1
            return
          endif
        elseif(tfcomplexq(k,cx))then
          ky=tfeexpr(cx%dbody(1),k1,mtfplus)
          ke=tfeexpr(ky,cx%dbody(2),mtfcomplex)
          return
        elseif(tfcomplexq(k1,cx))then
          ky=tfeexpr(cx%dbody(2),k,mtfplus)
          ke=tfeexpr(cx%dbody(1),ky,mtfcomplex)
          return
        endif
      case (mtffun)
        if(k%k == ktfoper+mtfnull)then
          ke=kxpfaloc(k1)
          return
        endif
      end select
 4900 isp=isp+3
      ktastk(isp-2)=ktfoper+iopc
      dtastk(isp-1)=k1
      dtastk(isp  )=ky
      ke=kxcompose(isp-2)
      isp=isp-3
      return
 5000 if(ktfrealq(ky))then
        ke=kxavaloc(-1,1,kle)
        kle%dbody(1)=ky
      else
        ke=kxadaloc(-1,1,klx)
        call descr_rlist(ke,kle)
        klx%dbody(1)=dtfcopy(ky)
      endif
      kle%head%k=ktfoper+iopc
      return
      end function

      function tfearray(k1,k,iopc1,irtc) result(kx)
      use complex
      use dot,only:tfdot
      implicit none
      type (sad_descriptor) ,intent(in):: k,k1
      type (sad_descriptor) kx,ky
      type (sad_dlist), pointer :: kl,kl1
      integer*4 ,intent(out):: irtc
      integer*8 ka1
      integer*4 ne,ne1,i,iopc1,isp0
      logical*4 list1,list
c     begin initialize for preventing compiler warning
c     end   initialize for preventing compiler warning
c$$$      if(tfmatrixqd(k1,kl1))then
c$$$        if(tfmatrixqd(k,kl2))then
c$$$          call tfematrix(kl1,kl2,kx,iopc1,irtc)
c$$$        else
c$$$c     call tfematrix1(kl1,k,kx,iopc1,irtc)
c$$$        endif
c$$$        if(irtc /= 0)then
c$$$          go to 101
c$$$        endif
c$$$        return
c$$$      elseif(tfmatrixqd(k,kl2))then
c$$$c     call tfmatrix2(k1,kl2,kx,iopc1,irtc)
c$$$        if(irtc /= 0)then
c$$$          exit
c$$$        endif
c$$$        return
c$$$      endif
      irtc=0
      do
        select case(iopc1)
        case (mtfplus:mtfless,mtfand:mtfnot,mtfcomplex)
          kx=tfecmplxl(k1,k,iopc1)
          return
        case (mtfsame:mtfunsame)
          exit
        end select
        if(ktflistq(k1,kl1))then
          if(tfcomplexq(k1))then
            ne1=0
            list1=.false.
          else
            if(tfexprq(k1))then
              exit
            endif
            ne1=kl1%nl
            list1=.true.
          endif
        else
          ne1=0
          list1=.false.
        endif
        if(ktflistq(k,kl))then
          if(tfcomplexq(k))then
            list=.false.
            ne=0
          else
            if(tfexprq(k))then
              exit
            endif
            if(iopc1 == mtfdot)then
              kx=tfdot(k1,k,irtc)
              return
            endif
            ne=kl%nl
            list=.true.
            if(list1)then
              if(iopc1 == mtfequal)then
                if(ne /= ne1)then
                  kx%k=0
                else
                  kx=tfecmplxl(k1,k,iopc1)
                endif
                return
              elseif(iopc1 == mtfunequal)then
                if(ne /= ne1)then
                  kx%k=ktftrue
                else
                  kx=tfecmplxl(k1,k,iopc1)
                endif
                return
              elseif(ne /= ne1)then
                exit
              endif
            endif
          endif
        else
          list=.false.
          ne=0
        endif
        if(list1)then
          if(list)then
            if(iopc1 == mtfequal)then
              kx%k=ktftrue
              do i=1,ne
                ky=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc /= 0)then
                  return
                endif
                if(ky%k == 0 .or. ky%k == ktfmzero)then
                  kx%k=0
                  return
                elseif(.not. ktfrealq(ky))then
                  exit
                endif
              enddo
              return
            elseif(iopc1 == mtfunequal)then
              kx%k=0
              do i=1,ne
                ky=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc /= 0)then
                  return
                endif
                if(ky%k == ktftrue)then
                  kx%k=ktftrue
                  return
                elseif(.not. ktfrealq(ky))then
                  exit
                endif
              enddo
              return
            else
              isp0=isp
              do i=1,ne
                isp=isp+1
                dtastk(isp)=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc /= 0)then
                  kx=dxnullo
                  isp=isp0
                  return
                endif
              enddo
              kx=kxmakelist(isp0)
              isp=isp0
            endif
          else
            select case (iopc1)
            case (mtfset,mtfsetdelayed)
              do i=1,ne1
                if(ktfrefq(kl1%dbody(i),ka1))then
                  call tflocald(dlist(ka1))
                  dlist(ka1)=dtfcopy(k)
                endif
              enddo
              kx=k
              return
            case (mtfequal,mtfunequal)
              exit
            end select
            isp0=isp
            do i=1,ne
              isp=isp+1
              dtastk(isp)=tfcmplx(kl1%dbody(i),k,iopc1,irtc)
              if(irtc /= 0)then
                kx=dxnullo
                isp=isp0
                return
              endif
            enddo
            kx=kxmakelist(isp0)
            isp=isp0
          endif
        else
          if(iopc1 == mtfequal .or. iopc1 == mtfunequal)then
            exit
          endif
          isp0=isp
          do i=1,ne
            isp=isp+1
            dtastk(isp)=tfcmplx(k1,kl%dbody(i),iopc1,irtc)
            if(irtc /= 0)then
              kx=dxnullo
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
        endif
        return
      enddo
      kx=tfeexpr(k1,k,iopc1)
      return
      end

      function tfset1(k10,k20,mopc,irtc) result(kx)
      use tfstk
      use eexpr,only:tfsetpart,tftagset
      use tfmessage,only:tfnewset
      use mackw
      implicit none
      type (sad_descriptor) k1,k2,kx,ks,ka,tfearray
      type (sad_descriptor) ,intent(in):: k10,k20
      type (sad_dlist),pointer :: list,kls1,kla,kls
      type (sad_symbol), pointer ::sym
      type (sad_symdef),pointer :: symd
      type (sad_namtbl), pointer :: loc
      integer*8 ka1,kaa,kas,kar
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mopc
      integer*4 itfmessage,itfmessageexp,itfmessagestr,i
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
          if(symd%sym%override == -2)then
            if(ktflistq(symd%value,kls1) .and.
     $           kls1%head%k == ktfoper+mtflist)then
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
 12     if(symd%sym%override == -2)then
          if(ktfprotectedqo(symd%sym) .and.
     $         symd%sym%gen /= -3)then
            irtc=itfmessage(999,'General::protect','""')
            return
          endif
          call tfdset(k2,symd%downval,kx,k1)
        else
          go to 9900
        endif
      elseif(ktfsymbolq(k1,sym))then
        if(sym%override == 0)then
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
        if(k2%k == ktfref)then
          ks=dtfcopy1(sad_descr(symd%sym))
          kx%k=ktfoper+mtfnull
        else
          ks=dtfcopy(k2)
          kx=ks
        endif
        if(symd%value%k == kfromd(sad_descr(symd%sym)) .and. tfnewset(.false.)
     $       .and. symd%sym%gen <= 0)then
          call loc_namtbl(symd%sym%loc,loc)
          if(loc%kind ==0)then
            irtc=itfmessagestr(9,'General::newset',loc%str%str(1:loc%str%nch))
            call tferrorhandle(dlist(loc%cont),irtc)
            irtc=0
          endif
        endif
        call tflocald(symd%value)
        symd%value=ks
      elseif(ktfrefq(k1,ka1))then
c        if(ka1 > 0 .and. ktfrealq(k2))then
        if(ka1 > 0)then
          call tflocald(dlist(ka1))
          dlist(ka1)=dtfcopy(k2)
        else
          irtc=itfmessage(999,'General::invset','"Illegal Memory Location"')
          return
        endif
      else
        go to 9900
      endif
      return
 9900 irtc=itfmessageexp(999,'General::invset',k10)
      return
      end

      module context
      use tfstk

      contains
      logical*4 function tfcontextqk(k)
      implicit none
      type (sad_descriptor) k
      type (sad_symbol), pointer :: sym
      tfcontextqk=ktfsymbolq(k,sym) .and. sym%gen == -3
      return
      end

      subroutine tfassigncont(kp,name)
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
      if(klist(kp) == 0)then
        klist(kp)=ktfsymbol+ktfcopy1(ka)
      endif
      return
      end

      subroutine tfsetcontext(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) kx
      type (sad_symdef), pointer :: symd
      integer*4 isp1,irtc,itfmessage
      if(isp1+1 /= isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      kx=dtastk(isp)
      if(.not. tfcontextqk(kx))then
        irtc=itfmessage(9,'General::wrongtype','"Context (Symbol ending `)"')
        return
      endif
      call loc_sad(ktfaddr(kx),symd)
      itfcontext=ktfaddrd(symd%value)
      irtc=0
      return
      end

      subroutine tfsetcontextpath(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ki
      type (sad_dlist), pointer :: kl
      type (sad_symdef), pointer :: symd
      integer*8 ka1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,m,i,isp0
      if(isp1+1 /= isp)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      kx=dtastk(isp)
      if(.not. tflistq(kx,kl) .or. ktfreallistq(kl))then
        go to 9000
      endif
      m=kl%nl
      isp0=isp
      do i=1,m
        ki=kl%dbody(i)
        if(.not. tfcontextqk(ki))then
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
      isp=isp0
      call tfree(itfcontextpath)
      itfcontextpath=ka1
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"List of Contexts"')
      return
      end

      end module context

      module efun
      use tfstk

      contains
      recursive function tfefunref(isp1,upvalue,irtc) result(kx)
      use tfmem
      use tfshare
      use tfcsi,only:cssetlfno,icslfno,icslfnm
      use findr,only:tffindroot,tffit
      use tfcx,only:tfatt,tfsolvemember
      use mathfun
      use eexpr
      use funs
      use eeval
      use eval1
      use context,only:tfsetcontext,tfsetcontextpath
      use modul,only:tfmodule,tfwith
      use repl, only:tfreplace1,tfreplacerepeated1
      use dset,only:tfupset
      use part,only:tfpartition
      use table,only:tftable,tfrange
      use take,only:tfdifference,tfprotect,tfreverse,tfrotateright1,tftake,tfrest
      use attrib,only:tfattributes,tfsetattributes,tfreleasehold
      use dot
      use sameq
      use convstr
      use readbuf
      use complex
      use gammaf
      implicit none
      type (sad_descriptor) kx,k1,k,kh,
     $     tfefun1,tfeval1,tfeintf,tfpuref,
     $     tfeintf2,tfget,tfmap,tfminmax,
     $     tfgetcommandline,tfreplacepart,tfpart,tfwrite,
     $     tftemporaryname,tfmapfile,tfunmapfile,
     $     tfgaussiancoulomb,tfgaussiancoulombu
      type (sad_dlist), pointer :: kl,kl1,klx,klh
      type (sad_symbol), pointer :: sym1
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      integer*8 ka1,ka,kax,kop
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,id,narg,nc,isp2,itfdepth,
     $     itfgetrecl,ltr0,iaf,itfopenwrite,itfmessage,
     $     itfopenappend,itfopcode,isp0,nsize
      real*8 vx,v1,v,f
      complex*16 c1
      character*3 opcx
      character*8 char8
      logical*4 ,intent(in):: upvalue
      logical*4 euv,rep
      type (sad_descriptor), save :: ilog2
      data ilog2%k /0/

c     DOUBLE specific math intrinsic function
c     from Fortran77    77   77   77
      intrinsic dabs,dsqrt,dexp,dlog
c     from Fortran77   77   77    77    77    77     77
      intrinsic dsin,dcos,dtan,dasin,dacos,datan
c     from Fortran 77    77    77
      intrinsic dsinh,dcosh,dtanh
c     from Fortran 08    08  <-- somehow cannot pass as an argument @ gfortran 5
c      intrinsic derf,derfc

c     DOUBLE COMPLEX specific math intrinsic function
c     from Fortran EX    EX     EX    EX    EX
      intrinsic cdsin,cdcos,cdsqrt,cdexp,cdlog

c     DOUBLE specific `t-prefix' math function proxy
c     for Fortran2008 generic function
      real*8   tasinh,tacosh,tatanh
      external tasinh,tacosh,tatanh

c     DOUBLE COMPLEX specific `t-prefix' math function proxy
c     for vendor extended math intrinsic function
      complex*16 tctan
      external   tctan

c     DOUBLE COMPLEX specific math function implemented by SAD
c      real*8 , external :: dlgama

c      real*8 , external:: aloggamma1,factorial,gammaq,gammap,inverseerf,
c     $     productlog,gamma0,ferf,ferfc
c      complex*16, external:: cloggamma1,cfactorial,cerfc,cerf,
c     $     cproductlog
      if(upvalue)then
        LOOP_I: do i=isp1+1,isp
          k1=dtastk(i)
 12       if(ktflistq(k1,kl1))then
            klh=>kl1
            k1=kl1%head
            do while(ktflistq(k1,kl1))
              k1=kl1%head
            enddo
            if(k1%k == ktfoper+mtfatt .or.
     $           k1%k == ktfoper+mtfslot)then
              k1=tfsolvemember(klh,rep,irtc)
              if(irtc == 0)then
                dtastk(i)=k1
                go to 12
              elseif(irtc > 0)then
                return
              endif
            endif
          endif
          if(ktfsymbolq(k1,sym1))then
            if(sym1%override == 0)then
              sym1=>tfsydef(sym1)
            endif
            call sym_symdef(sym1,symd)
            if(symd%upval /= 0)then
              call tfdeval(isp1,dfromk(sad_loc(symd%sym%loc)),kx,0,.false.,euv,irtc)
c              call tfdeval(isp1,sad_loc(symd%sym%loc),kx,0,.false.,euv,irtc)
              if(euv)then
                return
              endif
            endif
          endif
        enddo LOOP_I
      endif
      k1=dtastk(isp1)
      narg=isp-isp1
      if(ktfoperq(k1,ka1))then
        if(ka1 .le. mtfend)then
          go to 6000
        endif
        if(narg == 0)then
          narg=1
          isp=isp+1
          ktastk(isp)=ktfoper+mtfnull
        endif
        k=dtastk(isp)
        id=iget_fun_id(ka1)
c        write(*,*)'tfefunref-1 ',id
        irtc=-1
        go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 
     $         210, 220, 230, 240, 240, 260, 270, 280, 290, 300,
     $         310, 320, 330, 340, 350, 360, 370, 380, 380, 400,
     $         410, 420, 430, 440, 440, 460, 470, 480, 490, 500,
     $         510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
c              Sin  Cos  Tan  Sinh Cosh Tanh Exp  Log  Atan Det
c              Sqrt Flor Ceil Min  Max  Mod  StrL Arg  Sign Leng
c              Dims RplP ASin ACos ASh  ACh  ATh  Tabl Do   Attr
c              Peek Abs  Revs Modl Blck StrR SwiC Fltn If   Take
c              Selc Whil Join Apnd Ppnd Clr  Prot Unpr Drop MpAt
c
     $         610, 620, 630, 640, 650, 660, 670, 680, 690, 700,
     $         700, 720, 730, 730, 730, 760, 770, 780, 790, 800,
     $         810, 820, 830, 840, 850, 860, 870, 880, 890, 900,
     $         910, 920, 930, 940, 950, 960, 970, 980, 990,1000,
     $        1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
c              Innr Trps SglV DiaM LinS IdtN Eigs Oper Posi Sum
c              Prod Rang Re   Im   Conj ToSt Dpth Levl Writ Get
c              OpnW OpnA Clos Flsh Prnt WStr Retn Head RLiQ Pttn
c              Thrw Ctch Thrd SetA MpId FrCh ToCh CmpQ Tr   SvSM
c              Stch Sort Uni1 Ordr MChk Scan Iden TimU NumQ VecQ
c
     $        1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,
     $        1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,
     $        1310,1320,1330,1340,1350,1360,1370,1380,90,90,
     $        1290,1420,1430,1440,1450,1460,1470,1480,1490,1500,
     $        1510,1520,1530,1540,1540,1560,1570,1580,1590,1600,
c              AtmQ Outr MatQ TrcP Defi READ Ints Cmpl Roun IErf
c              StIO PolG ToIS ReaS OpnR ToEx StrM StrP ToUp Brek
c              Cont Goto Four IFou Chek Whic MapF UnmF GetU GetG
c              ToLo Unev Case DelC Vect LogG Nams GbCl LgG1 Fact
c              With WhiC Ovrr AppT PreT FndR GamR GmRP Erf  Erfc
c
     $        1610,1620, 760,1640,1650,1660,1670,1680,1690,1700,
     $        1710,1720,1730,1740,1750,1760,1770,1770,1770,1770,
     $        1810,1820,1830,1840,1850,1860,1870,1870,1870,1900,
     $        1910,1920,1930,1940,1950,1960,1970,1980,1980,1980,
     $        2010,2020,2030,2040,2050,2060,2070,2080,2080,2100,
c              Fit  Symb SyNm Extr Read Skip TmpN Exit StrF Rstr
c              MM   Shrt $SOT Dir  SDir Wait BesJ BesY BesI BesK
c              BasF StrT S2S  Even OddQ DatS Inst Delt FlAt Repl
c              SetE Spl$ FInd SCnt SCnP ToCt Ctxt BAnd BOr  BXor
c              RepM MSca StdF Abrt ChkA RelH NaNQ MapT ScaT Last
c
     $        2100,2100,2100,2140,2150,2160,2170,2180,2190,2200,
     $        2210,2220,2230,2240,2250,2260,2270,2280,2290,2300,
     $        2310,2320,2330,2340,2350,2360,2370,2380,2390,2400,
     $        2410,2420,2430,2440,2450,2460,2470,2480,2490,2500,
     $        2510,2520,2530,2540,2550,2560,2570,2580,2590
c              Frst Scnd Thrd ObjS PrdL GauC ClMO MAll Dupl GCLn
c              Seek DigQ LetQ ReaQ NSmQ OpSh RdSh WrSh ShSz FBuQ
c              GaCU GaCF Rest RRt1 Diff Gam0 XSIn Poch Hg21 H21R 
c              Hg11 H11R Hg01 H01R Zeta PGmM DZet HZet DHZt PGmH
c              Gamm HgU  HgPQ HPQR PLog Beta HLcP LchP BerB
     $       ),id
c
        if(id > 0)then
          if(id .le. 2000)then
            kx=tfefun1(isp1,id,.true.,irtc)
            go to 6900
          elseif(id .le. 3000) then
            call tfefun2(isp1,id,k,kx,irtc)
            go to 6900
          elseif(id .le. 4000) then
            call tfefun3ep(isp1,id,kx,irtc)
            go to 6900
          elseif(id .le. 5000) then
            call tfefunctbl8(isp1,id,kx,irtc)
            go to 6900
          endif
        endif
        go to 100
 90     irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 100    write(*,*)'Function implementation error: ',id
        irtc=0
        go to 7000
 110    if(narg == 1)then
          kx=tfeintf(dsin,cdsin,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 120    if(narg == 1)then
          kx=tfeintf(dcos,cdcos,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 130    if(narg == 1)then
          kx=tfeintf(dtan,tctan,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 140    if(narg == 1)then
          kx=tfeintf(dsinh,tcsinh,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 150    if(narg == 1)then
          kx=tfeintf(dcosh,tccosh,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 160    if(narg == 1)then
          kx=tfeintf(dtanh,tctanh,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 170    if(narg == 1)then
          kx=tfeintf(dexp,cdexp,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 180    if(narg == 1)then
          kx=tfeintf(dlog,cdlog,k,.true.,0.d0,dinfinity,irtc)
        elseif(narg == 2)then
          if(ilog2%k == 0)then
            ilog2=kxsymbolf('Log2$',5,.true.)
          endif
          dtastk(isp1)=ilog2
          call tfefun(isp1,kx,.true.,.false.,irtc)
          isp=isp1+2
          dtastk(isp1)=k1
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
        endif
        go to 6900
 190    if(narg == 1)then
          kx=tfeintf(datan,tcatan,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfeintf2(tatan2,tcatan2,dtastk(isp-1),k,irtc)
        else
          go to 6812
        endif
        go to 6900
 200    call tfdet(isp1,kx,irtc)
        go to 6900
 210    if(narg == 1)then
          kx=tfeintf(dsqrt,cdsqrt,k,.true.,0.d0,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 220    if(narg == 1)then
          kx=tfeintf(tfloor,tcfloor,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 230    if(narg == 1)then
          kx=tfeintf(tceiling,tcceiling,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 240    kx=tfminmax(isp1,id-13,irtc)
        go to 6900
 260    call tfmod(isp1,kx,0,irtc)
        go to 6900
 270    if(narg == 1)then
          if(ktfstringq(k,str))then
            kx=dfromr(dble(str%nch))
            go to 8000
          else
            irtc=itfmessage(9,'General::wrongtype','"Character-string"')
          endif
        else
          go to 6811
        endif
        go to 6900
 280    if(narg == 1)then
          kx=tfeintf(tfarg,tfcarg,k,.false.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 290    if(narg == 1)then
          kx=tfeintf(tfsign,tfcsign,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 300    if(narg /= 1)then
          go to 6811
        endif
        if(ktflistq(k,kl))then
          kx%x(1)=dble(kl%nl)
        else
          kx%x(1)=0.d0
        endif
        go to 8000
 310    call tfdimensions(isp1,kx,irtc)
        go to 6900
 320    kx=tfreplacepart(isp1,0,irtc)
        go to 6900
 330    if(narg == 1)then
          kx=tfeintf(dasin,tcasin,k,.true.,-1.d0,1.d0,irtc)
        else
          go to 6811
        endif
        go to 6900
 340    if(narg == 1)then
          kx=tfeintf(dacos,tcacos,k,.true.,-1.d0,1.d0,irtc)
        else
          go to 6811
        endif
        go to 6900
 350    if(narg == 1)then
          kx=tfeintf(tasinh,tcasinh,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 360    if(narg == 1)then
          kx=tfeintf(tacosh,tcacosh,k,.true.,1.d0,dinfinity,
     $         irtc)
        else
          go to 6811
        endif
        go to 6900
 370    if(narg == 1)then
          kx=tfeintf(tatanh,tcatanh,k,.true.,-1.d0,1.d0,
     $         irtc)
        else
          go to 6811
        endif
        go to 6900
 380    isp2=isp
        kx=tftable(isp1,isp1+2,isp2,29-id,irtc)
        go to 6900
 400    call tfattributes(isp1,kx,irtc)
        go to 6900
 410    if(narg == 2)then
          if(ktfnonrealq(dtastk(isp),f))then
            go to 6812
          endif
        elseif(narg /= 1)then
          go to 6811
        else
          f=0.d0
        endif
        if(ktfnonrealq(k))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
        else
          ka=int8(rtastk(isp1+1))
          if(f == 0.d0 .and. .not. tfchecklastp(ka))then
            irtc=itfmessage(9,'General::wrongnum',
     $           '"within allocated block"')
          else
            kx=kxadaloc(-1,4,klx)
            klx%rbody(1)=rlist(ka)
            klx%rbody(2)=dble(klist(ka))
            klx%rbody(3)=dble(ilist(1,ka))
            klx%rbody(4)=dble(ilist(2,ka))
            klx%dbody(5)=kxsalocb(0,transfer(klist(ka),char8),8)
            irtc=0
          endif
        endif
        go to 6900
 420    if(narg == 1)then
          kx=tfeintf(dabs,ccdabs,k,
     $         .false.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 430    call tfreverse(isp1,kx,irtc)
        go to 6900
 440    kx=tfmodule(isp1,id == nfunmodule,.true.,irtc)
        go to 6900
 460    call tfstringreplace(isp1,kx,irtc)
        go to 6900
 470    call tfswitchcases(isp1,kx,0,irtc)
        go to 6900
 480    call tfflatten(isp1,kx,irtc)
        go to 6900
 490    call tfiff(isp1,kx,irtc)
        go to 6900
 500    if(narg /= 2)then
          go to 6812
        else
          kx=tftake(dtastk(isp-1),k,.true.,.true.,irtc)
        endif
        go to 6900
 510    call tfselect(isp1,kx,irtc)
        go to 6900
 520    call tfwhile(isp1,kx,irtc)
        go to 6900
 530    kx=tfjoin(isp1,.true.,irtc)
        go to 6900
 540    if(narg == 2)then
          kx=tfappend(dtastk(isp1+1),k,.true.,0,irtc)
        elseif(narg /= 1)then
          go to 6812
        endif
        go to 6900
 550    if(narg == 2)then
          kx=tfappend(dtastk(isp1+1),k,.true.,1,irtc)
        elseif(narg /= 1)then
          go to 6812
        endif
        go to 6900
 560    call tfclear(isp1,kx,irtc)
        go to 6900
 570    call tfprotect(isp1,kx,.true.,irtc)
        go to 6900
 580    call tfprotect(isp1,kx,.false.,irtc)
        go to 6900
 590    if(narg /= 2)then
          go to 6812
        else
          kx=tftake(dtastk(isp-1),k,.false.,.true.,irtc)
        endif
        go to 6900
 600    kx=tfreplacepart(isp1,1,irtc)
        go to 6900
 610    if(narg /= 4)then
          irtc=itfmessage(9,'General::wrongnum','"4"')
        else
          call tfinner(dtastk(isp-2),dtastk(isp-1),
     $         kx,k,dtastk(isp-3),irtc)
        endif
        go to 6900
 620    if(narg /= 1)then
          go to 6811
        else
          call tftranspose(k,kx,irtc)
        endif
        go to 6900
 630    call tfsingularvalues(isp1,kx,irtc)
        go to 6900
 640    if(narg /= 1)then
          go to 6811
        else
          call tfdiagonalmatrix(k,kx,irtc)
        endif
        go to 6900
 650    call tflinearsolve(isp1,kx,irtc)
        go to 6900
 660    if(narg /= 1)then
          go to 6811
        else
          call tfidentitymatrix(k,kx,irtc)
        endif
        go to 6900
 670    if(narg /= 1)then
          go to 6811
        else
          call tfeigensystem(k,kx,irtc)
        endif
        go to 6900
 680    if(narg /= 1)then
          go to 6811
        elseif(.not. ktfstringq(k))then
          irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        else
          ka=ktfaddr(k)
          nc=ilist(1,ka)
          if(nc == 0 .or. nc > 3)then
            irtc=itfmessage(9,'General::invop',' ')
            go to 6900
          endif
          opcx=' '
          call tmovb(ilist(1,ka+1),opcx,nc)
          kax=itfopcode(opcx)
          if(kax .lt. 0)then
            irtc=itfmessage(9,'General::invop',' ')
          else
            kx%k=ktfoper+kax
            irtc=0
          endif
        endif
        go to 6900
 690    call tfposition(isp1,kx,0,irtc)
        go to 6900
 700    isp2=isp
        kx=tftable(isp1,isp1+2,isp2,id-58,irtc)
        go to 6900
 720    kx=tfrange(isp1,irtc)
        go to 6900
 730    if(narg == 1)then
          kx=tfcmplxf(k,id-62,int(ka1),irtc)
        else
          go to 6811
        endif
        go to 6900
 760    call tftostring(isp1,kx,id == 153,irtc)
        go to 6900
 770    if(narg == 1)then
          kx=dfromr(dble(itfdepth(k)))
          go to 8000
        else
          go to 6811
        endif
 780    if(narg == 2)then
          call tflevel(ktastk(isp1+1),k,kx,irtc)
        else
          go to 6812
        endif
        go to 6900
 790    kx=tfwrite(isp1,irtc)
        go to 6900
 800    if(narg /= 1)then
          go to 6811
        else
          kx=tfget(k,irtc)
        endif
        go to 6900
 810    if(narg /= 1)then
          go to 6811
        else
          vx=itfopenwrite(k,irtc)
        endif
        go to 829
 820    if(narg /= 1)then
          go to 6811
        else
          vx=itfopenappend(k,irtc)
        endif
 829    if(irtc /= 0)then
          if(vx == -2.d0)then
            call tfaddmessage(' ',0,icslfnm())
            irtc=0
            kx%k=kxfailed
            return
          endif
          go to 6900
        endif
        kx=dfromr(vx)
        go to 8000
 830    if(narg /= 1)then
          go to 6811
        else
          call tfclosef(k,irtc)
        endif
        go to 8200
 840    call tfflush(isp1,kx,irtc)
        go to 6900
 850    call tfprintf(isp1,kx,irtc)
        go to 6900
 860    call tfwritestring(isp1,kx,irtc)
        go to 6900
 870    if(narg == 1)then
          call tfthrow(irtcret,k,irtc)
          return
        else
          go to 6811
        endif
        go to 6900
 880    call tfhead(k,kx)
        go to 8100
 890    if(narg == 1)then
          kx%k=merge(ktftrue,i00,tfreallistq(k))
          irtc=0
        else
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        go to 6900
 900    call tfpartition(isp1,kx,irtc)
        go to 6900
 910    if(narg == 1)then
          call tfthrow(irtcthrow,k,irtc)
          return
        else
          go to 6811
        endif
        go to 6900
 920    call tfcatch(isp1,kx,irtc)
        go to 6900
 930    call tfthread(isp1,kx,0,irtc)
        go to 6900
 940    call tfsetattributes(isp1,kx,irtc)
        go to 6900
 950    kx=tfmap(isp1,4,1,irtc)
        go to 6900
 960    call tffromcharactercode(isp1,kx,irtc)
        go to 6900
 970    call tftocharactercode(isp1,kx,irtc)
        go to 6900
 980    call tfcomplexlistqkf(isp1,kx,irtc)
        go to 6900
 990    call tftr(isp1,kx,irtc)
        go to 6900
c 990    irtc=itfmessage(999,'General::unregister',' ')
c        go to 6900
 1000   call tfsavesharedmap()
        irtc=0
        go to 6900
c 1000   irtc=itfmessage(999,'General::unregister',' ')
c        go to 6900
 1010   call tfswitch(isp1,kx,irtc)
        go to 6900
 1020   call tfsort(isp1,kx,0,irtc)
        go to 6900
 1030   call tfsort(isp1,kx,1,irtc)
        go to 6900
 1040   call tforder(isp1,kx,irtc)
        go to 6900
 1050   call tfmemcheck(isp1,kx,irtc)
        go to 6900
 1060   kx=tfmap(isp1,1,1,irtc)
        go to 6900
 1070   if(narg == 1)then
          kx=k
          irtc=0
          return
        else
          go to 6811
        endif
 1080   if(narg == 1)then
          call cputime(vx,irtc)
          kx=dfromr(vx*1.d-6)
          go to 8000
        else
          go to 6811
        endif
        go to 6900
 1090   if(narg == 1)then
          kx%k=merge(ktftrue,ktffalse,tfnumberq(k))
          go to 8000
        else
          go to 6811
        endif
        go to 6900
 1100   call tfvectorqf(isp1,kx,irtc)
        go to 6900
 1110   if(narg == 1)then
          kx%k=merge(merge(ktftrue,ktffalse,
     $         kl%head%k == ktfoper+mtfcomplex),
     $         ktftrue,ktflistq(k,kl))
          go to 8000
        else
          go to 6811
        endif
        go to 6900
 1120   call tfouter(isp1,kx,irtc)
        go to 6900
 1130   call tfmatchqf(isp1,kx,irtc)
        go to 6900
 1140   if(narg == 1)then
          if(ktflistq(k,kl))then
            if(kl%head%k /= ktfoper+mtfcomp)then
              call tfprint1(k,
     $             6,-itfgetrecl(),4,.true.,.true.,irtc)
            endif
          else
            call tfprint1(k,6,-itfgetrecl(),4,.true.,.true.,irtc)
          endif
          ltr0=ltrace
          ltrace=6
          kx=tfeevalref(k,irtc)
          ltrace=ltr0
        else
          go to 6811
        endif
        go to 6900
 1150   call tfdefinition(isp1,kx,irtc)
        go to 6900
 1160   if(narg /= 1)then
          go to 6811
        else
          call nfread(k,kx,irtc)
        endif
        go to 6900
 1170   call tfintersection(isp1,kx,0,irtc)
        go to 6900
 1180   call tfintersection(isp1,kx,1,irtc)
        go to 6900
 1190   if(narg == 1)then
          kx=tfeintf(tround,tcround,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfeintf2(tround2,tcround2,dtastk(isp-1),k,irtc)
        else
          go to 6811
        endif
        go to 6900
 1200   if(narg == 1)then
          kx=tfeintf(inverseerf,cinverseerf,k,.true.,
     $         -1.d0,1.d0,irtc)
        else
          go to 6811
        endif
        go to 6900
 1210   call tfsetioch(isp1,kx,irtc)
        go to 6900
 1220   if(narg == 1)then
          kx=tfeintf(polygamma,cpolygamma,k,.true.,
     $         -dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfbessel(isp1,4,irtc)
        else
          go to 6812
        endif
        go to 6900
 1230   call tftoinputstring(isp1,kx,irtc)
        go to 6900
 1240   call tfreadstring(isp1,kx,.false.,.false.,irtc)
        go to 6900
 1250   call tfopenread(isp1,kx,irtc)
        go to 6900
 1260   call tftoexpression(isp1,kx,irtc)
        go to 6900
 1270   call tfstringmatchq(isp1,kx,irtc)
        go to 6900
 1280   call tfstringposition(isp1,kx,irtc)
        go to 6900
 1290   call tftouppercase(isp1,kx,id-119,irtc)
        go to 6900
 1300   continue
 1310   call tfbreak(id-123,narg,kx,irtc)
        go to 6900
 1320   if(narg == 1)then
          call tfthrow(irtcgoto,k,irtc)
          return
        else
          go to 6811
        endif
        go to 6900
 1330   continue
 1340   if(narg == 1)then
          call tffourier(id == 124,k,kx,irtc)
        else
          go to 6811
        endif
        go to 6900
 1350   kx=tfcheck(isp1,irtc)
        go to 6900
 1360   call tfwhich(isp1,kx,irtc)
        go to 6900
 1370   kx=tfmapfile(isp1,irtc)
        go to 6900
 1380   kx=tfunmapfile(isp1,irtc)
        go to 6900
 1420   kx=tfsequence(isp1,isp)
        irtc=0
        return
 1430   call tfcases(isp1,kx,irtc)
        go to 6900
 1440   call tfposition(isp1,kx,2,irtc)
        go to 6900
 1450   call tfvectorize(isp1,kx,irtc)
        go to 6900
 1460   if(narg == 1)then
          kx=tfeintf(aloggamma,cloggamma,k,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1470   call tfnames(isp1,kx,irtc)
        go to 6900
 1480   call tfgarbagecollect(isp1,kx,irtc)
        go to 6900
 1490   if(narg == 1)then
          kx=tfeintf(aloggamma1,cloggamma1,k,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1500   if(narg == 1)then
          kx=tfeintf(factorial,cfactorial,k,
     $         .true.,-1.d0+5.56d-17,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1510   call tfwith(isp1,kx,.true.,irtc)
        go to 6900
 1520   call tfswitchcases(isp1,kx,1,irtc)
        go to 6900
 1530   kx=tfoverride(isp1,irtc)
        go to 6900
 1540   call tfappendto(isp1,kx,id-143,irtc)
        go to 6900
 1560   call tffindroot(isp1,kx,irtc)
        go to 6900
 1570   kx=tfbessel(isp1,9,irtc)
        go to 6900
 1580   kx=tfbessel(isp1,10,irtc)
        go to 6900
 1590   if(narg == 1)then
          kx=tfeintf(ferf,cerf,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1600   if(narg == 1)then
          kx=tfeintf(ferfc,cerfc,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1610   call tffit(isp1,kx,irtc)
        go to 6900
 1620   call tfsymbol(isp1,kx,irtc)
        go to 6900
 1640   call tfextract(isp1,kx,irtc)
        go to 6900
 1650   call tfread(isp1,kx,irtc)
        go to 6900
 1660   call tfskip(isp1,kx,irtc)
        go to 6900
 1670   kx=tftemporaryname(isp1,irtc)
        go to 6900
 1680   call tfexit(isp1,kx,irtc)
        go to 6900
 1690   call tfstringfill(isp1,kx,irtc)
        go to 6900
 1700   call tfrestrict(isp1,kx,irtc)
        go to 6900
 1710   kx=tfminmax(isp1,0,irtc)
        go to 6900
 1720   call tfshort(isp1,kx,irtc)
        go to 6900
 1730   call setompnumthreads(isp1,kx,irtc)
        go to 6900
 1740   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1750   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1760   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1770   kx=tfbessel(isp1,id-167,irtc)
        go to 6900
 1810   call tfbaseform(isp1,kx,irtc)
        go to 6900
 1820   call tfstringtrim(isp1,kx,irtc)
        go to 6900
 1830   call tfstringtostream(isp1,kx,irtc)
        go to 6900
 1840   if(narg == 1)then
          kx=tfeintf(tfevenq,tfcevenq,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1850   if(narg == 1)then
          kx=tfeintf(tfoddq,tfcoddq,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 1860   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1870   kx=tfreplacepart(isp1,id-175,irtc)
        go to 6900
 1900   if(narg /= 2)then
          go to 6812
        endif
        call tfreplace(dtastk(isp1+1),dtastk(isp),
     $       kx,.false.,.true.,.false.,irtc)
        go to 6900
 1910   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1920   call tfspline(isp1,kx,irtc)
        go to 6900
 1930   call tffindindex(isp1,kx,irtc)
        go to 6900
 1940   call tfsetcontext(isp1,kx,irtc)
        go to 6900
 1950   call tfsetcontextpath(isp1,kx,irtc)
        go to 6900
 1960   call tftocontext(isp1,kx,irtc)
        go to 6900
 1970   call tfcontext(isp1,kx,irtc)
        go to 6900
 1980   call tfmod(isp1,kx,id-187,irtc)
        go to 6900
 2010   call tfreplacemember(isp1,kx,irtc)
        go to 6900
 2020   call tfmemberscan(isp1,kx,irtc)
        go to 6900
 2030   call tfstandardform(isp1,kx,irtc)
        go to 6900
 2040   irtc=-7
        if(narg == 1)then
          if(ktfrealq(k,v))then
            irtc=int(min(max(-3.d0,v),-1.d0)-3.d0)
          endif
        endif
        go to 6900
 2050   go to 6900
 2060   call tfreleasehold(isp1,kx,irtc)
        go to 6900
 2070   call tfnanqk(isp1,kx,irtc)
        go to 6900
 2080   call tfthread(isp1,kx,id-197,irtc)
        go to 6900
 2100   call tffirst(isp1,kx,id-201,irtc)
        go to 6900
 2140   call tfobjectsymbol(isp1,kx,irtc)
        go to 6900
 2150   if(narg == 1)then
          kx=tfeintf(productlog,cproductlog,k,
     $         .true.,-exp(-1.d0),dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 2160   kx=tfgaussiancoulomb(isp1,irtc)
        go to 6900
 2170   call tfclearmemberobject(isp1,kx,irtc)
        go to 6900
 2180   call tfmalloc(isp1,kx,irtc)
        go to 6900
 2190   if(narg == 1)then
          kx=k
          if(ktflistq(k))then
            kx=kxcopylist(k)
          endif
          irtc=0
          return
        else
          go to 6811
        endif
 2200   kx=tfgetcommandline(isp1,irtc)
        go to 6900
 2210   call tfseek(isp1,kx,irtc)
        go to 6900
 2220   call tfdigitQ(isp1,kx,irtc)
        go to 6900
 2230   call tfletterQ(isp1,kx,irtc)
        go to 6900
 2240   call tfrealqk(isp1,kx,irtc)
        go to 6900
 2250   if(narg /= 4)then
          go to 6814
        endif
        if(ktfnonrealq(dtastk(isp-1)) .or.
     $       ktfnonrealq(dtastk(isp)))then
          irtc=itfmessage(9,'General::wrongarg',
     $         '"$NearlySameQ[a,b,relthre,absthre]"')
          go to 7000
        endif
        irtc=0
        kx%k=merge(ktftrue,ktffalse,
     $       tfnearlysameqf(dtastk(isp1+1),dtastk(isp1+2),rtastk(isp1+3),rtastk(isp)))
        go to 6900
 2260   kx=tfopenshared(isp1,irtc)
        go to 6900
 2270   kx=tfreadshared(isp1,irtc)
        go to 6900
 2280   kx=tfwriteshared(isp1,irtc)
        go to 6900
 2290   if(narg /= 1)then
          irtc=itfmessage(9,'General::narg','"1"')
          go to 7000
        endif
        isp0=isp1+1
        call tfsharedsize(isp0,k%k,nsize,irtc)
        isp=isp1+1
        if(irtc /= 0)then
          go to 7000
        endif
        kx=dfromr(dble(max(0,nsize-2)*8))
        go to 6900
 2300   if(narg /= 1)then
           irtc=itfmessage(9,'General::narg','"1"')
           go to 7000
        endif
        irtc=0
        kx%k=merge(ktftrue,ktffalse,ktfoperq(k))
        go to 6900
 2310   kx=tfgaussiancoulombu(isp1,irtc)
        go to 6900
 2320   call tfgaussiancoulombfitted(isp1,kx,irtc)
        go to 6900
 2330   call tfrest(isp1,kx,irtc)
        go to 6900
 2340   call tfrotateright1(isp1,kx,irtc)
        go to 6900
 2350   call tfdifference(isp1,kx,irtc)
        go to 6900
 2360   if(narg == 1)then
          kx=tfeintf(gamma0,cgamma0,k,.false.,0.d0,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 2370   if(narg == 1)then
          kx=tfeintf(xsin,tcxsin,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6811
        endif
        go to 6900
 2380   if(narg == 2)then
          kx=tfbessel(isp1,12,irtc)
        else
          go to 6812
        endif
        go to 6900
 2390   if(narg == 4)then
          kx=tfhg(isp1,irtc)
        else
          go to 6814
        endif
        go to 6900
 2400   if(narg == 4)then
          isp=isp+1
          kx=tfhg(isp1,irtc)
        else
          go to 6814
        endif
        go to 6900
 2410   if(narg == 3)then
          kx=tfhg(isp1,irtc)
        else
          go to 6813
        endif
        go to 6900
 2420   if(narg == 3)then
          dtastk(isp+1)=dxnullo
          isp=isp+2
          kx=tfhg(isp1,irtc)
        else
          go to 6813
        endif
        go to 6900
 2430   if(narg == 2)then
          kx=tfhg(isp1,irtc)
        else
          go to 6812
        endif
        go to 6900
 2440   if(narg == 2)then
          dtastk(isp+1)=dxnullo
          isp=isp+3
          kx=tfhg(isp1,irtc)
        else
          go to 6812
        endif
        go to 6900
 2450   if(narg == 1)then
          kx=tfeintf(zeta,czeta,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfeintf2(zeta2,czeta2,dtastk(isp-1),k,irtc)
        else
          go to 6812
        endif
        go to 6900
 2460   if(narg == 1)then
          kx=tfeintf(polygamma,cpolygamma,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfbessel(isp1,5,irtc)
        else
          go to 6812
        endif
        go to 6900
 2470   if(narg == 1)then
          kx=tfeintf(dzeta,dczeta,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfeintf2(dzeta2,dczeta2,dtastk(isp-1),k,irtc)
        else
          go to 6812
        endif
        go to 6900
 2480   if(narg == 1)then
          kx=tfeintf(zeta,czeta,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfbessel(isp1,6,irtc)
        else
          go to 6812
        endif
        go to 6900
 2490   if(narg == 1)then
          kx=tfeintf(dzeta,dczeta,k,.true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfbessel(isp1,7,irtc)
        else
          go to 6812
        endif
        go to 6900
 2500   kx=tfbessel(isp1,8,irtc)
        go to 6900
 2510   if(narg == 1)then
          kx=tfeintf(rgamma,cgamma,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        elseif(narg == 2)then
          kx=tfbessel(isp1,11,irtc)
        else
          go to 6812
        endif
        go to 6900
 2520   if(narg == 3)then
          kx=kxhg(dtastk(isp1+1),dxnullo,dtastk(isp1+2),
     $         dtastk(isp),.true.,3,irtc)
        else
          go to 6812
        endif
        go to 6900
 2530   if(narg == 3)then
          kx=kxhgpq(isp1,.false.,irtc)
        else
          go to 6812
        endif
        go to 6900
 2540   if(narg == 3)then
          kx=kxhgpq(isp1,.true.,irtc)
        else
          go to 6812
        endif
        go to 6900
 2550   if(narg == 2)then
          kx=tfbessel(isp1,13,irtc)
        else
          go to 6812
        endif
        go to 6900
 2560   if(narg == 2)then
          kx=tfbessel(isp1,14,irtc)
        else
          go to 6812
        endif
        go to 6900
 2570   if(narg == 3)then
          kx=kxhg(dtastk(isp1+2),dtastk(isp),dtastk(isp),dtastk(isp1+1),
     $         .false.,4,irtc)
        else
          go to 6812
        endif
        go to 6900
 2580   if(narg == 3)then
          kx=kxhg(dtastk(isp1+2),dtastk(isp),dtastk(isp),dtastk(isp1+1),
     $         .false.,5,irtc)
        else
          go to 6812
        endif
        go to 6900
 2590   if(narg == 1)then
          if(ktfrealq(dtastk(isp),v) .and. anint(v) == v)then
            kx%x(1)=bernbf(nint(v))*factorial(v)
            irtc=0
          endif
        elseif(narg == 2)then
          irtc=-1
          if(ktfrealq(dtastk(isp1+1),v) .and. anint(v) == v .and.
     $         tfnumberq(dtastk(isp),c1))then
            kx=kxcalocc(-1,berpol(nint(v),c1))
            irtc=0
          endif
        else
          go to 6812
        endif
        go to 6900
 8000   continue
 8100   irtc=0
        return
 8200   if(irtc /= 0)then
          go to 6900
        endif
        kx%k=ktfoper+mtfnull
        return
 6000   iaf=int(ka1)
c        go to (
c     $       6010,
c     $       7000,7000,6600,7000,6600,7000,6620,6630,6690,6690,
c     $       6680,6680,6680,6680,6700,6700,6610,6100,6100,6510,
c     $       7000,7000,6010,7000,6640,6650,6750,7000,7000,7000,
c     $       7000,6002,6010,6010,6010,6660,6670,6560,6560,6710,
c     $       6010,6740,7000,7000,6200,6090,6520,6530,6540,6010,
c     $       6010,6550,6570,6570,6570,6570,6580,6580,6590,6720,
c     $       6730,6010,7000,7000,6010,7000),
c     $       iaf+1
c            null
c            m    i    +    -    *    /    v    ^    e    n    
c            >    <    g    l    E    N    ~    &&   o    c
c            [    ]    {    }    s    =    C    (    )    ,
c            ;    &    :    r    d    RepA RepR u    U    S    
c            ?    f    #    ##   .    |    M    MA   A    rept 
c            repn ineq AT   SF   TB   DB   INc  Dec  Part @
c            msgn TagS (*   *)   Hold z
c        go to 6001
        select case(iaf)
        case (mtfnull,mtflist,mtfcolon,mtfhold,mtftagset,
     $       mtfrepeated,mtfrepeatednull,mtfpattest,
     $       mtfrule,mtfruledelayed)
          kx=kxcompose(isp1)
          irtc=0
          return
        case (mtfset)
          kx=tfset(isp1,.true.,irtc)
          go to 6900
        case (mtfplus,mtftimes)
          kx=tfplus(isp1,iaf,irtc)
          go to 6900
        case (mtfrevpower)
          kx=tfrevpower(isp1,irtc)
          go to 6900
        case (mtfpower)
          kx=tfpower(isp1,irtc)
          go to 6900
        case (mtfgreater,mtfless,mtfgeq,mtfleq)
          kx=tfrelation(isp1,iaf,irtc)
          go to 6900
        case (mtfinequality)
          call tfinequality(isp1,kx,irtc)
          go to 6900
        case (mtfpart)
          if(ktflistq(ktastk(isp1+1)))then
            kx=tfpart(isp1+1,.true.,irtc)
            if(irtc == 0)then
              kx=tfeevalref(kx,irtc)
            endif
          elseif(ktfsymbolq(ktastk(isp1+1)))then
            irtc=-1
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"List or composition"')
          endif
          go to 6900
        case (mtfnot)
          kx=tfnot(isp1,iaf,irtc)
          go to 6900
        case (mtfupset,mtfupsetdelayed)
          if(narg /= 2)then
            go to 6812
          endif
          kx=tfupset(dtastk(isp1+1),dtastk(isp),i00,irtc)
          go to 6900
        case (mtffun)
          if(isp == isp1+2)then
            if(ktastk(isp) == ktfoper+mtfnull)then
              kx=kxpfaloc(dtastk(isp-1))
              irtc=0
              return
            endif
          endif
        case (mtfalt,mtfand,mtfor)
          if(narg .lt. 2 .and. iaf == mtfalt)then
            irtc=-1
            go to 6900
          endif
          if(narg == 2)then
            kx=tfeval1(dtastk(isp1+1),dtastk(isp),iaf,irtc)
            go to 6900
          endif
          if(narg == 0)then
            kx=dfromr(dble(mtfor-iaf))
            irtc=0
            return
          endif
          kx=dtastk(isp1+1)
          if(narg == 1)then
            irtc=0
            return
          endif
          do i=isp1+2,isp
            k1=kx
            kx=tfeval1(k1,dtastk(i),iaf,irtc)
            if(irtc /= 0)then
              go to 6900
            endif
          enddo
          return
        case (mtfdot)
          if(narg == 2)then
            kx=tfeval1(dtastk(isp1+1),dtastk(isp),iaf,irtc)
            go to 6900
          endif
          if(narg == 0)then
            irtc=itfmessage(9,'General::narg','"1 or more"')
            go to 6900
          endif
          kx=dtastk(isp)
          if(narg == 1)then
            irtc=0
            return
          endif
          do i=isp-1,isp1+1,-1
            k1=kx
            kx=tfeval1(dtastk(i),k1,iaf,irtc)
            if(irtc /= 0)then
              go to 6900
            endif
          enddo
          return
        case (mtfconcat)
          call tfstringjoin(isp1,kx,irtc)
          go to 6900
        case (mtfmap)
          kx=tfmap(isp1,3,1,irtc)
          go to 6900
        case (mtfmapall)
          call tfmapall(isp1,kx,irtc)
          go to 6900
        case (mtfapply)
          call tfapply(isp1,kx,irtc)
          go to 6900
        case (mtfaddto,mtfsubtractfrom,mtftimesby,mtfdivideby)
          if(narg /= 2)then
            go to 6812
          endif
          kx=tfeval1to(dtastk(isp1+1),dtastk(isp),iaf,.false.,irtc)
          go to 6900
        case (mtfincrement,mtfdecrement)
          v1=merge(1.d0,-1.d0,iaf == mtfincrement)
          if(narg == 1)then
            kx=tfeval1to(dtastk(isp1+1),dfromr(v1),mtfaddto,.true.,irtc)
          elseif(narg == 2 .and.
     $           ktastk(isp1+1) == ktfoper+mtfnull)then
            kx=tfeval1to(dtastk(isp),dfromr(v1),mtfaddto,.false.,irtc)
          else
            go to 6811
          endif
          go to 6900
        case (mtfsetdelayed)
          if(narg /= 2)then
            go to 6812
          endif
          kx=tfset(isp1,.true.,irtc)
          go to 6900
        case (mtfreplace)
          kx=tfreplace1(isp1,irtc)
          go to 6900
        case (mtfreplacerepeated)
          kx=tfreplacerepeated1(isp1,irtc)
          go to 6900
        case (mtfequal,mtfunequal)
          kx=tfequal(isp1,iaf,irtc)
          go to 6900
        case (mtfsame,mtfunsame)
          kx=tfsameq1(isp1,iaf,irtc)
          go to 6900
        case (mtfunset)
          if(narg /= 1)then
            go to 6811
          endif
          isp=isp1+2
          ktastk(isp)=ktfref
          kx=tfset(isp1,.true.,irtc)
          isp=isp1+1
          kx%k=ktfoper+mtfnull
          go to 6900
        case (mtfatt)
          call tfatt(isp1,kx,.true.,irtc)
          go to 6900
        case (mtfmessagename)
          call tfdeval(isp1,dlist(ifunbase+mtfmessagename),kx,1,.false.,euv,irtc)
          go to 6900
        case (mtfflag)
          call tfflagordef(isp1,kx,irtc)
          go to 6900
        case (mtfcomplex)
          if(narg /= 2)then
            go to 6812
          endif
          kx=tfeval1(dtastk(isp1+1),dtastk(isp),mtfcomplex,irtc)
          go to 6900
        case default
          if(iaf .le. mtfend)then
            go to 7000
          endif
          irtc=itfmessage(999,'General::invop',' ')
          return
        end select
      elseif(ktfsymbolqdef(k1%k,symd))then
        if(symd%sym%override /= 0)then
          if(symd%downval /= 0)then
            call tfdeval(isp1,k1,kx,1,.false.,euv,irtc)
            go to 6900
          else
            go to 6800
          endif
        else
          go to 6800
        endif
      elseif(ktflistq(k1,kl1))then
        kx=k1
        if(ktfoperq(kl1%head,kop))then
          id=iget_fun_id(kop)
          select case (id)
          case (-mtffun)
            kx=tfpuref(isp1,kl1,irtc)
            go to 6900
          case (-mtfnull)
            if(kl1%nl == 0)then
              kx=tfsequence(isp1,isp)
              if(ktflistq(kx,klx))then
                kx=tfleval(klx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                kx=tfeevalref(kx,irtc)
              endif
              go to 6900
            elseif(kl1%nl == 1 .and.
     $             ktfnonreallistqo(kl1))then
              kx=kxmakelist(isp1,klx)
              kh=klx%dbody(1)
              klx%head=dtfcopy(kh)
              irtc=0
              return
            endif
          case (-mtflist)
            kx=tfpart(isp1,.true.,irtc)
            if(irtc == 0)then
              if(ktflistq(kx,klx))then
                kx=tfleval(klx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                kx=tfeevalref(kx,irtc)
              endif
            endif
            go to 6900
          case (-mtfmap,-mtfapply)
            if(kl1%nl == 1)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=kl1%dbody(1)
              dtastk(isp2+2:isp2+isp-isp1+1)=dtastk(isp1+1:isp)
              isp=isp+isp2-isp1+1
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          case (nfunappend,nfunprepend,nfuncases,nfundelcases,
     $           nfunselcases,nfundelete,nfunposition,nfunselect,
     $           nfunreppart,nfunextract,nfunswicases)
            if(kl1%nl == 1)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=dtastk(isp1+1)
              dtastk(isp2+2)=kl1%dbody(1)
              if(isp > isp1+1)then
                dtastk(isp2+3:isp2+isp-isp1+1)=dtastk(isp1+2:isp)
              endif
              isp=isp2+isp-isp1+1
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          case (nfuninsert)
            if(kl1%nl == 2)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=dtastk(isp1+1)
              dtastk(isp2+2)=kl1%dbody(1)
              dtastk(isp2+3)=kl1%dbody(2)
              isp=isp2+3
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          end select
          go to 6800
        endif
        do while(ktflistq(kx,klx))
          kx=klx%head
        enddo
        if(ktfsymbolqdef(kx%k,symd))then
          if(symd%sym%override /= 0 .and. symd%downval /= 0)then
            call tfdeval(isp1,kx,kx,1,.false.,euv,irtc)
            go to 6900
          endif
        endif
        go to 6800
      elseif(ktfstringq(k1))then
        if(narg .le. 0 .or. narg > 2)then
          go to 6800
        elseif(ktfnonrealq(ktastk(isp)) .or.
     $         ktfnonrealq(ktastk(isp1+1)))then
          go to 6800
        endif
        kx=kxsubstring(k1,isp1+1,isp)
        irtc=0
        return
      else
        go to 6800
      endif
      go to 6900
 6800 irtc=0
      go to 7000
 6811 irtc=itfmessage(9,'General::narg','"1"')
      go to 7000
 6812 irtc=itfmessage(9,'General::narg','"2"')
      go to 7000
 6813 irtc=itfmessage(9,'General::narg','"3"')
      go to 7000
 6814 irtc=itfmessage(9,'General::narg','"4"')
      go to 7000
 6900 if(irtc == 0 .or. irtc .lt. -1)then
        return
      elseif(irtc == -1)then
        irtc=0
      endif
 7000 isp=isp1+narg
      kx=kxcrelistm(narg,ktastk(isp1+1:isp),dtastk(isp1))
      if(irtc > 0)then
        if(ierrorprint /= 0)then
          call tferrorhandle(kx,irtc)
        else
          call tfdebugprint(kx,'... in',3)
        endif
      elseif(irtc == -1)then
        irtc=0
      endif
      return
      end

      function tfcheck(isp1,irtc) result(kx)
      use eeval
      use tfcsi
      implicit none
      type (sad_descriptor) kx,kf,kxcheckmessage,kxmessagelist
      type (sad_symdef), pointer,save :: symd
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,itgetfpe,itfmessage,narg,irtc1
      data kxmessagelist%k,kxcheckmessage%k /0,0/
      narg=isp-isp1
      if(narg .lt. 2)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      isp0=isp
      rlist(ierrorgen)=0.d0
      kx=tfeevalref(dtastk(isp1+1),irtc)
      isp=isp0
      if(irtc == 0)then
        if(itgetfpe() /= 0)then
          call tclrfpe
          irtc=itfmessage(9,'General::fpe','""')
        endif
      endif
      if(irtc > 0)then
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',0,icslfnm())
        endif
      elseif(irtc == irtcabort .and. rlist(ierrorgen) == 0.d0)then
        irtc=itfmessage(999,'General::abort',' ')
      elseif(irtc .le. irtcret)then
        rlist(ierrorgen)=1.d0
      endif
      if(rlist(ierrorgen) /= 0.d0)then
        if(irtc .le. irtcret)then
          call tfcatchreturn(modethrow,kx,irtc)
c          write(*,*)'tfcheck-1 ',modethrow,irtc
        endif
        if(narg > 2)then
          if(kxcheckmessage%k == 0)then
            kxcheckmessage=kxsymbolz('Check$Message',13)
          endif
          dtastk(isp1+1)=kxcheckmessage
          kf=tfefunref(isp1+1,.false.,irtc1)
          if(irtc1 == 0 .and. ktfrealq(kf) .and. kf%k /= 0)then
            rlist(ierrorgen)=0.d0
            kx=tfeevalref(dtastk(isp1+2),irtc)
          endif
        else
          if(kxmessagelist%k == 0)then
            kxmessagelist=kxsymbolz('$MessageList',12)
            call descr_sad(kxmessagelist,symd)
          endif
          call tflocald(symd%value)
          symd%value=dtfcopy1(dxnulll)
          rlist(ierrorgen)=0.d0
          kx=tfeevalref(dtastk(isp1+2),irtc)
        endif
      endif
      return
      end

      function tfbessel(isp1,mode,irtc) result(kx)
      use gammaf
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        kx=dxnullo
        return
      endif
      select case (mode)
      case (0)
c        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesselj,irtc)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesj,irtc)
      case (1)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesy,irtc)
      case (2)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesi,irtc)
      case (3)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesk,irtc)
      case (4)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cgpolygamma2,irtc)
      case (5)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cpolygamma2,irtc)
      case (6)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,chzeta2,irtc)
      case (7)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,dchzeta2,irtc)
      case (8)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,chpolygamma2,irtc)
      case (9)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cgammaq,irtc)
      case (10)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cgammap,irtc)
      case (11)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cgamma2,irtc)
      case (12)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cpochh,irtc)
      case (13)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cpolylog,irtc)
      case (14)
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbeta2,irtc)
      end select
      return
      end

      recursive subroutine tfbesself(k1,k2,kx,cfun,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: kl1,kl2
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,n1,n2,isp0,isp2,i
      real*8 x1,x2
      complex*16 cx,c1,c2
      complex*16 ,external:: cfun
      if(tfcomplexnumlistqk(k1%k,kl1))then
        n1=kl1%nl
        if(tfcomplexnumlistqk(k2%k,kl2))then
          if(n1 /= kl2%nl)then
            irtc=itfmessage(9,'General::equalleng','"#1 and #2"')
            return
          endif
          isp0=isp
          call tfgetllstkall(kl1)
          call tfgetllstkall(kl2)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfbesself(dtastk(isp0+i),dtastk(isp0+n1+i),
     $           dtastk(isp),cfun,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        elseif(tfnumberq(k2))then
          isp0=isp
          call tfgetllstkall(kl1)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfbesself(dtastk(isp0+i),k2,dtastk(isp),cfun,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        else
          irtc=-1
          return
        endif
      elseif(tfcomplexnumlistqk(k2%k,kl2))then
        n2=kl2%nl
        isp0=isp
        call tfgetllstkall(kl2)
        isp2=isp
        do i=1,n2
          isp=isp+1
          call tfbesself(k1,dtastk(isp0+i),dtastk(isp),cfun,irtc)
          if(irtc /= 0)then
            isp=isp0
            return
          endif
        enddo
        kx=kxmakelist(isp2)
        isp=isp0
      elseif(ktfrealq(k1,x1))then
        if(ktfrealq(k2,x2))then
          cx=cfun(dcmplx(x1,0.d0),dcmplx(x2,0.d0))
        elseif(tfcomplexq(k2,c2))then
          cx=cfun(dcmplx(x1,0.d0),c2)
        else
          irtc=-1
          return
        endif
        go to 10
      elseif(tfcomplexq(k1,c1))then
        if(ktfrealq(k2,x2))then
          cx=cfun(c1,dcmplx(x2,0.d0))
        elseif(tfcomplexq(k2,c2))then
          cx=cfun(c1,c2)
        else
          irtc=-1
          return
        endif
        go to 10
      else
        irtc=-1
        return
      endif
      irtc=0
      return
 10   kx=kxcalocc(-1,cx)
      irtc=0
      return
      end

      subroutine tfvectorqf(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) kxi
      type (sad_dlist), pointer ::kl
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 m,i,ispf,narg,itfmessage
      narg=isp-isp1
      if(narg == 1)then
        kx%k=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp),kl))then
          kx%k=0
        else
          m=kl%nl
          if(m /= 0)then
            if(ktfnonreallistqo(kl))then
              do i=1,m
                if(tflistq(kl%dbody(i)))then
                  kx%k=0
                  return
                endif
              enddo
            endif
          endif
        endif
      elseif(narg == 2)then
        kx%k=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp-1),kl))then
          kx%k=0
        else
          m=kl%nl
          if(m /= 0)then
            ispf=isp+1
            ktastk(ispf)=ktfcopy(ktastk(isp))
            do i=1,m
              isp=ispf+1
              dtastk(isp)=kl%dbody(i)
              if(tflistq(ktastk(ispf+1)))then
                isp=ispf-1
                kx%k=0
                go to 100
              endif
              kxi=tfefunref(ispf,.true.,irtc)
              isp=ispf-1
              if(irtc /= 0)then
                go to 100
              endif
              if(ktfnonrealq(kxi) .or. kxi%k == 0)then
                kx%k=0
                go to 100
              endif
            enddo
 100        call tflocal(ktastk(ispf))
          endif
        endif
      else
        irtc=itfmessage(9,'General::narg','"1 or 2"')
      endif
      return
      end

      subroutine tfextract(isp1,kx,irtc)
      use tfstk
      use eeval,only:tfleval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kh,k,kind
      type (sad_dlist), pointer :: klind,kli
      type (sad_dlist), pointer :: kl
      integer*4 narg,i,isp0,itfmessage,m
      narg=isp-isp1
      if(narg == 2)then
        kh%k=ktfref
      elseif(narg == 3)then
        kh=dtastk(isp)
      elseif(narg == 1)then
        irtc=-1
        return
      else
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
      endif
      k=dtastk(isp1+1)
      kind=dtastk(isp1+2)
      if(.not. ktflistq(k,kl) .or. .not. tflistq(kind,klind))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1, List for #2"')
        return
      endif
      kind=dtfcopy1(kind)
      m=klind%nl
      isp0=isp
      if(ktfnonreallistqo(klind))then
        do i=1,m
          if(ktfnonlistq(klind%dbody(i)))then
            call tfextract1(kl,klind,kh,irtc)
            if(irtc /= 0)then
              go to 9000
            endif
            go to 100
          endif
        enddo
        do i=1,m
          call descr_sad(klind%dbody(i),kli)
          call tfextract1(kl,kli,kh,irtc)
          if(irtc /= 0)then
            go to 9000
          endif
        enddo
 100    kx=kxmakelist(isp0)
      elseif(m == 0)then
        if(kh%k /= ktfref)then
          isp=isp+1
          dtastk(isp)=kh
          isp=isp+1
          dtastk(isp)=k
          kx=tfefunref(isp-1,.true.,irtc)
        else
          kx=tfleval(kl,.true.,irtc)
        endif
      else
        call tfextract1(kl,klind,kh,irtc)
        kx=dtastk(isp)
      endif
 9000 isp=isp0
      call tflocal1d(kind)
      return
      end

      subroutine tfextract1(kl,kll,kh,irtc)
      use tfstk
      use part,only:tfpartrstk
      use eeval,only:tfeevalref
      implicit none
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(in):: kh
      type (sad_dlist) ,intent(inout):: kl
      type (sad_dlist) ,intent(inout):: kll
      integer*4 i,isp0,isp1,isp2,isp3
      logical*4 list,eval
      isp1=isp
      call tfgetllstkall(kll)
      isp2=isp
      call tfpartrstk(kl,isp1,isp2,list,
     $     .false.,.false.,eval,.true.,irtc)
      if(irtc /= 0)then
        isp=isp1
        return
      endif
      isp3=isp1+isp-isp2
      if(kh%k == ktfref)then
        do i=1,isp-isp2
          dtastk(isp1+i)=tfeevalref(dtastk(isp2+i),irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
        enddo
      elseif(kh%k == ktfoper+nfununeval)then
        ktastk(isp1+1:isp1+isp-isp2)=ktastk(isp2+1:isp)
c        do i=1,isp-isp2
c          ktastk(isp1+i)=ktastk(isp2+i)
c        enddo
        irtc=0
      else
        isp0=isp
        do i=1,isp-isp2
          isp=isp0+2
          dtastk(isp-1)=kh
          dtastk(isp)=dtastk(isp2+i)
          dtastk(isp1+i)=tfefunref(isp-1,.true.,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
        enddo
      endif
      isp=isp3
      return
      end

      subroutine tfthread(isp1,kx,mode,irtc)
      use sameq,only:tfconstheadqk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k,kh,kf,ki,kj
      type (sad_dlist), pointer :: list,klx,kli
      type (sad_rlist), pointer :: klj
      integer*8 kai
      integer*4 narg,i,j,n,m,itfmessage,isp0,isp2,isp3,m1,m2
      logical*4 allv,ch
      integer*4 ,parameter :: mth=2**16
      narg=isp-isp1
      if(mode == 0)then
        if(narg == 1)then
          kh%k=ktfoper+mtflist
        elseif(narg == 2)then
          kh=dtastk(isp)
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
          return
        endif
        k=dtastk(isp1+1)
        if(ktfnonlistq(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition for #1"')
          return
        endif
        kf=list%head
      else
        if(narg /= 2)then
          irtc=itfmessage(9,'General::narg','"2"')
          return
        endif
        k=dtastk(isp)
        if(ktfnonlistq(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition for #2"')
          return
        endif
        kh%k=ktfoper+mtflist
        kf=dtastk(isp1+1)
      endif
      m=list%nl
      n=-1
      ch=tfconstheadqk(kf)
      allv=ch .and. kh%k == ktfoper+mtflist
      isp3=isp
      if(ktfnonreallistqo(list))then
        do i=1,m
          isp=isp+1
          ktastk(isp)=ktfref
          ki=list%dbody(i)
          if(ktflistq(ki,kli) .and. kli%head%k == kh%k)then
            ktastk(isp)=ktfaddr(ki)
            if(n .lt. 0)then
              n=kli%nl
            elseif(n /= kli%nl)then
              irtc=itfmessage(9,'General::equalleng',
     $             '"elements in Thread"')
              return
            endif
            allv=allv .and. ktfreallistq(kli)
          else
            allv=.false.
          endif
        enddo
      endif
      if(n .lt. 0)then
        kx=k
        irtc=0
        return
      endif
      isp=isp+1
      isp2=isp
      dtastk(isp)=kh
      if(ch)then
        if(mode /= 2)then
          if(allv)then
            kx=kxadaloc(-1,n,klx)
            if(m < mth)then
              do j=1,n
                klx%dbody(j)=kxavaloc(0,m,klj)
                klj%rbody(1:m)=rlist(ktastk(isp3+1:isp3+m)+j)
                klj%head=dtfcopy(kf)
                klj%attr=ior(klj%attr,lconstlist)
              enddo
            else
              do j=1,n
                klx%dbody(j)=kxavaloc(0,m,klj)
                m1=1
                m2=mth
                do while (m1 < m)
                  klj%rbody(m1:m2)=rlist(ktastk(isp3+m1:isp3+m2)+j)
                  m1=m1+mth
                  m2=min(m2+mth,m)
                enddo
c                do i=1,m
c                  klj%rbody(i)=rlist(ktastk(isp3+i)+j)
c                enddo
                klj%head=dtfcopy(kf)
                klj%attr=ior(klj%attr,lconstlist)
              enddo
            endif
            klx%attr=ior(klx%attr,lconstlist)
            irtc=0
            isp=isp2-1
            return
          else
            do j=1,n
              isp=isp+1
              isp0=isp
              dtastk(isp0)=kf
              do concurrent (i=1:m)
                kai=ktastk(isp3+i)
                if(kai == ktfref)then
                  dtastk(isp0+i)=list%dbody(i)
                else
                  dtastk(isp0+i)=dlist(kai+j)
                endif
c                  dtastk(isp0+i)=merge(list%dbody(i),dlist(kai+j),kai == ktfref)
              enddo
              isp=isp0+m
              dtastk(isp0)=kxcompose(isp0)
              isp=isp0
            enddo
          endif
        endif
      elseif(mode /= 2)then
        do j=1,n
          isp=isp+1
          isp0=isp
          dtastk(isp0)=kf
          do i=1,m
            isp=isp+1
            kai=ktastk(isp3+i)
            if(kai == ktfref)then
              dtastk(isp)=list%dbody(i)
            else
              dtastk(isp)=dlist(kai+j)
            endif
c            dtastk(isp)=merge(list%dbody(i),dlist(kai+j),kai == ktfref)
          enddo
          call tfefunrefstk(isp0,isp0,irtc)
          if(irtc /= 0)then
            isp=isp2-1
            return
          endif
        enddo
      else
        do j=1,n
          isp=isp+1
          isp0=isp
          dtastk(isp0)=kf
          do concurrent (i=1:m)
            kai=ktastk(isp3+i)
            if(kai == ktfref)then
              dtastk(isp0+i)=list%dbody(i)
            else
              dtastk(isp0+i)=dlist(kai+j)
            endif
c            dtastk(isp0+i)=merge(list%dbody(i),dlist(kai+j),kai == ktfref)
          enddo
          isp=isp0+m
          kj=tfefunref(isp0,.true.,irtc)
          if(irtc /= 0)then
            if(irtc == -3)then
              go to 1000
            elseif(irtc == -2)then
              irtc=0
            else
              isp=isp2-1
              return
            endif
          endif
          isp=isp0-1
        enddo
      endif
 1000 if(mode /= 2)then
        kx=tfefunref(isp2,.true.,irtc)
      else
        kx%k=ktfoper+mtfnull
        irtc=0
      endif
      isp=isp2-1
      return
      end

      subroutine tfcontext(isp1,kx,irtc)
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx,k
      type (sad_symbol),pointer :: sym
      type (sad_namtbl),pointer :: loc
      integer*8 ka
      integer*4 isp1,irtc,itfmessage
      if(isp1+1 /= isp)then
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

      subroutine tfsetioch(isp1,kx,irtc)
      use tfcsi
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      real*8 ch,v
      integer*4 ich
      irtc=0
      if(isp .eq. isp1+2 .and. ktfrealq(dtastk(isp1+1),ch) .and. ktfrealq(dtastk(isp1+2),v))then
        ich=int(ch)
        select case (ich)
        case (0)
          if(v .eq. -1.d0)then
            lfni=5
          else
            lfni=int(v)
          endif
        case (1)
          if(v .eq. -1.d0)then
            lfno=6
          else
            lfno=int(v)
          endif
        case (2)
          if(v .eq. -1.d0)then
            lfnm=6
          else
            lfnm=int(v)
          endif
        case default
          irtc=-1
        end select
        kx=dfromr(v)
      elseif(isp .eq. isp1+1 .and. ktfrealq(dtastk(isp1+1),ch))then
        ich=int(ch)
        select case (ich)
        case (0)
          if(lfni .eq. 5)then
            kx=dfromr(-1.d0)
          else
            kx=dfromr(dble(lfni))
          endif
        case (1)
          if(lfno .eq. 6)then
            kx=dfromr(-1.d0)
          else
            kx=dfromr(dble(lfno))
          endif
        case (2)
          if(lfnm .eq. 6)then
            kx=dfromr(-1.d0)
          else
            kx=dfromr(dble(lfnm))
          endif
        case default
          irtc=-1
        end select
      else
        irtc=-1
      endif
      return
      end subroutine tfsetioch

      end module efun

      function tfefunrefu(isp1,irtc) result(kx)
      use efun
      implicit none
      type (sad_descriptor) :: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      kx=tfefunref(isp1,.true.,irtc)
      return
      end function

      function tfefunrefd(isp1,irtc) result(kx)
      use efun
      implicit none
      type (sad_descriptor) :: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      kx=tfefunref(isp1,.false.,irtc)
      return
      end function

      function tfefunreff(isp1,up,irtc) result(kx)
      use efun
      implicit none
      type (sad_descriptor) :: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: up
      kx=tfefunref(isp1,up,irtc)
      return
      end function

      subroutine tfefun(isp1,kx,ref,upvalue,irtc)
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      logical*4 ref,upvalue
      if(ref)then
        kx=tfefunref(isp1,upvalue,irtc)
      else
        call tfefundef(isp1,kx,irtc)
      endif
      return
      end

      subroutine tfefunrefc(isp1,kx,irtc)
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      levele=levele+1
      kx=tfefunref(isp1,.true.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfefunrefstk(isp1,isp2,irtc)
      use efun
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl
      integer*4 isp1,irtc,isp2
      kx=tfefunref(isp1,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      if(ktflistq(kx,kl))then
        if(kl%head%k == ktfoper+mtfnull)then
          isp=isp2-1
          call tfgetllstkall(kl)
          return
        endif
      endif
      isp=isp2
      dtastk(isp)=kx
      return
      end

      subroutine tfefundef(isp1,kx,irtc)
      use tfstk
      use funs
      use eeval
      use tfcx,only:tfatt,tfsolvemember
      use sameq,only:tfrefq
      implicit none
      type (sad_descriptor) kx,k1,tfefun1,tfpuref
      type (sad_dlist), pointer :: kl,kli,klx
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      integer*8 ka1
      integer*4 isp1,irtc,narg,id,itfmessageexp
      logical*4 rep
      irtc=0
      narg=isp-isp1
 1    k1=dtastk(isp1)
      if(ktfoperq(k1,ka1))then
        id=iget_fun_id(ka1)
        if(ka1 > mtfend .and.
     $       id > 999 .and. id .le. 2000)then
          kx=tfefun1(isp1,id,.false.,irtc)
          if(irtc == 0 .and. tfrefq(kx))then
          else
            irtc=itfmessageexp(999,'General::invset',k1)
          endif
          go to 6900
        elseif(ka1 == mtfatt)then
          call tfatt(isp1,kx,.false.,irtc)
          go to 6900
        else
          kx=kxcompose(isp1)
          return
        endif
      elseif(ktfsymbolqdef(k1%k,symd))then
        if(symd%sym%override /= 0)then
          if(ktfoperq(symd%value))then
            dtastk(isp1)=symd%value
            go to 1
          endif
        endif
        go to 7000
      elseif(ktflistq(k1,kl))then
        kx=k1
        if(kl%head%k == ktfoper+mtffun)then
          kx=tfpuref(isp1,kl,irtc)
          go to 6900
        elseif(kl%head%k == ktfoper+mtfnull)then
          if(kl%nl == 0)then
            kx=tfsequence(isp1,isp)
            if(ktflistq(kx,klx))then
              kx=tfleval(klx,.false.,irtc)
            elseif(ktfsymbolq(kx,sym))then
              if(sym%override == 0)then
                sym=>tfsydef(sym)
                kx=sad_descr(sym)
              endif
            endif
            go to 6900
          elseif(kl%nl == 1 .and.
     $           ktfnonreallistqo(kl))then
            kx=kxmakelist(isp1,klx)
            klx%head=dtfcopy(kl%head)
            return
          endif
        elseif(kl%head%k == ktfoper+mtflist)then
          go to 6900
        elseif(ktflistq(kl%head,kli))then
          k1=tfsolvemember(kli,rep,irtc)
          if(irtc == 0)then
            dtastk(isp1)=k1
          endif
        endif
        go to 7000
      elseif(ktfstringq(k1))then
        irtc=itfmessageexp(999,'General::invset',k1)
      else
        go to 7000
      endif
 6900 if(irtc == 0 .or. irtc .lt. -1)then
        return
      elseif(irtc == -1)then
        irtc=0
      endif
 7000 isp=isp1+narg
      kx=kxcrelistm(narg,ktastk(isp1+1:isp1+narg),dtastk(isp1))
      if(irtc > 0)then
        call tferrorhandle(kx,irtc)
      elseif(irtc == -1)then
        irtc=0
      endif
      return
      end

c     Type fixed function proxy for generic math function/vendor extension
c     DOUBLE instance of ASINH in Fortran2008
      real*8 function tasinh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tasinh=asinh(x)
      return
      end
c     DOUBLE instance of ACOSH in Fortran2008
      real*8 function tacosh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tacosh=acosh(x)
      return
      end
c     DOUBLE instance of ATANH in Fortran2008
      real*8 function tatanh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tatanh=atanh(x)
      return
      end
c     DOUBLE COMPLEX proxy of ZTAN vendor extension
      complex*16 function tctan(z)
      implicit none
      complex*16 z,ztan
      tctan=ztan(z)
      return
      end

      subroutine tfefundummy
      use mackw
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

      subroutine tftocontext(isp1,kx,irtc)
      use efun
      use context
      use funs
      use eeval
      use repl,only:tfgetoptionstk
      implicit none
      type (sad_descriptor) kx,k,kc,ki,ks,ic
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      type (sad_dlist), pointer :: kl
      integer*8 kai,ktsydefc,ktfsymbolc,ktcontaloc
      integer*4 isp1,irtc,itfmessage,i,isp0,nc,ispopt,isp2
      character*4 optname(1)
      type (sad_descriptor) kaopt(1)
      save kaopt,optname
      data kaopt%k /0/
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
          if(kai == mtfnull)then
            if(i == isp1+1)then
              ic%k=0
            endif
            cycle LOOP_I
          endif
          call loc_symstr(klist(klist(ifunbase+kai)),str)
        elseif(ktflistq(ki,kl))then
          if(kl%head%k == ktfoper+mtfslot)then
            k=tfslot(mtfslot,kl,.false.,irtc)
            if(irtc /= 0)then
              return
            endif
            if(ktfsymbolq(k,sym))then
              call sym_symstr(sym,str)
            elseif(ktfoperq(ki,kai))then
              if(kai == mtfnull)then
                if(i == isp1+1)then
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
        if(i /= isp2 .and. ic%k /= 0)then
          if(str%str(nc:nc) /= '`')then
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
      if(ispopt > isp0)then
        isp=isp0
        kx=tfsyeval(ks,irtc)
      else
        isp=isp+1
        dtastk(isp)=ks
        kx=tfefunref(isp0+1,.true.,irtc)
        isp=isp0
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Symbol or Character-string"')
      return
      end

      integer*8 function ktfcopy(k)
      use tfstk, ktfc=>ktfcopy
      implicit none
      integer*8 ,intent(in):: k
      ktfcopy=ktfc(k)
      return
      end
