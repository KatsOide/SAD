      module take
      use tfstk

      contains
      function tftake(k,kn,take0,eval,irtc) result(kx)
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k,kn
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: klr,kln
      integer*4 ,intent(out):: irtc
      integer*4 m,iv,n1,n2,mn,mx,itfmessage,i
      logical*4 ,intent(in):: take0,eval
      logical*4 take,list,d
      kx%k=ktfoper+mtfnull
      if(ktfnonlistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      irtc=0
      m=kl%nl
      if(m == 0)then
        kx=k
        return
      endif
      take=take0
      if(ktfrealq(kn,iv))then
        if(iv .lt. 0)then
          n1=m+iv+1
          n2=m
        else
          n1=1
          n2=iv
        endif
      elseif(tfreallistq(kn%k,kln))then
        mn=kln%nl
        if(mn .gt. 2)then
          irtc=itfmessage(9,'General::wrongleng',
     $         '"range-list","less than 3"')
          return
        endif
        if(mn == 0)then
          n1=1
          n2=0
        elseif(mn == 1)then
          n1=int(kln%rbody(1))
          if(n1 .lt. 0)then
            n1=m+n1+1
          endif
          n2=n1
        else
          n1=int(kln%rbody(1))
          if(n1 .lt. 0)then
            n1=m+n1+1
          endif
          n2=int(kln%rbody(2))
          if(n2 .lt. 0)then
            n2=m+n2+1
          endif
        endif
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real number or a list of reals for #2"')
        return
      endif
      list=kl%head%k == ktfoper+mtflist
      if(n1 .lt. 1 .or. n2 .lt. 1 .or.
     $     n1 .gt. m .or. n2 .gt. m)then
        if(take)then
          if(list)then
            kx=dxnulll
          else
            kx=kxaaloc(-1,0,klx)
            go to 1
          endif
          return
        else
          if(n2 .lt. 1 .or. n1 .gt. m)then
            kx=k
            return
          else
            n1=max(1,n1)
            n2=min(m,n2)
          endif
        endif
      endif
      if(n2 .lt. n1)then
        if(take)then
          if(list)then
            kx=dxnulll
            return
          else
            kx=kxaaloc(-1,0,klx)
            go to 1
          endif
        else
          kx=k
        endif
        return
      endif
      if(take)then
        mx=n2-n1+1
      else
        if(n1 == 1)then
          n1=n2+1
          n2=m
          mx=m-n1+1
          take=.true.
        elseif(n2 == m)then
          n2=n1-1
          n1=1
          mx=n2
          take=.true.
        else
          mx=m-(n2-n1+1)
        endif
      endif
      if(mx == m)then
        kx=k
        return
      elseif(mx .gt. 0)then
        kx=kxavaloc(-1,mx,klr)
        call descr_sad(kx,klx)
        if(take)then
          if(ktfreallistq(kl))then
            klr%rbody(1:mx)=kl%rbody(n1:n1+mx-1)
            klr%attr=ior(kl%attr,kconstarg)
          else
            d=.false.
            do i=1,mx
              klr%dbody(i)=dtfcopy(kl%dbody(n1+i-1))
              d=d .or. ktfnonrealq(klr%dbody(i))
            enddo
            if(d)then
              klr%attr=ior(iand(kl%attr,kconstarg+lconstlist),
     $             lnonreallist)
            endif
          endif
        else
          if(ktfreallistq(kl))then
            klr%rbody(1:n1-1)=kl%rbody(1:n1-1)
            klr%dbody(n1:m+n1-n2-1)=kl%dbody(n2+1:m)
            klr%attr=ior(kl%attr,kconstarg)
          else
            d=.false.
            do i=1,n1-1
              klx%dbody(i)=dtfcopy(kl%dbody(i))
              d=d .or. ktfnonrealq(klx%dbody(i))
            enddo
            do i=1,m-n2
              klx%dbody(n1+i-1)=dtfcopy(kl%dbody(n2+i))
              d=d .or. ktfnonrealq(kl%dbody(n2+i))
            enddo
            if(d)then
              klx%attr=ior(iand(kl%attr,kconstarg+lconstlist),
     $             lnonreallist)
            endif
          endif
        endif
      else
        if(list)then
          kx=dxnulll
          return
        endif
        kx=kxaaloc(-1,0,klx)
      endif
      if(list)then
        return
      endif
 1    klx%head=dtfcopy(kl%head)
      if(eval)then
        kx=tfleval(klx,.true.,irtc)
      endif      
      return
      end

      subroutine tfreverse(isp1,kx,irtc)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: list,listx
      integer*4 m,itfmessage
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonlistq(ktastk(isp),list))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      m=list%nl
      if(m .le. 1)then
        kx=dtastk(isp)
        irtc=0
        return
      endif
      if(iand(list%attr,lnonreallist) == 0)then
        call loc_sad(ktavaloc(-1,m),listx)
c        do i=1,m
        listx%dbody(m:1:-1)=list%dbody(1:m)
c        enddo
      else
        call loc_sad(ktadaloc(-1,m),listx)
c        do i=1,m
        call ktfcopym(list%body(1:m))
        listx%dbody(m:1:-1)=list%dbody(1:m)
c        enddo
      endif
      listx%head=dtfcopy(list%head)
      listx%attr=list%attr
      irtc=0
      kx=tfleval(listx,.true.,irtc)
      return
      end

      subroutine tfrotateright1(isp1,kx,irtc)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: klr
      integer*4 i,m,itfmessage,n
      if(isp /= isp1+1 .and. isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      elseif(ktfnonlistq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      call loc_sad(ktfaddr(ktastk(isp1+1)),kl)
      m=kl%nl
      if(m .le. 1)then
        go to 8000
      endif
      if(isp == isp1+1)then
        n=1
      elseif(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real for #2"')
        return
      else
        n=int(rtastk(isp))
      endif
      if(n .lt. 0)then
        n=m-(mod(-n-1,m)+1)
      else
        n=mod(n-1,m)+1
      endif
      if(n == 0)then
        go to 8000
      endif
      if(ktfreallistq(kl))then
        kx=kxavaloc(-1,m,klr)
        call descr_sad(kx,klx)
        klr%rbody(1:n)=kl%rbody(m-n+1:m)
        klr%rbody(n+1:m)=kl%rbody(1:m-n)
      else
        kx=kxadaloc(-1,m,klx)
        do i=1,n
          klx%dbody(i)=dtfcopy(kl%dbody(m-n+i))
        enddo
        do i=n+1,m
          klx%dbody(i)=dtfcopy(kl%dbody(i-n))
        enddo
      endif
      klx%head=dtfcopy(kl%head)
      klx%attr=kl%attr
      kx=tfleval(klx,.true.,irtc)
      return
 8000 kx=dtastk(isp1+1)
      irtc=0
      return
      end

      subroutine tfdifference(isp1,kx,irtc)
      use complex
      use eeval
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k0,k1,ks,krv
      type (sad_dlist), pointer :: klx
      type (sad_dlist), pointer :: kl
      type (sad_rlist), pointer :: klr
      integer*4 i,m,itfmessage,isp0
      krv%x(1)=-1.d0
      if(isp == isp1+2)then
        if(.not. ktfrealq(dtastk(isp),krv%x(1)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real for #2"')
          return
        endif
      elseif(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      if(.not. ktflistq(dtastk(isp1+1),kl))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      m=kl%nl
      if(m .le. 1)then
        irtc=itfmessage(9,'General::wrongleng',
     $       '"#1","longer than 1"')
        return
      endif
      if(ktfreallistq(kl))then
        kx=kxavaloc(-1,m-1,klr)
        call c_f_pointer(c_loc(klr),klx)
        if(krv%x(1) == -1.d0)then
          klr%rbody(1:m-1)=kl%rbody(2:m)-kl%rbody(1:m-1)
        elseif(krv%x(1) == 1.d0)then
          klr%rbody(1:m-1)=kl%rbody(2:m)+kl%rbody(1:m-1)
        else
          klr%rbody(1:m-1)=kl%rbody(2:m)+krv%x(1)*kl%rbody(1:m-1)
        endif
      else
        isp0=isp
        k0=kl%dbody(1)
        do i=2,m
          isp=isp+1
          k1=kl%dbody(i)
          ks=tfecmplxl(krv,k0,mtfmult)
          dtastk(isp)=tfecmplxl(ks,k1,mtfplus)
          k0=k1
        enddo
        kx=kxmakelist(isp0,klx)
        isp=isp0
      endif
      klx%head=dtfcopy(kl%head)
      kx=tfleval(klx,.true.,irtc)
      return
      isp=isp0
      return
      end

      subroutine tfrest(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonlistq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      kx=tftake(dtastk(isp),sad_descr(1.d0),.false.,.true.,irtc)
      return
      end

      subroutine tfprotect(isp1,kx,protect,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 i,itfmessageexp
      logical*4 ,intent(in):: protect
      kx%k=ktfoper+mtfnull
      irtc=0
      LOOP_I: do i=isp1+1,isp
        k=dtastk(i)
        if(ktfoperq(k,ka))then
          k%k=ktfsymbol+klist(ifunbase+ka)
        endif
        if(ktfsymbolq(k,sym))then
          if(sym%override == 0)then
            sym=>tfsydef(sym)
          endif
          if(sym%gen .le. 0)then
            if(protect)then
              sym%attr=ior(sym%attr,iattrprotected)
            else
              if(iand(sym%attr,iattrconstant) /= 0)then
                irtc=itfmessageexp(9,'General::unprotconst',k)
              endif
              sym%attr=ior(sym%attr,iattrprotected+iattrconstant)
     $             -iattrprotected-iattrconstant
            endif
          endif
        endif
      enddo LOOP_I
      return
      end

      end module take

      module attrib
      use tfstk
      implicit none
      type (sad_descriptor), save :: iaxnone,
     $     iaxholdall,iaxholdfirst,iaxholdrest,iaxholdnone,
     $     iaxconstant,iaximmediate1,iaxdynamic,iaxordless,
     $     iaxprotected,iaxnumeric
      data iaxnone%k /0/

      contains
        subroutine tfattrinit
        iaxnone      =kxsymbolf('None',4,.true.)
        iaxholdall   =kxsymbolf('HoldAll',7,.true.)
        iaxholdfirst =kxsymbolf('HoldFirst',9,.true.)
        iaxholdrest  =kxsymbolf('HoldRest',8,.true.)
        iaxholdnone  =kxsymbolf('HoldNone',8,.true.)
        iaxconstant  =kxsymbolf('Constant',8,.true.)
        iaximmediate1=kxsymbolf('Immediate',9,.true.)
        iaxnumeric   =kxsymbolf('Numeric',7,.true.)
        iaxdynamic   =kxsymbolf('Dynamic',7,.true.)
        iaxordless   =kxsymbolf('OrderLess',9,.true.)
        iaxprotected =kxsymbolf('Protected',9,.true.)
        return
        end subroutine tfattrinit

      subroutine tfattributes(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 iat,isp0,itfmessage
      if(iaxnone%k == 0)then
        call tfattrinit
      endif
      if(isp /= isp1+1)then
        irtc=itfmessage(999,'General::narg','"1"')
        return
      endif
      k=dtastk(isp)
      If(.not. ktfsymbolq(k,sym))then
        if(ktfoperq(k,ka))then
          call loc_sad(klist(ifunbase+ka),sym)
        else
          irtc=itfmessage(9,'General::wrongtype','"Symbol or Operator"')
          return
        endif
      endif
      iat=sym%attr
      isp0=isp
      if(iat == 0)then
        isp=isp+1
        dtastk(isp)=iaxnone
      else
        if(iand(iat,iattrholdfirst) /= 0)then
          isp=isp+1
          dtastk(isp)=merge(iaxholdall,iaxholdfirst,
     $         iand(iat,iattrholdrest) /= 0)
        else
          isp=isp+1
          dtastk(isp)=merge(iaxholdrest,iaxholdnone,
     $         iand(iat,iattrholdrest) /= 0)
        endif
        if(iand(iat,iattrconstant) /= 0)then
          isp=isp+1
          dtastk(isp)=iaxconstant
        endif
        if(iand(iat,iattrprotected) /= 0)then
          isp=isp+1
          dtastk(isp)=iaxprotected
        endif
        if(iand(iat,iattrimmediate) /= 0)then
          isp=isp+1
          dtastk(isp)=iaximmediate1
        endif
        if(iand(iat,iattrdynamic) /= 0)then
          isp=isp+1
          dtastk(isp)=iaxdynamic
        endif
        if(iand(iat,iattrorderless) /= 0)then
          isp=isp+1
          dtastk(isp)=iaxordless
        endif
      endif
      kx=kxmakelist(isp0)
      irtc=0
      return
      end

      recursive subroutine tfsetattributes(isp1,kx,irtc)
      use tfcode
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k,kv
      type (sad_dlist), pointer :: kl
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 narg,iattrib,i,itfmessage,isp0,
     $     itfmessageexp
      logical*4 prot
      if(iaxnone%k == 0)then
        call tfattrinit
      endif
      narg=isp-isp1
      if(narg /= 2)then
        irtc=itfmessage(999,'General::narg','"2"')
        return
      endif
      k=dtastk(isp1+1)
      irtc=0
      if(ktfoperq(k,ka))then
        k%k=ktfsymbol+klist(ifunbase+ka)
      endif
      if(ktfsymbolq(k,sym))then
        if(sym%override == 0 .or. sym%override == 1)then
          sym=>tfsydef(sym)
          ka=ksad_loc(sym%loc)
        endif
        kv=dtastk(isp)
        if(tflistq(kv,kl))then
          isp0=isp
          do i=1,kl%nl
            isp=isp0+1
            dtastk(isp)=k
            isp=isp+1
            dtastk(isp)=kl%dbody(i)
            call tfsetattributes(isp0,kx,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          enddo
          isp=isp0
          return
        endif
        if(.not. ktfsymbolq(kv))then
          go to 9000
        endif
        prot=ktfprotectedq(ktfaddrd(k)) .and. sym%override /= -3
        iattrib=sym%attr
        if(tfsamesymbolq(kv,iaxnone))then
          iattrib=0
        elseif(tfsamesymbolq(kv,iaxholdall))then
          iattrib=ior(iattrib,iattrholdall)
        elseif(tfsamesymbolq(kv,iaxholdfirst))then
          iattrib=iand(-iattrholdall-1,iattrib)+iattrholdfirst
        elseif(tfsamesymbolq(kv,iaxholdrest))then
          iattrib=iand(-iattrholdall-1,iattrib)+iattrholdrest
        elseif(tfsamesymbolq(kv,iaxholdnone))then
          iattrib=iand(-iattrholdall-1,iattrib)
        elseif(tfsamesymbolq(kv,iaxconstant))then
          if(sym%gen .gt. 0)then
            irtc=itfmessage(9,'General::localconst',' ')
            return
          endif
          prot=prot .and. iand(iattrib,iattrconstant) == 0
          iattrib=ior(iattrib,iattrconstant+iattrprotected)
        elseif(tfsamesymbolq(kv,iaxprotected))then
          prot=prot .and. iand(iattrib,iattrprotected) == 0
          iattrib=ior(iattrib,iattrprotected)
        elseif(tfsamesymbolq(kv,iaximmediate1))then
          iattrib=ior(iattrib,iattrimmediate)
        elseif(tfsamesymbolq(kv,iaxnumeric))then
          iattrib=ior(iattrib,iattrnumeric)
        elseif(tfsamesymbolq(kv,iaxdynamic))then
          iattrib=ior(iattrib,iattrdynamic)
        elseif(tfsamesymbolq(kv,iaxordless))then
          iattrib=ior(iattrib,iattrorderless)
        else
          go to 9000
        endif
        if(prot)then
          go to 9100
        endif
        sym%attr=iattrib
      elseif(tflistq(k,kl))then
        isp0=isp
        do i=1,kl%nl
          isp=isp+1
          dtastk(isp)=kl%dbody(i)
          isp=isp+1
          ktastk(isp)=ktastk(isp0)
          call tfsetattributes(isp0,kx,irtc)
          isp=isp0
          if(irtc /= 0)then
            return
          endif
        enddo
      endif
      kx%k=ktfoper+mtfnull
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $         '"Attribute",'//
     $'" to be None, HoldAll, HoldFirst, HoldRest, "'//
     $'"HoldNone, Immediate, Constant, OrderLess, or Dynamic"')
      return
 9100 irtc=itfmessageexp(999,'General::protect',k)
      return
      end

      subroutine tfreleasehold(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer ::klx
      integer*4 isp0
      irtc=-1
      if(isp .le. isp1)then
        return
      endif
      isp0=isp
      call tfreleaseholdstk(isp1,isp0,irtc)
      if(irtc == 0)then
        if(isp == isp0+1)then
          kx=dtastk(isp)
        else
          kx=kxcompose(isp0,klx)
          klx%head%k=ktfoper+mtfnull
        endif
      endif
      isp=isp0
      return
      end

      recursive subroutine tfreleaseholdstk(isp1,isp2,irtc)
      use eeval
      implicit none
      type (sad_dlist), pointer :: kli
      integer*4 isp1,irtc,i,isp0,isp2,isp3,j,isp4
      if(isp .le. isp1)then
        irtc=0
        return
      endif
      isp0=isp
      do i=isp1+1,isp2
        isp=isp+1
        dtastk(isp)=dtastk(i)
        if(ktflistq(dtastk(i),kli))then
          isp=isp-1
          if(kli%head%k == ktfoper+mtfhold)then
            call tfevallstkall(kli,.true.,.true.,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          else
            isp4=isp
            call tfgetllstkall(kli)
            isp=isp+1
            dtastk(isp)=kli%head
            isp3=isp
            call tfreleaseholdstk(isp4,isp-1,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
            if(kli%head%k == ktfoper+mtfnull)then
              do j=isp3+1,isp
                isp=isp4+j-isp3
                ktastk(isp)=ktastk(j)
              enddo
            else
              call tfefunrefstk(isp3,isp4+1,irtc)
              if(irtc /= 0)then
                isp=isp0
                return
              endif
            endif
          endif
        endif
      enddo
      return
      end

      end module attrib
