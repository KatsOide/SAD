      recursive logical*4 function tfsameqk(ka,kp) result(lx)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_list), pointer :: kla,klp
      type (sad_symbol), pointer :: syma,symp
      type (sad_string), pointer :: stra,strp
      type (sad_pat), pointer :: pata,patp
      integer*8 nc,ka,kp
      logical*4 tfsamelistqo,tfsamesymbolqk,tfsamesymbolqo
      if(ka .eq. kp)then
        lx=.true.
        return
      endif
      lx=.false.
      if(ktfrealq(ka) .or. ktfoperq(ka) .or.
     $     ktftype(ka) .ne. ktftype(kp))then
        return
      endif
      if(ktfsymbolq(ka,syma))then
        call loc_sad(ktfaddr(kp),symp)
        lx=tfsamesymbolqo(syma,symp)
      elseif(ktfstringq(ka,stra))then
        call loc_sad(ktfaddr(kp),strp)
        nc=stra%nch
        if(nc .eq. strp%nch)then
          lx=stra%str(1:nc) .eq. strp%str(1:nc)
        endif
      elseif(ktflistq(ka,kla))then
        call loc_sad(ktfaddr(kp),klp)
        if(kla%nl .eq. klp%nl)then
          lx=tfsamelistqo(kla,klp)
        endif
      elseif(ktfpatq(ka,pata))then
        call loc_sad(ktfaddr(kp),patp)
        if(.not. tfsameqk(pata%expr%k,patp%expr%k))then
          return
        endif
        if(.not. tfsameqk(pata%head%k,patp%head%k))then
          return
        endif
        if(.not. tfsameqk(pata%default%k,patp%default%k))then
          return
        endif
        lx=tfsamesymbolqk(ktfaddr(pata%sym%alloc),
     $       ktfaddr(patp%sym%alloc))
      endif
      return
      end

      logical*4 function tfsameqd(ka,kp)
      use tfstk
      implicit none
      type (sad_descriptor) ka,kp
      logical*4 tfsameqk
      tfsameqd=tfsameqk(ka%k,kp%k)
      return
      end

      logical*4 function tfsamestringqk(ka1,kp1)
      use tfstk
      implicit none
      type (sad_string), pointer :: stra,strp
      integer*8 ka1,kp1
      logical*4 tfsamestringqo
      if(ktfaddr(ka1) .eq. ktfaddr(kp1))then
        tfsamestringqk=.true.
      else
        call loc_sad(ktfaddr(ka1),stra)
        call loc_sad(ktfaddr(kp1),strp)
        tfsamestringqk=tfsamestringqo(stra,strp)
      endif
      return
      end

      logical*4 function tfsamesymbolqk(ka1,kp1)
      use tfstk
      implicit none
      type (sad_symbol) ,pointer :: syma,symp
      integer*8 ka1,kp1,ka,kp
      logical*4 tfsamesymbolqo
      ka=ktfaddr(ka1)
      kp=ktfaddr(kp1)
      if(ka .eq. kp)then
        tfsamesymbolqk=.true.
      elseif(ka .eq. 0 .or. kp .eq. 0)then
        tfsamesymbolqk=.false.
      else
        call loc_sad(ka,syma)
        call loc_sad(kp,symp)
        tfsamesymbolqk=tfsamesymbolqo(syma,symp)
      endif
      return
      end

      logical*4 function tfsamesymbolqd(k1,k2)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2
      logical*4 tfsamesymbolqk
      tfsamesymbolqd=tfsamesymbolqk(k1%k,k2%k)
      return
      end

      logical*4 function tfsamestringqo(sa,sp)
      use tfcode
      implicit none
      type(sad_string) sa,sp
      tfsamestringqo=sa%nch .eq. sp%nch .and.
     $     sa%str(:sa%nch) .eq. sp%str(:sp%nch)
      return
      end

      logical*4 function tfsamesymbolqo(sa,sp)
      use tfcode
      implicit none
      type(sad_symbol) sa,sp
      tfsamesymbolqo=sa%loc .eq. sp%loc .and.
     $     max(0,sa%gen) .eq. max(0,sp%gen)
      return
      end

      logical*4 function tfsamelistqo(lista,listp)
      use tfstk
      implicit none
      type (sad_list) lista,listp
      integer*8 kai,kpi,kaai,kapi
      integer*4 i,m
      logical*4 tfsameqk
      tfsamelistqo=.false.
      m=lista%nl
      if(m .ne. listp%nl)then
        return
      endif
      do i=0,m
        kai=lista%body(i)
        kpi=listp%body(i)
        if(kai .ne. kpi)then
          if(ktfobjq(kai))then
            if(tfsameqk(kai,kpi))then
              kaai=ktfaddr(kai)
              kapi=ktfaddr(kpi)
              if(ilist(1,kapi-1) .ge. ilist(1,kaai-1))then
                call tflocal1(kai)
                lista%body(i)=ktfcopy1(kpi)
              else
                call tflocal1(kpi)
                listp%body(i)=ktfcopy1(kai)
              endif
              cycle
            endif
          endif
          return
        endif
      enddo
      tfsamelistqo=.true.
      return
      end

      recursive logical*4 function tfnearlysameqf(k1,k2,re,ae)
     $     result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2
      type (sad_list), pointer :: kl1,kl2
      integer*4 i,n1
      real*8 re,ae,s,d
      complex*16 cx1,cx2
      logical*4 tfsameqd
      lx=.false.
      if(tfnumberqd(k1,cx1))then
        if(tfnumberqd(k2,cx2))then
          s=abs(cx1)+abs(cx2)
          d=abs(cx2-cx1)
          if(d .le. s*re .or. d .le. ae)then
            lx=.true.
          endif
        endif
      elseif(tfnumberqd(k2))then
      elseif(k1%k .ne. k2%k)then
      elseif(ktfnonlistqd(k1,kl1))then
        lx=tfsameqd(k1,k2)
      elseif(ktfnonlistqd(k2,kl2))then
      else
        n1=kl1%nl
        if(n1 .ne. kl2%nl)then
          return
        endif
        if(.not. tfnearlysameqf(kl1%dbody(0),kl2%dbody(0),re,ae))then
          return
        endif
        do i=1,n1
          if(.not. tfnearlysameqf(kl1%dbody(i),kl2%dbody(i),re,ae))then
            return
          endif
        enddo
        lx=.true.
      endif
      return
      end

      subroutine tfrealqk(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      if(ktfrealq(ktastk(isp)))then
        kx%k=ktftrue
      else
        kx%k=0
      endif
      return
      end

      logical*4 function tfexprqk(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      tfexprqk=ktflistqd(k,kl) .and. kl%head .ne. ktfoper+mtflist
      return
      end

      logical*4 function tfsameheadqk(k1,k2)
      use tfstk
      implicit none
      integer*8 k1,k2
      logical*4 tfsameqk
      tfsameheadqk=tfsameqk(klist(ktfaddr(k1)),klist(ktfaddr(k2)))
      return
      end

      recursive logical*4 function tfconstqk(k) result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      type (sad_symdef), pointer ::symd
      type (sad_pat), pointer :: pat
      logical*4 tfconstlistqo
      lx=.true.
      if(ktfsymbolqdef(k%k,symd))then
        lx=ktfconstantsymq(symd%sym) .and.
     $       symd%value%k .eq. ktfsymbol+ktfaddrd(k)
     $       .and. symd%upval .eq. 0
      elseif(ktflistqd(k,kl))then
        lx=tfconstlistqo(kl)
      elseif(ktfpatqd(k,pat))then
        lx=tfconstqk(pat%expr)
      endif
      return
      end

      logical*4 function tfconstlistqo(list)
      use tfstk
      implicit none
      type (sad_list) list
      type (sad_descriptor) kh
      integer*4 i
      logical*4 tfconstqk,tfconstheadqk,tfseqqo,nr
      if(iand(lnoconstlist,list%attr) .ne. 0)then
        tfconstlistqo=.false.
        return
      endif
      tfconstlistqo=.true.
      if(iand(lconstlist,list%attr) .ne. 0)then
        return
      endif
      kh=list%dbody(0)
      nr=ktfnonreallistqo(list)
      if(kh%k .eq. ktfoper+mtfhold .or.
     $     kh%k .eq. ktfoper+mtffun)then
        tfconstlistqo=.not. tfseqqo(list)
      elseif(kh%k .eq. ktfoper+mtfcomplex)then
        tfconstlistqo=.not. nr
      elseif(.not. tfconstheadqk(kh))then
        tfconstlistqo=.false.
      elseif(nr)then
        if(iand(kconstarg,list%attr) .eq. 0)then
          do i=1,list%nl
            if(.not. tfconstqk(list%body(i)))then
              tfconstlistqo=.false.
              list%attr=ior(lnoconstlist,list%attr)
              return
            endif
          enddo
          list%attr=ior(list%attr,kconstarg)
        endif
      endif
      if(tfconstlistqo)then
c        call tfdebugprint(ktflist+ka,'constlist',3)
        list%attr=ior(lconstlist+kconstarg,list%attr)
      else
c        call tfdebugprint(ktflist+ka,'nonconstlist',3)
        list%attr=ior(lnoconstlist,list%attr)
      endif
      return
      end

      logical*4 function tfconstheadqk(k)
      use ophash
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: list
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symd
      integer*8 ka
      logical*4 tfconstlistqo,tfconstqk
      tfconstheadqk=.true.
      if(ktfoperqd(k,ka))then
        if(ka .le. mtfend)then
          tfconstheadqk=constop(ka)
        else
          tfconstheadqk=.false.
        endif
      elseif(ktfstringqd(k))then
        tfconstheadqk=.false.
      elseif(ktfsymbolqdef(k%k,symd))then
        tfconstheadqk=ktfconstantsymq(symd%sym) .and.
     $       symd%value%k .eq. ktfsymbol+ka
     $       .and. symd%downval .eq. 0
      elseif(ktflistqd(k,list))then
        if(list%head .eq. ktfoper+mtflist .or.
     $       list%head .eq. ktfoper+mtffun)then
          tfconstheadqk=.false.
        elseif(iand(lconstlist,list%attr) .eq. 0)then
          if(iand(lnoconstlist,list%attr) .ne. 0)then
            tfconstheadqk=.false.
          elseif(.not. tfconstlistqo(list))then
            tfconstheadqk=.false.
          endif
        endif
      elseif(ktfpatqd(k,pat))then
        tfconstheadqk=tfconstqk(pat%expr%k)
      endif
      return
      end

      recursive logical*4 function tfseqqo(list) result(lx)
      use tfstk
      use tfcode
      implicit none
      type (sad_list) list
      type (sad_list), pointer :: kli
      integer*4 iadv,i
      logical*4 tfconstlistqo
      lx=.false.
      iadv=list%attr
      if(iand(lnonreallist,iadv) .eq. 0)then
        return
      endif
      if(iand(iadv,kconstarg+knoseqarg) .ne. 0)then
        return
      elseif(iand(iadv,kseqarg) .ne. 0)then
        lx=.true.
        return
      endif
      do i=1,list%nl
        if(ktflistq(list%body(i),kli))then
          if(kli%head .eq. ktfoper+mtfnull)then
            list%attr=ior(iadv,kseqarg)
            lx=.true.
            return
          endif
          if(tfconstlistqo(kli))then
          elseif(tfseqqo(kli))then
            list%attr=ior(iadv,kseqarg)
            lx=.true.
            return
          endif
        endif
      enddo
      list%attr=ior(iadv,knoseqarg)
      return
      end


      subroutine tfcomplexlistqkf(isp1,vx,irtc)
      use tfstk
      implicit none
      integer*4 irtc,isp1,itfmessage
      real*8 vx
      logical*4 tfcomplexlistqk
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(tfcomplexlistqk(dtastk(isp)))then
        vx=1.d0
      else
        vx=0.d0
      endif
      irtc=0
      return
      end

      recursive logical*4 function tfcomplexlistqk(k)
     $     result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      integer*8 kh
      integer*4 i
      lx=.false.
      if(ktflistqd(k,kl))then
        kh=kl%head
        if(iand(kl%attr,lnonreallist) .eq. 0)then
          if(kh .eq. ktfoper+mtfcomplex .and.
     $         kl%nl .eq. 2)then
            lx=.true.
          endif
        elseif(kh .eq. ktfoper+mtflist)then
          do i=1,kl%nl
            if(tfcomplexlistqk(kl%dbody(i)))then
              lx=.true.
              return
            endif
          enddo
        endif
      endif
      return
      end

      recursive logical*4 function tfsymbollistqo(kl) result(lx)
      use tfstk
      implicit none
      type (sad_list) kl
      type (sad_list), pointer :: kli
      type (sad_pat), pointer :: kpi
      integer*4 i
      if(iand(ksymbollist,kl%attr) .ne. 0)then
        lx=.true.
        return
      elseif(iand(knosymbollist,kl%attr) .ne. 0)then
        lx=.false.
        return
      endif
      lx=.false.
      do i=0,kl%nl
        if(ktfsymbolq(kl%body(i)))then
          lx=.true.
          kl%attr=ior(kl%attr,ksymbollist)
          return
        elseif(ktflistq(kl%body(i),kli))then
          lx=tfsymbollistqo(kli)
          if(lx)then
            kl%attr=ior(kl%attr,ksymbollist)
            return
          endif
        elseif(ktfpatq(kl%body(i),kpi))then
          if(kpi%sym%loc .ne. 0 .or. ktftype(kpi%expr%k) .ne. ktfref
     $         .or. ktftype(kpi%default%k) .ne. ktfref)then
            lx=.true.
            kl%attr=ior(kl%attr,ksymbollist)
            return
          endif
        endif
        if(i .eq. 0 .and. ktfreallistqo(kl))then
          exit
        endif
      enddo
      kl%attr=ior(kl%attr,knosymbollist)
      return
      end
        
      subroutine tfvectorqf(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kxi
      type (sad_list), pointer ::kl
      integer*4 isp1,irtc,m,i,ispf,narg,itfmessage
      narg=isp-isp1
      if(narg .eq. 1)then
        kx%k=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp),kl))then
          kx%k=0
        else
          m=kl%nl
          if(m .ne. 0)then
            if(ktfnonreallistqo(kl))then
              do i=1,m
                if(tflistqk(kl%body(i)))then
                  kx%k=0
                  return
                endif
              enddo
            endif
          endif
        endif
      elseif(narg .eq. 2)then
        kx%k=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp-1),kl))then
          kx%k=0
        else
          m=kl%nl
          if(m .ne. 0)then
            ispf=isp+1
            ktastk(ispf)=ktfcopy(ktastk(isp))
            do i=1,m
              isp=ispf+1
              ktastk(isp)=kl%body(i)
              if(tflistqk(ktastk(ispf+1)))then
                isp=ispf-1
                kx%k=0
                go to 100
              endif
              call tfefunref(ispf,kxi,.true.,irtc)
              isp=ispf-1
              if(irtc .ne. 0)then
                go to 100
              endif
              if(ktfnonrealqd(kxi) .or. kxi%k .ne. 0)then
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

      subroutine tfmatchqf(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfpmatc,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      irtc=0
      if(itfpmatc(ktastk(isp-1),ktastk(isp)) .ge. 0)then
        kx%k=ktftrue
      else
        kx%k=0
      endif
      return
      end

      logical function tfrepeatedqk(k,null)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list) , pointer :: kl
      logical*4 null
      tfrepeatedqk=.false.
      if(ktflistqd(k,kl))then
        if(kl%head .eq. ktfoper+mtfrepeated)then
          tfrepeatedqk=.true.
          null=.false.
        elseif(kl%head .eq. ktfoper+mtfrepeatednull)then
          tfrepeatedqk=.true.
          null=.true.
        endif
      endif
      return
      end

      logical function tfinequalityqk(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      tfinequalityqk=ktflistqd(k,kl) .and.
     $     kl%head .eq. ktfoper+mtfinequality
      return
      end

      logical*4 function tfcontextqk(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_symbol), pointer :: sym
      tfcontextqk=ktfsymbolqd(k,sym) .and. sym%gen .eq. -3
      return
      end

      recursive logical*4 function tfrefq(k) result(l)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer ::kl
      integer*4 i
      if(ktfrefq(k%k))then
        l=.true.
      else
        l=.false.
        if(tflistqd(k,kl))then
          if(kl%nl .gt. 0)then
            do i=1,kl%nl
              if(.not. tfrefq(kl%dbody(i)))then
                return
              endif
            enddo
            l=.true.
          endif
        endif
      endif
      return
      end

      logical*4 function tfruleqk(k)
      use tfstk, tf => tfruleqk
      integer*8 k
      tfruleqk=tf(k)
      return
      end

      subroutine tfsameqdummy
      use mackw
      return
      end
