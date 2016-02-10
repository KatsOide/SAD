      recursive logical*4 function tfsameqk(ka,kp) result(lx)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_pat), pointer :: pata,patp
      integer*8 kaa,kap,kta,i,nc,kaa1,kap1
      integer*8 ka,kp
      logical*4 tfsamelistqo,tfsamesymbolqo
      if(ka .eq. kp)then
        lx=.true.
        return
      endif
      lx=.false.
      kta=iand(ktfmask,ka)
      if(kta .ne. iand(ktfmask,kp))then
        return
      elseif(kta .eq. ktfoper)then
        return
      endif
      kaa=ktfaddr(ka)
      kap=ktfaddr(kp)
      if(kta .eq. ktfsymbol)then
        lx=tfsamesymbolqo(klist(kaa-3),klist(kap-3))
      elseif(kta .eq. ktfstring)then
        nc=ilist(1,kaa)
        if(nc .eq. ilist(1,kap))then
          do i=1,nc/4+1
            if(ilist(i,kaa+1) .ne. ilist(i,kap+1))then
              return
            endif
          enddo
          lx=.true.
        endif
      elseif(kta .eq. ktflist)then
        if(ilist(2,kaa-1) .eq. ilist(2,kap-1))then
          lx=tfsamelistqo(klist(kaa-3),klist(kap-3))
c          call tfdebugprint(ktflist+kaa,'tfsameqk',3)
c          call tfdebugprint(ktflist+kap,'===',3)
c          write(*,*)'==> ',lx
        endif
      elseif(kta .eq. ktfpat)then
        call pat_loc(kaa,pata)
        call pat_loc(kap,patp)
        if(.not. tfsameqk(pata%expr,patp%expr))then
          return
        endif
        if(.not. tfsameqk(pata%head,patp%head))then
          return
        endif
        kaa1=ktfaddr(pata%sym%alloc)
        kap1=ktfaddr(patp%sym%alloc)
        if(kaa1 .ne. 0 .and. kap1 .ne. 0)then
          lx=tfsamesymbolqo(klist(kaa1-3),klist(kap1-3))
        else
          lx=kaa1 .eq. 0 .and. kap1 .eq. 0
        endif
      endif
      return
      end

      logical*4 function tfsamestringqk(ka1,kp1)
      use tfstk
      implicit none
      integer*8 ka1,kp1,ka,kp
      logical*4 tfsamestringqo
      ka=ktfaddr(ka1)
      kp=ktfaddr(kp1)
      tfsamestringqk=ka .eq. kp .or.
     $     tfsamestringqo(klist(ka-3),klist(kp-3))
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

      logical*4 function tfsamesymbolqk(ka1,kp1)
      use tfstk
      implicit none
      integer*8 ka,kp,ka1,kp1
      logical*4 tfsamesymbolqo
      ka=ktfaddr(ka1)
      kp=ktfaddr(kp1)
      tfsamesymbolqk=ka .eq. kp .or.
     $     tfsamesymbolqo(klist(ka-3),klist(kp-3))
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
      integer*8 k1,k2,ka1,ka2
      integer*4 i,n1
      real*8 re,ae,s,d
      complex*16 cx1,cx2,tfgetnumber
      logical*4 tfnumberqk,tfsameqk
      lx=.false.
      if(tfnumberqk())then
        if(tfnumberqk(k2))then
          cx1=tfgetnumber(k1)
          cx2=tfgetnumber(k2)
          s=abs(cx1)+abs(cx2)
          d=abs(cx2-cx1)
          if(d .le. s*re .or. d .le. ae)then
            lx=.true.
          endif
        endif
      elseif(tfnumberqk(k2))then
      elseif(k1 .ne. k2)then
      elseif(ktfnonlistq(k1))then
        lx=tfsameqk(k1,k2)
      elseif(ktfnonlistq(k2))then
      else
        ka1=ktfaddr(k1)
        ka2=ktfaddr(k2)
        n1=ilist(2,ka1-1)
        if(n1 .ne. ilist(2,ka2-1))then
          return
        endif
        if(.not. tfnearlysameqf(klist(ka1),klist(ka2),re,ae))then
          return
        endif
        do i=1,n1
          if(.not. tfnearlysameqf(klist(ka1+i),klist(ka2+i),
     $         re,ae))then
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
      integer*8 kx
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      if(ktfrealq(ktastk(isp)))then
        kx=ktftrue
      else
        kx=0
      endif
      return
      end

      logical*4 function tfcomplexqk(k,n)
      use tfstk
      implicit none
      type (sad_list), pointer :: kl,kl1,kl2
      integer*8 k,k1,k2
      integer*4 n
      logical*4 tflistqk
      tfcomplexqk=.false.
      if(ktflistq(k))then
        call list_loc(ktfaddr(k),kl)
        if(kl%head .eq. ktfoper+mtfcomplex .and. kl%nl .eq. 2)then
          if(ktfreallistqo(kl))then
            tfcomplexqk=.true.
            n=0
          else
 1          k1=kl%body(1)
            k2=kl%body(2)
            if(tflistqk(k1) .and. tflistqk(k2))then
              call list_loc(k1,kl1)
              call list_loc(k2,kl2)
              if(kl1%nl .and. kl2%nl .and.
     $             ktfreallistqo(kl1) .and. ktfreallistqo(kl2))then
                n=kl1%nl
                tfcomplexqk=.true.
                n=0
              endif
            endif
          endif
        endif
      endif
      return
      end

      logical*4 function tflistqk(k)
      use tfstk
      implicit none
      integer*8 k
      tflistqk=ktflistq(k) .and. klist(ktfaddr(k)) .eq. ktfoper+mtflist
      return
      end

      logical*4 function tfnonlistqk(k)
      use tfstk
      implicit none
      integer*8 k
      tfnonlistqk=ktfnonlistq(k) .or.
     $     klist(ktfaddr(k)) .ne. ktfoper+mtflist
      return
      end

      logical*4 function tfreallistq(k)
      use tfstk
      implicit none
      integer*8 k,ka
      if(ktflistq(k))then
        ka=ktfaddr(k)
        tfreallistq=klist(ka) .eq. ktfoper+mtflist
     $       .and. ktfreallistq(ka)
      else
        tfreallistq=.false.
      endif
      return
      end

      logical*4 function tfnumlistqnk(k,n)
      use tfstk
      implicit none
      integer*8 k,ka
      integer*4 n 
      if(ktflistq(k))then
        ka=ktfaddr(k)
        tfnumlistqnk=klist(ka) .eq. ktfoper+mtflist
     $       .and. ktfreallistq(ka) .and. ilist(2,ka-1) .eq. n
      else
        tfnumlistqnk=.false.
      endif
      return
      end

      logical*4 function tfexprqk(k)
      use tfstk
      implicit none
      integer*8 k
      tfexprqk=ktflistq(k) .and. klist(ktfaddr(k)) .ne. ktfoper+mtflist
      return
      end

      logical*4 function tfsequenceqk(k)
      use tfstk
      implicit none
      integer*8 k
      if(ktftype(k) .eq. ktflist)then
        tfsequenceqk=klist(ktfaddr(k)) .eq. ktfoper+mtfnull
      else
        tfsequenceqk=.false.
      endif
      return
      end

      logical*4 function tfsameheadqk(k1,k2)
      use tfstk
      implicit none
      integer*8 k1,k2
      logical*4 tfsameqk
      tfsameheadqk=tfsameqk(klist(ktfaddr(k1)),
     $     klist(ktfaddr(k2)))
      return
      end

      recursive logical*4 function tfruleqk(k) result(lx)
      use tfstk
      implicit none
      type (sad_list), pointer :: list
      integer*8 k
      integer*4 i
      lx=.false.
      if(ktflistq(k))then
        call list_loc(ktfaddr(k),list)
        if(list%head .eq. ktfoper+mtflist)then
          if(list%nl .gt. 0)then
            if(iand(list%attr,lnonreallist) .ne. 0)then
              do i=1,list%nl
                if(.not. tfruleqk(list%body(i)))then
                  return
                endif
              enddo
              lx=.true.
            endif
          endif
        elseif(list%head .eq. ktfoper+mtfrule .or.
     $         list%head .eq. ktfoper+mtfruledelayed)then
          lx=list%nl .eq. 2
        endif
      endif
      return
      end

      recursive logical*4 function tfconstqk(k) result(lx)
      use tfstk
      implicit none
      integer*8 k,kt,ka
      logical*4 tfconstlistqo
      lx=.true.
      ka=ktfaddr(k)
      kt=k-ka
      if(kt .eq. ktfsymbol)then
        lx=ktfconstantq(ka) .and. klist(ka-4) .eq. ktfsymbol+ka
     $       .and. klist(ka-6) .eq. 0
      elseif(kt .eq. ktflist)then
        lx=tfconstlistqo(klist(ka-3))
      elseif(kt .eq. ktfpat)then
        lx=tfconstqk(klist(ka))
      endif
      return
      end

      logical*4 function tfconstlistqo(list)
      use tfstk
      implicit none
      type (sad_list) list
      integer*8 i,kh
      logical*4 tfconstqk,tfconstheadqk,tfseqqo,nr
      if(iand(lnoconstlist,list%attr) .ne. 0)then
        tfconstlistqo=.false.
        return
      endif
      tfconstlistqo=.true.
      if(iand(lconstlist,list%attr) .ne. 0)then
        return
      endif
      kh=list%head
      nr=ktfnonreallistqo(list)
      if(kh .eq. ktfoper+mtfhold .or.
     $     kh .eq. ktfoper+mtffun)then
        tfconstlistqo=.not. tfseqqo(list)
      elseif(kh .eq. ktfoper+mtfcomplex)then
        tfconstlistqo=.not. nr
      elseif(.not. tfconstheadqk(kh))then
        tfconstlistqo=.false.
      elseif(nr)then
        do i=1,list%nl
          if(.not. tfconstqk(list%body(i)))then
            tfconstlistqo=.false.
            list%attr=ior(lnoconstlist,list%attr)
            return
          endif
        enddo
      endif
      if(tfconstlistqo)then
c        call tfdebugprint(ktflist+ka,'constlist',3)
        list%attr=ior(lconstlist,list%attr)
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
      type (sad_list), pointer :: list
      integer*8 k,kt,ka
      logical*4 tfconstlistqo,tfconstqk
      tfconstheadqk=.true.
      ka=ktfaddr(k)
      kt=k-ka
      if(kt .eq. ktfoper)then
        if(ka .le. mtfend)then
          tfconstheadqk=constop(ka)
        else
          tfconstheadqk=.false.
        endif
      elseif(kt .eq. ktfstring)then
        tfconstheadqk=.false.
      elseif(kt .eq. ktfsymbol)then
        tfconstheadqk=ktfconstantq(ka) .and.
     $       klist(ka-4) .eq. ktfsymbol+ka
     $       .and. klist(ka-5) .eq. 0
      elseif(kt .eq. ktflist)then
        call list_loc(ka,list)
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
      elseif(kt .eq. ktfpat)then
        tfconstheadqk=tfconstqk(klist(ka))
      endif
      return
      end

      logical*4 function tfheldqk(k)
      use tfstk
      implicit none
      integer*8 k
      tfheldqk=ktflistq(k) .and.
     $     klist(ktfaddr(k)) .eq. ktfoper+mtfhold
      return
      end

      recursive logical*4 function tfseqqo(list) result(lx)
      use tfstk
      use tfcode
      implicit none
      type (sad_list) list
      type (sad_list), pointer :: listi
      integer*8 i,kai
      integer*4 iadv
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
        if(ktftype(list%body(i)) .eq. ktflist)then
          kai=ktfaddr(list%body(i))
          if(klist(kai) .eq. ktfoper+mtfnull)then
            list%attr=ior(iadv,kseqarg)
            lx=.true.
            return
          endif
          call list_loc(kai,listi)
          if(tfconstlistqo(listi))then
          elseif(tfseqqo(listi))then
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
      if(tfcomplexlistqk(ktastk(isp)))then
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
      integer*8 k,ka,kh
      integer*4 i
      lx=.false.
      if(ktflistq(k))then
        ka=ktfaddr(k)
        kh=klist(ka)
        if(iand(ilist(2,ka-3),lnonreallist) .eq. 0)then
          if(kh .eq. ktfoper+mtfcomplex .and.
     $         ilist(2,ka-1) .eq. 2)then
            lx=.true.
          endif
        elseif(kh .eq. ktfoper+mtflist)then
          do i=1,ilist(2,ka-1)
            if(.not. tfcomplexlistqk(klist(ka+i)))then
              return
            endif
          enddo
          lx=.true.
        endif
      endif
      return
      end

      logical*4 function tfcomplexnumlistqk(k)
      use tfstk
      implicit none
      integer*8 k,ka
      integer*4 i
      logical*4 tfnumberqk,tflistqk
      tfcomplexnumlistqk=.false.
      if(tflistqk(k))then
        ka=ktfaddr(k)
        if(iand(ilist(2,ka-3),lnoconstlist) .eq. 0)then
          tfcomplexnumlistqk=.true.
        else
          do i=1,ilist(2,ka-1)
            if(.not. tfnumberqk(klist(ka+i)))then
              return
            endif
          enddo
          tfcomplexnumlistqk=.true.
        endif
      endif
      return
      end
        
      subroutine tfvectorqf(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,ktfcopy,ka,kxi
      integer*4 isp1,irtc,m,i,ispf,narg,itfmessage
      logical*4 tflistqk
      narg=isp-isp1
      if(narg .eq. 1)then
        kx=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp)))then
          kx=0
        else
          ka=ktfaddr(ktastk(isp))
          m=ilist(2,ka-1)
          if(m .ne. 0)then
            if(ktfnonreallistq(ka))then
              do i=1,m
                if(tflistqk(klist(ka+i)))then
                  kx=0
                  return
                endif
              enddo
            endif
          endif
        endif
      elseif(narg .eq. 2)then
        kx=ktftrue
        irtc=0
        if(ktfnonlistq(ktastk(isp-1)))then
          kx=0
        else
          ka=ktfaddr(ktastk(isp-1))
          m=ilist(2,ka-1)
          if(m .ne. 0)then
            ispf=isp+1
            ktastk(ispf)=ktfcopy(ktastk(isp))
            do i=1,m
              isp=ispf+1
              ktastk(isp)=klist(ka+i)
              if(tflistqk(ktastk(ispf+1)))then
                isp=ispf-1
                kx=0
                go to 100
              endif
              call tfefunref(ispf,kxi,.true.,irtc)
              isp=ispf-1
              if(irtc .ne. 0)then
                go to 100
              endif
              if(ktfnonrealq(kxi) .or. kxi .ne. 0)then
                kx=0
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
      integer*8 kx
      integer*4 isp1,irtc,itfpmatc,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      irtc=0
      if(itfpmatc(ktastk(isp-1),ktastk(isp)) .ge. 0)then
        kx=ktftrue
      else
        kx=0
      endif
      return
      end

      logical function tfrepeatedqk(k,null)
      use tfstk
      implicit none
      integer*8 k,ka
      logical*4 null
      tfrepeatedqk=.false.
      if(iand(ktfmask,k) .eq. ktflist)then
        ka=ktfaddr(k)
        if(klist(ka) .eq. ktfoper+mtfrepeated)then
          tfrepeatedqk=.true.
          null=.false.
        elseif(klist(ka).eq. ktfoper+mtfrepeatednull)then
          tfrepeatedqk=.true.
          null=.true.
        endif
      endif
      return
      end

      logical function tfinequalityqk(k)
      use tfstk
      implicit none
      integer*8 k
      tfinequalityqk=ktflistq(k) .and.
     $     klist(ktfaddr(k)) .eq. ktfoper+mtfinequality
      return
      end

      logical*4 function tfcontextqk(k)
      use tfstk
      implicit none
      integer*8 k
      tfcontextqk=ktfsymbolq(k) .and. ilist(2,ktfaddr(k)-1) .eq. -3
      return
      end

      logical*4 function tfrefq(k)
      use tfstk
      implicit none
      integer*8 k,ka
      integer*4 i
      logical*4 tflistqk
      if(ktfrefq(k))then
        tfrefq=.true.
      else
        tfrefq=.false.
        if(tflistqk(k))then
          ka=ktfaddr(k)
          if(ilist(2,ka-1) .gt. 0)then
            do i=1,ilist(2,ka-1)
              if(ktfnonrefq(klist(ka+i)))then
                return
              endif
            enddo
            tfrefq=.true.
          endif
        endif
      endif
      return
      end
