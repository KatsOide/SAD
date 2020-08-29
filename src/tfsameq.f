      recursive logical*4 function tfnearlysameqf(k1,k2,re,ae)
     $     result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_dlist), pointer :: kl1,kl2
      integer*4 i,n1
      real*8 ,intent(in):: re,ae
      real*8 s,d
      complex*16 cx1,cx2
      lx=.false.
      if(tfnumberq(k1,cx1))then
        if(tfnumberq(k2,cx2))then
          s=abs(cx1)+abs(cx2)
          d=abs(cx2-cx1)
          if(d .le. s*re .or. d .le. ae)then
            lx=.true.
          endif
        endif
      elseif(tfnumberq(k2))then
      elseif(k1%k .ne. k2%k)then
      elseif(ktfnonlistq(k1,kl1))then
        lx=tfsameq(k1,k2)
      elseif(ktfnonlistq(k2,kl2))then
      else
        n1=kl1%nl
        if(n1 .ne. kl2%nl)then
          return
        endif
        if(.not. tfnearlysameqf(kl1%head,kl2%head,re,ae))then
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
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,ktfrealq(ktastk(isp)))
      return
      end

      subroutine tfnanqk(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,
     $     ktfrealq(ktastk(isp)) .and. ktfenanq(rtastk(isp)))
      return
      end

      logical*4 function tfconstlistqo(list)
      use tfstk
      implicit none
      type (sad_dlist) ,intent(inout):: list
      type (sad_descriptor) kh
      integer*4 i
      logical*4 tfconstheadqk,tfseqqo,nr
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
            if(.not. tfconstq(list%dbody(i)%k))then
              tfconstlistqo=.false.
              list%attr=ior(lnoconstlist,list%attr)
              return
            endif
          enddo
          list%attr=ior(list%attr,kconstarg)
        endif
      endif
      list%attr=ior(merge(lconstlist+kconstarg,lnoconstlist,
     $     tfconstlistqo),list%attr)
      return
      end

      logical*4 function tfconstheadqk(k)
      use ophash
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: list
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symd
      integer*8 ka
      logical*4 tfconstlistqo
      tfconstheadqk=.true.
      if(ktfoperq(k,ka))then
        tfconstheadqk=merge(constop(ka),.false.,ka .le. mtfend)
      elseif(ktfstringq(k))then
        tfconstheadqk=.false.
      elseif(ktfsymbolqdef(k%k,symd))then
        tfconstheadqk=ktfconstantsymq(symd%sym) .and.
     $       symd%value%k .eq. ktfsymbol+ka
     $       .and. symd%downval .eq. 0
      elseif(ktflistq(k,list))then
        if(list%head%k .eq. ktfoper+mtflist .or.
     $       list%head%k .eq. ktfoper+mtffun)then
          tfconstheadqk=.false.
        elseif(iand(lconstlist,list%attr) .eq. 0)then
          if(iand(lnoconstlist,list%attr) .ne. 0)then
            tfconstheadqk=.false.
          elseif(.not. tfconstlistqo(list))then
            tfconstheadqk=.false.
          endif
        endif
      elseif(ktfpatq(k,pat))then
        tfconstheadqk=tfconstq(pat%expr%k)
      endif
      return
      end

      recursive logical*4 function tfseqqo(list) result(lx)
      use tfstk
      use tfcode
      implicit none
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist), pointer :: kli
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
        if(ktflistq(list%dbody(i),kli))then
          if(kli%head%k .eq. ktfoper+mtfnull)then
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
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      real*8 ,intent(out):: vx
      logical*4 tfcomplexlistqk
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      vx=merge(1.d0,0.d0,tfcomplexlistqk(dtastk(isp)))
      irtc=0
      return
      end

      recursive logical*4 function tfcomplexlistqk(k)
     $     result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      integer*8 kh
      integer*4 i
      lx=.false.
      if(ktflistq(k,kl))then
        kh=kl%head%k
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
      type (sad_dlist) ,intent(inout):: kl
      type (sad_dlist), pointer :: kli
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
        if(ktfsymbolq(kl%dbody(i)))then
          lx=.true.
          kl%attr=ior(kl%attr,ksymbollist)
          return
        elseif(ktflistq(kl%dbody(i),kli))then
          lx=tfsymbollistqo(kli)
          if(lx)then
            kl%attr=ior(kl%attr,ksymbollist)
            return
          endif
        elseif(ktfpatq(kl%dbody(i),kpi))then
          if(kpi%sym%loc .ne. 0 .or. ktftype(kpi%expr%k) .ne. ktfref
     $         .or. ktftype(kpi%default%k) .ne. ktfref)then
            lx=.true.
            kl%attr=ior(kl%attr,ksymbollist)
            return
          endif
        endif
        if(i .eq. 0 .and. ktfreallistq(kl))then
          exit
        endif
      enddo
      kl%attr=ior(kl%attr,knosymbollist)
      return
      end
        
      subroutine tfvectorqf(isp1,kx,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kxi
      type (sad_dlist), pointer ::kl
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 m,i,ispf,narg,itfmessage
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
                if(tflistq(kl%dbody(i)))then
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
              dtastk(isp)=kl%dbody(i)
              if(tflistq(ktastk(ispf+1)))then
                isp=ispf-1
                kx%k=0
                go to 100
              endif
              kxi=tfefunref(ispf,.true.,irtc)
              isp=ispf-1
              if(irtc .ne. 0)then
                go to 100
              endif
              if(ktfnonrealq(kxi) .or. kxi%k .eq. 0)then
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
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfpmatc,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,
     $     itfpmatc(ktastk(isp-1),ktastk(isp)) .ge. 0)
      return
      end

      logical function tfrepeatedqk(k,null)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) , pointer :: kl
      logical*4 ,intent(out):: null
      tfrepeatedqk=.false.
      if(ktflistq(k,kl))then
        if(kl%head%k .eq. ktfoper+mtfrepeated)then
          tfrepeatedqk=.true.
          null=.false.
        elseif(kl%head%k .eq. ktfoper+mtfrepeatednull)then
          tfrepeatedqk=.true.
          null=.true.
        endif
      endif
      return
      end

      logical*4 function tfcontextqk(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_symbol), pointer :: sym
      tfcontextqk=ktfsymbolq(k,sym) .and. sym%gen .eq. -3
      return
      end

      recursive logical*4 function tfrefq(k) result(l)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer ::kl
      integer*4 i
      if(ktfrefq(k))then
        l=.true.
      else
        l=.false.
        if(tflistq(k,kl) .and. kl%nl .gt. 0)then
          do i=1,kl%nl
            if(.not. tfrefq(kl%dbody(i)))then
              return
            endif
          enddo
          l=.true.
        endif
      endif
      return
      end

      logical*4 function tfsamesymbolqk(k1,k2)
      use tfstk, tf => tfsamesymbolqk
      integer*8 ,intent(in):: k1,k2
      tfsamesymbolqk=tf(k1,k2)
      return
      end

      logical*4 function tfruleqk(k)
      use tfstk, tf => tfruleqk_dlist
      integer*8 ,intent(in):: k
      tfruleqk=tf(k)
      return
      end

      subroutine tfsameqdummy
      use mackw
      return
      end
