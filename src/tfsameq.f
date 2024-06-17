      module sameq
      use tfstk

      contains
      recursive logical*4 function tfnearlysameqf(k1,k2,re,ae) result(lx)
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
          if(d <= s*re .or. d <= ae)then
            lx=.true.
          endif
        endif
      elseif(tfnumberq(k2))then
      elseif(k1%k /= k2%k)then
      elseif(ktfnonlistq(k1,kl1))then
        lx=tfsameq(k1,k2)
      elseif(ktfnonlistq(k2,kl2))then
      else
        n1=kl1%nl
        if(n1 /= kl2%nl)then
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
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,ktfrealq(ktastk(isp)))
      return
      end

      subroutine tfnanqk(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,
     $     ktfrealq(ktastk(isp)) .and. ktfenanq(rtastk(isp)))
      return
      end

      subroutine tfcomplexlistqkf(isp1,kx,irtc)
      implicit none
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      type (sad_descriptor) ,intent(out):: kx
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      kx%x(1)=merge(1.d0,0.d0,tfcomplexlistqk(dtastk(isp)))
      irtc=0
      return
      end

      recursive logical*4 function tfcomplexlistqk(k) result(lx)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      integer*8 kh
      integer*4 i
      lx=.false.
      if(ktflistq(k,kl))then
        kh=kl%head%k
        if(iand(kl%attr,lnonreallist) == 0)then
          if(kh == ktfoper+mtfcomplex .and. kl%nl == 2)then
            lx=.true.
          endif
        elseif(kh == ktfoper+mtflist)then
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
      implicit none
      type (sad_dlist) ,intent(inout):: kl
      type (sad_dlist), pointer :: kli
      type (sad_pat), pointer :: kpi
      integer*4 i
      if(iand(ksymbollist,kl%attr) /= 0)then
        lx=.true.
        return
      elseif(iand(knosymbollist,kl%attr) /= 0)then
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
          if(kpi%sym%loc /= 0 .or. ktftype(kpi%expr%k) /= ktfref
     $         .or. ktftype(kpi%default%k) /= ktfref)then
            lx=.true.
            kl%attr=ior(kl%attr,ksymbollist)
            return
          endif
        endif
        if(i == 0 .and. ktfreallistq(kl))then
          exit
        endif
      enddo
      kl%attr=ior(kl%attr,knosymbollist)
      return
      end

      subroutine tfmatchqf(isp1,kx,irtc)
      use pmat,only:itfpmatc
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      irtc=0
      kx%k=merge(ktftrue,ktffalse,itfpmatc(dtastk(isp-1),dtastk(isp)) >= 0)
      return
      end

      logical function tfrepeatedqk(k,null)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) , pointer :: kl
      logical*4 ,intent(out):: null
      tfrepeatedqk=.false.
      if(ktflistq(k,kl))then
        if(kl%head%k == ktfoper+mtfrepeated)then
          tfrepeatedqk=.true.
          null=.false.
        elseif(kl%head%k == ktfoper+mtfrepeatednull)then
          tfrepeatedqk=.true.
          null=.true.
        endif
      endif
      return
      end

      recursive logical*4 function tfrefq(k) result(l)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer ::kl
      integer*4 i
      if(ktfrefq(k))then
        l=.true.
      else
        l=.false.
        if(tflistq(k,kl) .and. kl%nl > 0)then
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

      subroutine tfsameqdummy
      use mackw
      return
      end

      end module

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
