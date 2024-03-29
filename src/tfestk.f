      subroutine tfestk(isp0,iprior,lastfirst,irtc)
      use tfstk
      use ophash
      use eexpr
      implicit none
      type (sad_descriptor) kx,kh,tfeval1
      type (sad_dlist), pointer :: klx
      integer*4 ,intent(in):: isp0,iprior(0:mtfnopc)
      integer*4 ,intent(out):: irtc
      integer*4 iop,iop1,isp1,i,itgetfpe,itfmessage
      logical*4 ,intent(in):: lastfirst(0:mtfnopc)
      logical*4 ineq
      ineq(iop1)=iop1 >= mtfgreater .and. iop1 <= mtfunequal .or.
     $     iop1 == mtfsame .or. iop1 == mtfunsame
      irtc=0
      do while(isp > isp0)
        iop=itastk2(1,isp)
c        call tfdebugprint(dtastk(isp),'estk',1)
        if(iop == mtfnull .or. iop == mtfleftparen .or. iop == mtflist)then
          return
        endif
        isp1=isp-1
        iop1=itastk2(1,isp1)
c        write(*,*)'with ',isp,isp0,iop,iop1
        if(iop == mtfrightbra)then
          if(iop1 == mtfrightbra)then
            if(ktastk(isp) == ktfoper+mtfnull)then
              isp=isp-1
              do i=isp,isp0,-1
                if(itastk2(1,i) == mtfpart)then
                  isp1=i
                  kx=kxmakelist(isp1-1,klx)
                  klx%head%k=ktfoper+mtfpart
                  iop=mtfnull
                  go to 1010
                endif
              enddo
            endif
            go to 9010
          elseif(iop1 == mtfcomma .or. iop1 == mtfleftbra)then
            do i=isp1,isp0,-1
              if(itastk2(1,i) == mtfleftparen
     $             .or. itastk2(1,i) == mtflist)then
                go to 9010
              elseif(itastk2(1,i) == mtfpart)then
                return
              elseif(itastk2(1,i) == mtfleftbra)then
                isp1=i
                if(isp == isp1+1 .and. ktastk(isp) == ktfoper+mtfnull)then
                  isp=isp1
                endif
                kh=dtastk(i)
                call tfcomposefull(i,kh,kx,irtc)
                if(irtc /= 0)then
                  return
                endif
                iop=mtfnull
                go to 1010
              endif
            enddo
            go to 9010
          elseif(iop1 == mtfpart)then
            return
          endif
        elseif(iop == mtfrightbrace .and. (iop1 == mtfcomma .or. iop1 == mtflist))then
          do i=isp1,isp0,-1
            if(itastk2(1,i) == mtfleftparen
     $           .or. itastk2(1,i) == mtfleftbra
     $           .or. itastk2(1,i) == mtfpart)then
              go to 9010
            elseif(itastk2(1,i) == mtflist)then
              isp1=i
              kx=kxmakelist(isp1)
              iop=mtfnull
              go to 1010
            endif
          enddo
          go to 9010
        elseif(iop == mtfrightparen .and. iop1 == mtfleftparen)then
          ktastk(isp1)=ktastk(isp)
          itastk2(1,isp1)=mtfnull
          isp=isp1
          return
        elseif(iop1 == mtffun .and. ktastk(isp) == ktfoper+mtfnull)then
          kx=kxpfaloc(dtastk(isp-1))
          go to 1010
        endif
        if(iprior(iop) < iprior(iop1))then
          return
        elseif(iprior(iop) == iprior(iop1) .and. lastfirst(iop))then
          return
        endif
        if(iop1 == mtfcomma)then
          if(iop /= mtfcomma)then
            go to 9020
          endif
          return
        elseif(iop1 == mtfpower)then
          if(ktastk(isp) == ktfoper+mtfnull)then
            return
          endif
        elseif(iop1 == mtfleftbra .or. iop1 == mtfpart)then
          return
        endif
        if(ineq(iop1))then
          call tfeinequal(dtastk(isp1),dtastk(isp),kx,iop1,ineq(iop))
        else
          if(iop1 == mtfneg)then
            iop1=mtftimes
          elseif(iop1 == mtfnot)then
            if(ktastk(isp1) /= ktfoper+mtfnull)then
              go to 9030
            endif
          endif
          if(tfconstq(ktastk(isp1)) .and.
     $         .not. tfheldqd(dtastk(isp1)) .and.
     $         tfconstq(ktastk(isp)) .and.
     $         .not. tfheldqd(dtastk(isp)) .and.
     $         ktfimmediateq(klist(ifunbase+iop1)))then
            kx=tfeval1(dtastk(isp1),dtastk(isp),iop1,irtc)
            if(irtc /= 0)then
              return
            endif
            if(itgetfpe() > 0)then
              irtc=itfmessage(9,'General::fpe','""')
              return
            endif
          else
c            call tfmemcheckprint1('estk-eexpr',1,.true.)
c            call tfdebugprint(dtastk(isp1),'estk-eexpr',1)
c            write(*,*)iop1
c            call tfdebugprint(dtastk(isp),':',1)
            kx=tfeexpr(dtastk(isp1),dtastk(isp),iop1)
c            call tfdebugprint(kx,'estk-eexpr-end',1)
          endif
        endif
 1010   isp=isp1
        itastk2(1,isp)=iop
        dtastk(isp)=kx
      enddo
      return
 9010 irtc=itfmessage(9999,'General::mismatch',
     $     '"'//opcode(iop)//'"')
      return
 9020 irtc=itfmessage(9999,'General::incomplete','""')
      return
 9030 irtc=itfmessage(9,'General::narg','"1 for ~"')
      return
      end

      subroutine tfeinequal(k1,k2,kx,iop1,nextrel)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: kl1,kl2
      type (sad_dlist), pointer :: klx
      integer*4 ,intent(in):: iop1
      integer*4 i,m1,m2,irtc
      logical*4 ,intent(in):: nextrel
c      call tfdebugprint(k1,'einequal',1)
c      call tfdebugprint(k2,'and',1)
c      write(*,*)'with ',iop1,nextrel
      if(tfinequalityq(k1))then
        call loc_sad(ktfaddrd(k1),kl1)
        m1=kl1%nl
        if(tfinequalityq(k2))then
          call loc_sad(ktfaddrd(k2),kl2)
          m2=kl2%nl
          kx=kxadaloc(-1,m1+m2+1,klx)
          klx%head%k=ktfoper+mtfinequality
          do i=1,m1
            klx%dbody(i)=dtfcopy(kl1%dbody(i))
          enddo
          klx%dbody(m1+1)%k=ktfoper+iop1
          do i=1,m2
            klx%dbody(i+m1+1)=dtfcopy(kl2%dbody(i))
          enddo
        else
          kx=kxadaloc(-1,m1+2,klx)
          klx%head%k=ktfoper+mtfinequality
          do i=1,m1
            klx%dbody(i)=dtfcopy(kl1%dbody(i))
          enddo
          klx%dbody(m1+1)%k=ktfoper+iop1
          klx%dbody(m1+2)=dtfcopy(k2)
        endif
      else
        if(tfinequalityq(k2))then
          call loc_sad(ktfaddrd(k2),kl2)
          m2=kl2%nl
          kx=kxadaloc(-1,m2+2,klx)
          klx%head%k=ktfoper+mtfinequality
          klx%dbody(1)=dtfcopy(k1)
          klx%dbody(2)%k=ktfoper+iop1
          do i=1,m2
            klx%dbody(i+2)=dtfcopy(kl2%dbody(i))
          enddo
        else
          if(nextrel)then
            kx=kxadaloc(-1,3,klx)
            klx%head%k=ktfoper+mtfinequality
            klx%dbody(1)=dtfcopy(k1)
            klx%dbody(2)%k=ktfoper+iop1
            klx%dbody(3)=dtfcopy(k2)
          else
            kx=kxadaloc(-1,2,klx)
            klx%head%k=ktfoper+iop1
            klx%dbody(1)=dtfcopy(k1)
            klx%dbody(2)=dtfcopy(k2)
          endif
        endif
      endif
      if(.not. nextrel)then
        if(tfconstlistqo(klx))then
          kx=tfleval(klx,.true.,irtc)
        endif
      endif
      return
      end

      subroutine tfinequality(isp1,kx,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) tfeval1,kxi
      integer*8 ka
      integer*4 i,narg,itfmessage,iop1
      logical*4 ineq
      ineq(iop1)=iop1 >= mtfgreater .and. iop1 <= mtfunequal .or.
     $     iop1 == mtfsame .or. iop1 == mtfunsame
      narg=isp-isp1
      if(narg < 3 .or. narg /= (narg/2)*2+1)then
        irtc=itfmessage(9,'General::narg','"3 or more"')
        return
      endif
      kx%k=ktftrue
      do i=isp1+1,isp,2
        dtastk(i)=tfeevalref(dtastk(i),irtc)
        if(irtc /= 0)then
          return
        endif
        if(i > isp1+1)then
          if(ktfoperq(dtastk(i-1)) .and.
     $         ineq(int(ktfaddr(ktastk(i-1)))))then
            ka=ktfaddr(ktastk(i-1))
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"==, <>, <, >, <=, >=, ===, <=>"')
            return
          endif
          kxi=tfeval1(dtastk(i-2),dtastk(i),int(ka),irtc)
          if(irtc /= 0)then
            return
          endif
          if(kxi%k == 0)then
            kx=kxi
            return
          elseif(ktfnonrealq(kxi))then
            kx=tfeval1(kx,kxi,mtfand,irtc)
            if(irtc /= 0)then
              return
            endif
          endif
        endif
      enddo
      return
      end

      function tfcompose(isp1,kh,irtc) result (kx)
      use tfstk
      implicit none
c      include 'DEBUG.inc'
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(in):: kh
      type (sad_descriptor) tfcomposefun,tfcomposeoper
      integer*8 kah
      integer*4 iah,isp0
      if(rlist(iaximmediate) /= 0.d0)then
        if(ktfoperq(kh,kah))then
          iah=int(kah)
          if(iah > mtfend)then
            kx=tfcomposefun(isp1,iah,.false.,irtc)
          else
            kx=tfcomposeoper(isp1,iah,.true.,isp0,irtc)
          endif
          if(irtc == 0)then
            return
          endif
        elseif(ktfstringq(kh))then
          if(isp == isp1+1 .or. isp == isp1+2)then
            if(ktfrealq(dtastk(isp1+1)) .and.
     $           ktfrealq(dtastk(isp)))then
              kx=kxsubstring(kh,isp1+1,isp)
              irtc=0
              return
            endif
          endif
        endif
      endif
      kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),kh)
      irtc=0
      return
      end

      subroutine tfcomposefull(isp1,kh,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(in):: kh
      type (sad_descriptor) tfcomposefun,tfcomposeoper,tfefunrefu
      integer*8 kah
      integer*4 i,isp0,iah
      if(isp /= isp1 .and. rlist(iaximmediate) /= 0.d0)then
        if(ktfoperq(kh,kah))then
          iah=int(kah)
          if(iah > mtfend)then
            kx=tfcomposefun(isp1,iah,.true.,irtc)
            if(irtc == 0)then
              return
            elseif(irtc > 0)then
              go to 10
            endif
            if(.not. ktfimmediateq(klist(ifunbase+iah)))then
              go to 10
            endif
          else
            kx=tfcomposeoper(isp1,iah,.true.,isp0,irtc)
            if(irtc == 0)then
              return
            elseif(irtc > 0)then
              if(iah == mtfcomplex)then
                go to 1
              endif
              go to 10
            endif
          endif
        elseif(ktfstringq(kh))then
          if(isp == isp1+1 .or. isp == isp1+2)then
            if(ktfrealq(ktastk(isp1+1)) .and. ktfrealq(ktastk(isp)))then
              kx=kxsubstring(kh,isp1+1,isp)
              irtc=0
              return
            endif
          endif
          go to 10
        else
          go to 10
        endif
 1      isp0=isp
        isp=isp+1
        dtastk(isp)=kh
        do i=isp1+1,isp0
          if(tfconstq(ktastk(i)))then
            if(ktastk(i) == ktfoper+mtfnull)then
              isp=isp0
              go to 10
            elseif(tfheldqd(dtastk(i)))then
              isp=isp0
              go to 10
            endif
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            isp=isp0
            go to 10
          endif
        enddo
        kx=tfefunrefu(isp0+1,irtc)
        isp=isp0
        return
      endif
 10   kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),kh)
      irtc=0
      return
      end

      function  tfcomposeoper(isp1,iah,comp,isp0,irtc) result(kx)
      use tfstk
      use tfcx,only:tfclassmember
      use convstr,only:tfstringjoin
      use ophash
      use eexpr
      implicit none
      type (sad_descriptor) kx,tfpart
      integer*4 ,intent(in):: isp1,iah
      integer*4 ,intent(out):: irtc,isp0
      type (sad_dlist), pointer :: klx
      integer*8 iee0
      integer*4 i,isp00,iv
      logical*4 ,intent(in):: comp
      real*8 vx
      irtc=1
      kx=dxnullo
      if(constop(iah))then
        return
      elseif(isp1+2 == isp)then
        select case (iah)
          case (mtfset,mtfsetdelayed,mtfmap,mtfapply,
     $         mtfmapall,mtfunset,mtfmessagename)
          return
        case (mtfplus,mtftimes)
          if(tfconstq(ktastk(isp1+1)))then
            if(tfconstq(ktastk(isp)))then
              kx=tfplus(isp1,iah,irtc)
              return
            endif
          endif
        case (mtfpower,mtfrevpower)
          if(tfconstq(ktastk(isp1+1)))then
            if(tfconstq(ktastk(isp)))then
              kx=tfecmplxl(dtastk(isp1+1),dtastk(isp),mtfpower)
              return
            endif
          endif
        case (mtfequal)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktftrue,ktffalse,rtastk(isp1+1) == rtastk(isp))
            irtc=0
            return
          elseif(ktftype(ktastk(isp1+1)) == ktfstring .and.
     $           ktftype(ktastk(isp)) == ktfstring)then
            kx%k=merge(ktftrue,ktffalse,
     $           tfsamestringq(ktastk(isp1+1),ktastk(isp)))
            irtc=0
            return
          endif
        case (mtfunequal)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktffalse,ktftrue,rtastk(isp1+1) == rtastk(isp))
            irtc=0
            return
          elseif(ktftype(ktastk(isp1+1)) == ktfstring .and.
     $           ktftype(ktastk(isp)) == ktfstring)then
            kx%k=merge(ktffalse,ktftrue,
     $           tfsamestringq(ktastk(isp1+1),ktastk(isp)))
            irtc=0
            return
          endif
        case (mtfsame)
          if(tfconstq(ktastk(isp1+1)) .and. tfconstq(ktastk(isp)))then
            kx%k=merge(ktftrue,ktffalse,
     $           tfsameq(ktastk(isp1+1),ktastk(isp)))
            irtc=0
            return
          endif
        case (mtfunsame)
          if(tfconstq(ktastk(isp1+1)) .and. tfconstq(ktastk(isp)))then
            kx%k=merge(ktffalse,ktftrue,
     $           tfsameq(ktastk(isp1+1),ktastk(isp)))
            irtc=0
            return
          endif
        case (mtfgreater)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktftrue,ktffalse,
     $           rtastk(isp1+1) > rtastk(isp))
            irtc=0
            return
          endif
        case (mtfgeq)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktftrue,ktffalse,
     $           rtastk(isp1+1) >= rtastk(isp))
            irtc=0
            return
          endif
        case(mtfless)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktftrue,ktffalse,
     $           rtastk(isp1+1) < rtastk(isp))
            irtc=0
            return
          endif
        case (mtfleq)
          if(iand(ktrmask,ktastk(isp1+1)) /= ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) /= ktfnr )then
            kx%k=merge(ktftrue,ktffalse,
     $           rtastk(isp1+1) <= rtastk(isp))
            irtc=0
            return
          endif
        case (mtfpart)
          if(ktastk(isp1+1) == klist(iaxslotnull))then
            if(ktfrealq(ktastk(isp)))then
              iv=int(rtastk(isp))
              if(iv >= -2 .and. iv <= nslots)then
                kx=dlist(iaxslotpart+iv+2)
                irtc=0
                return
              endif
            endif
            return
          endif
        case (mtfatt)
          call tfclassmember(dtastk(isp1+1),dtastk(isp),kx,.false.,irtc)
          return
        case (mtffun)
          if(ktastk(isp) == ktfoper+mtfnull)then
            kx=kxpfaloc(dtastk(isp-1))
            irtc=0
          endif
          return
        end select
      endif
      select case (iah)
      case (mtfpart)
        if(tfconstq(ktastk(isp1+1)))then
          if(tflistq(ktastk(isp1+1)))then
            do i=isp1+2,isp
              if(ktfnonrealq(ktastk(i)))then
                if(ktastk(i) /= ktfoper+mtfnull)then
                  return
                endif
              endif
            enddo
            isp0=isp
            iee0=ierrorexp
            ierrorexp=1
            kx=tfpart(isp1+1,.false.,irtc)
            ierrorexp=iee0
            isp=isp0
            if(irtc /= 0)then
              irtc=1
              return
            endif
          endif
        endif
        return
      case (mtfplus)
        vx=0.d0
        isp0=isp
        isp=isp+1
        do i=isp1+1,isp0
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            vx=vx+rtastk(i)
          endif
        enddo
        if(isp == isp0+1)then
          kx=dfromr(vx)
          irtc=0
          return
        elseif(abs(vx) == 0.d0)then
          if(isp-isp0-1 /= isp0-isp1)then
            isp=isp0
            return
          endif
          isp00=isp0
          isp0=isp0+1
          go to 110
        else
          rtastk(isp0+1)=vx
          isp00=isp0
          go to 110
        endif
      case (mtftimes)
        vx=1.d0
        isp0=isp
        isp=isp+1
        do i=isp1+1,isp0
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            vx=vx*rtastk(i)
          endif
        enddo
        if(isp == isp0+1)then
          kx=dfromr(vx)
          isp=isp0
          irtc=0
          return
        elseif(vx == 1.d0)then
          if(isp-isp0-1 /= isp0-isp1)then
            isp=isp0
            return
          endif
          isp00=isp0
          isp0=isp0+1
          go to 110
        else
          rtastk(isp0+1)=vx
          isp00=isp0
          go to 110
        endif
      case (mtfand)
        isp0=isp
        do i=isp1+1,isp
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          elseif(abs(rtastk(i)) == 0.d0)then
            kx%k=0
            isp=isp0
            irtc=0
            return
          endif
        enddo
        if(isp == isp0)then
          kx%k=ktftrue
          irtc=0
          return
        else
          isp00=isp0
          go to 110
        endif
      case (mtfor)
        isp0=isp
        do i=isp1+1,isp
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          elseif(abs(rtastk(i)) /= 0.d0)then
            kx%k=ktftrue
            isp=isp0
            irtc=0
            return
          endif
        enddo
        if(isp == isp0)then
          kx%k=0
          irtc=0
          return
        else
          isp00=isp0
          go to 110
        endif
      case (mtfconcat)
        do i=isp1+1,isp
          if(ktfnonstringq(ktastk(i)))then
            go to 1
          endif
        enddo
        call tfstringjoin(isp1,kx,irtc)
        return
      end select
 1    irtc=-1
      return
 110  irtc=0
      if(isp == isp0+1)then
        kx=dtastk(isp)
        if(ktfrealq(kx))then
          isp=isp00
          return
        elseif(ktflistq(kx,klx))then
          if(klx%head%k == ktfoper+mtfcomplex)then
            isp=isp00
            return
          endif
        endif
      endif
      if(comp)then
        kx=kxmakelist(isp0,klx)
        klx%head%k=ktfoper+iah
        isp=isp00
      else
        kx%k=ktfref
      endif
      return
      end

      recursive function tfcomposefun(isp1,iah,full,irtc)
     $     result(kx)
      use tfstk
      use modul,only:tfmodule
      use pmat,only:itfpmatc
      implicit none
      type (sad_descriptor) kx,dh,tfefunrefd,tfefunrefu,tfefunreff
      integer*8 ka,kti,kai,i
      integer*4 ,intent(in):: isp1,iah
      integer*4 ,intent(out):: irtc
      integer*4 narg,id,isp2,j
      real*8 rimmediate0,v
      logical*4 ,intent(in):: full
      logical*4 re
      irtc=1
      kx=dxnullo
      id=iget_fun_id(int8(iah))
      select case (id)
      case (nfunif)
        narg=isp-isp1
        if(narg >= 2 .and. narg <= 4)then
          if(ktastk(isp1+1) == 0)then
            kx=merge(dtastk(isp1+3),dxnullo,narg >= 3)
          elseif(ktfrealq(ktastk(isp1+1)))then
            kx=dtastk(isp1+2)
          else
            kx%k=ktfoper+mtfnull
            return
          endif
          irtc=0
          return
        endif
        return

      case (nfunwhich)
        if(mod(isp-isp1,2) == 0)then
          isp2=isp1+1
          do j=isp1+1,isp-1,2
            if(ktfrealq(ktastk(j),v))then
              if(abs(v) /= 0.d0)then
                kx=dtastk(j+1)
                irtc=0
                return
              else
                isp2=j+2
              endif
            else
              isp2=j
              exit
            endif
          enddo
          if(isp2 > isp1+1)then
            dh%k=ktfoper+int8(iah)
            kx=kxcrelistm(isp-isp2+1,ktastk(isp2:isp),dh)
            irtc=0
          endif
        endif
        return

      case (nfunswitch)
        if(mod(isp-isp1,2) == 1)then
          if(tfconstq(dtastk(isp1+1)))then
            isp2=isp1+1
            do j=isp1+2,isp-1,2
              if(tfconstq(dtastk(j)))then
                if(itfpmatc(dtastk(isp1+1),dtastk(j)) >= 0)then
                  kx=dtastk(j+1)
                  irtc=0
                  return
                else
                  isp2=j+1
                endif
              else
                isp2=j-1
                exit
              endif
            enddo
            if(isp2 > isp1+1)then
              dh%k=ktfoper+int8(iah)
              dtastk(isp2)=dtastk(isp1+1)
              kx=kxcrelistm(isp-isp2+1,ktastk(isp2:isp),dh)
              irtc=0
            endif
          endif
        endif
        return

      case (nfunmodule)
        if(full)then
 11       if(isp == isp1+2)then
            kx=tfmodule(isp1,.true.,.false.,irtc)
            if(irtc /= 0)then
              if(irtc > 0 .and. ierrorprint /= 0)then
                call tfreseterror
              endif
              irtc=1
            endif
            isp=isp1+2
          elseif(isp > isp1+2)then
            dh=dtastk(isp1+1)
            dtastk(isp1+1)=dtastk(isp1)
            kx=tfcomposefun(isp1+1,iah,.true.,irtc)
            dtastk(isp1+1)=dh
            if(irtc /= 0)then
              if(irtc > 0 .and. ierrorprint /= 0)then
                call tfreseterror
              endif
              irtc=1
              return
            endif
            isp=isp1+2
            dtastk(isp)=kx
            go to 11
          endif
        endif
        return

      case (nfunwith)
        if(full)then
 21       if(isp == isp1+2)then
            kx=kxcompose(isp1)
            irtc=0
          elseif(isp > isp1+2)then
            dh=dtastk(isp1+1)
            dtastk(isp1+1)=dtastk(isp1)
            kx=tfcomposefun(isp1+1,iah,.true.,irtc)
            dtastk(isp1+1)=dh
            if(irtc /= 0)then
              if(irtc > 0 .and. ierrorprint /= 0)then
                call tfreseterror
              endif
              irtc=1
              return
            endif
            isp=isp1+2
            dtastk(isp)=kx
            go to 21
          endif
          return
        endif

      case (nfunlength)
        if(isp == isp1+1)then
          if(ktftype(ktastk(isp)) == ktflist)then
            ka=iand(ktamask,ktastk(isp))
            if(klist(ka) == ktfoper+mtflist)then
              if(iand(ilist(2,ka-3),lnonreallist) /= 0)then
                do i=ka+1,ka+ilist(2,ka-1)
                  kti=ktftype(klist(i))
                  if(kti == ktfsymbol)then
                    return
                  elseif(kti == ktflist)then
                    kai=iand(ktamask,klist(i))
                    if(klist(kai) == ktfoper+mtfnull .or.
     $                   klist(kai) == ktfoper+mtfslotseq)then
                      return
                    endif
                  endif
                enddo
              endif
              kx=dfromr(dble(ilist(2,ka-1)))
              irtc=0
            endif
          endif
        endif
        return
      end select
      irtc=-1
      kx%k=ktfoper+mtfnull
      if(ktfnumericq(klist(ifunbase+iah)))then
        rimmediate0=rlist(iaximmediate)
        rlist(iaximmediate)=-1.d0
        if(isp == isp1+1)then
          if(ktfrealq(ktastk(isp)))then
            kx=tfefunrefd(isp1,irtc)
          elseif(tfconstq(ktastk(isp)))then
            if(ktflistq(ktastk(isp)))then
              if(ilist(2,ktfaddr(ktastk(isp))-1) <= 15)then
                kx=tfefunrefu(isp1,irtc)
              endif
            else
              kx=tfefunrefu(isp1,irtc)
            endif
          endif
        else
          re=.false.
          do i=isp1+1,isp
            if(ktfrealq(ktastk(i)))then
            elseif(tfconstq(ktastk(i)))then
              re=.true.
            else
              rlist(iaximmediate)=rimmediate0
              return
            endif
          enddo
          kx=tfefunreff(isp1,re,irtc)
        endif
        rlist(iaximmediate)=rimmediate0
      endif
      return
      end

      subroutine tfevalb(string,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tfeval
      integer*4 ,intent(out):: irtc
      integer*4 istop
      character*(*) string
      levele=levele+1
      kx=tfeval(string,1,istop,.false.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfevalc(string)
      use tfstk
      implicit none
      type (sad_descriptor) kx,tfeval
      integer*4 irtc,istop,l
      character*(*) ,intent(in):: string
      levele=levele+1
      kx=tfeval(string,1,istop,.false.,irtc)
      if(irtc > 0 .and. ierrorprint /= 0)then
        call tfreseterror
      endif
      l=itfdownlevel()
      return
      end

      subroutine tfevals(string,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tfeval
      integer*4 ,intent(out):: irtc
      integer*4 istop
      character*(*) ,intent(in):: string
      levele=levele+1
      kx=tfeval(string,1,istop,.false.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      logical*4 function tfreadevalbuf(istart,istop,l,ipr)
      use tfcsi
      implicit none
      integer*4 , intent (in) :: ipr
      integer*4 , intent(out) :: l,istart
      integer*4 , intent(inout) :: istop
      integer*4 ip1
      ip1=ipoint
      ipoint=istop
      call tprmptget(ipr,.true.)
      if(rep)then
        ipoint=ip1
      else
        istop=ipoint
      endif
      tfreadevalbuf=ios == 0
      if(tfreadevalbuf)then
        istart=istop
        l=lrecl
      endif
      return
      end

      character*(*) function tfkname(k)
      use tfstk
      implicit none
      integer*8 ,intent(in):: k
      if(ktfrealq(k))then
        tfkname='Real'
      else
        select case(ktftype(k))
        case(ktfoper)
          tfkname='Function'
        case(ktflist)
          tfkname='List'
        case(ktfstring)
          tfkname='String'
        case(ktfsymbol)
          tfkname='Symbol'
        case(ktfpat)
          tfkname='Pattern'
        case(ktfref)
          tfkname='Reference'
        case default
          tfkname='Unknown'
        end select
      endif
      return
      end
