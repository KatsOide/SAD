      module eeval
      use tfstk
      
      contains
      function tfeeval(k,ref,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      logical*4 ,intent(in):: ref
      integer*4 ,intent(out):: irtc
      if(ref)then
        kx=tfeevalref(k,irtc)
      else
        kx=tfeevaldef(k,irtc)
      endif
      return
      end

      function tfeevalref(k,irtc) result(kx)
      use mackw
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: list
      integer*8 ka
      integer*4 ,intent(out):: irtc
      kx=k
      irtc=0
      if(ktflistq(k,list))then
        if(iand(lconstlist,list%attr) == 0)then
          kx=tfleval(list,.true.,irtc)
        endif
      elseif(ktfsymbolq(k))then
        kx=tfsyeval(k,irtc)
      elseif(ktfpatq(k))then
        kx=tfpateval(k,irtc)
      elseif(ktfrefq(k,ka))then
        kx=dlist(ka)
      endif
      return
      end

      function tfleval(list,ref,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist) ,intent(inout):: list
      integer*8 kaa
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: ref
      kaa=ksad_loc(list%head%k)
      kx%k=ktflist+kaa
      irtc=0
      if(iand(lconstlist,list%attr) /= 0)then
        return
      endif
      if(iand(lnoconstlist,list%attr) == 0)then
        if(tfconstlistqo(list))then
          return
        endif
      endif
      list%ref=list%ref+1
      call tfseval(kaa,list%nl,list%head,kx,.false.,
     $     ktfreallistq(list),ref,irtc)
      call tflocal1(kaa)
      return
      end

      subroutine tflevalstk(list,ref,irtc)
      implicit none
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist) , pointer :: kl
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: ref
      if(list%head%k == ktfoper+mtfnull)then
        call tfevallstkall(list,ref,ref,irtc)
      else
        isp=isp+1
        if(iand(lconstlist,list%attr) /= 0)then
          dtastk(isp)=sad_descr(list)
          irtc=0
        else
          dtastk(isp)=tfleval(list,ref,irtc)
          if(irtc /= 0)then
            return
          endif
          if(ktfsequenceq(ktastk(isp),kl))then
            isp=isp-1
            call tfgetllstkall(kl)
          endif
        endif
      endif
      return
      end

      recursive subroutine tfargevalstk(isp0,list,m,iaat,mf,av,irtc)
      implicit none
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist), pointer :: kli
      type (sad_descriptor) ki
      integer*8 ,intent(in):: iaat
      integer*4 ,intent(in):: m,isp0
      integer*4 ,intent(out):: irtc
      integer*4 mf,i
      logical*4 ,intent(in):: av
      irtc=0
      if(m <= 0)then
        return
      endif
      if(av)then
        call tfcopyarray(list%dbody(1:m),dtastk(isp+1:isp+m),m)
c        dtastk(isp+1:isp+m)=list%dbody(1:m)
        isp=isp+m
      else
        select case (iaat)
        case (1:)
          do i=1,m
            ki=list%dbody(i)
            if(ktflistq(ki,kli))then
              if(kli%head%k == ktfoper+mtfnull)then
                call tfargevalstk(isp0,kli,kli%nl,iaat,mf,
     $               ktfreallistq(kli),irtc)
                if(irtc /= 0)then
                  return
                endif
              elseif(ilist(2,iaat+min(isp-isp0+1,mf+1)) == 0)then
                call tflevalstk(kli,.true.,irtc)
                if(irtc /= 0)then
                  return
                endif
              else
                isp=isp+1
                dtastk(isp)=ki
              endif
            elseif((ktfsymbolq(ki) .or. ktfpatq(ki)) .and.
     $             ilist(2,iaat+min(isp-isp0+1,mf+1)) == 0)then
              call tfevalstk(ki,.true.,irtc)
              if(irtc /= 0)then
                return
              endif
            else
              isp=isp+1
              dtastk(isp)=ki
            endif
          enddo
        case (-iattrholdfirst)
          do i=1,m
            ki=list%dbody(i)
            if(ktfsequenceq(ki,kli))then
              call tfevallstkall(kli,isp > isp0,.true.,irtc)
              if(irtc /= 0)then
                return
              endif
            elseif(isp > isp0)then
              isp=isp+1
              dtastk(isp)=tfeevalref(ki,irtc)
              if(irtc /= 0)then
                return
              endif
            else
              isp=isp+1
              dtastk(isp)=ki
            endif
          enddo
        case (-iattrholdrest)
          do i=1,m
            ki=list%dbody(i)
            if(ktfsequenceq(ki,kli))then
              call tfevallstkall(kli,isp == isp0,.false.,irtc)
              if(irtc /= 0)then
                return
              endif
            elseif(isp == isp0)then
              isp=isp+1
              dtastk(isp)=tfeevalref(ki,irtc)
              if(irtc /= 0)then
                return
              endif
            else
              isp=isp+1
              dtastk(isp)=ki
            endif
          enddo
        case (-iattrholdall)
          call tfgetllstkall(list)
        case default
          call tfevallstkall(list,.true.,.true.,irtc)
        end select
      endif
      return
      end

      subroutine tfseqevalstkall(ks,m,av,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: ks(m)
      integer*4 ,intent(in):: m
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0
      logical*4 ,intent(in):: av
      irtc=0
      if(m <= 0)then
        return
      endif
      if(av)then
        isp0=isp
        dtastk(isp0+1:isp0+m)=ks
        isp=isp0+m
      else
        do i=1,m
          call tfevalstk(ks(i),.true.,irtc)
          if(irtc /= 0)then
            return
          endif
        enddo
      endif
      return
      end

      subroutine tfseqevalstk(ks,m,i,av,irtc)
      implicit none
      type (sad_dlist), pointer :: kli
      type (sad_descriptor) ,intent(in):: ks(m)
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: m,i
      logical*4 ,intent(in):: av
      irtc=0
      if(i <= m)then
        if(av)then
          isp=isp+1
          dtastk(isp)=ks(i)
        else
          if(ktflistq(ks(i),kli))then
            call tflevalstk(kli,.true.,irtc)
          elseif(ktfsymbolq(ks(i)) .or. ktfpatq(ks(i)))then
            call tfevalstk(ks(i),.true.,irtc)
          else
            isp=isp+1
            dtastk(isp)=ks(i)
          endif
        endif
      endif
      return
      end

      function tfevallev(list,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ki0,ki
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: listi
      integer*8 kai
      integer*4 ,intent(out):: irtc
      integer*4 i,isp1
      logical*4 ev
      if(ktfreallistq(list) .or. iand(list%attr,kconstarg) /= 0)then
        kx=sad_descr(list)
        irtc=0
        return
      endif
      ev=.false.
      isp=isp+1
      isp1=isp
      dtastk(isp)=list%head
      do i=1,list%nl
        ki=list%dbody(i)
        kai=ktfaddr(ki)
        select case(ki%k-kai)
        case (ktflist)
          call loc_dlist(kai,listi)
          if(iand(lconstlist,listi%attr) == 0)then
            ki0=ki
            ki=tfleval(listi,.true.,irtc)
            if(irtc /= 0)then
              go to 9000
            endif
            ev=ev .or. ki%k /= ki0%k
          endif
        case (ktfsymbol)
          ki0=ki
          ki=tfsyeval(ki0,irtc)
          if(irtc /= 0)then
            go to 9000
          endif
          ev=ev .or. ki%k /= ki0%k
        case (ktfpat)
          ki0=ki
          ki=tfpateval(ki0,irtc)
          if(irtc /= 0)then
            go to 9000
          endif
          ev=ev .or. ki%k /= ki0%k
        end select
        isp=isp+1
        dtastk(isp)=ki
        if(ktfsequenceq(ki%k,listi))then
          ev=.true.
          isp=isp-1
          call tfgetllstkall(listi)
        endif
      enddo
      if(ev)then
        kx=kxcompose(isp1)
      else
        kx=sad_descr(list)
      endif
      irtc=0
 9000 isp=isp1-1
      return
      end

      subroutine tfevalstk(k,ref,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: ref
      if(ktfsequenceq(k,kl))then
        call tfevallstkall(kl,ref,ref,irtc)
      else
        isp=isp+1
        dtastk(isp)=tfeeval(k,ref,irtc)
        if(irtc /= 0)then
          return
        endif
        if(ktfsequenceq(ktastk(isp),kl))then
          isp=isp-1
          call tfgetllstkall(kl)
        endif
      endif
      return
      end

      subroutine tfevallstkall(list,ref1,ref,irtc)
      implicit none
      type (sad_dlist) ,intent(in):: list
      type (sad_descriptor) ki
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0,m
      logical*4 ,intent(in):: ref,ref1
      m=list%nl
      if(ktfreallistq(list))then
        isp0=isp
        call tfcopyarray(list%dbody(1:m),dtastk(isp0+1:isp0+m),m)
c        dtastk(isp0+1:isp0+m)=list%dbody(1:m)
        isp=isp0+m
      elseif(m > 0)then
        ki=list%dbody(1)
        if(ktfrealq(ki))then
          isp=isp+1
          dtastk(isp)=ki
        else
          call tfevalstk(ki,ref1,irtc)
          if(irtc /= 0)then
            return
          endif
        endif
        do i=2,m
          ki=list%dbody(i)
          if(ktfrealq(ki))then
            isp=isp+1
            dtastk(isp)=ki
          else
            call tfevalstk(ki,ref,irtc)
            if(irtc /= 0)then
              return
            endif
          endif
        enddo
      endif
      irtc=0
      return
      end

      function tfeevaldef(k,irtc) result(kx)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: list
      type (sad_symbol), pointer :: sym
      integer*4 ,intent(out):: irtc
      if(ktflistq(k,list))then
        if(iand(lconstlist,list%attr) == 0)then
          kx=tfleval(list,.false.,irtc)
          return
        endif
      elseif(ktfsymbolq(k,sym) .and. sym%override == 0)then
        sym=>tfsydef(sym)
        kx=sad_descr(sym)
        irtc=0
        return
      endif
      irtc=0
      kx=k
      return
      end

      function tfpateval(k,irtc) result(kx)
      use tfcode
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ke
      type (sad_pat), pointer :: pat
      integer*8 ka,kae,kte
      integer*4 ,intent(out):: irtc
      ka=ktfaddrd(k)
      call loc_pat(ka,pat)
      kx%k=ktfpat+ka
      irtc=0
      ke=pat%expr
      kae=ktfaddr(ke)
      kte=ke%k-kae
      if(kte /= ktfref)then
        ke=tfeevalref(ke,irtc)
        if(irtc /= 0)then
          return
        endif
        kae=ktfaddr(ke)
        kte=ke%k-kae
      endif
      if(ke%k /= pat%expr%k)then
        kx=kxpcopyss(ke,pat%head,
     $       pat%sym%alloc,pat%default)
      endif
      return
      end

      function tfsyeval(ka,irtc) result(kx)
      use tfcode
      use iso_c_binding
      use mackw
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: ka
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      type (sad_object), pointer :: obj
      integer*4 maxlevel
      parameter (maxlevel=4096)
      integer*8 kas,j
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,level1
      integer*4 level,ig,ig1
      data level/0/
      kas=ktfaddrd(ka)
      call loc_sym(kas,sym)
      irtc=0
      ig=sym%gen
      if(ig == -3)then
        kx%k=ktfsymbol+kas
        return
      endif
      if(sym%override /= 0)then
        call loc_symdef(kas,symd)
      else
        call loc_namtbl(sym%loc,loc)
        kas=loc%symdef
        ig=max(0,ig)
        do while(kas > 0)
          call loc1_symdef(kas,symd)
          ig1=max(0,symd%sym%gen)
          if(ig == ig1)then
            exit
          elseif(ig < ig1)then
            kas=symd%next
          else
            exit
          endif
        enddo
      endif
      kx=symd%value
      if(kx%k /= ktfsymbol+kas .and. ktfnonrealq(kx%k))then
        level=level+1
        if(level > maxlevel-64)then
          level1=level
          level=0
          irtc=itfmessage(999,'General::deep','""')
          level=level1
          return
        endif
        kx=tfeevalref(kx,irtc)
        level=level-1
        if(ktfobjq(kx,obj))then
          if(ktfaddr(obj%alloc) == 0)then
            j=itflocal+levele
            obj%alloc%k=ktftype(obj%alloc%k)+klist(j)
            klist(j)=ksad_loc(obj%alloc%k)
          endif
        endif
      endif
      return
      end

      subroutine tfseval(ks,ns,kh,kx,stk,av,ref,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      use funs
      use tfcsi,only:icslfno
      implicit none
      type (sad_descriptor) kx,kf,kh,kl,tfefunrefu
      type (sad_funtbl), pointer :: fun
      type (sad_symbol), pointer :: sym
      type (sad_dlist), pointer :: list,klf,kls,kls1
      integer*4 maxlevel
      parameter (maxlevel=2**12)
      integer*8 ks,iaat,kaf
      integer*4 irtc,ns,iaf,isp10,isp11,j
      integer*4 isp1,level,itfgetrecl,
     $     isp0,mf,i1,level1,i,itfmessage,lpw,l
      real*8 v
      logical*4 ref,ev,evalh,rep,stk,av
      data level/0/
      level=level+1
      if(level .ge. maxlevel-64)then
        level1=level
        level=0
        irtc=itfmessage(999,'General::deep','""')
        level=level1
        go to 9000
      endif
      call loc_sad(ks,kls)
      isp0=isp
      kf=kh
      evalh=.false.
      kaf=ktfaddr(kf%k)
      select case(kf%k-kaf)
      case (ktfoper)
        call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
        iaf=-fun%id
        if(fun%narg < 0)then
          go to 100
        endif
      case (ktfsymbol)
        if(ref)then
          kf=tfsyeval(kf,irtc)
          if(irtc /= 0)then
            go to 8000
          endif
        else
          call loc_sym(kaf,sym)
          sym=>tfsydef(sym)
          kf=sad_descr(sym)
        endif
      case (ktflist)
        call loc_dlist(kaf,list)
        if(iand(lconstlist,list%attr) == 0)then
          kf=tfleval(list,ref,irtc)
          if(irtc /= 0)then
            go to 8000
          endif
        endif
      case (ktfpat)
        if(ref)then
          kf=tfpateval(kf,irtc)
          if(irtc /= 0)then
            go to 8000
          endif
        endif
      end select
      ev=.true.
      iaat=0

      do
        if(ktfoperq(kf%k))then
          kaf=ktfaddr(kf)
          call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
          iaat=klist(ifunbase+kaf)+1
          mf=fun%narg
          if(mf < 0)then
            iaf=-fun%id
            evalh=.true.
            exit
          endif
          ev=fun%mapeval(2,1) == -2
          if(fun%mapeval(2,1) == -1)then
            iaat=0
          endif
        elseif(ktfsymbolq(kf%k,sym))then
          if(sym%override /= 0)then
            iaat=iand(iattrholdall,sym%attr)
            if(iaat /= 0)then
              ev=.false.
              iaat=merge(i00,-iaat,iaat == iattrholdall)
            endif
          endif
        elseif(ktflistq(kf,klf))then
          if(klf%head%k == ktfoper+mtfnull)then
            if(klf%nl == 1)then
              kf=klf%dbody(1)
              cycle
            elseif(klf%nl == 0)then
              kf%k=ktfoper+mtfnull
              evalh=.true.
              exit
            endif
          elseif(ktfsymbolq(klf%head,sym) .and. .not. ref)then
            if(sym%gen == -3 .and. ktfreallistq(klf))then
              ev=.false.
              iaat=0
            endif
          endif
        endif
        isp=isp0+1
        isp1=isp
        go to 3000
      enddo

 100  isp=isp0+1
      isp1=isp
      select case(iaf)
      case(mtfnull,mtflist,mtfrule,mtfrepeated,mtfrepeatednull)
        if(stk)then
          ev=.true.
          go to 3000
        endif
        if(evalh .or. .not. ref)then
          if(tfonstackq(ks) .or. kls%ref > 0)then
            kls1=>tfduplist(kls)
            kls=>kls1
          endif
          kls%head%k=ktfoper+iaf
          if(ref .and. tfconstlistqo(kls))then
            kx=sad_descr(kls)
            irtc=0
            go to 8000
          endif
        endif
        if(ref)then
          kx=tfevallev(kls,irtc)
        else
          call tfevallstkall(kls,.false.,.false.,irtc)
          if(irtc == 0)then
            levele=levele+1
            go to 6000
          endif
        endif
      case(mtfset)
        rep=tfgetseqstk(ks,ns)
        if(isp > isp1)then
          dtastk(isp)=tfeevalref(dtastk(isp),irtc)
          if(irtc /= 0)then
            go to 8000
          endif
        endif
        levele=levele+1
        go to 6000
      case (mtfand)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%dbody(1),ns,i,av,irtc)
          if(irtc /= 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j),v))then
              if(abs(v) == 0.d0)then
                kx%k=0
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp == isp1)then
          kx%k=ktftrue
        elseif(isp == isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfor)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%dbody(1),ns,i,av,irtc)
          if(irtc /= 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j),v))then
              if(v/=0.d0)then
                kx%k=ktftrue
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp == isp1)then
          kx%k=ktffalse
        elseif(isp == isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfpart)
        if(ref .or. ns == 0)then
          ev=.true.
          go to 3000
        endif
        isp=isp+1
        dtastk(isp)=tfeevaldef(kls%dbody(1),irtc)
        if(irtc /= 0)then
          go to 8000
        endif
        call tfseqevalstkall(kls%dbody(2),ns-1,av,irtc)
        if(irtc == 0)then
          levele=levele+1
          go to 6000
        endif
      case (mtfslot,mtfslotseq)
        kx=tfslot(iaf,kls,ref,irtc)
      case (mtfcomp)
        if(ns == 0)then
          kx%k=ktfoper+mtfnull
          go to 8000
        endif
        i1=1
        do
          if(ltrace > 0)then
            levele=levele+1
            lpw=min(131,itfgetrecl())
            do i=i1,ns
              call tfprint1(kls%dbody(i),icslfno(),-lpw,4,.true.,.true.,irtc)
              kx=tfeevalref(kls%dbody(i),irtc)
              if(irtc /= 0)then
                go to 1320
              endif
            enddo
          else
            levele=levele+1
            irtc=0
            do i=i1,ns
              kx=kls%dbody(i)
              if(ktflistq(kx,list))then
                kx=tfleval(list,.true.,irtc)
              elseif(ktfsymbolq(kx))then
                kx=tfsyeval(kx,irtc)
              elseif(ktfpatq(kx))then
                kx=tfpateval(kx,irtc)
              endif
              if(irtc /= 0)then
                go to 1320
              endif
            enddo
          endif
          go to 7000
 1320     if(irtc > irtcret)then
            go to 7000
          endif
          call tfcatchreturn(irtcgoto,kl,irtc)
          l=itfdownlevel()
          if(irtc /= 0)then
            exit
          endif
          call tffindlabel(kls,ns,i1,kl)
          if(i1 .le. 0)then
            call tfthrow(irtcgoto,kl,irtc)
            exit
          endif
          i1=i1+1
        enddo
      case (mtffun,mtfpattest,mtftagset,mtfhold)
        rep=tfgetseqstk(ks,ns)
        if(rep .or. stk .or. evalh)then
          levele=levele+1
          go to 6000
        endif
        kx%k=ktflist+ks
      case default
        write(*,*)'tfleval-implementation error: ',kaf,iaf
        call abort
      end select
      go to 8000

 3000 if(levele .ge. maxlevele-32)then
        irtc=itfmessage(999,'General::deep','""')
        go to 8000
      endif
      levele=levele+1
      if(ev)then
        call tfseqevalstkall(kls%dbody(1),ns,av,irtc)
        if(irtc /= 0)then
          go to 7000
        endif
      elseif(iaat == 0 .or. av)then
        rep=tfgetseqstk(ks,ns)
      else
        call tfargevalstk(isp1,kls,ns,iaat,mf,.false.,irtc)
        if(irtc /= 0)then
          go to 7000
        endif
      endif
 6000 dtastk(isp1)=kf
      if(ref)then
        kx=tfefunrefu(isp1,irtc)
      else
        call tfefundef(isp1,kx,irtc)
      endif
 7000 continue
      call tfconnect(kx,irtc)
 8000 isp=isp0
 9000 level=max(0,level-1)
      return
      end

      function tfslot(iopc,kls,ref,irtc) result(kx)
      use funs,only:tfsequence
      implicit none
      type (sad_descriptor) kx,ka
      type (sad_dlist) ,intent(in):: kls
      type (sad_symbol), pointer :: sym
      type (sad_namtbl), pointer :: nam
      integer*8 kaa
      integer*4 ,intent(in):: iopc
      integer*4 ,intent(out):: irtc
      integer*4 ind,nc,isp1,isps,
     $     itfmessage,ns,ipf0,naf0,ls,isp2,itfmessagestr
      real*8 ffval,vx
      character*256 name
      character*12 inds
      logical*4 ,intent(in):: ref
      logical*4 exist
      kx=dxnullo
      ns=kls%nl
      if(ns > 1)then
        irtc=itfmessage(9,'General::narg','"0 or 1"')
        return
      endif
      if(ns == 0)then
        ind=1
      else
        ka=kls%dbody(1)
        if(ktfoperq(ka,kaa))then
          if(kaa == mtfnull)then
            ind=1
          else
            irtc=itfmessage(999,'General::invop',' ')
            return
          endif
        elseif(ktfrealq(ka,ind))then
        elseif(ktfsymbolq(ka,sym) .and. iopc == mtfslot)then
          call sym_namtbl(sym,nam)
          nc=nam%str%nch+1
          name(2:nc)=nam%str%str(1:nc-1)
          name(1:1)='#'
          call capita(name(1:nc))
          vx=ffval(name(1:nc),exist)
          if(exist)then
            kx=dfromr(vx)
            irtc=0
          else
            irtc=itfmessage(999,'FFS::undef','"element"')
          endif
          return
        else
          irtc=itfmessage(999,'General::wrongtype','"Number or symbol"')
          return
        endif
        if(ind < 0)then
          ind=napuref+ind+1
        endif
      endif
      isps=ipurefp+ind
      if(iopc == mtfslot)then
        if(ipurefp == 0 .or. ind <= 0 .or. ind > napuref)then
          call strfromil(ind,inds,ls)
          irtc=itfmessagestr(999,'General::slot',
     $         '#'//inds(:ls))
          return
        endif
        kx=dtastk(isps)
      else
        if(ipurefp == 0 .or. ind <= 0)then
          call strfromil(ind,inds,ls)
          irtc=itfmessagestr(999,'General::slot',
     $         '##'//inds(:ls))
          return
        endif
        isp1=isp
        isp2=ipurefp+napuref
        kx=tfsequence(isps-1,isp2)
      endif
      ipf0=ipurefp
      naf0=napuref
      ipurefp=itastk(1,ipf0+naf0+1)
      napuref=itastk(2,ipf0+naf0+1)
      kx=tfeeval(kx,ref,irtc)
      ipurefp=ipf0
      napuref=naf0
      return
      end

      subroutine tffindlabel(list,m,i,kr)
      use tfstk
      implicit none
      type (sad_descriptor) kr
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: listj
      integer*4 ,intent(out):: i
      integer*4 ,intent(in):: m
      integer*4 j
      type (sad_descriptor), save :: kxlabel
      data kxlabel%k /0/
      if(kxlabel%k == 0)then
        kxlabel=kxsymbolf('Label',5,.true.)
      endif
      do j=1,m
        if(ktflistq(list%dbody(j),listj))then
          if(listj%nl == 1)then
            if(tfsameq(kxlabel,listj%head))then
              if(tfsameq(kr,listj%dbody(1)))then
                i=j
                return
              endif
            endif
          endif
        endif
      enddo
      i=0
      return
      end

      logical*4 function tfgetseqstk(ks,ns)
      use tfstk
      implicit none
      type (sad_dlist), pointer :: kl
      integer*8 ,intent(in):: ks
      integer*8 ki
      integer*4 ,intent(in):: ns
      integer*4 i
      tfgetseqstk=.false.
      if(ns > 0)then
        do i=1,ns
          ki=klist(ks+i)
          if(ktfsequenceq(ki,kl))then
            tfgetseqstk=.true.
            call tfgetllstkall(kl)
            cycle
          endif
          isp=isp+1
          ktastk(isp)=ki
        enddo
      endif
      return
      end

      subroutine tfstringreplace(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx,kr
      type (sad_strbuf), pointer :: strb
      type (sad_string), pointer :: str,strs,stri
      integer*4 isp1,irtc,isp0,i,ir,ls,isp2,
     $     j,imin,ii,nr,indexb,itfmessage
      logical*4 full
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. ktfstringq(dtastk(isp1+1),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      isp0=isp
      call tfreplace(dxnullo,dtastk(isp1+2),kx,
     $     .false.,.false.,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      if(isp == isp0)then
        kx=dtastk(isp1+1)
        irtc=0
        return
      endif
      do i=isp0+1,isp,2
        if(.not. ktfstringq(dtastk(i)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"String -> String"')
          return
        endif
      enddo
      nullify(strb)
      ir=1
      ls=str%nch
      isp2=isp
 1    j=0
      imin=ls+1
      do i=isp0+1,isp2,2
        call descr_sad(dtastk(i),stri)
        ii=indexb(str%str,ls,stri%str,stri%nch,ir)
        if(ii > 0 .and. ii < imin)then
          imin=ii
          j=i
          itastk2(1,i)=stri%nch
        endif
      enddo
      if(j /= 0)then
        if(.not. associated(strb))then
          call getstringbuf(strb,0,.true.)
        endif
        if(imin == ir+1)then
          call putstringbufb1(strb,str%str(ir:ir))
        elseif(imin > ir)then
          call putstringbufb(strb,str%str(ir:imin-1),imin-ir,full)
        endif
        if(.not. ktfstringq(dtastk(j+1)))then
          kr=tfeevalref(dtastk(j+1),irtc)
          if(irtc /= 0)then
            go to 9000
          endif
          if(.not. ktfstringq(kr))then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"List of (String -> String)"')
            go to 9000
          endif
          dtastk(j+1)=kr
        endif
        call descr_sad(dtastk(j+1),strs)
        nr=strs%nch
        if(nr == 1)then
          call putstringbufb1(strb,strs%str)
        elseif(nr > 0)then
          call putstringbufb(strb,strs%str,nr,full)
        endif
        ir=imin+itastk2(1,j)
        if(ir <= ls)then
          go to 1
        endif
      endif
      if(.not. associated(strb))then
        kx=dtastk(isp1+1)
      else
        if(ls == ir)then
          call putstringbufb1(strb,str%str(ir:ir))
        elseif(ls > ir)then
          call putstringbufb(strb,str%str(ir:ls),ls-ir+1,full)
        endif
        kx=kxstringbuftostring(strb)
      endif
      isp=isp0
      irtc=0
      return
 9000 call tfreestringbuf(strb)
      isp=isp0
      return
      end

      end module eeval
