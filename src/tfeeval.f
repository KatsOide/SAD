      subroutine tfeeval(k,kx,ref,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx
      integer*4 irtc
      logical*4 ref
      if(ref)then
        call tfeevalref(k,kx,irtc)
      else
        call tfeevaldef(k,kx,irtc)
      endif
      return
      end

      subroutine tfeevalref(k,kx,irtc)
      use tfstk
      use mackw
      implicit none
      type (sad_descriptor) k,kx
      type (sad_list), pointer :: list
      integer*8 ka
      integer*4 irtc
      kx=k
      irtc=0
      if(ktflistqd(k,list))then
        if(iand(lconstlist,list%attr) .eq. 0)then
          call tfleval(list,kx,.true.,irtc)
        endif
      elseif(ktfsymbolqd(k))then
        call tfsyeval(k,kx,irtc)
      elseif(ktfpatqd(k))then
        call tfpateval(k,kx,irtc)
      elseif(ktfrefqd(k,ka))then
        kx=dlist(ka)
      endif
      return
      end

      subroutine tfleval(list,kx,ref,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_list) list
      integer*8 kaa
      integer*4 irtc
      logical*4 ref,tfconstlistqo
      kaa=ksad_loc(list%head)
      kx%k=ktflist+kaa
      irtc=0
      if(iand(lconstlist,list%attr) .ne. 0)then
        return
      endif
      if(iand(lnoconstlist,list%attr) .eq. 0)then
        if(tfconstlistqo(list))then
          return
        endif
      endif
      list%ref=list%ref+1
      call tfseval(kaa,list%nl,list%head,kx,.false.,
     $     ktfreallistqo(list),ref,irtc)
      call tflocal1(kaa)
      return
      end

      subroutine tfseval(ks,ns,kh,kx,stk,av,ref,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_funtbl), pointer :: fun
      type (sad_symbol), pointer :: sym
      type (sad_list), pointer :: list,klf,kls
      integer*4 maxlevel
      parameter (maxlevel=2**12)
      integer*8 ks,kh,kaf,kf,kl,iaat
      integer*4 irtc,ns
      integer*4 isp1,isp10,isp11,level,itfgetrecl,
     $     j,isp0,mf,i1,level1,i,itfmessage,lpw,l,itfdownlevel
      logical*4 ref,ev,evalh,rep,stk,tfconstlistqo,tfgetseqstk,av
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
      kaf=ktfaddr(kf)
      select case(kf-kaf)
      case (ktfoper)
        call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
        if(fun%narg .lt. 0)then
          go to 100
        endif
      case (ktfsymbol)
        if(ref)then
          call tfsyeval(kaf,kf,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        else
          call loc_sym(kaf,sym)
          call tfsydef(sym,sym)
          kf=ktfsymbol+ksad_loc(sym%loc)
        endif
      case (ktflist)
        call loc_list(kaf,list)
        if(iand(lconstlist,list%attr) .eq. 0)then
          call tfleval(list,kf,ref,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
      case (ktfpat)
        if(ref)then
          call tfpateval(kaf,kf,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
      end select
      ev=.true.
      iaat=0

      do while(.true.)
        if(ktfoperq(kf))then
          kaf=ktfaddr(kf)
          call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
          iaat=klist(ifunbase+kaf)+1
          mf=fun%narg
          if(mf .lt. 0)then
            evalh=.true.
            exit
          endif
          ev=fun%mapeval(2,1) .eq. -2
          if(fun%mapeval(2,1) .eq. -1)then
            iaat=0
          endif
        elseif(ktfsymbolq(kf,sym))then
          if(sym%override .ne. 0)then
            iaat=iand(iattrholdall,sym%attr)
            if(iaat .ne. 0)then
              ev=.false.
              if(iaat .eq. iattrholdall)then
                iaat=0
              else
                iaat=-iaat
              endif
            endif
          endif
        elseif(ktflistq(kf,klf))then
          if(klf%head .eq. ktfoper+mtfnull)then
            if(klf%nl .eq. 1)then
              kf=klf%body(1)
              cycle
            elseif(klf%nl .eq. 0)then
              kf=ktfoper+mtfnull
              evalh=.true.
              exit
            endif
          elseif(ktfsymbolq(klf%head,sym) .and. .not. ref)then
            if(sym%gen .eq. -3 .and. ktfreallistqo(klf))then
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
      select case(kaf)
c      go to (
c     $     1000,
c     $     1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,
c     $     1900,1900,1900,1900,1900,1900,1900,1100,1110,1900,
c     $     1900,1900,1000,1900,1900,1200,1900,1900,1900,1900,
c     $     1300,1010,1900,1000,1900,1900,1900,1900,1900,1900,
c     $     1010,1900,1400,1400,1900,1900,1900,1900,1900,1000,
c     $     1000,1900,1900,1900,1900,1900,1900,1900,1500,1900,
c     $     1900,1010,1900,1900,1010,1900),
c     $     kaf+1
c          null
c          m    i    +    -    *    /    v    ^    e    n    
c          >    <    g    l    E    N    ~    &&    o    c
c          [    ]    {    }    s    =    C    (    )    ,
c          ;    &    :    r    d    RepA RepR u    U    S    
c          ?    f    #    ##   .    |    M    MA   A    rept 
c          repn ineq AT   SF   TB   DB   Inc  Dec  Part @
c          msgn TagS (*   *)   Hold z
      case(mtfnull,mtflist,mtfrule,mtfrepeated,mtfrepeatednull)
        if(stk)then
          ev=.true.
          go to 3000
        endif
        if(evalh .or. .not. ref)then
          if(tfonstackq(ks) .or. kls%ref .gt. 0)then
            call tfduplist(kls,kls)
          endif
          kls%head=ktfoper+kaf
c          call tfloadlstk(lista,lista)
c          lista%head=ktfoper+kaf
c          call tfstk2l(lista,lista)
          if(ref .and. tfconstlistqo(kls))then
            kx%k=ktflist+ksad_loc(kls%head)
            irtc=0
            go to 8000
          endif
        endif
        if(ref)then
          call tfevallev(kls,kx,irtc)
        else
          call tfevallstkall(kls,.false.,.false.,irtc)
          if(irtc .eq. 0)then
            levele=levele+1
            go to 6000
          endif
        endif
      case(mtfset)
        rep=tfgetseqstk(ks,ns)
        if(isp .gt. isp1)then
          call tfeevalref(ktastk(isp),ktastk(isp),irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
        levele=levele+1
        go to 6000
      case (mtfpart)
        if(ref .or. ns .eq. 0)then
          ev=.true.
          go to 3000
        endif
        isp=isp+1
        call tfeevaldef(kls%body(1),ktastk(isp),irtc)
        if(irtc .ne. 0)then
          go to 8000
        endif
        call tfseqevalstkall(kls%body(2),ns-1,av,irtc)
        if(irtc .eq. 0)then
          levele=levele+1
          go to 6000
        endif
      case (mtfslot,mtfslotseq)
        call tfslot(kaf,kls,kx,ref,irtc)
      case (mtfcomp)
        if(ns .eq. 0)then
          kx%k=ktfoper+mtfnull
          go to 8000
        endif
        i1=1
        do while(.true.)
          if(ltrace .gt. 0)then
            levele=levele+1
            lpw=min(131,itfgetrecl())
            do i=i1,ns
              call tfprint1(kls%body(i),6,-lpw,4,.true.,
     $             .true.,irtc)
              call tfeevalref(kls%body(i),kx,irtc)
              if(irtc .ne. 0)then
                go to 1320
              endif
            enddo
          else
            levele=levele+1
            irtc=0
            do i=i1,ns
              kx%k=kls%body(i)
              if(ktflistqd(kx,list))then
                call tfleval(list,kx,.true.,irtc)
              elseif(ktfsymbolqd(kx))then
                call tfsyeval(kx,kx,irtc)
              elseif(ktfpatqd(kx))then
                call tfpateval(kx,kx,irtc)
              endif
              if(irtc .ne. 0)then
                go to 1320
              endif
            enddo
          endif
          go to 7000
 1320     if(irtc .ge. -3)then
            go to 7000
          endif
          call tfcatchreturn(2,kl,irtc)
          l=itfdownlevel()
          if(irtc .ne. 0)then
            exit
          endif
          call tffindlabel(kls,ns,i1,kl)
          if(i1 .le. 0)then
            call tfthrow(2,kl,irtc)
            exit
          endif
          i1=i1+1
        enddo
      case (mtfand)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%body(1),ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j)))then
              if(ktastk(j) .eq. 0)then
                kx%k=0
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx%k=ktftrue
        elseif(isp .eq. isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfor)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%body(1),ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j)))then
              if(ktastk(j) .ne. 0)then
                kx%k=ktftrue
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx%k=0
        elseif(isp .eq. isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtffun,mtfpattest,mtftagset,mtfhold)
        rep=tfgetseqstk(ks,ns)
        if(rep .or. stk .or. evalh)then
          levele=levele+1
          go to 6000
        endif
        kx%k=ktflist+ks
      case default
        write(*,*)'tfleval-implementation error: ',kaf
        call abort
      end select
      go to 8000

 3000 if(levele .ge. maxlevele-32)then
        irtc=itfmessage(999,'General::deep','""')
        go to 8000
      endif
      levele=levele+1
      if(ev)then
        call tfseqevalstkall(kls%body(1),ns,av,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      elseif(iaat .eq. 0 .or. av)then
        rep=tfgetseqstk(ks,ns)
      else
        call tfargevalstk(isp1,kls,ns,iaat,mf,.false.,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      endif
 6000 ktastk(isp1)=kf
c      call tfmemcheckprint('seval-efun',.false.,irtc)
c      if(irtc .ne. 0)then
c        call tfdebugprint(kf,'tfseval-efunref-in',3)
c        call tfdebugprint(ktastk(isp1+1),' ',3)
c        call tfdebugprint(ktastk(isp),' ',3)
c      endif
      if(ref)then
        call tfefunref(isp1,kx,.true.,irtc)
c          call tfdebugprint(kf,'tfseval-efunref-out',3)
c          call tfdebugprint(ktastk(isp1+1),' ',3)
c          call tfdebugprint(ktastk(isp),' ',3)
      else
        call tfefundef(isp1,kx,irtc)
      endif
 7000 continue
c      call tfdebugprint(kx,'tfseval-connect',3)
      call tfconnect(kx,irtc)
 8000 isp=isp0
 9000 level=max(0,level-1)
      return
      end

      logical*4 function tfgetseqstk(ks,ns)
      use tfstk
      implicit none
      type (sad_list), pointer :: kl
      integer*8 ks,ki
      integer*4 ns,i
      tfgetseqstk=.false.
      if(ns .gt. 0)then
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

      integer*4 function itfdownlevel()
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_list), pointer :: list
      type (sad_pat), pointer :: pat
      integer*8 j,kt,ka,i,lastfree
      integer*4 kk
      lastfree=0
      j=itflocal+levele
      i=klist(j)
      do while(i .ne. j)
        ka=ktfaddr(klist(i))
        if(ilist(1,i+1) .le. 0)then
          kt=klist(i)-ka
          if(kt .eq. ktflist)then
            call loc_list(i+2,list)
            if(iand(list%attr,lnonreallist) .ne. 0)then
              do kk=list%nl,1,-1
                call tflocalrel(list%dbody(kk),ka)
              enddo
            endif
            call tflocalrel(list%dbody(0),ka)
            i=i-list%lenp
            ilist(1,i-1)=list%lenp+list%lena+list%nl+4
          elseif(kt .eq. ktfpat)then
            call loc_pat(i+2,pat)
            if(pat%sym%ref .gt. 1)then
              call tflocalrel(ktfsymbol+ksad_loc(pat%sym%loc),ka)
              pat%sym%attr=pat%len-7
              pat%sym%gen=0
              pat%len=7
            endif
            if(pat%sym%loc .ne. 0)then
              call tflocalrel(pat%sym%alloc,ka)
            endif
            pat%sym%alloc%k=ktfsymbol
            call tflocalrel(pat%default,ka)
            call tflocalrel(pat%head,ka)
            call tflocalrel(pat%expr,ka)
          endif
          call tfreel(i,lastfree)
        else
          klist(i)=klist(i)-ka
        endif
        i=ka
      enddo
      call tfreel(int8(0),lastfree)
      klist(j)=j
      levele=max(0,levele-1)
      itfdownlevel=levele
      return
      end

      subroutine tflocalrel(k,kan)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_object), pointer :: obj
      integer*8 kan
      if(ktfobjqd(k,obj))then
        obj%ref=obj%ref-1
        if(obj%ref .le. 0)then
          obj%ref=0
          if(ktfaddr(obj%alloc) .eq. 0)then
            obj%alloc=ktftype(obj%alloc)+kan
            kan=ktfaddr(k%k)-2
          endif
        endif
      endif
      return
      end

      subroutine tfreel(ip,last)
      use tfstk
      implicit none
      integer*8 ip,last
      if(ip .eq. 0)then
        if(last .ne. 0)then
          call tfree(last)
          last=0
        endif
      else
        ilist(1,ip-1)=max(4,ilist(1,ip-1))
        if(last .eq. 0)then
          last=ip
        else
          if(ilist(1,last-1)+last .eq. ip)then
            ilist(1,last-1)=ilist(1,last-1)+ilist(1,ip-1)
          elseif(ilist(1,ip-1)+ip .eq. last)then
            ilist(1,ip-1)=ilist(1,ip-1)+ilist(1,last-1)
            last=ip
          else
            call tfree(last)
            last=ip
          endif
        endif
      endif
      return
      end

      subroutine tflocal(k)
      use tfstk
      implicit none
      type (sad_object), pointer :: obj
      integer*8 k,ka,itfroot
      if(ktfobjq(k))then
        ka=ktfaddr(k)
        call loc_obj(ka,obj)
        obj%ref=obj%ref-1
        if(obj%ref .le. 0)then
          obj%ref=0
          if(ktfaddr(obj%alloc) .eq. 0)then
            itfroot=itflocal+levele
            obj%alloc=obj%alloc+ktfaddr(klist(itfroot))
            klist(itfroot)=ka-2
          endif
        endif
      endif
      return
      end

      subroutine tflocal1(k)
      use tfstk
      implicit none
      type (sad_object), pointer :: obj
      integer*8 ka,k,itfroot
      ka=ktfaddr(k)
      call loc_obj(ka,obj)
      obj%ref=obj%ref-1
      if(obj%ref .le. 0)then
        obj%ref=0
        if(ktfaddr(obj%alloc) .eq. 0)then
          itfroot=itflocal+levele
          obj%alloc=obj%alloc+ktfaddr(klist(itfroot))
          klist(itfroot)=ka-2
        endif
      endif
      return
      end

      integer*4 function itfuplevel()
      use tfstk
      implicit none
      levele=min(maxlevele,levele+1)
      klist(itflocal+levele)=itflocal+levele
      itfuplevel=levele
      return
      end

      subroutine tflevalstk(list,ref,irtc)
      use tfstk
      implicit none
      type (sad_list) list
      type (sad_list) , pointer :: kl
      integer*4 irtc
      logical*4 ref
      if(list%head .eq. ktfoper+mtfnull)then
        call tfevallstkall(list,ref,ref,irtc)
      else
        isp=isp+1
        if(iand(lconstlist,list%attr) .ne. 0)then
          ktastk(isp)=ktflist+ksad_loc(list%head)
          irtc=0
        else
          call tfleval(list,ktastk(isp),ref,irtc)     
          if(irtc .ne. 0)then
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
      use tfstk
      implicit none
      type (sad_list) :: list
      type (sad_list), pointer :: kli
      integer*8 ki,iaat
      integer*4 isp0,mf,irtc,i,m
      logical*4 av
      irtc=0
      if(m .le. 0)then
        return
      endif
      if(av)then
        ktastk(isp+1:isp+m)=list%body(1:m)
        isp=isp+m
      else
        select case (iaat)
        case (1:)
          do i=1,m
            ki=list%body(i)
            if(ktflistq(ki,kli))then
              if(kli%head .eq. ktfoper+mtfnull)then
                call tfargevalstk(isp0,kli,kli%nl,iaat,mf,
     $               ktfreallistqo(kli),irtc)
                if(irtc .ne. 0)then
                  return
                endif
              elseif(ilist(2,iaat+min(isp-isp0+1,mf+1)) .eq. 0)then
c                call tfmemcheckprint('aevstk',.false.,irtc)
c                if(irtc .ne. 0)then
c                  call tfdebugprint(ki,'argevstk',1)
c                endif
                call tflevalstk(kli,.true.,irtc)
                if(irtc .ne. 0)then
                  return
                endif
              else
                isp=isp+1
                ktastk(isp)=ki
              endif
            elseif((ktfsymbolq(ki) .or. ktfpatq(ki)) .and.
     $             ilist(2,iaat+min(isp-isp0+1,mf+1)) .eq. 0)then
              call tfevalstk(ki,.true.,irtc)
              if(irtc .ne. 0)then
                return
              endif
            else
              isp=isp+1
              ktastk(isp)=ki
            endif
          enddo
        case (-iattrholdfirst)
          do i=1,m
            ki=list%body(i)
            if(ktfsequenceq(ki,kli))then
              call tfevallstkall(kli,isp .gt. isp0,.true.,irtc)
              if(irtc .ne. 0)then
                return
              endif
            elseif(isp .gt. isp0)then
              isp=isp+1
              call tfeevalref(ki,ktastk(isp),irtc)
              if(irtc .ne. 0)then
                return
              endif
            else
              isp=isp+1
              ktastk(isp)=ki
            endif
          enddo
        case (-iattrholdrest)
          do i=1,m
            ki=list%body(i)
            if(ktfsequenceq(ki,kli))then
              call tfevallstkall(kli,isp .eq. isp0,.false.,irtc)
              if(irtc .ne. 0)then
                return
              endif
            elseif(isp .eq. isp0)then
              isp=isp+1
              call tfeevalref(ki,ktastk(isp),irtc)
              if(irtc .ne. 0)then
                return
              endif
            else
              isp=isp+1
              ktastk(isp)=ki
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
      use tfstk
      implicit none
      integer*8 ks(m)
      integer*4 irtc,i,isp0,m
      logical*4 av
      irtc=0
      if(m .le. 0)then
        return
      endif
      if(av)then
        isp0=isp
        ktastk(isp0+1:isp0+m)=ks
        isp=isp0+m
      else
        do i=1,m
          call tfevalstk(ks(i),.true.,irtc)
          if(irtc .ne. 0)then
            return
          endif
        enddo
      endif
      return
      end

      subroutine tfseqevalstk(ks,m,i,av,irtc)
      use tfstk
      implicit none
      type (sad_list), pointer :: kli
      integer*8 ks(m)
      integer*4 irtc,m,i
      logical*4 av
      irtc=0
      if(i .le. m)then
        if(av)then
          isp=isp+1
          ktastk(isp)=ks(i)
        else
          if(ktflistq(ks(i),kli))then
            call tflevalstk(kli,.true.,irtc)
          elseif(ktfsymbolq(ks(i)) .or. ktfpatq(ks(i)))then
            call tfevalstk(ks(i),.true.,irtc)
          else
            isp=isp+1
            ktastk(isp)=ks(i)
          endif
        endif
      endif
      return
      end

      subroutine tfevallev(list,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_list) list
      type (sad_list), pointer :: listi
      integer*8 ki,kai,ki0
      integer*4 irtc,i,isp1
      logical*4 ev
      if(ktfreallistqo(list) .or. iand(list%attr,kconstarg) .ne. 0)then
        kx%k=ktflist+ksad_loc(list%head)
        irtc=0
        return
      endif
      ev=.false.
      isp=isp+1
      isp1=isp
      ktastk(isp)=list%head
      do i=1,list%nl
        ki=list%body(i)
        kai=ktfaddr(ki)
        select case(ki-kai)
        case (ktflist)
          call loc_list(kai,listi)
          if(iand(lconstlist,listi%attr) .eq. 0)then
            ki0=ki
            call tfleval(listi,ki,.true.,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            ev=ev .or. ki .ne. ki0
          endif
        case (ktfsymbol)
          ki0=ki
          call tfsyeval(ki0,ki,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          ev=ev .or. ki .ne. ki0
        case (ktfpat)
          ki0=ki
          call tfpateval(ki0,ki,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          ev=ev .or. ki .ne. ki0
        end select
        isp=isp+1
        ktastk(isp)=ki
        if(ktfsequenceq(ki,listi))then
          ev=.true.
          isp=isp-1
          call tfgetllstkall(listi)
        endif
      enddo
      if(ev)then
        kx=kxcompose(isp1)
      else
        kx%k=ktflist+ksad_loc(list%head)
      endif
      irtc=0
 9000 isp=isp1-1
      return
      end

      subroutine tfevalstk(k,ref,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      integer*4 irtc
      logical*4 ref
      if(ktfsequenceqd(k,kl))then
        call tfevallstkall(kl,ref,ref,irtc)
      else
        isp=isp+1
        call tfeeval(k,dtastk(isp),ref,irtc)
        if(irtc .ne. 0)then
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
      use tfstk
      implicit none
      type (sad_list) list
      integer*8 ki
      integer*4 irtc,i,isp0,m
      logical*4 ref,ref1
      m=list%nl
      if(ktfreallistqo(list))then
        isp0=isp
        ktastk(isp0+1:isp0+m)=list%body(1:m)
        isp=isp0+m
      elseif(m .gt. 0)then
        ki=list%body(1)
        if(ktfrealq(ki))then
          isp=isp+1
          ktastk(isp)=ki
        else
          call tfevalstk(ki,ref1,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
        do i=2,m
          ki=list%body(i)
          if(ktfrealq(ki))then
            isp=isp+1
            ktastk(isp)=ki
          else
            call tfevalstk(ki,ref,irtc)
            if(irtc .ne. 0)then
              return
            endif
          endif
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfgetllstk(list,i1,i2)
      use tfstk
      implicit none
      type (sad_list) list
      type (sad_list), pointer :: kl
      integer*4 i,m,i1,i2
      if(i2 .ge. 0)then
        m=min(i2,list%nl)
      else
        m=list%nl+i2+1
      endif
      if(i1 .gt. m)then
        return
      endif
      do i=max(0,i1),m
        isp=isp+1
        ktastk(isp)=list%body(i)
        if(ktfsequenceq(ktastk(isp),kl))then
          isp=isp-1
          call tfgetllstkall(kl)
        endif
      enddo
      return
      end

      recursive subroutine tfgetllstkall(list)
      use tfstk
      implicit none
      type (sad_list) list
      type (sad_list),pointer :: listi
      integer*4 i,m
      logical*4 noseq
      m=list%nl
      if(iand(list%attr,lnoseqlist) .ne. 0)then
        ktastk(isp+1:isp+m)=list%body(1:m)
c        call tmov(klist(ka+1),ktastk(isp+1),m)
        isp=isp+m
        return
      endif
      noseq=.true.
      if(ktfreallistqo(list))then
        ktastk(isp+1:isp+m)=list%body(1:m)
        isp=isp+m
      else
        do i=1,m
          isp=isp+1
          ktastk(isp)=list%body(i)
          if(ktfsequenceq(ktastk(isp)))then
            noseq=.false.
            isp=isp-1
            call loc_list(ktfaddr(list%body(i)),listi)
            call tfgetllstkall(listi)
          endif
        enddo
      endif
      if(noseq)then
        list%attr=ior(list%attr,lnoseqlist)
      endif
      return
      end

      subroutine tfsyeval(ka,kx,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      use mackw
      implicit none
      type (sad_descriptor) ka,kx
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      type (sad_object), pointer :: obj
      integer*4 maxlevel
      parameter (maxlevel=4096)
      integer*8 kas,j
      integer*4 irtc,itfmessage,level1
      integer*4 level,ig,ig1
      data level/0/
      kas=ktfaddrd(ka)
      call loc_sym(kas,sym)
      irtc=0
      ig=sym%gen
      if(ig .eq. -3)then
        kx%k=ktfsymbol+kas
        return
      endif
      if(sym%override .ne. 0)then
        call loc_symdef(kas,symd)
      else
        call loc_namtbl(sym%loc,loc)
        kas=loc%symdef
        ig=max(0,ig)
        do while(kas .gt. 0)
          call loc1_symdef(kas,symd)
          ig1=max(0,symd%sym%gen)
          if(ig .eq. ig1)then
            exit
          elseif(ig .lt. ig1)then
            kas=symd%next
          else
            exit
          endif
        enddo
      endif
      kx=symd%value
      if(kx%k .ne. ktfsymbol+kas .and. ktfnonrealq(kx%k))then
        level=level+1
        if(level .gt. maxlevel-64)then
          level1=level
          level=0
          irtc=itfmessage(999,'General::deep','""')
          level=level1
          return
        endif
        call tfeevalref(kx,kx,irtc)
        level=level-1
        if(ktfobjqd(kx,obj))then
          if(ktfaddr(obj%alloc) .eq. 0)then
            j=itflocal+levele
            obj%alloc=ktftype(obj%alloc)+klist(j)
            klist(j)=ksad_loc(obj%alloc)
          endif
        endif
      endif
      return
      end

      subroutine tfeevaldef(k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx
      type (sad_list), pointer :: list
      type (sad_symbol), pointer :: sym
      integer*4 irtc
      if(ktflistqd(k,list))then
        if(iand(lconstlist,list%attr) .eq. 0)then
          call tfleval(list,kx,.false.,irtc)
          return
        endif
      elseif(ktfsymbolqd(k,sym) .and. sym%override .eq. 0)then
        call tfsydef(sym,sym)
        kx%k=ktfsymbol+ksad_loc(sym%loc)
        irtc=0
        return
      endif
      irtc=0
      kx=k
      return
      end

      subroutine tfslot(kopc,kls,kx,ref,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ka
      type (sad_list) kls
      type (sad_symbol), pointer :: sym
      type (sad_namtbl), pointer :: nam
      integer*8 kaa,kopc,kax
      integer*4 ind,irtc,nc,isp1,isps,
     $     itfmessage,ns,ipf0,naf0,ls,isp2
      real*8 ffval,vx
      character*256 name
      character*12 inds
      logical*4 exist,ref
      ns=kls%nl
      if(ns .gt. 1)then
        irtc=itfmessage(9,'General::narg','"0 or 1"')
        return
      endif
      if(ns .eq. 0)then
        ind=1
      else
        ka=kls%dbody(1)
        if(ktfoperqd(ka,kaa))then
          if(kaa .eq. mtfnull)then
            ind=1
          else
            irtc=itfmessage(999,'General::invop',' ')
            return
          endif
        elseif(ktfrealqdi(ka,ind))then
        elseif(ktfsymbolqd(ka,sym) .and. kopc .eq. mtfslot)then
          call sym_namtbl(sym,nam)
          nc=nam%str%nch+1
          name(2:nc)=nam%str%str(1:nc-1)
          name(1:1)='#'
          call capita(name(1:nc))
          vx=ffval(name(1:nc),kax,exist)
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
        if(ind .lt. 0)then
          ind=napuref+ind+1
        endif
      endif
      isps=ipurefp+ind
      if(kopc .eq. mtfslot)then
        if(ind .le. 0 .or. ind .gt. napuref)then
          call strfromil(ind,inds,ls)
          irtc=itfmessage(999,'General::slot',
     $         '"#'//inds(:ls)//'"')
          return
        endif
        kx=dtastk(isps)
      else
        if(ipurefp .eq. 0 .or. ind .le. 0)then
          call strfromil(ind,inds,ls)
          irtc=itfmessage(999,'General::slot',
     $         '"##'//inds(:ls)//'"')
          return
        endif
        isp1=isp
        isp2=ipurefp+napuref
        call tfsequence(isps-1,isp2,kx)
      endif
      ipf0=ipurefp
      naf0=napuref
      ipurefp=itastk(1,ipf0+naf0+1)
      napuref=itastk(2,ipf0+naf0+1)
      call tfeeval(kx,kx,ref,irtc)
      ipurefp=ipf0
      napuref=naf0
      return
      end

      subroutine tfpateval(k,kx,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) k,kx,ke
      type (sad_pat), pointer :: pat
      integer*8 ka,kae,kte
      integer*4 irtc
      ka=ktfaddrd(k)
      call loc_pat(ka,pat)
      kx%k=ktfpat+ka
      irtc=0
      ke=pat%expr
      kae=ktfaddrd(ke)
      kte=ke%k-kae
      if(kte .ne. ktfref)then
c        call tfdebugprint(ke,'pateval',1)
c        write(*,*)'at ',kae
        call tfeevalref(ke,ke,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kae=ktfaddrd(ke)
        kte=ke%k-kae
      endif
      if(ke%k .ne. pat%expr%k)then
        kx=kxpcopyss(ke,pat%head,
     $       pat%sym%alloc,pat%default)
      endif
      return
      end

      subroutine tffindlabel(list,m,i,kr)
      use tfstk
      implicit none
      type (sad_descriptor) kr
      type (sad_list) list
      type (sad_list), pointer :: listj
      integer*4 i,j,m
      logical*4 tfsameqd,tfsameqk
      type (sad_descriptor), save :: kxlabel
      data kxlabel%k /0/
      if(kxlabel%k .eq. 0)then
        kxlabel=kxsymbolf('Label',5,.true.)
      endif
      do j=1,m
        if(ktflistq(list%body(j),listj))then
          if(listj%nl .eq. 1)then
            if(tfsameqk(kxlabel,listj%head))then
              if(tfsameqd(kr,listj%dbody(1)))then
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
