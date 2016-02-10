c$Header: /SAD/cvsroot/oldsad/src/tfeeval.f,v 1.151.2.21 2012/09/06 23:15:00 oide Exp $
      subroutine tfeeval(k,kx,ref,irtc)
      use tfstk
      implicit none
      integer*8 k,kx
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
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*8 k,kx,kt,ka
      integer*4 irtc
      ka=ktfaddr(k)
      kt=k-ka
      if(kt .eq. ktfsymbol)then
        call tfsyeval(ka,kx,irtc)
        return
      elseif(kt .eq. ktfpat)then
        call tfpateval(ka,kx,irtc)
        return
      elseif(kt .eq. ktflist)then
        ka=ktfaddr(k)
        if(iand(lconstlist,ilist(2,ka-3)) .eq. 0)then
          call tfleval(ka,kx,.true.,irtc)
          return
        endif
      elseif(kt .eq. ktfref)then
        kx=klist(ka)
        return
      endif
      irtc=0
      kx=k
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfleval(k,kx,ref,irtc)
      use tfstk
      implicit none
      integer*8 k,kx,kaa,ka
      integer*4 irtc
      logical*4 ref,tfconstlistqk
c      include 'DEBUG.inc'
      ka=ktfaddr(k)
      kx=ktflist+ka
      irtc=0
      if(iand(lconstlist,ilist(2,ka-3)) .ne. 0)then
        return
      endif
      if(iand(lnoconstlist,ilist(2,ka-3)) .eq. 0)then
        if(tfconstlistqk(ka))then
          return
        endif
      endif
      kaa=ka
      ilist(1,kaa-1)=ilist(1,kaa-1)+1
      call tfseval(kaa,ilist(2,kaa-1),klist(kaa),kx,.false.,
     $     ktfreallistq(ka),ref,irtc)
      call tflocal1(kaa)
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfseval(ks,ns,kh,kx,stk,av,ref,irtc)
      use tfstk
      implicit none
      integer*4 maxlevel
      parameter (maxlevel=4096)
      integer*8 ks,kh,kx,kaa,kaf,ktf,kf,kl,ktx,iaat,
     $     ktfloadlstk,ktfstk2l
      integer*4 irtc,ns
      integer*4 isp1,isp10,isp11,level,itfgetrecl,
     $     j,isp0,mf,i1,level1,i,itfmessage,lpw,l,itfdownlevel
      logical*4 ref,ev,evalh,rep,stk,tfconstlistqk,tfgetseqstk,av
      data level/0/
c      include 'DEBUG.inc'
      level=level+1
      if(level .ge. maxlevel-64)then
        level1=level
        level=0
        irtc=itfmessage(999,'General::deep','""')
        level=level1
        go to 9000
      endif
      isp0=isp
      kf=kh
      evalh=.false.
      kaf=ktfaddr(kf)
      ktf=kf-kaf
      if(ktf .eq. ktfoper)then
        if(ilist(1,klist(ifunbase+kaf)+1) .lt. 0)then
          go to 100
        endif
      elseif(ktf .eq. ktfsymbol)then
        if(ref)then
          call tfsyeval(kaf,kf,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        else
          call tfsydef(kaf,kf)
        endif
      elseif(ktf .eq. ktflist)then
        if(iand(lconstlist,ilist(2,kaf-3)) .eq. 0)then
          call tfleval(kaf,kf,ref,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
      elseif(ktf .eq. ktfpat .and. ref)then
        call tfpateval(kaf,kf,irtc)
        if(irtc .ne. 0)then
          go to 8000
        endif
      endif
      ev=.true.
      iaat=0

      do while(.true.)
        kaf=ktfaddr(kf)
        ktf=kf-kaf
        if(ktf .eq. ktfoper)then
          iaat=klist(ifunbase+kaf)+1
          mf=ilist(1,iaat)
          if(mf .lt. 0)then
            evalh=.true.
            exit
          endif
          ev=ilist(2,iaat+1) .eq. -2
          if(ilist(2,iaat+1) .eq. -1)then
            iaat=0
          endif
        elseif(ktf .eq. ktfsymbol)then
          if(ilist(2,kaf-3) .ne. 0)then
            iaat=iand(iattrholdall,ilist(1,kaf-3))
            if(iaat .ne. 0)then
              ev=.false.
              if(iaat .eq. iattrholdall)then
                iaat=0
              else
                iaat=-iaat
              endif
            endif
          endif
        elseif(ktf .eq. ktflist)then
          if(klist(kaf) .eq. ktfoper+mtfnull)then
            if(ilist(2,kaf-1) .eq. 1)then
              kf=klist(kaf+1)
              cycle
            elseif(ilist(2,kaf-1) .eq. 0)then
              kf=ktfoper+mtfnull
              evalh=.true.
              exit
            endif
          elseif(ktfsymbolq(klist(kaf)) .and. .not. ref)then
            if(ilist(2,iand(ktamask,klist(kaf))) .eq. -3 .and.
     $           iand(ilist(2,kaf-3),lnonreallist) .eq. 0)then
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
          kaa=ktfloadlstk(ks)
          klist(kaa)=ktfoper+kaf
          kaa=ktfstk2l(kaa)
          if(ref .and. tfconstlistqk(kaa))then
            kx=ktflist+kaa
            irtc=0
            go to 8000
          endif
        else
          kaa=ks
        endif
        if(ref)then
          call tfevallev(kaa,kx,irtc)
        else
          call tfevallstkall(kaa,.false.,.false.,irtc)
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
        call tfeevaldef(klist(ks+1),ktastk(isp),irtc)
        if(irtc .ne. 0)then
          go to 8000
        endif
        call tfseqevalstkall(ks+1,ns-1,av,irtc)
        if(irtc .eq. 0)then
          levele=levele+1
          go to 6000
        endif
      case (mtfslot,mtfslotseq)
        call tfslot(kaf,ks,ns,kx,ref,irtc)
      case (mtfcomp)
        if(ns .eq. 0)then
          kx=ktfoper+mtfnull
          go to 8000
        endif
        i1=1
        do while(.true.)
          if(ltrace .gt. 0)then
            levele=levele+1
            lpw=min(131,itfgetrecl())
            do i=i1,ns
              call tfprint1(klist(ks+i),6,-lpw,4,.true.,
     $             .true.,irtc)
              call tfeevalref(klist(ks+i),kx,irtc)
              if(irtc .ne. 0)then
                go to 1320
              endif
            enddo
          else
            levele=levele+1
            irtc=0
            do i=i1,ns
              kx=klist(ks+i)
              ktx=ktftype(kx)
              if(ktx .eq. ktflist)then
                call tfleval(kx,kx,.true.,irtc)
              elseif(ktx .eq. ktfsymbol)then
                call tfsyeval(kx,kx,irtc)
              elseif(ktx .eq. ktfpat)then
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
          call tffindlabel(ks,ns,i1,kl)
          if(i1 .le. 0)then
            call tfthrow(2,kl,irtc)
            exit
          endif
          i1=i1+1
        enddo
      case (mtfand)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(ks,ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j)))then
              if(ktastk(j) .eq. 0)then
                kx=0
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx=ktftrue
        elseif(isp .eq. isp1+1)then
          kx=ktastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfor)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(ks,ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j)))then
              if(ktastk(j) .ne. 0)then
                kx=ktftrue
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx=0
        elseif(isp .eq. isp1+1)then
          kx=ktastk(isp)
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
        kx=ktflist+ks
      case default
        write(*,*)'tfleval-implementation error: ',kaf
        stop
      end select
      go to 8000

 3000 if(levele .ge. maxlevele-32)then
        irtc=itfmessage(999,'General::deep','""')
        go to 8000
      endif
      levele=levele+1
      if(ev)then
        call tfseqevalstkall(ks,ns,av,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      elseif(iaat .eq. 0 .or. av)then
        rep=tfgetseqstk(ks,ns)
      else
        call tfargevalstk(isp1,ks,ns,iaat,mf,.false.,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      endif
 6000 ktastk(isp1)=kf
      if(ref)then
c        if(alo .and. abs(klist(idb)) .gt. 100000000 .or.
c     $       idb2+klist(idb2+6) .gt. idb1)then
c          call tfdebugprint(kf,'tfseval-efunref-in',3)
c          call tfdebugprint(ktastk(isp1+1),' ',3)
c          call tfdebugprint(ktastk(isp),' ',3)
c        endif
        call tfefunref(isp1,kx,.true.,irtc)
c        if(alo .and. abs(klist(idb)) .gt. 100000000 .or.
c     $       idb2+klist(idb2+6) .gt. idb1)then
c          call tfdebugprint(kf,'tfseval-efunref-out',3)
c          call tfdebugprint(ktastk(isp1+1),' ',3)
c          call tfdebugprint(ktastk(isp),' ',3)
c        endif
      else
        call tfefundef(isp1,kx,irtc)
      endif
 7000 continue
      call tfconnectk(kx,irtc)
 8000 isp=isp0
 9000 level=max(0,level-1)
      return
      include 'inc/TFSF.inc'
      end

      logical*4 function tfgetseqstk(ks,ns)
      use tfstk
      implicit none
      integer*8 ks,ki
      integer*4 ns,i
      tfgetseqstk=.false.
      if(ns .gt. 0)then
        do i=1,ns
          ki=klist(ks+i)
          if(ktfsequenceq(ki))then
            tfgetseqstk=.true.
            call tfgetllstkall(ki)
            cycle
          endif
          isp=isp+1
          ktastk(isp)=ki
        enddo
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfconnectk(k,irtc)
      use tfstk
      implicit none
      integer*4 l,itfdownlevel,irtc
      integer*8 k,ka
      if(levele .gt. 0)then
        if(irtc .ne. 0)then
          l=itfdownlevel()
        elseif(ktfobjq(k))then
          ka=ktfaddr(k)
          ilist(1,ka-1)=ilist(1,ka-1)+1
          l=itfdownlevel()
          call tflocal1(ka)
c          call tfdebugprint(ktftype(klist(ka-2))+ka,'tfconnectk',1)
c          write(*,*)'with ',ilist(1,ka-1),ktfaddr(klist(ka-2))
        else
          l=itfdownlevel()
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      integer*4 function itfdownlevel()
      use tfstk
      implicit none
      integer*8 j,kt,ka,i,lastfree
      integer*4 n,kk
      lastfree=0
      j=itflocal+levele
      i=klist(j)
      do while(i .ne. j)
        ka=ktfaddr(klist(i))
        if(ilist(1,i+1) .le. 0)then
          kt=klist(i)-ka
          if(kt .eq. ktflist)then
            if(iand(ilist(2,i-1),lnonreallist) .ne. 0)then
              n=ilist(2,i+1)
              do kk=n+2,3,-1
                call tflocalrel(klist(i+kk),ka)
              enddo
            endif
            call tflocalrel(klist(i+2),ka)
          elseif(kt .eq. ktfpat)then
            call tflocalrel(klist(i+10),ka)
            if(ilist(1,i+8) .gt. 1)then
              call tflocalrel(ktfsymbol+i+9,ka)
              ilist(1,i+6)=ilist(1,i-1)-7
              ilist(2,i+6)=0
              ilist(1,i-1)=7
            endif
            if(klist(i+9) .ne. 0)then
              call tflocalrel(klist(i+7),ka)
            endif
            klist(i+7)=ktfsymbol
            call tflocalrel(klist(i+3),ka)
            call tflocalrel(klist(i+2),ka)
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
      include 'inc/TFSF.inc'
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

      subroutine tflocalrel(k,kan)
      use tfstk
      implicit none
      integer*8 k,ka,kan
c      include 'DEBUG.inc'
      if(ktfobjq(k))then
        ka=ktfaddr(k)
        ilist(1,ka-1)=max(0,ilist(1,ka-1)-1)
        if(ilist(1,ka-1) .eq. 0)then
          if(ktfaddr(klist(ka-2)) .eq. 0)then
            klist(ka-2)=ktftype(klist(ka-2))+kan
            kan=ka-2
          endif
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tflocal(k)
      use tfstk
      implicit none
c      include 'DEBUG.inc'
      integer*8 k,ka,itfroot
      if(ktfobjq(k))then
        ka=ktfaddr(k)
        ilist(1,ka-1)=max(0,ilist(1,ka-1)-1)
        if(ilist(1,ka-1) .eq. 0)then
          if(ktfaddr(klist(ka-2)) .eq. 0)then
            itfroot=itflocal+levele
            klist(ka-2)=klist(ka-2)+ktfaddr(klist(itfroot))
            klist(itfroot)=ka-2
          endif
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tflocal1(k)
      use tfstk
      implicit none
      integer*8 ka,k,itfroot
c      include 'DEBUG.inc'
      ka=ktfaddr(k)
      ilist(1,ka-1)=max(0,ilist(1,ka-1)-1)
      if(ilist(1,ka-1) .eq. 0)then
        if(ktfaddr(klist(ka-2)) .eq. 0)then
          itfroot=itflocal+levele
          klist(ka-2)=klist(ka-2)+ktfaddr(klist(itfroot))
          klist(itfroot)=ka-2
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      integer*4 function itfuplevel()
      use tfstk
      implicit none
      levele=min(maxlevele,levele+1)
      klist(itflocal+levele)=itflocal+levele
      itfuplevel=levele
      return
      end

      subroutine tflevalstk(ka,ref,irtc)
      use tfstk
      implicit none
      integer*8 ka
      integer*4 irtc
      logical*4 ref
      if(klist(ka) .eq. ktfoper+mtfnull)then
        call tfevallstkall(ka,ref,ref,irtc)
      else
        isp=isp+1
        if(iand(lconstlist,ilist(2,ka-3)) .ne. 0)then
          ktastk(isp)=ktflist+ka
          irtc=0
        else
          call tfleval(ka,ktastk(isp),ref,irtc)          
          if(irtc .ne. 0)then
            return
          endif
          if(ktfsequenceq(ktastk(isp)))then
            isp=isp-1
            call tfgetllstkall(ktastk(isp+1))
          endif
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      recursive subroutine tfargevalstk(isp0,ka,m,iaat,mf,av,irtc)
      use tfstk
      implicit none
c      include 'DEBUG.inc'
      integer*8 ka,ki,kti,kai,iaat
      integer*4 isp0,mf,irtc,i,m
      logical*4 av
      irtc=0
      if(m .le. 0)then
        return
      endif
      if(av)then
        ktastk(isp+1:isp+m)=klist(ka+1:ka+m)
        isp=isp+m
      elseif(iaat .gt. 0)then
        do i=1,m
          ki=klist(ka+i)
          kti=ktftype(ki)
          if(kti .eq. ktflist)then
            kai=ki-kti
            if(klist(kai) .eq. ktfoper+mtfnull)then
              call tfargevalstk(isp0,kai,ilist(2,kai-1),iaat,mf,
     $             ktfreallistq(kai),irtc)
              if(irtc .ne. 0)then
                return
              endif
            elseif(ilist(2,iaat+min(isp-isp0+1,mf+1)) .eq. 0)then
              call tflevalstk(kai,.true.,irtc)
              if(irtc .ne. 0)then
                return
              endif
            else
              isp=isp+1
              ktastk(isp)=ki
            endif
          elseif((kti .eq. ktfsymbol .or. kti .eq. ktfpat) .and.
     $           ilist(2,iaat+min(isp-isp0+1,mf+1)) .eq. 0)then
            call tfevalstk(ki,.true.,irtc)
            if(irtc .ne. 0)then
              return
            endif
          else
            isp=isp+1
            ktastk(isp)=ki
          endif
        enddo
      elseif(iaat .eq. -iattrholdfirst)then
        do i=1,m
          ki=klist(ka+i)
          kti=ktftype(ki)
          kai=ki-kti
          if(kti .eq. ktflist .and. klist(kai) .eq. ktfoper+mtfnull)then
            call tfevallstkall(kai,isp .gt. isp0,.true.,irtc)
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
      elseif(iaat .eq. -iattrholdrest)then
        do i=1,m
          ki=klist(ka+i)
          kti=ktftype(ki)
          kai=ki-kti
          if(kti .eq. ktflist .and. klist(kai) .eq. ktfoper+mtfnull)then
            call tfevallstkall(kai,isp .eq. isp0,.false.,irtc)
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
      elseif(iaat .eq. -iattrholdall)then
        call tfgetllstkall(ka)
      else
        call tfevallstkall(ka,.true.,.true.,irtc)
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfseqevalstkall(ka,m,av,irtc)
      use tfstk
      implicit none
      integer*8 ka
      integer*4 irtc,i,isp0,m
      logical*4 av
      irtc=0
      if(m .le. 0)then
        return
      endif
      if(av)then
        isp0=isp
        ktastk(isp0+1:isp0+m)=klist(ka+1:ka+m)
        isp=isp0+m
      else
        do i=1,m
          call tfevalstk(klist(ka+i),.true.,irtc)
          if(irtc .ne. 0)then
            return
          endif
        enddo
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfseqevalstk(ka,m,i,ev,irtc)
      use tfstk
      implicit none
      integer*8 ka,kti
      integer*4 irtc,m,i
      logical*4 ev
      irtc=0
      if(i .le. m)then
        if(ev)then
          isp=isp+1
          rtastk(isp)=rlist(ka+i)
        else
          kti=ktftype(klist(ka+i))
          if(kti .eq. ktflist)then
            call tflevalstk(iand(ktamask,klist(ka+i)),.true.,irtc)
          elseif(kti .eq. ktfsymbol .or. kti .eq. ktfpat)then
            call tfevalstk(klist(ka+i),.true.,irtc)
          else
            isp=isp+1
            ktastk(isp)=klist(ka+i)
          endif
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfevallev(ka0,kx,irtc)
      use tfstk
      implicit none
      integer*8 ka,kx,ki,kti,kai,kax,ki0,kai0,ktfcompose,ka0
      integer*4 irtc,i,isp1
      logical*4 ev
      ka=ktfaddr(ka0)
      if(ktfreallistq(ka) .or. iand(ilist(2,ka-3),kconstarg) .ne. 0)then
        kx=ktflist+ka
        irtc=0
        return
      endif
      ev=.false.
      isp=isp+1
      isp1=isp
      ktastk(isp)=klist(ka)
      do i=1,ilist(2,ka-1)
        ki=klist(ka+i)
        kai=ktfaddr(ki)
        kti=ki-kai
        if(kti .eq. ktflist)then
          ki0=ki
          if(iand(lconstlist,ilist(2,kai-3)) .eq. 0)then
            kai0=kai
            call tfleval(kai0,ki,.true.,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            ev=ev .or. ki .ne. ki0
          endif
        elseif(kti .eq. ktfsymbol)then
          ki0=ki
          call tfsyeval(ki0,ki,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          ev=ev .or. ki .ne. ki0
        elseif(kti .eq. ktfpat)then
          ki0=ki
          call tfpateval(ki0,ki,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          ev=ev .or. ki .ne. ki0
        elseif(iand(ktrmask,ki) .ne. ktfnr)then
          ki=klist(ka+i)
        endif
        isp=isp+1
        ktastk(isp)=ki
        if(ktfsequenceq(ki))then
          ev=.true.
          isp=isp-1
          call tfgetllstkall(ki)
        endif
      enddo
      if(ev)then
        kax=ktfcompose(isp1)
      else
        kax=ka
      endif
      kx=ktflist+kax
      irtc=0
 9000 isp=isp1-1
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfevalstk(k,ref,irtc)
      use tfstk
      implicit none
      integer*8 k
      integer*4 irtc
      logical*4 ref
      if(ktfsequenceq(k))then
        call tfevallstkall(ktfaddr(k),ref,ref,irtc)
      else
        isp=isp+1
        call tfeeval(k,ktastk(isp),ref,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfsequenceq(ktastk(isp)))then
          isp=isp-1
          call tfgetllstkall(ktastk(isp+1))
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfevallstkall(ka,ref1,ref,irtc)
      use tfstk
      implicit none
      integer*8 ka,ki
      integer*4 irtc,i,isp0,m
      logical*4 ref,ref1
      m=ilist(2,ka-1)
      if(ktfreallistq(ka))then
        isp0=isp
        ktastk(isp0+1:isp0+m)=klist(ka+1:ka+m)
        isp=isp0+m
      elseif(m .gt. 0)then
        ki=klist(ka+1)
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
          ki=klist(ka+i)
          if(iand(ktrmask,ki) .ne. ktfnr)then
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
      include 'inc/TFSF.inc'
      end

      subroutine tfgetlstk(ka,i)
      use tfstk
      implicit none
      integer*8 ka
      integer*4 i
      isp=isp+1
      if(i .ge. 0 .and. i .le. ilist(2,ka-1))then
        ktastk(isp)=klist(ka+i)
      else
        ktastk(isp)=kxnull
      endif
      if(ktfsequenceq(ktastk(isp)))then
        isp=isp-1
        call tfgetllstkall(ktastk(isp+1))
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfgetl(ka,i,kx)
      use tfstk
      implicit none
      integer*8 ka,kx
      integer*4 i
      kx=klist(iand(ktamask,ka)+i)
      return
      end

      subroutine tfgetllstk(k,i1,i2)
      use tfstk
      implicit none
      integer*8 k,ka
      integer*4 i,m,i1,i2
      ka=ktfaddr(k)
      if(i2 .ge. 0)then
        m=min(i2,ilist(2,ka-1))
      else
        m=ilist(2,ka-1)+i2+1
      endif
      if(i1 .gt. m)then
        return
      endif
      do i=max(0,i1),m
        call tfgetlstk(ka,i)
      enddo
      return
      include 'inc/TFSF.inc'
      end

      recursive subroutine tfgetllstkall(ka0)
      use tfstk
      implicit none
      integer*8 ka0,ka
      integer*4 i,m
      logical*4 noseq
      ka=ktfaddr(ka0)
      m=ilist(2,ka-1)
      if(iand(ilist(2,ka-3),lnoseqlist) .ne. 0)then
        klist(isp+1:isp+m)=klist(ka+1:ka+m)
c        call tmov(klist(ka+1),ktastk(isp+1),m)
        isp=isp+m
        return
      endif
      noseq=.true.
      if(iand(ilist(2,ka-3),lnonreallist) .eq. 0)then
        ktastk(isp+1:isp+m)=klist(ka+1:ka+m)
        isp=isp+m
      else
        do i=1,m
          isp=isp+1
          ktastk(isp)=klist(ka+i)
          if(ktfsequenceq(ktastk(isp)))then
            noseq=.false.
            isp=isp-1
            call tfgetllstkall(ktastk(isp+1))
          endif
        enddo
      endif
      if(noseq)then
        ilist(2,ka-3)=ior(ilist(2,ka-3),lnoseqlist)
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfsyeval(ka,kx,irtc)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 maxlevel
      parameter (maxlevel=4096)
      integer*8 ka,kx,kas,kas1,loc,kp
      integer*4 irtc,itfmessage,level1
      integer*4 level,ig,ig1
      data level/0/
      kas=ktfaddr(ka)
      irtc=0
      ig=ilist(2,kas-1)
      if(ig .eq. -3)then
        kx=ktfsymbol+kas
        return
      endif
      if(ilist(2,kas-3) .ne. 0)then
        kas1=kas
      else
        loc=klist(kas)-5
        kp=klist(loc+1)
        ig=max(0,ig)
        do while(kp .gt. 0)
          ig1=max(0,ilist(2,kp+7))
          if(ig .eq. ig1)then
            kas1=kp+8
            exit
          elseif(ig .lt. ig1)then
            kp=klist(kp)
          else
            kas1=kas
            exit
          endif
        enddo
      endif
      kx=klist(kas1-4)
c      call tfdebugprint(ktfstring+klist(kas),'tfsyeval',1)
c      call tfdebugprint(kx,'==> ',1)
      if(kx .ne. ktfsymbol+kas1 .and. ktfnonrealq(kx))then
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
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfeevaldef(k,kx,irtc)
      use tfstk
      implicit none
      integer*8 k,kx,ka,kt
      integer*4 irtc
      ka=ktfaddr(k)
      kt=k-ka
      if(kt .eq. ktflist)then
        if(iand(lconstlist,ilist(2,ka-3)) .eq. 0)then
          call tfleval(ka,kx,.false.,irtc)
          return
        endif
      elseif(kt .eq. ktfsymbol)then
        if(ilist(2,ka-3) .eq. 0)then
          call tfsydef(ka,kx)
          irtc=0
          return
        endif
      endif
      irtc=0
      kx=k
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfslot(kopc,k1,ns,kx,ref,irtc)
      use tfstk
      implicit none
      integer*8 ka,kx,k1,kaa,kta,kopc,kfromr,kax,kas
      integer*4 ind,irtc,nc,isp1,isps,
     $     itfmessage,ns,ipf0,naf0,ls,isp2
      real*8 ffval,vx,rfromk
      character*256 name
      character*12 inds
      logical*4 exist,ref
      if(ns .gt. 1)then
        irtc=itfmessage(9,'General::narg','"0 or 1"')
        return
      endif
      if(ns .eq. 0)then
        ind=1
      else
        ka=klist(k1+1)
        kaa=ktfaddr(ka)
        kta=ka-kaa
        if(kta .eq. ktfoper)then
          if(kaa .eq. mtfnull)then
            ind=1
          else
            irtc=itfmessage(999,'General::invop',' ')
            return
          endif
        elseif(ktfrealq(ka))then
          ind=int(rfromk(ka))
        elseif(kta .eq. ktfsymbol .and. kopc .eq. mtfslot)then
          kas=klist(kaa)
          nc=ilist(1,kaa)+1
          call tmovb(ilist(1,kas+1),name(2:),nc-1)
          name(1:1)='#'
          call capita(name(1:nc))
          vx=ffval(name(1:nc),kax,exist)
          if(exist)then
            kx=kfromr(vx)
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
        kx=ktastk(isps)
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
      include 'inc/TFSF.inc'
      end

      subroutine tfpateval(k,kx,irtc)
      use tfstk
      implicit none
      integer*8 k,kx,ka,ke,kae,kte,kd,ktpcopyss
      integer*4 irtc
      ka=ktfaddr(k)
      kx=ktfpat+ka
      irtc=0
      ke=klist(ka)
      kae=ktfaddr(ke)
      kte=ke-kae
      if(kte .ne. ktfref)then
        call tfeevalref(ke,ke,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kae=ktfaddr(ke)
        kte=ke-kae
      endif
      if(ke .ne. klist(ka))then
        kd=klist(ka+8)
        kx=ktfpat+ktpcopyss(ke,klist(ka+1),ktfaddr(klist(ka+5)),kd)
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tffindlabel(ka,m,i,kr)
      use tfstk
      implicit none
      integer*8 ka,kr,kl,kaj,ktfsymbolf
      integer*4 i,j,m
      logical*4 tfsameqk
      integer*8 itflabel
      data itflabel /0/
      if(itflabel .eq. 0)then
        itflabel=ktfsymbolf('Label',5,.true.)
      endif
      if(ka .ne. 0)then
        do j=1,m
          if(ktftype(klist(ka+j)) .eq. ktflist)then
            kaj=iand(ktamask,klist(ka+j))
            if(ilist(2,kaj-1) .eq. 1)then
              if(tfsameqk(ktfsymbol+itflabel,klist(kaj)))then
                kl=klist(kaj+1)
                if(tfsameqk(kr,kl))then
                  i=j
                  return
                endif
              endif
            endif
          endif
        enddo
      endif
      i=0
      return
      include 'inc/TFSF.inc'
      end
