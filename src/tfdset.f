      module dset
      use tfstk

      contains
      subroutine tfcleardaloc(kad)
      use tfcode
      implicit none
      type (sad_deftbl), pointer :: dtbl
      integer*8 ,intent(in):: kad
      call loc_deftbl(kad,dtbl)
      call tflocal1d(dtbl%arg)
      call tflocal1d(dtbl%argc)
      call tflocald(dtbl%body)
      call tflocald(dtbl%bodyc)
      return
      end

      subroutine tfreplacedef(k,karg,kr,kargr,irtc)
      use repl,only:tfreplacesymbolstk1,tfreplacesymbolstk
      use pmat,only:tfinitpat
      implicit none
      type (sad_descriptor) ,intent(in):: k,karg
      type (sad_descriptor) ,intent(inout):: kr,kargr
      type (sad_descriptor) kx1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,ispr,i,ipsbase,nrule
      logical*4 rep,same
      isp0=isp
      call tfinitpat(isp0,karg)
      if(isp == isp0)then
        kr=k
        kargr=karg
        irtc=0
        return
      endif
      ispr=isp
      ipsbase=0
      same=.true.
      do i=isp0+1,ispr-1,2
        ipsbase=ipsbase+1
        isp=isp+2
        ktastk(isp-1)=ktfsymbol+ktastk(i+1)
        ivstk2(2,isp-1)=0
        dtastk(isp  )=kxargsym(ipsbase)
        same=same .and. ktastk(i+1) == klist(ktfaddr(ktastk(isp)))
      enddo
c      call tfresetpat(karg)
      if(same)then
        isp=isp0
        kr=k
        kargr=karg
        irtc=0
        return
      endif
      nrule=(isp-ispr)/2
      kx1=tfreplacesymbolstk(karg,ispr,nrule,.false.,rep,irtc)
      if(irtc /= 0)then
        isp=isp0
        return
      endif
      kargr=kx1
      kr=tfreplacesymbolstk1(k,ispr,nrule,.false.,rep,irtc)
      isp=isp0
      return
      end

      integer*8 function ktdaloc(kan0,kb0,kn0,k,karg,kr,kargr,tbl)
      use tfcode
      use pmat,only:tfinitpatlist
      implicit none
      type (sad_descriptor) ,intent(in):: karg,k
      type (sad_descriptor) ,intent(in):: kr,kargr
      type (sad_descriptor) kx,karg1
      type (sad_dlist), pointer :: klx
      type (sad_deftbl), pointer :: dtbl
      integer*8 ,intent(in):: kan0,kb0,kn0
      integer*8 kan,ktalocr,kb,kn
      integer*4 npat,isp0
      logical*4 ,intent(in):: tbl
      logical*4 rep,sym
      isp0=isp
      if(tbl)then
        kx=tfdefsymbol(kargr,rep,sym)
        karg1=kx
        call descr_sad(karg1,klx)
        call tfinitpatlist(isp0,klx)
        npat=(isp-isp0)/2
      else
        karg1=karg
        npat=0
      endif
      kan=kan0
      kn=kn0
      kb=kb0
      if(kan == 0)then
        kan=ktalocr(npat+8)
        call loc_deftbl(kan,dtbl)
        dtbl%next=kn
        dtbl%prev=kb
        klist(kb)=kan
        if(kn /= 0)then
          klist(kn+1)=kan
        endif
      else
        call loc_deftbl(kan,dtbl)
        call tfcleardaloc(kan)
        if(npat > dtbl%npat)then
          write(*,*)'ktdaloc-implementation error ',npat,dtbl%npat
          call abort
        endif
      endif
      dtbl%npat=npat
      dtbl%attr=0
      dtbl%arg=dtfcopy1(karg)   ! orifinal arg
      dtbl%argc=dtfcopy1(karg1)   ! replaced arg
      dtbl%body=dtfcopy(k)       ! original body
      dtbl%pat=-1
      if(.not. ktfobjq(kr))then
        dtbl%bodyc=kr
      else
        dtbl%bodyc=dtfcopy1(kr)    ! replaced body
        if(.not. tfconstq(kr%k))then
          dtbl%pat=0
        endif
      endif
      dtbl%compile=0.d0   ! compile level
      dtbl%pattbl(1:npat)=dtastk(isp0+1:isp0+npat*2-1:2)
      isp=isp0
      ktdaloc=kan
      return
      end
      
      recursive function tfdefsymbol(k,rep,sym) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: list
      type (sad_rlist), pointer :: klr
      type (sad_symbol), pointer :: ks,ks1
      integer*4 i,isp0,m
      logical*4 ,intent(out):: rep,sym
      logical*4 rep1
      rep=.false.
      sym=.false.
      kx=k
      if(ktflistq(k,list))then
        if(iand(lnodefsymbol,list%attr) /= 0)then
          return
        endif
        isp0=isp
        isp=isp+1
        dtastk(isp)=tfdefsymbol(list%head,rep,sym)
        m=list%nl
        if(ktfreallistq(list))then
          if(rep)then
            kx=kxavaloc(-1,m,klr)
            klr%head=dtfcopy(dtastk(isp))
            call tfcopyarray(list%rbody(1:m),klr%rbody(1:m),m)
c            klr%rbody(1:m)=list%rbody(1:m)
          elseif(.not. sym)then
            list%attr=ior(lnodefsymbol,list%attr)
          endif
          isp=isp0
          return
        endif
        do i=1,list%nl
          isp=isp+1
          dtastk(isp)=tfdefsymbol(list%dbody(i),rep1,sym)
          rep=rep .or. rep1
        enddo
        if(rep)then
          kx=kxcomposer(isp0+1)
        elseif(.not. sym)then
          list%attr=ior(lnodefsymbol,list%attr)
        endif
        isp=isp0
        return
      elseif(ktfsymbolq(k,ks))then
        if(ks%override == 0)then
          sym=.true.
          if(ks%gen /= maxgeneration)then
            ks1=>tfsydef(ks)
            kx=sad_descr(ks1)
            rep=.true.
          endif
        endif
      endif
      return
      end

      integer*8 function ktdhtaloc(kad1,kad0,nhash)
      use tfcode
      implicit none
      type (sad_defhash), pointer :: dhash
      integer*8 ,intent(in):: kad1,kad0
      integer*8 kan,kad,ktalocr
      integer*4 ,intent(in):: nhash
      kad=kad0
      kan=ktalocr(nhash+5)
      call loc_defhash(kan,dhash)
      dhash%attr=0
      dhash%next=kad
      dhash%prev=kad1
      dhash%gen=maxgeneration
      dhash%nhash=nhash
      dhash%dhash(0:nhash)%k=0
      klist(kad1 )=kan
      if(kad /= 0)then
        klist(kad+1)=kan
      endif
      ktdhtaloc=kan
      return
      end

      integer*4 function itfhasharg(k,nh)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      integer*8 ka
      integer*4 ,intent(in):: nh
      integer*4 ih,m,ix(2),ih1
      integer*2 h(2),h1(2)
      real*8 x,ha
      parameter (ha=7.d0**5+1.d0/7.d0**5)
      equivalence (h,ih),(x,ix),(h1,ih1)
      ih=0
      if(ktflistq(k,kl))then
        m=kl%nl
        ih=m
        if(m == 1)then
          ih=itfhash(kl%dbody(1),nh)+1
        elseif(m /= 0)then
          ih=ih+itfhash(kl%dbody(1),nh)+itfhash(kl%dbody(m),nh)
        endif
      elseif(ktfrefq(k,ka))then
        if(ka == 0)then
          ih=0
        else
          ka=ka-ispbase
          m=itastk2(1,ka)-int(ka)
          if(m == 0)then
            ih=m
          elseif(m == 1)then
            ih=itfhash(dtastk(ka+1),nh)+m
          else
            ih=itfhash(dtastk(ka+1),nh)+itfhash(dtastk(ka+m),nh)+m
          endif
        endif
      else
        call tfdebugprint(k,'itfhasharg-implementation error: ',3)
        rlist(7)=0.d0
      endif
      ih=(h(1)+h(2))
      itfhasharg=iand(ih,nh)
      return
      end

      recursive integer*4 function itfhash(k,nh) result(ih)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 ,intent(in):: nh
      integer*4 m,ix(2),ih1
      integer*2 h(2)
      real*8 x,ha,v
      parameter (ha=7.d0**5+1.d0/7.d0**5)
      equivalence (h,ih1),(x,ix)
      if(ktfrealq(k,v))then
        ih=int(v*ha)
      elseif(ktfoperq(k,ka))then
        ih=int(ka)
      elseif(ktflistq(k,kl))then
        m=kl%nl
        ih=itfhash(kl%head,nh)+m
        if(m == 1)then
          ih=ih+itfhash1(kl%dbody(1))
        elseif(m /= 0)then
          ih=ih+itfhash1(kl%dbody(1))+itfhash1(kl%dbody(m))
        endif
      elseif(ktfrefq(k,ka))then
        if(ka == 0)then
          ih=0
        else
          ka=ka-ispbase
          m=itastk2(1,ka)-int(ka)
          ih=itfhash(dtastk(ka),nh)+m
          if(m == 1)then
            ih=ih+itfhash1(dtastk(ka+1))
          elseif(m /= 0)then
            ih=ih+itfhash1(dtastk(ka+1))+itfhash1(dtastk(ka+m))
          endif
        endif
      elseif(ktfstringq(k))then
        ka=ktfaddrd(k)
        ih=ilist(1,ka)+ilist(1,ka+1)+ilist(2,ka+1)
      elseif(ktfsymbolq(k,sym))then
        ih=int(sym%loc)
      else
        ih=0
      endif
      ih1=ih
      ih=(h(1)+h(2))
      ih=iand(ih,nh)
      return
      end

      recursive integer*4 function itfhash1(k) result(ih)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 ix(2)
      real*8 x,v
      real*8 ,parameter ::ha=7.d0**5+1.d0/7.d0**5
      equivalence (x,ix)
      if(ktfrealq(k,v))then
        ih=int(v*ha)
      elseif(ktfoperq(k,ka))then
        ih=int(ka)
      elseif(ktflistq(k,kl))then
        ih=itfhash1(kl%head)+kl%nl
      elseif(ktfrefq(k,ka))then
        if(ka == 0)then
          ih=0
        else
          ka=ka-ispbase
          ih=itfhash1(dtastk(ka))+itastk2(1,ka)-int(ka)
        endif
      elseif(ktfstringq(k))then
        ka=ktfaddrd(k)
        ih=ilist(1,ka)+ilist(1,ka+1)+ilist(2,ka+1)
      elseif(ktfsymbolq(k,sym))then
        ih=int(sym%loc)
      else
        ih=0
      endif
      return
      end

      subroutine tfdeval1(isp1,kad0,kx,ev,ordless,ns,irtc)
      use tfcode
      use eeval
      use pmat,only:itfpmat,itfseqmatstk,itfseqmatstk1
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kh,kal
      type (sad_dlist), pointer :: larg
      type (sad_deftbl), pointer :: dtbl
      integer*8 ,intent(in):: kad0
      integer*8 kad,kap
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc,ns
      integer*4 mstk0,mat,im,isp0,iop0,irs
      logical*4 ,intent(out):: ev
      logical*4 ,intent(in):: ordless
      isp0=isp
      mstk0=mstk
      iop0=iordless
      iordless=0
      im=isp1+1
      ns=0
      kad=kad0
      do while(kad /= 0)
        mat=-1
        call loc_deftbl(kad,dtbl)
        kap=ktfaddr(dtbl%argc)
        call loc_sad(kap,larg)
        kh=larg%head
        mat=itfpmat(dtastk(isp1),kh)
        if(mat >= 0)then
          if(larg%nl == 1)then
            if(itfseqmatstk1(im,isp0,larg%dbody(1)) < 0)then
              go to 30
            endif
          else
            irs=itfseqmatstk(im,isp0,larg%dbody(1),
     $           larg%nl,1,ktfreallistq(larg),merge(kap,i00,ordless))
            if(irs < 0)then
              go to 30
            endif
          endif
          ev=.true.
          if(dtbl%compile < rlist(levelcompile))then
            call tfsetarg(dtbl,irtc)
            if(irtc /= 0)then
              call tfunsetpattbl(dtbl)
              isp=isp0
              iordless=iop0
              mstk=mstk0
              return
            endif
          endif
          iordless=iop0
          mstk=mstk0
          if(dtbl%npat /= 0)then
            kal=dtfcopy1(dtbl%argc)
            kx=tfevalwitharg(dtbl,dtbl%bodyc,irtc)
            call tflocal1d(kal)
          else
            kx=tfeevalref(dtbl%bodyc,irtc)
          endif
          if(irtc /= 0)then
            call tfcatchreturn(irtcret,kx,irtc)
          endif
          isp=isp0
          return
 30       call tfunsetpattbl(dtbl)
        endif
        iordless=iop0
        mstk=mstk0
        isp=isp0
        kad=dtbl%next
        ns=ns+1
      enddo
      irtc=-1
      return
      end

      function tfevalwitharg(dtbl,k,irtc) result(kx)
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kh
      type (sad_deftbl) ,intent(in):: dtbl
      integer*8 ks
      integer*4 ,intent(out):: irtc
      integer*4 isp0,m
      logical*4 av
      isp0=isp
      kx=tfevalwitharg1(k,ks,m,kh,av,irtc)
      call tfunsetpattbl(dtbl)
      if(irtc /= 0)then
        return
      endif
      if(kx%k == ktfref)then
        call tfseval(ks,m,kh,kx,.true.,av,.true.,irtc)
      else
        kx=tfeevalref(kx,irtc)
      endif
      call tfcatchreturn(irtcret,kx,irtc)
      isp=isp0
      return
      end

      function tfevalwitharg1(k,ks,m,kh,av,irtc) result(kx)
      use tfcode
      use funs,only:tfgetstkstk
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(out):: kh
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) :: kp(0:2),ksy,tfcomposefun,tfcomposeoper
      type (sad_dlist), pointer :: list,klh
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symdef
      integer*8 ,intent(out):: ks
      integer*8 i,kts,kah,ks1
      integer*4 ,intent(out):: irtc,m
      integer*4 isp1,itfmessageexp,isp0
      logical*4 rep,rep1
      logical*4 ,intent(out):: av
      irtc=0
      kx=k
      if(ktflistq(k,list))then
        kx%k=ktfref
        m=list%nl
        isp1=isp
        kh=list%head
        if(ktfreallistq(list))then
          if(ktflistq(kh,klh))then
            if(iand(larglist,klh%attr) == 0)then
              kx=k
              return
            endif
          elseif(ktfnonsymbolq(kh) .and. .not. ktfpatq(kh))then
            kx=k
            return
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc /= 0)then
            return
          endif
          ks=ktfaddrd(k)
          av=.true.
          return
        else
          if(ktflistq(kh,klh))then
            if(iand(larglist,klh%attr) == 0)then
              go to 20
            endif
          elseif(ktfnonsymbolq(kh) .and. .not. ktfpatq(kh))then
            go to 20
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc /= 0)then
            return
          endif
 20       if(ktfoperq(kh))then
            kah=ktfaddrd(kh)
            if(iget_fun_id(ktfaddrd(kh)) == nfunif)then
              call tfreplaceifarg(list,kx,rep,irtc)
              return
            endif
          endif
          if(iand(list%attr,lnopatlist) /= 0)then
            ks=ktfaddrd(k)
            av=.false.
            return
          endif
          isp=isp1+1
          dtastk(isp)=kh
          ks=isp+ispbase
          do i=1,m
            rep1=tfreplaceargstk(list%dbody(i),irtc)
            if(irtc /= 0)then
              return
            endif
          enddo
          if(rlist(iaximmediate) /= 0.d0)then
            if(ktfoperq(kh))then
              if(kah > mtfend)then
                ks1=isp+ispbase
                kx=tfcomposefun(isp1+1,kah,.false.,irtc)
              else
                kx=tfcomposeoper(isp1+1,kah,.false.,isp0,irtc)
                ks1=isp0+ispbase
              endif
              if(irtc == 0)then
                if(kx%k /= ktfref)then
                  return
                endif
                ks=ks1
              elseif(irtc > 1)then
                return
              else
                kx%k=ktfref
              endif
              irtc=0
            endif
          endif
          m=int(isp+ispbase-ks)
          av=.true.
          do i=ks-ispbase+1,isp
            if(ktfnonrealq(ktastk(i)))then
              av=.false.
              exit
            endif
          enddo
          return
        endif
      elseif(ktfpatq(k,pat))then
        rep=.false.
        if(pat%sym%loc /= 0)then
          rep=tfreplacearg(pat%sym%alloc,ksy,irtc)
          if(irtc /= 0)then
            return
          endif
          if(rep .and. ktfnonsymbolq(ksy))then
            irtc=itfmessageexp(999,'General::reppat',ksy)
            return
          endif
        else
          ksy%k=0
        endif
        rep1=tfreplacearg(pat%expr,kp(0),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%head,kp(1),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%default,kp(2),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        if(rep)then
          kx=kxpcopyss(kp(0),kp(1),ksy,kp(2))
        endif
      elseif(ktfsymbolqdef(k%k,symdef))then
        if(symdef%sym%override == 1)then
          kx=symdef%value
          kts=ktftype(kx%k)
          if(kts == ktflist .or. kts == ktfref)then
            call tfgetstkstk(kx,rep)
          endif
        endif
      endif
      return
      end

      logical*4 function tfreplacearg(k,kx,irtc)
      use funs
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: irtc
      integer*4 isp0
      isp0=isp
      tfreplacearg=tfreplaceargstk(k,irtc)
      if(irtc /= 0)then
        isp=isp0
        return
      endif
      if(isp == isp0)then
        kx=dxnull
      elseif(isp == isp0+1)then
        kx=dtastk(isp)
      else
        kx=tfsequence(isp0,isp)
      endif
      isp=isp0
      return
      end

      recursive logical*4 function tfreplaceargstk(k,irtc)
     $     result(lx)
      use tfcode
      use funs,only:tfgetstkstk
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ks,kx,kh,kp(0:2),tfcomposefun,tfcomposeoper
      type (sad_dlist), pointer :: list,klx,klh
      type (sad_rlist), pointer :: klr
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symd
      integer*8 kah
      integer*4 ,intent(out):: irtc
      integer*4 i,m,isp1,itfmessageexp,isp0
      logical*4 rep,rep1
      irtc=0
      rep=.false.
      lx=.false.
      if(ktflistq(k,list))then
        kx=k
        isp1=isp
        kh=list%head
        if(ktfreallistq(list))then
          if(ktflistq(kh,klh))then
            if(iand(larglist,klh%attr) == 0)then
              go to 10
            endif
          elseif(ktfnonsymbolq(kh) .and. .not. ktfpatq(kh))then
            go to 10
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc /= 0)then
            return
          endif
          if(rep)then
            m=list%nl
            kx=kxavaloc(-1,m,klr)
            call tfcopyarray(list%rbody(1:m),klr%rbody(1:m),m)
c            klr%rbody(1:m)=list%rbody(1:m)
            klr%head=dtfcopy(kh)
          endif
 10       isp=isp1+1
          dtastk(isp)=kx
          return
        else
          if(ktflistq(kh,klh))then
            if(iand(larglist,klh%attr) == 0)then
              go to 20
            endif
          elseif(ktfnonsymbolq(kh) .and. .not. ktfpatq(kh))then
            go to 20
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc /= 0)then
            return
          endif
 20       isp=isp1+1
          dtastk(isp)=kh
          if(ktfoperq(kh,kah))then
            if(iget_fun_id(kah) == nfunif)then
              call tfreplaceifarg(list,kx,rep,irtc)
              isp=isp1+1
              dtastk(isp)=kx
              lx=rep
              return
            elseif(kah == mtfnull)then
              isp=isp1
            endif
          endif
          if(iand(list%attr,lnopatlist) /= 0)then
            if(rep)then
              call tfgetllstkall(list)
            elseif(kh%k == ktfoper+mtfnull)then
              call tfgetllstkall(list)
              lx=.true.
              return
            else
              isp=isp1+1
              dtastk(isp)=k
              return
            endif
          else
            do i=1,list%nl
              rep1=tfreplaceargstk(list%dbody(i),irtc)
              if(irtc /= 0)then
                return
              endif
              rep=rep .or. rep1
            enddo
          endif
          lx=rep
          if(kh%k == ktfoper+mtfnull)then
            return
          endif
          if(rlist(iaximmediate) /= 0.d0)then
            if(ktfoperq(kh))then
              if(kah > mtfend)then
                kx=tfcomposefun(isp1+1,kah,.false.,irtc)
              else
                kx=tfcomposeoper(isp1+1,kah,.true.,isp0,irtc)
              endif
              if(irtc == 0)then
                lx=.true.
                if(ktflistq(kx,klx))then
                  if(klx%head%k == ktfoper+mtfnull)then
                    isp=isp1
                    call tfgetllstkall(klx)
                    return
                  endif
                endif
                isp=isp1+1
                dtastk(isp)=kx
                return
              elseif(irtc > 1)then
                return
              endif
            elseif(ktfstringq(kh))then
              if(isp == isp1+1 .or. isp == isp1+2)then
                if(ktfrealq(ktastk(isp1+1)) .and.
     $               ktfrealq(ktastk(isp)))then
                  kx=kxsubstring(kh,isp1+1,isp)
                  lx=.true.
                  irtc=0
                  return
                endif
              endif
            endif
          endif
          dtastk(isp1+1)=kxcrelistm(isp-isp1-1,
     $         ktastk(isp1+2:isp),kh)
          isp=isp1+1
          irtc=0
          lx=.true.
          return
        endif
      elseif(ktfpatq(k,pat))then
        rep=.false.
        kx=k
        if(pat%sym%loc /= 0)then
          rep=tfreplacearg(pat%sym%alloc,ks,irtc)
          if(irtc /= 0)then
            return
          endif
          if(rep .and. ktfnonsymbolq(ks))then
            irtc=itfmessageexp(999,'General::reppat',ks)
            return
          endif
        else
          ks%k=0
        endif
        rep1=tfreplacearg(pat%expr,kp(0),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%head,kp(1),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%default,kp(2),irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
c          if(ktftype(kp(i)) == ktfref)then
c            call tfsequence(kp(i),kp(i))
c            rep=.true.
c          endif
        lx=rep
        if(rep)then
          kx=kxpcopyss(kp(0),kp(1),ks,kp(2))
        endif
        isp=isp+1
        dtastk(isp)=kx
      elseif(ktfsymbolqdef(k%k,symd))then
        if(symd%sym%override == 1)then
          ks=symd%value
          if(ktflistq(ks) .or. ktfrefq(ks))then
            call tfgetstkstk(ks,rep)
          else
            isp=isp+1
            dtastk(isp)=ks
          endif
          rep=.true.
        else
          isp=isp+1
          dtastk(isp)=k
        endif
        lx=rep
        return
      else
        isp=isp+1
        dtastk(isp)=k
      endif
      return
      end

      subroutine tfreplaceifarg(list,kx,rep,irtc)
      implicit none
      type (sad_descriptor) kx,ki,tfcompose
      type (sad_dlist) list
      type (sad_dlist), pointer :: klx
      integer*4 irtc,i,isp1,j
      logical*4 rep,rep1
      irtc=0
      isp1=isp
      if(list%nl == 0)then
        if(rep)then
          kx=kxaaloc(-1,0,klx)
          klx%head=dtastk(isp1)
        else
          kx=sad_descr(list)
        endif
        return
      endif
      ki=list%dbody(1)
      select case(ktftype(ki%k))
      case (ktflist)
        if(iand(larglist,list%attr) /= 0)then
          rep1=tfreplaceargstk(ki,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
          rep=rep .or. rep1
        else
          isp=isp+1
          dtastk(isp)=ki
        endif
      case (ktfsymbol,ktfpat)
        rep1=tfreplaceargstk(ki,irtc)
        if(irtc /= 0)then
          isp=isp1
          return
        endif
        rep=rep .or. rep1
      case default
        isp=isp+1
        dtastk(isp)=ki
      end select
      i=0
      if(isp == isp1+1)then
        if(ktfrealq(ktastk(isp1+1)))then
          if(ktastk(isp1+1) /= 0)then
            i=2
          else
            i=3
          endif
        endif
      else
        rep=.true.
      endif
      if(i /= 0)then
        rep=.true.
        j=isp+1
        if(i .le. list%nl)then
          ki=list%dbody(i)
          select case (ktftype(ki%k))
          case (ktflist)
            if(iand(larglist,list%attr) /= 0)then
              rep1=tfreplaceargstk(ki,irtc)
              if(irtc /= 0)then
                isp=isp1
                return
              endif
              kx=merge(dtastk(j),dxnullo,j .le. isp)
            else
              kx=ki
            endif
          case (ktfsymbol,ktfpat)
            rep1=tfreplaceargstk(ki,irtc)
            if(irtc /= 0)then
              isp=isp1
              return
            endif
            kx=merge(dtastk(j),dxnullo,j .le. isp)
          case default
            kx=ki
          end select
        else
          kx=dxnullo
        endif
        return
      endif
      do i=2,list%nl
        ki=list%dbody(i)
        select case(ktftype(ki%k))
        case (ktflist)
          if(iand(larglist,list%attr) /= 0)then
            rep1=tfreplaceargstk(ki,irtc)
            if(irtc /= 0)then
              isp=isp1
              return
            endif
            rep=rep .or. rep1
          else
            isp=isp+1
            dtastk(isp)=ki
          endif
        case (ktfsymbol,ktfpat)
          rep1=tfreplaceargstk(ki,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
          rep=rep .or. rep1
        case default
          isp=isp+1
          dtastk(isp)=ki
        end select
      enddo
      if(rep)then
        kx=tfcompose(isp1,list%head%k,irtc)
      else
        kx=sad_descr(list)
      endif
      return
      end

      subroutine tfunsetpattbl(dtbl)
      use tfcode
      implicit none
      type (sad_deftbl) dtbl
      type (sad_pat), pointer :: pat
      integer*4 np,i
      np=dtbl%npat
      if(np /= maxgeneration)then
        do i=1,np
          call descr_sad(dtbl%pattbl(i),pat)
          pat%mat=0
          pat%value%k=ktfref
        enddo
      endif
      return
      end

      subroutine tfsetarg(dtbl,irtc)
      use tfcode
      use repl,only:tfsortsymbolstk
      use sameq,only:tfsymbollistqo
      implicit none
      type (sad_descriptor) kx
      type (sad_deftbl) dtbl
      type (sad_dlist), pointer :: klx
      type (sad_pat), pointer :: pat
      integer*4 isp0,i,irtc,nrule1
      logical*4 rep,member
      irtc=0
      kx=dtbl%bodyc
      if(ktflistq(kx,klx))then
        if(.not. tfsymbollistqo(klx))then
          dtbl%compile=1.d100
          return
        endif
      elseif(tfconstq(kx%k))then
        dtbl%compile=1.d100
        return
      endif
      if(dtbl%compile == 0.d0)then
        if(dtbl%npat /= 0)then
          isp0=isp
          do i=1,dtbl%npat
            call descr_pat(dtbl%pattbl(i),pat)
            do while(associated(pat%equiv))
              pat=>pat%equiv
            enddo
            isp=isp+2
            ktastk(isp-1)=ktfaddr(pat%sym%alloc%k)
            dtastk(isp  )=sad_descr(pat%sym)
          enddo
          call tfsortsymbolstk(isp0,dtbl%npat,nrule1)
          kx=tfsetuparg(dtbl%bodyc,isp0,nrule1*2,rep,member,irtc)
          isp=isp0
          if(irtc /= 0)then
            return
          endif
        else
          member=.true.
        endif
        if(member)then
          dtbl%compile=rlist(levelcompile)
          if(ktflistq(kx,klx))then
            if(iand(lmemberlist,klx%attr) /= 0)then
              kx=tfrecompilearg(kx,rep,irtc)
              if(irtc /= 0)then
                return
              endif
            endif
          elseif(ktfpatq(kx))then
            kx=tfrecompilearg(kx,rep,irtc)
            if(irtc /= 0)then
              return
            endif
          endif
        else
          dtbl%compile=1.d100
        endif
      else
        if(ktflistq(kx,klx))then
          if(iand(lmemberlist,klx%attr) /= 0)then
            kx=tfrecompilearg(kx,rep,irtc)
            if(irtc /= 0)then
              return
            endif
            dtbl%compile=rlist(levelcompile)
          endif
        elseif(ktfpatq(kx))then
          kx=tfrecompilearg(kx,rep,irtc)
          if(irtc /= 0)then
            return
          endif
          dtbl%compile=rlist(levelcompile)
        endif
      endif
      call tflocald(dtbl%bodyc)
      dtbl%bodyc=dtfcopy(kx)
      return
      end

      recursive function tfsetuparg(k,ispr,nrule2,rep,member,irtc)
     $     result(kx)
      use tfcode
      use repl,only:tfmatchsymstk
      implicit none
      type (sad_descriptor) k,kx,ki,k1,ks,kd
      type (sad_dlist), pointer :: list,klx
      type (sad_rlist), pointer :: klr
      type (sad_pat), pointer :: pat
      integer*8 ka1
      integer*4 irtc,i,m,isp1,ispr,nrule2,itfmessageexp,j
      logical*4 rep,rep1,rep2,member,member1
      irtc=0
      member=.false.
      rep=.false.
      kx=k
      if(ktfsymbolq(k%k))then
        if(tfmatchsymstk(ktfaddr(k%k),ispr,nrule2/2,j))then
          kx=dtastk(ivstk2(1,j)+1)
          rep=.true.
        endif
        return
      elseif(ktflistq(k,list))then
        k1=tfsetuparg(list%head,ispr,nrule2,rep,member,irtc)
        if(irtc /= 0)then
          return
        endif
        m=list%nl
        if(rep .and. ktfreallistq(list))then
          kx=kxavaloc(-1,m,klr)
          call tfcopyarray(list%rbody(1:m),klr%rbody(1:m),m)
c          klr%rbody(1:m)=list%rbody(1:m)
          klr%attr=ior(larglist,list%attr)
          klr%head=dtfcopy(k1)
          if(member)then
            klr%attr=ior(lmemberlist,klr%attr)
          else
            klr%attr=iand(-lmemberlist-1,klr%attr)
          endif
c          klr%attr=merge(ior(lmemberlist,klr%attr),
c     $         iand(-lmemberlist-1,klr%attr),member)
          return
        endif
        isp1=isp
        isp=isp+1
        dtastk(isp)=k1
        rep2=.false.
        LOOP_I: do i=1,m
          ki=tfsetuparg(list%dbody(i),
     $       ispr,nrule2,rep1,member1,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
          rep2=rep2 .or. rep1
          member=member .or. member1
          isp=isp+1
          dtastk(isp)=ki
        enddo LOOP_I
        if(rep2)then
          kx=kxcomposev(isp1+1,klx)
          klx%attr=ior(larglist,list%attr)
          rep=.true.
        else
          list%attr=ior(list%attr,kpatarg)
          if(rep)then
            kx=kxcomposev(isp1+1,klx)
            klx%attr=ior(larglist,list%attr)
          else
            kx=k
            klx=>list
            if(list%head%k == ktfoper+mtfatt)then
              member=.true.
            endif
          endif
        endif
        if(member)then
          klx%attr=ior(lmemberlist,klx%attr)
        else
          klx%attr=iand(-lmemberlist-1,klx%attr)
        endif
c        klx%attr=merge(ior(lmemberlist,klx%attr),
c     $       iand(-lmemberlist-1,klx%attr),member)
        isp=isp1
      elseif(ktfpatq(k,pat))then
        rep=.false.
        if(pat%sym%loc /= 0)then
          ks=tfsetuparg(pat%sym%alloc,ispr,nrule2,rep,member,irtc)
          if(irtc /= 0)then
            return
          endif
          if(rep)then
            if(ktftype(ks%k) /= ktfsymbol)then
              irtc=itfmessageexp(999,'General::reppat',k)
              return
            endif
          endif
        else
          ks%k=0
        endif
        k1=pat%expr
        if(.not. ktfrefq(k1,ka1) .or. ka1 > 3)then
          k1=tfsetuparg(k1,ispr,nrule2,rep1,member1,irtc)
          member=member .or. member1
          if(irtc /= 0)then
            return
          endif
          rep=rep .or. rep1
        endif
        kd=pat%default
        if(ktftype(kd%k) /= ktfref)then
          kd=tfsetuparg(kd,ispr,nrule2,rep1,member1,irtc)
          member=member .or. member1
          if(irtc /= 0)then
            return
          endif
          rep=rep .or. rep1
        endif
        if(rep)then
          kx=kxpcopyss(k1,pat%head,ks,kd)
        endif
      endif
      return
      end

      recursive function tfrecompilearg(k,rep,irtc) result(kx)
      use tfcode
      implicit none
      type (sad_descriptor) kx,k,k1,k2,kd
      type (sad_dlist), pointer :: list,klx
      type (sad_rlist), pointer :: klr
      type (sad_pat), pointer :: pat
      type (sad_symbol), pointer :: sym2
      integer*8 ka1
      integer*4 irtc,i,m,isp1
      logical*4 rep,rep1,rep2
      irtc=0
      rep=.false.
      kx=k
      if(ktflistq(k,list))then
        if(iand(list%attr,lmemberlist) == 0)then
          return
        endif
        k1=list%head
        if(k1%k == ktfoper+mtfhold)then
          return
        endif
        k1=tfrecompilearg(k1,rep,irtc)
        if(ktfreallistq(list))then
          if(rep)then
            m=list%nl
            kx=kxavaloc(-1,m,klr)
            call tfcopyarray(list%rbody(1:m),klr%rbody(1:m),m)
c            klr%rbody(1:m)=list%rbody(1:m)
c            call tmov(rlist(ka+1),rlist(kax+1),m)
            klr%attr=ior(larglist,list%attr)
            klr%head=dtfcopy(k1)
          endif
          return
        endif
        isp1=isp
        isp=isp+1
        dtastk(isp)=k1
        rep2=.false.
        do i=1,list%nl
          isp=isp+1
          dtastk(isp)=tfrecompilearg(list%dbody(i),rep1,irtc)
          if(irtc /= 0)then
            isp=isp1
            return
          endif
          rep2=rep2 .or. rep1
        enddo
        if(list%head%k == ktfoper+mtfatt)then
          if(isp == isp1+3)then
            k2=dtastk(isp1+2)
            if(ktfsymbolq(k2,sym2))then
              if(sym2%override /= 0)then
                if(iand(sym2%attr,iattrdynamic) /= 0)then
                  go to 120
                endif
              endif
c              call tfdebugprint(ktastk(isp1+1),'rcmparg',3)
              call tfatt(isp1+1,kx,.false.,irtc)
              if(irtc > 0)then
                isp=isp1
                return
              elseif(irtc == 0)then
c                call tfdebugprint(kx,'==>',3)
                isp=isp1
                rep=.true.
                return
              endif
              irtc=0
            endif
          endif
        endif
 120    if(rep .or. rep2)then
          kx%k=ktflist+ktfcompose(isp1+1,klx)
          klx%attr=ior(larglist,list%attr)
          rep=.true.
        endif
        isp=isp1
      elseif(ktfpatq(k,pat))then
        k1=pat%expr
        if(ktfrefq(k1,ka1) .and. ka1 > 3)then
          k1=tfrecompilearg(k1,rep,irtc)
          if(irtc /= 0)then
            return
          endif
        endif
        kd=pat%default
        kd=tfrecompilearg(kd,rep1,irtc)
        if(irtc /= 0)then
          return
        endif
        rep=rep .or. rep1
        if(rep)then
          kx=kxpcopyss(k1,pat%head,pat%sym%alloc,kd)
        endif
      endif
      return
      end function

      integer*8 function ktfrehash(kad0,iup)
      use tfcode
      implicit none
      type (sad_descriptor) karg,ka
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash0,dhash
      integer*8 kad,kad0,ka0,kas,i
      integer*4 iup,n,n0,ih
      integer*4 maxhash
      parameter (maxhash=32767)
      call loc_defhash(kad0,dhash0)
      n0=dhash0%nhash
      if(n0 >= maxhash)then
        ktfrehash=kad0
        return
      endif
      n=n0*2+1
      kas=dhash0%prev-iup+6
      kad=ktdhtaloc(dhash0%prev,dhash0%next,n)
      call loc_defhash(kad,dhash)
      dhash%attr=dhash0%attr
      do i=0,dhash0%nhash
        ka=dhash0%dhash(i)
        do while(ka%k /= 0)
          call loc_deftbl(ka%k,dtbl)
          karg=dtbl%arg
          ih=itfhasharg(karg,n)
          ka0=dtbl%next
          dtbl%next=dhash%dhash(ih)%k
          if(dtbl%next /= 0)then
            dlist(dtbl%next+1)=ka
          endif
          dtbl%prev=sad_loc(dhash%dhash(ih))
          dhash%dhash(ih)=ka
          ka%k=ka0
        enddo
      enddo
      call tfree(kad0)
      ktfrehash=kad
      return
      end

      logical*4 function tfmaxgenerationq(kad)
      implicit none
      integer*8 ,intent(in):: kad
      if(kad /= 0)then
        tfmaxgenerationq=ilist(1,kad+2) == maxgeneration
      else
        tfmaxgenerationq=.false.
      endif
      return
      end

      function tfupset(k1,k2,kas,irtc) result(kx)
      use funs,only:tfgetllstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,ki,karg
      type (sad_dlist), pointer :: kl,kli
      type (sad_symbol), pointer :: symi
      type (sad_symdef), pointer :: symd
      integer*8 ,intent(in):: kas
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0,isp1,m,itfmessage
      kx=dxnullo
      if(ktfnonlistq(k1,kl))then
        irtc=itfmessage(999,'General::wrongtype','"Expression"')
        return
      endif
      m=kl%nl
      if(m .le. 0)then
        irtc=itfmessage(999,'General::wrongleng',
     $       '"Expression","longer than 0"')
        return
      endif
      isp0=isp
      isp1=isp0+1
      call tfgetllstk(kl,0,-1)
      karg=kxcompose(isp1)
      LOOP_I: do i=isp1+1,isp
        ki=dtastk(i)
        do while(ktflistq(ki,kli))
          ki=kli%head
        enddo
        if(ktfsymbolqdef(ki%k,symd))then
          if(symd%sym%override /= 0)then
            if(symd%sym%gen .lt. 0 .and. symd%sym%gen /= -3)then
              cycle LOOP_I
            endif
            if(kas == 0 .or. kas == ktfaddr(ki))then
              call tfdset(k2,symd%upval,kx,karg)
              if(kas /= 0)then
                cycle LOOP_I
              endif
            endif
          else
            symi=>tfsydef(symd%sym)
            if(symi%gen .lt. 0 .and. symi%gen /= -3)then
              cycle LOOP_I
            endif
            if(kas == 0 .or. kas == ksad_loc(symi%loc))then
              call sym_symdef(symi,symd)
              call tfdset(k2,symd%upval,kx,karg)
              if(kas /= 0)then
                cycle LOOP_I
              endif
            endif
          endif
        endif
      enddo LOOP_I
      kx=k2
      isp=isp0
      irtc=0
      return
      end

      end module dset

      subroutine tfdset(k,kadi,kx,karg)
      use dset
      use pmat,only:itflistmat
      use tfcode
      implicit none
      type (sad_descriptor) ,intent(in):: k,karg
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) :: kr,kargr
      type (sad_dlist), pointer :: larg,largl,largd,largr
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash
      integer*8 ,intent(in):: kadi
      integer*8 kap,kad,kad1,kad0,kad10,kan,kih
      integer*4 m,ihash,isp0,mstk0,iop0,irtc,minhash,maxhash
      parameter (minhash=7,maxhash=32767)

c     Initialize to avoid compiler warning
      kad0=0
c
      call descr_sad(karg,larg)
      kad1=ksad_loc(kadi)
      kad=kadi
      if(tfconstpatternq(karg%k))then
        if(.not. tfmaxgenerationq(kad))then
          if(k%k /= ktfref)then
            ihash=itfhasharg(karg,minhash)
            kih=ktdhtaloc(kad1,kad,minhash)
            call loc_defhash(kih,dhash)
            kad1=sad_loc(dhash%dhash(ihash))
            kap=ktdaloc(i00,kad1,i00,sad_descr(ktfref),karg,k,karg,.false.)
            dhash%attr=ior(dhash%attr,1)
            kx=k
          else
            kx=dxnullo
          endif
          return
        else
          ihash=itfhasharg(karg,ilist(2,kad+2))
          call loc_defhash(kad,dhash)
          dhash%attr=ior(dhash%attr,1)
          kad1=sad_loc(dhash%dhash(ihash))
          kad=klist(kad1)
          do while(kad /= 0)
            call loc_deftbl(kad,dtbl)
            call descr_sad(dtbl%arg,largd)
            if(larg%nl == largd%nl)then
              if(tfsamelistqo(larg,largd))then
                if(k%k == ktfref)then
                  go to 8000
                else
                  kan=ktdaloc(kad,i00,i00,
     $                 sad_descr(ktfref),karg,k,karg,.false.)
                  kx=k
                  return
                endif
              endif
            endif
            kad1=kad
            kad=dtbl%next
          enddo
          if(k%k /= ktfref)then
            kad=klist(kad1)
            kan=ktdaloc(i00,
     $           kad1,kad,sad_descr(ktfref),karg,k,karg,.false.)
            kx=k
          else
            kx=dxnullo
          endif
          return
        endif
      endif
      if(tfmaxgenerationq(kad))then
        kad1=kad
        kad=klist(kad)
      endif
      kad10=0
      isp0=isp
      call tfreplacedef(k,karg,kr,kargr,irtc)
      if(irtc /= 0)then
        return
      endif
      call descr_sad(kargr,largr)
      if(rlist(iaxpriority) == 0.d0)then
        mstk0=mstk
        iop0=iordless
        iordless=0
        do while(kad > 0)
          call loc_deftbl(kad,dtbl)
          call descr_sad(dtbl%argc,largl)
          m=itflistmat(kargr,largl)
          isp=isp0
          if(m >= 0)then
            call tfunsetpattbl(dtbl)
          endif
 1        if(m == 1)then
            if(k%k == ktfref)then
              mstk=mstk0
              iordless=iop0
              go to 8000
            else
              kan=ktdaloc(kad,i00,i00,k,karg,kr,kargr,.true.)
            endif
            iordless=iop0
            mstk=mstk0
            kx=k
            return
          else
            if(largr%nl == largl%nl)then
              if(tfsamelistqo(largr,largl))then
                m=1
                go to 1
              endif
            endif
            if(m == 0)then
              if(kad10 == 0)then
                kad10=kad1
                kad0=kad
              endif
              exit
            endif
            kad1=kad
            kad=dtbl%next
          endif
c$$$          npat=dtbl%npat
c$$$          do i=1,npat
c$$$            kpa=dtbl%pattbl(i)
c$$$            call descr_sad(kpa,pat)
c$$$            call tfdebugprint(kpa,'tfdset-tbl-1',1)
c$$$            write(*,*)': ',i,m,ktfaddr(kpa),pat%ref
c$$$            if(associated(pat%equiv))then
c$$$              call tfdebugprint(sad_descr(pat%equiv),'ktdaloc-equiv',1)
c$$$            endif
c$$$          enddo
        enddo
        mstk=mstk0
        iordless=iop0
      endif
      if(k%k /= ktfref)then
        if(kad10 /= 0)then
          kad1=kad10
          kad=kad0
        endif
        kan=ktdaloc(i00,kad1,kad,k,karg,kr,kargr,.true.)
        kx=k
      else
        kx=dxnullo
      endif
      return
 8000 klist(kad1)=klist(kad)
      if(klist(kad) /= 0)then
        klist(klist(kad)+1)=kad1
      endif
      call tfcleardaloc(kad)
      call tfree(kad)
      kx%k=ktfoper+mtfnull
      return
      end

      subroutine tfdeval(isp1,kad00,kx,iup,def,ev,irtc)
      use dset
      use tfcode
      use eeval
      use pmat,only:itfseqmatstk,tfordlessq
      implicit none
      type (sad_descriptor) kx,kad00
      type (sad_dlist), pointer :: larg,klx
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash
      integer*8 kad0,kad,kadv,kadv0,kap,khash
      integer*4 ,intent(in):: isp1,iup
      integer*4 ,intent(out):: irtc
      integer*4 is,m,mstk0,im,i,isp0,ns,nh,iord0
      logical*4 ,intent(in):: def
      logical*4 ,intent(out):: ev
      logical*4 ordless
      ev=.false.
      im=isp1+1
      ordless=iup == 0 .and. tfordlessq(dtastk(isp1)) .and. isp > im
      kad0=ktfaddr(kad00)
      kad=klist(iup+kad0-6)
      isp0=isp
      if(tfmaxgenerationq(kad))then
        call loc_defhash(kad,dhash)
        if(iand(dhash%attr,1) == 0)then
          go to 20
        endif
        nh=dhash%nhash
        itastk2(1,isp1)=isp
        khash=sad_loc(dhash%dhash(
     $       itfhasharg(dfromk(ktfref+isp1+ispbase),nh)))
        kadv=klist(khash)
        if(kadv /= 0)then
          is=1
          kadv0=kadv
          ns=0
          do while(kadv /= 0)
            call loc_deftbl(kadv,dtbl)
            kap=ktfaddr(dtbl%arg)
            call loc_sad(kap,larg)
            m=larg%nl
            if(m == isp-isp1)then
              if(.not. tfsameq(dtastk(isp1),larg%head))then
                go to 10
              endif
              if(m /= 0)then
                if(ordless)then
                  iord0=iordless
                  iordless=0
                  mstk0=mstk
                  if(.not. itfseqmatstk(im,isp0,larg%dbody(1),
     $                 m,is,ktftype(larg%dbody(1)%k) == ktfoper,
     $                 kap) >= 0)then
                    iordless=iord0
                    mstk=mstk0
                    isp=isp0
                    go to 10
                  endif
                  iordless=iord0
                  mstk=mstk0
                else
                  if(ktfreallistq(larg))then
                    do i=1,m
                      if(ktastk(isp1+i) /= larg%dbody(i)%k)then
                        go to 10
                      endif
                    enddo
                  else
                    do i=1,m
                      if(.not. tfsameq(ktastk(isp1+i),
     $                     larg%body(i)))then
                        go to 10
                      endif
                    enddo
                  endif
                endif
              endif
              if(kadv /= kadv0)then
                if((nh < 127 .and. ns > 3)
     $               .or. (nh < 1023 .and. ns > 7)
     $               .or. (ns > 15 .and. nh .le. 32767))then
                  kad=ktfrehash(kad,iup)
                else
                  klist(dtbl%prev)=dtbl%next
                  if(dtbl%next /= 0)then
                    klist(dtbl%next+1)=dtbl%prev
                  endif
                  dtbl%next=kadv0
                  klist(kadv0+1)=kadv
                  dtbl%prev=khash
                  klist(khash)=kadv
                endif
              endif
              ev=.true.
              if(dtbl%compile < rlist(levelcompile))then
                call tfsetarg(dtbl,irtc)
                if(irtc /= 0)then
                  return
                endif
              endif
              irtc=0
              if(def)then
                kx%k=ktfref+ksad_loc(dtbl%bodyc%k)
                return
              endif
              kx=dtbl%bodyc
              if(dtbl%pat /= -1)then
                if(ktflistq(kx,klx))then
                  kx=tfleval(klx,.true.,irtc)
                elseif(ktfsymbolq(kx))then
                  kx=tfsyeval(kx,irtc)
                elseif(ktfpatq(kx))then
                  kx=tfpateval(kx,irtc)
                endif
                if(irtc /= 0)then
                  call tfcatchreturn(irtcret,kx,irtc)
                endif
              endif
              return
            endif
 10         kadv=dtbl%next
            ns=ns+1
          enddo
        endif
 20     if(def)then
          irtc=-1
          return
        endif
        kad=dhash%next
      endif
      if(kad /= 0 .and. .not. def)then
c$$$        call loc_deftbl(kad,dtbl)
c$$$        do i=1,dtbl%npat
c$$$          kpa=dtbl%pattbl(i)
c$$$          call descr_sad(kpa,pat)
c$$$          call tfdebugprint(kpa,'deval-dtbl',1)
c$$$          write(*,*)': ',i,ktfaddr(kpa),pat%ref
c$$$          if(associated(pat%equiv))then
c$$$            call tfdebugprint(sad_descr(pat%equiv),'ktdaloc-equiv',1)
c$$$          endif
c$$$        enddo
        call tfdeval1(isp1,kad,kx,ev,ordless,ns,irtc)
      else
        irtc=-1
      endif
      return
      end

        subroutine tfgetdefargp(kl,kas,kp,ev,irtc)
        use tfstk
        implicit none
        type (sad_dlist) , intent(inout)::kl
        type (sad_descriptor) kv
        integer*8 , intent(in)::kas
        integer*8 kp
        integer*4 isp0
        integer*4 , intent(out)::irtc
        logical*4 , intent(out)::ev
        isp=isp+1
        isp0=isp
        dtastk(isp)=kl%head
        call tfgetllstkall(kl)
        call tfdeval(isp0,dfromk(kas),kv,1,.true.,ev,irtc)
        kp=ktfaddrd(kv)
        isp=isp0-1
        return
        end subroutine
