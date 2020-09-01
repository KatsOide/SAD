      subroutine tfreplace(k,kr,kx,all,eval,rule,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k,kr
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ki,tfreplacestk
      type (sad_dlist), pointer :: lri,lr,klr
      integer*4 ,intent(out):: irtc
      integer*4 i,nrule,isp1,j,itfmessage
      logical*4 ,intent(in):: all,eval,rule
      logical*4 symbol,rep
      irtc=0
      isp1=isp
      symbol=.true.
      if(tflistq(kr,klr))then
        call tfflattenstk(klr,-1,ktfoper+mtflist,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        nrule=isp-isp1
        do i=isp,isp1+1,-1
          ki=dtastk(i)
          if(tfruleq(ki%k,lri))then
            j=i+i-isp1
            dtastk(j-1:j)=lri%dbody(1:2)
            ivstk2(2,j-1)=merge(1,0,.not. tfconstpatternq(ktastk(j-1)))
            symbol=symbol .and.
     $           iand(ktfmask,ktastk(j-1)) .eq. ktfsymbol
          elseif(ki%k .eq. ktfoper+mtfnull)then
            ktastk(i+i-isp1-1:isp+nrule-2)=ktastk(i+i-isp1+1:isp+nrule)
            nrule=nrule-1
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Rule or List of rules"')
            go to 9000
          endif
        enddo
        isp=isp+nrule
      elseif(tfruleq(kr%k,lr))then
        call tfgetllstkall(lr)
        ivstk2(2,isp-1)=merge(1,0,.not. tfconstpatternq(ktastk(isp-1)))
        symbol=iand(ktfmask,ktastk(isp-1)) .eq. ktfsymbol
        nrule=1
      else
        irtc=itfmessage(9,'General::wrongtype',
     $           '"Rule or List of rules"')
        go to 9000
      endif
      if(rule)then
        return
      endif
      if(nrule .le. 0)then
        kx=k
      else
        if(symbol)then
          call tfreplacesymbolstk(k,isp1,nrule,kx,.false.,rep,irtc)
        else
          call tfinitrule(isp1,nrule)
          kx=tfreplacestk(k,isp1,nrule,all,rep,irtc)
          call tfresetrule(isp1,nrule)
        endif
      endif
      if(eval .and. rep .and. irtc .eq. 0)then
        kx=tfeevalref(kx,irtc)
      endif
 9000 isp=isp1
      return
      end

      subroutine tfinitrule(ispr,nrule)
      use tfstk
      implicit none
      type (sad_descriptor) kp
      integer*4 ,intent(in):: ispr,nrule
      integer*4 i,isp0
      do i=ispr+1,ispr+nrule*2,2
        kp=dtastk(i)
        if(ktfnonrealq(kp) .and. ivstk2(2,i) .eq. 1)then
          isp0=isp
          call tfinitpat(isp0,kp)
          ivstk2(1,i)=isp0
          ivstk2(2,i)=isp
        endif
      enddo
      return
      end

      subroutine tfresetrule(ispr,nrule)
      use tfstk
      implicit none
      integer*4 ,intent(in):: ispr,nrule
      integer*4 i
      do i=ispr+1,ispr+nrule*2,2
        if(ktfpatq(dtastk(i)) .or. ktflistq(dtastk(i)))then
          call tfresetpat(dtastk(i))
        endif
      enddo
      return
      end

      recursive function tfreplacestk(k,ispr,nrule,all,rep,irtc)
     $     result(kx)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kp,ki,k1,kir,ks,kd,tfcompose
      type (sad_dlist), pointer :: klir,kl
      type (sad_rlist), pointer :: klr
      type (sad_pat), pointer :: pat
      integer*8 kair
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: irtc
      integer*4 i,m,isp1,isp0,isp2,
     $     itfmessageexp,mstk0,itfpmat,iop
      logical*4 ,intent(out):: rep
      logical*4 ,intent(in):: all
      logical*4 noreal,rep1
      irtc=0
      noreal=.true.
      mstk0=mstk
      isp0=isp
      LOOP_I: do i=ispr+1,ispr+nrule*2,2
        kp=dtastk(i)
        noreal=noreal .and. ktfnonrealq(kp%k) .and.
     $       ivstk2(2,i) .eq. 0
        if(ktfnonrealq(kp) .and. ivstk2(2,i) .ne. 0)then
          iop=iordless
          iordless=0
          m=itfpmat(k,kp)
          iordless=iop
          if(m .ge. 0)then
            kx=dtastk(i+1)
            isp2=isp
            call tfpvrulestk(ivstk2(1,i),ivstk2(2,i))
            if(isp .gt. isp2)then
              call tfreplacesymbolstk(kx,isp2,(isp-isp2)/2,kx,
     $             .false.,rep1,irtc)
            endif
            call tfresetpat(kp)
            isp=isp0
            mstk=mstk0
            rep=.true.
            return
          endif
          isp=isp0
          mstk=mstk0
        else
          if(tfsameq(k,kp))then
            kx=dtastk(i+1)
            rep=.true.
            return
          endif
        endif
      enddo LOOP_I
      rep=.false.
      if(.not. all)then
        kx=k
        return
      endif
      if(ktflistq(k,kl))then
        if(noreal .and. ktfreallistq(kl))then
          ki=kl%head
          k1=tfreplacestk(ki,ispr,nrule,.true.,rep,irtc)
          if(irtc .ne. 0 .or. .not. rep)then
            kx=k
            return
          endif
          m=kl%nl
          kx=kxavaloc(-1,m,klr)
          klr%head=dtfcopy(k1)
          klr%rbody(1:m)=kl%rbody(1:m)
        else
          isp1=isp
          isp=isp+1
          dtastk(isp)=tfreplacestk(kl%head,ispr,nrule,.true.,rep,irtc)
          if(irtc .ne. 0)then
            kx=k
            isp=isp1
            return
          endif
          do i=1,kl%nl
            isp=isp+1
            kir=tfreplacestk(kl%dbody(i),ispr,nrule,.true.,rep1,irtc)
            if(irtc .ne. 0)then
              kx=k
              isp=isp1
              return
            endif
            rep=rep .or. rep1
            dtastk(isp)=kir
            if(ktflistq(kir,klir))then
              kair=ktfaddrd(kir)
              if(klir%head%k .eq. ktfoper+mtfnull)then
                if(ktastk(isp1+1) .ne. ktfoper+mtffun)then
                  rep=.true.
                  isp=isp-1
                  call tfgetllstkall(klir)
                endif
              endif
            endif
          enddo
          kx=merge(tfcompose(isp1+1,ktastk(isp1+1),irtc),k,rep)
          isp=isp1
        endif
      elseif(ktfpatq(k,pat))then
        rep=.false.
        if(pat%sym%loc .ne. 0)then
          ks=tfreplacestk(pat%sym%alloc,ispr,nrule,
     $         .false.,rep,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(rep)then
            if(ktfnonsymbolq(ks))then
              kx=k
              irtc=itfmessageexp(999,'General::reppat',k)
              return
            endif
          endif
        else
          ks%k=0
        endif
        if(ktftype(pat%expr%k) .ne. ktfref)then
          k1=tfreplacestk(pat%expr,ispr,nrule,.true.,rep1,irtc)
          if(irtc .ne. 0)then
            return
          endif
          rep=rep .or. rep1
        else
          k1=pat%expr
        endif
        kd=pat%default
        if(ktftype(kd%k) .ne. ktfref)then
          kd=tfreplacestk(kd,ispr,nrule,.true.,rep1,irtc)
          if(irtc .ne. 0)then
            return
          endif
          rep=rep .or. rep1
        endif
        kx=merge(kxpcopyss(k1,pat%head,ks,kd),k,rep)
      else
        kx=k
      endif
      return
      end

      subroutine tfreplacesymbolstk(k,ispr,nrule,kx,scope,rep,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tfreplacesymbolstk1
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: irtc
      integer*4 nrule1
      logical*4 ,intent(in):: scope
      logical*4 ,intent(out):: rep
      call tfsortsymbolstk(ispr,nrule,nrule1)
      kx=tfreplacesymbolstk1(k,ispr,nrule1,scope,rep,irtc)
      return
      end

      recursive function tfreplacesymbolstk1(k,ispr,nrule,
     $     scope,rep,irtc) result(kx)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kd,k1,ki,ks,tfcompose
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: kl,kli
      type (sad_rlist), pointer :: klr
      integer*8 ka1,kas
      integer*4 ,intent(out):: irtc
      integer*4 i,m,isp1,j,itfmessageexp, id
      integer*4 ,intent(in):: ispr,nrule
      logical*4 ,intent(out):: rep
      logical*4 ,intent(in):: scope
      logical*4 rep1,tfmatchsymstk,tfsymbollistqo
      irtc=0
      rep=.false.
      kx=k
      if(ktflistq(k,kl))then
        if(.not. tfsymbollistqo(kl))then
c     call tfdebugprint(ktflist+ktfaddr(k),'repsymstk',3)
          return
        endif
        m=kl%nl
        k1=tfreplacesymbolstk1(kl%head,ispr,nrule,scope,rep,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfreallistq(kl))then
          if(rep)then
            kx=kxavaloc(-1,m,klr)
            klr%rbody(1:m)=kl%rbody(1:m)
            klr%head=dtfcopy(k1)
          endif
        else
          isp1=isp
          isp=isp+1
          dtastk(isp)=k1
          if(ktfoperq(k1))then
            ka1=ktfaddrd(k1)
            if(scope .and. ka1 .eq. mtffun)then
              if(m .eq. 2)then
                call tfreplacefunstk(kl,ispr,nrule,kx,rep1,irtc)
                rep=rep .or. rep1
                isp=isp1
                return
              endif
            else
              id=ilist(2,klist(ifunbase+ka1))
              if(id .eq. nfunif)then
                call tfreplaceifstk(kl,
     $               ispr,nrule,kx,scope,rep1,irtc)
                rep=rep .or. rep1
                isp=isp1
                return
              elseif(scope .and. id .eq. nfunwith)then
                if(kl%nl .eq. 2)then
                  call tfreplacewithstk(kl,
     $                 ispr,nrule,kx,rep1,irtc)
                  rep=rep .or. rep1
                  isp=isp1
                  return
                endif
              endif
            endif
          endif
          do i=1,m
c            call tfdebugprint(kl%dbody(i),'repsymstk',1)
            ki=tfreplacesymbolstk1(kl%dbody(i),ispr,nrule,
     $           scope,rep1,irtc)
            if(irtc .ne. 0)then
              return
            endif
c            call tfdebugprint(ki,'==> ',1)
            rep=rep .or. rep1
            isp=isp+1
            dtastk(isp)=ki
            if(ktfsequenceq(ki%k,kli))then
              isp=isp-1
              call tfgetllstkall(kli)
            endif
          enddo
          if(rep)then
            kx=tfcompose(isp1+1,ktastk(isp1+1),irtc)
          endif
          isp=isp1
        endif
        return
      elseif(ktfpatq(k,pat))then
        if(pat%sym%loc .ne. 0)then
          ks=pat%sym%alloc
          kas=ktfaddr(ks)
          if(tfmatchsymstk(kas,ispr,nrule,j))then
            if(ktfnonsymbolq(ktastk(ivstk2(1,j)+1)))then
              irtc=itfmessageexp(999,'General::reppat',k)
              return
            endif
            ks=dtastk(ivstk2(1,j)+1)
            rep=.true.
          endif
        endif
        if(ktftype(pat%expr%k).ne. ktfref)then
          k1=tfreplacesymbolstk1(pat%expr,ispr,nrule,
     $         scope,rep1,irtc)
          if(irtc .ne. 0)then
            return
          endif
          rep=rep .or. rep1
        else
          k1=pat%expr
        endif
        kd=pat%default
        if(ktftype(kd%k) .ne. ktfref)then
          kd=tfreplacesymbolstk1(kd,ispr,nrule,scope,rep1,irtc)
          if(irtc .ne. 0)then
            return
          endif
          rep=rep .or. rep1
        endif
        if(rep)then
          kx=kxpcopyss(k1,pat%head,ks,kd)
        endif
        return
      elseif(ktfsymbolq(k))then
c        call tfdebugprint(k,'repsymstk-symbol',1)
        if(tfmatchsymstk(ktfaddr(k),ispr,nrule,j))then
          kx=dtastk(ivstk2(1,j)+1)
c          if(ktfsymbolq(kx))then
c            call tfdebugprint(kx,'==> ',1)
c            write(*,*)'with ',ilist(2,ktfaddr(kx)-3)
c          endif
          rep=.true.
          return
        endif
      endif
      return
      end

      logical*4 function tfmatchsymstk(ka,ispr,nrule,j)
      use tfstk
      implicit none
      type (sad_symbol), pointer :: sym
      integer*8 ,intent(in):: ka
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: j
      logical*4 tfmatchsymstk1
      call loc_sym(ka,sym)
      tfmatchsymstk=tfmatchsymstk1(sym%loc,max(0,sym%gen),ispr,nrule,j)
      return
      end

      logical*4 function tfmatchsymstk1(loc,iag,ispr,nrule,j)
      use tfstk
      implicit none
      integer*8 ,intent(in):: loc
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: j
      integer*4 jm,jl,jh,iag
      jl=1
      jh=nrule
      do while (jh .ge. jl)
        jm=(jl+jh)/2
        j=ispr+jm*2-1
        if(loc .eq. ktastk2(j+1))then
          if(iag .eq. ivstk2(2,j))then
            tfmatchsymstk1=.true.
            return
          elseif(iag .lt. ivstk2(2,j))then
            jh=jm-1
          else
            jl=jm+1
          endif
        elseif(loc .lt. ktastk2(j+1))then
          jh=jm-1
        else
          jl=jm+1
        endif
      enddo
      tfmatchsymstk1=.false.
      return
      end

      subroutine tfsortsymbolstk(ispr,n,n1)
      use tfstk
      implicit none
      integer*4 ,intent(in):: ispr,n
      integer*4 ,intent(out):: n1
      integer*4 i,isp1,j,ig0,ig1
      integer*8 kai,kz0
      integer*8, parameter:: k32=2**32
      integer*8, allocatable :: kz(:),kg(:)
      integer*4, allocatable :: itab(:)
      if(n .eq. 1)then
        kai=ktfaddr(ktastk(ispr+1))
        ivstk2(1,ispr+1)=ispr+1
        ivstk2(2,ispr+1)=
     $       max(0,ilist(2,kai-1))
        ktastk2(ispr+2)=klist(kai)
        n1=1
      else
        allocate (kz(n),itab(n),kg(n))
        do i=1,n
          itab(i)=i
          kai=ktfaddr(ktastk(ispr+i*2-1))
          kz(i)=klist(kai)
          kg(i)=max(0,ilist(2,kai-1))*k32+i
        enddo
        call tfsorti(itab,kz,kg,n)
        isp1=ispr-1
        kz0=0
        ig0=0
        do i=1,n
          j=itab(i)
          ig1=int(kg(j)/k32)
          if(kz(j) .ne. kz0 .or. ig1 .ne. ig0)then
            kz0=kz(j)
            ig0=ig1
            isp1=isp1+2
            ivstk2(1,isp1)=ispr+j*2-1
            ivstk2(2,isp1)=ig0
            ktastk2(isp1+1)=kz0
          endif
        enddo
        n1=(isp1-ispr+1)/2
        deallocate(kz,itab,kg)
      endif
      return
      end

      recursive subroutine tfsorti(itab,iz,kg,n)
      implicit none
      integer*8 ,intent(in):: iz(n),kg(n)
      integer*4 ,intent(in):: n
      integer*4 ,intent(out):: itab(n)
      integer*4 m,i1,i2,is,im,ip1,ip2
      if(n .le. 1)then
        return
      endif
      i1=itab(1)
      i2=itab(n)
      if(iz(i1) .gt. iz(i2) .or.
     $     iz(i1) .eq. iz(i2) .and. kg(i1) .gt. kg(i2))then
        is=i1
        i1=i2
        i2=is
      endif
      if(n .eq. 2)then
        itab(1)=i1
        itab(2)=i2
        return
      endif
      m=(n+1)/2
      im=itab(m)
      if(iz(i1) .lt. iz(im) .or.
     $     iz(i1) .eq. iz(im) .and. kg(i1) .lt. kg(im))then
        if(iz(im) .gt. iz(i2) .or.
     $     iz(im) .eq. iz(i2) .and. kg(im) .gt. kg(i2))then
          is=im
          im=i2
          i2=is
        endif
      else
        is=im
        im=i1
        i1=is
      endif
      itab(1)=i1
      itab(m)=im
      itab(n)=i2
      if(n .eq. 3)then
        return
      endif
      itab(m)=itab(2)
      itab(2)=im
      ip1=3
      ip2=n-1
      do while(ip1 .le. ip2)
        do while((iz(itab(ip1)) .lt. iz(im) .or.
     $       iz(itab(ip1)) .eq. iz(im) .and. kg(itab(ip1)) .lt. kg(im))
     $       .and. ip1 .le. ip2)
          ip1=ip1+1
        enddo
        do while((iz(im) .lt. iz(itab(ip2)) .or.
     $       iz(itab(ip2)) .eq. iz(im) .and. kg(im) .lt. kg(itab(ip2)))
     $       .and. ip1 .le. ip2)
          ip2=ip2-1
        enddo
        if(ip2 .gt. ip1)then
          is=itab(ip1)
          itab(ip1)=itab(ip2)
          itab(ip2)=is
          ip1=ip1+1
          ip2=ip2-1
        endif
      enddo
      ip1=ip1-1
      is=itab(ip1)
      itab(ip1)=im
      itab(2)=is
      call tfsorti(itab,iz,kg,ip1-1)
      call tfsorti(itab(ip1+1),iz,kg,n-ip1)
      return
      end

      subroutine tfpvrulestk(isp1,isp2)
      use tfstk
      use tfcode
      use funs
      implicit none
      type (sad_descriptor) kx
      type (sad_pat), pointer :: pat
      integer*8 kp,kap
      integer*4 ,intent(in):: isp1,isp2
      integer*4 i,ispb,ispe
      logical*4 rep
      do i=isp1+1,isp2-1,2
        kp=ktastk(i)
        kap=ktfaddr(kp)
        call loc_pat(kap,pat)
        do while(associated(pat%equiv))
          pat=>pat%equiv
        enddo
        ispb=isp
        call tfgetstkstk(pat%value,rep)
        ispe=isp
        kx=tfsequence(ispb,ispe)
        isp=ispb+1
        dtastk(isp)=pat%sym%alloc
        isp=isp+1
        dtastk(isp)=kx
      enddo
      return
      end

      subroutine tfreplaceifstk(list,ispr,nrule,kx,
     $     scope,rep,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ki,kr
      type (sad_descriptor) tfreplacesymbolstk1,tfcompose
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: klx
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: irtc
      integer*4 i,isp1,j
      logical*4 ,intent(out):: rep
      logical*4 ,intent(in):: scope
      logical*4 rep1,rep2
      irtc=0
      isp1=isp
      if(list%nl .eq. 0)then
        if(rep)then
          kx=kxaaloc(-1,0,klx)
          klx%head%k=ktfoper+ktfaddr(ktastk(isp1))
        else
          kx=sad_descr(list)
        endif
        return
      endif
      ki=list%dbody(1)
      if(ktfrealq(ki) .or. ktfstringq(ki) .or. ktfoperq(ki))then
        isp=isp+1
        dtastk(isp)=ki
      else
        kr=tfreplacesymbolstk1(ki,ispr,nrule,scope,rep1,irtc)
        if(irtc .ne. 0)then
          isp=isp1
          return
        endif
        call tfgetstkstk(kr,rep2) 
        rep=rep .or. rep2 .or. rep1
      endif
      if(isp .eq. isp1+1)then
        i=merge(merge(2,3,ktastk(isp1+1) .ne. 0),0,
     $       ktfrealq(ktastk(isp1+1)))
      else
        rep=.true.
        i=0
      endif
      if(i .ne. 0)then
        rep=.true.
        j=isp+1
        if(i .le. list%nl)then
          ki=list%dbody(i)
          if(ktfrealq(ki) .or. ktfstringq(ki) .or. ktfoperq(ki))then
            kx=ki
          else
            kr=tfreplacesymbolstk1(ki,ispr,nrule,scope,rep1,irtc)
            if(irtc .ne. 0)then
              isp=isp1
              return
            endif
            call tfgetstkstk(kr,rep2) 
            rep=rep .or. rep1 .or. rep2
            kx=merge(dtastk(j),dxnullo,j .le. isp)
          endif
        else
          kx%k=ktfoper+mtfnull
        endif
        return
      endif
      do i=2,list%nl
        ki=list%dbody(i)
        if(ktfrealq(ki) .or. ktfstringq(ki) .or. ktfoperq(ki))then
          isp=isp+1
          dtastk(isp)=ki
        else
          kr=tfreplacesymbolstk1(ki,ispr,nrule,scope,rep1,irtc)
          if(irtc .ne. 0)then
            isp=isp1
            return
          endif
          call tfgetstkstk(kr,rep2)
          rep=rep2 .or. rep .or. rep1
        endif
      enddo
      kx=merge(tfcompose(isp1,ktastk(isp1),irtc),sad_descr(list),rep)
      return
      end

      subroutine tfreplacewithstk(list,ispr,nrule,kx,rep,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) k1,k2,ki,ki2
      type (sad_descriptor) tfreplacesymbolstk1
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist), pointer :: kl1,kli,klx,klx1
      integer*8 kai,ki1,ka1,kai1,ksave
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: irtc
      integer*4 i,ispj,ispa,ispb
      logical*4 ,intent(out):: rep
      logical*4 rep1,tfmatchsymstk
      irtc=0
      ispa=isp
      k1=list%dbody(1)
      rep=.false.
      if(tflistq(k1,kl1))then
        ktastk(ispa+1:ispa+nrule*2)=ktastk(ispr+1:ispr+nrule*2)
        ktastk2(ispa+1:ispa+nrule*2)=ktastk2(ispr+1:ispr+nrule*2)
c        do i=1,nrule*2
c          ktastk(ispa+i)=ktastk(ispr+i)
c          ktastk2(ispa+i)=ktastk(ivstkoffset+ispr+i)
c        enddo
        isp=isp+nrule*2
        ispb=isp
        ka1=ktfaddrd(k1)
        do i=1,kl1%nl
          ki=kl1%dbody(i)
          if(ktflistq(ki,kli))then
            kai=ktfaddr(ki)
            if((kli%head%k .eq. ktfoper+mtfset .or.
     $           kli%head%k .eq. ktfoper+mtfsetdelayed) .and.
     $           kli%nl .eq. 2)then
              ki1=kli%dbody(1)%k
              if(ktfsymbolq(ki1))then
                kai1=ktfaddr(ki1)
                if(tfmatchsymstk(kai1,ispa,nrule,ispj))then
                  ivstk2(2,ispj)=min(-ivstk2(2,ispj)-1,ivstk2(2,ispj))
                  ki2=tfreplacesymbolstk1(kli%dbody(2),ispr,nrule,
     $                 .true.,rep1,irtc)
                  if(irtc .ne. 0)then
                    isp=ispa
                    return
                  endif
                  isp=isp+1
                  if(rep1)then
                    dtastk(isp)=kxadaloc(-1,2,klx1)
                    klx1%head=dlist(kai)
                    klx1%dbody(1)%k=ktfcopy1(ki1)
                    klx1%dbody(2)=dtfcopy(ki2)
                    rep=.true.
                  else
                    ktastk(isp)=ki1
                  endif
                  cycle
                endif
              endif
            endif
          endif
          isp=isp+1
          dtastk(isp)=tfreplacesymbolstk1(ki,ispr,nrule,
     $         .true.,rep1,irtc)
          if(irtc .ne. 0)then
            isp=ispa
            return
          endif
          rep=rep .or. rep1
        enddo
        if(rep)then
          k1=kxmakelist(ispb)
        endif
        k2=tfreplacesymbolstk1(list%dbody(2),ispa,nrule,
     $       .true.,rep1,irtc)
        if(irtc .ne. 0)then
          isp=ispa
          return
        endif
c        ilist(2,ktfaddr(k2)-3)=ior(ilist(2,ktfaddr(k2)-3),kmodsymbol)
        rep=rep .or. rep1
        if(rep)then
          kx=kxadaloc(-1,2,klx)
          klx%head=list%head
          klx%dbody(1)=dtfcopy1(k1)
          klx%dbody(2)=dtfcopy(k2)
        else
          kx=sad_descr(list)
        endif
      else
        ksave=list%head%k
        list%head%k=ktfoper+mtfhold
        kx=tfreplacesymbolstk1(sad_descr(list),
     $       ispa,nrule,.true.,rep,irtc)
        if(irtc .eq. 0 .and. ktflistq(kx,klx))then
          klx%head%k=ksave
        endif
        list%head%k=ksave
      endif
      isp=ispa
      return
      end

      subroutine tfreplacefunstk(list,ispr,nrule,kx,rep,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,ki
      type (sad_descriptor) tfreplacesymbolstk1
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist), pointer :: kl1,klx
      integer*8 kai,ka1,ksave
      integer*4 ,intent(in):: ispr,nrule
      integer*4 ,intent(out):: irtc
      integer*4 i,isp1,j,ispi
      logical*4 ,intent(out):: rep
      logical*4 rep1,rej,tfmatchsymstk1
      irtc=0
      isp1=isp
      k1=list%dbody(1)
      if(tflistq(k1,kl1))then
        do i=1,kl1%nl
          ki=kl1%dbody(i)
          if(ktfsymbolq(ki))then
            kai=ktfaddr(ki)
            isp=isp+1
            ktastk(isp)=klist(kai)
            ivstk2(2,isp)=max(0,ilist(2,kai-1))
          endif
        enddo
      elseif(ktfsymbolq(k1))then
        isp=isp+1
        ka1=ktfaddrd(k1)
        ktastk(isp)=klist(ka1)
        ivstk2(2,isp)=max(0,ilist(2,ka1-1))
      endif
      rej=.false.
      if(isp .gt. isp1)then
        r1: do j=isp1+1,isp
          if(tfmatchsymstk1(ktastk(j),ivstk(2,j),ispr,nrule,ispi))then
            ivstk2(2,ispi)=Min(-ivstk2(2,ispi)-1,ivstk2(2,ispi))
            rej=.true.
            cycle r1
          endif
        enddo r1
      endif
      ksave=list%head%k
      list%head%k=ktfoper+mtfhold
      kx=tfreplacesymbolstk1(sad_descr(list),ispr,nrule,
     $     .true.,rep1,irtc)
      if(irtc .eq. 0 .and. ktflistq(kx,klx))then
        klx%head%k=ksave
      endif
      list%head%k=ksave
      rep=rep .or. rep1
      if(rej)then
        do i=1,nrule
          ispi=ispr+i*2-1
          ivstk2(2,ispi)=max(-(ivstk2(2,ispi)+1),ivstk2(2,ispi))
        enddo
      endif
      isp=isp1
      return
      end

      subroutine tfreplacerepeated(k,kr,kx,all,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k,kr
      type (sad_descriptor) ,intent(out)::kx
      type (sad_descriptor) k1
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: all,eval
      kx=k
      k1%k=ktfref
      irtc=0
      do
        if(tfsameq(k1,kx))then
          return
        endif
        k1=kx
        call tfreplace(k1,kr,kx,all,eval,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      enddo
      end

      subroutine tfgetoption(symbol,kr,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kr,kx,tfgetoption1
      type (sad_dlist), pointer :: lr
      integer*4 irtc
      character*(*) symbol
      logical*4 rep
      if(tfruleq(kr%k,lr))then
        kx=tfgetoption1(ktfsymbolz(symbol,len(symbol)),lr,rep)
        irtc=0
        if(.not. rep)then
          kx%k=ktfref
        endif
      else
        irtc=-1
      endif
      return
      end

      recursive function tfgetoption1(ka,list,rep) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist) list
      type (sad_dlist), pointer :: listi
      integer*8 ka
      integer*4 i
      logical*4 rep
      rep=.false.
      if(list%head%k .eq. ktfoper+mtflist)then
        do i=1,list%nl
          call descr_sad(list%dbody(i),listi)
          kx=tfgetoption1(ka,listi,rep)
          if(rep)then
            return
          endif
        enddo
      elseif(ktfnonreallistqo(list))then
        if(ktfsymbolq(list%dbody(1)))then
          rep=tfsamesymbolq(ka,list%dbody(1)%k)
          if(rep)then
            kx=list%dbody(2)
          endif
        endif
      endif
      return
      end

      subroutine tfgetoptionstk(isp1,kaopt,optname,nopt,ispopt,irtc)
      use tfstk
      implicit none
      type (sad_dlist), pointer :: lri
      type (sad_descriptor) kaopt(nopt),tfgetoption1
      integer*4 isp1,nopt,ispopt,i,j,isp0,irtc,lenw
      logical*4 rep
      character*(*) optname(nopt)
      isp0=isp
      if(kaopt(1)%k .eq. 0)then
        do i=1,nopt
          kaopt(i)=kxsymbolz(optname(i),lenw(optname(i)))
        enddo
      endif
      do i=isp0,isp1,-1
        if(.not. tfruleq(ktastk(i)))then
          ispopt=i+1
          go to 1
        endif
      enddo
      ispopt=isp1
 1    irtc=0
      isp=isp0+nopt
      if(ispopt .gt. isp)then
        ktastk(isp0+1:isp0+nopt)=ktfref
        return
      endif
      LOOP_J: do j=1,nopt
        do i=ispopt,isp0
          call loc_sad(ktfaddr(ktastk(i)),lri)
          dtastk(isp0+j)=tfgetoption1(kaopt(j)%k,lri,rep)
          if(rep)then
            cycle LOOP_J
          endif
        enddo
        ktastk(isp0+j)=ktfref
      enddo LOOP_J
      return
      end
