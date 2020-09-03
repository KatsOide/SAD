      function tfmodule(isp1,module,eval,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx,ke,kxl1,tfredefslist
      type (sad_dlist), pointer :: lvlist,kle
      type (sad_symdef), pointer :: symdi
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: isp1
      integer*4 i,isp0,isp2,itfmessage
      logical*4 ,intent(in):: module,eval
      logical*4 rep
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        kx=dxnullo
        return
      endif
      if(ktfnonlistq(dtastk(isp1+1),lvlist))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        kx=dxnullo
        return
      endif
      if(module .and. eval)then
        lgeneration=lgeneration+1
      endif
      isp0=isp
      call tfmlocalv(lvlist,module,eval,irtc)
      if(irtc .ne. 0)then
        if(eval)then
          go to 9200
        else
          kx=dxnullo
          isp=isp0
          return
        endif
      endif
      isp2=isp
      rep=.false.
      ke=dtastk(isp0)
      if(isp2 .eq. isp0)then
        if(.not. eval)then
          kx=dtastk(isp0)
          return
        endif
      elseif(module)then
        call tfreplacesymbolstk(dtastk(isp0),isp0,(isp2-isp0)/2,ke,
     $     .false.,rep,irtc)
        if(irtc .ne. 0)then
          isp=isp2
          go to 9200
        endif
      endif
      if(eval)then
        if(ktflistq(ke,kle))then
c          call tfdebugprint(ke,'tfmodule',1)
          kx=tfleval(kle,.true.,irtc)
c          call tfdebugprint(kx,'tfmodule-1',1)
        else
          kx=tfeevalref(ke,irtc)
        endif
        isp=isp2
      elseif(rep)then
        if(isp2 .gt. isp0)then
c          call tfdebugprint(ktflist+ksad_loc(lvlist%head),'module-1',2)
          kxl1=tfredefslist(isp0,isp2,lvlist)
c          call tfdebugprint(ktflist+kal1,'==>',2)
        else
          kxl1=sad_descr(lvlist)
        endif
        isp=isp0+1
        ktastk(isp)=ktastk(isp1)
        isp=isp+1
        dtastk(isp)=kxl1
        isp=isp+1
        dtastk(isp)=ke
        kx=kxcompose(isp0+1)
        isp=isp0
        return
      else
        kx=ke
        isp=isp0
        return
      endif
 9200 do i=isp,isp0+2,-2
c        call tfdebugprint(dtastk(i),'tfmodule-delete',1)
        call descr_sad(dtastk(i),symdi)
        call tfdelete(symdi,.true.,.true.)
      enddo
      isp=isp0
      return
      end

      subroutine tfmlocalv(list,module,eval,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ki,ki1i,ki1,ki2
      type (sad_dlist) list
      type (sad_dlist), pointer :: kli1,kli2,kli
      type (sad_symbol), pointer :: symi1i,symi,symi1
      type (sad_symdef), pointer :: symdi1,symdv
      integer*4 m,mi,lgi,ii,i,irtc,itfmessage,lg,isps
      logical*4 module,eval
      if(list%head%k .ne. ktfoper+mtflist)then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      m=list%nl
      if(m .eq. 0)then
        irtc=0
        return
      endif
      if(ktfreallistq(list))then
        go to 100
      endif
      if(eval)then
        do i=1,m
          ki=list%dbody(i)
          if(ktfsymbolq(ki,symi))then
            lgi=symi%gen
            if(module)then
              lg=lgeneration
            else
              lg=max(0,lgi)
            endif
            isp=isp+2
            dtastk(isp-1)=ki
            dtastk(isp  )=kxnaloc1(lg,symi%loc)
          elseif(ktflistq(ki,kli))then
            if(kli%head%k .ne. ktfoper+mtfset .and.
     $           kli%head%k .ne. ktfoper+mtfsetdelayed)then
              go to 100
            endif
            if(kli%nl .ne. 2)then
              go to 100
            endif
            ki1=kli%dbody(1)
            ki2=kli%dbody(2)
            if(kli%head%k .eq. ktfoper+mtfset)then
              if(ktflistq(ki2,kli2))then
                isps=isp
                ki2=tfleval(kli2,.true.,irtc)
                if(irtc .ne. 0)then
                  go to 200
                endif
              elseif(ktfsymbolq(ki2) .or. ktfpatq(ki2))then
                isps=isp
                ki2=tfeevalref(ki2,irtc)
                if(irtc .ne. 0)then
                  go to 200
                endif
              endif
            endif
            if(ktflistq(ki1,kli1))then
              if(kli1%head%k .eq. ktfoper+mtflist)then
                mi=kli1%nl
                if(tflistq(ki2,kli2))then
                  if(kli2%nl .ne. mi)then
                    irtc=itfmessage(9,'General::equalleng',
     $                   'lhs and rhs of Set')
                    return
                  endif
                  if(ktfreallistq(kli1))then
                    go to 100
                  endif
                  do ii=1,mi
                    ki1i=kli1%dbody(ii)
                    if(ktfsymbolq(ki1i,symi1i))then
                      lgi=symi1i%gen
                      if(module)then
                        lg=lgeneration
                      else
                        lg=max(0,lgi)
                      endif
                      isp=isp+2
                      dtastk(isp-1)=ki1i
                      dtastk(isp  )=kxnaloc1(lg,symi1i%loc)
                      call descr_sad(dtastk(isp),symdv)
                      isps=isp
                      symdv%value=dtfcopy(kli2%dbody(ii))
                      symdv%sym%ref=1
                    else
                      go to 100
                    endif
                  enddo
                else
                  go to 100
                endif
              endif
            elseif(ktfsymbolqdef(ki1%k,symdi1))then
              lgi=symdi1%sym%gen
              if(module)then
                lg=lgeneration
              else
                lg=max(0,lgi)
              endif
              isp=isp+2
              dtastk(isp-1)=ki1
              dtastk(isp  )=kxnaloc1(lg,symdi1%sym%loc)
              call descr_sad(dtastk(isp),symdv)
              isps=isp
              symdv%value=dtfcopy(ki2)
              symdv%sym%ref=1
            else
              go to 100
            endif
          else
            go to 100
          endif
        enddo
      else
        do i=1,m
          ki=list%dbody(i)
          if(ktfsymbolq(ki,symi))then
            isp=isp+2
            dtastk(isp-1)=ki
            dtastk(isp  )=dxsycopy(symi)
          elseif(ktflistq(ki,kli))then
            if(kli%head%k .ne. ktfoper+mtfset .and.
     $           kli%head%k .ne. ktfoper+mtfsetdelayed)then
              go to 100
            endif
            if(kli%nl .ne. 2)then
              go to 100
            endif
            ki1=kli%dbody(1)
            if(tflistq(ki1,kli1))then
              if(ktfreallistq(kli1))then
                go to 100
              endif
              mi=kli1%nl
              do ii=1,mi
                ki1i=kli1%dbody(ii)
                if(ktfsymbolq(ki1i,symi1i))then
                  isp=isp+2
                  dtastk(isp-1)=ki1i
                  dtastk(isp  )=dxsycopy(symi1i)
                else
                  go to 100
                endif
              enddo
            elseif(ktfsymbolq(ki1,symi1))then
              isp=isp+2
              dtastk(isp-1)=ki1
              dtastk(isp  )=dxsycopy(symi1)
            else
              go to 100
            endif
          else
            go to 100
          endif
        enddo
      endif
      irtc=0
      return
 100  irtc=itfmessage(9,'General::wrongtype',
     $'"Symbol, Symbol (:)= value, {Symbol1, ..} (:)= {Value1, ..}"')
      return
 200  isp=isps
      return
      end

      recursive function tfredefslist(isp1,isp2,lvlist)
     $     result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k11,kj,k1
      type (sad_dlist) lvlist
      type (sad_dlist), pointer :: kl1,klj
      integer*4 isp1,isp2,i,isp0,m,j
      m=lvlist%nl
      if(m .eq. 0)then
        kx=sad_descr(lvlist)
        return
      endif
      isp=isp+1
      isp0=isp
      do j=1,m
        kj=lvlist%dbody(j)
        if(ktfsymbolq(kj))then
          do i=isp1+1,isp2,2
            if(tfsamesymbolq(dtastk(i),kj))then
              isp=isp+1
              dtastk(isp)=dtastk(i+1)
              exit
            endif
          enddo
        elseif(ktflistq(kj,klj))then
          isp=isp+1
          dtastk(isp)=klj%head
          k1=klj%dbody(1)
          if(ktflistq(k1,kl1))then
            k11=tfredefslist(isp1,isp2,kl1)
            isp=isp+1
            dtastk(isp)=k11
          else
            do i=isp1+1,isp2,2
              if(tfsamesymbolq(dtastk(i),k1))then
                isp=isp+1
                dtastk(isp)=dtastk(i+1)
                exit
              endif
            enddo
          endif
          isp=isp+1
          dtastk(isp)=klj%dbody(2)
          dtastk(isp-2)=kxcomposev(isp-2)
          isp=isp-2
        endif
      enddo
      ktastk(isp0)=ktfoper+mtflist
      kx=kxcomposev(isp0)
      isp=isp0-1
      return
      end

      subroutine tfwith(isp1,kx,eval,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx,ki,k1,kxi
      type (sad_descriptor) tfreplacestk
      type (sad_descriptor) tfreplacesymbolstk1
      type (sad_dlist), pointer :: list,listi,klx
      type (sad_symbol), pointer :: symi
      integer*4 isp1,irtc,isp2,i,itfmessage,lg0,n
      logical*4 rep,symbol,eval
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. tflistq(ktastk(isp1+1),list))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      irtc=0
      isp2=isp
      lg0=lgeneration
      if(list%nl .gt. 0)then
        symbol=.true.
        do i=1,list%nl
          isp=isp+2
          ki=list%dbody(i)
          if(ktflistq(ki,listi))then
            if(listi%head%k .eq. ktfoper+mtfset .or.
     $           listi%head%k .eq. ktfoper+mtfsetdelayed)then
              k1=listi%dbody(1)
              dtastk(isp-1)=k1
              if(eval)then
                if(.not. tfconstpatternq(k1))then
                  ivstk2(2,isp-1)=1
                elseif(.not. ktfrealq(k1))then
                  ivstk2(2,isp-1)=0
                endif
                symbol=symbol .and. ktfsymbolq(k1)
                if(listi%head%k .eq. ktfoper+mtfsetdelayed)then
                  dtastk(isp)=listi%dbody(2)
                  cycle
                else
                  ki=listi%dbody(2)
                endif
              elseif(ktfsymbolq(k1,symi))then
                ivstk2(2,isp-1)=0
              else
                isp=isp-2
                cycle
              endif
              go to 30
            endif
          endif
          if(eval)then
            symbol=symbol .and. ktfsymbolq(ki)
            if(.not. tfconstpatternq(ki))then
              ivstk2(2,isp-1)=1
            elseif(.not. ktfrealq(ki))then
              ivstk2(2,isp-1)=0
            endif
            dtastk(isp-1)=ki
          else
            isp=isp-2
            cycle
          endif
 30       if(eval)then
            if(ktflistq(ki,listi))then
              kxi=tfleval(listi,.true.,irtc)
            else
              kxi=tfeevalref(ki,irtc)
            endif
            if(irtc .ne. 0)then
              isp=isp2
              return
            endif
            dtastk(isp)=kxi
          else
            lgeneration=min(lg0+1,lgeneration+1)
            dtastk(isp)=dxsycopy(symi)
            symi%gen=lgeneration
          endif
        enddo
        if(eval)then
          if(symbol)then
            call tfreplacesymbolstk(dtastk(isp2),isp2,list%nl,kx,
     $           .true.,rep,irtc)
          else
            call tfinitrule(isp2,list%nl)
            kx=tfreplacestk(dtastk(isp2),
     $           isp2,list%nl,.true.,rep,irtc)
            call tfresetrule(isp2,list%nl)
          endif
          if(irtc .ne. 0)then
            isp=isp2
            return
          endif
        else
          if(isp .gt. isp2)then
            n=(isp-isp2)/2;
            isp=isp+1
            call tfreplacesymbolstk(dtastk(isp1+1),isp2,n,dtastk(isp),
     $           .true.,rep,irtc)
            isp=isp+1
            dtastk(isp)=tfreplacesymbolstk1(dtastk(isp1+2),isp2,n,
     $           .true.,rep,irtc)
            kx=kxmakelist(isp-2,klx);
            klx%head=dtastk(isp1)
          else
            irtc=-1
          endif
        endif
      else
        if(eval)then
          kx=dtastk(isp2)
        else
          irtc=-1
        endif
      endif
      if(eval)then
        if(ktflistq(kx,klx))then
          kx=tfleval(klx,.true.,irtc)
        elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
          kx=tfeevalref(kx,irtc)
        endif
      endif
      isp=isp2
      return
      end

      recursive function tfredefsymbol(isp1,isp2,k,rep)
     $     result(kx)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx,k,k1,k2,kd,ks
      type (sad_dlist), pointer :: list
      type (sad_rlist), pointer :: klr
      type (sad_symbol), pointer :: sym
      type (sad_pat), pointer :: pat
      integer*8 kas
      integer*4 isp1,isp2,i,isp0
      logical*4 rep,rep1
      rep=.false.
      kx%k=ktfoper+mtfnull
      if(ktflistq(k,list))then
        isp=isp+1
        isp0=isp
        dtastk(isp)=tfredefsymbol(isp1,isp2,list%head,rep)
        if(ktfnonreallistqo(list))then
          do i=1,list%nl
            isp=isp+1
            dtastk(isp)=tfredefsymbol(isp1,isp2,list%dbody(i),rep1)
            rep=rep .or. rep1
          enddo
          if(rep)then
            kx=kxcomposev(isp0)
          endif
        elseif(rep)then
          kx=kxavaloc(-1,list%nl,klr)
          klr%head=dtfcopy(dtastk(isp0))
          klr%dbody(1:list%nl)=list%dbody(1:list%nl)
        endif
        isp=isp0-1
      elseif(ktfpatq(k,pat))then
        k1=tfredefsymbol(isp1,isp2,pat%expr,rep)
        k2=tfredefsymbol(isp1,isp2,pat%head,rep1)
        rep=rep .or. rep1
        kd=tfredefsymbol(isp1,isp2,pat%default,rep1)
        rep=rep .or. rep1
        if(pat%sym%loc .ne. 0)then
          ks=pat%sym%alloc
          kas=ktfaddr(ks)
          do i=isp1,isp2-1,2
            if(ilist(2,ktastk(i+1)-1) .eq. pat%sym%gen .and.
     $           klist(ktastk(i+1)) .eq. klist(kas))then
              kx=kxpcopyss(k1,k2,k_descr(ktfsymbol+ktastk(i)),kd)
              rep=.true.
              return
            endif
          enddo
        else
          ks%k=0
        endif
        if(rep)then
          kx=kxpcopyss(k1,k2,ks,kd)
        endif
      elseif(ktfsymbolq(k,sym))then
        do i=isp1,isp2-1,2
          if(ilist(2,ktastk(i+1)-1) .eq. sym%gen .and.
     $         klist(ktastk(i+1)) .eq. sym%loc)then
            kx%k=ktfsymbol+ktastk(i)
            rep=.true.
            return
          endif
        enddo
      endif
      return
      end
