      subroutine tfmodule(isp1,kx,module,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ke,kxl1
      type (sad_list), pointer :: lvlist,kle
      type (sad_symdef), pointer :: symdi
      integer*4 irtc,i, isp0,isp1,isp2, itfmessage
      logical*4 module,rep,eval
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+1),lvlist))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
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
        call tfreplacesymbolstk(ktastk(isp0),isp0,(isp2-isp0)/2,ke,
     $     .false.,rep,irtc)
        if(irtc .ne. 0)then
          isp=isp2
          go to 9200
        endif
      endif
      if(eval)then
        if(ktflistq(ke,kle))then
c          call tfdebugprint(ke,'tfmodule',1)
          call tfleval(kle,kx,.true.,irtc)
c          call tfdebugprint(kx,'tfmodule-1',1)
        else
          call tfeevalref(ke,kx,irtc)
        endif
        isp=isp2
      elseif(rep)then
        if(isp2 .gt. isp0)then
c          call tfdebugprint(ktflist+ksad_loc(lvlist%head),'module-1',2)
          call tfredefslist(isp0,isp2,lvlist,kxl1)
c          call tfdebugprint(ktflist+kal1,'==>',2)
        else
          kxl1%k=ktflist+sad_loc(lvlist%head%k)
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
        call descr_sad(dtastk(i),symdi)
        call tfdelete(symdi,.true.,.true.)
      enddo
      isp=isp0
      return
      end

      subroutine tfmlocalv(list,module,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ki,ki1i,ki1,ki2
      type (sad_list) list
      type (sad_list), pointer :: kli1,kli2,kli
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
          if(ktfsymbolqd(ki,symi))then
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
                call tfleval(kli2,ki2,.true.,irtc)
                if(irtc .ne. 0)then
                  go to 200
                endif
              elseif(ktfsymbolqd(ki2) .or. ktfpatqd(ki2))then
                isps=isp
                call tfeevalref(ki2,ki2,irtc)
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
                    if(ktfsymbolqd(ki1i,symi1i))then
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
                      symdv%value%k=ktfcopy(kli2%body(ii))
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
          if(ktfsymbolqd(ki,symi))then
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
                if(ktfsymbolqd(ki1i,symi1i))then
                  isp=isp+2
                  dtastk(isp-1)=ki1i
                  dtastk(isp  )=dxsycopy(symi1i)
                else
                  go to 100
                endif
              enddo
            elseif(ktfsymbolqd(ki1,symi1))then
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

      recursive subroutine tfredefslist(isp1,isp2,lvlist,kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k11,kj,k1
      type (sad_list) lvlist
      type (sad_list), pointer :: kl1,klj
      integer*4 isp1,isp2,i,isp0,m,j
      logical*4 tfsamesymbolqd
      m=lvlist%nl
      if(m .eq. 0)then
        kx%k=ktflist+sad_loc(lvlist%head%k)
        return
      endif
      isp=isp+1
      isp0=isp
      do j=1,m
        kj=lvlist%dbody(j)
        if(ktfsymbolqd(kj))then
          do i=isp1+1,isp2,2
            if(tfsamesymbolqd(dtastk(i),kj))then
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
            call tfredefslist(isp1,isp2,kl1,k11)
            isp=isp+1
            dtastk(isp)=k11
          else
            do i=isp1+1,isp2,2
              if(tfsamesymbolqd(dtastk(i),k1))then
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
      implicit none
      type (sad_descriptor) kx,ki,k1,kxi
      type (sad_list), pointer :: list,listi,klx
      type (sad_symbol), pointer :: symi
      integer*4 isp1,irtc,isp2,i,itfmessage,lg0,n
      logical*4 rep,symbol,tfconstpatternqk,eval
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
                if(.not. tfconstpatternqk(k1%k))then
                  ivstk2(2,isp-1)=1
                elseif(.not. ktfrealq(k1))then
                  ivstk2(2,isp-1)=0
                endif
                symbol=symbol .and. ktfsymbolqd(k1)
                if(listi%head%k .eq. ktfoper+mtfsetdelayed)then
                  ktastk(isp)=listi%body(2)
                  cycle
                else
                  ki=listi%dbody(2)
                endif
              elseif(ktfsymbolqd(k1,symi))then
                ivstk2(2,isp-1)=0
              else
                isp=isp-2
                cycle
              endif
              go to 30
            endif
          endif
          if(eval)then
            symbol=symbol .and. ktfsymbolqd(ki)
            if(.not. tfconstpatternqk(ki%k))then
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
              call tfleval(listi,kxi,.true.,irtc)
            else
              call tfeevalref(ki,kxi,irtc)
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
            call tfreplacesymbolstk(ktastk(isp2),isp2,list%nl,kx,
     $           .true.,rep,irtc)
          else
            call tfinitrule(isp2,list%nl)
            call tfreplacestk(ktastk(isp2),
     $           isp2,list%nl,kx,.true.,rep,irtc)
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
            call tfreplacesymbolstk(ktastk(isp1+1),isp2,n,ktastk(isp),
     $           .true.,rep,irtc)
            isp=isp+1
            call tfreplacesymbolstk1(ktastk(isp1+2),isp2,n,ktastk(isp),
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
          call tfleval(klx,kx,.true.,irtc)
        elseif(ktfsymbolqd(kx) .or. ktfpatqd(kx))then
          call tfeevalref(kx,kx,irtc)
        endif
      endif
      isp=isp2
      return
      end

      recursive subroutine tfredefsymbol(isp1,isp2,k,kx,rep)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx,k,k1,k2,kd,ks
      type (sad_list), pointer :: list
      type (sad_rlist), pointer :: klr
      type (sad_symbol), pointer :: sym
      type (sad_pat), pointer :: pat
      integer*8 kas
      integer*4 isp1,isp2,i,isp0
      logical*4 rep,rep1
      rep=.false.
      if(ktflistq(k,list))then
        isp=isp+1
        isp0=isp
        call tfredefsymbol(isp1,isp2,list%head,dtastk(isp),rep)
        if(ktfnonreallistqo(list))then
          do i=1,list%nl
            isp=isp+1
            call tfredefsymbol(isp1,isp2,list%dbody(i),dtastk(isp),rep1)
            rep=rep .or. rep1
          enddo
          if(rep)then
            kx=kxcomposev(isp0)
          endif
        elseif(rep)then
          kx=kxavaloc(-1,list%nl,klr)
          klr%head=dtfcopy(dtastk(isp0))
          klr%body(1:list%nl)=list%body(1:list%nl)
        endif
        isp=isp0-1
      elseif(ktfpatqd(k,pat))then
        call tfredefsymbol(isp1,isp2,pat%expr,k1,rep)
        call tfredefsymbol(isp1,isp2,pat%head,k2,rep1)
        rep=rep .or. rep1
        call tfredefsymbol(isp1,isp2,pat%default,kd,rep1)
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
      elseif(ktfsymbolqd(k,sym))then
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
