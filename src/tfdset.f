      subroutine tfdset(k,kadi,kx,karg)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) k,kx,kr,karg,kargr
      type (sad_list), pointer :: larg,largl,largd,largr
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash
      integer*8 ktdaloc,kap,kadi,
     $     ktdhtaloc,kad,kad1,kad0,kad10,kan,kih
      integer*4 m,ihash,itfhasharg,isp0,
     $     mstk0,iop0,itflistmat,irtc,minhash,maxhash
      parameter (minhash=7,maxhash=32767)
      logical*4 tfconstpatternqk,tfsamelistqo,tfmaxgenerationq

c     Initialize to avoid compiler warning
      kad0=0
c
      call descr_list(karg,larg)
      kad1=ksad_loc(kadi)
      kad=kadi
      if(tfconstpatternqk(karg%k))then
        if(.not. tfmaxgenerationq(kad))then
          if(k%k .ne. ktfref)then
            ihash=itfhasharg(karg,minhash)
            kih=ktdhtaloc(kad1,kad,minhash)
            call loc_defhash(kih,dhash)
            kad1=ksad_loc(dhash%hash(ihash))
            kap=ktdaloc(int8(0),kad1,int8(0),
     $           ktfref,karg,k,karg,.false.)
            dhash%attr=ior(dhash%attr,1)
            kx=k
          else
            kx%k=ktfoper+mtfnull
          endif
          return
        else
          ihash=itfhasharg(karg,ilist(2,kad+2))
          call loc_defhash(kad,dhash)
          dhash%attr=ior(dhash%attr,1)
          kad1=ksad_loc(dhash%hash(ihash))
          kad=klist(kad1)
          do while(kad .ne. 0)
            call loc_deftbl(kad,dtbl)
            call descr_list(dtbl%arg,largd)
            if(larg%nl .eq.largd%nl)then
              if(tfsamelistqo(larg,largd))then
                if(k%k .eq. ktfref)then
                  go to 8000
                else
                  kan=ktdaloc(kad,int8(0),int8(0),
     $                 ktfref,karg,k,karg,.false.)
                  kx=k
                  return
                endif
              endif
            endif
            kad1=kad
            kad=dtbl%next
          enddo
          if(k%k .ne. ktfref)then
            kad=klist(kad1)
            kan=ktdaloc(int8(0),kad1,kad,ktfref,karg,k,karg,.false.)
            kx=k
          else
            kx%k=ktfoper+mtfnull
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
      if(irtc .ne. 0)then
        return
      endif
      call descr_list(kargr,largr)
      if(rlist(iaxpriority) .eq. 0.d0)then
        mstk0=mstk
        iop0=iordless
        iordless=0
        do while(kad .gt. 0)
          call loc_deftbl(kad,dtbl)
          call descr_list(dtbl%argc,largl)
          m=itflistmat(kargr,largl)
          isp=isp0
          if(m .ge. 0)then
            call tfunsetpattbl(dtbl)
          endif
c          call tfdebugprint(ktflist+kargr,'tfdset',1)
c          call tfdebugprint(ktflist+kal  ,'=?=   ',1)
 1        if(m .eq. 1)then
            if(k%k .eq. ktfref)then
              mstk=mstk0
              iordless=iop0
              go to 8000
            else
              kan=ktdaloc(kad,int8(0),int8(0),
     $             k,karg,kr,kargr,.true.)
            endif
            iordless=iop0
            mstk=mstk0
            kx=k
            return
          else
            if(largr%nl .eq. largl%nl)then
              if(tfsamelistqo(largr,largl))then
                m=1
                go to 1
              endif
            endif
            if(m .eq. 0)then
              if(kad10 .eq. 0)then
                kad10=kad1
                kad0=kad
              endif
              exit
            endif
            kad1=kad
            kad=dtbl%next
          endif
        enddo
        mstk=mstk0
        iordless=iop0
      endif
      if(k%k .ne. ktfref)then
        if(kad10 .ne. 0)then
          kad1=kad10
          kad=kad0
        endif
c        call tfdebugprint(ktflist+kargr,'tfdset-7',3)
c        call tfdebugprint(kr,':= ',3)
        kan=ktdaloc(int8(0),kad1,kad,k,karg,kr,kargr,.true.)
        kx=k
      else
        kx%k=ktfoper+mtfnull
      endif
      return
 8000 klist(kad1)=klist(kad)
      if(klist(kad) .ne. 0)then
        klist(klist(kad)+1)=kad1
      endif
      call tfcleardaloc(kad)
      call tfree(kad)
      kx%k=ktfoper+mtfnull
      return
      end

      subroutine tfcleardaloc(kad)
      use tfstk
      use tfcode
      implicit none
      type (sad_deftbl), pointer :: dtbl
      integer*8 kad
      call loc_deftbl(kad,dtbl)
      call tflocal1d(dtbl%arg)
      call tflocal1d(dtbl%argc)
      call tflocal(dtbl%body%k)
      call tflocal(dtbl%bodyc%k)
      return
      end

      subroutine tfreplacedef(k,karg,kr,kargr,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kr,kx1,karg,kargr
      integer*4 irtc,isp0,ispr,i,ipsbase,nrule
      logical*4 rep,same
      isp0=isp
      call tfinitpat(isp0,karg)
      if(isp .eq. isp0)then
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
        same=same .and.
     $       ktastk(i+1) .eq. klist(ktfaddr(ktastk(isp)))
      enddo
      call tfresetpat(karg)
      if(same)then
        isp=isp0
        kr=k
        kargr=karg
        irtc=0
        return
      endif
      nrule=(isp-ispr)/2
      call tfreplacesymbolstk(karg,ispr,nrule,kx1,
     $     .false.,rep,irtc)
      if(irtc .ne. 0)then
        isp=isp0
        return
      endif
      kargr=kx1
      call tfreplacesymbolstk1(k,ispr,nrule,kr,.false.,rep,irtc)
      isp=isp0
      return
      end

      integer*8 function ktdaloc(kan0,kb0,kn0,
     $     k,karg,kr,kargr,tbl)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) k,kr,kx,karg,kargr,karg1
      type (sad_list), pointer :: klx
      type (sad_deftbl), pointer :: dtbl
      integer*8 kb0,kn0,kan,kan0,ktalocr,kb,kn
      integer*4 npat,isp0,i
      logical*4 tbl,rep,tfconstqk,sym
      isp0=isp
      if(tbl)then
        call tfdefsymbol(kargr,kx,rep,sym)
        karg1=kx
        call descr_list(karg1,klx)
        call tfinitpatlist(isp0,klx)
        npat=(isp-isp0)/2
      else
        karg1=karg
        npat=0
      endif
      kan=kan0
      kn=kn0
      kb=kb0
      if(kan .eq. 0)then
        kan=ktalocr(npat+8)
        call loc_deftbl(kan,dtbl)
        dtbl%next=kn
        dtbl%prev=kb
        klist(kb)=kan
        if(kn .ne. 0)then
          klist(kn+1)=kan
        endif
      else
        call loc_deftbl(kan,dtbl)
        call tfcleardaloc(kan)
        if(npat .gt. dtbl%npat)then
          write(*,*)'ktdaloc-implementation error ',npat,dtbl%npat
          call forcesf()
        endif
      endif
      dtbl%npat=npat
      dtbl%attr=0
      dtbl%arg=dtfcopy1(karg)   ! orifinal arg
      dtbl%argc=dtfcopy1(karg1)   ! replaced arg
      dtbl%body=dtfcopy(k)       ! original body
      dtbl%pat=-1
      if(.not. ktfobjqd(kr))then
        dtbl%bodyc=kr
      else
        dtbl%bodyc=dtfcopy1(kr)    ! replaced body
        if(.not. tfconstqk(kr%k))then
          dtbl%pat=0
        endif
      endif
      dtbl%compile=0.d0   ! compile level
      do i=1,npat
        dtbl%pattbl(i)=dtastk(isp0+i*2-1)
c        write(*,*)'loc.cont ',klist(klist(ktfaddr(klist(kan+7+i))+7)-3)
      enddo
      isp=isp0
      ktdaloc=kan
      return
      end
      
      recursive subroutine tfdefsymbol(k,kx,rep,sym)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k
      type (sad_list), pointer :: list,klx
      type (sad_symbol), pointer :: ks,ks1
      integer*4 i,isp0,m
      logical*4 rep,rep1,sym
      rep=.false.
      sym=.false.
      kx=k
      if(ktflistqd(k,list))then
        if(iand(lnodefsymbol,list%attr) .ne. 0)then
          return
        endif
        isp0=isp
        isp=isp+1
        call tfdefsymbol(list%dbody(0),dtastk(isp),rep,sym)
        m=list%nl
        if(ktfreallistqo(list))then
          if(rep)then
            kx=kxavaloc(-1,m,klx)
            klx%head=ktfcopy(ktastk(isp))
            klx%body(1:m)=list%body(1:m)
          elseif(.not. sym)then
            list%attr=ior(lnodefsymbol,list%attr)
          endif
          isp=isp0
          return
        endif
        do i=1,list%nl
          isp=isp+1
          call tfdefsymbol(list%dbody(i),dtastk(isp),rep1,sym)
          rep=rep .or. rep1
        enddo
        if(rep)then
          kx=kxcomposer(isp0+1)
        elseif(.not. sym)then
          list%attr=ior(lnodefsymbol,list%attr)
        endif
        isp=isp0
        return
      elseif(ktfsymbolqd(k,ks))then
        if(ks%override .eq. 0)then
          sym=.true.
          if(ks%gen .ne. maxgeneration)then
            call tfsydef(ks,ks1)
            kx%k=ktfsymbol+ksad_loc(ks1%loc)
            rep=.true.
          endif
        endif
      endif
      return
      end

      integer*8 function ktdhtaloc(kad1,kad0,nhash)
      use tfstk
      use tfcode
      implicit none
      type (sad_defhash), pointer :: dhash
      integer*8 kan,kad1,kad0,kad,ktalocr
      integer*4 nhash
      kad=kad0
      kan=ktalocr(nhash+5)
      call loc_defhash(kan,dhash)
      dhash%attr=0
      dhash%next=kad
      dhash%prev=kad1
      dhash%gen=maxgeneration
      dhash%nhash=nhash
      dhash%hash(0:nhash)=0
      klist(kad1 )=kan
      if(kad .ne. 0)then
        klist(kad+1)=kan
      endif
      ktdhtaloc=kan
      return
      end

      subroutine tfdeval(isp1,kad00,kx,iup,def,ev,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx,kad00
      type (sad_list), pointer :: larg,klx
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash
      integer*8 kad0,kad,kadv,kadv0,kap,ktfrehash,khash
      integer*4 isp1,iup,irtc,is,itfhasharg,
     $     m,mstk0,im,i,itfseqmatstk,isp0,ns,nh,iord0
      logical*4 ev,ordless,tfsameqk,def,tfmaxgenerationq,tfordlessq
      ev=.false.
      im=isp1+1
      ordless=iup .eq. 0 .and. tfordlessq(ktastk(isp1))
     $     .and. isp .gt. im
      kad0=ktfaddr(kad00)
      kad=klist(iup+kad0-6)
      isp0=isp
      if(tfmaxgenerationq(kad))then
        call loc_defhash(kad,dhash)
        if(iand(dhash%attr,1) .eq. 0)then
          go to 20
        endif
        nh=dhash%nhash
        itastk2(1,isp1)=isp
        khash=ksad_loc(dhash%hash(itfhasharg(ktfref+isp1,nh)))
        kadv=klist(khash)
        if(kadv .ne. 0)then
          is=1
          kadv0=kadv
          ns=0
          do while(kadv .ne. 0)
            call loc_deftbl(kadv,dtbl)
            kap=ktfaddr(dtbl%arg)
            call loc_list(kap,larg)
            m=larg%nl
            if(m .eq. isp-isp1)then
              if(.not. tfsameqk(ktastk(isp1),larg%head))then
                go to 10
              endif
              if(m .ne. 0)then
                if(ordless)then
                  iord0=iordless
                  iordless=0
                  mstk0=mstk
                  if(.not. itfseqmatstk(im,isp0,larg%body(1),
     $                 m,is,ktftype(larg%body(1)) .eq. ktfoper,
     $                 kap) .ge. 0)then
                    iordless=iord0
                    mstk=mstk0
                    isp=isp0
                    go to 10
                  endif
                  iordless=iord0
                  mstk=mstk0
                else
                  if(ktfreallistqo(larg))then
                    do i=1,m
                      if(ktastk(isp1+i) .ne. larg%body(i))then
                        go to 10
                      endif
                    enddo
                  else
                    do i=1,m
                      if(.not. tfsameqk(ktastk(isp1+i),
     $                     larg%body(i)))then
                        go to 10
                      endif
                    enddo
                  endif
                endif
              endif
              if(kadv .ne. kadv0)then
                if((nh .lt. 127 .and. ns .gt. 3)
     $               .or. (nh .lt. 1023 .and. ns .gt. 7)
     $               .or. (ns .gt. 15 .and. nh .le. 32767))then
                  kad=ktfrehash(kad,iup)
                else
                  klist(dtbl%prev)=dtbl%next
                  if(dtbl%next .ne. 0)then
                    klist(dtbl%next+1)=dtbl%prev
                  endif
                  dtbl%next=kadv0
                  klist(kadv0+1)=kadv
                  dtbl%prev=khash
                  klist(khash)=kadv
                endif
              endif
              ev=.true.
              if(dtbl%compile .lt. rlist(levelcompile))then
                call tfsetarg(dtbl,irtc)
                if(irtc .ne. 0)then
                  return
                endif
              endif
              irtc=0
              if(def)then
                kx%k=ktfref+ksad_loc(dtbl%bodyc%k)
                return
              endif
              kx=dtbl%bodyc
              if(dtbl%pat .ne. -1)then
                if(ktflistqd(kx,klx))then
                  call tfleval(klx,kx,.true.,irtc)
                elseif(ktfsymbolqd(kx))then
                  call tfsyeval(kx,kx,irtc)
                elseif(ktfpatqd(kx))then
                  call tfpateval(kx,kx,irtc)
                endif
                if(irtc .ne. 0)then
                  call tfcatchreturn(0,kx,irtc)
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
      if(kad .ne. 0 .and. .not. def)then
        call tfdeval1(isp1,kad,kx,ev,ordless,ns,irtc)
      else
        irtc=-1
      endif
      return
      end

      integer*4 function itfhasharg(k,nh)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      integer*8 ka
      integer*4 ih,m,itfhash,nh,ix(2),ih1
      integer*2 h(2),h1(2)
      real*8 x,ha
      parameter (ha=7.d0**5+1.d0/7.d0**5)
      equivalence (h,ih),(x,ix),(h1,ih1)
      if(ktflistqd(k,kl))then
        m=kl%nl
        if(m .eq. 0)then
          ih=0
        elseif(m .eq. 1)then
          ih=itfhash(kl%dbody(1),nh)
        else
          ih=itfhash(kl%dbody(1),nh)+itfhash(kl%dbody(m),nh)
        endif
        ih=ih+m
      elseif(ktfrefqd(k,ka))then
        if(ka .eq. 0)then
          ih=0
        else
          m=itastk2(1,ka)-int(ka)
          if(m .eq. 0)then
            ih=0
          elseif(m .eq. 1)then
            ih=itfhash(dtastk(ka+1),nh)
          else
            ih=itfhash(dtastk(ka+1),nh)+itfhash(dtastk(ka+m),nh)
          endif
          ih=ih+m
        endif
      else
        call tfdebugprint(k,
     $       'itfhasharg-implementation error: ',3)
        rlist(7)=0.d0
      endif
      ih=(h(1)+h(2))
      itfhasharg=iand(ih,nh)
      return
      end

      recursive integer*4 function itfhash(k,nh) result(ih)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 m,itfhash1,nh,ix(2),ih1
      integer*2 h(2)
      real*8 x,ha,v
      parameter (ha=7.d0**5+1.d0/7.d0**5)
      equivalence (h,ih1),(x,ix)
      if(ktfrealqd(k,v))then
        ih=int(v*ha)
      elseif(ktfoperqd(k,ka))then
        ih=int(ka)
      elseif(ktflistqd(k,kl))then
        m=kl%nl
        if(m .eq. 0)then
          ih=0
        elseif(m .eq. 1)then
          ih=itfhash1(kl%dbody(1))
        else
          ih=itfhash1(kl%dbody(1))+itfhash1(kl%dbody(m))
        endif
        ih=ih+itfhash(kl%dbody(0),nh)+m
      elseif(ktfrefqd(k,ka))then
        if(ka .eq. 0)then
          ih=0
        else
          m=itastk2(1,ka)-int(ka)
          if(m .eq. 0)then
            ih=0
          elseif(m .eq. 1)then
            ih=itfhash1(dtastk(ka+1))
          else
            ih=itfhash1(dtastk(ka+1))+itfhash1(dtastk(ka+m))
          endif
          ih=ih+itfhash(dtastk(ka),nh)+m
        endif
      elseif(ktfsymbolqd(k,sym))then
        ih=int(sym%loc)
      elseif(ktfstringqd(k))then
        ka=ktfaddrd(k)
        ih=ilist(1,ka)+ilist(1,ka+1)+ilist(2,ka+1)
      else
        ih=0
      endif
      ih1=ih
      ih=(h(1)+h(2))
      ih=iand(ih,nh)
      return
      end

      recursive integer*4 function itfhash1(k) result(ih)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      type (sad_symbol), pointer :: sym
      integer*8 ka
      integer*4 ix(2)
      real*8 x,ha,v
      parameter (ha=7.d0**5+1.d0/7.d0**5)
      equivalence (x,ix)
      if(ktfrealqd(k,v))then
        ih=int(v*ha)
      elseif(ktfoperqd(k,ka))then
        ih=int(ka)
      elseif(ktflistqd(k,kl))then
        ih=itfhash1(kl%dbody(0))+kl%nl
      elseif(ktfrefqd(k,ka))then
        if(ka .eq. 0)then
          ih=0
        else
          ih=itfhash1(dtastk(ka))+itastk2(1,ka)-int(ka)
        endif
      elseif(ktfsymbolqd(k,sym))then
        ih=int(sym%loc)
      elseif(ktfstringqd(k))then
        ka=ktfaddrd(k)
        ih=ilist(1,ka)+ilist(1,ka+1)+ilist(2,ka+1)
      else
        ih=0
      endif
      return
      end

      subroutine tfdeval1(isp1,kad0,kx,ev,ordless,ns,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx,kh,kal
      type (sad_list), pointer :: larg
      type (sad_deftbl), pointer :: dtbl
      integer*8 kad0,kad,kap,kpp
      integer*4 isp1,irtc,mstk0,mat,im,itfpmat,itfseqmatstk,isp0,
     $     iop0,itfseqmatstk1,ns
      logical*4 ev,ordless
      isp0=isp
      mstk0=mstk
      iop0=iordless
      iordless=0
      im=isp1+1
      ns=0
      kad=kad0
      do while(kad .ne. 0)
        call loc_deftbl(kad,dtbl)
        mat=-1
        kap=ktfaddr(dtbl%argc)
        call loc_list(kap,larg)
        kh=larg%dbody(0)
        mat=itfpmat(ktastk(isp1),kh)
        if(mat .ge. 0)then
          if(larg%nl .eq. 1)then
            if(itfseqmatstk1(im,isp0,larg%body(1)) .lt. 0)then
              go to 30
            endif
          else
            if(ordless)then
              kpp=kap
            else
              kpp=0
            endif
            if(itfseqmatstk(im,isp0,larg%body(1),
     $           larg%nl,1,ktfreallistqo(larg),kpp)
     $           .lt. 0)then
              go to 30
            endif
          endif
          ev=.true.
          if(dtbl%compile .lt. rlist(levelcompile))then
            call tfsetarg(dtbl,irtc)
            if(irtc .ne. 0)then
              call tfunsetpattbl(dtbl)
              isp=isp0
              iordless=iop0
              mstk=mstk0
              return
            endif
          endif
          iordless=iop0
          mstk=mstk0
          if(dtbl%npat .ne. 0)then
            kal=dtfcopy1(dtbl%argc)
            call tfevalwitharg(dtbl,dtbl%bodyc,kx,irtc)
            call tflocal1d(kal)
          else
            call tfeevalref(dtbl%bodyc%k,kx,irtc)
          endif
          if(irtc .ne. 0)then
            call tfcatchreturn(0,kx,irtc)
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

      subroutine tfevalwitharg(dtbl,k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,kh
      type (sad_deftbl) dtbl
      integer*8 ks
      integer*4 irtc,isp0,m
      logical*4 av
      isp0=isp
      call tfevalwitharg1(k,ks,m,kh,kx,av,irtc)
      call tfunsetpattbl(dtbl)
      if(irtc .ne. 0)then
        return
      endif
      if(kx%k .eq. ktfref)then
        call tfseval(ks,m,kh,kx,.true.,av,.true.,irtc)
      else
        call tfeevalref(kx,kx,irtc)
      endif
      call tfcatchreturn(0,kx,irtc)
      isp=isp0
      return
      end

      subroutine tfevalwitharg1(k,ks,m,kh,kx,av,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx,k,kh,kp(0:2),ksy
      type (sad_list), pointer :: list,klh
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symdef
      integer*8 i,kts,kah,ks1,ks
      integer*4 irtc,m,isp1,itfmessageexp,isp0
      logical*4 rep,rep1,tfreplacearg,tfreplaceargstk, av
      irtc=0
      kx=k
      if(ktflistqd(k,list))then
        kx%k=ktfref
        m=list%nl
        isp1=isp
        kh=list%dbody(0)
        if(ktfreallistqo(list))then
          if(ktflistqd(kh,klh))then
            if(iand(larglist,klh%attr) .eq. 0)then
              kx=k
              return
            endif
          elseif(ktfnonsymbolqd(kh) .and. .not. ktfpatqd(kh))then
            kx=k
            return
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc .ne. 0)then
            return
          endif
          ks=ktfaddrd(k)
          av=.true.
          return
        else
          if(ktflistqd(kh,klh))then
            if(iand(larglist,klh%attr) .eq. 0)then
              go to 20
            endif
          elseif(ktfnonsymbolqd(kh) .and. .not. ktfpatqd(kh))then
            go to 20
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc .ne. 0)then
            return
          endif
 20       if(ktfoperqd(kh))then
            kah=ktfaddrd(kh)
            if(iget_fun_id(ktfaddrd(kh)) .eq. nfunif)then
              call tfreplaceifarg(list,kx,rep,irtc)
              return
            endif
          endif
          if(iand(list%attr,lnopatlist) .ne. 0)then
            ks=ktfaddrd(k)
            av=.false.
            return
          endif
          isp=isp1+1
          dtastk(isp)=kh
          ks=isp
          do i=1,m
c            call tfdebugprint(list%body(i),'evwarg1',3)
            rep1=tfreplaceargstk(list%dbody(i),irtc)
c            call tfdebugprint(list%body(i),'==> ',3)
            if(irtc .ne. 0)then
              return
            endif
          enddo
          if(rlist(iaximmediate) .ne. 0.d0)then
            if(ktfoperqd(kh))then
              if(kah .gt. mtfend)then
                ks1=isp
                call tfcomposefun(isp1+1,kah,kx,.false.,irtc)
              else
                call tfcomposeoper(isp1+1,kah,kx,.false.,isp0,irtc)
                ks1=isp0
              endif
              if(irtc .eq. 0)then
                if(kx%k .ne. ktfref)then
                  return
                endif
                ks=ks1
              elseif(irtc .gt. 1)then
                return
              else
                kx%k=ktfref
              endif
              irtc=0
            endif
          endif
          m=int(isp-ks)
          av=.true.
          do i=ks+1,isp
            if(ktfnonrealq(ktastk(i)))then
              av=.false.
              exit
            endif
          enddo
          return
        endif
      elseif(ktfpatqd(k,pat))then
        rep=.false.
        if(pat%sym%loc .ne. 0)then
          rep=tfreplacearg(pat%sym%alloc,ksy,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(rep .and. ktfnonsymbolqd(ksy))then
            irtc=itfmessageexp(999,'General::reppat',ksy)
            return
          endif
        else
          ksy%k=0
        endif
        rep1=tfreplacearg(pat%expr,kp(0),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%head,kp(1),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%default,kp(2),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        if(rep)then
          kx=kxpcopyss(kp(0),kp(1),ksy,kp(2))
        endif
      elseif(ktfsymbolqdef(k%k,symdef))then
        if(symdef%sym%override .eq. 1)then
          kx=symdef%value
          kts=ktftype(kx%k)
          if(kts .eq. ktflist .or. kts .eq. ktfref)then
            call tfgetstkstk(kx,rep)
          endif
        endif
      endif
      return
      end

      logical*4 function tfreplacearg(k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx
      integer*4 irtc,isp0
      logical*4 tfreplaceargstk
      isp0=isp
c      call tfdebugprint(k,'tfreparg-in',1)
      tfreplacearg=tfreplaceargstk(k,irtc)
      if(irtc .ne. 0)then
        isp=isp0
        return
      endif
      if(isp .eq. isp0)then
        kx=dxnull
      elseif(isp .eq. isp0+1)then
        kx=dtastk(isp)
      else
        call tfsequence(isp0,isp,kx)
      endif
c      call tfdebugprint(kx,'tfreparg-out',1)
      isp=isp0
      return
      end

      recursive logical*4 function tfreplaceargstk(k,irtc)
     $     result(lx)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) k,ks,kx,kh,kp(0:2)
      type (sad_list), pointer :: list,klx,klh
      type (sad_pat), pointer :: pat
      type (sad_symdef), pointer :: symd
      integer*8 kah
      integer*4 irtc,i,m,isp1,itfmessageexp,isp0
      logical*4 rep,rep1,tfreplacearg
      irtc=0
      rep=.false.
      lx=.false.
      if(ktflistqd(k,list))then
        kx=k
        isp1=isp
        kh=list%dbody(0)
        if(ktfreallistqo(list))then
          if(ktflistqd(kh,klh))then
            if(iand(larglist,klh%attr) .eq. 0)then
              go to 10
            endif
          elseif(ktfnonsymbolqd(kh) .and. .not. ktfpatqd(kh))then
            go to 10
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(rep)then
            m=list%nl
            kx=kxavaloc(-1,m,klx)
            klx%body(1:m)=list%body(1:m)
            klx%dbody(0)=dtfcopy(kh)
          endif
 10       isp=isp1+1
          dtastk(isp)=kx
          return
        else
          if(ktflistqd(kh,klh))then
            if(iand(larglist,klh%attr) .eq. 0)then
              go to 20
            endif
          elseif(ktfnonsymbolqd(kh) .and. .not. ktfpatqd(kh))then
            go to 20
          endif
          rep=tfreplacearg(kh,kh,irtc)
          if(irtc .ne. 0)then
            return
          endif
 20       isp=isp1+1
          dtastk(isp)=kh
          if(ktfoperqd(kh,kah))then
            if(iget_fun_id(kah) .eq. nfunif)then
              call tfreplaceifarg(list,kx,rep,irtc)
              isp=isp1+1
              dtastk(isp)=kx
              lx=rep
              return
            elseif(kah .eq. mtfnull)then
              isp=isp1
            endif
          endif
          if(iand(list%attr,lnopatlist) .ne. 0)then
            if(rep)then
              call tfgetllstkall(list)
            elseif(kh%k .eq. ktfoper+mtfnull)then
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
c              call tfdebugprint(list%body(i),'repargstk',3)
              rep1=tfreplaceargstk(list%dbody(i),irtc)
              if(irtc .ne. 0)then
                return
              endif
c              call tfdebugprint(list%body(i),'==>',3)
              rep=rep .or. rep1
            enddo
          endif
          lx=rep
          if(kh%k .eq. ktfoper+mtfnull)then
            return
          endif
          if(rlist(iaximmediate) .ne. 0.d0)then
            if(ktfoperqd(kh))then
              if(kah .gt. mtfend)then
                call tfcomposefun(isp1+1,kah,kx,.false.,irtc)
              else
                call tfcomposeoper(isp1+1,kah,kx,.true.,isp0,irtc)
              endif
              if(irtc .eq. 0)then
                lx=.true.
                if(ktflistqd(kx,klx))then
                  if(klx%head .eq. ktfoper+mtfnull)then
                    isp=isp1
                    call tfgetllstkall(klx)
                    return
                  endif
                endif
                isp=isp1+1
                dtastk(isp)=kx
                return
              elseif(irtc .gt. 1)then
                return
              endif
            elseif(ktfstringqd(kh))then
              if(isp .eq. isp1+1 .or. isp .eq. isp1+2)then
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
      elseif(ktfpatqd(k,pat))then
        rep=.false.
        kx=k
        if(pat%sym%loc .ne. 0)then
          rep=tfreplacearg(pat%sym%alloc,ks,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(rep .and. ktfnonsymbolqd(ks))then
            irtc=itfmessageexp(999,'General::reppat',ks)
            return
          endif
        else
          ks%k=0
        endif
        rep1=tfreplacearg(pat%expr,kp(0),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%head,kp(1),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        rep1=tfreplacearg(pat%default,kp(2),irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
c          if(ktftype(kp(i)) .eq. ktfref)then
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
        if(symd%sym%override .eq. 1)then
          ks=symd%value
c          call tfdebugprint(k,'tfrepargstk-symbol',1)
c          call tfdebugprint(ks,':= ',1)
          if(ktflistqd(ks) .or. ktfrefqd(ks))then
            call tfgetstkstk(ks,rep)
          else
            isp=isp+1
            dtastk(isp)=ks
          endif
          rep=.true.
        else
c          call tfdebugprint(k,'tfrepargstk-symbol-reject',1)
c          write(*,*)'with ',symd%sym%override
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

      recursive subroutine tfgetstkstk(ks,rep)
      use tfstk
      implicit none
      type (sad_descriptor) ks,ki
      type (sad_list), pointer :: kl,kli
      integer*8 i,ka
      logical*4 rep,rep1
      if(ktfrefqd(ks,ka))then
        rep=.true.
        if(ka .gt. 3)then
          do i=ka+1,itastk2(1,ka)
            ki=dtastk(i)
            if(ktfrefqd(ki))then
              call tfgetstkstk(ki,rep1)
            elseif(ktflistqd(ki,kli))then
              if(kli%head .eq. ktfoper+mtfnull)then
                call tfgetllstkall(kli)
              else
                isp=isp+1
                dtastk(isp)=ki
              endif
            else
              isp=isp+1
              dtastk(isp)=ki
            endif
          enddo
          return
        endif
        isp=isp+1
        dtastk(isp)=dxnull
        return
      elseif(ktflistqd(ks,kl))then
        if(kl%head .eq. ktfoper+mtfnull)then
          call tfgetllstkall(kl)
          rep=.true.
          return
        endif
      endif        
      rep=.false.
      isp=isp+1
      dtastk(isp)=ks
      return
      end

      subroutine tfreplaceifarg(list,kx,rep,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ki
      type (sad_list) list
      type (sad_list), pointer :: klx
      integer*4 irtc,i,isp1,j
      logical*4 rep,rep1,tfreplaceargstk
      irtc=0
      isp1=isp
      if(list%nl .eq. 0)then
        if(rep)then
          kx=kxaaloc(-1,0,klx)
          klx%head=ktastk(isp1)
        else
          kx%k=ktflist+ksad_loc(list%head)
        endif
        return
      endif
      ki=list%dbody(1)
      select case(ktftype(ki%k))
      case (ktflist)
        if(iand(larglist,list%attr) .ne. 0)then
          rep1=tfreplaceargstk(ki,irtc)
          if(irtc .ne. 0)then
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
        if(irtc .ne. 0)then
          isp=isp1
          return
        endif
        rep=rep .or. rep1
      case default
        isp=isp+1
        dtastk(isp)=ki
      end select
      i=0
      if(isp .eq. isp1+1)then
        if(ktfrealq(ktastk(isp1+1)))then
          if(ktastk(isp1+1) .ne. 0)then
            i=2
          else
            i=3
          endif
        endif
      else
        rep=.true.
      endif
      if(i .ne. 0)then
        rep=.true.
        j=isp+1
        if(i .le. list%nl)then
          ki=list%dbody(i)
          select case (ktftype(ki%k))
          case (ktflist)
            if(iand(larglist,list%attr) .ne. 0)then
              rep1=tfreplaceargstk(ki,irtc)
              if(irtc .ne. 0)then
                isp=isp1
                return
              endif
              if(j .le. isp)then
                kx=dtastk(j)
              else
                kx%k=ktfoper+mtfnull
              endif
            else
              kx=ki
            endif
          case (ktfsymbol,ktfpat)
            rep1=tfreplaceargstk(ki,irtc)
            if(irtc .ne. 0)then
              isp=isp1
              return
            endif
            if(j .le. isp)then
              kx=dtastk(j)
            else
              kx%k=ktfoper+mtfnull
            endif
          case default
            kx=ki
          end select
        else
          kx%k=ktfoper+mtfnull
        endif
        return
      endif
      do i=2,list%nl
        ki=list%dbody(i)
        select case(ktftype(ki%k))
        case (ktflist)
          if(iand(larglist,list%attr) .ne. 0)then
            rep1=tfreplaceargstk(ki,irtc)
            if(irtc .ne. 0)then
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
          if(irtc .ne. 0)then
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
        call tfcompose(isp1,list%head,kx,irtc)
      else
        kx%k=ktflist+ksad_loc(list%head)
      endif
      return
      end

      subroutine tfunsetpattbl(dtbl)
      use tfstk
      use tfcode
      implicit none
      type (sad_deftbl) dtbl
      integer*4 np,i
      np=dtbl%npat
      if(np .ne. maxgeneration)then
        do i=1,np
          call tfunsetpat(dtbl%pattbl(i))
        enddo
      endif
      return
      end

      subroutine tfsetarg(dtbl,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx
      type (sad_deftbl) dtbl
      type (sad_list), pointer :: klx
      type (sad_pat), pointer :: pat
      integer*4 isp0,i,irtc,nrule1
      logical*4 rep,member,tfconstqk,tfsymbollistqo
      irtc=0
      kx=dtbl%bodyc
      if(ktflistqd(kx,klx))then
        if(.not. tfsymbollistqo(klx))then
          dtbl%compile=1.d100
          return
        endif
      elseif(tfconstqk(kx%k))then
c        call tfdebugprint(kx,'setarg-const',3)
        dtbl%compile=1.d100
        return
      endif
      if(dtbl%compile .eq. 0.d0)then
        if(dtbl%npat .ne. 0)then
          isp0=isp
          do i=1,dtbl%npat
            call descr_pat(dtbl%pattbl(i),pat)
            isp=isp+2
            ktastk(isp-1)=ktfaddr(pat%sym%alloc%k)
            ktastk(isp  )=ktfsymbol+ksad_loc(pat%sym%loc)
          enddo
          call tfsortsymbolstk(isp0,dtbl%npat,nrule1)
          call tfsetuparg(dtbl%bodyc,isp0,nrule1*2,
     $         kx,rep,member,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            return
          endif
        else
          member=.true.
        endif
        if(member)then
          dtbl%compile=rlist(levelcompile)
          if(ktflistqd(kx,klx))then
            if(iand(lmemberlist,klx%attr) .ne. 0)then
              call tfrecompilearg(kx,kx,rep,irtc)
              if(irtc .ne. 0)then
                return
              endif
            endif
          elseif(ktfpatqd(kx))then
            call tfrecompilearg(kx,kx,rep,irtc)
            if(irtc .ne. 0)then
              return
            endif
          endif
        else
          dtbl%compile=1.d100
        endif
      else
        if(ktflistqd(kx,klx))then
          if(iand(lmemberlist,klx%attr) .ne. 0)then
            call tfrecompilearg(kx,kx,rep,irtc)
            if(irtc .ne. 0)then
              return
            endif
            dtbl%compile=rlist(levelcompile)
          endif
        elseif(ktfpatqd(kx))then
          call tfrecompilearg(kx,kx,rep,irtc)
          if(irtc .ne. 0)then
            return
          endif
          dtbl%compile=rlist(levelcompile)
        endif
      endif
      call tflocal(dtbl%bodyc%k)
      dtbl%bodyc=dtfcopy(kx)
      return
      end

      recursive subroutine tfsetuparg(k,ispr,nrule2,kx,rep,member,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) k,kx,ki,k1,ks,kd
      type (sad_list), pointer :: list,klx
      type (sad_pat), pointer :: pat
      integer*8 ka1
      integer*4 irtc,i,m,isp1,ispr,nrule2,itfmessageexp,j
      logical*4 rep,rep1,rep2,member,member1,tfmatchsymstk
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
      elseif(ktflistqd(k,list))then
        call tfsetuparg(list%dbody(0),ispr,nrule2,k1,rep,member,irtc)
        if(irtc .ne. 0)then
          return
        endif
        m=list%nl
        if(rep .and. ktfreallistqo(list))then
          kx=kxavaloc(-1,m,klx)
          klx%body(1:m)=list%body(1:m)
          klx%attr=ior(larglist,list%attr)
          klx%dbody(0)=dtfcopy(k1)
          if(member)then
            klx%attr=ior(lmemberlist,klx%attr)
          else
            klx%attr=iand(-lmemberlist-1,klx%attr)
          endif
          return
        endif
        isp1=isp
        isp=isp+1
        dtastk(isp)=k1
        rep2=.false.
        LOOP_I: do i=1,m
          call tfsetuparg(list%dbody(i),
     $       ispr,nrule2,ki,rep1,member1,irtc)
          if(irtc .ne. 0)then
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
            if(list%head .eq. ktfoper+mtfatt)then
              member=.true.
            endif
          endif
        endif
        if(member)then
          klx%attr=ior(lmemberlist,klx%attr)
        else
          klx%attr=iand(-lmemberlist-1,klx%attr)
        endif
        isp=isp1
      elseif(ktfpatqd(k,pat))then
        rep=.false.
        if(pat%sym%loc .ne. 0)then
          call tfsetuparg(pat%sym%alloc,ispr,nrule2,ks,rep,member,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(rep)then
            if(ktftype(ks%k) .ne. ktfsymbol)then
              irtc=itfmessageexp(999,'General::reppat',k)
              return
            endif
          endif
        else
          ks%k=0
        endif
        k1=pat%expr
        if(.not. ktfrefqd(k1,ka1) .or. ka1 .gt. 3)then
          call tfsetuparg(k1,ispr,nrule2,k1,rep1,member1,irtc)
          member=member .or. member1
          if(irtc .ne. 0)then
            return
          endif
          rep=rep .or. rep1
        endif
        kd=pat%default
        if(ktftype(kd%k) .ne. ktfref)then
          call tfsetuparg(kd,ispr,nrule2,kd,rep1,member1,irtc)
          member=member .or. member1
          if(irtc .ne. 0)then
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

      recursive subroutine tfrecompilearg(k,kx,rep,irtc)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) kx,k,k1,k2,kd
      type (sad_list), pointer :: list,klx
      type (sad_pat), pointer :: pat
      type (sad_symbol), pointer :: sym2
      integer*8 ka1
      integer*4 irtc,i,m,isp1
      logical*4 rep,rep1,rep2
      irtc=0
      rep=.false.
      kx=k
      if(ktflistqd(k,list))then
        if(iand(list%attr,lmemberlist) .eq. 0)then
          return
        endif
        k1=list%dbody(0)
        if(k1%k .eq. ktfoper+mtfhold)then
          return
        endif
        call tfrecompilearg(k1,k1,rep,irtc)
        if(ktfreallistqo(list))then
          if(rep)then
            m=list%nl
            kx=kxavaloc(-1,m,klx)
            klx%body(1:m)=list%body(1:m)
c            call tmov(rlist(ka+1),rlist(kax+1),m)
            klx%attr=ior(larglist,list%attr)
            klx%dbody(0)=dtfcopy(k1)
          endif
          return
        endif
        isp1=isp
        isp=isp+1
        dtastk(isp)=k1
        rep2=.false.
        do i=1,list%nl
          isp=isp+1
          call tfrecompilearg(list%dbody(i),dtastk(isp),rep1,irtc)
          if(irtc .ne. 0)then
            isp=isp1
            return
          endif
          rep2=rep2 .or. rep1
        enddo
        if(list%head .eq. ktfoper+mtfatt)then
          if(isp .eq. isp1+3)then
            k2=dtastk(isp1+2)
            if(ktfsymbolqd(k2,sym2))then
              if(sym2%override .ne. 0)then
                if(iand(sym2%attr,iattrdynamic) .ne. 0)then
                  go to 120
                endif
              endif
c              call tfdebugprint(ktastk(isp1+1),'rcmparg',3)
              call tfatt(isp1+1,kx,.false.,irtc)
              if(irtc .gt. 0)then
                isp=isp1
                return
              elseif(irtc .eq. 0)then
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
      elseif(ktfpatqd(k,pat))then
        k1=pat%expr
        if(ktfrefqd(k1,ka1) .and. ka1 .gt. 3)then
          call tfrecompilearg(k1,k1,rep,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
        kd=pat%default
        call tfrecompilearg(kd,kd,rep1,irtc)
        if(irtc .ne. 0)then
          return
        endif
        rep=rep .or. rep1
        if(rep)then
          kx=kxpcopyss(k1,pat%head,pat%sym%alloc,kd)
        endif
      endif
      return
      end

      integer*8 function ktfrehash(kad0,iup)
      use tfstk
      use tfcode
      implicit none
      type (sad_descriptor) karg
      type (sad_deftbl), pointer :: dtbl
      type (sad_defhash), pointer :: dhash0,dhash
      integer*8 ktdhtaloc,kad,kad0,ka,ka0,kas,i
      integer*4 iup,n,n0,ih,itfhasharg
      integer*4 maxhash
      parameter (maxhash=32767)
      call loc_defhash(kad0,dhash0)
      n0=dhash0%nhash
      if(n0 .ge. maxhash)then
        ktfrehash=kad0
        return
      endif
      n=n0*2+1
      kas=dhash0%prev-iup+6
      kad=ktdhtaloc(dhash0%prev,dhash0%next,n)
      call loc_defhash(kad,dhash)
      dhash%attr=dhash0%attr
      do i=0,dhash0%nhash
        ka=dhash0%hash(i)
        do while(ka .ne. 0)
          call loc_deftbl(ka,dtbl)
          karg=dtbl%arg
          ih=itfhasharg(karg,n)
          ka0=dtbl%next
          dtbl%next=dhash%hash(ih)
          if(dtbl%next .ne. 0)then
            klist(dtbl%next+1)=ka
          endif
          dtbl%prev=ksad_loc(dhash%hash(ih))
          dhash%hash(ih)=ka
          ka=ka0
        enddo
      enddo
      call tfree(kad0)
      ktfrehash=kad
      return
      end

      logical*4 function tfmaxgenerationq(kad)
      use tfstk
      implicit none
      integer*8 kad
      if(kad .ne. 0)then
        tfmaxgenerationq=ilist(1,kad+2) .eq. maxgeneration
      else
        tfmaxgenerationq=.false.
      endif
      return
      end
