      module pmat
      use tfstk

      contains
      subroutine tflinkedpat(pat0,pat)
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_pat), target :: pat0
      type (sad_pat), pointer, intent(out):: pat
c      logical*4 rep
c      rep=.false.
      pat=>pat0
      do while(associated(pat%equiv))
        pat=>pat%equiv
c        rep=.true.
      enddo
c      if(rep .and. pat%mat == 0)then
c        pat%equiv=>pat0
c        pat=>pat0
c        pat%mat=0
c        pat%value%k=ktfref
c        nullify(pat%equiv)
c      endif
      return
      end subroutine

      integer*4 function itfpmat(k,kp)
      implicit none
      type (sad_descriptor) ,intent(in):: k,kp
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: list
      itfpmat=-1
      if(ktfpatq(kp,pat))then
        itfpmat=itfpatmat(k,pat)
      elseif(ktflistq(kp,list))then
        if(ktfnonlistq(k) .and. iand(list%attr,lsimplepat) /= 0)then
          return
        endif
        itfpmat=itflistmat(k,list)
      elseif(tfsameq(k,kp))then
        itfpmat=1
      endif
      return
      end function

      integer*4 function itfpmatc(k,kp)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(in):: kp
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: klp
      integer*4 mstk0,iop,isp0
      if(ktflistq(kp,klp))then
        itfpmatc=itfpmatcl(k,klp)
      elseif(ktfpatq(kp,pat))then
        iop=iordless
        iordless=0
        isp0=isp
        mstk0=mstk
        call tfinitpat(isp0,kp)
        isp=isp0
        itfpmatc=itfpatmat(k,pat)
        call tfresetpat(kp)
        mstk=mstk0
        isp=isp0
        iordless=iop
      else
        itfpmatc=merge(1,-1,tfsameq(k,kp))
      endif
      return
      end

      integer*4 function itfpmatcl(k,klp)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,intent(inout):: klp
      integer*4 mstk0,iop,isp0
      if(ktfnonlistq(k) .and. iand(klp%attr,lsimplepat) /= 0)then
        itfpmatcl=-1
        return
      endif
      iop=iordless
      iordless=0
      mstk0=mstk
      isp0=isp
      if(iand(klp%attr,lnopatlist) == 0)then
        call tfinitpatlist(isp0,klp)
        isp=isp0
        itfpmatcl=itflistmat(k,klp)
        call tfresetpatlist(klp)
      else
        itfpmatcl=itflistmat(k,klp)
      endif
      mstk=mstk0
      isp=isp0
      iordless=iop
      return
      end

      integer*4 function itflistmat1(k,listp,f)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,intent(inout):: listp
      logical*4 ,intent(in):: f
      itflistmat1=itflistmat(k,listp)
      return
      end

      integer*4 function itfconstpart(list,m,n,kx)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ki
      type (sad_dlist) ,intent(inout):: list
      integer*4 ,intent(in):: n,m
      integer*4 i,is,n1
      itfconstpart=0
      n1=n+1
      is=list%attr
      if(iand(kallnofixedarg,is) /= 0)then
        return
      else
        do i=n1,m
          ki=list%dbody(i)
          if(ktfpatq(ki))then
            list%attr=ior(knofixedarg,list%attr)
          elseif(ktflistq(ki))then
            if(tfconstpatternq(ki))then
              itfconstpart=i
              kx=ki
              return
            endif
            list%attr=ior(knoconstlist,list%attr)
          else
            itfconstpart=i
            kx=ki
            return
          endif
        enddo
        if(n1 == 1)then
          list%attr=ior(kallnofixedarg,list%attr)
        endif
      endif
      return
      end

      integer*4 function itfpatmat(k,pat0)
      use tfcode
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kv1
      type (sad_pat) ,intent(in):: pat0
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: list
      integer*8 kav1
      integer*4 isp1,isp2,i,np,iss
      itfpatmat=-1
      call tflinkedpat(pat0,pat)
      kv1=pat%value
      if(pat%mat == 0)then
        itfpatmat=itfsinglepat(k,pat)
      elseif(ktfrefq(kv1,kav1) .and. kav1 > 3)then
        if(ktflistq(k,list))then
          if(list%head%k == ktfoper+mtfnull)then
            isp1=isp
            call tfgetllstkall(list)
            iss=int(kav1-ispbase)
            isp2=isp
            np=itastk2(1,iss)-iss
            if(isp2-isp1 .ge. np)then
              do i=1,np
                if(.not. tfsameq(dtastk(isp1+i),dtastk(iss+i)))then
                  isp=isp1
                  return
                endif
              enddo
              itfpatmat=1
            endif
          endif
        endif
      elseif(tfsameq(k,pat%value))then
        itfpatmat=1
      endif
      return
      end

      integer*4 function itfseqmatseq(isp1,isp2,kl,mp1)
      implicit none
      type (sad_dlist) ,intent(inout):: kl
      integer*4 ,intent(in):: isp1,isp2,mp1
      if(isp1 > isp2 .and. mp1 > kl%nl)then
        itfseqmatseq=0
      else
        itfseqmatseq=itfseqmatstk21(isp1,isp2,kl%dbody(1),
     $       kl%nl,mp1,ktfreallistq(kl),i00,.false.)
      endif        
c      itfseqmatseq=merge(0,
c     $     itfseqmatstk21(isp1,isp2,kl%dbody(1),
c     $     kl%nl,mp1,ktfreallistq(kl),i00,
c     $     isp1 > isp2 .and. mp1 > kl%nl),
c     $     isp1 > isp2 .and. mp1 > kl%nl)
      return
      end

c$$$      integer*4 function itfseqmat(isp10,isp20,kp,mp,mp1)
c$$$      use tfstk
c$$$      implicit none
c$$$      type (sad_descriptor) ,intent(inout):: kp(1)
c$$$      integer*4 ,intent(in):: isp10,isp20,mp,mp1
c$$$      integer*4 isp1,isp2,itfseqmatstk,itfseqmatstk1
c$$$      isp1=isp10
c$$$      isp2=isp20
c$$$      write(*,*)'seqmat ',mp,mp1,isp1,isp2
c$$$      itfseqmat=merge(0,merge(itfseqmatstk1(isp1,isp2,kp(1)),
c$$$     $     itfseqmatstk(isp1,isp2,kp,mp,mp1,.false.,i00),
c$$$     $     mp == 1 .and. mp1 == 1),
c$$$     $     isp1 > isp2 .and. mp1 > mp)    
c$$$      return
c$$$      end

      integer*4 function itfseqmatstk1(isp10,isp20,kp)
      implicit none
      type (sad_descriptor) ,intent(inout):: kp
      integer*4 ,intent(in):: isp10,isp20
      integer*4 isps,i,ispm,iop0,ispm0,mstk0
      mstk0=mstk
      ispm=isp10-1
      iop0=iordless
      ispm0=ispm
      do
        i=itfseqm(isp10,isp20,kp,ispm,isps,isp20)
        if(i .ge. 0)then
          if(isp20 .lt. isps)then
            itfseqmatstk1=min(1,i)
          else
            itfseqmatstk1=-1
          endif
          if(itfseqmatstk1 .ge. 0)then
            return
          endif
          call tfresetpat(kp)
          if(iordless /= iop0)then
            ispm=ispm0
            cycle
          endif
        else
          iordless=iop0
        endif
        exit
      enddo
      mstk=mstk0
      itfseqmatstk1=-1
      return
      end

      recursive integer*4 function
     $     itfseqmatstk21(isp10,isp20,kp0,map,mp1,realp,kpp,f)
     $      result(ix)
      implicit none
      type (sad_descriptor) ,intent(inout):: kp0(1)
      integer*8 ,intent(in):: kpp
      integer*4 ,intent(in):: isp10,isp20,map,mp1
      logical*4 ,intent(in):: realp,f
      ix=itfseqmatstk(isp10,isp20,kp0,map,mp1,realp,kpp)
      return
      end

      recursive integer*4 function
     $     itfseqmatstk(isp10,isp20,kp0,map,mp1,realp,kpp)
     $      result(ix)
      implicit none
      type (sad_descriptor) ,intent(inout):: kp0(1)
      type (sad_descriptor) kp,k1
      integer*8 ,intent(in):: kpp
      integer*4 ,intent(in):: isp10,isp20,map,mp1
      integer*4 isps,i,ispf,np,ispp,mop,iop,iop0,np0,npd,npdi,
     $     isp1a,isp2a,ispt,mstk0
      logical*4 ,intent(in):: realp
      mop=0
      if(mp1 > map)then
        ix=merge(1,-1,isp20 .lt. isp10)
        return
      endif
      if(realp)then
        kp%k=i00
      else
        kp=kp0(mp1)
      endif
c      kp%k=merge(i00,kp0(mp1)%k,realp)
      mstk0=mstk
      ispf=isp20
      if(kpp == 0)then
        isp1a=isp10
        isp2a=isp20
        if(map == mp1+1)then
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
            do
              i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
              if(i .ge. 0)then
                ix=min(i,itfseqmatstk1(isps,isp2a,kp0(map)))
                if(ix .ge. 0)then
                  return
                endif
                call tfresetpat(kp)
                if(iordless /= iop0)then
                  ispf=ispt
                  cycle
                endif
              else
                iordless=iop0
              endif
              exit
            enddo
          enddo
        elseif(map /= mp1)then
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
            do
              i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
              if(i .ge. 0)then
                ix=min(i,itfseqmatstk(isps,isp2a,
     $               kp0,map,mp1+1,realp,i00))
                if(ix .ge. 0)then
                  return
                endif
                call tfresetpat(kp)
                if(iordless /= iop0)then
                  ispf=ispt
                  cycle
                endif
              else
                iordless=iop0
              endif
              exit
            enddo
          enddo
        else
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
            do
              i=itfseqm(isp1a,isp2a,kp,ispf,isps,isp2a)
              if(i .ge. 0)then
                ix=merge(1,-1,isp2a .lt. isps)
                if(ix .ge. 0)then
                  return
                endif
                call tfresetpat(kp)
                if(iordless /= iop0)then
                  ispf=ispt
                  cycle
                endif
              else
                iordless=iop0
              endif
              exit
            enddo
          enddo
        endif
        mstk=mstk0
        ix=-1
        return
      endif
      iop=iordless
      do
        if(iop /= 0)then
          if(kpp > 0)then
            do while(iop /= 0 .and. ktastk(iop) /= kpp)
              iop=itastk2(1,iop)
            enddo
            if(iop == 0)then
              cycle
            endif
            iordless=iop
          endif
          ispf=itastk2(2,iop)
          isp1a=itastk(1,iop+1)
          isp2a=itastk(2,iop+1)
          np=itastk2(1,iop+1)
        else
          mstk=mstk-2
          iop=mstk+1
          ktastk(iop)=kpp
          itastk2(1,iop)=iordless
          itastk2(2,iop)=ispf
          itastk(1,iop+1)=isp10
          itastk(2,iop+1)=isp20
          itastk2(1,iop+1)=0
          np=0
          iordless=iop
          isp1a=isp10
          isp2a=isp20
          ispf=isp20
        endif
        exit
      enddo
      do
        do while(ispf .ge. isp1a-1)
          iop0=iordless
          ispt=ispf
          ix=-1
          do while (ix .lt. 0)
            ispt=merge(isp2a,ispf,map == mp1)
            i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
            if(i .ge. 0)then
              if(mp1 .ge. map)then
                if(isp2a .lt. isps)then
                  ix=min(1,i)
                endif
              else
c                ix=min(i,merge(itfseqmatstk1(isps,isp2a,kp0(map)),
c     $               itfseqmatstk(isps,isp2a,kp0,map,mp1+1,
c     $               realp,i00),map == mp1+1))
                if(map == mp1+1)then
                  ix=min(i,itfseqmatstk1(isps,isp2a,kp0(map)))
                else
                  ix=min(i,itfseqmatstk(isps,isp2a,kp0,map,mp1+1,realp,i00))
                endif
              endif
              if(ix .ge. 0)then
                itastk2(2,iop)=ispf
                return
              endif
              call tfresetpat(kp)
              if(iordless /= iop0)then
                ispf=ispt
              else
                exit
              endif
            else
              iordless=iop0
              exit
            endif
          enddo
        enddo
        iordless=iop
        np=np+1
        npd=np
        np0=isp2a-isp1a+1
        do i=2,np0+1
          npdi=npd/i
          if(npdi*i /= npd)then
            mop=np0-i
            exit
          endif
          npd=npdi
        enddo
        if(mop .lt. 0)then
          mstk=mstk0
          ix=-1
          iordless=itastk2(1,iordless)
          return
        endif
        itastk2(1,iop+1)=np
        ispp=isp1a+mop
        k1=dtastk(ispp)
        dtastk(ispp)=dtastk(ispp+1)
        dtastk(ispp+1)=k1
        ispf=isp20
      enddo
      end

      subroutine tfstkpat(isp1,isp2,pat)
      use tfcode
      implicit none
      type (sad_pat) ,intent(inout):: pat
      integer*4 ,intent(in):: isp1,isp2
      if(pat%sym%loc /= 0)then
        pat%mat=1
        if(isp1 .ge. isp2)then
          pat%value=dxnull
        elseif(isp2 == isp1+1)then
          pat%value=dtastk(isp2)
        else
          pat%value%k=ktfref+isp1+ispbase
          itastk2(1,isp1)=isp2
        endif
      endif
      return
      end

      recursive subroutine tfunsetpat(k)
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: klp
      if(ktflistq(k,klp))then
        call tfunsetpatlist(klp)
      elseif(ktfpatq(k,pat))then
c        do while(associated(pat%equiv))
c          pat=>pat%equiv
c          pat%mat=0
c          pat%value%k=ktfref
c          call tfunsetpat(sad_descr(pat%equiv))
c        enddo
        pat%mat=0
        pat%value%k=ktfref
        call tfunsetpat(pat%expr)
      endif
      return
      end

      recursive subroutine tfunsetpatlist(klp)
      implicit none
      type (sad_dlist) ,intent(inout):: klp
      integer*4 i
      if(iand(klp%attr,lnopatlist) /= 0)then
        return
      endif
      call tfunsetpat(klp%head)
      if(ktfnonreallistqo(klp))then
        if(iand(klp%attr,knopatarg) == 0)then
          do i=1,klp%nl
            call tfunsetpat(klp%dbody(i))
          enddo
        endif
      endif
      return
      end

      recursive subroutine tfresetpat(kp)
      implicit none
      type (sad_descriptor) ,intent(in):: kp
      type (sad_pat), pointer :: pat
      type (sad_dlist), pointer :: list
      if(ktflistq(kp,list))then
        call tfresetpatlist(list)
      elseif(ktfpatq(kp,pat))then
        call tfresetpatpat(pat)
      endif
      return
      end

      recursive subroutine tfresetpatpat(pat)
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_pat) ,intent(inout):: pat
      call tfresetpat(pat%expr)
      pat%mat=0
      pat%value%k=ktfref
      if(associated(pat%equiv))then
c        pat=pat%equiv
c        nullify(pat%equiv)
      endif
      return
      end

      recursive subroutine tfresetpatlist(list)
      implicit none
      type (sad_dlist) ,intent(inout):: list
      integer*4 i
      if(iand(list%attr,lnopatlist) /= 0)then
        return
      endif
      call tfresetpat(list%head)
      if(ktfnonreallistqo(list))then
        if(iand(list%attr,knopatarg) == 0)then
          do i=1,list%nl
            call tfresetpat(list%dbody(i))
          enddo
        endif
      endif
      return
      end

      recursive subroutine tfinitpat(isp0,k)
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      type (sad_pat), pointer :: pat
      integer*8 kas
      integer*4 ,intent(in):: isp0
      integer*4 i
      if(ktflistq(k,kl))then
        call tfinitpatlist(isp0,kl)
      elseif(ktfpatq(k,pat))then
        loop10: do
          nullify(pat%equiv)
          if(pat%sym%loc /= 0)then
            kas=ktfaddr(pat%sym%alloc)
            do i=isp0+1,isp-1,2
              if(kas == ktastk(i+1) .and.
     $             ktastk(i) /= k%k)then
                call loc_pat(ktfaddr(ktastk(i)),pat%equiv)
                exit loop10
              endif
            enddo
            isp=isp+2
            dtastk(isp-1)=k
            ktastk(isp  )=kas
          endif
          pat%value%k=ktfref
          exit
        enddo loop10
        call tfinitpat(isp0,pat%expr)
        pat%mat=0
      endif
      return
      end

      recursive subroutine tfinitpatlist(isp0,kl)
      implicit none
      type (sad_dlist) ,intent(inout):: kl
      integer*4 ,intent(in):: isp0
      integer*4 isp1,m,i
      if(iand(kl%attr,lnopatlist) /= 0)then
        return
      endif
      isp1=isp
      call tfinitpat(isp0,kl%head)
      if(ktfnonreallistqo(kl))then
        if(iand(kl%attr,knopatarg) == 0)then
          m=kl%nl
          if(kl%head%k == ktfoper+mtfpattest)then
            m=1
          endif
          do i=1,m
            call tfinitpat(isp0,kl%dbody(i))
          enddo
          if(isp == isp1)then
            kl%attr=ior(kl%attr,knopatarg)
          endif
        endif
      endif
      if(isp == isp1)then
        kl%attr=ior(kl%attr,lnopatlist)
      endif
      return
      end

      logical*4 function tfordlessq(k)
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_funtbl), pointer ::fun
      type (sad_symbol), pointer ::sym
      if(ktfoperq(k))then
        call c_f_pointer(c_loc(klist(klist(ifunbase+ktfaddrd(k))-9)),
     $       fun)
        tfordlessq=iand(iattrorderless,fun%def%sym%attr) /= 0
      elseif(ktfsymbolq(k,sym))then
        tfordlessq=iand(iattrorderless,sym%attr) /= 0
      else
        tfordlessq=.false.
      endif
      return
      end

      subroutine itfpmatdummy
      use mackw
      return
      end

      integer*4 function itflistmat(k,listp)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kph,ki,kc,kh,kf,k1,kp1,tfefunrefu
      type (sad_dlist) ,intent(inout):: listp
      type (sad_dlist), pointer :: lista,kl1
      integer*8 kpp,ierrorth0
      integer*4 i,isp3,irtc,isp1,isp0,iop,isp2,ic,icm,np
      logical*4 ra
      itflistmat=-1
      kph=listp%head
      main: do
        if(iand(lsimplepat,listp%attr) /= 0)then
          if(ktfnonlistq(k,lista))then
            return
          endif
          isp0=isp
          kh=lista%head
          itflistmat=itfpmat(kh,listp%head)
          if(itflistmat .lt. 0)then
            return
          endif
          isp1=isp
          np=listp%nl
          if(tfordlessq(kh) .and. lista%nl > 1)then
            iop=iordless
            do while(iop /= 0)
              if(ktastk(iop) /= ksad_loc(listp%head%k))then
                iop=itastk2(1,iop)
              else
                exit
              endif
            enddo
            if(iop /= 0)then
              iordless=iop
              itflistmat=min(itflistmat,
     $             itfseqmatstk(itastk(1,iop+1),itastk(2,iop+1),
     $             listp%dbody(1),np,1,ktfreallistq(listp),int8(-1)))
              if(itflistmat .lt. 0)then
                call tfresetpat(kph)
              endif
              return
            else
              kpp=sad_loc(listp%head)
            endif
            icm=isp1
            call tfgetllstkall(lista)
            if(np == 0)then
              if(isp1 .lt. isp)then
                itflistmat=-1
                call tfresetpat(kph)
                isp=isp0
              else
                isp=isp1
              endif
              return
            endif
          else
            call tfgetllstkall(lista)
            if(np == 0)then
              if(isp1 .lt. isp)then
                exit main
              else
                isp=isp1
              endif
              return
            endif
            if(ktfreallistq(listp))then
              if(np /= isp-isp1)then
                exit main
              endif
              do i=1,np
                if(ktfnonrealq(ktastk(isp1+i)) .or.
     $               ktastk(isp1+i) /= listp%dbody(i)%k)then
                  exit main
                endif
              enddo
              isp=isp1
              return
            elseif(tfconstpatternlistbodyqo(listp))then
              if(np /= isp-isp1)then
                exit main
              endif
              do i=1,np
                if(.not. tfsameq(ktastk(isp1+i),listp%dbody(i)%k))then
                  exit main
                endif
              enddo
              return
            endif
            icm=isp1
            ic=itfconstpart(listp,np,0,kc)
            loop20: do while(ic /= 0)
              LOOP_I: do i=icm+1,isp
                if(tfsameq(dtastk(i),kc))then
                  if(i == icm+1 .and. i-isp1 == ic)then
                    icm=i
                  endif
                  exit loop20
                endif
              enddo LOOP_I
              exit main
            enddo loop20
            ic=itfconstpart(listp,np,ic,kc)
          endif
          kpp=0
          isp2=isp
          if(itflistmat .ge. 0)then
            if(np == icm-isp1+1)then
              itflistmat=min(itflistmat,
     $             itfseqmatstk1(icm+1,isp2,listp%dbody(np)))
            else
              itflistmat=min(0,
     $             itfseqmatstk(icm+1,isp2,listp%dbody(1),np,icm-isp1+1,
     $             ktfreallistq(listp),kpp))
            endif
          endif
          if(itflistmat .lt. 0)then
            call tfresetpat(kph)
            isp=isp0
          endif
        else
          select case (kph%k)
          case (ktfoper+mtfalt)
            if(ktfreallistq(listp))then
              if(ktfrealq(k))then
                do i=1,listp%nl
                  if(k%k == listp%dbody(i)%k)then
                    itflistmat=1
                    return
                  endif
                enddo
              endif
            else
              ra=.true.
              do i=1,listp%nl
                ki=listp%dbody(i)
                if(ktfpatq(ki) .or. ktflistq(ki))then
                  ra=.false.
                endif
                itflistmat=itfpmat(k,ki)
                if(itflistmat .ge. 0)then
                  if(ra)then
                    if(i > 1)then
                      listp%dbody(2:i)=listp%dbody(1:i-1)
                      listp%dbody(1)=ki
                    endif
                  endif
                  return
                endif
              enddo
            endif
          case (ktfoper+mtfpattest)
            if(listp%nl == 2)then
              itflistmat=itfpmat(k,listp%dbody(1))
              if(itflistmat .ge. 0)then
                isp3=isp
                isp=isp+1
                dtastk(isp)=tfeevalref(listp%dbody(2),irtc)
                if(irtc /= 0)then
                  if(irtc > 0 .and. ierrorprint /= 0)then
                    call tfreseterror
                  endif
                  itflistmat=-1
                  call tfresetpat(listp%dbody(1))
                  isp=isp3
                else
                  isp=isp+1
                  dtastk(isp)=k
                  ierrorth0=ierrorth
                  ierrorth=10
                  kf=tfefunrefu(isp3+1,irtc)
                  ierrorth=ierrorth0
                  isp=isp3
                  do
                    if(irtc /= 0)then
                      if(irtc > 0 .and. ierrorprint /= 0)then
                        call tfreseterror
                      endif
                      exit
                    elseif(ktfnonrealq(kf) .or. kf%k == 0)then
                      exit
                    endif
                    return
                  enddo
                  itflistmat=-1
                  call tfresetpat(listp%dbody(1))
                  return
                endif
              endif
            endif
          case (ktfoper+mtfrepeated,ktfoper+mtfrepeatednull)
            if(ktfsequenceq(k%k))then
              isp1=isp
              isp=isp+1
              dtastk(isp)=lista%dbody(1)
              isp2=isp
              itflistmat=itfseqmatseq(isp1+1,isp2,listp,1)
              if(itflistmat .lt. 0)then
                isp=isp1
              endif
            else
              if(k%k == ktfoper+mtfnull)then
                if(kph%k == ktfoper+mtfrepeatednull)then
                  itflistmat=1
                endif
              else
                kp1=listp%dbody(1)
                do
                  k1=k
                  itflistmat=itfpmat(k1,kp1)
                  if(itflistmat .lt. 0 .and. ktflistq(k1,kl1))then
                    k1=kl1%head
                    cycle
                  endif
                  exit
                enddo
              endif
            endif
          case default
            if(kph%k == kxliteral)then
              if(listp%nl == 0)then
                if(k%k == ktfoper+mtfnull)then
                  itflistmat=1
                endif
              else
                kp1=listp%dbody(1)
                itflistmat=itfpmat(k,kp1)
              endif
              return
            endif
            listp%attr=ior(lsimplepat,listp%attr)
            cycle
          end select
        endif
        return
      enddo main
      itflistmat=-1
      call tfresetpat(kph)
      isp=isp0
      return
      end

      recursive integer*4 function
     $     itfseqm(isp1,isp2,kp0,ispf,isps,ispt) result(ix)
      use tfstk
      use tfcode
      use eeval
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(in):: kp0
      type (sad_descriptor) kp1,kp,kd,k2,kf,tfefunrefu
      type (sad_pat), pointer :: pat,pat0
      type (sad_dlist), pointer :: listp
      integer*8 ka2,kad,ierrorth0
      integer*4 ,intent(out):: ispf,isps
      integer*4 ,intent(in):: isp1,isp2,ispt
      integer*4 m,iss,isp0,irtc,ireppat,np,i
      ix=-1
      isps=0
      kp=kp0
      if(ktfpatq(kp,pat0))then
        pat=>pat0
c        call tfdebugprint(sad_descr(pat),'itfseqm-0',1)
        do while(associated(pat%equiv))
c          write(*,*)': ',ktfaddr(sad_descr(pat%equiv)),
c     $         ktftyped(sad_descr(pat%equiv))-ktfpat
c          call tfdebugprint(sad_descr(pat%equiv),'itfseqm-e',1)
          pat=>pat%equiv
        enddo
c        call tfdebugprint(sad_descr(pat),'itfseqm-1',1)
c        call tflinkedpat(pat0,pat)
        k2=pat%value
c        call tfdebugprint(sad_descr(pat),'seqm-pat',1)
c        call tfdebugprint(k2,':=',1)
        if(k2%k /= ktfref)then
          ispf=isp1-2
          if(ktfrefq(k2,ka2) .and. ka2 > 3)then
            iss=int(ka2-ispbase)
            np=itastk2(1,iss)-iss
            m=isp1+np-1
            if(m .le. ispt)then
              do i=1,np
                if(.not. tfsameq(dtastk(iss+i),dtastk(isp1+i-1)))then
                  return
                endif
              enddo
              ix=0
              isps=m+1
            endif
          elseif(isp1 .le. isp2)then
            if(tfsameq(k2,dtastk(isp1)))then
              ix=1
              isps=isp1+1
            endif
          endif
          return
        endif
        if(pat%default%k /= ktfref .and. ispt == isp1-1)then
          pat%value=pat%default
          ix=0
          isps=isp1
          ispf=isp1-2
          return
        endif
        kd=pat%expr
        if(.not. ktfrefq(kd,kad) .or. kad > 3)then
          ix=itfseqm(isp1,isp2,kd,ispf,isps,ispt)
          if(ix .ge. 0)then
            call tfstkpat(isp1-1,isps-1,pat)
          endif
        elseif(kad == 1)then
          ispf=min(ispt,isp1)-1
          if(ispt .ge. isp1 .and. isp2 .ge. isp1)then
            ix=itfsinglepat(dtastk(isp1),pat)
            isps=isp1+1
          endif
        else
          ispf=ispt-1
          isps=ispt+1
          if(ispt == isp1)then
            ix=itfsinglepat(dtastk(isp1),pat)
          elseif(ispt .ge. isp1-int(kad)+2)then
            if(pat%head%k /= ktfref)then
              itastk2(1,isp1-1)=ispt
              ix=itfsinglepat(sad_descr(ktfref+isp1-1+ispbase),pat)
              if(ix .lt. 0)then
                return
              endif
            else
              ix=0
            endif
            call tfstkpat(isp1-1,ispt,pat)
          endif
        endif
        return
      endif
      if(ktfnonlistq(kp,listp))then
        if(isp1 .le. isp2)then
          isps=isp1+1
          if(tfsameq(dtastk(isp1),kp))then
            ix=1
          endif
c          call tfdebugprint(ktastk(isp1),'seqm',1)
c          call tfdebugprint(kp,'=?=',1)
c          write(*,*)'==> ',ix
        endif
        ispf=min(ispt,isp1)-1
        return
      else
        if(iand(lsimplepatlist,listp%attr) == 0)then
          if(listp%head%k == ktfoper+mtfrepeated .or.
     $         listp%head%k == ktfoper+mtfrepeatednull)then
            kp1=listp%dbody(1)
            ireppat=merge(-3,-2,
     $             listp%head%k == ktfoper+mtfrepeatednull)
            if(ispt .ge. isp1+ireppat+2)then
              ix=1
              ispf=ispt-1
              do i=isp1,ispt
                ix=min(ix,itfpmat(dtastk(i),kp1))
                if(ix .lt. 0)then
                  return
                endif
              enddo
              isps=ispt+1
            else
              ispf=isp1-2
            endif
            return
          elseif(listp%head%k == ktfoper+mtfpattest)then
            kp1=listp%dbody(1)
            ix=itfseqm(isp1,isp2,kp1,ispf,isps,ispt)
            if(ix .lt. 0)then
              return
            endif
            isp0=isp
            call tfevalstk(listp%dbody(2),.true.,irtc)
            do
              if(irtc /= 0)then
                if(irtc > 0 .and. ierrorprint /= 0)then
                  call tfreseterror
                endif
                isp=isp0
                exit
              endif
              ktastk(isp+1:isp+isps-isp1)=ktastk(isp1:isps-1)
              isp=isp+isps-isp1
              ierrorth0=ierrorth
              ierrorth=10
              kf=tfefunrefu(isp0+1,irtc)
              ierrorth=ierrorth0
              isp=isp0
              if(irtc /= 0)then
                if(irtc > 0 .and. ierrorprint /= 0)then
                  call tfreseterror
                endif
                exit
              elseif(ktfnonrealq(kf) .or. kf%k == 0)then
                exit
              endif
              return
            enddo
            ix=-1
            ispf=ispt-1
            call tfresetpat(kp1)
            return
          else
            listp%attr=ior(lsimplepatlist,listp%attr)
          endif
        else
          listp%attr=ior(lsimplepatlist,listp%attr)
        endif
        if(isp1 .le. isp2)then
          isps=isp1+1
          if(listp%head%k == ktfoper+mtfnull)then
            m=isp1+listp%nl-1
            if(m .le. ispt)then
              ix=itfseqmatseq(isp1,m,listp,1)
              isps=m+1
            endif
          else
c            ix=merge(-1,itflistmat(dtastk(isp1),listp),
c     $           ktfnonlistq(ktastk(isp1)) .and.
c     $           iand(lsimplepat,listp%attr) /= 0)
            if(ktfnonlistq(ktastk(isp1)) .and. iand(lsimplepat,listp%attr) /= 0)then
              ix=-1
            else
              ix=itflistmat(dtastk(isp1),listp)
            endif
          endif
        elseif(listp%head%k == ktfoper+mtfnull .and.
     $         listp%nl == 0)then
          ix=0
          isps=isp1
          return
        endif
        ispf=isp1-2
        return
      endif
      end

      integer*4 function itfsinglepat(k,pat)
      use tfstk
      use tfcode
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kx,ke,k1,kv,kh
      type (sad_pat) ,intent(inout):: pat
      type (sad_pat), pointer :: pata
      integer*8 kax,ka,ka1
      integer*4 isp3,irtc,isp1,i
      itfsinglepat=-1
      kx=pat%expr
      if(ktfrefq(kx,kax))then
        if(pat%head%k /= ktfref)then
          isp3=isp
          kv=tfeevalref(pat%head,irtc)
          isp=isp3
          if(irtc /= 0)then
            if(irtc > 0 .and. ierrorprint /= 0)then
              call tfreseterror
            endif
            return
          endif
          if(ktfrefq(k,ka) .and. kax > 1)then
            ka=ka-ispbase
            isp1=itastk2(1,ka)
            do i=int(ka)+1,isp1
              call tfhead(dtastk(i),kh)
              if(.not. tfsameq(kh,kv))then
                return
              endif
            enddo
          else
            call tfhead(k,kh)
            if(.not. tfsameq(kh,kv))then
              return
            endif
          endif
        endif
        if(ktfpatq(k,pata))then
          k1=pata%expr
          if(ktfrefq(k1,ka1) .and. ka1 .le. 3 .and.
     $         ka1 > kax)then
            return
          endif
        endif
      elseif(ktfpatq(k))then
        itfsinglepat=itfpmat(k,kx)
        if(itfsinglepat .lt. 0)then
          return
        endif
      else
        isp3=isp
        ke=tfeevalref(kx,irtc)
        isp=isp3
        if(irtc /= 0)then
          if(irtc > 0 .and. ierrorprint /= 0)then
            call tfreseterror
          endif
          return
        endif
        itfsinglepat=itfpmat(k,ke)
        if(itfsinglepat .lt. 0)then
          return
        endif
      endif
      itfsinglepat=0
      if(pat%sym%loc /= 0)then
        pat%mat=1
        pat%value=k
c        call tfdebugprint(sad_descr(pat),'singlepat',1)
c        call tfdebugprint(k,':=',1)
c        write(*,*)'at ',sad_loc(pat%value)
      endif
      return
      end

      end module pmat
