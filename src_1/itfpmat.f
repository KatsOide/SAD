      module tfpmat

      contains
      subroutine tflinkedpat(pat0,pat)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_pat), target :: pat0
      type (sad_pat), pointer, intent(out):: pat
      logical*4 rep
      rep=.false.
      pat=>pat0
      do while(associated(pat%equiv))
        pat=>pat%equiv
        rep=.true.
      enddo
      if(rep .and. pat%mat .eq. 0)then
        pat%equiv=>pat0
        pat=>pat0
        pat%mat=0
        pat%value%k=ktfref
        nullify(pat%equiv)
      endif
      return
      end subroutine

      end module

      integer*4 function itfpmat(k,kp)
      use tfstk
      implicit none
      type (sad_descriptor) k,kp
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: list
      integer*4 itfpatmat,itflistmat
      logical*4 tfsameqd
      itfpmat=-1
      if(ktfpatqd(kp,pat))then
        itfpmat=itfpatmat(k,pat)
      elseif(ktflistqd(kp,list))then
        if(ktfnonlistqd(k) .and. iand(list%attr,lsimplepat) .ne. 0)then
          return
        endif
        itfpmat=itflistmat(k,list)
      elseif(tfsameqd(k,kp))then
        itfpmat=1
      endif
      return
      end

      integer*4 function itflistmat(k,listp)
      use tfstk
      implicit none
      type (sad_descriptor) kph,k,ki,kc,kh,kf,k1,kp1
      type (sad_list) listp
      type (sad_list), pointer :: lista,kl1
      integer*8 kpp,ierrorth0,isp0,isp1,isp3,isp2,iop,icm,i
      integer*4 itfpmat,irtc,itfseqmatseq,itfseqmatstk,ic,
     $     itfconstpart,np,itfseqmatstk1
      logical*4 tfsameqk,tfconstpatternlistbodyqo,ra,tfordlessq,
     $     tfsameqd
      itflistmat=-1
      kph=listp%dbody(0)
 1    if(iand(lsimplepat,listp%attr) .ne. 0)then
        if(ktfnonlistqd(k,lista))then
          return
        endif
        isp0=isp
        kh=lista%dbody(0)
        itflistmat=itfpmat(kh,listp%head)
        if(itflistmat .lt. 0)then
          return
        endif
        isp1=isp
        np=listp%nl
        if(tfordlessq(kh) .and. lista%nl .gt. 1)then
          iop=iordless
          do while(iop .ne. 0)
            if(klist(iop) .ne. ksad_loc(listp%head))then
              iop=ktastk2(iop)
            else
              exit
            endif
          enddo
          if(iop .ne. 0)then
            iordless=iop
            itflistmat=min(itflistmat,
     $           itfseqmatstk(klist(iop+1),klist(iop+2),
     $           listp%dbody(1),np,1,ktfreallistqo(listp),int8(-1)))
            if(itflistmat .lt. 0)then
              call tfresetpat(kph)
            endif
            return
          else
            kpp=ksad_loc(listp%head)
          endif
          icm=isp1
          call tfgetllstkall(lista)
          if(np .eq. 0)then
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
          if(np .eq. 0)then
            if(isp1 .lt. isp)then
              go to 9000
            else
              isp=isp1
            endif
            return
          endif
          if(ktfreallistqo(listp))then
            if(np .ne. isp-isp1)then
              go to 9000
            endif
            do i=1,np
              if(ktfnonrealq(ktastk(isp1+i)))then
                go to 9000
              elseif(ktastk(isp1+i) .ne. listp%body(i))then
                go to 9000
              endif
            enddo
            isp=isp1
            return
          elseif(tfconstpatternlistbodyqo(listp))then
            if(np .ne. isp-isp1)then
              go to 9000
            endif
            do i=1,np
              if(.not. tfsameqk(ktastk(isp1+i),listp%body(i)))then
                go to 9000
              endif
            enddo
            return
          endif
          icm=isp1
          ic=itfconstpart(listp,np,0,kc)
          do while(ic .ne. 0)
            LOOP_I: do i=icm+1,isp
              if(tfsameqd(dtastk(i),kc))then
                go to 20
              endif
            enddo LOOP_I
            go to 9000
 20         if(i .eq. icm+1 .and. i-isp1 .eq. ic)then
              icm=i
            endif
            ic=itfconstpart(listp,np,ic,kc)
          enddo
          kpp=0
        endif
        isp2=isp
        if(np .eq. icm-isp1+1)then
          itflistmat=min(itflistmat,
     $         itfseqmatstk1(icm+1,isp2,listp%dbody(np)))
        else
          itflistmat=min(itflistmat,
     $         itfseqmatstk(icm+1,isp2,listp%dbody(1),
     $         np,int(icm-isp1+1),ktfreallistqo(listp),kpp),0)
        endif
        if(itflistmat .lt. 0)then
          call tfresetpat(kph)
          isp=isp0
        endif
        return
      else
        select case (kph%k)
        case (ktfoper+mtfalt)
          if(ktfreallistqo(listp))then
            if(ktfrealqd(k))then
              do i=1,listp%nl
                if(k%k .eq. listp%body(i))then
                  itflistmat=1
                  return
                endif
              enddo
            endif
          else
            ra=.true.
            do i=1,listp%nl
              ki=listp%dbody(i)
              if(ktfpatqd(ki) .or. ktflistqd(ki))then
                ra=.false.
              endif
              itflistmat=itfpmat(k,ki)
              if(itflistmat .ge. 0)then
                if(ra)then
                  if(i .gt. 1)then
                    listp%dbody(2:i)=listp%dbody(1:i-1)
                    listp%dbody(1)=ki
                  endif
                endif
                return
              endif
            enddo
          endif
        case (ktfoper+mtfpattest)
          if(listp%nl .eq. 2)then
            itflistmat=itfpmat(k,listp%body(1))
            if(itflistmat .ge. 0)then
              isp3=isp
              isp=isp+1
              call tfeevalref(listp%body(2),ktastk(isp),irtc)
              if(irtc .ne. 0)then
                if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
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
                call tfefunref(isp3+1,kf,.true.,irtc)
                ierrorth=ierrorth0
                isp=isp3
                if(irtc .ne. 0)then
                  if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                    call tfreseterror
                  endif
                  go to 2010
                elseif(ktfnonrealqd(kf))then
                  go to 2010
                elseif(kf%k .eq. 0)then
                  go to 2010
                endif
                return
 2010           itflistmat=-1
                call tfresetpat(listp%dbody(1))
                return
              endif
            endif
          endif
        case (ktfoper+mtfrepeated,ktfoper+mtfrepeatednull)
          if(ktfsequenceq(k%k))then
            isp1=isp
            isp=isp+1
            ktastk(isp)=lista%body(1)
            isp2=isp
            itflistmat=itfseqmatseq(isp1+1,isp2,listp,1)
            if(itflistmat .lt. 0)then
              isp=isp1
            endif
          else
            if(k%k .eq. ktfoper+mtfnull)then
              if(kph%k .eq. ktfoper+mtfrepeatednull)then
                itflistmat=1
              endif
            else
              kp1=listp%dbody(1)
 101          k1=k
              itflistmat=itfpmat(k1,kp1)
              if(itflistmat .lt. 0 .and. ktflistqd(k1,kl1))then
                k1=kl1%dbody(0)
                go to 101
              endif
            endif
          endif
        case default
          if(kph%k .eq. kxliteral)then
            if(listp%nl .eq. 0)then
              if(k%k .eq. ktfoper+mtfnull)then
                itflistmat=1
              endif
            else
              kp1=listp%dbody(1)
              itflistmat=itfpmat(k,kp1)
            endif
            return
          endif
          listp%attr=ior(lsimplepat,listp%attr)
          go to 1
        end select
        return
      endif
 9000 itflistmat=-1
      call tfresetpat(kph)
      isp=isp0
      return
      end

      integer*4 function itfconstpart(list,m,n,kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ki
      type (sad_list) list
      integer*4 n,m,i,is,n1
      logical*4 tfconstpatternqk
      itfconstpart=0
      n1=n+1
      is=list%attr
      if(iand(kallnofixedarg,is) .ne. 0)then
        return
      else
        do i=n1,m
          ki=list%dbody(i)
          if(ktfpatqd(ki))then
            list%attr=ior(knofixedarg,list%attr)
          elseif(ktflistqd(ki))then
            if(tfconstpatternqk(ki))then
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
        if(n1 .eq. 1)then
          list%attr=ior(kallnofixedarg,list%attr)
        endif
      endif
      return
      end

      integer*4 function itfpatmat(k,pat0)
      use tfstk
      use tfcode
      use tfpmat
      implicit none
      type (sad_descriptor) k,kv1
      type (sad_pat) pat0
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: list
      integer*8 kav1,isp1,isp2
      integer*4 itfsinglepat,i,np,iss
      logical*4 tfsameqd
      itfpatmat=-1
      call tflinkedpat(pat0,pat)
      kv1=pat%value
      if(pat%mat .eq. 0)then
        itfpatmat=itfsinglepat(k,pat)
      elseif(ktfrefqd(kv1,kav1) .and. kav1 .gt. 3)then
        if(ktflistqd(k,list))then
          if(list%head .eq. ktfoper+mtfnull)then
            isp1=isp
            call tfgetllstkall(list)
            iss=int(kav1)
            isp2=isp
            np=ktastk2(iss)-iss
            if(isp2-isp1 .ge. np)then
              do i=1,np
                if(.not. tfsameqd(dtastk(isp1+i),dtastk(iss+i)))then
                  isp=isp1
                  return
                endif
              enddo
              itfpatmat=1
            endif
          endif
        endif
      elseif(tfsameqd(k,pat%value))then
        itfpatmat=1
      endif
      return
      end

      integer*4 function itfseqmatseq(isp1,isp2,kl,mp1)
      use tfstk
      implicit none
      type (sad_list) kl
      integer*8 isp1,isp2
      integer*4 mp1,itfseqmatstk
      if(isp1 .gt. isp2 .and. mp1 .gt. kl%nl)then
        itfseqmatseq=0
      else
        itfseqmatseq=itfseqmatstk(isp1,isp2,kl%body(1),
     $       kl%nl,mp1,ktfreallistqo(kl),int8(0))
      endif
      return
      end

      integer*4 function itfseqmat(isp10,isp20,kp,mp,mp1)
      use tfstk
      implicit none
      type (sad_descriptor) kp(1)
      integer*8 isp1,isp2,isp10,isp20
      integer*4 mp,mp1,itfseqmatstk,itfseqmatstk1
      isp1=isp10
      isp2=isp20
      if(isp1 .gt. isp2 .and. mp1 .gt. mp)then
        itfseqmat=0
      elseif(mp .eq. 1 .and. mp1 .eq. 1)then
        itfseqmat=itfseqmatstk1(isp1,isp2,kp(1))
      else
        itfseqmat=itfseqmatstk(isp1,isp2,kp,mp,mp1,.false.,int8(0))
      endif
      return
      end

      integer*4 function itfseqmatstk1(isp10,isp20,kp)
      use tfstk
      implicit none
      type (sad_descriptor) kp
      integer*4 itfseqm,i
      integer*8 isp10,isp20,iop0,ispm0,mstk0,ispm,isps
      mstk0=mstk
      ispm=isp10-1
      iop0=iordless
      ispm0=ispm
 1020 i=itfseqm(isp10,isp20,kp,ispm,isps,isp20)
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
        if(iordless .ne. iop0)then
          ispm=ispm0
          go to 1020
        endif
      else
        iordless=iop0
      endif
      mstk=mstk0
      itfseqmatstk1=-1
      return
      end

      recursive integer*4 function
     $     itfseqmatstk(isp10,isp20,kp0,map,mp1,realp,kpp)
     $      result(ix)
      use tfstk
      implicit none
      type (sad_descriptor) kp,kp0(1),k1
      integer*8 kpp,isps,ispf,isp10,isp20,isp1a,isp2a,ispt,mstk0,
     $     iop,iop0,ispp
      integer*4 mp1,itfseqm,i,map,np,mop,np0,npd,npdi,itfseqmatstk1
      logical*4 realp
      mop=0
      if(mp1 .gt. map)then
        if(isp20 .lt. isp10)then
          ix=1
        else
          ix=-1
        endif
        return
      endif
      if(realp)then
        kp%k=0
      else
        kp=kp0(mp1)
      endif
      mstk0=mstk
      if(kpp .eq. 0)then
        isp1a=isp10
        isp2a=isp20
        ispf=isp20
        if(map .eq. mp1+1)then
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
 1011       i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
            if(i .ge. 0)then
              ix=min(i,itfseqmatstk1(isps,isp2a,kp0(map)))
              if(ix .ge. 0)then
                return
              endif
              call tfresetpat(kp)
              if(iordless .ne. iop0)then
                ispf=ispt
                go to 1011
              endif
            else
              iordless=iop0
            endif
          enddo
        elseif(map .ne. mp1)then
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
 1010       i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
            if(i .ge. 0)then
              ix=min(i,itfseqmatstk(isps,isp2a,
     $             kp0,map,mp1+1,realp,int8(0)))
              if(ix .ge. 0)then
                return
              endif
              call tfresetpat(kp)
              if(iordless .ne. iop0)then
                ispf=ispt
                go to 1010
              endif
            else
              iordless=iop0
            endif
          enddo
        else
          do while(ispf .ge. isp1a-1)
            iop0=iordless
            ispt=ispf
 1020       i=itfseqm(isp1a,isp2a,kp,ispf,isps,isp2a)
            if(i .ge. 0)then
              if(isp2a .lt. isps)then
                ix=min(1,i)
              else
                ix=-1
              endif
              if(ix .ge. 0)then
                return
              endif
              call tfresetpat(kp)
              if(iordless .ne. iop0)then
                ispf=ispt
                go to 1020
              endif
            else
              iordless=iop0
            endif
          enddo
        endif
        mstk=mstk0
        ix=-1
        return
      endif
      iop=iordless
 3    if(iop .ne. 0)then
        if(kpp .gt. 0)then
          do while(iop .ne. 0 .and. klist(iop) .ne. kpp)
            iop=ktastk2(iop)
          enddo
          if(iop .eq. 0)then
            go to 3
          endif
          iordless=iop
        endif
        ispf=ktastk2(iop)
        isp1a=klist(iop+1)
        isp2a=klist(iop+2)
        np=itastk2(1,iop+2)
      else
        mstk=mstk-3
        iop=mstk+1
        klist(iop)=kpp
        klist(iop+1)=isp10
        klist(iop+2)=isp20
        ktastk2(iop)=iordless
        ktastk2(iop+1)=ispf
        itastk2(1,iop+2)=0
        np=0
        iordless=iop
        isp1a=isp10
        isp2a=isp20
        ispf=isp20
      endif
      do while(.true.)
        do while(ispf .ge. isp1a-1)
          iop0=iordless
          ispt=ispf
          ix=-1
          do while (ix .lt. 0)
            if(map .eq. mp1)then
              ispt=isp2a
            else
              ispt=ispf
            endif
            i=itfseqm(isp1a,isp2a,kp,ispf,isps,ispt)
            if(i .ge. 0)then
              if(mp1 .ge. map)then
                if(isp2a .lt. isps)then
                  ix=min(1,i)
                endif
              else
                if(map .eq. mp1+1)then
                  ix=min(i,itfseqmatstk1(isps,isp2a,kp0(map)))
                else
                  ix=min(i,
     $                 itfseqmatstk(isps,isp2a,kp0,map,mp1+1,
     $                 realp,int8(0)))
                endif
              endif
              if(ix .ge. 0)then
                ktastk2(iop+1)=ispf
                return
              endif
              call tfresetpat(kp)
              if(iordless .ne. iop0)then
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
        np0=int(isp2a-isp1a+1)
        do i=2,np0+1
          npdi=npd/i
          if(npdi*i .ne. npd)then
            mop=np0-i
            exit
          endif
          npd=npdi
        enddo
        if(mop .lt. 0)then
          mstk=mstk0
          ix=-1
          iordless=ktastk2(iordless)
          return
        endif
        itastk2(1,iop+2)=np
        ispp=isp1a+mop
        k1=dtastk(ispp)
        dtastk(ispp)=dtastk(ispp+1)
        dtastk(ispp+1)=k1
c     call tfdebugprint(k1,'seqmatstk',1)
c     call tfdebugprint(ktastk(ispp),'<->',1)
c     write(*,*)'at ',ispp,' with ',mop,np
        ispf=isp20
      enddo
      end

      recursive integer*4 function
     $     itfseqm(isp1,isp2,kp0,ispf,isps,ispt) result(ix)
      use tfstk
      use tfcode
      use tfpmat
      use iso_c_binding
      implicit none
      type (sad_descriptor) kp1,kp0,kp,kd,k2,kf
      type (sad_pat), pointer :: pat,pat0
      type (sad_list), pointer :: listp
      integer*8 ka2,kad,ierrorth0,isp1,isp2,ispf,isp0,ispt,isps,m,i,
     $     iss
      integer*4 itfpmat,irtc,itfseqmatseq,itflistmat,
     $     itfsinglepat,ireppat,np
      logical*4 tfsameqd
      ix=-1
      kp=kp0
      if(ktfpatqd(kp,pat0))then
        call tflinkedpat(pat0,pat)
        k2=pat%value
        if(k2%k .ne. ktfref)then
          ispf=isp1-2
          if(ktfrefqd(k2,ka2) .and. ka2 .gt. 3)then
            iss=ka2
            np=int(ktastk2(iss)-iss)
            m=isp1+np-1
            if(m .le. ispt)then
              do i=1,np
                if(.not. tfsameqd(dtastk(iss+i),dtastk(isp1+i-1)))then
                  return
                endif
              enddo
              ix=0
              isps=m+1
            endif
          elseif(isp1 .le. isp2)then
            if(tfsameqd(k2,dtastk(isp1)))then
              ix=1
              isps=isp1+1
            endif
          endif
          return
        endif
        if(pat%default%k .ne. ktfref .and. ispt .eq. isp1-1)then
          pat%value=pat%default
          ix=0
          isps=isp1
          ispf=isp1-2
          return
        endif
        kd=pat%expr
        if(.not. ktfrefqd(kd,kad) .or. kad .gt. 3)then
          ix=itfseqm(isp1,isp2,kd,ispf,isps,ispt)
          if(ix .ge. 0)then
            call tfstkstk(isp1-1,isps-1,pat)
          endif
        elseif(kad .eq. 1)then
          ispf=min(ispt,isp1)-1
          if(ispt .ge. isp1 .and. isp2 .ge. isp1)then
            ix=itfsinglepat(dtastk(isp1),pat)
            isps=isp1+1
          endif
        else
          ispf=ispt-1
          isps=ispt+1
          if(ispt .eq. isp1)then
            ix=itfsinglepat(dtastk(isp1),pat)
          elseif(ispt .ge. isp1-kad+2)then
            if(pat%head%k .ne. ktfref)then
              ktastk2(isp1-1)=ispt
              ix=itfsinglepat(ktfref+isp1-1,pat)
              if(ix .lt. 0)then
                return
              endif
            else
              ix=0
            endif
            call tfstkstk(isp1-1,ispt,pat)
          endif
        endif
        return
      endif
      if(ktfnonlistqd(kp,listp))then
        if(isp1 .le. isp2)then
          isps=isp1+1
          if(tfsameqd(dtastk(isp1),kp))then
            ix=1
          endif
c          call tfdebugprint(ktastk(isp1),'seqm',1)
c          call tfdebugprint(kp,'=?=',1)
c          write(*,*)'==> ',ix
        endif
        ispf=min(ispt,isp1)-1
        return
      else
        if(iand(lsimplepatlist,listp%attr) .eq. 0)then
          if(listp%head .eq. ktfoper+mtfrepeated .or.
     $         listp%head .eq. ktfoper+mtfrepeatednull)then
            kp1=listp%dbody(1)
            if(listp%head .eq. ktfoper+mtfrepeatednull)then
              ireppat=-3
            else
              ireppat=-2
            endif
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
          elseif(listp%head .eq. ktfoper+mtfpattest)then
            kp1=listp%dbody(1)
            ix=itfseqm(isp1,isp2,kp1,ispf,isps,ispt)
            if(ix .lt. 0)then
              return
            endif
            isp0=isp
            call tfevalstk(listp%body(2),.true.,irtc)
            if(irtc .ne. 0)then
              if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                call tfreseterror
              endif
              isp=isp0
              go to 301
            endif
            ktastk(isp+1:isp+isps-isp1)=ktastk(isp1:isps-1)
            isp=isp+isps-isp1
            ierrorth0=ierrorth
            ierrorth=10
            call tfefunref(isp0+1,kf,.true.,irtc)
            ierrorth=ierrorth0
            isp=isp0
            if(irtc .ne. 0)then
              if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                call tfreseterror
              endif
              go to 301
            elseif(ktfnonrealqd(kf))then
              go to 301
            elseif(kf%k .eq. 0)then
              go to 301
            endif
            return
 301        ix=-1
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
          if(listp%head .eq. ktfoper+mtfnull)then
            m=isp1+listp%nl-1
            if(m .le. ispt)then
              ix=itfseqmatseq(isp1,m,listp,1)
              isps=m+1
            endif
          else
            if(ktfnonlistq(ktastk(isp1)) .and.
     $           iand(lsimplepat,listp%attr) .ne. 0)then
              ix=-1
            else
              ix=itflistmat(ktastk(isp1),listp)
            endif
          endif
        elseif(listp%head .eq. ktfoper+mtfnull .and.
     $         listp%nl .eq. 0)then
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
      implicit none
      type (sad_descriptor) k,kx,ke,k1,kv,kh
      type (sad_pat) pat
      type (sad_pat), pointer :: pata
      integer*8 kax,ka,ka1,i,isp1,isp3
      integer*4 irtc,itfpmat
      logical*4 tfsameqd
      itfsinglepat=-1
      kx=pat%expr
      if(ktfrefqd(kx,kax))then
        if(pat%head%k .ne. ktfref)then
          isp3=isp
          call tfeevalref(pat%head,kv,irtc)
          isp=isp3
          if(irtc .ne. 0)then
            if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
              call tfreseterror
            endif
            return
          endif
          if(ktfrefqd(k,ka) .and. kax .gt. 1)then
            isp1=ktastk2(ka)
            do i=ka+1,isp1
              call tfhead(dtastk(i),kh)
              if(.not. tfsameqd(kh,kv))then
                return
              endif
            enddo
          else
            call tfhead(k,kh)
            if(.not. tfsameqd(kh,kv))then
              return
            endif
          endif
        endif
        if(ktfpatqd(k,pata))then
          k1=pata%expr
          if(ktfrefqd(k1,ka1) .and. ka1 .le. 3 .and.
     $         ka1 .gt. kax)then
            return
          endif
        endif
      elseif(ktfpatqd(k))then
        itfsinglepat=itfpmat(k,kx)
        if(itfsinglepat .lt. 0)then
          return
        endif
      else
        isp3=isp
        call tfeevalref(kx,ke,irtc)
        isp=isp3
        if(irtc .ne. 0)then
          if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
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
      if(pat%sym%loc .ne. 0)then
        pat%mat=1
        pat%value=k
      endif
      return
      end

      subroutine tfstkstk(isp1,isp2,pat)
      use tfstk
      use tfcode
      implicit none
      type (sad_pat) pat
      integer*8 isp1,isp2
      if(pat%sym%loc .ne. 0)then
        pat%mat=1
        if(isp1 .ge. isp2)then
          pat%value=dxnull
        elseif(isp2 .eq. isp1+1)then
          pat%value%k=ktastk(isp2)
        else
          pat%value%k=ktfref+isp1
          ktastk2(isp1)=isp2
        endif
      endif
      return
      end

      integer*4 function itfpmatc(k,kp)
      use tfstk
      implicit none
      type (sad_descriptor) kp,k
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: klp
      integer*8 mstk0,iop,isp0
      integer*4 itfpatmat,itfpmatcl
      logical*4 tfsameqd
      if(ktflistqd(kp,klp))then
        itfpmatc=itfpmatcl(k,klp)
      elseif(ktfpatqd(kp,pat))then
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
      elseif(tfsameqd(k,kp))then
        itfpmatc=1
      else
        itfpmatc=-1
      endif
      return
      end

      integer*4 function itfpmatcl(k,klp)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list) klp
      integer*8 mstk0,iop,isp0
      integer*4 itflistmat
      if(ktfnonlistqd(k) .and. iand(klp%attr,lsimplepat) .ne. 0)then
        itfpmatcl=-1
        return
      endif
      iop=iordless
      iordless=0
      mstk0=mstk
      isp0=isp
      if(iand(klp%attr,lnopatlist) .eq. 0)then
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

      recursive subroutine tfunsetpat(k)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) k
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: klp
      if(ktflistqd(k,klp))then
        call tfunsetpatlist(klp)
      elseif(ktfpatqd(k,pat))then
        pat%mat=0
        pat%value%k=ktfref
        call tfunsetpat(pat%expr)
      endif
      return
      end

      subroutine tfunsetpatlist(klp)
      use tfstk
      implicit none
      type (sad_list) klp
      integer*4 i
      if(iand(klp%attr,lnopatlist) .ne. 0)then
        return
      endif
      call tfunsetpat(klp%dbody(0))
      if(ktfnonreallistqo(klp))then
        if(iand(klp%attr,knopatarg) .eq. 0)then
          do i=1,klp%nl
            call tfunsetpat(klp%dbody(i))
          enddo
        endif
      endif
      return
      end

      subroutine tfresetpat(kp)
      use tfstk
      implicit none
      type (sad_descriptor) kp
      type (sad_pat), pointer :: pat
      type (sad_list), pointer :: list
      if(ktflistqd(kp,list))then
        call tfresetpatlist(list)
      elseif(ktfpatqd(kp,pat))then
        call tfresetpatpat(pat)
      endif
      return
      end

      subroutine tfresetpatpat(pat)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_pat) pat
      if(pat%value%k .ne. ktfref)then
        pat%mat=0
        pat%value%k=ktfref
        nullify(pat%equiv)
        call tfresetpat(pat%expr%k)
      endif
      return
      end

      subroutine tfresetpatlist(list)
      use tfstk
      implicit none
      type (sad_list) list
      integer*4 i
      if(iand(list%attr,lnopatlist) .ne. 0)then
        return
      endif
      call tfresetpat(list%dbody(0))
      if(ktfnonreallistqo(list))then
        if(iand(list%attr,knopatarg) .eq. 0)then
          do i=1,list%nl
            call tfresetpat(list%dbody(i))
          enddo
        endif
      endif
      return
      end

      recursive subroutine tfinitpat(isp0,k)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      type (sad_pat), pointer :: pat
      integer*8 kas,i,isp0
      if(ktflistqd(k,kl))then
        call tfinitpatlist(isp0,kl)
      elseif(ktfpatqd(k,pat))then
        if(.not. associated(pat%equiv))then
          if(pat%sym%loc .ne. 0)then
            kas=ktfaddr(pat%sym%alloc)
            do i=isp0+1,isp-1,2
              if(kas .eq. ktastk(i+1))then
                call loc_pat(ktfaddr(ktastk(i)),pat%equiv)
                go to 10
              endif
            enddo
            isp=isp+2
            dtastk(isp-1)=k
            ktastk(isp  )=kas
          endif
          pat%value%k=ktfref
          nullify(pat%equiv)
        endif
 10     call tfinitpat(isp0,pat%expr)
        pat%mat=0
      endif
      return
      end

      subroutine tfinitpatlist(isp0,kl)
      use tfstk
      implicit none
      type (sad_list) kl
      integer*8 isp1,isp0
      integer*4 m,i
      if(iand(kl%attr,lnopatlist) .ne. 0)then
        return
      endif
      isp1=isp
      call tfinitpat(isp0,kl%dbody(0))
      if(ktfnonreallistqo(kl))then
        if(iand(kl%attr,knopatarg) .eq. 0)then
          m=kl%nl
          if(kl%head .eq. ktfoper+mtfpattest)then
            m=1
          endif
          do i=1,m
            call tfinitpat(isp0,kl%dbody(i))
          enddo
          if(isp .eq. isp1)then
            kl%attr=ior(kl%attr,knopatarg)
          endif
        endif
      endif
      if(isp .eq. isp1)then
        kl%attr=ior(kl%attr,lnopatlist)
      endif
      return
      end

      logical*4 function tfconstpatternqk(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      logical*4 tfconstpatternheadqk,tfconstpatternlistbodyqo
      tfconstpatternqk=.true.
      if(ktfpatqd(k))then
        tfconstpatternqk=.false.
      elseif(ktflistqd(k,kl))then
        tfconstpatternqk=tfconstpatternheadqk(kl%dbody(0)) .and.
     $       tfconstpatternlistbodyqo(kl)
      endif
      return
      end

      recursive logical*4 function tfconstpatternheadqk(k)
     $     result(lx)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      integer*8 ka
      logical*4 tfconstpatternlistbodyqo,tfsamesymbolqk
      if(ktfpatqd(k))then
        lx=.false.
      elseif(ktfoperqd(k,ka))then
        lx=ka .ne. mtfalt .and. ka .ne. mtfrepeated .and.
     $       ka .ne. mtfrepeatednull
      elseif(ktflistqd(k,kl))then
        lx=tfconstpatternheadqk(kl%dbody(0)) .and.
     $       tfconstpatternlistbodyqo(kl)
      elseif(ktfsymbolqd(k))then
        lx=.not. tfsamesymbolqk(k%k,kxliteral)
      else
        lx=.true.
      endif
      return
      end

      logical function tfconstpatternlistbodyqo(list)
      use tfstk
      implicit none
      type (sad_descriptor) ki
      type (sad_list) list
      integer*4 i
      logical*4 tfconstpatternqk
      tfconstpatternlistbodyqo=.true.
      if(ktfnonreallistqo(list))then
        if(iand(list%attr,kconstlist) .ne. 0)then
          return
        elseif(iand(list%attr,knoconstlist) .ne. 0)then
          tfconstpatternlistbodyqo=.false.
          return
        endif
        do i=1,list%nl
          ki=list%dbody(i)
          if(ktfpatqd(ki))then
            tfconstpatternlistbodyqo=.false.
          elseif(ktflistqd(ki))then
            tfconstpatternlistbodyqo=tfconstpatternqk(ki)
          endif
          if(.not. tfconstpatternlistbodyqo)then
            list%attr=ior(list%attr,knoconstlist)
            return
          endif
        enddo
        list%attr=ior(list%attr,kconstlist)
      endif
      return
      end

      logical*4 function tfordlessq(k)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) k
      type (sad_funtbl), pointer ::fun
      type (sad_symbol), pointer ::sym
      if(ktfoperqd(k))then
        call c_f_pointer(c_loc(klist(klist(ifunbase+ktfaddrd(k))-9)),
     $       fun)
        tfordlessq=iand(iattrorderless,fun%def%sym%attr) .ne. 0
      elseif(ktfsymbolqd(k,sym))then
        tfordlessq=iand(iattrorderless,sym%attr) .ne. 0
      else
        tfordlessq=.false.
      endif
      return
      end

      subroutine itfpmatdummy
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      return
      end
