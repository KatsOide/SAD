      subroutine tfpart(isp1,kx,err,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer ::kl
      integer*4 isp1,irtc,isp2
      logical*4 err
      if(isp .le. isp1)then
        kx=dtastk(isp1)
        irtc=0
        return
      endif
      isp2=isp
      call loc_sad(ktfaddr(ktastk(isp1)),kl)
      call tfpart1(kl,isp1,isp2,kx,err,irtc)
      return
      end

      recursive subroutine tfpart1(list,isp1,isp2,kx,err,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kl,ki
      type (sad_list) list
      type (sad_list), pointer :: listi,listl
      integer*4 isp1,irtc,narg,ivi,iv,ma,m,i,
     $     isp0,itfmessage,itfmessageexp,isp2
      logical*4 err
      narg=isp2-isp1
      ma=list%nl
      kl=dtastk(isp1+1)
      isp0=isp
      irtc=0
      if(ktfrealqdi(kl,iv))then
        if(iv .lt. 0)then
          iv=ma+iv+1
          if(iv .lt. 0)then
            ivi=iv
            go to 9030
          endif
        elseif(iv .gt. ma)then
          ivi=iv
          go to 9030
        endif
        kx=list%dbody(iv)
        if(narg .ne. 1)then
          if(ktflistqd(kx,listi))then
            call tfpart1(listi,isp1+1,isp2,kx,err,irtc)
          else
            go to 9040
          endif
        endif
        return
      endif
      if(tflistqd(kl,listl))then
        if(ktfnonreallistqo(listl))then
          if(err)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"List of reals as index"')
          else
            irtc=-1
          endif
          return
        endif
        m=listl%nl
        if(narg .gt. 1)then
          do i=1,m
            ivi=int(listl%rbody(i))
            if(ivi .lt. 0)then
              ivi=ma+ivi+1
              if(ivi .lt. 0)then
                go to 9030
              endif
            elseif(ivi .gt. ma)then
              go to 9030
            endif
            ki=list%dbody(ivi)
            if(ktflistqd(ki,listi))then
              isp=isp+1
              call tfpart1(listi,isp1+1,isp2,dtastk(isp),err,irtc)
              if(irtc .ne. 0)then
                isp=isp0
                return
              endif
            else
              go to 9040
            endif
          enddo
        else
          do i=1,m
            ivi=int(listl%rbody(i))
            if(ivi .lt. 0)then
              ivi=ma+ivi+1
              if(ivi .lt. 0)then
                go to 9030
              endif
            elseif(ivi .gt. ma)then
              go to 9030
            endif
            ktastk(isp0+i)=list%body(ivi)
          enddo
          isp=isp0+m
        endif
        kx=kxmakelist(isp0)
        isp=isp0
      elseif(kl%k .eq. ktfoper+mtfnull)then
        if(narg .eq. 1)then
          kx%k=ktflist+ksad_loc(list%head)
        else
          if(ktfnonreallistqo(list))then
            do i=1,ma
              if(ktflistq(list%body(i),listi))then
                isp=isp+1
                call tfpart1(listi,isp1+1,isp2,dtastk(isp),err,irtc)
                if(irtc .ne. 0)then
                  isp=isp0
                  return
                endif
              else
                go to 9040
              endif
            enddo
          endif
          kx=kxmakelist(isp0)
          isp=isp0
        endif
      elseif(ierrorexp .gt. 0)then
        irtc=-1
      elseif(err)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real, List of Reals, or Null as index"')
      else
        irtc=-1
      endif
      return
 9030 if(err)then
        irtc=itfmessageexp(9,'General::index',dble(ivi))
      else
        irtc=-1
      endif
      isp=isp0
      return
 9040 if(err)then
        irtc=itfmessage(9,'General::toomany','"indices"')
      else
        irtc=-1
      endif
      isp=isp0
      return
      end

      recursive subroutine tfpartrstk(lar,isp1,isp2,
     $     list,last,write,eval,err,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ki,kl,k1
      type (sad_list) lar
      type (sad_list),pointer :: lari,listl
      integer*4 isp1,irtc,narg,ivi,iv,ma,m,i,isp0,itfmessage,
     $     itfmessageexp,isp2
      logical*4 err,list,list1,last,write,eval,eval1
      irtc=0
      eval=.false.
      narg=isp2-isp1
      ma=lar%nl
      if(narg .le. 0)then
        return
      endif
      kl=dtastk(isp1+1)
      if(ktfrealqdi(kl,iv))then
        if(iv .lt. 0)then
          iv=ma+iv+1
          if(iv .lt. 0)then
            ivi=iv
            go to 9030
          endif
        elseif(last .and. iv .eq. ma+1)then
        elseif(iv .gt. ma)then
          ivi=iv
          go to 9030
        endif
        if(narg .gt. 1)then
          k1=lar%dbody(iv)
          if(ktflistqd(k1,lari))then
            if(write)then
              call tfclonelist(lari,lari)
              call tfreplist(lar,iv,
     $             ktflist+ksad_loc(lari%head),eval1)
              eval=eval .or. eval1
            endif
            call tfpartrstk(lari,isp1+1,isp2,
     $           list,last,write,eval1,err,irtc)
            eval=eval .or. eval1
          else
            go to 9040
          endif
        elseif(write)then
          isp=isp+1
          ktastk(isp)=ksad_loc(lar%head)
          itastk2(1,isp)=iv
          list=.false.
        else
          isp=isp+1
          ktastk(isp)=lar%body(iv)
        endif
      else
        list=.true.
        isp0=isp
        if(tflistqd(kl,listl))then
          if(ktfnonreallistqo(listl))then
            if(err)then
              irtc=itfmessage(9,'General::wrongtype',
     $             '"List of reals as index"')
            else
              irtc=-1
            endif
            return
          endif
          m=listl%nl
          if(narg .gt. 1)then
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              ki=lar%dbody(ivi)
              if(ktflistqd(ki,lari))then
                if(write)then
                  call tfclonelist(lari,lari)
                  call tfreplist(lar,ivi,
     $                 ktflist+ksad_loc(lari%head),eval1)
                  eval=eval .or. eval1
                endif
                call tfpartrstk(lari,isp1+1,isp2,
     $               list1,last,write,eval1,err,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                eval=eval .or. eval1
              else
                isp=isp0
                go to 9040
              endif
            enddo
          elseif(write)then
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(last .and. ivi .eq. ma+1)then
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              isp=isp+1
              ktastk(isp)=ksad_loc(lar%head)
              itastk2(1,isp)=ivi
            enddo
          else
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(last .and. ivi .eq. ma+1)then
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              isp=isp+1
              ktastk(isp)=lar%body(ivi)
            enddo
          endif
        elseif(kl%k .eq. ktfoper+mtfnull)then
          if(narg .eq. 1)then
            if(write)then
              do i=1,ma
                isp=isp+1
                ktastk(isp)=ksad_loc(lar%head)
                itastk2(1,isp)=i
              enddo
            else
              call tfgetllstkall(lar)
            endif
          else
            do i=1,ma
              ki=lar%dbody(i)
              if(ktflistqd(ki,lari))then
                if(write)then
                  call tfclonelist(lari,lari)
                  call tfreplist(lar,i,
     $                 ktflist+ksad_loc(lari%head),eval1)
                  eval=eval .or. eval1
                endif
                call tfpartrstk(lari,isp1+1,isp2,list1,
     $               last,write,eval1,err,irtc)
                if(irtc .ne. 0)then
                  isp=isp0
                  return
                endif
                eval=eval .or. eval1
              else
                isp=isp0
                go to 9040
              endif
            enddo
          endif
        elseif(err)then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real, List of Reals, or Null as index"')
        else
          irtc=-1
        endif
      endif
      return
 9030 if(err)then
        irtc=itfmessageexp(9,'General::index',dble(ivi))
      else
        irtc=-1
      endif
      return
 9040 if(err)then
        irtc=itfmessage(9,'General::toomany','"indices"')
      else
        irtc=-1
      endif
      return
      end

      subroutine tffirst(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k
      type (sad_list), pointer :: kl
      integer*4 isp1,irtc,itfmessage,mode,m
      if(isp .ne. isp1+1)then
        if(rlist(iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        return
      endif
      k=dtastk(isp)
      if(ktfnonlistqd(k,kl))then
        if(rlist(iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::wrongtype','"List"')
        endif
        return
      endif
      m=kl%nl
      if(m .le. max(0,mode))then
        if(rlist(iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::index','""')
        endif
        return
      endif
      if(mode .lt. 0)then
        kx=kl%dbody(m)
      else
        kx=kl%dbody(mode+1)
      endif
      irtc=0
      call tfeevalref(kx,kx,irtc)
      return
      end

      subroutine tfextract(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kh,k,kind
      type (sad_list), pointer :: klind,kl,kli
      integer*4 irtc,isp1,narg,i,isp0,itfmessage,m
      narg=isp-isp1
      if(narg .eq. 2)then
        kh%k=ktfref
      elseif(narg .eq. 3)then
        kh=dtastk(isp)
      else
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
      endif
      k=dtastk(isp1+1)
      kind=dtastk(isp1+2)
      if(.not. ktflistqd(k,kl) .or. .not. tflistqd(kind,klind))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1, List for #2"')
        return
      endif
      kind=dtfcopy1(kind)
      m=klind%nl
      isp0=isp
      if(ktfnonreallistqo(klind))then
        do i=1,m
          if(ktfnonlistq(klind%body(i)))then
            call tfextract1(kl,klind,kh,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            go to 100
          endif
        enddo
        do i=1,m
          call loc_sad(ktfaddr(klind%body(i)),kli)
          call tfextract1(kl,kli,kh,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
        enddo
 100    kx=kxmakelist(isp0)
      elseif(m .eq. 0)then
        if(kh%k .ne. ktfref)then
          isp=isp+1
          dtastk(isp)=kh
          isp=isp+1
          dtastk(isp)=k
          call tfefunref(isp-1,kx,.true.,irtc)
        else
          call tfleval(kl,kx,.true.,irtc)
        endif
      else
        call tfextract1(kl,klind,kh,irtc)
        kx=dtastk(isp)
      endif
 9000 isp=isp0
      call tflocal1d(kind)
      return
      end

      subroutine tfextract1(kl,kll,kh,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kh
      type (sad_list) kl,kll
      integer*4 irtc,i,isp0,isp1,isp2,isp3
      logical*4 list,eval
      isp1=isp
      call tfgetllstkall(kll)
      isp2=isp
      call tfpartrstk(kl,isp1,isp2,list,
     $     .false.,.false.,eval,.true.,irtc)
      if(irtc .ne. 0)then
        isp=isp1
        return
      endif
      isp3=isp1+isp-isp2
      if(kh%k .eq. ktfref)then
        do i=1,isp-isp2
          call tfeevalref(ktastk(isp2+i),ktastk(isp1+i),irtc)
          if(irtc .ne. 0)then
            isp=isp1
            return
          endif
        enddo
      elseif(kh%k .eq. ktfoper+nfununeval)then
        ktastk(isp1+1:isp1+isp-isp2)=ktastk(isp2+1:isp)
c        do i=1,isp-isp2
c          ktastk(isp1+i)=ktastk(isp2+i)
c        enddo
        irtc=0
      else
        isp0=isp
        do i=1,isp-isp2
          isp=isp0+2
          dtastk(isp-1)=kh
          dtastk(isp)=dtastk(isp2+i)
          call tfefunref(isp-1,ktastk(isp1+i),.true.,irtc)
          if(irtc .ne. 0)then
            isp=isp1
            return
          endif
        enddo
      endif
      isp=isp3
      return
      end

      subroutine tfreplacepart(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,kn,k1,kf
      type (sad_list), pointer :: list,kln,klx
      integer*4 irtc,isp1,narg,i,mode,itfmessage
      logical*4 single,seq,rep
      narg=isp-isp1
      if(mode .lt. 3)then
        if(narg .ne. 3)then
          irtc=itfmessage(9,'General::narg','"3"')
          return
        endif
      else
        if(narg .ne. 2)then
          irtc=itfmessage(9,'General::narg','"2"')
          return
        endif
      endif
      kn=dtastk(isp)
      if(mode .eq. 1)then
        k=dtastk(isp-1)
        if(ktfnonlistqd(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
          return
        endif
        kf=dtfcopy(dtastk(isp-2))
      else
        k=dtastk(isp1+1)
        if(ktfnonlistqd(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
          return
        endif
      endif
      call tfclonelist(list,list)
      if(ktfrealqd(kn))then
        single=.true.
      elseif(tflistqd(kn,kln))then
        k1=kln%dbody(1)
        single=ktfnonlistqd(k1)
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
        go to 9000
      endif
      if(single)then
        call tfreplacepart1(mode,list,kn,kf,seq,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      else
        do i=1,kln%nl
          call tfreplacepart1(mode,list,kln%dbody(i),kf,seq,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
        enddo
      endif
      if(seq)then
        call tfrebuildl(list,klx,rep)
      else
        list%attr=ior(list%attr,ktoberebuilt)-ktoberebuilt
        klx=>list
      endif
      call tfleval(klx,kx,.true.,irtc)
 9000 if(mode .eq. 1)then
        call tflocald(kf)
      endif
      return
      end

      subroutine tfreplacepart1(mode,kln,kn,kf,seq,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ki,kn,kf,kxi
      type (sad_list) kln
      type (sad_list), pointer :: kl,kli,klxi
      integer*4 mode,irtc,isp0,ivi,itfmessage,isp2,isp3,i
      logical*4 list,seq,seq1
      isp0=isp
      if(ktfrealqd(kn))then
        isp=isp+1
        dtastk(isp)=kn
      elseif(tflistqd(kn,kl))then
        call tfgetllstkall(kl)
        if(isp .eq. isp0)then
          irtc=0
          return
        endif
      else
        isp=isp0
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
        return
      endif
      isp2=isp
      call tfpartrstk(kln,isp0,isp2,list,
     $     mode .eq. 2,.true.,seq,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      if(mode .eq. 0)then
        do i=isp2+1,isp
          call loc_sad(ktastk(i),kli)
          call tfreplist(kli,itastk2(1,i),dtastk(isp0-1),seq1)
          seq=seq .or. seq1
        enddo
      elseif(mode .eq. 1)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_list(ktastk(i),kli)
          ivi=itastk2(1,i)
          isp=isp3+1
          dtastk(isp)=kf
          isp=isp+1
          dtastk(isp)=kli%dbody(ivi)
          call tfefunref(isp3+1,ki,.true.,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          call tfreplist(kli,ivi,ki,seq1)
          seq=seq .or. seq1
          isp=isp3
        enddo
      elseif(mode .eq. 2)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_list(ktastk(i),kli)
          ivi=itastk2(1,i)
          if(ivi .gt. kli%nl)then
            if(ivi .eq. 1)then
              kli%nl=1
              call tfreplist(kli,ivi,dtastk(isp0-1),seq1)
            else
              isp=isp3
              isp=isp+1
              ktastk(isp)=kli%body(ivi-1)
              isp=isp+1
              ktastk(isp)=ktastk(isp0-1)
              kxi=kxmakelist(isp3,klxi)
              klxi%head=ktfoper+mtfnull
              call tfreplist(kli,ivi,kxi,seq1)
            endif
          else
            isp=isp3+1
            ktastk(isp)=ktastk(isp0-1)
            isp=isp+1
            ktastk(isp)=kli%body(ivi)
            kxi=kxmakelist(isp3,klxi)
            klxi%head=ktfoper+mtfnull
            call tfreplist(kli,ivi,kxi,seq1)
          endif
          seq=seq .or. seq1
        enddo
      elseif(mode .eq. 3)then
        do i=isp2+1,isp
          call loc_list(ktastk(i),kli)
          call tfreplist(kli,itastk2(1,i),dxnull,seq1)
          seq=seq .or. seq1
        enddo
      elseif(mode .eq. 4)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_list(ktastk(i),kli)
          ivi=itastk2(1,i)
          ki=kli%dbody(ivi)
          if(ktflistqd(ki,kl))then
            call tfgetllstkall(kl)
            kxi=kxmakelist(isp3,klxi)
            klxi%head=ktfoper+mtfnull
            isp=isp3
            call tfreplist(kli,ivi,kxi,seq1)
            seq=seq .or. seq1
          endif
        enddo
      endif
 9000 isp=isp0
      return
      end

      subroutine tfpartition(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kp,k,ks
      type (sad_list), pointer :: klp,kl,kls
      integer*4 isp1,irtc,narg,i,m,id,itfmessage
      real*8 vp,vs
      narg=isp-isp1
      k=dtastk(isp1+1)
      if(ktfnonlistqd(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      kp=dtastk(isp1+2)
      if(narg .eq. 2)then
        if(ktfrealqdi(kp,itastk(1,isp+1)))then
          isp=isp+1
          itastk(2,isp)=itastk(1,isp)
          call tfpartitionstk(isp,isp,kl,kx,irtc)
          isp=isp1+2
        elseif(tflistqd(kp,klp) .and. ktfreallistqo(klp))then
          do i=1,klp%nl
            isp=isp+1
            itastk(1,isp)=int(klp%rbody(i))
            itastk(2,isp)=itastk(1,isp)
          enddo
          call tfpartitionstk(isp1+3,isp,kl,kx,irtc)
          isp=isp1+2
        else
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real or List of Reals for index"')
        endif
        return
      elseif(narg .eq. 3)then
        ks=dtastk(isp)
        if(ktfrealqd(kp,vp) .and. ktfrealqd(ks,vs))then
          isp=isp+1
          itastk(1,isp)=vp
          itastk(2,isp)=vs
          call tfpartitionstk(isp,isp,kl,kx,irtc)
          isp=isp1+3
          return
        elseif(tflistqd(kp,klp) .and. ktfreallistqo(klp))then
          m=klp%nl
          if(tflistqd(ks,kls) .and. ktfreallistqo(kls) .and.
     $         kls%nl .eq. m)then
            do i=1,m
              isp=isp+1
              itastk(1,isp)=int(klp%rbody(i))
              itastk(2,isp)=int(kls%rbody(i))
            enddo
          elseif(ktfrealqdi(ks,id))then
            do i=1,m
              isp=isp+1
              itastk(1,isp)=int(klp%rbody(i))
              itastk(2,isp)=id
            enddo
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Real or List of Reals for index"')
          endif
          call tfpartitionstk(isp1+4,isp,kl,kx,irtc)
          isp=isp1+3
          return
        endif
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
      else
        irtc=itfmessage(9,'General::narg','"2 or 3"')
      endif
      return
      end

      recursive subroutine tfpartitionstk(isp10,isp20,kl,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kk,kxthread
      type (sad_list) kl
      type (sad_list), pointer :: klj
      integer*8 kai
      integer*4 isp1,ma,ne,no,i,mp,isp10,isp20,isp2,isp3,isp0,irtc,j,
     $     itfmessage
      data kxthread%k /0/
      if(kxthread%k .eq. 0)then
        kxthread=kxsymbolz('`System`Thread',14)
      endif
      isp1=isp10
      isp2=isp20
      ma=kl%nl
      ne=itastk(1,isp1)
      if(ne .le. 0)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"#2","positive number or list of them"')
        return
      endif
      no=itastk(2,isp1)
      mp=(ma-ne+no)/no
      isp0=isp
      isp=isp+1
      ktastk(isp)=kl%head
      if(mp .gt. 0)then
        if(isp1 .eq. isp2)then
          do i=1,mp
            kai=(i-1)*no
            isp=isp+1
            isp3=isp
            ktastk(isp)=kl%head
            ktastk(isp+1:isp+ne)=kl%body(kai+1:kai+ne)
            isp=isp+ne
            call tfefunrefstk(isp3,isp3,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
        else
          do i=1,mp
            kai=(i-1)*no
            isp=isp+1
            isp3=isp
            dtastk(isp)=kxthread
            isp=isp+1
            ktastk(isp)=kl%head
            do j=1,ne
              isp=isp+1
              ktastk(isp)=kl%body(kai+j)
              if(ktflistq(ktastk(isp),klj))then
                call tfpartitionstk(isp1+1,isp2,klj,kk,irtc)
                if(irtc .ne. 0)then
                  isp=isp0
                  return
                endif
                dtastk(isp)=kk
              endif
            enddo
            dtastk(isp3+1)=kxcompose(isp3+1)
            isp=isp3+1
            call tfefunrefstk(isp3,isp3,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
        endif
      endif
      call tfefunref(isp0+1,kx,.true.,irtc)
      isp=isp0
      return
      end

      subroutine tfthread(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,kh,kf,ki,kj
      type (sad_list), pointer :: list,klj,klx,kli
      integer*8 kai
      integer*4 isp1,irtc,narg,mode,
     $     i,j,n,m,itfmessage,isp0,isp2,isp3
      logical*4 tfconstheadqk,allv,ch
      narg=isp-isp1
      if(mode .eq. 0)then
        if(narg .eq. 1)then
          kh%k=ktfoper+mtflist
        elseif(narg .eq. 2)then
          kh=dtastk(isp)
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
          return
        endif
        k=dtastk(isp1+1)
        if(ktfnonlistqd(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition for #1"')
          return
        endif
        kf=list%dbody(0)
      else
        if(narg .ne. 2)then
          irtc=itfmessage(9,'General::narg','"2"')
          return
        endif
        k=dtastk(isp)
        if(ktfnonlistqd(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition for #2"')
          return
        endif
        kh%k=ktfoper+mtflist
        kf=dtastk(isp1+1)
      endif
      m=list%nl
      n=-1
      ch=tfconstheadqk(kf%k)
      allv=ch .and. kh%k .eq. ktfoper+mtflist
      isp3=isp
      if(ktfnonreallistqo(list))then
        do i=1,m
          isp=isp+1
          ktastk(isp)=ktfref
          ki=list%dbody(i)
          if(ktflistqd(ki,kli) .and. kli%head .eq. kh%k)then
            ktastk(isp)=ktfaddr(ki)
            if(n .lt. 0)then
              n=kli%nl
            elseif(n .ne. kli%nl)then
              irtc=itfmessage(9,'General::equalleng',
     $             '"elements in Thread"')
              return
            endif
            allv=allv .and. ktfreallistqo(kli)
          else
            allv=.false.
          endif
        enddo
      endif
      if(n .lt. 0)then
        kx=k
        irtc=0
        return
      endif
      isp=isp+1
      isp2=isp
      dtastk(isp)=kh
      if(ch)then
        if(mode .ne. 2)then
          if(allv)then
            kx=kxadaloc(-1,n,klx)
            do j=1,n
              klx%dbody(j)=kxavaloc(0,m,klj)
              do i=1,m
                klj%rbody(i)=rlist(ktastk(isp3+i)+j)
              enddo
              klj%dbody(0)=dtfcopy(kf)
              klj%attr=ior(klj%attr,lconstlist)
            enddo
            klx%attr=ior(klx%attr,lconstlist)
            irtc=0
            isp=isp2-1
            return
          else
            do j=1,n
              isp=isp+1
              isp0=isp
              dtastk(isp0)=kf
              do i=1,m
                isp=isp+1
                kai=ktastk(isp3+i)
                if(kai .eq. ktfref)then
                  ktastk(isp)=list%body(i)
                else
                  ktastk(isp)=klist(kai+j)
                endif
              enddo
              dtastk(isp0)=kxcompose(isp0)
              isp=isp0
            enddo
          endif
        endif
      elseif(mode .ne. 2)then
        do j=1,n
          isp=isp+1
          isp0=isp
          dtastk(isp0)=kf
          do i=1,m
            isp=isp+1
            kai=ktastk(isp3+i)
            if(kai .eq. ktfref)then
              ktastk(isp)=list%body(i)
            else
              ktastk(isp)=klist(kai+j)
            endif
          enddo
          call tfefunrefstk(isp0,isp0,irtc)
          if(irtc .ne. 0)then
            isp=isp2-1
            return
          endif
        enddo
      else
        do j=1,n
          isp=isp+1
          isp0=isp
          dtastk(isp0)=kf
          do i=1,m
            isp=isp+1
            kai=ktastk(isp3+i)
            if(kai .eq. ktfref)then
              ktastk(isp)=list%body(i)
            else
              ktastk(isp)=klist(kai+j)
            endif
          enddo
          call tfefunref(isp0,kj,.true.,irtc)
          if(irtc .ne. 0)then
            if(irtc .eq. -3)then
              go to 1000
            elseif(irtc .eq. -2)then
              irtc=0
            else
              isp=isp2-1
              return
            endif
          endif
          isp=isp0-1
        enddo
      endif
 1000 if(mode .ne. 2)then
        call tfefunref(isp2,kx,.true.,irtc)
      else
        kx%k=ktfoper+mtfnull
        irtc=0
      endif
      isp=isp2-1
      return
      end

      subroutine tfsetpart(kln,k2,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k2,kr,k1
      type (sad_list) kln
      type (sad_list), pointer :: listl,list2,kli
      integer*8 kap
      integer*4 irtc,i,isp1,isp2,n,itfmessageexp,mode
      logical*4 list,def,eval,eval1,tfgetstoredp,app
      if(tfgetstoredp(kln%body(1),kap,def,irtc))then
        k1=dlist(kap)
        if(ktfnonlistqd(k1,listl))then
          irtc=itfmessageexp(999,'General::invset',k1)
          return
        endif
      else
        if(irtc .eq. 0)then
          irtc=itfmessageexp(999,'General::invset',
     $         ktflist+ksad_loc(kln%head))
        endif
        return
      endif
      isp1=isp
      call tfgetllstk(kln,2,-1)
      isp2=isp
      call tfclonelist(listl,listl)
      call tfpartrstk(listl,isp1,isp2,list,
     $     .false.,.true.,eval,.true.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(mode .eq. 0)then
        if(list)then
          if(tflistqd(k2,list2))then
            n=list2%nl
            if(n .eq. isp-isp2)then
              do i=1,n
                call loc_sad(ktastk(isp2+i),kli)
                call tfreplist(kli,
     $               itastk2(1,isp2+i),list2%dbody(i),eval1)
                eval=eval .or. eval1
              enddo
              go to 100
            endif
          endif
          do i=isp2+1,isp
            call loc_sad(ktastk(i),kli)
            call tfreplist(kli,itastk2(1,i),k2,eval1)
            eval=eval .or. eval1
          enddo
        elseif(isp .eq. isp2+1)then
          call loc_sad(ktastk(isp),kli)
          call tfreplist(kli,itastk2(1,isp),k2,eval1)
          eval=eval .or. eval1
        endif
      else
        app=mode .eq. 1
        if(list)then
          if(tflistqd(k2,list2))then
            n=list2%nl
            if(n .eq. isp-isp2)then
              do i=1,n
                call tfappendto1(ktastk(isp2+i)+itastk2(1,isp2+i),
     $               list2%body(i),kr,app,eval1,irtc)
                if(irtc .ne. 0)then
                  go to 1000
                endif
                if(eval1)then
                  call loc_sad(ktastk(isp2+i),kli)
                  call tfreplist(kli,itastk2(1,isp2+i),kr,eval1)
                  eval=eval .or. eval1
                endif
              enddo
              go to 100
            endif
          endif
          do i=1,n
            call tfappendto1(ktastk(isp2+i)+itastk2(1,isp2+i),
     $           k2,kr,app,eval1,irtc)
            if(irtc .ne. 0)then
              go to 1000
            endif
            if(eval1)then
              call loc_sad(ktastk(isp2+i),kli)
              call tfreplist(kli,itastk2(1,isp2+i),kr,eval1)
              eval=eval .or. eval1
            endif
          enddo
        elseif(isp .eq. isp2+1)then
          call tfappendto1(ktastk(isp)+itastk2(1,isp),
     $         k2,kr,app,eval1,irtc)
          if(irtc .ne. 0)then
            go to 1000
          endif
          if(eval1)then
            call loc_sad(ktastk(isp),kli)
            call tfreplist(kli,itastk2(1,isp),kr,eval1)
            eval=eval .or. eval1
          endif
        endif
      endif
 100  if(eval)then
        call tfleval(listl,kx,.true.,irtc)
        if(irtc .ne. 0)then
          go to 1000
        endif
      else
        kx%k=ktflist+ksad_loc(listl%head)
      endif
      call tflocal(klist(kap))
      dlist(kap)=dtfcopy(kx)
 1000 isp=isp1
      kx=k2
      return
      end

      subroutine tfreplist(list,iv,k,eval)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list) list
      integer*4 iv,i,lattr
      logical*4 eval,tfconstqk,tfconstlistqo
      eval=.false.
      if(iv .eq. 0)then
        call tflocal(list%head)
        list%dbody(0)=dtfcopy(k)
        list%attr=iand(list%attr,lnonreallist+ktoberebuilt)
      elseif(ktfreallistqo(list))then
        if(ktfrealqd(k))then
          list%dbody(iv)=k
          list%attr=ior(list%attr,ktoberebuilt)-ktoberebuilt
        else
          list%dbody(iv)=dtfcopy(k)
          eval=ktfsequenceq(k%k)
          if(eval)then
            list%attr=ktoberebuilt+lnonreallist
          else
            list%attr=lnonreallist
          endif
        endif
      else
        call tflocal(list%body(iv))
        if(ktfnonrealqd(k))then
          list%dbody(iv)=dtfcopy(k)
          eval=ktfsequenceq(k%k)
          lattr=iand(list%attr,ktoberebuilt)+lnonreallist
          if(tfconstlistqo(list))then
            if(.not. tfconstqk(k%k) .or. eval)then
              list%attr=lattr
            endif
          else
            list%attr=lattr
          endif
        else
          list%dbody(iv)=k
          do i=1,list%nl
            if(ktfnonrealq(list%body(i)))then
              list%attr=iand(list%attr,ktoberebuilt)+lnonreallist
              return
            endif
          enddo
          list%attr=iand(list%attr,lconstlist+lnoconstlist)
        endif
      endif
      return
      end
