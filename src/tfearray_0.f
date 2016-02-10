      subroutine tfearray(k1,k,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_list), pointer :: l,l1
      integer*8 k1,k,kx,kt,kt1,ka,ka1,ktaaloc,kax,ky,ktfmakelist
      integer*4 irtc,ne,ne1,i,iopc1,isp0
      logical*4 tfexprqk,tfcomplexqk,list1,list,clist,clist1
      call tfecmplxl(k1,k,kx,iopc1,irtc)
      if(irtc .ne. -1)then
        return
      endif
      irtc=0
      clist=.false.
      clist1=.false.
      if(iopc1 .ge. mtfgreater .and. iopc1 .le. mtfnot)then
        go to 101
      endif
      ka1=ktfaddr(k1)
      kt1=k1-ka1
      ka=ktfaddr(k)
      kt =k-ka
      if(kt1 .eq. ktflist)then
        if(tfcomplexqk(k1,ne1))then
          list1=ne1 .gt. 0
          clist1=list1
        else
          if(tfexprqk(k1))then
            go to 101
          endif
          call list_loc(ka1,l1)
          ne1=l1%nl
          list1=.true.
        endif
      else
        ne1=0
        list1=.false.
      endif
      if(kt .eq. ktflist)then
        if(tfcomplexqk(k,ne))then
          list=ne .gt. 0
          clist=list
        else
          if(tfexprqk(k))then
            go to 101
          endif
          if(iopc1 .eq. mtfdot)then
            call tfdot(k1,k,kx,irtc)
            return
          endif
          call list_loc(ka,l)
          ne=l%nl
          list=.true.
          if(list1 .and. ne .ne. ne1)then
            if(iopc1 .eq. mtfequal)then
              kx=0
            elseif(iopc1 .eq. mtfunequal)then
              kx=ktftrue
            else
              go to 101
            endif
            return
          endif
        endif
      else
        list=.false.
        ne=0
      endif
      if(clist1)then
        call list_loc(ktfaddr(l1%body(1)),l1r)
        call list_loc(ktfaddr(l1%body(2)),l1i)
      endif
      if(clist)then
        call list_loc(ktfaddr(l%body(1)),lr)
        call list_loc(ktfaddr(l%body(2)),li)
      endif
      if(clist1)then
        if(clist)then
          if(iopc1 .eq. mtfequal)then
            kx=ktftrue
            do i=1,ne
              call tfcmplx(l1r%body(i),lr%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. 0)then
                kx=0
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
              call tfcmplx(l1i%body(i),li%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. 0)then
                kx=0
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
            enddo
            return
          elseif(iopc1 .eq. mtfunequal)then
            kx=0
            do i=1,ne
              call tfcmplx(l1r%body(i),lr%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. ktftrue)then
                kx=ktftrue
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
              call tfcmplx(l1i%body(i),li%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. ktftrue)then
                kx=ktftrue
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
            enddo
            return
          else
ccccc
            isp0=isp
            kax=ktaaloc(-1,ne)
            do i=1,ne
              isp=isp+1
              call tfcmplx(l1%body(i),l%body(i),ktastk(isp),
     $             iopc1,irtc)
              if(irtc .ne. 0)then
                isp=isp0
                return
              endif
            enddo
            kax=ktfmakelist(isp0)
            isp=isp0
          endif
        else
          if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
            go to 101
          endif
          isp0=isp
          do i=1,ne
            isp=isp+1
            call tfcmplx(l1%body(i),k,ktastk(isp),iopc1,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kax=ktfmakelist(isp0)
          isp=isp0
        endif


      if(list1)then
        if(list)then
          if(iopc1 .eq. mtfequal)then
            kx=ktftrue
            do i=1,ne
              call tfcmplx(l1%body(i),l%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. 0)then
                kx=0
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
            enddo
            return
          elseif(iopc1 .eq. mtfunequal)then
            kx=0
            do i=1,ne
              call tfcmplx(l1%body(i),l%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky .eq. ktftrue)then
                kx=ktftrue
                return
              elseif(ktfnonrealq(ky))then
                go to 101
              endif
            enddo
            return
          else
            isp0=isp
            kax=ktaaloc(-1,ne)
            do i=1,ne
              isp=isp+1
              call tfcmplx(l1%body(i),l%body(i),ktastk(isp),
     $             iopc1,irtc)
              if(irtc .ne. 0)then
                isp=isp0
                return
              endif
            enddo
            kax=ktfmakelist(isp0)
            isp=isp0
          endif
        else
          if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
            go to 101
          endif
          isp0=isp
          do i=1,ne
            isp=isp+1
            call tfcmplx(l1%body(i),k,ktastk(isp),iopc1,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kax=ktfmakelist(isp0)
          isp=isp0
        endif
      else
        if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
          go to 101
        endif
        isp0=isp
        do i=1,ne
          isp=isp+1
          call tfcmplx(k1,l%body(i),ktastk(isp),iopc1,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
        enddo
        kax=ktfmakelist(isp0)
        isp=isp0
      endif
      kx=ktflist+kax
      return
 101  call tfeexpr(k1,k,kx,iopc1)
      return
      end

      recursive subroutine tfecmplxl(k1,k2,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_list), pointer :: kl1,kl2,klx
      integer*8 k1,k2,kx,ktfmakelist,ktavaloc,kax,ka1,ka2,ka10,
     $     k1i,k2i,kxi,ktcalocv,kfromr,ir,ix1,ix2
      integer*4 irtc,m1,m2,i,iopc1,isp0,iopc2
      real*8 v1,v2i,rfromk
      complex*16 c1,cx,tfgetnumber,cx1,cx2,tfcmplxmathv
      logical*4 tfcomplexqk,tflistqk,tfnumberqk,d,c
      if(iopc1 .le. mtfnull .or.
     $     (iopc1 .gt. mtfpower .and. iopc1 .ne. mtfcomplex))then
        irtc=-1
        return
      endif
      irtc=0
      if(ktfrealq(k1))then
        v1=rfromk(k1)
        if(tflistqk(k2))then
          ka2=ktfaddr(k2)
          call list_loc(ka2,kl2)
          m2=kl2%nl
          if(m2 .eq. 0)then
            kx=kxnulll
            return
          endif
          if(ktfreallistqo(kl2))then
            if(iopc1 .eq. mtfcomplex)then
              isp0=isp
              d=.false.
              c=.true.
              do i=1,m2
                isp=isp+1
                if(kl2%body(i) .eq. 0)then
                  ktastk(isp)=k1
                else
                  ktastk(isp)=ktflist+
     $                 ktcalocv(-1,v1,kl2%rbody(i))
                  d=.true.
                endif
              enddo
              go to 8000
            endif
            kax=ktavaloc(-1,m2)
            call list_loc(kax,klx)
            kx=ktflist+kax
            klx%attr=ior(klx%attr,lconstlist)
            select case(iopc1)
            case (mtfplus)
              klx%rbody(1:m2)=kl2%rbody(1:m2)+v1
            case (mtftimes)
              klx%rbody(1:m2)=kl2%rbody(1:m2)*v1
            case (mtfrevpower)
              ir=int8(v1)
              if(dble(ir) .eq. v1)then
                if(ir .eq. -1)then
                  klx%rbody(1:m2)=1.d0/kl2%rbody(1:m2)
                else
                  klx%rbody(1:m2)=kl2%rbody(1:m2)**ir
                endif
              else
                klx%rbody(1:m2)=kl2%rbody(1:m2)**v1
              endif
            case (mtfpower)
              do i=1,m2
                ir=int8(kl2%rbody(i))
                if(dble(ir) .eq. kl2%rbody(i))then
                  if(ir .eq. -1)then
                    klx%rbody(i)=1.d0/v1
                  else
                    klx%rbody(i)=v1**ir
                  endif
                else
                  klx%rbody(i)=v1**kl2%rbody(i)
                endif
              enddo
            case default
              irtc=-1
            end select
            return
          else
            isp0=isp
            d=.false.
            c=.true.
            c1=v1
            select case(iopc1)
c            go to (1900,1900,3030,1900,3050,1900,3070,3080),iopc1
c                  m    i    +    -    *    /    v    ^
            case (mtfplus)
              do i=1,m2
                k2i=kl2%body(i)
                isp=isp+1
                if(tfnumberqk(k2i))then
                  cx=c1+tfgetnumber(k2i)
                  if(imag(cx) .eq. 0.d0)then
                    ktastk(isp)=kfromr(dble(cx))
                  else
                    ktastk(isp)=ktflist+
     $                   ktcalocv(-1,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  d=d .or. ktfnonrealq(ktastk(isp))
                endif
              enddo
              go to 8000
            case (mtftimes)
              do i=1,m2
                k2i=kl2%body(i)
                isp=isp+1
                if(tfnumberqk(k2i))then
                  cx=c1*tfgetnumber(k2i)
                  if(imag(cx) .eq. 0.d0)then
                    ktastk(isp)=kfromr(dble(cx))
                  else
                    ktastk(isp)=ktflist+
     $                   ktcalocv(-1,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  d=d .or. ktfnonrealq(ktastk(isp))
                endif
              enddo
              go to 8000
            case (mtfrevpower)
              ix1=int8(v1)
              if(dble(ix1) .eq. v1)then
                if(ix1 .eq. -1)then
                  do i=1,m2
                    k2i=kl2%body(i)
                    isp=isp+1
                    if(tfnumberqk(k2i))then
                      cx=1.d0/tfgetnumber(k2i)
                      if(imag(cx) .eq. 0.d0)then
                        ktastk(isp)=kfromr(dble(cx))
                      else
                        ktastk(isp)=ktflist+
     $                       ktcalocv(-1,dble(cx),imag(cx))
                        d=.true.
                      endif
                    else
                      call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                      if(irtc .ne. 0)then
                        return
                      endif
                      d=d .or. ktfnonrealq(ktastk(isp))
                    endif
                  enddo
                  go to 8000
                else
                  do i=1,m2
                    k2i=kl2%body(i)
                    isp=isp+1
                    if(tfnumberqk(k2i))then
                      cx=tfgetnumber(k2i)**ix1
                      if(imag(cx) .eq. 0.d0)then
                        ktastk(isp)=kfromr(dble(cx))
                      else
                        ktastk(isp)=ktflist+
     $                       ktcalocv(-1,dble(cx),imag(cx))
                        d=.true.
                      endif
                    else
                      call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                      if(irtc .ne. 0)then
                        return
                      endif
                      d=d .or. ktfnonrealq(ktastk(isp))
                    endif
                  enddo
                  go to 8000
                endif
              else
                do i=1,m2
                  k2i=kl2%body(i)
                  isp=isp+1
                  if(tfnumberqk(k2i))then
                    cx=tfgetnumber(k2i)**v1
                    if(imag(cx) .eq. 0.d0)then
                      ktastk(isp)=kfromr(dble(cx))
                    else
                      ktastk(isp)=ktflist+
     $                     ktcalocv(-1,dble(cx),imag(cx))
                      d=.true.
                    endif
                  else
                    call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                    if(irtc .ne. 0)then
                      return
                    endif
                    d=d .or. ktfnonrealq(ktastk(isp))
                  endif
                enddo
                go to 8000
              endif
            case (mtfpower)
              do i=1,m2
                k2i=kl2%body(i)
                isp=isp+1
                if(tfnumberqk(k2i))then
                  cx2=tfgetnumber(k2i)
                  if(imag(cx2) .eq. 0.d0)then
                    ix2=int8(dble(cx2))
                    if(ix2 .eq. dble(cx2))then
                      if(ix2 .eq. -1)then
                        cx=1.d0/v1
                      else
                        cx=v1**ix2
                      endif
                    else
                      cx=v1**dble(cx2)
                    endif
                  else
                    cx=v1**cx2
                  endif
                  if(imag(cx) .eq. 0.d0)then
                    ktastk(isp)=kfromr(dble(cx))
                  else
                    ktastk(isp)=ktflist+
     $                   ktcalocv(-1,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  d=d .or. ktfnonrealq(ktastk(isp))
                endif
              enddo
              go to 8000
            case default
              if(iopc1 .eq. mtfcomplex)then
                do i=1,m2
                  k2i=kl2%body(i)
                  isp=isp+1
                  if(tfnumberqk(k2i))then
                    cx2=tfgetnumber(k2i)
                    cx=dcmplx(v1-imag(cx2),dble(cx2))
                    if(imag(cx) .eq. 0.d0)then
                      ktastk(isp)=kfromr(dble(cx))
                    else
                      ktastk(isp)=ktflist+
     $                     ktcalocv(-1,dble(cx),imag(cx))
                      d=.true.
                    endif
                  else
                    call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                    if(irtc .ne. 0)then
                      return
                    endif
                    d=d .or. ktfnonrealq(ktastk(isp))
                  endif
                enddo
                go to 8000
              else
                irtc=-1
                isp=isp0
                return
              endif
            end select
          endif
        else
          call tfeexpr(k1,k2,kx,iopc1)
        endif
      elseif(tfcomplexqk(k1))then
        ka1=ktfaddr(k1)
        call list_loc(ka1,kl1)
        if(tflistqk(k2))then
          ka2=ktfaddr(k2)
          call list_loc(ka2,kl2)
          m2=ilist(2,ka2-1)
          c1=cfromr(kl1%rbody(1:2))
          isp0=isp
          d=.false.
          c=.true.
          if(ktfreallistqo(kl2))then
            call list_loc(ka2,kl2)
            do i=1,m2
              isp=isp+1
              v2i=kl2%rbody(i)
              cx=tfcmplxmathv(c1,dcmplx(v2i,0.d0),iopc1)
              if(imag(cx) .eq. 0.d0)then
                rtastk(isp)=dble(cx)
              else
                ktastk(isp)=ktflist+
     $               ktcalocv(-1,dble(cx),imag(cx))
                d=.true.
              endif
            enddo
          else
            do i=1,m2
              isp=isp+1
              k2i=kl2%body(i)
              if(tfnumberqk(k2i))then
                cx=tfcmplxmathv(c1,tfgetnumber(k2i),iopc1)
                if(imag(cx) .eq. 0.d0)then
                  rtastk(isp)=dble(cx)
                else
                  ktastk(isp)=ktflist+
     $                 ktcalocv(-1,dble(cx),imag(cx))
                  d=.true.
                endif
              else
                call tfeval1(k1,k2i,ktastk(isp),iopc1,irtc)
                if(ktfnonrealq(ktastk(isp)))then
                  c=.false.
                  d=.true.
                endif
                if(irtc .ne. 0)then
                  return
                endif
              endif
            enddo
          endif
          go to 8000
        else
          call tfeexpr(k1,k2,kx,iopc1)
          return
        endif
      elseif(tflistqk(k1))then
        if(.not. tflistqk(k2))then
          if(iopc1 .eq. mtfpower)then
            call tfecmplxl(k2,k1,kx,mtfrevpower,irtc)
          elseif(iopc1 .eq. mtfrevpower)then
            call tfecmplxl(k2,k1,kx,mtfpower,irtc)
          elseif(iopc1 .eq. mtfcomplex)then
            call tfeval1(0.d0,k2,kxi,mtfcomplex,irtc)
            if(irtc .ne. 0)then
              return
            endif
            call tfecmplxl(kxi,k1,kx,mtfplus,irtc)
          else
            call tfecmplxl(k2,k1,kx,iopc1,irtc)
          endif
          return
        endif
        ka1=ktfaddr(k1)
        ka2=ktfaddr(k2)
        call list_loc(ka1,kl1)
        call list_loc(ka2,kl2)
        m1=kl1%nl
        m2=kl2%nl
        if(m1 .ne. m2)then
          call tfeexpr(k1,k2,kx,iopc1)
          return
        endif
        if(ktfreallistqo(kl1) .and. ktfreallistqo(kl2))then
          if(iopc1 .eq. mtfcomplex)then
            isp0=isp
            d=.false.
            c=.true.
            do i=1,m2
              isp=isp+1
              if(kl2%rbody(i) .eq. 0.d0)then
                rtastk(isp)=kl1%rbody(i)
              else
                ktastk(isp)=ktflist+
     $               ktcalocv(-1,kl1%rbody(i),kl2%rbody(i))
                d=.true.
              endif
            enddo
            go to 8000
          endif
          kax=ktavaloc(-1,m1)
          kx=ktflist+kax
          call list_loc(kax,klx)
          klx%attr=lconstlist
          select case (iopc1)
          case (mtfplus)
            klx%rbody(1:m2)=kl2%rbody(1:m2)+kl1%rbody(1:m2)
          case (mtftimes)
            klx%rbody(1:m2)=kl2%rbody(1:m2)*kl1%rbody(1:m2)
          case (mtfrevpower)
            do i=1,m2
              ir=int8(kl1%rbody(i))
              if(ir .eq. kl1%rbody(i))then
                if(ir .eq. -1)then
                  klx%rbody(i)=1.d0/kl2%rbody(i)
                else
                  klx%rbody(i)=kl2%rbody(i)**ir
                endif
              else
                klx%rbody(i)=kl2%rbody(i)**kl1%rbody(i)
              endif
            enddo
          case (mtfpower)
            do i=1,m2
              ir=int8(kl2%rbody(i))
              if(ir .eq. kl2%rbody(i))then
                if(ir .eq. -1)then
                  klx%rbody(i)=1.d0/kl1%rbody(i)
                else
                  klx%rbody(i)=kl1%rbody(i)**ir
                endif
              else
                klx%rbody(i)=kl1%rbody(i)**kl2%rbody(i)
              endif
            enddo
          case default
            irtc=-1
          end select
          return
        else
          isp0=isp
          c=.true.
          iopc2=iopc1
 4000     select case (iopc2)
          case (mtfneg,mtfinv,mtfminus,mtfdiv)
            irtc=-1
            return
          case (mtfplus)
            do i=1,m2
              isp=isp+1
              k1i=kl1%body(i)
              k2i=kl2%body(i)
              if(tfnumberqk(k1i) .and. tfnumberqk(k2i))then
                cx=tfgetnumber(k1i)+tfgetnumber(k2i)
                if(imag(cx) .eq. 0.d0)then
                  ktastk(isp)=kfromr(dble(cx))
                else
                  ktastk(isp)=ktflist+
     $                 ktcalocv(-1,dble(cx),imag(cx))
                  d=.true.
                endif
              else
                call tfeval1(k1i,k2i,ktastk(isp),iopc1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ktfnonrealq(ktastk(isp)))then
                  d=.true.
                  c=.false.
                endif
              endif
            enddo
            go to 8000
          case (mtftimes)
            do i=1,m2
              isp=isp+1
              k1i=kl1%body(i)
              k2i=kl2%body(i)
              if(tfnumberqk(k1i) .and. tfnumberqk(k2i))then
                cx=tfgetnumber(k1i)*tfgetnumber(k2i)
                if(imag(cx) .eq. 0.d0)then
                  ktastk(isp)=kfromr(dble(cx))
                else
                  ktastk(isp)=ktflist+
     $                 ktcalocv(-1,dble(cx),imag(cx))
                  d=.true.
                endif
              else
                call tfeval1(k1i,k2i,ktastk(isp),iopc1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ktfnonrealq(ktastk(isp)))then
                  d=.true.
                  c=.false.
                endif
              endif
            enddo
            go to 8000
          case (mtfrevpower)
            do i=1,m2
              isp=isp+1
              k1i=kl1%body(i)
              k2i=kl2%body(i)
              if(tfnumberqk(k1i) .and. tfnumberqk(k2i))then
                cx1=tfgetnumber(k1i)
                if(imag(cx1) .eq. 0.d0)then
                  ix1=int8(cx1)
                  if(dble(ix1) .eq. dble(cx1))then
                    if(ix1 .eq. -1)then
                      cx=1.d0/tfgetnumber(k2i)
                    else
                      cx=tfgetnumber(k2i)**ix1
                    endif
                  else
                    cx=tfgetnumber(k2i)**dble(cx1)
                  endif
                else
                  cx=tfgetnumber(k2i)**cx1
                endif
                if(imag(cx) .eq. 0.d0)then
                  ktastk(isp)=kfromr(dble(cx))
                else
                  ktastk(isp)=ktflist+
     $                 ktcalocv(-1,dble(cx),imag(cx))
                  d=.true.
                endif
              else
                call tfeval1(k1i,k2i,ktastk(isp),iopc1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ktfnonrealq(ktastk(isp)))then
                  d=.true.
                  c=.false.
                endif
              endif
            enddo
            go to 8000
          case (mtfpower)
            iopc2=mtfrevpower
            ka10=ka1
            ka1=ka2
            ka2=ka10
            go to 4000
          case default
            if(iopc1 .eq. mtfcomplex)then
              do i=1,m2
                isp=isp+1
                k1i=kl1%body(i)
                k2i=kl2%body(i)
                if(tfnumberqk(k1i) .and. tfnumberqk(k2i))then
                  cx2=tfgetnumber(k2i)
                  cx=tfgetnumber(k1i)+dcmplx(-imag(cx2),dble(cx2))
                  if(imag(cx) .eq. 0.d0)then
                    ktastk(isp)=kfromr(dble(cx))
                  else
                    ktastk(isp)=ktflist+
     $                   ktcalocv(-1,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  call tfeval1(k1i,k2i,ktastk(isp),iopc1,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  if(ktfnonrealq(ktastk(isp)))then
                    d=.true.
                    c=.false.
                  endif
                endif
              enddo
              go to 8000
            else
              irtc=-1
              isp=isp0
              return
            endif
          end select
        endif
      else
        call tfeexpr(k1,k2,kx,iopc1)
        return
      endif
      return
 8000 kax=ktfmakelist(isp0)
      call list_loc(kax,klx)
      isp=isp0
      if(c)then
        klx%attr=ior(lconstlist,klx%attr)
      endif
      kx=ktflist+kax
      return
      end

      logical*4 function tfnumberqk(k)
      use tfstk
      implicit none
      integer*8 k,ka
      if(ktfrealq(k))then
        tfnumberqk=.true.
      elseif(ktflistq(k))then
        ka=ktfaddr(k)
        tfnumberqk=klist(ka) .eq. ktfoper+mtfcomplex .and.
     $       ktfreallistq(ka) .and. ilist(2,ka-1) .eq. 2 
      else
        tfnumberqk=.false.
      endif
      return
      end

      complex*16 function tfgetnumber(k)
      use tfstk
      implicit none
      type (sad_list), pointer :: kl
      integer*8 k
      real*8 rfromk
      if(ktfrealq(k))then
        tfgetnumber=rfromk(k)
      else
        call list_loc(ktfaddr(k),kl)
        tfgetnumber=dcmplx(kl%rbody(1),kl%rbody(2))
      endif
      return
      end

      subroutine tfcmplxmath(c1,c2,kx,iopc1,irtc)
      use tfstk
      implicit none
      integer*8 kx,ktcalocv,kfromr
      integer*4 iopc1,irtc
      complex*16 c1,c2,cx,tfcmplxmathv
      if(iopc1 .gt. mtfunequal .and. iopc1 .ne. mtfcomplex)then
        irtc=-1
      else
        cx=tfcmplxmathv(c1,c2,iopc1)
        irtc=0
        if(imag(cx) .eq. 0.d0)then
          kx=kfromr(cx)
        else
          kx=ktflist+ktcalocv(-1,dble(cx),imag(cx))
        endif
      endif
      return
      end

      complex*16 function tfcmplxmathv(c1,c2,iopc1)
      use tfstk
      implicit none
      integer*4 iopc1
      integer*8 i1,i2
      complex*16 c1,c2
      integer*8 itredmath,ktfsymbolz
      data itredmath /0/
      if(itredmath .eq. 0)then
        itredmath=ktfsymbolz('`System`$ReduceMath',19)-4
      endif
      select case(iopc1)
      case (mtfneg)
        tfcmplxmathv=-c2
      case (mtfinv)
        tfcmplxmathv=1.d0/c2
      case (mtfplus)
        tfcmplxmathv=c1+c2
      case (mtftimes)
        tfcmplxmathv=c1*c2
      case (mtfrevpower)
        if(imag(c1) .eq. 0.d0)then
          i1=int8(c1)
          if(i1 .eq. dble(c1))then
            if(i1 .eq. -1)then
              tfcmplxmathv=1.d0/c2
            else
              tfcmplxmathv=c2**i1
            endif
          else
            tfcmplxmathv=c2**dble(c1)
          endif
        else
          tfcmplxmathv=c2**c1
        endif
      case(mtfpower)
        if(imag(c2) .eq. 0.d0)then
          i2=int8(c2)
          if(i2 .eq. dble(c2))then
            if(i2 .eq. -1)then
              tfcmplxmathv=1.d0/c1
            elseif(i2 .eq. 0 .and.
     $             rlist(itredmath) .ne. 0.d0)then
              tfcmplxmathv=1.d0
            else
              tfcmplxmathv=c1**i2
            endif
          else
            tfcmplxmathv=c1**dble(c2)
          endif
        else
          tfcmplxmathv=c1**c2
        endif
      case (mtfequal)
        if(c1 .eq. c2)then
          tfcmplxmathv=1.d0
        else
          tfcmplxmathv=0.d0
        endif
      case (mtfunequal)
        if(c1 .ne. c2)then
          tfcmplxmathv=1.d0
        else
          tfcmplxmathv=0.d0
        endif
      case default
        if(iopc1 .eq. mtfcomplex)then
          tfcmplxmathv=c1+dcmplx(-imag(c2),dble(c2))
        else
          tfcmplxmathv=0.d0
        endif
      end select
      return
      end
