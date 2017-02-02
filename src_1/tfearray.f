      subroutine tfearray(k1,k,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,k1,kx,ky
      type (sad_list), pointer :: kl,kl1
      integer*8 isp0
      integer*4 irtc,ne,ne1,i,iopc1
      logical*4 tfexprqk,list1,list
c     begin initialize for preventing compiler warning
c     end   initialize for preventing compiler warning
c$$$      if(tfmatrixqd(k1,kl1))then
c$$$        if(tfmatrixqd(k,kl2))then
c$$$          call tfematrix(kl1,kl2,kx,iopc1,irtc)
c$$$        else
c$$$c     call tfematrix1(kl1,k,kx,iopc1,irtc)
c$$$        endif
c$$$        if(irtc .ne. 0)then
c$$$          go to 101
c$$$        endif
c$$$        return
c$$$      elseif(tfmatrixqd(k,kl2))then
c$$$c     call tfmatrix2(k1,kl2,kx,iopc1,irtc)
c$$$        if(irtc .ne. 0)then
c$$$          go to 101
c$$$        endif
c$$$        return
c$$$      endif
      if(iopc1 .ge. mtfplus .and. iopc1 .le. mtfpower
     $     .or. iopc1 .eq. mtfcomplex)then
        call tfecmplxl(k1,k,kx,iopc1)
        irtc=0
        return
      endif
      irtc=0
      if(iopc1 .ge. mtfgreater .and. iopc1 .le. mtfnot)then
        go to 101
      endif
      if(ktflistqd(k1,kl1))then
        if(tfcomplexqd(k1))then
          ne1=0
          list1=.false.
        else
          if(tfexprqk(k1%k))then
            go to 101
          endif
          ne1=kl1%nl
          list1=.true.
        endif
      else
        ne1=0
        list1=.false.
      endif
      if(ktflistqd(k,kl))then
        if(tfcomplexqd(k))then
          list=.false.
          ne=0
        else
          if(tfexprqk(k%k))then
            go to 101
          endif
          if(iopc1 .eq. mtfdot)then
            call tfdot(k1,k,kx,irtc)
            return
          endif
          ne=kl%nl
          list=.true.
          if(list1 .and. ne .ne. ne1)then
            if(iopc1 .eq. mtfequal)then
              kx%k=0
            elseif(iopc1 .eq. mtfunequal)then
              kx%k=ktftrue
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
      if(list1)then
        if(list)then
          if(iopc1 .eq. mtfequal)then
            kx%k=ktftrue
            do i=1,ne
              call tfcmplx(kl1%body(i),kl%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky%k .eq. 0)then
                kx%k=0
                return
              elseif(.not. ktfrealqd(ky))then
                go to 101
              endif
            enddo
            return
          elseif(iopc1 .eq. mtfunequal)then
            kx%k=0
            do i=1,ne
              call tfcmplx(kl1%body(i),kl%body(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky%k .eq. ktftrue)then
                kx%k=ktftrue
                return
              elseif(.not. ktfrealqd(ky))then
                go to 101
              endif
            enddo
            return
          else
            isp0=isp
            do i=1,ne
              isp=isp+1
              call tfcmplx(kl1%body(i),kl%body(i),dtastk(isp),
     $             iopc1,irtc)
              if(irtc .ne. 0)then
                isp=isp0
                return
              endif
            enddo
            kx=kxmakelist(isp0)
            isp=isp0
          endif
        else
          if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
            go to 101
          endif
          isp0=isp
          do i=1,ne
            isp=isp+1
            call tfcmplx(kl1%body(i),k,dtastk(isp),iopc1,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
        endif
      else
        if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
          go to 101
        endif
        isp0=isp
        do i=1,ne
          isp=isp+1
          call tfcmplx(k1,kl%dbody(i),dtastk(isp),iopc1,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
        enddo
        kx=kxmakelist(isp0)
        isp=isp0
      endif
      return
 101  call tfeexpr(k1,k,kx,iopc1)
      return
      end

      recursive subroutine tfecmplxl(k1,k2,kx,iopc1)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,kxi,k2i,k1i
      type (sad_list), pointer :: kl1,kl2,klx,kl10
      integer*8 ir,ix1,ix2
      integer*4 irtc,m1,m2,i,iopc1,iopc2
      real*8 v1,v2i
      complex*16 c1,cx,cx1,cx2,tfcmplxmathv
      logical*4 d,c,tfconstqk
      if(ktfrealqd(k1,v1))then
        if(tflistqd(k2,kl2))then
          m2=kl2%nl
          if(m2 .eq. 0)then
            kx=dxnulll
            return
          endif
          if(ktfreallistqo(kl2))then
            kx=kxavaloc(-1,m2,klx)
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
            case (mtfcomplex)
              d=.false.
              do i=1,m2
                if(kl2%body(i) .eq. 0)then
                  klx%dbody(i)=k1
                else
                  klx%dbody(i)=kxcalocv(0,v1,kl2%rbody(i))
                  d=.true.
                endif
              enddo
              if(d)then
                klx%attr=ior(klx%attr,lnonreallist)
              endif
            case default
              call tfeexpr(k1,k2,kx,iopc1)
            end select
            return
          else
            c=.true.
            c1=v1
            kx=kxadaloc(-1,m2,klx)
            d=.false.
            select case(iopc1)
c            go to (1900,1900,3030,1900,3050,1900,3070,3080),iopc1
c                  m    i    +    -    *    /    v    ^
            case (mtfplus)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberqd(k2i,cx2))then
                  cx=c1+cx2
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                elseif(tflistqd(k2i))then
                  call tfecmplxl(k1,k2i,kxi,mtfplus)
                  d=.true.
                  c=c .and. tfconstqk(kxi%k)
                  klx%dbody(i)=dtfcopy(kxi)
                else
                  call tfeexpr(k1,k2i,kxi,mtfplus)
                  if(ktfnonrealqd(kxi))then
                    c=c .and. tfconstqk(kxi%k)
                    d=.true.
                    kxi=dtfcopy(kxi)                  
                  endif
                  klx%dbody(i)=kxi
                endif
              enddo
              go to 9000
            case (mtftimes)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberqd(k2i,cx2))then
                  cx=c1*cx2
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                elseif(tflistqd(k2i))then
                  call tfecmplxl(k1,k2i,kxi,mtftimes)
                  d=.true.
                  c=c .and. tfconstqk(kxi%k)
                  klx%dbody(i)=dtfcopy(kxi)
                else
                  call tfeexpr(k1,k2i,kxi,mtftimes)
                  if(ktfnonrealqd(kxi))then
                    c=c .and. tfconstqk(kxi%k)
                    d=.true.
                    kxi=dtfcopy(kxi)                  
                  endif
                  klx%dbody(i)=kxi
                endif
              enddo
              go to 9000
            case (mtfrevpower)
              ix1=int8(v1)
              if(dble(ix1) .eq. v1)then
                if(ix1 .eq. -1)then
                  do i=1,m2
                    k2i=kl2%dbody(i)
                    if(tfnumberqd(k2i,cx2))then
                      cx=1.d0/cx2
                      if(imag(cx) .eq. 0.d0)then
                        klx%rbody(i)=dble(cx)
                      else
                        klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                        d=.true.
                      endif
                    elseif(tflistqd(k2i))then
                      call tfecmplxl(k1,k2i,kxi,mtfrevpower)
                      d=.true.
                      c=c .and. tfconstqk(kxi%k)
                      klx%dbody(i)=dtfcopy(kxi)
                    else
                      call tfeexpr(k1,k2i,kxi,mtfrevpower)
                      if(ktfnonrealqd(kxi))then
                        c=c .and. tfconstqk(kxi%k)
                        d=.true.
                        kxi=dtfcopy(kxi)                  
                      endif
                      klx%dbody(i)=kxi
                    endif
                  enddo
                  go to 9000
                else
                  do i=1,m2
                    k2i=kl2%dbody(i)
                    if(tfnumberqd(k2i,cx2))then
                      cx=cx2**ix1
                      if(imag(cx) .eq. 0.d0)then
                        klx%rbody(i)=dble(cx)
                      else
                        klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                        d=.true.
                      endif
                    elseif(tflistqd(k2i))then
                      call tfecmplxl(k1,k2i,kxi,mtfrevpower)
                      d=.true.
                      c=c .and. tfconstqk(kxi%k)
                      klx%dbody(i)=dtfcopy(kxi)
                    else
                      call tfeexpr(k1,k2i,kxi,mtfrevpower)
                      if(ktfnonrealqd(kxi))then
                        c=c .and. tfconstqk(kxi%k)
                        d=.true.
                        kxi=dtfcopy(kxi)                  
                      endif
                      klx%dbody(i)=kxi
                    endif
                  enddo
                  go to 9000
                endif
              else
                do i=1,m2
                  k2i=kl2%dbody(i)
                  if(tfnumberqd(k2i,cx2))then
                    cx=cx2**cx1
                    if(imag(cx) .eq. 0.d0)then
                      klx%rbody(i)=dble(cx)
                    else
                      klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                      d=.true.
                    endif
                  elseif(tflistqd(k2i))then
                    call tfecmplxl(k1,k2i,kxi,mtfrevpower)
                    d=.true.
                    c=c .and. tfconstqk(kxi%k)
                    klx%dbody(i)=dtfcopy(kxi)
                  else
                    call tfeexpr(k1,k2i,kxi,mtfrevpower)
                    if(ktfnonrealqd(kxi))then
                      c=c .and. tfconstqk(kxi%k)
                      d=.true.
                      kxi=dtfcopy(kxi)                  
                    endif
                    klx%dbody(i)=kxi
                  endif
                enddo
                go to 9000
              endif
            case (mtfpower)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberqd(k2i,cx2))then
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
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                elseif(tflistqd(k2i))then
                  call tfecmplxl(k1,k2i,kxi,mtfpower)
                  d=.true.
                  c=c .and. tfconstqk(kxi%k)
                  klx%dbody(i)=dtfcopy(kxi)
                else
                  call tfeexpr(k1,k2i,kxi,mtfpower)
                  if(ktfnonrealqd(kxi))then
                    c=c .and. tfconstqk(kxi%k)
                    d=.true.
                    kxi=dtfcopy(kxi)                  
                  endif
                  klx%dbody(i)=kxi
                endif
              enddo
              go to 9000
            case (mtfcomplex)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberqd(k2i,cx2))then
                  cx=dcmplx(v1-imag(cx2),dble(cx2))
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                elseif(tflistqd(k2i))then
                  call tfecmplxl(k1,k2i,kxi,mtfcomplex)
                  d=.true.
                  c=c .and. tfconstqk(kxi%k)
                  klx%dbody(i)=dtfcopy(kxi)
                else
                  call tfeexpr(k1,k2i,kxi,mtfcomplex)
                  if(ktfnonrealqd(kxi))then
                    c=c .and. tfconstqk(kxi%k)
                    d=.true.
                    kxi=dtfcopy(kxi)                  
                  endif
                  klx%dbody(i)=kxi
                endif
              enddo
              go to 9000
            case default  
              call tfeexpr(k1,k2,kx,iopc1)
              return
            end select
          endif
        else
          call tfeexpr(k1,k2,kx,iopc1)
        endif
      elseif(tfcomplexqk(k1%k,kl1))then
        if(tflistqd(k2,kl2))then
          m2=kl2%nl
          c1=kl1%cbody(1)
          kx=kxadaloc(-1,m2,klx)
          c=.true.
          d=.false.
          if(ktfreallistqo(kl2))then
            do i=1,m2
              v2i=kl2%rbody(i)
              cx=tfcmplxmathv(c1,dcmplx(v2i,0.d0),iopc1)
              if(imag(cx) .eq. 0.d0)then
                klx%rbody(i)=dble(cx)
              else
                klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                d=.true.
              endif
            enddo
          else
            do i=1,m2
              k2i=kl2%dbody(i)
              if(tfnumberqd(k2i,cx2))then
                cx=tfcmplxmathv(c1,cx2,iopc1)
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
              elseif(tflistqd(k2i))then
                call tfecmplxl(k1,k2i,kxi,iopc1)
                d=.true.
                c=c .and. tfconstqk(kxi%k)
                klx%dbody(i)=dtfcopy(kxi)
              else
                call tfeexpr(k1,k2i,kxi,iopc1)
                if(ktfnonrealqd(kxi))then
                  c=c .and. tfconstqk(kxi%k)
                  d=.true.
                  kxi=dtfcopy(kxi)                  
                endif
                klx%dbody(i)=kxi
              endif
            enddo
          endif
          go to 9000
        else
          call tfeexpr(k1,k2,kx,iopc1)
          return
        endif
      elseif(tflistqd(k1,kl1))then
        if(.not. tflistqd(k2,kl2))then
          if(iopc1 .eq. mtfpower)then
            call tfecmplxl(k2,k1,kx,mtfrevpower)
          elseif(iopc1 .eq. mtfrevpower)then
            call tfecmplxl(k2,k1,kx,mtfpower)
          elseif(iopc1 .eq. mtfcomplex)then
            call tfeval1(0.d0,k2,kxi,mtfcomplex,irtc)
            if(irtc .ne. 0)then
              call tfreseterror
              call tfeexpr(k1,k2,kx,mtfcomplex)
            else
              call tfecmplxl(kxi,k1,kx,mtfplus)
            endif
          else
            call tfecmplxl(k2,k1,kx,iopc1)
          endif
          return
        endif
        m1=kl1%nl
        if(m1 .ne. kl2%nl)then
          call tfeexpr(k1,k2,kx,iopc1)
          return
        endif
        if(ktfreallistqo(kl1) .and. ktfreallistqo(kl2))then
          kx=kxavaloc(-1,m1,klx)
          klx%attr=lconstlist
          select case (iopc1)
          case (mtfplus)
            klx%rbody(1:m1)=kl2%rbody(1:m1)+kl1%rbody(1:m1)
          case (mtftimes)
            klx%rbody(1:m1)=kl2%rbody(1:m1)*kl1%rbody(1:m1)
          case (mtfrevpower)
            do i=1,m1
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
            do i=1,m1
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
          case (mtfcomplex)
            d=.false.
            do i=1,m1
              if(kl2%rbody(i) .eq. 0.d0)then
                klx%rbody(i)=kl1%rbody(i)
              else
                klx%dbody(i)=kxcalocv(0,kl1%rbody(i),kl2%rbody(i))
                d=.true.
              endif
            enddo
            if(d)then
              klx%attr=ior(klx%attr,lnonreallist+lconstlist)
            endif
          case default
            
          end select
          return
        else
          c=.true.
          iopc2=iopc1
 4000     select case (iopc2)
          case (mtfplus)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(tfnumberqd(k1i,cx1) .and. tfnumberqd(k2i,cx2))then
                cx=cx1+cx2
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
                cycle
              elseif(tflistqd(k1i) .or. tflistqd(k2i))then
                call tfecmplxl(k1i,k2i,kxi,mtfplus)
              else
                call tfeexpr(k1i,k2i,kxi,mtfplus)
              endif
              if(ktfnonrealqd(kxi))then
                c=c .and. tfconstqk(kxi%k)
                d=.true.
                kxi=dtfcopy(kxi)                  
              endif
              klx%dbody(i)=kxi
            enddo
            go to 9000
          case (mtftimes)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(tfnumberqd(k1i,cx1) .and. tfnumberqd(k2i,cx2))then
                cx=cx1*cx2
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
                cycle
              elseif(tflistqd(k1i) .or. tflistqd(k2i))then
                call tfecmplxl(k1i,k2i,kxi,mtftimes)
              else
                call tfeexpr(k1i,k2i,kxi,mtftimes)
              endif
              if(ktfnonrealqd(kxi))then
                c=c .and. tfconstqk(kxi%k)
                d=.true.
                kxi=dtfcopy(kxi)                  
              endif
              klx%dbody(i)=kxi
            enddo
            go to 9000
          case (mtfrevpower)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(tfnumberqd(k1i,cx1) .and. tfnumberqd(k2i,cx2))then
                if(imag(cx1) .eq. 0.d0)then
                  ix1=int8(cx1)
                  if(dble(ix1) .eq. dble(cx1))then
                    if(ix1 .eq. -1)then
                      cx=1.d0/cx2
                    else
                      cx=cx2**ix1
                    endif
                  else
                    cx=cx2**dble(cx1)
                  endif
                else
                  cx=cx2**cx1
                endif
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
                cycle
              elseif(tflistqd(k1i) .or. tflistqd(k2i))then
                call tfecmplxl(k1i,k2i,kxi,mtfrevpower)
              else
                call tfeexpr(k1i,k2i,kxi,mtfrevpower)
              endif
              if(ktfnonrealqd(kxi))then
                c=c .and. tfconstqk(kxi%k)
                d=.true.
                kxi=dtfcopy(kxi)                  
              endif
              klx%dbody(i)=kxi
            enddo
            go to 9000
          case (mtfpower)
            iopc2=mtfrevpower
            kl10=>kl1
            kl1=>kl2
            kl2=>kl10
            go to 4000
          case (mtfcomplex)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(tfnumberqd(k1i,cx1) .and. tfnumberqd(k2i,cx2))then
                cx=cx1+dcmplx(-imag(cx2),dble(cx2))
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                endif
                cycle
              elseif(tflistqd(k1i) .or. tflistqd(k2i))then
                call tfecmplxl(k1i,k2i,kxi,mtfcomplex)
              else
                call tfeexpr(k1i,k2i,kxi,mtfcomplex)
              endif
              if(ktfnonrealqd(kxi))then
                c=c .and. tfconstqk(kxi%k)
                d=.true.
                kxi=dtfcopy(kxi)                  
              endif
              klx%dbody(i)=kxi
            enddo
            go to 9000
          case default
            call tfeexpr(k1,k2,kx,iopc1)
          end select
        endif
      else
        call tfeexpr(k1,k2,kx,iopc1)
        return
      endif
      return
 9000 if(.not. d)then
        klx%attr=0
      endif
      if(c)then
        klx%attr=ior(klx%attr,lconstlist)
      endif
      return
      end

      logical*4 function tfgetcomplex(k,c)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_complex), pointer :: cx
      complex*16 c
      tfgetcomplex=ktflistqx(k%k,cx) .and.
     $     cx%head .eq. ktfoper+mtfcomplex .and.
     $     iand(lnonreallist,cx%attr) .eq. 0 .and. cx%nl .eq. 2 
      if(tfgetcomplex)then
        c=cx%cx(1)
      endif
      return
      end

      subroutine tfcmplxmath(c1,c2,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 iopc1,irtc
      complex*16 c1,c2,cx,tfcmplxmathv
      if(iopc1 .gt. mtfunequal .and. iopc1 .ne. mtfcomplex)then
        irtc=-1
      else
        cx=tfcmplxmathv(c1,c2,iopc1)
        irtc=0
        if(imag(cx) .eq. 0.d0)then
          kx=dfromr(dble(cx))
        else
          kx=kxcalocv(-1,dble(cx),imag(cx))
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
            elseif(i2 .eq. 0 .and. redmath%value%k .ne. 0)then
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
      case (mtfcomplex)
        tfcmplxmathv=c1+dcmplx(-imag(c2),dble(c2))
      case default
        tfcmplxmathv=0.d0
      end select
      return
      end
