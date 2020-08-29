      module earray

      contains
      function tfecmplx1(k1,k2i,iopc1,c,d) result(kxi)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kxi,tfecmplxl
      type (sad_descriptor) ,intent(in):: k1,k2i
      integer*4 ,intent(in):: iopc1
      logical*4 ,intent(inout):: c,d
      if(tflistq(k2i))then
        kxi=dtfcopy(tfecmplxl(k1,k2i,iopc1))
        d=.true.
        c=c .and. tfconstq(kxi%k)
        kxi=dtfcopy(kxi)
      else
        kxi=tfeexpr(k1,k2i,iopc1)
        if(ktfnonrealq(kxi))then
          c=c .and. tfconstq(kxi%k)
          d=.true.
          kxi=dtfcopy(kxi)                  
        endif
      endif
      return
      end function

      end module

      recursive function tfecmplxl(k1,k2,iopc1) result(kx)
      use tfstk
      use eexpr
      use earray
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kxi,k2i,k1i,tfeval1
      type (sad_dlist), pointer :: kl1,kl2,klx,kl10
      type (sad_rlist), pointer :: klr,klr1,klr2
      integer*4 ,parameter::lseg=1000000
      integer*8 ir,ix1,ix2
      integer*4 ,intent(in):: iopc1
      integer*4 irtc,m1,m2,i,iopc2
      real*8 v1,v1i,v2i
      complex*16 c1,cx,cx1,cx2,tfcmplxmathv
      logical*4 d,c
c      call tfdebugprint(k1,'cmplxl',1)
c      call tfdebugprint(k2,'and',1)
c      write(*,*)'with ',iopc1
      if(ktfrealq(k1,v1))then
        if(tflistq(k2,kl2))then
          m2=kl2%nl
          if(m2 .eq. 0)then
            kx=dxnulll
            return
          endif
          if(ktfreallistq(k2,klr2))then
            kx=kxavaloc(-1,m2,klr)
            klr%attr=ior(klr%attr,lconstlist)
            select case(iopc1)
            case (mtfplus)
              klr%rbody(1:m2)=klr2%rbody(1:m2)+v1
            case (mtftimes)
              klr%rbody(1:m2)=klr2%rbody(1:m2)*v1
c              nseg=(m2-1)/lseg+1
c              do i=0,nseg-1
c                klr%rbody(i*lseg+1:min((i+1)*lseg,m2))=
c     $               kl2%rbody(i*lseg+1:min((i+1)*lseg,m2))*v1
c              enddo
c              write(*,*)'tfecmplx-rl*2 '
            case (mtfrevpower)
              ir=int8(v1)
              klr%rbody(1:m2)=merge(merge(1.d0/klr2%rbody(1:m2),
     $             klr2%rbody(1:m2)**ir,ir .eq. -1),
     $             klr2%rbody(1:m2)**v1,dble(ir) .eq. v1)
            case (mtfpower)
              do concurrent (i=1:m2)
                ir=int8(klr2%rbody(i))
                klr%rbody(i)=merge(merge(1.d0/v1,v1**ir,ir .eq. -1),
     $               v1**klr2%rbody(i),dble(ir) .eq. klr2%rbody(i))
              enddo
            case (mtfgreater)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1>klr2%rbody(1:m2))
            case (mtfless)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1<klr2%rbody(1:m2))
            case (mtfgeq)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1>=klr2%rbody(1:m2))
            case (mtfleq)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1<=klr2%rbody(1:m2))
            case (mtfand)
              klr%rbody(1:m2)=merge(1.d0,0.d0,
     $             v1 .ne. 0.d0 .and. klr2%rbody(1:m2) .ne. 0.d0)
            case (mtfor)
              klr%rbody(1:m2)=merge(1.d0,0.d0,
     $             v1 .ne. 0.d0 .or. klr2%rbody(1:m2) .ne. 0.d0)
            case (mtfequal)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1 == klr2%rbody(1:m2))
            case (mtfunequal)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1 /= klr2%rbody(1:m2))
            case (mtfcomplex)
              d=.false.
              do i=1,m2
                if(klr2%dbody(i)%k .eq. 0)then
                  klr%dbody(i)=k1
                else
                  klr%dbody(i)=kxcalocv(0,v1,klr2%rbody(i))
                  d=.true.
                endif
              enddo
              if(d)then
                klr%attr=ior(klr%attr,lnonreallist)
              endif
            case default
              kx=tfeexpr(k1,k2,iopc1)
            end select
            return
          else
            c=.true.
            c1=v1
            kx=kxadaloc(-1,m2,klx)
            d=.false.
            select case(iopc1)
            case (mtfplus)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberq(k2i,cx2))then
                  cx=c1+cx2
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtftimes)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberq(k2i,cx2))then
                  cx=c1*cx2
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfrevpower)
              ix1=int8(v1)
              if(dble(ix1) .eq. v1)then
                if(ix1 .eq. -1)then
                  do i=1,m2
                    k2i=kl2%dbody(i)
                    if(tfnumberq(k2i,cx2))then
                      cx=1.d0/cx2
                      if(imag(cx) .eq. 0.d0)then
                        klx%rbody(i)=dble(cx)
                      else
                        klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                        d=.true.
                      endif
                    else
                      klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                    endif
                  enddo
                  go to 9000
                else
                  do i=1,m2
                    k2i=kl2%dbody(i)
                    if(tfnumberq(k2i,cx2))then
                      cx=cx2**ix1
                      if(imag(cx) .eq. 0.d0)then
                        klx%rbody(i)=dble(cx)
                      else
                        klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                        d=.true.
                      endif
                    else
                      klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                    endif
                  enddo
                  go to 9000
                endif
              else
                do i=1,m2
                  k2i=kl2%dbody(i)
                  if(tfnumberq(k2i,cx2))then
                    cx=cx2**cx1
                    if(imag(cx) .eq. 0.d0)then
                      klx%rbody(i)=dble(cx)
                    else
                      klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                      d=.true.
                    endif
                  else
                    klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                  endif
                enddo
                go to 9000
              endif
            case (mtfpower)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberq(k2i,cx2))then
                  if(imag(cx2) .eq. 0.d0)then
                    ix2=int8(dble(cx2))
                    cx=merge(merge(1.d0/v1,v1**ix2,ix2 .eq. -1),
     $                   v1**dble(cx2),ix2 .eq. dble(cx2))
                  else
                    cx=v1**cx2
                  endif
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfcomplex)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(tfnumberq(k2i,cx2))then
                  cx=dcmplx(v1-imag(cx2),dble(cx2))
                  if(imag(cx) .eq. 0.d0)then
                    klx%rbody(i)=dble(cx)
                  else
                    klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                    d=.true.
                  endif
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfgreater)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1>v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfless)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1<v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfgeq)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1>=v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfleq)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1<=v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfand)
              if(v1 .eq. 0.d0)then
                klx%rbody(1:m2)=0.d0
              else
                do i=1,m2
                  if(ktfrealq(k2i,v2i))then
                    klx%rbody(i)=merge(1.d0,0.d0,v2i .ne. 0.d0)
                  else
                    klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                  endif
                enddo
              endif
              go to 9000
            case (mtfor)
              if(v1 .ne. 0.d0)then
                klx%rbody(1:m2)=1.d0
              else
                do i=1,m2
                  k2i=kl2%dbody(i)
                  if(ktfrealq(k2i,v2i))then
                    klx%rbody(i)=merge(1.d0,0.d0,v2i .ne. 0.d0)
                  else
                    klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                  endif
                enddo
              endif
              go to 9000
            case (mtfequal)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1 == v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfunequal)
              do i=1,m2
                k2i=kl2%dbody(i)
                if(ktfrealq(k2i,v2i))then
                  klx%rbody(i)=merge(1.d0,0.d0,v1 /= v2i)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case default
              kx=tfeexpr(k1,k2,iopc1)
              return
            end select
          endif
        else
          kx=tfeexpr(k1,k2,iopc1)
        endif
      elseif(tfcomplexq(k1,c1))then
        if(tflistq(k2,kl2))then
          m2=kl2%nl
          kx=kxadaloc(-1,m2,klx)
          c=.true.
          d=.false.
          if(ktfreallistq(k2,klr2))then
            do i=1,m2
              v2i=klr2%rbody(i)
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
              if(tfnumberq(k2i,cx2))then
                cx=tfcmplxmathv(c1,cx2,iopc1)
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
              elseif(tflistq(k2i))then
                kxi=tfecmplxl(k1,k2i,iopc1)
                d=.true.
                c=c .and. tfconstq(kxi%k)
                klx%dbody(i)=dtfcopy(kxi)
              else
                kxi=tfeexpr(k1,k2i,iopc1)
                if(ktfnonrealq(kxi))then
                  c=c .and. tfconstq(kxi%k)
                  d=.true.
                  kxi=dtfcopy(kxi)                  
                endif
                klx%dbody(i)=kxi
              endif
            enddo
          endif
          go to 9000
        else
          kx=tfeexpr(k1,k2,iopc1)
          return
        endif
      elseif(tflistq(k1,kl1))then
        if(.not. tflistq(k2,kl2))then
          select case(iopc1)
          case(mtfpower,mtfrevpower)
            kx=tfecmplxl(k2,k1,mtfrevpower+mtfpower-iopc1)
          case(mtfcomplex)
            kxi=tfeval1(dfromr(0.d0),k2,mtfcomplex,irtc)
            if(irtc .ne. 0)then
              call tfreseterror
              kx=tfeexpr(k1,k2,mtfcomplex)
            else
              kx=tfecmplxl(kxi,k1,mtfplus)
            endif
          case(mtfgreater,mtfless,mtfgeq,mtfleq)
            kx=tfecmplxl(k2,k1,mtfgreater+mtfless-iopc1)
          case default
            kx=tfecmplxl(k2,k1,iopc1)
          end select
          return
        endif
        m1=kl1%nl
        if(m1 .ne. kl2%nl)then
          kx=tfeexpr(k1,k2,iopc1)
          return
        endif
        if(ktfreallistq(k1,klr1) .and. ktfreallistq(k2,klr2))then
          kx=kxavaloc(-1,m1,klr)
          klr%attr=lconstlist
          select case (iopc1)
          case (mtfplus)
            klr%rbody(1:m1)=klr2%rbody(1:m1)+klr1%rbody(1:m1)
          case (mtftimes)
            klr%rbody(1:m1)=klr2%rbody(1:m1)*klr1%rbody(1:m1)
          case (mtfrevpower)
            do concurrent (i=1:m1)
              ir=int8(klr1%rbody(i))
              klr%rbody(i)=merge(merge(1.d0/klr2%rbody(i),
     $             klr2%rbody(i)**ir,ir .eq. -1),
     $             klr2%rbody(i)**klr1%rbody(i),ir .eq. klr1%rbody(i))
            enddo
          case (mtfpower)
            do concurrent (i=1:m1)
              ir=int8(klr2%rbody(i))
              klr%rbody(i)=merge(merge(1.d0/klr1%rbody(i),
     $             klr1%rbody(i)**ir,ir .eq. -1),
     $             klr1%rbody(i)**klr2%rbody(i),
     $             ir .eq. klr2%rbody(i))
            enddo
          case (mtfcomplex)
            d=.false.
            do i=1,m1
              if(klr2%rbody(i) .eq. 0.d0)then
                klr%rbody(i)=klr1%rbody(i)
              else
                klr%dbody(i)=kxcalocv(0,klr1%rbody(i),klr2%rbody(i))
                d=.true.
              endif
            enddo
            if(d)then
              klr%attr=ior(klr%attr,lnonreallist+lconstlist)
            endif
          case (mtfgreater)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)>klr2%rbody(1:m1))
          case (mtfless)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)<klr2%rbody(1:m1))
          case (mtfgeq)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)>=klr2%rbody(1:m1))
          case (mtfleq)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)<=klr2%rbody(1:m1))
          case (mtfand)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)/=0.d0 .and. klr2%rbody(1:m1)/=0.d0)
          case (mtfor)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)/=0.d0 .or. klr2%rbody(1:m1)/=0.d0)
          case (mtfequal)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)==klr2%rbody(1:m1))
          case (mtfunequal)
            klr%rbody(1:m1)=merge(1.d0,0.d0,
     $           klr1%rbody(1:m1)/=klr2%rbody(1:m1))
          case default
            kx=tfeexpr(k1,k2,iopc1)
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
              if(tfnumberq(k1i,cx1) .and. tfnumberq(k2i,cx2))then
                cx=cx1+cx2
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
                cycle
              else
                kxi=merge(tfecmplxl(k1i,k2i,mtfplus),
     $                 tfeexpr(k1i,k2i,mtfplus),
     $                 tflistq(k1i) .or. tflistq(k2i))
              endif
              if(ktfnonrealq(kxi))then
                c=c .and. tfconstq(kxi%k)
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
              if(tfnumberq(k1i,cx1) .and. tfnumberq(k2i,cx2))then
                cx=cx1*cx2
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                  d=.true.
                endif
                cycle
              else
                kxi=merge(tfecmplxl(k1i,k2i,mtftimes),
     $               tfeexpr(k1i,k2i,mtftimes),
     $               tflistq(k1i) .or. tflistq(k2i))
              endif
              if(ktfnonrealq(kxi))then
                c=c .and. tfconstq(kxi%k)
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
              if(tfnumberq(k1i,cx1) .and. tfnumberq(k2i,cx2))then
                if(imag(cx1) .eq. 0.d0)then
                  ix1=int8(cx1)
                  cx=merge(merge(1.d0/cx2,cx2**ix1,ix1 .eq. -1),
     $                 cx2**dble(cx1),dble(ix1) .eq. dble(cx1))
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
              else
                kxi=merge(tfecmplxl(k1i,k2i,mtfrevpower),
     $               tfeexpr(k1i,k2i,mtfrevpower),
     $               tflistq(k1i) .or. tflistq(k2i))
              endif
              if(ktfnonrealq(kxi))then
                c=c .and. tfconstq(kxi%k)
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
              if(tfnumberq(k1i,cx1) .and. tfnumberq(k2i,cx2))then
                cx=cx1+dcmplx(-imag(cx2),dble(cx2))
                if(imag(cx) .eq. 0.d0)then
                  klx%rbody(i)=dble(cx)
                else
                  klx%dbody(i)=kxcalocv(0,dble(cx),imag(cx))
                endif
                cycle
              else
                kxi=merge(tfecmplxl(k1i,k2i,mtfcomplex),
     $               tfeexpr(k1i,k2i,mtfcomplex),
     $               tflistq(k1i) .or. tflistq(k2i))    
              endif
              if(ktfnonrealq(kxi))then
                c=c .and. tfconstq(kxi%k)
                d=.true.
                kxi=dtfcopy(kxi)                  
              endif
              klx%dbody(i)=kxi
            enddo
            go to 9000
          case (mtfgreater)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i>v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfless)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i<v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfgeq)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i>=v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfleq)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i<=v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfequal)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i==v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfunequal)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i/=v2i)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfand)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i/=0.d0 .and. v2i/=0.d0)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case (mtfor)
            kx=kxadaloc(-1,m1,klx)
            d=.false.
            do i=1,m1
              k1i=kl1%dbody(i)
              k2i=kl2%dbody(i)
              if(ktfrealq(k1i,v1i) .and. ktfrealq(k2i,v2i))then
                klx%rbody(i)=merge(1.d0,0.d0,v1i/=0.d0 .or. v2i/=0.d0)
              else
                klx%dbody(i)=tfecmplx1(k1i,k2i,iopc1,c,d)
              endif
            enddo
            go to 9000
          case default
            kx=tfeexpr(k1,k2,iopc1)
          end select
        endif
      else
        kx=tfeexpr(k1,k2,iopc1)
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
      end function

      recursive function tfcmplxf(k,mode,iaf) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: klr
      integer*4 ,intent(in):: mode,iaf
      integer*4 m,i,isp0
      if(tfnumberq(k))then
        if(ktfrealq(k))then
          if(mode .eq. 1)then
            kx=k
          elseif(mode .eq. 2)then
            kx%k=0
          else
            kx=k
          endif
        else
          call loc_sad(ktfaddrd(k),kl)
          if(mode .eq. 1)then
            kx=kl%dbody(1)
          elseif(mode .eq. 2)then
            kx=kl%dbody(2)
          else
            kx=kxcalocv(-1,kl%rbody(1),-kl%rbody(2))
          endif
        endif
        return
      elseif(tfreallistq(k,klr))then
        if(mode .eq. 1 .or. mode .eq. 3)then
          kx=k
        else
          kx%k=ktflist+ktraaloc(-1,klr%nl)
        endif
        return
      elseif(ktflistq(k,kl))then
        m=kl%nl
        if(mode .eq. 3 .or. kl%head%k .eq. ktfoper+mtflist)then
          isp0=isp
          do i=1,m
            isp=isp+1
            dtastk(isp)=tfcmplxf(kl%dbody(i),mode,iaf)
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
          return
        endif
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head%k=ktfoper+iaf
      klx%dbody(1)=dtfcopy(k)
      return
      end function

      subroutine tfcmplx(k1,k2,kx,iopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tfenum,tfeval1
      type (sad_complex), pointer :: cx1,cx2
      integer*8 ki1,ki2
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: iopc
      real*8 v1,v2
      complex*16 c1
      irtc=0
      if(iopc .eq. mtfnot)then
        if(ktfrealq(k2))then
          if(k2%k .eq. 0)then
            kx%k=ktftrue
          else
            kx%k=0
          endif
        else
          go to 8000
        endif
        return
      endif
      if(ktfrealq(k1,v1))then
        if(ktfrealq(k2,v2))then
          if(iopc .eq. mtfcomplex)then
            if(k2%k .ne. 0)then
              kx=kxcalocv(-1,v1,v2)
            else
              kx=k1
            endif
          elseif(iopc .eq. mtfpower)then
            if(k1%k .lt. 0)then
              ki2=int(v2)
              if(ki2 .ne. v2)then 
                c1=dcmplx(v1,0.d0)**v2
                if(imag(c1) .eq. 0.d0)then
                  kx=sad_descr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=dfromr(v1**ki2)
              endif
              return
            endif
            kx=tfenum(k1%x(1),k2%x(1),mtfpower,irtc)
          elseif(iopc .eq. mtfrevpower)then
            if(k2%k .lt. 0)then
              ki1=int(v1)
              if(ki1 .ne. v1)then 
                c1=dcmplx(v2,0.d0)**v1
                if(imag(c1) .eq. 0.d0)then
                  kx=sad_descr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=sad_descr(v2**ki1)
              endif
              return
            endif
            kx=tfenum(v2,v1,mtfpower,irtc)
          else
            kx=tfenum(k1%x(1),k2%x(1),iopc,irtc)
          endif
          return
         elseif(tfcomplexq(k2,cx2))then
          call tfcmplxmath(dcmplx(v1,0.d0),cx2%cx(1),kx,iopc,irtc)
          return
        else
          go to 8000
        endif
      elseif(tfcomplexq(k1,cx1))then
        if(ktfrealq(k2,v2))then
          call tfcmplxmath(cx1%cx(1),dcmplx(v2,0.d0),
     $         kx,iopc,irtc)
          return
        elseif(tfcomplexq(k2,cx2))then
          call tfcmplxmath(cx1%cx(1),cx2%cx(1),kx,iopc,irtc)
          return
        endif
      endif
 8000 kx=tfeval1(k1,k2,iopc,irtc)
      return
      end

      type (sad_descriptor) function tfenum(v1,v,iopc1,irtc)
      use tfstk
      implicit none
      real*8 ,intent(in):: v1,v
      real*8 x1
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      integer*4 ix,itfmessage
      irtc=0
      select case (iopc1)
      case (mtfneg)
        x1=-v
      case (mtfinv)
        x1=1.d0/v
      case (mtfplus)
        x1=v1+v
      case (mtftimes)
        x1=v1*v
c        if(abs(v) .eq. dinfinity .or. abs(v1) .eq. dinfinity)then
c          x1=merge(dnotanumber,merge(dinfinity,-dinfinity,
c     $         v .gt. 0.d0 .and. v1 .gt. 0.d0 .or.
c     $         v .lt. 0.d0 .and. v1 .lt. 0.d0),
c     $         v .eq. 0.d0 .or. v1 .eq. 0.d0)
c        else
c          x1=v1*v
c        endif
      case (mtfpower)
        if(v .eq. -1.d0)then
          x1=merge(0.d0,1.d0/v1,abs(v1) .eq. dinfinity)
        elseif(v .eq. 2.d0)then
          x1=v1**2
        elseif(v .eq. .5d0)then
          x1=sqrt(v1)
        elseif(v .eq. 0.d0 .and. redmath%value%k .ne. 0)then
          x1=1.d0
        else
          ix=int(v)
          x1=merge(v1**ix,v1**v,ix .eq. v)
        endif
      case (mtfrevpower)
        if(v1 .eq. -1.d0)then
          x1=1.d0/v
        elseif(v1 .eq. 2.d0)then
          x1=v**2
        elseif(v1 .eq. .5d0)then
          x1=sqrt(v)
        elseif(v1 .eq. 0.d0 .and. redmath%value%k .ne. 0)then
          x1=1.d0
        else
          ix=int(v1)
          x1=merge(v**ix,v**v1,ix .eq. v1)
        endif
      case (mtfequal)
        x1=merge(1.d0,0.d0,v1 .eq. v)
      case (mtfunequal)
        x1=merge(1.d0,0.d0,v1 .ne. v)
      case (mtfgreater)
        x1=merge(1.d0,0.d0,v1 > v)
      case (mtfless)
        x1=merge(1.d0,0.d0,v1 < v)
      case (mtfgeq)
        x1=merge(1.d0,0.d0,v1 >= v)
      case (mtfleq)
        x1=merge(1.d0,0.d0,v1 <= v)
      case (mtfnot)
        x1=merge(1.d0,0.d0,v .eq. 0.d0)
      case (mtfand)
        x1=merge(1.d0,0.d0,v1 .ne. 0.d0 .and. v .ne. 0.d0)
      case (mtfor)
        x1=merge(1.d0,0.d0,v1 .ne. 0.d0 .or. v .ne. 0.d0)
      case default
        write(*,*)'invop: ',iopc1
        irtc=itfmessage(999,'General::invop',' ')
        x1=0.d0
      end select
      tfenum%x(1)=x1
      return
      end

      subroutine tfjoine(k1,k2,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) k10,k20,ky1
      type (sad_dlist), pointer ::kl1,kl2
      integer*4 ,intent(out):: irtc
      integer*4 ma1,ma2,m,iopc,isp1
      if(ktfnonlistq(k1,kl1) .or. ktfnonlistq(k2,kl2))then
        irtc=-1
        return
      endif
      if( .not. tfsameheadq(k1,k2))then
        irtc=-1
        return
      endif
      irtc=0
      iopc=int(ktfaddr(kl1%head%k))
      if(iopc .eq. mtfplus .or. iopc .eq. mtftimes)then
        if(.not. tfnumberq(kl1%dbody(1)))then
          call tfjoin2(k2,k1,kx,.false.,irtc)
          go to 1000
        endif
        if(.not. tfnumberq(kl2%dbody(1)))then
          call tfjoin2(k1,k2,kx,.false.,irtc)
          go to 1000
        endif
        ma1=kl1%nl
        ma2=kl2%nl
        m=ma1+ma2-1
        k10=kl1%dbody(1)
        k20=kl2%dbody(1)
        call tfcmplx(k10,k20,ky1,iopc,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(ky1))then
          if(iopc .eq. mtfplus)then
            if(ky1%k .eq. 0)then
              m=m-1
            endif
          else
            if(ky1%k .eq. 0)then
              kx%k=0
              return
            elseif(ky1%k .eq. ktftrue)then
              m=m-1
            endif
          endif
        endif
        if(m .eq. 1)then
          kx=kl2%dbody(2)
          return
        endif
        isp1=isp
        if(m .eq. ma1+ma2-2)then
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        else
          isp=isp+1
          dtastk(isp)=ky1
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        endif
        kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),k_descr(ktfoper+iopc))
        isp=isp1
      else
        call tfjoin2(k1,k2,kx,.false.,irtc)
        return
      endif
 1000 isp=isp+1
      dtastk(isp)=kx
      call tfsort(isp-1,kx,0,irtc)
      isp=isp-1
      return
      end

      subroutine tfjoin2(k1,k2,kx,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: irtc
      integer*4 isp1
      logical*4 ,intent(in):: eval
      isp1=isp
      isp=isp+1
      dtastk(isp)=k1
      isp=isp+1
      dtastk(isp)=k2
      call tfjoin(isp1,kx,eval,irtc)
      isp=isp1
      return
      end

      subroutine tfjoin(isp1,kx,eval,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kf
      type (sad_dlist), pointer :: kl1,kli
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i,narg,isp0
      logical*4 ,intent(in):: eval
      logical*4 ev
      narg=isp-isp1
      if(ktfnonlistq(ktastk(isp1+1),kl1))then
        go to 9010
      endif
      if(narg .le. 1)then
        if(narg .eq. 1)then
          kx=dtastk(isp1+1)
          irtc=0
        else
          irtc=itfmessage(9,'General::narg','"1 or more"')
        endif
        return
      endif
      kf=kl1%head
      isp=isp+1
      isp0=isp
      call tfgetllstkall(kl1)
      if(isp .ge. mstk)then
        irtc=itfmessage(9999,'General::stack','"Join"')
        return
      endif
      do i=isp1+2,isp0-1
        if(ktfnonlistq(ktastk(i),kli))then
          isp=isp0-1
          go to 9000
        endif
        if(.not. tfsameq(kli%head,kf))then
          go to 9100
        endif
        call tfgetllstkall(kli)
        if(isp .ge. mstk)then
          irtc=itfmessage(9999,'General::stack','"Join"')
          return
        endif
      enddo
      dtastk(isp0)=kf
      ev=eval
      if(ev)then
        if(kf%k .eq. ktfoper+mtflist .or.
     $       kf%k .eq. ktfoper+mtfalt .or.
     $       kf%k .eq. ktfoper+mtfnull)then
          ev=.false.
        endif
      endif
      if(ev)then
        kx=tfefunref(isp0,.true.,irtc)
      else
        kx=kxcompose(isp0)
        irtc=0
      endif
      isp=isp0-1
      return
 9000 isp=isp0-1
 9010 irtc=itfmessage(9,'General::wrongtype',
     $     '"List or composition for all args"')
      return
 9100 irtc=itfmessage(9,'General::samehead',' ')
      isp=isp0-1
      return
      end
