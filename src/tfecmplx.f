      module complex
      use tfstk
      type (sad_descriptor),external::tfeexpr

      contains
      function tfcmplx(k1,k2,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) tfeval1
      type (sad_complex), pointer :: cx1,cx2
      integer*8 ki1,ki2
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: iopc
      real*8 v1,v2
      complex*16 c1
      irtc=0
      if(iopc == mtfnot)then
        if(ktfrealq(k2))then
          if(k2%k == 0)then
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
          if(iopc == mtfcomplex)then
            kx=kxcalocv(-1,v1,v2)
          elseif(iopc == mtfpower)then
            if(k1%k .lt. 0)then
              ki2=int(v2)
              if(ki2 /= v2)then 
                c1=dcmplx(v1,0.d0)**v2
                kx=kxcalocc(-1,c1)
              else
                kx=dfromr(v1**ki2)
              endif
              return
            endif
            kx=tfenum(k1%x(1),k2%x(1),mtfpower,irtc)
          elseif(iopc == mtfrevpower)then
            if(k2%k .lt. 0)then
              ki1=int(v1)
              if(ki1 /= v1)then 
                c1=dcmplx(v2,0.d0)**v1
                kx=kxcalocc(-1,c1)
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
          kx=tfcmplxmath(dcmplx(v1,0.d0),cx2%cx(1),iopc,irtc)
          return
        else
          go to 8000
        endif
      elseif(tfcomplexq(k1,cx1))then
        if(ktfrealq(k2,v2))then
          kx=tfcmplxmath(cx1%cx(1),dcmplx(v2,0.d0),iopc,irtc)
          return
        elseif(tfcomplexq(k2,cx2))then
          kx=tfcmplxmath(cx1%cx(1),cx2%cx(1),iopc,irtc)
          return
        endif
      endif
 8000 kx=tfeval1(k1,k2,iopc,irtc)
      return
      end

      function tfcmplxmath(c1,c2,iopc1,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      complex*16 ,intent(in):: c1,c2
      complex*16 cx
      if(iopc1 .gt. mtfunequal .and. iopc1 /= mtfcomplex)then
        kx=dxnullo
        irtc=-1
      else
        cx=tfcmplxmathv(c1,c2,iopc1)
        kx=kxcalocc(-1,cx)
        irtc=0
      endif
      return
      end

      function tfecmplx1(k1,k2i,iopc1,c,d) result(kxi)
      implicit none
      type (sad_descriptor) kxi
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

      recursive function tfecmplxl(k1,k2,iopc1) result(kx)
c      use omp_lib
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
      real*8 v1,v2,v1i,v2i
      complex*16 c1,cx,cx1,cx2
      logical*4 d,c
c      call tfdebugprint(k1,'cmplxl',1)
c      call tfdebugprint(k2,'and:',1)
c      write(*,*)'with ',iopc1
      if(iopc1 == mtfnot)then
        if(ktfrealq(k2,v2))then
          kx%k=merge(ktftrue,ktffalse,v2 /= 0.d0)
        elseif(tfreallistq(k2,klr2) .and. klr2%nl .gt. 0)then
          kx=kxavaloc(-1,klr2%nl,klr)
          klr%rbody(1:klr2%nl)=merge(1.d0,0.d0,
     $             klr2%rbody(1:klr2%nl) == 0.d0)
        else
          kx=tfeexpr(k1,k2,iopc1)
        endif
        return
      endif
      if(ktfrealq(k1,v1))then
        if(tflistq(k2,kl2))then
          m2=kl2%nl
          if(m2 == 0)then
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
              if(dble(ir) == v1)then
                if(ir == -1)then
                  klr%rbody(1:m2)=1.d0/klr2%rbody(1:m2)
                else
                  klr%rbody(1:m2)=klr2%rbody(1:m2)**ir
                endif
              else
                klr%rbody(1:m2)=klr2%rbody(1:m2)**v1
              endif
c              klr%rbody(1:m2)=merge(merge(1.d0/klr2%rbody(1:m2),
c     $             klr2%rbody(1:m2)**ir,ir == -1),
c     $             klr2%rbody(1:m2)**v1,dble(ir) == v1)
            case (mtfpower)
              do concurrent (i=1:m2)
                ir=int8(klr2%rbody(i))
                klr%rbody(i)=merge(merge(1.d0/v1,v1**ir,ir == -1),
     $               v1**klr2%rbody(i),dble(ir) == klr2%rbody(i))
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
     $             v1 /= 0.d0 .and. klr2%rbody(1:m2) /= 0.d0)
            case (mtfor)
              klr%rbody(1:m2)=merge(1.d0,0.d0,
     $             v1 /= 0.d0 .or. klr2%rbody(1:m2) /= 0.d0)
            case (mtfequal)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1 == klr2%rbody(1:m2))
            case (mtfunequal)
              klr%rbody(1:m2)=merge(1.d0,0.d0,v1 /= klr2%rbody(1:m2))
            case (mtfcomplex)
              d=.false.
              do i=1,m2
                klr%dbody(i)=kxcalocv(0,v1,klr2%rbody(i),d)
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
                  klx%dbody(i)=kxcalocc(0,cx,d)
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
                  klx%dbody(i)=kxcalocc(0,cx,d)
                else
                  klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                endif
              enddo
              go to 9000
            case (mtfrevpower)
              ix1=int8(v1)
              if(dble(ix1) == v1)then
                if(ix1 == -1)then
                  do i=1,m2
                    k2i=kl2%dbody(i)
                    if(tfnumberq(k2i,cx2))then
                      cx=1.d0/cx2
                      klx%dbody(i)=kxcalocc(0,cx,d)
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
                      klx%dbody(i)=kxcalocc(0,cx,d)
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
                    klx%dbody(i)=kxcalocc(0,cx,d)
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
                  if(imag(cx2) == 0.d0)then
                    ix2=int8(dble(cx2))
                    cx=merge(merge(1.d0/v1,v1**ix2,ix2 == -1),
     $                   v1**dble(cx2),ix2 == dble(cx2))
                  else
                    cx=v1**cx2
                  endif
                  klx%dbody(i)=kxcalocc(0,cx,d)
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
                  klx%dbody(i)=kxcalocc(0,cx,d)
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
              if(v1 == 0.d0)then
                klx%rbody(1:m2)=0.d0
              else
                do i=1,m2
                  if(ktfrealq(k2i,v2i))then
                    klx%rbody(i)=merge(1.d0,0.d0,v2i /= 0.d0)
                  else
                    klx%dbody(i)=tfecmplx1(k1,k2i,iopc1,c,d)
                  endif
                enddo
              endif
              go to 9000
            case (mtfor)
              if(v1 /= 0.d0)then
                klx%rbody(1:m2)=1.d0
              else
                do i=1,m2
                  k2i=kl2%dbody(i)
                  if(ktfrealq(k2i,v2i))then
                    klx%rbody(i)=merge(1.d0,0.d0,v2i /= 0.d0)
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
              klx%dbody(i)=kxcalocc(0,cx,d)
            enddo
          else
            do i=1,m2
              k2i=kl2%dbody(i)
              if(tfnumberq(k2i,cx2))then
                cx=tfcmplxmathv(c1,cx2,iopc1)
                klx%dbody(i)=kxcalocc(0,cx,d)
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
            if(irtc /= 0)then
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
        if(m1 /= kl2%nl)then
          kx=tfeexpr(k1,k2,iopc1)
          return
        endif
        if(ktfreallistq(k1,klr1) .and. ktfreallistq(k2,klr2))then
          kx=kxavaloc(-1,m1,klr)
          klr%attr=lconstlist
          select case (iopc1)
          case (mtfplus)
c            call omp_set_dynamic(.true.)
c            call omp_set_num_threads(int(1+m1/400))
c!$OMP PARALLEL WORKSHARE shared(klr,klr1,klr2,m1)
            klr%rbody(1:m1)=klr2%rbody(1:m1)+klr1%rbody(1:m1)
c!$OMP END PARALLEL WORKSHARE
          case (mtftimes)
c            call omp_set_dynamic(.true.)
c            call omp_set_num_threads(int(1+m1/400))
c!$OMP PARALLEL WORKSHARE shared(klr,klr1,klr2,m1)
            klr%rbody(1:m1)=klr2%rbody(1:m1)*klr1%rbody(1:m1)
c!$OMP END PARALLEL WORKSHARE
          case (mtfrevpower)
            do concurrent (i=1:m1)
              ir=int8(klr1%rbody(i))
              if(ir == klr1%rbody(i))then
                if(ir == -1)then
                  klr%rbody(i)=1.d0/klr2%rbody(i)
                else
                  klr%rbody(i)=klr2%rbody(i)**ir
                endif
              else
                klr%rbody(i)=klr2%rbody(i)**klr1%rbody(i)
              endif
c              klr%rbody(i)=merge(merge(1.d0/klr2%rbody(i),
c     $             klr2%rbody(i)**ir,ir == -1),
c     $             klr2%rbody(i)**klr1%rbody(i),ir == klr1%rbody(i))
            enddo
          case (mtfpower)
            do concurrent (i=1:m1)
              ir=int8(klr2%rbody(i))
              if(ir == klr2%rbody(i))then
                if(ir == -1)then
                  klr%rbody(i)=1.d0/klr1%rbody(i)
                else
                  klr%rbody(i)=klr1%rbody(i)**ir
                endif
              else
                klr%rbody(i)=klr1%rbody(i)**klr2%rbody(i)
              endif
c              klr%rbody(i)=merge(merge(1.d0/klr1%rbody(i),
c     $             klr1%rbody(i)**ir,ir == -1),
c     $             klr1%rbody(i)**klr2%rbody(i),
c     $             ir == klr2%rbody(i))
            enddo
          case (mtfcomplex)
            d=.false.
            do i=1,m1
              klr%dbody(i)=kxcalocv(0,klr1%rbody(i),klr2%rbody(i),d)
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
                klx%dbody(i)=kxcalocc(0,cx,d)
                cycle
              elseif(tflistq(k1i) .or. tflistq(k2i))then
                kxi=tfecmplxl(k1i,k2i,mtfplus)
              else
                kxi=tfeexpr(k1i,k2i,mtfplus)
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
                klx%dbody(i)=kxcalocc(0,cx,d)
                cycle
              elseif(tflistq(k1i) .or. tflistq(k2i))then
                kxi=tfecmplxl(k1i,k2i,mtftimes)
              else
                kxi=tfeexpr(k1i,k2i,mtftimes)
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
                if(imag(cx1) == 0.d0)then
                  ix1=int8(cx1)
                  cx=merge(merge(1.d0/cx2,cx2**ix1,ix1 == -1),
     $                 cx2**dble(cx1),dble(ix1) == dble(cx1))
                else
                  cx=cx2**cx1
                endif
                klx%dbody(i)=kxcalocc(0,cx,d)
                cycle
              elseif(tflistq(k1i) .or. tflistq(k2i))then
                kxi=tfecmplxl(k1i,k2i,mtfrevpower)
              else
                kxi=tfeexpr(k1i,k2i,mtfrevpower)
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
                klx%dbody(i)=kxcalocc(0,cx,d)
                cycle
              elseif(tflistq(k1i) .or. tflistq(k2i))then
                kxi=tfecmplxl(k1i,k2i,mtfcomplex)
              else
                kxi=tfeexpr(k1i,k2i,mtfcomplex)
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

      complex*16 function tfcmplxmathv(c1,c2,iopc1)
      implicit none
      integer*4 ,intent(in):: iopc1
      integer*8 i1,i2
      complex*16 ,intent(in):: c1,c2
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
        if(imag(c1) == 0.d0 .or. imag(c1) == -0.d0)then
          i1=int8(c1)
          tfcmplxmathv=merge(merge(1.d0/c2,c2**i1,i1 == -1),
     $         c2**dble(c1),i1 == dble(c1))
        else
          tfcmplxmathv=c2**c1
        endif
      case(mtfpower)
        if(imag(c2) == 0.d0 .or. imag(c2) == -0.d0)then
          i2=int8(c2)
          tfcmplxmathv=merge(merge(1.d0/c1,merge((1.d0,0.d0),c1**i2,
     $         i2 == 0 .and. redmath%value%k /= 0),
     $         i2 == -1),c1**dble(c2),i2 == dble(c2))
        else
          tfcmplxmathv=c1**c2
        endif
      case (mtfequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 == c2)
      case (mtfunequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 /= c2)
      case (mtfcomplex)
        tfcmplxmathv=c1+dcmplx(-imag(c2),dble(c2))
      case default
        tfcmplxmathv=0.d0
      end select
      return
      end

      recursive function tfcmplxf(k,mode,iaf,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      type (sad_rlist), pointer :: klr
      integer*4 ,intent(in):: mode,iaf
      integer*4 ,intent(out):: irtc
      integer*4 m,i,isp0
      if(tfnumberq(k))then
        if(ktfrealq(k))then
          kx=merge(dxzero,k,mode == 2)
        else
          call loc_sad(ktfaddrd(k),kl)
          if(mode ==3)then
            kx=kxcalocv(-1,kl%rbody(1),-kl%rbody(2))
          else
            kx=kl%dbody(mode)
          endif
        endif
        irtc=0
        return
      elseif(tfreallistq(k,klr))then
        if(mode == 1 .or. mode == 3)then
          kx=k
        else
          kx=kxraaloc(-1,klr%nl)
        endif
        irtc=0
        return
      elseif(tflistq(k,kl))then
        m=kl%nl
        isp0=isp
        do i=1,m
          isp=isp+1
          dtastk(isp)=tfcmplxf(kl%dbody(i),mode,iaf,irtc)
          if(irtc /= 0)then
            dtastK(isp)=kl%dbody(i)
            irtc=0
          endif
        enddo
        kx=kxmakelist(isp0)
        isp=isp0
        irtc=0
        return
      endif
      irtc=-1
      return
      end function

      type (sad_descriptor) function tfenum(v1,v,iopc1,irtc)
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
c        if(abs(v) == dinfinity .or. abs(v1) == dinfinity)then
c          x1=merge(dnotanumber,merge(dinfinity,-dinfinity,
c     $         v .gt. 0.d0 .and. v1 .gt. 0.d0 .or.
c     $         v .lt. 0.d0 .and. v1 .lt. 0.d0),
c     $         v == 0.d0 .or. v1 == 0.d0)
c        else
c          x1=v1*v
c        endif
      case (mtfpower)
        if(v == -1.d0)then
          x1=merge(0.d0,1.d0/v1,abs(v1) == dinfinity)
        elseif(v == 2.d0)then
          x1=v1**2
        elseif(v == .5d0)then
          x1=sqrt(v1)
        elseif(v == 0.d0 .and. redmath%value%k /= 0)then
          x1=1.d0
        else
          ix=int(v)
          x1=merge(v1**ix,v1**v,ix == v)
        endif
      case (mtfrevpower)
        if(v1 == -1.d0)then
          x1=1.d0/v
        elseif(v1 == 2.d0)then
          x1=v**2
        elseif(v1 == .5d0)then
          x1=sqrt(v)
        elseif(v1 == 0.d0 .and. redmath%value%k /= 0)then
          x1=1.d0
        else
          ix=int(v1)
          x1=merge(v**ix,v**v1,ix == v1)
        endif
      case (mtfequal)
        x1=merge(1.d0,0.d0,v1 == v)
      case (mtfunequal)
        x1=merge(1.d0,0.d0,v1 /= v)
      case (mtfgreater)
        x1=merge(1.d0,0.d0,v1 > v)
      case (mtfless)
        x1=merge(1.d0,0.d0,v1 < v)
      case (mtfgeq)
        x1=merge(1.d0,0.d0,v1 >= v)
      case (mtfleq)
        x1=merge(1.d0,0.d0,v1 <= v)
      case (mtfnot)
        x1=merge(1.d0,0.d0,v == 0.d0)
      case (mtfand)
        x1=merge(1.d0,0.d0,v1 /= 0.d0 .and. v /= 0.d0)
      case (mtfor)
        x1=merge(1.d0,0.d0,v1 /= 0.d0 .or. v /= 0.d0)
      case default
        write(*,*)'invop: ',iopc1
        irtc=itfmessage(999,'General::invop',' ')
        x1=0.d0
      end select
      tfenum%x(1)=x1
      return
      end

      end module complex
