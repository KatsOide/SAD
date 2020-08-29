      subroutine tfearray(k1,k,kx,iopc1,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) ,intent(in):: k,k1
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ky,tfdot,tfecmplxl
      type (sad_dlist), pointer :: kl,kl1
      integer*4 ,intent(out):: irtc
      integer*4 ne,ne1,i,iopc1,isp0
      logical*4 list1,list
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
c      call tfdebugprint(k1,'earray',1)
c      call tfdebugprint(k,'and',1)
c      write(*,*)'with ',iopc1
      irtc=0
      select case(iopc1)
      case (mtfplus:mtfunequal,mtfand,mtfor,mtfcomplex)
        kx=tfecmplxl(k1,k,iopc1)
        return
      case (mtfnot:mtfsame)
        go to 101
      end select
      if(ktflistq(k1,kl1))then
        if(tfcomplexq(k1))then
          ne1=0
          list1=.false.
        else
          if(tfexprq(k1))then
            go to 101
          endif
          ne1=kl1%nl
          list1=.true.
        endif
      else
        ne1=0
        list1=.false.
      endif
      if(ktflistq(k,kl))then
        if(tfcomplexq(k))then
          list=.false.
          ne=0
        else
          if(tfexprq(k))then
            go to 101
          endif
          if(iopc1 .eq. mtfdot)then
            kx=tfdot(k1,k,irtc)
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
              call tfcmplx(kl1%dbody(i),kl%dbody(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky%k .eq. 0)then
                kx%k=0
                return
              elseif(.not. ktfrealq(ky))then
                go to 101
              endif
            enddo
            return
          elseif(iopc1 .eq. mtfunequal)then
            kx%k=0
            do i=1,ne
              call tfcmplx(kl1%dbody(i),kl%dbody(i),ky,iopc1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ky%k .eq. ktftrue)then
                kx%k=ktftrue
                return
              elseif(.not. ktfrealq(ky))then
                go to 101
              endif
            enddo
            return
          else
            isp0=isp
            do i=1,ne
              isp=isp+1
              call tfcmplx(kl1%dbody(i),kl%dbody(i),dtastk(isp),
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
            call tfcmplx(kl1%dbody(i),k,dtastk(isp),iopc1,irtc)
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
 101  kx= tfeexpr(k1,k,iopc1)
      return
      end

      subroutine tfcmplxmath(c1,c2,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      complex*16 ,intent(in):: c1,c2
      complex*16 cx,tfcmplxmathv
      if(iopc1 .gt. mtfnot .and. iopc1 .ne. mtfcomplex)then
        irtc=-1
      else
        cx=tfcmplxmathv(c1,c2,iopc1)
        kx=merge(dfromr(dble(cx)),
     $       kxcalocv(-1,dble(cx),imag(cx)),
     $       imag(cx) .eq. 0.d0)
        irtc=0
      endif
      return
      end

      complex*16 function tfcmplxmathv(c1,c2,iopc1)
      use tfstk
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
        if(imag(c1) .eq. 0.d0)then
          i1=int8(c1)
          tfcmplxmathv=merge(merge(1.d0/c2,c2**i1,i1 .eq. -1),
     $         c2**dble(c1),i1 .eq. dble(c1))
        else
          tfcmplxmathv=c2**c1
        endif
      case(mtfpower)
        if(imag(c2) .eq. 0.d0)then
          i2=int8(c2)
          tfcmplxmathv=merge(merge(1.d0/c1,merge((1.d0,0.d0),c1**i2,
     $         i2 .eq. 0 .and. redmath%value%k .ne. 0),
     $         i2 .eq. -1),c1**dble(c2),i2 .eq. dble(c2))
        else
          tfcmplxmathv=c1**c2
        endif
      case (mtfequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 .eq. c2)
      case (mtfunequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 .ne. c2)
      case (mtfcomplex)
        tfcmplxmathv=c1+dcmplx(-imag(c2),dble(c2))
      case default
        tfcmplxmathv=0.d0
      end select
      return
      end
