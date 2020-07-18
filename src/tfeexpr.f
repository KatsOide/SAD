      subroutine tfcmplx(k1,k2,kx,iopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,tfenum,tfeval1
      type (sad_complex), pointer :: cx1,cx2
      integer*8 ki1,ki2
      integer*4 irtc,iopc
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
      real*8 v1,v,x1
      integer*4 ix,iopc1,irtc,itfmessage
      irtc=0
c      go to (
c     $     2000,
c     $     1010,1020,1030,2000,1050,2000,1070,1080,1090,1100,
c     $     1110,1120,1130,1140,1090,1100,1170,1180,1190),iopc1+1
c          m    i    +    -    *    /    v    ^    e    n    
c          >    <    g    l    E    N    ~    &    o    c
c      go to 2000
      select case (iopc1)
      case (mtfneg)
        x1=-v
      case (mtfinv)
        x1=1.d0/v
      case (mtfplus)
        x1=v1+v
      case (mtftimes)
        if(abs(v) .eq. dinfinity .or. abs(v1) .eq. dinfinity)then
          if(v .eq. 0.d0 .or. v1 .eq. 0.d0)then
            x1=dnotanumber
          elseif(v .gt. 0.d0 .and. v1 .gt. 0.d0 .or.
     $           v .lt. 0.d0 .and. v1 .lt. 0.d0)then
            x1=dinfinity
          else
            x1=-dinfinity
          endif
        else
          x1=v1*v
        endif
      case (mtfpower)
        if(v .eq. -1.d0)then
          if(abs(v1) .eq. dinfinity)then
            x1=0.d0
          else
            x1=1.d0/v1
          endif
        elseif(v .eq. 2.d0)then
          x1=v1**2
        elseif(v .eq. .5d0)then
          x1=sqrt(v1)
        elseif(v .eq. 0.d0 .and. redmath%value%k .ne. 0)then
          x1=1.d0
        else
          ix=int(v)
          if(ix .eq. v)then
            x1=v1**ix
          else
            x1=v1**v
          endif
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
          if(ix .eq. v1)then
            x1=v**ix
          else
            x1=v**v1
          endif
        endif
      case (mtfequal)
        if(v1 .eq. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfunequal)
        if(v1 .ne. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfgreater)
        if(v1 .gt. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfless)
        if(v1 .lt. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfgeq)
        if(v1 .ge. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfleq)
        if(v1 .le. v)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfnot)
        if(v .eq. 0.d0)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfand)
        if(v1 .ne. 0.d0 .and. v .ne. 0.d0)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case (mtfor)
        if(v1 .ne. 0.d0 .or. v .ne. 0.d0)then
          x1=1.d0
        else
          x1=0.d0
        endif
      case default
        irtc=itfmessage(999,'General::invop',' ')
        x1=0.d0
      end select
      tfenum%x(1)=x1
      return
      end

      subroutine tfjoine(k1,k2,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,k10,k20,ky1
      type (sad_dlist), pointer ::kl1,kl2
      integer*4 irtc,ma1,ma2,m,iopc,isp1
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
      type (sad_descriptor) k1,k2,kx
      integer*4 irtc,isp1
      logical*4 eval
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
      type (sad_descriptor) kx,kf
      type (sad_dlist), pointer :: kl1,kli
      integer*4 isp1,irtc,itfmessage,i,narg,isp0
      logical*4 eval,ev
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
