      real*8 function tfloor(x)
      implicit none
      real*8 x
      tfloor=aint(x+1.d-15)
      if(x .lt. 0.d0 .and. x .ne. tfloor)then
        tfloor=tfloor-1.d0
      endif
      return
      end

      complex*16 function tcfloor(z)
      implicit none
      real*8 tfloor
      complex*16 z
      tcfloor=dcmplx(tfloor(dble(z)),tfloor(imag(z)))
      return
      end

      real*8 function tceiling(x)
      implicit none
      real*8 x
      tceiling=-aint(-x+1.d-15)
      if(x .gt. 0.d0 .and. x .ne. tceiling)then
        tceiling=1.d0+tceiling
      endif
      return
      end

      complex*16 function tcceiling(z)
      implicit none
      real*8 tceiling
      complex*16 z
      tcceiling=dcmplx(tceiling(dble(z)),tceiling(imag(z)))
      return
      end

      real*8 function tround(x)
      implicit none
      real*8 x
      tround=sign(aint(abs(x)+0.5d0),x)
      return
      end

      complex*16 function tcround(z)
      implicit none
      real*8 tround
      complex*16 z
      tcround=dcmplx(tround(dble(z)),tround(imag(z)))
      return
      end

      real*8 function tfevenq(x)
      implicit none
      real*8 x
      if(x*.5d0 .eq. aint(x*.5d0))then
        tfevenq=1.d0
      else
        tfevenq=0.d0
      endif
      return
      end

      real*8 function tfoddq(x)
      implicit none
      real*8 x,tfevenq
      tfoddq=tfevenq(x+1.d0)
      return
      end

      complex*16 function tfcevenq(z)
      implicit none
      complex*16 z
      real*8 tfevenq
      tfcevenq=dcmplx(tfevenq(dble(z)),tfevenq(imag(z)))
      return
      end

      complex*16 function tfcoddq(z)
      implicit none
      complex*16 z
      real*8 tfoddq
      tfcoddq=dcmplx(tfoddq(dble(z)),tfoddq(imag(z)))
      return
      end

      complex*16 function tcsinh(z)
      implicit none
      complex*16 z,z1
      z1=dcmplx(-imag(z),dble(z))
      tcsinh=sin(z1)
      tcsinh=dcmplx(imag(tcsinh),-dble(tcsinh))
      return
      end

      complex*16 function tccosh(z)
      implicit none
      complex*16 z
      tccosh=cos(dcmplx(-imag(z),dble(z)))
      return
      end

      complex*16 function tctanh(z)
      implicit none
      complex*16 z,z1
      z1=dcmplx(-imag(z),dble(z))
      tctanh=sin(z1)/cos(z1)
      tctanh=dcmplx(imag(tctanh),-dble(tctanh))
      return
      end

      complex*16 function tcatan(z)
      implicit none
      complex*16 z
      if(dble(z) .eq. 0.d0)then
        tcatan=dcmplx(0.d0,atanh(dimag(z)))
      else
        tcatan=.5d0*dcmplx(atan2(dble(z),1.d0+imag(z))
     $       -atan2(-dble(z),1.d0-imag(z)),
     $       log(((1.d0+imag(z))**2+dble(z)**2)/
     $       ((1.d0-imag(z))**2+dble(z)**2))*.5d0)
      endif
      return
      end

      complex*16 function tcatan2(z1,z)
      implicit none
      real*8 pih
      parameter (pih=3.1415926535897932385d0*0.5d0)
      complex*16 z,z1,tcatan
      if(z .eq. (0.d0,0.d0))then
        if(dble(z1) .eq. 0.d0)then
          if(imag(z1) .eq. 0.d0)then
            tcatan2=(0.d0,0.d0)
          else
            tcatan2=sign(pih,imag(z1))
          endif
        else
          tcatan2=sign(pih,dble(z1))
        endif
      else
        tcatan2=tcatan(z1/z)
      endif
      return
      end

      complex*16 function tcasin(z)
      implicit none
      complex*16 z,z1
      z1=sqrt(1.d0-z**2)+(0.d0,1.d0)*z
      tcasin=dcmplx(atan2(imag(z1),dble(z1)),
     $     -.5d0*log(dble(z1)**2+imag(z1)**2))
      return
      end

      complex*16 function tcacos(z)
      implicit none
      complex*16 z,z1
      z1=sqrt(1.d0-z**2)+(0.d0,1.d0)*z
      tcacos=dcmplx(atan2(dble(z1),imag(z1)),
     $     .5d0*log(dble(z1)**2+imag(z1)**2))
      return
      end

      complex*16 function tcasinh(z)
      implicit none
      complex*16 z,z1,tcasin
      z1=dcmplx(-imag(z),dble(z))
      z1=tcasin(z1)
      tcasinh=dcmplx(imag(z1),-dble(z1))
      return
      end

      complex*16 function tcacosh(z)
      implicit none
      complex*16 z,z1,tcacos
      z1=tcacos(z)
      tcacosh=dcmplx(imag(z1),-dble(z1))
      return
      end

      complex*16 function tcatanh(z)
      implicit none
      complex*16 z,z1,tcatan
      z1=dcmplx(-imag(z),dble(z))
      z1=tcatan(z1)
      tcatanh=dcmplx(imag(z1),-dble(z1))
      return
      end

      complex*16 function tcxsin(x)
      implicit none
      complex*16 x,x2
      if(abs(x) .gt. .1d0)then
        tcxsin=x-sin(x)
      else
        x2=x**2
        tcxsin=x*x2/6.d0*(1.d0-x2/20.d0*(1.d0-x2/42.d0*(
     1       1.d0-x2/72.d0*(1.d0-x2/110.d0))))
      endif
      return
      end

      subroutine tfmod(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage,mode,i
      if(isp .le. isp1+1)then
        irtc=itfmessage(9,'General::narg','"2"')
      elseif(isp .ne. isp1+2)then
        if(mode .eq. 0)then
          irtc=itfmessage(9,'General::narg','"2"')
        else
          kx=dtastk(isp1+1)
          do i=isp1+2,isp
            call tfmodf(kx,dtastk(i),kx,mode,irtc)
            if(irtc .ne. 0)then
              return
            endif
          enddo
        endif
      else
        call tfmodf(dtastk(isp1+1),dtastk(isp),kx,mode,irtc)
      endif
      return
      end

      recursive subroutine tfmodf(k1,k2,kx,mode,irtc)
      use iso_c_binding
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      type (sad_list), pointer :: klx,kl1,kl2
      type (sad_rlist), pointer :: klr
      integer*8 ka1,ka2
      integer*4 irtc,i,n1,n2,itfmessage,isp0,isp2,mode
      real*8 v1,v2,vx,tfloor
      complex*16 c1,c2,cx,tcfloor
c     begin initialize for preventing compiler warning
      cx=0.d0
      vx=0.d0
c     end   initialize for preventing compiler warning
      ka1=ktfaddrd(k1)
      ka2=ktfaddrd(k2)
      if(tfcomplexnumlistqk(k1%k,kl1))then
        n1=kl1%nl
        if(tfcomplexnumlistqk(k2%k,kl2))then
          if(n1 .ne. kl2%nl)then
            irtc=itfmessage(9,'General::equalleng','"#1 and #2"')
            return
          endif
          isp0=isp
          if(ktfreallistq(kl1) .and. ktfreallistq(kl2))then
            kx=kxavaloc(-1,n1,klr)
            call descr_list(kx,klx)
            klr%attr=ior(klr%attr,lconstlist)
            if(mode .eq. 0)then
              do i=1,n1
                klr%rbody(i)=kl1%rbody(i)
     $               -tfloor(kl1%rbody(i)/kl2%rbody(i))*kl2%rbody(i)
              enddo
            elseif(mode .eq. 1)then
              do i=1,n1
                klr%rbody(i)=iand(int8(kl1%rbody(i)),int8(kl2%rbody(i)))
              enddo
            elseif(mode .eq. 2)then
              do i=1,n1
                klr%rbody(i)=ior(int8(kl1%rbody(i)),int8(kl2%rbody(i)))
              enddo
            elseif(mode .eq. 3)then
              do i=1,n1
                klr%rbody(i)=ieor(int8(kl1%rbody(i)),int8(kl2%rbody(i)))
              enddo
            endif
          else
            call tfgetllstkall(kl1)
            call tfgetllstkall(kl2)
            isp2=isp
            do i=1,n1
              isp=isp+1
              call tfmodf(dtastk(isp0+i),dtastk(isp0+n1+i),
     $             dtastk(isp),mode,irtc)
              if(irtc .ne. 0)then
                isp=isp0
                return
              endif
            enddo
            kx=kxmakelist(isp2)
            isp=isp0
          endif
        elseif(tfnumberqd(k2))then
          isp0=isp
          call tfgetllstkall(kl1)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfmodf(dtastk(isp0+i),k2,dtastk(isp),mode,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        else
          irtc=-1
          return
        endif
      elseif(tfcomplexnumlistqk(k2%k,kl2) .and. (tfcomplexq(k1) .or.
     $       tfnonlistq(k1)))then
        n2=kl2%nl
        isp0=isp
        call tfgetllstkall(kl2)
        isp2=isp
        do i=1,n2
          isp=isp+1
          call tfmodf(k1,dtastk(isp0+i),dtastk(isp),mode,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
        enddo
        kx=kxmakelist(isp2)
        isp=isp0
      elseif(ktfrealq(k1,v1))then
        if(ktfrealq(k2,v2))then
          if(mode .eq. 0)then
            if(k2%k .eq. 0)then
              irtc=itfmessage(9,'General::wrongval','"#2","nonzero"')
              return
            endif
            vx=v1-tfloor(v1/v2)*v2
          elseif(mode .eq. 1)then
            vx=iand(int8(v1),int8(v2))
          elseif(mode .eq. 2)then
            vx=ior(int8(v1),int8(v2))
          elseif(mode .eq. 3)then
            vx=ieor(int8(v1),int8(v2))
          endif
          kx=dfromr(vx)
        elseif(tfcomplexq(k2,c2))then
          if(mode .eq. 0)then
            cx=v1-tcfloor(v1/c2)*c2
          elseif(mode .eq. 1)then
            cx=dcmplx(dble(iand(int(v1),int(dble(c2)))),
     $           dble(iand(int(v1),int(imag(c2)))))
          elseif(mode .eq. 2)then
            cx=dcmplx(dble(ior(int(v1),int(dble(c2)))),
     $           dble(ior(int(v1),int(imag(c2)))))
          elseif(mode .eq. 3)then
            cx=dcmplx(dble(ieor(int(v1),int(dble(c2)))),
     $           dble(ieor(int(v1),int(imag(c2)))))
          endif
          go to 10
        else
          irtc=-1
          return
        endif
      elseif(tfcomplexq(k1,c1))then
        if(ktfrealq(k2,v2))then
          if(mode .eq. 0)then
            if(v2 .eq. 0.d0)then
              irtc=itfmessage(9,'General::wrongval','"#2","nonzero"')
              return
            endif
            cx=c1-tcfloor(c1/v2)*v2
          elseif(mode .eq. 1)then
            cx=dcmplx(dble(iand(int8(dble(c1)),int8(v2))),
     $           dble(iand(int8(imag(c1)),int8(v2))))
          elseif(mode .eq. 2)then
            cx=dcmplx(dble(ior(int8(dble(c1)),int8(v2))),
     $           dble(ior(int8(imag(c1)),int8(v2))))
          elseif(mode .eq. 3)then
            cx=dcmplx(dble(ieor(int8(dble(c1)),int8(v2))),
     $           dble(ieor(int8(imag(c1)),int8(v2))))
          endif
        elseif(tfcomplexq(k2,c2))then
          if(mode .eq. 0)then
            cx=c1-tcfloor(c1/c2)*c2
          elseif(mode .eq. 1)then
            cx=dcmplx(dble(iand(int8(dble(c1)),int8(dble(c2)))),
     $           dble(iand(int8(imag(c1)),int8(imag(c2)))))
          elseif(mode .eq. 2)then
            cx=dcmplx(dble(ior(int8(dble(c1)),int8(dble(c2)))),
     $           dble(ior(int8(imag(c1)),int8(imag(c2)))))
          elseif(mode .eq. 3)then
            cx=dcmplx(dble(ieor(int8(dble(c1)),int8(dble(c2)))),
     $           dble(ieor(int8(imag(c1)),int8(imag(c2)))))
          endif
        else
          irtc=-1
          return
        endif
        go to 10
      elseif(tflistq(k1%k,kl1))then
        n1=kl1%nl
        if(tflistq(k2%k,kl2))then
          if(n1 .ne. kl2%nl)then
            irtc=itfmessage(9,'General::equalleng','"#1 and #2"')
            return
          endif
          isp0=isp
          call tfgetllstkall(kl1)
          call tfgetllstkall(kl2)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfmodf(dtastk(isp0+i),dtastk(isp0+n1+i),
     $           dtastk(isp),mode,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        else
          isp0=isp
          call tfgetllstkall(kl1)
          isp2=isp
          do i=1,isp2-isp0
            isp=isp+1
            call tfmodf(dtastk(isp0+i),k2,dtastk(isp),mode,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        endif
      else
        irtc=-1
        return
      endif
      irtc=0
      return
 10   if(imag(cx) .ne. 0.d0)then
        kx=kxcalocv(-1,dble(cx),imag(cx))
      else
        kx=dfromr(dble(cx))
      endif
      irtc=0
      return
      end

      subroutine tfminmax(isp1,kx,mode,irtc)
      use tfstk
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: kli,klx
      type (sad_rlist), pointer :: klr,klir
      integer*4 isp1,irtc,mode,i
      real*8 xmin,xmax,ymin,ymax
      logical*4 cpx
      xmin=dinfinity
      xmax=-xmin
      cpx=.false.
      ymin=dinfinity
      ymax=-ymin
      if(isp .eq. isp1+1)then
        if(ktastk(isp) .eq. ktfoper+mtfnull)then
          go to 10
        endif
      endif
      do i=isp1+1,isp
        if(ktfrealq(ktastk(i)))then
          xmin=min(xmin,rtastk(i))
          xmax=max(xmax,rtastk(i))
        elseif(ktflistq(ktastk(i),kli))then
          if(kli%head%k .eq. ktfoper+mtfcomplex)then
            if(ktfnonreallistqo(kli) .or. kli%nl .ne. 2)then
              go to 9000
            endif
            call c_f_pointer(c_loc(kli),klir)
            cpx=.true.
            xmin=min(xmin,klir%rbody(1))
            xmax=max(xmax,klir%rbody(1))
            ymin=min(ymin,klir%rbody(2))
            ymax=max(ymax,klir%rbody(2))
          elseif(kli%head%k .eq. ktfoper+mtflist)then
            call tfminandmaxl(ksad_loc(kli%head%k),
     $           xmin,xmax,ymin,ymax,cpx,irtc)
            if(irtc .ne. 0)then
              return
            endif
          else
            go to 9000
          endif
        else
          go to 9000
        endif
      enddo
 10   if(mode .eq. 0)then
        if(cpx)then
          kx=kxadaloc(-1,2,klx)
          klx%dbody(1)=kxcalocv(0,xmin,ymin)
          klx%dbody(2)=kxcalocv(0,xmax,ymax)
        else
          kx=kxavaloc(-1,2,klr)
          call descr_list(kx,klx)
          klr%rbody(1)=xmin
          klr%rbody(2)=xmax
        endif
        klx%attr=ior(klx%attr,lconstlist)
      elseif(mode .eq. 1)then
        if(cpx)then
          kx=kxcalocv(-1,xmin,ymin)
        else
          kx=dfromr(xmin)
        endif
      else
        if(cpx)then
          kx=kxcalocv(-1,xmax,ymax)
        else
          kx=dfromr(xmax)
        endif
      endif
      irtc=0
      return
 9000 irtc=-1
      return
      end

      recursive subroutine tfminandmaxl(ka,xmin,xmax,ymin,ymax,cpx,irtc)
      use tfstk
      use iso_c_binding
      implicit none
      type (sad_descriptor) ki
      type (sad_list), pointer :: kl,kli
      type (sad_complex), pointer :: klic
      integer*8 ka,i
      integer*4 irtc,n
      real*8 xmin,xmax,ymin,ymax
      logical*4 cpx
      n=ilist(2,ka-1)
      if(n .ne. 0)then
        call loc_list(ka,kl)
        if(ktfreallistq(kl))then
          xmin=min(xmin,minval(kl%rbody(1:n)))
          xmax=max(xmax,maxval(kl%rbody(1:n)))
        else
          do i=1,n
            ki=kl%dbody(i)
            if(ktfrealq(ki))then
              xmin=min(xmin,kl%rbody(i))
              xmax=max(xmax,kl%rbody(i))
            elseif(ktflistq(ki,kli))then
              if(tfcomplexq(ki,klic))then
                cpx=.true.
                xmin=min(xmin,klic%re)
                xmax=max(xmax,klic%re)
                ymin=min(ymin,klic%im)
                ymax=max(ymax,klic%im)
              elseif(kli%head%k .eq. ktfoper+mtflist)then
                call tfminandmaxl(ksad_loc(kli%head%k),
     $               xmin,xmax,ymin,ymax,cpx,irtc)
                if(irtc .ne. 0)then
                  return
                endif
              else
                go to 9000
              endif
            else
              go to 9000
            endif
          enddo
        endif
      endif
      irtc=0
      return
 9000 irtc=-1
      return
      end

      subroutine tfrestrict(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"3"')
      elseif(ktfnonrealq(ktastk(isp-1)) .or.
     $       ktfnonrealq(ktastk(isp)))then
        irtc=-1
c        irtc=itfmessage(9,'General::wrongtype','"x, min, max"')
      else
        call tfrestrictl(dtastk(isp1+1),kx,
     $       rtastk(isp-1),rtastk(isp),irtc)
      endif
      return
      end

      recursive subroutine tfrestrictl(k,kx,x1,x2,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx
      type (sad_list), pointer :: kl,klx
      type (sad_rlist), pointer :: klr
      integer*4 irtc,i,isp0,n
      real*8 x1,x2,v
      irtc=0
      if(ktfrealq(k,v))then
        kx=dfromr(min(x2,max(x1,v)))
      elseif(tflistq(k,kl))then
        n=kl%nl
        if(n .eq. 0)then
          kx=k
          return
        endif
        if(ktfreallistq(kl))then
          kx=kxavaloc(-1,n,klr)
          call descr_list(kx,klx)
          klr%attr=ior(klr%attr,lconstlist)
          klr%rbody(1:n)=min(x2,max(x1,kl%rbody(1:n)))
c          do i=1,n
c            rlist(kax+i)=min(x2,max(x1,rlist(ka+i)))
c          enddo
        else
          isp0=isp
          do i=1,n
            isp=isp+1
            call tfrestrictl(kl%dbody(i),dtastk(isp),x1,x2,irtc)
            if(irtc .ne. 0)then
              return
            endif
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
        endif
      else
        irtc=-1
c        irtc=itfmessage(9,'General::wrongtype',
c     $       '"Real or List of Reals"')
      endif
      return
      end

      complex*16 function ccdabs(c)
      implicit none
      complex*16 c
      ccdabs=dcmplx(abs(c),0.d0)
      return
      end

      real*8 function tfsign(x)
      implicit none
      real*8 x
      if(x .gt. 0.d0)then
        tfsign=1.d0
      elseif(x .eq. 0.d0)then
        tfsign=0.d0
      else
        tfsign=-1.d0
      endif
      return
      end

      complex*16 function tfcsign(cx)
      implicit none
      complex*16 cx
      if(cx .eq. (0.d0,0.d0))then
        tfcsign=0.d0
      else
        tfcsign=cx/abs(cx)
      endif
      return
      end

      real*8 function tfarg(x)
      implicit none
      real*8 x
      if(x .ge. 0.d0)then
        tfarg=0.d0
      else
        tfarg=asin(1.d0)*2.d0
      endif
      return
      end

      complex*16 function tfcarg(z)
      implicit none
      complex*16 z
      tfcarg=imag(log(z))
      return
      end

      integer*4 function itfsyserr(level)
      implicit none
      character*132 string
      integer*4 level,itfmessage,i,l
      call gerror(string)
      do i=1,132
        if(string(i:i) .eq. char(0))then
          string(i:i)=' '
          l=i-1
          go to 10
        endif
      enddo
      l=len_trim(string)
 10   itfsyserr=itfmessage(level,'System::error',
     $     '"'//string(1:l)//'"')
      return
      end

      module cfunc
      implicit none
      integer*4 :: icrtc=0
      end module

      complex*16 function tfdummy(c)
      use cfunc
      implicit none
      complex*16 c
      icrtc=-1
      tfdummy=(0.d0,0.d0)
      return
      end

      recursive subroutine tfeintf(fun,cfun,k,kx,cmpl,rmin,rmax,ir)
      use tfstk
      use cfunc
      use tmacro
      implicit none
      type (sad_descriptor) k,kx
      type (sad_list), pointer ::klx,kl
      type (sad_rlist), pointer ::klr
      integer*4 ir,i,m,isp0
      real*8 fun,rmin,rmax
      complex*16 cfun,cv
      real*8 v
      logical*4 cmpl
      external fun,cfun
      ir=0
      if(ktfrealq(k,v))then
        if(v .lt. rmin .or. v .gt. rmax)then
          cv=cfun(dcmplx(v,0.d0))
          if(imag(cv) .ne. 0.d0)then
            kx=kxcalocv(-1,dble(cv),imag(cv))
          else
            kx=dfromr(dble(cv))
          endif
        else
          kx=dfromr(fun(v))
        endif
      elseif(tfnumberqd(k,cv))then
        if(cmpl)then
          cv=cfun(cv)
          if(imag(cv) .ne. 0.d0)then
            kx=kxcalocv(-1,dble(cv),imag(cv))
          else
            kx=dfromr(dble(cv))
          endif
          return
        else
          icrtc=0
          kx=dfromr(dble(cfun(cv)))
          ir=icrtc
        endif
      elseif(ktflistq(k,kl))then
        if(kl%head%k .eq. ktfoper+mtflist)then
          m=kl%nl
          if(ktfreallistq(kl))then
            isp0=isp
            if(rmin .ne. -dinfinity .or. rmax .ne. dinfinity)then
              do i=1,m
                isp=isp+1
                if(kl%rbody(i) .lt. rmin .or. kl%rbody(i) .gt. rmax)then
                  cv=cfun(dcmplx(kl%rbody(i),0.d0))
                  if(imag(cv) .eq. 0.d0)then
                    rtastk(isp)=dble(cv)
                  else
                    dtastk(isp)=kxcalocv(-1,dble(cv),imag(cv))
                  endif
                else
                  rtastk(isp)=fun(kl%rbody(i))
                endif
              enddo
              kx=kxmakelist(isp0)
              isp=isp0
            else
              kx=kxavaloc(-1,m,klr)
              call descr_list(kx,klx)
              klr%attr=ior(klr%attr,lconstlist)
              do i=1,m
                klr%rbody(i)=fun(kl%rbody(i))
              enddo
            endif
          else
            isp0=isp
            do i=1,m
              isp=isp+1
              if(ktfrealq(kl%body(i)))then
                if(kl%rbody(i) .lt. rmin .or. kl%rbody(i) .gt. rmax)then
                  cv=cfun(dcmplx(kl%rbody(i),0.d0))
                  if(imag(cv) .eq. 0.d0)then
                    rtastk(isp)=dble(cv)
                  else
                    dtastk(isp)=kxcalocv(-1,dble(cv),imag(cv))
                  endif
                else
                  rtastk(isp)=fun(kl%rbody(i))
                endif
              elseif(tfnumberqd(kl%dbody(i),cv))then
                if(cmpl)then
                  cv=cfun(cv)
                  if(imag(cv) .ne. 0.d0)then
                    dtastk(isp)=kxcalocv(-1,dble(cv),imag(cv))
                  else
                    rtastk(isp)=dble(cv)
                  endif
                else
                  rtastk(isp)=dble(cfun(cv))
                endif
              elseif(tflistq(kl%body(i)))then
                call tfeintf(fun,cfun,kl%dbody(i),
     $               dtastk(isp),cmpl,rmin,rmax,ir)
                if(ir .ne. 0.d0)then
                  isp=isp0
                  return
                endif
              else
                isp=isp0
                ir=-1
                return
              endif
            enddo
            kx=kxmakelist(isp0)
            isp=isp0
          endif
        else
          ir=-1
        endif
      else
        ir=-1
      endif
      return
      end

      recursive subroutine tfeintf2(fun,cfun,k,k1,cmpl,kx,ir)
      use tfstk
      use cfunc
      use tmacro
      implicit none
      type (sad_descriptor) kx,k,k1,ki,k1i
      type (sad_list), pointer ::klx,kl,kl1
      type (sad_rlist), pointer ::klr
      integer*4 ir,i,m,m1,isp0
      logical*4 cmpl
      external fun,cfun
      real*8 fun,v,v1
      complex*16 cv,cv1,cfun
      
      ir=0
      if(ktfrealq(k,v) .and. ktfrealq(k1,v1))then
        kx=dfromr(fun(v,v1))
      elseif(tfnumberqd(k,cv) .and. tfnumberqd(k1,cv1))then
        if(cmpl)then
          cv=cfun(cv,cv1)
          if(imag(cv) .ne. 0.d0)then
            kx=kxcalocv(-1,dble(cv),imag(cv))
          else
            kx=dfromr(dble(cv))
          endif
          return
        else
          icrtc=0
          kx=dfromr(dble(cfun(cv,cv1)))
          ir=icrtc
        endif
      elseif(tflistq(k,kl) .and. tflistq(k1,kl1))then
        m=kl%nl
        m1=kl1%nl
        if(m .ne. m1)then
          ir=8
          return
        endif
        if(iand(kl%attr,lnonreallist) .eq. 0 .and.
     $       iand(kl1%attr,lnonreallist) .eq. 0)then
          kx=kxavaloc(-1,m,klr)
          call descr_list(kx,klx)
          klr%attr=ior(klr%attr,lconstlist)
          do i=1,m
            klr%rbody(i)=fun(kl%rbody(i),kl1%rbody(i))
          enddo
        else
          isp0=isp
          do i=1,m
            ki=kl%dbody(i)
            k1i=kl1%dbody(i)
            if(ktfrealq(ki) .and.
     $           ktfrealq(k1i))then
              rtastk(isp)=fun(kl%rbody(i),kl1%rbody(i))
            elseif(tflistq(ki) .and. tflistq(k1i))then
              call tfeintf2(fun,cfun,ki,k1i,cmpl,dtastk(isp),ir)
              if(ir .ne. 0.d0)then
                return
              endif
            else
              ir=-1
              return
            endif
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
        endif
      else
        ir=-1
      endif
      return
      end

      subroutine tfmain(isp1,kx,irtc)
      use tfstk
      use trackbypass, only: bypasstrack
      use tfrbuf
      use tmacro
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,infl0,itfmessage,lfn,ierrfl,ierr,nc
      ierr=0
      if(isp .eq. isp1+2)then
        ierr=int(rtastk(isp))
      elseif(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
c      lfn=itfopenread(ktastk(isp1+1),.false.,irtc)
      call tfopenread(isp1,kx,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(.not. ktfrealq(kx,lfn))then
        irtc=itfmessage(9,'General::wonrgtype','"Real"')
        return
      endif
      infl0=infl
      infl=lfn
      ierrfl=errfl
      errfl=ierr
      bypasstrack=.true.
c      write(*,*)'tfmain-1 '
      call toplvl
c      write(*,*)'tfmain-2 '
      bypasstrack=.false.
      errfl=ierrfl
      call tfreadbuf(irbclose,lfn,int8(0),int8(0),nc,' ')
      infl=infl0
      kx%k=ktfoper+mtfnull
      return
      end

      recursive subroutine tfrgbcolor(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kxi
      type (sad_rlist), pointer ::kl
      integer*8 ka
      integer*4 isp1,irtc,ir,ig,ib,isp2,isp0,itfmessage,i,n
      character*2 r256(0:255)
      data r256/
     $     '00','01','02','03','04','05','06','07',
     $     '08','09','0A','0B','0C','0D','0E','0F',
     $     '10','11','12','13','14','15','16','17',
     $     '18','19','1A','1B','1C','1D','1E','1F',
     $     '20','21','22','23','24','25','26','27',
     $     '28','29','2A','2B','2C','2D','2E','2F',
     $     '30','31','32','33','34','35','36','37',
     $     '38','39','3A','3B','3C','3D','3E','3F',
     $     '40','41','42','43','44','45','46','47',
     $     '48','49','4A','4B','4C','4D','4E','4F',
     $     '50','51','52','53','54','55','56','57',
     $     '58','59','5A','5B','5C','5D','5E','5F',
     $     '60','61','62','63','64','65','66','67',
     $     '68','69','6A','6B','6C','6D','6E','6F',
     $     '70','71','72','73','74','75','76','77',
     $     '78','79','7A','7B','7C','7D','7E','7F',
     $     '80','81','82','83','84','85','86','87',
     $     '88','89','8A','8B','8C','8D','8E','8F',
     $     '90','91','92','93','94','95','96','97',
     $     '98','99','9A','9B','9C','9D','9E','9F',
     $     'A0','A1','A2','A3','A4','A5','A6','A7',
     $     'A8','A9','AA','AB','AC','AD','AE','AF',
     $     'B0','B1','B2','B3','B4','B5','B6','B7',
     $     'B8','B9','BA','BB','BC','BD','BE','BF',
     $     'C0','C1','C2','C3','C4','C5','C6','C7',
     $     'C8','C9','CA','CB','CC','CD','CE','CF',
     $     'D0','D1','D2','D3','D4','D5','D6','D7',
     $     'D8','D9','DA','DB','DC','DD','DE','DF',
     $     'E0','E1','E2','E3','E4','E5','E6','E7',
     $     'E8','E9','EA','EB','EC','ED','EE','EF',
     $     'F0','F1','F2','F3','F4','F5','F6','F7',
     $     'F8','F9','FA','FB','FC','FD','FE','FF'/
      if(isp .eq. isp1+1)then
        if(.not. tflistq(ktastk(isp)))then
          go to 9000
        endif
        ka=iand(ktamask,ktastk(isp))
        n=ilist(2,ka-1)
        if(iand(ilist(2,ka-3),lnonreallist) .eq. 0)then
          go to 9000
        endif
        isp0=isp
        do i=1,n
          if(.not. tfnumlistqn(dlist(ka+i),3,kl))then
            go to 9100
          endif
          isp2=isp
          call tfgetllstkall(kl)
          call tfrgbcolor(isp2,kxi,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
          isp=isp2+1
          dtastk(isp)=kxi
        enddo
        kx=kxmakelist(isp0)
        isp=isp0
      elseif(isp .ne. isp1+3)then
        go to 9000
      else
        if(iand(ktrmask,ktastk(isp1+1)) .eq. ktfnr)then
          go to 9000
        endif
        if(iand(ktrmask,ktastk(isp1+2)) .eq. ktfnr)then
          go to 9000
        endif
        if(iand(ktrmask,ktastk(isp)) .eq. ktfnr)then
          go to 9000
        endif
        ir=min(255,max(0,nint(256.d0*rtastk(isp1+1))))
        ig=min(255,max(0,nint(256.d0*rtastk(isp1+2))))
        ib=min(255,max(0,nint(256.d0*rtastk(isp))))
        kx=kxsalocb(-1,'#'//r256(ir)//r256(ig)//r256(ib),7)
        irtc=0
      endif
      return
 9100 isp=isp0
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"r,g,b or {{r1, g1, b1}, ...}"')
      return
      end
