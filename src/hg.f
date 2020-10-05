      module hg

      contains
      real*8 recursive function hgrr(a,b,c,x) result(f)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      if(b .lt. a)then
        f=hgrr(b,a,c,x)
        return
      endif
      if(anint(a) .eq. a .and. a .le. 0.d0)then
        f=hgrp(a,b,c,x,.true.)
      elseif(x .le. -1.d0)then
        f=hgrr1(a,b,c,x)
      elseif(x .lt. 0.d0)then
        f=(1.d0-x)**(-a)*hgrr(a,c-b,c,x/(x-1.d0))
      elseif(x .eq. 0.d0)then
        f=1.d0/gamma(c)
      elseif(x .le. 0.5d0)then
        f=hgrr3(a,b,c,x)
      elseif(x .le. 1.d0)then
        f=hgrr4(a,b,c,x)
      else
        f=0.d0
      endif
      return
      end function

      real*8 function hgrp(a,b,c,x,reg) result(f)
      use gammaf
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 u,s,s1
      integer*4 i
      logical*4 ,intent(in):: reg
      if(reg)then
        u=1.d0/gamma(c)
      else
        u=1.d0
      endif
      f=u
      s=0.d0
      do i=1,nint(-a)
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        if(c+s .eq. 0.d0)then
          exit
        endif
        f=f+u
        s=s1
      enddo
      return
      end function

      real*8 function hgrr1(a,b,c,x) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 sab,x1,m,k,k1,f,u,lx,ba
      integer*4 i
      x1=1.d0/(1.d0-x)
      ba=b-a
      m=anint(ba)
      if(ba .ne. m)then
        sab=sin(m_pi*(ba))
        if(sab .ne. 0.d0)then
          f1=(x1**a*hgrr(a,c-b,1.d0-ba,x1)/gamma(b)/gamma(c-a)
     $         -x1**b*hgrr(b,c-a,ba+1.d0,x1)/gamma(a)/gamma(c-b)
     $         )/sab*m_pi
          return
        endif
      endif
      if(m .ne. 0.d0)then
        f=gamma(m)
        u=f
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(c-b+(i-1))/i/(m-i)*x1
          f=f+u
        enddo
        f=f/gamma(b)/gamma(c-a)
      else
        f=0.d0
      endif
      lx=-log(x1)
      u=(-x1)**m/gamma(m+1.d0)/gamma(a)/gamma(c-b)
      f=f+u*(lx+polygamma(1.d0)+polygamma(m+1.d0)
     $     -polygamma(b)-polygamma(c-a))
      k=1.d0
      do
        k1=k+1.d0
        u=u*(b+k-1.d0)*(c-a+k-1.d0)/k/(k+m)*x1
        f1=f+u*(lx+polygamma(k1)+polygamma(m+k1)
     $       -polygamma(b+k)-polygamma(c-a+k))
        if(f1 .eq. f)then
          f1=f1*x1**a
          return
        endif
        k=k1
        f=f1
      enddo
      return
      end function

      real*8 function hgrr2(a,b,c,x)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      hgrr2=(1.d0-x)**(-a)*hgrr(a,c-b,c,x/(x-1.d0))
      return
      end

      real*8 function hgr3(a,b,c,x) result(f1)
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 f,s,s1,u
      if(c .eq. 0.d0)then
        f1=1.d0/0.d0
        return
      endif
      u=1.d0
      s=0.d0
      f=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f+u
        if(f .eq. f1)then
          return
        endif
        f=f1
        s=s1
      enddo
      end function

      real*8 function hgrr3(a,b,c,x) result(f1)
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 f,s,s1,u
      if(anint(c) .eq. c .and. c .le. 0.d0)then
        s=-c+1.d0
        u=gamma(a+s)/gamma(a)*gamma(b+s)/gamma(b)
     $       /gamma(s+1.d0)*x**s
      else
        u=1.d0/gamma(c)
        s=0.d0
      endif
      f=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f+u
        if(f .eq. f1)then
          return
        endif
        f=f1
        s=s1
      enddo
      end function

      real*8 recursive function hgrr4(a,b,c,x) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 cab,sabc,x1,m,k,k1,lx,u,f
      integer*4 i
      x1=1.d0-x
      cab=c-a-b
      m=anint(cab)
      if(m .ne. cab)then
        sabc=sin(m_pi*(cab))
        if(sabc .ne. 0.d0)then
          f1=(hgrr(a,b,1.d0-cab,x1)
     $         /gamma(c-a)/gamma(c-b)
     $         -hgrr(c-a,c-b,cab+1.d0,x1)*x1**cab
     $         /gamma(a)/gamma(b))/sabc*m_pi
          return
        endif
      endif
      if(m .lt. 0.d0)then
        f1=x1**m*hgrr4(a+m,b+m,c,x)
        return
      endif
      if(m .ne. 0.d0)then
        u=gamma(c)
        f=u
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(b+(i-1))/i/(m-i)*x1
          f=f+u
        enddo
        f=f/gamma(a+m)/gamma(b+m)
      else
        f=0.d0
      endif
      lx=log(x1)
      u=-(-x1)**m/gamma(a)/gamma(b)/gamma(m+1.d0)
      f=f+u*(lx-polygamma(1.d0)-polygamma(m+1.d0)
     $     +polygamma(a+m)+polygamma(b+m))
      k=0.d0
      do
        k1=k+1.d0
        u=u*(a+m+k)*(b+m+k)/k1/(k1+m)*x1
        f1=f+u*(lx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $       +polygamma(a+m+k1)+polygamma(b+m+k1))
        if(f .eq. f1)then
          return
        endif
        k=k1
        f=f1
      enddo
      return
      end function

      complex*16 function chgrr5(a,b,c,x) result(cf1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 cab,sabc,x1,m,k,k1,u,f
      complex*16 clx,cf
      integer*4 i
      x1=1.d0-1.d0/x
      cab=c-a-b
      m=anint(cab)
      if(m .ne. cab)then
        sabc=sin(m_pi*(cab))
        if(sabc .ne. 0.d0)then
          cf1=(hgrr(a,a-c+1.d0,1.d0-cab,x1)
     $         /gamma(c-a)/gamma(c-b)*x**(-a)
     $         -hgrr(c-a,1.d0-a,cab+1.d0,x1)
     $         *dcmplx(1.d0-x,0.d0)**(cab)*x**(a-c)
     $         /gamma(a)/gamma(b))/sabc*m_pi
          return
        endif
      endif
      if(m .ne. 0.d0)then
        u=gamma(m)/gamma(b+m)
        f=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(b+m-(i-1))/i/(m-i)*x1
          f=f+u
        enddo
        f=f/gamma(a+m)
      else
        f=0.d0
      endif
      clx=log(dcmplx(-x1,0.d0))
      u=-(x1)**m/gamma(m+1.d0)/gamma(b)/gamma(a)
      cf=f+u*(clx-polygamma(1.d0)-polygamma(m+1.d0)
     $     +polygamma(a+m)+polygamma(b))
      k=0.d0
      do
        k1=k+1.d0
        u=-u*(a+m+k)/k1/(k1+m)*x1
        if(b .gt. k1 .or. b-k1 .ne. anint(b-k1))then
          cf1=cf+u*(clx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +polygamma(a+m+k1)+polygamma(b-k1))/gamma(b-k1)
        else
          cf1=cf+u*(-1.d0)**(k-b)*gamma(k1-b+1.d0)
        endif
        if(cf .eq. cf1)then
          cf1=cf1*x**(-a)
          return
        endif
        k=k1
        cf=cf1
      enddo
      return
      end function

      complex*16 function chgrr6(a,b,c,x) result(cf1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 ba,sab,x1,m,k,k1,u,f,d
      complex*16 clx,cf,cx
      integer*4 i
      x1=1.d0/x
      ba=b-a
      m=anint(ba)
      if(m .ne. ba)then
        sab=sin(m_pi*ba)
        if(sab .ne. 0.d0)then
          x1=1.d0/x
          cx=dcmplx(-x,0.d0)
          cf1=(cx**(-a)*hgrr(a,a-c+1.d0,1.d0-ba,x1)
     $         /gamma(b)/gamma(c-a)-
     $         cx**(-b)*hgrr(b,b-c+1.d0,ba+1.d0,x1)
     $         /gamma(a)/gamma(c-b))/sab*m_pi
          return
        endif
      endif
      if(m .ne. 0.d0)then
        u=gamma(m)/gamma(c-a)
        f=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(c-a-(i-1))/i/(m-i)*x1
          f=f+u
        enddo
        f=f/gamma(b)
      else
        f=0.d0
      endif
      clx=log(dcmplx(-x,0.d0))
      u=x1**(m)/gamma(a)/gamma(m+1.d0)
      d=c-b
      if(anint(d) .eq. d .and. d .le. 0.d0)then
        cf=f+u*(-1.d0)**(d)*gamma(1.d0-d)
      else
        cf=f+u*(clx+polygamma(1.d0)+polygamma(m+1.d0)
     $       -polygamma(b)-polygamma(d))/gamma(d)
      endif
      k=0.d0
      do
        k1=k+1.d0
        u=-u*(b+k)*x1/k1/(k1+m)
        d=c-b-k1
        if(anint(d) .eq. d .and. d .le. 0.d0)then
          cf1=cf+u*(-1.d0)**d*gamma(1.d0-d)
        else
          cf1=cf+u*(clx+polygamma(k1+1.d0)+polygamma(m+k1+1.d0)
     $         -polygamma(b+k1)-polygamma(d))/gamma(d)
        endif
        if(cf1 .eq. cf)then
          cf1=cf1*(-x)**(-a)
          return
        endif
        k=k1
        cf=cf1
      enddo
      return
      end function

      recursive function tfhg(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      real*8 a,b,c
      if(isp .ne. isp1+4)then
        go to 9000
      endif
      if(ktfnonrealq(dtastk(isp1+1),a))then
        go to 9000
      endif
      if(ktfnonrealq(dtastk(isp1+2),b))then
        go to 9000
      endif
      if(ktfnonrealq(dtastk(isp1+3),c))then
        go to 9000
      endif
      kx=kxhg(a,b,c,dtastk(isp),irtc)
      return
 9000 irtc=-1
      kx=dxnullo
      return
      end function

      recursive function kxhg(a,b,c,k,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist),pointer::kxl,kl
      type (sad_rlist),pointer::klv
      integer*4 ,intent(out):: irtc
      real*8 ,intent(in):: a,b,c
      real*8 x
      integer*4 i
      logical*4 d
      if(ktfrealq(k,x))then
        kx=kxhgr(a,b,c,x)
        irtc=0
      elseif(tfreallistq(k,klv))then
        kx=kxadaloc(-1,klv%nl,kxl)
        d=.false.
        do i=1,klv%nl
          kxl%dbody(i)=dtfcopyd(kxhgr(a,b,c,klv%rbody(i)),d)
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
        irtc=0
      elseif(tflistq(k,kl))then
        kx=kxadaloc(-1,kl%nl,kxl)
        d=.false.
        do i=1,kl%nl
          kxl%dbody(i)=dtfcopyd(kxhg(a,b,c,kl%dbody(i),irtc),d)
          if(irtc .ne. 0)then
            kxl%dbody(i:kl%nl)=dxnullo
            exit
          endif
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      else
        irtc=-1
      endif
      return
      end

      function kxhgr(a,b,c,x) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      real*8 ,intent(in):: a,b,c,x
      complex*16 cx
      if(x .eq. 0.d0)then
        kx%x(1)=1.d0
      elseif(c .eq. 0.d0)then
        kx%x(1)=1.d0/0.d0
      elseif(a .le. 0.d0 .and. anint(a) .eq. a)then
        kx%x(1)=hgrp(a,b,c,x,.false.)
      elseif(x .le. 1.d0)then
        kx%x(1)=hgrr(a,b,c,x)*gamma(c)
      elseif(x .le. 2.d0)then
        cx=chgrr5(a,b,c,x)*gamma(c)
        if(imag(cx) .ne. 0.d0)then
          kx=kxcalocv(-1,dble(cx),imag(cx))
        else
          kx%x(1)=dble(cx)
        endif
      else
        cx=chgrr6(a,b,c,x)*gamma(c)
        if(imag(cx) .ne. 0.d0)then
          kx=kxcalocv(-1,dble(cx),imag(cx))
        else
          kx%x(1)=dble(cx)
        endif
      endif
      return
      end function

      end module
