      module hg
      use mathfun
      use macmath
      real*8 ,parameter ::sconf=2.d0**55,xmp=30.d0,xmth=0.96d0,
     $     xmth1=1.1d0,zt2=m_pi**2/6.d0,dz0=-.5d0*log(m_2pi),
     $     xmeps=1.d-6,
     $     dpmg1=-1.4132139976024971836d0,
     $     ddpmg1=-0.10668467023117844470d0,
     $     dpzmg1=1.1217084897192276354d0,
     $     dpmg2=-0.295789424393121511d0,
     $     ddpmg2=1.58338508513914499d0,
     $     dpzmg2=0.853223513334086317d0,
     $     epso=5e-17**2/4.d0,epsba=1.d-10,epsba1=1.d-7
      complex*16 ,parameter :: one=(1.d0,1.d-18),onec=(1.d0,-1.d-18)
      complex*16 ,parameter :: xs1=dcmplx(0.5d0,sqrt(.75d0)),
     $     xs2=dcmplx(0.5d0,-sqrt(.75d0))
      real*8 ,parameter :: epsx=0.3d0
      complex*16 cpgn,cpgz
      real*8 pgnc,pgn,pgz

      contains
      complex*16 recursive function chg(a,b,c,x,reg) result(f)
      use gammaf
      implicit none
      complex*16,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 ax,ax1
      if(dble(b) .lt. dble(a))then
        f=chg(b,a,c,x,reg)
        return
      elseif(imag(a) == 0.d0)then
        if(imag(b) == 0.d0 .and. imag(c) == 0.d0
     $       .and. imag(x) == 0.d0 .and. dble(x) <= 1.d0)then
          f=dcmplx(hgrr(dble(a),dble(b),dble(c),dble(x),reg),0.d0)
          return
        elseif(anint(dble(a)) == dble(a) .and. dble(a) <= 0.d0)then
          f=chgp(dble(a),b,c,x,reg)
          return
        endif
      endif
      if(abs(x-xs1) .lt. epsx .or. abs(x-xs2) .lt. epsx)then
        f=chg1s(a,b,c,x,reg)
      elseif(dble(x) <= 0.5d0)then
        ax=abs(x)
        if(ax .ge. 1.d0)then
          f=chg1(a,b,c,x,reg)
c          write(*,'(a,1p10g12.4)')'chg-1 ',a,b,c,x,f
        else
          ax1=abs(x-1.d0)
          if(ax1 > 1.d0)then
            f=zeroim(1.d0-x)**(-a)*chg(a,c-b,c,zeroim(x/(x-1.d0)),reg)
c            write(*,'(a,1p10g12.4)')'chg-2 ',a,b,c,x,f
          else
            f=chg3(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-3 ',a,b,c,x,f
          endif
        endif
      else
        ax1=abs(x-1.d0)
        if(ax1 > 1.d0)then
          f=chg6(a,b,c,x,reg)
c          write(*,'(a,1p10g12.4)')'chg-6 ',a,b,c,x,f
        else
          ax=abs(x)
          if(ax > 1.d0)then
            f=chg5(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-5 ',a,b,c,x,f
          else
            f=chg4(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-4 ',a,b,c,x,f
          endif
        endif
      endif
      return
      end function

      real*8 recursive function hgrr(a,b,c,x,reg) result(f)
      use gammaf
      use mathfun
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 ,parameter ::bth=2.d0**32
      if(b .lt. a)then
        f=hgrr(b,a,c,x,reg)
        return
      endif
      if(anint(a) == a .and. a <= 0.d0)then
        f=hgrp(a,b,c,x,reg)
      elseif(x <= -1.d0)then
        f=hgrr1(a,b,c,x,reg)
      elseif(x .lt. 0.d0)then
        if(b .lt. bth .or. x <= -0.5d0)then
          f=(1.d0-x)**(-a)*hgrr(a,c-b,c,x/(x-1.d0),reg)
        else
          f=hgrr3(a,b,c,x,reg)
        endif
      elseif(x == 0.d0)then
        if(reg)then
          f=gammai(c)
        else
          f=1.d0
        endif
      elseif(x <= 0.5d0)then
        f=hgrr3(a,b,c,x,reg)
      elseif(x <= 1.d0)then
        f=hgrr4(a,b,c,x,reg)
      else
        f=0.d0
      endif
      return
      end function

      complex*16 function chgp(a,b,c,x,reg) result(f)
      use gammaf
      use mathfun
      implicit none
      complex*16 ,intent(in):: b,c,x
      real*8 ,intent(in):: a
      logical*4 ,intent(in):: reg
      complex*16 g0,g1
      real*8 k
      integer*4 i
      if(reg)then
        f=cgammai(c)
      else
        f=(1.d0,0.d0)
      endif
      if(a == 0.d0)then
        return
      endif
      g0=f
      f=f-f*b/c*x
      if(a == -1.d0)then
        return
      endif
      g1=f
      k=-1.d0
      do i=2,-nint(a)
        f=(k*(x-1.d0)*g0-(c-2.d0*k+(k-b)*x)*g1)/(k-c)
        g0=g1
        g1=f
        k=k-1.d0
      enddo
      return
      end function

      real*8 function hgrp(a,b,c,x,reg) result(f)
      use gammaf
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 g0,g1,k
      integer*4 i
      logical*4 ,intent(in):: reg
      if(reg)then
        f=gammai(c)
      else
        f=1.d0
      endif
      if(a == 0.d0)then
        return
      endif
      g0=f
      f=f-f*b/c*x
      if(a == -1.d0)then
        return
      endif
      g1=f
      k=-1.d0
      do i=2,-nint(a)
        f=(k*(x-1.d0)*g0-(c-2.d0*k+(k-b)*x)*g1)/(k-c)
        g0=g1
        g1=f
        k=k-1.d0
      enddo
      return
      end function

      complex*16 function chg1s(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 f,u,x1,df,b2c,g,g0,g1
      real*8 k,k1
      integer*4 no
      u=zeroim(1.d0-.5d0*x)**(-a)
      if(reg)then
        g0=cgammai(c)
      else
        g0=1.d0
      endif
      g1=g0*(1.d0-2.d0*b/c)
      x1=x/(x-2.d0)
      f=u*g0
      u=u*a*x1
      f=f+u*g1
      k=1.d0
      b2c=b*2.d0-c
      no=nogam
      do
        k1=k+1.d0
        u=u*(a+k)/k1*x1
        g=(g0*k-b2c*g1)/(c+k)
        df=u*g
        f1=f+df
        no=no+10
c        write(*,'(a,i10,1p10g12.4)')'chg1s ',no,k1,u,f1,df
        if(abs(df)**2 <= no*abs(f1)**2*epso)then
          return
        endif
        f=f1
        k=k1
        g0=g1
        g1=g
      enddo
      end function

      complex*16 function chg1(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 sba,x1,df,u,lx,ba,b1,lc
      real*8 m,k,k1
      integer*4 i,no
      x1=zeroim(1.d0/(1.d0-x))
      ba=b-a
      m=anint(dble(ba))
      if(dble(ba) /= m .or. imag(ba) /= 0.d0)then
        sba=csinp(ba)/m_pi
        if(abs(sba) > epsba1)then
          lx=log(x1)
          if(reg)then
            f1=(chg(a,c-b,1.d0-ba,x1,.true.)
     $           *exp(a*lx-cloggamma(b)-cloggamma(c-a))
     $           -chg(b,c-a,ba+1.d0,x1,.true.)
     $           *exp(b*lx-cloggamma(a)-cloggamma(c-b))
     $           )/sba
          else
            lc=cloggamma(c)
            f1=(chg(a,c-b,1.d0-ba,x1,.true.)
     $           *exp(a*lx-cloggamma(b)-cloggamma(c-a)+lc)
     $           -chg(b,c-a,ba+1.d0,x1,.true.)
     $           *exp(b*lx-cloggamma(a)-cloggamma(c-b)+lc)
     $           )/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        f1=exp(log_gamma(m)-cloggamma(b1)-cloggamma(c-a))
        u=f1
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(c-b1+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nogam*3
      else
        f1=(0.d0,0.d0)
        no=0
      endif
      lx=-log(x1)
      u=exp(m*log(zeroim(-x1))-log_gamma(m+1.d0)
     $     -cloggamma(a)-cloggamma(c-b1))
c      u=zeroim(-x1)**m*gammai(m+1.d0)*cgammai(a)*cgammai(c-b1)
      f1=f1+u*(lx-m_euler+polygamma(m+1.d0)
     $     -cpolygamma(b1)-cpolygamma(c-a))
c      write(*,'(a,1p10g12.4)')'chg1 ',x,u,f1
      no=no+nogam*3+nopg*3
      k=1.d0
      do
        k1=k+1.d0
        u=u*(b1+k-1.d0)*(c-a+k-1.d0)/k/(k+m)*x1
        df=u*(lx+polygamma(k1)+polygamma(m+k1)
     $       -cpolygamma(b1+k)-cpolygamma(c-a+k))
        no=no+nopg*4
        f1=f1+df
        if(abs(df)**2 <= no*abs(f1)**2*epso)then
          f1=f1*x1**a
          f1=f1+(ba-m)*(cpolygamma(b1)*(1.d0-f1)
     $         +cpolygamma(b1+1.d0)*a*b1/c*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      real*8 function hgrr1(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 sba,x1,m,k,k1,df,u,lx,ba,b1
      integer*4 i,no
      x1=1.d0/(1.d0-x)
      ba=b-a
      m=anint(ba)
      if(ba /= m)then
        sba=sinp(ba)/m_pi
        if(abs(sba) > epsba)then
          if(reg)then
            f1=(hgrr(a,c-b,1.d0-ba,x1,.true.)
     $           *x1**a*gammai(b)*gammai(c-a)
     $           -hgrr(b,c-a,ba+1.d0,x1,.true.)
     $           *x1**b*gammai(a)*gammai(c-b))/sba
          else
            f1=(hgrr(a,c-b,1.d0-ba,x1,.true.)
     $           *x1**a*gammai(b)*pochh(c-a,a)
     $           -hgrr(b,c-a,ba+1.d0,x1,.true.)
     $           *x1**b*gammai(a)*pochh(c-b,b))/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        f1=pochh(b1,-a)*gammai(c-a)
c        f1=exp(log_gamma(m)-log_gamma(b1)-log_gamma(c-a))
        u=f1
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(c-b1+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*5+nogam
      else
        f1=(0.d0,0.d0)
        no=0.d0
      endif
      lx=-log(x1)
      u=(-x1)**nint(m)*gammai(m+1.d0)*gammai(a)*gammai(c-b1)
c      u=(-x1)**nint(m)*exp(-log_gamma(m+1.d0)-log_gamma(a)
c     $     -log_gamma(c-b1))
      f1=f1+u*(lx-m_euler+polygamma(m+1.d0)
     $     -polygamma(b1)-polygamma(c-a))
      no=no+nogam*3+nopg*3
      k=1.d0
      do
        k1=k+1.d0
        u=u*(b1+k-1.d0)*(c-a+k-1.d0)/k/(k+m)*x1
        df=u*(lx+polygamma(k1)+polygamma(m+k1)
     $       -polygamma(b1+k)-polygamma(c-a+k))
        f1=f1+df
        no=no+nopg*4
        if(df**2 <= no*f1**2*epso)then
          f1=f1*x1**a
          f1=f1+(ba-m)*polygamma(b1)*(1.d0-f1)
          if(.not. reg)then
            f1=f1*gamma(c)
          endif
          return
        endif
        k=k1
      enddo
      return
      end function

      complex*16 function chg3(a,b,c,x,reg) result(f1)
      use gammaf
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 u
      real*8 s,s1
      integer*4 no
      if(imag(c) == 0.d0 .and.
     $     anint(dble(c)) == dble(c) .and. dble(c) <= 0.d0)then
        s=-dble(c)+1.d0
        if(reg)then
          u=exp(cloggamma(a+s)-cloggamma(a)+cloggamma(b+s)-cloggamma(b)
     $         -log_gamma(s+1.d0)+s*log(zeroim(x))-cloggamma(c))
          no=nogam*6
        else
          u=exp(cloggamma(a+s)-cloggamma(a)+cloggamma(b+s)-cloggamma(b)
     $         -log_gamma(s+1.d0)+s*log(zeroim(x)))
          no=nogam*5
        endif
c        u=cpochh(a,dcmplx(s,0.d0))*cpochh(b,dcmplx(s,0.d0))
c     $       *gammai(s+1.d0)*zeroim(x)**s
      else
        if(reg)then
          u=cgammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
      endif
      f1=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f1+u
        no=no+9
        if(abs(u)**2 <= no*abs(f1)**2*epso)then
          exit
        endif
        s=s1
      enddo
      end function

      real*8 function hgrr3(a,b,c,x,reg) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 s,s1,u,x1
      real*8 ,parameter ::xthb=-20.d0
      integer*4 no
      x1=x*b
      if(anint(c) == c .and. c <= 0.d0)then
        s=-c+1.d0
        if(reg)then
          u=x**s*pochh(a,s)*pochh(b,s)*gammai(s+1.d0)*gammai(c)
c          u=exp(log_gamma(a+s)-log_gamma(a)+log_gamma(b+s)-log_gamma(b)
c     $         -log_gamma(s+1.d0)+s*log(x)-log_gamma(c))
          no=nogam*4+nolog
        else
          u=x**s*pochh(a,s)*pochh(b,s)*gammai(s+1.d0)
c          u=exp(log_gamma(a+s)-log_gamma(a)+log_gamma(b+s)-log_gamma(b)
c     $         -log_gamma(s+1.d0)+s*log(x))
          no=nogam*3+nolog
        endif
c        u=pochh(a,s)*pochh(b,s)*gammai(s+1.d0)*x**s
      else
        if(reg)then
          u=gammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
      endif
      f1=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f1+u
        no=no+4
        if(u**2 <= no*f1**2*epso)then
          return
        endif
        s=s1
      enddo
      end function

      complex*16 recursive function chg4(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 cab,scab,x1,lx,u,df,c1,lc
      real*8 m,k,k1
      integer*4 i,no
      x1=zeroim(1.d0-x)
      cab=c-a-b
      m=anint(dble(cab))
      if(m /= dble(cab) .or. imag(cab) /= 0.d0)then
        scab=csinp(cab)/m_pi
        if(abs(scab) > epsba1)then
          lx=log(x1)
          if(reg)then
            f1=(chg(a,b,1.d0-cab,x1,.true.)
     $           *exp(-cloggamma(c-a)-cloggamma(c-b))
     $           -chg(c-a,c-b,cab+1.d0,x1,.true.)
     $           *exp(lx*cab-cloggamma(a)-cloggamma(b)))/scab
          else
            lc=cloggamma(c)
            f1=(chg(a,b,1.d0-cab,x1,.true.)
     $           *exp(-cloggamma(c-a)-cloggamma(c-b)+lc)
     $           -chg(c-a,c-b,cab+1.d0,x1,.true.)
     $           *exp(lx*cab-cloggamma(a)-cloggamma(b)+lc))/scab
          endif
          return
        endif
      endif
      if(dble(m) .lt. 0.d0)then
        f1=x1**m*chg4(a+m,b+m,c,x,reg)
        return
      endif
      c1=a+b+m
      if(m /= 0.d0)then
c        u=cgammai(a+m)*cgammai(b+m)*gamma(m)
c        u=cgammai(a+m)*cpochh(b+m,-b)
        u=exp(-cloggamma(a+m)-cloggamma(b+m)+log_gamma(m))
        f1=u
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(b+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*3+nogam*3
      else
        f1=(0.d0,0.d0)
        no=0
      endif
      if(x1 /= (0.d0,0.d0))then
        lx=log(x1)
c        u=-zeroim(-x1)**m*cgammai(a)*cgammai(b)*gammai(m+1.d0)
        u=-exp(log(zeroim(-x1))*m-cloggamma(a)-cloggamma(b)
     $       -log_gamma(m+1.d0))
        f1=f1+u*(lx+m_euler-polygamma(m+1.d0)
     $       +cpolygamma(a+m)+cpolygamma(b+m))
        no=no+nogam*3+nopg*3
        k=0.d0
        do
          k1=k+1.d0
          u=u*(a+m+k)*(b+m+k)/k1/(k1+m)*x1
          df=u*(lx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +cpolygamma(a+m+k1)+cpolygamma(b+m+k1))
          f1=f1+df
          no=no+nopg*4
          if(abs(df)**2 <= no*abs(f1)**2*epso)then
            f1=f1-(cab-m)*(cpolygamma(c1)*(1.d0-f1)
     $           +cpolygamma(c1+1.d0)*a*b/c1*x)
            if(.not. reg)then
              f1=f1*cgamma(c)
            endif
            exit
          endif
          k=k1
        enddo
      endif
      return
      end function

      real*8 recursive function hgrr4(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 cab,scab,x1,m,k,k1,lx,u,df,c1
      integer*4 i,no
      x1=1.d0-x
      cab=c-a-b
      m=anint(cab)
      if(m /= cab)then
        scab=sinp(cab)/m_pi
        if(abs(scab) > epsba)then
          if(reg)then
            f1=(hgrr(a,b,1.d0-cab,x1,.true.)
     $           *gammai(c-a)*gammai(c-b)
     $           -hgrr(c-a,c-b,cab+1.d0,x1,.true.)
     $           *gammai(a)*gammai(b)*x1**cab)/scab
          else
            f1=(hgrr(a,b,1.d0-cab,x1,.true.)
     $           *gammai(c-a)*pochh(c-b,b)
     $           -hgrr(c-a,c-b,cab+1.d0,x1,.true.)
     $           *gammai(a)*pochh(b,c-b)*x1**cab)/scab
          endif
          return
        endif
      endif
      if(m .lt. 0.d0)then
        f1=x1**m*hgrr4(a+m,b+m,c,x,reg)
        return
      endif
      c1=a+b+m
      if(m /= 0.d0)then
        u=gammai(a+m)*pochh(b+m,-b)
        f1=u
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(b+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*3+nogam*3
      else
        f1=0.d0
        no=0
      endif
      if(x1 /= 0.d0)then
        lx=log(x1)
        u=-(-x1)**m*gammai(a)*gammai(b)*gammai(m+1.d0)
        f1=f1+u*(lx+m_euler-polygamma(m+1.d0)
     $       +polygamma(a+m)+polygamma(b+m))
        no=no+nogam*3+nopg*3
        k=0.d0
        do
          k1=k+1.d0
          u=u*(a+m+k)*(b+m+k)/k1/(k1+m)*x1
          df=u*(lx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +polygamma(a+m+k1)+polygamma(b+m+k1))
          no=no+nopg*4
          if(df**2 <= no*f1**2*epso)then
            f1=f1-(cab-m)*(polygamma(c1)*(1.d0-f1)
     $           +polygamma(c1+1.d0)*a*b/c1*x)
            if(.not. reg)then
              f1=f1*gamma(c)
            endif
            exit
          endif
          k=k1
        enddo
      endif
      return
      end function

      complex*16 function chg5(a,b,c,x,reg) result(f1)
      use gammaf
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 cab,scab,x1,u,df,c1
      complex*16 clx,lc,lx,lx1
      real*8 k,k1,m
      integer*4 i,no
      x1=1.d0-1.d0/x
      cab=c-a-b
      m=anint(dble(cab))
      if(m /= cab .or. imag(cab) /= 0.d0)then
        scab=csinp(cab)/m_pi
        if(abs(scab) >= epsba1)then
          lx=log(x)
          lx1=log(zeroim(1.d0-x))
          if(reg)then
            f1=( chg(a,a-c+1.d0,1.d0-cab,x1,.true.)
     $           *exp(-a*lx-cloggamma(c-a)-cloggamma(c-b))
     $           -chg(c-a,1.d0-a,1.d0+cab,x1,.true.)
     $           *exp((a-c)*lx+lx1*cab-cloggamma(a)-cloggamma(b)))/scab
          else
            lc=cloggamma(c)
            f1=( chg(a,a-c+1.d0,1.d0-cab,x1,.true.)
     $           *exp(-a*lx-cloggamma(c-a)-cloggamma(c-b)+lc)
     $           -chg(c-a,1.d0-a,1.d0+cab,x1,.true.)
     $           *exp((a-c)*lx+lx1*cab-cloggamma(a)-cloggamma(b)+lc)
     $           )/scab
          endif
          return
        endif
      endif
      if(m .lt. 0.d0)then
        f1=zeroim(1.d0-x)**m*chg5(a+m,b+m,c,x,reg)
        return
      endif
      c1=m+a+b
      if(m /= 0.d0)then
        u=cpochh(b+m,-b)*cgammai(a+m)
        f1=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(b+m-i)/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nogam*2
      else
        f1=(0.d0,0.d0)
        no=0
      endif
      clx=log(zeroim(-x1))
      u=-x1**m*gammai(m+1.d0)*cgammai(a)
      if(imag(b) == 0.d0 .and. dble(b) <= 0.d0 .and.
     $     dble(b) == anint(dble(b)))then
        f1=f1+u*(-1.d0)**b*cgamma(1.d0-b)
      else
        f1=f1+u*(clx+m_euler-polygamma(m+1.d0)
     $       +cpolygamma(a+m)+cpolygamma(b))*cgammai(b)
        no=no+nopg*4
      endif
      no=no+nogam*3+nopg*3
      k=0.d0
      do
        k1=k+1.d0
        u=-u*(a+m+k)/k1/(k1+m)*x1
        if(imag(b) == 0.d0 .and. dble(b) <= k1 .and.
     $       dble(b-k1) == anint(dble(b-k1)))then
          df=u*(-1.d0)**(k-b)*cgamma(k1-b+1.d0)
        else
          df=u*(clx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +cpolygamma(a+m+k1)+cpolygamma(b-k1))
     $         *cgammai(b-k1)
          no=no+nopg*4
        endif
        no=no+nogam
        f1=f1+df
        if(abs(df)**2 <= no*abs(f1)**2*epso)then
          f1=f1*zeroim(x)**(-a)
          f1=f1-(cab-m)*(cpolygamma(c1)*(1.d0-f1)
     $         +cpolygamma(c1+1.d0)*a*b/c1*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      complex*16 function chg6(a,b,c,x,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 ba,sba,x1,u,df,d,clx,b1,lc,lx
      real*8 k,k1,m
      integer*4 i,no
      x1=1.d0/x
      ba=b-a
      m=anint(dble(ba))
      if(imag(ba) /= 0.d0 .or. m /= dble(ba))then
        sba=csinp(ba)/m_pi
        if(abs(sba) > epsba1)then
          lx=log(zeroim(-x))
          if(reg)then
            f1=(chg(a,a-c+1.d0,1.d0-ba,x1,.true.)
     $           *exp(-a*lx-cloggamma(b)-cloggamma(c-a))
     $           -chg(b,b-c+1.d0,ba+1.d0,x1,.true.)
     $           *exp(-b*lx-cloggamma(a)-cloggamma(c-b)))/sba
          else
            lc=cloggamma(c)
            f1=(chg(a,a-c+1.d0,1.d0-ba,x1,.true.)
     $           *exp(-a*lx-cloggamma(b)-cloggamma(c-a)+lc)
     $           -chg(b,b-c+1.d0,ba+1.d0,x1,.true.)
     $           *exp(-b*lx-cloggamma(a)-cloggamma(c-b)+lc))/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        u=exp(log_gamma(m)-cloggamma(b1)-cloggamma(c-a))
c        u=cpochh(b1,m-b1)*cgammai(c-a)
        f1=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(c-a-i)/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nogam*2
      else
        f1=(0.d0,0.d0)
        no=0
      endif
      clx=log(zeroim(-x))
c      u=x**(-m)*cgammai(a)*gammai(m+1.d0)
      u=exp(-m*log(x)-cloggamma(a)-log_gamma(m+1.d0))
      d=c-b1
      if(imag(d) == 0.d0 .and.
     $     anint(dble(d)) == dble(d) .and. dble(d) <= 0.d0)then
        f1=f1+u*(-1.d0)**dble(d)*gamma(1.d0-dble(d))
        no=no+nogam*3
      else
        f1=f1+u*(clx-m_euler+polygamma(m+1.d0)
     $       -cpolygamma(b1)-cpolygamma(d))*cgammai(d)
        no=no+nogam*3+nopg*3
      endif
      k=0.d0
      do
        k1=k+1.d0
        d=c-b1-k1
        u=-u*(b1+k)*x1/k1/(k1+m)
        if(imag(d) == 0.d0 .and.
     $       anint(dble(d)) == dble(d) .and. dble(d) <= 0.d0)then
          df=u*(-1.d0)**dble(d)*gamma(1.d0-dble(d))
          no=no+nogam
        else
          df=u*(clx+polygamma(k1+1.d0)+polygamma(m+k1+1.d0)
     $         -cpolygamma(b1+k1)-cpolygamma(d))*cgammai(d)
          no=no+nogam+nopg*3
        endif
        f1=f1+df
        if(abs(df)**2 <= no*abs(f1)**2*epso)then
          f1=f1*zeroim(-x)**(-a)
          f1=f1+(ba-m)*(cpolygamma(b1)*(1.d0-f1)
     $         +cpolygamma(b1+1.d0)*a*b1/c*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      complex*16 recursive function confhg0(c,x0) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: c,x0
      complex*16 f,u,x
      real*8 s,s1
      integer*4 no
      x=zeroim(x0)
      if(imag(c) == 0.d0 .and. imag(x) == 0.d0)then
        f1=confhgrr0(dble(c),dble(x))
      else
        if(imag(c) == 0.d0 .and.
     $       anint(dble(c)) == dble(c) .and. dble(c) <= 0.d0)then
          s=-dble(c)+1.d0
          u=exp(log(x)*s-log_gamma(s+1.d0))
c          u=gammai(s+1.d0)*x**s
        else
          u=cgammai(c)
          s=0.d0
        endif
        no=nogam
        f=u
        do
          s1=s+1.d0
          u=u/(c+s)/s1*x
          f1=f+u
          no=no+5
          if(abs(u)**2 <= no*abs(f1)**2*epso)then
            exit
          endif
          f=f1
          s=s1
        enddo
      endif
      return
      end function

      real*8 recursive function confhgrr0(c,x) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: c,x
      real*8 f,s,s1,u
      integer*4 no
      if(anint(c) == c .and. c <= 0.d0)then
        s=-c+1.d0
c        u=gammai(s+1.d0)*x**s
        u=exp(s*log(x)-log_gamma(s+1.d0))
      else
        u=gammai(c)
        s=0.d0
      endif
      no=nogam
      f=u
      do
        s1=s+1.d0
        u=u/(c+s)/s1*x
        f1=f+u
        if(u**2 <= no*f1**2*epso)then
          return
        endif
        f=f1
        s=s1
      enddo
      return
      end function

      complex*16 recursive function confhg1(a,c,x0,reg) result(f1)
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: a,c,x0
      complex*16 u,x
      real*8 s,s1
      logical*4 ,intent(in) , optional:: reg
      logical*4 reg1
      integer*4 no
      if(present(reg))then
        reg1=reg
      else
        reg1=.true.
      endif
      x=zeroim(x0)
      if(imag(a) == 0.d0 .and. imag(c) == 0.d0 .and.
     $     imag(x) == 0.d0)then
        f1=confhgrr1(dble(a),dble(c),dble(x),reg1)
      elseif(imag(c) == 0.d0 .and.
     $       dble(c) == anint(dble(c)) .and. dble(c) <= 0.d0)then
        if(reg1)then
          f1=cpochh(a,1.d0-c)*x**(1.d0-c)
     $         *confhg1(a-c+1.d0,2.d0-c,x,.true.)
        else
          f1=dcmplx(1.0/0.d0,0.d0)
        endif
      elseif(dble(x) .ge. 0.d0)then
        if(reg1)then
          u=cgammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
        f1=u
        do
          s1=s+1.d0
          if(a+s == (0.d0,0.d0))then
            exit
          else
            u=u*(a+s)/(c+s)/s1*x
          endif
          f1=f1+u
          no=no+7
          if(abs(u)**2 <= no*abs(f1)**2*epso)then
            exit
          endif
          s=s1
        enddo
      else
        f1=exp(x)*confhg1(c-a,c,-x,reg1)
      endif
      return
      end function

      real*8 recursive function confhgrr1(a,c,x,reg) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,c,x
      real*8 s,s1,u
      integer*4 no
      logical*4 ,intent(in), optional:: reg
      logical*4 reg1
      if(present(reg))then
        reg1=reg
      else
        reg1=.true.
      endif
      if(c == anint(c) .and. c <= 0.d0)then
        if(reg1)then
          f1=pochh(a,1.d0-c)*x**(1.d0-c)
     $         *confhgrr1(a-c+1.d0,2.d0-c,x,.true.)
        else
          f1=1.d0/0.d0
        endif
      elseif(x .ge. 0.d0)then
        if(reg1)then
          u=gammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
c        write(*,'(a,1p10g12.4)')'chgr ',a,c,x,s,u
        f1=u
        do
          if(a+s == 0.d0)then
            exit
          else
            s1=s+1.d0
            u=u*(a+s)/(c+s)/s1*x
          endif
          f1=f1+u
          no=no+3
          if(u**2 <= no*f1**2*epso)then
            exit
          endif
          s=s1
        enddo
      else
        f1=exp(x)*confhgrr1(c-a,c,-x,reg1)
      endif
      return
      end function

      recursive function tfhg(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 iconf
      logical*4 reg
      reg=.false.
      iconf=0
      if(isp == isp1+5)then
        reg=.true.
        if(ktastk(isp1+3) == dxnullo%k)then
          iconf=2
        elseif(ktastk(isp1+4) == dxnullo%k)then
          iconf=1
        endif
      elseif(isp == isp1+2)then
        iconf=2
      elseif(isp == isp1+3)then
        iconf=1
      elseif(isp /= isp1+4)then
        go to 9000
      endif
      kx=kxhg(dtastk(isp1+1),dtastk(isp1+2),dtastk(isp1+3-iconf),
     $     dtastk(isp1+4-iconf),reg,iconf,irtc)
      return
 9000 irtc=-1
      kx=dxnullo
      return
      end function

      recursive function kxhg(ka,kb,kc,k,reg,iconf,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k,ka,kb,kc
      type (sad_dlist),pointer::kxl,kl
      type (sad_rlist),pointer::klv
      integer*4 ,intent(out):: irtc
      complex*16 ca,cb,cc,cx
      integer*4 ,intent(in):: iconf
      integer*4 i
      logical*4 ,intent(in):: reg
      logical*4 d
      irtc=0
      if(.not. tfnumberq(kc,cc) .or.
     $     (iconf /= 2 .and. .not. tfnumberq(ka,ca))
     $       .or. (iconf == 0 .and. .not. tfnumberq(kb,cb)))then
        kx=dxnullo
        irtc=-1
        return
      endif
      if(tfnumberq(k,cx))then
        kx=kxhgc(ca,cb,cc,cx,iconf,reg)
      elseif(tfreallistq(k,klv))then
        kx=kxadaloc(-1,klv%nl,kxl)
        d=.false.
        do i=1,klv%nl
          kxl%dbody(i)=dtfcopyd(kxhgc(ca,cb,cc,
     $         dcmplx(klv%rbody(i),0.d0),iconf,reg),d)
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      elseif(tflistq(k,kl))then
        kx=kxadaloc(-1,kl%nl,kxl)
        d=.false.
        do i=1,kl%nl
          kxl%dbody(i)=dtfcopyd(kxhg(ka,kb,kc,kl%dbody(i),
     $         reg,iconf,irtc),d)
          if(irtc /= 0)then
            kxl%dbody(i:kl%nl)=dxnullo
            exit
          endif
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      else
        kx=dxnullo
        irtc=-1
      endif
      return
      end function 

      recursive function kxhgc(a,b,c,x,iconf,reg) result(kx)
      use gammaf
      use tfstk
      implicit none
      type (sad_descriptor) kx
      complex*16 ,intent(in):: a,b,c,x
      integer*4 ,intent(in):: iconf
      complex*16 cx
      logical*4 ,intent(in):: reg
      select case (iconf)
      case (0)
        cx=chg(a,b,c,x,reg)
      case (1)
        cx=confhg1(a,c,x,reg)
      case (2)
        cx=confhg0(c,x)
        if(.not. reg)then
          cx=cx*cgamma(c)
        endif
      case (3)
        cx=chgu(a,c,x)
      case (4)
        cx=chlerch(x,a,c)
      case (5)
        cx=clerch(x,a,c)
      case default
        kx=dxnullo
        return
      end select
      kx=kxcalocc(-1,cx)
      return
      end function

      function kxhgpq(isp1,reg,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist),pointer ::kla,klb
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: reg
      complex*16 ,allocatable::a(:),b(:)
      integer*4 na,nb,m
      logical*4 cmplm,realm,veca,vecb
      irtc=-1
      kx=dxnullo
      call tfmatrixmaybeq(dtastk(isp1+1),cmplm,realm,veca,na,m,kla)
      if(m /= 0)then
        return
      endif
      call tfmatrixmaybeq(dtastk(isp1+2),cmplm,realm,vecb,nb,m,klb)
      if(m /= 0)then
        return
      endif
      allocate(a(0:na))
      call tfl2cm(kla,a(1:na),na,0,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      a(0)=merge((1.d0,0.d0),(0.d0,0.d0),veca)
      allocate(b(0:nb))
      call tfl2cm(klb,b(1:nb),nb,0,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      b(0)=merge((1.d0,0.d0),(0.d0,0.d0),vecb)
      kx=kxhgpqa(a,b,dtastk(isp),reg,irtc)
      return
      end function

      recursive function kxhgpqa(a,b,k,reg,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,pointer ::kl
      type (sad_dlist) ,pointer ::kxl
      complex*16 ,intent(in):: a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 ,intent(out):: irtc
      complex*16 z
      integer*4 i
      logical*4 d
      irtc=0
      if(tfnumberq(k,z))then
        kx=kxcalocc(-1,chgpq(a,b,z,reg))
      elseif(tflistq(k,kl))then
        kx=kxadaloc(-1,kl%nl,kxl)
        d=.false.
        do i=1,kl%nl
          kxl%dbody(i)=dtfcopyd(kxhgpqa(a,b,kl%dbody(i),reg,irtc),d)
          if(irtc /= 0)then
            kxl%dbody(i:kl%nl)=dxnullo
            return
          endif
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      endif
      return
      end function

      complex*16  function chgpq(a,b,z,reg) result(f)
      use gammaf
      implicit none
      complex*16 ,intent(in):: z,a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 na,nb
      complex*16 u
      real*8 k
      integer*4 i,no
      if(imag(z) == 0.d0 .and.
     $     a(0) /= (0.d0,0.d0) .and. b(0) /= (0.d0,0.d0))then
        f=hgpq(a,b,dble(z),reg)
        return
      endif
      na=size(a)-1
      nb=size(b)-1
c      write(*,'(a,2i5,1p10g12.4)')'chgpq ',na,nb,b(1)
      if(nb == 0 .and. na == 0)then
        f=exp(z)
      elseif(nb == 1)then
        select case (na)
        case (0)
          if(reg)then
            f=confhg0(b(1),z)
          else
            f=confhg0(b(1),z)*cgammai(b(1))
          endif
          return
        case (1)
          f=confhg1(a(1),b(1),z,reg)
          return
        case (2)
          f=chg(a(1),a(2),b(1),z,reg)
          return
        end select
      endif
      u=(0.d0,0.d0)
      if(reg)then
        do i=1,nb
          u=u+cloggamma(b(i))
        enddo
        no=nb*nogam
      endif
      u=exp(-u)
      f=u
      k=0.d0
      no=0
      main: do
        do i=1,na
          u=u*(a(i)+k)
          if(u == (0.d0,0.d0))then
c            write(*,'(a,3i5,1p10g12.4)')'hgpq-u0 ',i,na,nb,k,a(i),u
            return
          endif
        enddo
        do i=1,nb
          u=u/(b(i)+k)
        enddo
        k=k+1.d0
        u=u*z/k
        no=no+2*(na+nb+1)
        f=f+u
c        write(*,'(a,3i5,1p10g12.4)')'hgpq-0 ',i,na,nb,z,u,f
        if(abs(u)**2 <= no*abs(f)**2*epso)then
          return
        endif
      enddo main
      end function 

      complex*16  function hgpq(a,b,z,reg) result(f)
      use gammaf
      implicit none
      real*8 ,intent(in):: z
      complex*16 ,intent(in):: a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 na,nb
      real*8 u,k,f1
      integer*4 i,no
      na=size(a)-1
      nb=size(b)-1
      if(nb == 0 .and. na == 0)then
        f=exp(z)
      elseif(nb == 1)then
        select case (na)
        case (0)
          if(reg)then
            f=confhgrr0(dble(b(1)),z)*gammai(dble(b(1)))
          else
            f=confhgrr0(dble(b(1)),z)
          endif
          return
        case (1)
          if(reg)then
            f=confhgrr1(dble(a(1)),dble(b(1)),z)*gammai(dble(b(1)))
          else
            f=confhgrr1(dble(a(1)),dble(b(1)),z)
          endif
          return
        case (2)
          f=chg(a(1),a(2),b(1),dcmplx(z,0.d0),reg)
          return
        end select
      endif
      u=1.d0
      if(reg)then
        do i=1,nb
          u=u*gamma(dble(b(i)))
        enddo
        no=nb*nogam
      endif
      f1=u
      k=0.d0
      no=0
      main: do
        do i=1,na
          u=u*(dble(a(i))+k)
          if(u == 0.d0)then
c            write(*,'(a,3i5,1p10g12.4)')'hgpq-u0 ',i,na,nb,k,a(i),u
            f=f1
            return
          endif
        enddo
        do i=1,nb
          u=u/(dble(b(i))+k)
        enddo
        k=k+1.d0
        u=u*z/k
        no=no+2*(na+nb+1)
        f1=f1+u
c        write(*,'(a,3i5,1p10g12.4)')'hgpq-0 ',i,na,nb,z,u,f
        if(u**2 <= no*f1**2*epso)then
          f=f1
          return
        endif
      enddo main
      end function 

      real*8  function pgin(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=polygamma(pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      real*8  function pgin1(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=gpolygamma2(2.d0,pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      real*8  function pgip(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=gpolygamma2(pgnc,pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      complex*16  function cpgin(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=cpolygamma(cpgz*t+1.d0)*(1.d0-t)**cpgn
      return
      end

      complex*16  function cpgin1(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=cgpolygamma2((2.d0,0.d0),cpgz*t+1.d0)*(1.d0-t)**cpgn
      return
      end

      complex*16  function cpgip(t) result(f1)
      use gammaf
      implicit none
      real*8 ,intent(in):: t
      f1=cgpolygamma2(dcmplx(pgnc,0.d0),cpgz*t+1.d0)*
     $     (1.d0-t)**cpgn
      return
      end

      complex*16 function cpolygamma2(n,x0) result(f1)
      use tfstk, only:ktfenanq
      use gammaf
      use mathfun
      use macmath
      implicit none
      complex*16 ,intent(in):: n,x0
      complex*16 u,f,lx,df,x
      complex*16 ,external:: cbint
      integer*4 i
      integer*4 ,parameter :: nmax=1000
      complex*16 xk,xn,xk1,sabcp,cab,f10,xkk,xkk1,j1n,sabc
      complex*16 ,save :: ns=(-1.d100,0.d0),cgi(-2:nmax),
     $     chk1(0:nmax),chk3(0:nmax),chk2(0:nmax),
     $     chk4(0:nmax),chk5(0:nmax),chk6(0:nmax)
      real*8 j,j1,axk,m,k,dfa
      integer*4 ,save :: im1,im2,im3,im4,im5,im6
      integer*4 jj,icg,no
      x=zeroim(x0)
      if(x == (0.d0,0.d0))then
        if(dble(n) .ge. -1.d0)then
          f1=dcmplx(1.d0/0.d0,0.d0)
        else
          f1=(0.d0,0.d0)
        endif
        return
      elseif(imag(n) == 0.d0)then
        if(dble(n) == 0.d0)then
          f1=cpolygamma(x)
          return
        elseif(imag(x) == 0.d0 .and. dble(x) .ge. 0.d0)then
          f1=dcmplx(polygamma2(dble(n),dble(x)),0.d0)
          return
        endif
        if(dble(n) == anint(dble(n)) .and.
     $         dble(n) > 0.d0)then
          f1=(-1.d0)**(dble(n)+1.d0)*factorial(dble(n))
     $         *chzeta2(n+1.d0,x)
          return
        endif
      endif
      if(abs(x) .lt. xmth)then
        u=x*cgammai(2.d0-n)
        f=(m_euler-log(x)+cpolygamma(-n))*cgammai(-n)/x
     $       -m_euler*cgammai(1.d0-n)+zt2*u
        k=2.d0
        no=nogam
        do
          u=-u*x*k/(k-n)
          df=zeta(k+1.d0)*u
          f=f+df
          if(abs(df)**2 <= no*abs(f)**2*epso)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
          no=no+nozt+10
        enddo
      elseif(abs(n+1.d0)+abs(x-1.d0) .lt. xmeps)then
        f1=(n+1.d0)*(dpmg1+ddpmg1*.5d0*(n+1.d0)+
     $       dpzmg1*(x-1.d0))
     $       +(x-1.d0)*(-m_euler+(x-1.d0)*.5d0*zt2)
      elseif(abs(n+1.d0)+abs(x-2.d0) .lt. xmeps)then
        f1=(n+1.d0)*(dpmg2+ddpmg2*.5d0*(n+1.d0)+
     $       dpzmg2*(x-2.d0))
     $       +(x-2.d0)*((1.d0-m_euler)+(x-1.d0)*.5d0*(zt2-1.d0))
c$$$      elseif(abs(x) > xmp)then
c$$$        u=-cgammai(-n)/x**(n+1.d0)
c$$$        f=cgpolygamma2(n,x)+dz0*u
c$$$        k=2.d0
c$$$        do
c$$$          u=-u*(n+k)/x
c$$$          f1=f+u*gpolygamma2(-k,1.d0)
c$$$          if(f1 == f)then
c$$$            exit
c$$$          endif
c$$$          f=f1
c$$$          k=k+1.d0
c$$$        enddo
      else
        if(ns /= n)then
          ns=n
          im1=0
          im2=0
          im3=0
          im4=0
          im5=0
          im6=0
          icg=1
          cgi(-2)=cgammai(2.d0-n)
          cgi(-1)=cgammai(1.d0-n)
          cgi(0)=cgammai(-n)
          cgi(1)=cgammai(-1.d0-n)
          chk1(0)=-m_euler-cpolygamma(1.d0-n)
          chk2(0)=cgi(-2)
          chk3(0)=chk2(0)
          chk4(0)=cgammai(2.d0+n)
          chk5(0)=chk4(0)
          chk6(0)=(-m_euler-cpolygamma(-n))
        endif
        f=(n*log(x)-m_euler*(x+n)-n*cpolygamma(-n))*cgi(-1)/x
     $       +cgi(-2)*zt2*x
        no=nogam*2+nopg*2
        k=1.d0
        f1=f
        do
          f10=f
          dfa=abs(f)**2*epso
          xk=zeroim(-x/k)
          xkk=xk/k
          axk=abs(xk)
c          write(*,'(a,1p10g12.4)')'cp2-k ',k,xk,n,f
c          if(.true.)then
c            f=f+(cgi(-2)-chg((1.d0,0.d0),(2.d0,0.d0),
c     $           (2.d0,0.d0)-n,xk))*xkk
          if(abs(xk-xs1) .lt. epsx .or. abs(xk-xs2) .lt. epsx)then
            f=f+(cgi(-2)-chg1s((1.d0,0.d0),(2.d0,0.d0),
     $           (2.d0,0.d0)-n,xk,.true.))*xkk
          elseif(dble(xk) <= 0.5d0)then
            if(axk .ge. 1.d0)then
              xk1=zeroim(1.d0/(1.d0-xk))
              xkk1=xkk*xk1
              lx=-log(xk1)
              u=-cgi(0)*xkk1*xk1
              f=f-cgi(-1)*xkk1-u*(lx+chk1(0))+cgi(-2)*xkk
              j=1.d0
              do jj=1,nmax
                j1=j+1.d0
                u=u*(j-n)/j*xk1
                if(jj > im1)then
                  chk1(jj)=polygamma(j1)-cpolygamma(j1-n)
                  im1=jj
                  no=no+nopg*2
                endif
                df=-u*(lx+chk1(jj))
                f=f+df
                if(abs(df)**2 <= no*dfa)then
                  no=no+16*3+4+jj*(3*3+1)
c                  write(*,'(a,1p10g12.4)')'pgm-1 ',xk,f1,f10
                  exit
                endif
                j=j1
              enddo
            elseif(abs(xk-1.d0) > 1.d0)then
c$$$              f=(one-x)**(-a)*chg(a,c-b,c,x/(x-1.d0))
              xk1=zeroim(xk/(xk-1.d0))
              xkk1=xk1/k
              f=f+cgi(-2)*(xkk1+xkk)
              xn=xk1*xkk1
              do i=1,nmax
                if(i > im2)then
                  chk2(i)=chk2(i-1)*((i-1)-n)/((i+1)-n)
                  im2=i
                endif
                df=chk2(i)*xn
                f=f+df
                if(abs(df)**2 <= no*dfa)then
c                  write(*,'(a,1p10g12.4)')'pgm-2 ',xk,f1,f10
                  no=no+6*3+i*(6+1)
                  exit
                endif
                xn=xn*xk1
              enddo
            else
              xn=-xk*xkk
              do i=1,nmax
                if(i > im3)then
                  chk3(i)=chk3(i-1)*(i+1)/((1+i)-n)
                  im3=i
                endif
                df=chk3(i)*xn
                f=f+df
                if(abs(df)**2 <= no*dfa)then
                  no=no+3*3+i*(3+1)
c                  if(k <= 3.d0 .or. mod(k,100000.d0) == 0.d0)then
c                    write(*,'(a,1p10g12.4)')'pgm-3 ',k,xk,f10,f-f10
c                  endif
                  exit
                endif
                xn=xn*xk
              enddo
            endif
          else
c            if(abs(xk-1.d0) .ge. 1.d0)then
c              f=f+(cgi(-2)-chg((1.d0,0.d0),(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-3 ',xk,f,f-f10
            if(abs(xk-1.d0) .ge. 1.d0)then
              xk1=1.d0/xk
              xkk1=xk1*xkk
              lx=log(zeroim(-xk))
              u=xk1*xkk1
              f=f+cgi(-1)*xkk1+u*(lx+chk6(0))*cgi(0)
     $             +cgi(-2)*xkk
              j=0.d0
              icg=0
              do jj=1,nmax
                j1=j+1.d0
                u=-u*xk1/j1
                j1n=-j1-n
                if(jj > icg)then
                  cgi(jj)=cgammai(j1n)
                  icg=jj
                  no=no+nogam
                endif
                if(jj > im6)then
                  chk6(jj)=polygamma(j1+1.d0)-cpolygamma(j1n)
                  im6=jj
                  no=no+nopg*2
                endif
                if(imag(j1n) == 0.d0 .and.
     $               dble(j1n) == anint(dble(j1n)) .and.
     $               dble(j1n) <= 0.d0)then
                  df=-u*(-1.d0)**(1.d0-dble(j1n))*gamma(1.d0-dble(j1n))
                  no=no+nogam
                else
                  df=u*(lx+chk6(jj))*cgi(jj)
                endif
                f=f+df
                if(abs(df)**2 <= no*dfa)then
                  no=no+19*3+jj*(4*3+2)
c                  write(*,'(a,i8,1p10g12.4)')'cpg-6 ',no,j1+n,xk,df,f
                  exit
                endif
                j=j1
              enddo
c            elseif(axk > 1.d0)then
c              f=f+(cgi(-2)-chg((1.d0,0.d0),(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-5 ',xk,f,f-f10
            elseif(axk > 1.d0)then
              xk1=zeroim(1.d0-1.d0/xk)
              cab=-1.d0-n
              m=dble(cab)
              f=f+cgi(-2)*xkk
              if(imag(n) == 0.d0 .and. m == anint(m)
     $             .or. csinp(cab) == 0.d0)then
                if(m /= 0.d0)then
c                  u=-gamma(m)*gammai(m+2.d0)*gammai(m+1.d0)/k
                  u=-gammai(m+2.d0)/m/k
                  f=f+u
                  do i=1,int(m)-1
                    u=u*(3.d0+m-i)/(m-i)*xk1
                    f=f+u
                  enddo
                  no=no+nogam+int(m)*3
                endif
                lx=log(zeroim(-xk1))
                u=xk1**m*gammai(m+1.d0)/k
c                write(*,'(a,1p10g12.4)')'cpg2-n50 ',xk,m,u,f-f10
                f=f+u*((lx+1.d0)-xk1*(lx-1.d0))
                u=u*xk1
                j=2.d0
                no=no+nogam
c                write(*,'(a,1p10g12.4)')'cpg2-n51 ',xk,m,u,f-f10
                do
                  u=u/j*max(1.d0,j-2.d0)*xk1
                  f=f-u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+20*3+7+int(j)*(3*3+1)
c                    write(*,'(a,1p10g12.4)')'cpg2-n5 ',xk,m,j,u,f1-f10
                    exit
                  endif
                  j=j+1.d0
                enddo
              else
                sabcp=m_pi/csinp(cab)
                xkk1=-cgi(-1)*cgi(0)/k*sabcp
c     f1=( chg(a,a-c+1.d0,1.d0-cab,x1)
c     $             *cgammai(c-a)*cgammai(c-b)*x**(-a)
c     $             -chg(c-a,1.d0-a,1.d0+cab,x1)
c     $             *(one-x)**cab*x**(a-c)
c     $             *cgammai(a)*cgammai(b))/sabc*m_pi
                f=f+chk5(0)*xkk1
     $               +cgi(0)*zeroim(1.d0-xk)**cab*xk**(n-1.d0)*sabcp*xkk
                xn=xk1*xkk1
                do i=1,nmax
                  if(i > im5)then
                    chk5(i)=chk5(i-1)*(n+i-1)/(i+1+n)
                    im5=i
                  endif
                  df=chk5(i)*xn
                  f=f+df
                  if(abs(df)**2 <= no*dfa)then
c                    write(*,'(a,i8,1p10g12.4)')'pgm-5 ',no,xk,f,df
                    no=no+(23*3+2)+i*(3+1)
                    exit
                  endif
                  xn=xn*xk1
                enddo
              endif
c            elseif(.true.)then
c              f=f+(cgi(-2)-chg((1.d0,0.d0),(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-4 ',xk,f,f-f10
            else
              xk1=zeroim(1.d0-xk)
              cab=-1.d0-n
              m=dble(cab)
              f=f+cgi(-2)*xkk
              sabc=csinp(cab)
              if(imag(n) == 0.d0 .and. m == anint(m)
     $             .or.  sabc == 0.d0)then
                if(m /= 0.d0)then
c                   u=-exp(-log_gamma(m+1.d0)-log_gamma(m+2.d0))*xkk
                  u=-gammai(1.d0+m)*gammai(2.d0+m)*xkk
                  f=f+u
                  do i=1,int(m)-1
                    u=-u*(1.d0+i)/(m-i)*xk1
                    f=f+u
                  enddo
                  no=no+nogam*2
                endif
                if(xk1 /= (0.d0,0.d0))then
                  lx=log(xk1)
                  u=zeroim(-xk1)**m*gammai(m+1.d0)*xkk
                  f=f+u*(lx+m_euler+polygamma(2.d0+m))
                  j=0.d0
                  no=no+nogam+nopg
                  do
                    j1=j+1.d0
                    u=u*(j1+m+1.d0)/j1*xk1
                    df=u*(lx-polygamma(j1+1.d0)+polygamma(2.d0+m+j1))
                    f=f+df
                    no=no+nopg*2
                    if(abs(df)**2 <= no*dfa)then
c                      write(*,'(a,1p10g12.4)')'cpg2-n4 ',
c     $                     xk,m,j1,f1,f1-f10
                      exit
                    endif
                    j=j1
                  enddo
                endif
c              elseif(xk == (1.d0,0.d0))then
c                f1=f+dcmplx(1.d0/0.d0,0.d0)
c                return
              else
                if(xk1 == 0.d0 .and.
     $               (imag(cab) /= 0.d0 .or. dble(cab) .lt. 0.d0))then
                  f1=dcmplx(-1.d0/0.d0,-1.d0/0.d0)
                  return
                endif
                sabcp=m_pi/sabc
c$$$  f1=(chg(a,b,1.d0-cab,x1)
c$$$  $             *cgammai(c-a)*cgammai(c-b)
c$$$  $             -chg(c-a,c-b,cab+1.d0,x1)*x1**cab
c$$$  $             *cgammai(a)*cgammai(b))/sabc*m_pi
                xkk1=-cgi(-1)*cgi(0)*xkk*sabcp
                f=f+chk4(0)*xkk1
                xn=xk1*xkk1
                do i=1,nmax
                  if(i > im4)then
                    chk4(i)=chk4(i-1)*(i+1)/(i+1+n)
                    im4=i
                  endif
                  u=chk4(i)*xn
                  f=f+u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+25+i*4
                    exit
                  endif
                  xn=xn*xk1
                enddo
                u=cgi(0)*xk1**cab*xkk*sabcp
c                write(*,'(a,1p10g12.4)')'pgm-41 ',xk,xk1,cab,sabcp,u
                f=f+u
                j=0.d0
                do
                  j1=j+1.d0
                  u=u*(j1-n)/j1*xk1
                  f=f+u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+(3*4+1)+int(j1)*(3*3+1)
c                    write(*,'(a,i8,1p10g12.4)')'pgm-4 ',no,xk,f-f10
                    exit
                  endif
                  j=j1
                enddo
              endif
            endif
          endif
          if(abs(f-f10)**2 <= no*dfa)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
        enddo
c      real*8 ,parameter ::eps=1.d-9
c$$$      elseif(dble(n) > 0.d0)then
c$$$        cpgz=x
c$$$        in=ceiling(dble(n))
c$$$        pgnc=in+2
c$$$        cpgn=pgnc-n-1.d0
c$$$        u=-cgammai(1.d0-n)/x**n
c$$$        f1=((cbint(cpgip,0.d0,1.d0,eps,eps)*x
c$$$     $       +gpolygamma2(pgnc-1.d0,1.d0))*x/cpgn
c$$$     $       +gpolygamma2(pgnc-2.d0,1.d0))/(cpgn-1.d0)
c$$$     $       *cgammai(cpgn-1.d0)*x**(cpgn-1.d0)
c$$$     $       -((log(x)-cpolygamma(-n)-m_euler)*n/x-m_euler)*u
c$$$c        write(*,'(a,1p10g12.4)')'cp2-p ',pgnc,f1,
c$$$c     $      ( cbint(pgip,0.d0,1.d0,1.d-2,1.d-11)*x
c$$$c     $       +gpolygamma2(pgnc,1.d0))/pgn
c$$$        do i=1,in-1
c$$$          u=-u*x*i/(i-n)
c$$$          f1=f1+zeta(i+1.d0)*u
c$$$        enddo
c$$$      elseif(dble(n) .ge. -2.d0)then
c$$$        cpgz=x
c$$$        cpgn=1.d0-n
c$$$        f1=(((cbint(cpgin1,0.d0,1.d0,eps,eps)*x
c$$$     $       +m_pi**2/6.d0)*x/cpgn
c$$$     $       -m_euler)/(-n)
c$$$     $       -(log(x)-cpolygamma(-n)-m_euler)/x)/x**n*cgammai(-n)
c$$$      else
c$$$        cpgz=x
c$$$        cpgn=-1.d0-n
c$$$        f1=(cbint(cpgin,0.d0,1.d0,eps,eps)
c$$$     $       -(log(x)-cpolygamma(-n)-m_euler)/x)/x**n*cgammai(-n)
      endif
      return
      end function

      real*8 recursive function polygamma2(n,x) result(f1)
      use tfstk, only:ktfenanq
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: n,x
      integer*4 ,parameter :: nmax=1000
      real*8 ,parameter ::eps=1.d-9
      real*8 f,xk,xn,u,xk1,lx,f10,xkk,xkk1,df
      real*8 ,save :: ns=(-1.d100,0.d0),cgi(-2:nmax),
     $     chk1(0:nmax),chk2(0:nmax)
      real*8 k,j,j1,axk,dfa
      integer*4 ,save :: im1,im2
      integer*4 i,jj,icg,no
      if(x == 0.d0)then
        if(n .ge. -1.d0)then
          f1=1.d0/0.d0
        else
          f1=0.d0
        endif
      elseif(n == anint(n) .and. n .ge. 0.d0)then
        f1=gpolygamma2(n,x)
      elseif(abs(x) .lt. xmth)then
        u=gammai(2.d0-n)*x
        f=(m_euler-log(x)+polygamma(-n))*gammai(-n)/x
     $       -m_euler*gammai(1.d0-n)+zt2*u
        no=nogam*2
        k=2.d0
        do
          u=-u*x*k/(k-n)
          df=zeta(k+1.d0)*u
          f=f+df
          if(df**2 <= no*f**2*epso)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
          no=no+nozt
        enddo
      elseif(abs(n+1.d0)+abs(x-1.d0) .lt. xmeps)then
        f1=(n+1.d0)*(dpmg1+ddpmg1*.5d0*(n+1.d0)+
     $       dpzmg1*(x-1.d0))
     $       +(x-1.d0)*(-m_euler+(x-1.d0)*.5d0*zt2)
      elseif(abs(n+1.d0)+abs(x-2.d0) .lt. xmeps)then
        f1=(n+1.d0)*(dpmg2+ddpmg2*.5d0*(n+1.d0)+
     $       dpzmg2*(x-2.d0))
     $       +(x-2.d0)*(1.d0-m_euler+(x-2.d0)*.5d0*(-1.d0+zt2))
c$$$      elseif(abs(x) > xmp)then
c$$$        u=-gammai(-n)/x**(n+1.d0)
c$$$        f=gpolygamma2(n,x)+dz0*u
c$$$        k=2.d0
c$$$        do
c$$$          u=-u*(n+k)/x
c$$$          f1=f+u*gpolygamma2(-k,1.d0)
c$$$          if(f1 == f)then
c$$$            exit
c$$$          endif
c$$$          f=f1
c$$$          k=k+1.d0
c$$$        enddo
      else
        if(ns /= n)then
          ns=n
          im1=0
          im2=0
          icg=1
          cgi(-2)=gammai(2.d0-n)
          cgi(-1)=gammai(1.d0-n)
          cgi(0)=gammai(-n)
          cgi(1)=gammai(-1.d0-n)
          chk1(0)=-m_euler-polygamma(1.d0-n)
          chk2(0)=cgi(-2)
        endif
        f=(n*log(x)-m_euler*(x+n)-n*polygamma(-n))*cgi(-1)/x
     $       +cgi(-2)*zt2*x
c        f=(m_euler-log(x)+polygamma(-n))*cgi(0)/x
c     $       -m_euler*cgi(-1)+cgi(-2)*zt2*x
        k=1.d0
        no=nogam*2
        do
          f10=f
          dfa=f**2*epso
          xk=-x/k
          xkk=xk/k
          axk=abs(xk)
          if(xk <= 0.5d0)then
            if(axk .ge. 1.d0)then
              xk1=1.d0/(1.d0-xk)
              xkk1=xkk*xk1
              lx=-log(xk1)
              u=-cgi(0)*xkk1*xk1
              f=f-cgi(-1)*xkk1-u*(lx+chk1(0))+cgi(-2)*xkk
              j=1.d0
              do jj=1,nmax
                j1=j+1.d0
                u=u*(j-n)/j*xk1
                if(jj > im1)then
                  chk1(jj)=polygamma(j1)-polygamma(j1-n)
                  im1=jj
                endif
                df=-u*(lx+chk1(jj))
                f=f+df
                if(df**2 <= no*dfa)then
                  no=no+16*1+4+jj*(3*1+1)
c                  write(*,'(a,i5,1p10g12.4)')'pgm-1 ',jj,xk,f10-f1
                  exit
                endif
                j=j1
              enddo
            elseif(abs(xk-1.d0) > 1.d0)then
c$$$              f=(one-x)**(-a)*chg(a,c-b,c,x/(x-1.d0))
              xk1=xk/(xk-1.d0)
              xkk1=xk1/k
              f=f+cgi(-2)*(xkk1+xkk)
              xn=xk1*xkk1
              do i=1,nmax
                if(i > im2)then
                  chk2(i)=chk2(i-1)*((i-1)-n)/((i+1)-n)
                  im2=i
                endif
                df=chk2(i)*xn
                f=f+df
                if(df**2 <= no*dfa)then
c                  write(*,'(a,1p10g12.4)')'pgm-2 ',xk,f1,f10
                  no=no+6*1+i*(6+1)
                  exit
                endif
                xn=xn*xk1
              enddo
            else
              go to 9000
            endif
          else
            go to 9000
          endif
          if(abs(f-f10)**2 <= no*dfa)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
        enddo
      endif
      return
 9000 write(*,'(a,1p8g15.7)')
     $     'polygamma2-implementation error ',n,x
      f1=0.d0
      return
      end function

      complex*16 function cgammaq(a,x) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x
      if(a == (0.d0,0.d0) .and. x == (0.d0,0.d0))then
        f1=(0.d0,0.d0)
      else
        f1=1.d0-confhg1(a,a+1.d0,-x)*zeroim(x)**a
      endif
      return
      end

      real*8 function gammaq(a,x) result(f1)
      implicit none
      real*8 ,intent(in):: a,x
      if(a == 0.d0 .and. x == 0.d0)then
        f1=0.d0
      else
        f1=1.d0-confhgrr1(a,a+1.d0,-x)*x**a
      endif
      return
      end

      complex*16   function cgammap(a,x) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x
      if(a == (0.d0,0.d0) .and. x == (0.d0,0.d0))then
        f1=(1.d0,0.d0)
      else
        f1=confhg1(a,a+1.d0,-x)*zeroim(x)**a
      endif
      return
      end

      complex*16 recursive  function cgamma2(a,x) result(f1)
      use gammaf
      use macmath
      implicit none
      complex*16 ,intent(in):: a,x
      f1=exp(-x)*chgu(1.d0-a,1.d0-a,x)
      return
      end function

      real*8 recursive   function gamma2(a,x) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,x
      f1=exp(-x)*hgu(1.d0-a,1.d0-a,x)
      return
      end function

      complex*16 recursive function chgu(a,b,x0) result(f1)
      use gammaf
      use macmath
      implicit none
      complex*16 ,intent(in):: a,b,x0
      complex*16 lx,u,ab1,x,df
      real*8 k,k1,rb,ra,rab1
      integer*4 n,i,no
      ra=dble(a)
      rb=dble(b)
      x=zeroim(x0)
      if(imag(a) == 0.d0 .and. imag(b) == 0.d0
     $     .and. imag(x) == 0.d0 .and. dble(x) .ge. 0.d0)then
        f1=dcmplx(hgu(ra,rb,dble(x)),0.d0)
        return
      endif
      ab1=a-b+1.d0
      rab1=dble(ab1)
      if(x == (0.d0,0.d0) .and. rb .lt. 1.d0)then
        f1=cpochh(ab1,-a)
      elseif(imag(a) == 0.d0 .and.
     $       ra == anint(ra) .and. ra <= 0.d0)then
        n=nint(-ra)
        u=x**n
        k=0.d0
        f1=u
        do i=1,n
          k1=k+1.d0
          u=-u*(-ra-k)/k1*(rb-ra-k1)/x
          f1=f1+u
          k=k1
        enddo
      elseif(imag(ab1) == 0.d0 .and.
     $       rab1 == anint(rab1) .and. rab1 <= 0.d0)then
c        f1=(-1.d0)**nint(-rab1)/x**(b-1.d0)*cpochh(2.d0-b,-ab1)
c     $       *confhg1(ab1,2.d0-b,x,.false.)
        f1=(-1.d0)**nint(-rab1)*exp(-(b-1.d0)*log(x)-cloggamma(2.d0-b)
     $       +cloggamma(1.d0-a))
     $       *confhg1(ab1,2.d0-b,x,.false.)
      elseif(imag(b) == 0.d0 .and. anint(rb) == rb)then
        n=nint(rb-1.d0)
        if(n <= -1.d0)then
          f1=x**(-n)*chgu(ab1,2.d0-b,x)
        else
          lx=log(x)
          if(n > 0)then
            u=exp((1.d0-rb)*lx+log_gamma(rb-1.d0)-cloggamma(a))
            f1=u
            k=2.d0
            do i=2,n
              u=u/max(rb-1.d0-k,1.d0)*(rb-a-k)/k*x
              f1=f1+u
              k=k+1.d0
            enddo
          else
            f1=(0.d0,0.d0)
          endif
c          u=-(-1.d0)**n*cgammai(a-rb+1.d0)*gammai(rb)
          u=-(-1.d0)**n*exp(-cloggamma(a-rb+1.d0)-log_gamma(rb))
          f1=f1+u*(lx+cpolygamma(a)+m_euler-polygamma(rb))
          k=0.d0
          no=nolog*2+nogam*3+n*5
          do
            k1=k+1.d0
            u=u*(a+k)/k1/(rb+k)*x
            df=u*(lx+cpolygamma(a+k1)-polygamma(k1+1.d0)
     $           -polygamma(rb+k1))
            f1=f1+df
            no=no+nogam*3
            if(abs(df)**2 <= no*abs(f1)**2*epso)then
              exit
            endif
            k=k1
          enddo
        endif
      else
        f1=(cgammai(ab1)*confhg1(a,b,x,.true.)
     $       -exp((1.d0-b)*log(x)-cloggamma(a))
     $       *confhg1(ab1,2.d0-b,x,.true.))*m_pi/csinp(b)
      endif
      return
      end function

      real*8 recursive function hgu(a,b,x) result(f1)
      use gammaf
      use macmath
      implicit none
      real*8 ,intent(in):: a,b,x
      real*8 lx,u,f,k,k1,ab1
      integer*4 n,i
      ab1=a-b+1.d0
c      write(*,'(a,1p10g12.4)')'hgu-0 ',a,b,x,ab1
      if(x == 0.d0 .and. b .lt. 1.d0)then
        f1=pochh(ab1,-a)
      elseif(a == anint(a) .and. a <= 0.d0)then
        n=nint(-a)
        u=x**n
        k=0.d0
        f1=u
        do i=1,n
          k1=k+1.d0
          u=-u*(-a-k)/k1*(b-a-k1)/x
          f1=f1+u
          k=k1
        enddo
      elseif(ab1 == anint(ab1) .and. ab1 <= 0.d0)then
        f1=(-1.d0)**nint(-ab1)*pochh(2.d0-b,-ab1)*
     $       confhgrr1(ab1,2.d0-b,x,.false.)/x**(b-1.d0)
c        write(*,'(a,1p10g12.4)')'hgu-ab1 ',ab1,f1,
c     $       confhgrr1(ab1,2.d0-b,x,.false.),pochh(2.d0-b,-ab1)
      elseif(b == anint(b))then
        n=nint(b-1.d0)
        if(n .ge. 0)then
          if(n == 0)then
            f=0.d0
          else
c     u=exp(-log_gamma(a)+log_gamma(b-a)-log_gamma(2.d0-a)
c     $         -log_gamma(b-1.d0))/x
            u=pochh(2.d0-a,b-2.d0)*gammai(a)*gammai(b-1.d0)/x
            f=u
            k=1.d0
            do i=2,n
              k1=k+1.d0
              u=u*k*(b-k1)/(k1-a)/x
              f=f+u
              k=k1
            enddo
          endif
          lx=log(x)
c     u=-(-1.d0)**n*exp(-log_gamma(ab1)-log_gamma(b))
          u=-(-1.d0)**n*gammai(ab1)*gammai(b)
c     write(*,'(a,1p10g12.4)')'hgu ',n,u
          f=f+u*(lx+polygamma(a)+m_euler-polygamma(b))
          k=0.d0
          do
            k1=k+1.d0
            u=u*(a+k)/k1/(b+k)*x
            f1=f+u*(lx+polygamma(a+k1)-polygamma(k1+1.d0)
     $           -polygamma(b+k1))
            if(f1 == f)then
              exit
            endif
            f=f1
            k=k1
          enddo
        else
          f1=x**(-n)*hgu(ab1,2.d0-b,x)
        endif
      else
        f1=(gammai(ab1)*confhgrr1(a,b,x,.true.)
     $       -x**(1.d0-b)*gammai(a)
     $       *confhgrr1(ab1,2.d0-b,x,.true.))*m_pi/sinp(b)
      endif
      return
      end function

      complex*16 function cerf(z) result(f1)
      implicit none
      complex*16 ,intent(in):: z
      f1=2.d0*z/m_sqrtpi
     $     *confhg1((.5d0,0.d0),(1.5d0,0.d0),zeroim(-z**2),.false.)
      return
      end function

      complex*16 function cerfc(z) result(f1)
      use macmath
      implicit none
      complex*16 ,intent(in):: z
      f1=1.d0-2.d0*z/m_sqrtpi*confhg1((.5d0,0.d0),(1.5d0,0.d0),
     $     zeroim(-z**2),.false.)
      return
      end function

      complex*16 recursive function cinverseerf(x) result(f1)
      use gammaf,only:nolog
      use macmath
      implicit none
      complex*16 ,intent(in):: x
      complex*16 df,f,x2
      real*8, parameter:: c=sqrt(m_pi_2),eps=5.d-17**2,
     $     dimax=2.d0,drmax=0.5d0,dmin=0.1d0
      integer*4 ,parameter :: imax=1000
      real*8 adf,adf0,fact
      integer*4 no,i
      if(dble(x) .lt. 0.d0)then
        f1=-cinverseerf(zeroim(-x))
c        write(*,'(a,1p10g12.4)')'cierf-1 ',x,f1
        return
      elseif(sign(1.d0,imag(x)) .lt. 0.d0)then
        f1=conjz(cinverseerf(conjz(x)))
c        write(*,'(a,1p10g12.4)')'cierf-2 ',x,f1
        return
      endif
      If(abs(x-1.d0) + abs(x-1.8d0) .lt. 1.2d0)then
        f1=log(sqrt(2.d0/m_pi)/(1.d0-x))
        f1=sqrt(zeroim(f1-.5d0*log(2.d0*f1)))
        no=nolog*2
      elseif(abs(x) .lt. 1.d0)then
        x2=m_pi*x**2
        f1=.5d0*m_sqrtpi*x*(1.d0+x2*(1.d0/12.d0
     $       +x2*(7.d0/480.d0+x2*(127.d0/40320.d0
     $       +x2*4369.d0/5806080.d0))))
        no=20+nolog
      else
        f1=sqrt(-4.d0/m_pi*log(0.5d0-(0.d0,1.d0)*x))
        no=nolog+10
      endif
      f=f1
      adf0=1.d100
      fact=1.d0
      do i=1,imax
        df=(x-cerf(f))*c*exp(f**2)
        adf=abs(df)**2
        fact=max(1.d0,sqrt(adf/adf0))
        if(fact > 1.d4)then
          exit
        endif
        f1=f+df/fact
c        if(fact /= 1.d0)then
c          write(*,'(a,1p10g12.4)')'cierf-i ',x,fact,f1,df
c        endif
        if(adf/fact**2 <= no*abs(f1)**2*eps)then
          exit
        endif
        f=f1
        adf0=adf*4.d0
        no=no+40
      enddo
c      write(*,'(a,1p10g12.4)')'cierf-3 ',x,f1
      return
      end

      real*8  recursive function inverseerf(x) result(f1)
      use gammaf,only:nolog
      use macmath
      implicit none
      real*8 ,intent(in):: x
      real*8 df,x2
      integer*4 no,i
      real*8, parameter:: c=sqrt(m_pi_2),eps=5.d-17
      integer*4 ,parameter :: imax=1000
      if(x .lt. 0.d0)then
        f1=-inverseerf(-x)
        return
      endif
      if(abs(x) == 1.d0)then
        f1=x/0.d0
        return
      elseif(x > 0.8d0)then
        f1=log(sqrt(2.d0/m_pi)/(1.d0-x))
        f1=sqrt(f1-.5d0*log(2.d0*f1))
        no=nolog*2+10
      else
        x2=m_pi*x**2
        f1=.5d0*m_sqrtpi*x*(1.d0+x2*(1.d0/12.d0
     $       +x2*(7.d0/480.d0+x2*(127.d0/40320.d0
     $       +x2*4369.d0/5806080.d0))))
        no=11
      endif
      do i=1,imax
        df=(x-erf(f1))*c*exp(f1**2)
        f1=f1+df
        if(df**2 <= no*f1**2*eps)then
          exit
        endif
        no=no+nolog
      enddo
      return
      end

      complex*16 function cbesj(n,z) result(f)
      implicit none
      complex*16 ,intent(in):: n,z
      f=zeroim(.5d0*z)**n*confhg0(n+1.d0,-.25d0*z**2)
      return
      end

      complex*16 function cbesi(n,z) result(f)
      implicit none
      complex*16 ,intent(in):: n,z
      f=zeroim(.5d0*z)**n*confhg0(n+1.d0,.25d0*z**2)
      return
      end

      complex*16 function cbesy(n,z0) result(f)
      use macmath
      implicit none
      complex*16 ,intent(in):: n,z0
      complex*16 z
      z=zeroim(z0)
      f=((.5d0*z)**n*ccosp(n)*confhg0(n+1.d0,-.25d0*z**2)
     $     -(2.d0/z)**n*confhg0(1.d0-n,-.25d0*z**2))/csinp(n)
      return
      end

      complex*16 function cbesk(n,z0) result(f)
      use macmath
      implicit none
      complex*16 ,intent(in):: n,z0
      complex*16 z
      z=zeroim(z0)
      f=.5d0*m_pi/csinp(n)*((2.d0/z)**n*confhg0(1.d0-n,.25d0*z**2)
     $     -(.5d0*z)**n*confhg0(1.d0+n,.25d0*z**2))
      return
      end
      end module
