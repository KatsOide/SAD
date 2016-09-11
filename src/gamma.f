c$$$      real*8 function gamma(x)
c$$$      implicit none
c$$$      real*8 x,factorial
c$$$      gamma=factorial(x)/x
c$$$      return
c$$$      end

      real*8 function factorial(x)
      use macmath
      implicit none
c     Including pi = m_pi
      integer*4 ix,lt
      parameter (lt=23)
      real*8 x,aloggamma1,fac(0:lt)
      data fac/
     $1d0,
     $1d0,
     $2d0,
     $6d0,
     $24d0,
     $120d0,
     $720d0,
     $5040d0,
     $40320d0,
     $362880d0,
     $3628800d0,
     $39916800d0,
     $479001600d0,
     $6227020800d0,
     $87178291200d0,
     $1307674368000d0,
     $20922789888000d0,
     $355687428096000d0,
     $6402373705728000d0,
     $121645100408832000d0,
     $2432902008176640000d0,
     $51090942171709440000d0,
     $1124000727777607680000d0,
     $25852016738884976640000d0/
      if(x .lt. 0.d0)then
        factorial=exp(-aloggamma1(-x))*pi*x/sin(pi*x)
      else
        ix=int(x)
        if(ix .le. lt .and. x-ix .eq. 0.d0)then
          factorial=fac(ix)
        else
          factorial=exp(aloggamma1(x))
        endif
      endif
      return
      end

      real*8 function aloggamma1(x)
      use macmath
      implicit none
c     Including pi = m_pi
      real*8 x,x1
      real*8 c0,c1,c2,c3,c4,c5,c6,c7,c8
      real*8 gamma,alogsqrt2pi,alogpi
      parameter (gamma=5.d0,alogsqrt2pi=0.918938533204672742d0,
     $     alogpi=1.144729885849400174d0)
      parameter (c0=0.999999999999997524d0,
     $      c1=76.1800917309326077d0,
     $      c2=-86.5053204008552d0,
     $      c3=24.01409906379226027d0,
     $      c4=-1.231743354454618365d0,
     $      c5=0.001216872118636531519d0,
     $      c6=-0.00001408915744128554778d0,
     $      c7=4.00491935010387864d-6,
     $      c8=-4.93961492826482964d-7)
      if(x .lt. 0.d0)then
        x1=-x
        aloggamma1=log(pi*x/sin(pi*x))-
     $       ((x1+.5d0)*log(x1+gamma+.5d0)-(x1+gamma+.5d0)
     $       +alogsqrt2pi+log(
     $       c0+c1/(1.d0+x1)+c2/(2.d0+x1)+c3/(3.d0+x1)
     $       +c4/(4.d0+x1)+c5/(5.d0+x1)+c6/(6.d0+x1)
     $       +c7/(7.d0+x1)+c8/(8.d0+x1)))
      else
        aloggamma1=(x+.5d0)*log(x+gamma+.5d0)-(x+gamma+.5d0)
     $       +alogsqrt2pi+log(
     $       c0+c1/(1.d0+x)+c2/(2.d0+x)+c3/(3.d0+x)
     $       +c4/(4.d0+x)+c5/(5.d0+x)+c6/(6.d0+x)
     $       +c7/(7.d0+x)+c8/(8.d0+x))
      endif
      return
      end

      complex*16 function cgamma(x)
      implicit none
      complex*16 x,cloggamma1
      cgamma=exp(cloggamma1(x))/x
      return
      end

      complex*16 function cfactorial(x)
      implicit none
      complex*16 x,cloggamma1
      cfactorial=exp(cloggamma1(x))
      return
      end

      complex*16 function cloggamma1(z)
      use macmath
      implicit none
c     Including pi = m_pi
      complex*16 z,z1
      real*8 c0,c1,c2,c3,c4,c5,c6,c7,c8
      real*8 gamma,alogsqrt2pi,alogpi
      parameter (gamma=5.d0,alogsqrt2pi=0.918938533204672742d0,
     $     alogpi=1.144729885849400174d0)
      parameter (c0=0.999999999999997524d0,
     $      c1=76.1800917309326077d0,
     $      c2=-86.5053204008552d0,
     $      c3=24.01409906379226027d0,
     $      c4=-1.231743354454618365d0,
     $      c5=0.001216872118636531519d0,
     $      c6=-0.00001408915744128554778d0,
     $      c7=4.00491935010387864d-6,
     $      c8=-4.93961492826482964d-7)
      if(dble(z) .lt. 0.d0)then
        z1=-z
        cloggamma1=log(pi*z/sin(pi*z))-
     $       ((z1+.5d0)*log(z1+gamma+.5d0)-(z1+gamma+.5d0)
     $       +alogsqrt2pi+log(
     $       c0+c1/(1.d0+z1)+c2/(2.d0+z1)+c3/(3.d0+z1)
     $       +c4/(4.d0+z1)+c5/(5.d0+z1)+c6/(6.d0+z1)
     $       +c7/(7.d0+z1)+c8/(8.d0+z1)))
      else
        cloggamma1=(z+.5d0)*log(z+gamma+.5d0)-(z+gamma+.5d0)
     $       +alogsqrt2pi+log(
     $       c0+c1/(1.d0+z)+c2/(2.d0+z)+c3/(3.d0+z)
     $       +c4/(4.d0+z)+c5/(5.d0+z)+c6/(6.d0+z)
     $       +c7/(7.d0+z)+c8/(8.d0+z))
      endif
      return
      end

      real*8 function erf(x)
      implicit none
      real*8 x,gammap
      erf=sign(gammap(.5d0,x**2),x)
      return
      end

      real*8 function erfc(x)
      implicit none
      real*8 x,gammaq,gammap
      if(x .lt. 0.d0)then
        erfc=1.d0+gammap(.5d0,x**2)
      else
        erfc=gammaq(.5d0,x**2)
      endif
      return
      end
      
      real*8 function gammap(a,x)
      implicit none
      real*8 a,x,gammaser,gammacf
      if(x .lt. 0.d0 .or. a .lt. 0.d0)then
        gammap=x/0.d0
        return
      endif
      if(a .eq. 0.d0)then
        gammap=1.d0
      elseif(x .lt. a+1.d0)then
        gammap=gammaser(a,x)
      else
        gammap=1.d0-gammacf(a,x)
      endif
      return
      end

      real*8 function gammaq(a,x)
      implicit none
      real*8 a,x,gammaser,gammacf
      if(x .lt. 0.d0 .or. a .lt. 0.d0)then
        gammaq=x/0.d0
        return
      endif
      if(a .eq. 0.d0)then
        gammaq=0.d0
      elseif(x .lt. a+1.d0)then
        gammaq=1.d0-gammaser(a,x)
      else
        gammaq=gammacf(a,x)
      endif
      return
      end

      real*8 function gamma0(x)
      use macmath
      implicit none
c     Including Euler's gamma(euler)
      real*8 x,gamma0ser,gamma0cf
      if(x .lt. 0.d0)then
        gamma0=x/0.d0
        return
      endif
      if(x .lt. 1.d0)then
        gamma0=gamma0ser(x)-log(x)-euler
      else
        gamma0=gamma0cf(x)
      endif
      return
      end

      real*8 function gamma0log(x)
      use macmath
      implicit none
c     Including Euler's gamma(euler)
      real*8 x,gamma0ser,gamma0cf
      if(x .lt. 0.d0)then
        gamma0log=x/0.d0
        return
      endif
      if(x .lt. 1.d0)then
        gamma0log=gamma0ser(x)
      else
        gamma0log=gamma0cf(x)+log(x)+euler
      endif
      return
      end

      real*8 function gammaser(a,x)
      implicit none
      integer*4 itmax,i
      real*8 a,x,eps,gln,ap,sum,del,aloggamma1
      parameter (itmax=300,eps=1.d-13)
      if(x .lt. 0.d0)then
        gammaser=x/0.d0
        return
      elseif(x .eq. 0.d0)then
        gammaser=0.d0
        return
      elseif(a .eq. 0.d0)then
        gammaser=1.d0
        return
      endif
      ap=a
      sum=1.d0
      del=sum
      do i=1,itmax
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del) .lt. abs(sum)*eps)then
          go to 1
        endif
      enddo
 1    gln=aloggamma1(a-1.d0)
      gammaser=sum*exp(-x+a*log(x)-gln)/a
      return
      end
      
      real*8 function gammacf(a,x)
      implicit none
      integer*4 itmax,i
      real*8 a,x,eps,fpmin,a1,gln,b,c,d,h,an,del,aloggamma1
      parameter (itmax=300,eps=1.d-13,fpmin=1.d-30)
      if(a .eq. 0.d0)then
        gammacf=0.d0
        return
      endif
      a1=a-1.d0
      b=x-a1
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d) .lt. fpmin)then
          d=fpmin
        endif
        c=b+an/c
        if(abs(c) .lt. fpmin)then
          c=fpmin
        endif
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0) .lt. eps)then
          go to 1
        endif
      enddo
 1    gln=aloggamma1(a1)
      gammacf=exp(-x+a*log(x)-gln)*h
      return
      end

      real*8 function gamma0ser(x)
      implicit none
      integer*4 itmax,i
      real*8 x,eps,ap,sum,del
      parameter (itmax=300,eps=1.d-13)
      if(x .le. 0.d0)then
        gamma0ser=x/0.d0
        return
      endif
      ap=1.d0
      sum=x
      del=x
      do i=1,itmax
        del=-del*x*ap
        ap=ap+1.d0
        del=del/ap**2
        sum=sum+del
        if(abs(del) .lt. abs(sum)*eps)then
          go to 1
        endif
      enddo
 1    gamma0ser=sum
      return
      end
      
      real*8 function gamma0cf(x)
      implicit none
      integer*4 itmax,i
      real*8 x,eps,fpmin,a1,b,c,d,h,an,del
      parameter (itmax=300,eps=1.d-13,fpmin=1.d-30)
      a1=-1.d0
      b=x-a1
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do i=1,itmax
        an=-i**2
        b=b+2.d0
        d=an*d+b
        if(abs(d) .lt. fpmin)then
          d=fpmin
        endif
        c=b+an/c
        if(abs(c) .lt. fpmin)then
          c=fpmin
        endif
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0) .lt. eps)then
          go to 1
        endif
      enddo
 1    gamma0cf=exp(-x)*h
      return
      end

      complex*16 function cerf(z)
      implicit none
      real*8 a
      complex*16 z,cerfcf,cerfs,cerfcd
      a=dble(z)**2+imag(z)**2
      if(dble(z) .lt. 0.d0)then
        if(dble(-z) .gt. 1.d0)then
          cerf=-1.d0+cerfcf(-z)
        elseif(a .gt. 0.01d0)then
          cerf=-1.d0+cerfcd(-z)
        else
          cerf=-cerfs(-z)
        endif
      else
        if(dble(z) .gt. 1.d0)then
          cerf=1.d0-cerfcf(z)
        elseif(a .gt. 0.01d0)then
          cerf=1.d0-cerfcd(z)
        else
          cerf=cerfs(z)
        endif
      endif
      return
      end

      complex*16 function cerfc(z)
      implicit none
      real*8 a
      complex*16 z,cerfcf,cerfs,cerfcd
      a=dble(z)**2+imag(z)**2
      if(dble(z) .lt. 0.d0)then
        if(dble(z) .gt. 1.d0)then
          cerfc=2.d0-cerfcf(-z)
        elseif(a .gt. 0.01d0)then
          cerfc=2.d0-cerfcd(-z)
        else
          cerfc=1.d0+cerfs(-z)
        endif
      else
        if(dble(z) .gt. 1.d0)then
          cerfc=cerfcf(z)
        elseif(a .gt. 0.01d0)then
          cerfc=cerfcd(z)
        else
          cerfc=1.d0-cerfs(z)
        endif
      endif
      return
      end

      complex*16 function cerfs(z)
      use macmath
      implicit none
      real*8 r
c     Including m_2_sqrtpi:	2 / Sqrt[Pi]
      parameter (r=m_2_sqrtpi)
      complex*16 z,z2
      z2=z**2
      cerfs=r*z*(1.d0-z2*(1.d0/3.d0-z2*(1.d0/10.d0-
     $     z2*(1.d0/42.d0-z2*(1.d0/216.d0-z2/1320.d0)))))
      return
      end

      complex*16 function cerfcd(z)
      use macmath
      implicit none
c     Including m_2_sqrtpi:	2 / Pi
      real*8 r,h,an,an1,an2,an0
      parameter (r=m_2_pi,h=0.2d0)
      complex*16 z,cs,cs1,ca,ca2,cd,z1,z2
      integer*4 i,itmax
      parameter (itmax=66)
      real*8 a(itmax)
      data a /itmax*0.d0/
      if(a(1) .eq. 0.d0)then
        an=-1.d0
        do i=1,itmax
          an=an+2.d0
          a(i)=exp(-(an*h)**2)
        enddo
      endif
      if(imag(z) .gt. 0.d0)then
        z1=dcmplx(imag(z),dble(z))
        an0=anint(dble(z1)/2.d0/h)*2.d0
        z2=z1-an0*h
        an1=an0+1.d0
        an2=an0-1.d0
        ca=exp(2.d0*h*z2)
        ca2=ca**2
        cs=a(1)*(ca/an1+1.d0/an2/ca)
        do i=2,itmax
          an1=an1+2.d0
          an2=an2-2.d0
          ca=ca*ca2
          cd=a(i)*(ca/an1+1.d0/an2/ca)
          cs1=cs+cd
          if(cs1 .eq. cs)then
            go to 11
          endif
          cs=cs1
        enddo
        write(*,*)'cerfcd convergence error'
 11     cs=exp(z1**2-z2**2)*cs
        cerfcd=dcmplx(1.d0-r*imag(cs),-r*dble(cs))
      else
        z1=dcmplx(-imag(z),dble(z))
        an0=anint(dble(z1)/2.d0/h)*2.d0
        z2=z1-an0*h
        an1=an0+1.d0
        an2=an0-1.d0
        ca=exp(2.d0*h*z2)
        ca2=ca**2
        cs=a(1)*(ca/an1+1.d0/an2/ca)
        do i=2,itmax
          an1=an1+2.d0
          an2=an2-2.d0
          ca=ca*ca2
          cd=a(i)*(ca/an1+1.d0/an2/ca)
          cs1=cs+cd
          if(cs1 .eq. cs)then
            go to 10
          endif
          cs=cs1
        enddo
        write(*,*)'cerfcd convergence error'
 10     cs=exp(z1**2-z2**2)*cs
        cerfcd=dcmplx(1.d0-r*imag(cs),r*dble(cs))
      endif
      return
      end

      complex*16 function cerfcf(z)
      use macmath
      implicit none
c     Including m_2_sqrtpi:	2 / Sqrt[Pi]
      real*8 fpmin,eps,r
      parameter (fpmin=1.d-300,eps=1.d-15,r=m_2_sqrtpi)
      integer*4 i,itmax
      parameter (itmax=1000)
      real*8 an,a
      complex*16 z,cb,cd,ch,cdel,cc,z2
      z2=z**2
      cb=2.d0*z2+1.d0
      cc=1.d0/fpmin
      cd=1.d0/cb
      ch=cd
      an=-1.d0
      do i=1,itmax
        an=an+2.d0
        a=-an*(an+1.d0)
        cb=cb+4.d0
        cd=1.d0/(a*cd+cb)
        cc=cb+a/cc
        cdel=cc*cd
        ch=ch*cdel
        if(abs(cdel-1.d0) .lt. eps)then
          go to 2
        endif
      enddo
      write(*,*)'erfc convergence error'
 2    cerfcf=ch*z*r*exp(-z2)
      return
      end

      real*8 function productlog(x)
      use tfstk
      use macmath
      implicit none
c     Including m_e(Napier's constant: Exp[1])
      real*8 x,w,f1,d,eps,en
      parameter (eps=3.d-16,en=m_e)
      if(x .le. -1.d0/en)then
        productlog=-1.d0
        return
      endif
      if(x .eq. 0.d0)then
        productlog=0.d0
        return
      endif
      if(x .lt. -0.25d0)then
        w=sqrt1(2.d0*en*x+1.d0)
c        w=sqrt(2.d0*(en*x+1.d0))-1.d0
      elseif(x .lt. 2.d0)then
        w=.5d0*sqrt1(4.d0*x)
c        w=.5d0*(sqrt(4.d0*x+1.d0)-1.d0)
      else
        w=log(x/log(x/log(x)))
      endif
      f1=x*exp(-w)
      d=(f1-w)/(f1+1.d0)
      do while(abs(d) .gt. eps*abs(w))
        w=w+d
        f1=x*exp(-w)
        d=(f1-w)/(f1+1.d0)
      enddo
      productlog=w+d
      return
      end

      complex*16 function cproductlog(z)
      use macmath
      implicit none
c     Including m_e(Napier's constant: Exp[1])
      complex*16 z,w,f1,d
      real*8 productlog,eps,en
      parameter (eps=3.d-16,en=m_e)
      if(imag(z) .eq. 0.d0)then
        if(dble(z) .ge. -1.d0/en)then
          cproductlog=productlog(dble(z))
          return
        endif
      endif
      if(abs(z+1.d0/en) .lt. 0.2d0)then
        w=sqrt(2.d0*(en*z+1.d0))-1.d0
      elseif(abs(z) .lt. 3.d0)then
        w=.5d0*(sqrt(4.d0*z+1.d0)-1.d0)
      else
        w=log(z/log(z/log(z)))
      endif
      f1=z*exp(-w)
      d=(f1-w)/(f1+1.d0)
      do while(abs(d) .gt. abs(w)*eps)
        w=w+d
        f1=z*exp(-w)
        d=(f1-w)/(f1+1.d0)
      enddo
      cproductlog=w+d
      return
      end

      real*8 function inverseerf(x)
      implicit none
      real*8 x,erf,y,dy,ady0,y1
      real*8, parameter:: c=1.12837916709551257d0, a=0.140012d0,
     $     pa=0.461814d0
c      if(abs(x) .gt. 0.6d0)then
c        xl2=log((1.d0-x)*(1.d0+x))
c        xp=2.d0/pa+xl2*.5d0
c        y=sign(sqrt(-xl2/(sqrt(xp**2-xl2/a)+xp)/a),x)
c      else
        y=x
c      endif
      ady0=1.d100
      do while(.true.)
        dy=(x-erf(y))/c/exp(-y**2)
        if(abs(dy) .ge. ady0)then
          exit
        endif
        y1=y+dy
        if(y .eq. y1)then
          exit
        endif
        y=y1
        ady0=abs(dy)
      enddo
      inverseerf=y
      return
      end
