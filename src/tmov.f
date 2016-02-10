      subroutine tmov(a,b,n)
      implicit none
      integer*4 n
      real*8 a(n),b(n)
      b=a
      return
      end

      subroutine tmovi(ia,ib,n)
      implicit none
      integer*4 n
      integer*4 ia(n),ib(n)
      ib=ia
      return
      end

      subroutine tadd(a,b,c,n)
      implicit none
      integer*4 n
      real*8 a(n),b(n),c(n)
      c=a+b
      return
      end

      subroutine tsub(a,b,c,n)
      real*8 a(n),b(n),c(n)
      c=a-b
      return
      end

      subroutine tfill(a,x,n)
      implicit none
      integer*4 n
      real*8 a(n),x
      a=x
      return
      end

      subroutine tclr(a,n)
      implicit none
      integer*4 n
      real*8 a(n)
      a=0.d0
      return
      end

      subroutine tmul(a,c,b,n)
      implicit none
      integer*4 n
      real*8 a(n),c,b(n)
      b=c*a
      return
      end
