      real*8 function xsin(x)
      implicit none
      real*8 x,x2
      if(abs(x) .gt. .1d0)then
        xsin=x-sin(x)
      else
        x2=x**2
        xsin=x*x2/6.d0*(1.d0-x2/20.d0*(1.d0-x2/42.d0*(
     1       1.d0-x2/72.d0*(1.d0-x2/110.d0))))
      endif
      return
      end

      real*8 function xsinh(x)
      implicit none
      real*8 x,x2
      if(abs(x) .gt. .1d0)then
        xsinh=x-sinh(x)
      else
        x2=x**2
        xsinh=-x*x2/6.d0*(1.d0+x2/20.d0*(1.d0+x2/42.d0*(
     1       1.d0+x2/72.d0*(1.d0+x2/110.d0))))
      endif
      return
      end

      real*8 function xlog(x)
      implicit none
      real*8 x
      if(abs(x) .gt. 1.d-2)then
        xlog=log(1.d0+x)
      else
        xlog=x*(1.d0-x*(.5d0-x*(1.d0/3.d0
     1      -x*(.25d0-x*(.2d0-x*(1.d0/6.d0
     1      -x*(1.d0/7.d0-x*(.125d0-x/9.d0))))))))
      endif
      return
      end

      real*8 function sinc(x)
      implicit none
      real*8 x,x2
      if(abs(x) .gt. .1d0)then
        sinc=x*cos(x)-sin(x)
      else
        x2=x**2
        sinc=-x*x2/6.d0*(2.d0-x2/20.d0*(4.d0-x2/42.d0*(
     1       6.d0-x2/72.d0*(8.d0-x2/11.d0))))
      endif
      return
      end

      real*8 function sinhc(x)
      implicit none
      real*8 x,x2
      if(abs(x) .gt. .1d0)then
        sinhc=x*cosh(x)-sinh(x)
      else
        x2=x**2
        sinhc=x*x2/6.d0*(2.d0+x2/20.d0*(4.d0+x2/42.d0*(
     1       6.d0+x2/72.d0*(8.d0+x2/11.d0))))
      endif
      return
      end
