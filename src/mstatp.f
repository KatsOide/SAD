      subroutine mstatp(a,n,xmax,xmin,xrms,xmean,imax)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      if(n.eq.0) then
        xmax=0.d0
        xmin=0.d0
        xmean=0.d0
        xrms=0.d0
        return
      endif

      imax=1
      imin=1
      xmax=a(1)
      xmin=a(1)
      do i=2,n
        if(a(i).ge.xmax) then
          xmax=a(i)
          imax=i
        endif
        if(a(i).le.xmin) then
          xmin=a(i)
          imin=i
        endif
      enddo
      
      xmean=0.d0
      xrms=0.d0
      do i=1,n
        xmean=xmean+a(i)
        xrms=xrms+a(i)*a(i)
      enddo

      xrms=sqrt(xrms/float(n))
      xmean=xmean/float(n)
      if(xmax.lt.-xmin)then
        imax=imin
      endif
      return
      end
