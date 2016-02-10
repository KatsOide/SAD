      subroutine mstat(a,n,xmax,xmin,xrms,xmean)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      xmax=-1d36
      xmin= 1d36
      xmean=0d0
      xrms=0d0
      do 10 i=1,n
        xmax=max(xmax,a(i))
        xmin=min(xmin,a(i))
   10 continue
      do 20 i=1,n
        xmean=xmean+a(i)
        xrms=xrms+a(i)*a(i)
   20 continue
      if(n.ne.0) then
        xrms=sqrt(xrms/n)
        xmean=xmean/n
      else
        xmax=0d0
        xmin=0d0
      endif
      return
      end
