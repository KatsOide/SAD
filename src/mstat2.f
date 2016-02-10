      subroutine mstat2(a,la,n,xmax,xmin,xrms,xmean)
      implicit real*8 (a-h,o-z)
      dimension a(n),la(n)
      xmax=-1d36
      xmin= 1d36
      xmean=0d0
      xrms=0d0
      do 10 i=1,n
        j=la(i)
        xmax=max(xmax,a(j))
        xmin=min(xmin,a(j))
   10 continue
      do 20 i=1,n
        j=la(i)
        xmean=xmean+a(j)
        xrms=xrms+a(j)*a(j)
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
