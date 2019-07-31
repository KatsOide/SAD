      subroutine  pack(ipoint,iflag,nlive,ndim)
      parameter (large=2**19)
      dimension ipoint(ndim),iflag(*)
      n=0
      do 11 i=1,ndim
        j=ipoint(i)
        if(iflag(j).eq.0) then
          n=n+1
        else
          ipoint(i)=ipoint(i)+large
        endif
   11 continue
      call msort(ipoint,ndim)
      do 12 i=n+1,ndim
        ipoint(i)=ipoint(i)-large
   12 continue
      nlive=n
      return
      end
