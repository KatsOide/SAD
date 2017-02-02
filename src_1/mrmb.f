      subroutine  mrmb(is,xs,b,ncb,nbp,a,idim)
      implicit real*8 (a-h,o-z)
      dimension is(ncb+4,*),xs(ncb+2,2,*),a(idim,*)
     z,         b(idim,*)
      call pclr(b,2*idim*nbp)
      do 12 j=1,nbp
        ncor=is(2,j)
        j1=2*j-1
        do 11 i=1,idim
          do 10 k=1,ncor
            b(i,j1  )=b(i,j1  )+a(i,is(k+2,j))*xs(k,1,j)
            b(i,j1+1)=b(i,j1+1)+a(i,is(k+2,j))*xs(k,2,j)
   10     continue
   11   continue
   12 continue
      return
      end
