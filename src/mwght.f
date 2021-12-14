      subroutine mwght(twiss,ip,a,b,n,m,nd,imon,nmon,cor)
      use ffs
      use tffitcode
      logical cor(8)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),imon(nmona,4),
     1          a(nd,m),b(n)
      include 'inc/common.inc'
c
      kcx=0
      kcy=0
      kdx=0
      kdy=0
      if( cor(1).or.cor(3) ) then
        kcx=nmon
      elseif( cor(2).or.cor(4) ) then
        kcy=nmon
      endif
      kc=kcx+kcy
      if( cor(5).or.cor(7) ) then
        kdx=nmon
      elseif( cor(6).or.cor(8) ) then
        kdy=nmon
      endif
      do 10 i=1,kcx
        l=imon(i,2)
        sqrbxi=sqrt(twiss(imon(l,1),ip,2))
        b(i)=b(i)*sqrbxi
        do 10 j=1,m
          a(i,j) = a(i,j)*sqrbxi
   10 continue
      do 20 i=kcx+1,kc
        l=imon(i-kcx,2)
        sqrbyi=sqrt(twiss(imon(l,1),ip,5))
        b(i)=b(i)*sqrbyi
        do 20 j=1,m
          a(i,j) = a(i,j)*sqrbyi
   20 continue
      do 30 i=kc+1,kc+kdx
        l=imon(i-kc,2)
        sqrbxi=sqrt(twiss(imon(l,1),ip,2))
        b(i)=b(i)*sqrbxi
        do 30 j=1,m
          a(i,j) = a(i,j)*sqrbxi
   30 continue
      k=kc+kdx
      do 40 i=k+1,k+kdy
        l=imon(i-k,2)
        sqrbyi=sqrt(twiss(imon(l,1),ip,5))
        b(i)=b(i)*sqrbyi
        do 40 j=1,m
          a(i,j) = a(i,j)*sqrbyi
   40 continue
      return
      end
