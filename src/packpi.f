      subroutine  packpi(ipair,npair,ia,ipoint,iflag,nlive,ndim,packP)
      logical mhogan,packP
      dimension ipair(2,npair),ia(ndim),ipoint(ndim),iflag(ndim)
      external pack
      if(packP) then
        do 11 i=1,ndim
          j=ipoint(i)
          do 10 jp=1,npair
            if(mhogan(ipair(1,jp),ipair(2,jp)+1,ia(j))) then
              iflag(j)=iflag(j)+1
              goto 11
            endif
   10     continue
   11   continue
      else
        do 21 i=1,ndim
          j=ipoint(i)
          do 20 jp=1,npair
            if(mhogan(ipair(1,jp),ipair(2,jp)+1,ia(j))) then
              iflag(j)=iflag(j)-1
              goto 21
            endif
   20     continue
   21   continue
      endif
      call pack(ipoint,iflag,nlive,ndim)
      return
      end
