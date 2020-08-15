      subroutine tflag(word,exist)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer
      implicit none
      integer*4 i
      logical ,intent(out):: exist
      character*(*) ,intent(in):: word
      do i=1,nflag
        if(fname(i) .eq. ' ')then
          cycle
        endif
        if(word .eq. fname(i))then
          flags(i)=.true.
          exist=.true.
          return
        elseif(word .eq. 'NO'//fname(i))then
          flags(i)=.false.
          exist=.true.
          return
        endif
        if(sino(i) .eq. ' ')then
          cycle
        endif
        if(word .eq. sino(i))then
          flags(i)=.false.
          exist=.true.
          return
        elseif(word .eq. 'NO'//sino(i))then
          flags(i)=.true.
          exist=.true.
          return
        endif
      enddo
      exist=.false.
      return
      end
