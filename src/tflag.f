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
        if(fname(i) == ' ')then
          cycle
        endif
        if(word == fname(i))then
          flags(i)=.true.
          exist=.true.
          return
        elseif(word == 'NO'//fname(i))then
          flags(i)=.false.
          exist=.true.
          return
        endif
        if(sino(i) == ' ')then
          cycle
        endif
        if(word == sino(i))then
          flags(i)=.false.
          exist=.true.
          return
        elseif(word == 'NO'//sino(i))then
          flags(i)=.true.
          exist=.true.
          return
        endif
      enddo
      exist=.false.
      return
      end
