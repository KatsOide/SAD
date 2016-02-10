      subroutine tfgeti(ival,scale,word,lfno,exist)
      implicit none
      real*8 val,scale
      integer*4 lfno,ival
      character*(*) word
      logical exist
      call tfgetr(val,scale,word,lfno,exist)
      if(exist)then
        ival=val
      endif
      return
      end
