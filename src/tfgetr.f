      subroutine tfgetr(val,scale,word,lfno,exist)
      implicit none
      real*8 val,getva,scale,v
      integer*4 lfno
      character*(*) word
      logical exist
      v=getva(exist)*scale
      if(exist)then
        val=v
      else
        write(lfno,*)'?Missing parameter for ',word(1:len_trim(word))
        call skipline
      endif
      return
      end
