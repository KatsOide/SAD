      subroutine  mdpmax(word,lfno)
      implicit real*8 (a-h,o-z)
      logical exist
      character*(*) word
      include 'inc/common.inc'
c     maximum dp/p shift attainable
      dps=getva(exist)
      if(exist) then
        dpshft=dps
      endif
      call getwdl(word)
      return
      end
