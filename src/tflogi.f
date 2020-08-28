      logical*4 function tflogi(word,exist)
      use ffs
      implicit none
      integer*4 i
      character*(*) ,intent(in):: word
      character*32 fname1
      logical*4 ,intent(out):: exist
      tflogi=.false.
      exist=.true.
      if(word .eq. ' ')then
        exist=.false.
        tflogi=.false.
      endif
      do i=1,nflag
        if(word .eq. fname(i))then
          tflogi=flags(i)
          return
        else
          fname1='NO'//fname(i)
          if(word .eq. fname1)then
            tflogi=.not. flags(i)
            return
          elseif(word .eq. sino(i))then
            tflogi=.not. flags(i)
            return
          else
            fname1='NO'//sino(i)
            if(word .eq. fname1)then
              tflogi=flags(i)
              return
            endif
          endif
        endif
      enddo
      exist=.false.
      return
      end
