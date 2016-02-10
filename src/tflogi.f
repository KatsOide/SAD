      logical*4 function tflogi(word,flg,fname,sino,nflag,exist)
      implicit none
      integer*4 nflag,i
      character*(*) word,fname(nflag),sino(nflag)
      character*32 fname1
      logical*4 flg(nflag),exist
      tflogi=.false.
      exist=.true.
      if(word .eq. ' ')then
        exist=.false.
        tflogi=.false.
      endif
      do i=1,nflag
        if(word .eq. fname(i))then
          tflogi=flg(i)
          return
        else
          fname1='NO'//fname(i)
          if(word .eq. fname1)then
            tflogi=.not. flg(i)
            return
          elseif(word .eq. sino(i))then
            tflogi=.not. flg(i)
            return
          else
            fname1='NO'//sino(i)
            if(word .eq. fname1)then
              tflogi=flg(i)
              return
            endif
          endif
        endif
      enddo
      exist=.false.
      return
      end
