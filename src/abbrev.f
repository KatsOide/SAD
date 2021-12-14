      logical*4 function abbrev(word,ref,delim)
      implicit none
      character*(*) word,ref
      character delim
      integer*4 id,lenw,lw,lr
      id=index(ref,delim)
      if(id .le. 0)then
        abbrev=word .eq. ref
        return
      endif
      abbrev=.false.
      if(word(1:id-1) .eq. ref(1:id-1))then
        lw=lenw(word)
        lr=lenw(ref)
        if(lw .gt. lr-1)then
          return
        endif
        if(lw .ge. id)then
          if(word(id:lw) .ne. ref(id+1:lw+1))then
            return
          endif
        endif
      else
        return
      endif
      abbrev=.true.
      return
      end
