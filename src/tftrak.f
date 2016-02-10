      subroutine tftrak(word,trdtbl,trval,lfno,exist)
      use ffs
      use tffitcode
      real*8 trdtbl(3,6)
      logical*4 exist,abbrev
      character*(*) word
c     begin initialize for preventing compiler warning
      j=0
c     end   initialize for preventing compiler warning
      exist=.false.
      if(abbrev(word,'TRA_CK','_'))then
        if(dapert)then
          if(np0 .lt. 3)then
            call termes(lfno,'NP must be larger than 2.',' ')
            call skipline
            exist=.false.
            return
          endif
          call tfsetparam
          call tftrad(trdtbl,trval,lfno )
          exist=.true.
        else
          write(*,*)' Only DPAERT mode has been implemented.'
          return
        endif
      elseif(word .eq. 'NXR' .or. word .eq. 'NYR' .or.
     $       word .eq. 'NZR')then
        if(word .eq. 'NXR')then
          j=1
        elseif(word .eq. 'NYR')then
          j=3
        elseif(word .eq. 'NZR')then
          j=6
        endif
        do i=1,3
          call tfgetr(trdtbl(i,j),1.d0,word,lfno,exist)
          if(.not. exist)then
            return
          endif
        enddo
      endif
      return
      end
