      real*8 function getvad(exist,word,def)
      implicit none
      real*8 def,getva
      logical*4 exist
      character*(*) word
      getvad=getva(exist)
      if(.not. exist)then
        call getwdl(word)
        if(word .eq. '*')then
          getvad=def
          exist=.true.
        endif
      endif
      return
      end

      real*8 function getva(exist)
      use tfcsi
      implicit none
      logical*4 exist
      real*8 getval
      getva=getval()
      exist=ios .eq. 0
      if(ios .eq. -2)then
        ios=0
      endif
      return
      end

      real*8 function getval()
      use tfstk
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      integer*4 irtype,itfpeeko,next
      irtype=itfpeeko(kx,next)
      if(irtype .eq. 0 .and. ktfrealqd(kx,getval))then
        call cssetp(next)
      else
        getval=0.d0
        call cssets(-2)
      endif
      return
      end
