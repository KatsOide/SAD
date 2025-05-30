      real*8 function getvad(exist,word,def)
      implicit none
      real*8 ,intent(in):: def
      real*8 getva
      logical*4 ,intent(out):: exist
      character*(*) ,intent(out):: word
      getvad=getva(exist)
      if(.not. exist)then
        call getwdl(word)
        if(word == '*')then
          getvad=def
          exist=.true.
        endif
      endif
      return
      end

      real*8 function getva(exist)
      use tfcsi
      implicit none
      logical*4 ,intent(out):: exist
      real*8 getval
      getva=getval()
      exist=ios == 0
      if(ios == -2)then
        ios=0
      endif
      return
      end

      real*8 function getval()
      use geto
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      integer*4 irtype,next
      irtype=itfpeeko(kx,next)
      if(irtype == 0 .and. ktfrealq(kx,getval))then
        ipoint=next
      else
        getval=0.d0
        ios=-2
      endif
      return
      end
