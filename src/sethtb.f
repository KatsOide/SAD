      function sethtb(token,type,ival)
      use maccbk
      implicit none
      character*(*) token
      integer*4 sethtb,type,ival,decloc
      integer idx,hsrch,lenw
      include 'inc/MACCODE.inc'
      sethtb=0
c
      idx= hsrch(token(:lenw(token)))
c       print *,"sethtb: ",token(:lenw(token)),lenw(token),
c     $      type,ival,idx
       if(idx .le. 0 .or. idx .gt. HTMAX) then
          call errmsg('sethtb'
     &               ,'illegal index value for sethashtble'
     &               , 0,16)
       else
         idtype(idx)=type
         if(type .eq. icRSVD) then
            idval(idx)=decloc(ival)
         else
            idval(idx)=ival
         endif
         sethtb=idx
       endif

       return

       end
