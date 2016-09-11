      character*(*) function tfkwrd(lt,ioff)
      use tfstk
      use mackw
      implicit none
      integer*4 lt,ioff,i
      if(ioff .eq. 0)then
        tfkwrd=pname(kytbl(0,lt))(2:)
      else
        if(ioff .ge. kytbl(kwMAX,lt))then
          tfkwrd=' '
        else
          i=kyindex(ioff,lt)
          if(i .ne. 0)then
            tfkwrd=pname(kytbl(i,0))(2:)
          else
            tfkwrd='-'
          endif
        endif
      endif
      return
      end
