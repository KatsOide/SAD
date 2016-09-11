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
          do i=1,kwMAX-1
            if(kytbl(i,lt) .eq. ioff)then
              tfkwrd=pname(kytbl(i,0))(2:)
              return
            endif
          enddo
          tfkwrd='-'
        endif
      endif
      return
      end

      character*(*) function tfkwrd1(lt,ioff,ii)
      use tfstk
      use mackw
      implicit none
      integer*4 lt,ioff,i,ii
      ii=0
      if(ioff .eq. 0)then
        tfkwrd1=pname(kytbl(0,lt))(2:)
      else
        if(ioff .ge. kytbl(kwMAX,lt))then
          tfkwrd1=' '
        else
          do i=1,kwMAX-1
            if(kytbl(i,lt) .eq. ioff)then
              tfkwrd1=pname(kytbl(i,0))(2:)
              ii=i
              return
            endif
          enddo
          tfkwrd1='-'
        endif
      endif
      return
      end
