      subroutine tday(day)
      implicit none
c      external idate
      integer*4 month,iday,iy,iw,iv(8)
      character*9 week(0:6)
      character*(*) day
      data week /'Sunday   ','Monday   ','Tuesday  ',
     1           'Wednesday','Thursday ','Friday   ',
     1           'Saturday '/
c      call idate(iarray)
c      iday=iarray(1)
c      month=iarray(2)
c      iy=iarray(3)
      call datetime1(iv)
      iy=iv(1)
      month=iv(2)
      iday=iv(3)
      if(month .lt. 3)then
        iy=iy-1
        month=month+12
      endif
      iw=int(iy*1.25)-int(iy*.01)+int(iy/400.)
     1   +int(2.6*month+1.601)+iday
      iw=mod(iw,7)
      day=week(iw)
      return
      end

      character*2 function ord(m)
      implicit none
      integer*4 m,mm
      mm=mod(m,100)
      if(mm .gt. 3 .and. mm .lt. 21)then
        ord='th'
      else
        mm=mod(m,10)
        if(mm .eq. 1)then
          ord='st'
        elseif(mm .eq. 2)then
          ord='nd'
        elseif(mm .eq. 3)then
          ord='rd'
        else
          ord='th'
        endif
      endif
      return
      end

      character function rad62(i)
      implicit none
      integer*4 ,intent(in):: i
      character*62 ,parameter ::c62=
     1'0123456789AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz'
      if(i .ge. 0 .and. i .le. 61)then
        rad62=c62(i+1:i+1)
      else
        rad62='*'
      endif
      return
      end
