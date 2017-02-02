      subroutine close_lun(lfn)
      implicit none
      integer*4 lfn
      close(lfn)
      return
      end

      integer*4 function open_read(fn, nc)
      implicit none
      character*(*) fn
      integer*4 nc
      integer*4 in, lun, nextfn
      
      write(*,*)'open_read ',fn(:nc)
      lun = -2
      in = nextfn(0)
      if(in .ne. 0)then
         open(in,file=fn(1:nc),status='OLD',
     $        err=9000)
         lun = in
      else
         lun = -1
      endif
 9000 continue
      open_read = lun
      return
      end

      integer*4 function open_write(fn, nc)
      implicit none
      character*(*) fn
      integer*4 nc
      integer*4 in, lun, nextfn

      lun = -2
      in = nextfn(0)
      if(in .ne. 0)then
         open(in,file=fn(1:nc),status='UNKNOWN',
     $        err=9000)
         lun = in
      else
         lun = -1
      endif
 9000 continue
      open_write = lun
      return
      end

      integer*4 function open_append(fn, nc)
      implicit none
      character*(*) fn
      integer*4 nc
      integer*4 in, lun, nextfn

      lun = -2
      in = nextfn(0)
      if(in .ne. 0)then
         open(in,file=fn(1:nc),status='UNKNOWN',access='APPEND',
     $        err=9000)
         lun = in
      else
         lun = -1
      endif
 9000 continue
      open_append = lun
      return
      end
