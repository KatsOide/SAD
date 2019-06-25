      integer*4 function getfd(lfn)
      implicit none
      integer*4 lfn, fd, ierr
      external PXFFILENO
      call PXFFILENO(lfn, fd, ierr)
      getfd = fd
      return
      end
