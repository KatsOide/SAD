      integer*4 function getfd(lfn)
      implicit none
      integer*4 lfn
      getfd = -1
      write(*,*) 'GETFD: UNIX file descriptor retrieving feature',
     $     ' is not supported!'
      stop 1
      return
      end
