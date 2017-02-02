      integer*4 function fseek_subroutine(lun, offset, whence)
      implicit none
      integer*4 lun, offset, whence
      call fseek(lun, offset, whence, fseek_subroutine)
      end
