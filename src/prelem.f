      subroutine prelem(idxe,head)
      use maccbk
      implicit real*8 (a-h,o-z)
      integer idxe
      character*(*) head
c
      call prelm0(idxe,idval(idxe),head)
      return
      end
