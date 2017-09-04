      subroutine tkptblini(kptbl)
      use tfstk
      use tmacro
      implicit none
      integer, parameter :: nkptbl = 6
      integer*4 kptbl(np0,nkptbl),i

      kptbl(1:np0,1:nkptbl)=0

      do i=1,np0
        kptbl(i,1)=i
        kptbl(i,2)=i
      enddo

      return
      end
