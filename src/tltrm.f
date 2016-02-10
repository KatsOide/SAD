      subroutine tltrm(latt,kptbl)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 i
      integer*4 latt(2,nlat),kptbl(np0,6)
      call tclrpara(latt,nlat-1)
      do i=1,np0
c        write(*,*)'tltrm ',i,np0,kptbl(i,3)
c     $       ilist(1,abs(kptbl(i,3)))
        if(kptbl(i,3) .ne. 0)then
          call tfree(int8(abs(kptbl(i,3))))
          kptbl(i,3)=0
        endif
      enddo
      return
      end
