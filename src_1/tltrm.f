      subroutine tltrm(latt,kptbl)
      use tfstk
      use tmacro
      implicit none
      integer*4 i
      integer*8 latt(nlat)
      integer*4 kptbl(np0,6)
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
