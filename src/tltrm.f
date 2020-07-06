      subroutine tltrm(kptbl)
      use tfstk
      use tmacro
      use sad_main
      use tparastat
      implicit none
      integer*4 i
      integer*4 kptbl(np0,6)
      call tclrpara
      do i=1,np0
c         write(*,*)'tltrm ',i,np0,kptbl(i,3)
c     $       ilist(1,abs(kptbl(i,3)))
        if(kptbl(i,3) .ne. 0)then
          call tfree(int8(abs(kptbl(i,3))))
          kptbl(i,3)=0
        endif
      enddo
      return
      end
