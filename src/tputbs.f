      subroutine tputbs(a,label,lfno)
      implicit none
      real*8 ,intent(in):: a(21)
      integer*4 ,intent(in):: lfno
      integer*4 i,j
      character*(*) ,intent(in):: label(6)
      character*60 vout
      character*10 autofg
      write(lfno,'(10X,6A)')label
      do 10 i=1,6
        vout=' '
        do j=1,i
          vout((j-1)*10+1:j*10)=autofg(a(i*(i-1)/2+j),'10.7')
        enddo
        write(lfno,9015)label(i),vout(1:i*10)
 9015   format(1X,A,a)
 10   continue
      write(lfno,*)
      return
      end
