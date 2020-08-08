      subroutine tput(a,label,label1,form,lcol,lfno)
c      use tfstk, only:ktfenanq
      implicit none
      integer*4 lfno,i,j,lcol,ll,ll1
      real*8 a(lcol,6)
      character*(*) label(6),label1(lcol),form
      character*12 autofg
      character*80 buf
      ll1=len(label1(1))
      ll=len(label(1))
      buf(1:ll1)=' '
      if(lcol .eq. 6)then
        write(lfno,'(7A)')buf(1:ll1),label
      endif
      buf(ll1+1:ll1+1)=':'
      do 10 i=1,lcol
        buf(1:ll1)=label1(i)
        do 20 j=1,6
          buf(ll1+(j-1)*ll+2:ll1+j*ll+1)=' '//autofg(a(i,j),form)
c          write(*,*)autofg(a(i,j),form),' ',a(i,j),ktfenanq(a(i,j))
20      continue
        write(lfno,'(1X,A)')buf(1:ll1+ll*6+1)
10    continue
      if(lcol .eq. 6)then
        write(lfno,*)
      endif
      return
      end
