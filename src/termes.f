      subroutine termes(mess,mess1)
      use tfcsi,only:icslfnm
      implicit none
      character*(*) mess
      character*(*) mess1
      character*20 head
      integer*4 l1,l2,l
      l=icslfnm()
      if(l .eq. 0)then
        return
      endif
      l1=min(69,len_trim(mess))
      l2=min(70-l1,len_trim(mess1))
      if(mess(1:4) .eq. 'Info' .or. mess(1:5) .eq. 'Usage')then
        if(l .le. 0)then
          return
        endif
        head=' ???-FFS-'
      else
        head=' ???-FFS-Error-'
        call setffserr(1)
      endif
1     write(l,'(A,A,1X,A)',err=900)head(1:len_trim(head)),
     1       mess(1:l1),mess1(1:l2)
      if(l .ne. 6)then
        l=6
        go to 1
      endif
      return
900   write(*,*)'???-FFS-Error during message output to',l
      write(*,'(A,A,1X,A)')head(1:len_trim(head)),
     1       mess(1:l1),mess1(1:l2)
      return
      end

      subroutine setffserr(n)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 n
      iffserr=n
      return
      end
