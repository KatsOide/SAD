      subroutine tmast
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 i,ie,id,j
      do 10 i=1,nlat
        master(i)=0
10    continue
      ie=0
      do 1010 i=1,nlat-1
        if(i .le. ie)then
          go to 1010
        endif
        id=idtypec(i)
        if(id .eq. icSOL)then
          go to 1010
        elseif(id .eq. 1 .or. id .gt. icDODECA .and.
     $         id .ne. icMULT .and. id .ne. icCAVI .and.
     $         id .ne. icTCAV)then
          go to 1010
        else
          ie=i
          do 1020 j=i+1,nlat-1
            if(iele1(j) .eq. iele1(i))then
              ie=j
              master(j)=-i
            elseif(pos(j) .ne. pos(j+1))then
              go to 1021
            endif
 1020     continue
 1021     master(i)=ie
        endif
 1010 continue
      return
      end
