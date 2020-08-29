      subroutine tfrenumber_element(isp1, kx, irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit none
      integer*8 kx, ia, ki
      integer(4), intent(in)  :: isp1
      integer(4), intent(out) :: irtc
c
      integer(4) :: itfmessage, ifany1
c
      character(len=MAXPNAME) :: name
      integer(4) :: namel
      real(8) :: rfromk
      integer :: index, count, n, last, status, i, j, ivi
      integer, allocatable :: forward_map(:)
c
      kx = ktfoper + mtfnull
      irtc = 0
c
c     Reset LINE$P table.  -- 2012/1/8 By KO
      call tfclearlinep()
c
c
      BODY: do
c     Checking number of arguments

         index = 0
         if(isp1 + 2 .ne. isp)exit BODY

c     Checking 1st argument
         index = -1
         call tfgetstrns(ktastk(isp1 + 1), name, namel)
         if(namel .le. 0)exit BODY
         if(ifany1(name(1:namel), namel, '*%{<|.', 1) .gt. 0)exit BODY

         count = 0
         do i=1,nlat
            if(name(1:namel) .eq. pnamec(i))then
               count = count + 1
            endif
         enddo
         if(count .lt. 1)exit BODY
         if(count .lt. 2)return

c     Checking 2nd argument
         index = -2
         if(ktfnonlistq(ktastk(isp1 + 2)))exit BODY

         ia = ktfaddr(ktastk(isp1 + 2))
         n = ilist(2, ia - 1)
         if(n .lt. 1)return

         allocate(forward_map(min(n,count)), STAT=status)
         if(status .ne. 0)then
            irtc = itfmessage(9, 'General::memoryfull', '""')
            return
         endif

c     Load sequence number list and check conflicts
         forward_map(1:min(n,count)) = 0
         last = 0
         do i=1,min(n, count)
           ki = klist(ia + i)
            if(ktfnonrealq(ki))exit BODY
            ivi = int(rfromk(ki))
            if(ivi .lt. 1)exit BODY
            if(ivi .gt. last)then
               last = ivi
            else
               do j=1,i-1
                  if(forward_map(j) .eq. ivi)exit BODY
               enddo
            endif
            forward_map(i) = ivi
         enddo

c     Execute renumbering
         j = 0
         do i=1,nlat
            if(name(1:namel) .eq. pnamec(i))then
               if(j .lt. min(n,count))then
                  j = j + 1 
                  ilist(i, ifmult) = forward_map(j)
               else
                  last = last + 1
                  ilist(i, ifmult) = last
               endif
            endif
         enddo
         return
      enddo BODY

      select case(index)
      case(0)
         irtc = itfmessage(9, 'General::narg', '"2"')

      case(-1)
         irtc = itfmessage(9, 'General::wrongtype',
     $        '"Name of component family for #1"')

      case(-2)
         irtc = itfmessage(9, 'General::wrongtype',
     $        '"List of positive sequence numbers for #2"')

      case default
         stop 'tfrenumber_element is reached to invalid error handler'

      end select
      return
      end subroutine tfrenumber_element
