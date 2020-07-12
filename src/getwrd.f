c     Peek/Get word from input buffer with case preserving
      integer*4 function peekwd0(outstr,next)
      use ffs_flag
      use tfcsi
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 ,intent(out):: next
      integer*4 is,notany,ifany
      next=ipoint
      peekwd0=0
c      write(*,*)'peekwd0 ',lfni,ipoint,lrecl,ios
      if(ipoint .gt. lrecl .or. ipoint .lt. 1)then
        outstr=' '
        return
      endif
      is=notany(buffer(1:lrecl),delim(3:ldel),ipoint)
      if(is .ge. ipoint)then
        if(buffer(is:is) .eq. char(10) .or.
     $       buffer(is:is) .eq. ';')then
          next=is
          outstr=' '
          return
        elseif(is .eq. lrecl)then
          next=lrecl+1
          outstr=buffer(ipoint:is)
          peekwd0=is-ipoint+1
          return
        endif
      else
        next=lrecl+1
        outstr=' '
        return
      endif
      next=ifany(buffer(1:lrecl),delim(1:ldel),is+1)
      if(next .le. 0)then
        next=lrecl+1
      endif
      outstr=buffer(is:next-1)
      peekwd0=next-is
      return
      end

      integer*4 function getwdl0(outstr)
      use tfcsi
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 next,peekwd0
      getwdl0=peekwd0(outstr,next)
c      write(*,*)'getwdl0 ',next,ipoint,lrecl
      call cssetp(next)
      return
      end

      integer*4 function getwrd0(outstr)
      use tfcsi
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 getwdl0
 1    getwrd0=getwdl0(outstr)
      if(outstr(1:1) .eq. ' ')then
        call skiplnget
        if(icsstat() .ne. 0)then
          outstr=' '
          return
        endif
        go to 1
      endif
      return
      end

c     Peek/Get word from input buffer with case normalization
      subroutine peekwd(outstr,next)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 ,intent(out):: next
      integer*4 l,peekwd0
      l=peekwd0(outstr,next)
      if(convcase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
      return
      end

      subroutine getwdl(outstr)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 l,getwdl0
      l=getwdl0(outstr)
      if(convcase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
      return
      end

      subroutine getwrd(outstr)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 l,getwrd0
      l=getwrd0(outstr)
      if(convcase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
c      write(*,*)'getwrd ',outstr(1:min(len(outstr),l))
      return
      end

c     Peek/Get word from input buffer with case preserving by preservecase
      subroutine peekwdp(outstr,next)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 ,intent(out):: next
      integer*4 l,peekwd0
      l=peekwd0(outstr,next)
      if(convcase .and. .not. preservecase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
      return
      end

      subroutine getwdlp(outstr)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 l,getwdl0
      l=getwdl0(outstr)
      if(convcase .and. .not. preservecase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
      return
      end

      subroutine getwrdp(outstr)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr
      integer*4 l,getwrd0
      l=getwrd0(outstr)
      if(convcase .and. .not. preservecase .and. l .gt. 0)then
        call capita1(outstr(1:min(len(outstr),l)))
      endif
      return
      end

c     Peek/Get word from input buffer with case normalization/preserved
c     1st arg(outstr)	case normalized word
c     2nd arg(outstrp)	case preserved word by preservecase flag
      subroutine peekwd2(outstr,outstrp,next)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr,outstrp
      integer*4 ,intent(out):: next
      integer*4 l,peekwd0
      l=peekwd0(outstrp,next)
      if(l .lt. 1)then
        outstr=' '
        return
      endif
      outstr=outstrp(1:min(len(outstrp),l))
      if(convcase)then
        call capita1(outstr(1:min(len(outstr),l)))
        if(.not. preservecase)then
          if(len(outstr) .ge. l)then
            outstrp=outstr(1:l)
          else
            call capita1(outstrp(1:min(len(outstrp),l)))
          endif
        endif
      endif
      return
      end

      subroutine getwdl2(outstr,outstrp)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr,outstrp
      integer*4 l,getwdl0
      l=getwdl0(outstrp)
      if(l .lt. 1)then
        outstr=' '
        return
      endif
      outstr=outstrp(1:min(len(outstrp),l))
      if(convcase)then
        call capita1(outstr(1:min(len(outstr),l)))
        if(.not. preservecase)then
          if(len(outstr) .ge. l)then
            outstrp=outstr(1:l)
          else
            call capita1(outstrp(1:min(len(outstrp),l)))
          endif
        endif
      endif
      return
      end

      subroutine getwrd2(outstr,outstrp)
      use ffs_flag
      use tmacro
      implicit none
      character*(*) ,intent(out):: outstr,outstrp
      integer*4 l,getwrd0
      l=getwrd0(outstrp)
      if(l .lt. 1)then
        outstr=' '
        return
      endif
      outstr=outstrp(1:min(len(outstrp),l))
      if(convcase)then
        call capita1(outstr(1:min(len(outstr),l)))
        if(.not. preservecase)then
          if(len(outstr) .ge. l)then
            outstrp=outstr(1:l)
          else
            call capita1(outstrp(1:min(len(outstrp),l)))
          endif
        endif
      endif
      return
      end

c     Peek character from input buffer with case normalization
      character function peekch(next)
      use ffs_flag
      use tmacro
      use tfcsi
      implicit none
      integer*4 ,intent(out):: next
      integer*4 notspace
      next=ipoint
      if(ipoint .ge. lrecl)then
        peekch=' '
        return
      endif
      next=notspace(buffer(1:lrecl),ipoint)
      if(next .eq. 0)then
        peekch=' '
        next=lrecl
      else
        peekch=buffer(next:next)
        if(convcase)then
          call capita1(peekch)
        endif
        if(peekch .eq. char(10))then
          peekch=' '
        endif
        next=next+1
      endif
      return
      end

      subroutine skiplnget
      use tfcsi
      implicit none
      integer*4 ifany,is
      if(ipoint .ge. ipbase .and. ipoint .le. lrecl)then
        is=ifany(buffer(1:lrecl),delim(1:2),ipoint)
        if(is .ge. ipoint .and. is .lt. lrecl)then
          ipoint=is+1
        else
          ipoint=lrecl+1
          call getbuf
        endif
      else
        ipoint=lrecl+1
        call getbuf
      endif
      return
      end

      subroutine skipline
      use tfcsi
      implicit none
      integer*4 ipoint1,ifchar
      if(ipoint .ge. lrecl)then
        ipoint=lrecl+1
      else
        ipoint1=ifchar(buffer(1:lrecl),char(10),ipoint)+1
        if(ipoint1 .le. 1)then
          ipoint=lrecl+1
        else
          ipoint=ipoint1
        endif
      endif
      return
      end

      subroutine csrst(lfn0)
      use tfcsi
      implicit none
      integer*4 lfn0
      lfn1=lfn0
      call skipline
      ios=0
      return
      end
