      subroutine prelm0(idxe,idxl,head)
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      implicit real*8 (a-h,o-z)
      integer idxe,idxl
      character*(*) head
c
c for debug
c     print *,'>enter into prelem',pname(idxe)
c end debug
      if (idtype(idxe) .eq. icNULL) then
        write(outfl,*)
     &                 pname(idxe),':not defined yet'
      else if (idtype(idxe) .ge. icMXEL) then
        call errmsg('prelem'
     &       ,'invalid ID name '//pname(idxe)//'!',0,0)
      else
        call prkwnm(
     &        pname(kytbl(0,idtype(idxe)))(2:),
     &        pname(idxe))
        do 1000 i=1,kwMAX-2
        if (kytbl(i,idtype(idxe)) .gt. 0) then
          call prkwdv(head,pname(kytbl(i,0))(2:),
     &           rlist(idxl+kytbl(i,idtype(idxe))))
        endif
 1000   continue
c flush print buffer
        call prkwot(');')
c print error list (if exist)
        idxerr=ilist(2,idxl)
        nerr=0
        if (idxerr .ne. 0) then
          call prkwnm(
     &          pname(kytbl(0,idtype(idxe)))(2:),
     &          pname(idxe))
 2100     continue
            if (idxerr .eq. 0) go to 2000
            if (kytbl(ilist(1,idxerr),idtype(idxe)) .gt. 0) then
              nerr=nerr+1
c             write(cerr,'(I8.8)')nerr
c             call prkw('!'//cerr//'-th error list for ')
              call prkw(pname(kytbl(ilist(1,idxerr),0))(2:)//'=(')
              if (ilist(2,idxerr+1) .le. 0) then
                call Lprnum(ilist(1,idxerr+1),idxerr+2,head)
              else
c               write(outfl,*)
c    &                       pname(ilist(2,idxerr+1))(2:)
                call prkw(   pname(ilist(2,idxerr+1))(2:))
                call Lprnum(ilist(1,idxerr+1),idxerr+2,head)
              endif
              call prkw(')')
            endif
            idxerr=ilist(2,idxerr)
            go to 2100
 2000     continue
          call prkwot(');')
        endif
      endif
c for debug
c     print *,'<leave from prelem',pname(idxe)
c end debug
      return
      end
