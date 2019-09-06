      module prmpt
        integer*4 , save:: ipr=0
      end module

      subroutine getbuf
      use tfrbuf
      use tfcsi
      use iso_c_binding
      implicit none
      integer*4 lrecl0,lrecl00,nc
      logical*4 unmapped
      if(lfni .le. 0)then
        ios=99999
        return
      endif
      unmapped=itbuf(lfni) .le. moderead
2     if(ipoint .gt. lrecl)then
        call tprmpt(lfni,lfno,lfn1)
        if(unmapped)then
          if(.not. rec)then
            lrecl0=linep+1
          else
            lrecl0=lrecl+1
          endif
          ipoint=lrecl0
          lrecl0=max(lrecl0,1)
          lrecl00=lrecl0
 1        ios=-999
          if(lrecl0 .gt. nbmax-256)then
            ios=999999
            go to 10
          endif
          call tfreadbuf(lfni,lrecl0,nc)
          if(nc .eq. -99)then
            go to 20
          elseif(nc .eq.  -999)then
            go to 10
          elseif(nc .le. 0)then
            lrecl=max(lrecl0-1,0)
          else
            lrecl=len_trim(buffer(1:lrecl0+nc-1))
          endif
          if(lfn1 .gt. 0)then
            write(lfn1,'(1x,a)')buffer(lrecl0:lrecl)
          endif
          ipoint=lrecl00
          ios=0
          if(lrecl .lt. lrecl0 .or. buffer(lrecl0:lrecl) .eq. ' ')then
            lrecl=max(lrecl0-1,0)
          else
            if(buffer(lrecl:lrecl) .eq. '\\')then
              lrecl0=lrecl
              go to 1
            endif
            lrecl=lrecl+1
          endif
          if(lrecl .gt. 0)then
            buffer(lrecl:lrecl)=char(10)
          endif
        else
          ios=0
          call tfreadbuf(lfni,0,nc)
          if(nc .eq. -99)then
            go to 20
          elseif(nc .eq.  -999)then
            go to 10
          endif
          if(lfn1 .gt. 0)then
            if(buffer(lrecl:lrecl) .eq. char(10))then
              write(lfn1,'(1x,a)')buffer(ipoint:lrecl-1)
            else
              write(lfn1,'(1x,a)')buffer(ipoint:lrecl)
            endif
          endif
        endif
      else
        call skipline
        if(ipoint .gt. lrecl)then
          go to 2
        endif
      endif
      return
 10   if(ios .le. 0)then
        ios=9999
      endif
      if(lfni .eq. 5)then
        write(*,*)'???-getbuf-buffer overfolw for input stream'
        stop
      endif
      return
 20   if(ios .le. 0)then
        ios=99999
      endif
c      write(*,*)'getbuf-99999'
      if(lfni .eq. 5)then
        write(*,*)'???-getbuf-end of input stream'
        stop
      endif
      return
      end

      subroutine tprmpt(lfni,lfno,lfn1)
      use tfstk
      use ffs
      use ophash, only: opcode
      use tffitcode
      use prmpt
      implicit none
      character*80 pr
      character*10 n,autofg
      integer*4 lfni,lfno,lfn1,l,nc
      if(lfni .eq. 5 .and. lfno .eq. 6 .and. lfn1 .eq. 0)then
        if(ipr .eq. 0)then
          if(ffsprmpt)then
            call elname(mfpnt,pr)
            l=len_trim(pr)
            if(mfpnt .ne. mfpnt1)then
              pr(l+1:l+1)=':'
              call elname(mfpnt1,pr(l+2:80))
              l=len_trim(pr)
            endif
            pr(l+1:l+1)='/'
            call elname(id1,pr(l+2:80))
            l=len_trim(pr)
            if(id1 .ne. id2)then
              pr(l+1:l+1)=':'
              call elname(id2,pr(l+2:80))
              l=len_trim(pr)
            endif
            pr(l+1:l+1)='>'
            write(lfno,'(a,$)')pr(1:l+1)
          else
            n=autofg(rlist(iaxline)+1.d0,'S10.0')
            write(lfno,'('' In['',a,'']:= '',$)')n(1:len_trim(n))
          endif
        elseif(ipr .gt. 0)then
          nc=len_trim(opcode(ipr))
          pr(1:9)=' ...'//opcode(ipr)(1:nc)//'    '
          write(lfno,'(a,$)')pr(1:9)
        endif
      elseif(lfno .eq. -1)then
        ipr=lfni
      endif
      return
      end

      subroutine tprmptget(ipr1)
      use prmpt
      implicit none
      integer*4 , intent(in)::ipr1
      ipr=ipr1
      call getbuf
      ipr=0
      return
      end
