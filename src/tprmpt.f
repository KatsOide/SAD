      subroutine tprmpt(lfni,lfno,lfn1)
      use tfstk
      use ffs
      use ophash, only: opcode
      use tffitcode
      implicit none
      character*80 pr
      character*10 n,autofg
      integer*4 lfni,lfno,lfn1,l,nc
      integer*4 , save:: ipr=0
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

      subroutine tprmptget(ipr)
      implicit none
      integer*4 , intent(in)::ipr
      call tprmpt(ipr,-1,0)
      call getbuf
      call tprmpt(0,-1,0)
      return
      end
