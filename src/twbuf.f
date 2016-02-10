      subroutine twbuf(word,lfno,ilmgn,irmgn,itab,icmd)
      character*(*) word
      character*255 buff
      integer*4 ip
      save buff,ip
      if(icmd .eq. 0)then
        ip=0
      elseif(icmd .eq. -1)then
        if(ip .gt. ilmgn-1)then
          write(lfno,'(A)')buff(1:ip)
        endif
        ip=0
      else
10      ip1=(max(ip,ilmgn-1)+itab-1)/itab*itab
        if(ip1 .ge. irmgn)then
          write(lfno,'(A)')buff(1:ip)
          ip=0
          go to 10
        endif
        l=len_trim(word)+1
        if(l .gt. 1)then
          if(ip1+l .gt. irmgn)then
            write(lfno,'(A)')buff(1:ip)
            ip=0
            go to 10
          endif
          buff(ip+1:ip1+1)=' '
          buff(ip1+2:)=word(1:l-1)
          ip=ip1+l
        endif
      endif
      return
      end
