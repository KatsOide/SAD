      subroutine twbuf(word,unit,lfno,ilmgn,irmgn,itab,icmd)
      implicit none
      integer*4 , intent(in)::lfno,ilmgn,irmgn,itab,icmd
      character*(*) , intent(in)::word,unit
      character*256 buff
      integer*4 ip,l,ip1,lw
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
        lw=len_trim(word)
        l=lw+len(unit)+1
        if(l .gt. 1)then
          if(ip1+l .gt. irmgn)then
            write(lfno,'(A)')buff(1:ip)
            ip=0
            go to 10
          endif
          buff(ip+1:ip1+1)=' '
          ip=ip1+l
          buff(ip1+2:ip)=word(1:lw)//unit
        endif
      endif
      return
      end

      subroutine twelm(lfno,idisp1,idisp2,title,irmgn,itab)
      use ffs
      use tffitcode
      use ffs_pointer
      implicit none
      integer*4 ,intent(in)::lfno,idisp1,idisp2,irmgn,itab
      character*(*) , intent(in)::title
      character*16 name,name1,title1
      call elname(idisp1,name)
      if(idisp2 .gt. 0)then
        call elname(idisp2,name1)
      else
        name1=' '
      endif
      title1=title
      call twbuf(
     1    title1(1:len_trim(title1))//' '//
     $     name(1:len_trim(name))//' '//name1,'',
     1    lfno,1,irmgn,itab,1)
      return
      end
