      subroutine tdinit(lfnd1,cmd1,style1)
      implicit none
      integer*4 n,icol,lfnd,lfnd1
      character*(*) style1,cmd1
      character*80 buf
      character*12 cmd,pat
      common /td/ n,icol,lfnd
      common /tdb/ buf,cmd,pat
      lfnd=lfnd1
      cmd=cmd1
      n=0
      icol=0
      buf=' '
      write(lfnd,*)' SET CARD 80'
      if(style1 .ne. ' ')then
        pat='PATTERNED'
        write(lfnd,*)' SET PATTERN ',style1
      else
        pat=' '
      endif
      return
      end

      subroutine tdterm
      implicit none
      integer*4 n,icol,lfnd
      character*80 buf
      character*12 cmd,pat
      common /td/ n,icol,lfnd
      common /tdb/ buf,cmd,pat
      if(icol .ne. 0)then
        write(lfnd,'(A)')buf(1:icol)
        write(lfnd,9001)cmd,pat
9001    format(1x,a,2x,a)
      endif
      return
      end

      subroutine tdput(x0,y0)
      implicit none
      integer*4 l,n,icol,lfnd,lene,ifany
      real*8 x,y,x0,y0
      character*11 autofg,bx,by
      character*23 b
      character*80 buf
      character*12 cmd,pat
      external trim
      common /td/ n,icol,lfnd
      common /tdb/ buf,cmd,pat
      x=x0
      y=y0
      if(abs(x) .lt. 1.d-30)then
        x=0.d0
      endif
      if(abs(y) .lt. 1.d-30)then
        y=0.d0
      endif
      bx=autofg(x,'S9.6')
      if(ifany(bx,'E',1) .ne. 0)then
        bx=autofg(x,'S11.8')
      endif
      by=autofg(y,'S9.6')
      if(ifany(by,'E',1) .ne. 0)then
        by=autofg(y,'S11.8')
      endif
      b=bx//' '//by
      call trim(b)
      l=lene(b)
1     if(icol+l+1 .le. 79)then
        if(icol .eq. 0)then
          buf(1:l)=b
          icol=l
        else
          buf(icol+1:icol+l+1)=';'//b
          icol=icol+l+1
        endif
      else
        write(lfnd,'(A)')buf(1:icol)
        icol=0
        go to 1
      endif
      n=n+1
      if(n .gt. 1000 .and. icol .ne. 0)then
        write(lfnd,'(A)')buf(1:icol)
        write(lfnd,9001)cmd,pat
9001    format(1x,a,2x,a)
        n=0
      endif
      return
      end
