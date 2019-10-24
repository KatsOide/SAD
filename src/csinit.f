      subroutine csinit(lfn0,iconv1,cmnt1)
      use tfstk
      use tfrbuf
      use tfcsi
      use iso_c_binding
      implicit none
      integer*4 lfn0,iconv1
      character*(*) cmnt1
c      allocate(buffer0)
      buffer=>buffer0
      ibcloc=transfer(c_loc(buffer0(1:1)),ibcloc)/8
c      write(*,*)'csinit-0 ',ibcloc*8,
c     $     transfer(c_loc(jlist(1,ibcloc)),ibcloc)
      ibuf(lfni)=ibcloc
      lfn1=lfn0
      ipoint=>mbuf(lfni)
      lrecl=>lbuf(lfni)
      lrecl=0
      buffer(1:1)=char(10)
      ipoint=1
      iconv=iconv1
c      call tfsetconvcase(iconv1 .eq. 1)
      delim=';'//char(10)//' ,'//char(9)
      ldel=len_trim(delim)
      cmnt=cmnt1
      lcmnt=len_trim(cmnt)
      rep=.false.
      lastln=0
      ibegt=0
      lastt=0
      call csrst(lfn0)
c      write(*,*)'csinit ',lfni,ipoint,lrecl,ibcloc,buffer(1:1)
      return
      end

      integer*4 function lene(string)
      implicit none
      character*(*) string
      lene=len_trim(string)
      return
      end

      integer*4 function lene1(string,l)
      implicit none
      integer*4 l
      character*(*) string
      lene1=len_trim(string(:l))
      return
      end

      integer*4 function lenw(string)
      implicit none
      character*(*) string
      lenw=index(string,' ')-1
      if(lenw .lt. 0)then
        lenw=len(string)
      endif
      return
      end

      integer*4 function ifany(string,pat,istart)
      implicit none
      character*(*) string,pat
      integer*4 istart,ifany1
      ifany=ifany1(string,len_trim(string),
     $     pat(:len_trim(pat)),istart)
      return
      end

      integer*4 function ifany1(string,l,pat,istart)
      implicit none
      integer*4 istart,i,l
      character*(*) pat
      character string(l)
      if(istart .gt. 0)then
        do i=istart,l
          if(index(pat,string(i)) .gt. 0)then
            ifany1=i
            return
          endif
        enddo
      endif
      ifany1=0
      return
      end

      integer*4 function notany(string,pat,istart)
      implicit none
      integer*4 istart,notany1
      character*(*) string,pat
      notany=notany1(string,len_trim(string),
     $     pat(:len_trim(pat)),istart)
      return
      end

      integer*4 function notany1(string,l,pat,istart)
      implicit none
      integer*4 i,l,istart
      character*(*) string,pat
      if(istart .gt. 0)then
        do i=istart,l
          if(index(pat,string(i:i)) .le. 0)then
            notany1=i
            return
          endif
        enddo
      endif
      notany1=0
      return
      end

      integer*4 function indexs(string,pat,istart)
      implicit none
      character*(*) string,pat
      integer*4 istart,l,m
      l=len_trim(string)
      m=len_trim(pat)
      indexs=0
      if(istart .gt. 0 .and. istart .le. l)then
        indexs=index(string(istart:l),pat(1:m))
        if(indexs .gt. 0)then
          indexs=indexs+istart-1
        endif
      endif
      return
      end

      integer*4 function indexb(string,ns,pat,np,istart)
      implicit none
      integer*1 string(*),pat(*)
      integer*4 ns,np,istart
      integer*4 i,j

      do indexb=istart,ns-np+1
         i=indexb
         j=1
         do while(j .le. np .and. string(i) .eq. pat(j))
            i=i+1
            j=j+1
         enddo
         if(j .gt. np) then
            return
         endif
      enddo
      indexb=0
      return
      end

      subroutine trim(string)
      implicit none
      character*(*) string
      logical tr
      integer*4 l,len,j,i
      l=len(string)
      j=1
      tr=.true.
      do i=1,l
        if(tr)then
          tr=string(i:i) .eq. ' '
          if(.not. tr)then
            string(j:j)=string(i:i)
            j=j+1
          endif
        else
          tr=string(i:i) .eq. ' '
          string(j:j)=string(i:i)
          j=j+1
        endif
      enddo
      string(j:l)=' '
      return
      end

      subroutine capita(string)
      implicit none
      character*(*) string
      integer*4 i,ic,ioff
c      parameter (ioff=ichar('A')-ichar('a'))
      parameter (ioff=-32)
      do i=1,len(string)
        ic=ichar(string(i:i))
        if(ic .ge. ichar('a') .and. ic .le. ichar('z'))then
          string(i:i)=char(ic+ioff)
        endif
      enddo
      return
      end

      subroutine small(string)
      implicit none
      character*(*) string
      integer*4 i,ic,ioff
c      parameter (ioff=ichar('A')-ichar('a'))
      parameter (ioff=-32)
      do i=1,len(string)
        ic=ichar(string(i:i))
        if(ic .ge. ichar('A') .and. ic .le. ichar('Z'))then
          string(i:i)=char(ic-ioff)
        endif
      enddo
      return
      end

      integer*4 function notchar(string,ch,istart)
      implicit none
      character*(*) string
      character ch
      integer*4 istart,notchar1
      notchar=notchar1(string,len_trim(string),ch,istart)
      return
      end

      integer*4 function notchar1(string,l,ch,istart)
      implicit none
      character*(*) string
      character ch
      integer*4 istart,i,l
      if(istart .gt. 0)then
        do i=istart,l
          if(string(i:i) .ne. ch)then
            notchar1=i
            return
          endif
        enddo
      endif
      notchar1=0
      return
      end

      integer*4 function notspace(string,istart)
      implicit none
      character*(*) string
      integer*4 istart,i
      if(istart .gt. 0)then
        do i=istart,len(string)
          if(string(i:i) .ne. ' ')then
            notspace=i
            return
          endif
        enddo
      endif
      notspace=0
      return
      end

      integer*4 function ifchar(string,ch,istart)
      implicit none
      character*(*) string
      character ch
      integer*4 istart,l
      if(istart .gt. 0)then
        l=len_trim(string)
        ifchar=index(string(istart:l),ch)
        if(ifchar .gt. 0)then
          ifchar=ifchar+istart-1
        endif
      else
        ifchar=0
      endif
      return
      end

      subroutine capita1(string)
      character*(*) string
      character c,qc
      logical quote
      l=len(string)
      ioff=ichar('A')-ichar('a')
      quote=.false.
      qc=' '
      do 10 i=1,l
        c=string(i:i)
        if(quote)then
          if(c .eq. qc)then
            quote=.false.
          endif
        else        
          if(c .eq. '''' .or. c .eq. '"')then
            quote=.true.
            qc=c
          elseif(c .ge. 'a' .and. c .le. 'z')then
            string(i:i)=char(ichar(c)+ioff)
          endif
        endif
 10   continue
      return
      end
