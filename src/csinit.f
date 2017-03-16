      module tfcsi
        use tfcbk, only:maxlbuf
        integer*4, parameter :: nbmax=maxlbuf
        integer*4 ipoint,iconv,ios,lfn1,lfni,lfno,lrecl,ldel,
     $       lcmnt,linep,lastln,ibegt,lastt
        logical*4 rec
        character*16 delim,cmnt
        character*(nbmax) buffer

        contains
        subroutine cssetp(ip)
        implicit none
        integer*4 ip
        ipoint=ip
        return
        end subroutine

        subroutine cssets(ip)
        implicit none
        integer*4 ip
        ios=ip
        return
        end subroutine

        subroutine cssetl(ip)
        implicit none
        integer*4 ip
        lrecl=ip
        return
        end subroutine

        subroutine cssetlfni(ip)
        implicit none
        integer*4 ip
        lfni=ip
        return
        end subroutine

        subroutine cssetlfno(ip)
        implicit none
        integer*4 ip
        lfno=ip
        return
        end subroutine

        subroutine cssetlfn1(ip)
        implicit none
        integer*4 ip
        lfn1=ip
        return
        end subroutine

        subroutine cssetrec(f)
        implicit none
        logical*4 f
        rec=f
        return
        end subroutine

        subroutine cssetlinep(ip)
        implicit none
        integer*4 ip
        linep=ip
        return
        end subroutine

        integer*4 function icsmrk()
        implicit none
        icsmrk=ipoint
        return
        end function

        integer*4 function icsstat()
        implicit none
        icsstat=ios
        return
        end function

        integer*4 function icslrecl()
        implicit none
        icslrecl=lrecl
        return
        end function

        integer*4 function icslfni()
        implicit none
        icslfni=lfni
        return
        end function

        integer*4 function icslfno()
        implicit none
        icslfno=lfno
        return
        end function

        integer*4 function icslfn1()
        implicit none
        icslfn1=lfn1
        return
        end function

        integer*4 function icslinep()
        implicit none
        icslinep=linep
        return
        end function

        logical*4 function csrec()
        implicit none
        csrec=rec
        return
        end function

      end module

      subroutine csinit(lfn0,iconv1,cmnt1,rec0)
      use tfcsi
      implicit none
      integer*4 lfn0,iconv1
      character*(*) cmnt1
      logical*4 rec0
      lfn1=lfn0
      lrecl=1
      buffer(1:1)=char(10)
      ipoint=1
      iconv=iconv1
      call tfsetconvcase(iconv1 .eq. 1)
      delim=';'//char(10)//' ,'
      ldel=len_trim(delim)
      cmnt=cmnt1
      lcmnt=len_trim(cmnt)
      rec=rec0
      linep=1
      lastln=0
      ibegt=0
      lastt=0
      call csrst(lfn0)
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
