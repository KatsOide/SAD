      subroutine getbuf
      use tfrbuf
      use tfcsi
      implicit none
      integer*4 lrecl0,lr
      if(lfni .le. 0)then
        ios=99999
        return
      endif
2     if(ipoint .ge. lrecl)then
        call tprmpt(lfni,lfno,lfn1)
        if(.not. rec)then
          lrecl0=linep+1
        else
          lrecl0=lrecl+1
        endif
        ipoint=lrecl0
        if(lrecl0 .gt. 1 .and.
     $       index(delim(1:ldel),buffer(lrecl0-1:lrecl0-1)) .le. 0)then
          write(*,*)'getbuf-unusual record ',rec,lrecl0,
     $         buffer(max(1,lrecl0-32):lrecl0)
          buffer(lrecl0:lrecl0)=char(10)
          lrecl0=lrecl0+1
        endif
 1      ios=-999
        if(lrecl0 .gt. nbmax-256)then
          ios=999999
          go to 10
        endif
        call tfreadbuf(irbreadrecord,lfni,int8(0),int8(0),
     $       lr,buffer(lrecl0:))
        if(lr .eq. -99)then
          go to 20
        elseif(lr .eq.  -999)then
          go to 10
        endif
        if(lr .le. 0)then
          lrecl=lrecl0-1
        else
          lrecl=len_trim(buffer(1:lrecl0+lr-1))
        endif
c        write(*,*)'getbuf ',lfni,lrecl0,lr,buffer(lrecl0:lrecl)
        if(lfn1 .gt. 0)then
          write(lfn1,'(1X,A)')buffer(lrecl0:lrecl)
        endif
        ios=0
        if(lrecl .ge. lrecl0)then
          call removetab(buffer(lrecl0:lrecl))
          call removecomment(buffer(lrecl0:lrecl),cmnt(1:lcmnt),'''"')
        endif
        if(lrecl .lt. lrecl0 .or. buffer(lrecl0:lrecl) .eq. ' ')then
          lrecl=lrecl0
        else
          lrecl=len_trim(buffer(1:lrecl))
          if(buffer(lrecl:lrecl) .eq. '\')then !'
            lrecl0=lrecl
            go to 1
          endif
          lrecl=lrecl+1
        endif
        buffer(lrecl:lrecl)=char(10)
      else
        call skipline
        if(ipoint .ge. lrecl)then
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
      if(lfni .eq. 5)then
        write(*,*)'???-getbuf-end of input stream'
        stop
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

      subroutine setbuf(string,nc)
      use tfcsi
      implicit none
      integer*4 nc
      character*(nc) string
      ipoint=lrecl+1
      linep=ipoint
      lrecl=ipoint+nc
      buffer(ipoint:lrecl-1)=string
      call removetab(buffer(ipoint:lrecl))
      call removecomment(buffer(ipoint:lrecl),cmnt(1:lcmnt),'''"')
      buffer(lrecl:lrecl)=char(10)
c      write(*,*)'setbuf ',ipoint,lrecl,string
      return
      end

      subroutine savebuf(string,nc)
      use tfcsi
      implicit none
      integer*4 nc,i
      character*(*) string
      nc=min(len(string),max(lrecl-ipoint,0))
      if(nc .gt. 0)then
        do i=ipoint,ipoint+nc-1
          if(buffer(i:i) .eq. char(10))then
            nc=i-ipoint+1
            exit
          endif
        enddo
        string(:nc)=buffer(ipoint:lrecl-1)
      endif
      return
      end

      subroutine removecomment(buf,comment,quote)
      character*(*) buf,quote
      character comment
      integer*4 i1,iq,ic,iq1,ifany1
      i1=1
      ic=index(buf,comment)
      do while(ic .gt. 0)
        iq=ifany1(buf,len(buf),quote,i1)
        if(iq .eq. 0 .or. iq .gt. ic)then
          buf(ic:)=' '
          return
        endif
        iq1=index(buf(iq+1:),buf(iq:iq))
        if(iq1 .le. 0)then
          return
        else
          i1=iq+iq1+1
          ic=index(buf(i1:),comment)
          if(ic .le. 0)then
            return
          endif
          ic=i1+ic-1
        endif
      enddo
      return
      end

      subroutine removetab(buf)
      implicit none
      character*(*) buf
      integer*4 i1,i
      i1=index(buf,char(9))
      do while(i1 .gt. 0)
        buf(i1:i1)=' '
        i=index(buf(i1+1:),char(9))
        if(i .gt. 0)then
          i1=i1+i
        else
          return
        endif
      enddo
      return
      end
