      integer*4 function itfgeto(kx)
      use tfstk
      use tfcsi
      use tfrbuf
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 irtc, m
      if(ipoint .gt. lrecl .or. ipoint .le. 0)then
        itfgeto=-1
        kx%k=ktfoper+mtfnull
        return
      endif
      m=0
c      write(*,*)'itfgeto ',lrecl,ipoint,buffer(max(ipoint-16,1):ipoint)
      call tfeval(buffer(1:lrecl),ipoint,m,kx,.true.,irtc)
c      if(lfni .lt. 100)then
c        call tfdebugprint(kx,'itfgeto',3)
c        write(*,'(a,10i12)')'with ',lfni,irtc,ipoint,lrecl,m,
c     $       kerror,ktfaddr(kerror)
c      endif
      ipoint=m
      if(irtc .eq. 0)then
        itfgeto=0
      elseif(irtc .gt. 0 .and. kerror .ne. 0)then
       if(rlist(ktfaddr(kerror)+1) .ge. 1000.d0)then
          itfgeto=-3
        else
          itfgeto=-2
        endif
        call tfreseterror
      else
        itfgeto=max(-3,min(-1,irtc+3))
      endif
      return
      end

      integer*4 function itfpeeko(kx,next)
      use tfstk
      use tfcsi, only:ipoint
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: next
      integer*4 ip0,itfgeto
      ip0=ipoint
      itfpeeko=itfgeto(kx)
      next=ipoint
      ipoint=ip0
      return
      end

      character*(*) function tfgetstrs(k,nc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_symbol), pointer :: sym
      type (sad_namtbl), pointer :: nam
      type (sad_string), pointer :: str
      integer*4 ,intent(inout):: nc
      if(ktfsymbolq(k,sym))then
        call sym_namtbl(sym,nam)
        str=>nam%str
      elseif(.not. ktfstringq(k,str))then
        tfgetstrs=' '
        nc=-1
        return
      endif
      nc=min(str%nch,len(tfgetstrs))
      if(nc .gt. 0)then
        tfgetstrs(1:nc)=str%str(1:nc)
        tfgetstrs(nc+1:)=' '
      else
        tfgetstrs=' '
      endif
      return
      end

      subroutine tfgetstrns(k,str,nc)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_string), pointer :: ks
      type (sad_symbol), pointer :: sym
      type (sad_namtbl), pointer :: nam
      integer*4 nc
      character*(*) str
      if(ktfstringq(k,ks))then
      elseif(ktfsymbolq(k,sym))then
        call loc_namtbl(sym%loc,nam)
        ks=>nam%str
      else
        nc=-1
        return
      endif
      nc=min(ks%nch,len(str))
      str(1:nc)=ks%str(1:nc)
      return
      end

      character*(*) function tfgetstr(k,nc)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_string), pointer :: str
      integer*4 nc
      if(ktfaddr(k) .eq. 0)then
        tfgetstr=' '
        nc=0
      else
        call descr_sad(k,str)
        nc=min(str%nch,len(tfgetstr))
        if(nc .gt. 0)then
          tfgetstr(1:nc)=str%str(1:nc)
          tfgetstr(nc+1:)=' '
        else
          tfgetstr=' '
        endif
      endif
      return
      end

      character*(*) function tfgetstrv(name)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_string), pointer :: str
      character*(*) name
      k=kxsymbolv(name,len_trim(name))
      if(ktfstringq(k,str))then
        tfgetstrv(1:str%nch)=str%str(1:str%nch)
        tfgetstrv(str%nch+1:)=' '
      else
        tfgetstrv=' '
      endif
      return
      end

      subroutine tfgettok(word)
      use tfstk
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 notspace,istop,nc,irt,i
      character*(*) word
 1    if(ipoint .ge. lrecl)then
        call getbuf
        if(ios .ne. 0)then
          word=' '
          return
        endif
      endif
      i=notspace(buffer(1:lrecl),ipoint)
      if(i .le. 0)then
        ipoint=lrecl
        go to 1
      endif
c      write(*,*)'itfgeto-tftok ',i,lrecl
      call tfetok(buffer(i:lrecl),istop,kx,-1,irt)
      nc=istop-1
      istop=istop+i-1
      ipoint=istop
c      write(*,*)'tfgettok ',ipoint,i,istop,nc,
c     $     buffer(istop-nc:istop-1)
      if(nc .le. 0)then
        go to 1
      endif
      if(ktfstringq(kx,str))then
        nc=min(str%nch,len(word))
        word(1:nc)=str%str(1:nc)
      else
        word(1:nc)=buffer(istop-nc:istop-1)
      endif
      word(nc+1:)=' '
      call capita(word(1:nc))
      return
      end
