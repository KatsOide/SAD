      subroutine tfstandardform(isp1,kx,irtc)
      use tfstk
      use tfform
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      type (sad_symbol), pointer :: symf,symp
      type (sad_symdef), pointer :: symfd,sympd
      integer*8 kf1,kp1
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(iaxform .eq. 0)then
        call tfforminit
      endif
      symf=>tfsydef(symform)
      symp=>tfsydef(sympw)
      call sym_symdef(symf,symfd)
      call sym_symdef(symp,sympd)
      kf1=symfd%value%k
      kp1=sympd%value%k
      symfd%value=dtfcopy1(dxnulls)
      sympd%value=dfromr(dble(nbmax-256))
c      call tfdebugprint(sympd%value,'standardform',1)
      call tfeeval(ktastk(isp),kx,.true.,irtc)
      call tflocald(symfd%value)
      symfd%value%k=kf1
      call tflocald(sympd%value)
      sympd%value%k=kp1
      return
      end

      character*(*) function tfgetform()
      use tfstk
      use tfform
      implicit none
      type (sad_descriptor) ks
      type (sad_string) , pointer :: str
      integer*4 nc,irtc
      if(iaxform .eq. 0)then
        call tfforminit
      endif
      call tfsyeval(iaxform,ks,irtc)
      if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
        tfgetform=' '
      else
        if(ktfstringq(ks,str))then
          if(str%nch .le. 0)then
            tfgetform=' '
          else
            nc=min(str%nch,len(tfgetform))
            tfgetform(1:nc)=str%str(1:nc)
            tfgetform(nc+1:)=' '
          endif
        else
          tfgetform=' '
        endif
      endif
      return
      end

      integer*4 function itfgetrecl()
      use tfstk
      use tfform
      implicit none
      type (sad_descriptor) ks
      integer*4 irtc
      if(iaxform .eq. 0)then
        call tfforminit
      endif
      call tfsyeval(iaxpagewidth,ks,irtc)
      if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
        itfgetrecl=131
      elseif(ktfnonrealq(ks,itfgetrecl))then
        itfgetrecl=131
      else
        if(itfgetrecl .le. 0)then
          itfgetrecl=131
        endif
      endif
      return
      end

      subroutine tfreestringbuf(strb)
      use tfstk
      use strbuf
      implicit none
      type (sad_strbuf) strb
      integer*8 kbuf
      kbuf=sad_loc(strb%nch)
      if(.not. tfonstackq(kbuf))then
        strb%indw=strb%maxnch/8+5
        call tfree(kbuf-2)
      endif
      return
      end

      subroutine flushstringbuf(strb,indent,cr,lfno,irtc)
      use strbuf
      implicit none
      type (sad_strbuf) strb
      integer*4 lcr,icr,m,irtc,lfno,maxint
      parameter (maxint=1073741823)
      logical*4 cr
      logical*4 indent
      icr=strb%column
      lcr=strb%indw
      if(icr .gt. 0 .and.
     $     (.not. indent .or. min(32,lcr) .lt. icr))then
        m=strb%nch
        strb%nch=icr
        call writestringbuf(strb,cr,lfno)
      else
        icr=0
        m=0
        lcr=strb%llevel
        call writestringbuf(strb,cr,lfno)
      endif
      strb%indw=maxint
      strb%maxllevel=lcr
      strb%column=0
      strb%remlines=strb%remlines-1
      if(strb%remlines .le. 0)then
        if(cr)then
          write(lfno,*,ERR=100)'...'
        else
          write(lfno,'(a,$)',ERR=100)'...'
        endif
 100    strb%nch=0
        irtc=-1
c        irtc=itfmessage(9,'General::longstr',' ')
        return
      endif
      if(indent)then
        strb%nch=min(lcr,32)
        strb%str(1:strb%nch)=' '
      endif
      if(icr .gt. 0)then
        strb%str(strb%nch+1:strb%nch+m-icr)=strb%str(icr+1:icr+m)
        strb%nch=strb%nch+m-icr
        if(strb%str(strb%nch:strb%nch) .eq. '\')then
          strb%nch=strb%nch-1
        endif
      endif
      irtc=0
      return
      end

      subroutine tmovb(from,to,n)
      implicit none 
      integer*4, intent(in):: n
      character*(*) , intent(in)::from
      character*(*) , intent(out)::to
      to(1:n)=from(1:n)
      return
      end

      subroutine writeb(string,l,nl,lfno)
      implicit none
      integer*4 l,lfno,i,j,maxlen
      integer*1 string(l)
      logical*4 nl
      parameter (maxlen=1024)
      if(l .gt. 0)then
c
c       I hope the below works for all compilers... but not for DEC
c      write(lfno,'(a,$)')string(1:l)
c
c        write(form,*)l
c        do i=12,1,-1
c          if(form(i:i) .eq. ' ')then
c            form(i:i)='a'
c            exit
c          endif
c        enddo
c        write(lfno,'($,'//form//')',ERR=100)string
c     parameter maxlen is a maximum array size for input of WRITE statement
c     maxlen = 1024 for Interl Fortran
         j=1
         do while(j .le. l)
            write(lfno,'($,1024a1)',ERR=100)(string(i),
     $           i=j,min(l,j+maxlen-1))
            j=j+maxlen
         enddo
      endif
      if(nl)then
         write(lfno,*,ERR=100)
       endif         
      return
 100  write(*,*)'Write error to file ',lfno
      return
      end
