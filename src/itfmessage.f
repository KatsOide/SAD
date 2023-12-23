      module tfmessage
      use tfstk
      logical*4 newsym,newset
      data newsym,newset/.true.,.true./

      contains
      type (sad_descriptor) function kxmessagename(mess,comp)
      implicit none
      type (sad_descriptor) k1
      type (sad_symdef), pointer :: symd
      integer*4 i,l,isp1
      character*(*) ,intent(in):: mess
      logical*4 ,optional,intent(in)::comp
      i=index(mess,'::')
      if(i <= 0)then
        write(*,*)'itfmessagename implementation error: ',mess
        call abort
      endif
      isp1=isp+1
      ktastk(isp1)=ktfoper+mtfmessagename
      k1=kxsymbolz(mess(1:i-1),i-1,symd)
      if(ktfoperq(symd%value))then
        dtastk(isp1+1)=symd%value
      else
        dtastk(isp1+1)=k1
      endif
      l=len(mess)
      isp=isp1+2
      dtastk(isp)=kxsymbolz(mess(i+2:l),l-i-1)
      if(.not. present(comp) .or. comp)then
        kxmessagename=kxcomposev(isp1)
        isp=isp1-1
      else
        kxmessagename%k=0
      endif
      return
      end function 

      logical*4 function tfnewsym(ini) result(r)
      use dset
      implicit none
      logical*4 ,intent(in):: ini
      type (sad_descriptor) km,kx
      type (sad_symdef) ,pointer ::symd
      integer*4 isp0,irtc
      logical*4 ev
      logical*4 ,save::init=.false.
      integer*8,save :: ka=0
      if(ini)then
        isp0=isp
        km=kxmessagename('General::newsym',.false.)
        call loc_symdef(klist(ifunbase+mtfmessagename),symd)
        call tfdeval(isp0+1,dfromk(ksad_loc(symd%sym%loc)),kx,1,.true.,ev,irtc)
        ka=ktfaddr(kx)
        isp=isp0
        init=.true.
      endif
      r=newsym .and. init .and. ktfstringq(dlist(ka))
      return
      end function

      logical*4 function tfnewset(ini) result(r)
      use dset
      implicit none
      logical*4 ,intent(in):: ini
      type (sad_descriptor) km,kx
      type (sad_symdef) ,pointer ::symd
      integer*4 isp0,irtc
      logical*4 ev
      logical*4 ,save::init=.false.
      integer*8,save :: ka=0
      if(ini)then
        isp0=isp
        km=kxmessagename('General::newset',.false.)
        call loc_symdef(klist(ifunbase+mtfmessagename),symd)
        call tfdeval(isp0+1,dfromk(ksad_loc(symd%sym%loc)),kx,1,.true.,ev,irtc)
        ka=ktfaddr(kx)
        isp=isp0
        init=.true.
      endif
      r=newset .and. init .and. ktfstringq(dlist(ka))
      return
      end function

      end module

      integer*4 function itfmessage(level,mess0,arg0)
      use tfstk
      use tfmessage
      use tfcsi
      use eeval
      implicit none
      type (sad_descriptor) dm,ks,mn
      type (sad_dlist), pointer :: klx,klm,klhm,klhms
      integer*4 ,intent(in):: level
      integer*4 irtc,ltr0,l,la,isp0
      character*(*) ,intent(in):: mess0,arg0
      character*256 mess,arg
      logical*4 ,save::iter=.false.
      data mn%k /0/
      if(iter .or. level .lt. ierrorth .or.
     $     rlist(iaximmediate) .lt. 0.d0)then
        itfmessage=-1
        return
      endif
      if(mn%k == 0)then
        mn=kxsymbolz('MessageString',13)
      endif
      iter=.true.
      ltr0=ltrace
      ltrace=0
      l=len_trim(mess0)
      mess=mess0(:l)
      la=len_trim(arg0)
      arg=arg0(:la)
      dm=dtfcopy1(kxmessagename(mess(:l)))
      call descr_sad(dm,klm)
      ks=tfleval(klm,.true.,irtc)
      if(irtc /= 0)then
        if(irtc > 0)then
          call tfaddmessage(' ',0,icslfnm())
        endif
        itfmessage=-1
        call tflocal1(dm%k)
      elseif(ktfstringq(ks))then
        call tflocal(kerror)
        kerror=ktflist+ktadaloc(0,4,klx)
        klx%rbody(1)=dble(level)
        klx%dbody(2)=kxadaloc(0,1,klhm)
        klx%dbody(3)=kxsalocb(0,mess,l)
        klx%dbody(4)=kxadaloc(0,1,klhms)
        klhm%head%k=ktfoper+mtfhold
        klhm%dbody(1)=dm
        klhms%head%k=ktfoper+mtfhold
        isp0=isp
        isp=isp+1
        dtastk(isp)=mn
        isp=isp+1
        dtastk(isp)=ks
        call tfstringliststk(arg(:la))
        klhms%dbody(1)=dtfcopy1(kxcomposev(isp0+1))
        isp=isp0
        rlist(ierrorgen)=rlist(ierrorgen)+1.d0
        ierrorprint=sad_loc(klx%head%k)
        itfmessage=1
      else
        call tflocal1(dm%k)
        itfmessage=-1
      endif
      iter=.false.
      ltrace=ltr0
      return
      end

c     tfstringliststk destroy `string' to decode backslash escape
      subroutine tfstringliststk(string)
      use tfstk
      implicit none
      character*(*) ,intent(inout):: string
      integer*4 l,is,is1,in
      l=len_trim(string)
c      write(*,*)'tfstringliststk ',l,' ',string(1:l)
      is=1
      do while(is <= l)
c     Search '"' from string(is:l) with backslash escape
        is1=is
        do while(is <= l .and. string(is1:is1) /= '"')
           if(string(is1:is1) == '\\') then
              string(is1:l-1)=string(is1+1:l)
              l=l-1
           endif
           is1=is1+1
        enddo
        if(is1 > l) then
           is1=0
        endif
        if(is1 <= 0)then
          isp=isp+1
          dtastk(isp)=kxsalocb(-1,string(is:),l-is+1)
          return
        endif
        is=is1+1
c     Search '"' from string(is+1:l) with backslash escape
        in=is
        do while(in <= l .and. string(in:in) /= '"')
           if(string(in:in) == '\\') then
              string(in:l-1)=string(in+1:l)
              l=l-1
           endif
           in=in+1
        enddo
        in=in-is
        if(in > l) then
           in=0
        endif
        if(in <= 0)then
          isp=isp+1
          dtastk(isp)=kxsalocb(-1,string(is:),l-is+1)
          return
        endif
        isp=isp+1
        dtastk(isp)=kxsalocb(-1,string(is:),in)
        is=is+in+1
      enddo
      return
      end

      integer*4 function itfmessageexp(level,mess,k)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(in):: level
      integer*4 nc,itfmessage
      character*(*) ,intent(in):: mess
      character*256 buf
      call tfconvstrs(buf(2:),k,nc,.true.,' ')
      nc=max(0,min(254,nc))
      buf(1:1)='"'
      buf(nc+2:nc+2)='"'
      itfmessageexp=itfmessage(level,mess,buf(:nc+2))
      return
      end

      integer*4 function itfmessagestr(level,mess,str)
      use tfstk
      use strbuf
      implicit none
      integer*4 ,parameter :: lbuf=128
      type (sad_strbuf) , pointer :: strb
      character*(*) ,intent(in):: str
      integer*4 ,intent(in):: level
      character*(*) ,intent(in):: mess
      integer*4 irtc,isp0,itfmessage
      isp0=isp
      call getstringbuf(strb,lbuf*2,.true.)
      call tfquotestring(strb,str,min(lbuf-1,len(str)),0,irtc)
      itfmessagestr=itfmessage(level,mess,strb%str(1:strb%nch))
      isp=isp0
      return
      end

      subroutine tferrorhandle(kx,irtc)
      use tfstk
      use tfcsi,only:icslfnm
      implicit none
      type (sad_descriptor) ,intent(in):: kx
      integer*4 ,intent(in):: irtc
      integer*4 lf
      lf=icslfnm()
      if(lf /= 0 .and. ierrorprint /= 0)then
        ierrorf=kx%k
        if(irtc > 0  .and. tflistq(kerror))then
          call tfaddmessage(' ',0,lf)
        endif
        ierrorprint=0
      endif
      return
      end
      
      subroutine tferrorhandles(kx,str,irtc)
      use tfstk
      use tfcsi
      implicit none
      type (sad_descriptor) ,intent(in):: kx
      character*(*) ,intent(in):: str
      integer*4 ,intent(in):: irtc
      if(ierrorprint /= 0)then
        ierrorf=kx%k
        if(irtc > 0)then
          if(tflistq(kerror))then
            call tfaddmessage(str(1:len_trim(str)),
     $           len_trim(str),icslfnm())
          endif
        endif
        ierrorprint=0
      endif
      return
      end
      
      subroutine tfaddmessage(str,ip,lfn)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kxaddmess,tfefunrefd
      type (sad_dlist), pointer :: kle
      integer*4 ,intent(in):: ip,lfn
      integer*4 irtc,isp1,ltr0
      character*(*) ,intent(in):: str
      logical*4 iter
      data kxaddmess%k,iter /0,.false./
      if(iter)then
        call tfdebugprint(kerror,'???Error in error handling.',10)
        return
      endif
      if(.not. tflistq(kerror,kle))then
        return
      endif
      iter=.true.
      ltr0=ltrace
      ltrace=0
      if(kxaddmess%k == 0)then
        kxaddmess=kxsymbolz('Add$Message',11)
      endif
      isp1=isp+1
      if(lfn /= 0)then
        isp=isp1+1
        dtastk(isp1)=kxaddmess
        ktastk(isp)=kerror
        kx=tfefunrefd(isp1,irtc)
        isp=isp1-1
        if(irtc == 0)then
          if(ktflistq(kx) .and. lfn /= 0)then
            call tfemes(kx,str,ip,lfn)
          else
            call tflocal(kerror)
            kerror=0
            call tclrfpe
          endif
        endif
      endif
      ierrorf=0
      ierrorprint=0
      iter=.false.
      ltrace=ltr0
      return
      end

      subroutine tfemes(k,string,ip,lfn)
      use tfstk
      implicit none
      integer*4 ,parameter :: lstr=512,lenbuf=2**23
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl
      character*(*) ,intent(in):: string
      integer*4 ,intent(in):: ip,lfn
      integer*4 ls,ifchar,is,i,lb,nc,itfgetrecl,nc1
      character*10 autofg
      character*(lstr) tfgetstr,buff
      if(ktfaddr(k) == 0 .or. ierrorprint == 0)then
        ierrorf=0
        return
      elseif(tflistq(k,kl))then
        buff(1:3)='???'
        buff(4:)=tfgetstr(kl%dbody(3),nc)
        buff(nc+4:nc+6)=': '
        buff(nc+7:)=tfgetstr(kl%dbody(4),nc1)
        call capita(buff(nc+7:nc+7))
      else
        buff='???Unexpected error: '//autofg(dble(k%k),'S10')
      endif
      ierrorprint=0
      if(lfn /= 0)then
        lb=min(len_trim(buff),lstr-1)
        if(ierrorf /= 0)then
          buff(lb+1:lb+4)=' in '
          call tfconvstrs0(buff(lb+5:),lenbuf,ierrorf,nc,.true.,' ')
c     str1=tfconvstr(ierrorf,nc,'*')
          lb=min(lb+nc+4,lstr)
        else
          lb=lb+1
          buff(lb:lb)=':'
        endif
        call tftruncprint(buff(1:lb),itfgetrecl()-1,' ,})]'//char(10),.false.,lfn)
        if(lfn /= 6)then
          call tftruncprint(buff(1:lb),itfgetrecl()-1,' ,})]'//char(10),.false.,6)
        endif
        if(ip > 0)then
          do i=ip-2,1,-1
            if(string(i:i) == char(10))then
              is=i+1
              go to 1
            endif
          enddo
          is=max(1,ip-32)
 1        ls=ifchar(string(1:ip),char(10),is)-1
          if(ls <= 0)then
            ls=len_trim(string(1:ip))-1
          endif
          if(ls .ge. is)then
            write(lfn,'(a)')string(is:ls)
            if(lfn /= 6)then
              write(6,'(a)')string(is:ls)
            endif
            if(ierrorf == 0)then
              ls=ip-is
              if(ls <= 0)then
                ls=1
              endif
              buff(1:ls)=' '
              buff(ls:ls)='^'
              write(lfn,'(a)')buff(1:ls)
              if(lfn /= 6)then
                write(6,'(a)')buff(1:ls)
              endif
            endif
          endif
        endif
      endif
      ierrorf=0
      call tflocal(kerror)
      kerror=0
      call tclrfpe
      return
      end

      subroutine tfreseterror
      use tfstk
      implicit none
      if(kerror /= 0)then
        call tflocal(kerror)
        kerror=0
      endif
      if(ierrorprint /= 0)then
        ierrorprint=0
      endif
      return
      end

      subroutine tffserrorhandle(l,irtc)
      use tfstk
      implicit none
      integer*4 l,irtc
      call tferrorhandle(sad_descr(dble(l)),irtc)
      call tfreseterror
      return
      end
