      module readopt
        type ropt
        sequence
        character*64 delim
        integer*4 ndel
        logical*4 new,null,del,opt
        end type
      end module

      subroutine tfwrite(isp1,kx,irtc)
      use tfstk
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_strbuf), pointer :: strb
      integer*4 isp1,irtc,itfgetrecl,narg,lfn,isp11,
     $     j,i,nc,lpw,itfmessage,isp2,itfgetlfn
      logical*4 exist
      narg=isp-isp1
      if(narg .le. 1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      lfn=itfgetlfn(isp1,.false.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      lpw=itfgetrecl()
      isp11=isp
      call getstringbuf(strb,-lpw,.true.)
      exist=.false.
      isp2=isp
      do j=isp1+2,isp11
        call tfevalstk(dtastk(j),.true.,irtc)
        if(irtc .ne. 0)then
          go to 10
        endif
        do i=isp2+1,isp
          call tfconvstrb(strb,dtastk(i),nc,
     $         .false.,.false.,lfn,'*',irtc)
          if(irtc .ne. 0)then
            go to 10
          endif
          exist=exist .or. nc .gt. 0
        enddo
        isp=isp2
      enddo
      call writestringbufn(strb,.true.,lfn)
      if(.not. exist)then
        write(lfn,*,ERR=9010)
      endif
 9010 kx%k=ktfoper+mtfnull
      irtc=0
 10   call tfreestringbuf(strb)
      isp=isp11
      return
      end

      integer*4 function itfgetlfn(isp1,read,irtc)
      use tfstk
      use tfcsi
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      integer*4 isp1,irtc,itfmessage
      logical*4 read
      itfgetlfn=0
      if(isp .le. isp1)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      if(ktfrealq(ktastk(isp1+1)))then
        itfgetlfn=int(rtastk(isp1+1))
        go to 100
      elseif(tflistqd(dtastk(isp1+1),kl))then
        if(read)then
          k=kl%dbody(1)
        else
          k=kl%dbody(2)
        endif
        if(ktfrealqdi(k,itfgetlfn))then
          go to 100
        endif
      endif
      irtc=itfmessage(9,'General::wrongtype',
     $     '"Real number or List of two Reals"')
      return
 100  irtc=0
      if(itfgetlfn .eq. -1)then
        if(read)then
          itfgetlfn=icslfni()
        else
          itfgetlfn=icslfno()
        endif
      elseif(itfgetlfn .lt. 0)then
        irtc=itfmessage(9,'General::wrongnum','"positive, 0 or -1"')
      endif
      return
      end

      subroutine tfwritestring(isp1,kx,irtc)
      use tfstk
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_strbuf), pointer :: strb
      integer*4 mmax
      parameter (mmax=1000000)
      integer*4 isp1,irtc,narg,lfn,isp11,j,i,itfgetlfn,
     $             nc,itfmessage,isp2
      narg=isp-isp1
      if(narg .le. 1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      lfn=itfgetlfn(isp1,.false.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      isp11=isp
      do j=isp1+2,isp11
        call tfevalstk(dtastk(j),.true.,irtc)
        if(irtc .ne. 0)then
          isp=isp11
          return
        endif
        isp2=isp
        call getstringbuf(strb,0,.true.)
        do i=isp11+1,isp2
          call tfconvstrb(strb,dtastk(i),
     $         nc,.false.,.false.,-1,'*',irtc)
          if(irtc .ne. 0)then
            go to 10
          endif
          call writestringbuf(strb,.false.,lfn)
        enddo
        call tfreestringbuf(strb)
        isp=isp11
      enddo
      kx%k=ktfoper+mtfnull
      return
 10   call tfreestringbuf(strb)
      isp=isp11
      return
      end

      subroutine tfprintf(isp1,kx,irtc)
      use tfstk
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,isp0
      isp0=isp
      isp=isp+1
      rtastk(isp)=icslfno()
      ktastk(isp+1:isp+isp0-isp1)=ktastk(isp1+1:isp0)
      isp=isp+isp0-isp1
c      do i=1,isp0-isp1
c        isp=isp+1
c        ktastk(isp)=ktastk(isp1+i)
c      enddo
      call tfwrite(isp0,kx,irtc)
      isp=isp0
      return
      end

      subroutine tfdebugprint(k,pr1,nline)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) k
      integer*4 irtc,nline,l,itfdownlevel
      character*(*) pr1
      prolog=pr1
      ncprolog=min(len_trim(pr1)+1,len(prolog))
      prolog(ncprolog:ncprolog)=' '
      levele=levele+1
      call tfprint1(k,6,79,nline,.true.,.true.,irtc)
      if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
      l=itfdownlevel()
      return
      end

      subroutine tfshort(isp1,kx,irtc)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfgetrecl,nret,itfmessage
      if(isp .eq. isp1+1)then
        call tfprint1(dtastk(isp),
     $       icslfno(),-itfgetrecl(),1,.true.,.true.,irtc)
      elseif(isp .eq. isp1+2)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype','"real number"')
          return
        endif
        nret=int(rtastk(isp))
        call tfprint1(dtastk(isp1+1),
     $       icslfno(),-itfgetrecl(),abs(nret),nret .ge. 0,
     $       .true.,irtc)
      else
        irtc=itfmessage(9,'General::narg','"1 or 2"')
      endif
      kx%k=ktfoper+mtfnull
      return
      end

      subroutine tfprint1(k,lfno,lrec,nret,cr,str,irtc)
      use tfstk
      use tfrbuf
      use strbuf
      implicit none
      type (sad_descriptor) k,kxlongstr,ks,kr
      type (sad_strbuf), pointer :: strb
      integer*4 lfno,irtc,nc,lrec,nret,irtc1,isp0
      logical*4 cr,str,tfsameqd
      save kxlongstr
      data kxlongstr%k /0/
      isp0=isp
      call getstringbuf(strb,lrec,.true.)
      if(nret .gt. 0)then
        strb%remlines=nret
      endif
      if(ncprolog .gt. 0)then
        call putstringbufp(strb,prolog(1:ncprolog),lfno,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      endif
      call tfconvstrb(strb,k,nc,str,.false.,lfno,' ',irtc)
      if(irtc .eq. 0)then
        call writestringbufn(strb,cr,lfno)
      elseif(irtc .eq. -1)then
        irtc=0
      elseif(irtc .gt. 0)then
        kr=dlist(ktfaddr(kerror)+2)
        if(kxlongstr%k .eq. 0)then
          call tfevals('Hold[General::longstr]',ks,irtc1)
          if(irtc1 .ne. 0)then
            irtc=irtc1
            go to 9000
          endif
          kxlongstr=dtfcopy1(ks)
        endif
        if(tfsameqd(kr,kxlongstr))then
          call tfreseterror
          irtc=0
        endif
      endif
 9000 call tfreestringbuf(strb)
      ncprolog=0
      isp=isp0
      return
      end

      subroutine tfdefinition(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ki,k
      type (sad_symbol), pointer ::sym
      type (sad_symdef), pointer ::symd
      type (sad_list), pointer :: klh,kli,klx
      integer*8 ii,ka0,kadi,kad,ka
      integer*4 isp1,irtc,i,isp00,isp0,kk,iop,ispv,icomp
      icomp=0
      call tfgetoption('Compiled',ktastk(isp),kx,irtc)
      if(irtc .eq. -1)then
        ispv=isp
      elseif(irtc .ne. 0)then
        return
      elseif(ktfrealqd(kx))then
        if(kx%k .ne. 0)then
          icomp=1
        endif
        ispv=isp-1
      else
        ispv=isp
      endif
      isp00=isp
      do i=isp1+1,ispv
        isp0=isp
        k=dtastk(i)
        if(ktfoperqd(k,ka))then
          ka=klist(ifunbase+ka)
          k%k=ktfsymbol+ka
        endif
        if(ktfsymbolqd(k,sym))then
          call tfsydef(sym,sym)
          call sym_symdef(sym,symd)
          ka0=ksad_loc(symd%upval)
          if(ktfrefqd(symd%value,ka))then
            k=dlist(ka)
          else
            k=symd%value
          endif
        else
          ka0=0
        endif
        isp=isp+1
        if(ka0 .ne. 0)then
          dtastk(isp)=kxadaloc(-1,2,klh)
          klh%head=ktfoper+mtfsetdelayed
          klh%body(1)=ktfsymbol+ktfcopy1(ka0+6)
          klh%dbody(2)=dtfcopy(k)
          do kk=1,0,-1
            if(kk .eq. 1)then
              iop=mtfsetdelayed
            else
              iop=mtfupsetdelayed
            endif
            kad=klist(ka0+kk)
            do while(kad .ne. 0)
              if(ilist(1,kad+2) .eq. maxgeneration)then
                do ii=kad+3,kad+ilist(2,kad+2)+3
                  kadi=klist(ii)
                  do while(kadi .ne. 0)
                    isp=isp+1
                    dtastk(isp)=kxadaloc(-1,2,kli)
                    kli%head=ktfoper+iop
                    kli%body(1)=ktfcopy1(klist(kadi+3+icomp))
                    kli%body(2)=ktfcopy(klist(kadi+5+icomp))
                    if(kli%body(2) .eq. ktfref)then
                      kli%body(2)=ktfcopy(klist(kadi+6))
                    endif
                    kadi=klist(kadi)
                  enddo
                enddo
              else
                isp=isp+1
                dtastk(isp)=kxadaloc(-1,2,kli)
                kli%head=ktfoper+iop
                kli%body(1)=ktfcopy1(klist(kad+3+icomp))
                kli%body(2)=ktfcopy(klist(kad+5+icomp))
                if(kli%body(2) .eq. ktfref)then
                  kli%body(2)=ktfcopy(klist(kad+6))
                endif
              endif
              kad=klist(kad)
            enddo
          enddo
        else
          dtastk(isp)=k
        endif
        ki=kxmakelist(isp0)
        isp=isp0+1
        dtastk(isp)=ki
      enddo
      kx=kxmakelist(isp00,klx)
      klx%head=ktfoper+mtfhold
      isp=isp00
      irtc=0
      return
      end

      subroutine tfgetf(fname)
      use tfstk
      use tfcode
      implicit none
      character*(*) fname
      type (sad_descriptor) kx,k
      integer*4 irtc
      k=kxsalocb(-1,fname,len_trim(fname))
      call tfget(k,kx,irtc)
      return
      end

      subroutine tfget(k,kx,irtc)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) k,kx,kf,kfn
      integer*4 irtc,itfgeto,lfni0,lfn10,ip0,lr0,
     $     lfn,isp0,itf,nc
      isp0=isp
      isp=isp+1
      dtastk(isp)=k
      call tfopenread(isp0,kfn,irtc)
      isp=isp0
      if(irtc .ne. 0)then
        return
      endif
      lfn=int(rfromd(kfn))
      lfni0=icslfni()
      lfn10=icslfn1()
      ip0=icsmrk()
      lr0=icslrecl()
      rec=csrec()
      linep=icslinep()
      call cssetrec(.false.)
      call cssetlinep(lr0)
      call cssetp(lr0)
      call cssetlfni(lfn)
      call cssetlfn1(0)
      levele=levele+1
      kx%k=ktfoper+mtfnull
      itf=0
      do while(itf .ge. 0)
        itf=itfgeto(kf)
        if(itf .ge. 0)then
          kx=kf
        endif
        if(itf .eq. -1)then
          itf=0
          call skipln
        endif
        if(icsstat() .ne. 0)then
          call cssets(0)
          itf=min(-1,itf)
        endif
      enddo
      call tfconnect(kx,0)
      close(lfn)
      call tfreadbuf(irbclose,lfn,int8(0),int8(0),nc,' ')
      call cssetlfni(lfni0)
      call cssetlfn1(lfn10)
      call cssetl(lr0)
      call cssetp(ip0)
      call cssetlinep(linep)
      call cssetrec(rec)
      if(itf .eq. -3)then
        irtc=irtcabort
      else
        irtc=0
      endif
      return
      end

      subroutine tfread1(isp1,lfn,kx,irtc)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) kx,kf
      integer*4 isp1,irtc,itfgeto,nc,
     $     lfni0,lfn10,ip0,lr0,itf,
     $     lfn,linep0,itfmessage
      logical*4 openf
      if(isp .gt. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      lfni0=icslfni()
      openf=lfn .gt. 0 .and. lfn .ne. lfni0
      if(openf)then
        lfn10=icslfn1()
        ip0=icsmrk()
        lr0=icslrecl()
        linep0=icslinep()
        call cssetlinep(lr0)
        nc=0
        call tfreadbuf(irbreadbuf,lfn,int8(0),int8(0),nc,buffer)
        if(nc .gt. 0)then
          call setbuf(buffer(1:nc),nc)
        else
          call cssetp(lr0)
        endif
        call cssetlfni(lfn)
        call cssetlfn1(0)
      endif
      levele=levele+1
 1    itf=itfgeto(kf)
      if(itf .ge. 0)then
        kx=kf
      else
        kx%k=ktfoper+mtfnull
      endif
      if(itf .eq. -1)then
        call tprmpt(-1,-1,0)
        call getbuf
        call tprmpt(0,-1,0)
        if(icsstat() .eq. 0)then
          go to 1
        endif
      endif
      if(icsstat() .ne. 0)then
        call cssets(0)
        kx%k=kxeof
        if(openf)then
          call tfreadbuf(irbreset,lfn,int8(0),int8(0),nc,' ')
        endif
      elseif(openf)then
        call savebuf(buffer,nc)
        if(nc .gt. 0)then
          call tfreadbuf(irbsetbuf,lfn,int8(0),int8(0),nc,buffer)
        else
          call tfreadbuf(irbreset,lfn,int8(0),int8(0),nc,' ')
        endif
      endif
      call tfconnect(kx,0)
      if(openf)then
        call cssetlfni(lfni0)
        call cssetl(lr0)
        call cssetlinep(linep0)
        call cssetp(ip0)
        call cssetlfn1(lfn10)
      endif
      irtc=0
      return
      end

      recursive subroutine tfskip(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,i,narg,isp0,nrpt
      narg=isp-isp1
      if(narg .ge. 3 .and.
     $     ktfrealqdi(dtastk(isp1+3),nrpt))then
        isp0=isp
        do i=1,nrpt
          ktastk(isp0+1)=ktastk(isp1+1)
          ktastk(isp0+2)=ktastk(isp1+2)
          ktastk(isp0+3:isp0+narg-1)=ktastk(isp1+4:isp1+narg)
c          do j=4,narg
c            ktastk(isp0+j-1)=ktastk(isp1+j)
c          enddo
          isp=isp0+narg-1
          call tfskip(isp0,kx,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            return
          endif
          if(kx%k .eq. kxeof)then
            return
          endif
        enddo
      else
        call tfread(isp1,kx,irtc)
      endif
      return
      end

      subroutine tfread(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,lfn,itfgetlfn,irtc
      lfn=itfgetlfn(isp1,.true.,irtc)
      if(irtc .eq. 0)then
        call tfreadf(isp1,lfn,kx,irtc)
      endif
      return
      end

      subroutine tfreadf(isp1,lfn,kx,irtc)
      use tfstk
      use readopt
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,lfn
      type(ropt) opts
      call tfreadoptions(isp1,opts,.true.,irtc)
      if(irtc .eq. 0)then
        isp=min(isp,isp1+2)
        call tfreadfs(isp1,lfn,opts,kx,irtc)
      endif
      return
      end

      subroutine tfreadfs(isp1,lfn,opts,kx,irtc)
      use tfstk
      use readopt
      implicit none
      type (sad_descriptor) kx,k1,k2
      type (sad_list), pointer :: list
      integer*4 isp1,irtc,isp0,n1,narg,itfmessage,lfn,i1
      logical*4 tfsamesymbolqd
      type(ropt) opts
      type (sad_descriptor)
     $     itfexprs,itfstrs,itfwords,itfreals,itfchars,
     $     itfbyte,itfbinint,itfbinsint,itfbinreal,itfbindble
      data itfexprs%k,itfstrs%k,itfwords%k,itfreals%k,itfchars%k,
     $     itfbyte%k,itfbinint%k,itfbinsint%k,itfbinreal%k,itfbindble%k
     $     /0,0,0,0,0, 0,0,0,0,0/
      narg=isp-isp1
      if(narg .eq. 0)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
      elseif(narg .eq. 1)then
        call tfread1(isp1,lfn,kx,irtc)
      else
        if(itfexprs%k .eq. 0)then
          itfexprs=kxsymbolf('Expression',10,.true.)
          itfstrs=kxsymbolf('String',6,.true.)
          itfwords=kxsymbolf('Word',4,.true.)
          itfreals=kxsymbolf('Real',4,.true.)
          itfchars=kxsymbolf('Character',9,.true.)
          itfbyte=kxsymbolf('Byte',4,.true.)
          itfbinsint=kxsymbolf('BinaryShortInteger',18,.true.)
          itfbinint=kxsymbolf('BinaryInteger',13,.true.)
          itfbinreal=kxsymbolf('BinaryReal',10,.true.)
          itfbindble=kxsymbolf('BinaryDouble',12,.true.)
        endif
        k2=dtastk(isp1+2)
        if(ktfsymbolqd(k2))then
          if(tfsamesymbolqd(k2,itfexprs))then
            isp0=isp
            isp=isp+1
            ktastk(isp)=ktastk(isp1+1)
            call tfread1(isp0,lfn,kx,irtc)
            isp=isp0
          elseif(tfsamesymbolqd(k2,itfstrs))then
            isp0=isp
            isp=isp+1
            ktastk(isp)=ktastk(isp1+1)
            n1=opts%ndel
            opts%ndel=0
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
            opts%ndel=n1
            isp=isp0
          elseif(tfsamesymbolqd(k2,itfwords))then
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
          elseif(tfsamesymbolqd(k2,itfchars))then
            call tfreadstringf(lfn,kx,.true.,opts,irtc)
          elseif(tfsamesymbolqd(k2,itfreals))then
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(ktfstringqd(kx))then
              isp0=isp
              isp=isp+1
              dtastk(isp)=kx
              call tftoexpression(isp0,kx,irtc)
              isp=isp0
            endif
          elseif(tfsamesymbolqd(k2,itfbyte))then
            call tfreadbyte(isp1,lfn,kx,1,irtc)
          elseif(tfsamesymbolqd(k2,itfbinsint))then
            call tfreadbyte(isp1,lfn,kx,2,irtc)
          elseif(tfsamesymbolqd(k2,itfbinint))then
            call tfreadbyte(isp1,lfn,kx,3,irtc)
          elseif(tfsamesymbolqd(k2,itfbinreal))then
            call tfreadbyte(isp1,lfn,kx,4,irtc)
          elseif(tfsamesymbolqd(k2,itfbindble))then
            call tfreadbyte(isp1,lfn,kx,5,irtc)
          else
            go to 9000
          endif
        elseif(ktflistqd(k2,list))then
          if(list%head .eq. ktfoper+mtflist)then
            call tfreadfm(lfn,list,list%nl,opts,.false.,kx,irtc)
          elseif(list%head .eq. ktfoper+mtftimes .and.
     $           list%nl .eq. 2)then
            k1=list%dbody(1)
            if(ktfnonrealqdi(k1,i1))then
              go to 9000
            endif
            call tfreadfm(lfn,list,i1,opts,.true.,kx,irtc)
          else
            go to 9000
          endif
        else
          go to 9000
        endif            
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $     '"(n*) BinaryDouble, BinaryInteger, BinaryReal, '//
     $     ' BinaryShortInteger, Byte, Character, Expression, '//
     $     ' BinaryFloat, Real, String, or Word is","for 2nd arg"')
      return
      end

      subroutine tfreadfm(lfn,list,m,opts,mult,kx,irtc)
      use tfstk
      use readopt
      implicit none
      type (sad_descriptor) kx
      type (sad_list) list
      type (sad_list), pointer ::kl,klx
      integer*8 kxi
      integer*4 lfn,irtc,isp0,isp2,kk,m
      logical*4 mult
      type(ropt) opts
      isp2=isp
      do kk=1,m
        isp0=isp
        isp=isp0+2
        if(mult)then
          ktastk(isp)=list%body(2)
        else
          ktastk(isp)=list%body(kk)
        endif
        call tfreadfs(isp0,lfn,opts,kxi,irtc)
        if(irtc .ne. 0)then
          isp=isp2
          return
        endif
        isp=isp0
        if(ktfsequenceq(kxi,kl))then
          call tfgetllstkall(kl)
        else
          isp=isp+1
          ktastk(isp)=kxi
        endif
      enddo
      kx=kxmakelist0(isp2,klx)
      if(mult)then
        klx%head=ktfoper+mtfnull
      endif
      irtc=0
      return
      end

      subroutine tfreadstring(isp1,kx,char1,del,irtc)
      use tfstk
      use tfrbuf
      use readopt
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfgetlfn,lfn
      logical*4 char1,del
      type(ropt)  opts
      lfn=itfgetlfn(isp1,.true.,irtc)
      call tfreadoptions(isp1,opts,del,irtc)
      if(irtc .eq. 0)then
        call tfreadstringf(lfn,kx,char1,opts,irtc)
      endif
      return
      end

      subroutine tfreadoptions(isp1,opts,del,irtc)
      use tfstk
      use tfrbuf
      use readopt
      implicit none
      integer*4 isp1,irtc,isp0,ispopt,itfmessage
      logical*4 del
      character*64 tfgetstr
      integer*4 nopt
      parameter (nopt=3)
      integer*8 kaopt(nopt)
      character*14 optname(nopt)
      data kaopt/nopt*0/
      data optname/
     $     'WordSeparators',
     $     'ReadNewRecord ',
     $     'NullWords     '/
      type(ropt)  opts
      opts%null=.false.
      opts%new=.true.
      if(del)then
        opts%ndel=3
        opts%delim=' ,'//char(9)
      else
        opts%ndel=0
      endif
      if(isp .gt. isp1+2)then
        isp0=isp
        call tfgetoptionstk(isp1+3,kaopt,optname,nopt,ispopt,irtc)
        isp=isp0
        if(irtc .ne. 0)then
          return
        endif
        if(ispopt .ne. isp1+3)then
          irtc=itfmessage(9,'General::narg','"2 (+ options)"')
          return
        endif
        opts%opt=.false.
        if(ktfstringq(ktastk(isp0+1)))then
          opts%delim=tfgetstr(ktastk(isp0+1),opts%ndel)
          opts%opt=.true.
        elseif(ktastk(isp0+1) .ne. ktfref)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"Character-string is","for WordSeparators ->"')
          return
        endif
        if(ktfrealq(ktastk(isp0+2)))then
          opts%new=rtastk(isp0+2) .ne. 0.d0
          opts%opt=.true.
        elseif(ktastk(isp0+2) .ne. ktfref)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"True or False is","for ReadNewRecord ->"')
          return
        endif
        if(ktfrealq(ktastk(isp0+3)))then
          opts%null=rtastk(isp0+3) .ne. 0.d0
          opts%opt=.true.
        elseif(ktastk(isp0+3) .ne. ktfref)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"True or False is","for NullWords ->"')
          return
        endif
        if(.not. opts%opt)then
          irtc=itfmessage(9,'General::wrongopt',' ')
        endif
      endif
      return
      end

      subroutine tfreadstringf(lfn,kx,char1,opts,irtc)
      use tfstk
      use tfrbuf
      use tfcsi
      use readopt
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*8 ib,is
      integer*4 irtc,lfn,isw,next,nc1,nc
      logical*4 char1
      type(ropt) opts
      irtc=0
      isw=1
      if(lfn .le. 0 .or. lfn .eq. icslfni())then
        call tfreadstringfb(lfn,kx,char1,opts,irtc)
        return
      elseif(itbuf(lfn) .le. 2)then
        call tfreadstringfb(lfn,kx,char1,opts,irtc)
        return
      endif
 20   call tfreadbuf(irbgetpoint,lfn,ib,is,nc,' ')
 10   if(nc .lt. 0)then
        call tfreadbuf(irbreadrecordbuf,lfn,ib,is,nc,' ')
        if(nc .ge. 0)then
          go to 20
        else
          go to 101
        endif
      else
        if(nc .eq. 0)then
          if(lfn .eq. 0)then
            kx%k=kxeof
            return
          endif
          if(opts%new)then
            if((opts%ndel .gt. 0 .and. .not. opts%null) .or. char1)then
              nc=-1
              go to 10
            endif
            call tfreadbuf(irbeor2bor,lfn,int8(0),int8(0),nc,' ')
          endif
          kx=dxnulls
          return
        else
          isw=1
          if(char1)then
            nc1=1
            next=2
          elseif(opts%ndel .gt. 0)then
            call tfword(jlist(is,ib),nc,opts%delim(1:opts%ndel),
     $           opts%ndel,isw,nc1,next,opts%null)
            if(nc1 .le. 0 .and. .not. opts%null)then
              if(opts%new)then
                nc=-1
                go to 10
              else
                nc1=0
                next=nc+1
              endif
            endif
          else
            nc1=nc
            next=nc+1
          endif
          if(opts%new .and. next .gt. nc)then
            call tfreadbuf(irbbor,lfn,int8(0),int8(0),nc,' ')
          else
            call tfreadbuf(irbmovepoint,lfn,int8(0),int8(0),next-1,' ')
          endif
          nc=nc1
        endif
        call loc_sad(ib-1,str)
        kx=kxsalocb(-1,str%str(is+isw-1:is+isw-2+nc),nc)
        return
      endif
      return
 101  kx%k=kxeof
      return
      end

      subroutine tfreadstringfb(lfn,kx,char1,opts,irtc)
      use tfstk
      use tfrbuf
      use tfcsi
      use readopt
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*8 ib,is
      integer*4 irtc,lfn,isw,next,nc1,nc
      character*(maxlbuf) buff
      logical*4 char1
      type(ropt) opts
      irtc=0
      isw=1
 20   call tfreadbuf(irbgetpoint,lfn,ib,is,nc,' ')
 10   if(nc .lt. 0)then
        if(lfn .gt. 0 .and. lfn .ne. icslfni())then
          call tfreadbuf(irbreadrecordbuf,lfn,ib,is,nc,buff)
          if(nc .ge. 0)then
            go to 20
          else
            go to 101
          endif
        else
          call savebuf(buff,nc)
          if(nc .le. 0)then
            if(lfn .eq. 0)then
              go to 101
            endif
            if(opts%new)then
              call tprmpt(-1,-1,0)
              call getbuf
              if(icsstat() .ne. 0)then
                go to 101
              endif
              call tprmpt(0,-1,0)
              call savebuf(buff,nc)
            else
              nc=0
              go to 1
            endif
          endif
          if(char1)then
            nc=1
            call cssetp(icsmrk()+1)
          elseif(opts%ndel .gt. 0)then
c            write(*,*)'tfreadstring-1 ',opts%delim(1:opts%ndel),'$ ',
c     $           opts%ndel,opts%null
            call tfword(buff(1:nc),nc,opts%delim(1:opts%ndel),
     $           opts%ndel,isw,nc1,next,opts%null)
            if(next .gt. 0)then
              call cssetp(icsmrk()+next-1)
            else
              call cssetp(icsmrk()+nc)
            endif
            nc=nc1
          else
            call cssetp(icsmrk()+nc)
          endif
        endif
        go to 1
 101    kx%k=kxeof
        return
      else
        if(nc .eq. 0)then
          if(lfn .eq. 0)then
            kx%k=kxeof
            return
          endif
          if(opts%new)then
            if((opts%ndel .gt. 0 .and. .not. opts%null) .or. char1)then
              nc=-1
              go to 10
            endif
            call tfreadbuf(irbeor2bor,lfn,int8(0),int8(0),nc,' ')
          endif
          kx=dxnulls
          return
        else
          isw=1
          if(char1)then
            nc1=1
            next=2
          elseif(opts%ndel .gt. 0)then
            call tfword(jlist(is,ib),nc,opts%delim(1:opts%ndel),
     $           opts%ndel,isw,nc1,next,opts%null)
            if(nc1 .le. 0 .and. .not. opts%null)then
              if(opts%new)then
                nc=-1
                go to 10
              else
                nc1=0
                next=nc+1
              endif
            endif
          else
            nc1=nc
            next=nc+1
          endif
          if(opts%new .and. next .gt. nc)then
            call tfreadbuf(irbbor,lfn,int8(0),int8(0),nc,' ')
          else
            call tfreadbuf(irbmovepoint,lfn,int8(0),int8(0),next-1,' ')
          endif
          nc=nc1
        endif
        call loc_sad(ib-1,str)
        kx=kxsalocb(-1,str%str(is+isw-1:is+isw-2+nc),nc)
c        kx=kxsalocb(-1,jlist(is+isw-1,ib),nc)
        return
      endif
 1    kx=kxsalocb(-1,buff(isw:),nc)
      return
      end

      subroutine tfreadbyte(isp1,lfn,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,lfn,itfmessage,fgetc,mode,narg,
     $     isp0,ibuf(2),jbuf(4),i,nb,          ispopt
      real*8 rbuf,vx
      real*4 sbuf(2)
      character bbuf(8)
      equivalence (bbuf,rbuf),(rbuf,ibuf),(rbuf,jbuf),(rbuf,sbuf)
      logical*4 opt,little,tfsamesymbolqd
      integer*4 nopt
      parameter (nopt=1)
      integer*8 kaopt(nopt)
      character*9 optname(nopt)
      data kaopt/nopt*0/
      type (sad_descriptor), save :: itflittle,itfbig
      data itflittle%k,itfbig%k /0,0/
      data optname/
     $     'ByteOrder'/
      narg=isp-isp1
      irtc=0
      if(lfn .eq. 0)then
        kx%k=kxeof
        return
      endif
      if(itflittle%k .eq. 0)then
        itflittle=kxsymbolf('LittleEndian',12,.true.)
        itfbig=kxsymbolf('BigEndian',9,.true.)
      endif
      little=.false.
      if(narg .gt. 2)then
        isp0=isp
        call tfgetoptionstk(isp1+3,kaopt,optname,nopt,ispopt,irtc)
        isp=isp0
       if(irtc .ne. 0)then
          return
        endif
        if(ispopt .ne. isp1+3)then
          irtc=itfmessage(9,'General::narg','"2 (+ options)"')
          return
        endif
        opt=.false.
        if(iand(ktfmask,ktastk(isp0+1)) .eq. ktfsymbol)then
          if(tfsamesymbolqd(dtastk(isp0+1),itflittle))then
            little=.true.
            opt=.true.
          elseif(tfsamesymbolqd(dtastk(isp0+1),itfbig))then
            little=.false.
            opt=.true.
          else
            irtc=itfmessage(9,'General::wrongval',
     $           '"LittleEndian or BigEndian","for ByteOrder ->"')
          endif
        elseif(iand(ktfmask,ktastk(isp0+1)) .ne. ktfoper)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"LittleEndian or BigEndian","for ByteOrder ->"')
          return
        endif
        if(.not. opt)then
          irtc=itfmessage(9,'General::wrongopt',' ')
          return
        endif
      endif
      nb=mode
      if(mode .eq. 3)then
        nb=4
      elseif(mode .eq. 5)then
        nb=8
      endif
      rbuf=0.d0
      if(little)then
        do i=nb,1,-1
          irtc=fgetc(lfn,bbuf(i))
          if(irtc .ne. 0)then
            go to 101
          endif
        enddo
      else
        do i=1,nb
          irtc=fgetc(lfn,bbuf(i))
          if(irtc .ne. 0)then
            go to 101
          endif
        enddo
      endif
      go to (10,20,30,40,50),mode
 10   vx=ichar(bbuf(1))
      go to 1000
 20   vx=jbuf(1)
      go to 1000
 30   vx=ibuf(1)
      go to 1000
 40   vx=sbuf(1)
      go to 1000
 50   vx=rbuf
 1000 kx=dfromr(vx)
      return
 101  kx%k=kxeof
      irtc=0
      return
      end

      subroutine tfword(str,nc,del,ndel,is,nw,next,null)
      implicit none
      integer*4 nw,is,nc,i,next,ndel
      character*(ndel) del
      logical*4 null
      integer*1 str(nc)
c      write(*,*)'tfword ',len(del)
      if(null)then
        is=1
        if(index(del,char(str(1))) .gt. 0)then
          next=2
          nw=0
          return
        endif
      else
        do i=1,nc
          if(index(del,char(str(i))) .le. 0)then
            is=i
            go to 10
          endif
        enddo
        is=nc+1
        nw=0
        next=nc+1
        return
      endif
 10   do i=is+1,nc
        if(index(del,char(str(i))) .gt. 0)then
c          write(*,*)'tfword ',i,' $',del,'$',
c     $         char(str(i)),'$'
          nw=i-is
          next=i+1
          return
        endif
      enddo
      next=nc+1
      nw=next-is
      return
      end

      subroutine tftoexpression(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,nc,itfmessage
      if(isp .le. isp1)then
        kx%k=ktfoper+mtfnull
        irtc=0
        return
      endif
      if(isp .gt. isp1+1 .or.
     $     .not. ktfstringqd(dtastk(isp),str))then
        kx=dtastk(isp)
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      nc=str%nch
      irtc=0
      if(nc .le. 0)then
        kx%k=ktfoper+mtfnull
        return
      endif
      call tfevalb(str%str(1:nc),nc,kx,irtc)
      return
      end

      subroutine tfopenread(isp1,kx,irtc)
      use tfstk
      use tfcsi
      use tfrbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*8 ka,kfromr,kfile,mapallocfile,ksize
      integer*4 irtc,isp1,ifile,lfn,itfmessage,nc
      logical*4 disp
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfstringq(ktastk(isp)))then
        ka=ktfaddr(ktastk(isp))
        call loc_string(ka,str)
        nc=str%nch
        if(nc .le. 0)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"filename or \"!command\"","argument"')
          return
        endif
        disp=.false.
        if(str%str(1:1) .eq. '!')then
          call tfsystemcommand(str%str,nc,kx,irtc)
          if(irtc .ne. 0)then
            return
          endif
          disp=.true.
c          if(str%str(nc:nc) .eq. '&')then
c            lfn=itfopenread(kx,disp,irtc)
c            kx%k=kfromr(dble(lfn))
c            return
c          endif
          call loc_string(ktfaddr(kx),str)
        else
          if(nc .gt. 2) then
             if (str%str(nc-1:nc) .eq. '.z' .or.
     $           str%str(nc-1:nc) .eq. '.Z')then
               call tfuncompress(str%str,nc,kx,irtc)
               if(irtc .ne. 0)then
                 return
               endif
               call loc_string(ktfaddr(kx),str)
               disp=.true.
             endif
          endif
          if(nc .gt. 3) then
             if (str%str(nc-2:nc) .eq. '.gz')then
               call tfungzip(str%str,nc,kx,irtc)
               if(irtc .ne. 0)then
                 return
               endif
               call loc_string(ktfaddr(kx),str)
               disp=.true.
             endif
          endif
          if(nc .gt. 4) then
             if (str%str(nc-3:nc) .eq. '.bz2')then
               call tfunbzip2(str%str,nc,kx,irtc)
               if(irtc .ne. 0)then
                 return
               endif
               call loc_string(ktfaddr(kx),str)
               disp=.true.
             endif
          endif
        endif
        kfile=mapallocfile(str%str,ifile,ksize,irtc)
        if(irtc .ne. 0)then
          irtc=itfmessage(999,'General::fileopen',str%str(1:nc))
          kx=dxfailed
        else
          call tfreadbuf(irbopen,lfn,kfile/8,ksize+modemapped,ifile,' ')
c          write(*,*)'tfopenread ',
c     $         kfile,ifile,ksize,lfn,irtc,' ',
c     $         str%str(1:nc)
          kx%k=kfromr(dble(lfn))
        endif
      endif
      return
      end

      subroutine tfuncompress(cmd,nc,kx,irtc)
      use tfstk
      implicit none 
      type (sad_descriptor) kx
      integer*4 irtc,nc
      character*2048 cmd
c      call tfsystemcommand('!uncompress -c '//cmd(:nc),
c     $     nc+15,itx,iax,vx,irtc)
      call tfsystemcommand('!cat '//cmd(:nc)//'|uncompress -c',
     $     nc+19,kx,irtc)
      return
      end

      subroutine tfungzip(cmd,nc,kx,irtc)
      use tfstk
      implicit none 
      type (sad_descriptor) kx
      integer*4 irtc,nc
      character*2048 cmd
      call tfsystemcommand('!gzip -dc '//cmd(:nc),
     $     nc+10,kx,irtc)
      return
      end

      subroutine tfunbzip2(cmd,nc,kx,irtc)
      use tfstk
      implicit none 
      type (sad_descriptor) kx
      integer*4 irtc,nc
      character*2048 cmd
      call tfsystemcommand('!bzip2 -dc '//cmd(:nc),
     $     nc+11,kx,irtc)
      return
      end

      subroutine tfsystemcommand(cmd,nc,kx,irtc)
      use tfrbuf
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor), save :: ntable(1000)
      integer*4 irtc,nc,system,itfsyserr,lw,ir,i,m
      integer*4 l
      character post
      character*2048 cmd
      character*2048 buff,tfgetstr
      data ntable(:)%k /1000*0/
      l=nextfn(moderead)
      if(ntable(l)%k .eq. 0)then
        call tftemporaryname(isp,kx,irtc)
        ntable(l)=dtfcopy(kx)
      else
        kx=ntable(l)
        irtc=0
      endif
      buff=tfgetstr(kx,lw)
      post=' '
      m=nc
      do i=nc,1,-1
        m=i
        if(cmd(i:i) .eq. '&')then
          post='&'
          m=i-1
          exit
        elseif(cmd(i:i) .ne. ' ')then
          exit
        endif
      enddo
      ir=system(cmd(2:m)//' > '//buff(:lw)//' '//post)
c      write(*,*)'tfsyscmd ',ir,' ',cmd(2:m),' ',buff(:lw)
      if(ir .lt. 0)then
        irtc=itfsyserr(999)
        return
      endif
      call tpause(10000)
      if(post .eq. '&')then
        call tpause(10000)
      endif
c      ir=system('chmod 777 '//buff(:lw)//char(0))
c      if(ir .lt. 0)then
c        irtc=itfsyserr(999)
c        return
c      endif
      return
      end

      subroutine tfstringtostream(isp1,kx,irtc)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) kx
      integer*8 kfromr
      integer*4 isp1,irtc,itfmessage,iu,nc
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonstringq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Character-string"')
        return
      endif
      call tfreadbuf(irbopen,iu,ktfaddr(ktastk(isp)),
     $     modestring,nc,' ')
      if(iu .le. 0)then
        irtc=itfmessage(9,'General::fileopen','"(String)"')
      else
        kx%k=kfromr(dble(iu))
        irtc=0
      endif
      irtc=0
      return
      end

      subroutine tfclosef(k,irtc)
      use tfstk, only:ktfnonrealqdi
      use tfcode
      use tfrbuf
      implicit none
      type (sad_descriptor) k
      integer*4 irtc,iu,itfmessage,nc
      if(ktfnonrealqdi(k,iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
        return
      endif
      call tfreadbuf(irbclose,iu,int8(0),int8(0),nc,' ')
      irtc=0
      return
      end

      subroutine tftemporaryname(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,lw,i,tftmpnam,getpid,lm,itfmessage
      character*4096 buff
      character*256 machine
      data lm /0/
      save machine
      type (sad_descriptor) ixmachine
      data ixmachine%k /0/
      if(isp .eq. isp1+1)then
        if(ktastk(isp) .ne. kxnull)then
          irtc=itfmessage(9,'General::wrongval','Null')
          return
        endif
      elseif(isp1 .ne. isp)then
        irtc=itfmessage(9,'General::narg','0')
        return
      endif
      if(ixmachine%k .eq. 0)then
        ixmachine=kxsymbolz('System`$MachineName',19)
      endif
      if(lm .eq. 0)then
        call descr_sad(ixmachine,symd)
        call descr_sad(symd%value,str)
        lm=min(256,str%nch)
        machine(1:lm)=str%str(1:lm)
        do i=1,lm
          if(machine(i:i) .eq. '.')then
            lm=i-1
            exit
          endif
        enddo
      endif
      buff(1:9)='/tmp/tmp_'
      buff(10:10+lm)=machine(:lm)//"_"
      write(buff(11+lm:18+lm),'(I8.8)')getpid()
      buff(19+lm:25+lm)='.XXXXXX'
      buff(26+lm:26+lm)=char(0)
      lw=tftmpnam(buff)
      kx=kxsalocb(-1,buff,lw)
      irtc=0
      return
      end

      subroutine tfseek(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      integer*4 isp1,irtc,lfn,ioff,ls,   ierr
      logical*4 tfsamesymbolqd,tfsamesymbolqk
      integer*4 fseek,itfmessage,itfsyserr
      external fseek
      type (sad_descriptor), save :: iabof,iacp
      data iabof%k,iacp%k/0,0/
      if(iabof%k .eq. 0)then
        iabof=kxsymbolf('BeginningOfFile',15,.true.)
        iacp= kxsymbolf('CurrentPosition',15,.true.)
      endif
      if(isp .eq. isp1+2)then
        ls=1
      elseif(isp .ne. isp1+3)then
        go to 9000
      else
        k1=dtastk(isp)
        if(.not. ktfsymbolqd(k1))then
          go to 9000
        endif
        if(tfsamesymbolqk(k1,kxeof))then
          ls=2
        elseif(tfsamesymbolqd(k1,iabof))then
          ls=0
        elseif(tfsamesymbolqd(k1,iacp))then
          ls=1
        else
          go to 9000
        endif
      endif
      if(.not. ktfrealqd(dtastk(isp1+1)) .or.
     $     .not. ktfrealqd(dtastk(isp1+2)))then
        go to 9000
      endif
      lfn=int(rtastk(isp1+1))
      ioff=int(rtastk(isp1+2))
      ierr=fseek(lfn,ioff,ls)
      if(ierr .ne. 0)then
        irtc=itfsyserr(9)
      endif
      kx%k=ktfoper+mtfnull
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"Seek[fileno,offset,'//
     $       'EndOfFile|BeginningOfFile|CurrentPosition]"')
      return
      end

      subroutine tfflush(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,narg,lfn,itfmessage,itfgetlfn
      narg=isp-isp1
      if(narg .ne. 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      lfn=itfgetlfn(isp1,.false.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      endfile(lfn,err=1)
 1    call flush(lfn)
      backspace(lfn,err=10)
 10   irtc=0
      kx%k=ktfoper+mtfnull
      return
      end
