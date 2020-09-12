      module readopt
        type ropt
        sequence
        character*64 delim
        integer*4 ndel
        logical*4 new,null,del,opt
        end type
      end module

      function tfwrite(isp1,irtc) result(kx)
      use tfstk
      use strbuf
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_strbuf), pointer :: strb
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfgetrecl,narg,lfn,isp11,
     $     j,i,nc,lpw,itfmessage,isp2,itfgetlfn
      logical*4 exist
      kx=dxnullo
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

      function tfmapfile(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_string), pointer :: fname
      type (sad_rlist), pointer :: kl
      integer*8 map,ksize,maprwfile
      integer*4 ifd,itfmessage
      kx=dxnullo;
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
c      call tfdebugprint(dtastk(isp1+1),'mapfile',1)
      if(.not. ktfstringq(dtastk(isp1+1),fname))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Filename for #1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp),ksize))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Bytes/8 for #2"')
        return
      endif
      ksize=ksize*8
      map=maprwfile(fname%str,ifd,ksize,irtc)
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'General::mmap','""')        
        return
      endif
      kx=kxraaloc(-1,2,kl)
      kl%rbody(1)=dble(map)
      kl%rbody(2)=dble(ifd)
      return
      end

      function tfunmapfile(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 map,ksize
      integer*4 ifd,itfmessage,unmap
      kx=dxnullo;
      if(isp .ne. isp1+3)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+1),map))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Address for #1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+2),ksize))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Bytes/8 for #2"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+1),ifd))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"File descriptor for #3"')
        return
      endif
      ksize=ksize*8
      irtc=unmap(map,ksize,ifd)
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'General::mmap','"(Unmap)"')
        return
      endif
      return
      end

      integer*4 function itfgetlfn(isp1,read,irtc) result(iv)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) k
      type (sad_dlist), pointer :: kl
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      logical*4 ,intent(in):: read
      iv=0
      if(isp .le. isp1)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      if(ktfrealq(ktastk(isp1+1),iv))then
        go to 100
      elseif(tflistq(dtastk(isp1+1),kl))then
        k=kl%dbody(merge(1,2,read))
        if(ktfrealq(k,iv))then
          go to 100
        endif
      endif
      irtc=itfmessage(9,'General::wrongtype',
     $     '"Real number or List of two Reals"')
      return
 100  irtc=0
      select case (iv)
      case (-1)
        iv=merge(lfni,lfno,read)
      case (0)
      case default
        if(iv .lt. 0)then
          irtc=itfmessage(9,'General::wrongnum','"positive, 0 or -1"')
        endif
      end select
      return
      end

      subroutine tfwritestring(isp1,kx,irtc)
      use tfstk
      use strbuf
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_strbuf), pointer :: strb
      integer*4 ,parameter ::mmax=1000000
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 narg,lfn,isp11,j,i,itfgetlfn,nc,itfmessage,isp2
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
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tfwrite
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0
      isp0=isp
      isp=isp+1
      rtastk(isp)=lfno
      ktastk(isp+1:isp+isp0-isp1)=ktastk(isp1+1:isp0)
      isp=isp+isp0-isp1
c      do i=1,isp0-isp1
c        isp=isp+1
c        ktastk(isp)=ktastk(isp1+i)
c      enddo
      kx=tfwrite(isp0,irtc)
      isp=isp0
      return
      end

      subroutine tfdebugprint(k,pr1,nline)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(in):: nline
      integer*4 l,itfdownlevel,irtc
      character*(*) ,intent(in):: pr1
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
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfgetrecl,nret,itfmessage
      if(isp .eq. isp1+1)then
        call tfprint1(dtastk(isp),
     $       lfno,-itfgetrecl(),1,.true.,.true.,irtc)
      elseif(isp .eq. isp1+2)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype','"real number"')
          return
        endif
        nret=int(rtastk(isp))
        call tfprint1(dtastk(isp1+1),
     $       lfno,-itfgetrecl(),abs(nret),nret .ge. 0,
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
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kxlongstr,ks,kr
      type (sad_strbuf), pointer :: strb
      integer*4 ,intent(in):: lfno,nret
      integer*4 ,intent(out):: irtc
      integer*4 nc,lrec,irtc1,isp0
      logical*4 ,intent(in):: cr,str
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
        if(tfsameq(kr,kxlongstr))then
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
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ki,k
      type (sad_symbol), pointer ::sym,sym1
      type (sad_symdef), pointer ::symd
      type (sad_dlist), pointer :: klh,kli,klx
      integer*8 ii,ka0,kadi,kad,ka
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,isp00,isp0,kk,iop,ispv,icmpl
      icmpl=0
      call tfgetoption('Compiled',ktastk(isp),kx,irtc)
      if(irtc .eq. -1)then
        ispv=isp
      elseif(irtc .ne. 0)then
        return
      elseif(ktfrealq(kx))then
        if(kx%k .ne. 0)then
          icmpl=1
        endif
        ispv=isp-1
      else
        ispv=isp
      endif
      isp00=isp
      do i=isp1+1,ispv
        isp0=isp
        k=dtastk(i)
        if(ktfoperq(k,ka))then
          ka=klist(ifunbase+ka)
          k%k=ktfsymbol+ka
        endif
        if(ktfsymbolq(k,sym1))then
          sym=>tfsydef(sym1)
          call sym_symdef(sym,symd)
          ka0=ksad_loc(symd%upval)
          if(ktfrefq(symd%value,ka))then
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
          klh%head%k=ktfoper+mtfsetdelayed
          klh%dbody(1)%k=ktfsymbol+ktfcopy1(ka0+6)
          klh%dbody(2)=dtfcopy(k)
          do kk=1,0,-1
            iop=merge(mtfsetdelayed,mtfupsetdelayed, kk .eq. 1)
            kad=klist(ka0+kk)
            do while(kad .ne. 0)
              if(ilist(1,kad+2) .eq. maxgeneration)then
                do ii=kad+3,kad+ilist(2,kad+2)+3
                  kadi=klist(ii)
                  do while(kadi .ne. 0)
                    isp=isp+1
                    dtastk(isp)=kxadaloc(-1,2,kli)
                    kli%head%k=ktfoper+iop
                    kli%dbody(1)=dtfcopy1(dlist(kadi+3+icmpl))
                    kli%dbody(2)=dtfcopy(dlist(kadi+5+icmpl))
                    if(kli%dbody(2)%k .eq. ktfref)then
                      kli%dbody(2)=dtfcopy(dlist(kadi+6))
                    endif
                    kadi=klist(kadi)
                  enddo
                enddo
              else
                isp=isp+1
                dtastk(isp)=kxadaloc(-1,2,kli)
                kli%head%k=ktfoper+iop
                kli%dbody(1)=dtfcopy1(dlist(kad+3+icmpl))
                kli%dbody(2)=dtfcopy(dlist(kad+5+icmpl))
                if(kli%dbody(2)%k .eq. ktfref)then
                  kli%dbody(2)=dtfcopy(dlist(kad+6))
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
      klx%head%k=ktfoper+mtfhold
      isp=isp00
      irtc=0
      return
      end

      subroutine tfgetf(fname)
      use tfstk
      use tfcode
      implicit none
      character*(*) ,intent(in):: fname
      type (sad_descriptor) kx,k,tfget
      integer*4 irtc
      k=kxsalocb(-1,fname,len_trim(fname))
      kx=tfget(k,irtc)
      return
      end

      recursive function tfget(k,irtc) result(kx)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) kx,kf,kfn
      type (csiparam) sav
      integer*4 ,intent(out):: irtc
      integer*4 itfgeto,lfn,isp0,itf
      isp0=isp
      isp=isp+1
      dtastk(isp)=k
c      call tfdebugprint(k,'tfget',1)
      call tfopenread(isp0,kfn,irtc)
      isp=isp0
      if(irtc .ne. 0)then
        kx%k=ktfoper+mtfnull
        return
      endif
      lfn=int(rfromd(kfn))
      sav=savep
      rep=.false.
      levele=levele+1
      kx%k=ktfoper+mtfnull
      call trbassign(lfn)
      lfn1=0
      call skiplnget
      itf=0
      do while(itf .ge. 0)
        itf=itfgeto(kf)
        if(itf .ge. 0)then
          kx=kf
        elseif(itf .eq. -1)then
          itf=0
          call skiplnget
        endif
        if(ios .ne. 0)then
          ios=0
          itf=min(-1,itf)
        endif
      enddo
      call tfconnect(kx,0)
      call trbclose(lfn)
      savep=sav
      call trbassign(lfni)
      irtc=merge(irtcabort,0,itf .eq. -3)
      return
      end

      subroutine tfread1(lfn,kx)
      use tfstk
      use tfrbuf
      use tfcsi
      implicit none
      type (sad_descriptor) , intent(out)::kx
      type (csiparam) sav
      integer*4 , intent(in)::lfn
      integer*4 itfgeto,itf
      logical*4 openf
      if(lfn .le. 0)then
        kx%k=kxeof
        return
      endif
      sav=savep
      openf=lfn .ne. sav%lfni
      if(openf)then
        call trbassign(lfn)
        lfn1=0
      endif
      levele=levele+1
      itf=-1
      ios=0
      do while (itf .eq. -1 .and. ios .eq. 0)
        itf=itfgeto(kx)
        if(itf .eq. -1)then
          call tprmptget(-1,.true.)
        endif
      enddo
      if(itf .lt. 0)then
        kx=dxnullo
      endif
      if(ios .ne. 0)then
        ios=0
        kx%k=kxeof
        if(openf)then
          call trbreset(lfn)
        endif
      endif
      call tfconnect(kx,0)
      if(openf)then
        savep=sav
        call trbassign(lfni)
      endif
      return
      end

      recursive subroutine tfskip(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,narg,isp0,nrpt
      narg=isp-isp1
      if(narg .ge. 3 .and.
     $     ktfrealq(dtastk(isp1+3),nrpt))then
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
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 lfn,itfgetlfn
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
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1,lfn
      integer*4 ,intent(out):: irtc
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
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) k1,k2
      type (sad_dlist), pointer :: list
      integer*4 ,intent(in):: isp1,lfn
      integer*4 ,intent(out):: irtc
      integer*4 isp0,n1,narg,itfmessage,i1
      type(ropt) opts
      type (sad_descriptor)
     $     itfexprs,itfstrs,itfwords,itfreals,itfchars,
     $     itfbyte,itfbinint,itfbinsint,itfbinreal,itfbindble
      data itfexprs%k,itfstrs%k,itfwords%k,itfreals%k,itfchars%k,
     $     itfbyte%k,itfbinint%k,itfbinsint%k,itfbinreal%k,itfbindble%k
     $     /0,0,0,0,0, 0,0,0,0,0/
      narg=isp-isp1
      irtc=0
      if(narg .eq. 0)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
      elseif(narg .eq. 1)then
        call tfread1(lfn,kx)
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
        if(ktfsymbolq(k2))then
          if(tfsamesymbolq(k2,itfexprs))then
            isp0=isp
            isp=isp+1
            ktastk(isp)=ktastk(isp1+1)
            call tfread1(lfn,kx)
            isp=isp0
          elseif(tfsamesymbolq(k2,itfstrs))then
            isp0=isp
            isp=isp+1
            ktastk(isp)=ktastk(isp1+1)
            n1=opts%ndel
            opts%ndel=0
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
            opts%ndel=n1
            isp=isp0
          elseif(tfsamesymbolq(k2,itfwords))then
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
          elseif(tfsamesymbolq(k2,itfchars))then
            call tfreadstringf(lfn,kx,.true.,opts,irtc)
          elseif(tfsamesymbolq(k2,itfreals))then
            call tfreadstringf(lfn,kx,.false.,opts,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(ktfstringq(kx))then
              isp0=isp
              isp=isp+1
              dtastk(isp)=kx
              call tftoexpression(isp0,kx,irtc)
              isp=isp0
            endif
          elseif(tfsamesymbolq(k2,itfbyte))then
            call tfreadbyte(isp1,lfn,kx,1,irtc)
          elseif(tfsamesymbolq(k2,itfbinsint))then
            call tfreadbyte(isp1,lfn,kx,2,irtc)
          elseif(tfsamesymbolq(k2,itfbinint))then
            call tfreadbyte(isp1,lfn,kx,3,irtc)
          elseif(tfsamesymbolq(k2,itfbinreal))then
            call tfreadbyte(isp1,lfn,kx,4,irtc)
          elseif(tfsamesymbolq(k2,itfbindble))then
            call tfreadbyte(isp1,lfn,kx,5,irtc)
          else
            go to 9000
          endif
        elseif(ktflistq(k2,list))then
          if(list%head%k .eq. ktfoper+mtflist)then
            call tfreadfm(lfn,list,list%nl,opts,.false.,kx,irtc)
          elseif(list%head%k .eq. ktfoper+mtftimes .and.
     $           list%nl .eq. 2)then
            k1=list%dbody(1)
            if(ktfnonrealq(k1,i1))then
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
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kxi
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer ::kl,klx
      integer*4 ,intent(in):: lfn,m
      integer*4 ,intent(out):: irtc
      integer*4 isp0,isp2,kk
      logical*4 ,intent(in):: mult
      type(ropt) ,intent(in):: opts
      isp2=isp
      do kk=1,m
        isp0=isp
        isp=isp0+2
        dtastk(isp)=list%dbody(merge(2,kk,mult))
        call tfreadfs(isp0,lfn,opts,kxi,irtc)
        if(irtc .ne. 0)then
          isp=isp2
          return
        endif
        isp=isp0
        if(ktfsequenceq(kxi%k,kl))then
          call tfgetllstkall(kl)
        else
          isp=isp+1
          dtastk(isp)=kxi
        endif
      enddo
      kx=kxmakelist0(isp2,klx)
      if(mult)then
        klx%head%k=ktfoper+mtfnull
      endif
      irtc=0
      return
      end

      subroutine tfreadstring(isp1,kx,char1,del,irtc)
      use tfstk
      use tfrbuf
      use readopt
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfgetlfn,lfn
      logical*4 ,intent(in):: char1,del
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
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,ispopt,itfmessage
      logical*4 ,intent(in):: del
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
      type(ropt)  ,intent(out):: opts
      opts%null=.false.
      opts%new=.true.
      if(del)then
        opts%ndel=5
        opts%delim=' ,'//char(10)//char(9)//char(13)
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
        if(ktfstringq(dtastk(isp0+1)))then
          opts%delim=tfgetstr(dtastk(isp0+1),opts%ndel)
          opts%opt=.true.
        elseif(ktastk(isp0+1) .ne. ktfref)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"Character-string is","for WordSeparators ->"')
          return
        endif
        if(ktfrealq(dtastk(isp0+2)))then
          opts%new=rtastk(isp0+2) .ne. 0.d0
          opts%opt=.true.
        elseif(ktastk(isp0+2) .ne. ktfref)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"True or False is","for ReadNewRecord ->"')
          return
        endif
        if(ktfrealq(dtastk(isp0+3)))then
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
      use readbuf
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (csiparam) sav
      integer*4 ,intent(in):: lfn
      integer*4 ,intent(out):: irtc
      integer*4 ie,is
      integer*4 isw,next,nc1,nc
      logical*4 ,intent(in):: char1
      logical*4 fb
      type (ropt) opts
      irtc=0
      if(lfn .le. 0 .or. ibuf(lfn) .eq. 0)then
        kx%k=kxeof
        return
      endif
      isw=1
      sav=savep
      if(lfn .ne. lfni)then
        call trbassign(lfn)
      endif
      fb=itbuf(lfn) .lt. modestring
      nc=sign(abs(lrecl-ipoint)+1,lrecl-ipoint)
      is=ipoint
      do while (.true.)
        if(nc .lt. 0)then
          if(.not. fb .or. lfn .ne. lfni)then
            call tfreadbuf(lfn,1,nc)
            if(nc .lt. 0)then
              go to 101
            endif
          else
            if(opts%new)then
              do while (nc .lt. 0)
                call tprmptget(-1,.false.)
                if(ios .ne. 0)then
                  go to 101
                endif
                nc=max(0,lrecl-ipoint)
              enddo
            else
              nc=0
              go to 1
            endif
          endif
        endif
        if(nc .eq. 0)then
          if(opts%new)then
            if((opts%ndel .gt. 0 .and. .not. opts%null) .or. char1)then
              nc=-1
              cycle
            endif
            call trbeor2bor(lfn)
          endif
          kx=dxnulls
          go to 1000
        endif
        if(char1)then
          nc1=1
          next=2
        elseif(opts%ndel .gt. 0)then
          call tfword(buffer(ipoint:ipoint+nc-1),
     $         opts%delim(1:opts%ndel),
     $         isw,nc1,next,opts%null)
          if(nc1 .le. 0 .and. .not. opts%null)then
            if(opts%new)then
              nc=-1
              cycle
            else
              nc1=0
              next=nc+1
            endif
          endif
        else
          nc1=nc
          next=nc+1
        endif
        exit
      enddo
      is=ipoint
      if(opts%new .and. next .gt. nc)then
        call trbnextl(lfn)
      else
        call trbmovepoint(lfn,next-1)
      endif
      nc=nc1
      ie=is+isw-2+nc
      if(buffer(ie:ie) .eq. char(10))then
        nc=max(nc-1,0)
      endif
 1    kx=kxsalocb(-1,buffer(is+isw-1:),nc)
      go to 1000
 101  kx%k=kxeof
 1000 if(lfn .ne. sav%lfni)then
        savep=sav
        call trbassign(sav%lfni)
      endif
      return
      end

      subroutine tfreadbyte(isp1,lfn,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1,lfn,mode
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,fgetc,narg,
     $     isp0,ibuf(2),jbuf(4),i,nb,ispopt
      real*8 rbuf,vx
      real*4 sbuf(2)
      character bbuf(8)
      equivalence (bbuf,rbuf),(rbuf,ibuf),(rbuf,jbuf),(rbuf,sbuf)
      logical*4 opt,little
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
          if(tfsamesymbolq(dtastk(isp0+1),itflittle))then
            little=.true.
            opt=.true.
          elseif(tfsamesymbolq(dtastk(isp0+1),itfbig))then
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

      subroutine tfword(str,del,is,nw,next,null)
      implicit none
      integer*4 ,intent(out):: is,nw,next
      integer*4 nc,i
      character*(*) ,intent(in):: str,del
      logical*4 ,intent(in):: null
      nc=len(str)
      if(null)then
        is=1
        if(index(del,str(1:1)) .gt. 0)then
          next=2
          nw=0
          return
        endif
      else
        do i=1,nc
          if(index(del,str(i:i)) .le. 0)then
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
        if(index(del,str(i:i)) .gt. 0)then
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
      type (sad_descriptor) ,intent(out):: kx
      type (sad_string), pointer :: str
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 nc,itfmessage
      if(isp .le. isp1)then
        kx%k=ktfoper+mtfnull
        irtc=0
        return
      endif
      if(isp .gt. isp1+1 .or.
     $     .not. ktfstringq(dtastk(isp),str))then
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
      call tfevalb(str%str(1:nc),kx,irtc)
      return
      end

      subroutine tfopenread(isp1,kx,irtc)
      use tfstk
      use tfcsi
      use readbuf, only:trbopenmap
      use tfrbuf
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_string), pointer :: str
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,nc,itfmessagestr
      logical*4 disp
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfstringq(dtastk(isp),str))then
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
        call trbopenmap(str%str(1:nc),kx,irtc)
        if(Irtc .ne. 0)then
          irtc=itfmessagestr(999,'General::fileopen',str%str(1:nc))
          kx=dxfailed
        endif
      endif
      return
      end

      subroutine tfuncompress(cmd,nc,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: nc
      integer*4 ,intent(out):: irtc
      character*(*) ,intent(in):: cmd
c      call tfsystemcommand('!uncompress -c '//cmd(:nc),
c     $     nc+15,itx,iax,vx,irtc)
      call tfsystemcommand('!cat '//cmd(:nc)//'|uncompress -c',
     $     nc+19,kx,irtc)
      return
      end

      subroutine tfungzip(cmd,nc,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: nc
      integer*4 ,intent(out):: irtc
      character*(*) ,intent(in):: cmd
      call tfsystemcommand('!gzip -dc '//cmd(:nc),
     $     nc+10,kx,irtc)
      return
      end

      subroutine tfunbzip2(cmd,nc,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: nc
      integer*4 ,intent(out):: irtc
      character*(*) ,intent(in):: cmd
      call tfsystemcommand('!bzip2 -dc '//cmd(:nc),
     $     nc+11,kx,irtc)
      return
      end

      subroutine tfsystemcommand(cmd,nc,kx,irtc)
      use tfrbuf
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) tftemporaryname
      integer*4 ,intent(in):: nc
      integer*4 ,intent(out):: irtc
      character*(*) ,intent(in):: cmd
      integer*4 , parameter :: llbuf=1024
      integer*4 system,itfsyserr,lw,ir,i,m
      integer*4 l
      character post
      character*(llbuf) buff,tfgetstr
      l=nextfn(moderead)
      if(ntable(l)%k .eq. 0)then
        kx=tftemporaryname(isp,irtc)
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
      use readbuf, only:trbopen
      use tfrbuf, only:modestring
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,iu,nc
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonstringq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Character-string"')
        return
      endif
      call trbopen(iu,ktfaddr(ktastk(isp)),int8(modestring),nc)
      if(iu .le. 0)then
        irtc=itfmessage(9,'General::fileopen','"(String)"')
      else
        kx=dfromr(dble(iu))
        irtc=0
      endif
      return
      end

      subroutine tfclosef(k,irtc)
      use tfstk, only:ktfnonrealq
      use tfcode
      use tfrbuf
      implicit none
      type (sad_descriptor) k
      integer*4 irtc,iu,itfmessage
      if(ktfnonrealq(k,iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
        return
      endif
      call trbclose(iu)
      irtc=0
      return
      end

      function tftemporaryname(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      integer*4 lw,i,tftmpnam,getpid,lm,itfmessage
      character*1024 buff
      character*256 machine
      data lm /0/
      save machine
      type (sad_descriptor) ixmachine
      data ixmachine%k /0/
      kx=dxnullo
      if(isp .eq. isp1+1)then
        if(ktastk(isp) .ne. ktfoper+mtfnull)then
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
c      write(*,*)'temporaryName ',buff(:lw)
      kx=kxsalocb(-1,buff,lw)
      irtc=0
      return
      end

      subroutine tfseek(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k1
      integer*4 lfn,ioff,ls,   ierr
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
        if(.not. ktfsymbolq(k1))then
          go to 9000
        endif
        if(tfsamesymbolq(k1%k,kxeof))then
          ls=2
        elseif(tfsamesymbolq(k1,iabof))then
          ls=0
        elseif(tfsamesymbolq(k1,iacp))then
          ls=1
        else
          go to 9000
        endif
      endif
      if(.not. ktfrealq(dtastk(isp1+1)) .or.
     $     .not. ktfrealq(dtastk(isp1+2)))then
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
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 narg,lfn,itfmessage,itfgetlfn
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
