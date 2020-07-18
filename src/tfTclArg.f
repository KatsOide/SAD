      module tclstrbuf

      contains

      subroutine tftclstringbufu(strb,str,nc,full)
      use tfstk
      use strbuf
      implicit none
      type (sad_strbuf), pointer :: strb
      integer*4 nc,m
      character*(nc) str
      character*(nc*3) buf
      logical*4 tfxftq,full
      if(tfxftq())then
        call tfuniconvs(str,nc,buf,m)
c        write(*,*)'unicode ',buf(1:m)
        call putstringbufb(strb,buf,m,full)
      else
        call putstringbufb(strb,str,nc,full)
      endif
      return
      end subroutine

      subroutine tftclstringbuf(strb,str,nc,full)
      use tfstk
      use strbuf
      implicit none
      type (sad_strbuf), pointer :: strb
      integer*4 nc,i,j
      character*(*) str
      character*(nc*3) buf
      character ch
      logical*4 full,tfxftq
      logical uniconv
      call putstringbufb1(strb,'"')
      uniconv=tfxftq()
      if(uniconv)then
        call tfuniconvs(str,nc,buf,j)
        j=j+1
        i=1
        do while(i .lt. j)
          ch=buf(i:i)
          i=i+1
          if(ch .eq. '"' .or. ch .eq. '$' .or. ch .eq. '\'
     $         .or. ch .eq. '[' .or. ch .eq. ']')then
            buf(i:j)=buf(i-1:j-1)
            buf(i-1:i-1)='\'
            j=j+1
            i=i+1
          endif
        enddo
      else
        j=1
        do i=1,nc
          ch=str(i:i)
          if(ch .eq. '"' .or. ch .eq. '$' .or. ch .eq. '\'
     $         .or. ch .eq. '[' .or. ch .eq. ']'
     $         .or. (ch .gt. char(159) .and. ch .lt. char(192)))then
            buf(j:j)='\'
            j=j+1
          endif
          buf(j:j)=ch
          j=j+1
        enddo
      endif
      call putstringbufb(strb,buf,j-1,full)
      call putstringbufb1(strb,'"')
      return
      end subroutine

      recursive function tftclarggen(isp1,strb,single,irtc)
     $     result(kx)
      use tfstk
      use strbuf
      implicit none
      type (sad_descriptor) kx,ki
      type (sad_dlist), pointer :: kl
      type (sad_strbuf), pointer :: strb
      type (sad_string), pointer :: str,stri,strj
      integer*8 kai,kti,kih,kaih,kaj
      integer*4 isp1,irtc,i,isp0,nc,j,isp3,isp00,isp2
      logical*4 single,full,ev,list,ini
      save isp00
      integer*8 ifwidget,iftkpathname,iflabel,ifstring
      data ifwidget,iftkpathname,iflabel,ifstring/0,0,0,0/

      if(ifwidget .eq. 0)then
        ifwidget=ktfsymbolz('Widget',6)
        iftkpathname=ktfsymbolz('TkPathName',10)
        iflabel=ktfsymbolz('TkOptionLabel',13)
        ifstring=ktfsymbolz('String',6)
      endif
      irtc=0
      if(isp1 .ge. isp)then
        kx=dxnulls
        return
      endif
      isp2=isp
      ini=.false.
      if(.not. associated(strb))then
        ini=.true.
        isp00=isp
        call getstringbuf(strb,0,.true.)
        if(single)then
          call putstringbufb1(strb,'{')
        endif
      endif
      LOOP_I: do i=isp1+1,isp2
        ki=dtastk(i)
        kai=ktfaddrd(ki)
        kti=ki%k-kai
        if(kti .eq. ktflist)then
          kih=klist(kai)
          kaih=ktfaddr(kih)
          if(kih .eq. ktfoper+mtfrule .or.
     $         kih .eq. ktfoper+mtfruledelayed)then
            isp3=isp
            isp=isp+1
            ktastk(isp)=ktfsymbol+iflabel
            isp=isp+1
            ktastk(isp)=klist(kai+1)
            call tfdeval(isp-1,iflabel,kx,1,.false.,ev,irtc)
            isp=isp3
            if(irtc .ne. 0)then
              go to 9000
            elseif(.not. ktfstringq(kx,str))then
              irtc=-1
              go to 9000
            endif
            call putstringbufb(strb,str%str,str%nch,full)
            isp0=isp
            isp=isp+1
            ktastk(isp)=klist(kai+2)
            list=tflistq(ktastk(isp))
            if(list)then
              call putstringbufb1(strb,'{')
            endif
            kx=tftclarggen(isp0,strb,single,irtc)
            if(irtc .gt. 0)then
              go to 9000
            endif
            if(list)then
              call putstringbufb1(strb,'}')
            endif
            isp=isp0
          elseif(kih .eq. ktfoper+mtflist)then
            if(ktfreallistq(kai))then
              do j=1,ilist(2,kai-1)
                call tfconvreal(strb,rlist(kai+j))
              enddo
            else
              isp0=isp
              call loc_dlist(kai,kl)
              call tfgetllstkall(kl)
c              call tfgetllstkall(klist(kai-3))
              kx=tftclarggen(isp0,strb,single,irtc)
              if(irtc .gt. 0)then
                go to 9000
              endif
              isp=isp0
            endif
          elseif(ktfsymbolq(kih))then
            if(tfsamesymbolq(kih,ifwidget))then
              isp0=isp
              isp=isp+1
              ktastk(isp)=ktfsymbol+iftkpathname
              isp=isp+1
              ktastk(isp)=klist(kai+1)
              call tfdeval(isp-1,iftkpathname,kx,1,.false.,ev,irtc)
              if(irtc .ne. 0)then
                go to 9000
              endif
              isp=isp0+1
              dtastk(isp)=kx
              kx=tftclarggen(isp0,strb,single,irtc)
              if(irtc .ne. 0)then
                go to 9000
              endif
              isp=isp0
            elseif(tfsamesymbolq(kih,ifstring))then
              isp0=isp
              call loc_dlist(kai,kl)
              call tfgetllstkall(kl)
c              call tfgetllstkall(klist(kai-3))
              do j=isp0+1,isp
                if(ktfstringq(ktastk(j)))then
                  kaj=ktfaddr(ktastk(j))
                  call loc_string(kaj,strj)
                  call tftclstringbufu(strb,strj%str,strj%nch,full)
                else                  
                  call tfconvstrb(strb,dtastk(j),nc,
     $                 .false.,.false.,-1,' ',irtc)
                  if(irtc .ne. 0)then
                    go to 9000
                  endif
                endif
              enddo
              isp=isp0
            endif
          endif
        elseif(kti .eq. ktfstring)then
          call loc_string(kai,stri)
          if(stri%nch .gt. 0)then
            if(single)then
              call tftclstringbufu(strb,stri%str,stri%nch,full)
            else
              call tftclstringbuf(strb,stri%str,stri%nch,full)
            endif
          elseif(.not. single)then
            call putstringbufb(strb,'{}',2,full)
          endif
          if(irtc .ne. 0)then
            go to 9000
          endif
        elseif(ktfrealq(ki))then
          call tfconvreal(strb,rtastk(i))
        elseif(ki%k .eq. ktfoper+mtfnull)then
          cycle LOOP_I
        else
          call tfconvstrb(strb,ki,nc,.false.,.true.,-1,' ',irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
        endif
        if(i .ne. isp2)then
          call putstringbufb1(strb,' ')
          if(irtc .ne. 0)then
            go to 9000
          endif
        endif
      enddo LOOP_I
      if(ini)then
        if(single)then
          call putstringbufb(strb,'} ',2,full)
        else
          call putstringbufb1(strb,' ')
        endif
        kx=kxstringbuftostring(strb)
        isp=isp00
      else
        if(.not. single)then
          call putstringbufb1(strb,' ')
        endif
      endif
      return
 9000 if(ini)then
        call tfreestringbuf(strb)
        isp=isp00
      endif
      return
      end function

      end module

      subroutine tftclarg(isp1,kx,irtc)
      use tfstk
      use strbuf
      use tclstrbuf
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_strbuf), pointer :: strb
      integer*4 isp1,irtc
      nullify(strb)
      kx=tftclarggen(isp1,strb,.false.,irtc)
      end

      subroutine tftclarg1(isp1,kx,irtc)
      use tfstk
      use strbuf
      use tclstrbuf
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_strbuf), pointer :: strb
      integer*4 isp1,irtc
      nullify(strb)
      kx=tftclarggen(isp1,strb,.true.,irtc)
      end

      recursive subroutine tfgetcanvasargstk(ka,irtc)
      use tfstk
      use strbuf
      use tclstrbuf
      implicit none
      type (sad_descriptor) :: k1,k2,kx,k2i
      type (sad_strbuf), pointer :: strb
      type (sad_string), pointer :: str
      type (sad_dlist), pointer ::kl,k2l
      integer*8 ktrsaloc,ka,ki,kai
      integer*4 n,i,isp3,irtc,nc
      logical*4 ev,full
      real*8 x
      integer*8 iflabel
      data iflabel/0/
      if(iflabel .eq. 0)then
        iflabel=ktfsymbolz('TkOptionLabel',13)
      endif
      call loc_sad(ka,kl)
      n=kl%nl
      if(kl%head%k .eq. ktfoper+mtflist)then
        if(ktfreallistq(kl))then
          rtastk(isp+1:isp+n)=kl%rbody(1:n)
          isp=isp+n
        else
          do i=1,n
            ki=kl%dbody(i)%k
            if(ktfrealq(ki))then
              isp=isp+1
              ktastk(isp)=ki
            elseif(ktflistq(ki))then
              kai=ktfaddr(ki)
              if(ktfoperq(klist(kai)))then
                call tfgetcanvasargstk(kai,irtc)
                if(irtc .ne. 0)then
                  return
                endif
              endif
            endif
          enddo
        endif
      elseif((kl%head%k .eq. ktfoper+mtfrule .or.
     $       kl%head%k .eq. ktfoper+mtfruledelayed) .and. 
     $       n .eq. 2)then
        k1=kl%dbody(1)
        if(ktfsymbolq(k1) .or. ktfoperq(k1))then
          k2=kl%dbody(2)
          if(ktfrealq(k2) .or. ktfstringq(k2) .or.
     $         ktflistq(k2,k2l) .and.
     $         k2l%head%k .eq. ktfoper+mtflist)then
            isp3=isp
            isp=isp+1
            ktastk(isp)=ktfsymbol+iflabel
            isp=isp+1
            dtastk(isp)=k1
            call tfdeval(isp-1,iflabel,kx,1,.false.,ev,irtc)
            if(irtc .ne. 0)then
              return
            endif
c            call tfdebugprint(kx,'TkOptionLabel',2)
            if(ktfstringq(kx,str))then
              isp=isp3+1
              call getstringbuf(strb,0,.true.)
              call putstringbufb(strb,str%str(2:str%nch-1),
     $             str%nch-2,full)
              dtastk(isp3+1)=kxstringbuftostring(strb)
              if(ktfrealq(k2,x))then
                k2%k=ktfstring+ktrsaloc(-1,x)
              elseif(ktflistq(k2,k2l))then
                call getstringbuf(strb,0,.true.)
                do i=1,k2l%nl
                  k2i=k2l%dbody(i)
                  if(ktfstringq(k2i,str))then
                    call tftclstringbuf(strb,str%str,str%nch,full)
                  else
                    call tfconvstrb(strb,k2i,nc,
     $                   .false.,.false.,-1,' ',irtc)
                  endif
                  call putstringbufb1(strb,' ')
                enddo
                k2=kxstringbuftostring(strb)
              endif
c              call tfdebugprint(k2,'== ',2)
              isp=isp3+2
              dtastk(isp)=k2
            endif
          endif
        endif
      endif
      irtc=0
      return
      end

      subroutine tfuniconv(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_string) , pointer :: str0
      integer*4 isp1,irtc,itfmessage,m,m1
      character*32768 str
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonstringq(ktastk(isp),str0))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      m=str0%nch
      call tfuniconvs(str0%str,m,str,m1)
      kx=kxsalocb(-1,str,m1)
      irtc=0
      return
      end

      logical*4 function tfxftq()
      use tfstk
      implicit none
      logical*4 xft
      integer*8 kaxfs,kxxft
      save kaxfs,kxxft,xft
      data kaxfs /0/
      if(kaxfs .eq. 0)then
        kaxfs=ktfsymbolz('$FontSystem',11)-4
        kxxft=ktfstring+ktsalocb(0,'Xft',3)
        xft=tfsameq(klist(kaxfs),kxxft)
      endif
      tfxftq=xft
      return
      end

      subroutine tfuniconvs(in,m,out,m1)
      implicit none
      integer*4 m,m1,i,j
      character*(m) in
      character*(*) out
      character(len=8) :: ubuf
      integer :: ulen
      logical*4 alt
      integer, external :: Unicode2UTF8
      integer :: kgcode(0:255)=(/
c     00h - 1Fh
     $     int2(z'0000'),int2(z'0001'),int2(z'0002'),int2(z'0003'),
     $     int2(z'0004'),int2(z'0005'),int2(z'0006'),int2(z'0007'),
     $     int2(z'0008'),int2(z'0009'),int2(z'000a'),int2(z'000b'),
     $     int2(z'000c'),int2(z'000d'),int2(z'000e'),int2(z'000f'),
     $     int2(z'0010'),int2(z'0011'),int2(z'0012'),int2(z'0013'),
     $     int2(z'0014'),int2(z'0015'),int2(z'0016'),int2(z'0017'),
     $     int2(z'0018'),int2(z'0019'),int2(z'001a'),int2(z'001b'),
     $     int2(z'001c'),int2(z'001d'),int2(z'001e'),int2(z'001f'),

c     20h - 3Fh
     $     int2(z'0020'),int2(z'0021'),int2(z'2200'),int2(z'0023'),
     $     int2(z'2203'),int2(z'0025'),int2(z'0026'),int2(z'220b'),
     $     int2(z'0028'),int2(z'0029'),int2(z'2217'),int2(z'002b'),
     $     int2(z'002c'),int2(z'2212'),int2(z'002e'),int2(z'002F'),
     $     int2(z'0030'),int2(z'0031'),int2(z'0032'),int2(z'0033'),
     $     int2(z'0034'),int2(z'0035'),int2(z'0036'),int2(z'0037'),
     $     int2(z'0038'),int2(z'0039'),int2(z'003a'),int2(z'003b'),
     $     int2(z'003c'),int2(z'003d'),int2(z'003e'),int2(z'003f'),

c     40h - 5Fh
     $     int2(z'2245'),int2(z'0391'),int2(z'0392'),int2(z'03a7'),
     $     int2(z'0394'),int2(z'0395'),int2(z'03a6'),int2(z'0393'),
     $     int2(z'0397'),int2(z'0399'),int2(z'03d1'),int2(z'039a'),
     $     int2(z'039b'),int2(z'039c'),int2(z'039d'),int2(z'039f'),
     $     int2(z'03a0'),int2(z'0398'),int2(z'03a1'),int2(z'03a3'),
     $     int2(z'03a4'),int2(z'03a5'),int2(z'03c2'),int2(z'03a9'),
     $     int2(z'039e'),int2(z'03a8'),int2(z'0396'),int2(z'005b'),
     $     int2(z'2234'),int2(z'005d'),int2(z'22a5'),int2(z'005f'),

c     60h - 7Fh
     $     int2(z'203e'),int2(z'03b1'),int2(z'03b2'),int2(z'03c7'),
     $     int2(z'03b4'),int2(z'03b5'),int2(z'03c6'),int2(z'03b3'),
     $     int2(z'03b7'),int2(z'03b9'),int2(z'03d5'),int2(z'03ba'),
     $     int2(z'03bb'),int2(z'03bc'),int2(z'03bd'),int2(z'03bf'),
     $     int2(z'03c0'),int2(z'03b8'),int2(z'03c1'),int2(z'03c3'),
     $     int2(z'03c4'),int2(z'03c5'),int2(z'03d6'),int2(z'03c9'),
     $     int2(z'03be'),int2(z'03c8'),int2(z'03b6'),int2(z'007b'),
     $     int2(z'007c'),int2(z'007d'),int2(z'223c'),int2(z'0020'),

c     80h - 9Fh
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),

c     A0h - BFh
     $     int2(z'0020'),int2(z'03d2'),int2(z'2032'),int2(z'2264'),
     $     int2(z'2044'),int2(z'221e'),int2(z'0192'),int2(z'2663'),
     $     int2(z'2666'),int2(z'2665'),int2(z'2660'),int2(z'2194'),
     $     int2(z'2190'),int2(z'2191'),int2(z'2192'),int2(z'2193'),
     $     int2(z'00b0'),int2(z'00b1'),int2(z'2033'),int2(z'2265'),
     $     int2(z'00d7'),int2(z'221d'),int2(z'2202'),int2(z'2022'),
     $     int2(z'00f7'),int2(z'2260'),int2(z'2261'),int2(z'2248'),
     $     int2(z'2026'),int2(z'0020'),int2(z'0020'),int2(z'21b5'),

c     C0h - DFh
     $     int2(z'2135'),int2(z'2111'),int2(z'211C'),int2(z'2118'),
     $     int2(z'2297'),int2(z'2295'),int2(z'2205'),int2(z'2229'),
     $     int2(z'222a'),int2(z'2283'),int2(z'2287'),int2(z'2284'),
     $     int2(z'2282'),int2(z'2286'),int2(z'2208'),int2(z'2209'),
     $     int2(z'2220'),int2(z'2207'),int2(z'00ae'),int2(z'00a9'),
     $     int2(z'2122'),int2(z'220f'),int2(z'221a'),int2(z'22c5'),
     $     int2(z'00ac'),int2(z'2227'),int2(z'2228'),int2(z'21d4'),
     $     int2(z'21d0'),int2(z'21d1'),int2(z'21d2'),int2(z'21d3'),

c     E0h -FFh
     $     int2(z'22c4'),int2(z'2329'),int2(z'00ae'),int2(z'00a9'),
     $     int2(z'2122'),int2(z'2211'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'232a'),int2(z'222b'),int2(z'2320'),
     $     int2(z'0020'),int2(z'2321'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020'),
     $     int2(z'0020'),int2(z'0020'),int2(z'0020'),int2(z'0020')/)
      alt=.false.
      i=1
      j=0
      do while(i .le. m)
        if(i .lt. m)then
          if(in(i:i+1) .eq. '`f')then
            alt=.true.
            j=j+2
            out(j-1:j)='`f'
            i=i+2
            cycle
          elseif(in(i:i+1) .eq. '`n')then
            alt=.false.
            j=j+2
            out(j-1:j)='`n'
            i=i+2
            cycle
          elseif(in(i:i+1) .eq. '``')then
          elseif(in(i:i) .eq. '`')then
            j=j+2
            out(j-1:j)=in(i:i+1)
            i=i+2
            cycle
          endif
        endif
        if(alt)then
          ulen=Unicode2UTF8(kgcode(ichar(in(i:i))),ubuf)
          out(j+1:j+ulen)=ubuf(1:ulen)
          j=j+ulen
          i=i+1
          cycle
        elseif(ichar(in(i:i)) .gt. 127)then
          ulen=Unicode2UTF8(ichar(in(i:i)),ubuf)
          out(j+1:j+ulen)=ubuf(1:ulen)
          j=j+ulen
          i=i+1
          cycle
        endif
        j=j+1
        out(j:j)=in(i:i)
        i=i+1
      enddo
      m1=j
      return
      end
