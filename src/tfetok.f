      module tftok
      implicit none
      character*29 oper
      character*24 oper1
      character*15 oper2
      character*4 oper3
      character*11 opern
      parameter (oper=' +-*/=(){}[],;@#:^&<>|~.?"'''//char(9))
      parameter (oper1='+-*/=(){}[],;@#:^&<>|~.?')
      parameter (oper2='=><&|/.+-:#@[*)')
      parameter (oper3='=>@.')
      parameter (opern='0123456789.')
      integer*4 levelop(0:255),ichorder(0:255)
      data levelop /256*0/

      contains
        type (sad_descriptor) function
     $       kxmakestring(string,l,del,ii)
        use tfstk
        use strbuf
        use ISO_C_BINDING
        implicit none
        type (sad_strbuf), pointer :: strb
        logical*4 full
        integer*4 l,i,j,jmax,ii
        integer :: ucode
        character(len=8) :: ubuf
        character*(*) string
        character ch,del
        integer(4), external :: Unicode2UTF8
        call getstringbuf(strb,l,.true.)
        i=1
        do while(i .le. l)
          ch=string(i:i)
          i=i+1
          if(ch .eq. '\\' .and. i .le. l)then
c '
            ch=string(i:i)
            if('0' .le. ch .and. ch .le. '7')then
c     Octal character literal:		\#[##]
              ucode=iachar(ch)-iachar('0')
              do j=1,min(2,l-i)
                ch=string(i+j:i+j)
                if('0' .le. ch .and. ch .le. '7')then
                  ucode=ucode*8+(iachar(ch)-iachar('0'))
                else
                  exit
                endif
              enddo
              ch=char(ucode)
              i=i+j
            elseif(i+1 .le. l .and. (ch .eq. 'x' .or. ch .eq. 'X'))then
c     Hexadecimal character literal:	\x#[#] or \X#[#]
              ucode=0
              do j=1,min(2,l-i)
                ch=string(i+j:i+j)
                if(    '0' .le. ch .and. ch .le. '9')then
                  ucode=ucode*16+(iachar(ch)-iachar('0'))
                elseif('a' .le. ch .and. ch .le. 'f')then
                  ucode=ucode*16+(iachar(ch)-iachar('a')+10)
                elseif('A' .le. ch .and. ch .le. 'F')then
                  ucode=ucode*16+(iachar(ch)-iachar('A')+10)
                else
                  exit
                endif
              enddo
              ch=char(ucode)
              i=i+j
            elseif(i+1 .le. l .and. (ch .eq. 'u' .or. ch .eq. 'U'))then
c     Unicode character literal:	\u#[###] or \U#[#######]
              if(ch .eq. 'U')then
                jmax=min(8,l-i)
              else
                jmax=min(4,l-i)
              endif
              ucode=0
              do j=1,jmax
                ch=string(i+j:i+j)
                if(    '0' .le. ch .and. ch .le. '9')then
                  ucode=ucode*16+(iachar(ch)-iachar('0'))
                elseif('a' .le. ch .and. ch .le. 'f')then
                  ucode=ucode*16+(iachar(ch)-iachar('a')+10)
                elseif('A' .le. ch .and. ch .le. 'F')then
                  ucode=ucode*16+(iachar(ch)-iachar('A')+10)
                else
                  exit
                endif
              enddo
              if(j .gt. 1)then
                jmax=Unicode2UTF8(ucode,ubuf)
                if(jmax .gt. 0)then
                  call putstringbufb(strb,ubuf,jmax,full)
                  i=i+j
                  cycle
                endif
              endif
              ch=string(i:i)
              i=i+1
            elseif(ch .eq. 'n')then
              ch=C_NEW_LINE
              i=i+1
            elseif(ch .eq. 'r')then
              ch=C_CARRIAGE_RETURN
              i=i+1
            elseif(ch .eq. 't')then
              ch=C_HORIZONTAL_TAB
              i=i+1
            elseif(ch .eq. 'f')then
              ch=C_FORM_FEED
              i=i+1
            elseif(ch .eq. 'e')then
              ch=char(27)
              i=i+1
            elseif(ch .eq. C_NEW_LINE)then
              i=i+1
              cycle
            else
              i=i+1
            endif
          elseif(notabspace(del) .and. ch .eq. del)then
            exit
          endif
          call putstringbufb1(strb,ch)
        enddo
        ii=i
        kxmakestring=sad_descr(strb%string(1))
        return
        end function 

        logical*4 function notabspace(ch)
        implicit none
        character ch
        notabspace=ch .ne. ' ' .and. ch .ne. char(9)
        return
        end function

      end module

      subroutine tfetok(string,istop,kx,icont,irt)
      use ISO_C_BINDING
      use tfstk
      use tftok
      use tfcode
      use strbuf
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klr
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer ::loc
      type (sad_strbuf), pointer :: strb
      integer*4 istop,nc,is2,isp0,lm,ich,irt
      integer*8 kax,ka,icont
      character*(*) string
      real*8 vx,eval1
      integer*4 m,l,i1,is1,notspace,ifromstr,
     $     i,jc,itfopcode,it1,it2,notany1
      character*1 ch
      kx%k=ktfoper+mtfnull
      l=len(string)
      irt=-2
      if(l .le. 0)then
        istop=2
        return
      endif
      istop=l+1
      if(notabspace(string(1:1)))then
        is1=1
        go to 1
      endif
      do i=2,l
        if(notabspace(string(i:i)))then
          is1=i
          go to 1
        endif
      enddo
      return
 1    ch=string(is1:is1)
      ich=ichar(ch)
      if(levelop(ich) .ne. 0)then
        if(ch .eq. C_NEW_LINE)then
          irt=-2
          istop=is1+1
          return
        elseif(ch .eq. '''' .or. ch .eq. '"')then
          isp0=isp
          kx=kxmakestring(string(is1+1:l),l-is1,ch,i1)
          call descr_strbuf(kx,strb)
          if(is1+i1 .eq. l+1 .and. string(l:l) .ne. ch)then
            istop=l
            irt=-4
            call tfreestringbuf(strb)
          else
            istop=is1+i1
            irt=0
            kx=kxstringbuftostring(strb)
          endif
          isp=isp0
          return
        endif
        if(is1 .eq. l)then
          lm=is1
        elseif(is1 .eq. l-1)then
          if(levelop(ichar(string(l:l))) .gt. 1)then
            lm=l
          else
            lm=is1
          endif
        elseif(levelop(ichar(string(is1+1:is1+1))) .gt. 2)then
          if(levelop(ichar(string(is1+2:is1+2))) .gt. 3)then
            lm=is1+2
          else
            lm=is1+1
          endif
        else
          lm=is1
        endif
        do i=lm,is1,-1
          jc=itfopcode(string(is1:i))
          if(jc .ge. 0)then
            irt=-1
            kax=jc
            istop=i+1
            if(jc .eq. mtfdot)then
              vx=eval1(string,l,is1,m)
              if(m .gt. 0)then
                istop=is1+m
                irt=0
                kx=dfromr(vx)
                return
              endif
            elseif(jc .eq. mtfreplace .or. jc .eq. mtfunset
     $           .or. jc .eq. mtfreplacerepeated
     $           .or. jc .eq. mtfincrement
     $           .or. jc .eq. mtfdecrement)then
              if(i .lt. l .and. string(istop:istop) .le. '9'
     $             .and. string(istop:istop) .ge. '0')then
                istop=i
                if(jc .eq. mtfreplace)then
                  kax=mtfdiv
                elseif(jc .eq. mtfreplacerepeated)then
                  kax=mtfconcat
                elseif(jc .eq. mtfincrement)then
                  kax=mtfplus
                elseif(jc .eq. mtfdecrement)then
                  kax=mtfminus
                else
                  kax=mtfset
                endif
              endif
            elseif(jc .eq. mtfatt)then
              if(is1 .gt. 1)then
                irt=-2
                istop=is1
                return
              endif
            elseif(jc .eq. mtfpower)then
              if(is1 .ne. 1 .and. is1 .lt. l-1
     $             .and. string(is1:is1+2) .eq. '^^^')then
                irt=-2
                istop=is1
                return
              endif
            elseif(jc .eq. mtfleftcomment)then
              istop=l+1
              if(is1 .lt. l-2)then
                is2=index(string(is1+2:l),'*)')
                if(is2 .gt. 0)then
                  is2=is2+is1+1
                  if(is2 .lt. l-1)then
                    is1=notspace(string(1:l),is2+2)
                    if(is1 .gt. 0)then
                      go to 1
                    endif
                  endif
                  irt=-2
                  return
                endif
              endif
              irt=-3
              return
            endif
            kx%k=ktfoper+kax
            return
          endif
        enddo
        istop=is1+1
      elseif(ch .ge. '0' .and. ch .le. '9')then
        vx=eval1(string,l,is1,m)
        if(m .gt. 0)then
          istop=is1+m
          irt=0
          kx=dfromr(vx)
        endif
        return
      elseif(ch .eq. '!')then
        return
      else
        irt=0
        nc=l-is1+1
        do i=is1+1,l
          if(levelop(ichar(string(i:i))) .ne. 0)then
            nc=i-is1
            istop=i
            exit
          endif
        enddo
        if(icont .ge. 0)then
          it1=is1
          it2=it1+nc-1
          if(index(string(it1:it2),'_') .ne. 0)then
            irt=0
            kx=kxpaloc(string(it1:it2))
            return
          elseif(string(it1:it1) .eq. '%')then
            if(nc .eq. 1)then
              irt=0
              kx=kxadaloc(-1,1,klx)
              klx%head%k=ktfsymbol+ktfcopy1(iaxout)
              klx%dbody(1)%k=ktfsymbol+ktfcopy1(iaxline+4)
              go to 9000
            elseif(notany1(string(it1+1:it2),max(0,it2-it1),
     $             '0123456789',1) .le. 0)then
              irt=0
              kx=kxavaloc(-1,1,klr)
              call descr_sad(kx,klx)
              klx%head%k=ktfsymbol+ktfcopy1(iaxout)
              klr%rbody(1)=ifromstr(string(it1+1:it2))
              go to 9000
            endif
          endif
          kx=kxsymbolz(string(it1:it2),nc,symd)
c          write(*,*)'tfetok ',string(it1:it2),kax,klist(kax)
c          call tfdebugprint(kx,'etok',1)
          if(rlist(iaximmediate) .ne. 0.d0)then
            call loc_namtbl(symd%sym%loc,loc)
            ka=loc%symdef
            if(ka .ne. 0)then
              call loc1_symdef(ka,symd)
              if(ktfconstantsymq(symd%sym))then
                kx=symd%value
              endif
            endif
          endif
c          call tfdebugprint(kx,'tfetok-symbol',1)
c          write(*,*)'kax ',kax
        endif
      endif
 9000 return
      end

      integer function itfopcode(oper)
      use ophash
      implicit none
      integer*4 iop1,ih,j
      character*(*) oper
      character*4 oper1
      if(opini)then
        call tfopcodehash
      endif
      oper1=oper
      iop1=transfer(oper1,0)
      ih=ichar(oper1(1:1))+ichar(oper1(2:2))
     $     +ichar(oper1(3:3))+ichar(oper1(4:4))
      ih=iand(ih,63)
      j=iophash(1,ih)
      if(j .lt. 0)then
        go to 10
      elseif(iop1 .eq. transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(2,ih)
      if(j .lt. 0)then
        go to 10
      elseif(iop1 .eq. transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(3,ih)
      if(j .lt. 0)then
        go to 10
      elseif(iop1 .eq. transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(4,ih)
      if(j .ge. 0)then
        if(iop1 .eq. transfer(opcode(j),0))then
          itfopcode=j
          return
        endif
      endif
 10   if(oper1 .eq. '=>  ')then
        itfopcode=mtfgeq
      elseif(oper1 .eq. '=<  ')then
        itfopcode=mtfleq
      elseif(oper1 .eq. '><  ')then
        itfopcode=mtfunequal
      else
        itfopcode=-1
      endif
      return
      end

      real*8 function fflogi(name,exist)
      use tfstk
      implicit none
      character*(*) name
      logical*4 exist,v,tflogi
      call capita(name)
      v=tflogi(name,exist)
      if(exist)then
        if(v)then
          fflogi=1.d0
        else
          fflogi=0.d0
        endif
      else
        fflogi=0.d0
      endif
      return
      end

      subroutine tfflags(isp1,kx,irtc)
      use tfstk
      use ffs, only:fff,nflag,fname
      use tmacro
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,klxi
      type (sad_string), pointer :: str
      integer*4 irtc,i,lenw,isp1,itfmessage
      if(isp .gt. isp1+1)then
        irtc=itfmessage(9,'General::narg','"0"')
        return
      endif
      kx=kxadaloc(-1,nflag,klx)
      do i=1,nflag
        klx%dbody(i)=kxadaloc(0,2,klxi)
        klxi%dbody(1)=kxsalocb(0,fname(i),lenw(fname(i)),str)
        if(fff%flags(i))then
          klxi%dbody(2)%k=ktftrue
        else
          klxi%rbody(2)=0.d0
        endif
      enddo
      irtc=0
      return
      end
      
      real*8 function ffval(name,exist)
      use tfstk
      use sad_main
      use ffs
      use ffs_pointer, only:klp,errk,tfvalvar,ival,pnamec
      implicit none
      logical*4 exist
      integer*4 i,lenw
      character*(*) name
      character*(MAXPNAME) name1
      ffval=0.d0
      if(name(1:1) .ne. '#')then
        exist=.false.
        return
      endif
      name1=name(2:lenw(name))
      do i=1,nele
        if(pnamec(klp(i)) .eq. name1)then
          exist=.true.
          ffval=tfvalvar(klp(i),ival(i))/errk(1,klp(i))
          return
        endif
      enddo
      exist=.false.
      return
      end

      subroutine tfsymbol(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. ktfstringq(dtastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      kx=kxsymbolf(str%str,str%nch,.false.)
      irtc=0
      return
      end

      subroutine tftokinit
      use tftok
      implicit none
      integer*4 i,ic1
      character ch
      do i=0,255
        ch=char(i)
        if(index(oper3,ch) .ne. 0)then
          levelop(i)=4
        elseif(index(oper2,ch) .ne. 0)then
          levelop(i)=3
        elseif(index(oper1,ch) .ne. 0)then
          levelop(i)=2
        elseif(index(oper//char(10),ch) .ne. 0)then
          levelop(i)=1
        endif
        if(ch .ge. 'a' .and. ch .le. 'z')then
          ic1=100000+(i-ichar('a'))*2
        elseif(ch .ge. 'A' .and. ch .le. 'Z')then
          ic1=100000+(i-ichar('A'))*2+1
        else
          ic1=i
        endif
        ichorder(i)=ic1
      enddo
      return
      end
