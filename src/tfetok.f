      module tftok
      implicit none
      character*24 oper1
      character*15 oper2
      character*4 oper3
      character*11 opern
      character*29, parameter ::
     $     oper=' +-*/=(){}[],;@#:^&<>|~.?"''\t\r'
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
        use convstr,only:Unicode2UTF8
        use ISO_C_BINDING
        implicit none
        type (sad_strbuf), pointer :: strb
        logical*4 full
        integer*4 ,intent(in):: l
        integer*4 ,intent(out):: ii
        integer*4 i,j,jmax
        integer :: ucode
        character(len=8) :: ubuf
        character*(*) ,intent(in):: string
        character ,intent(in):: del
        character ch
        call getstringbuf(strb,l,.true.)
        i=1
        do while(i <= l)
          ch=string(i:i)
          i=i+1
          if(ch == '\\' .and. i <= l)then
            ch=string(i:i)
            if('0' <= ch .and. ch <= '7')then
c     Octal character literal:		\#[##]
              ucode=iachar(ch)-iachar('0')
              do j=1,min(2,l-i)
                ch=string(i+j:i+j)
                if('0' <= ch .and. ch <= '7')then
                  ucode=ucode*8+(iachar(ch)-iachar('0'))
                else
                  exit
                endif
              enddo
              ch=char(ucode)
              i=i+j
            elseif(i+1 <= l .and. (ch == 'x' .or. ch == 'X'))then
c     Hexadecimal character literal:	\x#[#] or \X#[#]
              ucode=0
              do j=1,min(2,l-i)
                ch=string(i+j:i+j)
                if(    '0' <= ch .and. ch <= '9')then
                  ucode=ucode*16+(iachar(ch)-iachar('0'))
                elseif('a' <= ch .and. ch <= 'f')then
                  ucode=ucode*16+(iachar(ch)-iachar('a')+10)
                elseif('A' <= ch .and. ch <= 'F')then
                  ucode=ucode*16+(iachar(ch)-iachar('A')+10)
                else
                  exit
                endif
              enddo
              ch=char(ucode)
              i=i+j
            elseif(i+1 <= l .and. (ch == 'u' .or. ch == 'U'))then
c     Unicode character literal:	\u#[###] or \U#[#######]
              if(ch == 'U')then
                jmax=min(8,l-i)
              else
                jmax=min(4,l-i)
              endif
              ucode=0
              do j=1,jmax
                ch=string(i+j:i+j)
                if(    '0' <= ch .and. ch <= '9')then
                  ucode=ucode*16+(iachar(ch)-iachar('0'))
                elseif('a' <= ch .and. ch <= 'f')then
                  ucode=ucode*16+(iachar(ch)-iachar('a')+10)
                elseif('A' <= ch .and. ch <= 'F')then
                  ucode=ucode*16+(iachar(ch)-iachar('A')+10)
                else
                  exit
                endif
              enddo
              if(j > 1)then
                jmax=Unicode2UTF8(ucode,ubuf)
                if(jmax > 0)then
                  call putstringbufb(strb,ubuf,jmax,full)
                  i=i+j
                  cycle
                endif
              endif
              ch=string(i:i)
              i=i+1
            else
              select case (ch)
              case ('n')
                ch=C_NEW_LINE
                i=i+1
              case ('r')
                ch=C_CARRIAGE_RETURN
                i=i+1
              case ('t')
                ch=C_HORIZONTAL_TAB
                i=i+1
              case ('f')
                ch=C_FORM_FEED
                i=i+1
              case ('e')
                ch=char(27)
                i=i+1
              case (C_NEW_LINE)
                i=i+1
                cycle
              case default
                i=i+1
              end select
            endif
          elseif(notabspace(del) .and. ch == del)then
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
        character ,intent(in):: ch
        notabspace=ch /= ' ' .and. ch /= char(9)
        return
        end function

        subroutine tremovecont(string,buf)
        implicit none
        character*(*), intent(in)::string
        character*(*), intent(out)::buf
        integer*4 i,j
        i=0
        j=0
        do while(i < len(string))
          i=i+1
          if(i < len(string) .and. string(i:i+1) == '\\\n')then
            i=i+1
            cycle
          elseif(i < len(string)-1 .and.
     $           string(i:i+2) == '\\\r\n')then
            i=i+2
            cycle
          else
            j=j+1
            buf(j:j)=string(i:i)
          endif
        enddo
        return
        end subroutine

        real*8 function eval2(string,is1,istop,buf) result(vx)
        implicit none
        character*(*), intent(in) :: string
        character*(*), intent(out):: buf
        integer*4 , intent (in) :: is1
        integer*4 , intent (out) :: istop
        integer*4 l,m,j,nnl,is2
        real*8 eval1
        l=len(string)
        vx=eval1(string,l,is1,m)
        istop=is1+m
        if(m > 0 .or. string(is1:is1) == '.')then
          is2=is1+max(m,1)
          if(is2 < l .and. (string(is2:is2+1) == '\\\n' .or.
     $           index('eEdDxX',string(is2:is2)) /= 0))then
            call tedigicopy(string(is1:l),buf,j,nnl)
            vx=eval1(buf,j,1,m)
            istop=is1+m+nnl
          endif
        endif
        return
        end function

        subroutine tedigicopy(string,buf,j,nnl)
        implicit none
        character*(*), intent(in)::string
        character*(*), intent(out)::buf
        integer*4 , intent(out)::j,nnl
        integer*4 i
        character ch
        j=0
        i=0
        nnl=0
        do while (i < len(string))
          i=i+1
          ch=string(i:i)
          if(ch == '\\')then
            if(i < len(string) .and. string(i+1:i+1) == '\n')then
              nnl=nnl+2
              i=i+1
              cycle
            else
              exit
            endif
          elseif(ch >= '0' .and. ch <= '9'
     $           .or. index('eEdDxX-+.',ch) /= 0)then
            j=j+1
            buf(j:j)=ch
          else
            exit
          endif
        enddo
        return
        end subroutine

        subroutine tftokinit
        implicit none
        integer*4 i,ic1
        character ch
        do i=0,255
          ch=char(i)
          if(index(oper3,ch) /= 0)then
            levelop(i)=4
          elseif(index(oper2,ch) /= 0)then
            levelop(i)=3
          elseif(index(oper1,ch) /= 0)then
            levelop(i)=2
          elseif(index(oper//char(10),ch) /= 0)then
            levelop(i)=1
          endif
          if(ch >= 'a' .and. ch <= 'z')then
            ic1=100000+(i-ichar('a'))*2
          elseif(ch >= 'A' .and. ch <= 'Z')then
            ic1=100000+(i-ichar('A'))*2+1
          else
            ic1=i
          endif
          ichorder(i)=ic1
        enddo
        return
        end subroutine

      end module

      subroutine tfetok(string,istop,kx,icont,irt)
      use ISO_C_BINDING
      use tfstk
      use tftok
      use tfcode
      use strbuf
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klr
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer ::loc
      type (sad_strbuf), pointer :: strb
      integer*4 ,intent(out):: istop,irt
      integer*8 ,intent(in):: icont
      integer*4 nc,is2,isp0,lm,ich,nnl
      integer*8 kax,ka
      character*(*) ,intent(in):: string
      real*8 vx
      integer*4 l,i1,is1,notspace,ifromstr,
     $     i,jc,itfopcode,it1,it2,notany1
      character*1 ch
      character*2048 buf
      kx%k=ktfoper+mtfnull
      l=len(string)
      irt=-2
      if(l <= 0)then
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
      if(levelop(ich) /= 0)then
        select case (ch)
        case (C_NEW_LINE)
          if(is1 == 1 .or. string(is1-1:is1-1) /= '\\')then
            irt=-2
            istop=is1+1
            return
          endif
        case (C_CARRIAGE_RETURN)
          if(is1 == 1 .or. string(is1-1:is1-1) /= '\\')then
            irt=-2
            if(is1 < l .and.
     $           string(is1+1:is1+1) == C_NEW_LINE)then
              istop=is1+2
            else
              istop=is1+1
            endif
            return
          endif
        case ('''','"')
          isp0=isp
          kx=kxmakestring(string(is1+1:l),l-is1,ch,i1)
          call descr_strbuf(kx,strb)
          if(is1+i1 == l+1 .and. string(l:l) /= ch)then
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
        end select
        if(is1 == l)then
          lm=is1
        elseif(is1 == l-1)then
          if(levelop(ichar(string(l:l))) > 1)then
            lm=l
          else
            lm=is1
          endif
        elseif(levelop(ichar(string(is1+1:is1+1))) > 2)then
          if(levelop(ichar(string(is1+2:is1+2))) > 3)then
            lm=is1+2
          else
            lm=is1+1
          endif
        else
          lm=is1
        endif
        do i=lm,is1,-1
          jc=itfopcode(string(is1:i))
          if(jc >= 0)then
            irt=-1
            kax=jc
            istop=i+1
            select case (jc)
            case (mtfdot)
              vx=eval2(string,i,istop,buf)
              if(istop > i)then
                irt=0
                kx=dfromr(vx)
                return
              endif
              istop=i+1
            case (mtfreplace,mtfunset,mtfreplacerepeated,
     $             mtfincrement,mtfdecrement)
              if(i < l .and. string(istop:istop) <= '9'
     $             .and. string(istop:istop) >= '0')then
                istop=i
                select case (jc)
                case(mtfreplace)
                  kax=mtfdiv
                case(mtfreplacerepeated)
                  kax=mtfconcat
                case(mtfincrement)
                  kax=mtfplus
                case(mtfdecrement)
                  kax=mtfminus
                case default
                  kax=mtfset
                end select
              endif
            case(mtfatt)
              if(is1 > 1)then
                irt=-2
                istop=is1
                return
              endif
            case(mtfpower)
              if(is1 /= 1 .and. is1 < l-1
     $             .and. string(is1:is1+2) == '^^^')then
                irt=-2
                istop=is1
                return
              endif
            case(mtfleftcomment)
              istop=l+1
              if(is1 < l-2)then
                is2=index(string(is1+2:l),'*)')
                if(is2 > 0)then
                  is2=is2+is1+1
                  if(is2 < l-1)then
                    is1=notspace(string(1:l),is2+2)
                    if(is1 > 0)then
                      go to 1
                    endif
                  endif
                  irt=-2
                  return
                endif
              endif
              irt=-3
              return
            end select
            kx%k=ktfoper+kax
            return
          endif
        enddo
        istop=is1+1
      elseif(ch >= '0' .and. ch <= '9')then
        vx=eval2(string,is1,istop,buf)
        kx=dfromr(vx)
        irt=0
        return
      elseif(ch == '!')then
        return
      else
        irt=0
        nc=l-is1+1
        nnl=0
        i=is1+1
        if(levelop(ichar(string(i:i))) /= 0)then
          nc=1
          istop=i
        else
          do i=is1+2,l
            if(levelop(ichar(string(i:i))) /= 0)then
              if(string(i-1:i) == '\\\n')then
                nnl=nnl+2
                cycle
              elseif(string(i-1:i) == '\\\r')then
                nnl=nnl+merge(3,2,i < l .and. string(i+1:i+1) == '\n')
                cycle
              endif
              nc=i-is1
              istop=i
              exit
            endif
          enddo
        endif
        if(icont >= 0)then
          it1=is1
          it2=it1+nc-1
          if(nnl > 0)then
            call tremovecont(string(it1:it2),buf)
          endif
          if(index(string(it1:it2),'_') /= 0)then
            irt=0
            if(nnl == 0)then
              kx=kxpaloc(string(it1:it2))
            else
              kx=kxpaloc(buf(1:nc-nnl))
            endif
            return
          elseif(string(it1:it1) == '%')then
            if(nc == 1)then
              irt=0
              kx=kxadaloc(-1,1,klx)
              klx%head%k=ktfsymbol+ktfcopy1(iaxout)
              klx%dbody(1)%k=ktfsymbol+ktfcopy1(iaxline+4)
              go to 9000
            elseif(notany1(string(it1+1:it2),max(0,it2-it1),
     $             '0123456789',1) <= 0)then
              irt=0
              kx=kxavaloc(-1,1,klr)
              call descr_sad(kx,klx)
              klx%head%k=ktfsymbol+ktfcopy1(iaxout)
              if(nnl == 0)then
                klr%rbody(1)=ifromstr(string(it1+1:it2))
              else
                klr%rbody(1)=ifromstr(buf(1:nc-nnl))
              endif
              go to 9000
            endif
          endif
          if(nnl == 0)then
            kx=kxsymbolz(string(it1:it2),nc,symd)
          else
            kx=kxsymbolz(buf,nc-nnl,symd)
          endif
          if(rlist(iaximmediate) /= 0.d0)then
            call loc_namtbl(symd%sym%loc,loc)
            ka=loc%symdef
            if(ka /= 0)then
              call loc1_symdef(ka,symd)
              if(ktfconstantsymq(symd%sym))then
                kx=symd%value
              endif
            endif
          endif
        endif
      endif
 9000 return
      end

      integer function itfopcode(oper)
      use ophash
      implicit none
      integer*4 iop1,ih,j
      character*(*) ,intent(in):: oper
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
      if(j < 0)then
        go to 10
      elseif(iop1 == transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(2,ih)
      if(j < 0)then
        go to 10
      elseif(iop1 == transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(3,ih)
      if(j < 0)then
        go to 10
      elseif(iop1 == transfer(opcode(j),0))then
        itfopcode=j
        return
      endif
      j=iophash(4,ih)
      if(j >= 0)then
        if(iop1 == transfer(opcode(j),0))then
          itfopcode=j
          return
        endif
      endif
 10   if(oper1 == '=>  ')then
        itfopcode=mtfgeq
      elseif(oper1 == '=<  ')then
        itfopcode=mtfleq
      elseif(oper1 == '><  ')then
        itfopcode=mtfunequal
      else
        itfopcode=-1
      endif
      return
      end

      real*8 function fflogi(name,exist)
      use tfstk
      implicit none
      character*(*) ,intent(inout):: name
      logical*4 ,intent(out):: exist
      logical*4 tflogi
      call capita(name)
      fflogi=merge(1.d0,0.d0,tflogi(name,exist))
      return
      end

      subroutine tfflags(isp1,kx,irtc)
      use tfstk
      use ffs, only:fff,nflag,fname
      use tmacro
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: klx,klxi
      type (sad_string), pointer :: str
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: isp1
      integer*4 i,lenw,itfmessage
      if(isp > isp1+1)then
        irtc=itfmessage(9,'General::narg','"0"')
        return
      endif
      kx=kxadaloc(-1,nflag,klx)
      do i=1,nflag
        klx%dbody(i)=kxadaloc(0,2,klxi)
        klxi%dbody(1)=kxsalocb(0,fname(i),lenw(fname(i)),str)
        klxi%dbody(2)%k=merge(ktftrue,ktffalse,fff%flags(i))
      enddo
      irtc=0
      return
      end
      
      real*8 function ffval(name,exist)
      use tfstk
      use sad_main
      use ffs
      use ffs_pointer, only:errk,tfvalvar,pnamec
      implicit none
      logical*4 ,intent(out):: exist
      integer*4 i,lenw
      character*(*) ,intent(in):: name
      character*(MAXPNAME) name1
      ffval=0.d0
      if(name(1:1) /= '#')then
        exist=.false.
        return
      endif
      name1=name(2:lenw(name))
      do i=1,nele
        if(pnamec(nelvx(i)%klp) == name1)then
          exist=.true.
          ffval=tfvalvar(nelvx(i)%klp,
     $         nelvx(i)%ival)/errk(1,nelvx(i)%klp)
          return
        endif
      enddo
      exist=.false.
      return
      end

      subroutine tfsymbol(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_string), pointer :: str
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp /= isp1+1)then
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
