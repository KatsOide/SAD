      module convstr
      use tfstk

      contains
      subroutine tfstringfill(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_strbuf), pointer :: strb
      type (sad_string), pointer :: str1,str2
      integer*4 isp1,irtc,n1,n2,m,n,itfmessage,isp0,na
      logical*4 full
      if(isp1+3 /= isp)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(.not. ktfstringq(dtastk(isp1+1),str1) .or.
     $     .not. ktfstringq(dtastk(isp-1),str2) .or.
     $     .not. ktfrealq(dtastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"[String, String, Real]"')
        return
      endif
      irtc=0
      n=int(rtastk(isp))
      na=abs(n)
      n1=str1%nch
      if(n == 0)then
        kx=dxnulls
        return
      elseif(na == n1)then
        kx=dtastk(isp1+1)
        return
      endif
      isp0=isp
      m=min(n1,na)
      call getstringbuf(strb,na,.true.)
      n2=str2%nch
      if(n > 0)then
        call putstringbufb(strb,str1%str(1:1),m,full)
        do while(m < na)
          call putstringbufb(strb,str2%str,min(n2,na-m),full)
          m=m+n2
        enddo
      else
        if(m < na)then
          n=mod(na-m,n2)
          if(n > 0)then
            call putstringbufb(strb,str2%str(1+(n2-n):n2),n,full)
            na=na-n
          endif
          do while(m < na)
            call putstringbufb(strb,str2%str(1:n2),n2,full)
            na=na-n2
          enddo
        endif
        call putstringbufb(strb,str1%str(1+(n1-m):n1),m,full)
      endif
      kx=kxstringbuftostring(strb)
      isp=isp0
      return
      end

      subroutine tfstringtrim(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,i,i1,i2,nc,itfmessage
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. ktfstringq(dtastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Character-string"')
        return
      endif
      irtc=0
      nc=str%nch
      do i=1,nc
        if(str%str(i:i) /= ' ' .and. str%str(i:i) /= '\t')then
          i1=i
          go to 10
        endif
      enddo
      kx=dxnulls
      return
 10   do i=nc,1,-1
        if(str%str(i:i) /= ' ' .and. str%str(i:i) /= '\t')then
          i2=i
          go to 20
        endif
      enddo
      i2=nc
 20   kx=kxsalocb(-1,str%str(i1:i2),i2-i1+1)
      return
      end

      subroutine tfbaseform(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_strbuf), pointer :: strb
      integer*4 isp1,irtc,itfmessage,i,nc1,isp0
      real*8 base
      if(isp < isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(iand(ktrmask,ktastk(isp)) == ktfnr)then
        irtc=itfmessage(9,'General::wrongtype','"real number"')
        return
      endif
      base=aint(rtastk(isp))
      if(base < 2.d0 .or. base > 62.d0)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Base","between 2 and 62"')
        return
      endif
      isp0=isp
      irtc=0
      call getstringbuf(strb,0,.true.)
      do i=isp1+1,isp0-1
        if(iand(ktrmask,ktastk(i)) /= ktfnr)then
          call tfconvbase(strb,rtastk(i),base)
        else
          call tfconvstrb(strb,dtastk(i),
     $         nc1,.false.,.false.,-1,'*',irtc)
          if(irtc /= 0)then
            go to 10
          endif
        endif
      enddo
 10   kx=kxstringbuftostring(strb)
      isp=isp0
      return
      end

      subroutine tftoinputstring(isp1,kx,irtc)
      use tfcode
      use strbuf
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx,k
      type (sad_strbuf), pointer :: strb
      type (sad_namtbl), pointer :: loc
      type (sad_symbol), pointer :: sym
      integer*8 ka,ic,i
      integer*4 isp1,irtc,nc1,narg,isp0,itfmessage
      narg=isp-isp1
      if(narg /= 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      k=dtastk(isp)
      irtc=0
      if(ktfsymbolq(k,sym))then
        if(sym%gen <= 0
     $       .or. sym%gen == maxgeneration)then
          call sym_namtbl(sym,loc)
          ic=loc%cont
          do i=itfcontextpath,
     $         itfcontextpath+ilist(2,itfcontextpath-1)-1
            if(klist(i) == ic)then
              kx=loc%str%alloc
              return
            endif
          enddo
        endif  
      elseif(ktfoperq(k,ka))then
        call loc_namtbl(klist(ifunbase+ka),loc)
        kx=loc%str%alloc
        return
      endif
      isp0=isp
      call getstringbuf(strb,0,.true.)
      call tfconvstrb(strb,k,nc1,.true.,.false.,-1,'*',irtc)
      kx=kxstringbuftostring(strb)
      isp=isp0
      return
      end

      subroutine tffromcharactercode(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: kl
      type (sad_string), pointer :: str
      integer*4 isp1,narg,irtc,i,m,itfmessage
      narg=isp-isp1
      if(narg /= 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfrealq(ktastk(isp)))then
        kx=kxsalocb(-1,char(int(rtastk(isp))),1)
      elseif(tfreallistq(ktastk(isp),kl))then
        m=kl%nl
        kx=kxsalocbb(-1,m,str)
        do i=1,m
          str%str(i:i)=char(int(kl%rbody(i)))
        enddo
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals"')
        return
      endif
      irtc=0
      return
      end

      subroutine tftocharactercode(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      type (sad_rlist), pointer :: klr
      integer*4 isp1,narg,irtc,i,m,itfmessage
      narg=isp-isp1
      if(narg /= 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(ktfnonstringq(ktastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      m=str%nch
      kx=kxavaloc(-1,m,klr)
      if(m > 0)then
        do concurrent (i=1:m)
          klr%rbody(i)=iachar(str%str(i:i))
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfstringmatchq(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: stra,strp
      integer*4 isp1,irtc,nc,nc1,itfmessage
      logical*4 tmatchl
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(.not. ktfstringq(dtastk(isp1+1),stra)
     $     .or. .not. ktfstringq(dtastk(isp),strp))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"[String, Pattern-string]"')
        return
      endif
      nc=stra%nch
      irtc=0
      kx%k=0
      if(nc > 0)then
        nc1=strp%nch
        if(nc1 > 0)then
          if(tmatchl(stra%str,nc,strp%str,nc1))then
            kx%k=ktftrue
          endif
        endif
      endif
      return
      end

      subroutine tfstringposition(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx,ki
      type (sad_dlist), pointer ::kl,klx
      type (sad_rlist), pointer ::kli
      type (sad_string), pointer :: str,strp,stri
      integer*4 isp1,irtc,nc,nc1,i,j,
     $     i1,ip,isp0,l,itfmessage,indexb
      if(isp == isp1+1 .or. isp > isp1+3)then
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
      elseif(.not. ktfstringq(dtastk(isp1+1),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      if(isp == isp1+3)then
        if(.not. ktfrealq(dtastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real for #3"')
          return
        endif
        l=int(rtastk(isp))
      else
        l=str%nch
      endif
      nc=str%nch
      isp0=isp
      if(ktfstringq(dtastk(isp1+2), strp))then
        nc1=strp%nch
        if(nc1 > 0)then
          i1=1
          do while(i1 > 0 .and. i1 <= nc
     $         .and. isp < isp0+l)
            ip=indexb(str%str,nc,strp%str,nc1,i1)
            if(ip <= 0)then
              i1=0
            else
              isp=isp+1
              itastk(1,isp)=ip
              itastk(2,isp)=ip+nc1-1
              i1=ip+1
            endif
          enddo
        endif
      elseif(ktflistq(dtastk(isp1+2),kl))then
        LOOP_I: do i=1,kl%nl
          ki=kl%dbody(i)
          if(.not. ktfstringq(ki,stri))then
            isp=isp0
            return
          endif
          nc1=stri%nch
          if(nc1 > 0)then
            i1=1
            do while(i1 > 0 .and. i1 <= nc)
              if(isp >= isp0+l)then
                go to 100
              endif
              ip=indexb(str%str,nc,stri%str,nc1,i1)
              if(ip <= 0)then
                i1=0
              else
                if(i1 == 1)then
                  do j=isp0+1,isp
                    if(itastk(1,j) == ip
     $                   .and. itastk(2,j) == ip+nc1-1)then
                      cycle LOOP_I
                    endif
                  enddo
                endif
                isp=isp+1
                itastk(1,isp)=ip
                itastk(2,isp)=ip+nc1-1
                i1=ip+1
              endif
            enddo
          endif
        enddo LOOP_I
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"String or List of Strings for #2')
        return
      endif
 100  kx=kxadaloc(-1,isp-isp0,klx)
      do i=isp0+1,isp
        klx%dbody(i-isp0)=kxavaloc(0,2,kli)
        kli%rbody(1)=dble(itastk(1,i))
        kli%rbody(2)=dble(itastk(2,i))
      enddo
      isp=isp0
      irtc=0
      return
      end

      subroutine tftouppercase(isp1,kx,mode,irtc)
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str,strx
      integer*4 isp1,irtc,nc,mode,itfmessage,i,kk
      logical*4 rep
      integer*1 jdif
c      parameter (jdif=ichar('a')-ichar('A'))
      parameter (jdif=32)
c      include 'DEBUG.inc'
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(.not. ktfstringq(dtastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      nc=str%nch
      if(nc == 0)then
        kx=dxnulls
        irtc=0
        return
      elseif(nc == 1)then
        if(mode == 0)then
          if(str%str(1:1) >= 'a' .and.
     $         str%str(1:1) <= 'z')then
            kk=ichar(str%str(1:1))-jdif
            kx%k=ktfstring+iaxschar+kk*5+3
          else
            kx=dtastk(isp)
          endif
        else
          if(str%str(1:1) >= 'A' .and.
     $         str%str(1:1) <= 'Z')then
            kk=ichar(str%str(1:1))+jdif
            kx%k=ktfstring+iaxschar+kk*5+3
          else
            kx=dtastk(isp)
          endif
        endif
      else
        kx=kxsalocbb(-1,nc,strx)
        rep=.false.
        if(mode == 0)then
          do i=1,nc
            if(str%str(i:i) >= 'a' .and.
     $           str%str(i:i) <= 'z')then
              strx%str(i:i)=char(ichar(str%str(i:i))-jdif)
              rep=.true.
            else
              strx%str(i:i)=str%str(i:i)
            endif
          enddo
        else
          do i=1,nc
            if(str%str(i:i) >= 'A' .and.
     $           str%str(i:i) <= 'Z')then
              strx%str(i:i)=char(ichar(str%str(i:i))+jdif)
              rep=.true.
            else
              strx%str(i:i)=str%str(i:i)
            endif
          enddo
        endif
        if(.not. rep)then
          kx=dtastk(isp)
        endif
      endif
      irtc=0
      return
      end

      subroutine tfconvround(strb,x)
      use strbuf
      implicit none
      type (sad_strbuf), pointer :: strb
      integer*4 ir,l,ich
      real*8 x
      character*11 str
      logical*4 full
      if(strb%nch /= 0)then
        ich=ichar(strb%str(strb%nch:strb%nch))
        if(ich >= ichar('0') .and. ich <= ichar('9'))then
          call putstringbufb1(strb,' ')
        endif
      endif
      ir=int(sign(aint(abs(x)+0.5d0),x))
      call strfromil(ir,str,l)
      call putstringbufb(strb,str,l,full)
      return
      end

      integer*8 function ktrsaloc(mode,x)
      implicit none
      integer*4 l,lenw,mode
      real*8 x
      character*22 str1,autos1
      str1=autos1(x)
      l=lenw(str1)
      ktrsaloc=ktsalocb(mode,str1,l)
      return
      end

      subroutine tfdigitq(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,itfmessage,i
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(.not. ktfstringq(dtastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      irtc=0
      do i=1,str%nch
        if(str%istr(i) < ichar('0') .or.
     $       str%istr(i) > ichar('9'))then
          kx%k=0
          return
        endif
      enddo
      kx%k=ktftrue
      return
      end

      subroutine tfletterq(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 isp1,irtc,itfmessage,i
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(.not. ktfstringq(dtastk(isp),str))then
        irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        return
      endif
      irtc=0
      do i=1,str%nch
        if(str%istr(i) < ichar('A') .or.
     $       str%istr(i) > ichar('Z') .and.
     $       str%istr(i) < ichar('a') .or.
     $       str%istr(i) > ichar('z'))then
          kx%k=0
          return
        endif
      enddo
      kx%k=ktftrue
      return
      end

      integer function Unicode2UTF8(ucode, buffer)
      integer,          intent(in)  :: ucode
      character(len=*), intent(out) :: buffer
      if(ucode < 0)then
        Unicode2UTF8 = 0
      elseif(len(buffer) >= 1 .and. ucode <= ior(Z'0000007F',0))then
        buffer(1:1)  = char(            iand(Z'07F', ucode)         )
        Unicode2UTF8 = 1
      elseif(len(buffer) >= 2 .and. ucode <= ior(Z'000007FF',0))then
        buffer(1:1)  = char(ior(Z'0C0', iand(Z'01F', ucode / 2**6 )))
        buffer(2:2)  = char(ior(Z'080', iand(Z'03F', ucode        )))
        Unicode2UTF8 = 2
      elseif(len(buffer) >= 3 .and. ucode <= ior(Z'0000FFFF',0))then
        buffer(1:1)  = char(ior(Z'0E0', iand(Z'00F', ucode / 2**12)))
        buffer(2:2)  = char(ior(Z'080', iand(Z'03F', ucode / 2**6 )))
        buffer(3:3)  = char(ior(Z'080', iand(Z'03F', ucode        )))
        Unicode2UTF8 = 3
      elseif(len(buffer) >= 4 .and. ucode <= ior(Z'001FFFFF',0))then
        buffer(1:1)  = char(ior(Z'0F0', iand(Z'007', ucode / 2**18)))
        buffer(2:2)  = char(ior(Z'080', iand(Z'03F', ucode / 2**12)))
        buffer(3:3)  = char(ior(Z'080', iand(Z'03F', ucode / 2**6 )))
        buffer(4:4)  = char(ior(Z'080', iand(Z'03F', ucode        )))
        Unicode2UTF8 = 4
      elseif(len(buffer) >= 5 .and. ucode <= ior(Z'03FFFFFF',0))then
        buffer(1:1)  = char(ior(Z'0F8', iand(Z'003', ucode / 2**24)))
        buffer(2:2)  = char(ior(Z'080', iand(Z'03F', ucode / 2**18)))
        buffer(3:3)  = char(ior(Z'080', iand(Z'03F', ucode / 2**12)))
        buffer(4:4)  = char(ior(Z'080', iand(Z'03F', ucode / 2**6 )))
        buffer(5:5)  = char(ior(Z'080', iand(Z'03F', ucode        )))
        Unicode2UTF8 = 5
      elseif(len(buffer) >= 6 .and. ucode <= ior(Z'7FFFFFFF',0))then
        buffer(1:1)  = char(ior(Z'0FC', iand(Z'001', ucode / 2**30)))
        buffer(2:2)  = char(ior(Z'080', iand(Z'03F', ucode / 2**24)))
        buffer(3:3)  = char(ior(Z'080', iand(Z'03F', ucode / 2**18)))
        buffer(4:4)  = char(ior(Z'080', iand(Z'03F', ucode / 2**12)))
        buffer(5:5)  = char(ior(Z'080', iand(Z'03F', ucode / 2**6 )))
        buffer(6:6)  = char(ior(Z'080', iand(Z'03F', ucode        )))
        Unicode2UTF8 = 6
      else
        Unicode2UTF8 = 0
      endif
      return
      end

      subroutine tfstringjoin(isp1,kx,irtc)
      use strbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: ksi,ksx
      type (sad_strbuf), pointer :: strb
      integer*4 isp1,irtc,i,nc1,l,is,isp0
      if(isp1+1 == isp)then
        call tftostring(isp1,kx,.false.,irtc)
        return
      elseif(isp1 == isp)then
        kx=dxnulls
        irtc=0
        return
      endif
      l=0
      do i=isp1+1,isp
        if(ktfstringq(dtastk(i),ksi))then
          l=l+ksi%nch
        else
          l=0
          go to 1
        endif
      enddo
      if(l == 0)then
        kx=dxnulls
      else
        kx=kxscopy(ktfaddr(ktastk(isp1+1)),l,ksx)
        call loc_sad(ktfaddr(ktastk(isp1+1)),ksi)
        is=ksi%nch
        do i=isp1+2,isp
          call loc_sad(ktfaddr(ktastk(i)),ksi)
          if(ksi%nch /= 0)then
            ksx%str(is+1:is+ksi%nch)=ksi%str(1:ksi%nch)
            is=is+ksi%nch
          endif
        enddo
      endif
      irtc=0
      return
 1    isp0=isp
      call getstringbuf(strb,l,.true.)
      do i=isp1+1,isp0
        call tfconvstrb(strb,dtastk(i),
     $       nc1,.false.,.false.,-1,'*',irtc)
        if(irtc /= 0)then
          go to 9000
        endif
      enddo
      kx=kxstringbuftostring(strb)
      isp=isp0
      return
 9000 call tfreestringbuf(strb)
      isp=isp0
      return
      end subroutine

      end module

      character*(*) function tfconvstr(k,nc,form)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(out):: nc
      character*(*) ,intent(in):: form
      character*(len(tfconvstr)) buff
      call tfconvstrs(buff,k,nc,.false.,form)
      nc=min(nc,len(tfconvstr))
      tfconvstr=buff(1:nc)
      return
      end

      subroutine tfconvstrs(buff,k,nc,str,form)
      use tfstk
      use strbuf
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(out):: nc
      logical*4 ,intent(in):: str
      character*(*) ,intent(in):: form
      character*(*) ,intent(out):: buff
      call tfconvstrs0(buff,len(buff),k,nc,str,form)
      return
      end

      subroutine tfconvstrs0(buff,l,k,nc,str,form)
      use tfstk
      use strbuf
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_strbuf), pointer :: strb
      integer*4 ,intent(in):: l
      integer*4 ,intent(out):: nc
      character*(*) ,intent(in):: form
      character*(*) ,intent(out):: buff
      logical*4 ,intent(in):: str
      integer*4 irtc,isp0,lh,lb,lx
      isp0=isp
      lb=len(buff)
      call getstringbuf(strb,-l,.true.)
      call tfconvstrb(strb,k,nc,str,.false.,-1,form,irtc)
      if(irtc /= 0)then
        if(irtc > 0 .and. ierrorprint /= 0)then
          call tfreseterror
        endif
        nc=l+1
      endif
      if(nc > lb)then
        lh=lb/2+index(strb%str(lb/2+1:nc),',')
        buff(:lh)=strb%str(:lh)
        buff(lh+1:lh+5)=' ... '
        lx=max(index(strb%str(nc-(lb-lh-6):nc),','),1)
        lb=lb-lx+1
        buff(lh+6:lb)=strb%str(nc-(lb-lh-6):nc)
c        write(*,*)'convstrs ',l,lb,lh,nc,strb%str(lh-3:lh+3),' ',strb%str(nc-(lb-lh-6)-3:nc-(lb-lh-6)+3)
        nc=lb
      else
        buff(1:nc)=strb%str(1:nc)
      endif
      call tfreestringbuf(strb)
      isp=isp0
      return
      end

      subroutine tftostring(isp1,kx,symb,irtc)
      use tfstk
      use tfcode
      use strbuf
      use eeval
      use repl, only:tfgetoption
      use iso_c_binding
      implicit none
      type (sad_descriptor) k,kx,k1
      type (sad_namtbl), pointer :: loc
      type (sad_strbuf), pointer :: strb
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*8 ka,ic,icp
      integer*4 isp1,irtc,nc1,narg,isp0,itfmessage,i
      character*32 form
      logical*4 infm,symb,hold,gens
      type (sad_descriptor), save :: inputf,iholdf,igenf,istandf
      data inputf%k /0/
c      include 'DEBUG.inc'
      infm=.false.
      hold=symb
      gens=.false.
      narg=isp-isp1
      form='*'
      if(symb .and. narg /= 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      else
        if(narg > 1)then
          if(inputf%k == 0)then
            inputf =kxsymbolz('InputForm',9)
            iholdf =kxsymbolz('HoldForm',8)
            igenf  =kxsymbolz('GenericSymbolForm',17)
            istandf=kxsymbolz('StandardForm',12)
          endif
          do i=isp,isp1+2,-1
            if(ktflistq(dtastk(i)))then
              call tfgetoption('FormatType',dtastk(i),k1,irtc)
            else
              k1=dtastk(i)
            endif
            if(ktfsymbolq(k1))then
              if(tfsamesymbolq(k1,inputf))then
                infm=.true.
              elseif(tfsamesymbolq(k1,iholdf))then
                hold=.true.
              elseif(tfsamesymbolq(k1,igenf))then
                gens=.true.
              elseif(tfsamesymbolq(k1,istandf))then
                form=' '
              endif
            elseif(ktfstringq(k1,str))then
              form=str%str(1:str%nch)
            endif
          enddo
        endif
      endif
      if(hold)then
        k=dtastk(isp1+1)
        irtc=0
      else
        k=tfeevalref(dtastk(isp1+1),irtc)
        if(irtc /= 0)then
          return
        endif
      endif
      if(ktfsymbolq(k,sym))then
        if(gens .or. sym%gen <= 0
     $       .or. sym%gen == maxgeneration)then
          call sym_namtbl(sym,loc)
          ic=loc%cont
          do icp=itfcontextpath,
     $         itfcontextpath+ilist(2,itfcontextpath-1)-1
            if(klist(icp) == ic)then
              kx=sad_descr(loc%str)
              return
            endif
          enddo
        endif
      elseif(ktfoperq(k,ka))then
        call loc_namtbl(klist(klist(ifunbase+ka)),loc)
        kx=loc%str%alloc
        return
      elseif(symb)then
        irtc=itfmessage(9,'General::wrongtype','"Symbol"')
        return
      elseif(ktfstringq(k) .and. .not. infm)then
        kx=k
        return
      endif
      isp0=isp
      call getstringbuf(strb,0,.true.)
      call tfconvstrb(strb,k,nc1,infm,gens,-1,form,irtc)
      kx=kxstringbuftostring(strb)
      isp=isp0
      return
      end
