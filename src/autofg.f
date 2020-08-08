      character*(*) function autofg(x,form1)
      use tfcbk
      use tfstk, only:ktfenanq
      implicit none
      character*(*) form1
      character*16 form
      character*5 expstr
      character*32 buff,buf1,autos1
      character*32 zero
      parameter (zero='0000000000000000000000000000000')
      character*32 ovfl
      parameter (ovfl='$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
      logical*4 tzero,math,canv
      integer*4 idot,lf1,notspace,lenw,lc2,
     $     lf,lc,lc1,lm,i,li1,is,iexp,ifromstr,isign
      real*8 x,ax
      logical*4 shift
      if(form1 .eq. ' ' .or. form1 .eq. 'S')then
        autofg=autos1(x)
        return
      elseif(form1(1:1) .eq. 'F' .or. form1(1:1) .eq. 'G')then
        form='('//form1(:lenw(form1))//')'
        write(autofg,form)x
        if(autofg(1:1) .eq. '*')then
          autofg(1:len_trim(autofg))=ovfl
        endif
        return
      elseif(form1(1:1) .eq. 'E')then
        form='(1PD'//form1(2:lenw(form1))//')'
        write(autofg,form)x
        return
      endif
      buff=' '
      expstr=' '
      is=len(buf1)
      math=form1(1:1) .eq. 'M'
      canv=form1(1:1) .eq. 'C'
      tzero=form1(1:1) .eq. 'S' .or. math .or. canv
      if(tzero)then
        form=form1(2:)
      else
        form=form1
      endif
      if(form .eq. ' ')then
        tzero=.true.
        lc=18
        lf=15
      else
        idot=index(form,'.')
        lc=ifromstr(form(1:idot-1))
        lf=ifromstr(form(idot+1:))
      endif
      lc=min(18,lc)
      lf=min(15,lf)
      if(ktfenanq(x))then
        if(tzero)then
          autofg='NaN'
        else
          autofg(:lc-3)=' '
          autofg(lc-2:)='NaN'
        endif
        return
      elseif(x .eq. 0.d0)then
        if(tzero)then
          autofg='0'
        elseif(lf .eq. 0)then
          autofg(:lc-1)=' '
          autofg(lc:)='0'
        else
          autofg(1:lc-lf-2)=' '
          autofg(lc-lf-1:lc)=' .'//zero
          autofg(lc+1:)=' '
        endif
        return
      elseif(x .eq. dinfinity)then
        if(tzero)then
          autofg='INF'
        else
          autofg(:lc-3)=' '
          autofg(lc-2:)='INF'
        endif
        return
      elseif(x .eq. -dinfinity)then
        if(tzero)then
          autofg='-INF'
        else
          autofg(:lc-4)=' '
          autofg(lc-3:)='-INF'
        endif
        return
      endif
      ax=abs(x)
      call strfromd(ax,buff(2:19),isign,iexp)
      if(tzero)then
        if(x .ge. 0.d0)then
          lc1=lc
          lc2=lc1
        else
          lc1=lc-1
          lc2=lc1
        endif
      else
        if(x .ge. 0.d0)then
          lc1=lc-min(lf,1)
          lc2=lc-1
        else
          lc1=lc-1-min(lf,1)
          lc2=lc-2
        endif
      endif
 10   if(iexp .ge. 101)then
        lm=lc2-5
      elseif(iexp .ge. 11)then
        lm=lc2-4
      elseif(iexp .ge. 1)then
        lm=lc2-3
      elseif(iexp .ge. -8)then
        lm=lc2-4
      elseif(iexp .ge. -98)then
        lm=lc2-5
      else
        lm=lc2-6
      endif
      if(lm .lt. 0)then
        if(iexp .ge. 1)then
          buf1='***'
          is=3
        else
          buf1='0.'
          is=2
        endif
        go to 1000
      endif
      buff(1:1)='0'
      if(iexp .le. min(10,lc2) .and. iexp .gt. -min(3,lf))then
        call roundnumstr(buff(1:19),min(lc1,iexp+lf)+2,shift)
        if(shift)then
          iexp=iexp+1
          go to 10
        endif
        lf1=min(lf,lc1-1)
        li1=lc1-lf1-1
        is=lc1
        if(iexp .eq. lc1)then
          buf1=buff(2:lc1+1)
        elseif(iexp .ge. li1)then
          if(lf .eq. 0)then
            buf1=buff(2:iexp+1)
          else
            buf1=buff(2:iexp+1)//'.'//buff(iexp+2:lc1)
          endif
        elseif(iexp .ge. 1)then
          buf1(1:li1-iexp)=' '
          if(lf .eq. 0)then
            buf1(li1-iexp+1:)=buff(2:iexp+1)
          else
            buf1(li1-iexp+1:)=buff(2:iexp+1)//'.'//
     $           buff(iexp+2:min(lc1,iexp+lf1+1))
          endif
        else
          buf1(1:li1)=' '
          buf1(li1+1:li1+1)='.'
          if(iexp .eq. 0)then
            buf1(li1+2:)=buff(2:lf1+1)
          else
            buf1(li1+2:li1-iexp+1)=zero
            buf1(li1-iexp+2:)=buff(2:lf1+iexp+1)
          endif
        endif
      else
        call roundnumstr(buff(1:19),max(3,lm+2),shift)
        if(shift)then
          iexp=iexp+1
          go to 10
        endif
        call strfromi(iexp-1,expstr(2:))
        expstr(1:1)='E'
        if(lm .eq. 0)then
          buf1=buff(2:2)//expstr
          is=1
        elseif(lm .eq. 1)then
          buf1=buff(2:2)//'.'//expstr
          is=2
        else
          buf1=buff(2:2)//'.'//buff(3:lm+1)//expstr
          is=lm+1
        endif
      endif
 1000 if(tzero)then
        if(index(buf1,'.') .gt. 0)then
          do i=is,2,-1
            if(buf1(i:i) .eq. '.')then
              buf1(i:)=expstr
              go to 20
            elseif(buf1(i:i) .ne. ' ' .and.
     $             buf1(i:i) .ne. '0')then
              buf1(i+1:)=expstr
              go to 20
            endif
          enddo
        endif
      endif
 20   if(tzero)then
        if(x .ge. 0.d0)then
          autofg=buf1
        else
          i=notspace(buf1,1)
          autofg=' '//buf1
          autofg(i:i)='-'
        endif
      else
        if(x .ge. 0.d0)then
          autofg=' '//buf1
        else
          autofg='  '//buf1
          i=notspace(autofg,3)
          autofg(i-1:i-1)='-'
        endif
      endif
      if(math)then
        i=index(autofg,'E')
        if(i .gt. 0)then
          autofg=autofg(1:i-1)//' 10^'//autofg(i+1:)
        endif
      elseif(canv)then
        i=index(autofg,'E')
        if(i .gt. 0)then
          autofg=autofg(1:i-1)//'x10`u'//
     $         autofg(i+1:len_trim(autofg))//'`n'
        endif
      elseif(tzero)then
        do i=1,len(autofg)
          if(autofg(i:i) .ne. ' ')then
            autofg=autofg(i:)
            return
          endif
        enddo
      endif
      return
      end

      character*(*) function strfromis(n)
      implicit none
      integer*4 n,l
      call strfromil(n,strfromis,l)
      return
      end

      subroutine strfromi(n,string)
      implicit none
      integer*4 n,l
      character*(*) string
      call strfromil(n,string,l)
      return
      end

      subroutine strfromil(n,string,leng)
      implicit none
      integer*4 n,n1,n2,l,ich0,leng
c      parameter (ich0=ichar('0'))
      parameter (ich0=48)
      character*(*) string
      character*2,parameter :: ch(0:99)=[
     $     '00','01','02','03','04','05','06','07','08','09',
     $     '10','11','12','13','14','15','16','17','18','19',
     $     '20','21','22','23','24','25','26','27','28','29',
     $     '30','31','32','33','34','35','36','37','38','39',
     $     '40','41','42','43','44','45','46','47','48','49',
     $     '50','51','52','53','54','55','56','57','58','59',
     $     '60','61','62','63','64','65','66','67','68','69',
     $     '70','71','72','73','74','75','76','77','78','79',
     $     '80','81','82','83','84','85','86','87','88','89',
     $     '90','91','92','93','94','95','96','97','98','99']
      if(n .eq. 0)then
        string='0'
        leng=1
        return
      endif
      n1=abs(n)
      l=len(string)
      leng=l
      string=' '
      do while(n1 .ne. 0)
        n2=n1/100
        string(l-1:l)=ch(n1-n2*100)
        l=l-2
        n1=n2
      enddo
      if(string(l+1:l+1) .eq. '0')then
        l=l+1
      endif
      if(n .lt. 0)then
        string(l:l)='-'
        l=l-1
      endif
      if(l .ne. 0)then
        if(l .lt. 0)then
          string='***'
          leng=3
        else
          leng=leng-l
          string=string(l+1:)
        endif
      endif
      return
      end

      subroutine strfromd(x,string,isign,iexp)
      implicit none
      real*8 x
      real*8 xm(-9:9),xl,xh,a,a1,ai
      parameter (xl=1.d8,xh=1.d9)
      integer*4 mexp(-9:9),iexp,iae,i,iam,iaf,isign,ix
      character*(*) string
      character*32 zero
      parameter (zero='0000000000000000000000000000000')
      data xm/1.d-256,1.d-128,1.d-64,1.d-32,
     $        1.d-16, 1.d-8,  1.d-4, 1.d-2,
     $        1.d-1,  1.d0,   1.d1,  1.d2,
     $        1.d4,   1.d8,   1.d16, 1.d32,
     $        1.d64,  1.d128, 1.d256/
      data mexp/ -256,   -128,   -64,   -32,
     $           -16,    -8,     -4,    -2,
     $           -1,     0,      1,     2,
     $           4,      8,      16,    32,
     $           64,     128,    256/
      if(x .eq. 0.d0)then
        string=zero
        isign=0
        iexp=0
        return
      endif
      a=abs(x)
      if(a .lt. 2.d0**31)then
        ai=aint(x)
        if(ai .eq. x)then
          ix=int(x)
          call strfromifixed(ix,string(1:10))
          do i=1,10
            if(string(i:i) .ne. '0')then
              iexp=11-i
              string(1:iexp)=string(i:i+iexp-1)
              string(iexp+1:)=zero
              go to 1
            endif
          enddo
        endif
      endif
      iae=9
      if(a .ge. xh)then
        do i=-9,-1
          a1=a*xm(i)
          if(a1 .ge. xl)then
            a=a1
            iae=iae-mexp(i)
          endif
        enddo
        if(a .ge. xh)then
          a=a*.1d0
          iae=iae+1
        endif
        if(a .ge. xh)then
          string='INF'
          return
        endif
      elseif(a .lt. xl)then
        do i=9,1,-1
          a1=a*xm(i)
          if(a1 .lt. xh)then
            a=a1
            iae=iae-mexp(i)
          endif
        enddo
        if(a .lt. xl)then
          a=a*10.d0
          iae=iae-1
        endif
      endif
      iexp=iae
      iam=int(a)
      iaf=int((a-iam)*1.d9+.5d0)
      if(iaf .gt. 999999999)then
        iaf=0
        iam=iam+1
        if(iam .gt. 999999999)then
          iam=iam/10
          iexp=iexp+1
        endif
      endif
      call strfromifixed(iam,string(1:9))
      if(iaf .eq. 0)then
        string(10:18)=zero
      else
        call strfromifixed(iaf,string(10:18))
      endif
 1    if(x .lt. 0)then
        isign=-1
      else
        isign=0
      endif
      return
      end
      
      subroutine strfromifixed(n,string)
      implicit none
      integer*4 n,n1,n2,n3,i,icharzero
c      parameter (icharzero=ichar('0'))
      parameter (icharzero=48)
      character*(*) string
      character*32 zero
      parameter (zero='0000000000000000000000000000000')
      n1=n
      do i=len(string),1,-1
        if(n1 .eq. 0)then
          string(:i)=zero
          return
        endif
        n2=n1/10
        n3=n1-n2*10
        string(i:i)=char(icharzero+n3)
        n1=n2
      enddo
      return
      end

      subroutine roundnumstr(string,icol,shift)
      implicit none
      integer*4 icol,i,j,l
      character*(*) string
      character ch
      character*32 zero
      parameter (zero='0000000000000000000000000000000')
      logical*4 shift,inc
      l=len(string)
      shift=.false.
      if(icol .gt. l)then
        return
      endif
      ch=string(icol:icol)
      string(icol:l)=zero
      j=icol-1
      inc=ch .ge. '5'
      string(1:1)='0'
      do while(inc)
        string(j:j)=char(ichar(string(j:j))+1)
        inc=string(j:j) .gt. '9'
        if(inc)then
          string(j:j)='0'
          j=j-1
        endif
      enddo
      if(string(1:1) .eq. '1')then
        shift=.true.
        do i=icol,2,-1
          string(i:i)=string(i-1:i-1)
        enddo
        string(1:1)='0'
      endif
      return
      end

      integer function ifromstr(string)
      implicit none
      integer*4 n,is,i,notspace,i1,k
      character*(*) string
      i1=notspace(string,1)
      n=0
      is=1
      do i=i1,len(string)
        if(string(i:i) .eq. '-')then
          is=-1
        elseif(string(i:i) .eq. '+')then
        else
          k=ichar(string(i:i))-ichar('0')
          if(k .lt. 0 .or. k .gt. 9)then
            ifromstr=is*n
            return
          else
            n=n*10+k
          endif
        endif
      enddo
      ifromstr=is*n
      return
      end

      integer function ifromstrb(string,l)
      implicit none
      integer*4 n,is,i,notspace,i1,k,l
      character*(*) string
      i1=notspace(string(1:l),1)
      n=0
      is=1
      do i=i1,l
        if(string(i:i) .eq. '-')then
          is=-1
        elseif(string(i:i) .eq. '+')then
        else
          k=ichar(string(i:i))-ichar('0')
          if(k .lt. 0 .or. k .gt. 9)then
            ifromstrb=is*n
            return
          else
            n=n*10+k
          endif
        endif
      enddo
      ifromstrb=is*n
      return
      end

c
c A wrapper to amend a bug? of gdtoa, up to the version 20180730
c https://github.com/10110111/gdtoa-desktop
c
      character*(*) function autos1(x) result(s1)
      implicit none
      real*8 ,intent(in):: x
      integer*4 ic,i,l
      character*32 autos
      s1=autos(x)
      l=len_trim(s1)
      ic=index(s1(1:l),':')
      if(ic .le. 0)then
        return
      endif
c      write(*,*)'autos1 ',s1(1:l)
      s1(ic:l)=s1(ic+1:l)
      l=l-1
      do i=ic-1,1,-1
        if(s1(i:i) .eq. '.')then
          cycle
        elseif(s1(i:i) .eq. '-')then
          s1='-1'//s1(i+1:l)
          return
        elseif(ichar(s1(i:i)) .le. ichar('8'))then
          s1(i:i)=char(ichar(s1(i:i))+1)
          return
        else
          s1(i:i)='0'
          cycle
        endif
      enddo
      s1='1'//s1(1:l)
      return
      end
