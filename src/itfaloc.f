      integer*8 function ktfaloc(mode,ktype,nw)
      use tfstk
      implicit none
      integer*4 nw,mode
      integer*8 ktype,k1,itfroot
      k1=ktaloc(nw+2)
      ilist(2,k1-1)=0
      if(mode .eq. 0)then
        klist(k1)=ktype
      else
        itfroot=itflocal+levele
        klist(k1)=ktype+klist(itfroot)
        klist(itfroot)=k1
      endif
      ilist(1,k1+1)=mode+1
      ktfaloc=k1+2
      return
      end

      integer*8 function ktfalocr(mode,ktype,nw)
      use tfstk
      implicit none
      integer*4 nw,mode
      integer*8 ktype,ktalocr,k1,itfroot
      k1=ktalocr(nw+2)
      ilist(2,k1-1)=0
      if(mode .eq. 0)then
        klist(k1)=ktype
      else
        itfroot=itflocal+levele
        klist(k1)=ktype+klist(itfroot)
        klist(itfroot)=k1
      endif
      ilist(1,k1+1)=mode+1
      ktfalocr=k1+2
      return
      end

      integer*8 function ktlookupc(name,nc,icont,cre)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_namtbl), pointer :: loc
      integer*4 nc
      logical*4 cre
      character*(nc) name
      integer*8 i,i1,i0,itfroot,ktalocr,icont
      integer*4 ithash,nw
      itfroot=icont+ithash(name,nc)+1
      i=klist(itfroot)
      i1=i
      i0=i
      do while(i .ne. itfroot)
        call loc_namtbl(i,loc)
        if(nc .eq. loc%str%nch)then
          if(loc%str%str(1:nc) .eq. name)then
            if(i .ne. i1)then
              klist(i0)=loc%next
              loc%next=klist(itfroot)
              klist(itfroot)=i
            endif
            ktlookupc=i
            return
          endif
        endif
        i0=i
        i=loc%next
      enddo
      if(cre)then
        nw=nc/8+8
        i=ktalocr(nw)
        call loc_namtbl(i,loc)
        loc%next=klist(itfroot)
        klist(itfroot)=i
        loc%symdef=0
        loc%cont=icont
        loc%str%len=nc/8+5
        loc%str%alloc=sad_descr(loc%str)
        loc%str%ref=1
        loc%str%gen=0
        loc%str%nch=nc
        loc%str%nc=0
        loc%str%str(1:nc)=name
        ktlookupc=i
      else
        ktlookupc=0
      endif
      return
      end

      integer*8 function ktalocr(n)
      use tfstk
      implicit none
      integer*8 kresv
      integer*4 nresv,n,nsize
      parameter (nsize=2**18-1)
      save kresv
      data kresv,nresv /0,0/
      if(n .gt. nsize-5)then
        ktalocr=ktaloc(n)
        return
      endif
 1    if(nresv .eq. 0)then
        kresv=ktaloc(nsize)
        nresv=nsize
      endif
      if(n .gt. nresv-5)then
        ilist(1,kresv-1)=nresv+1
        call tfree(kresv)
        nresv=0
        go to 1
      endif
      ktalocr=kresv
      ilist(1,kresv-1)=n+1
      kresv=kresv+n+1
      nresv=nresv-n-1
      return
      end

      integer*8 function ktfsymbolf(name,l,const)
      use tfstk
      implicit none
      integer*4 l
      character name(l)
      logical*4 const
      ktfsymbolf=ktfaddr(kxsymbolf(name,l,const))
      return
      end

      recursive integer*8 function ktfsymbolc(name,l,icont)
     $     result(kres)
      use tfstk
      implicit none
      type (sad_symdef), pointer :: contd
      character*(*) name
      integer*8 ktsydefc,k1,ktcontaloc,icont,ic,ic1,j
      integer*4 l,i
      i=index(name(1:l), '`')
      if(i .le. 0)then
        if(icont .eq. 0)then
          k1=ktsydefc(name,l,itfcontext,.false.)
          if(k1 .ne. 0)then
            kres=k1
            return
          else
            do j=itfcontextpath,itfcontextpath
     $           +ilist(2,itfcontextpath-1)-1
              if(klist(j) .ne. itfcontext)then
                k1=ktsydefc(name,l,klist(j),.false.)
                if(k1 .ne. 0)then
                  kres=k1
                  return
                endif
              endif
            enddo
            kres=ktsydefc(name,l,itfcontext,.true.)
          endif
        else
          kres=ktsydefc(name,l,icont,.true.)
        endif
      elseif(i .eq. 1)then
        if(l .eq. 1)then
          kres=ktsydefc('`',1,itfcontroot,.true.)
        else
          kres=ktfsymbolc(name(2:l),l-1,itfcontroot)
        endif
      else
        if(icont .eq. 0)then
          ic=ktsydefc(name(1:i),i,itfcontroot,.true.)
        else
          ic=ktsydefc(name(1:i),i,icont,.true.)
        endif
        call loc_sad(ic,contd)
        ic1=contd%value%k
        if(contd%sym%gen .ne. -3)then
          call tflocal(ic1)
          ic1=ktcontaloc(ic)
        else
          ic1=ktfaddr(ic1)
        endif
        if(i .eq. l)then
          kres=ic
        else
          kres=ktfsymbolc(name(i+1:l),l-i,ic1)
        endif
      endif
      return
      end

      integer*8 function ktcontaloc(ic)
      use tfstk
      implicit none
      type (sad_symdef), pointer :: contd
      integer*8 itf,ic,i,ktavalocr
      ktcontaloc=ktavalocr(0,nsymhash+1)
      if(ic .ne. 0)then
        call loc_symdef(ic,contd)
        contd%value%k=ktflist+ktcontaloc
        contd%sym%gen=-3
        klist(ktcontaloc)=ktfsymbol+ic
      else
        klist(ktcontaloc)=0
      endif
      ilist(2,ktcontaloc-1)=0
      do i=0,nsymhash
        itf=ktcontaloc+i+1
        klist(itf)=itf
      enddo
      return
      end

      integer*8 function ktsydefc(string,leng,icont,cre)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      character*(*) string
      integer*8 icont,loc,ktlookupc
      integer*4 leng
      logical*4 cre
      integer*4 lg,ls,i,j,is
      lg=0
      ls=leng
      is=index(string(:ls),'$')
      if(is .gt. 0)then
        do1: do i=leng-1,is,-1
          if(string(i:i) .eq. '$')then
            do j=i+1,leng
              if(string(j:j) .lt. '0' .or.
     $           string(j:j) .gt. '9')then
                lg=0
                exit do1
              endif
              lg=lg*10+ichar(string(j:j))-ichar('0')
            enddo
            ls=i-1
            exit
          endif
        enddo do1
      endif
      loc=ktlookupc(string,ls,icont,cre)
      if(loc .ne. 0)then
        call tfsydefg(loc,kx,lg)
        ktsydefc=ktfaddrd(kx)
      else
        ktsydefc=0
      endif
      return
      end

      subroutine tfsydefg(locp,kx,ig0)
      use tfstk
      use tfcode
      use mackw
      implicit none
      type (sad_descriptor) kx
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      integer*8 kas1,locp
      integer*4 nc,idx
      character*(MAXPNAME) name
      integer*4 ig,ig1,ig0
      integer*4 hsrchz1
      call loc_namtbl(locp,loc)
      kas1=loc%symdef
      ig=max(0,ig0)
      do while(kas1 .ne. 0)
        ig1=max(0,ilist(2,kas1+7))
c        call tfdebugprint(ktfsymbol+kas1+8,'sydefg',1)
c        write(*,*)ig,ig1,kas1,klist(loc+1)
        if(ig .eq. ig1)then
          kx%k=ktfsymbol+kas1+8
          return
        elseif(ig .lt. ig1)then
          kas1=klist(kas1)
        else
          kx=kxnaloc1(ig,locp)
          return
        endif
      enddo
      if(kas1 .le. 0)then
        nc=loc%str%nch
        if(nc .le. MAXPNAME)then
          if(nc .le. 0)then
            write(*,*)'tfsydef ',loc
            call abort
          else
            name(1:nc)=loc%str%str(1:nc)
c            write(*,*)'tfsydefg ',name(1:nc)
            idx=hsrchz1(name(1:nc))
c            write(*,*)' : ',idx
            if(idx .ne. 0 .and. idtype(idx) .eq. icGLR)then
              kx=kxnaloc1(ig,locp)
              call descr_sad(kx,symd)
              symd%value%k=ktfref+idval(idx)
              return
            endif
          endif
        endif
        kx=kxnaloc1(ig,locp)
      else
        kx%k=ktfsymbol+kas1+8
      endif
      return
      end

      integer*4 function ithash(name,nc)
      use tfstk
      implicit none
      integer*4 nc,i,ih
      character name(nc)
      ih=ichar(name(1))
      if(nc .eq. 2)then
        ih=ih+ichar(name(2))
      elseif(nc .eq. 3)then
        ih=ih+ichar(name(2))+ichar(name(3))
      elseif(nc .eq. 4)then
        ih=ih+ichar(name(2))+ichar(name(3))+ichar(name(4))
      elseif(nc .eq. 5)then
        ih=ih+ichar(name(2))+ichar(name(3))+ichar(name(4))
     $       +ichar(name(5))
      elseif(nc .eq. 6)then
        ih=ih+ichar(name(2))+ichar(name(3))+ichar(name(4))
     $       +ichar(name(5))+ichar(name(6))
      else
        do i=2,nc
          ih=ih+ichar(name(i))
        enddo
      endif
      ithash=iand(ih,nsymhash)
      return
      end

      integer*8 function ktavalocr(mode,nd)
      use tfstk
      implicit none
      integer*4 mode,nd
      integer*8 k1,ktalocr,itfroot
      k1=ktalocr(nd+3)
      ilist(1,k1-1)=0
      ilist(2,k1-1)=kconstarg
      if(mode .eq. 0)then
        klist(k1)=ktflist
      else
        itfroot=itflocal+levele
        klist(k1)=ktflist+klist(itfroot)
        klist(itfroot)=k1
      endif
      ilist(1,k1+1)=mode+1
      ilist(2,k1+1)=nd
      klist(k1+2)=ktfoper+mtflist
      ktavalocr=k1+2
      return
      end

      integer*8 function ktcalocm(n)
      use tfstk
      implicit none
      type(sad_dlist), pointer :: kl
      integer*4 n,i
      integer*8 k
      k=ktaloc(n*6-1)+2
      do i=1,n
        call loc_sad(k+(i-1)*6,kl)
        kl%lenp=int2(0)
        kl%lena=int2(0)
        kl%attr=0
        kl%alloc%k=ktflist
        kl%ref=1
        kl%nl=2
        kl%head%k=ktfoper+mtfcomplex
      enddo
      ktcalocm=k
      return
      end
      
      integer*8 function ktrvaloc(name,x)
      use tfstk
      implicit none
      integer*4 lenw
      real*8 x
      character*(*) name
      ktrvaloc=ktfsymbolz(name,lenw(name))-4
      call tflocal(klist(ktrvaloc))
      rlist(ktrvaloc)=x
      return
      end
      
      integer*8 function ktcvaloc(name,x,y)
      use tfstk
      implicit none
      integer*4 lenw
      real*8 x,y
      character*(*) name
      ktcvaloc=ktfsymbolz(name,lenw(name))-4
      call tflocal(klist(ktcvaloc))
      dlist(ktcvaloc)=kxcalocv(0,x,y)
      return
      end
      
      subroutine tsvaloc(name,str)
      use tfstk
      implicit none
      integer*8 kax,kas
      integer*4 lenw,loc
      character*(*) name,str
      loc=-1
      kax=ktfsymbolz(name,lenw(name))-4
      call tflocal(klist(kax))
      kas=ktsalocb(0,str,len_trim(str))
      klist(kax)=ktfstring+kas
      return
      end
      
      integer*8 function ktsalocbi(mode,string,i,leng)
      use tfstk
      implicit none
      integer*4 mode,i,leng
      character string(i+leng-1)
      ktsalocbi=ktsalocb(mode,string(i),leng)
      return
      end

      subroutine tfpadstr(string,kp,leng)
      use tfstk
      implicit none
      integer*8 kp
      integer*4 leng
      character string(leng)
      klist(leng/8+kp)=0
c     Terminate string buffer by NULL character
      call tmovb(string,jlist(1,kp),leng)
      return
      end

      integer*8 function ktaaloc(mode,n)
      use tfstk, kf => ktaaloc
      implicit none
      integer*4 mode,n
      ktaaloc=kf(mode,n)
      return
      end

      integer*8 function ktadaloc(mode,n)
      use tfstk, kf => ktadaloc
      implicit none
      integer*4 mode,n
      ktadaloc=kf(mode,n)
      return
      end

      integer*8 function ktavaloc(mode,n)
      use tfstk, kf => ktavaloc
      implicit none
      integer*4 mode,n
      ktavaloc=kf(mode,n)
      return
      end

      integer*8 function ktfsymbolz(name,l)
      use tfstk, kf => ktfsymbolz
      implicit none
      integer*4 l
      character name(l)
      ktfsymbolz=kf(name,l)
      return
      end

      integer*8 function ktsalocb(mode,str,leng)
      use tfstk, kf => ktsalocb
      implicit none
      integer*4 mode,leng
      character str(leng)
      ktsalocb=kf(mode,str,leng)
      return
      end
