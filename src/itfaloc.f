      integer*8 function ktfaloc(mode,ktype,nw)
      use tfstk
      implicit none
      integer*4 ,intent(in):: nw,mode
      integer*8 ,intent(in):: ktype
      integer*8 k1,itfroot
      k1=ktaloc(nw+2)
      ilist(2,k1-1)=0
      if(mode == 0)then
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
      integer*4 ,intent(in):: nw,mode
      integer*8 ,intent(in):: ktype
      integer*8 ktalocr,k1,itfroot
      k1=ktalocr(nw+2)
      ilist(2,k1-1)=0
      if(mode == 0)then
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
      integer*4 ,intent(in):: nc
      logical*4 ,intent(in):: cre
      character*(nc),intent(in):: name
      integer*8 ,intent(in):: icont
      integer*8 i,i1,i0,itfroot,ktalocr
      integer*4 ithash,nw
      itfroot=icont+ithash(name,nc)+1
      i=klist(itfroot)
      i1=i
      i0=i
      do while(i /= itfroot)
        call loc_namtbl(i,loc)
        if(nc == loc%str%nch)then
          if(loc%str%str(1:nc) == name)then
            if(i /= i1)then
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
      integer*4 ,intent(in):: n
      integer*4 nresv,nsize
      parameter (nsize=2**18-1)
      save kresv
      data kresv,nresv /0,0/
      if(n > nsize-5)then
        ktalocr=ktaloc(n)
        return
      endif
 1    if(nresv == 0)then
        kresv=ktaloc(nsize)
        nresv=nsize
      endif
      if(n > nresv-5)then
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
      integer*4 ,intent(in):: l
      character ,intent(in):: name(l)
      logical*4 ,intent(in):: const
      ktfsymbolf=ktfaddr(kxsymbolf(name,l,const))
      return
      end

      recursive integer*8 function ktfsymbolc(name,l,icont)
     $     result(kres)
      use tfstk
      implicit none
      type (sad_symdef), pointer :: contd
      character*(*) ,intent(in):: name
      integer*8 ktsydefc,k1,ktcontaloc,icont,ic,ic1,j
      integer*4 ,intent(in):: l
      integer*4 i
      i=index(name(1:l), '`')
      if(i <= 0)then
        if(icont == 0)then
          k1=ktsydefc(name,l,itfcontext,.false.)
          if(k1 /= 0)then
            kres=k1
            return
          else
            do j=itfcontextpath,itfcontextpath
     $           +ilist(2,itfcontextpath-1)-1
              if(klist(j) /= itfcontext)then
                k1=ktsydefc(name,l,klist(j),.false.)
                if(k1 /= 0)then
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
      elseif(i == 1)then
        if(l == 1)then
          kres=ktsydefc('`',1,itfcontroot,.true.)
        else
          kres=ktfsymbolc(name(2:l),l-1,itfcontroot)
        endif
      else
        ic=ktsydefc(name(1:i),i,
     $       merge(itfcontroot,icont,icont == 0),.true.)
        call loc_sad(ic,contd)
        ic1=contd%value%k
        if(contd%sym%gen /= -3)then
          call tflocal(ic1)
          ic1=ktcontaloc(ic)
        else
          ic1=ktfaddr(ic1)
        endif
        if(i == l)then
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
      integer*8 ,intent(in):: ic
      integer*8 itf,i,ktavalocr
      ktcontaloc=ktavalocr(0,nsymhash+1)
      if(ic /= 0)then
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
      character*(*) ,intent(in):: string
      integer*8 ,intent(in):: icont
      integer*8 loc,ktlookupc
      integer*4 ,intent(in):: leng
      logical*4 ,intent(in):: cre
      integer*4 lg,ls,i,j,is
      lg=0
      ls=leng
      is=index(string(:ls),'$')
      if(is > 0)then
        do1: do i=leng-1,is,-1
          if(string(i:i) == '$')then
            do j=i+1,leng
              if(string(j:j) < '0' .or.
     $           string(j:j) > '9')then
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
      if(loc /= 0)then
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
      use tfmessage,only:tfnewsym
      use mackw
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_symdef), pointer :: symd
      type (sad_namtbl), pointer :: loc
      integer*8 ,intent(in):: locp
      integer*8 kas1
      integer*4 nc,idx
      character*(MAXPNAME) name
      integer*4 ,intent(in):: ig0
      integer*4 ig,ig1,itfmessagestr,irtc
      integer*4 hsrchz1
      call loc_namtbl(locp,loc)
      kas1=loc%symdef
      ig=max(0,ig0)
      do while(kas1 /= 0)
        ig1=max(0,ilist(2,kas1+7))
c        call tfdebugprint(dfromk(ktfsymbol+kas1+8),'tfsydefg',1)
c        write(*,*)ig,ig1,kas1
        if(ig == ig1)then
          kx%k=ktfsymbol+kas1+8
          return
        elseif(ig < ig1)then
          kas1=klist(kas1)
        else
          kx=kxnaloc1(ig,locp)
          return
        endif
      enddo
      if(kas1 <= 0)then
        nc=loc%str%nch
        if(nc <= MAXPNAME)then
          if(nc <= 0)then
            write(*,*)'tfsydef ',loc
            call abort
          else
            name(1:nc)=loc%str%str(1:nc)
            idx=hsrchz1(name(1:nc))
            if(idx /= 0 .and. idtype(idx) == icGLR)then
              kx=kxnaloc1(ig,locp)
              call descr_sad(kx,symd)
              symd%value%k=ktfref+idval(idx)
              return
            endif
          endif
        endif
        kx=kxnaloc1(ig,locp)
        if(tfnewsym(.false.))then
          irtc=itfmessagestr(9,'General::newsym',loc%str%str(1:loc%str%nch))
          call tferrorhandle(dlist(loc%cont),irtc)
        endif
      else
        kx%k=ktfsymbol+kas1+8
      endif
      return
      end

      integer*4 function ithash(name,nc)
      use tfstk
      implicit none
      integer*4 ,intent(in):: nc
      integer*4 i,ih
      character ,intent(in):: name(nc)
      ih=ichar(name(1))
      if(nc == 2)then
        ih=ih+ichar(name(2))
      elseif(nc == 3)then
        ih=ih+ichar(name(2))+ichar(name(3))
      elseif(nc == 4)then
        ih=ih+ichar(name(2))+ichar(name(3))+ichar(name(4))
      elseif(nc == 5)then
        ih=ih+ichar(name(2))+ichar(name(3))+ichar(name(4))
     $       +ichar(name(5))
      elseif(nc == 6)then
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
      integer*4 ,intent(in):: mode,nd
      integer*8 k1,ktalocr,itfroot
      k1=ktalocr(nd+3)
      ilist(1,k1-1)=0
      ilist(2,k1-1)=kconstarg
      if(mode == 0)then
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
      
      subroutine tsvaloc(name,str)
      use tfstk
      implicit none
      integer*8 kax,kas
      integer*4 lenw,loc
      character*(*) ,intent(in):: name,str
      loc=-1
      kax=ktfsymbolz(name,lenw(name))-4
      call tflocal(klist(kax))
      kas=ktsalocb(0,str,len_trim(str))
      klist(kax)=ktfstring+kas
      return
      end
      
      integer*8 function ktsalocbi(mode,string,i,n)
      use tfstk, kf => ktsalocbi
      implicit none
      integer*4 ,intent(in):: mode,n,i
      character,intent(in):: string(n)
      ktsalocbi=kf(mode,string,i,n)
      return
      end

      integer*8 function ktaaloc(mode,n)
      use tfstk, kf => ktaaloc
      implicit none
      integer*4 ,intent(in):: mode,n
      ktaaloc=kf(mode,n)
      return
      end

      integer*8 function ktadaloc(mode,n)
      use tfstk, kf => ktadaloc
      implicit none
      integer*4 ,intent(in):: mode,n
      ktadaloc=kf(mode,n)
      return
      end

      integer*8 function ktavaloc(mode,n)
      use tfstk, kf => ktavaloc
      implicit none
      integer*4 ,intent(in):: mode,n
      ktavaloc=kf(mode,n)
      return
      end

      integer*8 function ktfsymbolz(name,l)
      use tfstk, kf => ktfsymbolz
      implicit none
      integer*4 ,intent(in):: l
      character ,intent(in):: name(l)
      ktfsymbolz=kf(name,l)
      return
      end

      integer*8 function ktsalocb(mode,str,leng)
      use tfstk, kf => ktsalocb
      implicit none
      integer*4 ,intent(in):: mode,leng
      character ,intent(in):: str(leng)
      ktsalocb=kf(mode,str,leng)
      return
      end
