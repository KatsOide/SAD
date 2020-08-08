      subroutine tfmemcheck(isp1,kx,irtc)
      use tfstk
      use tfmem
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx5,kld
      type (sad_rlist), pointer :: klri,klr
      integer*8 ip1,ip0,ip,ix,i,j, i1, ki,m,ih,mf
      integer*4 isp1,irtc,l,nf,isp0,itfmessage
      real*8 xnmem,xnnet
      character*12 tftypename
      logical*4 tfcheckelement,freel,check
      if(isp1+1 .eq. isp .and. ktfrealq(ktastk(isp)))then
        check=.true.
        freel=rtastk(isp) .eq. 2.d0
      else
        check=.false.
        freel=.false.
      endif
      l=-1
      ih=-1
      if(check)then
        call tfcheckcont(itfcontext,irtc)
        if(irtc .ne. 0)then
          return
        endif
        do l=1,levele
          j=itflocal+l
          i1=j
          i=klist(j)
          do while(i .ne. j)
            if(i .le. 0)then
              write(*,*)'Wrong temporary element:',i
              go to 9100
            endif
            if(klist(i) .le. 0)then
              write(*,*)'Wrong temporary element pointed:',i,klist(i)
              go to 9100
            endif
            ki=iand(ktfmask,klist(i))+i+2
            if(.not. tfcheckelement(ki,.true.))then
              write(*,*)' Error in a temporary element:'
              call tfdebugprint(ki,' ',3)
              go to 9100
            endif
            i1=i
            i=ktfaddr(klist(i))
          enddo
        enddo
      endif
      nf=0
      mf=0
      isp0=isp
      do ix=0,nindex
        ip0=icp+ix*2
        ip=klist(ip0)
        do while(ip .ne. ip0)
          ip1=klist(ip)
          if(klist(ip1+1) .ne. ip)then
            write(*,*)
     $           'Free area (size-list) is inconsistent: area =',ip,
     $           ' next area =',ip1,' previous(next) area =',
     $           klist(ip1+1),' index =',ix
            go to 9200
          endif
          nf=nf+1
          if(ix .lt. nindex)then
            m=ix+1
            if(ilist(1,ip-1) .ne. m)then
              write(*,*)
     $             'Free area size is inconsistent: area =',ip,
     $             ' index =',m,' stored size =',ilist(1,ip-1)
              go to 9200
            endif
          else
            m=ilist(1,ip-1)
          endif
          if(.not. tfchecklastp(ip) .or.
     $         .not. tfchecklastp(ip+m-1))then
            write(*,*)
     $       'Free area is out of allocated block: area =',ip,
     $           'size =',m
            go to 9200
          endif
          mf=mf+m
          if(freel)then
            isp=isp+2
            ktastk(isp-1)=ip-1
            ktastk(isp)=m
          endif
          ip=ip1
        enddo
      enddo
      xnnet=nnet
      xnmem=nmem
      if(freel)then
        kx=kxadaloc(-1,5,kld)
        call descr_rlist(kx,klr)
        kld%dbody(5)=kxadaloc(0,nf,klx5)
        do i=1,nf*2,2
          klx5%dbody((i+1)/2)=kxavaloc(0,2,klri)
          klri%rbody(1)=ktastk(isp0+i)
          klri%rbody(2)=ktastk(isp0+i+1)
        enddo
        isp=isp0
      else
        kx=kxavaloc(-1,4,klr)
      endif
      klr%rbody(1)=xnnet
      klr%rbody(2)=xnmem
      klr%rbody(3)=dble(nf)
      klr%rbody(4)=mf+xnnet-xnmem
      irtc=0
      return
 9100 write(*,*)'Level = ',l,
     $     ', root =',j,
     $     ', temporary element =',i,
     $     ', previous =',i1
      write(*,*)', next =',ktfaddr(klist(i)),
     $     ', type = ',tftypename(klist(i)),
     $     ', size =',ilist(1,i-1)
 9200 irtc=itfmessage(999,'General::memory',' ')
      kx%k=ktfoper+mtfnull
      return
      end

      character*(*) function tftypename(k)
      use tfstk
      implicit none
      integer*8 k
      select case (ktftype(k))
      case (ktflist)
        tftypename='List'
      case (ktfsymbol)
        tftypename='Symbol'
      case (ktfpat)
        tftypename='Pattern'
      case (ktfstring)
        tftypename='String'
      case (ktfoper)
        tftypename='Function'
      case (ktfref)
        tftypename='Reference'
      case default
        if(ktfrealq(k))then
          tftypename='Real'
        else
          tftypename='Undefined'
        endif
      end select
      return
      end

      recursive subroutine tfcheckcont(icont,irtc)
      use tfstk
      use tfmem
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_namtbl), pointer :: loc
      integer*8 k0,kad,kad0,kadi,kadi0,ks,ii,kt,ka,icont,ih
      integer*8 i,j,k,i1,k1
      integer*4 irtc,ig,nc,l,m,itfmessage
      character*16 name
      logical*4 tfcheckelement
      irtc=0
      do m=0,nsymhash
        j=icont+m+1
        i=klist(j)
        do while(i .ne. j)
          call loc_namtbl(i,loc)
          nc=loc%str%nch
          if(nc .le. 0)then
            write(*,*)'Zero or negative length for symbol name: ',nc
            go to 9000
          endif
          if(loc%len .ne. nc/8+9)then
            write(*,*)'Inconsistent loc table size: ',loc%len,
     $           ' <> ',nc/8+9
            go to 9000
          endif
c          write(*,*)'checkcont ',loc%str%str(1:nc)
          k=loc%symdef
          k0=ksad_loc(loc%symdef)
          do while(k .gt. 0)
            if(klist(k+1) .ne. k0)then
              write(*,*)'Wrong previous name table: ',klist(k+1),
     $             ' <> ',k0
              go to 9000
            endif
            ks=klist(k+4)
            if(ilist(2,k+7) .eq. -3)then
              if(ktfaddr(ks) .ne. icont)then
                call tfcheckcont(ktfaddr(ks),irtc)
                if(irtc .ne. 0)then
                  return
                endif
              endif
            elseif(iand(ktfmask,ks) .ne. ktfref)then
              if(.not. tfcheckelement(ks,.false.))then
                go to 9000
              endif
            endif
            do l=0,1
              kad0=k+2+l
              kad=klist(kad0)
              do while(kad .ne. 0)
                if(kad .lt. 0)then
                  write(*,*)
     $                 'Negative definition pointer: ',
     $                 kad,' previous pointer =',kad0
                  go to 9000
                endif
                if(.not. tfchecklastp(kad))then
                  write(*,*)
     $                 'Out of memory block ',kad,
     $                 ' block = ',itfcbk(kad)
                  go to 9000
                endif
                if(klist(kad+1) .ne. kad0)then
                  write(*,*)'Wrong previous definition:',klist(kad+1),
     $                 ' <>',kad0
                  go to 9000
                endif
                if(ilist(1,kad+2) .eq. maxgeneration)then
                  do ii=kad+3,kad+ilist(2,kad+2)+3
                    ih=ii-kad-3
                    kadi0=ii
                    kadi=klist(ii)
                    do while(kadi .ne. 0)
                      call tfcheckdef(kadi,kadi0,irtc)
                      if(irtc .ne. 0)then
                        go to 9000
                      endif
                      kadi0=kadi
                      kadi=klist(kadi)
                    enddo
                  enddo
                else
                  call tfcheckdef(kad,kad0,irtc)
                  if(irtc .ne. 0)then
                    go to 9000
                  endif
                endif
                kad0=kad
                kad=klist(kad)
              enddo
            enddo
            if(ilist(1,k-1) .ne. 10)then
              if(ktfoperq(klist(k+4)))then
              else
                write(*,*)'Inconsistent size of generation table:',
     $               ilist(1,k-1)
                go to 9000
              endif
            endif
            k1=klist(k)
            if(k1 .lt. 0)then
              write(*,*)
     $             'Negative pointer to next generation table:',k1
              go to 9000
            endif
            if(klist(k+8) .ne. i)then
              write(*,*)
     $             'Inconsistent pointer to loc table: ',
     $             klist(k+8),' name table =',i+5
              go to 9000
            endif
            k0=k
            k=k1
          enddo
          i1=loc%next
          if(i1 .le. 0)then
            write(*,*)'Zero or negative pointer to next name:',i1
            go to 9000
          endif
          i=i1
        enddo
      enddo
      return
 9000 nc=loc%str%nch
      name=loc%str%str(1:nc)
      write(*,*)'Name = ',name(1:nc),
     $     ', context = ',icont,
     $     ', symbol hash =',m,
     $     ', arg hash =',ih,
     $     ', upvalue =',l,
     $     ', generation table =',k,
     $     ', symbol location =',i,
     $     ', element type =',kt,
     $     ', element pointer =',ka,
     $     ', generation =',ig,
     $     ', stack pointer =',isp,
     $     ', current generation =',lgeneration
      irtc=itfmessage(999,'General::memory',' ')
      return
      end

      subroutine tfcheckdef(kad,kad0,irtc)
      use tfstk
      implicit none
      integer*8 kad,kad0
      integer*4 irtc
      logical*4 tfcheckelement
      irtc=1
      if(kad .lt. 0)then
        write(*,*)
     $       'Negative definition pointer:',kad
        return
      endif
      if(klist(kad+1) .ne. kad0)then
        write(*,*)
     $       'Inconsistent definition table:',
     $       ' at ',kad,' previous table =',kad0,
     $       ' <> stored previous table',klist(kad+1)
        return
      endif
      if(.not. tfcheckelement(ktflist+klist(kad+3),.false.))then
        write(*,*)'Error in original arg',kad
        return
      endif
      if(.not. tfcheckelement(ktflist+klist(kad+4),.false.))then
        write(*,*)'Error in replaced arg',kad
        return
      endif
      if(.not. tfcheckelement(klist(kad+4),.false.))then
        write(*,*)'Error in original value',kad
        return
      endif
      if(.not. tfcheckelement(klist(kad+4),.false.))then
        write(*,*)'Error in replaced value',kad
        return
      endif
      irtc=0
      return
      end

      recursive logical*4 function tfcheckelement(k,temp)
     $     result(lx)
      use tfstk
      implicit none
      integer*8 k,ka,kt, kt1,loc
      integer*4 i,m,nc,nw,lg,nrecmax
      parameter (nrecmax=1024)
      integer*4 nrec
      logical*4 temp,nr
      character*10 tfkname
      save nrec
      data nrec/0/
      lx=.false.
      if(ktfrealq(k))then
        lx=.true.
        return
      endif
      ka=ktfaddr(k)
      kt=k-ka
      if(ka .lt. 0)then
        write(*,*)'Negative pointer:'
        return
      endif
      if(kt .eq. ktfoper)then
        if(ka .lt. 0 .or. ka .gt. 16384)then
          write(*,*)'Illegal function number: ',ka
          return
        endif
      else
        if(ka .eq. 0)then
          write(*,*)'Zero pointer:'
          return
        endif
        if(.not. tfchecklastp(ka))then
          write(*,*)
     $         'Element is out of allocated block: element =',ka
          return
        endif
        kt1=iand(ktfmask,klist(ka-2))
        if(kt .eq. ktfsymbol)then
          if(ilist(2,ka-3) .ne. 0)then
            if(kt1 .ne. ktfsymbol)then
              write(*,*)
     $             'Inconsistent element type for definition:',
     $             tfkname(kt),tfkname(kt1)
              return
            endif
          else
            if(kt1 .eq. ktfstring)then
              if(ilist(2,ka-1) .ne. ka)then
                write(*,*)
     $               'Inconsistent element type (string as symbol):',
     $               ka,ilist(2,ka-1)
                return
              endif
            elseif(kt1 .ne. ktfsymbol)then
              write(*,*)
     $             'Inconsistent element type, should be Symbol, but ',
     $             tfkname(kt1)
              return
            endif
          endif
        elseif(kt1 .ne. kt)then
          write(*,*)'Inconsistent element type, should be ',
     $         tfkname(kt),', but ',tfkname(kt1)
          return
        endif
        if(ilist(1,ka-1) .le. 0 .and. .not. temp)then
          write(*,*)'Wrong degeneration count:',ilist(1,ka-1)
          return
        endif
        if(kt .eq. ktflist)then
          m=ilist(2,ka-1)
          if(m .lt. 0)then
            write(*,*)'Negative length of list:',m
            return
          endif
          nrec=nrec+1
          if(nrec .lt. nrecmax)then
            if(.not. tfcheckelement(klist(ka),.false.))then
              write(*,*)' Error in head of list.'
              call tfdebugprint(klist(ka),
     $             'Error in head of list:',3)
              return
            endif
            if(ktfreallistq(ka))then
              do i=1,m
                if(ktfnonrealq(klist(ka+i)))then
                  write(*,*)' Non-real in a real list, index =',i,
     $                 ' size =',m
c                  call tfdebugprint(klist(ka),'Head:',3)
c                  call tfdebugprint(klist(ka+i),
c     $             'Non-real in a real list:',3)
                  return
                endif
              enddo
            elseif(m .ne. 0)then
              nr=.false.
              do i=1,m
                nr=nr .or. ktfnonrealq(klist(ka+i))
                if(.not. tfcheckelement(klist(ka+i),.false.))then
                  return
                endif
              enddo
              if(.not. nr)then
                write(*,*)'All reals in a non-real list'
c                call tfdebugprint(k,' ',1)
                return
              endif
            endif
          endif
          nrec=nrec-1
        elseif(kt .eq. ktfstring)then
          nc=ilist(1,ka)
          nw=ilist(1,ka-3)-4
          if(nw .lt. nc/8+1)then
            write(*,*)'Inconsistent length of string:',
     $           nc,' size of allocated memory:',nw
            return
          endif
        elseif(kt .eq. ktfsymbol)then
          nw=ilist(1,ka-3)
c          if((nw .lt. 4 .and. nw .ne. 0) .or. nw .gt. 9)then
c            if(it .ne. ntfsymbol .or. ilist(2,ia-1) .ne. ka)then
c              write(*,*)'Inconsistent length of symbol: ',nw
c              return
c            endif
c          endif
          lg=ilist(2,ka-1)
          if(lg .lt. -3)then
            write(*,*)'Illegal generation of symbol:',lg
            return
          endif
          loc=klist(ka)
          if(loc .le. 0)then
            write(*,*)'Illegal loc table of symbol:',loc
            return
          endif
        elseif(kt .eq. ktfpat)then
          nw=ilist(1,ka-3)
c          if(nw .lt. 13 .or. nw .gt. 16)then
c            write(*,*)
c     $           'Inconsistent size of pattern: length =',nw
c            return
c          endif
        else
          write(*,*)'Undefined type',kt
          return
        endif
      endif
      lx=.true.
      return
      end

      subroutine tfgarbagecollect(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc
      return
      end

      subroutine tfmemcheckprint1(tag,k,pri)
      use tfstk
      implicit none
      integer*4 irtc,k
      character*(*) tag
      logical*4 pri
      call tfmemcheckprint(tag,k,pri,irtc)
      return
      end

      subroutine tfmemcheckprint(tag,k,pri,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1
      integer*4 isp0,irtc,itfdownlevel,i,k
      character*(*) tag
      logical*4 pri
      isp0=isp
      isp=isp+1
      rtastk(isp)=1.d0
      levele=levele+1
      call tfmemcheck(isp0,k1,irtc)
      if(irtc .ne. 0)then
        write(*,*)'memcheck-error: ',tag,' ',k
        call tfreseterror
      elseif(pri)then
        write(*,*)'memcheck-ok: ',tag,' ',k
      endif
      i=itfdownlevel()
      isp=isp0
      return
      end

      subroutine tfprintnan
      use tfstk,only:dnotanumber
      implicit none
      integer*8 kfromr
      real*8 x0,x1
      x1=-1.d0
      x0=0.d0
      write(*,'(a,8('' '',z16))')'NaN: ',
     $     kfromr(dnotanumber),kfromr(x0/(x1-x1)),
     $     kfromr(sqrt(x1)),kfromr(log(x1)),kfromr(x0*(x1/x0)),
     $     kfromr((-x0)*((-x1)/x0)),kfromr(x1/x0)
      write(*,*)x0/x0,x1/x0,x0/(x1-x1)
      return
      end
