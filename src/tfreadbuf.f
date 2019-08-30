      subroutine tfreadbuf(icmd,lfn,ib,is,nc)
      use tfstk
      use tfrbuf
      use tfshare
      use tfcsi
      use iso_c_binding
      implicit none
      integer*8 ib,is,ia,ls1,mapresizefile,lenfile,ib1
      integer*4 icmd,lfn,nc,j,itfgetbuf,irtc,ls,ie,i
      if(lfn .le. 0)then
        select case (icmd)
        case (irbgetpoint)
          nc=-1
          return
        case (irbinit)
          go to 9100
        case (irbopen)
          do j=nbuf,11,-1
            if(itbuf(j) .eq. modeclose)then
              lfn=j
              call irbopen1(lfn,ib,is,nc)
              return
            endif
          enddo
          lfn=0
          go to 9100
        case (irbreadrecordbuf,irbreadrecord)
          go to 9000
        case (irbreadbuf)
          go to 9100
        case (irbibuf)
          if(itbuf(lfn) .eq. is)then
            ib=ibuf(lfn)
          else
            ib=0
          endif
        end select
        return
      endif
      select case (icmd)
      case (irbmovepoint)
        mbuf(lfn)=min(lbuf(lfn)+1,max(1,mbuf(lfn)+nc))
      case (irbbor)
        mbuf(lfn)=lbuf(lfn)+1
      case (irbeor2bor)
        if(mbuf(lfn) .eq. lbuf(lfn))then
          mbuf(lfn)=lbuf(lfn)+1
        endif
      case (irbgetpoint)
c        write(*,*)'irbgetpoint ',lfn,mbuf(lfn),lbuf(lfn)
        if(ibuf(lfn) .eq. 0)then
          go to 9000
        elseif(mbuf(lfn) .gt. lbuf(lfn)
     $         .or. mbuf(lfn) .le. 0)then
          nc=-1
        else
          nc=lbuf(lfn)-mbuf(lfn)+1
          ib=ibuf(lfn)
          is=mbuf(lfn)
        endif
      case (irbinit)
        itbuf(lfn)=int(ib)
        if(ib .le. modewrite)then
          if(ibuf(lfn) .eq. 0)then
            ibuf(lfn)=ktaloc(maxlbuf/8)
            lenbuf(lfn)=maxlbuf
          endif
          ilist(2,ibuf(lfn)-1)=0
        endif
        lbuf(lfn)=0
        mbuf(lfn)=1
      case (irbreset)
        if(itbuf(lfn) .le. modewrite)then
          lbuf(lfn)=0
          mbuf(lfn)=1
        else
          mbuf(lfn)=lbuf(lfn)+1
        endif
      case (irbclose)
        select case (itbuf(lfn))
        case (modewrite)
          close(lfn)
          if(ibuf(lfn) .gt. 0)then
            if(ilist(2,ibuf(lfn)-1) .ne. 0)then
              call unixclose(ilist(2,ibuf(lfn)-1))
              ilist(2,ibuf(lfn)-1)=0
            endif
            call tfree(ibuf(lfn))
          endif
        case(modeclose,moderead)
          close(lfn)
        case (modestring)
          call tflocal1(ibuf(lfn)-1)
          ibuf(lfn)=0
        case (modeshared)
          ia=ibuf(lfn)
          call tfreeshared(ia)
          ibuf(lfn)=0
        case default
          call unmapfile(klist(ibuf(lfn)),int8(lenbuf(lfn)))
          call unixclose(ifd(lfn))
        end select
        lbuf(lfn)=0
        mbuf(lfn)=1
        itbuf(lfn)=modeclose
      case (irbopen)
        do j=nbuf,11,-1
          if(itbuf(j) .eq. modeclose)then
            lfn=j
            call irbopen1(lfn,ib,is,nc)
            return
          endif
        enddo
        lfn=0
        go to 9100
      case (irbreadrecordbuf,irbreadrecord)
        if(itbuf(lfn) .le. modewrite)then
          if(icmd .eq. irbreadrecordbuf)then
            nc=itfgetbuf(lfn,ilist(1,ibuf(lfn)),maxlbuf,irtc)
            lenbuf(lfn)=nc
            lbuf(lfn)=0
            mbuf(lfn)=1
          else
            nc=itfgetbuf(lfn,jlist(ib,ibuf(lfn)),
     $           maxlbuf-int(ib)-256,irtc)
            lbuf(lfn)=int(ib)-1
            lenbuf(lfn)=lbuf(lfn)+nc
            mbuf(lfn)=int(ib)
          endif
          if(irtc .ne. 0)then
            return
          endif
        else
 11       ls=lenbuf(lfn)
          if(lbuf(lfn) .lt. ls .and.
     $         jlist(lbuf(lfn)+1,ibuf(lfn)) .eq. 10)then
            if(lfn1 .gt. 0)then
              write(lfn1,*)
            endif
            lbuf(lfn)=lbuf(lfn)+1
          endif
c          if(ls .eq. 1)then
c            write(*,*)'irbrr ',lbuf(lfn),mbuf(lfn)
c          endif
          if(lbuf(lfn) .lt. ls)then
            ie=ls
            do i=lbuf(lfn)+1,ls
              if(jlist(i,ibuf(lfn)) .eq. 10)then
                if(i .eq. 1 .or. jlist(i-1,ibuf(lfn)) .ne.
     $               ichar('\\'))then
                  ie=i
                  exit
                endif
              endif
            enddo
            nc=int(ie-lbuf(lfn))
            mbuf(lfn)=lbuf(lfn)+1
            lbuf(lfn)=ie
            if(icmd .eq. irbreadrecord)then
c              if(ib .eq. 0)then
c                nc=min(nc,len(buff))
c                call tmovb(jlist(mbuf(lfn),ibuf(lfn)),buff,nc)
c              else
                is=int8(min(mbuf(lfn),ie))
c              endif
            endif
          else
            if(itbuf(lfn) .ge. modemapped)then
              ls1=lenfile(ifd(lfn))
              if(ls .lt. ls1)then
                ib1=mapresizefile(klist(ibuf(lfn)),ifd(lfn),
     $               ls,ls1)/8
                if(ib1 .lt. 0)then
                  go to 9000
                endif
                mbuf(lfn)=ls
                ibuf(lfn)=ib1
                lenbuf(lfn)=int(ls1)
                go to 11
              endif
            elseif(mbuf(lfn) .lt. lbuf(lfn))then
              nc=1
              mbuf(lfn)=lbuf(lfn)
              is=lbuf(lfn)
              return
            endif
            mbuf(lfn)=ls+1
            go to 9000
          endif
        endif
      case (irbassign)
        if(ibuf(lfn) .ne. 0)then
          call c_f_pointer(c_loc(jlist(1,ibuf(lfn))),buffer)
        else
          buffer=>buffer0
          ibuf(lfn)=transfer(c_loc(buffer0),ib)/8
        endif
        ipoint=>mbuf(lfn)
        lrecl=>lbuf(lfn)
        lfni=lfn
      case (irbreadbuf)
        if(ibuf(lfn) .eq. 0 .or. mbuf(lfn) .le. 0 .or.
     $       mbuf(lfn) .gt. lbuf(lfn))then
          nc=0
        else
          nc=lbuf(lfn)-mbuf(lfn)
          if(jlist(lbuf(lfn),ibuf(lfn)) .ne. 10)then
            nc=nc+1
          endif
          nc=max(nc,1)
          ipoint=min(mbuf(lfn),lenbuf(lfn))
        endif
      case (irbsetbuf)
      case (irbcloseinp)
        if(itbuf(lfn) .le. 2)then
          itbuf(lfn)=modeclose
          if(ibuf(lfn) .gt. 0)then
            if(ilist(2,ibuf(lfn)-1) .ne. 0)then
              call unixclose(ilist(2,ibuf(lfn)-1))
              ilist(2,ibuf(lfn)-1)=0
            endif
          endif
        endif
      case (irbsetinp)
        if(itbuf(lfn) .le. modewrite)then
          if(ibuf(lfn) .gt. 0)then
            ilist(2,ibuf(lfn)-1)=int(ib)
          endif
        endif
      case (irbibuf)
        if(itbuf(lfn) .eq. is)then
          ib=ibuf(lfn)
        else
          ib=0
        endif
      end select
      return
 9000 nc=-99
      return
 9100 nc=-999
      return
      end

      subroutine irbopen1(j,ib,is,nc)
      use tfrbuf
      use tfstk
      implicit none
      integer*4 j,nc
      integer*8 ib,is
      if(itbuf(j) .eq. 0)then
        itbuf(j)=int(is)
        lenbuf(j)=0
        ifd(J)=0
        select case (int(is))
        case (moderead,modewrite)
        case (modestring)
          ibuf(j)=ktfcopy1(ib)+1
          lenbuf(j)=ilist(1,ib)
        case (modeshared)
          ibuf(j)=ib
          lenbuf(j)=ilist(1,ib-1)
        case default
          ibuf(j)=ib
          lenbuf(j)=int(is-modemapped)
          ifd(j)=nc
        end select
        lbuf(j)=0
        mbuf(j)=1
c        write(*,*)'irbopen1 ',j,ibuf(j)
        return
      endif
      return
      end
 
      subroutine readstr(in,str,irtc)
      use tfrbuf
      use tfcsi, only:buffer
      implicit none
      integer*4 in,irtc,nc
      integer*8 is
      character*(*) str
c      write(*,*)'reststr ',in
      call tfreadbuf(irbreadrecord,in,int8(lbuf(in)+1),is,nc)
c      write(*,*)': ',nc,'''',buffer(mbuf(in):mbuf(in)+nc-1),''''
      irtc=0
      if(nc .gt. 0)then
        str=buffer(mbuf(in):mbuf(in)+nc-1)
        mbuf(in)=mbuf(in)+nc
        if(str(nc:nc) .eq. char(10))then
          str(nc:)=' '
        else
          str(nc+1:)=' '
        endif
      else
        str=' '
        if(nc .lt. 0)then
          irtc=-1
        endif
      endif
      return
      end

      subroutine tfopenshared(isp1,kx,irtc)
      use tfstk
      use tfrbuf
      use tfshare
      implicit none
      type (sad_descriptor) kx
      integer*8 ia
      integer*4 isp1,irtc,itfmessage,n,m,iu,nc
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      n=int(rtastk(isp))
      if(n .lt. 0)then
        irtc=itfmessage(9,'General::wrongval','"Non-negative"')
        kx%k=kxfailed
        irtc=0
        return
      endif
      m=(n+7)/8+2
      irtc=0
      ia=ktfallocshared(m+1)
c      ia=mapallocfixed8(rlist(0), m+1, 8, irtc)
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'General::mmap','""')
        kx%k=kxfailed
        irtc=0
        return
      endif
      call tfreadbuf(irbopen,iu,ia,int8(modeshared),nc)
      if(iu .le. 0)then
        call tfreeshared(ia)
        irtc=itfmessage(9,'General::fileopen','"(Shared)"')
        return
      endif
      ilist(1,ia)=m
      ilist(2,ia)=0
      klist(ia+3)=ktfoper+mtfnull
      kx=dfromr(dble(iu))
      return
      end

      subroutine tfreadshared(isp1,kx,irtc)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*8 ia
      integer*4 isp1,irtc,itfmessage,isp0,iu,nc,ist
      logical*4 tfcheckelement
c      call tfdebugprint(ktastk(isp),'readshard',3)
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp),iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      irtc=0
      call tfreadbuf(irbibuf,iu,ia,int8(4),nc)
      if(ia .eq. 0)then
        kx%k=kxeof
        return
      endif
      ist=ilist(2,ia)
      do while(ist .ne. 0)
        if(ist .ne. 1 .and. ist .ne. -1)then
          write(*,*)'tfreadshared shared memory destructed ',ist,ia
          call abort
        endif
        call tpause(10000)
        ist=ilist(2,ia)
      enddo
      ilist(2,ia)=1
      kx=dlist(ia+1)
c      write(*,*)'readshared ',ia,kx%k
      if(ktfobjq(kx))then
c        write(*,*)'readshared-obj '
        if(ktfsymbolq(kx))then
c          write(*,*)'readshared-symbol '
          if( .not. tfcheckelement(kx,.false.))then
            irtc=itfmessage(99,'General::wrongtype',
     $           '"undefined symbol returned in Shared"')
            kx%k=ktfoper+mtfnull
            ilist(2,ia)=0
            return
          endif
        elseif(ktfstringq(kx,str))then
c          call tfdebugprint(kx,'readshared-string',1)
c          write(*,*)'at ',ia
          kx=kxsalocb(-1,str%str(1:str%nch),str%nch)
        elseif(ktflistq(kx))then
c          write(*,*)'readshared-list '
          isp0=isp
          call tfrecallshared(isp0,ktflist+ia+3,kx,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            kx%k=ktfoper+mtfnull
          endif
        else
c          write(*,*)'readshared-other '
          kx%k=ktfoper+mtfnull
        endif
      endif
      ilist(2,ia)=0
      return
      end

      subroutine tfwriteshared(isp1,kx,irtc)
      use tfstk
      use tfrbuf
      implicit none
      integer*8 kx,kas,ka,kt,kap,k
      integer*4 isp1,irtc,itfmessage,itfmessageexp,isp0,n,i,iu,nc,ist
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      iu=int(rtastk(isp1+1))
      call tfreadbuf(irbibuf,iu,kas,int8(4),nc)
      if(kas .eq. 0)then
        irtc=itfmessage(99,'Shared::notopen','""')
        return
      endif
      ist=ilist(2,kas)
      do while(ist .ne. 0)
        if(ist .ne. 1 .and. ist .ne. -1)then
          write(*,*)'writeshared shared memory destructed: ',ist,kas
          call abort
        endif
        call tpause(10000)
        ist=ilist(2,kas)
      enddo
      ilist(2,kas)=-1
      k=ktastk(isp)
      if(ktfnonobjq(k) .or. ktfsymbolq(k))then
        klist(kas+1)=k
      else
        ka=ktfaddr(k)
        kt=k-ka
        if(kt .eq. ktfstring)then
          if(ilist(1,ka) .gt. ilist(1,kas)*8)then
            ilist(2,kas)=0
            irtc=itfmessageexp(9,'Shared::toolarge',
     $         dble(ilist(1,ka)-ilist(1,kas)*8))
            return
          endif
          call tmov(ilist(1,ka+1),ilist(1,kas+3),
     $         (ilist(1,ka)+7)/8)
          ilist(1,kas+2)=ilist(1,ka)
          klist(kas+1)=ktfstring+kas+2
c          write(*,*)'writeshared-string ',kas,klist(kas+1)
        elseif(kt .eq. ktflist)then
          isp0=isp
          call tfsharedsize(isp0,k,n,irtc)
          if(irtc .ne. 0)then
            ilist(2,kas)=0
            isp=isp1+2
            return
          endif
          if(n .gt. ilist(1,kas))then
            ilist(2,kas)=0
            irtc=itfmessageexp(9,'Shared::toolarge',
     $           dble((n-ilist(1,kas))*8))
            return
          endif
          kap=kas+2
          do i=isp1+3,isp
            n=itastk2(1,i)
            ktastk2(i)=kap
            kap=kap+n
          enddo
          do i=isp1+3,isp
            call tfstoreshared(isp1+3,ktastk(i),ktastk2(i))
          enddo
          isp=isp1+2
          klist(kas+1)=ktflist+kas+3
        else
          ilist(2,kas)=0
          irtc=itfmessage(9,'General::wrongtype',
     $         '"No Patterns"')
          return
        endif
      endif
      ilist(2,kas)=0
      irtc=0
      kx=ktfoper+mtfnull
      return
      end
