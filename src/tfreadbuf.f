      subroutine tfreadbuf(lfn,ib,nc)
      use tfstk
      use tfshare
      use tfcsi
      use tfrbuf
c      use iso_c_binding
      implicit none
      integer*8 ls1,mapresizefile,lenfile,ib1
      integer*4 lfn,itfgetbuf,irtc,ls,ie,i,nc,ib
      if(lfn .le. 0 .or. ibuf(lfn) .eq. 0)then
        go to 9000
      endif
      if(itbuf(lfn) .le. modewrite)then
c        write(*,*)'tfreadbuf ',lfn,itbuf(lfn),modewrite
        nc=itfgetbuf(lfn,jlist(ib,ibuf(lfn)),
     $       maxlbuf-ib-256,irtc)
        if(irtc .ne. 0)then
          return
        endif
        lbuf(lfn)=ib-1
        lenbuf(lfn)=lbuf(lfn)+nc
c        write(*,*)'tfreadbuf ',lfn,nc,lbuf(lfn),mbuf(lfn),ib
        mbuf(lfn)=ib
      else
 11     ls=lenbuf(lfn)
c        if(lbuf(lfn) .lt. ls .and.
c     $       jlist(lbuf(lfn)+1,ibuf(lfn)) .eq. 10)then
c          lbuf(lfn)=lbuf(lfn)+1
c        endif
        if(lbuf(lfn) .lt. ls)then
          ie=ls
          do i=lbuf(lfn)+1,ls
            if(jlist(i,ibuf(lfn)) .eq. 10)then
              if(i .eq. 1 .or. jlist(i-1,ibuf(lfn)) .ne.
     $             ichar('\\'))then
                ie=i
                exit
              endif
            endif
          enddo
          nc=ie-lbuf(lfn)
          mbuf(lfn)=lbuf(lfn)+1
          lbuf(lfn)=ie
        else
          if(itbuf(lfn) .ge. modemapped)then
            ls1=lenfile(ifd(lfn))
            if(ls .lt. ls1)then
              ib1=mapresizefile(klist(ibuf(lfn)),ifd(lfn),
     $             ls,ls1)/8
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
            return
          endif
          mbuf(lfn)=ls+1
          go to 9000
        endif
      endif
      return
 9000 nc=irbeof
      return
      end subroutine

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
      character*(*) str
c      write(*,*)'reststr ',in
      call tfreadbuf(in,lbuf(in)+1,nc)
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
      call trbopen(iu,ia,int8(modeshared),nc)
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
      integer*4 isp1,irtc,itfmessage,isp0,iu,ist
      logical*4 tfcheckelement
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp),iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      irtc=0
      ia=itrbibuf(iu,modeshared)
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
      integer*4 isp1,irtc,itfmessage,itfmessageexp,isp0,n,i,iu,ist
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      iu=int(rtastk(isp1+1))
      kas=itrbibuf(iu,modeshared)
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
