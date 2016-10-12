      subroutine tfreadbuf(icmd,lfn,ib,is,nc,buff)
      use tfstk
      use tfrbuf
      use tfshare
      implicit none
      integer*8 ib,is,ia
      integer*4 icmd,lfn,nc,i,j,ls,ie,itfgetbuf,irtc
      character*(*) buff
      select case (icmd)
      case (irbmovepoint)
        if(lfn .gt. 0)then
          mbuf(lfn)=min(lbuf(lfn),max(0,mbuf(lfn)+nc))
        endif
      case (irbbor)
        if(lfn .gt. 0)then
          mbuf(lfn)=lbuf(lfn)+1
        endif
      case (irbeor2bor)
        if(lfn .gt. 0)then
          if(mbuf(lfn) .eq. lbuf(lfn))then
            mbuf(lfn)=lbuf(lfn)+1
          endif
        endif
      case (irbgetpoint)
        if(lfn .eq. 0)then
          nc=-1
          return
        elseif(ibuf(lfn) .eq. 0 .and. itbuf(lfn) .le. 1)then
          ibuf(lfn)=ktaloc(maxlbuf/8)
          ilist(2,ibuf(lfn)-1)=0
          lbuf(lfn)=0
          mbuf(lfn)=-1
        endif
        if(ibuf(lfn) .eq. 0)then
          go to 9000
        elseif(mbuf(lfn) .gt. lbuf(lfn)
     $         .or. mbuf(lfn) .lt. 0)then
          nc=-1
        else
          nc=lbuf(lfn)-mbuf(lfn)
          ib=ibuf(lfn)
          is=mbuf(lfn)+1
        endif
      case (irbinit)
        if(lfn .gt. 0)then
          itbuf(lfn)=ib
          if(ib .le. 2)then
            if(ibuf(lfn) .eq. 0)then
              ibuf(lfn)=ktaloc(maxlbuf/8)
            endif
            ilist(2,ibuf(lfn)-1)=0
          endif
          lbuf(lfn)=0
          mbuf(lfn)=-1
        else
          go to 9100
        endif
      case (irbreset)
        if(lfn .gt. 0)then
          if(itbuf(lfn) .le. 2)then
            lbuf(lfn)=0
            mbuf(lfn)=-1
          else
            mbuf(lfn)=lbuf(lfn)
          endif
        endif
      case (irbclose)
        if(lfn .gt. 0)then
          if(itbuf(lfn) .le. 2)then
            close(lfn)
            if(ibuf(lfn) .gt. 0)then
              if(ilist(2,ibuf(lfn)-1) .ne. 0)then
                call unixclose(ilist(2,ibuf(lfn)-1))
                ilist(2,ibuf(lfn)-1)=0
              endif
            endif
          elseif(itbuf(lfn) .eq. 3)then
            call tflocal1(ibuf(lfn)-1)
            ibuf(lfn)=0
          elseif(itbuf(lfn) .eq. 4)then
            ia=ibuf(lfn)
            call tfreeshared(ia)
c            irtc=mapfree(rlist(ia))
            ibuf(lfn)=0
c            if(irtc .ne. 0)then
c              go to 9000
c            endif
          endif
          lbuf(lfn)=0
          mbuf(lfn)=-1
          itbuf(lfn)=0
        endif
      case (irbopen)
        do j=nbuf,11,-1
          if(itbuf(j) .eq. 0)then
            lfn=j
            itbuf(j)=is
            if(is .eq. 3)then
              ibuf(j)=ktfcopy1(ib)+1
            elseif(is .eq. 4)then
              ibuf(j)=ib
            endif
            lbuf(j)=0
            mbuf(j)=-1
            return
          endif
        enddo
        lfn=0
        go to 9100
      case (irbreadrecordbuf)
        if(itbuf(lfn) .le. 2)then
          if(lfn .gt. 0)then
            nc=itfgetbuf(lfn,ilist(1,ibuf(lfn)),maxlbuf,irtc)
            if(irtc .ne. 0)then
              return
            endif
            nc=min(nc,len(buff))
            lbuf(lfn)=nc
            mbuf(lfn)=0
          else
            nc=itfgetbuf(lfn,buff,len(buff),irtc)
            if(irtc .ne. 0)then
              return
            endif
            nc=min(nc,len(buff))
          endif
        else
          ls=ilist(1,ibuf(lfn)-1)
          if(jlist(lbuf(lfn)+1,ibuf(lfn)) .eq. 10)then
            lbuf(lfn)=lbuf(lfn)+1
          endif
          if(lbuf(lfn) .lt. ls)then
            do i=lbuf(lfn)+1,ls
              if(jlist(i,ibuf(lfn)) .eq. 10)then
                ie=i-1
                go to 10
              endif
            enddo
            ie=ls
 10         nc=ie-lbuf(lfn)
            mbuf(lfn)=lbuf(lfn)
            lbuf(lfn)=ie
          else
            go to 9000
          endif
        endif
      case (irbreadrecord)
        if(itbuf(lfn) .le. 2)then
          nc=itfgetbuf(lfn,buff,len(buff),irtc)
c          write(*,*)'tfreadbuf ',lfn,getfd(lfn),nc,irtc,buff(:nc)
          if(irtc .ne. 0)then
            return
          endif
c          read(lfn,'(q,a)',end=9000,err=9100)nc,
c     $         buff(:min(nc,len(buff)))
          nc=min(nc,len(buff))
          if(lfn .gt. 0)then
            lbuf(lfn)=0
            mbuf(lfn)=-1
          endif
        else
          ls=ilist(1,ibuf(lfn)-1)
          if(jlist(lbuf(lfn)+1,ibuf(lfn)) .eq. 10)then
            lbuf(lfn)=lbuf(lfn)+1
          endif
          if(lbuf(lfn) .lt. ls)then
            do i=lbuf(lfn)+1,ls
              if(jlist(i,ibuf(lfn)) .eq. 10)then
                ie=i-1
                go to 20
              endif
            enddo
            ie=ls
 20         nc=ie-lbuf(lfn)
            mbuf(lfn)=lbuf(lfn)
            lbuf(lfn)=ie
            call tmovb(jlist(mbuf(lfn)+1,ibuf(lfn)),buff,nc)
          else
            go to 9000
          endif
        endif
      case (irbreadbuf)
        if(lfn .gt. 0)then
          if(ibuf(lfn) .eq. 0 .or. mbuf(lfn) .lt. 0 .or.
     $         mbuf(lfn) .ge. lbuf(lfn))then
            nc=0
          else
            nc=lbuf(lfn)-mbuf(lfn)
            call tmovb(jlist(mbuf(lfn)+1,ibuf(lfn)),buff,nc)
          endif
        else
          go to 9100
        endif
      case (irbsetbuf)
        if(lfn .gt. 0)then
          if(itbuf(lfn) .le. 2)then
            call tmovb(buff,jlist(1,ibuf(lfn)),nc)
            lbuf(lfn)=nc
            mbuf(lfn)=0
          elseif(itbuf(lfn) .eq. 3)then
            mbuf(lfn)=max(mbuf(lfn),lbuf(lfn)-nc)
          endif
        endif
      case (irbcloseinp)
        if(lfn .gt. 0)then
          if(itbuf(lfn) .le. 2)then
            itbuf(lfn)=0
            if(ibuf(lfn) .gt. 0)then
              if(ilist(2,ibuf(lfn)-1) .ne. 0)then
                call unixclose(ilist(2,ibuf(lfn)-1))
                ilist(2,ibuf(lfn)-1)=0
              endif
            endif
          endif
        endif
      case (irbsetinp)
        if(lfn .gt. 0)then
          if(itbuf(lfn) .le. 2)then
            if(ibuf(lfn) .gt. 0)then
              ilist(2,ibuf(lfn)-1)=int(ib)
            endif
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

      subroutine readstr(in,str,irtc)
      implicit none
      integer*4 in,irtc,nc,itfgetbuf
      character*(*) str
      nc=itfgetbuf(in,str,len(str),irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(nc .gt. 0 .and. str(nc:nc) .eq. char(10))then
        str(nc:)=' '
      else
        str(nc+1:)=' '
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
c      write(*,*)'openshared ',ia,m
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'General::mmap','""')
        kx%k=kxfailed
        irtc=0
        return
      endif
      call tfreadbuf(irbopen,iu,ia,int8(4),nc,' ')
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
      integer*4 isp1,irtc,itfmessage,isp0,iu,nc
      logical*4 tfcheckelement
c      call tfdebugprint(ktastk(isp),'readshard',3)
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      irtc=0
      iu=int(rtastk(isp))
      call tfreadbuf(irbibuf,iu,ia,int8(4),nc,' ')
      if(ia .eq. 0)then
        kx%k=kxeof
        return
      endif
      do while(ilist(2,ia) .ne. 0)
c        write(*,*)'readshared ',ilist(2,ia)
        call tpause(10000)
      enddo
      ilist(2,ia)=1
      kx=dlist(ia+1)
c      write(*,*)'readshared ',ia,kx%k
      if(ktfobjqd(kx))then
c        write(*,*)'readshared-obj '
        if(ktfsymbolqd(kx))then
c          write(*,*)'readshared-symbol '
          if( .not. tfcheckelement(kx,.false.))then
            irtc=itfmessage(99,'General::wrongtype',
     $           '"undefined symbol returned in Shared"')
            kx%k=ktfoper+mtfnull
            ilist(2,ia)=0
            return
          endif
        elseif(ktfstringqd(kx,str))then
c          call tfdebugprint(kx,'readshared-string',1)
c          write(*,*)'at ',ia
          kx=kxsalocb(-1,str%str(1:str%nch),str%nch)
        elseif(ktflistqd(kx))then
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
      integer*4 isp1,irtc,itfmessage,itfmessageexp,isp0,n,i,iu,nc
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(ktfnonrealq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      iu=int(rtastk(isp1+1))
      call tfreadbuf(irbibuf,iu,kas,int8(4),nc,' ')
      if(kas .eq. 0)then
        irtc=itfmessage(99,'Shared::notopen','""')
        return
      endif
      do while(ilist(2,kas) .ne. 0)
c        write(*,*)'writeshared ',ilist(2,kas)
        call tpause(10000)
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
