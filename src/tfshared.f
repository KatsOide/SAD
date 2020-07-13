      module tfshare
      integer*4, parameter:: nshmax=1024
      integer*8, save:: kashare(nshmax)=0,lshare(nshmax)=0
      integer*4, save :: ishared(nshmax)=0,kstshare(0:nshmax)=0
      
      contains
      integer*8 function ktfallocshared(n)
      use tfstk
      use iso_c_binding
      implicit none
      integer*4, save :: lps=0
      integer*4 irtc,n,getpagesize,nsh1,i,na
      integer*8 k,kpb,kcp
      if(lps .eq. 0)then
        lps=getpagesize()/8
      endif
      na=((n+3)/lps+2)*lps
      nsh1=0
      do i=1,kstshare(0)
        if(kstshare(i) .eq. 0 .and. lshare(i) .ge. na)then
          kstshare(i)=1
          ktfallocshared=kashare(i)
          return
        elseif(kstshare(i) .ge. 0)then
          nsh1=i
        endif
      enddo
      kstshare(0)=nsh1
      nsh1=0
      if(kstshare(0) .eq. nshmax)then
        do i=1,kstshare(0)
          if(kstshare(i) .le. 0)then
            if(lshare(i) .gt. 0)then
              kstshare(i)=1
              if(nsh1 .ne. 0)then
                kstshare(nsh1)=1-kstshare(i)
              endif
              call tfreleaseshared(kashare(i))
              nsh1=i
              exit
            elseif(nsh1 .eq. 0)then
              kstshare(i)=1-kstshare(i)
              nsh1=i
            endif
          endif
        enddo
      else
        nsh1=kstshare(0)+1
        kstshare(nsh1)=1
        kstshare(0)=nsh1
      endif
      k=ktaloc(na)
      kcp=transfer(c_loc(klist(k)),k)/8
      kpb=k+((kcp+1+lps)/lps)*lps-kcp
c      write(*,*)'ktfallocshared-mmap ',nsh1,kpb,na-lps
      call mapallocshared8(klist(kpb),int8(na-lps),8,irtc)
      if(irtc .ne. 0)then
        write(*,*)'ktfallocshared failed: ',kpb,na-lps
        call abort
      endif
      ktfallocshared=kpb+2
      klist(kpb)=k
      klist(kpb+1)=na-lps
      if(nsh1 .ne. 0)then
        ishared(nsh1)=1
        kstshare(nsh1)=1
        kashare(nsh1)=kpb+2
        lshare(nsh1)=na
      endif
c      write(*,*)'allocshared ',kpb,na,lps,
c     $     transfer(c_loc(klist(kpb)),k)/8
      return
      end function

      subroutine tfreeshared(kpb,ist)
      use tfstk
      implicit none
      integer*4, optional :: ist
      integer*4 i,is
      integer*8 kpb
      is=0
      if(present(ist))then
        is=ist
      endif
      do i=1,kstshare(0)
        if(kashare(i) .eq. kpb)then
          kstshare(i)=is
          if(is .lt. 0)then
            lshare(i)=0
          endif
          return
        endif
      enddo
      call tfreleaseshared(kpb)
      return
      end subroutine

      subroutine tfreleaseshared(kpb)
      use tfstk
      implicit none
      integer*8 kpb,k
      integer*4 irtc
      k=klist(kpb-2)
      call mapallocfixed8(klist(kpb-2),klist(kpb-1),8,irtc)
      if(irtc .ne. 0)then
        write(*,*)'tffreecshared failed: ',kpb,klist(kpb-1)
        call abort
      endif
c      write(*,*)'tfreeshared ',kpb,klist(kpb-1),irtc
      if(itfcbk(k) .eq. 0)then
        call tfentercbk(k,klist(kpb-1)+1)
      endif
      call tfree(k)
      return
      end subroutine

      end module

      subroutine tfsavesharedmap()
      use tfstk
      use tfshare
      implicit none
      ishared=0
      kstshare(0)=0
      return
      end subroutine

      subroutine tfresetsharedmap()
      use tfstk
      use tfshare
      implicit none
      integer*4 i
      do i=1,kstshare(0)
        if(ishared(i) .ne. 0)then
          kstshare(i)=-1
        endif
      enddo
      return
      end subroutine

      subroutine ktfinitshare
      use tfstk
      use iso_c_binding
      use tfshare
      implicit none
      integer*4, save :: lps=0
      integer*4 getpagesize
      if(lps .eq. 0)then
        lps=getpagesize()/8
      endif
      return
      end

      recursive subroutine tfrecallshared(isp0,k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,k0,ki
      type (sad_string), pointer :: str
      type (sad_dlist), pointer :: kl
      integer*8 kax
      integer*4 isp0,m,irtc,itfmessage,i
      logical*4 tfcheckelement
      do i=isp0+1,isp
        if(ktastk(i) .eq. k%k)then
          kx=dtastk2(i)
          return
        endif
      enddo
      if(ktfstringq(k,str))then
        kx=kxsalocb(-1,str%str,str%nch)
      elseif(ktfsymbolq(k))then
        if( .not. tfcheckelement(k,.false.))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"defined symbol returned in Shared"')
          return
        endif
        kx=k
      elseif(ktflistq(k,kl))then
c        call tfdebugprint(k,'recallshared',3)
        k0=kl%head
        if(ktfobjq(k0))then
          call tfrecallshared(isp0,k0,k0,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
        m=kl%nl
        if(kl%ref .eq. 0)then
          kax=ktavaloc(-1,m)
          dlist(kax+1:kax+m)=kl%dbody(1:m)
        else
          kax=ktadaloc(-1,m)
          do i=1,m
            ki=kl%dbody(i)
            if(ktfobjq(ki))then
              call tfrecallshared(isp0,ki,ki,irtc)
              if(irtc .ne. 0)then
                klist(kax+1:kax+m)=ktfoper+mtfnull
                exit
              endif
              dlist(kax+i)=dtfcopy1(ki)
            else
              dlist(kax+i)=ki
            endif
          enddo
        endif
        dlist(kax)=dtfcopy(k0)
        kx%k=ktflist+kax
      else
        kx%k=ktfoper+mtfnull
      endif
      isp=isp+1
      dtastk(isp)=k
      dtastk2(isp)=kx
      irtc=0
      return
      end

      subroutine tfstoreshared(isp0,k,kap)
      use tfstk
      implicit none
      integer*8 k,ka,kap,kh,ki,kt
      integer*4 isp0,i,j,m
      ka=ktfaddr(k)
      kt=k-ka
      if(kt .eq. ktfstring)then
        klist(kap)=klist(ka)
        do i=1,(ilist(1,ka)+7)/8
          rlist(kap+i)=rlist(ka+i)
        enddo
      elseif(kt .eq. ktflist)then
        kh=klist(ka)
        klist(kap+1)=kh
        if(ktfobjq(kh) .and. ktfnonsymbolq(kh))then
c          call tfdebugprint(kh,'storeshared-head',1)
          do j=isp0,isp
            if(ktastk(j) .eq. kh)then
              if(ktfstringq(kh))then
                klist(kap+1)=ktfstring+ktastk2(j)
              else
                klist(kap+1)=ktflist+ktastk2(j)+1
              endif
c              klist(kap+1)=iand(ktfmask,kh)+ktastk2(j)
              exit
            endif
          enddo
        endif
        m=ilist(2,ka-1)
        ilist(2,kap)=m
        if(ktfreallistq(ka))then
          ilist(1,kap)=0
          rlist(kap+2:kap+m+1)=rlist(ka+1:ka+m)
c          do i=1,m
c            rlist(kap+i+1)=rlist(ka+i)
c          enddo
        else
          ilist(1,kap)=1
          do i=1,m
            ki=klist(ka+i)
            if(ktfnonobjq(ki) .or. ktfsymbolq(ki))then
              klist(kap+i+1)=ki
            else
              do j=isp0,isp
                if(ktastk(j) .eq. ki)then
                  if(ktfstringq(ki))then
                    klist(kap+i+1)=ktfstring+ktastk2(j)
                  else
                    klist(kap+i+1)=ktflist+ktastk2(j)+1
                  endif
                  exit
                endif
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      recursive subroutine tfsharedsize(isp0,k,n,irtc)
      use tfstk
      implicit none
      integer*8 k,ka,kt,ki,kh
      integer*4 isp0,i,irtc,itfmessage,ni,n
      irtc=0
      if(ktfnonobjq(k) .or. ktfsymbolq(k))then
        n=1
      else
        ka=ktfaddr(k)
        kt=k-ka
        if(kt .eq. ktfstring)then
          n=3+(ilist(1,ka)+7)/8
          do i=isp0+1,isp
            if(ktastk(i) .eq. k)then
              n=1
              exit
            endif
          enddo
          isp=isp+1
          ktastk(isp)=k
          itastk2(1,isp)=n
        elseif(kt .eq. ktflist)then
          do i=isp0+1,isp
            if(ktastk(i) .eq. k)then
              n=1
              return
            endif
          enddo
          n=3+ilist(2,ka-1)
          isp=isp+1
          ktastk(isp)=k
          itastk2(1,isp)=n
          kh=klist(ka)
          if(ktfobjq(kh) .and. ktfnonsymbolq(kh))then
            call tfsharedsize(isp0,klist(ka),ni,irtc)
            if(irtc .ne. 0)then
              return
            endif
            n=n+ni-1
          endif
          if(ktfnonreallistq(ka))then
            do i=1,ilist(2,ka-1)
              ki=klist(ka+i)
              if(ktfobjq(ki) .and. ktfnonsymbolq(ki))then
                call tfsharedsize(isp0,ki,ni,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                n=n+ni-1
              endif
            enddo
          endif
        else
          go to 9000
        endif
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"No Symbols nor Patterns"')
      return
      end
