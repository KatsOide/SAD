      integer*4 function italoc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 n,n1,m,m1,na,l,l1
      integer*8 i,ic,i1,j,ixl(16)
      do l=1,16
        ixl(l)=ktaloc(n)
        if(ixl(l) .ge. 2**31-n .or. ixl(l) .le. 0)then
          if(l .ne. 1)then
            cycle
          else
            n1=max(n,3)
            m=n1+1
            do na=n1,nindex
              ic=icp+na*2
              i1=ic
              i=klist(i1)
              do while(i .ne. ic)
                if(i .gt. 0 .and. i .lt. 2**31)then
                  m1=ilist(1,i-1)
                  if(m1 .eq. m)then
                    klist(i1)=klist(i)
                    klist(klist(i)+1)=i1
                    j=ich+iand(i+m+2,mhash)
                    do while(klist(j) .ne. i+2)
                      j=klist(j)
                    enddo
                    klist(j)=klist(i+2)
                    nnet=nnet+m
                    italoc=int(i)
                    if(ic .eq. maxic .and. klist(ic) .eq. ic)then
                      maxic=maxic-2
                    endif
                    call tfree(ixl(1))
                    return
                  elseif(m1 .ge. m+minseg2)then
                    klist(i1)=klist(i)
                    klist(klist(i)+1)=i1
                    j=ich+iand(i+m1+2,mhash)
                    do while(klist(j) .ne. i+2)
                      j=klist(j)
                    enddo
                    klist(j)=klist(i+2)
                    call tsetindexhash(i+m,m1-m)
                    ilist(1,i-1)=m
                    nnet=nnet+m
                    italoc=int(i)
                    if(ic .eq. maxic .and. klist(ic) .eq. ic)then
                      maxic=maxic-2
                    endif
                    call tfree(ixl(1))
                    return
                  endif
                endif
                i1=i
                i=klist(i)
              enddo
            enddo
          endif
        else
          do l1=l-1,1,-1
            call tfree(ixl(l))
          enddo
          italoc=int(ixl(l))
          return
        endif
      enddo
      write(*,*)'italoc memory allocation error ',ixl(16),n
      call meminfo()
      call abort
      italoc=-1
      return
      end

      integer*4 function itcaloc(n)
      use tfstk
      implicit none
      integer*4 n,italoc,i
      i=italoc(n)
      klist(i:i+n-1)=0
      itcaloc=i
      return
      end

      integer*8 function ktcaloc(n)
      use tfstk
      implicit none
      integer*4 n
      integer*8 k
      k=ktaloc(n)
      klist(k:k+n-1)=0
      ktcaloc=k
      return
      end

      subroutine talocp(m,ip1)
      use tfstk
      use tfmem
      implicit none
      integer*4 np,m
      integer*8 ip1,j,j1,ics,na
      np=int(1+(m+3)/mpsize)
      na=mpsize*np
      ip1=ktfsadalloc(na)+1
c      write(*,*)'talocp ',ip1,ip1+na
      nmem=nmem+na
      nnet=nnet+na
      j=klist(icsep)
      j1=icsep
      do while(j .ne. icsep)
        if(j+4 .eq. ip1)then
          klist(j1)=klist(j)
          ip1=j
          na=na+4
          exit
        endif
        ilist(2,j-1)=ilist(2,j-1)-1
        if(ilist(2,j-1) .le. 0)then
          klist(j1)=klist(j)
        endif
        j1=j
        j=klist(j)
      enddo
      ilist(1,ip1-1)=int(na)-4
      call tfree(ip1)
      ics=ip1+na-4
      klist(ics)=klist(icsep)
      klist(icsep)=ics
      ilist(:,ics-1)=4
      return
      end

      integer*4 function mtaloc(n)
      use tfstk
      implicit none
      integer*4 n,italoc
      mtaloc=italoc(max(3,n-1))-1
      return
      end

      integer*4 function mctaloc(n)
      use tfstk
      implicit none
      integer*4 n,i,italoc
      i=italoc(max(3,n-1))-1
      klist(i:i+n-1)=0
      mctaloc=i
      return
      end

      logical*4 function tfaloc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 n,n1
      n1=max(n,3)
      if(n1 .lt. nindex)then
        tfaloc=klist(icp+n1*2) .ne. icp+n1*2
      else
        tfaloc=.false.
      endif
      return
      end

      subroutine tfmalloc(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*8 k
      integer*4 isp1,irtc,i,itfmessage,m,ir
      kx%k=ktfoper+mtfnull
      do i=isp1+1,isp
        if(ktfnonrealq(dtastk(i),ir))then
          irtc=itfmessage(9,'General::wrongtype','"Real"')
          return
        endif
        m=max(3,(ir+7)/8)
        k=ktaloc(m)
        if(k .le. 0)then
          irtc=itfmessage(9999,'Memory::alloc','""')
          return
        endif
        call tfree(k)
      enddo
      irtc=0
      return
      end

      subroutine tprintaloc(tag)
      use tfstk
      use tfmem
      implicit none
      character*(*) tag
      write(*,*)'aloc-',tag,nnet,nmem
      return
      end

      subroutine tfreem(ip,n)
      use tfstk
      implicit none
      integer*8 ip
      integer*4 n
      ilist(1,ip)=max(n,4)
      call tfree(ip+1)
      return
      end

      subroutine meminfo()
      use tfmem
      implicit none
      write(*,*)' klist(0): ',kcpklist0,
     $     ', nitaloc: ',nitaloc,', icp: ',icp,
     $     ', firstalloc: ',kfirstalloc
      return
      end
