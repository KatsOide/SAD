c$$$      integer*4 function italoc(n)
c$$$      use tfstk
c$$$      use tfmem
c$$$      implicit none
c$$$      integer*4 n,n1,m,m1,na,l,l1
c$$$      integer*8 i,ic,i1,j,ixl(16)
c$$$      do l=1,16
c$$$        ixl(l)=ktaloc(n)
c$$$        if(ixl(l) .ge. 2**31-n .or. ixl(l) .le. 0)then
c$$$          if(l .ne. 1)then
c$$$            cycle
c$$$          else
c$$$            n1=max(n,3)
c$$$            m=n1+1
c$$$            do na=n1,nindex
c$$$              ic=icp+na*2
c$$$              i1=ic
c$$$              i=klist(i1)
c$$$              do while(i .ne. ic)
c$$$                if(i .gt. 0 .and. i .lt. 2**31)then
c$$$                  m1=ilist(1,i-1)
c$$$                  if(m1 .eq. m)then
c$$$                    klist(i1)=klist(i)
c$$$                    klist(klist(i)+1)=i1
c$$$                    j=ich+iand(i+m+2,mhash)
c$$$                    do while(klist(j) .ne. i+2)
c$$$                      j=klist(j)
c$$$                    enddo
c$$$                    klist(j)=klist(i+2)
c$$$                    nnet=nnet+m
c$$$                    italoc=int(i)
c$$$                    if(ic .eq. maxic .and. klist(ic) .eq. ic)then
c$$$                      maxic=maxic-2
c$$$                    endif
c$$$                    call tfree(ixl(1))
c$$$                    return
c$$$                  elseif(m1 .ge. m+minseg2)then
c$$$                    klist(i1)=klist(i)
c$$$                    klist(klist(i)+1)=i1
c$$$                    j=ich+iand(i+m1+2,mhash)
c$$$                    do while(klist(j) .ne. i+2)
c$$$                      j=klist(j)
c$$$                    enddo
c$$$                    klist(j)=klist(i+2)
c$$$                    call tsetindexhash(i+m,m1-m)
c$$$                    ilist(1,i-1)=m
c$$$                    nnet=nnet+m
c$$$                    italoc=int(i)
c$$$                    if(ic .eq. maxic .and. klist(ic) .eq. ic)then
c$$$                      maxic=maxic-2
c$$$                    endif
c$$$                    call tfree(ixl(1))
c$$$                    return
c$$$                  endif
c$$$                endif
c$$$                i1=i
c$$$                i=klist(i)
c$$$              enddo
c$$$            enddo
c$$$          endif
c$$$        else
c$$$          do l1=l-1,1,-1
c$$$            call tfree(ixl(l))
c$$$          enddo
c$$$          italoc=int(ixl(l))
c$$$          return
c$$$        endif
c$$$      enddo
c$$$      write(*,*)'italoc memory allocation error ',ixl(16),n
c$$$      call meminfo()
c$$$      call abort
c$$$      italoc=-1
c$$$      return
c$$$      end
c$$$
c$$$      integer*4 function itcaloc(n)
c$$$      use tfstk
c$$$      implicit none
c$$$      integer*4 n,italoc,i
c$$$      i=italoc(n)
c$$$      klist(i:i+n-1)=0
c$$$      itcaloc=i
c$$$      return
c$$$      end
c$$$
      integer*8 function ktcaloc(n)
      use tfstk
      implicit none
      integer*4 ,intent(in):: n
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
      integer*4 ,intent(in):: m
      integer*4 np
      integer*8 ,intent(out):: ip1
      integer*8 j,j1,ics,na
      np=int(1+(m+3)/mpsize)
      na=mpsize*np
      ip1=ktfsadalloc(na)+1
c      write(*,*)'talocp ',ip1,ip1+na
      if(ip1 .le. 0)then
        return
      endif
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

c$$$      integer*4 function mtaloc(n)
c$$$      use tfstk
c$$$      implicit none
c$$$      integer*4 n,italoc
c$$$      mtaloc=italoc(max(3,n-1))-1
c$$$      return
c$$$      end
c$$$
c$$$      integer*4 function mctaloc(n)
c$$$      use tfstk
c$$$      implicit none
c$$$      integer*4 n,i,italoc
c$$$      i=italoc(max(3,n-1))-1
c$$$      klist(i:i+n-1)=0
c$$$      mctaloc=i
c$$$      return
c$$$      end

      logical*4 function tfaloc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 ,intent(in):: n
      integer*4 n1
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
      type (sad_descriptor) ,intent(out):: kx
      integer*8 k
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage,m,ir
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
      character*(*) ,intent(in):: tag
      write(*,*)'aloc-',tag,nnet,nmem
      return
      end

      subroutine tfreem(ip,n)
      use tfstk
      implicit none
      integer*8 ,intent(in):: ip
      integer*4 ,intent(in):: n
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
