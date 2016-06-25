      integer*4 function italoc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 n,n1,m,m1,na
      integer*8 ix,ktaloc,i,ic,i1,j
      ix=ktaloc(n)
      if(ix .ge. 2**31 .or. ix .le. 0)then
c        write(*,*)'italoc ',ix
        call tfree(ix)
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
                return
              endif
            endif
            i1=i
            i=klist(i)
          enddo
        enddo
        write(*,*)'italoc memory allocation error ',ix,n
        call meminfo()
        call forcesf()
      endif
      italoc=int(ix)
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

      integer*8 function ktaloc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 n,m,n1,m1
      integer*8 ic1,i,ic,ic2,ip1,j,i1
      n1=max(n,3)
      m=n1+1
      if(n1 .lt. nindex)then
        ic=icp+n1*2
        i=klist(ic)
        if(i .ne. ic)then
          m=ilist(1,i-1)
          klist(ic)=klist(i)
          klist(klist(i)+1)=ic
          j=ich+iand(i+m+2,mhash)
          do while(klist(j) .ne. i+2)
            j=klist(j)
          enddo
          klist(j)=klist(i+2)
          nnet=nnet+m
          ktaloc=i
          return
        endif
        ic1=ic+min(m,minseg1)*2
        if(ic1 .le. maxic)then
          ic2=min(maxic-2,ic1+10)
          do ic=ic1,ic2,2
            i=klist(ic)
            if(ic .ne. i)then
              m1=ilist(1,i-1)
              klist(ic)=klist(i)
              klist(klist(i)+1)=ic
              j=ich+iand(i+m1+2,mhash)
              do while(klist(j) .ne. i+2)
                j=klist(j)
              enddo
              klist(j)=klist(i+2)
              call tsetindexhash(i+m,m1-m)
              ilist(1,i-1)=m
              nnet=nnet+m
              ktaloc=i
              return
            endif
          enddo
          do ic=maxic,ic2+2,-2
            i=klist(ic)
            if(ic .ne. i)then
              m1=ilist(1,i-1)
              klist(ic)=klist(i)
              klist(klist(i)+1)=ic
              j=ich+iand(i+m1+2,mhash)
              do while(klist(j) .ne. i+2)
                j=klist(j)
              enddo
              klist(j)=klist(i+2)
              call tsetindexhash(i+m,m1-m)
              ilist(1,i-1)=m
              nnet=nnet+m
              if(klist(ic) .eq. ic)then
                maxic=ic-2
              else
                maxic=ic
              endif
              ktaloc=i
              return
            endif
          enddo
          maxic=ic1-2
        endif
      endif
      ic=icp+nindex*2
 1000 i1=ic
      i=klist(i1)
      do while(i .ne. ic)
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
          ktaloc=i
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
          ktaloc=i
          return
        endif
        i1=i
        i=klist(i)
      enddo
      call talocp(m,ip1)
      if(ip1 .gt. 0)then
        go to 1000
      endif
      ktaloc=-1
      return
      end

      subroutine tsetindexhash(ip,m)
      use tfstk
      use tfmem
      implicit none
      integer*4 m
      integer*8 ip,ic,ic1,ia
      if(m .gt. nindex)then
        ic=icp+nindex*2
      else
        ic=icp+(m-1)*2
        maxic=max(ic,maxic)
      endif
      klist(ip  )=klist(ic)
      klist(ip+1)=ic
      klist(klist(ic)+1)=ip
      klist(ic)=ip
      ia=ip+2
      ic1=ich+iand(ia+m,mhash)
      klist(ia  )=klist(ic1)
      klist(ic1 )=ia
      ilist(1,ip-1)=m
      lastpend=max(lastpend,ip)
      return
      end

      subroutine tfree(ka)
      use tfstk
      use tfmem
      implicit none
      integer*8 ka,ix,ik,ik0,ip,ix1
      integer*4 m,mx
      m=ilist(1,ka-1)
      if(m .lt. 4)then
        if(m .ne. 0)then
          write(*,*)'tfree-too small segment: ',ka,m
          call forcesf()
        endif
        return
      endif
      nnet=nnet-m
      ix=ka+2
      ik=iand(ix,mhash)+ich
      ik0=ik
      ip=klist(ik)
      do while(ip .ne. ik0)
        if(ip .lt. ix)then
          mx=ilist(1,ip-3)
          if(ip+mx .eq. ix)then
            klist(ik  )=klist(ip)
            klist(klist(ip-2)+1)=klist(ip-1)
            klist(klist(ip-1)  )=klist(ip-2)
            m=m+mx
            ix=ip
            exit
          endif
        endif
        ik=ip
        ip=klist(ik)
      enddo
      ix1=ix+m
      if(ix1 .le. lastpend)then
        ik=iand(ix1+ilist(1,ix1-3),mhash)+ich
        ik0=ik
        ip=klist(ik)
        do while(ip .ne. ik0)
          if(ip .eq. ix1)then
            klist(ik  )=klist(ip)
            klist(klist(ip-2)+1)=klist(ip-1)
            klist(klist(ip-1)  )=klist(ip-2)
            m=m+ilist(1,ix1-3)
            exit
          endif
          ik=ip
          ip=klist(ik)
        enddo
      endif
      call tsetindexhash(ix-2,m)
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
      ilist(1,ics-1)=4
      ilist(2,ics-1)=4
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
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
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
      integer*8 ktaloc,k
      integer*4 isp1,irtc,i,itfmessage,m,ir
      kx%k=ktfoper+mtfnull
      do i=isp1+1,isp
        if(ktfnonrealqdi(dtastk(i),ir))then
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
      integer*4 ip,n
      ilist(1,ip)=max(n,4)
      call tfree(int8(ip+1))
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
