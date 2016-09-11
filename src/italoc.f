      integer*4 function italoc(n)
      use tfstk
      use tfmem
      implicit none
      integer*4 n,n1,m,m1,na,l,l1
      integer*8 ktaloc,i,ic,i1,j,ixl(16)
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
      call forcesf()
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
      integer*8 ktaloc,k
      k=ktaloc(n)
      klist(k:k+n-1)=0
      ktcaloc=k
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
c      call tfsetlastp(ip+m-1)
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
c      if(tfchecklastp(ix1))then
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
c      endif
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

      subroutine talocinit
      use tfstk
      use tfmem
      use iso_c_binding
      implicit none
      integer*8 ka,ic
      integer*8, parameter :: kcpthre1=(ktamask+1)/2,
     $     nitalocmin=2**20
      allocate(sadalloc(1)%ca(nindex*2+mhash+16))
      ka=transfer(c_loc(sadalloc(1)%ca(1)),int8(0))
      kfirstalloc=ka
c     write(*,*)'talocinit ',ka,kcpthre,kcpthre1
c$$$      if(ka .ge. kcpthre)then
c$$$        deallocate(sadalloc(1)%ca)
c$$$c        if(ka .ge. kcpthre1)then
c$$$c          kcpklist0=kcpthre1
c$$$c        else
c$$$          kcpklist0=kcpthre
c$$$c        endif
c$$$        call tfcbkinit
c$$$        irtc=-1
c$$$        nitaloc=nitalocmin*8
c$$$        do while(irtc .ne. 0)
c$$$          nitaloc=nitaloc/2
c$$$          if(nitaloc .lt. nitalocmin)then
c$$$            write(*,*)'talocinit-error ',kcpklist0
c$$$            call forcesf()
c$$$          endif
c$$$          mm=mapallocfixed8(klist(0),nitaloc,8,irtc)
c$$$        enddo
c$$$        write(*,*)'talocinit-1 ',mm,nitaloc
c$$$        icbk=1
c$$$        jcbk=1
c$$$        kcbk(1,1)=0
c$$$        kcbk(2,1)=nitaloc-1
c$$$        kcbk(3,1)=0
c$$$        icp=1
c$$$      else
        kcpklist0=0
        call tfcbkinit
        icp=ksad_loc(sadalloc(1)%ca(1))
        icbk=1
        jcbk=1
        kcbk(1,1)=icp
        kcbk(2,1)=icp+nindex*2+mhash+15
        kcbk(3,1)=kcbk(2,1)
c$$$      endif
c      write(*,*)'talocinit ',icp,kcpklist0
      do ic=icp,icp+nindex*2,2
        klist(ic)=ic
        klist(ic+1)=ic
      enddo
      ich=icp+nindex*2+4
      do ic=ich,ich+mhash
        klist(ic)=ic
      enddo
      icsep=ich+mhash+1
      klist(icsep)=icsep
c$$$      if(kcpklist0 .ne. 0)then
c$$$        kb=icsep+2
c$$$        nmem=int(nitaloc-1-(icsep+1-icp))
c$$$        nnet=nmem
c$$$        ilist(1,kb-1)=int(nmem)
c$$$c      write(*,*)'talocinit ',ilist(1,kb-1)+kb-1
c$$$        call tfree(kb)
c$$$c     ics=icp+nitaloc-3
c$$$c     klist(ics)=klist(icsep)
c$$$c     klist(icsep)=ics
c$$$c     ilist(1,ics-1)=4
c$$$c     ilist(2,ics-1)=4
c$$$      else
        nmem=0
        nnet=0
c$$$      endif
      return
      end subroutine

c$$$
c$$$      subroutine talocinit
c$$$      use tfstk
c$$$      use tfmem
c$$$      use iso_c_binding
c$$$      implicit none
c$$$      integer*8 ka,kb,ic
c$$$      integer*4 irtc
c$$$      integer*8, parameter ::
c$$$     $     kcpthre1=(ktamask+1)/2,nitalocmin=2**20,katry0=2**23
c$$$      integer*8 mapallocfixed8,mm
c$$$      integer*8, pointer, dimension(:) :: kla
c$$$      type (c_ptr) cp
c$$$      mm=0
c$$$      call c_f_pointer(transfer(int8(8),cp),kla,[16])
c$$$      irtc=1
c$$$      nitaloc=nitalocmin
c$$$      ka=katry0
c$$$      do while(irtc .ne. 0)
c$$$c        write(*,*)'talocinit ',ka,nitaloc
c$$$        mm= mapallocfixed8(kla(ka),nitaloc,8,irtc)
c$$$c        write(*,*)'talocinit-mmap ',irtc,mm
c$$$        ka=ka+nitaloc
c$$$        if(ka .gt. ktamask)then
c$$$          write(*,*)'Initial memory allocation failed. ',ka
c$$$          call forcesf()
c$$$        endif
c$$$      enddo
c$$$      write(*,*)'talocinit ',ka,mm
c$$$      kcpklist0=mm-8
c$$$      call tfcbkinit
c$$$      icp=1
c$$$      icbk=0
c$$$c     allocate(sadalloc(1)%ca(nindex*2+mhash+16))
c$$$c     ka=transfer(c_loc(sadalloc(1)%ca(1)),int8(0))
c$$$c     kfirstalloc=ka
c$$$c     write(*,*)'talocinit ',ka,kcpthre,kcpthre1
c$$$c$$$  if(ka .ge. kcpthre)then
c$$$c$$$  deallocate(sadalloc(1)%ca)
c$$$c$$$  if(ka .ge. kcpthre1)then
c$$$c$$$  kcpklist0=kcpthre1
c$$$c$$$  else
c$$$c$$$  kcpklist0=kcpthre
c$$$c$$$  endif
c$$$c$$$  call tfcbkinit
c$$$c$$$  irtc=-1
c$$$c$$$  nitaloc=nitaloc0*2
c$$$c$$$  do while(irtc .ne. 0)
c$$$c$$$  nitaloc=nitaloc/2
c$$$c$$$  if(nitaloc .lt. nitalocmin)then
c$$$c$$$  write(*,*)'talocinit-error ',kcpklist0
c$$$c$$$  call forcesf()
c$$$c$$$  endif
c$$$c$$$  mm= mapallocfixed8(klist(0),nitaloc,8,irtc)
c$$$c$$$  enddo
c$$$c$$$  c          write(*,*)'talocinit-1 ',nitaloc
c$$$c$$$  icbk=0
c$$$c$$$  icp=1
c$$$c$$$  else
c$$$c$$$  kcpklist0=0
c$$$c$$$  call tfcbkinit
c$$$c$$$  icbk=1;
c$$$c$$$  icp=ksad_loc(sadalloc(1)%ca(1))
c$$$c$$$  endif
c$$$c     write(*,*)'talocinit ',icp,kcpklist0
c$$$      do ic=icp,icp+nindex*2,2
c$$$        klist(ic)=ic
c$$$        klist(ic+1)=ic
c$$$      enddo
c$$$      ich=icp+nindex*2+4
c$$$      do ic=ich,ich+mhash
c$$$        klist(ic)=ic
c$$$      enddo
c$$$      icsep=ich+mhash+1
c$$$      klist(icsep)=icsep
c$$$c     if(ka .ge. kcpthre)then
c$$$      kb=icsep+2
c$$$      nmem=int(nitaloc-(icsep+1-icp))
c$$$      nnet=nmem
c$$$      ilist(1,kb-1)=int(nmem)
c$$$c      write(*,*)'talocinit ',nmem,icsep,icp,kb,kb-1+nmem-1
c$$$c      write(*,*)klist(kb-1+nmem-1)
c$$$      call tfree(kb)
c$$$c     ics=icp+nitaloc-3
c$$$c     klist(ics)=klist(icsep)
c$$$c     klist(icsep)=ics
c$$$c     ilist(1,ics-1)=4
c$$$c     ilist(2,ics-1)=4
c$$$c     endif
c$$$      return
c$$$      end subroutine
c$$$
