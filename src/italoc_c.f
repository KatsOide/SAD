c$Header: /SAD/cvsroot/oldsad/src/italoc.f,v 1.59.2.7 2012/08/21 15:47:05 oide Exp $
      integer*4 function italoc(n)
      implicit none
      include 'inc/TFCBK.inc'
      integer*4 n
      integer*8 ktaloc,ix
      ix=ktaloc(n)
      if(ix .ge. 2**31 .or. ix .le. 0)then
        write(*,*)'italoc memory allocation error ',ix
        stop
      endif
      italoc=ix
      return
      end

      integer*8 function ktaloc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,m,n1,m1
      integer*8 ic1,i,ic,ip0,ic2,ip1
c      include 'DEBUG.inc'
      n1=max(n,3)
      m=n1+1
      if(n1 .lt. nindex)then
        i=icp+n1*2
        ic=klist(i)
        if(ic .ne. i)then
          klist(i)=klist(ic)
          klist(klist(ic)+1)=i
          klist(klist(ic+3))  =klist(ic+2)
          klist(klist(ic+2)+1)=klist(ic+3)
          if(m .gt. 6)then
            klist(klist(ic+5))  =klist(ic+4)
            klist(klist(ic+4)+1)=klist(ic+5)
          endif
          ilist(1,ic)=m
          nnet=nnet+m
          ktaloc=ic+1
          return
        endif
        ic1=i+min(m,minseg1)*2
        if(ic1 .le. maxic)then
          ic2=min(maxic-2,ic1+10)
          do i=ic1,ic2,2
            ic=klist(i)
            if(ic .ne. i)then
              klist(i)=klist(ic)
              klist(klist(ic)+1)=i
              klist(klist(ic+3))  =klist(ic+2)
              klist(klist(ic+2)+1)=klist(ic+3)
              klist(klist(ic+5))  =klist(ic+4)
              klist(klist(ic+4)+1)=klist(ic+5)
              m1=klist(ic+6)
              call tsetindexhash(ic+m,m1-m)
              ilist(1,ic)=m
              nnet=nnet+m
              ktaloc=ic+1
              return
            endif
          enddo
          do i=maxic,ic2+2,-2
            ic=klist(i)
            if(ic .ne. i)then
              klist(i)=klist(ic)
              klist(klist(ic)+1)=i
              klist(klist(ic+3))  =klist(ic+2)
              klist(klist(ic+2)+1)=klist(ic+3)
              klist(klist(ic+5))  =klist(ic+4)
              klist(klist(ic+4)+1)=klist(ic+5)
              m1=klist(ic+6)
              call tsetindexhash(ic+m,m1-m)
              ilist(1,ic)=m
              nnet=nnet+m
              if(klist(i) .eq. ic)then
                maxic=i-2
              else
                maxic=i
              endif
              ktaloc=ic+1
              return
            endif
          enddo
          maxic=ic1-2
        endif
      endif
      ip0=icp+nindex*2
 1000 ic=klist(ip0)
      do while(ic .ne. ip0)
        m1=klist(ic+6)
        if(m1 .eq. m)then
          klist(klist(ic+1))  =klist(ic)
          klist(klist(ic  )+1)=klist(ic+1)
          klist(klist(ic+3))  =klist(ic+2)
          klist(klist(ic+2)+1)=klist(ic+3)
          klist(klist(ic+5))  =klist(ic+4)
          klist(klist(ic+4)+1)=klist(ic+5)
          ilist(1,ic)=m
          nnet=nnet+m
          ktaloc=ic+1
          return
        elseif(m1 .ge. m+minseg2)then
          klist(klist(ic+1))  =klist(ic)
          klist(klist(ic  )+1)=klist(ic+1)
          klist(klist(ic+3))  =klist(ic+2)
          klist(klist(ic+2)+1)=klist(ic+3)
          klist(klist(ic+5))  =klist(ic+4)
          klist(klist(ic+4)+1)=klist(ic+5)
          call tsetindexhash(ic+m,m1-m)
          ilist(1,ic)=m
          nnet=nnet+m
          ktaloc=ic+1
          return
        endif
        ic=klist(ic)
      enddo
      call talocp(m,ip1)
      if(ip1 .gt. 0)then
        go to 1000
      endif
      ktaloc=-1
      return
      end

      integer*4 function mctaloc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      integer*4 n,i,j,italoc
      i=italoc(n-1)-1
      do j=i,i+n-1
        rlist(j)=0.d0
      enddo
      mctaloc=i
      return
      end

      subroutine tsetindexhash(ip,m)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
c      include 'DEBUG.inc'
      integer*4 m
      integer*8 ip,ic,ic1,ia
      ia=ip
      if(m .gt. nindex)then
        ic=icp+nindex*2
      else
        ic=icp+(m-1)*2
        maxic=max(ic,maxic)
      endif
      klist(ia  )=klist(ic)
      klist(ia+1)=ic
      klist(klist(ic)+1)=ia
      klist(ic)=ia
      ia=ia+2
      ic1=ich+iand(ia+m,mhash)*4
      klist(ia  )=klist(ic1)
      klist(ia+1)=ic1
      klist(klist(ic1)+1)=ia
      klist(ic1 )=ia
      if(m .gt. 6)then
        ia=ia+2
        ic1=ich+iand(ia,mhash)*4+2
        klist(ia  )=klist(ic1)
        klist(ia+1)=ic1
        klist(klist(ic1)+1)=ia
        klist(ic1 )=ia
        klist(ia+2)=m
      endif
      return
      end

      subroutine tfree(ka)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
c      include 'DEBUG.inc'
      integer*8 ka,ic,ix,ik,ik0,ip,ip0,ix1,ix2
      integer*4 m,n,mmx,mx
      m=ilist(1,ka-1)
      if(m .lt. 4)then
        if(m .ne. 0)then
          write(*,*)'tfree-too small segment: ',ka,m
        endif
        return
      endif
      lastpend=max(lastpend,ka+m-1)
      nnet=nnet-m
      ix=ka+1
      ik=iand(ix,mhash)*4+ich
      ik0=ik
      ip=klist(ik)
      do while(ip .ne. ik0)
        if(ip .ge. ix)then
          ik=ip
          ip=klist(ik)
          cycle
        elseif(ix .ge. ip+4 .and. ix .le. ip+6)then
          mx=ix-ip
        else
          mx=klist(ip+4)
          if(ip+mx .ne. ix)then
            ik=ip
            ip=klist(ik)
            cycle
          endif
          mmx=mod(mx,mhash)
          if(mmx .ge. 4 .and. mmx .le. 6)then
            ip0=icp+nindex*2
            ic=klist(ip0)
            do while(ic .ne. ip0)
              if(ic .eq. ip-2)then
                go to 10
              endif
            enddo
            ik=ip
            ip=klist(ik)
            cycle
          endif
        endif
 10     klist(klist(ip)+1)=ik
        klist(ik  )=klist(ip)
        klist(klist(ip-2)+1)=klist(ip-1)
        klist(klist(ip-1)  )=klist(ip-2)
        if(mx .gt. 6)then
          klist(klist(ip+2)+1)=klist(ip+3)
          klist(klist(ip+3)  )=klist(ip+2)
        endif
        m=m+mx
        ix=ip
        exit
      enddo
      ix1=ix+m
      ix2=ix1+2
      ik=iand(ix2,mhash)*4+ich+2
      ik0=ik
      ip=klist(ik)
      do while(ip .ne. ik0)
        if(ip .eq. ix2)then
          klist(klist(ip)+1)=ik
          klist(ik  )=klist(ip)
          klist(klist(ip-2)+1)=klist(ip-1)
          klist(klist(ip-1)  )=klist(ip-2)
          klist(klist(ip-4)+1)=klist(ip-3)
          klist(klist(ip-3)  )=klist(ip-4)
          m=m+klist(ip+2)
          go to 100
        endif
        ik=ip
        ip=klist(ik)
      enddo
      do n=4,6
        ik=iand(ix1+n,mhash)*4+ich
        ik0=ik
        ip=klist(ik)
        do while(ip .ne. ik0)
          if(ip .eq. ix1)then
            klist(klist(ip)+1)=ik
            klist(ik  )=klist(ip)
            klist(klist(ip-2)+1)=klist(ip-1)
            klist(klist(ip-1)  )=klist(ip-2)
            m=m+n
            exit
          endif
          ik=ip
          ip=klist(ik)
        enddo  
      enddo
 100  call tsetindexhash(ix-2,m)
      return
      end

      logical*4 function tfaloc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,n1
      n1=max(n,3)
      if(n1 .lt. nindex)then
        tfaloc=klist(icp+n1*2) .ne. icp+n1*2
      else
        tfaloc=.false.
      endif
      return
      end

      subroutine talocp(m,ip1)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 np,na,m,irtc
      integer*8 ip1,lmalloc8
      np=1+(m-1)/mpsize
      na=mpsize*np
      ip1=lmalloc8(na,irtc)+1
      if(irtc .ne. 0 .or. ip1 .le. 1)then
        na=mpsize
        ip1=lmalloc8(na)+1
        if(irtc .ne. 0 .or. ip1 .le. 1)then
          write(*,*)'Memory Full, request: ',m,ip1,lastpend
          stop 'talocp@italoc.f'
          ip1=-1
          return
        endif
      endif
      lastpend=max(lastpend,ip1+na-1)
      nmem=nmem+na
      nnet=nnet+na
      ilist(1,ip1-1)=na
      call tfree(ip1)
      return
      end

      subroutine talocinit
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 mfalloc
      integer*8 i
      data icp,nmem,nnet,maxic,minic/0,0,0,0,3/
      if(icp .eq. 0)then
        icp=mfalloc(nindex*2+mhash*4+10)+1
        lastpend=icp+nindex*2+mhash*4+9
        do i=icp,icp+nindex*2,2
          klist(i)=i
          klist(i+1)=i
        enddo
        ich=icp+nindex*2+4
        do i=ich,ich+mhash*4,4
          klist(i)=i
          klist(i+1)=i
          klist(i+2)=i+2
          klist(i+3)=i+2
        enddo
      endif
      return
      end

      subroutine tfreeme(ifree)
      implicit none
      include 'inc/TFCBK.inc'
      integer*4 ifree
      call tfree(int8(ifree+1))
      return
      end

      subroutine tfmalloc(isp1,kx,irtc)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,ktaloc,k
      integer*4 isp1,irtc,i,itfmessage,m
      kx=ktfoper+mtfnull
      do i=isp1+1,isp
        if(ktfnonrealq(ktastk(i)))then
          irtc=itfmessage(9,'General::wrongtype','"Real"')
          return
        endif
        m=max(3,(int(rtastk(i))+7)/8)
        k=ktaloc(m)
        if(k .le. 0)then
          irtc=itfmessage(9999,'Memory::alloc','""')
          return
        endif
        call tfree(k)
      enddo
      irtc=0
      return
      include 'inc/TFSF.inc'
      end

      subroutine tprintaloc(tag)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      character*(*) tag
      write(*,*)'aloc-',tag,nnet,nmem
      return
      end
