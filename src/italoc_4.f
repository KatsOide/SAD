c$Header: /SAD/cvsroot/oldsad/src/italoc.f,v 1.59 2011/07/25 23:38:14 oide Exp $
      integer*4 function italoc(n)
      implicit none
      integer*4 n
      integer*8 ktaloc
      italoc=ktaloc(n)
      return
      end

      integer*4 function italoca(n)
      implicit none 
      integer*4 n
      integer*8 ktaloc
      italoca=ktaloc(n)
      return
      end

      integer*8 function ktaloc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,irtc
      integer*8 m,ic1,i,n1,ic,m1,ip0,ip1,ic2,ia
      n1=max(n,3)
      m=n1+1
c      call tfreecheck('ktaloc-0 ',1,0,dble(n),irtc)
      if(n1 .lt. nindex)then
        i=icp+n1*3
        ia=klist(i+1)
        if(ia .ne. i)then
          ic=klist(ia)
          klist(i+1)=klist(ic)
          klist(klist(klist(ic))+1)=i
          klist(klist(klist(ic+3))+2)=klist(ic+2)
          klist(klist(klist(ic+2))+3)=klist(ic+3)
          klist(klist(klist(ic+5))+4)=klist(ic+4)
          klist(klist(klist(ic+4))+5)=klist(ic+5)
          if(ic .ne. ia+1)then
c            write(*,*)'ktaloc-1 ',ia,ic,klist(ic+6)
            call tfree(int(ic))
          endif
          ilist(1,ia-1)=m
          nnet=nnet+m
          ktaloc=ia
          return
        endif
        ic1=i+min(max(m,minseg0),minseg1)*3
        if(ic1 .le. maxic)then
          ic2=min(maxic-3,ic1+15)
          do i=ic1,ic2,3
            ia=klist(i+1)
            if(ia .ne. i)then
              ic=klist(ia)
              klist(i+1)=klist(ic)
              klist(klist(klist(ic))+1)=i
              klist(klist(klist(ic+3))+2)=klist(ic+2)
              klist(klist(klist(ic+2))+3)=klist(ic+3)
              klist(klist(klist(ic+5))+4)=klist(ic+4)
              klist(klist(klist(ic+4))+5)=klist(ic+5)
              m1=klist(ic+6)
              call tsetindexhash(ia+m,m1-m)
              ilist(1,ia-1)=m
              nnet=nnet+m
              ktaloc=ia
              return
            endif
          enddo
          do i=maxic,ic2+3,-3
            ia=klist(i+1)
            if(ia .ne. i)then
              ic=klist(ia)
              klist(i+1)=klist(ic)
              klist(klist(klist(ic  ))+1)=i
              klist(klist(klist(ic+3))+2)=klist(ic+2)
              klist(klist(klist(ic+2))+3)=klist(ic+3)
              klist(klist(klist(ic+5))+4)=klist(ic+4)
              klist(klist(klist(ic+4))+5)=klist(ic+5)
              m1=klist(ic+6)
              call tsetindexhash(ia+m,m1-m)
              ilist(1,ia-1)=m
              nnet=nnet+m
              ktaloc=ia
              if(klist(i) .eq. ic)then
                maxic=i-3
              else
                maxic=i
              endif
              return
            endif
          enddo
          maxic=ic1-3
        endif
      endif
      ip0=icp+nindex*3
 1000 ia=klist(ip0+1)
      do while(ia .ne. ip0)
        ic=klist(ia)
        m1=klist(ic+6)
        if(m1 .eq. m)then
          klist(klist(klist(ic+1))  )=klist(ic  )
          klist(klist(klist(ic  ))+1)=klist(ic+1)
          klist(klist(klist(ic+3))+2)=klist(ic+2)
          klist(klist(klist(ic+2))+3)=klist(ic+3)
          klist(klist(klist(ic+5))+4)=klist(ic+4)
          klist(klist(klist(ic+4))+5)=klist(ic+5)
          ilist(1,ia-1)=m
          nnet=nnet+m
          ktaloc=ia
          return
        elseif(m1 .ge. m+minseg2)then
          klist(klist(klist(ic+1))  )=klist(ic  )
          klist(klist(klist(ic  ))+1)=klist(ic+1)
          klist(klist(klist(ic+3))+2)=klist(ic+2)
          klist(klist(klist(ic+2))+3)=klist(ic+3)
          klist(klist(klist(ic+5))+4)=klist(ic+4)
          klist(klist(klist(ic+4))+5)=klist(ic+5)
          call tsetindexhash(ia+m,m1-m)
          ilist(1,ia-1)=m
          nnet=nnet+m
          ktaloc=ia
          return
        endif
        ia=klist(ic)
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
      integer*4 irtc
      integer*8 ip,m,ic,ic1,ktaloc,i
      if(m .gt. nindex)then
        i=icp+nindex*3
      else
        i=icp+(m-1)*3
        maxic=max(i,maxic)
      endif
      if(m .gt. 8)then
        ic=ip+1
      else
        ic=ktaloc(8)
      endif
c      call tfreecheck('tsih-0 ',1,0,dble(ip),irtc)
      klist(ip)=ic
      klist(ic  )=klist(i+1)
      klist(ic+1)=i
      klist(klist(klist(i+1))+1)=ip
      klist(i+1)=ip
      i=ich+iand(ip+m,mhash)*5
      klist(ic+2)=klist(i+1)
      klist(ic+3)=i
      klist(klist(klist(i+1))+3)=ip
      klist(i+1)=ip
      i=ich+iand(ip  ,mhash)*5
      klist(ic+4)=klist(i+3)
      klist(ic+5)=i
      klist(klist(klist(i+3))+5)=ip
      klist(i+3)=ip
      klist(ic+6)=m
c      call tfreecheck('tsih-1 ',1,0,dble(m),irtc)
c      if(irtc .ne. 0)then
c        write(*,*)'tsih-1 ',i,ip,ic
c      endif
      return
      end

      subroutine tfree(ia)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 ia,irtc
      integer*8 m,ix,ik,ik0,ix2,mx,ip,nfree,ifree(3),i,ic
      if(ia .eq. 0)then
        write(*,*)'tfree: detect null pointer'
        return
      endif
      m=ilist(1,ia-1)
      if(m .lt. 4)then
        if(m .ne. 0)then
          write(*,*)'tfree-too small segment: ',ia,m
        endif
        return
      endif
      nnet=nnet-m
      nfree=0
      ix=ia
      ik=iand(ix,mhash)*5+ich
      ip=klist(ik+1)
      ik0=ik
      do while(ip .ne. ik0)
        ic=klist(ip)
        mx=klist(ic+6)
        if(ip+mx .eq. ix)then
          klist(klist(klist(ic  ))+1)=klist(ic+1)
          klist(klist(klist(ic+1))  )=klist(ic)
          klist(klist(klist(ic+2))+3)=ik
          klist(klist(klist(ic+3))+2)=klist(ic+2)
          klist(klist(klist(ic+4))+5)=klist(ic+5)
          klist(klist(klist(ic+5))+4)=klist(ic+4)
          m=m+mx
          ix=ip
          if(ic .ne. ip+1)then
            nfree=nfree+1
            ifree(nfree)=ic
          endif
          exit
        endif
        ik=ip
        ip=klist(ic+2)
      enddo
      ix2=ix+m
      ik=iand(ix2,mhash)*5+ich
      ip=klist(ik+3)
      ik0=ik
      do while(ip .ne. ik0)
        ic=klist(ip)
        if(ip .eq. ix2)then
          klist(klist(klist(ic  ))+1)=klist(ic+1)
          klist(klist(klist(ic+1))  )=klist(ic)
          klist(klist(klist(ic+2))+3)=klist(ic+3)
          klist(klist(klist(ic+3))+2)=klist(ic+2)
          klist(klist(klist(ic+4))+5)=ik
          klist(klist(klist(ic+5))+4)=klist(ic+4)
          m=m+klist(ic+6)
          if(ic .ne. ip+1)then
            nfree=nfree+1
            ifree(nfree)=ic
          endif
          exit
        endif
        ik=ip
        ip=klist(ic+4)
      enddo
      call tsetindexhash(ix,m)
c      call tfreecheck('tfree-8',1,0,dble(ia),irtc)
c      if(irtc .ne. 0)then
c        write(*,*)'tfree ',ia,ix,ilist(1,ia-1),m
c      endif
      do i=1,nfree
        call tfree(int(ifree(i)))
      enddo
      return
      end

      subroutine taloc2(na,nb,iaa,iab)
      implicit none 
      include 'inc/TFCBK.inc'
      integer*4 na,nb,iaa,iab,italoc,na1,nb1
      logical*4 tfaloc
      if(tfaloc(na) .or. tfaloc(nb))then
        iaa=italoc(na)
        iab=italoc(nb)
      else
        na1=max(na,3)
        nb1=max(nb,3)
        iaa=italoc(na1+nb1+1)
        iab=iaa+na1+1
        ilist(1,iaa-1)=na1+1
        ilist(1,iab-1)=nb1+1
      endif
      return
      end

      subroutine taloc3(na,nb,nc,iaa,iab,iac)
      implicit none 
      include 'inc/TFCBK.inc'
      integer*4 na,nb,iaa,iab,italoc,nc,iac,na1,nb1,nc1
      logical*4 tfaloc
      if(tfaloc(na))then
        iaa=italoc(na)
        call taloc2(nb,nc,iab,iac)
      elseif(tfaloc(nb))then
        iab=italoc(nb)
        call taloc2(na,nc,iaa,iac)
      elseif(tfaloc(nc))then
        iac=italoc(nc)
        call taloc2(na,nb,iaa,iab)
      else
        na1=max(3,na)
        nb1=max(3,nb)
        nc1=max(3,nc)
        iaa=italoc(na1+nb1+nc1+2)
        iab=iaa+na1+1
        iac=iab+nb1+1
        ilist(1,iaa-1)=na1+1
        ilist(1,iab-1)=nb1+1
        ilist(1,iac-1)=nc1+1
      endif
      return
      end

      logical*4 function tfaloc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,ik
      ik=icp+max(n,3)*3
      if(n .le. nindex)then
        tfaloc=klist(ik+1) .ne. ik
      else
        tfaloc=.false.
      endif
      return
      end

      subroutine talocp(m,ip1)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 np,na,mfalloc
      integer*8 ip1,m
      np=1+(m+4)/mpsize
      na=mpsize*np
      ip1=mfalloc(na)+1
      if(ip1 .le. 1)then
        na=mpsize
        ip1=mfalloc(na)+1
        if(ip1 .le. 1)then
          write(*,*)'Memory Full, request: ',m,ip1,lastpend
          stop 'talocp@italoc.f'
          ip1=-1
          return
        endif
      endif
c      if(ip1 .ge. lastpend .and. ip1 .lt. lastpend+4)then
c        na=na+ip1-lastpend
c        ip1=lastpend
c      endif
      lastpend=ip1+na
      nmem=nmem+na
      nnet=nnet+na
      ilist(1,ip1-1)=na
      call tfree(int(ip1))
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
        icp=mfalloc(nindex*3+mhash*5+15)+1
        lastpend=icp+nindex*3+mhash*5+14
        do i=icp,icp+nindex*3,3
          klist(i)=i+1
          klist(i+1)=i
          klist(i+2)=i
        enddo
        ich=icp+nindex*3+3
        do i=ich,ich+mhash*5,5
          klist(i)=i-1
          klist(i+1)=i
          klist(i+2)=i
          klist(i+3)=i
          klist(i+4)=i
        enddo
      endif
      return
      end

      subroutine tfreeme(ifree)
      implicit none
      include 'inc/TFCBK.inc'
      integer*4 ifree
      ilist(1,ifree)=ilist(1,ifree)+1
      call tfree(ifree+1)
      return
      end

      subroutine tfmalloc(isp1,itx,iax,vx,irtc)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*4 isp1,itx,iax,irtc,i,itfmessage,m,k,italoc
      real*8 vx
      itx=ntfoper
      iax=mtfnull
      do i=isp1+1,isp
        if(itastk(1,i) .ne. ntfreal)then
          irtc=itfmessage(9,'General::wrongtype','"Real"')
          return
        endif
        m=max(1,(int(vstk(ivstkoffset+i))+7)/8)
        k=italoc(m)
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
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      character*(*) tag
      write(*,*)'aloc-',tag,nnet,nmem
      return
      end
