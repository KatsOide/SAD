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
      integer*4 n
      integer*8 m,ic1,i,n1,ic,m1,ip0,ip,ip1,ic2
      n1=max(n,3)
      m=n1+1
      call tfreecheck('ktaloc-0 ',1,0,dble(n))
      if(n1 .lt. nindex)then
        i=icp+n1*2
        ic=klist(i)
        if(ic .ne. i)then
          klist(i)=klist(ic)
          klist(klist(ic)+1)=i
          klist(klist(ic+3))  =klist(ic+2)
          klist(klist(ic+2)+1)=klist(ic+3)
          if(m .gt. 5)then
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
      integer*4 irtc
      integer*8 ip,m,ic,ic1,ia
      call tfreecheck('tsih-0 ',1,0,dble(ip))
      if(m .gt. nindex)then
        ic=icp+nindex*2
      else
        ic=icp+(m-1)*2
        maxic=max(ic,maxic)
      endif
      ia=ip
      klist(ia  )=klist(ic)
      klist(ia+1)=ic
      klist(klist(ic)+1)=ia
      klist(ic)=ia
      if(m .gt. 6)then
        klist(ia+6)=m
      endif
      call tfreecheck('tsih-1 ',1,0,dble(klist(ia)))
      ia=ia+2
      ic1=ich+iand(ia+m,mhash)*4
      klist(ia  )=klist(ic1)
      klist(ia+1)=ic1
      klist(klist(ic1)+1)=ia
      klist(ic1 )=ia
      call tfreecheck('tsih-2 ',1,0,dble(klist(ic1)))
      if(m .gt. 5)then
        ia=ia+2
        ic1=ich+iand(ia,mhash)*4+2
        klist(ia  )=klist(ic1)
        klist(ia+1)=ic1
        klist(klist(ic1)+1)=ia
        klist(ic1 )=ia
        call tfreecheck1('tsih-3 ',1,0,dble(klist(ich)),irtc)
        if(irtc .ne. 0)then
          write(*,*)'tsih-3 ',ich,ia,ic1
        endif
        if(m .gt. 6)then
          klist(ia+2)=m
        endif
      endif
      call tfreecheck('tsih ',1,0,dble(m))
      return
      end

      subroutine tfree(ia)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 ia
      integer*8 m,ix,ik,ip1,ik0,ix1,ix2,mx,n,ip
      if(ia .le. 0)then
        write(*,*)'tfree: zero or negative pointer',ia
        return
      endif
      call tfreecheck('tfree-0 ',1,0,dble(ia))
      m=ilist(1,ia-1)
      if(m .lt. 4)then
        if(m .ne. 0)then
          write(*,*)'tfree-too small segment: ',ia,m
        endif
        return
      endif
      nnet=nnet-m
      ix=ia+1
      ik=iand(ix,mhash)*4+ich
      ik0=ik
      ip=klist(ik)
      do while(ip .ne. ik0)
        if(ip+4 .eq. ix .or. ip+5 .eq. ix .or. ip+6 .eq. ix)then
          mx=ix-ip
        else
          mx=klist(ip+4)
          if(ip+mx .ne. ix)then
            ik=ip
            ip=klist(ik)
            cycle
          endif
        endif
        klist(klist(ip)+1)=ik
        klist(ik  )=klist(ip)
        klist(klist(ip-2)+1)=klist(ip-1)
        klist(klist(ip-1)  )=klist(ip-2)
        if(mx .gt. 5)then
          klist(klist(ip+2)+1)=klist(ip+3)
          klist(klist(ip+3)  )=klist(ip+2)
        endif
        m=m+mx
        ix=ip
        call tfreecheck('tfree-1 ',1,0,dble(ix))
        exit
      enddo
      ix1=ix+m
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
            if(n .eq. 6)then
              klist(klist(ip+2)+1)=klist(ip+3)
              klist(klist(ip+3)  )=klist(ip+2)
            endif
            m=m+n
            call tfreecheck('tfree-2 ',1,0,dble(m))
            go to 100
          endif
          ik=ip
          ip=klist(ik)
        enddo  
      enddo
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
          call tfreecheck('tfree-3 ',1,0,dble(m))
          exit
        endif
        ik=ip
        ip=klist(ik)
      enddo
 100  call tsetindexhash(ix-2,m)
      call tfreecheck('tfree-4 ',1,0,dble(ip))
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
        icp=mfalloc(nindex*2+mhash*4+10)+1
        lastpend=icp+nindex*2+mhash*4+9
        do i=icp,icp+nindex*2,2
          klist(i)=i
          klist(i+1)=i
        enddo
        ich=icp+nindex*2+2
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
        m=max(4,(int(vstk(ivstkoffset+i))+7)/8)
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
