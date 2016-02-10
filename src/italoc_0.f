c$Header: /SAD/cvsroot/oldsad/src/italoc.f,v 1.59 2011/07/25 23:38:14 oide Exp $
      integer*4 function italoc(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,m,ic1,i,n1,ic,m1,ip0,ip,ip1,ic2
      n1=max(n,3)
      m=n1+1
      if(n1 .lt. nindex)then
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          italoc=ic
          return
        endif
        ic1=i+min(m,minseg1)
        if(ic1 .le. maxic)then
          ic2=min(maxic-1,ic1+5)
          do i=ic1,ic2
            ic=ilist(1,i)
            if(ic .ne. i)then
              ilist(1,i)=ilist(1,ic)
              ilist(2,ilist(1,ic))=i
              ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
              if(ilist(2,ic+1) .ne. 0)then
                ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
              else
                ilistdummy(2)=ilist(2,ic-1)
              endif
              ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
              if(ilist(1,ic+1) .ne. 0)then
                ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
              else
                ilistdummy(1)=ilist(1,ic-1)
              endif
              call tsetindexhash1(ic+m,i-icp-m+1)
              ilist(1,ic-1)=m
              nnet=nnet+m
              italoc=ic
              return
            endif
          enddo
          do i=maxic,ic2+1,-1
            ic=ilist(1,i)
            if(ic .ne. i)then
              ilist(1,i)=ilist(1,ic)
              ilist(2,ilist(1,ic))=i
              ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
              if(ilist(2,ic+1) .ne. 0)then
                ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
              else
                ilistdummy(2)=ilist(2,ic-1)
              endif
              ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
              if(ilist(1,ic+1) .ne. 0)then
                ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
              else
                ilistdummy(1)=ilist(1,ic-1)
              endif
              call tsetindexhash1(ic+m,i-icp-m+1)
              ilist(1,ic-1)=m
              nnet=nnet+m
              if(ilist(1,i) .eq. ic)then
                maxic=i-1
              else
                maxic=i
              endif
              italoc=ic
              return
            endif
          enddo
          maxic=ic1-1
        endif
      endif
      ip0=icp+nindex
 1000 ip=ilist(1,ip0)
      do while(ip .ne. ip0)
        m1=ilist(1,ip+2)
        if(m1 .eq. m)then
          ilist(2,ilist(1,ip))=ilist(2,ip)
          ilist(1,ilist(2,ip))=ilist(1,ip)
          ilist(2,ilist(2,ip-1))=ilist(2,ip+1)
          if(ilist(2,ip+1) .ne. 0)then
            ilist(2,ilist(2,ip+1)-2)=ilist(2,ip-1)
          else
            ilistdummy(2)=ilist(2,ip-1)
          endif
          ilist(1,ilist(1,ip-1))=ilist(1,ip+1)
          if(ilist(1,ip+1) .ne. 0)then
            ilist(1,ilist(1,ip+1)-2)=ilist(1,ip-1)
          else
            ilistdummy(1)=ilist(1,ip-1)
          endif
          ilist(1,ip-1)=m
          nnet=nnet+m
          italoc=ip
          return
        elseif(m1 .ge. m+minseg2)then
          ilist(2,ilist(1,ip))=ilist(2,ip)
          ilist(1,ilist(2,ip))=ilist(1,ip)
          ilist(2,ilist(2,ip-1))=ilist(2,ip+1)
          if(ilist(2,ip+1) .ne. 0)then
            ilist(2,ilist(2,ip+1)-2)=ilist(2,ip-1)
          else
            ilistdummy(2)=ilist(2,ip-1)
          endif
          ilist(1,ilist(1,ip-1))=ilist(1,ip+1)
          if(ilist(1,ip+1) .ne. 0)then
            ilist(1,ilist(1,ip+1)-2)=ilist(1,ip-1)
          else
            ilistdummy(1)=ilist(1,ip-1)
          endif
          call tsetindexhash(ip+m,m1-m)
          ilist(1,ip-1)=m
          nnet=nnet+m
          italoc=ip
          return
        endif
        ip=ilist(1,ip)
      enddo
      call talocp(m,ip1)
      if(ip1 .gt. 0)then
        go to 1000
      endif
      italoc=-1
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

      integer*4 function italoca(n)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 n,m,ic1,i,n1,ic,ii,m1,ip0,ip,ip1,ic2
      n1=max(n,3)
      m=n1+1
      if(n1 .lt. nindex)then
        i=icp+n1
        if(i .le. maxic-3)then
          ic=ilist(1,i)
          if(ic .ne. i)then
            ilist(1,i)=ilist(1,ic)
            ilist(2,ilist(1,ic))=i
            ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
            if(ilist(2,ic+1) .ne. 0)then
              ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
            else
              ilistdummy(2)=ilist(2,ic-1)
            endif
            ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
            if(ilist(1,ic+1) .ne. 0)then
              ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
            else
              ilistdummy(1)=ilist(1,ic-1)
            endif
            ilist(1,ic-1)=i-i+m
            nnet=nnet+ilist(1,ic-1)
            italoca=ic
            return
          endif
          ii=i+1
          ic=ilist(1,ii)
          if(ic .ne. ii)then
            ilist(1,ii)=ilist(1,ic)
            ilist(2,ilist(1,ic))=ii
            ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
            if(ilist(2,ic+1) .ne. 0)then
              ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
            else
              ilistdummy(2)=ilist(2,ic-1)
            endif
            ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
            if(ilist(1,ic+1) .ne. 0)then
              ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
            else
              ilistdummy(1)=ilist(1,ic-1)
            endif
            ilist(1,ic-1)=ii-i+m
            nnet=nnet+ilist(1,ic-1)
            italoca=ic
            return
          endif
          ii=ii+1
          ic=ilist(1,ii)
          if(ic .ne. ii)then
            ilist(1,ii)=ilist(1,ic)
            ilist(2,ilist(1,ic))=ii
            ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
            if(ilist(2,ic+1) .ne. 0)then
              ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
            else
              ilistdummy(2)=ilist(2,ic-1)
            endif
            ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
            if(ilist(1,ic+1) .ne. 0)then
              ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
            else
              ilistdummy(1)=ilist(1,ic-1)
            endif
            ilist(1,ic-1)=ii-i+m
            nnet=nnet+ilist(1,ic-1)
            italoca=ic
            return
          endif
          ii=ii+1
          ic=ilist(1,ii)
          if(ic .ne. ii)then
            ilist(1,ii)=ilist(1,ic)
            ilist(2,ilist(1,ic))=ii
            ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
            if(ilist(2,ic+1) .ne. 0)then
              ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
            else
              ilistdummy(2)=ilist(2,ic-1)
            endif
            ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
            if(ilist(1,ic+1) .ne. 0)then
              ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
            else
              ilistdummy(1)=ilist(1,ic-1)
            endif
            ilist(1,ic-1)=ii-i+m
            nnet=nnet+ilist(1,ic-1)
            italoca=ic
            return
          endif
        endif
        ic1=i+min(m,minseg1)
        if(ic1 .le. maxic)then
          ic2=min(maxic-1,ic1+5)
          do i=ic1,ic2
            ic=ilist(1,i)
            if(ic .ne. i)then
              ilist(1,i)=ilist(1,ic)
              ilist(2,ilist(1,ic))=i
              ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
              if(ilist(2,ic+1) .ne. 0)then
                ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
              else
                ilistdummy(2)=ilist(2,ic-1)
              endif
              ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
              if(ilist(1,ic+1) .ne. 0)then
                ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
              else
                ilistdummy(1)=ilist(1,ic-1)
              endif
              call tsetindexhash1(ic+m,i-icp-m+1)
              ilist(1,ic-1)=m
              nnet=nnet+m
              italoca=ic
              return
            endif
          enddo
          do i=maxic,ic2+1,-1
            ic=ilist(1,i)
            if(ic .ne. i)then
              ilist(1,i)=ilist(1,ic)
              ilist(2,ilist(1,ic))=i
              ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
              if(ilist(2,ic+1) .ne. 0)then
                ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
              else
                ilistdummy(2)=ilist(2,ic-1)
              endif
              ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
              if(ilist(1,ic+1) .ne. 0)then
                ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
              else
                ilistdummy(1)=ilist(1,ic-1)
              endif
              call tsetindexhash1(ic+m,i-icp-m+1)
              ilist(1,ic-1)=m
              nnet=nnet+m
              if(ilist(1,i) .eq. ic)then
                maxic=i-1
              else
                maxic=i
              endif
              italoca=ic
              return
            endif
          enddo
          maxic=ic1-1
        endif
      endif
      ip0=icp+nindex
 1000 ip=ilist(1,ip0)
      do while(ip .ne. ip0)
        m1=ilist(1,ip+2)
        if(m1 .eq. m)then
          ilist(2,ilist(1,ip))=ilist(2,ip)
          ilist(1,ilist(2,ip))=ilist(1,ip)
          ilist(2,ilist(2,ip-1))=ilist(2,ip+1)
          if(ilist(2,ip+1) .ne. 0)then
            ilist(2,ilist(2,ip+1)-2)=ilist(2,ip-1)
          else
            ilistdummy(2)=ilist(2,ip-1)
          endif
          ilist(1,ilist(1,ip-1))=ilist(1,ip+1)
          if(ilist(1,ip+1) .ne. 0)then
            ilist(1,ilist(1,ip+1)-2)=ilist(1,ip-1)
          else
            ilistdummy(1)=ilist(1,ip-1)
          endif
          ilist(1,ip-1)=m
          nnet=nnet+m
          italoca=ip
          return
        elseif(m1 .ge. m+minseg2)then
          ilist(2,ilist(1,ip))=ilist(2,ip)
          ilist(1,ilist(2,ip))=ilist(1,ip)
          ilist(2,ilist(2,ip-1))=ilist(2,ip+1)
          if(ilist(2,ip+1) .ne. 0)then
            ilist(2,ilist(2,ip+1)-2)=ilist(2,ip-1)
          else
            ilistdummy(2)=ilist(2,ip-1)
          endif
          ilist(1,ilist(1,ip-1))=ilist(1,ip+1)
          if(ilist(1,ip+1) .ne. 0)then
            ilist(1,ilist(1,ip+1)-2)=ilist(1,ip-1)
          else
            ilistdummy(1)=ilist(1,ip-1)
          endif
          call tsetindexhash(ip+m,m1-m)
          ilist(1,ip-1)=m
          nnet=nnet+m
          italoca=ip
          return
        endif
        ip=ilist(1,ip)
      enddo
      call talocp(m,ip1)
      if(ip1 .gt. 0)then
        go to 1000
      endif
      italoca=-1
      return
      end

      subroutine tsetindexhash1(ip,m)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 ip,m,ic
      integer*4 ia,ip1
      ic=icp+m-1
      maxic=max(ic,maxic)
      ilist(1,ip)=ic
      ilist(1,ilist(2,ic))=ip
      ilist(2,ip)=ilist(2,ic)
      ilist(2,ic)=ip
      ilist(1,ip+2)=m
      ia=ip+1
      ip1=ich+iand(ia,mhash)
      ilist(2,ia)=ilist(2,ip1)
      if(ilist(2,ip1) .ne. 0)then
        ilist(2,ilist(2,ip1)-2)=ia
      else
        ilistdummy(2)=ia
      endif
      ilist(2,ip1)=ia
      ilist(2,ia-2)=ip1
      ip1=ich+iand(ia+m,mhash)
      ilist(1,ia)=ilist(1,ip1)
      if(ilist(1,ip1) .ne. 0)then
        ilist(1,ilist(1,ip1)-2)=ia
      else
        ilistdummy(1)=ia
      endif
      ilist(1,ip1)=ia
      ilist(1,ia-2)=ip1
      return
      end

      subroutine tsetindexhash(ip,m)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 ip,m,ic
      integer*4 ia,ip1
      if(m .gt. nindex)then
        ic=icp+nindex
        ip1=ilist(1,ic)
        do while(ip1 .lt. ip .and. ip1 .ne. ic)
          ip1=ilist(1,ip1)
        enddo
        ilist(1,ip)=ip1
        ilist(1,ilist(2,ip1))=ip
        ilist(2,ip)=ilist(2,ip1)
        ilist(2,ip1)=ip
      else
        ic=icp+m-1
        maxic=max(ic,maxic)
        ilist(1,ip)=ic
        ilist(1,ilist(2,ic))=ip
        ilist(2,ip)=ilist(2,ic)
        ilist(2,ic)=ip
      endif
      ilist(1,ip+2)=m
      ia=ip+1
      ip1=ich+iand(ia,mhash)
      ilist(2,ia)=ilist(2,ip1)
      if(ilist(2,ip1) .ne. 0)then
        ilist(2,ilist(2,ip1)-2)=ia
      else
        ilistdummy(2)=ia
      endif
      ilist(2,ip1)=ia
      ilist(2,ia-2)=ip1
      ip1=ich+iand(ia+m,mhash)
      ilist(1,ia)=ilist(1,ip1)
      if(ilist(1,ip1) .ne. 0)then
        ilist(1,ilist(1,ip1)-2)=ia
      else
        ilistdummy(1)=ia
      endif
      ilist(1,ip1)=ia
      ilist(1,ia-2)=ip1
      return
      end

      subroutine tfree(ia)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 ia,m,ix,ik,il,ih,ip,ip1,m1,ic
      if(ia .eq. 0)then
        write(*,*)'tfree: detect null pointer'
        return
      endif
      m=ilist(1,ia-1)
      if(m .lt. 4)then
        if(m .ne. 0)then
          write(*,*)'tfree ',ia,m
        endif
        return
      endif
      nnet=nnet-m
      ix=ia+m+1
      ik=iand(ix,mhash)+ich
      ip=ilist(2,ik)
      do while(ip .ne. 0)
        if(ip .eq. ix)then
          ilist(2,ilist(2,ip-2))=ilist(2,ip)
          if(ilist(2,ip) .ne. 0)then
            ilist(2,ilist(2,ip)-2)=ilist(2,ip-2)
          else
            ilistdummy(2)=ilist(2,ip-2)
          endif
          ilist(1,ilist(1,ip-2))=ilist(1,ip)
          if(ilist(1,ip) .ne. 0)then
            ilist(1,ilist(1,ip)-2)=ilist(1,ip-2)
          else
            ilistdummy(1)=ilist(1,ip-2)
          endif
          m=m+ilist(1,ip+1)
          ilist(1,ilist(2,ip-1))=ilist(1,ip-1)
          ilist(2,ilist(1,ip-1))=ilist(2,ip-1)
          ik=iand(ia+m+1,mhash)+ich
          exit
        endif
        ip=ilist(2,ip)
      enddo
      ix=ia+1
      ih=iand(ix,mhash)+ich
      ip=ilist(1,ih)
      do while(ip .ne. 0)
        m1=ilist(1,ip+1)
        if(ip+m1 .eq. ix)then
          ilist(1,ilist(2,ip-1))=ilist(1,ip-1)
          ilist(2,ilist(1,ip-1))=ilist(2,ip-1)
          ilist(1,ilist(1,ip-2))=ilist(1,ip)
          if(ilist(1,ip) .ne. 0)then
            ilist(1,ilist(1,ip)-2)=ilist(1,ip-2)
          else
            ilistdummy(1)=ilist(1,ip-2)
          endif
          m1=m1+m
          il=iand(ip+m1,mhash)+ich
          ilist(1,ip)=ilist(1,il)
          if(ilist(1,il) .ne. 0)then
            ilist(1,ilist(1,il)-2)=ip
          else
            ilistdummy(1)=ip
          endif
          ilist(1,il)=ip
          ilist(1,ip-2)=il
          ia=ip-1
          if(m1 .gt. nindex)then
            ic=icp+nindex
            ip1=ilist(1,ic)
            do while(ip1 .lt. ia .and. ip1 .ne. ic)
              ip1=ilist(1,ip1)
            enddo
            ilist(1,ia)=ip1
            ilist(1,ilist(2,ip1))=ia
            ilist(2,ia)=ilist(2,ip1)
            ilist(2,ip1)=ia
            ilist(1,ia+2)=m1
          else
            ic=icp+m1-1
            maxic=max(ic,maxic)
            ilist(1,ia)=ic
            ilist(1,ilist(2,ic))=ia
            ilist(2,ia)=ilist(2,ic)
            ilist(2,ic)=ia
            ilist(1,ia+2)=m1
          endif
          return
        endif
        ip=ilist(1,ip)
      enddo
      if(m .gt. nindex)then
        ic=icp+nindex
        ip1=ilist(1,ic)
        do while(ip1 .lt. ia .and. ip1 .ne. ic)
          ip1=ilist(1,ip1)
        enddo
        ilist(1,ia)=ip1
        ilist(1,ilist(2,ip1))=ia
        ilist(2,ia)=ilist(2,ip1)
        ilist(2,ip1)=ia
        ilist(1,ia+2)=m
      else
        ic=icp+m-1
        maxic=max(ic,maxic)
        ilist(1,ia)=ic
        ilist(1,ilist(2,ic))=ia
        ilist(2,ia)=ilist(2,ic)
        ilist(2,ic)=ia
        ilist(1,ia+2)=m
      endif
      ilist(2,ix)=ilist(2,ih)
      if(ilist(2,ih) .ne. 0)then
        ilist(2,ilist(2,ih)-2)=ix
      else
        ilistdummy(2)=ix
      endif
      ilist(2,ih)=ix
      ilist(2,ix-2)=ih
      ilist(1,ix)=ilist(1,ik)
      if(ilist(1,ik) .ne. 0)then
        ilist(1,ilist(1,ik)-2)=ix
      else
        ilistdummy(1)=ix
      endif
      ilist(1,ik)=ix
      ilist(1,ix-2)=ik
      return
      end

      subroutine taloc2(na,nb,iaa,iab)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 m,i,n1,ic,na,nb,iaa,iab,italoc,na1,nb1
      na1=max(na,3)
      n1=na1
      if(n1 .lt. nindex)then
        m=n1+1
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
c          call tremovehash(ic+1)
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          iaa=ic
          iab=italoc(nb)
          return
        endif
      endif
      nb1=max(nb,3)
      n1=nb1
      if(n1 .lt. nindex)then
        m=n1+1
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
c          call tremovehash(ic+1)
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          iab=ic
          iaa=italoc(na)
          return
        endif
      endif
      iaa=italoc(na1+nb1+1)
      iab=iaa+na1+1
      ilist(1,iaa-1)=na1+1
      ilist(1,iab-1)=nb1+1
      return
      end

      subroutine taloc3(na,nb,nc,iaa,iab,iac)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 m,i,n1,ic,
     $     na,nb,iaa,iab,italoc,na1,nb1,nc,iac,nc1
      na1=max(na,3)
      n1=na1
      if(n1 .lt. nindex)then
        m=n1+1
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
c          call tremovehash(ic+1)
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          iaa=ic
          call taloc2(nb,nc,iab,iac)
          return
        endif
      endif
      nb1=max(nb,3)
      n1=nb1
      if(n1 .lt. nindex)then
        m=n1+1
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
c          call tremovehash(ic+1)
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          iab=ic
          call taloc2(na,nc,iaa,iac)
          return
        endif
      endif
      nc1=max(nc,3)
      n1=nc1
      if(n1 .lt. nindex)then
        m=n1+1
        i=icp+n1
        ic=ilist(1,i)
        if(ic .ne. i)then
          ilist(1,i)=ilist(1,ic)
          ilist(2,ilist(1,ic))=i
c          call tremovehash(ic+1)
          ilist(2,ilist(2,ic-1))=ilist(2,ic+1)
          if(ilist(2,ic+1) .ne. 0)then
            ilist(2,ilist(2,ic+1)-2)=ilist(2,ic-1)
          else
            ilistdummy(2)=ilist(2,ic-1)
          endif
          ilist(1,ilist(1,ic-1))=ilist(1,ic+1)
          if(ilist(1,ic+1) .ne. 0)then
            ilist(1,ilist(1,ic+1)-2)=ilist(1,ic-1)
          else
            ilistdummy(1)=ilist(1,ic-1)
          endif
          ilist(1,ic-1)=m
          nnet=nnet+m
          iac=ic
          call taloc2(na,nb,iaa,iab)
          return
        endif
      endif
      iaa=italoc(na1+nb1+nc1+2)
      iab=iaa+na1+1
      iac=iab+nb1+1
      ilist(1,iaa-1)=na1+1
      ilist(1,iab-1)=nb1+1
      ilist(1,iac-1)=nc1+1
      return
      end

      subroutine talocp(m,ip1)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*4 m,np,na,ip1,mfalloc
      np=1+(m+4)/mpsize
      na=mpsize*np
      ip1=mfalloc(na)+1
      if(ip1 .le. 1)then
        na=mpsize
        ip1=mfalloc(na)+1
c     write(*,*)'italoc ',ip1,na
        if(ip1 .le. 1)then
          write(*,*)'Memory Full, request: ',m,ip1,lastpend
          stop 'talocp@italoc.f'
          ip1=-1
          return
        endif
      endif
      if(ip1 .ge. lastpend .and. ip1 .le. lastpend+4)then
        na=na+ip1-lastpend
        ip1=lastpend
      endif
      lastpend=ip1+na
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
      integer*4 i,mfalloc
      data icp,nmem,nnet,maxic,minic/0,0,0,0,3/
      if(icp .eq. 0)then
        icp=mfalloc(nindex+mhash+10)+1
        lastpend=icp+nindex+mhash+9
        do i=icp,icp+nindex
          ilist(1,i)=i
          ilist(2,i)=i
        enddo
        do i=icp+nindex+1,icp+nindex+1+mhash
          rlist(i)=0.d0
        enddo
        ilist(2,icp-1)=0
        ich=icp+nindex+1
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
        tfaloc=ilist(1,icp+n1) .ne. icp+n1
      else
        tfaloc=.false.
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
