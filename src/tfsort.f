      subroutine tfsort(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kf
      type (sad_dlist), pointer :: kl,klx
      integer*4 isp1,mode,irtc,narg,m,isp0,itfmessage
      narg=isp-isp1
      if(narg .gt. 2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+1),kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      irtc=0
      m=kl%nl
      if(m .le. 1)then
        kx=dtastk(isp1+1)
        return
      endif
      kf%k=merge(ktastk(isp),ktfref,narg .eq. 2)
      isp0=isp
      call tfsortl(kl,ktfreallistq(kl),m,mode,kf,.false.,irtc)
      if(irtc .ne. 0)then
        isp=isp0
        return
      endif
      kx%k=ktflist+ktfmakelist(isp0,klx)
      klx%head=dtfcopy(kl%head)
      isp=isp0
      return
      end

      subroutine tfsortl(kl,av,n,mode,kf,ins,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kf
      type (sad_dlist) kl
      integer*8 isort
      integer*4 n,mode,irtc,itab(n),i,j,isp0
      logical*4 av,ins
      integer*8 iafsort,ktfsymbolf
      save iafsort
      data iafsort/0/
      if(iafsort .eq. 0)then
        iafsort=ktfsymbolf('$SortMethodID',13,.false.)-4
      endif
c      do i=1,n
      itab(1:n)=[(i,i=1,n)]
c      enddo
      isort=int(rlist(iafsort))
      if(isort .eq. 1)then
        call tfsortml(itab,kl,av,n,mode,kf,irtc)
      else
        call tfsortql(itab,kl,av,n,mode,kf,irtc)
      endif
      isp0=isp
      if(mode .eq. 0)then
        dtastk(isp0+1:isp0+n)=kl%dbody(itab(1:n))
        isp=isp0+n
      elseif(ins)then
        do i=1,n
          j=itab(i)
          isp=isp+1
          ktastk(isp)=merge(j,0,j .ge. 0)
        enddo
      else
        do i=1,n
          j=itab(i)
          if(j .ge. 0)then
            isp=isp+1
            dtastk(isp)=kl%dbody(j)
          endif
        enddo
      endif
      return
      end

      recursive subroutine tfsortql(itab,kl,av,n,mode,kf,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kf
      type (sad_dlist) kl
      integer*4 n,itab(n),mode,irtc,itforderl,
     $     ip1,ip2,m,i1,i2,im,is,l1,l2,l
      real*8 v1,v2
      logical*4 av
      irtc=0
      if(n .le. 1)then
        return
      endif
      i1=itab(1)
      i2=itab(n)
      if(av .and. kf%k .eq. ktfref)then
        v1=kl%rbody(i1)
        v2=kl%rbody(i2)
        l=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
        if(mode .eq. 0 .or. l .ne. 0)then
          if(l .gt. 0)then
            is=i1
            i1=i2
            i2=is
          endif
        else
          itab(1)=min(i1,i2)
          itab(n)=itab(1)-i1-i2
          call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
          return
        endif
        if(n .eq. 2)then
          itab(1)=i1
          itab(2)=i2
          return
        endif
        m=(n+1)/2
        im=itab(m)
        v1=kl%rbody(i1)
        v2=kl%rbody(im)
        l=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
        if(mode .eq. 0 .or. l .ne. 0)then
          if(l .lt. 0)then
            v1=kl%rbody(im)
            v2=kl%rbody(i2)
            l=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
            if(mode .eq. 0 .or. l .ne. 0)then
              if(l .gt. 0)then
                is=im
                im=i2
                i2=is
              endif
            else
              itab(1)=i1
              itab(n)=-i2
              call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
              return
            endif
          elseif(l .ne. 0)then
            is=im
            im=i1
            i1=is
          endif
        else
          itab(1)=-i1
          itab(n)=i2
          call tfsortql(itab(2),kl,av,n-1,mode,kf,irtc)
          return
        endif
        itab(1)=i1
        itab(n)=i2
        itab(m)=im
        if(n .eq. 3)then
          return
        endif
        itab(m)=itab(2)
        itab(2)=im
        ip1=3
        ip2=n-1
        if(mode .eq. 0)then
          do while(ip1 .le. ip2)
            l1=0
            do while(l1 .le. 0 .and. ip1 .le. ip2)
              if(kl%rbody(itab(ip1)) .le. kl%rbody(im))then
                ip1=ip1+1
              else
                l1=1
              endif
            enddo
            l2=0
            do while(l2 .le. 0 .and. ip1 .le. ip2)
              if(kl%rbody(im) .le. kl%rbody(itab(ip2)))then
                ip2=ip2-1
              else
                l2=1
              endif
            enddo
            if(ip2 .gt. ip1)then
              is=itab(ip1)
              itab(ip1)=itab(ip2)
              itab(ip2)=is
              ip1=ip1+1
              ip2=ip2-1
            endif
          enddo
        else
          do while(ip1 .le. ip2)
            l1=0
            do while(l1 .le. 0 .and. ip1 .le. ip2)
              v1=kl%rbody(itab(ip1))
              v2=kl%rbody(im)
              l=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
              if(l1 .ne. 0)then
                if(l1 .le. 0)then
                  ip1=ip1+1
                endif
              else
                itab(n)=-itab(ip1)
                itab(ip1)=i2
                call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
                return
              endif
            enddo
            l2=0
            do while(l2 .le. 0 .and. ip1 .le. ip2)
              v1=kl%rbody(im)
              v2=kl%rbody(itab(ip2))
              l2=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
              if(l2 .ne. 0)then
                if(l2 .le. 0)then
                  ip2=ip2-1
                  if(ip2 .eq. 2)then
                    call tfsortql(itab(3),kl,av,n-2,mode,kf,irtc)
                    return
                  endif
                endif
              else
                itab(n)=-itab(ip2)
                itab(ip2)=i2
                call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
                return
              endif
            enddo
            if(ip2 .gt. ip1)then
              is=itab(ip1)
              itab(ip1)=itab(ip2)
              itab(ip2)=is
              ip1=ip1+1
              ip2=ip2-1
            endif
          enddo
        endif
      else
        l=itforderl(kl,av,i1,i2,kf,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(mode .eq. 0 .or. l .ne. 0)then
          if(l .gt. 0)then
            is=i1
            i1=i2
            i2=is
          endif
        else
          itab(1)=min(i1,i2)
          itab(n)=itab(1)-i1-i2
          call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
          return
        endif
        if(n .eq. 2)then
          itab(1)=i1
          itab(2)=i2
          return
        endif
        m=(n+1)/2
        im=itab(m)
        l=itforderl(kl,av,i1,im,kf,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(mode .eq. 0 .or. l .ne. 0)then
          if(l .lt. 0)then
            l=itforderl(kl,av,im,i2,kf,irtc)
            if(mode .eq. 0 .or. l .ne. 0)then
              if(l .gt. 0)then
                is=im
                im=i2
                i2=is
              endif
            else
              itab(1)=i1
              itab(m)=min(im,i2)
              itab(n)=itab(m)-im-i2
              call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
              return
            endif
          elseif(l .ne. 0)then
            is=im
            im=i1
            i1=is
          endif
        else
          itab(m)=min(i1,im)
          itab(1)=itab(m)-i1-im
          itab(n)=i2
          call tfsortql(itab(2),kl,av,n-1,mode,kf,irtc)
          return
        endif
        itab(1)=i1
        itab(n)=i2
        itab(m)=im
        if(n .eq. 3)then
          return
        endif
        itab(m)=itab(2)
        itab(2)=im
        ip1=3
        ip2=n-1
        if(mode .eq. 0)then
          do while(ip1 .le. ip2)
            l1=0
            do while(l1 .le. 0 .and. ip1 .le. ip2)
              l1=itforderl(kl,av,itab(ip1),im,kf,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(l1 .le. 0)then
                ip1=ip1+1
              endif
            enddo
            l2=0
            do while(l2 .le. 0 .and. ip1 .le. ip2)
              l2=itforderl(kl,av,im,itab(ip2),kf,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(l2 .le. 0)then
                ip2=ip2-1
              endif
            enddo
            if(ip2 .gt. ip1)then
              is=itab(ip1)
              itab(ip1)=itab(ip2)
              itab(ip2)=is
              ip1=ip1+1
              ip2=ip2-1
            endif
          enddo
        else
          do while(ip1 .le. ip2)
            l1=0
            do while(l1 .le. 0 .and. ip1 .le. ip2)
              l1=itforderl(kl,av,itab(ip1),im,kf,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(l1 .ne. 0)then
                if(l1 .le. 0)then
                  ip1=ip1+1
                endif
              else
                itab(2)=min(itab(ip1),im)
                itab(n)=itab(2)-itab(ip1)-im
                itab(ip1)=i2
                call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
                return
              endif
            enddo
            l2=0
            do while(l2 .le. 0 .and. ip1 .le. ip2)
              l2=itforderl(kl,av,im,itab(ip2),kf,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(l2 .ne. 0)then
                if(l2 .le. 0)then
                  ip2=ip2-1
                  if(ip2 .eq. 2)then
                    call tfsortql(itab(3),kl,av,n-2,mode,kf,irtc)
                    return
                  endif
                endif
              else
                itab(2)=min(im,itab(ip2))
                itab(n)=itab(2)-im-itab(ip2)
                itab(ip2)=i2
                call tfsortql(itab,kl,av,n-1,mode,kf,irtc)
                return
              endif
            enddo
            if(ip2 .gt. ip1)then
              is=itab(ip1)
              itab(ip1)=itab(ip2)
              itab(ip2)=is
              ip1=ip1+1
              ip2=ip2-1
            endif
          enddo
        endif
      endif
      ip1=ip1-1
      is=itab(ip1)
      itab(ip1)=im
      itab(2)=is
      call tfsortql(itab,kl,av,ip1-1,mode,kf,irtc)
      if(irtc .ne. 0)then
        return
      endif
      call tfsortql(itab(ip1+1),kl,av,n-ip1,mode,kf,irtc)
      return
      end

c     Merge sort core engine
c     mode = 0 Merge sort
c            1 Merge union
c            2 Merge union with dropped index list (for Override[])
      subroutine tfsortml(itab,kl,av,n,mode,kf,irtc)
      use tfstk
      use iso_c_binding
      implicit none
      type (sad_descriptor) kf
      type (sad_dlist) kl
      type (sad_rlist), pointer :: klr
      integer*4 n,mode,irtc
      integer*4 itab(n)
      integer*4 isp0,l,m,im,is,i1,i2,l1,l2,j1,j2,p0,p1,p2
      integer*4 itforderl
      logical*4 av

      irtc=0
      if(n .le. 1)then
        return
      endif

c     Special case: real list
      if(av .and. kf%k .eq. ktfref)then
        call dlist_rlist(kl,klr)
        call tfsortmrl(itab,klr,n,mode)
        return
      endif

c     Special case: n = 2
      if(n .eq. 2)then
        l=itforderl(kl,av,1,2,kf,irtc)
        if(irtc .ne. 0)return
        if(l .gt. 0)then
          itab(1)=2
          itab(2)=1
        elseif(mode .ne. 0 .and. l .eq. 0)then
          itab(2)=-2
        endif
        return
      endif

c     Initialize Array Index
      isp0=isp
      isp=isp+n

      p0=2
      p2=p0
      im=0
      i1=1
      do while(i1 .le. n)
        im=im+1
        itab(im)=p2-p0+1

        p1=p2
        p2=p2+1
        itastk(p2,isp0)=i1

        i2=i1+1
        do while(i2 .le. n)
          l=itforderl(kl,av,i2-1,i2,kf,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
          if(l .gt. 0)then
            if(p2-p1 .lt. 2)then
              itastk(p2,isp0)=i2
              p2=p2+1
              itastk(p2,isp0)=i1
              i1=i2+1
            else
              i1=i2
            endif
            exit
          elseif(l .lt. 0)then
            p2=p2+1
            itastk(p2,isp0)=i2
          elseif(mode .eq. 0)then
            p2=p2+1
            itastk(p2,isp0)=i2
          endif
          i2=i2+1
        enddo
        if(i2 .gt. n)exit
      enddo
      itab(im+1)=p2-p0+1
      m=im

c     Merge
      p1=2
      isp=isp0+p2-p0
      do while(m .gt. 1)
        p0=p2
        im=0
        is=1
        do while(is+1 .le. m)
          i1=itab(is)
          i2=itab(is+1)
          l1=i2
          l2=itab(is+2)
          is=is+2

          im=im+1
          itab(im)=p2-p0+1
          do
            j1=itastk(p1+i1,isp0)
            j2=itastk(p1+i2,isp0)

            l=itforderl(kl,av,j1,j2,kf,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
            if(l .lt. 0)then
              p2=p2+1
              itastk(p2,isp0)=j1
              i1=i1+1
              if(i1 .ge. l1)then
                do while(i2 .lt. l2)
                  p2=p2+1
                  itastk(p2,isp0)=itastk(p1+i2,isp0)
                  i2=i2+1
                enddo
                exit
              endif
            elseif(l .gt. 0)then
              p2=p2+1
              itastk(p2,isp0)=j2
              i2=i2+1
              if(i2 .ge. l2)then
                do while(i1 .lt. l1)
                  p2=p2+1
                  itastk(p2,isp0)=itastk(p1+i1,isp0)
                  i1=i1+1
                enddo
                exit
              endif
            else
              if(mode .eq. 0)then
                p2=p2+1
                itastk(p2,isp0)=j1
                i1=i1+1
                if(i1 .ge. l1)then
                  do while(i2 .lt. l2)
                    p2=p2+1
                    itastk(p2,isp0)=itastk(p1+i2,isp0)
                    i2=i2+1
                  enddo
                  exit
                endif
              else
                i2=i2+1
                if(i2 .ge. l2)then
                  do while(i1 .lt. l1)
                    p2=p2+1
                    itastk(p2,isp0)=itastk(p1+i1,isp0)
                    i1=i1+1
                  enddo
                  exit
                endif
              endif
            endif
          enddo
        enddo
        if(is .eq. m)then
          im=im+1
          itab(im)=p2-p0+1
          do i1=itab(is),itab(is+1)-1
            p2=p2+1
            itastk(p2,isp0)=itastk(p1+i1,isp0)
          enddo
        endif
        itab(im+1)=p2-p0+1
        m=im
        p2=p1
        p1=p0
      enddo

c     Copy index table to return area
      l1=itab(2)
c      do i1=1,l1-1
      itab(1:l1-1)=itastk(p1+1:p1+l1-1,isp0)
c      enddo
      if(mode .ne. 2)then
        itab(l1:n)=-1
c        do i1=l1,n
c          itab(i1)=-1
c        enddo
      else
c        do im=1,n
        itastk(1:n,isp0+1)=-[(im,im=1,n)]
c        enddo
c        do i1=1,l1-1
        itastk(itab(1:l1-1),isp0+1)=0
c        enddo
        i1=l1
        do im=1,n
          i2=itastk(im,isp0+1)
          if(i2 .lt. 0)then
            itab(i1)=i2
            i1=i1+1
          endif
        enddo
      endif

      isp=isp0
      return
      end

c     Special version of tfsortml for Real List
      subroutine tfsortmrl(itab,kl,n,mode)
      use tfstk
      implicit none
      type (sad_rlist) kl
      integer*4 n,mode
      integer*4 itab(n)
      integer*4 isp0,m,im,is,i1,i2,l1,l2,j1,j2,p0,p1,p2
      real*8 v1,v2

c     Special case: n = 2
      if(n .eq. 2)then
        v1=kl%rbody(1)
        v2=kl%rbody(2)
        if(v1 .gt. v2)then
          itab(1)=2
          itab(2)=1
        elseif(mode .ne. 0 .and. v1 .eq. v2)then
          itab(2)=-2
        endif
        return
      endif

c     Initialize Array Index
      isp0=isp

      p0=2
      p2=p0
      im=0
      i1=1
      do while(i1 .le. n)
        im=im+1
        itab(im)=p2-p0+1

        p1=p2
        p2=p2+1
        itastk(p2,isp0)=i1

        i2=i1+1
        do while(i2 .le. n)
          v1=kl%rbody(i2-1)
          v2=kl%rbody(i2)
          if(v1 .gt. v2)then
            if(p2-p1  .lt. 2)then
              itastk(p2,isp0)=i2
              p2=p2+1
              itastk(p2,isp0)=i1
              i1=i2+1
            else
              i1=i2
            endif
            exit
          elseif(v1 .eq. v2)then
            if(mode .eq. 0)then
              p2=p2+1
              itastk(p2,isp0)=i2
            endif
          else
            p2=p2+1
            itastk(p2,isp0)=i2
          endif
          i2=i2+1
        enddo
        if(i2 .gt. n)exit
      enddo
      itab(im+1)=p2-p0+1
      m=im

c     Merge
      p1=2
      do while(m .gt. 1)
        p0=p2
        im=0
        is=1
        do while(is+1 .le. m)
          i1=itab(is)
          i2=itab(is+1)
          l1=i2
          l2=itab(is+2)
          is=is+2

          im=im+1
          itab(im)=p2-p0+1
          do
            j1=itastk(p1+i1,isp0)
            j2=itastk(p1+i2,isp0)
            v1=kl%rbody(j1)
            v2=kl%rbody(j2)
            if(v1 .gt. v2)then
              p2=p2+1
              itastk(p2,isp0)=j2
              i2=i2+1
              if(i2 .ge. l2)then
                do while(i1 .lt. l1)
                  p2=p2+1
                  itastk(p2,isp0)=itastk(p1+i1,isp0)
                  i1=i1+1
                enddo
                exit
              endif
            elseif(v1 .eq. v2)then
              if(mode .eq. 0)then
                p2=p2+1
                itastk(p2,isp0)=j1
                i1=i1+1
                if(i1 .ge. l1)then
                  do while(i2 .lt. l2)
                    p2=p2+1
                    itastk(p2,isp0)=itastk(p1+i2,isp0)
                    i2=i2+1
                  enddo
                  exit
                endif
              else
                i2=i2+1
                if(i2 .ge. l2)then
                  do while(i1 .lt. l1)
                    p2=p2+1
                    itastk(p2,isp0)=itastk(p1+i1,isp0)
                    i1=i1+1
                  enddo
                  exit
                endif
              endif
            else
              p2=p2+1
              itastk(p2,isp0)=j1
              i1=i1+1
              if(i1 .ge. l1)then
                do while(i2 .lt. l2)
                  p2=p2+1
                  itastk(p2,isp0)=itastk(p1+i2,isp0)
                  i2=i2+1
                enddo
                exit
              endif
            endif
          enddo
        enddo
        if(is .eq. m)then
          im=im+1
          itab(im)=p2-p0+1
          do i1=itab(is),itab(is+1)-1
            p2=p2+1
            itastk(p2,isp0)=itastk(p1+i1,isp0)
          enddo
        endif
        itab(im+1)=p2-p0+1
        m=im
        p2=p1
        p1=p0
      enddo

c     Copy index table to return area
      l1=itab(2)
c      do i1=1,l1-1
      itab(1:l1-1)=itastk(p1+1:p1+l1-1,isp0)
c      enddo
      if(mode .ne. 2)then
        itab(l1:n)=-1
c        do i1=l1,n
c          itab(i1)=-1
c        enddo
      else
c        do im=1,n
        itastk(1:n,isp0+1)=-[(im,im=1,n)]
c        enddo
c        do i1=1,l1-1
        itastk(itab(1:l1-1),isp0+1)=0
c        enddo
        i1=l1
        do im=1,n
          i2=itastk(im,isp0+1)
          if(i2 .lt. 0)then
            itab(i1)=i2
            i1=i1+1
          endif
        enddo
      endif

      isp=isp0
      return
      end

      recursive subroutine tfintersection(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kj,kr
      type (sad_dlist), pointer :: kl1,klr
      integer*8 ka1,kax
      integer*4 isp1,mode,irtc,i,j,kk,isp0,
     $     isp2,itfcanonicalorder,il,ih,m,narg,itfmessage
      real*8 vj
      narg=isp-isp1
      if(narg .eq. 1)then
        call tfsort(isp1,kx,1,irtc)
        return
      endif
      do i=isp1+1,isp
        if(ktfnonlistq(ktastk(i)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
          return
        endif
      enddo
      ka1=ktfaddr(ktastk(isp1+1))
      call loc_sad(ka1,kl1)
      m=kl1%nl
      kx%k=ktflist+ka1
      if(m .eq. 0)then
        irtc=0
        return
      endif
      isp0=isp
      if(narg .gt. 2)then
        do i=isp1+2,isp
          isp=isp+1
          dtastk(isp)=kx
          isp=isp+1
          ktastk(isp)=ktastk(i)
          call tfintersection(isp0,kx,mode,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            return
          endif
        enddo
        return
      endif
      if(.not. tfsameheadq(ka1,ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::samehead',' ')
        return
      endif
      isp=isp+1
      ktastk(isp)=ktastk(isp0)
      call tfsort(isp0,kr,1,irtc)
      isp=isp0
      if(irtc .ne. 0)then
        return
      endif
      call loc_sad(ktfaddrd(kr),klr)
      call tfgetllstkall(klr)
      isp2=isp
      if(mode .eq. 0)then
        if(ktfreallistq(kl1))then
          LOOP_J_1: do j=1,m
            vj=kl1%rbody(j)
            il=isp0+1
            if(ktfrealq(ktastk(il)))then
              if(vj .eq. rtastk(il))then
                go to 110
              elseif(vj .lt. rtastk(il))then
                cycle LOOP_J_1
              endif
            else
              cycle LOOP_J_1
            endif
            ih=isp2
            if(ktfrealq(ktastk(ih)))then
              if(vj .eq. rtastk(ih))then
                go to 110
              elseif(vj .gt. rtastk(ih))then
                cycle LOOP_J_1
              endif
            endif
            do while(ih .gt. il+1)
              kk=il+(ih-il)/2
              if(ktfrealq(ktastk(kk)))then
                if(vj .eq. rtastk(kk))then
                  go to 110
                elseif(vj .lt. rtastk(kk))then
                  ih=kk
                else
                  il=kk
                endif
              else
                il=kk
              endif
            enddo
            cycle LOOP_J_1
 110        isp=isp+1
            rtastk(isp)=vj
          enddo LOOP_J_1
        else
          LOOP_J_2: do j=1,m
            kj=kl1%dbody(j)
            il=isp0+1
            i=itfcanonicalorder(kj,dtastk(il))
            if(i .eq. 0)then
              go to 210
            elseif(i .lt. 0)then
              cycle LOOP_J_2
            endif
            ih=isp2
            i=itfcanonicalorder(kj,dtastk(ih))
            if(i .eq. 0)then
              go to 210
            elseif(i .gt. 0)then
              cycle LOOP_J_2
            endif
            do while(ih .gt. il+1)
              kk=il+(ih-il)/2
              i=itfcanonicalorder(kj,dtastk(kk))
              if(i .eq. 0)then
                go to 210
              elseif(i .lt. 0)then
                ih=kk
              else
                il=kk
              endif
            enddo
            cycle LOOP_J_2
 210        isp=isp+1
            dtastk(isp)=kj
          enddo LOOP_J_2
        endif
      elseif(mode .eq. 1)then
        if(ktfreallistq(kl1))then
          LOOP_J_3: do j=1,m
            vj=kl1%rbody(j)
            il=isp0+1
            if(ktfrealq(ktastk(il)))then
              if(vj .eq. rtastk(il))then
                cycle LOOP_J_3
              elseif(vj .lt. rtastk(il))then
                go to 310
              endif
            else
              go to 310
            endif
            ih=isp2
            if(ktfrealq(ktastk(ih)))then
              if(vj .eq. rtastk(ih))then
                cycle LOOP_J_3
              elseif(vj .gt. rtastk(ih))then
                go to 310
              endif
            endif
            do while(ih .gt. il+1)
              kk=il+(ih-il)/2
              if(ktfrealq(ktastk(kk)))then
                if(vj .eq. rtastk(kk))then
                  cycle LOOP_J_3
                elseif(vj .lt. rtastk(kk))then
                  ih=kk
                else
                  il=kk
                endif
              else
                il=kk
              endif
            enddo
 310        isp=isp+1
            rtastk(isp)=vj
          enddo LOOP_J_3
        else
          LOOP_J_4: do j=1,m
            kj=kl1%dbody(j)
            il=isp0+1
            i=itfcanonicalorder(kj,dtastk(il))
            if(i .eq. 0)then
              cycle LOOP_J_4
            elseif(i .lt. 0)then
              go to 410
            endif
            ih=isp2
            i=itfcanonicalorder(kj,dtastk(ih))
            if(i .eq. 0)then
              cycle LOOP_J_4
            elseif(i .gt. 0)then
              go to 410
            endif
            do while(ih .gt. il+1)
              kk=il+(ih-il)/2
              i=itfcanonicalorder(kj,dtastk(kk))
              if(i .eq. 0)then
                cycle LOOP_J_4
              elseif(i .lt. 0)then
                ih=kk
              else
                il=kk
              endif
            enddo
 410        isp=isp+1
            dtastk(isp)=kj
          enddo LOOP_J_4
        endif
      endif
      kax=ktfmakelist(isp2)
      dlist(kax)=dtfcopy(kl1%head)
      kx%k=ktflist+kax
      isp=isp0+1
      dtastk(isp)=kx
      call tfsort(isp0,kx,1,irtc)
      isp=isp0
      return
      end

      integer function itforderl(kl,av,i1,i2,kf,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) k1,k2,kf,kx,kx1
      type (sad_dlist) kl
      integer*4 irtc,itfmessage,itfcanonicalorder,isp1,i1,i2
      real*8 v1,v2
      logical*4 av
      irtc=0
      if(av)then
        if(kf%k .eq. ktfref)then
          v1=kl%rbody(i1)
          v2=kl%rbody(i2)
          itforderl=merge(1,merge(0,-1,v1 .eq. v2),v1 .gt. v2)
          return
        endif
      endif
      k1=kl%dbody(i1)
      k2=kl%dbody(i2)
      if(kf%k .eq. ktfref)then
        itforderl=itfcanonicalorder(k1,k2)
        return
c        write(*,*)'itforderl ',i1,i2,itforderl
      else
        isp=isp+1
        isp1=isp
        dtastk(isp)=kf
        isp=isp+1
        dtastk(isp)=k1
        isp=isp+1
        dtastk(isp)=k2
        kx=tfefunref(isp1,.true.,irtc)
        if(irtc .ne. 0)then
          itforderl=-1
          return
        endif
        if(ktfnonrealq(kx))then
          itforderl=-1
          irtc=itfmessage(9,'General::wrongval',
     $         '"Real number","as the result of order-function"')
          return
        endif
        isp=isp1
        dtastk(isp)=kf
        isp=isp+1
        dtastk(isp)=k2
        isp=isp+1
        dtastk(isp)=k1
        kx1=tfefunref(isp1,.true.,irtc)
        isp=isp1-1
        if(irtc .ne. 0)then
          itforderl=-1
          return
        endif
        if(ktfnonrealq(kx1))then
          irtc=itfmessage(9,'General::wrongval',
     $         '"Real number","as the result of order-function"')
          itforderl=-1
          return
        endif
        itforderl=merge(0,merge(1,-1,kx%k .eq. 0),kx%k .eq. kx1%k)
      endif
      return
      end

      subroutine tforder(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfcanonicalorder,itfmessage
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      kx%x(1)=dble(itfcanonicalorder(dtastk(isp),dtastk(isp1+1)))
      irtc=0
      return
      end

      recursive integer*4 function itfcanonicalorder(k1,k2)
     $     result(ix)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) k1,k2,k1c,k2c
      type (sad_dlist), pointer :: kl1c,kl2c
      type (sad_symbol), pointer :: sym1c,sym2c
      type (sad_string), pointer :: str1,str2
      type (sad_namtbl), pointer :: loc1,loc2
      type (sad_pat), pointer :: pat1,pat2
      type (sad_complex), pointer :: cx1,cx2
      integer*8 icont1,icont2
      integer*4 m1,m2,l,i,itfstringorder,itfpatorder
      real*8 d,v1,v2
      ix=-1
      if(ktfrealq(k1,v1))then
        if(ktfrealq(k2,v2))then
          if(v1 .gt. v2)then
            ix=1
          elseif(v1 .eq. v2)then
            ix=0
          endif
        endif
        return
      elseif(ktfrealq(k2))then
        ix=1
        return
      endif
      if(k1%k .eq. k2%k)then
        ix=0
        return
      endif
      if(tfcomplexq(k1,cx1))then
        if(tfcomplexq(k2,cx2))then
          d=cx1%re-cx2%re
          if(d .gt. 0.d0)then
            ix=1
          elseif(d .eq. 0.d0)then
            d=cx1%im-cx2%im
            if(d .gt. 0.d0)then
              ix=1
            elseif(d .eq. 0.d0)then
              ix=0
            endif
          endif
        endif
        return
      elseif(tfcomplexq(k2))then
        ix=1
        return
      endif
      if(ktfstringq(k1,str1))then
        if(ktfstringq(k2,str2))then
          ix=itfstringorder(str1,str2)
        elseif(ktfrealq(k2))then
          ix=1
        endif
        return
      elseif(ktfstringq(k2))then
        ix=1
        return
      endif
      call tfcanonicalconv(k1,k1c)
      call tfcanonicalconv(k2,k2c)
      if(ktfsymbolq(k1c,sym1c))then
        call loc_namtbl(sym1c%loc,loc1)
        if(ktfsymbolq(k2c,sym2c))then
          call loc_namtbl(sym2c%loc,loc2)
          icont1=loc1%cont
          icont2=loc2%cont
          if(icont1 .eq. icont2)then
            ix=itfstringorder(loc1%str,loc2%str)
            if(ix .eq. 0)then
              if(sym1c%gen .gt. sym2c%gen)then
                ix=1
              elseif(sym1c%gen .lt. sym2c%gen)then
                ix=-1
              endif
            endif
          else
            ix=itfcanonicalorder(dlist(icont1),dlist(icont2))
          endif
        endif
      elseif(ktfsymbolq(k2c))then
        ix=1
      elseif(ktflistq(k1c,kl1c))then
        if(ktflistq(k2c,kl2c))then
          if(kl1c%head%k .ne. kl2c%head%k)then
            l=itfcanonicalorder(kl1c%head,kl2c%head)
            if(l .ne. 0)then
              ix=l
              return
            endif
          endif
          m1=kl1c%nl
          m2=kl2c%nl
          if(m1 .gt. m2)then
            ix=1
          elseif(m1 .eq. m2)then
            if(ktfreallistq(kl1c))then
              if(ktfreallistq(kl2c))then
                do i=1,m1
                  if(kl1c%rbody(i) .gt. kl2c%rbody(i))then
                    ix=1
                    return
                  elseif(kl1c%rbody(i) .lt. kl2c%rbody(i))then
                    ix=-1
                    return
                  endif
                enddo
                ix=0
              endif
            else
              if(ktfreallistq(kl2c))then
                ix=1
              else
                do i=1,m1
                  l=itfcanonicalorder(kl1c%dbody(i),kl2c%dbody(i))
                  if(l .ne. 0)then
                    ix=l
                    return
                  endif
                enddo
                ix=0
              endif
            endif
          endif
        else
          ix=1
        endif
      elseif(ktflistq(k2c))then
        ix=1
      elseif(ktfpatq(k1c,pat1))then
        if(ktfpatq(k2c,pat2))then
          ix=itfpatorder(pat1,pat2)
        else
          ix=1
        endif
      elseif(ktfpatq(k2c))then
        ix=1
      else
        ix=merge(1,-1,ktftype(k1c%k) .gt. ktftype(k2c%k))
      endif
      return
      end

      subroutine tfcanonicalconv(k,kx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      integer*8 ka
      kx%k=merge(ktfsymbol+klist(ifunbase+ka),k%k,ktfoperq(k,ka))
      return
      end

      integer*4 function itfstringorder(str1,str2)
      use tfstk
      use tftok
      implicit none
      type (sad_string) str1,str2
      integer*4 m1,m2,i,ic1,ic2
      m1=str1%nch
      m2=str2%nch
      do i=1,min(m1,m2)
        ic1=ichar(str1%str(i:i))
        ic2=ichar(str2%str(i:i))
        if(ic1 .ne. ic2)then
          itfstringorder=merge(1,-1,ichorder(ic1) .gt. ichorder(ic2))
          return
        endif
      enddo
      itfstringorder=min(1,max(-1,m1-m2))
      return
      end

      integer*4 function itfcharacterorder(ch1,ch2)
      use tftok
      implicit none
      character ch1,ch2
      integer*4 ic1,ic2
      ic1=ichorder(ichar(ch1))
      ic2=ichorder(ichar(ch2))
      itfcharacterorder=max(-1,min(1,ic1-ic2))
      return
      end

      integer*4 function itfpatorder(pat1,pat2)
      use tfstk
      use tfcode
      implicit none
      type (sad_pat) pat1,pat2
      type (sad_string), pointer :: str1,str2
      integer*8 ih1,ih2
      integer*4 itfstringorder,itfcanonicalorder
      ih1=pat1%sym%loc
      ih2=pat2%sym%loc
      if(ih1 .eq. 0)then
        if(ih2 .ne. 0)then
          itfpatorder=-1
          return
        endif
      elseif(ih2 .eq. 0)then
        itfpatorder=1
        return
      endif
      call loc_symstr(ih1,str1)
      call loc_symstr(ih2,str2)
      itfpatorder=itfstringorder(str1,str2)
      if(itfpatorder .ne. 0)then
        return
      endif
      itfpatorder=itfcanonicalorder(pat1%expr,pat2%expr)
      if(itfpatorder .ne. 0)then
        return
      endif
      itfpatorder=itfcanonicalorder(pat1%head,pat2%head)
      if(itfpatorder .ne. 0)then
        return
      endif
      itfpatorder=itfcanonicalorder(pat1%default,pat2%default)
      return
      end
