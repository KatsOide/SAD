      integer*8 function ktfmaloc(k,n,m,vec,trans,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*8 ktfmalocp
      integer*4 ,intent(inout):: n,m,irtc
      logical*4 ,intent(in):: vec,trans
      ktfmaloc=ktfmalocp(k,n,m,vec,trans,.false.,.true.,irtc)
      return
      end

      integer*8 function ktfmalocp(k,n,m,vec,trans,map,err,irtc)
      use tfstk
      use tfshare
      use tmacro
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer :: kl,kl1,kli
      integer*8 kap,i0,ip0
      integer*4 ,intent(out):: n,m,irtc
      integer*4 i,itfmessage
      logical*4 ,intent(in):: vec,trans,map,err
c
      ktfmalocp=-9999
c
      if(.not. tflistq(k,kl))then
        go to 9100
      endif
      if(vec)then
        if(ktfreallistq(kl))then
          n=kl%nl
          m=0
          if(map .and. nparallel .gt. 1)then
            irtc=1
            kap=ktfallocshared(n)
c            kap=mapalloc8(rlist(0), n, 8, irtc)
          else
            kap=ktaloc(n)
          endif
          rlist(kap:kap+n-1)=kl%rbody(1:n)
c          call tmov(rlist(ka+1),rlist(kap),n)
          ktfmalocp=kap
          irtc=0
          return
        elseif(tfnumberq(kl%dbody(1)))then
          go to 9000
        endif
      elseif(ktfreallistq(kl)
     $       .or. tfnonlistq(kl%dbody(1),kl1))then
        go to 9000
      endif
      irtc=0
      n=kl%nl
      m=kl1%nl
      do i=1,n
        if(tfnonlistq(kl%dbody(i),kli))then
          go to 9000
        endif
        if(kli%nl .ne. m)then
          go to 9000
        endif
        if(ktfnonreallistqo(kli))then
          irtc=merge(itfmessage(9,'General::wrongtype','"Real matrix"'),
     $           -1,err)
          return
        endif
      enddo
      if(map .and. nparallel .gt. 1)then
        irtc=1
        kap=ktfallocshared(m*n)
c        kap=mapalloc8(rlist(0), m*n, 8, irtc)
      else
        kap=ktaloc(m*n)
      endif
      if(kap .lt. 0)then
        irtc=itfmessage(999,'General::memoryfull',' ')
        return
      endif
      if(trans)then
        do i=1,n
          i0=ktfaddr(kl%dbody(i)%k)
          ip0=kap+(i-1)*m-1
          rlist(ip0+1:ip0+m)=rlist(i0+1:i0+m)
        enddo
      else
        do i=1,n
          i0=ktfaddr(kl%dbody(i)%k)
          rlist(kap+i-1:kap+(m-1)*n+i-1:n)=rlist(i0+1:i0+m)
c          do j=1,m
c            ip=kap+(j-1)*n+i-1
c            rlist(ip)=rlist(i0+j)
c          enddo
        enddo
      endif
      irtc=0
      ktfmalocp=kap
      return
 9000 irtc=merge(itfmessage(9,'General::wrongtype','"Matrix"'),
     $     -1,err)
      return
 9100 irtc=merge(itfmessage(9,'General::wrongtype',
     $     '"Numerical List or Matrix"'),-1,err)
      return
      end

      subroutine tfl2m(kl,a,n,m,trans)
      use tfstk
      implicit none
      type (sad_dlist) ,intent(in):: kl
      type (sad_dlist), pointer :: kli
      integer*4 ,intent(in):: n,m
      integer*4 i
      real*8, intent(out):: a(n,m)
      logical*4 ,intent(in):: trans
      if(trans)then
        do i=1,m
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          a(:,i)=kli%rbody(1:n)
        enddo
      else
        do i=1,n
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          a(i,:)=kli%rbody(1:m)
        enddo
      endif
      return
      end

      subroutine tfl2cm(kl,c,n,m,trans,irtc)
      use tfstk
      implicit none
      type (sad_dlist) ,intent(in):: kl
      type (sad_dlist), pointer :: kli
      type (sad_complex), pointer :: cx
      integer*8 kij
      integer*4 ,intent(in)::n,m
      integer*4 ,intent(out)::irtc
      integer*4 i,j
      complex*16, intent(out):: c(n,max(m,1))
      real*8 x
      logical*4 trans
      if(m .eq. 0)then
        do j=1,n
          kij=kl%body(j)
          if(ktfrealq(kij,x))then
            c(j,1)=dcmplx(x,0.d0)
          elseif(tfcomplexq(kij,cx))then
            c(j,1)=cx%cx(1)
          else
            irtc=-1
            return
          endif
        enddo
      elseif(trans)then
        do i=1,m
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          do j=1,n
            kij=kli%dbody(j)%k
            if(ktfrealq(kij,x))then
              c(j,i)=x
            elseif(tfcomplexq(kij,cx))then
              c(j,i)=cx%cx(1)
            else
              irtc=-1
              return
            endif
          enddo
        enddo
      else
        do i=1,n
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          do j=1,m
            kij=kli%dbody(j)%k
            if(ktfrealq(kij,x))then
              c(i,j)=x
            elseif(tfcomplexq(kij,cx))then
              c(i,j)=cx%cx(1)
            else
              irtc=-1
              return
            endif
          enddo
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfmsize(k,n,m,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist), pointer ::kl,kl1,kli
      integer*4 ,intent(inout):: n,m,irtc
      integer*4 i,itfmessage
      if(.not. tflistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      if(ktfreallistq(kl))then
        n=kl%nl
        m=-1
        irtc=0
        return
      endif
      if(tfnonlistq(kl%dbody(1),kl1))then
        go to 9000
      endif
      n=kl%nl
      m=kl1%nl
      do i=1,n
        if(tfnonlistq(kl%dbody(i),kli))then
          go to 9000
        endif
        if(kli%nl .ne. m)then
          go to 9000
        endif
        if(ktfnonreallistqo(kli))then
          irtc=itfmessage(9,'General::wrongtype','"Real matrix"')
          return
        endif
      enddo
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Matrix"')
      return
      end

      integer*8 function ktfcmaloc(k,n,m,vec,trans,err,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_complex), pointer :: cxi
      type (sad_dlist), pointer :: kl,kl1,kli
      type (sad_rlist), pointer :: klj
      integer*8 kap,ki,ip0,kj,kaj
      integer*4 ,intent(out):: n,m,irtc
      integer*4 i,j,itfmessage
      logical*4 ,intent(in):: vec,trans,err
      ktfcmaloc=0
      if(.not. tflistq(k,kl))then
        go to 9100
      endif
      irtc=0
      if(vec)then
        n=kl%nl
        m=0
        if(ktfreallistq(kl))then
          kap=ktaloc(n*2)
          rlist(kap:kap+2*n-2:2)=kl%rbody(1:n)
          rlist(kap+1:kap+2*n-1:2)=0.d0
c          do i=1,n
c            rlist(kap+i*2-2)=kl%rbody(i)
c            rlist(kap+i*2-1)=0.d0
c          enddo
          ktfcmaloc=kap
          return
        elseif(tfnumberq(kl%dbody(1)))then
          kap=ktaloc(n*2)
          do i=1,n
            ki=kl%dbody(i)%k
            if(ktfrealq(ki))then
              klist(kap+i*2-2)=ki
              klist(kap+i*2-1)=0
            else
              if(tfcomplexq(ki,cxi))then
                rlist(kap+i*2-2)=cxi%re
                rlist(kap+i*2-1)=cxi%im
              else
                go to 8900
              endif
            endif
          enddo
          ktfcmaloc=kap
          return
        endif
      elseif(ktfreallistq(kl) .or. tfnumberq(kl%dbody(1)))then
        go to 9000
      endif
      if(tfnonlistq(kl%dbody(1),kl1))then
        ktfcmaloc=-1
        go to 9000
      endif
      irtc=0
      n=kl%nl
      m=kl1%nl
      do i=1,n
        if(tfnonlistq(kl%dbody(i)%k,kli))then
          ktfcmaloc=-1
          go to 9000
        endif
        if(kli%nl .ne. m)then
          ktfcmaloc=-1
          go to 9000
        endif
      enddo
      kap=ktaloc(m*n*2)
      if(trans)then
        do i=1,n
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          ip0=kap+(i-1)*2*m-2
          if(ktfreallistq(kli))then
            rlist(ip0+2:ip0+2*m:2)=kli%rbody(1:m)
            rlist(ip0+3:ip0+2*m+1:2)=0.d0
c            do j=1,m
c              rlist(ip0+j*2  )=kli%rbody(j)
c              rlist(ip0+j*2+1)=0.d0
c            enddo
          else
            do j=1,m
              kj=kli%dbody(j)%k
              if(ktfrealq(kj))then
                klist(ip0+j*2  )=kj
                klist(ip0+j*2+1)=0
              else
                kaj=ktfaddr(kj)
                if(ktfreallistq(kj,klj) .and.
     $               klj%head%k .eq. ktfoper+mtfcomplex
     $               .and. klj%nl .eq. 2)then
                  rlist(ip0+j*2  )=klj%rbody(1)
                  rlist(ip0+j*2+1)=klj%rbody(2)
                else
                  call tfree(kap)
                  ktfcmaloc=-1
                  go to 9000
                endif
              endif
            enddo
          endif
        enddo
      else
        do i=1,n
          call loc_sad(ktfaddr(kl%dbody(i)%k),kli)
          if(ktfreallistq(kli))then
            rlist(kap+i*2-2:kap+2*((m-1)*n+i-1))=kli%rbody(1:m)
c            do j=1,m
c              ip0=kap+(j-1)*2*n+i*2-2
c              rlist(ip0  )=kli%rbody(j)
c              rlist(ip0+1)=0.d0
c            enddo
          else
            do j=1,m
              kj=kli%dbody(j)%k
              ip0=kap+(j-1)*2*n+i*2-2
              if(ktfrealq(kj))then
                klist(ip0  )=kj
                klist(ip0+1)=0
              else
                kaj=ktfaddr(kj)
                if(ktfreallistq(kj,klj)  .and.
     $               klj%head%k .eq. ktfoper+mtfcomplex
     $               .and. klj%nl .eq. 2)then
                  rlist(ip0  )=rlist(kaj+1)
                  rlist(ip0+1)=rlist(kaj+2)
                else
                  call tfree(kap)
                  ktfcmaloc=-1
                  go to 9000
                endif
              endif
            enddo
          endif
        enddo
      endif
      irtc=0
      ktfcmaloc=kap
      return
 8900 call tfree(kap)
      ktfcmaloc=-1
      if(err)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Numeric vector"')
      else
        irtc=-1
      endif
      return
 9000 irtc=merge(itfmessage(9,'General::wrongtype','"Matrix"'),
     $     -1,err)
      return
 9100 irtc=merge(itfmessage(9,'General::wrongtype','"List"'),
     $     -1,err)
      return
      end

      integer*8 function ktfc2l(cx,n)
      use tfstk
      implicit none
      integer*8 ktfcm2l
      integer*4 ,intent(in)::n
      complex*16 ,intent(in):: cx(n)
      ktfc2l=ktfcm2l(cx,0,n,1,.false.,.false.)
      return
      end

      subroutine nanchecka(a,n,m,ndim,tag,k)
      use tfstk, only:ktfenanq
      implicit none
      integer*4 n,m,ndim,i,j,k
      real*8, intent(in):: a(ndim,m)
      character*(*) tag
      do i=1,m
        do j=1,n
          if(ktfenanq(a(j,i)))then
            write(*,'(2a,3i6,Z18)')'nanchecka ',tag,k,j,i,a(j,i)
            return
          endif
        enddo
      enddo
      return
      end

      integer*8 function ktfcm2l(a,n,m,nd,trans,conj)
      use tfstk
      implicit none
      type (sad_dlist), pointer :: klx,klxi
      type (sad_rlist), pointer :: klj,kl
      integer*8 kax,kaxi,ktcalocm,kai,kc,kaj
      integer*4,intent(in):: n,m,nd
      integer*4 i,j
      logical*4 ,intent(in)::trans,conj
      logical*4 c
      real*8 imag_sign
      complex*16, intent(in):: a(nd,m)
      imag_sign=merge(-1.d0,1.d0,conj)
      kc=0
      if(n .eq. 0)then
        kax=ktaaloc(-1,m,klx)
        c=.false.
        do i=1,m
          if(imag(a(1,i)) .eq. 0.d0)then
            klx%rbody(i)=dble(a(1,i))
            if(kc .ne. 0)then
              kai=kc+(i-1)*6
              call tflocal1(kai)
            endif
          else
            if(kc .eq. 0)then
              kc=ktcalocm(m-i+1)-(i-1)*6
            endif
            kai=kc+(i-1)*6
            call loc_sad(kai,kl)
            kl%rbody(1)=dble(a(1,i))
            kl%rbody(2)=imag_sign*imag(a(1,i))
            klx%dbody(i)%k=ktflist+kai
c            klist(kax+i)=ktflist+ktcalocv(0,dble(a(1,i)),
c     $           imag_sign*imag(a(1,i)))
            c=.true.
          endif
        enddo
        klx%attr=merge(lconstlist+lnonreallist,lconstlist,c)
      else
        if(trans)then
          kax=ktadaloc(-1,m,klx)
          do i=1,m
            kaxi=ktaaloc(0,n,klxi)
            kc=0
            c=.false.
            do j=1,n
              if(imag(a(j,i)) .eq. 0.d0)then
                klxi%rbody(j)=dble(a(j,i))
                if(kc .ne. 0)then
                  kaj=kc+(j-1)*6
                  call tflocal1(kaj)
                endif
              else
                if(kc .eq. 0)then
                  kc=ktcalocm(n-j+1)-(j-1)*6
                endif
                kaj=kc+(j-1)*6
                call loc_sad(kaj,klj)
                klj%rbody(1)=dble(a(j,i))
                klj%rbody(2)=imag_sign*imag(a(j,i))
                klxi%dbody(j)%k=ktflist+kaj
                c=.true.
              endif
            enddo
            klxi%attr=merge(lconstlist+lnonreallist,lconstlist,c)
            klx%dbody(i)%k=ktflist+kaxi
          enddo
          klxi%attr=ior(klxi%attr,lconstlist)
        else
          kax=ktadaloc(-1,n,klx)
          do i=1,n
            kaxi=ktaaloc(0,m,klxi)
            kc=0
            c=.false.
            do j=1,m
              if(imag(a(i,j)) .eq. 0.d0)then
                klxi%rbody(j)=dble(a(i,j))
                if(kc .ne. 0)then
                  kaj=kc+(j-1)*6
                  call tflocal1(kaj)
                endif
              else
                if(kc .eq. 0)then
                  kc=ktcalocm(m-j+1)-(j-1)*6
                endif
                kaj=kc+(j-1)*6
                call loc_sad(kaj,klj)
                klj%rbody(1)=dble(a(i,j))
                klj%rbody(2)=imag_sign*imag(a(i,j))
                klxi%dbody(j)%k=ktflist+kaj
                c=.true.
              endif
            enddo
            klxi%attr=merge(lconstlist+lnonreallist,lconstlist,c)
            klx%dbody(i)%k=ktflist+kaxi
          enddo
          klxi%attr=ior(klxi%attr,lconstlist)
        endif
      endif
      ktfcm2l=kax
      return
      end
