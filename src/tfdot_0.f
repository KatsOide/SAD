c$Header: /SAD/cvsroot/oldsad/src/tfdot.f,v 1.44 2011/05/16 05:21:25 amorita Exp $
      subroutine tfdot(k1,k2,kx,irtc)
      implicit none 
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 k1,k2,kx,ktfcmaloc,ktfmalocp,ktcvaloc,ktfcm2l,
     $     ktaloc,ktfcopy,ktavaloc,ktaaloc,ka1,kax,k1i,
     $     kr,kx1,kx2,kxi,ktcalocv,ktfmakelist
      integer*4 irtc,ir,j,i,id1,n1,m1,n2,m2,itfmessage,isp0
      logical*4 c1,c2,tflistqk
      if(.not. tflistqk(k1) .or. .not. tflistqk(k2))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      ka1=ktfaddr(k1)
      n1=ilist(2,ka1-1)
      if(tflistqk(klist(ka1+1)))then
        isp0=isp
        do i=1,n1
          isp=isp+1
          call tfdot(klist(ka1+i),k2,ktastk(isp),irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
        enddo
        kx=ktflist+ktfmakelist(isp0)
        isp=isp0
        return
      endif
      if(ktfreallistq(ka1))then
        m1=0
        kx1=ka1+1
        c1=.false.
      else
        kx1=ktfcmaloc(k1,n1,m1,.true.,.false.,.false.,irtc)
        if(irtc .gt. 0)then
          return
        elseif(irtc .ne. 0)then
          call tfinner(k1,k2,kx,ktfoper+mtfplus,ktfoper+mtfmult,irtc)
          return
        endif
        c1=.true.
      endif
      kx2=ktfmalocp(k2,n2,m2,.true.,.false.,.false.,.false.,irtc)
      if(irtc .eq. 0)then
        c2=.false.
      elseif(irtc .gt. 0)then
        if(c1)then
          call tfree(kx1)
        endif
        return
      else
        kx2=ktfcmaloc(k2,n2,m2,.true.,.false.,.false.,irtc)
        if(irtc .ne. 0)then
          if(c1)then
            call tfree(kx1)
          endif
          if(irtc .lt. 0)then
            call tfinner(k1,k2,kx,ktfoper+mtfplus,ktfoper+mtfmult,irtc)
          endif
          return
        endif
        c2=.true.
      endif
      if(n1 .ne. n2)then
        irtc=itfmessage(9,'General::equalleng','"arguments"')
        go to 9000
      endif
      if(n1 .eq. 0)then
        kx=0
        go to 9000
      endif
      if(c1)then
        if(c2)then
          kr=ktaloc(2*max(m2,1))
          call tfdotcc(rlist(kx1),rlist(kx2),rlist(kr),n1,max(m2,1))
        else
          kr=ktaloc(2*max(m2,1))
          call tfdotcr(rlist(kx1),rlist(kx2),rlist(kr),n1,max(m2,1))
        endif
      else
        if(c2)then
          kr=ktaloc(2*max(m2,1))
          call tfdotrec(rlist(kx1),rlist(kx2),rlist(kr),n1,max(m2,1))
        elseif(m2 .eq. 0)then
          call tfdotrr(rlist(kx1),rlist(kx2),kx,n1,1)
        else
          kax=ktavaloc(-1,m2)
          call tfdotrr(rlist(kx1),rlist(kx2),rlist(kax+1),n1,m2)
          kx=ktflist+kax
        endif
      endif
      if(m2 .eq. 0)then
        if(c1 .or. c2)then
          if(rlist(kr+1) .eq. 0.d0)then
            kx=klist(kr)
          else
            kx=ktflist+ktcalocv(-1,rlist(kr),rlist(kr+1))
          endif
          call tfree(kr)
        endif
      else
        if(c1 .or. c2)then
          kx=ktflist+ktfcm2l(rlist(kr),0,m2,1,.false.,.false.)
          call tfree(kr)
        endif
      endif
 9000 if(c1)then
        call tfree(kx1)
      endif
      call tfree(kx2)
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfdotrr(a,b,c,n,m)
      implicit none
      integer*4 n,m,i,j
      real*8 a(n),b(n,m),c(m),s
      do i=1,m
        s=a(1)*b(1,i)
        do j=2,n
          s=s+a(j)*b(j,i)
        enddo
        c(i)=s
      enddo
      return
      end
        
      subroutine tfdotcr(a,b,c,n,m)
      implicit none
      integer*4 n,m,i,j
      complex*16 a(n),c(m),s
      real*8 b(n,m)
      do i=1,m
        s=a(1)*b(1,i)
        do j=2,n
          s=s+a(j)*b(j,i)
        enddo
        c(i)=s
      enddo
      return
      end
        
      subroutine tfdotrec(a,b,c,n,m)
      implicit none
      integer*4 n,m,i,j
      complex*16 b(n,m),c(m),s
      real*8 a(n)
      do i=1,m
        s=a(1)*b(1,i)
        do j=2,n
          s=s+a(j)*b(j,i)
        enddo
        c(i)=s
      enddo
      return
      end
        
      subroutine tfdotcc(a,b,c,n,m)
      implicit none
      integer*4 n,m,i,j
      complex*16 a(n),b(n,m),c(m),s
      do i=1,m
        s=a(1)*b(1,i)
        do j=2,n
          s=s+a(j)*b(j,i)
        enddo
        c(i)=s
      enddo
      return
      end

      subroutine tfinner(k1,k2,kx,ks,kp,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 k1,k2,kx,ks,kp,ktaaloc,kax,ki,k1i,ka1,ka2,
     $     ktfcopy,ktfcopy1
      integer*4 irtc,m1,i,j,m2,isp0,itfmessage
      logical*4 tflistqk
      if(ktfnonlistq(k1) .or. ktfnonlistq(k2))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      ka1=ktfaddr(k1)
      ka2=ktfaddr(k2)
      m1=ilist(2,ka1-1)
      if(ktfnonreallistq(k1) .and. tflistqk(klist(ka1)))then
        kax=ktaaloc(-1,m1)
        do i=1,m1
          k1i=klist(ka1+i)
          call tfinner(k1i,k2,ki,ks,kp,irtc)
          if(irtc .ne. 0)then
            if(ktfnonreallistq(kax))then
              do j=i,m1
                klist(kax+j)=ktfoper+mtfnull
              enddo
            else
              do j=i,m1
                klist(kax+j)=0
              enddo
            endif
            return
          endif
          if(ktfnonrealq(ki))then
            klist(kax+i)=ktfcopy1(ki)
            ilist(2,kax-3)=ior(ilist(2,kax-3),lnonreallist)
          else
            klist(kax+i)=ki
          endif
        enddo
        kx=ktflist+kax
        return
      endif
      m2=ilist(2,ka2-1) 
      if(m2 .ne. m1)then
        irtc=itfmessage(9,'General::equalleng','"arguments"')
        return
      endif
      isp=isp+3
      isp0=isp
      do i=1,m1
        ktastk(isp-2)=kp
        ktastk(isp-1)=klist(ka1+i)
        ktastk(isp  )=klist(ka2+i)
        call tfefunref(isp-2,ki,.true.,irtc)
        isp=isp0
        if(irtc .ne. 0)then
          isp=isp-3
          return
        endif
        if(i .eq. 1)then
          kx=ki
        else
          ktastk(isp-2)=ks
          ktastk(isp-1)=kx
          ktastk(isp  )=ki
          call tfefunref(isp-2,kx,.true.,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            isp=isp-3
            return
          endif
        endif
      enddo
      isp=isp-3
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfouter(isp1,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,kh,kai,ktfmakelist,ktfcopy
      integer*4 isp1,irtc,narg,i,mi,ispi,isp0,j,mn,itfmessage
      real*8 vx
      logical*4 tfsameheadqk
      narg=isp-isp1
      if(narg .lt. 3)then
        irtc=itfmessage(9,'General::narg','"3 or more"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      itastk(1,ivstkoffset+isp1+2)=1
      do i=isp1+3,isp
        if(ktfnonlistq(ktastk(i)) .or.
     $       .not. tfsameheadqk(ktastk(isp1+2),ktastk(i)))then
          irtc=itfmessage(9,'General::samehead',' ')
          return
        endif
        itastk(1,ivstkoffset+i)=1
      enddo
      kh=ktfcopy(klist(ktfaddr(ktastk(isp1+2))))
      mn=ilist(2,ktfaddr(ktastk(isp1+narg))-1)
 1    do i=1,mn
        isp0=isp+1
        ktastk(isp0)=ktastk(isp1+1)
        isp=isp0
        do j=2,narg-1
          isp=isp+1
          ktastk(isp)=klist(ktfaddr(ktastk(isp1+j))+
     $         ilist(1,ivstkoffset+isp1+j))
        enddo
        isp=isp+1
        ktastk(isp)=klist(ktfaddr(ktastk(isp1+narg))+i)
        call tfefunrefstk(isp0,isp0,irtc)
        if(irtc .ne. 0)then
          call tflocal(kh)
          return
        endif
      enddo
      itastk(1,ivstkoffset+isp1+narg)=mn
      do i=narg,2,-1
        mi=ilist(2,ktfaddr(ktastk(isp1+i))-1)
        if(itastk(1,ivstkoffset+isp1+i) .ge. mi)then
          ispi=isp-mi
          kai=ktfmakelist(ispi)
          isp=ispi+1
          ktastk(isp)=ktflist+kai
          klist(kai)=ktfcopy(kh)
          itastk(1,ivstkoffset+isp1+i)=1
        else
          itastk(1,ivstkoffset+isp1+i)=itastk(1,ivstkoffset+isp1+i)+1
          go to 1
        endif
      enddo
      kx=ktastk(isp)
      isp=isp1+narg
      call tflocal(kh)
      irtc=0
      return
      include 'inc/TFSF.inc'
      end

      subroutine tftranspose(k,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 k,kx,ktadaloc,ktfcopy,ktavaloc,ktaaloc,kij,kax,kaxi,
     $     k1,ka,ka1,ki
      integer*4 irtc,m,n,i,j,itfmessage
      logical*4 tfnonlistqk,d
      if(tfnonlistqk(k))then
        go to 9000
      endif
      ka=ktfaddr(k)
      k1=klist(ka+1)
      if(tfnonlistqk(k1))then
        go to 9000
      endif
      ka1=ktfaddr(k1)
      m=ilist(2,ka -1)
      n=ilist(2,ka1-1)
      do i=2,m
        ki=klist(ka+i)
        if(tfnonlistqk(ki))then
          go to 9000
        endif
        if(ilist(2,ktfaddr(ki)-1) .ne. n)then
          irtc=itfmessage(9,'General::equalleng','"elements"')
          return
        endif
      enddo
      kax=ktadaloc(-1,n)
      do i=1,n
        kaxi=ktaaloc(0,m)
        d=.false.
        do j=1,m
          kij=klist(ktfaddr(klist(ka+j))+i)
          if(ktfrealq(kij))then
            klist(kaxi+j)=kij
          else
            d=.true.
            klist(kaxi+j)=ktfcopy(kij)
          endif
        enddo
        if(d)then
          ilist(2,kaxi-3)=ior(ilist(2,kaxi-3),lnonreallist)
        else
          ilist(2,kaxi-3)=lconstlist
        endif
        klist(kax+i)=ktflist+kaxi
      enddo
      kx=ktflist+kax
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Matrix"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfdet(isp1,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,k
      integer*4 isp1,irtc,n,m,narg,itfmessage
      logical*4 cmplm,realm,vec
      narg=isp-isp1
      if(narg .ne. 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      k=ktastk(isp)
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m)
      if(n .ne. m .or. m .eq. 0)then
        go to 9000
        return
      endif
      if(cmplm)then
        call tfcdet(k,kx,n,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      elseif(realm)then
        call tfrdet(k,kx,n)
      else
        go to 9000
      endif
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Numerical Square Matrix"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfrdet(k,kx,n)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,k,kfromr
      integer*4 n
      real*8 a(n,n),tdet
      call tfl2m(k,a,n,n,.false.)
      kx=kfromr(tdet(a,n,n))
      return
      end

      subroutine tfcdet(k,kx,n,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,k,ktcalocv,kfromr
      integer*4 n,irtc
      complex*16 c(n,n),cx,tcdet
      call tfl2cm(k,c,n,n,.false.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      cx=tcdet(c,n,n)
      if(imag(cx) .eq. 0.d0)then
        kx=kfromr(dble(cx))
      else
        kx=ktflist+ktcalocv(-1,dble(cx),imag(cx))
      endif
      return
      end

      subroutine tfsingularvalues(isp1,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,kax,k,kap,kwp,kux,kvx,kwx,ktadaloc,ktavaloc,
     $     ktfcopy1,ktfcmaloc,ktfmaloc,ktfcm2l,ktaloc,ktfm2l
      integer*4 isp1,irtc,mn,n,m,i,narg,itfmessage
      real*8 eps
      logical*4 inv
      logical*4 realm,cmplm,vec
      narg=isp-isp1
      if(narg .ne. 3)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      k=ktastk(isp-2)
      if(ktfnonrealq(ktastk(isp)) .or. ktfnonrealq(ktastk(isp-1)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real number for #2 and #3"')
        return
      endif
      eps=rtastk(isp-1)
      inv=rtastk(isp) .ne. 0.d0
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m)
      if(cmplm)then
        kap=ktfcmaloc(k,n,m,.false.,.false.,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kwp=ktaloc(m)
        if(max(n,m) .gt. 200)then
          call tftclupdate(3)
        endif
        call tcsvdma(rlist(kap),rlist(kwp),n,m,eps,inv,kux)
        mn=min(m,n)
        kvx=ktfcopy1(ktfcm2l(rlist(kap),mn,m,n,.false.,.false.))
      elseif(realm)then
        kap=ktfmaloc(k,n,m,.false.,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(max(n,m) .gt. 200)then
          call tftclupdate(3)
        endif
        kwp=ktaloc(m)
        call tsvdma(rlist(kap),rlist(kwp),m,n,eps,inv,kux)
        mn=min(m,n)
        kvx=ktfcopy1(ktfm2l(rlist(kap),mn,m,n,.false.))
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Matrix"')
        return
      endif
      kwx=ktavaloc(0,mn)
      do i=1,mn
        rlist(kwx+i)=rlist(kwp+i-1)
      enddo
      call tfree(kwp)
      call tfree(kap)
      kax=ktadaloc(-1,3)
      klist(kax+1)=ktflist+ktfcopy1(kux)
      klist(kax+2)=ktflist+kwx
      klist(kax+3)=ktflist+kvx
      kx=ktflist+kax
      return
      include 'inc/TFSF.inc'
      end

      subroutine tcsvdma(ca,w,m,n,eps,inv,kux)
      implicit none
      integer*8 kux,ktfcm2l
      integer*4 n,m,mn
      real*8 w(m),eps
      complex*16 ca(n,m),cu(n,n)
      logical*4 inv
      call tcsvdm(ca,cu,w,n,m,n,n,eps,inv)
      mn=min(m,n)
      kux=ktfcm2l(cu,mn,n,n,.false.,.false.)
      return
      end

      subroutine tsvdma(a,w,m,n,eps,inv,kux)
      implicit none
      integer*8 kux,ktfm2l
      integer*4 n,m,mn
      real*8 a(n,m),w(m),eps,u(n,n)
      logical*4 inv
      call tsvdm(a,u,w,n,m,n,n,eps,inv)
      mn=min(m,n)
      kux=ktfm2l(u,mn,n,n,.false.)
      return
      end

      subroutine tflinearsolve(isp1,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,k,kb
      integer*4 isp1,irtc,narg,nb,mb,n,m,itfmessage,mx
      real*8 eps
      logical*4 cmplm,realm,vec
      narg=isp-isp1
      if(narg .ne. 2 .and. narg .ne. 3)then
        irtc=itfmessage(9,'General::narg','"2 (+ option)"')
        return
      elseif(narg .eq. 2)then
        k=ktastk(isp-1)
        kb=ktastk(isp)
        eps=1.d-8
      else
        k=ktastk(isp-2)
        kb=ktastk(isp-1)
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype','"real number"')
          return
        endif
        eps=rtastk(isp)
      endif
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m)
      if(.not. realm)then
        go to 9000
      endif
      call tfmatrixmaybeq(kb,cmplm,realm,vec,mb,nb)
      if(.not. realm .and. .not. vec)then
        go to 9100
      endif
      if(nb .eq. 0 .or. vec)then
        nb=mb
        mb=0
      endif
      if(n .ne. nb)then
        irtc=itfmessage(9,'General::equalleng','"matrix and vector"')
        return
      endif
      if(max(n,m) .gt. 200)then
        call tftclupdate(3)
      endif
      mx=max(mb,1)
      call tflsolve(k,kb,kx,n,m,mb,mx,eps)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Real Matrix for #1"')
      return
 9100 irtc=itfmessage(9,'General::wrongtype',
     $     '"Real Vector or Matrix for #2"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tflsolve(k,kb,kx,n,m,mb,mx,eps)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*8 kx,k,kb,ktfm2l
      integer*4 n,m,mb,mx
      real*8 eps,a(n,m),b(n,mx),x(m,mx)
      call tfl2m(k,a,n,m,.false.)
      if(mb .eq. 0)then
        call tmov(rlist(kb+1),b,n)
      else
        call tfl2m(kb,b,n,mb,.true.)
      endif
      call tsolvm(a,b,x,n,m,mx,n,n,m,eps,.true.)
      if(mb .eq. 0)then
        kx=ktflist+ktfm2l(x,0,m,1,.false.)
      else
        kx=ktflist+ktfm2l(x,m,mb,m,.true.)
      endif
      return
      end

      subroutine tfdiagonalmatrix(k,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 k,kx,ktraaloc,ktadaloc,kax,kaxi,ka,ki,ktfcopy
      integer*4 irtc,m,i,itfmessage
      logical*4 tfnonlistqk
      if(tfnonlistqk(k))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      ka=ktfaddr(k)
      m=ilist(2,ka-1)
      kax=ktadaloc(-1,m)
      do i=1,m
        kaxi=ktraaloc(0,m)
        ki=klist(ka+i)
        if(ktfrealq(ki))then
          klist(kaxi+i)=ki
        else
          ilist(2,kaxi-3)=ior(lnonreallist,ilist(2,kaxi-3))
          klist(kaxi+i)=ktfcopy(ki)
        endif
        klist(kax+i)=ktflist+kaxi
      enddo
      kx=ktflist+kax
      irtc=0
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfidentitymatrix(k,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 k,kx,ktraaloc,ktadaloc,kax,kaxi
      integer*4 irtc,m,i,itfmessage
      real*8 rfromk
      if(ktfnonrealq(k))then
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
        return
      endif
      m=rfromk(k)
      if(m .le. 0)then
        irtc=itfmessage(9,'General::wrongnum','"Positive"')
        return
      endif
      kax=ktadaloc(-1,m)
      do i=1,m
        kaxi=ktraaloc(0,m)
        rlist(kaxi+i)=1.d0
        klist(kax+i)=ktflist+kaxi
      enddo
      kx=ktflist+kax
      irtc=0
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfeigensystem(k,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 k,kx,kax,kex,kvx,ktadaloc,ktfcopy1
      integer*4 irtc,n,m,itfmessage
      logical*4 realm,cmplm,vec
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m)
      if(n .ne. m)then
        go to 9000
      endif
      if(cmplm)then
        call tfceigen(k,n,kvx,kex,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      elseif(realm)then
        call tfreigen(k,n,kvx,kex)
      else
        go to 9000
      endif
      kax=ktadaloc(-1,2)
      klist(kax+1)=ktflist+ktfcopy1(kex)
      klist(kax+2)=ktflist+ktfcopy1(kvx)
      kx=ktflist+kax
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Square Matrix"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfceigen(k,m,kvx,kex,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 k,kvx,kex,ktfcm2l,ktfc2l
      complex*16 c(m,m),cw(m,m),ce(m)
      integer*4 irtc,m
      call tfl2cm(k,c,m,m,.false.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(m .gt. 200)then
        call tftclupdate(3)
      endif
      call tceigen(c,cw,ce,m,m)
      kvx=ktfcm2l(c,m,m,m,.true.,.false.)
      kex=ktfc2l(ce,m)
      return
      end

      subroutine tfreigen(k,m,kvx,kex)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 kvx,kex,ktfm2l,ktfcm2l,kvxkvr,kvxkvi,
     $     kzi1,kzi2,ktaaloc,ktcalocv,ktfc2l,k
      integer*4 m,i,j
      real*8 a(m,m),w(m,m)
      complex*16 ce(m)
      logical*4 d,con
      call tfl2m(k,a,m,m,.false.)
      if(m .gt. 200)then
        call tftclupdate(3)
      endif
      call teigen(a,w,ce,m,m)
      kvx=ktfm2l(a,m,m,m,.true.)
      kex=ktfc2l(ce,m)
      if(ktfnonreallistq(kex))then
        con=.true.
        do i=1,m
          if(con)then
            if(ktfnonrealq(klist(kex+i)))then
              kvxkvr=ktfaddr(klist(kvx+i  ))
              kvxkvi=ktfaddr(klist(kvx+i+1))
              kzi1=ktaaloc(0,m)
              kzi2=ktaaloc(0,m)
              d=.false.
              do j=1,m
                if(rlist(kvxkvi+j) .ne. 0.d0)then
                  d=.true.
                  klist(kzi1+j)=ktflist+
     $                 ktcalocv(0,rlist(kvxkvr+j),rlist(kvxkvi+j))
                  klist(kzi2+j)=ktflist+
     $                 ktcalocv(0,rlist(kvxkvr+j),-rlist(kvxkvi+j))
                else
                  rlist(kzi1+j)=rlist(kvxkvr+j)
                  rlist(kzi2+j)=rlist(kvxkvr+j)
                endif
              enddo
              if(d)then
                ilist(2,kzi1-3)=lconstlist+lnonreallist
                ilist(2,kzi2-3)=lconstlist+lnonreallist
              else
                ilist(2,kzi1-3)=lconstlist
                ilist(2,kzi2-3)=lconstlist
              endif
              call tflocal1(klist(kvx+i  ))
              call tflocal1(klist(kvx+i+1))
              klist(kvx+i  )=ktflist+kzi1
              klist(kvx+i+1)=ktflist+kzi2
              con=.false.
            endif
          else
            con=.true.
          endif
        enddo
      endif
      return
      include 'inc/TFSF.inc'
      end

      subroutine tffourier(inv,k,kx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/MACMATH.inc'
      integer*8 k,kx,ka
      integer*4 irtc,i,m,n,itfmessage
      logical*4 inv,tfcomplexqk,tfnonlistqk,even,power2
      if(tfnonlistqk(k))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      ka=ktfaddr(k)
      m=ilist(2,ka-1)
      irtc=-1
      if(m .eq. 0)then
        return
      elseif(m .eq. 1)then
        kx=k
        irtc=0
        return
      endif
      call tffft(ka,m,kx,inv,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Square Matrix"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tffft(ka,m,kx,inv,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/MACMATH.inc'
      integer*8 k,kx,ka,ktfc2l,kai
      integer*4 m,n,irtc,i
      complex*16 cx(m),cy(m)
      real*8 c,s,f,w,a(m*2),ay(m*2)
      logical*4 inv,power2,even,tfcomplexqk
      irtc=-1
      n=2
      do while(n .lt. m)
        n=n*2
      enddo
      power2=n .eq. m
      even=iand(m,1) .eq. 0
      if(ktfnonreallistq(ka))then
        f=1.d0/sqrt(dble(m))
        do i=1,m
          if(ktfrealq(klist(ka+i)))then
            cx(i)=f*rlist(ka+i)
          elseif(tfcomplexqk(klist(ka+i)))then
            kai=ktfaddr(klist(ka+i))
            cx(i)=f*dcmplx(rlist(kai+1),rlist(kai+2))
          else
            return
          endif
        enddo
        if(power2)then
          call tcftr(cx,m,inv)
        else
          if(inv)then
            call zfftw(cx,m,-1,cy)
          else
            call zfftw(cx,m,1,cy)
          endif
        endif
        kx=ktflist+ktfc2l(cx,m)
        irtc=0
        return
      elseif(even)then
        call tmov(rlist(ka+1),a,m)
        if(power2)then
          if(m .ge. 4)then
            call trftr(a,m,inv)
            f=1.d0/sqrt(dble(m))
            do i=1,m
              a(i)=f*a(i)
            enddo
            do i=1,m/2-1
              a((m-i)*2+1)= a(i*2+1)
              a((m-i)*2+2)=-a(i*2+2)
            enddo
            a(m+1)=a(2)
            a(m+2)=0.d0
            a(2  )=0.d0
            go to 1000
          else
            call tcftr(a,m/2,inv)
          endif
        else
          if(inv)then
            call zfftw(a,m/2,-1,ay)
          else
            call zfftw(a,m/2, 1,ay)
          endif
        endif
        f=.5d0/sqrt(dble(m))
        if(inv)then
          w=-2.d0*pi/m
        else
          w= 2.d0*pi/m
        endif
        do i=1,m/2-1
          c=cos(w*i)
          s=sin(w*i)
          a((m-i)*2+1)= f*(a(i*2+1)+a(m-i*2+1)
     $         +c*(a(i*2+2)+a(m-i*2+2))
     $         +s*(a(i*2+1)-a(m-i*2+1)))
          a((m-i)*2+2)=-f*(a(i*2+2)-a(m-i*2+2)
     $         +s*(a(i*2+2)+a(m-i*2+2))
     $         -c*(a(i*2+1)-a(m-i*2+1)))
        enddo
        a(m+1)=f*2.d0*(a(1)-a(2))
        a(1  )=f*2.d0*(a(1)+a(2))
        a(m+2)=0.d0
        a(2  )=0.d0
        do i=1,m/2-1
          a(i*2+1)= a((m-i)*2+1)
          a(i*2+2)=-a((m-i)*2+2)
        enddo
      else
        f=1.d0/sqrt(dble(m))
        do i=1,m
          a(i*2-1)=f*rlist(ka+i)
          a(i*2  )=0.d0
        enddo
        if(inv)then
          call zfftw(a,m,-1,ay)
        else
          call zfftw(a,m, 1,ay)
        endif
      endif
 1000 kx=ktflist+ktfc2l(a,m)
      irtc=0
      return
      include 'inc/TFSF.inc'
      end
c
c
c  K. Yokoya's FFT routines.   Received 2/26/1998.
c
c
C**************************ZFFTW **************************************
      SUBROUTINE ZFFTW(A,NN,ISGN,WORK)
C  FFT routine for Complex*16 array with arbitrary number of elements.
C  Needs work area but no restriction for NN.
C      B(k)= Sum A(j)* exp[+-2*pi*i*k*j/NN]   over j=0 to NN-1
C                     for k=0 to NN-1
C Input:
C  A     Complex*16 array of dimension A(0:NN-1)
C  NN    Length of A. (>=1)
C  ISGN  Sign of the phase, +1 or -1.
C Output:
C  A     B is stored in the same area as A.
C Work area
C  W     Same size as A.
C
      IMPLICIT NONE
      INTEGER MPRIME0,MPRIME,MREST
      PARAMETER (MPRIME0=172,MPRIME=500,MREST=MPRIME-MPRIME0)
      INTEGER NN,ISGN
      COMPLEX*16 A(0:NN-1),WORK(0:NN-1)
      INTEGER I,J,KK,L,M,N,N1,K,KPRIME
      INTEGER PRIME(MPRIME)
      DATA PRIME/
     % 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,
     % 89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,
     % 173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,
     % 263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,
     % 359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,
     % 457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,
     % 569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,
     % 659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,
     % 769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
     % 881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,
     % 997,1009,1013,1019,1021,MREST*0/
      COMPLEX*16 W,W0,WP
      REAL*8 THETA
      REAL*8 PI2
      DATA PI2 /6.28318 53071 79586 477D0/
C
      IF(NN.LE.1) GOTO 900
      IF(ABS(ISGN).NE.1) GOTO 920
      M=NN
      KPRIME=1
      KK=2
 200  N1=M
 220  IF(MOD(N1,KK).EQ.0) THEN
        N1=N1/KK
        GOTO 220
      ELSE
        L=NN/M
        N=M/N1
        M=N1
        IF(N.NE.1) THEN
          CALL ZFFTW0(A,L*M,N,KK,ISGN,WORK)
          IF(M.NE.1) THEN
            W0=1
            THETA=PI2/(ISGN*M*N)
            WP=DCMPLX(COS(THETA),SIN(THETA))
c            WP=DCMPLX(-2D0*SIN(THETA/2D0)**2,SIN(THETA))
            DO 260 J=0,M-1
              W=1
              DO 240 K=0,N-1
                DO 230 I=0,L-1
 230            A(I+L*(J+M*K))=W*A(I+L*(J+M*K))
                W=W0*W
 240          CONTINUE
c              W0=W0*WP+W0
              W0=W0*WP
 260        CONTINUE
            CALL ZFTWTR(A,L,M,N,WORK)
          ENDIF
        ENDIF
        IF(M.NE.1) THEN
          KPRIME=KPRIME+1
          IF(KPRIME.LE.MPRIME) THEN
            IF(KPRIME.GT.MPRIME0) CALL NEXTPRIM(PRIME,KPRIME)
            KK=PRIME(KPRIME)
          ELSE
            KK=KK+2
          ENDIF
          GOTO 200
        ENDIF
      ENDIF
      RETURN
 900  IF(NN.LE.0) WRITE(6,910) NN
 910  FORMAT(' (SUBR.ZFFTW) Invalid 2nd argument n=',I4,
     %  '. Must be >=1.')
      RETURN
 920  WRITE(6,930) ISGN
 930  FORMAT(' (SUBR.ZFFTW) Invalid 3rd argument ISGN=',I4,
     %  '. Must be 1 or -1.')
      RETURN
      END
C-------------------------- NEXTPRIM ------------------------------
      SUBROUTINE NEXTPRIM(PRIME,KPRIME)
      IMPLICIT NONE
      INTEGER KPRIME,PRIME(KPRIME)
      INTEGER K,P
      P=PRIME(KPRIME-1)
 100  P=P+2
      K=0
 120  K=K+1
      IF(PRIME(K).GT.P/PRIME(K)) THEN
        PRIME(KPRIME)=P
        RETURN
      ENDIF
      IF(MOD(P,PRIME(K)).EQ.0) GOTO 100
      GOTO 120
      END
C-------------------------- ZFFTW0 --------------------------------------
      SUBROUTINE ZFFTW0(A,M,N,K,ISGN,WORK)
C Fast Fourier Transformation of arbitrary base to be called by ZFFTW.
C    B(l,k)= Sum A(l,j)*exp[+-2*pi*i*k*j/N]   over j=0 to N-1
C             for k=0 to N-1, l=1 to M.
C Input:
C   A    Complex*16 array of dimension A(M,0:N-1)
C   M    First dimension of A.
C   N    Length of Fourier transformation. 
C         N= K**(positive integer). (Not checked for.)
C   K    base.  >=2. (need not be a prime)
C   ISGN +1 or -1, sign of the phase.
C Work area
C   WORK Complex*16 array of length K.
C Output:
C   A    B(l,k) is stored in the same place as A.
C   
      IMPLICIT NONE
      INTEGER M,N,K,ISGN
      COMPLEX*16 A(M,0:N-1)
      INTEGER M1,N1,I,J,NN,IS,K1,L
      COMPLEX*16 T,W,WP,W0,WP0,W2,WORK(0:K-1)
      REAL*8 THETA,SGN
      REAL*8 PI2
      DATA PI2 /6.28318 53071 79586 477D0/
      REAL*8 C31,C51,C52,C53,C54
      DATA C31/0.86602 54037 84438 64676D0/,
     %     C51/0.30901 69943 74947 42410D0/,
     %     C52/0.95105 65162 95153 57212D0/,
     %     C53/-.80901 69943 74947 42410D0/,
     %     C54/0.58778 52522 92473 12917D0/
      COMPLEX*16 X0,X1,X2,X3,X4,X5,X6,X7,X8
C     begin initialize for preventing compiler warning
      WP0=(0.D0,0.D0)
C     end   initialize for preventing compiler warning
      J=0
      K1=K-1
      NN=N/K*K1
      DO 120 I=0,N-1
        IF(J.GT.I) THEN
          DO 110 M1=1,M
          T=A(M1,J)
          A(M1,J)=A(M1,I)
 110      A(M1,I)=T
        ENDIF
C  Bit inversion in K-scale
        N1=NN
        DO WHILE (N1.GE.K1.AND.J.GE.N1)
          J=J-N1
          N1=N1/K
        ENDDO
        J=J+N1/K1
 120  CONTINUE
C  Following statements apply for any K>=2, but in
C  order for the computing time, the cases K<=5 are
C  written separately.
      NN=1
      SGN=DFLOAT(ISGN)
      IF(K.GE.6) THEN
        THETA=PI2/(ISGN*K)
c        WP0=DCMPLX(-2D0*SIN(THETA/2D0)**2,SIN(THETA))
        WP0=DCMPLX(COS(THETA),SIN(THETA))
      ENDIF
      DO WHILE (NN.LT.N)
        IS=K*NN
        THETA=PI2/(ISGN*IS)
c        WP=DCMPLX(-2D0*SIN(THETA/2D0)**2,SIN(THETA))
        WP=DCMPLX(COS(THETA),SIN(THETA))
        W=1
        DO 300 N1=0,NN-1
          DO 280 M1=1,M
            DO 260 I=N1,N1+N-1,IS
            IF(K.EQ.2) THEN
              T=W*A(M1,I+NN)
              A(M1,I+NN)=A(M1,I)-T
              A(M1,I)=A(M1,I)+T
            ELSEIF(K.EQ.3) THEN
              X0=A(M1,I)
              X1=W*A(M1,I+NN)
              X2=W**2*A(M1,I+2*NN)
              X3=X1+X2
              X4=DCMPLX(0D0,SGN*C31)*(X1-X2)
              X5=X0-0.5D0*X3
              A(M1,I)=X0+X3
              A(M1,I+NN)=X5+X4
              A(M1,I+2*NN)=X5-X4
            ELSEIF(K.EQ.4) THEN
              W2=W**2
              X1=W*A(M1,I+NN)
              X2=W2*A(M1,I+2*NN)
              X3=W2*W*A(M1,I+3*NN)
              X4=A(M1,I)+X2
              X5=A(M1,I)-X2
              X6=X1+X3
              X7=DCMPLX(0D0,SGN)*(X1-X3)
              A(M1,I)=X4+X6
              A(M1,I+NN)=X5+X7
              A(M1,I+2*NN)=X4-X6
              A(M1,I+3*NN)=X5-X7
            ELSEIF(K.EQ.5) THEN
              W2=W**2
              X1=W*A(M1,I+NN)
              X2=W2*A(M1,I+2*NN)
              X3=W2*W*A(M1,I+3*NN)
              X4=W2*W2*A(M1,I+4*NN)
              X0=A(M1,I)
              X5=X1+X4
              X6=DCMPLX(0D0,SGN)*(X1-X4)
              X7=X2+X3
              X8=DCMPLX(0D0,SGN)*(X2-X3)
              A(M1,I)=X0+X5+X7
              X1=X0+C51*X5+C53*X7
              X2=X0+C53*X5+C51*X7
              X3=C52*X6+C54*X8
              X4=C54*X6-C52*X8
              A(M1,I+NN)=X1+X3
              A(M1,I+2*NN)=X2+X4
              A(M1,I+3*NN)=X2-X4
              A(M1,I+4*NN)=X1-X3
            ELSE
              W0=W
              DO 200 J=0,K1
 200          WORK(J)=A(M1,I+J*NN)
              DO 240 L=0,K1
                T=WORK(K1)
                DO 220 J=K1-1,0,-1
 220              T=T*W0+WORK(J)
                A(M1,I+L*NN)=T
c                W0=W0*WP0+W0
                W0=W0*WP0
 240          CONTINUE
            ENDIF
 260        CONTINUE
 280      CONTINUE
c          W=W*WP+W
          W=W*WP
 300    CONTINUE
        NN=IS
      ENDDO
      RETURN
      END
C--------------- ZFTWTR -----------------
      SUBROUTINE ZFTWTR(A,L,M,N,W)
      IMPLICIT NONE
      INTEGER L,M,N,I,J,K
      COMPLEX*16 A(L,0:M*N-1),W(0:M*N-1)
      DO 300 I=1,L
        DO 200 J=0,M-1
        DO 200 K=0,N-1
 200    W(J+M*K)=A(I,J+M*K)
        DO 240 J=0,M-1
        DO 240 K=0,N-1
 240    A(I,K+N*J)=W(J+M*K)
 300  CONTINUE
      RETURN
      END

      subroutine tfmatrixmaybeq(k,cmplm,realm,vec,n,m)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*8 k,ka,kai
      integer*4 i,n,m
      logical*4 cmplm,realm,vec
      n=0
      m=0
      realm=.false.
      cmplm=.false.
      vec=.false.
      if(ktflistq(k))then
        ka=ktfaddr(k)
        if(klist(ka) .eq. ktfoper+mtflist)then
          n=ilist(2,ka-1)
          if(ktfnonreallistq(ka))then
            do i=1,n
              if(ktfnonlistq(klist(ka+i)))then
                return
              endif
              kai=ktfaddr(klist(ka+i))
              if(klist(kai) .ne. ktfoper+mtflist)then
                return
              endif
              if(i .eq. 1)then
                m=ilist(2,kai-1)
              elseif(m .ne. ilist(2,kai-1))then
                m=0
                return
              endif
              if(ktfnonreallistq(kai))then
                cmplm=.true.
                return
              endif
            enddo
            realm=.true.
          else
            vec=.true.
          endif
        endif
      endif
      return
      include 'inc/TFSF.inc'
      end
