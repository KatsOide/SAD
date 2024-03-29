      module dot
      use tfstk
      use maloc

      contains
      recursive function tfdot(k1,k2,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_dlist), pointer :: kl1,kl2,kl2i
      integer*4 i,itfmessage,isp0,m,n
      real*8 xr,xi
      complex*16 cx,cx1,cx2
      kx%k=ktfoper+mtfnull
      if(.not. tflistq(k1,kl1) .or. .not. tflistq(k2,kl2))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      n=kl1%nl
      if(n == 0)then
        kx%k=0
        return
      endif
      isp0=isp
      if(tflistq(kl1%dbody(1)))then
        do i=1,n
          isp=isp+1
          dtastk(isp)=tfdot(kl1%dbody(i),k2,irtc)
          if(irtc /= 0)then
            isp=isp0
            return
          endif
        enddo
        kx=kxmakelist(isp0)
        isp=isp0
        return
      endif
      if(kl2%nl /= n)then
        irtc=itfmessage(9,'General::equalleng','"arguments"')
        return
      endif
      if(ktfreallistq(kl1))then
        if(ktfreallistq(kl2))then
          xr=dot_product(kl1%rbody(1:n),kl2%rbody(1:n))
          kx=sad_descr(xr)
          return
        else
          if(tflistq(kl2%dbody(1),kl2i))then
            m=kl2i%nl
            isp=isp+m
            call tfloadrcstk(isp0,kl2i,kl1%rbody(1),irtc)
            if(irtc /= 0)then
              go to 3000
            endif
            do i=2,n
              if(tflistq(kl2%dbody(i),kl2i))then
                if(kl2i%nl /= m)then
                  irtc=itfmessage(9,'General::equalleng',
     $                 '"elements of matrix"')
                  isp=isp0
                  return
                endif
                call tfaddrcstk(isp0,kl2i,kl1%rbody(i),irtc)
                if(irtc /= 0)then
                  go to 3000
                endif
              else
                go to 3000
              endif
            enddo
          else
            xr=0.d0
            xi=0.d0
            do i=1,n
              if(ktfrealq(kl2%dbody(i)))then
                xr=xr+kl1%rbody(i)*kl2%rbody(i)
              elseif(tfcomplexq(kl2%dbody(i),cx))then
                xr=xr+kl1%rbody(i)*dble(cx)
                xi=xi+kl1%rbody(i)*imag(cx)
              else
                go to 3000
              endif
            enddo
            kx=kxcalocv(-1,xr,xi)
            return
          endif
        endif
      else
        if(ktfreallistq(kl2))then
          xr=0.d0
          xi=0.d0
          do i=1,n
            if(ktfrealq(kl1%dbody(i)))then
              xr=xr+kl1%rbody(i)*kl2%rbody(i)
            elseif(tfcomplexq(kl1%dbody(i),cx))then
              xr=xr+dble(cx)*kl2%rbody(i)
              xi=xi+imag(cx)*kl2%rbody(i)
            else
              go to 3000
            endif
          enddo
          kx=kxcalocv(-1,xr,xi)
          return
        else
          if(tflistq(kl2%dbody(1),kl2i))then
            m=kl2i%nl
            isp=isp+m
            rtastk (isp0+1:isp0+m)=0.d0
            rtastk2(isp0+1:isp0+m)=0.d0
            do i=1,n
              call descr_sad(kl2%dbody(i),kl2i)
              if(kl2i%nl /= m)then
                irtc=itfmessage(9,'General::equalleng',
     $               '"elements of matrix"')
                isp=isp0
                return
              endif
              if(ktfrealq(kl1%dbody(i)))then
                call tfaddrcstk(isp0,kl2i,kl1%rbody(i),irtc)
                if(irtc /= 0)then
                  go to 3000
                endif
              elseif(tfcomplexq(kl1%dbody(i),cx))then
                call tfaddccstk(isp0,kl2i,cx,irtc)
                if(irtc /= 0)then
                  go to 3000
                endif
              else
                go to 3000
              endif
            enddo
          else
            cx=(0.d0,0.d0)
            do i=1,n
              if(ktfrealq(kl1%dbody(i)))then
                if(ktfrealq(kl2%dbody(i)))then
                  cx=cx+kl1%rbody(i)*kl2%rbody(i)
                elseif(tfcomplexq(kl2%dbody(i),cx1))then
                  cx=cx+kl1%rbody(i)*cx1
                else
                  go to 3000
                endif
              elseif(tfcomplexq(kl1%dbody(i),cx1))then
                if(ktfrealq(kl2%dbody(i)))then
                  cx=cx+cx1*kl2%rbody(i)
                elseif(tfcomplexq(kl2%dbody(i),cx2))then
                  cx=cx+cx1*cx2
                else
                  go to 3000
                endif
              endif
            enddo
            kx=kxcalocc(-1,cx)
            return
          endif
        endif
      endif
      do i=isp0,isp
        if(rtastk2(i) /= 0.d0)then
          dtastk(i)=kxcalocv(-1,rtastk(i),rtastk2(i))
        endif
      enddo
      kx=kxmakelist(isp0)
      isp=isp0
      irtc=0
      return
 3000 isp=isp0
      call tfinner(k1,k2,kx,sad_descr(ktfoper+mtfplus),
     $     sad_descr(ktfoper+mtftimes),irtc)
      return
      end

      subroutine tfloadrcstk(isp0,kl,x,irtc)
      implicit none
      type (sad_dlist) ,intent(in):: kl
      integer*4 ,intent(in):: isp0
      integer*4 ,intent(out):: irtc
      integer*4 i
      real*8 ,intent(in):: x
      complex*16 c
      if(ktfreallistq(kl))then
        rtastk (isp0+1:isp0+kl%nl)=x*kl%rbody(1:kl%nl)
        rtastk2(isp0+1:isp0+kl%nl)=0.d0
      else
        do i=1,kl%nl
          if(ktfrealq(kl%dbody(i)))then
            rtastk (isp0+i)=x*kl%rbody(i)
            rtastk2(isp0+i)=0.d0
          elseif(tfcomplexq(kl%dbody(i),c))then
            rtastk (isp0+i)=x*dble(c)
            rtastk2(isp0+i)=x*imag(c)
          else
            irtc=-1
            return
          endif
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfaddrcstk(isp0,kl,x,irtc)
      implicit none
      type (sad_dlist) ,intent(in):: kl
      integer*4 ,intent(in):: isp0
      integer*4 ,intent(out):: irtc
      integer*4 i
      real*8 ,intent(in):: x
      complex*16 c
      if(ktfreallistq(kl))then
        rtastk(isp0+1:isp0+kl%nl)=
     $       rtastk(isp0+1:isp0+kl%nl)+x*kl%rbody(1:kl%nl)
      else
        do i=1,kl%nl
          if(ktfrealq(kl%dbody(i)))then
            rtastk(isp0+i)=rtastk(isp0+i)+x*kl%rbody(i)
          elseif(tfcomplexq(kl%dbody(i),c))then
            rtastk (isp0+i)=rtastk(isp0+i) +x*dble(c)
            rtastk2(isp0+i)=rtastk2(isp0+i)+x*imag(c)
          else
            irtc=-1
            return
          endif
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfaddccstk(isp0,kl,cx,irtc)
      implicit none
      type (sad_dlist) ,intent(in):: kl
      integer*4 ,intent(in):: isp0
      integer*4 ,intent(out):: irtc
      integer*4 i
      complex*16 ,intent(in):: cx
      complex*16 c,cx1
      if(ktfreallistq(kl))then
        rtastk(isp0+1:isp0+kl%nl)=
     $       rtastk (isp0+1:isp0+kl%nl)+dble(cx)*kl%rbody(1:kl%nl)
        rtastk2(isp0+1:isp0+kl%nl)=
     $       rtastk2(isp0+1:isp0+kl%nl)+imag(cx)*kl%rbody(1:kl%nl)
      else
        do i=1,kl%nl
          if(ktfrealq(kl%dbody(i)))then
            rtastk (isp0+i)=rtastk (isp0+i)+dble(cx)*kl%rbody(i)
            rtastk2(isp0+i)=rtastk2(isp0+i)+imag(cx)*kl%rbody(i)
          elseif(tfcomplexq(kl%dbody(i),c))then
            cx1=cx*c
            rtastk (isp0+i)=rtastk (isp0+i)+dble(cx1)
            rtastk2(isp0+i)=rtastk2(isp0+i)+imag(cx1)
          else
            irtc=-1
            return
          endif
        enddo
      endif
      irtc=0
      return
      end

      subroutine tfouter(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 kh,kai
      integer*4 narg,i,mi,ispi,isp0,j,mn,itfmessage
      narg=isp-isp1
      if(narg .lt. 3)then
        irtc=itfmessage(9,'General::narg','"3 or more"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+2)))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      itastk2(1,isp1+2)=1
      do i=isp1+3,isp
        if(ktfnonlistq(ktastk(i)) .or.
     $       .not. tfsameheadq(ktastk(isp1+2),ktastk(i)))then
          irtc=itfmessage(9,'General::samehead',' ')
          return
        endif
        itastk2(1,i)=1
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
     $         itastk2(1,isp1+j))
        enddo
        isp=isp+1
        ktastk(isp)=klist(ktfaddr(ktastk(isp1+narg))+i)
        call tfefunrefstk(isp0,isp0,irtc)
        if(irtc /= 0)then
          call tflocal(kh)
          return
        endif
      enddo
      itastk2(1,isp1+narg)=mn
      do i=narg,2,-1
        mi=ilist(2,ktfaddr(ktastk(isp1+i))-1)
        if(itastk2(1,isp1+i) .ge. mi)then
          ispi=isp-mi
          kai=ktfmakelist(ispi)
          isp=ispi+1
          ktastk(isp)=ktflist+kai
          klist(kai)=ktfcopy(kh)
          itastk2(1,isp1+i)=1
        else
          itastk2(1,isp1+i)=itastk2(1,isp1+i)+1
          go to 1
        endif
      enddo
      kx=dtastk(isp)
      isp=isp1+narg
      call tflocal(kh)
      irtc=0
      return
      end

      subroutine tftranspose(k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kij
      type (sad_dlist), pointer :: kl,kl1,kli,klj,klx,klxi
      integer*4 m,n,i,j,itfmessage
      logical*4 d
      if(.not. tflistq(k,kl))then
        go to 9000
      endif
      if(tfnonlistq(kl%dbody(1),kl1))then
        go to 9000
      endif
      m=kl%nl
      n=kl1%nl
      do i=2,m
        if(tfnonlistq(kl%dbody(i),kli))then
          go to 9000
        endif
        if(kli%nl /= n)then
          irtc=itfmessage(9,'General::equalleng','"rows of matrix"')
          return
        endif
      enddo
      kx=kxadaloc(-1,n,klx)
      do i=1,n
        klx%dbody(i)=kxaaloc(0,m,klxi)
        d=.false.
        do j=1,m
          call loc_sad(ktfaddr(kl%dbody(j)),klj)
          kij=klj%dbody(i)
          if(ktfrealq(kij))then
            klxi%dbody(j)=kij
          else
            d=.true.
            klxi%dbody(j)=dtfcopy(kij)
          endif
        enddo
        if(d)then
          klxi%attr=ior(klxi%attr,lnonreallist)
        else
          klxi%attr=lconstlist
        endif
      enddo
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Matrix"')
      return
      end

      subroutine tfdet(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kl
      integer*4 n,m,narg,itfmessage
      logical*4 cmplm,realm,vec
      narg=isp-isp1
      if(narg /= 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      call tfmatrixmaybeq(dtastk(isp),cmplm,realm,vec,n,m,kl)
      if(n /= m .or. m == 0)then
        go to 9000
        return
      endif
      if(cmplm)then
        call tfcdet(kl,kx,n,irtc)
        if(irtc /= 0)then
          go to 9000
        endif
      elseif(realm)then
        call tfrdet(kl,kx,n)
      else
        go to 9000
      endif
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Numerical Square Matrix"')
      return
      end

      subroutine tfrdet(kl,kx,n)
      implicit none
      type (sad_descriptor) , intent(out)::kx
      type (sad_dlist) , intent(in):: kl
      integer*4 ,intent(in)::n
      real*8, allocatable :: a(:,:)
      allocate (a(n,n))
      call tfl2m(kl,a,n,n,.false.)
      kx=dfromr(tdet(a,n,n))
      deallocate (a)
      return
      end

      subroutine tfcdet(kl,kx,n,irtc)
      implicit none
      type (sad_descriptor) ,intent(out)::kx
      type (sad_dlist) ,intent(in)::kl
      integer*4 ,intent(in):: n
      integer*4 ,intent(out):: irtc
      complex*16 , allocatable ::c(:,:)
      complex*16 cx
      allocate (c(n,n))
      call tfl2cm(kl,c,n,n,.false.,irtc)
      if(irtc /= 0)then
        deallocate(c)
        return
      endif
      cx=tcdet(c,n,n)
      kx=kxcalocc(-1,cx)
      deallocate(c)
      return
      end

      subroutine tfsingularvalues(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k
      type (sad_dlist), pointer :: kl,klx
      integer*8 kux,kvx,kwx
      integer*4 n,m,narg,itfmessage
      real*8 eps
      logical*4 inv
      logical*4 realm,cmplm,vec
      narg=isp-isp1
      if(narg /= 3)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(ktfnonrealq(dtastk(isp)) .or. ktfnonrealq(dtastk(isp-1)))then
        irtc=itfmessage(9,'General::wrongtype','"Real number for #2 and #3"')
        return
      endif
      eps=rtastk(isp-1)
      inv=rtastk(isp) /= 0.d0
      k=dtastk(isp1+1)
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m,kl)
      if(cmplm)then
        if(max(n,m) > 200)then
          call tftclupdate(3)
        endif
        call tcsvdma(kl,kux,kvx,kwx,n,m,eps,inv,irtc)
        if(irtc /= 0)then
          return
        endif
      elseif(realm)then
        if(max(n,m) > 200)then
          call tftclupdate(3)
        endif
        call tsvdma(kl,kux,kvx,kwx,n,m,eps,inv)
        irtc=0
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Matrix"')
        return
      endif
      kx=kxadaloc(-1,3,klx)
      klx%dbody(1)%k=ktflist+ktfcopy1(kux)
      klx%dbody(2)%k=ktflist+ktfcopy1(kwx)
      klx%dbody(3)%k=ktflist+ktfcopy1(kvx)
      return
      end

      subroutine tcsvdma(kl,kux,kvx,kwx,n,m,eps,inv,irtc)
      implicit none
      type (sad_dlist) , intent(in)::kl
      integer*8 , intent(out)::kux,kvx,kwx
      integer*4 ,intent(in)::n,m
      integer*4 mn
      integer*4 ,intent(out)::irtc
      real*8 ,intent(in)::eps
      complex*16, allocatable:: ca(:,:),cu(:,:)
      real*8 , allocatable:: w(:)
      logical*4 ,intent(in)::inv
      allocate (ca(n,m),cu(n,n),w(m))
      call tfl2cm(kl,ca,n,m,.false.,irtc)
      if(irtc /= 0)then
        deallocate(ca,cu,w)
        return
      endif
      call tcsvdm(ca,cu,w,n,m,n,n,eps,inv)
      mn=min(m,n)
      kux=ktfcm2l(cu,mn,n,n,.false.,.false.)
      kvx=ktfcm2l(ca,mn,m,n,.false.,.false.)
      kwx=ktavaloc(-1,mn)
      rlist(kwx+1:kwx+mn)=w(1:mn)
      deallocate(ca,cu,w)
      return
      end

      subroutine tsvdma(kl,kux,kvx,kwx,n,m,eps,inv)
      implicit none
      type (sad_dlist) ,intent(in)::kl
      integer*8,intent(out):: kux,kvx,kwx
      integer*4 ,intent(in)::n,m
      integer*4 mn,ndim
      real*8 , allocatable:: a(:,:),u(:,:),w(:)
      real*8 , intent(in)::eps
      logical*4 ,intent(in)::inv
      ndim=max(n,m-1)
      allocate (a(ndim,m),w(m),u(n,n))
      call tfl2m1(kl,a,n,m,ndim,.false.)
      call tsvdm(a,u,w,n,m,ndim,n,eps,inv)
      mn=min(m,n)
      kux=ktfaddr(kxm2l(u,mn,n,n,.false.))
      kvx=ktfaddr(kxm2l(a,mn,m,ndim,.false.))
      kwx=ktavaloc(-1,mn)
      rlist(kwx+1:kwx+mn)=w(1:mn)
      deallocate(a,w,u)
      return
      end

      subroutine tflinearsolve(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) , intent(out) :: kx
      type (sad_descriptor) k,kb
      type (sad_dlist), pointer :: kl,klb
      integer*4 , intent(in) :: isp1
      integer*4 , intent(out) :: irtc
      integer*4 narg,nb,mb,n,m,itfmessage,mx
      real*8 eps
      logical*4 cmplm,realm,vec
      narg=isp-isp1
      if(narg /= 2 .and. narg /= 3)then
        irtc=itfmessage(9,'General::narg','"2 (+ option)"')
        return
      elseif(narg == 2)then
        k=dtastk(isp-1)
        kb=dtastk(isp)
        eps=1.d-8
      else
        k=dtastk(isp-2)
        kb=dtastk(isp-1)
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype','"real number"')
          return
        endif
        eps=rtastk(isp)
      endif
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m,kl)
      if(.not. realm)then
        go to 9000
      endif
      call tfmatrixmaybeq(kb,cmplm,realm,vec,mb,nb,klb)
      if(.not. realm .and. .not. vec)then
        go to 9100
      endif
      if(nb == 0 .or. vec)then
        nb=mb
        mb=0
      endif
      if(n /= nb)then
        irtc=itfmessage(9,'General::equalleng','"matrix and vector"')
        return
      endif
      if(max(n,m) > 200)then
        call tftclupdate(3)
      endif
      mx=max(mb,1)
      call tflsolve(kl,klb,kx,n,m,mb,mx,eps)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Real Matrix for #1"')
      return
 9100 irtc=itfmessage(9,'General::wrongtype',
     $     '"Real Vector or Matrix for #2"')
      return
      end

      subroutine tflsolve(kl,klb,kx,n,m,mb,mx,eps)
      implicit none
      type (sad_descriptor) , intent(out)::kx
      type (sad_dlist) , intent(in)::kl,klb
      integer*4 , intent(in):: n,m,mb,mx
      real*8 ,intent(in)::eps
      real*8 , allocatable:: a(:,:),b(:,:),x(:,:)
      allocate (a(n,m),b(n,mx),x(m,mx))
      call tfl2m(kl,a,n,m,.false.)
      if(mb == 0)then
        b(:,1)=klb%rbody(1:n)
      else
        call tfl2m(klb,b,n,mb,.true.)
      endif
      call tsolvm(a,b,x,n,m,mx,n,n,m,eps,.true.)
      if(mb == 0)then
        kx=kxm2l(x,0,m,1,.false.)
      else
        kx=kxm2l(x,m,mb,m,.true.)
      endif
      deallocate (a,b,x)
      return
      end

      subroutine tfdiagonalmatrix(k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ki
      type (sad_dlist), pointer :: kl,klx
      type (sad_rlist), pointer :: klri
      integer*4 m,i,itfmessage
      real*8 x
      if(.not. tflistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      m=kl%nl
      kx=kxadaloc(-1,m,klx)
      do i=1,m
        klx%dbody(i)%k=ktflist+ktraaloc(0,m,klri)
        ki=kl%dbody(i)
        if(ktfrealq(ki,x))then
          klri%rbody(i)=x
        else
          klri%attr=ior(lnonreallist,klri%attr)
          klri%dbody(i)=dtfcopy(ki)
        endif
      enddo
      irtc=0
      return
      end

      subroutine tfidentitymatrix(k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klri
      integer*4 ,intent(out):: irtc
      integer*4 m,i,itfmessage
      if(ktfnonrealq(k,m))then
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
        return
      endif
      if(m .le. 0)then
        irtc=itfmessage(9,'General::wrongnum','"Positive"')
        return
      endif
      kx=kxadaloc(-1,m,klx)
      do i=1,m
        klx%dbody(i)%k=ktflist+ktraaloc(0,m,klri)
        klri%rbody(i)=1.d0
      enddo
      irtc=0
      return
      end

      subroutine tfeigensystem(k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kl,klx
      integer*8 kex,kvx
      integer*4 n,m,itfmessage
      logical*4 realm,cmplm,vec
      call tfmatrixmaybeq(k,cmplm,realm,vec,n,m,kl)
      if(n /= m)then
        go to 9000
      endif
      if(cmplm)then
        call tfceigen(kl,n,kvx,kex,irtc)
        if(irtc /= 0)then
          go to 9000
        endif
      elseif(realm)then
        call tfreigen(kl,n,kvx,kex)
      else
        go to 9000
      endif
      kx=kxadaloc(-1,2,klx)
      klx%dbody(1)%k=ktflist+ktfcopy1(kex)
      klx%dbody(2)%k=ktflist+ktfcopy1(kvx)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Square Matrix"')
      return
      end

      subroutine tfceigen(kl,m,kvx,kex,irtc)
      implicit none
      type (sad_dlist) ,intent(in):: kl
      integer*8 ,intent(out):: kvx,kex
      complex*16 , allocatable ::c(:,:),cw(:,:),ce(:)
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: m
      allocate (c(m,m),cw(m,m),ce(m))
      call tfl2cm(kl,c,m,m,.false.,irtc)
      if(irtc /= 0)then
        deallocate (c,cw,ce)
        return
      endif
      if(m > 200)then
        call tftclupdate(3)
      endif
      call tceigen(c,cw,ce,m,m)
      kvx=ktfcm2l(c,m,m,m,.true.,.false.)
      kex=ktfc2l(ce,m)
      deallocate (c,cw,ce)
      return
      end

      subroutine tfreigen(kl,m,kvx,kex)
      implicit none
      type (sad_dlist) ,intent(in):: kl
      integer*8 ,intent(out):: kvx,kex
      integer*8 kvxkvr,kvxkvi,kzi1,kzi2
      integer*4 ,intent(in):: m
      integer*4 i,j
      real*8 , allocatable :: a(:,:),w(:,:)
      complex*16 ce(m)
      logical*4 d,con
      allocate (a(m,m),w(m,m))
      call tfl2m(kl,a,m,m,.false.)
      if(m > 200)then
        call tftclupdate(3)
      endif
      call teigen(a,w,ce,m,m)
      kvx=ktfaddr(kxm2l(a,m,m,m,.true.))
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
                if(rlist(kvxkvi+j) /= 0.d0)then
                  d=.true.
                  dlist(kzi1+j)=
     $                 kxcalocv(0,rlist(kvxkvr+j),rlist(kvxkvi+j))
                  dlist(kzi2+j)=
     $                 kxcalocv(0,rlist(kvxkvr+j),-rlist(kvxkvi+j))
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
      deallocate(a,w)
      return
      end

      subroutine tffourier(inv,k,kx,irtc)
      use macmath
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: kl
      integer*4 ,intent(out):: irtc
      integer*4 m,itfmessage
      logical*4 ,intent(in):: inv
      if(.not. tflistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      m=kl%nl
      irtc=-1
      if(m == 0)then
        return
      elseif(m == 1)then
        kx=k
        irtc=0
        return
      endif
      call tffft(kl,m,kx,inv,irtc)
      if(irtc /= 0)then
        go to 9000
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"Numerical Square Matrix"')
      return
      end

      subroutine tffft(kl,m,kx,inv,irtc)
      use macmath
      use iso_c_binding
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist),intent(in):: kl
      type (sad_complex), pointer :: klic
      integer*4 ,intent(in):: m
      integer*4 ,intent(out):: irtc
      integer*4 n,i
      complex*16 ,allocatable,dimension(:)::cx,cy,ay
      complex*16 , pointer ::ca(:)
      real*8 , allocatable,dimension(:),target::a
      real*8 c,s,f,w
      logical*4 ,intent(in):: inv
      logical*4 power2,even
      allocate(cx(m),cy(m),ay(m),a(m*2))
      irtc=-1
      n=2
      do while(n .lt. m)
        n=n*2
      enddo
      power2=n == m
      even=iand(m,1) == 0
      if(ktfnonreallistqo(kl))then
        f=1.d0/sqrt(dble(m))
        do i=1,m
          if(ktfrealq(kl%dbody(i)))then
            cx(i)=f*kl%rbody(i)
          elseif(tfcomplexq(kl%dbody(i),klic))then
            cx(i)=f*klic%cx(1)
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
        kx%k=ktflist+ktfc2l(cx,m)
        irtc=0
        return
      endif
      call c_f_pointer(c_loc(a),ca,[m/2])
      if(even)then
        a(1:m)=kl%rbody(1:m)
        if(power2)then
          if(m .ge. 4)then
            call trftr(a,m,inv)
            f=1.d0/sqrt(dble(m))
            a(1:m)=f*a(1:m)
            a(m*2-1:m+3:-2)= a(3:m-1:2)
            a(m*2:  m+4:-2)=-a(4:m:  2)
c            do i=1,m/2-1
c              a((m-i)*2+1)= a(i*2+1)
c              a((m-i)*2+2)=-a(i*2+2)
c            enddo
            a(m+1)=a(2)
            a(m+2)=0.d0
            a(2  )=0.d0
            go to 1000
          else
            call tcftr(ca,m/2,inv)
          endif
        else
          if(inv)then
            call zfftw(ca,m/2,-1,ay)
          else
            call zfftw(ca,m/2, 1,ay)
          endif
        endif
        f=.5d0/sqrt(dble(m))
        w=merge(-2.d0,2.d0,inv)*pi/m
        do concurrent (i=1:m/2-1)
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
        a(3:m-1:2)= a(m*2-1:m+3:-2)
        a(4:m:  2)=-a(m*2:m+4:-2)        
c        do i=1,m/2-1
c          a(i*2+1)= a((m-i)*2+1)
c          a(i*2+2)=-a((m-i)*2+2)
c        enddo
      else
        f=1.d0/sqrt(dble(m))
        a(1:2*m-1:2)=f*kl%rbody(1:m)
        a(2:2*m  :2)=0.d0
c        do i=1,m
c          a(i*2-1)=f*kl%rbody(i)
c          a(i*2  )=0.d0
c        enddo
        if(inv)then
          call zfftw(ca,m,-1,ay)
        else
          call zfftw(ca,m, 1,ay)
        endif
      endif
 1000 kx%k=ktflist+ktfc2l(ca,m)
      irtc=0
      return
      end

      recursive subroutine tfinner(k1,k2,kx,ks,kp,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ,intent(in):: k1,k2,ks,kp
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k1i,ki,tfefunrefu
      type (sad_dlist), pointer :: kl1,kl2,klx
      integer*4 m1,i,m2,isp0,itfmessage
      if(ktfnonlistq(k1,kl1) .or. ktfnonlistq(k2,kl2))then
        irtc=itfmessage(9,'General::wrongtype','"List"')
        return
      endif
      m1=kl1%nl
      if(ktfnonreallistqo(kl1) .and. tflistq(kl1%head))then
        kx=kxaaloc(-1,m1,klx)
        do i=1,m1
          k1i=kl1%dbody(i)
          call tfinner(k1i,k2,ki,ks,kp,irtc)
          if(irtc /= 0)then
            klx%dbody(1:m1)%k=merge(ktfoper+mtfnull,i00,
     $           ktfnonreallistqo(klx))
            return
          endif
          if(ktfnonrealq(ki))then
            klx%dbody(i)=dtfcopy1(ki)
            klx%attr=ior(klx%attr,lnonreallist)
          else
            klx%dbody(i)=ki
          endif
        enddo
        return
      endif
      m2=kl2%nl
      if(m2 /= m1)then
        irtc=itfmessage(9,'General::equalleng','"arguments"')
        return
      endif
      isp=isp+3
      isp0=isp
      do i=1,m1
        dtastk(isp-2)=kp
        dtastk(isp-1)=kl1%dbody(i)
        dtastk(isp  )=kl2%dbody(i)
        ki=tfefunrefu(isp-2,irtc)
        isp=isp0
        if(irtc /= 0)then
          isp=isp-3
          return
        endif
        if(i == 1)then
          kx=ki
        else
          dtastk(isp-2)=ks
          dtastk(isp-1)=kx
          dtastk(isp  )=ki
          kx=tfefunrefu(isp-2,irtc)
          isp=isp0
          if(irtc /= 0)then
            isp=isp-3
            return
          endif
        endif
      enddo
      isp=isp-3
      return
      end

      subroutine tftr(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) :: tfefunrefu
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kl,kli
      integer*4 n,m,narg,itfmessage,i,isp0,nm
      real*8 s
      logical*4 cmplm,realm,vec
      narg=isp-isp1
      if(narg /= 1 .and. narg /= 2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      call tfmatrixmaybeq(dtastk(isp1+1),cmplm,realm,vec,n,m,kl)
      if(m == 0)then
        go to 9000
        return
      endif
      nm=min(n,m)
      if(narg == 1)then
        if(realm)then
          s=0.d0
          do i=1,nm
            if(ktflistq(kl%dbody(i),kli))then
              s=s+kli%rbody(i)
            else
              go to 9000
            endif
          enddo
          kx=dfromr(s)
          irtc=0
        else
          isp0=isp
          isp=isp+1
          ktastk(isp)=ktfoper+mtfplus
          do i=1,nm
            if(ktflistq(kl%dbody(i),kli))then
              isp=isp+1
              dtastk(isp)=kli%dbody(i)
            else
              isp=isp0
              go to 9000
            endif
          enddo
          kx=tfefunrefu(isp0+1,irtc)
          isp=isp0
        endif
      else
        isp0=isp
        isp=isp+1
        dtastk(isp)=dtastk(isp1+2)
        do i=1,nm
          if(ktflistq(kl%dbody(i),kli))then
            isp=isp+1
            dtastk(isp)=kli%dbody(i)
          else
            isp=isp0
            go to 9000
          endif
        enddo
        kx=tfefunrefu(isp0+1,irtc)
        isp=isp0
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"Matrix"')
      return
      end

      real*8 function tdet(a,n,ndim)
      use mathfun, only:p2h
      implicit none
      integer*4 ,intent(in)::n,ndim
      integer*4 i,j,k
      real*8 ,intent(inout)::a(ndim,n)
      real*8 d,di,p,x,u
      d=1.d0
      do i=1,n-1
        di=a(i,i)
        do j=i+1,n
          if(a(j,i) .ne. 0.d0)then
            if(di .eq. 0.d0)then
              do k=i+1,n
                x=a(j,k)
                a(j,k)=a(i,k)
                a(i,k)=x
              enddo
              di=a(j,i)
            else
              p=a(j,i)/di
c              u=1.d0/sqrt(1.d0+p**2)
              u=1.d0/p2h(p)
              do k=i+1,n
                x=a(j,k)
                a(j,k)=(x-p*a(i,k))*u
                a(i,k)=(a(i,k)+p*x)*u
              enddo
              di=(di+p*a(j,i))*u
            endif
          endif
        enddo
        d=d*di
        if(d .eq. 0.d0)then
          tdet=0.d0
          return
        endif
      enddo
      tdet=d*a(n,n)
      return
      end function

      complex*16 function tcdet(a,n,ndim)
      implicit none
      integer*4 ,intent(in)::n,ndim
      integer*4 i,j,k
      real*8 u
      complex*16 ,intent(inout)::a(ndim,n)
      complex*16 d,di,p,x,pc
      d=(1.d0,0.d0)
      do i=1,n-1
        di=a(i,i)
        do j=i+1,n
          if(a(j,i) .ne. (0.d0,0.d0))then
            if(di .eq. (0.d0,0.d0))then
              do k=i+1,n
                x=a(j,k)
                a(j,k)=a(i,k)
                a(i,k)=x
              enddo
              di=a(j,i)
            else
              p=a(j,i)/di
              pc=conjg(p)
              u=1.d0/sqrt(1.d0+dble(p*pc))
              do k=i+1,n
                x=a(j,k)
                a(j,k)=(x-p*a(i,k))*u
                a(i,k)=(a(i,k)+pc*x)*u
              enddo
              di=(di+pc*a(j,i))*u
            endif
          endif
        enddo
        d=d*di
        if(d .eq. 0.d0)then
          tcdet=0.d0
          return
        endif
      enddo
      tcdet=d*a(n,n)
      return
      end function

      end module
c
c
c  K. Yokoya's FFT routines.   Received 2/26/1998.
c  Revised Jul 6, 2020 for FORTRAN 2018 compliance.
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
      use macmath ,MPRIME=>nsprime,PRIME=>smallprime,MPRIME0=>nsprime
      IMPLICIT NONE
      integer*4 ,PARAMETER ::MREST=0
      INTEGER*4 NN,ISGN
      COMPLEX*16 A(0:NN-1),WORK(0:NN-1)
      INTEGER*4  I,J,KK,L,M,N,N1,K,KPRIME
      COMPLEX*16 W,W0,WP
      REAL*8 THETA
C
      IF(NN.LE.1) GOTO 900
      IF(ABS(ISGN)/=1) GOTO 920
      M=NN
      KPRIME=1
      KK=2
 200  N1=M
 220  IF(MOD(N1,KK)==0) THEN
        N1=N1/KK
        GOTO 220
      ELSE
        L=NN/M
        N=M/N1
        M=N1
        IF(N/=1) THEN
          CALL ZFFTW0(A,L*M,N,KK,ISGN,WORK)
          IF(M/=1) THEN
            W0=1
            THETA=PI2/(ISGN*M*N)
            WP=DCMPLX(COS(THETA),SIN(THETA))
c            WP=DCMPLX(-2D0*SIN(THETA/2D0)**2,SIN(THETA))
            DO J=0,M-1
              W=1
              DO K=0,N-1
                DO  I=0,L-1
                  A(I+L*(J+M*K))=W*A(I+L*(J+M*K))
                enddo
                W=W0*W
              enddo
c              W0=W0*WP+W0
              W0=W0*WP
            enddo
            CALL ZFTWTR(A,L,M,N,WORK)
          ENDIF
        ENDIF
        IF(M/=1) THEN
          KPRIME=KPRIME+1
          IF(KPRIME.LE.MPRIME) THEN
c            IF(KPRIME>MPRIME0) CALL NEXTPRIM(PRIME,KPRIME)
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
      IF(PRIME(K)>P/PRIME(K)) THEN
        PRIME(KPRIME)=P
        RETURN
      ENDIF
      IF(MOD(P,PRIME(K))==0) GOTO 100
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
      use macmath
      IMPLICIT NONE
      INTEGER M,N,K,ISGN
      COMPLEX*16 A(M,0:N-1)
      INTEGER M1,N1,I,J,NN,IS,K1,L
      COMPLEX*16 T(1:N),TM(1:M),W,WP,W0,WP0,W2,WORK(0:K-1),U
      REAL*8 THETA,SGN
      REAL*8 , parameter ::
     $     C31=0.86602 54037 84438 64676D0,
     %     C51=0.30901 69943 74947 42410D0,
     %     C52=0.95105 65162 95153 57212D0,
     %     C53=-.80901 69943 74947 42410D0,
     %     C54=0.58778 52522 92473 12917D0
      COMPLEX*16 X0,X1,X2,X3,X4,X5,X6,X7,X8
C     begin initialize for preventing compiler warning
      WP0=(0.D0,0.D0)
C     end   initialize for preventing compiler warning
      J=0
      K1=K-1
      NN=N/K*K1
      DO I=0,N-1
        IF(J.GT.I) THEN
          TM=A(:,J)
          A(:,J)=A(:,I)
          A(:,I)=TM
        ENDIF
C  Bit inversion in K-scale
        N1=NN
        DO WHILE (N1.GE.K1.AND.J.GE.N1)
          J=J-N1
          N1=N1/K
        ENDDO
        J=J+N1/K1
      ENDDO
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
        DO N1=0,NN-1
          DO  M1=1,M
            SELECT CASE (K)
            CASE (2)
              T(1:N:IS)=W*A(M1,N1+NN:N1+N-1+NN:IS)
              A(M1,N1+NN:N1+N-1+NN:IS)=A(M1,N1:N1+N-1:IS)-T(1:N:IS)
              A(M1,N1:N1+N-1:IS)=A(M1,N1:N1+N-1:IS)+T(1:N:IS)
            CASE (3)
              DO concurrent (I=N1:N1+N-1:IS)
                X0=A(M1,I)
                X1=W*A(M1,I+NN)
                X2=W**2*A(M1,I+2*NN)
                X3=X1+X2
                X4=DCMPLX(0D0,SGN*C31)*(X1-X2)
                X5=X0-0.5D0*X3
                A(M1,I)=X0+X3
                A(M1,I+NN)=X5+X4
                A(M1,I+2*NN)=X5-X4
              enddo
            CASE (4)
              DO concurrent (I=N1:N1+N-1:IS)
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
              enddo
            CASE (5)
              DO concurrent (I=N1:N1+N-1:IS)
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
              enddo
            CASE DEFAULT
              DO concurrent (I=N1:N1+N-1:IS)
                W0=W
                WORK(0:K1)=A(M1,I:I+K1*NN:NN)
                DO L=0,K1
                  U=WORK(K1)
                  DO J=K1-1,0,-1
                    U=U*W0+WORK(J)
                  ENDDO
                  A(M1,I+L*NN)=U
                  W0=W0*WP0
                ENDDO
              ENDDO
            END SELECT
          ENDDO
c     W=W*WP+W
          W=W*WP
        ENDDO
        NN=IS
      ENDDO
      RETURN
      END
C---------------ZFTWTR -----------------
      SUBROUTINE ZFTWTR(A,L,M,N,W)
      IMPLICIT NONE
      INTEGER L,M,N,I,J
      COMPLEX*16 A(L,0:M*N-1),W(0:M*N-1)
      DO concurrent (I=1:L)
        DO J=0,M-1
          W(J:J+M*(N-1):M)=A(I,J:J+M*(N-1):M)
        ENDDO
        DO J=0,M-1
          A(I,N*J:N-1+N*J)=W(J:J+M*(N-1):M)
        ENDDO
      ENDDO
      RETURN
      END

