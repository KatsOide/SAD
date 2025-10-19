      subroutine pmicad(a,b,x,nx0,m,n,mdim,norm,eps,svd,solvf)
c --- solve A.x=b  using "Best Corrector Method" ----------
      use solv
      implicit real*8 (a-h,o-z)
c     parameter (nmax=1000,prec=1d-2)
      parameter (nmax=1000,prec=1d-6)
      dimension a(mdim,n),b(m),x(n),iw(nmax)
      logical norm,svd,solvf
c     begin begin initialize for preventing compiler warning
      tm=0.d0
      small=0.d0
      k=0
      sk=0.d0
c     end   begin initialize for preventing compiler warning

c     call cputime(ctime0,irtc)
      nx=min(nx0,m,n)
      if(n.gt.nmax) then
        print *,' Too many variables:',n
        return
      endif
      do 5 i=1,n
        iw(i)=i
 5    continue
      if(norm) then
        tm=0d0
        do 7 j=1,n
          s=0d0
          do 6 i=1,m
            s=s+a(i,j)**2
 6        continue
          if(s.gt.tm) tm=s
 7      continue
      endif
      do 100 l=1,nx
        if(norm) small=tm*prec**2
        sm=0d0
        tm=0d0
        do 10 j=l,n
          s=0d0
          sd=0d0
          do 9 i=l,m
c09MAY94    t=a(i,j)**2
c09MAY94    s=s+t*b(i)**2
            s=s+a(i,j)*b(i)
            sd=sd+a(i,j)**2
    9     continue
          s=s**2
          if(sd.gt.tm) tm=sd
          if(norm) s=s/(sd+small)
c         ... normalize the inner product (aj,b) by norm of aj
          if(s.ge.sm)then
            sm=s
            k=j
            sk=sd
          endif
   10   continue
        i=iw(l)
        iw(l)=iw(k)
        iw(k)=i
        do 11 i=1,m
          t=a(i,l)
          a(i,l)=a(i,k)
          a(i,k)=t
   11   continue
        sk=sign(sqrt(sk),a(l,l))
c09MAY94h=sk*(sk+a(l,l))
        h=1d0/(sk*(sk+a(l,l)))
        a(l,l)=a(l,l)+sk
        do 13 j=l+1,n
          s=0d0
          do 12 i=l,m
            s=s+a(i,l)*a(i,j)
   12     continue
          s=s*h
          do i=l,m
            a(i,j)=a(i,j)-a(i,l)*s
          enddo
c09MAY94  x(j)=s/h
   13   continue
c09MAY94 --->
c       do 15 j=l+1,n
c         do 14 i=l,m
c           a(i,j)=a(i,j)-a(i,l)*x(j)
c  14     continue
c  15   continue
c09MAY94 <---
        s=0d0
        do 16 i=l,m
          s=s+a(i,l)*b(i)
   16   continue
c09MAY94s=s/h
        s=s*h
        do 17 i=l,m
          b(i)=b(i)-a(i,l)*s
   17   continue
c       s=0d0
c       do 18 i=l+1,m
c         s=s+b(i)**2
c  18   continue
c       s=sqrt(s/real(m))
c       print *,'l=',l,' iw=',iw(l),' Residual rms: ',sngl(s)
        a(l,l)=-sk
        do 19 i=l+1,m
          a(i,l)=0d0
   19   continue
  100 continue
c     call cputime(ctime1,irtc)
c     call tsolva(a,b,x,nx,nx,mdim,eps)
      if(solvf)then
        call tsvd(a,b,x,eps,svd)
      endif
c     call cputime(ctime2,irtc)
c     print *,'ctime(Htrans)=',sngl(1d-6*(ctime1-ctime0))
c    1       ,'ctime(Tsolvg)=',sngl(1d-6*(ctime2-ctime1))
c
      if(svd .and. solvf) then
        do i=1,nx
          err=0d0
          do j=1,nx
            err=err+(b(j)*a(j,i))**2
          enddo
          a(1,i)=sqrt(err)
        enddo
        do i=1,nx
          b(i)=a(1,i)
        enddo
        do i=nx+1,m
          b(i)=0d0
        enddo
      endif
      do 101 i=nx+1,n
        x(i)=0d0
  101 continue
      if(.not.solvf)then
        do i=1,nx
          x(i)=1.d0
        enddo
      endif
c     call msort2(iw,x,nx)
      call msortn(iw,x,b,nx,1,nx)
      do 102 i=nx,1,-1
        if(iw(i).eq.i) return
        x(iw(i))=x(i)
        x(i)=0d0
        if(svd) then
          b(iw(i))=b(i)
          b(i)=0d0
        endif
  102 continue
      return
      end
