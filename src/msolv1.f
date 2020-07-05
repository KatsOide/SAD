      subroutine msolv1(a,b,x,na,m,nd,c,d,nc,ndc,cond,micado,nx,norm,
     1                  xc,t,u,c1,svd)
      implicit real*8 (a-h,o-z)
      parameter (epsc=1d-13,tol=1d-13)
      logical cond,micado,norm,svd
      dimension a(nd,m),b(na),x(m),c(ndc,m),d(nc),xc(m),t(nc),u(nc),
     1          c1(nc,m)
      include 'inc/common.inc'
c*** check **********************************************
c     call msolvcheck(a,b,d,x,nd,na,nc,m,.true.)
c********************************************************
      if(cond) then
c       do 10101 j=1,m
c         do 10102 i=1,nc
c           if(abs(c(i,j)).lt.1d-7) c(i,j)=0d0
c0102     continue
c0101   continue
        do 17 l=1,nc
          do 11 i=1,l-1
            s=0d0
            do 10 k=1,m
              s=s+c(l,k)*c(i,k)
   10       continue
            t(i)=s/u(i)
   11     continue
          do 13 k=1,m
            s=0d0
            do 12 i=1,l-1
              s=s+t(i)*c(i,k)
   12       continue
            c(l,k)=c(l,k)-s
   13     continue
c         s=0d0
c         do 14 i=1,l-1
c           s=s+t(i)*d(i)
c  14     continue
c         d(l)=d(l)-s
          s=0d0
          do 15 k=1,m
            s=s+c(l,k)**2
   15     continue
          if(s.lt.m*epsc**2) then
            do 16 k=1,m
              c(l,k)=0d0
   16       continue
            u(l)=1d0
          else
            u(l)=s
          endif
   17   continue
        do 19 j=1,m
          do 18 i=1,nc
            c1(i,j)=c(i,j)
   18     continue
   19   continue
        if(micado) then
          call tmov(d,t,nc)
          call pmicad(c1,d,xc,nx,nc,m,nc,norm,1.d-8,svd,.true.)
          rmin=1d19
          do 21 l=1,nc
            s=0d0
            do 20 j=1,m
              s=s+c(l,j)*xc(j)
   20       continue
            rmin=min(rmin,abs(t(l)-s))
   21     continue
          if(rmin.gt.tol) then
            print *,' ?? Too few variables. (MSOLVG)'
          endif
        else
          call tsolvg(c1,d,xc,nc,m,nc)
        endif
        do 23 l=1,na
          s=0d0
          do 22 k=1,m
            s=s+a(l,k)*xc(k)
   22     continue
          b(l)=b(l)-s
   23   continue
        do 28 l=1,na
          do 25 i=1,nc
            s=0d0
            do 24 k=1,m
              s=s+a(l,k)*c(i,k)
   24       continue
            t(i)=s/u(i)
   25     continue
          do 27 i=1,nc
            do 26 k=1,m
              a(l,k)=a(l,k)-t(i)*c(i,k)
   26       continue
   27     continue
c         b(l)=b(l)-t(i)*d(i)
   28   continue
        if(.not.micado) then
          call tsolva(a,b,x,na,m,nd,eptsol)
        endif
        do 29 i=1,m
          x(i)=x(i)+xc(i)
   29   continue
      elseif(micado) then
        call pmicad(a,b,x,nx,na,m,nd,norm,eptsol,svd,.true.)
      else
        call tsolva(a,b,x,na,m,nd,eptsol)
      endif
c*** check **********************************************
c     call msolvcheck(a,b,d,x,nd,na,nc,m,.false.)
c********************************************************
      return
      end
c
      subroutine msolcheck()
      use tfstk
      use ffs
      use tffitcode
      return
      end
c
      subroutine msolvcheck(a,b,d,x,nd,na,nc,m,save)
      use tfstk
      use ffs
      use tffitcode
      logical save
      real*8 a(nd,m),b(na),d(nc),x(m),resrms
      real*8,allocatable::  a1(:,:),b1(:),d1(:)
      if(save) then
        allocate(a1(nd,m))
        allocate(b1(na))
        allocate(d1(nc))
        a1=a
        b1=b
        d1=d
c        it=ktaloc(nd*m)
c        it1=ktaloc(na+nc)
c        call tmov(a,rlist(it),nd*m)
c        call tmov(b,rlist(it1),na)
c        call tmov(d,rlist(it1+na),nc)
      else
        print *,'x:'
        write(*,'(1p,10g11.3)')(x(i),i=1,m)
        print *,'b-a.x:'
        call mresdue(a1,x,b1,resrms,nd,na+nc,m)
c        write(*,'(1p,10g11.3)') (rlist(it1+i),i=0,na+nc-1)
        print *,' |b-a*x|/sqrt(n) :',resrms
        deallocate(a1,b1,d1)
      endif
      return
      end
c
      subroutine mresdue(a,x,b,resrms,nd,n,m)
      implicit real*8 (a-h,o-z)
      real*8 a(nd,m),b(n),x(m)
c
      resrms=0d0
      do i=1,n
        s=0d0
        do j=1,m
          s=s+a(i,j)*x(j)
        enddo
        b(i)=b(i)-s
        resrms=resrms+b(i)**2
      enddo
      resrms=sqrt(resrms/dble(n))
      return
      end
