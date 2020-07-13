      subroutine  msolvg(a,b,x,n,m,nd,nc,cond,micado,nx,norm,svd)
      use tfstk
      use ffs
      use tffitcode
      logical cond,micado,norm,svd
      real*8  a(nd,m), b(n), x(m)
      real*8 ,allocatable::  a1(:,:), b1(:), x1(:)
      integer*8 ixc,it,iu,ic1,im
c
      na=n-nc
      if(cond) then
        ixc=ktaloc(m)
        it=ktaloc(nc)
        iu=ktaloc(nc)
        ic1=ktaloc(nc*m)
        if(micado) then
          allocate(a1(nd,m))
          allocate(b1(n))
          allocate(x1(m))
c          ia=ktaloc(nd*m)
c          ib=ktaloc(n)
          im=ktaloc(m)
          a1=a
          b1=b
c          call tmov(a,rlist(ia),nd*m)
c          call tmov(b,rlist(ib),n)
          call msolv1(a,b,x,na,m,nd,a(na+1,1),b(na+1),nc,nd,.false.,
     1         micado,nx,norm,rlist(ixc),rlist(it),rlist(iu),
     1         rlist(ic1),svd)
          m1=0
          do i=1,m
            if(x(i).ne.0d0)then
              m1=m1+1
              ilist(1,im-1+m1)=i
              a=a1
c              call tmov(rlist(ia+(i-1)*nd),a(1,m1),nd)
            endif
          enddo
          b=b1
c          call tmov(rlist(ib),b,n)
          call msolv1(a,b,x,na,m1,nd,a(na+1,1),b(na+1),nc,nd,.true.,
     1         .false.,nx,norm,rlist(ixc),rlist(it),rlist(iu),
     1         rlist(ic1),svd)
          do i=m1,1,-1
            j=ilist(1,im-1+i)
            if(j.ne.i)then
              x(j)=x(i)
              x(i)=0d0
            endif
          enddo
          call tfree(im)
c          call tfree(ib)
c          call tfree(ia)
        else
          call msolv1(a,b,x,na,m,nd,a(na+1,1),b(na+1),nc,nd,cond,
     1         micado,nx,norm,rlist(ixc),rlist(it),rlist(iu),
     1         rlist(ic1),svd)
        endif
        call tfree(ic1)
        call tfree(iu)
        call tfree(it)
        call tfree(ixc)
      else
        call msolv1(a,b,x,n,m,nd,a,b,nc,nd,cond,
     1              micado,nx,norm,b,b,b,a,svd)
      endif
c
      return
      end
