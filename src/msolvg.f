      subroutine  msolvg(a,b,x,n,m,nd,nc,cond,micado,nx,norm,svd)
      use tfstk
      use ffs
      use tffitcode
      logical cond,micado,norm,svd
      dimension a(nd,m), b(n), x(m)
      integer*8 ixc,it,iu,ic1,ia,ib,im
c
      na=n-nc
      if(cond) then
        ixc=ktaloc(m)
        it=ktaloc(nc)
        iu=ktaloc(nc)
        ic1=ktaloc(nc*m)
        if(micado) then
          ia=ktaloc(nd*m)
          ib=ktaloc(n)
          im=ktaloc(m)
          call tmov(a,rlist(ia),nd*m)
          call tmov(b,rlist(ib),n)
          call msolv1(a,b,x,na,m,nd,a(na+1,1),b(na+1),nc,nd,.false.,
     1         micado,nx,norm,rlist(ixc),rlist(it),rlist(iu),
     1         rlist(ic1),svd)
          m1=0
          do i=1,m
            if(x(i).ne.0d0)then
              m1=m1+1
              ilist(1,im-1+m1)=i
              call tmov(rlist(ia+(i-1)*nd),a(1,m1),nd)
            endif
          enddo
          call tmov(rlist(ib),b,n)
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
          call tfree(ib)
          call tfree(ia)
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
