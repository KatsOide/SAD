      subroutine  msolvg(a,b,x,n,m,nd,nc,cond,micado,nx,norm,svd)
      use tfstk
      use ffs
      use tffitcode
      logical cond,micado,norm,svd
      dimension a(nd,m), b(n), x(m)
c
      na=n-nc
      if(cond) then
        ixc=italoc(m)
        it=italoc(nc)
        iu=italoc(nc)
        ic1=italoc(nc*m)
        if(micado) then
          ia=italoc(nd*m)
          ib=italoc(n)
          im=italoc(m)
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
          call tfree(int8(im))
          call tfree(int8(ib))
          call tfree(int8(ia))
        else
          call msolv1(a,b,x,na,m,nd,a(na+1,1),b(na+1),nc,nd,cond,
     1         micado,nx,norm,rlist(ixc),rlist(it),rlist(iu),
     1         rlist(ic1),svd)
        endif
        call tfree(int8(ic1))
        call tfree(int8(iu))
        call tfree(int8(it))
        call tfree(int8(ixc))
      else
        call msolv1(a,b,x,n,m,nd,a,b,nc,nd,cond,
     1              micado,nx,norm,b,b,b,a,svd)
      endif
c
      return
      end
