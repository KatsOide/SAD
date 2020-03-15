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
      end

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
      end
