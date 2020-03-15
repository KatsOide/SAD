      subroutine tbal(w,a,v,n,ndim)
      implicit none
c   radix is the machine radix of floating.
      real*8 , parameter ::radix=2.d0,sqrdx=radix**2
      integer*4 n,ndim,i,j
      real*8 a(n,n),w(ndim,n),v(n)
      real*8 c,r,g,f,s,u
      logical*4 last
      a=w(1:n,:)
      w(1:n,:)=0.d0
      v=1.d0
      do 1 i=1,n
        w(i,i)=1.d0
1     continue
      last=.true.
      do while(last)
        last=.false.
        do i=1,n
          c=0.d0
          r=0.d0
          do j=1,n
            if(j .ne. i)then
              u=v(j)/v(i)
              c=c+abs(a(j,i)*u)
              r=r+abs(a(i,j)/u)
            endif
          enddo
          if(c .ne. 0.d0 .and. r .ne. 0.d0)then
            g=r/radix
            f=1.d0
            s=c+r
            do while(c .lt. g)
              f=f*radix
              c=c*sqrdx
            enddo
            g=r*radix
            do while(c .gt. g)
              f=f/radix
              c=c/sqrdx
            enddo
            if((c+r)/f .lt. .95d0*s)then
              last=.true.
              v(i)=v(i)/f
            endif
          endif
        enddo
      enddo
      return
      end
