      subroutine tbal(w,a,v,n,ndim)
      implicit none
      real*8 radix,sqrdx
c   radix is the machine radix of floating.
      parameter (radix=2.d0,sqrdx=radix**2)
      integer*4 n,ndim,i,j
      real*8 a(n,n),w(ndim,n),v(n)
      real*8 c,r,g,f,s,u
      logical*4 last
      do 1 i=1,n
        do 2 j=1,n
          a(j,i)=w(j,i)
          w(j,i)=0.d0
2       continue
        w(i,i)=1.d0
        v(i)=1.d0
1     continue
100   last=.false.
      do 10 i=1,n
        c=0.d0
        r=0.d0
        do 20 j=1,n
          if(j .ne. i)then
            u=v(j)/v(i)
            c=c+abs(a(j,i)*u)
            r=r+abs(a(i,j)/u)
          endif
20      continue
        if(c .ne. 0.d0 .and. r .ne. 0.d0)then
          g=r/radix
          f=1.d0
          s=c+r
200       if(c .lt. g)then
            f=f*radix
            c=c*sqrdx
            go to 200
          endif
          g=r*radix
300       if(c .gt. g)then
            f=f/radix
            c=c/sqrdx
            go to 300
          endif
          if((c+r)/f .lt. .95d0*s)then
            last=.true.
            v(i)=v(i)/f
          endif
        endif
10    continue
      if(last)then
        go to 100
      endif
      return
      end
