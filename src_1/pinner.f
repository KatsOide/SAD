      subroutine pinner(a,b,c,n,m,ndim)
      implicit none
C*DEC added by Y.Tange 10-Jan-1995
      integer*4 n,m,ndim,i,j
      real*8 a(ndim,m),b(m),c(n),s
C*DEC End
C*HP
C     real*8 a(ndim,m),b(m),c(n),s
C     integer*4 n,m,ndim,i,j
C*HP
      do 11 i=1,n
        s=0d0
        do 10 j=1,m
          s=s+a(i,j)*b(j)
   10   continue
        c(i)=s
   11 continue
      return
      end
