      subroutine teigen(a,w,eig,n,ndim)
c
c   Subroutine to obtain eigen values and eigen vectors of a real
c   matrix a.
c                                7-Dec-1987   K. Oide
c   Modified                    27-Oct-1989   K. Oide
c   Modified                    17-Nov-1989   K. Oide
c   Added Balancing              1-Nov-1990   K. Oide
c   Modified Convergence        20-Nov-1990   K. Oide
c   Modified Shifting            4-Dec-1990   K. Oide
c   Fast Givens' Transformation 21-Jan-1992   K. Oide
c   Maximum size 14,000         29-May-2010   K. Oide
c   Maximum size 32768           5-Oct-2010   K. Oide
c
c   Usage:       call teigen(a,w,eig,n,ndim)
c
c                where
c
c                a     is the input n*n matrixand the eigen vectors
c                      are returned here.
c                w     is a work area of n**2 words.
c                eig   is a (2,n) real array where the eigen values
c                      are returned.   eig(1,i) and eig(2,i) contain
c                      the real and imaginary part of the i-th eigen-
c                      valeu respectively.   When the eigenvalue is
c                      a complex, the output are stored as
c
c                      eig(1,i  ) =      er
c                      eig(2,i  ) =     !ei!
c                      eig(1,i+1) =      er
c                      eig(2,i+1) =    -!ei!
c
c                      where er and ei are the real and imaginary
c                      part of the i-th eigenvalue respectively.
c                n     is the matrix dimension.
c                ndim  is the size of the first dimension of a.
c
c   Restriction: When the matrix cannot be diagonalized (Jordan
c                standard form), the returned eigen vectors are not
c                correct.
c
      implicit none
      integer*4 n,ndim,nmax
      parameter (nmax=32768)
      integer*4 ia(nmax)
      real*8 a(ndim,n),w(n,n),eig(2,n),vx(nmax)
      if(n .gt. nmax)then
        write(*,*)'TEIGEN-Too large matrix.',n
      else
        call tbal(a,w,vx,n,ndim)
        call thess(w,a,vx,n,ndim)
        call tqr(w,a,eig,ia,vx,n,ndim)
      endif
      return
      end
