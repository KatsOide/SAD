c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      integer*4 nmult
      parameter (nmult=21)
      logical*4 gknini
      real*8 gkn(0:nmult,0:nmult)
      common /tma/gkn,gknini
