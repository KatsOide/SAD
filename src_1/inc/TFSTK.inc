c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c This file has been obsolete 9/1/2012 by K.O.
c      integer*4  ivstk(2,RBASE:RBASE+MAXMEM0-1)
c      real*8      vstk(  RBASE:RBASE+MAXMEM0-1)
c      integer*4 itastk(2,RBASE:RBASE+MAXMEM0-1)
c      integer*2 jtastk(4,RBASE:RBASE+MAXMEM0-1)
c      integer*8 ktastk(  RBASE:RBASE+MAXMEM0-1)
c      real*8    rtastk(  RBASE:RBASE+MAXMEM0-1)
c      equivalence ( ivstk(1,RBASE),ilist(1,RBASE))
c      equivalence (  vstk(  RBASE),ilist(1,RBASE))
c      equivalence (itastk(1,RBASE),ilist(1,RBASE))
c      equivalence (jtastk(1,RBASE),ilist(1,RBASE))
c      equivalence (rtastk(  RBASE),ilist(1,RBASE))
c      equivalence (ktastk(  RBASE),ilist(1,RBASE))
c      integer*4 mstk,isp,ivstkoffset,ipurefp,napuref,isporg
c      common /tfstk/mstk,isp,ivstkoffset,ipurefp,napuref,isporg
