c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      integer*4 FLAGON,FLAGOF
      parameter (FLAGON=-1,FLAGOF=0)
c
      integer*4 HTMAX,MAXPNAME,LILISTDUMMY
      parameter(MAXPNAME=32,LILISTDUMMY=3)
      character*(MAXPNAME) NULSTR
      parameter(HTMAX=32767,NULSTR='        ')
c     define `RBASE': base index number in pfalloc.f
      integer*4 RBASE
      parameter (RBASE=1513*1024)
