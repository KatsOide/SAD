c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      integer STDERR,STDIN,STDOUT,STDPLT,STDLST
      integer MAXLLEN
      parameter (STDERR=6,STDIN=5,STDOUT=6,STDPLT=8, STDLST=21)
      parameter (MAXLLEN=255)
c
      common /bufbk/pbuf,buf(MAXLLEN)
      common /files/errfl,infl,outfl,pltfl,msglvl,lstfl,inflpt
      integer*4 errfl,infl,outfl,pltfl,msglvl,lstfl,inflpt
      integer*4  pbuf
      character*1 buf
