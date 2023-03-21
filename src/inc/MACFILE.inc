c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer STDERR,STDIN,STDOUT,STDPLT,STDLST
c$$$      integer MAXLLEN
c$$$      parameter (STDERR=6,STDIN=5,STDOUT=6,STDPLT=8, STDLST=21)
c$$$      parameter (MAXLLEN=255)
c$$$c
c$$$      common /bufbk/pbuf,buf(MAXLLEN)
c$$$      common /files/errfl,infl,outfl,pltfl,msglvl,lstfl,inflpt
c$$$      integer*4 errfl,infl,outfl,pltfl,msglvl,lstfl
c$$$      integer*8 inflpt
c$$$      integer*4  pbuf
c$$$      character*1 buf
