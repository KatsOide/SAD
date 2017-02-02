c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*4 irbinit,irbopen,irbclose,irbreadrecord,irbreadbuf,
c$$$     $     irbmovepoint,irbeor2bor,irbgetpoint,irbreset,
c$$$     $     irbreadrecordbuf,irbbor,irbsetinp,irbcloseinp,
c$$$     $     irbsetbuf,irbibuf
c$$$      parameter (irbinit=1,irbopen=2,irbclose=3,irbreadrecord=4,
c$$$     $     irbreadbuf=5,irbmovepoint=6,irbbor=7,irbgetpoint=8,
c$$$     $     irbreset=9,irbreadrecordbuf=10,irbeor2bor=11,
c$$$     $     irbsetinp=12,irbcloseinp=13,irbsetbuf=14,irbibuf=15)
c$$$      integer*4 ncprolog
c$$$      character*128 prolog
c$$$      common /shortlog/ ncprolog,prolog
c$$$      integer*4 nbuf
c$$$      parameter (nbuf=1024)
c$$$      integer*4 lbuf(nbuf),mbuf(nbuf)
c$$$      integer*8 ibuf(nbuf),itbuf(nbuf)
c$$$      common /rbufcom/ ibuf,itbuf,lbuf,mbuf
