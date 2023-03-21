c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*4 mpsize,nindex,mhash,minseg0,minseg1,minseg2
c$$$      parameter (mpsize=2**15,nindex=64,mhash=32767,
c$$$     $     minseg0=9,minseg1=16,minseg2=16)
c$$$      integer*8 icp,nmem,nnet,ich,maxic,minic,lastpend,
c$$$     $     icsep
c$$$      common /itcache/icp,nmem,nnet,ich,maxic,minic,lastpend,icsep
