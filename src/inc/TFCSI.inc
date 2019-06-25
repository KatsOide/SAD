c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*4 ipoint,iconv,ios,lfn1,lfni,lfno,lrecl,ldel,
c$$$     $     lcmnt,linep,lastln,ibegt,lastt,nbmax
c$$$      parameter (nbmax=200000)
c$$$      character*16 delim,cmnt
c$$$      character*(nbmax) buffer
c$$$      logical*4 rec
c$$$      common /csi/ ipoint,iconv,ios,lfn1,lfni,lfno
c$$$      common /csbuf/ lrecl,buffer
c$$$      common /csidel/ ldel,delim
c$$$      common /csicm/  lcmnt,cmnt
c$$$      common /csirec/ rec,linep,lastln,ibegt,lastt
