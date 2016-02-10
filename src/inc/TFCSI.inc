c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      integer*4 ipoint,iconv,ios,lfn1,lfni,lfno,lrecl,ldel,
     $     lcmnt,linep,lastln,ibegt,lastt,nbmax
      parameter (nbmax=200000)
      character*16 delim,cmnt
      character*(nbmax) buffer
      logical*4 rec
      common /csi/ ipoint,iconv,ios,lfn1,lfni,lfno
      common /csbuf/ lrecl,buffer
      common /csidel/ ldel,delim
      common /csicm/  lcmnt,cmnt
      common /csirec/ rec,linep,lastln,ibegt,lastt
