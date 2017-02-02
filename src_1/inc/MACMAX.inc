c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      integer*4 pagesz,inipage
      parameter(pagesz=4096/8,inipage=4)

      integer MAXSTR,MAXMEM0,MAXMEM
      parameter (MAXSTR=255)
      parameter (MAXMEM=2*inipage*pagesz)
      parameter (MAXMEM0=6*1024*pagesz)
