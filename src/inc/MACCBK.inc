c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c      include 'inc/MACPARAM.inc'
c
c      character*(MAXPNAME) pname(HTMAX)
c      integer*4 idtype(HTMAX),idval(HTMAX)
c      common /idlist/pname,idtype,idval
c
c      include 'inc/MACMAX.inc'
c      integer*4 ilistroot,ilistdummy(LILISTDUMMY)
c      real*8 rlist
c      integer*4 ilist(2,RBASE:RBASE+MAXMEM0-1)
c      integer*8 klist(RBASE:RBASE+MAXMEM0-1)
c      integer*1 jlist(8,RBASE:RBASE+MAXMEM0-1)
c      common /rmem/ilistroot,ilistdummy,rlist(RBASE:RBASE+MAXMEM0-1)
c      equivalence (rlist(RBASE),ilist(1,RBASE))
c      equivalence (rlist(RBASE),klist(RBASE))
c      equivalence (jlist(1,RBASE),ilist(1,RBASE))
