c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*4 maxcond,lblname,mcommon
c$$$      parameter (maxcond=4001)
c$$$      parameter (lblname=MAXCOND/4)
c$$$      parameter (mcommon=maxcond*14+lblname+4*2+7+1)
c$$$      integer*8 iut
c$$$      integer*4 mfitp(maxcond),ifitp(maxcond),ifitp1(maxcond),
c$$$     $     kdp(maxcond),kfitp(maxcond),kfit(maxcond),
c$$$     $     icalc(3,maxcond),iqcol(maxcond),lfp(2,maxcond)
c$$$      real*8 ctime2,ctime0,rsconv,fitval(maxcond)
c$$$      integer*4 nvar,ntouch,itmax,measp,nfc,ncalc,
c$$$     $     blname(lblname),pading
c$$$      common /ffslocal/ fitval,ctime2,ctime0,rsconv,iut,
c$$$     $     mfitp,ifitp,ifitp1,kdp,kfitp,kfit,
c$$$     $     icalc,iqcol,lfp,blname,nvar,ntouch,
c$$$     $     itmax,measp,nfc,ncalc,pading
