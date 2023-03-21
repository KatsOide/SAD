c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      include 'inc/TMACRO1.inc'
c$$$      integer*4 ndim,ndima,nele,nfit,marki,iorgx,iorgy,iorgr,
c$$$     $     mfpnt,mfpnt1,id1,id2,nve,ielmhash,iscal
c$$$      integer*8 ifaux,ifibzl,ifmult,ifklp,
c$$$     $     ifival,iftwissp,iftwis,ifpos,ifgeo,ifsize,
c$$$     $     ifgamm,ifele,ifcoup,iferrk,ifvarele,ifvvar,
c$$$     $     ifvalvar,ifele1,ifele2,ifmast,iftouchele,
c$$$     $     iftouchv,lfnp,iffserr,ifwakep,ifivcomp
c$$$      integer*4 ndim1,nscal,niscal
c$$$      parameter (ndim1=201)
c$$$      parameter (nscal=64,niscal=64)
c$$$      integer*4 nelmhash
c$$$      parameter (nelmhash=1023)
c$$$      real*8 twissi,emx,emy,dpmax,geo0,xixf,xiyf,scal
c$$$      common /ffs/ ifaux,ifibzl,ifmult,ifklp,
c$$$     $     ifival,iftwissp,iftwis,ifpos,ifgeo,ifsize,
c$$$     $     ifgamm,ifele,ifcoup,iferrk,ifvarele,ifvvar,
c$$$     $     ifvalvar,ifele1,ifele2,ifmast,iftouchele,
c$$$     $     iftouchv,lfnp,iffserr,ifwakep,ifivcomp,
c$$$     $     twissi(ntwissfun),emx,emy,dpmax,geo0(3,4),
c$$$     $     xixf,xiyf,scal(nscal),iscal(niscal),
c$$$     $     ndim,ndima,nele,nfit,marki,iorgx,iorgy,iorgr,
c$$$     $     mfpnt,mfpnt1,id1,id2,nve,ielmhash
