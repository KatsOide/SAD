c$$$c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*8 iwakepold
c$$$      integer*4 nflag0,nlat,np0,nturn,isynch,lspect,nspect,
c$$$     $     lplot,nplot,nuse,nclas,irad,novfl,
c$$$     $     npelm,ipelm,ipoltr,ipolb,ipoll,ipolid,ipolo,
c$$$     $     ilattp,nparallel,pspac_nx,pspac_ny,pspac_nz,
c$$$     $     pspac_nturn,pspac_nturncalc
c$$$      include 'inc/MACPHYS.inc'
c$$$      real*8 c,hp,e
c$$$      parameter (c=cveloc,hp=plankr,e=elemch)
c$$$      parameter (nflag0=64)
c$$$      include 'inc/MACCODE.inc'
c$$$      include 'inc/MACKW.inc'
c$$$      include 'inc/MACFILE.inc'
c$$$      logical*4 oldflagsdummy,calint,caltouck,tparaed
c$$$      real*8
c$$$     $     amass,charge,h0,p0,omega0,trf0,crad,erad,epsrad,
c$$$     $     codin(6),dleng,anrad,urad,u0,vc0,wrfeff,dp0,brho,
c$$$     $     ccintr,cintrb,pbunch,coumin,re0,pgev,emidiv,
c$$$     $     emidib,emidiq,emidis,ctouck,dvemit,h1emit,
c$$$     $     anbunch,tdummy(6),zlost,alost,
c$$$     $     taurdx,taurdy,taurdz,fridiv,beamin(21),
c$$$     $     vccos,vcsin,vcphic,vcalpha,vceff,
c$$$     $     vcacc,dvcacc,ddvcacc,
c$$$     $     pspac_dx,pspac_dy,pspac_dz,dvfs,rcratio,rclassic,brhoz,
c$$$     $     dummycom(111)
c$$$      common /tparm/
c$$$     $     amass,charge,h0,p0,omega0,trf0,crad,erad,epsrad,
c$$$     $     codin,dleng,anrad,urad,u0,vc0,wrfeff,dp0,brho,
c$$$     $     ccintr,cintrb,pbunch,coumin,re0,pgev,emidiv,
c$$$     $     emidib,emidiq,emidis,ctouck,dvemit,h1emit,
c$$$     $     anbunch,tdummy,zlost,alost,
c$$$     $     taurdx,taurdy,taurdz,fridiv,beamin,
c$$$     $     vccos,vcsin,vcphic,vcalpha,vceff,
c$$$     $     vcacc,dvcacc,ddvcacc,
c$$$     $     pspac_dx,pspac_dy,pspac_dz,dvfs,rcratio,rclassic,brhoz,
c$$$     $     dummycom
c$$$c      integer*4 ntouckl,ntouckx,ntouckz
c$$$c      parameter (ntouckl=120,ntouckx=34,ntouckz=25)
c$$$c      real*8
c$$$c     $     touckm(ntouckz,ntouckx,3),tampl(ntouckx,3),touckl(ntouckl),
c$$$c     $     toucke(ntouckl,6000)
c$$$c      common /intr/ touckm,tampl,touckl,toucke
c$$$      common /itparm/ iwakepold,
c$$$     $     nlat,np0,nturn,isynch,lspect,nspect,lplot,nplot,
c$$$     $     nuse,nclas,irad,oldflagsdummy(32),
c$$$     $     novfl,calint,
c$$$     $     npelm,ipelm,ipoltr,ipolb,ipoll,ipolid,ipolo,
c$$$     $     caltouck,ilattp,tparaed,nparallel,
c$$$     $     pspac_nx,pspac_ny,pspac_nz,pspac_nturn,pspac_nturncalc
