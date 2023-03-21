c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer*4
c$$$     $     mfitax,mfitbx,mfitnx,mfitay,mfitby,mfitny,
c$$$     $     mfitex,mfitepx,mfitey,mfitepy,
c$$$     $     mfitr1,mfitr2,mfitr3,mfitr4,
c$$$     $     mfitdx,mfitdpx,mfitdy,mfitdpy,mfitdz,mfitddp,
c$$$     $     mfitpex,mfitpepx,mfitpey,mfitpepy,
c$$$     $     mfittrx,mfittry,mfitleng,
c$$$     $     mfitgx,mfitgy,mfitgz,mfitchi1,mfitchi2,mfitchi3
c$$$      parameter (
c$$$     $     mfitax=1,mfitbx=mfitax+1,mfitnx=mfitbx+1,
c$$$     $     mfitay=mfitnx+1,mfitby=mfitay+1,mfitny=mfitby+1,
c$$$     $     mfitex=mfitny+1,mfitepx=mfitex+1,
c$$$     $     mfitey=mfitepx+1,mfitepy=mfitey+1,
c$$$     $     mfitr1=mfitepy+1,mfitr2=mfitr1+1,
c$$$     $     mfitr3=mfitr2+1,mfitr4=mfitr3+1,
c$$$     $     mfitdx=mfitr4+1,mfitdpx=mfitdx+1,
c$$$     $     mfitdy=mfitdpx+1,mfitdpy=mfitdy+1,
c$$$     $     mfitdz=mfitdpy+1,mfitddp=mfitdz+1,
c$$$     $     mfitpex=mfitddp+1,mfitpepx=mfitpex+1,
c$$$     $     mfitpey=mfitpepx+1,mfitpepy=mfitpey+1,
c$$$     $     mfittrx=mfitpepy+1,mfittry=mfittrx+1,mfitleng=mfittry+1,
c$$$     $     mfitgx=mfitleng+1,mfitgy=mfitgx+1,mfitgz=mfitgy+1,
c$$$     $     mfitchi1=mfitgz+1,mfitchi2=mfitchi1+1,mfitchi3=mfitchi2+1)
c$$$      integer*4 ntwissfun,mfito,mfit,mfit1
c$$$      parameter (ntwissfun=mfitddp,mfito=mfittry,mfit=mfitchi3,
c$$$     $     mfit1=mfit+12)
