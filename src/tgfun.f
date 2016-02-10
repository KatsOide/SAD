      real*8 function tgfun(kf,kp,idp,utwiss,itwissp,
     1           tracex,tracey,pos,geo,nfam)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 kf,kp,idp,itwissp(nlat),nfam
      real*8 pos(nlat),geo(3,4,nlat),utwiss(ntwissfun,-nfam:nfam,1),
     $     tracex(-nfam:nfam),tracey(-nfam:nfam),
     $     tfchi,tphysdisp
      if(kf .le. mfitddp)then
        tgfun=utwiss(kf,idp,itwissp(kp))
      elseif(kf .ge. mfitpex .and. kf .le. mfitpepy)then
        tgfun=tphysdisp(kf,utwiss(1,idp,itwissp(kp)))
      elseif(kf .eq. mfittrx)then
        tgfun=tracex(idp)
      elseif(kf .eq. mfittry)then
        tgfun=tracey(idp)
      elseif(kf .eq. mfitleng)then
        tgfun=pos(kp)
      elseif(kf .le. mfitgz)then
        tgfun=geo(kf-mfitgx+1,4,kp)
      else
        tgfun=tfchi(geo(1,1,kp),kf-mfitchi1+1)
      endif
      return
      end

      real*8 function tphysdisp(kf,utwiss1)
      use ffs
      use tffitcode
      implicit none
      integer*4 kf
      real*8 utwiss1(ntwissfun),cc,r1,r2,r3,r4,detr
c     begin initialize for preventing compiler warning
      tphysdisp=0.d0
c     end   initialize for preventing compiler warning
      r1=utwiss1(mfitr1)
      r2=utwiss1(mfitr2)
      r3=utwiss1(mfitr3)
      r4=utwiss1(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(utwiss1(mfitdetr) .lt. 1.d0)then
        if(kf .eq. mfitpex)then
          tphysdisp=cc*utwiss1(mfitex)
     $         +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
        elseif(kf .eq. mfitpepx)then
          tphysdisp=cc*utwiss1(mfitepx)
     $         -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
        elseif(kf .eq. mfitpey)then
          tphysdisp=cc*utwiss1(mfitey)
     $         -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
        elseif(kf .eq. mfitpepy)then
          tphysdisp=cc*utwiss1(mfitepy)
     $         -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
        endif
      else
        if(kf .eq. mfitpex)then
          tphysdisp=cc*utwiss1(mfitey)
     $         -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
        elseif(kf .eq. mfitpepx)then
          tphysdisp=cc*utwiss1(mfitepy)
     $         -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
        elseif(kf .eq. mfitpey)then
          tphysdisp=cc*utwiss1(mfitex)
     $         +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
        elseif(kf .eq. mfitpepy)then
          tphysdisp=cc*utwiss1(mfitepx)
     $         -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
        endif
      endif
      return
      end

      subroutine tgetphysdisp(twiss,l,pe)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 l
      real*8 twiss(nlat,-ndim:ndim,*),pe(4)
      real*8 cc,r1,r2,r3,r4,detr
      r1=twiss(l,0,mfitr1)
      r2=twiss(l,0,mfitr2)
      r3=twiss(l,0,mfitr3)
      r4=twiss(l,0,mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twiss(l,0,mfitdetr) .lt. 1.d0)then
        pe(1)=cc*twiss(l,0,mfitex)
     $       +r4*twiss(l,0,mfitey)-r2*twiss(l,0,mfitepy)
        pe(2)=cc*twiss(l,0,mfitepx)
     $       -r3*twiss(l,0,mfitey)+r1*twiss(l,0,mfitepy)
        pe(3)=cc*twiss(l,0,mfitey)
     $       -r1*twiss(l,0,mfitex)-r2*twiss(l,0,mfitepx)
        pe(4)=cc*twiss(l,0,mfitepy)
     $       -r3*twiss(l,0,mfitex)-r4*twiss(l,0,mfitepx)
      else
        pe(1)=cc*twiss(l,0,mfitey)
     $       -r1*twiss(l,0,mfitex)-r2*twiss(l,0,mfitepx)
        pe(2)=cc*twiss(l,0,mfitepy)
     $       -r3*twiss(l,0,mfitex)-r4*twiss(l,0,mfitepx)
        pe(3)=cc*twiss(l,0,mfitex)
     $       +r4*twiss(l,0,mfitey)-r2*twiss(l,0,mfitepy)
        pe(4)=cc*twiss(l,0,mfitepx)
     $       -r3*twiss(l,0,mfitey)+r1*twiss(l,0,mfitepy)
      endif
      return
      end

      real*8 function tfchi(geo,i)
      implicit none
      real*8 geo(3,4)
      integer*4 i
      if(i .eq. 2)then
        tfchi=asin(geo(3,3))
      elseif(i .eq. 1)then
        if(geo(2,3) .eq. 0.d0)then
          tfchi=0.d0
        else
          tfchi=atan2(geo(2,3),geo(1,3))
        endif
      else
        if(geo(3,1) .eq. 0.d0)then
          tfchi=0.d0
        else
          tfchi=atan2(geo(3,1),-geo(3,2))
        endif
      endif
      return
      end

      subroutine tfchi2geo(chi1,chi2,chi3,geo)
      implicit none
      real*8 chi1,chi2,chi3,geo(3,3),
     $     schi1,cchi1,
     $     schi2,cchi2,
     $     schi3,cchi3
      cchi1=cos(chi1)
      cchi2=cos(chi2)
      cchi3=cos(chi3)
      schi1=sin(chi1)
      schi2=sin(chi2)
      schi3=sin(chi3)
      geo(1,1)= schi1*cchi3-cchi1*schi2*schi3
      geo(2,1)=-cchi1*cchi3-schi1*schi2*schi3
      geo(3,1)= cchi2*schi3
      geo(1,2)= schi1*schi3+cchi1*schi2*cchi3
      geo(2,2)=-cchi1*schi3+schi1*schi2*cchi3
      geo(3,2)=-cchi2*cchi3
      geo(1,3)= cchi1*cchi2
      geo(2,3)= schi1*cchi2
      geo(3,3)= schi2
      return
      end
