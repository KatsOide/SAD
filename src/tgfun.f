      real*8 function tgfun(kf,kp,idp)
      use tfstk
      use ffs
      use ffs_fit
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 kf,kp,idp
      real*8 tfchi,tphysdisp
      if(kf .le. mfitzpy)then
        tgfun=utwiss(kf,idp,itwissp(kp))
      elseif(kf .ge. mfitpex .and. kf .le. mfitpepy)then
        tgfun=tphysdisp(kf,utwiss(1,idp,itwissp(kp)))
      elseif(kf .eq. mfittrx)then
        tgfun=optstat(idp)%tracex
      elseif(kf .eq. mfittry)then
        tgfun=optstat(idp)%tracey
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
        elseif(kf .eq. mfitpzx)then
          tphysdisp=cc*utwiss1(mfitzx)
     $         +r4*utwiss1(mfitzy)-r2*utwiss1(mfitzpy)
        elseif(kf .eq. mfitpzpx)then
          tphysdisp=cc*utwiss1(mfitzpx)
     $         -r3*utwiss1(mfitzy)+r1*utwiss1(mfitzpy)
        elseif(kf .eq. mfitpzy)then
          tphysdisp=cc*utwiss1(mfitzy)
     $         -r1*utwiss1(mfitzx)-r2*utwiss1(mfitzpx)
        elseif(kf .eq. mfitpzpy)then
          tphysdisp=cc*utwiss1(mfitzpy)
     $         -r3*utwiss1(mfitzx)-r4*utwiss1(mfitzpx)
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
        elseif(kf .eq. mfitpzx)then
          tphysdisp=cc*utwiss1(mfitzy)
     $         -r1*utwiss1(mfitzx)-r2*utwiss1(mfitzpx)
        elseif(kf .eq. mfitpzpx)then
          tphysdisp=cc*utwiss1(mfitzpy)
     $         -r3*utwiss1(mfitzx)-r4*utwiss1(mfitzpx)
        elseif(kf .eq. mfitpzy)then
          tphysdisp=cc*utwiss1(mfitzx)
     $         +r4*utwiss1(mfitzy)-r2*utwiss1(mfitzpy)
        elseif(kf .eq. mfitpzpy)then
          tphysdisp=cc*utwiss1(mfitzpx)
     $         -r3*utwiss1(mfitzy)+r1*utwiss1(mfitzpy)
        endif
      endif
      return
      end

      subroutine tgetphysdisp(l,pe)
      implicit none
      integer*4 l
      real*8 pe(4)
      call tgetphysdispi(l,0,pe)
      return
      end

      subroutine tgetphysdispi(l,icol,pe)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 l,icol
      real*8 pe(4)
      real*8 cc,r1,r2,r3,r4,detr
      r1=twiss(l,icol,mfitr1)
      r2=twiss(l,icol,mfitr2)
      r3=twiss(l,icol,mfitr3)
      r4=twiss(l,icol,mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twiss(l,icol,mfitdetr) .lt. 1.d0)then
        pe(1)=cc*twiss(l,icol,mfitex)
     $       +r4*twiss(l,icol,mfitey)-r2*twiss(l,icol,mfitepy)
        pe(2)=cc*twiss(l,icol,mfitepx)
     $       -r3*twiss(l,icol,mfitey)+r1*twiss(l,icol,mfitepy)
        pe(3)=cc*twiss(l,icol,mfitey)
     $       -r1*twiss(l,icol,mfitex)-r2*twiss(l,icol,mfitepx)
        pe(4)=cc*twiss(l,icol,mfitepy)
     $       -r3*twiss(l,icol,mfitex)-r4*twiss(l,icol,mfitepx)
      else
        pe(1)=cc*twiss(l,icol,mfitey)
     $       -r1*twiss(l,icol,mfitex)-r2*twiss(l,icol,mfitepx)
        pe(2)=cc*twiss(l,icol,mfitepy)
     $       -r3*twiss(l,icol,mfitex)-r4*twiss(l,icol,mfitepx)
        pe(3)=cc*twiss(l,icol,mfitex)
     $       +r4*twiss(l,icol,mfitey)-r2*twiss(l,icol,mfitepy)
        pe(4)=cc*twiss(l,icol,mfitepx)
     $       -r3*twiss(l,icol,mfitey)+r1*twiss(l,icol,mfitepy)
      endif
      return
      end

      subroutine tgetphysdispz(l,pe)
      implicit none
      integer*4 l
      real*8 pe(4)
      call tgetphysdispzi(l,0,pe)
      return
      end

      subroutine tgetphysdispzi(l,icol,pe)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 l,icol
      real*8 pe(4)
      real*8 cc,r1,r2,r3,r4,detr
      r1=twiss(l,icol,mfitr1)
      r2=twiss(l,icol,mfitr2)
      r3=twiss(l,icol,mfitr3)
      r4=twiss(l,icol,mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twiss(l,icol,mfitdetr) .lt. 1.d0)then
        pe(1)=cc*twiss(l,icol,mfitzx)
     $       +r4*twiss(l,icol,mfitzy)-r2*twiss(l,icol,mfitzpy)
        pe(2)=cc*twiss(l,icol,mfitzpx)
     $       -r3*twiss(l,icol,mfitzy)+r1*twiss(l,icol,mfitzpy)
        pe(3)=cc*twiss(l,icol,mfitzy)
     $       -r1*twiss(l,icol,mfitzx)-r2*twiss(l,icol,mfitzpx)
        pe(4)=cc*twiss(l,icol,mfitzpy)
     $       -r3*twiss(l,icol,mfitzx)-r4*twiss(l,icol,mfitzpx)
      else
        pe(1)=cc*twiss(l,icol,mfitzy)
     $       -r1*twiss(l,icol,mfitzx)-r2*twiss(l,icol,mfitzpx)
        pe(2)=cc*twiss(l,icol,mfitzpy)
     $       -r3*twiss(l,icol,mfitzx)-r4*twiss(l,icol,mfitzpx)
        pe(3)=cc*twiss(l,icol,mfitzx)
     $       +r4*twiss(l,icol,mfitzy)-r2*twiss(l,icol,mfitzpy)
        pe(4)=cc*twiss(l,icol,mfitzpx)
     $       -r3*twiss(l,icol,mfitzy)+r1*twiss(l,0,mfitzpy)
      endif
      return
      end

      subroutine tgetphysdispu(utwiss1,pe)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      real*8 utwiss1(ntwissfun),pe(4)
      real*8 cc,r1,r2,r3,r4,detr
      r1=utwiss1(mfitr1)
      r2=utwiss1(mfitr2)
      r3=utwiss1(mfitr3)
      r4=utwiss1(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(utwiss1(mfitdetr) .lt. 1.d0)then
        pe(1)=cc*utwiss1(mfitex)
     $       +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
        pe(2)=cc*utwiss1(mfitepx)
     $       -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
        pe(3)=cc*utwiss1(mfitey)
     $       -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
        pe(4)=cc*utwiss1(mfitepy)
     $       -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
      else
        pe(1)=cc*utwiss1(mfitey)
     $       -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
        pe(2)=cc*utwiss1(mfitepy)
     $       -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
        pe(3)=cc*utwiss1(mfitex)
     $       +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
        pe(4)=cc*utwiss1(mfitepx)
     $       -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
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
