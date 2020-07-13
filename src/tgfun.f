      real*8 function tgfun(kf,kp,idp)
      use tfstk
      use ffs
      use ffs_fit
      use ffs_pointer
      use tffitcode
      use geolib
      implicit none
      integer*4 kf,kp,idp
      real*8 tphysdisp
      if(kf .le. mfitzpy)then
        tgfun=utwiss(kf,idp,itwissp(kp))
      elseif(kf .ge. mfitpex .and. kf .le. mfitgmz)then
        tgfun=tphysdisp(kf,utwiss(1:ntwissfun,idp,itwissp(kp)))
      else
        select case (kf)
        case (mfittrx)
          tgfun=optstat(idp)%tracex
        case (mfittry)
          tgfun=optstat(idp)%tracey
        case (mfitleng)
          tgfun=pos(kp)
        case default
          if(kf .le. mfitgz)then
            tgfun=geo(kf-mfitgx+1,4,kp)
          else
            tgfun=tfchi(geo(:,:,kp),kf-mfitchi1+1)
          endif
        end select
      endif
      return
      end

      real*8 function tphysdisp(kf,utwiss1)
      use ffs
      use tffitcode
      implicit none
      integer*4 , intent(in)::kf
      real*8 ,intent(in):: utwiss1(ntwissfun)
      real*8 cc,r1,r2,r3,r4,detr
c     begin initialize for preventing compiler warning
c      tphysdisp=0.d0
c     end   initialize for preventing compiler warning
      if(kf .le. mfitpzpy)then
        r1=utwiss1(mfitr1)
        r2=utwiss1(mfitr2)
        r3=utwiss1(mfitr3)
        r4=utwiss1(mfitr4)
        detr=r1*r4-r2*r3
        cc=sqrt(1.d0-detr)
        if(utwiss1(mfitdetr) .lt. 1.d0)then
          select case (kf)
          case (mfitpex)
            tphysdisp=cc*utwiss1(mfitex)
     $           +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
          case (mfitpepx)
            tphysdisp=cc*utwiss1(mfitepx)
     $           -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
          case (mfitpey)
            tphysdisp=cc*utwiss1(mfitey)
     $           -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
          case (mfitpepy)
            tphysdisp=cc*utwiss1(mfitepy)
     $           -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
          case (mfitpzx)
            tphysdisp=cc*utwiss1(mfitzx)
     $           +r4*utwiss1(mfitzy)-r2*utwiss1(mfitzpy)
          case (mfitpzpx)
            tphysdisp=cc*utwiss1(mfitzpx)
     $           -r3*utwiss1(mfitzy)+r1*utwiss1(mfitzpy)
          case (mfitpzy)
            tphysdisp=cc*utwiss1(mfitzy)
     $           -r1*utwiss1(mfitzx)-r2*utwiss1(mfitzpx)
          case (mfitpzpy)
            tphysdisp=cc*utwiss1(mfitzpy)
     $           -r3*utwiss1(mfitzx)-r4*utwiss1(mfitzpx)
          case default
            tphysdisp=0.d0
          end select
        else
          select case (kf)
          case(mfitpex)
            tphysdisp=cc*utwiss1(mfitey)
     $           -r1*utwiss1(mfitex)-r2*utwiss1(mfitepx)
          case (mfitpepx)
            tphysdisp=cc*utwiss1(mfitepy)
     $           -r3*utwiss1(mfitex)-r4*utwiss1(mfitepx)
          case (mfitpey)
            tphysdisp=cc*utwiss1(mfitex)
     $           +r4*utwiss1(mfitey)-r2*utwiss1(mfitepy)
          case (mfitpepy)
            tphysdisp=cc*utwiss1(mfitepx)
     $           -r3*utwiss1(mfitey)+r1*utwiss1(mfitepy)
          case (mfitpzx)
            tphysdisp=cc*utwiss1(mfitzy)
     $           -r1*utwiss1(mfitzx)-r2*utwiss1(mfitzpx)
          case (mfitpzpx)
            tphysdisp=cc*utwiss1(mfitzpy)
     $           -r3*utwiss1(mfitzx)-r4*utwiss1(mfitzpx)
          case (mfitpzy)
            tphysdisp=cc*utwiss1(mfitzx)
     $           +r4*utwiss1(mfitzy)-r2*utwiss1(mfitzpy)
          case (mfitpzpy)
            tphysdisp=cc*utwiss1(mfitzpx)
     $           -r3*utwiss1(mfitzy)+r1*utwiss1(mfitzpy)
          case default
            tphysdisp=0.d0
          end select
        endif
      else
c        write(*,*)'tphysdisp ',kf
        select case (kf)
        case (mfitgmx)
          tphysdisp=(1.d0+utwiss1(mfitax)**2)/utwiss1(mfitbx)
        case (mfitgmy)
          tphysdisp=(1.d0+utwiss1(mfitay)**2)/utwiss1(mfitby)
        case (mfitgmz)
          tphysdisp=(1.d0+utwiss1(mfitaz)**2)/utwiss1(mfitbz)
        case default
          tphysdisp=0.d0
        end select
      endif
      return
      end

      real*8 function tgetgm(kf,l,idp)
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 , intent(in)::l,kf,idp
      select case (kf)
      case (mfitgmx)
        tgetgm=(1.d0+twiss(l,idp,mfitax)**2)/twiss(l,idp,mfitbx)
      case (mfitgmy)
        tgetgm=(1.d0+twiss(l,idp,mfitay)**2)/twiss(l,idp,mfitby)
      case (mfitgmz)
        tgetgm=(1.d0+twiss(l,idp,mfitaz)**2)/twiss(l,idp,mfitbz)
      case default
        tgetgm=0.d0
      end select
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
