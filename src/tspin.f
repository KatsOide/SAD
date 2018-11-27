      module tspin
      use macphys
      real*8, parameter :: gspin=finest/m_2pi

      type spin
      sequence
      real*8 sx,sy,sz
      end type


      subroutine textractb(px0,py0,px1,py1,dp,bxds,byds)
      use tmacro, only:p0
      implicit none
      real*8 px0,py0,px1,py1,dp,bxds,byds,pxc,pyc,dpx,dpy,p,v
      pxc=(px1+px0)*.5d0
      pyc=(py1+py0)*.5d0
      dpx=px1-px0
      dpy=py1-py0
      p=(1.d0+dp)*p0
      v=p/sqrt(1.d0+p^2)
      bxds=-dpy*p
      byds=dpx*p
%%%%%%%%


      end module



