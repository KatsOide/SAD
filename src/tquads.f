c  Obsolete Aug. 10, 2020
c
      subroutine tquads(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak0,bz,
     1     dx,dy,theta,theta2,radlvl,fringe,f1in,f2in,f1out,f2out,
     $     mfring,ibsi,eps0)
      use tfstk
      use ffs_flag
      use tmacro
      use photontable,only:tsetpcvt,pcvt
      use kradlib
      use sol,only:tsolrot
      use mathfun, only:akang
c      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 , intent(in)::np,mfring,ibsi
      integer*4 i
      real*8 , intent(in)::theta2,al,bz,ak0,dx,dy,theta,radlvl,
     $     f1in,f2in,f1out,f2out,eps0
      real*8 , intent(inout)::
     $     x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np)
      real*8 alr,bxs,bys,bzs,a,b,ea,p,pxf,pyf,ak
      logical*4 , intent(in)::fringe
      logical*4 krad
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       4,0.d0,ak0,
     $       dx,dy,theta,.false.,.false.)
        return
      endif
      krad=rad .and. radlvl .ne. 1.d0
c      theta2=theta+akang(dcmplx(ak0,0.d0),al,cr1)
      ak=sign(ak0,al)
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.true.)
      if(calpol)then
        bsi=0.d0
      endif
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
      endif
      if(fringe .and. mfring .ne. 2)then
        call ttfrin(np,x,px,y,py,z,g,4,ak,al,bz)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        do 2110 i=1,np
c          p=(1.d0+g(i))**2
          p=(1.d0+g(i))
          a=f1in/p
          ea=exp(a)
          b=f2in/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
2110    continue
      endif
      if(krad)then
        if(photons)then
          call tsetpcvt(l_track,dx,dy,theta2,0.d0,0.d0,al)
          pcvt%fr0=-0.5d0*f1in/al
        endif
        if(f1in .ne. 0.d0)then
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,f1in,0.d0)
        endif
        call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,ak,bz,0.d0,0.d0,eps0,
     $       alr)
      else
        call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,ak,bz,0.d0,0.d0,ibsi,eps0)
      endif
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        do 2120 i=1,np
c          p=(1.d0+g(i))**2
          p=(1.d0+g(i))
          a=-f1out/p
          ea=exp(a)
          b=f2out/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
2120    continue
      endif
      if(fringe .and. mfring .ne. 1)then
        call ttfrin(np,x,px,y,py,z,g,4,-ak,al,bz)
      endif
      if(krad .and. f1out .ne. 0.d0)then
        pcvt%fr0=1.d0-.5d0*f1out/al
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,f1out,0.d0)
      endif
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.false.)
      return
      end
