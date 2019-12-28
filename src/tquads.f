      subroutine tquads(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak0,bz,
     1     dx,dy,theta,radlvl,fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0)
      use tfstk
      use ffs_flag
      use tmacro
      use photontable,only:tsetphotongeo
      use tspin
      use sol,only:tsolrot
      use mathfun, only:akang
c      use ffs_pointer, only:inext,iprev
      implicit none
      integer*4 np,mfring,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     pxr0(np),pyr0(np),zr0(np),bsi(np),
     $     al,bz,ak0,ak,dx,dy,theta,radlvl,alr,
     $     f1in,f2in,f1out,f2out,eps0,theta1,theta2,bxs,bys,bzs,
     $     a,b,ea,p,pxf,pyf
      real*8 sx(np),sy(np),sz(np)
      complex*16 cr1
      logical*4 enarad,fringe
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       4,0.d0,ak0,
     $       dx,dy,theta, 1.d0,.false.)
        return
      endif
      enarad=rad .and. radlvl .ne. 1.d0
      call akang(dcmplx(ak0,0.d0),al,theta1,cr1)
      if(theta1 .ne. 0.d0)then
        ak=-ak0
      else
        ak=ak0
      endif
      theta2=theta+theta1
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.true.)
      if(enarad)then
        pxr0=px
        pyr0=py
        zr0=z
        bsi=0.d0
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
      if(enarad)then
        if(photons)then
          call tsetphotongeo(0.d0,0.d0,theta2,.true.)
        endif
        if(f1in .ne. 0.d0)then
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         pxr0,pyr0,zr0,bsi,f1in,0.d0)
        endif
        call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       pxr0,pyr0,zr0,bsi,
     $       al,ak,bz,0.d0,0.d0,eps0,
     $       alr)
      else
        call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       bsi,al,ak,bz,0.d0,0.d0,0,eps0)
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
      if(enarad .and. f1out .ne. 0.d0)then
        if(photons)then
          call tsetphotongeo(0.d0,0.d0,0.d0,.false.)
        endif
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       pxr0,pyr0,zr0,bsi,f1out,0.d0)
      endif
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.false.)
      return
      end
