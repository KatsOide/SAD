      subroutine tquads(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,bz,
     1     dx,dy,theta,cost,sint,radlvl,
     1     fringe,f1in,f2in,f1out,f2out,
     $     mfring,eps0,ld,forward)
      use tfstk
      use ffs_flag
      use tmacro
      use photontable,only:tsetphotongeo
      use tspin
c      use ffs_pointer, only:inext,iprev
      implicit none
      type (sad_rlist), pointer :: klr
      integer*4 np,ld,mfring,i,irtc,ld1,level,m,itfuplevel,
     $     itfdownlevel
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     pxr0(np),pyr0(np),zr0(np),bsi(np),
     $     al,bz,ak,dx,dy,theta,cost,sint,radlvl,alr,
     $     f1in,f2in,f1out,f2out,eps0,
     $     a,aki,akm,ali,alm,b,ea,fx,fy,p,pr,px0,pxf,pyf,rb,x0
      real*8 sx(np),sy(np),sz(np)
      logical*4 enarad,fringe,forward
      character*13 vname
      character*2 ord
      integer*8 ifv,ifvh,kx
      save ifv,ifvh
      data ifv/0/
      data vname/'SolenoidShape'/
      if(ifv .eq. 0)then
        ifv=ktavaloc(0,1)
        ifvh=ktfsymbolz(vname,len(vname))
        klist(ifv)=ktfsymbol+ktfcopy1(ifvh)
      endif
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       4,0.d0,ak,
     $       dx,dy,theta,cost,sint, 1.d0,.false.)
        return
      endif
      enarad=rad .and. radlvl .ne. 1.d0
      fx= bz*dy*.5d0
      fy=-bz*dx*.5d0
      do i=1,np
c        pr=(1.d0+g(i))**2
        pr=(1.d0+g(i))
        x(i)=x(i)-dx
        y(i)=y(i)-dy
        px(i)=px(i)+fx/pr
        py(i)=py(i)+fy/pr
      enddo
      if(theta .ne. 0.d0)then
        do i=1,np
          x0=x(i)
          x(i)=cost*x0-sint*y(i)
          y(i)=sint*x0+cost*y(i)
          px0=px(i)
          px(i)=cost*px0-sint*py(i)
          py(i)=sint*px0+cost*py(i)
        enddo
      endif
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
      if(ifv .eq. 0)then
        if(enarad)then
          if(photons)then
            call tsetphotongeo(0.d0,0.d0,theta,.true.)
          endif
          if(f1in .ne. 0.d0)then
            call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           pxr0,pyr0,zr0,1.d0,0.d0,bsi,f1in)
          endif
          pxr0=px
          pyr0=py
          zr0=z
          bsi=0.d0
          call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         pxr0,pyr0,zr0,bsi,
     $         al,ak,bz,0.d0,0.d0,eps0,
     $         alr)
        else
          call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         bsi,al,ak,bz,0.d0,0.d0,0,eps0)
        endif
      else
        level=itfuplevel()
        ld1=ld
 1      rlist(ifv+1)=ld1
        call tfleval(klist(ifv-3),kx,.true.,irtc)
        if(irtc .ne. 0)then
          level=itfdownlevel()
          if(ierrorprint .ne. 0)then
            call tfaddmessage(' ',2,6)
          endif
          write(*,*)' Error in '//vname//' at ',ld1,ord(ld1),
     $         ' element.'
          return
        elseif(tfreallistq(kx,klr))then
          m=klr%nl
          if(m .le. 2)then
            write(*,*)' '//vname//' must have more than 2 numbers at ',
     $           ld1,ord(ld1)//' element.'
            return
          endif
          alm=al/(m-1)
          akm=ak/(m-1)
          if(enarad .and. photons)then
            call tsetphotongeo(0.d0,0.d0,theta,.true.)
          endif
          do i=1,m
            if(forward)then
              rb=klr%rbody(i)
            else
              rb=klr%rbody(m-i+1)
            endif
            if(i .eq. 1 .or. i .eq. m)then
              ali=alm*.5d0
              aki=akm*.5d0
            else
              ali=alm
              aki=akm
            endif
            if(enarad)then
              call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             pxr0,pyr0,zr0,bsi,ali,aki,
     $             bz*rb,0.d0,0.d0,eps0,alr)
            else
              call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             bsi,ali,aki,
     $             bz*rb,0.d0,0.d0,0,eps0)
            endif
          enddo
          level=itfdownlevel()
        else
          if(ld .eq. ld1)then
            ld1=ld+1
            go to 1
          else
            level=itfdownlevel()
            if(enarad)then
              call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             pxr0,pyr0,zr0,bsi,al,ak,
     $             bz,0.d0,0.d0,eps0,alr)
            else
              call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             bsi,al,ak,
     $             bz,0.d0,0.d0,0,eps0)
            endif
          endif
        endif
      endif
c      if(enarad)then
c        call trad(np,x,px,y,py,g,dv,0.d0,
c     1       b1,0.d0,0.d0,.5d0*al,f1r,f2r,al,al,-1.d0)
c      endif
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
     $       pxr0,pyr0,zr0,1.d0,0.d0,bsi,f1out)
      endif
      if(theta .ne. 0.d0)then
        do i=1,np
          x0=x(i)
          x(i)= cost*x0+sint*y(i)
          y(i)=-sint*x0+cost*y(i)
          px0=px(i)
          px(i)= cost*px0+sint*py(i)
          py(i)=-sint*px0+cost*py(i)
        enddo
      endif
      do i=1,np
c        pr=(1.d0+g(i))**2
        pr=(1.d0+g(i))
        px(i)=px(i)-fx/pr
        py(i)=py(i)-fy/pr
        x(i)=x(i)+dx
        y(i)=y(i)+dy
      enddo
      return
      end
