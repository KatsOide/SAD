      subroutine tbendi(np,x,px,y,py,z,g,dv,pz,al,phib,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dphix,dphiy,cost,sint,
     1     fb1,fb2,mfring,
     1     enarad,fringe,eps0)
      include 'inc/TMACRO.inc'
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      dimension x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np)
      logical*4 enarad,fringe
      include 'inc/TENT.inc'
      if(dphiy .ne. 0.d0)then
        do 3510 i=1,np
c          pr=(1.d0+g(i))**2
          pr=(1.d0+g(i))
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
3510    continue
      endif
      tanp1=sinp1/cosp1
      tanp2=sinp2/cosp2
      rhob=al/phib
      rho0=al/phi0
      aind=rho0/phi0*ak
      if(rad .and. enarad)then
        b=brhoz/rhob
        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
     1             1.d0/rho0,-tanp1*2.d0/al,.5d0*al,
     $       fb1,fb2,0.d0,al,1.d0)
      endif
      if(fb1 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -1)then
          dxfr1=fb1**2/rhob/24.d0
          dyfr1=fb1/rhob**2/6.d0
          dyfra1=4.d0*dyfr1/fb1**2
          do i=1,np
c            dp=g(i)*(2.d0+g(i))
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)+dxfr1*dp/pr
            py(i)=py(i)+(dyfr1-dyfra1*y(i)**2)*y(i)/pr**2
            z(i)=z(i)+(dxfr1*px(i)+
     $           (.5d0*dyfr1-.25d0*dyfra1*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
      if(eps0 .le. 0.d0)then
        eps=1.d-5
      else
        eps=1.d-5*eps0
      endif
      ndiv=1+int(sqrt(abs(ak*al)/eps/12.d0))
      if(fringe)then
        af=1.d0
      else
        af=0.d0
      endif
      do 1100 i=1,np
c        p=(1.d0+g(i))**2
        p=(1.d0+g(i))
        rhoe=rhob*p
        f=y(i)/rhoe
        fpx=af*px(i)
        py(i)=py(i)-(tanp1+fpx)*f
        ff  =af*y(i)*f*.5d0
        x(i)=x(i)+ff
        z(i)=z(i)-ff*fpx
        px(i)=px(i)+tanp1*x(i)/rhoe
1100  continue
      phin=phi0/ndiv
      akn=ak/ndiv
      aln=al/ndiv
      csphi0=cos(phin)
      snphi0=sin(phin)
      if(csphi0 .ge. 0.d0)then
        sinsq0=snphi0**2/(1.d0+csphi0)
      else
        sinsq0=1.d0-csphi0
      endif
c      sinsq0=2.d0*sin(phin*.5d0)**2
      drhob=rhob-rho0
      do 1120 n=1,ndiv
        if(n .eq. 1)then
          ak1=akn*.5d0
        else
          ak1=akn
        endif
        do 1130 i=1,np
c          dp=g(i)*(2.d0+g(i))
          dp=g(i)
          p=1.d0+dp
          xi=x(i)
          xr=xi/rho0
          pxi=px(i)-ak1*(xi+xi*xr*(.5d0-xr*(2.d0-xr)/12.d0))/p
          yi=y(i)
          pyi=py(i)+ak1*yi/p
          rhoe=rhob*p
          s=min(.95d0,pxi**2+pyi**2)
          dpz1=s*(-.5d0-s*(.125d0+s*(.0625d0+s*.0390625d0)))
          if(s .gt. 1.d-4)then
            dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
            dpz1=(dpz1**2-s)/(2.d0+2.d0*dpz1)
          endif
          pz1=1.d0+dpz1
          drho=drhob+rhoe*dpz1+rhob*dp
          dpx=-(xi-drho)/rhoe*snphi0-sinsq0*pxi
          pxf=pxi+dpx
          s=min(.95d0,pxf**2+pyi**2)
          dpz2=s*(-.5d0-s*(.125d0+s*(.0625d0+s*.0390625d0)))
          if(s .gt. 1.d-4)then
            dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
            dpz2=(dpz2**2-s)/(2.d0+2.d0*dpz2)
          endif
          pz2=1.d0+dpz2
          d=pxf*pz1+pxi*pz2
          if(d .eq. 0.d0)then
            sinda=2.d0*pxf*pz2/(pxf**2+pz2**2)
          else
            sinda=dpx*(pxf+pxi)/d
          endif
          s=sinda**2
          if(s .gt. 2.d-4)then
            da=sinda
     1     *(1.d0+s*(a3+s*(a5+s*(a7+s*(a9+s*(a11+s*(a13+s*a15)))))))
          else
            da=sinda*(1.d0+s*(a3+s*(a5+a7*s)))
          endif
          x(i)=xi*csphi0+rhoe*(snphi0*pxi-dpx*(pxi+pxf)/(pz1+pz2))
     1         +drho*sinsq0
          px(i)=pxf
          y(i)=yi+pyi*rhoe*(phin-da)
          py(i)=pyi
          z(i)=z(i)-phin*(dp*rhob+drhob)+da*rhoe-dv(i)*aln
1130    continue
1120  continue
      do 1140 i=1,np
c        p=(1.d0+g(i))**2
        p=(1.d0+g(i))
        xr=x(i)/rho0
        px(i)=px(i)-
     1        akn*(x(i)+x(i)*xr*(.5d0-xr*(2.d0-xr)/12.d0))/p*.5d0
        py(i)=py(i)+akn*y(i)/p*.5d0
1140  continue
c     if(np .gt. 16)then
        do 1150 i=1,np
c          p=(1.d0+g(i))**2
          p=(1.d0+g(i))
          rhoe=rhob*p
          px(i)=px(i)+tanp2*x(i)/rhoe
          f=y(i)/rhoe
          fpx=af*px(i)
          py(i)=py(i)-(tanp2-fpx)*f
          ff  =af*y(i)*f*.5d0
          x(i)=x(i)-ff
          z(i)=z(i)+ff*fpx
1150    continue
c     else
c*voption novec
c       do 1151 i=1,np
c1151    continue
c     endif
      if(fb2 .ne. 0.d0)then
        if(mfring .gt. 0 .or. mfring .eq. -2)then
          dxfr2=fb2**2/rhob/24.d0
          dyfr2=fb2/rhob**2/6.d0
          dyfra2=4.d0*dyfr2/fb2**2
          do i=1,np
            dp=g(i)
            pr=1.d0+dp
            x(i)=x(i)-dxfr2*dp/pr
            py(i)=py(i)+(dyfr2-dyfra2*y(i)**2)*y(i)/pr**2
            z(i)=z(i)-(dxfr2*px(i)-
     $           (.5d0*dyfr2-.25d0*dyfra2*y(i)**2)*y(i)**2/pr)/pr
          enddo
        endif
      endif
      if(rad .and. enarad)then
        call trad(np,x,px,y,py,g,dv,b,0.d0,b*(aind/rho0),
     1             1.d0/rho0,-tanp2*2.d0/al,.5d0*al,
     $       fb1,fb2,al,al,-1.d0)
      endif
      if(dphiy .ne. 0.d0)then
        do 3520 i=1,np
c          pr=(1.d0+g(i))**2
          pr=(1.d0+g(i))
          px(i)=px(i)+dphix/pr
          py(i)=py(i)+dphiy/pr
3520    continue
      endif
      include 'inc/TEXIT.inc'
      return
      end
