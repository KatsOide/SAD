      subroutine tcoord(np,x,px,y,py,z,
     1                 dx,dy,dz,chi1,chi2,chi3,dir)
      use tmacro
      implicit real*8(a-h,o-z)
      logical dir
      dimension x(np),px(np),y(np),py(np),z(np)
c
c   xi        c3 -s3 0      1  0  0      c1 0  s1      x-dx
c   eta  =    s3  c3 0      0  c2 s2     0  1  0       y-dy
c   zeta      0   0  1      0 -s2 c2    -s1 0  c1      z-dz
c
      cchi1=cos(chi1)
      schi1=sin(chi1)
      cchi2=cos(chi2)
      schi2=sin(chi2)
      cchi3=cos(chi3)
      schi3=sin(chi3)
      if(dir)then
        dx0=dx
        dy0=dy
        dz0=dz
        dx1=0.d0
        dy1=0.d0
        dz1=0.d0
        r11= cchi1*cchi3+schi1*schi2*schi3
        r12=-cchi2*schi3
        r13= schi1*cchi3-cchi1*schi2*schi3
        r21=-schi1*schi2*cchi3+cchi1*schi3
        r22= cchi2*cchi3
        r23= cchi1*schi2*cchi3+schi1*schi3
        r31=-schi1*cchi2
        r32=-schi2
        r33= cchi1*cchi2
      else
        dx0=0.d0
        dy0=0.d0
        dz0=0.d0
        dx1=dx
        dy1=-dy
        dz1=-dz
        r11= cchi1*cchi3+schi1*schi2*schi3
        r21= cchi2*schi3
        r31=-schi1*cchi3+cchi1*schi2*schi3
        r12= schi1*schi2*cchi3-cchi1*schi3
        r22= cchi2*cchi3
        r32= cchi1*schi2*cchi3+schi1*schi3
        r13= schi1*cchi2
        r23=-schi2
        r33= cchi1*cchi2
      endif
      do 10 i=1,np
        pxi=px(i)
        pyi=py(i)
        a=min(.95d0,pxi**2+pyi**2)
        dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        pzi=1.d0+dpz
        zi=pzi*z(i)-dz0
        xi=x(i)+pxi*z(i)-dx0
        yi=y(i)+pyi*z(i)-dy0
        xf =r11*xi +r12*yi +r13*zi
        yf =r21*xi +r22*yi +r23*zi
        zf =r31*xi +r32*yi +r33*zi+dz1
        px(i)=r11*pxi+r12*pyi+r13*pzi
        py(i)=r21*pxi+r22*pyi+r23*pzi
        pzf  =r31*pxi+r32*pyi+r33*pzi
        z(i)=zf/pzf
        x(i)=xf-px(i)*z(i)+dx1
        y(i)=yf-py(i)*z(i)+dy1
10    continue
      return
      end
