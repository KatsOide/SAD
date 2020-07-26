      subroutine tcoorde(trans,cod,beam,
     1                 dx,dy,dz,chi1,chi2,chi3,dir)
      use ffs_flag
      use tmacro
      use mathfun, only: sqrtl
      use temw,only:tmulbs
      implicit none
      integer*4 i
      real*8 dx,dy,dz,chi1,chi2,chi3,cchi1,schi1,
     $     cchi2,schi2,cchi3,schi3,dx0,dy0,dz0,dx1,dy1,dz1,
     $     r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $     pxi,pyi,a,pzi,zi,yi,xi,xf,yf,zf,
     $     pzf,p,z,pxf,pyf,u,v,x,y
      real*8 trans(6,12),cod(6),beam(42)
      logical*4 dir
      real*8 trans1(6,6),trans2(6,6)
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
      p=1.d0+cod(6)
      z=cod(5)
      pxi=cod(2)
      pyi=cod(4)
      a=pxi**2+pyi**2
      pzi=p*sqrtl(1.d0-a/p**2)
      zi=pzi/p*z-dz0
      xi=cod(1)-dx0+pxi/p*z
      yi=cod(3)-dy0+pyi/p*z
      call tinitr(trans1)
      trans1(1,2)= z/p
      trans1(1,5)= pxi/p
      trans1(1,6)=-pxi*z/p**2
      trans1(3,4)= z/p
      trans1(3,5)= pyi/p
      trans1(3,6)=-pyi*z/p**2
      trans1(5,2)=-pxi/p/pzi*z
      trans1(5,4)=-pyi/p/pzi*z
      trans1(5,5)=pzi/p
      trans1(5,6)=(pxi**2+pyi**2)/p**2/pzi*z
      trans1(6,2)=-pxi/pzi
      trans1(6,4)=-pyi/pzi
      trans1(6,6)=p/pzi
      xf =r11*xi +r12*yi +r13*zi
      yf =r21*xi +r22*yi +r23*zi
      zf =r31*xi +r32*yi +r33*zi+dz1
      pxf=r11*pxi+r12*pyi+r13*pzi
      pyf=r21*pxi+r22*pyi+r23*pzi
      pzf=r31*pxi+r32*pyi+r33*pzi
      do 110 i=1,6
        x=trans1(1,i)
        y=trans1(3,i)
        u=trans1(2,i)
        v=trans1(4,i)
        trans1(1,i)
     1     =r11*x+r12*y+r13*trans1(5,i)
        trans1(3,i)
     1     =r21*x+r22*y+r23*trans1(5,i)
        trans1(5,i)
     1     =r31*x+r32*y+r33*trans1(5,i)
        trans1(2,i)
     1     =r11*u+r12*v+r13*trans1(6,i)
        trans1(4,i)
     1     =r21*u+r22*v+r23*trans1(6,i)
        trans1(6,i)
     1     =r31*u+r32*v+r33*trans1(6,i)
110   continue
      call tinitr(trans2)
      trans2(1,2)=-zf/pzf
      trans2(1,5)=-pxf/pzf
      trans2(1,6)= pxf/pzf**2*zf
      trans2(3,4)=-zf/pzf
      trans2(3,5)=-pyf/pzf
      trans2(3,6)= pyf/pzf**2*zf
      trans2(5,2)= pxf/p/pzf*zf
      trans2(5,4)= pyf/p/pzf*zf
      trans2(5,5)= p/pzf
      trans2(5,6)=-(pxf**2+pyf**2)/p/pzf**2*zf
      trans2(6,2)=pxf/p
      trans2(6,4)=pyf/p
      trans2(6,6)=pzf/p
      call tmultr(trans1,trans2,6)
      call tmultr(trans,trans1,irad)
      call tmulbs(beam ,trans1,.true.)
      cod(1)=xf-pxf/pzf*zf+dx1
      cod(3)=yf-pyf/pzf*zf+dy1
      cod(5)=p/pzf*zf
      cod(2)=pxf
      cod(4)=pyf
      return
      end
