      subroutine tsconv(trans1,cod,l1,ent)
      use tfstk
      use tmacro
      implicit none
      integer*4 l1
      real*8 trans1(6,6),cod(6)
      real*8 p,pxi,pyi,pzi,xi,yi,ds,
     $     cchi1,schi1,cchi2,schi2,cchi3,schi3,
     $     r11,r12,r13,r21,r22,r23,r31,r32,r33,xf,yf,zf,pxf,pyf,pzf,
     $     dpzfdpx,dpzfdpy,dpzfdp
      logical*4 ent
c     begin initialize for preventing compiler warning
      ds=0.d0
c     end   initialize for preventing compiler warning
      p=1.d0+cod(6)
      if(ent)then
        pxi=cod(2)
        pyi=cod(4)
        pzi=p*(1.d0+sqrt1(-(pxi**2+pyi**2)/p**2))
c        pzi=p*sqrt(1.d0-(pxi**2+pyi**2)/p**2)
c        pzi=p-(pxi**2+pyi**2)/(p+sqrt((p-pyi)*(p+pyi)-pxi**2))
        xi=cod(1)
        yi=cod(3)
      else
        ds=rlist(l1+5)
        pxi=cod(2)
        pyi=cod(4)
        pzi=p*(1.d0+sqrt1(-(pxi**2+pyi**2)/p**2))
c        pzi=p*sqrt(1.d0-(pxi**2+pyi**2)/p**2)
c        pzi=p-(pxi**2+pyi**2)/(p+sqrt((p-pyi)*(p+pyi)-pxi**2))
        xi=cod(1)-rlist(l1+3)
        yi=cod(3)-rlist(l1+4)
      endif
      call tinitr(trans1)
      cchi1=cos(rlist(l1+9 ))
      schi1=sin(rlist(l1+9 ))
      cchi2=cos(rlist(l1+10))
      schi2=sin(rlist(l1+10))
      cchi3=cos(rlist(l1+11))
      schi3=sin(rlist(l1+11))
      r11= cchi1*cchi3+schi1*schi2*schi3
      r12=-cchi2*schi3
      r13= schi1*cchi3-cchi1*schi2*schi3
      r21=-schi1*schi2*cchi3+cchi1*schi3
      r22= cchi2*cchi3
      r23= cchi1*schi2*cchi3+schi1*schi3
      r31=-schi1*cchi2
      r32=-schi2
      r33= cchi1*cchi2
      xf =r11*xi +r12*yi 
      yf =r21*xi +r22*yi 
      zf =r31*xi +r32*yi 
      pxf=r11*pxi+r12*pyi+r13*pzi
      pyf=r21*pxi+r22*pyi+r23*pzi
      pzf=r31*pxi+r32*pyi+r33*pzi
      dpzfdpx=r31-r33*pxi/pzi
      dpzfdpy=r32-r33*pyi/pzi
      dpzfdp =    r33*p  /pzi
      trans1(2,2)=r11-r13*pxi/pzi
      trans1(2,4)=r12-r13*pyi/pzi
      trans1(2,6)=    r13*p  /pzi
      trans1(4,2)=r21-r23*pxi/pzi
      trans1(4,4)=r22-r23*pyi/pzi
      trans1(4,6)=    r23*p  /pzi
      trans1(1,1)=r11-pxf/pzf*r31
      trans1(1,2)=-zf/pzf*(trans1(2,2)-pxf/pzf*dpzfdpx)
      trans1(1,3)=r12-pxf/pzf*r32
      trans1(1,4)=-zf/pzf*(trans1(2,4)-pxf/pzf*dpzfdpy)
      trans1(1,6)=-zf/pzf*(trans1(2,6)-pxf/pzf*dpzfdp )
      trans1(3,1)=r21-pyf/pzf*r31
      trans1(3,2)=-zf/pzf*(trans1(4,2)-pyf/pzf*dpzfdpx)
      trans1(3,3)=r22-pyf/pzf*r32
      trans1(3,4)=-zf/pzf*(trans1(4,4)-pyf/pzf*dpzfdpy)
      trans1(3,6)=-zf/pzf*(trans1(4,6)-pyf/pzf*dpzfdp )
      trans1(5,1)=    p  /pzf*r31
      trans1(5,2)= zf/pzf*(           -p  /pzf*dpzfdpx)
      trans1(5,3)=    p  /pzf*r32
      trans1(5,4)= zf/pzf*(           -p  /pzf*dpzfdpy)
      trans1(5,6)= zf/pzf*(1.d0       -p  /pzf*dpzfdp )
      if(ent)then
        cod(1)=xf-pxf/pzf*zf+rlist(l1+3)
        cod(3)=yf-pyf/pzf*zf+rlist(l1+4)
        cod(5)=cod(5)+p/pzf*zf+rlist(l1+5)
      else
        cod(1)=xf-pxf/pzf*zf
        cod(3)=yf-pyf/pzf*zf
        cod(5)=cod(5)+p/pzf*zf+ds-dvemit*ds
      endif
      cod(2)=pxf
      cod(4)=pyf
      return
      end
