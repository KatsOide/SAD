      subroutine tsconv(trans1,cod,l1,ent)
      use tfstk
      use kyparam
      use tmacro
      implicit none
      integer*8 l1
      real*8 trans1(6,6),cod(6)
      real*8 ds
      logical*4 ent
c     begin initialize for preventing compiler warning
      ds=0.d0
c     end   initialize for preventing compiler warning
      if(ent)then
        call tsrote(trans1,cod,rlist(l1+ky_CHI1_SOL),
     $       rlist(l1+ky_CHI2_SOL),rlist(l1+ky_CHI3_SOL))
        cod(1)=cod(1)+rlist(l1+ky_DX_SOL)
        cod(3)=cod(3)+rlist(l1+ky_DY_SOL)
        cod(5)=cod(5)+rlist(l1+ky_DZ_SOL)
      else
        ds=rlist(l1+ky_DZ_SOL)
        cod(1)=cod(1)-rlist(l1+ky_DX_SOL)
        cod(3)=cod(3)-rlist(l1+ky_DY_SOL)
        call tsrote(trans1,cod,rlist(l1+ky_CHI1_SOL),
     $       rlist(l1+ky_CHI2_SOL),rlist(l1+ky_CHI3_SOL))
        cod(5)=cod(5)+ds-dvemit*ds
      endif
      return
      end

      subroutine tsrote(trans,cod,chi1,chi2,chi3)
      use tfstk, only:sqrt1
      implicit none
      real*8 trans(6,6),cod(6),chi1,chi2,chi3,
     $     r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $     pxi,pyi,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf,
     $     cchi1,schi1,cchi2,schi2,cchi3,schi3,
     $     dpzfdpx,dpzfdpy,dpzfdp,p
      cchi1=cos(chi1)
      schi1=sin(chi1)
      cchi2=cos(chi2)
      schi2=sin(chi2)
      cchi3=cos(chi3)
      schi3=sin(chi3)
      r11= cchi1*cchi3+schi1*schi2*schi3
      r12=-cchi2*schi3
      r13= schi1*cchi3-cchi1*schi2*schi3
      r21=-schi1*schi2*cchi3+cchi1*schi3
      r22= cchi2*cchi3
      r23= cchi1*schi2*cchi3+schi1*schi3
      r31=-schi1*cchi2
      r32=-schi2
      r33= cchi1*cchi2
      xi=cod(1)
      yi=cod(3)
      p=1.d0+cod(6)
      pxi=cod(2)
      pyi=cod(4)
      pzi=p*(1.d0+sqrt1(-(pxi**2+pyi**2)/p**2))
      xf =r11*xi +r12*yi 
      yf =r21*xi +r22*yi 
      zf =r31*xi +r32*yi 
      pxf=r11*pxi+r12*pyi+r13*pzi
      pyf=r21*pxi+r22*pyi+r23*pzi
      pzf=r31*pxi+r32*pyi+r33*pzi
      dpzfdpx=r31-r33*pxi/pzi
      dpzfdpy=r32-r33*pyi/pzi
      dpzfdp =    r33*p  /pzi
      trans(2,1)=0.d0
      trans(2,2)=r11-r13*pxi/pzi
      trans(2,3)=0.d0
      trans(2,4)=r12-r13*pyi/pzi
      trans(2,5)=0.d0
      trans(2,6)=    r13*p  /pzi
      trans(4,1)=0.d0
      trans(4,2)=r21-r23*pxi/pzi
      trans(4,3)=0.d0
      trans(4,4)=r22-r23*pyi/pzi
      trans(4,5)=0.d0
      trans(4,6)=    r23*p  /pzi
      trans(1,1)=r11-pxf/pzf*r31
      trans(1,2)=-zf/pzf*(trans(2,2)-pxf/pzf*dpzfdpx)
      trans(1,3)=r12-pxf/pzf*r32
      trans(1,4)=-zf/pzf*(trans(2,4)-pxf/pzf*dpzfdpy)
      trans(1,5)=0.d0
      trans(1,6)=-zf/pzf*(trans(2,6)-pxf/pzf*dpzfdp )
      trans(3,1)=r21-pyf/pzf*r31
      trans(3,2)=-zf/pzf*(trans(4,2)-pyf/pzf*dpzfdpx)
      trans(3,3)=r22-pyf/pzf*r32
      trans(3,4)=-zf/pzf*(trans(4,4)-pyf/pzf*dpzfdpy)
      trans(3,5)=0.d0
      trans(3,6)=-zf/pzf*(trans(4,6)-pyf/pzf*dpzfdp )
      trans(5,1)=    p  /pzf*r31
      trans(5,2)= zf/pzf*(           -p  /pzf*dpzfdpx)
      trans(5,3)=    p  /pzf*r32
      trans(5,4)= zf/pzf*(           -p  /pzf*dpzfdpy)
      trans(5,5)=1.d0
      trans(5,6)= zf/pzf*(1.d0       -p  /pzf*dpzfdp )
      trans(6,1:5)=0.d0
      trans(6,6)=1.d0
      cod(1)=xf-pxf/pzf*zf
      cod(3)=yf-pyf/pzf*zf
      cod(5)=cod(5)+p/pzf*zf
      cod(2)=pxf
      cod(4)=pyf
      return
      end
