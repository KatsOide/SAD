      subroutine tsconv(trans1,cod,rr,l1,ent)
      use tfstk
      use kyparam
      use tmacro
      implicit none
      integer*8 ,intent(in):: l1
      real*8, intent(inout):: cod(6)
      real*8 , intent(out)::trans1(6,6),rr(3,3)
      real*8 ds
      logical*4 ,intent(in):: ent
c     begin initialize for preventing compiler warning
      ds=0.d0
c     end   initialize for preventing compiler warning
      if(ent)then
        call tsrote(trans1,cod,rr,rlist(l1+ky_CHI1_SOL),
     $       rlist(l1+ky_CHI2_SOL),rlist(l1+ky_CHI3_SOL))
        cod(1)=cod(1)+rlist(l1+ky_DX_SOL)
        cod(3)=cod(3)+rlist(l1+ky_DY_SOL)
        cod(5)=cod(5)+rlist(l1+ky_DZ_SOL)
      else
        ds=rlist(l1+ky_DZ_SOL)
        cod(1)=cod(1)-rlist(l1+ky_DX_SOL)
        cod(3)=cod(3)-rlist(l1+ky_DY_SOL)
        call tsrote(trans1,cod,rr,rlist(l1+ky_CHI1_SOL),
     $       rlist(l1+ky_CHI2_SOL),rlist(l1+ky_CHI3_SOL))
        cod(5)=cod(5)+ds-dvemit*ds
      endif
      return
      end

      subroutine tsrote(trans,cod,rr,chi1,chi2,chi3)
      use mathfun, only:sqrt1
      implicit none
      real*8 ,intent(inout):: cod(6)
      real*8 ,intent(in):: chi1,chi2,chi3
      real*8 ,intent(out):: trans(6,6),rr(3,3)
      real*8 pxi,pyi,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf,
     $     cchi1,schi1,cchi2,schi2,cchi3,schi3,
     $     dpzfdpx,dpzfdpy,dpzfdp,p
      cchi1=cos(chi1)
      schi1=sin(chi1)
      cchi2=cos(chi2)
      schi2=sin(chi2)
      cchi3=cos(chi3)
      schi3=sin(chi3)
      rr(1,1)= cchi1*cchi3+schi1*schi2*schi3
      rr(1,2)=-cchi2*schi3
      rr(1,3)= schi1*cchi3-cchi1*schi2*schi3
      rr(2,1)=-schi1*schi2*cchi3+cchi1*schi3
      rr(2,2)= cchi2*cchi3
      rr(2,3)= cchi1*schi2*cchi3+schi1*schi3
      rr(3,1)=-schi1*cchi2
      rr(3,2)=-schi2
      rr(3,3)= cchi1*cchi2
      xi=cod(1)
      yi=cod(3)
      p=1.d0+cod(6)
      pxi=cod(2)
      pyi=cod(4)
      pzi=p*(1.d0+sqrt1(-(pxi**2+pyi**2)/p**2))
      xf =rr(1,1)*xi +rr(1,2)*yi 
      yf =rr(2,1)*xi +rr(2,2)*yi 
      zf =rr(3,1)*xi +rr(3,2)*yi 
      pxf=rr(1,1)*pxi+rr(1,2)*pyi+rr(1,3)*pzi
      pyf=rr(2,1)*pxi+rr(2,2)*pyi+rr(2,3)*pzi
      pzf=rr(3,1)*pxi+rr(3,2)*pyi+rr(3,3)*pzi
      dpzfdpx=rr(3,1)-rr(3,3)*pxi/pzi
      dpzfdpy=rr(3,2)-rr(3,3)*pyi/pzi
      dpzfdp =    rr(3,3)*p  /pzi
      trans(2,1)=0.d0
      trans(2,2)=rr(1,1)-rr(1,3)*pxi/pzi
      trans(2,3)=0.d0
      trans(2,4)=rr(1,2)-rr(1,3)*pyi/pzi
      trans(2,5)=0.d0
      trans(2,6)=    rr(1,3)*p  /pzi
      trans(4,1)=0.d0
      trans(4,2)=rr(2,1)-rr(2,3)*pxi/pzi
      trans(4,3)=0.d0
      trans(4,4)=rr(2,2)-rr(2,3)*pyi/pzi
      trans(4,5)=0.d0
      trans(4,6)=    rr(2,3)*p  /pzi
      trans(1,1)=rr(1,1)-pxf/pzf*rr(3,1)
      trans(1,2)=-zf/pzf*(trans(2,2)-pxf/pzf*dpzfdpx)
      trans(1,3)=rr(1,2)-pxf/pzf*rr(3,2)
      trans(1,4)=-zf/pzf*(trans(2,4)-pxf/pzf*dpzfdpy)
      trans(1,5)=0.d0
      trans(1,6)=-zf/pzf*(trans(2,6)-pxf/pzf*dpzfdp )
      trans(3,1)=rr(2,1)-pyf/pzf*rr(3,1)
      trans(3,2)=-zf/pzf*(trans(4,2)-pyf/pzf*dpzfdpx)
      trans(3,3)=rr(2,2)-pyf/pzf*rr(3,2)
      trans(3,4)=-zf/pzf*(trans(4,4)-pyf/pzf*dpzfdpy)
      trans(3,5)=0.d0
      trans(3,6)=-zf/pzf*(trans(4,6)-pyf/pzf*dpzfdp )
      trans(5,1)=    p  /pzf*rr(3,1)
      trans(5,2)= zf/pzf*(           -p  /pzf*dpzfdpx)
      trans(5,3)=    p  /pzf*rr(3,2)
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
