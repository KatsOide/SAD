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
      use chg, only:tbgrote
      use mathfun, only:xsincos
      implicit none
      real*8 ,intent(inout):: cod(6)
      real*8 ,intent(in):: chi1,chi2,chi3
      real*8 ,intent(out):: trans(6,6),rr(3,3)
      real*8 s1,xs1,c1,dc1,s2,xs2,c2,dc2,s3,xs3,c3,dc3,dr(3,3)
      call xsincos(chi1,s1,xs1,c1,dc1)
      call xsincos(chi2,s2,xs2,c2,dc2)
      call xsincos(chi3,s3,xs3,c3,dc3)
      dr(1,1)= dc1+dc3+dc1*dc3+s1*s2*s3
      dr(1,2)=-c2*s3
      dr(1,3)= s1*c3          -c1*s2*s3
      dr(2,1)=-s1*s2*c3       +c1*s3
      dr(2,2)= dc2+dc3+dc2*dc3
      dr(2,3)= c1*s2*c3       +s1*s3
      dr(3,1)=-s1*c2
      dr(3,2)=-s2
      dr(3,3)= dc1+dc2+dc1*dc2
      call tbgrote(trans,cod,dr,0.d0,0.d0,rr)
      return
      end
c$$$      cchi1=cos(chi1)
c$$$      schi1=sin(chi1)
c$$$      cchi2=cos(chi2)
c$$$      schi2=sin(chi2)
c$$$      cchi3=cos(chi3)
c$$$      schi3=sin(chi3)
c$$$      rr(1,1)= cchi1*cchi3+schi1*schi2*schi3
c$$$      rr(1,2)=-cchi2*schi3
c$$$      rr(1,3)= schi1*cchi3-cchi1*schi2*schi3
c$$$      rr(2,1)=-schi1*schi2*cchi3+cchi1*schi3
c$$$      rr(2,2)= cchi2*cchi3
c$$$      rr(2,3)= cchi1*schi2*cchi3+schi1*schi3
c$$$      rr(3,1)=-schi1*cchi2
c$$$      rr(3,2)=-schi2
c$$$      rr(3,3)= cchi1*cchi2
c$$$      xi=cod(1)
c$$$      yi=cod(3)
c$$$      p=1.d0+cod(6)
c$$$      pxi=cod(2)
c$$$      pyi=cod(4)
c$$$      pzi=p*(1.d0+sqrt1(-(pxi**2+pyi**2)/p**2))
c$$$      xf =rr(1,1)*xi +rr(1,2)*yi 
c$$$      yf =rr(2,1)*xi +rr(2,2)*yi 
c$$$      zf =rr(3,1)*xi +rr(3,2)*yi 
c$$$      pxf=rr(1,1)*pxi+rr(1,2)*pyi+rr(1,3)*pzi
c$$$      pyf=rr(2,1)*pxi+rr(2,2)*pyi+rr(2,3)*pzi
c$$$      pzf=rr(3,1)*pxi+rr(3,2)*pyi+rr(3,3)*pzi
c$$$      dpzfdpx=rr(3,1)-rr(3,3)*pxi/pzi
c$$$      dpzfdpy=rr(3,2)-rr(3,3)*pyi/pzi
c$$$      dpzfdp =    rr(3,3)*p  /pzi
c$$$      trans(2,1)=0.d0
c$$$      trans(2,2)=rr(1,1)-rr(1,3)*pxi/pzi
c$$$      trans(2,3)=0.d0
c$$$      trans(2,4)=rr(1,2)-rr(1,3)*pyi/pzi
c$$$      trans(2,5)=0.d0
c$$$      trans(2,6)=    rr(1,3)*p  /pzi
c$$$      trans(4,1)=0.d0
c$$$      trans(4,2)=rr(2,1)-rr(2,3)*pxi/pzi
c$$$      trans(4,3)=0.d0
c$$$      trans(4,4)=rr(2,2)-rr(2,3)*pyi/pzi
c$$$      trans(4,5)=0.d0
c$$$      trans(4,6)=    rr(2,3)*p  /pzi
c$$$      trans(1,1)=rr(1,1)-pxf/pzf*rr(3,1)
c$$$      trans(1,2)=-zf/pzf*(trans(2,2)-pxf/pzf*dpzfdpx)
c$$$      trans(1,3)=rr(1,2)-pxf/pzf*rr(3,2)
c$$$      trans(1,4)=-zf/pzf*(trans(2,4)-pxf/pzf*dpzfdpy)
c$$$      trans(1,5)=0.d0
c$$$      trans(1,6)=-zf/pzf*(trans(2,6)-pxf/pzf*dpzfdp )
c$$$      trans(3,1)=rr(2,1)-pyf/pzf*rr(3,1)
c$$$      trans(3,2)=-zf/pzf*(trans(4,2)-pyf/pzf*dpzfdpx)
c$$$      trans(3,3)=rr(2,2)-pyf/pzf*rr(3,2)
c$$$      trans(3,4)=-zf/pzf*(trans(4,4)-pyf/pzf*dpzfdpy)
c$$$      trans(3,5)=0.d0
c$$$      trans(3,6)=-zf/pzf*(trans(4,6)-pyf/pzf*dpzfdp )
c$$$      trans(5,1)=    p  /pzf*rr(3,1)
c$$$      trans(5,2)= zf/pzf*(           -p  /pzf*dpzfdpx)
c$$$      trans(5,3)=    p  /pzf*rr(3,2)
c$$$      trans(5,4)= zf/pzf*(           -p  /pzf*dpzfdpy)
c$$$      trans(5,5)=1.d0
c$$$      trans(5,6)= zf/pzf*(1.d0       -p  /pzf*dpzfdp )
c$$$      trans(6,1:5)=0.d0
c$$$      trans(6,6)=1.d0
c$$$      cod(1)=xf-pxf/pzf*zf
c$$$      cod(3)=yf-pyf/pzf*zf
c$$$      cod(5)=cod(5)+p/pzf*zf
c$$$      cod(2)=pxf
c$$$      cod(4)=pyf
c$$$      return
c$$$      end
