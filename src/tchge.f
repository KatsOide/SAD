      module chg

      contains
      subroutine tchge(trans,cod,beam,srot,
     $     dx0,dy,dz,theta,dtheta,dchi2,alg,phig,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      use sad_basics
      use mathfun, only:pxy2dpz,xsincos
      implicit none
      integer*4 , parameter :: ms=9
      real*8, intent(inout):: trans(6,12),cod(6),beam(42),srot(3,ms)
      real*8 ,intent(in):: dx0,dy,dz,theta,dtheta,alg,dchi2
      real*8 , target ::trans1(6,6),trans2(6,6),trans3(6,6)
      real*8 sx(ms),cost,sint,phig,th,xi,pxi,dx
      real*8 ,pointer :: transa(:,:)
      logical*4 ,intent(in):: ent
      logical*4 mcal
      mcal=irad .ge. 6
      if(phig /= 0.d0)then
        th=theta
      else
        th=theta+dtheta
      endif
      nullify(transa)
      if(ent)then
        if(th /= 0.d0)then
          cost=cos(th)
          sint=sin(th)
          dx= cost*dx0-sint*dy
          if(mcal)then
            call tinitr(trans1)
            trans1(1,1)= cost
            trans1(1,3)=-sint
            trans1(3,1)= sint
            trans1(3,3)= cost
            trans1(2,2)= cost
            trans1(2,4)=-sint
            trans1(4,2)= sint
            trans1(4,4)= cost
            transa=>trans1
          endif
          xi=cod(1)
          cod(1)= cost*xi-sint*cod(3)
          cod(3)= sint*xi+cost*cod(3)
          pxi=cod(2)
          cod(2)= cost*pxi-sint*cod(4)
          cod(4)= sint*pxi+cost*cod(4)
          if(calpol .and. irad .gt. 6)then
            sx=srot(1,:)
            srot(1,:)= cost*sx-sint*srot(2,:)
            srot(2,:)= sint*sx+cost*srot(2,:)
          endif
          cod(3)= cod(3)-(sint*dx0+cost*dy)
        else
          dx=dx0
          cod(3)=cod(3)-dy
        endif
        if(dx /= 0.d0 .or. dz /= 0.d0 .or. dtheta /= 0.d0 .or. dchi2 /= 0.d0)then
          if(.not. associated(transa))then
            call tinitr(trans3)
            transa=>trans3
          endif
          call tbrote(transa,cod,srot,alg,phig,dx,dz,dtheta,dchi2,.true.)
        endif
      else
        if(th /= 0.d0)then
          cost=cos(th)
          sint=sin(th)
          dx= cost*dx0-sint*dy
        else
          cost=1.d0
          sint=0.d0
          dx=dx0
        endif
        if(dx /= 0.d0 .or. dz /= 0.d0 .or. dtheta /= 0.d0 .or. dchi2 /= 0.d0)then
          call tinitr(trans2)
          call tbrote(trans2,cod,srot,alg,phig,dx,dz,dtheta,dchi2,.false.)
          if(mcal)then
            transa=>trans2
          endif
        endif
        if(th /= 0.d0)then
          cod(3)= cod(3)+(sint*dx0+cost*dy)
          if(mcal)then
            call tinitr(trans1)
            trans1(1,1)= cost
            trans1(1,3)= sint
            trans1(3,1)=-sint
            trans1(3,3)= cost
            trans1(2,2)= cost
            trans1(2,4)= sint
            trans1(4,2)=-sint
            trans1(4,4)= cost
          endif
          xi=cod(1)
          cod(1)= cost*xi+sint*cod(3)
          cod(3)=-sint*xi+cost*cod(3)
          pxi=cod(2)
          cod(2)= cost*pxi+sint*cod(4)
          cod(4)=-sint*pxi+cost*cod(4)
          if(mcal)then
            if(calpol .and. irad .gt. 6)then
              sx=srot(1,:)
              srot(1,:)= cost*sx+sint*srot(2,:)
              srot(2,:)=-sint*sx+cost*srot(2,:)
            endif
            if(associated(transa))then
              call tmultr5(transa,trans1,6)
            else
              transa=>trans1
            endif
          endif
        else
          cod(3)=cod(3)+dy
        endif
      endif
      if(associated(transa))then
        trans(:,1:irad)=matmul(transa,trans(:,1:irad))
        if(irad .gt. 6)then
          call tmulbs(beam,transa,.true.)
        endif
      endif
      return
      end

      subroutine tbrote(trans,cod,srot,alg,phig,dx,dz,dtheta,dchi2,ent)
      use sad_basics
      use tmacro, only:irad
      use ffs_flag,only:calpol
      use mathfun, only:xsincos
      implicit none
      real*8 ,intent(inout):: trans(6,6),cod(6),srot(3,9)
      real*8 ,intent(in):: phig,dtheta,alg,dchi2,dx,dz
      logical*4 ,intent(in):: ent
      real*8 ,target :: trans2(6,6)
      real*8 rr(3,3),s2,xs2,c2,dc2,st,xst,ct,dct,dl,c0,s0,dx1,dz1
      integer*4 i
      if(phig == 0.d0)then
        c0=1.d0
        s0=0.d0
        dl=alg
      else
        c0=cos(phig)
        s0=sin(phig)
        dl=alg*s0/phig
      endif
      call xsincos(dchi2,s2,xs2,c2,dc2)
      call xsincos(dtheta,st,xst,ct,dct)
      rr(1,1)=dct*c0**2+dc2*s0**2+c0*s0*s2*st
      rr(2,2)=c2*dct+dc2
      rr(3,3)=dc2*c0**2+dct*s0**2-c0*s0*s2*st
      if(ent)then
        dx1=dx*c0-dz*s0
        dz1=dz*c0+dx*s0
        rr(1,2)= s0*s2-c0*c2*st
        rr(1,3)= c0*((dct-dc2)*s0-c0*s2*st)
        rr(2,1)=-ct*s0*s2+c0*st
        rr(2,3)= c0*ct*s2+s0*st
        rr(3,1)= s0*(c0*(dct-dc2)+s0*s2*st)
        rr(3,2)=-c0*s2-c2*s0*st
        cod(1)=cod(1)-dx1
        if(dz1 /= 0.d0)then
          call tdrife1(trans,cod,1.d0,0.d0,dz1)
          call tbgrote(trans2,cod,rr,dl*s0,-dl*c0)
          call tmultr5(trans,trans2,6)
        else
          call tbgrote(trans,cod,rr,dl*s0,-dl*c0)
        endif
      else
        rr(2,1)= s0*s2-c0*c2*st
        rr(3,1)= c0*((dct-dc2)*s0-c0*s2*st)
        rr(1,2)=-ct*s0*s2+c0*st
        rr(3,2)= c0*ct*s2+s0*st
        rr(1,3)= s0*(c0*(dct-dc2)+s0*s2*st)
        rr(2,3)=-c0*s2-c2*s0*st
        call tbgrote(trans,cod,rr,dl*s0,-dl*c0)
        dx1=dx*c0+dz*s0
        dz1=dz*c0-dx*s0
        if(dz1 /= 0.d0)then
          call tdrife1(trans2,cod,1.d0,0.d0,-dz1)
          call tmultr5(trans,trans2,6)
        endif
        cod(1)=cod(1)+dx1
      endif
      if(calpol .and. irad .gt. 6)then
        do concurrent (i=1:9)
          srot(1,i)=dot_product(rr(1,:),srot(:,i))+srot(1,i)
          srot(2,i)=dot_product(rr(2,:),srot(:,i))+srot(2,i)
          srot(3,i)=dot_product(rr(3,:),srot(:,i))+srot(3,i)
        enddo
      endif
      return
      end

      subroutine tbgrote(trans,cod,dr,dx,dz,rr)
      use mathfun,only:pxy2dpz
      implicit none
      real*8 ,intent(inout):: trans(6,6),cod(6)
      real*8 ,intent(in):: dr(3,3),dx,dz
      real*8 ,intent(out) , optional ::rr(3,3)
      real*8 pxi,pyi,pzi,xi,yi,zi,xf,yf,zf,pxf,pyf,pzf,dpzi,dpzf,
     $     dpzidpx,dpzidpy,dpzidp,pr,r11,r22,r33,pxfpzf,pyfpzf,zfpzf,prpzf
      xi=cod(1)+dx
      yi=cod(3)
      pxi=cod(2)
      pyi=cod(4)
      pr=1.d0+cod(6)
      dpzi=pxy2dpz(pxi/pr,pyi/pr)*pr
      pzi=pr+dpzi
      zi =dz
      xf =dr(1,1)*xi +dr(1,2)*yi +dr(1,3)*zi +cod(1)
      yf =dr(2,1)*xi +dr(2,2)*yi +dr(2,3)*zi +yi
      zf =dr(3,1)*xi +dr(3,2)*yi +dr(3,3)*zi
      r11=1.d0+dr(1,1)
      r22=1.d0+dr(2,2)
      r33=1.d0+dr(3,3)
      if(present(rr))then
        rr=dr
        rr(1,1)=r11
        rr(2,2)=r22
        rr(3,3)=r33
      endif
      pxf =    r11*pxi+dr(1,2)*pyi+dr(1,3)*pzi
      pyf =dr(2,1)*pxi+    r22*pyi+dr(2,3)*pzi
      dpzf=dr(3,1)*pxi+dr(3,2)*pyi+    r33*dpzi+dr(3,3)*pr
      pzf =pr+dpzf
      pxfpzf=pxf/pzf
      pyfpzf=pyf/pzf
      zfpzf=zf/pzf
      prpzf=pr/pzf
      cod(2)=pxf
      cod(4)=pyf
      cod(1)=xf-pxf*zfpzf
      cod(3)=yf-pyf*zfpzf
      cod(5)=cod(5)+pr*zfpzf

      dpzidpx=-pxi/pzi
      dpzidpy=-pyi/pzi
      dpzidp =  pr/pzi
      trans(2,1)=0.d0
      trans(2,2)=    r11+dr(1,3)*dpzidpx
      trans(2,3)=0.d0
      trans(2,4)=dr(1,2)+dr(1,3)*dpzidpy
      trans(1:4,5)=0.d0
      trans(2,6)=dr(1,3)*dpzidp

      trans(4,1)=0.d0
      trans(4,2)=dr(2,1)+dr(2,3)*dpzidpx
      trans(4,3)=0.d0
      trans(4,4)=    r22+dr(2,3)*dpzidpy
      trans(4,6)=dr(2,3)*dpzidp

      trans(1,1)=    r11-dr(3,1)*pxfpzf
      trans(1,2)=-zfpzf*(trans(2,2)-pxfpzf*(dr(3,1)+    r33*dpzidpx))
      trans(1,3)=dr(1,2)-dr(3,2)*pxfpzf
      trans(1,4)=-zfpzf*(trans(2,4)-pxfpzf*(dr(3,2)+    r33*dpzidpy))
      trans(1,6)=-zfpzf*(dr(1,3)-    r33*pxfpzf)*dpzidp

      trans(3,1)=dr(2,1)-dr(3,1)*pyfpzf
      trans(3,2)=-zfpzf*(trans(4,2)-pyfpzf*(dr(3,1)+    r33*dpzidpx))
      trans(3,3)=    r22-dr(3,2)*pyfpzf
      trans(3,4)=-zfpzf*(trans(4,4)-pyfpzf*(dr(3,2)+    r33*dpzidpy))
      trans(3,6)=-zfpzf*(dr(2,3)-    r33*pyfpzf)*dpzidp

c      trans(5,1)=-trans(1,1)*trans(2,6)-trans(3,1)*trans(4,6)
      trans(5,1)= dr(3,1)*prpzf
c      trans(5,2)=-trans(1,2)*trans(2,6)+trans(2,2)*trans(1,6)
c     $           -trans(3,2)*trans(4,6)+trans(4,2)*trans(3,6)
      trans(5,2)= -zfpzf*prpzf*(dr(3,1)+r33*dpzidpx)
c      trans(5,3)=-trans(1,3)*trans(2,6)-trans(3,3)*trans(4,6)
      trans(5,3)= dr(3,2)*prpzf
c      trans(5,4)=-trans(1,4)*trans(2,6)+trans(2,4)*trans(1,6)
c     $           -trans(3,4)*trans(4,6)+trans(4,4)*trans(3,6)
      trans(5,4)= -zfpzf*prpzf*(dr(3,2)+r33*dpzidpy)
      trans(5,5)=1.d0
c      trans(5,6)=zfpzf*(1.d0-pr/pzf*    r33*dpzidp)
      trans(5,6)=r33*(dpzf*prpzf+dpzi/pzi)-dr(3,3)

      trans(6,1:5)=0.d0
      trans(6,6)=1.d0
      return
      end subroutine

      end module
