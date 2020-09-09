c  The varibles px, py in this level are canonical momenta divided by p.
c 
      module sol
      real*8,save,private::cchi1,cchi2,cchi3,schi1,
     $     schi2,schi3,dcchi1,dcchi2,
     $     s0s,fx,fy

      contains
      subroutine tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,dz,
     $     chi1,chi2,chi3,bxs,bys,bzs,ent)
      use mathfun, only:sqrt1,xsincos,pxy2dpz
      use ffs_flag, only:calpol,rad
      implicit none
      integer*4 ,intent(in):: np
      real*8 ,intent(inout)::
     $     x(np),px(np),y(np),py(np),z(np),g(np),
     $     sx(np),sy(np),sz(np),bxs,bys,bzs
      real*8 , intent(in)::al,bz,dx,dy,dz,
     $     chi1,chi2,chi3
      logical*4 ,intent(in)::ent
      integer*4 i,j
      real*8 b,phix,phiy,phiz,a,a12,a14,a22,a24,cosphi,
     $     dphizsq,pr,ds1,ds2,pz0,pz1,bzp,alb,s,pl,ps2,
     $     dpz0,dpl,plx,ply,plz,ptx,pty,ptz,pbx,pby,pbz,
     $     phi,sinphi,dcosphi,xsinphi,dphi,phi0,bxs0,xs,
     $     px0,pxi,pyi,pz2,x0,x1,y1,sv(3),rr(3,3)
      integer*4 , parameter ::itmax=10
      real*8 ,parameter :: conv=3.d-16
      if(ent)then
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or. chi2 .ne. 0.d0)then
          s0s=-al*.5d0
          call xsincos(chi1,schi1,xs,cchi1,dcchi1)
          call xsincos(chi2,schi2,xs,cchi2,dcchi2)
          if(bz .ne. 0.d0)then
            b=abs(bz)
            phix=(-schi1)         *sign(1.d0,bz)
            phiy=( cchi1)*(-schi2)*sign(1.d0,bz)
            phiz=( cchi1)*( cchi2)*sign(1.d0,bz)
            dphizsq=phix**2+phiy**2
            bxs=phix*b
            bys=phiy*b
            bzs=phiz*b
            do i=1,np
              pr=(1.d0+g(i))
              px(i)=px(i)+bz*y(i)/pr*.5d0
              py(i)=py(i)-bz*x(i)/pr*.5d0
              x(i)=x(i)-dx
              y(i)=y(i)-dy
              ds1  = schi1*x(i)+dcchi1*s0s-cchi1*dz
              x(i) = cchi1*x(i)-schi1*(s0s-dz)
              ds2  =-(schi2*y(i)+cchi2*ds1+dcchi2*s0s)
              y(i) = cchi2*y(i)-schi2*(ds1+s0s)
c              pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
              pz0=1.d0+pxy2dpz(px(i),py(i))
              pz1  = schi1*px(i)+cchi1*pz0
              px(i)= cchi1*px(i)-schi1*pz0
              py(i)= cchi2*py(i)-schi2*pz1
              bzp=bzs/pr
              alb=pr/b
              s=px(i)**2+py(i)**2
              dpz0=sqrt1(-s)
c     dpz0=-s/(1.d0+sqrt(1.d0-s))
              pz0=1.d0+dpz0
              dpl=px(i)*phix+py(i)*phiy+dpz0*phiz
              pl=phiz+dpl
              plx=pl*phix
              ply=pl*phiy
              plz=pl*phiz
              ptx=px(i)-plx
              pty=py(i)-ply
              ptz=dpz0 -dpl*phiz+dphizsq
              pbx=pty*phiz-ptz*phiy
              pby=ptz*phix-ptx*phiz
              pbz=ptx*phiy-pty*phix
              if(ds2 .ne. 0.d0)then
                ps2=ds2/alb
                phi=asin(ps2/pz0)
                do j=1,itmax
                  call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
c                  call sxsin(phi,sinphi,xsinphi)
c                  dcosphi=2.d0*sin(.5d0*phi)**2
                  s=plz*xsinphi+pz0*sinphi-pbz*dcosphi
                  dphi=(ps2-s)/(pz0+ptz*dcosphi+pbz*sinphi)
                  phi0=phi
                  phi=phi+dphi
                  if(phi0 .eq. phi .or. abs(dphi/phi) .lt. conv)then
                    exit
                  endif
                enddo
              else
                phi=0.d0
                sinphi=0.d0
                dcosphi=0.d0
              endif
              x(i)=x(i)+(plx*phi+ptx*sinphi+pbx*dcosphi)*alb
              y(i)=y(i)+(ply*phi+pty*sinphi+pby*dcosphi)*alb
              z(i)=z(i)-phi*alb
              px(i)=px(i)-ptx*dcosphi+pbx*sinphi-bzp*y(i)*.5d0
              py(i)=py(i)-pty*dcosphi+pby*sinphi+bzp*x(i)*.5d0
            enddo        
          else
            bxs=0.d0
            bys=0.d0
            bzs=0.d0
            do i=1,np
              pr=(1.d0+g(i))
              x(i)=x(i)-dx
              y(i)=y(i)-dy
              pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c     pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
              pz1  = schi1*px(i)+cchi1*pz0
              px(i)= cchi1*px(i)-schi1*pz0
              py(i)= cchi2*py(i)-schi2*pz1
              ds1  = schi1*x(i)+dcchi1*s0s-cchi1*dz
              x(i) = cchi1*x(i)-schi1*(s0s-dz)
              ds2  = schi2*y(i)+cchi2*ds1+dcchi2*s0s
              y(i) = cchi2*y(i)-schi2*(ds1+s0s)
              pz2=1.d0+sqrt1(-px(i)**2-py(i)**2)
c     pz2=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
              a=ds2/pz2
              x(i) =x(i)-a*px(i)
              y(i) =y(i)-a*py(i)
              z(i) =z(i)+a
            enddo
          endif
        else
          cchi1=0.d0
          cchi2=0.d0
          schi1=0.d0
          schi2=0.d0
          dcchi1=0.d0
          dcchi2=0.d0
          s0s=0.d0
          bzs=bz
          bxs=0.d0
          bys=0.d0
          fx= bzs*dy*.5d0
          fy=-bzs*dx*.5d0
          x=x-dx
          y=y-dy
          px=px+fx/(1.d0+g)
          py=py+fy/(1.d0+g)
        endif
        if(chi3 .ne. 0.d0)then
          cchi3=cos(chi3)
          schi3=sin(chi3)
          do i=1,np
            x0=x(i)
            x(i)=cchi3*x0-schi3*y(i)
            y(i)=schi3*x0+cchi3*y(i)
            px0=px(i)
            px(i)=cchi3*px0-schi3*py(i)
            py(i)=schi3*px0+cchi3*py(i)
          enddo
          bxs0=bxs
          bxs=cchi3*bxs0-schi3*bys
          bys=schi3*bxs0+cchi3*bys
        endif
      else
        if(chi3 .ne. 0.d0)then
          do i=1,np
            x0=x(i)
            x(i)= cchi3*x0+schi3*y(i)
            y(i)=-schi3*x0+cchi3*y(i)
            px0=px(i)
            px(i)= cchi3*px0+schi3*py(i)
            py(i)=-schi3*px0+cchi3*py(i)
          enddo
        endif
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or. chi2 .ne. 0.d0)then
          if(bz .ne. 0.d0)then
            do i=1,np
              pr=(1.d0+g(i))
              px(i)=px(i)+bzs/pr*y(i)*.5d0
              py(i)=py(i)-bzs/pr*x(i)*.5d0
              pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
c     pz0=sqrt((1.d0-px(i))*(1.d0+px(i))-py(i)**2)
              pz1  =-schi2*py(i)+cchi2*pz0
              pyi  = cchi2*py(i)+schi2*pz0
              pz2  =-schi1*px(i)+cchi1*pz1
              pxi  = cchi1*px(i)+schi1*pz1
              y1   = cchi2*y(i)-schi2*s0s
              ds1  =-schi2*y(i)-dcchi2*s0s
              ds2  =-schi1*x(i)+cchi1*ds1-dcchi1*s0s+dz
              x1   = cchi1*x(i)+schi1*(ds1-s0s)
              bzp=bz/pr
              phi=-bzp*ds2/pz2
              call xsincos(phi,a24,a12,a22,a14)
c              a24=sin(phi)
c              a12=a24/bzp
c              a22=cos(phi)
c              if(a22 .ge. 0.d0)then
c                a14=a24**2/(1.d0+a22)/bzp
c              else
c                a14=(1.d0-a22)/bzp
c              endif
              x(i) =x1 +(a24*pxi-a14*pyi)/bzp+dx
              y(i) =y1 +(a14*pxi+a24*pyi)/bzp+dy
              px(i)=    a22*pxi+a24*pyi-bzp*y(i)*.5d0
              py(i)=   -a24*pxi+a22*pyi+bzp*x(i)*.5d0
              z(i) =z(i)+ds2/pz2
            enddo
          else
            do i=1,np
              pr=(1.d0+g(i))
              pz0=1.d0+sqrt1(-px(i)**2-py(i)**2)
              pz1  =-schi2*py(i)+cchi2*pz0
              pyi  = cchi2*py(i)+schi2*pz0
              pz2  =-schi1*px(i)+cchi1*pz1
              pxi  = cchi1*px(i)+schi1*pz1
              y1   = cchi2*y(i)-schi2*s0s
              ds1  =-schi2*y(i)-dcchi2*s0s
              ds2  =-schi1*x(i)+cchi1*ds1-dcchi1*s0s+dz
              x1   = cchi1*x(i)+schi1*(ds1-s0s)
              a=ds2/pz2
              px(i)=pxi
              py(i)=pyi
              x(i) =x1-a*pxi+dx
              y(i) =y1-a*pyi+dy
              z(i) =z(i)+a
            enddo
          endif
        else
          px=px-fx/(1.d0+g)
          py=py-fy/(1.d0+g)
          x=x+dx
          y=y+dy
        endif
      endif
      if(rad .and. calpol)then
        call trot33(rr,ent)
        do i=1,np
          sv=matmul(rr,(/sx(i),sy(i),sz(i)/))
          sx(i)=sv(1)
          sy(i)=sv(2)
          sz(i)=sv(3)
        enddo
      endif
      return
      end subroutine

      subroutine tsolrote(trans,cod,beam,srot,al,bz,dx,dy,dz,
     $     chi1,chi2,chi3,bxs,bys,bzs,ent)
      use tmacro, only:irad
      use ffs_flag, only:calpol
      use mathfun, only: sqrtl,xsincos
      use temw,only:tmulbs
      implicit none
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 trans1(6,6),trans2(6,6),tb(6)
      real*8 ,intent(in):: bz,dx,dy,chi3,dz,chi1,chi2,al
      real*8 ,intent(out):: bxs,bys,bzs
      real*8 s0,rr(3,3),xs,x0,px0,bzh,
     $     ds1,pr,pz0,dpz0dpx,dpz0dpy,dpz0dp,
     $     pz1,ds2,bxs0,bzh1,pzmin,a,ptmax
      logical*4 ,intent(in):: ent
      logical*4 mul
      parameter (pzmin=1.d-10,ptmax=0.9999d0)
      mul=.false.
c     cod in/outs canonical momenta!
      if(ent)then
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or.
     $       chi2 .ne. 0.d0)then
          call tinitr(trans1)
          mul=.true.
          bzh1=bz*.5d0
          cod(2)=cod(2)+bzh1*cod(3)
          cod(4)=cod(4)-bzh1*cod(1)
          s0=-.5d0*al
          call xsincos(chi1,schi1,xs,cchi1,dcchi1)
          call xsincos(chi2,schi2,xs,cchi2,dcchi2)
          bzs=cchi2*cchi1*bz
          bxs=-schi1*bz
          bys=-schi2*cchi1*bz
          cod(1)=cod(1)-dx
          cod(3)=cod(3)-dy
          pr=1.d0+cod(6)
          a=cod(2)**2+cod(4)**2
          pz0=pr*sqrtl(1.d0-a/pr**2)
          dpz0dpx=-cod(2)/pz0
          dpz0dpy=-cod(4)/pz0
          dpz0dp = pr/pz0
          pz1   = schi1*cod(2)+cchi1*pz0
          cod(2)= cchi1*cod(2)-schi1*pz0
          cod(4)= cchi2*cod(4)-schi2*pz1
          trans1(2,2)= cchi1      -schi1*dpz0dpx
          trans1(2,4)=            -schi1*dpz0dpy
          trans1(2,6)=            -schi1*dpz0dp
          trans1(4,2)=-schi2*schi1-schi2*cchi1*dpz0dpx
          trans1(4,4)= cchi2      -schi2*cchi1*dpz0dpy
          trans1(4,6)=            -schi2*cchi1*dpz0dp
          ds1   = schi1*cod(1)+dcchi1*s0-cchi1*dz
          cod(1)= cchi1*cod(1)-schi1*(s0-dz)
          ds2   = schi2*cod(3)+cchi2*ds1+dcchi2*s0
          cod(3)= cchi2*cod(3)-schi2*(ds1+s0)
          trans1(1,1)= cchi1
          trans1(3,1)=-schi2*schi1
          trans1(3,3)= cchi2
          trans1(5,1)=-cchi2*schi1
          trans1(5,3)=-schi2
          trans1(5,5)=0.d0
          call tsoldz(trans2,cod,-ds2,bxs,bys,bzs,.false.)
          trans1=matmul(trans2,trans1)
          trans1(5,5)=1.d0
          bzh=bzs*.5d0
          cod(2)=cod(2)-bzh*cod(3)
          cod(4)=cod(4)+bzh*cod(1)
          trans1(:,1)=trans1(:,1)-bzh1*trans1(:,4)
          trans1(:,3)=trans1(:,3)+bzh1*trans1(:,2)
          trans1(2,:)=trans1(2,:)-bzh*trans1(3,:)
          trans1(4,:)=trans1(4,:)+bzh*trans1(1,:)
        else
          bzh=bz*.5d0
          cod(1)=cod(1)-dx
          cod(3)=cod(3)-dy
          cod(2)=cod(2)+bzh*dy
          cod(4)=cod(4)-bzh*dx
          cchi1=1.d0
          schi1=0.d0
          cchi2=1.d0
          schi2=0.d0
          dcchi1=0.d0
          dcchi2=0.d0
          bxs=0.d0
          bys=0.d0
          bzs=bz
        endif
        if(chi3 .ne. 0.d0)then
          if(.not. mul)then
            call tinitr(trans1)
          endif
          mul=.true.
          cchi3=cos(chi3)
          schi3=sin(chi3)
          x0=cod(1)
          cod(1)=cchi3*x0-schi3*cod(3)
          cod(3)=schi3*x0+cchi3*cod(3)
          px0=cod(2)
          cod(2)=cchi3*px0-schi3*cod(4)
          cod(4)=schi3*px0+cchi3*cod(4)
          tb=trans1(1,:)
          trans1(1,:)=cchi3*tb-schi3*trans1(3,:)
          trans1(3,:)=schi3*tb+cchi3*trans1(3,:)
          tb=trans1(2,:)
          trans1(2,:)=cchi3*tb-schi3*trans1(4,:)
          trans1(4,:)=schi3*tb+cchi3*trans1(4,:)
          bxs0=bxs
          bxs=cchi3*bxs0-schi3*bys
          bys=schi3*bxs0+cchi3*bys
        else
          cchi3=1.d0
          schi3=0.d0
        endif
      else
        if(chi3 .ne. 0.d0)then
          call tinitr(trans1)
          mul=.true.
          x0=cod(1)
          cod(1)= cchi3*x0+schi3*cod(3)
          cod(3)=-schi3*x0+cchi3*cod(3)
          px0=cod(2)
          cod(2)= cchi3*px0+schi3*cod(4)
          cod(4)=-schi3*px0+cchi3*cod(4)
          tb=trans1(1,:)
          trans1(1,:)= cchi3*tb+schi3*trans1(3,:)
          trans1(3,:)=-schi3*tb+cchi3*trans1(3,:)
          tb=trans1(2,:)
          trans1(2,:)= cchi3*tb+schi3*trans1(4,:)
          trans1(4,:)=-schi3*tb+cchi3*trans1(4,:)
        else
          cchi3=1.d0
          schi3=0.d0
        endif
        if(dz .ne. 0.d0 .or. chi1 .ne. 0.d0 .or.
     $       chi2 .ne. 0.d0)then
          if(.not. mul)then
            call tinitr(trans1)
          endif
          mul=.true.
          s0=-.5d0*al
          call tinitr(trans2)
          pr=1.d0+cod(6)
          bzh=bzs*.5d0
          cod(2)=cod(2)+bzh*cod(3)
          cod(4)=cod(4)-bzh*cod(1)
          a=cod(2)**2+cod(4)**2
          pz0=pr*sqrtl(1.d0-a/pr**2)
c          pz0=sqrt(max(pzmin,(pr-cod(2))*(pr+cod(2))-cod(4)**2))
          dpz0dpx=-cod(2)/pz0
          dpz0dpy=-cod(4)/pz0
          dpz0dp = pr/pz0
          ds1   =-schi2*cod(3)-dcchi2*s0
          cod(3)= cchi2*cod(3)-schi2*s0
          ds2   =-schi1*cod(1)+cchi1*ds1-dcchi1*s0+dz
          cod(1)= cchi1*cod(1)+schi1*(ds1-s0)
          pz1   =-schi2*cod(4)+cchi2*pz0
          cod(2)= cchi1*cod(2)+schi1*pz1
          cod(4)= cchi2*cod(4)+schi2*pz0
          trans2(1,1)= cchi1
          trans2(1,3)=-schi1*schi2
          trans2(3,3)= cchi2
          trans2(5,1)= schi1
          trans2(5,3)= cchi1*schi2
          trans2(2,2)= cchi1         +schi1*cchi2*dpz0dpx
          trans2(2,4)=-schi1*schi2   +schi1*cchi2*dpz0dpy
          trans2(2,3)= trans2(2,2)*bzh
          trans2(2,1)=-trans2(2,4)*bzh
          trans2(2,6)=                schi1*cchi2*dpz0dp
          trans2(4,2)=                schi2*dpz0dpx
          trans2(4,4)= cchi2         +schi2*dpz0dpy
          trans2(4,6)=                schi2*dpz0dp
          trans2(4,3)= trans2(4,2)*bzh
          trans2(4,1)=-trans2(4,4)*bzh
          trans1(5,5)=0.d0
          call tmultr5(trans1,trans2,6)
          call tsoldz(trans2,cod,-ds2,0.d0,0.d0,bz,.false.)
          trans1=matmul(trans2,trans1)
c          call tmultr(trans1,trans2,6)
          trans1(5,5)=1.d0
          cod(1)=cod(1)+dx
          cod(3)=cod(3)+dy
          bzh1=bz*.5d0
          cod(2)=cod(2)-bzh1*cod(3)
          cod(4)=cod(4)+bzh1*cod(1)
          trans1(2,:)=trans1(2,:)-bzh1*trans1(3,:)
          trans1(4,:)=trans1(4,:)+bzh1*trans1(1,:)
        else
          bzh=bz*.5d0
          cod(2)=cod(2)-bzh*dy
          cod(4)=cod(4)+bzh*dx
          cod(1)=cod(1)+dx
          cod(3)=cod(3)+dy
        endif
      endif
      if(mul)then
        call tmultr5(trans,trans1,irad)
        if(irad .gt. 6)then
          call tmulbs(beam,trans1,.true.)
          if(calpol)then
            call trot33(rr,ent)
            srot=matmul(rr,srot)
          endif
        endif
      endif
      return
      end subroutine

      subroutine tsoldz(trans,cod,al,bxs0,bys0,bzs0,drift)
      use mathfun
      implicit none
      integer*4 ,save::ndiag=15
      integer*4 ,parameter ::itmax=15
      real*8 ,intent(inout):: trans(6,6),cod(6)
      real*8 ,intent(in):: al,bxs0,bys0,bzs0
      real*8 bxs,bys,bzs,pxi,pyi,pz0,
     $     dpz0dpx,dpz0dpy,dpz0dp,phi,
     $     dphidz,dphidpx,dphidpy,dphidp,a24,a12,a22,a14,
     $     da12,da14,pr,dpz0,dv,dvdp,phix,phiy,phiz,babs,
     $     alb,pbx,pby,pbz,pl,dpl,dphizsq,a,r,
     $     dpldpx,dpldpy,dpldp,dplz,plx,ply,plz,ptx,pty,ptz,
     $     cosphi,sinphi,dcosphi,dphi,
     $     xsinphi,ax,ay,az,cx,cy,albabs,ap,bpr,db
      real*8 ,parameter ::conv=1.d-15,bzthre=1.d-20,ptmax=0.9999d0
      logical*4 ,intent(in):: drift
      bxs=bxs0
      bys=bys0
      bzs=bzs0
      babs=norm2([bzs,bxs,bys])
      if(abs(babs) .lt. bzthre)then
        bxs=0.d0
        bys=0.d0
        bzs=0.d0
        babs=0.d0
      endif
      call tinitr(trans)
      pr=1.d0+cod(6)
c cod does NOT have canonical momenta!
      pxi=cod(2)
      pyi=cod(4)
      a=pxi**2+pyi**2
      dpz0=-a/pr/(1.d0+sqrtl(1.d0-a/pr**2))
      pz0=pr+dpz0
      r=al/pz0
      dpz0dpx= -pxi/pz0
      dpz0dpy= -pyi/pz0
      dpz0dp =   pr/pz0
      if(bxs .eq. 0.d0 .and. bys .eq. 0.d0)then
        phi=bzs*r
        dphidz  = 1.d0/pz0
        dphidpx = -r/pz0*dpz0dpx
        dphidpy = -r/pz0*dpz0dpy
        dphidp  = -r/pz0*dpz0dp
        if(bzs .eq. 0.d0)then
          a24=0.d0
          a12=r
          a22=1.d0
          a14=0.d0
          da12=1.d0
          da14=0.d0
        else
          a24=sin(phi)
          a12=a24/bzs
          a22=cos(phi)
          a14=merge(a24**2/(1.d0+a22),1.d0-a22,a22 .eq. 0.d0)/bzs
          da12=a22
          da14=a24
        endif
        cod(1)=cod(1)+a12*pxi+a14*pyi
        cod(3)=cod(3)-a14*pxi+a12*pyi
        cod(2)=       a22*pxi+a24*pyi
        cod(4)=      -a24*pxi+a22*pyi
        trans(1,2)=( da12*pxi+da14*pyi)*dphidpx+a12
        trans(1,4)=( da12*pxi+da14*pyi)*dphidpy+a14
        trans(1,6)=( da12*pxi+da14*pyi)*dphidp
        trans(3,2)=(-da14*pxi+da12*pyi)*dphidpx-a14
        trans(3,4)=(-da14*pxi+da12*pyi)*dphidpy+a12
        trans(3,6)=(-da14*pxi+da12*pyi)*dphidp
        trans(2,2)=( -a24*pxi+ a22*pyi)*bzs*dphidpx+a22
        trans(2,4)=( -a24*pxi+ a22*pyi)*bzs*dphidpy+a24
        trans(2,6)=( -a24*pxi+ a22*pyi)*bzs*dphidp
        trans(4,2)=( -a22*pxi- a24*pyi)*bzs*dphidpx-a24
        trans(4,4)=( -a22*pxi- a24*pyi)*bzs*dphidpy+a22
        trans(4,6)=( -a22*pxi- a24*pyi)*bzs*dphidp
        trans(5,2)=r*pr/pz0*dpz0dpx
        trans(5,4)=r*pr/pz0*dpz0dpy
        if(drift)then
          call tgetdv(cod(6),dv,dvdp)
          cod(5)=cod(5)+al*(dpz0/pz0-dv)
          trans(5,6)=al*(a/pz0**3+dvdp)
        else
          trans(1,5)=( da12*pxi+da14*pyi)*dphidz
          trans(3,5)=(-da14*pxi+da12*pyi)*dphidz
          trans(2,5)=( -a24*pxi+ a22*pyi)*bzs*dphidz
          trans(4,5)=( -a22*pxi- a24*pyi)*bzs*dphidz
          cod(5)=cod(5)-r*pr
          trans(5,5)=-pr/pz0
          trans(5,6)= r*a/pz0**2
        endif
      else
        phix=bxs/babs
        phiy=bys/babs
        phiz=bzs/babs
        alb=1.d0/babs
        albabs=al*babs
        dphizsq=phix**2+phiy**2
        dpl=pxi*phix+pyi*phiy+dpz0*phiz
        pl=pr*phiz+dpl
        dpldpx=phix+phiz*dpz0dpx
        dpldpy=phiy+phiz*dpz0dpy
        dpldp =     phiz*dpz0dp
        plx=pl*phix
        ply=pl*phiy
        plz=pl*phiz
        ptx=pxi-plx
        pty=pyi-ply
        ptz=dpz0 -dpl*phiz+pr*dphizsq
        pbx=pty*phiz-ptz*phiy
        pby=ptz*phix-ptx*phiz
        pbz=ptx*phiy-pty*phix
        bpr=albabs/pr
        db=bpr-pbz
        ap=hypot(pz0,pbz)
        dphi=-atan(pbz,pz0)
        phi=asin(min(1.d0,max(-1.d0,db/ap)))-dphi
        if(al .ne. 0.d0)then
          if(plz .eq. 0.d0)then
            call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
          else
            call tsolconv(pz0,plz,pbz,bpr,
     $           phi,sinphi,xsinphi,cosphi,dcosphi,ndiag)
          endif
        else
          phi=0.d0
          xsinphi=0.d0
          sinphi=0.d0
          dcosphi=0.d0
          cosphi=1.d0
        endif
        dplz=-pr*dphizsq+dpl*phiz
        ax=pxi+ptx*dcosphi+pbx*sinphi
        ay=pyi+pty*dcosphi+pby*sinphi
        az=pz0+ptz*dcosphi+pbz*sinphi
        dphidpx=-(dpldpx*phiz*xsinphi+dpz0dpx*sinphi
     $       -phiy*dcosphi)/az
        dphidpy=-(dpldpy*phiz*xsinphi+dpz0dpy*sinphi
     $       +phix*dcosphi)/az
        dphidp =-(dpldp *phiz*xsinphi+dpz0dp*sinphi)/az
        dphidz =babs/az
        cod(1)=cod(1)+(plx*phi+ptx*sinphi-pbx*dcosphi)*alb
        cod(3)=cod(3)+(ply*phi+pty*sinphi-pby*dcosphi)*alb
        cod(2)=pxi+ptx*dcosphi+pbx*sinphi
        cod(4)=pyi+pty*dcosphi+pby*sinphi
c cod does NOT have canonical momenta!
        trans(1,2)=alb*(dpldpx*phix*xsinphi+sinphi
     $       +dpz0dpx*phiy       *dcosphi+ax*dphidpx)
        trans(1,4)=alb*(dpldpy*phix*xsinphi
     $       -(phiz-dpz0dpy*phiy)*dcosphi+ax*dphidpy)
        trans(1,6)=alb*(dpldp *phix*xsinphi
     $       +dpz0dp *phiy       *dcosphi+ax*dphidp )
        trans(3,2)=alb*(dpldpx*phiy*xsinphi
     $       -(dpz0dpx*phix-phiz)*dcosphi+ay*dphidpx)
        trans(3,4)=alb*(dpldpy*phiy*xsinphi+sinphi
     $       -dpz0dpy*phix       *dcosphi+ay*dphidpy)
        trans(3,6)=alb*(dpldp *phiy*xsinphi
     $       -dpz0dp *phix       *dcosphi+ay*dphidp )
        cx=-ptx*sinphi+pbx*cosphi
        cy=-pty*sinphi+pby*cosphi
        trans(2,2)=cosphi-dpldpx*phix*dcosphi
     $       -dpz0dpx*phiy*sinphi+cx*dphidpx
        trans(2,4)=      -dpldpy*phix*dcosphi
     $       +(phiz-dpz0dpy*phiy)*sinphi+cx*dphidpy
        trans(2,6)=      -dpldp *phix*dcosphi
     $       -dpz0dp *phiy*sinphi+cx*dphidp
        trans(4,2)=      -dpldpx*phiy*dcosphi
     $       +(dpz0dpx*phix-phiz)*sinphi+cy*dphidpx
        trans(4,4)=cosphi-dpldpy*phiy*dcosphi
     $       +dpz0dpy*phix*sinphi+cy*dphidpy
        trans(4,6)=      -dpldp *phiy*dcosphi
     $       +dpz0dp *phix*sinphi+cy*dphidp
        trans(5,2)= -pr*alb*dphidpx
        trans(5,4)= -pr*alb*dphidpy
        if(drift)then
          call tgetdv(cod(6),dv,dvdp)
          cod(5)=cod(5)+((dpl*phiz-dphizsq*pr)*xsinphi
     $         +dpz0*sinphi-pbz*dcosphi)*alb-dv*al
          trans(5,6)= -alb*(phi+pr*dphidp)+al*dvdp
        else
          trans(1,5)= alb*ax*dphidz
          trans(3,5)= alb*ay*dphidz
          trans(2,5)= cx*dphidz
          trans(4,5)= cy*dphidz
          cod(5)=cod(5)-pr*phi*alb
          trans(5,6)= -alb*(phi+pr*dphidp)
          trans(5,5)= -pr*alb*dphidz
        endif
      endif
      return
      end subroutine

      subroutine trot33(rr,ent)
      implicit none
      real*8, intent(out)::rr(3,3)
      logical*4 , intent(in)::ent
      rr(1,1)= cchi1*cchi3+schi1*schi2*schi3
      rr(2,2)= cchi2*cchi3
      rr(3,3)= cchi1*cchi2
      if(ent)then
        rr(1,2)=-cchi2*schi3
        rr(1,3)= schi1*cchi3-cchi1*schi2*schi3
        rr(2,1)=-schi1*schi2*cchi3+cchi1*schi3
        rr(2,3)= cchi1*schi2*cchi3+schi1*schi3
        rr(3,1)=-schi1*cchi2
        rr(3,2)=-schi2
      else
        rr(2,1)=-cchi2*schi3
        rr(3,1)= schi1*cchi3-cchi1*schi2*schi3
        rr(1,2)=-schi1*schi2*cchi3+cchi1*schi3
        rr(3,2)= cchi1*schi2*cchi3+schi1*schi3
        rr(1,3)=-schi1*cchi2
        rr(2,3)=-schi2
      endif
      return
      end subroutine

      subroutine tsolconv(pz0,plz,pbz,bpr,
     $     phi,sinphi,xsinphi,cosphi,dcosphi,ndiag)
      use mathfun
      implicit none
      integer*4 ,parameter :: itmax=20
      real*8 ,parameter :: conv=1.d-15
      real*8 ,intent(in):: pz0,plz,pbz,bpr
      real*8 ,intent(out):: phi,sinphi,xsinphi,cosphi,dcosphi
      real*8 dphi,s,u,phi0
      integer*4 ,intent(inout):: ndiag
      integer*4 j
      phi=asin(min(1.d0,max(-1.d0,(bpr-pbz)/hypot(pz0,pbz))))
     $     +atan(pbz,pz0)
      if(plz .eq. 0.d0)then
        call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
      else
        do j=1,itmax
          call xsincos(phi,sinphi,xsinphi,cosphi,dcosphi)
          s=plz*xsinphi+pz0*sinphi-pbz*dcosphi
          u=-plz*dcosphi+pz0*cosphi+pbz*sinphi
          dphi=merge((bpr-s)/u,(bpr-s)/pz0,u .ne. 0.d0)
          phi0=phi
          phi=phi+dphi
          if(phi0 .eq. phi .or. abs(dphi) .le. conv*abs(phi))then
            return
          endif
        enddo
        if(ndiag .ge. 0)then
          ndiag=ndiag-1
          write(*,'(a,1p8g13.5)')'tsolconv convergence error: ',
     $         phi,dphi,bpr,s,u,plz,pz0,pbz
          if(ndiag .eq. -1)then
            write(*,*)
     $           'Further tsolconv messages will be suppressed.'
          endif
        endif
      endif
      return
      end

      subroutine tradks(np,x,px,y,py,z,g,dv,sx,sy,sz,bzs,al,bzph)
      use kradlib, only:tradk
      implicit none
      integer*4 , intent(in)::np
      real*8 , intent(in)::bzs,al
      real*8 ,intent(inout)::
     $     x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np)
      real*8 , intent(out)::bzph(np)
      bzph=.5d0*bzs/(1.d0+g)
      px=px+bzph*y
      py=py-bzph*x
      call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,0.d0)
      bzph=.5d0*bzs/(1.d0+g)
      px=px-bzph*y
      py=py+bzph*x
      return
      end

      end module

      subroutine tsol(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     k,kstop,ke,insol,kptbl,la,n,
     $     nwak,nextwake,out)
      use kyparam
      use tfstk
      use ffs
      use ffs_wake
      use ffs_pointer,
     $     only:idelc,direlc,idtypec,idvalc,pnamec,lpnamec,compelc
      use sad_main
      use tparastat
      use ffs_seg
      use kradlib
      use sol
      implicit none
      real*8 conv
      parameter (conv=3.d-16)
      type (sad_comp), pointer::cmp
      type (sad_dlist), pointer :: lsegp
      integer*4 la1,la
      parameter (la1=15)
      integer*4 k,kbz
      integer*4 , intent(inout)::np
      real*8 ,intent(inout)::
     $     x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np)
      real*8 tfbzs,fw,bzs,al,theta,phi,phix,phiy,
     $     bz1,dx,dy,rot,rtaper,ph,bzph(np)
      integer*4 kptbl(np0,6),nwak,nextwake,n,
     $     i,ke,l,lt,itab(np),izs(np),
     $     kdx,kdy,krot,kstop,kb,lwl,lwt,irtc
      integer*8 iwpl,iwpt
      logical*4 , intent(inout) :: insol
      logical*4 out,seg,autophi,krad
      logical*1,save::dofr(0:1)=[.false.,.false.]
      real*8 ,save::dummy(256)=0.d0
      kb=merge(k,k+1,insol)
      do i=kb,nlat
        if(idtypec(i) .eq. icSOL)then
          if(rlist(idvalc(i)+ky_BND_SOL)
     $         .ne. 0.d0)then
            ke=i
            go to 20
          endif
        endif
      enddo
      write(*,*)' ???-TRACK-?Missing end of solenoid ',
     $     pname(idelc(k))(1:lpnamec(k))
      ke=nlat
 20   bzs=tfbzs(k,kbz)
      if(.not. insol)then
        call compelc(k,cmp)
        if(.not. cmp%update)then
          call tpara(cmp)
        endif
        call trots(np,x,px,y,py,z,dv,
     $       cmp%value(p_R11_SOL),
     $       cmp%value(p_R12_SOL),
     $       cmp%value(p_R13_SOL),
     $       cmp%value(p_R21_SOL),
     $       cmp%value(p_R22_SOL),
     $       cmp%value(p_R23_SOL),
     $       cmp%value(p_R31_SOL),
     $       cmp%value(p_R32_SOL),
     $       cmp%value(p_R33_SOL),
     $       cmp%value(ky_DX_SOL),
     $       cmp%value(ky_DY_SOL),
     $       cmp%value(ky_DZ_SOL),
     $       .true.)
        krad=rad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0
     $       .and. cmp%value(ky_F1_SOL) .ne. 0.d0 .and. bzs .ne. 0.d0
        if(krad)then
          pxr0=px
          pyr0=py
          zr0=z
          bsi=0.d0
        endif
        if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
          call tsfrin(np,x,px,y,py,z,g,bzs)
        endif
        if(krad)then
          call tradks(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         bzs,cmp%value(ky_F1_SOL),bzph)
c          call tserad(np,x,px,y,py,g,dv,l1,rho)
        endif
        insol=.true.
      endif
      iwpt=0
      iwpl=0
      do l=kb,min(ke,kstop)
        l_track=l
        if(la .le. 0)then
          call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $         kptbl,np,n,
     $         0.d0,0.d0,0.d0,0.d0,
     $         -alost,-alost,alost,alost,0.d0,0.d0,0.d0,0.d0)
          if(np .le. 0)then
            return
          endif
          la=la1
        else
          la=la-1
        endif
        lt=idtypec(l)
        call compelc(l,cmp)
        seg=tcheckseg(cmp,lt,al,lsegp,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(l,irtc)
          return
        endif
        if(l .eq. nextwake)then
          iwpl=abs(kwaketbl(1,nwak))
          lwl=merge((ilist(1,iwpl-1)-2)/2,0,iwpl .ne. 0)
          iwpt=abs(kwaketbl(2,nwak))
          lwt=merge((ilist(1,iwpt-1)-2)/2,0,iwpt .ne. 0)
          fw=(abs(charge)*e*pbunch*anbunch/amass)/np0*.5d0
          kdx=kytbl(kwDX,lt)
          dx=merge(cmp%value(kdx),0.d0,kdx .ne. 0)
          kdy=kytbl(kwDY,lt)
          dy=merge(cmp%value(kdy),0.d0,kdy .ne. 0)
          krot=kytbl(kwROT,lt)
          rot=merge(cmp%value(krot),0.d0,krot .ne. 0)
          call txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         dx,dy,rot,
     $         int(anbunch),
     $         fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $         p0,h0,itab,izs,.true.)
        endif
        select case (lt)
        case (icDRFT)
          al=cmp%value(ky_L_DRFT)
          if(spac)then
            if(bzs .ne. 0.d0)then
              call spdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,al,bzs,
     $             cmp%value(ky_RADI_DRFT),n,kptbl)
            else
              call spdrift_free(np,x,px,y,py,z,g,dv,sx,sy,sz,al,
     $             cmp%value(ky_RADI_DRFT),n,kptbl)
            endif
          elseif(bzs .ne. 0.d0)then
            call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           al,bzs,rad .and. cmp%value(ky_RAD_DRFT) .eq. 0.d0)
          else
            call tdrift_free(np,x,px,y,py,z,dv,al)
          endif
        case (icBEND)
          if(.not. cmp%update)then
            call tpara(cmp)
          endif
          al=cmp%value(ky_L_BEND)
          theta=cmp%value(ky_ROT_BEND)+cmp%value(ky_DROT_BEND)
          phi=cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND)
          phiy= phi*cos(theta)
          phix= phi*sin(theta)
          call tdrift(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al,bzs,phiy,phix,rad .and. al .ne. 0.d0
     $         .and. cmp%value(ky_RAD_BEND) .eq. 0.d0)
        case(icQUAD)
          if(.not. cmp%update)then
            call tpara(cmp)
          endif
          al=cmp%value(ky_L_QUAD)
          rtaper=1.d0
          krad=rad
          if(rad)then
            if(radcod .and. radtaper)then
              rtaper=1.d0-dp0
     $             +(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
            endif
            krad=cmp%value(ky_RAD_QUAD) .eq. 0.d0 .and. al .ne. 0.d0
          endif
          call tquad(np,x,px,y,py,z,g,dv,sx,sy,sz,al,
     $         cmp%value(ky_K1_QUAD)*rtaper,bzs,
     $         cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1         cmp%value(ky_ROT_QUAD),
     1         cmp%value(p_THETA2_QUAD),
     1         krad,.true.,cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $         cmp%value(p_AKF1F_QUAD)*rtaper,
     $         cmp%value(p_AKF2F_QUAD)*rtaper,
     $         cmp%value(p_AKF1B_QUAD)*rtaper,
     $         cmp%value(p_AKF2B_QUAD)*rtaper,
     $         cmp%ivalue(1,p_FRMD_QUAD),cmp%value(ky_EPS_QUAD),
     $         cmp%value(ky_KIN_QUAD) .eq. 0.d0)

        case (icMULT)
          rtaper=1.d0
          if(rad .and. radcod .and. radtaper)then
            rtaper=1.d0-dp0
     $           +(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          endif
          if(seg)then
            call tmultiseg(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           cmp,lsegp,bzs,rtaper,n,kptbl)
          else
            call tmulti1(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           cmp,bzs,rtaper,n,kptbl)
          endif
        case(icCAVI)
          if(.not. cmp%update)then
            call tpara(cmp)
          endif
          autophi=cmp%value(ky_APHI_CAVI) .ne. 0.d0
          ph=cmp%value(ky_DPHI_CAVI)
          if(autophi)then
            ph=ph+gettwiss(mfitdz,l_track)*cmp%value(p_W_CAVI)
          endif
          call tmultiacc(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         cmp%value(ky_L_CAVI),dummy,dummy,
     $         bzs,0.d0,0.d0,0.d0,
     1         cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $         0.d0,0.d0,0.d0,
     $         cmp%value(ky_ROT_CAVI),0.d0,cmp%value(ky_ROT_CAVI),
     $         0.d0,.false.,
     $         cmp%value(ky_FRIN_CAVI) .eq. 0.d0,
     $         0.d0,0.d0,0.d0,0.d0,
     $         cmp%ivalue(1,p_FRMD_CAVI),
     $         0.d0,0.d0,dofr,
     $         cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI),
     $         cmp%value(p_W_CAVI),
     $         cmp%value(ky_PHI_CAVI),ph,cmp%value(p_VNOMINAL_CAVI),
     $         0.d0,1.d0,autophi,1,0,
     $         n,kptbl)
        case (icSOL)
          bz1=merge(0.d0,tfbzs(l,kbz),l .eq. ke)
          krad=rad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0
     $         .and. cmp%value(ky_F1_SOL) .ne. 0.d0 .and. bzs .ne. bz1          
          if(krad)then
            bzph=.5d0*bzs/(1.d0+g)
            pxr0=px+bzph*y
            pyr0=py-bzph*x
            zr0=z
          endif
          if(cmp%value(ky_FRIN_SOL) .eq. 0.d0)then
            call tsfrin(np,x,px,y,py,z,g,bz1-bzs)
          endif
          if(l .eq. ke)then
            if(.not. cmp%update)then
              call tpara(cmp)
            endif
            if(krad)then
              call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $             cmp%value(ky_F1_SOL),0.d0)
c              call tserad(np,x,px,y,py,g,dv,lp,rho)
            endif
            call trots(np,x,px,y,py,z,dv,
     $           cmp%value(p_R11_SOL),
     $           cmp%value(p_R12_SOL),
     $           cmp%value(p_R13_SOL),
     $           cmp%value(p_R21_SOL),
     $           cmp%value(p_R22_SOL),
     $           cmp%value(p_R23_SOL),
     $           cmp%value(p_R31_SOL),
     $           cmp%value(p_R32_SOL),
     $           cmp%value(p_R33_SOL),
     $           cmp%value(ky_DX_SOL),
     $           cmp%value(ky_DY_SOL),
     $           cmp%value(ky_DZ_SOL),
     $           .false.)
          elseif(krad)then
            call tradks(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $           bz1,cmp%value(ky_F1_SOL),bzph)
c     call tserad(np,x,px,y,py,g,dv,lp,rhoe)
          endif
          bzs=bz1
        case(icMAP)
          call temap(np,np0,x,px,y,py,z,g,dv,sx,sy,sz,l,n,kptbl)
        case(icAprt)
          call tapert1(x,px,y,py,z,g,dv,sx,sy,sz,
     1         kptbl,np,n)
          if(np .le. 0)then
            return
          endif
        case(icMARK, icMONI)

        case default
          if(out)then
            out=.false.
            write(*,*)
     $'Only DRIFT, BEND, QUAD, SOL, MULT, '//
     $'CAVI, MAP, APERT, MARK, MON '//
     $'are supported within SOL: ',
     $           lt
          endif
        end select

        if(l .eq. nextwake .and. l .ne. ke)then
          call txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         dx,dy,rot,
     $         int(anbunch),
     $         fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $         p0,h0,itab,izs,.false.)
          nwak=nwak+1
          nextwake=merge(0,iwakeelm(nwak),nwak .gt. nwakep)
        endif
      enddo
      return
      end

      subroutine trots(np,x,px,y,py,z,dv,
     $     r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $     dx,dy,dz,ent)
      use mathfun
      implicit none
      integer*4 np,i
      real*8,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np)
      real*8 ,intent(in)::dx,dy,dz,r11,r12,r13,r21,r22,r23,r31,r32,r33
      real*8 pxi,pyi,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf
      logical*4,intent(in):: ent
      if(ent)then
        do concurrent (i=1:np)
          pxi=px(i)
          pyi=py(i)
          pzi=1.d0+pxy2dpz(pxi,pyi)
          xi=x(i)
          yi=y(i)
          xf =r11*xi +r12*yi
          yf =r21*xi +r22*yi
          zf =r31*xi +r32*yi
          pxf=r11*pxi+r12*pyi+r13*pzi
          pyf=r21*pxi+r22*pyi+r23*pzi
          pzf=r31*pxi+r32*pyi+r33*pzi
          px(i)=pxf
          py(i)=pyf
          x(i)=xf-pxf/pzf*zf+dx
          y(i)=yf-pyf/pzf*zf+dy
          z(i)=z(i)+zf/pzf+dz
        enddo
      else
        do concurrent (i=1:np)
          pxi=px(i)
          pyi=py(i)
          pzi=1.d0+pxy2dpz(pxi,pyi)
          xi=x(i)-dx
          yi=y(i)-dy
          xf =r11*xi +r12*yi
          yf =r21*xi +r22*yi
          zf =r31*xi +r32*yi
          pxf=r11*pxi+r12*pyi+r13*pzi
          pyf=r21*pxi+r22*pyi+r23*pzi
          pzf=r31*pxi+r32*pyi+r33*pzi
          x(i)=xf-pxf/pzf*zf
          y(i)=yf-pyf/pzf*zf
          z(i)=z(i)+zf/pzf+dz-dv(i)*dz
          px(i)=pxf
          py(i)=pyf
        enddo
      endif
      return
      end

      subroutine tsfrin(np,x,px,y,py,z,g,dbz)
      implicit none
      integer*4 ,intent(in):: np
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ,intent(in):: dbz
      integer*4 i
      real*8 x0,y0,px0,py0,z0,p,bq,bp,pr,pphi,phi0,w,c,s
c     Apply identity map if dbz == 0
      if(dbz .eq. 0.d0) then
        return
      endif
      do i=1,np
c       Map definition
c         [position @ enter]
c          x0  :=  x(i)
c         px0  := px(i)  [Px / p]
c          y0  :=  y(i)
c         py0  := py(i)  [Py / p]
c          z0  :=  z(i)
c          p   := 1 + g(i)
c
c         r    := hypot(x0, y0)
c         b    := .5 * dbz / p
c         bq   := -.25 * b
c         bp   :=  .5  * b
c         pr   := x0 * px0 + y0 * py0
c         pphi := x0 * py0 - y0 * px0
c         r1   := r * exp(bq * pphi)
c         phi  := atan(y0, x0)
c         phi1 := phi + bq * pr
c
c         [position @ exit]
c          x1  := x(i)  <- r1 * cos(phi1)
c          y1  := y(i)  <- r1 * sin(phi1)
c         px1  := px(i) <- (x1 * pr - y1 * pphi) / r1**2
c         py1  := py(i) <- (y1 * pr + x1 * pphi) / r1**2
c          z1  := z(i)  <- z(i) + bp * pr * pphi
c
c       Map expansion
c         phi = atan2(y0, x0) -> x0 = r * cos(phi), y0 = r * sin(phi)
c         phi0 := bq * pr
c         cos(phi1) = cos(phi + phi0)
c                   = cos(phi) * cos(phi0) - sin(phi) * sin(phi0)
c                   = x0/r * cos(phi0) - y0/r * sin(phi0)
c                   = (x0 * cos(phi0) - y0 * sin(phi0)) / r
c         sin(phi1) = sin(phi + phi0)
c                   = sin(phi) * cos(phi0) + cos(phi) * sin(phi0)
c                   = y0/r * cos(phi0) + x0/r * sin(phi0)
c                   = (y0 * cos(phi0) + x0 * sin(phi0)) / r
c
c          x1 = r1 * cos(phi1)
c             = (cos(phi0) *  x0 - sin(phi0) *  y0) * (r1 / r)
c          y1 = r1 * sin(phi1)
c             = (cos(phi0) *  y0 + sin(phi0) *  x0) * (r1 / r)
c         px1 = (cos(phi1) *  pr - sin(phi1) * pphi) / r1
c             = (cos(phi0) * px0 - sin(phi0) * py0) / (r1 / r)
c         py1 = (sin(phi1) *  pr + cos(phi1) * pphi) / r1
c             = (cos(phi0) * py0 + sin(phi0) * px0) / (r1 / r)
c          z1 = z0 + bp * pr * pphi
c
c         w := (r1 / r) = exp(bq * pphi)
c         c := cos(phi0)
c         s := sin(phi0)
c
c       Optimized map code
        x0  = x(i)
        y0  = y(i)
        px0 = px(i)
        py0 = py(i)
        z0  = z(i)
        p   = 1.d0 + g(i)

c       b     = .5d0 * dbz / p
c       bq    = -.25d0 * b
c       bp    =  .50d0 * b
        bq    = -.125d0 * dbz / p
        bp    = -2.d0 * bq
        pr    = x0 * px0 + y0 * py0
        pphi  = x0 * py0 - y0 * px0
        phi0  = bq * pr
        w     = exp(bq * pphi)
        c     = cos(phi0)
        s     = sin(phi0)

        x(i)  = (c *  x0 - s *  y0) * w
        y(i)  = (c *  y0 + s *  x0) * w
        px(i) = (c * px0 - s * py0) / w
        py(i) = (c * py0 + s * px0) / w
        z(i)  = z0 + bp * pr * pphi
      enddo
      return
      end

      subroutine tsfrie(trans,cod,dbz)
      implicit none
      real*8 ,intent(out):: trans(6,6),cod(6)
      real*8 ,intent(in):: dbz
      real*8 x0,px0,y0,py0,z0,p,bq,bp,pr,pphi,phi0,w,c,s,
     $     x1,px1,y1,py1,z1

c     Apply identity map if dbz == 0
c      if(dbz .eq. 0.d0)then
c        call tinitr(trans)
c        return
c      endif

c     Map definition
c       [position @ enter]
c        x0  := cod(1)
c       px0  := cod(2)
c        y0  := cod(3)
c       py0  := cod(4)
c        z0  := cod(5)
c        p   := 1 + cod(6)
c
c       r    := hypot(x0, y0)
c       b    := .5 * dbz / p**2
c       bq   := -.25 * b
c       bp   :=  .5  * b / p
c       pr   := x0 * px0 + y0 * py0
c       pphi := x0 * py0 - y0 * px0
c       r1   := r * exp(bq * pphi)
c       phi  := atan(y0, x0)
c       phi1 := phi + bq * pr
c
c       [position @ exit]
c        x1  := cod(1) <- r1 * cos(phi1)
c        y1  := cod(3) <- r1 * sin(phi1)
c       px1  := cod(2) <- (x1 * pr - y1 * pphi) / r1**2
c       py1  := cod(4) <- (y1 * pr + x1 * pphi) / r1**2
c        z1  := cod(5) <- z0 + bp * pr * pphi
c
c     Map expansion
c       phi = atan2(y0, x0) -> x0 = r * cos(phi), y0 = r * sin(phi)
c       phi0 := bq * pr
c       cos(phi1) = cos(phi + phi0)
c                 = cos(phi) * cos(phi0) - sin(phi) * sin(phi0)
c                 = x0/r * cos(phi0) - y0/r * sin(phi0)
c                 = (x0 * cos(phi0) - y0 * sin(phi0)) / r
c       sin(phi1) = sin(phi + phi0)
c                 = sin(phi) * cos(phi0) + cos(phi) * sin(phi0)
c                 = y0/r * cos(phi0) + x0/r * sin(phi0)
c                 = (y0 * cos(phi0) + x0 * sin(phi0)) / r
c
c        x1 = r1 * cos(phi1)
c           = (cos(phi0) *  x0 - sin(phi0) *  y0) * (r1 / r)
c        y1 = r1 * sin(phi1)
c           = (cos(phi0) *  y0 + sin(phi0) *  x0) * (r1 / r)
c       px1 = (cos(phi1) *  pr - sin(phi1) * pphi) / r1
c           = (cos(phi0) * px0 - sin(phi0) * py0) / (r1 / r)
c       py1 = (sin(phi1) *  pr + cos(phi1) * pphi) / r1
c           = (cos(phi0) * py0 + sin(phi0) * px0) / (r1 / r)
c        z1 = z0 + bp * pr * pphi
c
c       w := (r1 / r) = exp(bq * pphi)
c       c := cos(phi0)
c       s := sin(phi0)
c
c     Derivatives
c       dbq    =  bq * (-2 / p) * dp = bp * dp
c       dbp    =  bp * (-3 / p) * dp
c       dpr    =  px0 * dx + x0 * dpx + py0 * dy + y0 * dpy
c       dpphi  =  py0 * dx - y0 * dpx - px0 * dx + x0 * dpy
c       dphi0  =  bq * dpr + pr * dbq
c              =  bq * (px0 * dx + x0 * dpx + py0 * dy + y0 * dpy)
c               + bp * pr * dp
c       dw     =  w * (bq * dpphi + pphi * bp * dp)
c       d(1/w) =  (1 / w) * (-1 / w) * dw
c              = -(1 / w) * (bq * dpphi + pphi * bp * dp)
c       dc     = -s * dphi0
c       ds     =  c * dphi0
c
c        dx1 =  ( c *  dx - s *  dy) * w + (-s *  x0 - c *  y0) * w * dphi0
c             + ( x1 / w) * dw
c            =  (c *  dx - s *  dy) * w
c             -  y1 * dphi0 +  x1 * (bq * dpphi + bp * pphi * dp)
c            =  c * w *  dx + bq * (  x1 * ( py0) -   y1 * px0) *  dx
c                           + bq * (  x1 * ( -y0) -   y1 *  x0) * dpx
c             - s * w *  dy + bq * (  x1 * (-px0) -   y1 * py0) *  dy
c                           + bq * (  x1 * (  x0) -   y1 *  y0) * dpy
c                           + bp * (  x1 *  pphi  -   y1 *  pr) *  dp
c
c       dpx1 =  ( c * dpx - s * dpy) / w + (-s * px0 - c * py0) / w * dphi0
c             + (px1 * w) * d(1/w)
c            =  (c * dpx - s * dpy) / w
c             - py1 * dphi0 - px1 * (bq * dpphi + pphi * bp * dp)
c            =                bq * (-px1 * ( py0) -  py1 * px0) *  dx
c             + c / w * dpx + bq * (-px1 * ( -y0) -  py1 *  x0) * dpx
c                           + bq * (-px1 * (-px0) -  py1 * py0) *  dy
c             - s / w * dpy + bq * (-px1 * (  x0) -  py1 *  y0) * dpy
c                           + bp * (-px1 *  pphi  -  py1 *  pr) *  dp
c
c        dy1 =  ( c *  dy + s *  dx) * w + (-s *  y0 + c *  x0) * w * dphi0
c             + ( y1 / w) * dw
c            =  ( s *  dx + c *  dy) * w
c             +  x1 * dphi0 +  y1 * (bq * dpphi + bp * pphi * dp)
c            =  s * w *  dx + bq * (  y1 * ( py0) +   x1 * px0) *  dx
c                           + bq * (  y1 * ( -y0) +   x1 *  x0) * dpx
c             + c * w *  dy + bq * (  y1 * (-px0) +   x1 * py0) *  dy
c                           + bq * (  y1 * (  x0) +   x1 *  y0) * dpy
c                           + bp * (  y1 *  pphi  +   x1 *  pr) *  dp
c
c       dpy1 =  ( c * dpy + s * dpx) / w + (-s * py0 + c * px0) / w * dphi0
c             + (py1 * w) * d(1/w)
c            =  (c * dpx - s * dpy) / w
c             + px1 * dphi0 - py1 * (bq * dpphi + pphi * bp * dp)
c            =                bq * (-py1 * ( py0) +  px1 * px0) *  dx
c             + s / w * dpx + bq * (-py1 * ( -y0) +  px1 *  x0) * dpx
c                           + bq * (-py1 * (-px0) +  px1 * py0) *  dy
c             + c / w * dpy + bq * (-py1 * (  x0) +  px1 *  y0) * dpy
c                           + bp * (-py1 *  pphi  +  px1 *  pr) *  dp
c
c       dz1 = dz + pr * pphi * dbp + bp * pr * dpphi + bp + pphi * dpr
c           =                 bp * (  pr * ( py0) + pphi * px0) *  dx
c                           + bp * (  pr * ( -y0) + pphi *  x0) * dpx
c                           + bp * (  pr * (-px0) + pphi * py0) *  dy
c                           + bp * (  pr * (  x0) + pphi *  y0) * dpy
c             +          dz
c                           - 3 * bp * pr * pphi / p            *  dp
c
c     Optimized map code
      x0  = cod(1)
      px0 = cod(2)
      y0  = cod(3)
      py0 = cod(4)
      z0  = cod(5)
      p   = 1.d0 + cod(6)

c     b     = .5d0 * dbz / p**2
c     bq    = -.25d0 * b
c     bp    =  .50d0 * b / p
      bq    = -.125d0 * dbz / p**2
      bp    = -2.d0 * bq / p
      pr    = x0 * px0 + y0 * py0
      pphi  = x0 * py0 - y0 * px0
      phi0  = bq * pr
      w     = exp(bq * pphi)
      c     = cos(phi0)
      s     = sin(phi0)

      x1  = (c *  x0 - s *  y0) * w
      y1  = (c *  y0 + s *  x0) * w
      px1 = (c * px0 - s * py0) / w
      py1 = (c * py0 + s * px0) / w
      z1  = z0 + bp * pr * pphi

      cod(1) =  x1
      cod(2) = px1
      cod(3) =  y1
      cod(4) = py1
      cod(5) =  z1

c     Transfer Matrix: trans(i,j) := dcod(i)@exit / dcod(j)@enter
      trans(1,1) =  c * w + bq * (  x1 * ( py0) -   y1 * px0)
      trans(1,2) =          bq * (  x1 * ( -y0) -   y1 *  x0)
      trans(1,3) = -s * w + bq * (  x1 * (-px0) -   y1 * py0)
      trans(1,4) =          bq * (  x1 * (  x0) -   y1 *  y0)
      trans(1,5) =  0.d0
      trans(1,6) =          bp * (  x1 *   pphi -   y1 *  pr)

      trans(2,1) =          bq * (-px1 * ( py0) -  py1 * px0)
      trans(2,2) =  c / w + bq * (-px1 * ( -y0) -  py1 *  x0)
      trans(2,3) =          bq * (-px1 * (-px0) -  py1 * py0)
      trans(2,4) = -s / w + bq * (-px1 * (  x0) -  py1 *  y0)
      trans(2,5) =  0.d0
      trans(2,6) =          bp * (-px1 *   pphi -  py1 *  pr)

      trans(3,1) =  s * w + bq * (  y1 * ( py0) +   x1 * px0)
      trans(3,2) =          bq * (  y1 * ( -y0) +   x1 *  x0)
      trans(3,3) =  c * w + bq * (  y1 * (-px0) +   x1 * py0)
      trans(3,4) =          bq * (  y1 * (  x0) +   x1 *  y0)
      trans(3,5) =  0.d0
      trans(3,6) =          bp * (  y1 *   pphi +   x1 *  pr)

      trans(4,1) =          bq * (-py1 * ( py0) +  px1 * px0)
      trans(4,2) =  s / w + bq * (-py1 * ( -y0) +  px1 *  x0)
      trans(4,3) =          bq * (-py1 * (-px0) +  px1 * py0)
      trans(4,4) =  c / w + bq * (-py1 * (  x0) +  px1 *  y0)
      trans(4,5) =  0.d0
      trans(4,6) =          bp * (-py1 *   pphi +  px1 *  pr)

      trans(5,1) =          bp * (  pr * ( py0) + pphi * px0)
      trans(5,2) =          bp * (  pr * ( -y0) + pphi *  x0)
      trans(5,3) =          bp * (  pr * (-px0) + pphi * py0)
      trans(5,4) =          bp * (  pr * (  x0) + pphi *  y0)
      trans(5,5) =  1.d0
      trans(5,6) = -3.d0 * bp * pr * pphi / p

      trans(6,1) =  0.d0
      trans(6,2) =  0.d0
      trans(6,3) =  0.d0
      trans(6,4) =  0.d0
      trans(6,5) =  0.d0
      trans(6,6) =  1.d0

      return
      end
