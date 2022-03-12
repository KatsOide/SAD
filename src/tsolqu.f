      subroutine tsolqu(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $     al,ak,bz0,ak0x,ak0y,ibsi,eps0)
      use tsolz
      use kradlib, only:bsi
      use mathfun
      implicit none
      type (tzparams) tz
      integer*4 ,intent(in):: np,ibsi
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np),
     $     gp(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak,eps0,bz0,ak0x,ak0y
      real*8, parameter::phieps=1.d-7,epsdef=0.2d0
      integer*4 i,n,ndiv
      real*8 bz,a,c,akk,eps,bw,dw,r,ap,dpz,aln0,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,pxi,pyi,
     $     awu,dwu,dz1,dz2
        associate (
     $     s1=>tz%tz0%s1,dc1=>tz%tz0%dc1,xs1=>tz%tz0%xs1,
     $     sh2=>tz%tz0%sh2,dch2=>tz%tz0%dch2,xsh2=>tz%tz0%xsh2,
     $     akkp=>tz%tz0%akkp,bzp=>tz%tz0%bzp,
     $     w1=>tz%tz0%w1,phi1=>tz%tz0%phi1,
     $     w2=>tz%tz1%w2,ws=>tz%tz1%ws,w12=>tz%tz1%w12,wd=>tz%tz1%wd,
     $     phi2=>tz%tz1%phi2,wss=>tz%tz1%wss,csw1=>tz%tz1%csw1,
     $     csws=>tz%tz1%csws,ca1=>tz%tz1%ca1,dcw1=>tz%tz1%dcw1,
     $     dcw2=>tz%tz1%dcw2,cr2=>tz%tz1%cr2,cr3=>tz%tz1%cr3,
     $     g1=>tz%tz1%g1,g2=>tz%tz1%g2,cxs1=>tz%tz1%cxs1,
     $     cxs2=>tz%tz1%cxs2,aw1=>tz%tz1%aw1,aw2=>tz%tz1%aw2,
     $     wr1=>tz%tz1%wr1,wr2=>tz%tz1%wr2,dxs=>tz%tz1%dxs)
      if(ak*al < 0.d0)then
        write(*,*)'tsolqu-implementation error ',al,ak,bz0
        stop
      elseif(ak == 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $       al,bz0,ak0x,ak0y,.false.)
        return
      endif
      bz=bz0
      eps=merge(epsdef,epsdef*eps0,eps0 == 0.d0)
      ndiv=1+int(abs(al*hypot(ak,bz)/eps))
c      ndiv=1+int(abs(al*dcmplx(ak,bz))/eps)
      aln0=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz == 0.d0)then
        do concurrent (i=1:np)
          call tzsetparam0(tz%tz0,gp(i),aln0,akk)
          if(ibsi == 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
          elseif(ibsi >= 0)then
            bsi(i)=0.d0
          endif
          do n=1,ndiv
            pxi=px(i)
            pyi=py(i)
            ap=pxi**2+pyi**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*merge(aln0*.5d0,aln0,n == 1)
            x(i)=x(i)+pxi*r
            y(i)=y(i)+pyi*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            u1 =   xi*dc1 +pxi*s1/w1
            u2 =-xi*w1*s1 +pxi*dc1
            v1 =   yi*dch2+pyi*sh2/w1
            v2 = yi*w1*sh2+pyi*dch2
            x(i) =x(i)+u1
            px(i)=pxi +u2 
            y(i) =y(i)+v1
            py(i)=pyi +v2 
            z(i) =z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
     $           -0.25d0*(
     $           w1*(xi**2*xs1-yi**2*xsh2)
     $           +(ap-dv(i))*aln0
     $           +u1*px(i)+xi*pxi*dc1
     $           +v1*py(i)+yi*pyi*dch2)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln0*.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi == 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
          endif
        enddo
      else
        do concurrent (i=1:np)
          call tzsetparams(tz,gp(i),aln0,akk,bz)
          if(ibsi == 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          elseif(ibsi >= 0)then
            bsi(i)=bzp*al
          endif
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          do n=1,ndiv
            ap=px(i)**2+py(i)**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*merge(aln0*.5d0,aln0,n == 1)
            call xsincos(r*bzp,a24,a12,a22,a14)
            pxi=px(i)
            x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
            y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
            pxi  =      a22*pxi+a24*py(i)
            pyi  =     -a24*pxi+a22*py(i)
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            a=  (w2*ws*xi-bzp*pyi  )*wss
            bw= (  ws*pxi-bzp*w2*yi)*wss
            c=  (w1*wd*xi+pyi  )*wss
            dw= ( -wd*pxi+w1*yi)*wss
            u1w= a*dc1 +bw*s1
            u2w=-a*s1  +bw*dc1
            v1w= c*dch2+dw*sh2
            v2w= c*sh2 +dw*dch2
            x(i) =x(i)+u1w      +v1w*bzp
            px(i)=pxi +u2w*w1   +v2w*w2*bzp
            y(i) =y(i)+wd*u2w   +ws*v2w
            py(i)=pyi -wd*u1w*w1+ws*v1w*w2
            awu=a/ws*w1
            dwu=dw*w2
            call tztafs(tz, awu, pxi,pyi,aw1, ws, w12,wss, g1, dz1)
            call tztafs(tz,-dwu,-pyi,pxi,aw2,-w12,ws, -wss,g2, dz2)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln0*dv(i)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln0*.5d0
          call xsincos(r*bzp,a24,a12,a22,a14)
          pxi=px(i)
          x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
          y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
          px(i)=      a22*pxi+a24*py(i)-bzp*y(i)*.5d0
          py(i)=     -a24*pxi+a22*py(i)+bzp*x(i)*.5d0
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi == 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          endif
        enddo
      endif
      return
      end associate
      end

      subroutine tsolqum(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $     al,ak,bz,ak0x,ak0y,ibsi,eps0,tzs,ini)
      use tsolz
      use kradlib, only:bsi
c      use tmacro,only:l_track
      use mathfun
      implicit none
      integer*4 ,intent(in):: np,ibsi
      type (tzparams) ,intent(inout):: tzs(np)
c      type (tzparams)  tzs(np),tz
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np),
     $     gp(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak,eps0,bz,ak0x,ak0y
      logical*4 ,intent(in):: ini
      real*8, parameter::phieps=1.d-7,epsdef=0.2d0
      integer*4 i,n,ndiv
      real*8 a,c,akk,eps,bw,dw,r,ap,dpz,aln0,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,pxi,pyi,
     $     awu,dwu,dz1,dz2
      if(ak*al < 0.d0)then
        write(*,*)'tsolqu-implementation error ',al,ak,bz
        stop
      elseif(ak == 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $       al,bz,ak0x,ak0y,.false.)
        return
      endif
      eps=merge(epsdef,epsdef*eps0,eps0 == 0.d0)
      ndiv=1+int(abs(al*hypot(ak,bz)/eps))
c      ndiv=1+int(abs(al*dcmplx(ak,bz))/eps)
      aln0=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz == 0.d0)then
        if(ini)then
          do concurrent (i=1:np)
            call tzsetparam0(tzs(i)%tz0,gp(i),aln0,akk)
          enddo
        endif
        do concurrent (i=1:np)
          associate (
     $         w1=>tzs(i)%tz0%w1,s1=>tzs(i)%tz0%s1,dc1=>tzs(i)%tz0%dc1,
     $         xs1=>tzs(i)%tz0%xs1,sh2=>tzs(i)%tz0%sh2,
     $         dch2=>tzs(i)%tz0%dch2,xsh2=>tzs(i)%tz0%xsh2)
          if(ibsi == 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
          elseif(ibsi >= 0)then
            bsi(i)=0.d0
          endif
          do n=1,ndiv
            pxi=px(i)
            pyi=py(i)
            ap=pxi**2+pyi**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*merge(aln0*.5d0,aln0,n == 1)
            x(i)=x(i)+pxi*r
            y(i)=y(i)+pyi*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            u1 =   xi*dc1 +pxi*s1/w1
            u2 =-xi*w1*s1 +pxi*dc1
            v1 =   yi*dch2+pyi*sh2/w1
            v2 = yi*w1*sh2+pyi*dch2
            x(i) =x(i)+u1
            px(i)=pxi +u2 
            y(i) =y(i)+v1
            py(i)=pyi +v2 
            z(i) =z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
     $           -0.25d0*(
     $           w1*(xi**2*xs1-yi**2*xsh2)
     $           +(ap-dv(i))*aln0
     $           +u1*px(i)+xi*pxi*dc1
     $           +v1*py(i)+yi*pyi*dch2)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln0*.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi == 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
          endif
          end associate
        enddo
      else
        if(ini)then
          do concurrent (i=1:np)
            call tzsetparams(tzs(i),gp(i),aln0,akk,bz)
          enddo
        endif
        do concurrent (i=1:np)
          associate (
     $         w1=>tzs(i)%tz0%w1,bzp=>tzs(i)%tz0%bzp,
     $         akkp=>tzs(i)%tz0%akkp,pr=>tzs(i)%tz0%pr,
     $         s1=>tzs(i)%tz0%s1,dc1=>tzs(i)%tz0%dc1,
     $         xs1=>tzs(i)%tz0%xs1,sh2=>tzs(i)%tz0%sh2,
     $         dch2=>tzs(i)%tz0%dch2,xsh2=>tzs(i)%tz0%xsh2,
     $         w2=>tzs(i)%tz1%w2,ws=>tzs(i)%tz1%ws,w12=>tzs(i)%tz1%w12,
     $         wd=>tzs(i)%tz1%wd,wss=>tzs(i)%tz1%wss,
     $         ca1=>tzs(i)%tz1%ca1,dxs=>tzs(i)%tz1%dxs,
     $         g1=>tzs(i)%tz1%g1,g2=>tzs(i)%tz1%g2,
     $         aw1=>tzs(i)%tz1%aw1,aw2=>tzs(i)%tz1%aw2)
          if(ibsi == 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          elseif(ibsi >= 0)then
            bsi(i)=bzp*al
          endif
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          do n=1,ndiv
            ap=px(i)**2+py(i)**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*merge(aln0*.5d0,aln0,n == 1)
            call xsincos(r*bzp,a24,a12,a22,a14)
            pxi=px(i)
            x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
            y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
            pxi  =      a22*pxi+a24*py(i)
            pyi  =     -a24*pxi+a22*py(i)
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            a=  (w2*ws*xi-bzp*pyi  )*wss
            bw= (  ws*pxi-bzp*w2*yi)*wss
            c=  (w1*wd*xi+pyi  )*wss
            dw= ( -wd*pxi+w1*yi)*wss
            u1w= a*dc1 +bw*s1
            u2w=-a*s1  +bw*dc1
            v1w= c*dch2+dw*sh2
            v2w= c*sh2 +dw*dch2
            x(i) =x(i)+u1w      +v1w*bzp
            px(i)=pxi +u2w*w1   +v2w*w2*bzp
            y(i) =y(i)+wd*u2w   +ws*v2w
            py(i)=pyi -wd*u1w*w1+ws*v1w*w2
            awu=a/ws*w1
            dwu=dw*w2
            call tztafs(tzs(i), awu, pxi,pyi,aw1, ws, w12,wss, g1, dz1)
            call tztafs(tzs(i),-dwu,-pyi,pxi,aw2,-w12,ws, -wss,g2, dz2)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln0*dv(i)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln0*.5d0
          call xsincos(r*bzp,a24,a12,a22,a14)
          pxi=px(i)
          x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
          y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
          px(i)=      a22*pxi+a24*py(i)-bzp*y(i)*.5d0
          py(i)=     -a24*pxi+a22*py(i)+bzp*x(i)*.5d0
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi == 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          endif
          end associate
        enddo
      endif
      return
      end

      subroutine tsolqur(np,x,px,y,py,z,gp,dv,sx,sy,sz,al,ak,
     $     bz0,ak0x,ak0y,eps0,alr)
      use tsolz
      use kradlib, only:bsi,tradk
      use ffs_flag, only:ndivrad
      use photontable,only:tgswap,pcvt
      use mathfun
      implicit none
      type (tzparams) tz
      integer*4 ,intent(in):: np
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np),
     $     gp(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak,eps0,bz0,ak0x,ak0y
      real*8 ,intent(out):: alr
      real*8, parameter::phieps=1.d-7,epsdef=0.2d0,arad=0.01d0
      integer*4 ,parameter ::ndivmax=1000
      integer*4 i,n,ndiv
      real*8 bz,a,c,akk,eps,bw,dw,r,ap,dpz,aka,b,d,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,pxi,pyi,
     $     awu,dwu,dz1,dz2
        associate (
     $       c1=>tz%tz0%c1,s1=>tz%tz0%s1,dc1=>tz%tz0%dc1,
     $       xs1=>tz%tz0%xs1,ch2=>tz%tz0%ch2,sh2=>tz%tz0%sh2,
     $       dch2=>tz%tz0%dch2,xsh2=>tz%tz0%xsh2,pr=>tz%tz0%pr,
     $       akkp=>tz%tz0%akkp,aln=>tz%tz0%aln,bzp=>tz%tz0%bzp,
     $       w1=>tz%tz0%w1,phi1=>tz%tz0%phi1,
     $       w2=>tz%tz1%w2,ws=>tz%tz1%ws,w12=>tz%tz1%w12,wd=>tz%tz1%wd,
     $       phi2=>tz%tz1%phi2,wss=>tz%tz1%wss,csw1=>tz%tz1%csw1,
     $       csws=>tz%tz1%csws,ca1=>tz%tz1%ca1,dcw1=>tz%tz1%dcw1,
     $       dcw2=>tz%tz1%dcw2,cr2=>tz%tz1%cr2,cr3=>tz%tz1%cr3,
     $       g1=>tz%tz1%g1,g2=>tz%tz1%g2,cxs1=>tz%tz1%cxs1,
     $       cxs2=>tz%tz1%cxs2,aw1=>tz%tz1%aw1,aw2=>tz%tz1%aw2,
     $       wr1=>tz%tz1%wr1,wr2=>tz%tz1%wr2,dxs=>tz%tz1%dxs)
      if(ak*al < 0.d0)then
        write(*,*)'tsolqur-implementation error ',al,ak,bz0
        stop
      elseif(ak == 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $       al,bz0,ak0x,ak0y,.false.)
        alr=al
        return
      endif
      bz=bz0
      eps=merge(epsdef,epsdef*eps0,eps0 == 0.d0)
      aka=hypot(ak,bz)
      ndiv=1+int(abs(al)*aka/eps)
      ndiv=min(ndivmax,
     $     max(ndiv,ndivrad(hypot(ak0x,ak0y),ak,bz,eps0)))
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz == 0.d0)then
        alr=aln*0.5d0
        do n=1,ndiv
c!$OMP PARALLEL
c!$OMP DO
          do i=1,np
            call tzsetparam0(tz%tz0,gp(i),aln,akk)
            bsi(i)=merge(akk*(x(i)+dx0)*(y(i)+dy0),0.d0,n == 1)
            ap=px(i)**2+py(i)**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*alr
            x(i)=x(i)+px(i)*r
            y(i)=y(i)+py(i)*r
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            b=px(i)
            d=py(i)
            u1 =   xi*dc1 +b*s1/w1
            u2 =-xi*w1*s1 +b*dc1
            v1 =   yi*dch2+d*sh2/w1
            v2 = yi*w1*sh2+d*dch2
            x(i) =x(i) +u1
            px(i)=px(i)+u2 
            y(i) =y(i) +v1
            py(i)=py(i)+v2 
            z(i) =z(i)-0.25d0*(
     $           w1*(xi**2*xs1-yi**2*xsh2)
     $           +(b**2+d**2)*aln
     $           +u1*(u2+b)+xi*b*dc1
     $           +v1*(v2+d)+yi*d*dch2)
     $           -dv(i)*aln
            if(n == ndiv)then
              bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
            endif
          enddo
c!$OMP END DO
c!$OMP END PARALLEL
          alr=aln
          call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,alr,0.d0)
          pcvt%fr0=pcvt%fr0+alr/al
        enddo
        do i=1,np
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*0.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
        enddo
      else
        do i=1,np
          bzp=bz/(1.d0+gp(i))
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
        enddo
        alr=aln*0.5d0
        do n=1,ndiv
          do i=1,np
            call tzsetparams(tz,gp(i),aln,akk,bz)
            bsi(i)=merge(akk*(x(i)+dx0)*(y(i)+dy0)+bzp*alr,
     $           bzp*alr,n == 1)
            ap=px(i)**2+py(i)**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*alr
            call xsincos(r*bzp,a24,a12,a22,a14)
            pxi=px(i)
            x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
            y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
            pxi  =     a22*pxi+a24*py(i)
            pyi  =    -a24*pxi+a22*py(i)
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            xi=x(i)+dx0
            yi=y(i)+dy0
            a=  (w2*ws*xi-bzp*pyi)*wss
            bw= (ws*pxi-bzp*w2*yi)*wss
            c=  (w1*wd*xi +pyi)*wss
            dw= (-wd*pxi+w1*yi)*wss
            u1w= a*dc1+bw*s1
            u1=u1w*w1
            u2w=-a*s1 +bw*dc1
            u2=u2w*w1
            v1w= c*dch2+dw*sh2
            v1=v1w*w2
            v2w= c*sh2 +dw*dch2
            v2=v2w*w2
            x(i) =x(i)+u1w+v1w*bzp
            px(i)=pxi +u2 + v2*bzp
            y(i) =y(i)+wd*u2w+ws*v2w
            py(i)=pyi -wd*u1 +ws*v1
            awu=a/ws*w1
            dwu= dw*w2
            call tztafs(tz,awu,  pxi,pyi,aw1,  ws,w12, wss,g1,dz1)
            call tztafs(tz,-dwu,-pyi,pxi,aw2,-w12,ws, -wss,g2,dz2)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln*dv(i)
            if(n == ndiv)then
              bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*aln*.5d0
            endif
          enddo
          alr=aln
          call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,alr,0.d0)
          pcvt%fr0=pcvt%fr0+alr/al
        enddo
        do i=1,np
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*.5d0
          call xsincos(r*bzp,a24,a12,a22,a14)
          pxi=px(i)
          x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
          y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
          px(i)=      a22*pxi+a24*py(i)-bzp*y(i)*.5d0
          py(i)=     -a24*pxi+a22*py(i)+bzp*x(i)*.5d0
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
        enddo
      endif
      return
      end associate
      end
