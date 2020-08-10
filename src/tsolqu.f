      subroutine tsolqu(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $     al,ak,bz0,ak0x,ak0y,ibsi,eps0)
      use tsolz
      use kradlib, only:bsi
      use mathfun
      use tmacro, only:l_track
      implicit none
      type (tzparam) tz
      integer*4 ,intent(in):: np,ibsi
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),dv(np),
     $     gp(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak,eps0,bz0,ak0x,ak0y
      real*8, parameter::phieps=1.d-7,epsdef=0.2d0
      integer*4 i,n,ndiv
      real*8 bz,a,c,akk,eps,bw,dw,r,ap,dpz,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,ra,pxi,pyi,
     $     awu,dwu,dz1,dz2
        associate (
     $       w1=>tz%tz0%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%tz0%phi1,phi2=>tz%phi2,wss=>tz%wss,
     $       bzp=>tz%tz0%bzp,akkp=>tz%tz0%akkp,aln=>tz%tz0%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,g1=>tz%g1,g2=>tz%g2,
     $       wr1=>tz%wr1,wr2=>tz%wr2,aw1=>tz%aw1,aw2=>tz%aw2,
     $     c1=>tz%tz0%c1,s1=>tz%tz0%s1,dc1=>tz%tz0%dc1,xs1=>tz%tz0%xs1,
     $       ch2=>tz%tz0%ch2,sh2=>tz%tz0%sh2,dch2=>tz%tz0%dch2,
     $       xsh2=>tz%tz0%xsh2,pr=>tz%tz0%pr,
     $     c1p=>tz%tzp%c1p,s1p=>tz%tzp%s1p,xs1p=>tz%tzp%xs1p,
     $       ch2p=>tz%tzp%ch2p,sh2p=>tz%tzp%sh2p,xsh2p=>tz%tzp%xsh2p,
     $     w1p=>tz%tzp%w1p,w2p=>tz%tzp%w2p,wsp=>tz%tzp%wsp,
     $       w12p=>tz%tzp%w12p,wdp=>tz%tzp%wdp,phi1p=>tz%tzp%phi1p,
     $       phi2p=>tz%tzp%phi2p,g1p=>tz%tzp%g1p,g2p=>tz%tzp%g2p,
     $       wr1p=>tz%tzp%wr1p,wr2p=>tz%tzp%wr2p,
     $     wssip=>tz%tzp%wssip, ca1p=>tz%tzp%ca1p,
     $       dcw1p=>tz%tzp%dcw1p, dcw2p=>tz%tzp%dcw2p,
     $       csw1p=>tz%tzp%csw1p,cswsp=>tz%tzp%cswsp,dxs=>tz%dxs,
     $     cr2p=>tz%tzp%cr2p, cr3p=>tz%tzp%cr3p,
     $       dxsp=>tz%tzp%dxsp,
     $       aw1p=>tz%tzp%aw1p,aw2p=>tz%tzp%aw2p,cxs1=>tz%cxs1,
     $       cxs2=>tz%cxs2,cxs1p=>tz%tzp%cxs1p,cxs2p=>tz%tzp%cxs2p)
      if(ak*al .lt. 0.d0)then
        write(*,*)'tsolqu-implementation error ',al,ak,bz0
        stop
      elseif(ak .eq. 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $       al,bz0,ak0x,ak0y,.false.)
        return
      endif
      bz=bz0
      if(eps0 .eq. 0.d0)then
        eps=epsdef
      else
        eps=epsdef*eps0
      endif
      ndiv=1+int(abs(al*hypot(ak,bz)/eps))
c      ndiv=1+int(abs(al*dcmplx(ak,bz))/eps)
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz .eq. 0.d0)then
        do concurrent (i=1:np)
          tz%tz0=tzsetparam0(gp(i),aln,akk)
          ra=aln*0.5d0
          if(ibsi .eq. 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
          elseif(ibsi .ge. 0)then
            bsi(i)=0.d0
          endif
          do n=1,ndiv
            pxi=px(i)
            pyi=py(i)
            ap=pxi**2+pyi**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*ra
            ra=aln
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
     $           +(ap-dv(i))*aln
     $           +u1*px(i)+xi*pxi*dc1
     $           +v1*py(i)+yi*pyi*dch2)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi .eq. 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
          endif
        enddo
      else
        do concurrent (i=1:np)
          tz=tzsetparam(gp(i),aln,akk,bz)
          if(ibsi .eq. 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          elseif(ibsi .ge. 0)then
            bsi(i)=bzp*al
          endif
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          ra=aln*0.5d0
          do n=1,ndiv
            ap=px(i)**2+py(i)**2
            dpz=sqrt1(-ap)
            r=-dpz/(1.d0+dpz)*ra
            ra=aln
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
            u1 =u1w*w1
            u2w=-a*s1  +bw*dc1
            u2 =u2w*w1
            v1w= c*dch2+dw*sh2
            v1 =v1w*w2
            v2w= c*sh2 +dw*dch2
            v2=v2w*w2
            x(i) =x(i)+u1w   +v1w*bzp
            px(i)=pxi +u2    + v2*bzp
            y(i) =y(i)+wd*u2w+ws*v2w
            py(i)=pyi -wd*u1 +ws*v1
            awu=a/ws*w1
            dwu=dw*w2
            call tztaf(tz, awu, pxi,pyi,aw1, ws, w12,wss, g1, dz1)
            call tztaf(tz,-dwu,-pyi,pxi,aw2,-w12,ws, -wss,g2, dz2)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln*dv(i)
          enddo
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*.5d0
          call xsincos(r*bzp,a24,a12,a22,a14)
          pxi=px(i)
          x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
          y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
          px(i)=      a22*pxi+a24*py(i)
          py(i)=     -a24*pxi+a22*py(i)
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          px(i)=px(i)-bzp*y(i)*.5d0
          py(i)=py(i)+bzp*x(i)*.5d0
          if(ibsi .eq. 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          endif
        enddo
      endif
      return
      end associate
      end

      subroutine tsolqur(np,x,px,y,py,z,gp,dv,sx,sy,sz,al,ak,
     $     bz0,ak0x,ak0y,eps0,alr)
      use tsolz
      use kradlib, only:bsi,tradk
      use ffs_flag, only:ndivrad
      use photontable,only:tgswap,pcvt
      use mathfun
      implicit none
      type (tzparam) tz
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
     $       w1=>tz%tz0%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%tz0%phi1,phi2=>tz%phi2,wss=>tz%wss,
     $       bzp=>tz%tz0%bzp,akkp=>tz%tz0%akkp,aln=>tz%tz0%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,g1=>tz%g1,g2=>tz%g2,
     $       wr1=>tz%wr1,wr2=>tz%wr2,dxs=>tz%dxs,
     $     c1=>tz%tz0%c1,s1=>tz%tz0%s1,dc1=>tz%tz0%dc1,xs1=>tz%tz0%xs1,
     $       ch2=>tz%tz0%ch2,sh2=>tz%tz0%sh2,dch2=>tz%tz0%dch2,
     $       xsh2=>tz%tz0%xsh2,pr=>tz%tz0%pr,
     $     c1p=>tz%tzp%c1p,s1p=>tz%tzp%s1p,xs1p=>tz%tzp%xs1p,
     $       ch2p=>tz%tzp%ch2p,sh2p=>tz%tzp%sh2p,xsh2p=>tz%tzp%xsh2p,
     $     w1p=>tz%tzp%w1p,w2p=>tz%tzp%w2p,wsp=>tz%tzp%wsp,
     $       w12p=>tz%tzp%w12p,wdp=>tz%tzp%wdp,phi1p=>tz%tzp%phi1p,
     $       phi2p=>tz%tzp%phi2p,g1p=>tz%tzp%g1p,g2p=>tz%tzp%g2p,
     $       wr1p=>tz%tzp%wr1p,wr2p=>tz%tzp%wr2p,
     $     wssip=>tz%tzp%wssip, ca1p=>tz%tzp%ca1p,
     $       dcw1p=>tz%tzp%dcw1p, dcw2p=>tz%tzp%dcw2p,
     $       csw1p=>tz%tzp%csw1p,cswsp=>tz%tzp%cswsp,
     $     cr2p=>tz%tzp%cr2p, cr3p=>tz%tzp%cr3p,
     $       dxsp=>tz%tzp%dxsp,aw1=>tz%aw1,aw2=>tz%aw2,
     $       aw1p=>tz%tzp%aw1p,aw2p=>tz%tzp%aw2p,cxs1=>tz%cxs1,
     $       cxs2=>tz%cxs2,cxs1p=>tz%tzp%cxs1p,cxs2p=>tz%tzp%cxs2p)
      if(ak*al .lt. 0.d0)then
        write(*,*)'tsolqur-implementation error ',al,ak,bz0
        stop
      elseif(ak .eq. 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $       al,bz0,ak0x,ak0y,.false.)
        alr=al
        return
      endif
      bz=bz0
      if(eps0 .eq. 0.d0)then
        eps=epsdef
      else
        eps=epsdef*eps0
      endif
      aka=hypot(ak,bz)
      ndiv=1+int(abs(al)*aka/eps)
      ndiv=min(ndivmax,
     $     max(ndiv,ndivrad(hypot(ak0x,ak0y),ak,bz,eps0)))
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz .eq. 0.d0)then
        alr=aln*0.5d0
        do n=1,ndiv
c!$OMP PARALLEL
c!$OMP DO
          do concurrent (i=1:np)
            tz%tz0=tzsetparam0(gp(i),aln,akk)
            if(n .eq. 1)then
              bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
            else
              bsi(i)=0.d0
            endif
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
            if(n .eq. ndiv)then
              bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
            endif
          enddo
c!$OMP END DO
c!$OMP END PARALLEL
          alr=aln
          call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,alr,0.d0)
          pcvt%fr0=pcvt%fr0+alr/al
        enddo
        do concurrent (i=1:np)
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*0.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
        enddo
      else
        do concurrent (i=1:np)
          bzp=bz/(1.d0+gp(i))
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
        enddo
        alr=aln*0.5d0
        do n=1,ndiv
          do concurrent (i=1:np)
            tz=tzsetparam(gp(i),aln,akk,bz)
            if(n .eq. 1)then
              bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*alr
            else
              bsi(i)=bzp*alr
            endif
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
            b=bw*w1
            c=  (w1*wd*xi +pyi)*wss
            dw= (-wd*pxi+w1*yi)*wss
            d=dw*w2
            u1w= a*dc1+bw*s1
            u1=u1w*w1
            u2w=-a*s1 +bw*dc1
            u2=u2w*w1
            v1w= c*dch2+dw*sh2
            v1=v1w*w2
            v2w= c*sh2 +dw*dch2
            v2=v2w*w2
            x(i) =x(i) +u1w+v1w*bzp
            px(i)=pxi+u2 + v2*bzp
            y(i) =y(i) +wd*u2w+ws*v2w
            py(i)=pyi-wd*u1 +ws*v1
            awu=a/ws*w1
            dwu=d
            call tztaf(tz,awu,  pxi,pyi,aw1,  ws,w12, wss,g1,dz1)
            call tztaf(tz,-dwu,-pyi,pxi,aw2,-w12,ws, -wss,g2,dz2)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln*dv(i)
            if(n .eq. ndiv)then
              bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*aln*.5d0
            endif
          enddo
          alr=aln
          call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,alr,0.d0)
          pcvt%fr0=pcvt%fr0+alr/al
        enddo
        do concurrent (i=1:np)
          ap=px(i)**2+py(i)**2
          dpz=sqrt1(-ap)
          r=-dpz/(1.d0+dpz)*aln*.5d0
          call xsincos(r*bzp,a24,a12,a22,a14)
          pxi=px(i)
          x(i) =x(i)+(a24*pxi-a14*py(i))/bzp
          y(i) =y(i)+(a14*pxi+a24*py(i))/bzp
          px(i)=     a22*pxi+a24*py(i)
          py(i)=    -a24*pxi+a22*py(i)
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          px(i)=px(i)-bzp*y(i)*.5d0
          py(i)=py(i)+bzp*x(i)*.5d0
        enddo
      endif
      return
      end associate
      end
