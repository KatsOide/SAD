      recursive subroutine tsolqu(np,x,px,y,py,z,gp,dv,bsi,al,ak,bz0,
     $     ak0x,ak0y,ibsi,eps0)
      use tsolz
      use tfstk
      implicit none
      type (tzparam) tz
      integer*4 np,i,n,ndiv,ibsi
      real*8, parameter::smax=0.99d0,phieps=1.d-7
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),gp(np),bsi(np),
     $     sx(np),sy(np),sz(np)
      real*8 al,ak,eps0,bz,a,b,c,d,akk,eps,
     $     bw,dw,r,ap,dpz,ak0x,ak0y,bz0,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,z00,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,ra,phi,pxi,pyi,
     $     awu,dwu,dz1,dz2
        associate (
     $       w1=>tz%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%phi1,phi2=>tz%phi2,
     $     wss=>tz%wss,bzp=>tz%bzp,akkp=>tz%akkp,aln=>tz%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,
     $     c1=>tz%c1,s1=>tz%s1,dc1=>tz%dc1,xs1=>tz%xs1,ch2=>tz%ch2,
     $       sh2=>tz%sh2,dch2=>tz%dch2,xsh2=>tz%xsh2,pr=>tz%pr,
     $     c1p=>tz%c1p,s1p=>tz%s1p,xs1p=>tz%xs1p,ch2p=>tz%ch2p,
     $       sh2p=>tz%sh2p,xsh2p=>tz%xsh2p,
     $     w1p=>tz%w1p,w2p=>tz%w2p,wsp=>tz%wsp,w12p=>tz%w12p,
     $       wdp=>tz%wdp,phi1p=>tz%phi1p,phi2p=>tz%phi2p,
     $       g1=>tz%g1,g2=>tz%g2,g1p=>tz%g1p,g2p=>tz%g2p,
     $       wr1=>tz%wr1,wr2=>tz%wr2,wr1p=>tz%wr1p,wr2p=>tz%wr2p,
     $     wssip=>tz%wssip, ca1p=>tz%ca1p,
     $       dcw1p=>tz%dcw1p, dcw2p=>tz%dcw2p, 
     $       csw1p=>tz%csw1p,cswsp=>tz%cswsp,
     $     cr2p=>tz%cr2p, cr3p=>tz%cr3p,dxs=>tz%dxs,dxsp=>tz%dxsp,
     $     aw1=>tz%aw1,aw2=>tz%aw2,aw1p=>tz%aw1p,aw2p=>tz%aw2p,
     $       cxs1=>tz%cxs1,cxs2=>tz%cxs2,
     $       cxs1p=>tz%cxs1p,cxs2p=>tz%cxs2p)
      if(ak*al .lt. 0.d0)then
        call tsolqu(np,y,py,x,px,z,gp,dv,bsi,al,-ak,
     $       -bz0,-ak0y,-ak0x,ibsi,eps0)
        return
      endif
      bz=bz0
      if(ak .eq. 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,bsi,
     $       al,bz,ak0x,ak0y,.false.)
        return
      endif
      z00=z(1)
      if(eps0 .eq. 0.d0)then
        eps=0.2d0
      else
        eps=0.2d0*eps0
      endif
      ndiv=1+int(abs(al*dcmplx(ak,bz))/eps)
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz .eq. 0.d0)then
        do i=1,np
c     pr=(1.d0+gp(i))**2
          call tzsetparam0(tz,gp(i),akk)
          ra=aln*0.5d0
          if(ibsi .eq. 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
          else
            bsi(i)=0.d0
          endif
          do n=1,ndiv
            ap=min(smax,px(i)**2+py(i)**2)
            dpz=sqrt1(-ap)
c             dpz=-ap/(1.d0+sqrt(1.d0-ap))
            r=-dpz/(1.d0+dpz)*ra
            ra=aln
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
          enddo
          ap=min(smax,px(i)**2+py(i)**2)
          dpz=sqrt1(-ap)
c          dpz=-ap/(1.d0+sqrt(1.d0-ap))
          r=-dpz/(1.d0+dpz)*aln*.5d0
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          if(ibsi .eq. 2)then
            bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)
          endif
        enddo
      else
        do i=1,np
          call tzsetparam(tz,gp(i),akk,bz)
          if(ibsi .eq. 1)then
            bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*al
          else
            bsi(i)=bzp*al
          endif
          px(i)=px(i)+bzp*y(i)*.5d0
          py(i)=py(i)-bzp*x(i)*.5d0
          ra=aln*0.5d0
          do n=1,ndiv
            ap=min(smax,px(i)**2+py(i)**2)
            dpz=sqrt1(-ap)
c            dpz=-ap/(1.d0+sqrt(1.d0-ap))
            r=-dpz/(1.d0+dpz)*ra
            ra=aln
            phi=r*bzp
            a24=sin(phi)
            a12=a24/bzp
            a22=cos(phi)
            if(a22 .ge. 0.d0)then
              a14=a12*a24/(1.d0+a22)
            else
              a14=(1.d0-a22)/bzp
            endif
            pxi=px(i)
            x(i) =x(i)+a12*pxi+a14*py(i)
            y(i) =y(i)-a14*pxi+a12*py(i)
            px(i)=     a22*pxi+a24*py(i)
            py(i)=    -a24*pxi+a22*py(i)
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            pxi=px(i)
            pyi=py(i)
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
            call tztaf(0,tz,awu,pxi,pyi,aw1,ws,w12,wss,g1,
     $           0.d0,0.d0,0.d0,0.d0,0.d0,dz1,0.d0)
            call tztaf(0,tz,-dwu,-pyi,pxi,aw2,-w12,ws,-wss,g2,
     $           0.d0,0.d0,0.d0,0.d0,0.d0,dz2,0.d0)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln*dv(i)
          enddo
          ap=min(smax,px(i)**2+py(i)**2)
          dpz=sqrt1(-ap)
c          dpz=-ap/(1.d0+sqrt(1.d0-ap))
          r=-dpz/(1.d0+dpz)*aln*.5d0
          phi=r*bzp
          a24=sin(phi)
          a12=a24/bzp
          a22=cos(phi)
          if(a22 .ge. 0.d0)then
            a14=a12*a24/(1.d0+a22)
          else
            a14=(1.d0-a22)/bzp
          endif
          pxi=px(i)
          x(i) =x(i)+a12*pxi+a14*py(i)
          y(i) =y(i)-a14*pxi+a12*py(i)
          px(i)=     a22*pxi+a24*py(i)
          py(i)=    -a24*pxi+a22*py(i)
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

      recursive subroutine tsolqur(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,al,ak,
     $     bz0,ak0x,ak0y,eps0,alr)
      use tsolz
      use tfstk
      use tspin
      implicit none
      type (tzparam) tz
      integer*4 np,i,n,ndiv
      real*8 smax
      parameter (smax=0.99d0)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),gp(np),
     $     px0(np),py0(np),zr0(np),bsi(np)
      real*8 sx(np),sy(np),sz(np)
      real*8 al,ak,eps0,bz,a,b,c,d,akk,eps,alr,
     $     bw,dw,r,ap,dpz,ak0x,ak0y,bz0,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     dx0,dy0,xi,yi,a12,a14,a22,a24,phi,pxi,pyi,
     $     phieps,awu,dwu,dz1,dz2
      parameter (phieps=1.d-7)
        associate (
     $       w1=>tz%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%phi1,phi2=>tz%phi2,
     $     wss=>tz%wss,bzp=>tz%bzp,akkp=>tz%akkp,aln=>tz%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,
     $     c1=>tz%c1,s1=>tz%s1,dc1=>tz%dc1,xs1=>tz%xs1,ch2=>tz%ch2,
     $       sh2=>tz%sh2,dch2=>tz%dch2,xsh2=>tz%xsh2,pr=>tz%pr,
     $     c1p=>tz%c1p,s1p=>tz%s1p,xs1p=>tz%xs1p,ch2p=>tz%ch2p,
     $       sh2p=>tz%sh2p,xsh2p=>tz%xsh2p,
     $     w1p=>tz%w1p,w2p=>tz%w2p,wsp=>tz%wsp,w12p=>tz%w12p,
     $       wdp=>tz%wdp,phi1p=>tz%phi1p,phi2p=>tz%phi2p,
     $       g1=>tz%g1,g2=>tz%g2,g1p=>tz%g1p,g2p=>tz%g2p,
     $       wr1=>tz%wr1,wr2=>tz%wr2,wr1p=>tz%wr1p,wr2p=>tz%wr2p,
     $     wssip=>tz%wssip, ca1p=>tz%ca1p,
     $       dcw1p=>tz%dcw1p, dcw2p=>tz%dcw2p, 
     $       csw1p=>tz%csw1p,cswsp=>tz%cswsp,
     $     cr2p=>tz%cr2p, cr3p=>tz%cr3p,dxs=>tz%dxs,dxsp=>tz%dxsp,
     $     aw1=>tz%aw1,aw2=>tz%aw2,aw1p=>tz%aw1p,aw2p=>tz%aw2p,
     $       cxs1=>tz%cxs1,cxs2=>tz%cxs2,
     $       cxs1p=>tz%cxs1p,cxs2p=>tz%cxs2p)
      if(ak*al .lt. 0.d0)then
        call tsolqur(np,y,py,x,px,z,gp,dv,sy,sx,sz,
     $       py0,px0,zr0,bsi,al,-ak,
     $       -bz0,-ak0y,-ak0x,eps0,alr)
        return
      endif
      bz=bz0
      if(ak .eq. 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,sx,sy,sz,bsi,
     $       al,bz,ak0x,ak0y,.false.)
        alr=al
        return
      endif
      if(eps0 .eq. 0.d0)then
        eps=0.2d0
      else
        eps=0.2d0*eps0
      endif
      ndiv=1+int(abs(al*dcmplx(ak,bz))/eps)
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      if(bz .eq. 0.d0)then
        alr=aln*0.5d0
        do n=1,ndiv
          do i=1,np
c     pr=(1.d0+gp(i))**2
            call tzsetparam0(tz,gp(i),akk)
            if(n .eq. 1)then
              bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)
            else
              bsi(i)=0.d0
            endif
            ap=min(smax,px(i)**2+py(i)**2)
            dpz=sqrt1(-ap)
c             dpz=-ap/(1.d0+sqrt(1.d0-ap))
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
          alr=aln
          if(n .ne. ndiv)then
            call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $           px0,py0,zr0,bsi,alr)
            px0=px
            py0=py
            zr0=z
          endif
        enddo
        do i=1,np
          ap=min(smax,px(i)**2+py(i)**2)
          dpz=sqrt1(-ap)
c          dpz=-ap/(1.d0+sqrt(1.d0-ap))
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
            call tzsetparam(tz,gp(i),akk,bz)
            if(n .eq. 1)then
              bsi(i)=akk*(x(i)+dx0)*(y(i)+dy0)+bzp*alr
            else
              bsi(i)=bzp*alr
            endif
            ap=min(smax,px(i)**2+py(i)**2)
            dpz=sqrt1(-ap)
c            dpz=-ap/(1.d0+sqrt(1.d0-ap))
            r=-dpz/(1.d0+dpz)*alr
            phi=r*bzp
            a24=sin(phi)
            a12=a24/bzp
            a22=cos(phi)
            if(a22 .ge. 0.d0)then
              a14=a12*a24/(1.d0+a22)
            else
              a14=(1.d0-a22)/bzp
            endif
            pxi=px(i)
            x(i) =x(i)+a12*pxi+a14*py(i)
            y(i) =y(i)-a14*pxi+a12*py(i)
            px(i)=     a22*pxi+a24*py(i)
            py(i)=    -a24*pxi+a22*py(i)
            z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
            pxi=px(i)
            pyi=py(i)
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
            call tztaf(0,tz,awu,pxi,pyi,aw1,ws,w12,wss,g1,
     $           0.d0,0.d0,0.d0,0.d0,0.d0,dz1,0.d0)
            call tztaf(0,tz,-dwu,-pyi,pxi,aw2,-w12,ws,-wss,g2,
     $           0.d0,0.d0,0.d0,0.d0,0.d0,dz2,0.d0)
            z(i)=z(i)+
     $           bzp*(-((awu*dwu*dxs**2)/akkp) +
     $           ca1*pxi*pyi*wss)
     $           +dz1+dz2-aln*dv(i)
            if(n .eq. ndiv)then
              bsi(i)=-akk*(x(i)+dx0)*(y(i)+dy0)+bzp*aln*.5d0
            endif
          enddo
          alr=aln
          if(n .ne. ndiv)then
            call tradk(np,x,px,y,py,z,gp,dv,sx,sy,sz,
     $           px0,py0,zr0,bsi,alr)
            px0=px
            py0=py
            zr0=z
          endif
        enddo
        do i=1,np
          ap=min(smax,px(i)**2+py(i)**2)
          dpz=sqrt1(-ap)
c     dpz=-ap/(1.d0+sqrt(1.d0-ap))
          r=-dpz/(1.d0+dpz)*aln*.5d0
          phi=r*bzp
          a24=sin(phi)
          a12=a24/bzp
          a22=cos(phi)
          if(a22 .ge. 0.d0)then
            a14=a12*a24/(1.d0+a22)
          else
            a14=(1.d0-a22)/bzp
          endif
          pxi=px(i)
          x(i) =x(i)+a12*pxi+a14*py(i)
          y(i) =y(i)-a14*pxi+a12*py(i)
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
