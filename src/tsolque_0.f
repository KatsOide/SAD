      module tsolz
      type tzparam
      real*8 w1,w2,ws,w12,wd,phi1,phi2,
     $     wss,bzp,akkp,aln,csws,ca1,cr1,dcw1,dcw2,cr2,cr3,
     $     c1,s1,dc1,xs1,ch2,sh2,dch2,xsh2,pr,
     $     c1p,s1p,xs1p,ch2p,sh2p,xsh2p,
     $     w1p,w2p,wsp,w12p,wdp,
     $     wssip, cr1p, ca1p, dcw1p, dcw2p, cswsp,
     $     cr2p, cr3p,dxs,dxsp,
     $     aw1,aw2,aw1p,aw2p,cxs1,cxs2,cxs1p,cxs2p
      end type

      contains
        subroutine tztaf(mode, tz,
     $     zwu, pxt, pyt, zw1, wst, w12t, swss, g1,
     $     zwup, zw1p, wstp, w12tp, g1p, dz, dzp)
        implicit none
        type (tzparam) tz
        real*8 zwu, pxt, pyt, zw1, wst, w12t, swss, g1,
     $     zwup, zw1p, g1p, dz, dzp, f1,f2,f3,f1p,f2p,f3p,
     $       w12xsh,wsxs,xwsd,wsm,ws1,ws2,cr1,
     $       w12xshp,wsxsp,xwsdp,wsmp,ws1p,ws2p,cr1p,
     $       wstp,w12tp,cswa,cswap,akkw
        integer*4 mode
        wst = w1t + w2t
        w12t =  w1t - w2t
        akkw=w1t*w2t
        w12xsh = w12t*xsh2
        wsxs = wst*xs1
        cwsd = c1*w12xsh - ch2*wsxs
        cswa = (csw1 + cwsd)*swss
        wsm = 4*swss*w1*w2
        ws1 = (-1 + wsm)*wst
        ws2 = w12t*(1 + wsm)
        cr1 = bzp**2*swss
        f1 = cswa + csws - zw1
        f2 = ca1*cr1 + g1
        f3 = 2*cr1*cwsd + cr3*w12xsh +
     $       aln*(dcw1*ws1 + dcw2*ws2) - cr2*wsxs
        dz=zwu*(zwu*f3/2.+(f2*pxt + tz%bzp*f1*pyt)) - (zw1*pyt**2)/2.
        if(mode .ne. 0)then
          w12xshp = w12tp*xsh2 + w12t*xsh2p
          wsxsp = wstp*xs1 + wst*xs1p
          cwsdp = c1p*w12xsh + c1*w12xshp - ch2p*wsxs - ch2*wsxsp
          cswap = swss*(cwsdp + (csw1 + cwsd)*wssip)
          wsmp = 4*swss*(w1p*w2 + w1*(w2p + w2*wssip))
          ws1p = wsmp*wst + (-1 + wsm)*wstp
          ws2p = w12tp + w12tp*wsm + w12t*wsmp
          cr1p = (bzp**2*swss*(-2 + pr*wssip))/pr
          f1p = cswap + cswsp - zw1p
          f2p = ca1p*cr1 + ca1*cr1p + g1p
          f3p = 2*cr1p*cwsd + 2*cr1*cwsdp + cr3p*w12xsh + cr3*w12xshp +
     $         aln*(dcw1p*ws1 + dcw1*ws1p + dcw2p*ws2 + dcw2*ws2p) -
     $         cr2p*wsxs - cr2*wsxsp
          dzp=(zwu**2*f3p)/2. - (zw1p*pyt**2)/2. + 
     -         zwup*(zwu*f3 + f2*pxt + tz%bzp*f1*pyt) + 
     -         zwu*(f2p*pxt + tz%bzp*(f1p - f1/tz%pr)*pyt)
        endif
        return
        end

        subroutine tzsetparam(tz,dp,akk,bz)
        implicit none
        type (tzparam) tz
        real*8 dp,bz,akk,xsin,xsinh,wa
        tz%pr=1.d0+dp
        tz%bzp=bz/tz%pr
        tz%akkp=akk/tz%pr
        wa=sqrt(tz%bzp**4+4.d0*tz%akkp**2)
        tz%w1=sqrt((tz%bzp**2+wa)*.5d0)
        tz%w2=tz%akkp/tz%w1
        tz%wss=1.d0/(tz%w1**2+tz%w2**2)
        tz%ws=tz%w1+tz%w2
        tz%wd=tz%bzp/tz%ws
        tz%phi1=tz%aln*tz%w1
        tz%c1=cos(tz%phi1)
        tz%s1=sin(tz%phi1)
        tz%xs1=xsin(tz%phi1)
        if(tz%c1 .ge. 0.d0)then
          tz%dc1=-tz%s1**2/(1.d0+tz%c1)
        else
          tz%dc1=tz%c1-1.d0
        endif
        tz%phi2=tz%aln*tz%w2
        tz%ch2=cosh(tz%phi2)
        tz%sh2=sinh(tz%phi2)
        tz%xsh2=xsinh(tz%phi2)
        tz%dch2=tz%sh2**2/(1.d0+tz%ch2)
        tz%cxs1=tz%dc1*tz%s1 - tz%xs1
        tz%cxs2=tz%dch2*tz%sh2 - tz%xsh2
        tz%aw1=tz%aln + tz%w1*tz%cxs2/(tz%w2*tz%ws)
        tz%aw2=tz%aln + tz%w2*tz%cxs1/(tz%w1*tz%ws)
        tz%ca1=-(tz%ch2*tz%dc1) - tz%dch2 -
     $       (tz%s1*tz%sh2*tz%w12)/tz%ws
        tz%cr1=tz%bzp**2*tz%wss
        tz%cr2=tz%c1*tz%w2/tz%w1
        tz%cr3=tz%ch2*tz%w1/tz%w2
        tz%dxs=tz%w2*tz%xs1-tz%w1*tz%xsh2
        tz%csws=(tz%ch2*tz%w1*tz%phi1 + tz%c1*tz%w2*tz%phi2)*tz%wss
        tz%dcw1=tz%dc1*tz%w2
        tz%dcw2=tz%dch2*tz%w1
        return
        end

        subroutine tzsetparam0(tz,dp,akk)
        implicit none
        type (tzparam) tz
        real*8 dp,akk,xsin,xsinh
        tz%pr=(1.d0+dp)
        tz%akkp=akk/tz%pr
        tz%w1=sqrt(tz%akkp)
        tz%phi1=tz%aln*tz%w1
        tz%c1=cos(tz%phi1)
        tz%s1=sin(tz%phi1)
        tz%xs1=xsin(tz%phi1)
        if(tz%c1 .ge. 0.d0)then
          tz%dc1=-tz%s1**2/(1.d0+tz%c1)
        else
          tz%dc1=tz%c1-1.d0
        endif
        tz%phi2=tz%aln*tz%w2
        tz%ch2=cosh(tz%phi1)
        tz%sh2=sinh(tz%phi1)
        tz%xsh2=xsinh(tz%phi1)
        tz%dch2=tz%sh2**2/(1.d0+tz%ch2)
        return
        end
        
        subroutine tzsetparamp(tz)
        implicit none
        type (tzparam) tz
        tz%w1p=-tz%w1**3/tz%pr*tz%wss
        tz%w2p=-tz%w2**3/tz%pr*tz%wss
        tz%wsp=tz%w1p+tz%w2p
        tz%wdp=-tz%wd*(1.d0/tz%pr+tz%wsp/tz%ws)
        tz%wssip=2.d0*(tz%w1**4+tz%w2**4)*tz%wss**2/tz%pr
        tz%w12p=tz%w1p-tz%w2p
        tz%c1p=-tz%s1*tz%w1p*tz%aln
        tz%s1p=tz%c1*tz%w1p*tz%aln
        tz%xs1p=-tz%dc1*tz%w1p*tz%aln
        tz%ch2p=tz%sh2*tz%w2p*tz%aln
        tz%sh2p=tz%ch2*tz%w2p*tz%aln
        tz%xsh2p=-tz%dch2*tz%w2p*tz%aln        
        tz%cxs1p=tz%aln*(tz%dc1 + tz%c1*tz%dc1 - tz%s1**2)*tz%w1p
        tz%cxs2p=tz%aln*(tz%dch2 + tz%ch2*tz%dch2 + tz%sh2**2)*tz%w2p
        tz%aw1p=(tz%cxs2*(tz%w1p*tz%w2*tz%ws - tz%w1*tz%w2p*tz%ws) +
     $       tz%akkp*(tz%cxs2p*tz%ws - tz%cxs2*tz%wsp))
     $       /(tz%w2**2*tz%ws**2)
        tz%aw2p=(-(tz%cxs1*tz%w1p*tz%w2*tz%ws) +
     $       tz%cxs1*tz%w1*tz%w2p*tz%ws + 
     -       tz%akkp*(tz%cxs1p*tz%ws - tz%cxs1*tz%wsp))
     $       /(tz%w1**2*tz%ws**2)
        tz%ca1p=tz%aln*(tz%ch2*tz%s1*tz%w1p - tz%c1*tz%sh2*tz%w2p) - 
     -       (tz%s1*tz%sh2*tz%w12p + tz%aln*tz%c1*tz%sh2*tz%w12*tz%w1p
     $       + tz%aln*tz%ch2*tz%s1*tz%w12*tz%w2p)/tz%ws +
     $       (tz%s1*tz%sh2*tz%w12*tz%wsp)/tz%ws**2
        tz%cr1p=tz%bzp**2*tz%wss*(-2./tz%pr + tz%wssip)
        tz%cr2p=(-((tz%c1 + tz%s1*tz%phi1)*tz%w1p*tz%w2) +
     $       tz%c1*tz%w1*tz%w2p)/tz%w1**2
        tz%cr3p=(tz%ch2*tz%w1p*tz%w2 - tz%ch2*tz%w1*tz%w2p +
     $       tz%sh2*tz%phi1*tz%w2*tz%w2p)/tz%w2**2
        tz%dxsp=tz%w2p*(tz%phi1*tz%dch2 + tz%xs1) -
     $       tz%w1p*(tz%phi2*tz%dc1 + tz%xsh2)
        tz%cswsp= tz%aln*tz%wss*((-(tz%s1*tz%w1p*tz%w2*tz%phi2) +
     $       tz%sh2*tz%w1*tz%phi1*tz%w2p) + 
     -       tz%ch2*tz%w1*(2.d0*tz%w1p + tz%w1*tz%wssip) +
     $       tz%c1*tz%w2*(2.d0*tz%w2p + tz%w2*tz%wssip))
        tz%dcw1p=-(tz%phi2*tz%s1*tz%w1p) + tz%dc1*tz%w2p
        tz%dcw2p= tz%dch2*tz%w1p + tz%phi1*tz%sh2*tz%w2p
        return
        end

      end module

      subroutine tsolque(trans,cod,beam,al,ak,
     $     bz0,ak0x,ak0y,eps0,enarad,radcod,calpol,irad,ld)
      use tsolz
      implicit none
      type(tzparam) tz
      integer*4 i,n,ndiv,ld,itgetqraddiv,irad
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6)
      real*8 al,ak,eps0,bz,a,b,c,d,akk,eps,bzh,
     $     bw,dw,ak0x,ak0y,dx0,dy0,
     $     xi0,yi0,dy,dpy,
     $     adx,adpy,adp,bdy,bdpx,bdp,bwdy,bwdpx,bwdp,
     $     cdx,cdpy,cdp,ddy,ddpx,ddp,dwdy,dwdpx,dwdp,cwdp,
     $     u1wx,u1wpx,u1wy,u1wpy,
     $     u2wx,u2wpx,u2wy,u2wpy,
     $     v1wx,v1wpx,v1wy,v1wpy,
     $     v2wx,v2wpx,v2wy,v2wpy,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     u1p,u1wp,u2p,u2wp,v1p,v1wp,v2p,v2wp,
     $     dv,dvdp,xi,yi,pxi,pyi,xf,yf,pxf,pyf,
     $     tbrho,bx,by,bxy,b1,br,bz0,cw,phieps,
     $     awu,dwu,awup,dwup,
     $     dz1,dz2,dz1p,dz2p
      logical*4 enarad,calpol,radcod
      parameter (phieps=1.d-2)
      if(ak*(1.d0+cod(6)) .lt. 0.d0)then
c        write(*,'(a,1p8g13.5)')'tsolque-in  ',ak,bz,cod
        call texchg(trans,cod,beam)
        call tsolquerc(trans,cod,beam,al,-ak,
     $       -bz0,-ak0y,-ak0x,eps0,enarad,radcod,calpol,irad,ld)
        call texchg(trans,cod,beam)
c        write(*,'(a,1p8g13.5)')'tsolque-out ',ak,bz,cod
        return
      endif
      bz=bz0
      if(ak .eq. 0.d0)then
        call tdrife(trans,cod,beam,al,
     $       bz,ak0x,ak0y,.true.,enarad,.false.,irad,ld)
        return
      endif
      if(eps0 .eq. 0.d0)then
        eps=0.1d0
      else
        eps=0.1d0*eps0
      endif
      ndiv=1+int(sqrt((ak*al)**2+(bz*al)**2)/eps)
      if(enarad)then
        ndiv=max(ndiv,itgetqraddiv(cod,ak,al))
      endif
      tz%aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      br=tbrho()
      b1=br*akk
      call tinitr(trans1)
c     end   initialize for tz%preventing compiler tz%warning
      do n=1,ndiv
        if(enarad .and. radcod .or. n .eq. 1)then
          call tzsetparam(tz,cod(6),akk,bz)
          call tgetdv(cod(6),dv,dvdp)
          call tzsetparamp(tz)
        endif
        if(enarad)then
          bx= b1*(cod(3)+dy0)
          by= b1*(cod(1)+dx0)
          bxy=b1
          if(n .eq. 1)then
            call trade(trans,beam,cod,bx,by,bz*br,bz,
     $           0.d0,bxy,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,.5d0*tz%aln)
          else
            call trade(trans,beam,cod,bx,by,bz*br,bz,
     $           0.d0,bxy,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,tz%aln)
          endif
        endif
        if(n .eq. 1)then
          call tqente(trans,cod,beam,tz%aln*.5d0,bz,calpol,irad,ld)
        else
          call tqente(trans,cod,beam,tz%aln,bz,calpol,irad,ld)
        endif
        xi0=cod(1)
        yi0=cod(3)
        xi=xi0+dx0
        yi=yi0+dy0
        bzh=bz*.5d0
        pxi=(cod(2)+yi0*bzh)/tz%pr
        pyi=(cod(4)-xi0*bzh)/tz%pr
        a = (tz%w2*tz%ws*xi-tz%bzp*pyi)*tz%wss
        bw= (tz%ws*pxi-tz%bzp*tz%w2*yi)*tz%wss
        b = bw*tz%w1
        c = (tz%w1*tz%wd*xi+pyi)*tz%wss
        cw=c*tz%w2
        dw= (-tz%wd*pxi+tz%w1*yi)*tz%wss
        d = dw*tz%w2
        u1w= a*tz%dc1+bw*tz%s1
        u1=u1w*tz%w1
        u2w=-a*tz%s1+bw*tz%dc1
        u2=u2w*tz%w1
        v1w= c*tz%dch2+dw*tz%sh2
        v1=v1w*tz%w2
        v2w= c*tz%sh2+dw*tz%dch2
        v2=v2w*tz%w2
        xf =xi0+u1w+v1w*tz%bzp
        pxf=pxi+u2 +v2*tz%bzp
        dy =tz%wd*u2w+tz%ws*v2w
        yf =yi0+dy
        dpy=-tz%wd*u1 +tz%ws*v1
        pyf=pyi+dpy
        adx=tz%w2*tz%ws*tz%wss
        adpy=-tz%bzp*tz%wss
        adp=a*tz%wssip+
     $       ((tz%w2p*tz%ws+tz%w2*tz%wsp)*xi+tz%bzp*pyi/tz%pr)*tz%wss
        bwdpx=tz%ws*tz%wss
        bwdy=-tz%bzp*tz%w2*tz%wss
        bwdp=bw*tz%wssip+
     $       (tz%wsp*pxi-tz%bzp*(tz%w2p-tz%w2/tz%pr)*yi)*tz%wss
        bdpx=bwdpx*tz%w1
        bdy=bwdy*tz%w1
        bdp=bwdp*tz%w1+bw*tz%w1p
        cdx=tz%w1*tz%wd*tz%wss
        cdpy=tz%wss
        cdp=c*tz%wssip+((tz%w1p*tz%wd+tz%w1*tz%wdp)*xi)*tz%wss
        cwdp=cdp*tz%w2+cw*tz%w2p
        dwdpx=-tz%wd*tz%wss
        dwdy=tz%w1*tz%wss
        dwdp=dw*tz%wssip+(-tz%wdp*pxi+tz%w1p*yi)*tz%wss
c        write(*,'(a,1p7g15.7)')': ',dw*tz%wssip,-tz%wdp*pxi*tz%wss,tz%w1p*tz%wss*yi,
c     $       cdp*tz%dch2*tz%bzp,c*tz%ch2p*tz%bzp,dwdp*tz%sh2*tz%bzp,dw*tz%sh2p*tz%bzp
        ddpx=dwdpx*tz%w2
        ddy=dwdy*tz%w2
        ddp=dwdp*tz%w2+dw*tz%w2p
        u1wx=adx*tz%dc1
        u1wpx=bwdpx*tz%s1
        u1wy=bwdy*tz%s1
        u1wpy=adpy*tz%dc1
        u2wx=-adx*tz%s1
        u2wpx=bwdpx*tz%dc1
        u2wy=bwdy*tz%dc1
        u2wpy=-adpy*tz%s1
        v1wx =cdx*tz%dch2
        v1wpx=dwdpx*tz%sh2
        v1wy =dwdy*tz%sh2
        v1wpy=cdpy*tz%dch2
        v2wx =cdx*tz%sh2
        v2wpx=dwdpx*tz%dch2
        v2wy =dwdy*tz%dch2
        v2wpy=cdpy*tz%sh2
        u1wp=adp*tz%dc1+a*tz%c1p+bwdp*tz%s1+bw*tz%s1p
        u1p=u1wp*tz%w1+u1w*tz%w1p
        u2wp=-adp*tz%s1-a*tz%s1p+bwdp*tz%dc1+bw*tz%c1p
        u2p=u2wp*tz%w1+u2w*tz%w1p
        v1wp=cdp*tz%dch2+c*tz%ch2p+dwdp*tz%sh2+dw*tz%sh2p
        v1p=v1wp*tz%w2+v1w*tz%w2p
        v2wp=cdp*tz%sh2+c*tz%sh2p+dwdp*tz%dch2+dw*tz%ch2p
        v2p=v2wp*tz%w2+v2w*tz%w2p
        trans1(1,1)=1.d0+u1wx+v1wx*tz%bzp
        trans1(1,2)=u1wpx+v1wpx*tz%bzp
        trans1(1,3)=u1wy+v1wy*tz%bzp
        trans1(1,4)=u1wpy+v1wpy*tz%bzp
        trans1(1,6)=u1wp+(v1wp-v1w/tz%pr)*tz%bzp
        trans1(2,1)=(tz%w1*u2wx+tz%w2*v2wx*tz%bzp)
        trans1(2,2)=1.d0+(tz%w1*u2wpx+tz%w2*v2wpx*tz%bzp)
        trans1(2,3)=(tz%w1*u2wy+tz%w2*v2wy*tz%bzp)
        trans1(2,4)=(tz%w1*u2wpy+tz%w2*v2wpy*tz%bzp)
        trans1(2,6)=u2p+(v2p-v2/tz%pr)*tz%bzp
        trans1(3,1)=tz%wd*u2wx+tz%ws*v2wx
        trans1(3,2)=tz%wd*u2wpx+tz%ws*v2wpx
        trans1(3,3)=1.d0+tz%wd*u2wy+tz%ws*v2wy
        trans1(3,4)=tz%wd*u2wpy+tz%ws*v2wpy
        trans1(3,6)=tz%wdp*u2w+tz%wd*u2wp+tz%wsp*v2w+tz%ws*v2wp
        trans1(4,1)=-tz%wd*tz%w1*u1wx+tz%ws*tz%w2*v1wx
        trans1(4,2)=-tz%wd*tz%w1*u1wpx+tz%ws*tz%w2*v1wpx
        trans1(4,3)=-tz%wd*tz%w1*u1wy+tz%ws*tz%w2*v1wy
        trans1(4,4)=1.d0-tz%wd*tz%w1*u1wpy+tz%ws*tz%w2*v1wpy
        trans1(4,6)=-tz%wdp*u1-tz%wd*u1p+tz%wsp*v1+tz%ws*v1p
        tz%w12=tz%w1-tz%w2
        awu=a/tz%ws*tz%w1
        dwu=d
        awup=adp/tz%ws*tz%w1-a/tz%ws*(tz%w2p-tz%wsp/tz%ws*tz%w2)
        dwup=ddp
        call tztaf(1,tz,awu,pxi,pyi,tz%aw1,tz%w1,tz%w2,1.d0,
     $       tz%s1**2,
     $       awup,tz%aw1p,tz%w1p,tz%w2p,2.d0*tz%s1*tz%s1p,dz1,dz1p)
        call tztaf(1,tz,dwu,-pyi,pxi,tz%aw2,tz%w2,-tz%w1,-1.d0,
     $       -tz%sh2**2,
     $       dwup,tz%aw2p,tz%w2p,-tz%w1p,-2.d0*tz%sh2*tz%sh2p,dz2,dz2p)
        cod(5)=cod(5)+
     $       tz%bzp*(-((awu*dwu*tz%dxs**2)/tz%akkp) +
     $       tz%ca1*pxi*pyi*tz%wss)
     $       +dz1+dz2-tz%aln*dv
        trans1(5,6)=
     $       -((tz%bzp*tz%dxs*(awup*dwu*tz%dxs +
     $       awu*dwup*tz%dxs + 2*awu*dwu*tz%dxsp))/
     -     tz%akkp) + tz%bzp*pxi*pyi*tz%wss*(tz%ca1p +
     $       tz%ca1*(-1./tz%pr + tz%wssip))
     $       +dz1p+dz2p+dvdp*tz%aln
        cod(1)=xf
        cod(2)=pxf*tz%pr-bzh*yf
        cod(3)=yf
        cod(4)=pyf*tz%pr+bzh*xf
        do i=1,4
          trans1(i,2)=trans1(i,2)/tz%pr
          trans1(i,4)=trans1(i,4)/tz%pr
          trans1(i,1)=trans1(i,1)-bzh*trans1(i,4)
          trans1(i,3)=trans1(i,3)+bzh*trans1(i,2)
          trans1(i,6)=trans1(i,6)
     $         -(pxi*trans1(i,2)+pyi*trans1(i,4))
        enddo
        do i=1,6
          trans1(2,i)=trans1(2,i)*tz%pr-bzh*trans1(3,i)
          trans1(4,i)=trans1(4,i)*tz%pr+bzh*trans1(1,i)
        enddo
        trans1(2,6)=trans1(2,6)+pxf
        trans1(4,6)=trans1(4,6)+pyf
        trans1(5,1)=
     $       -trans1(1,1)*trans1(2,6)
     $       +trans1(2,1)*trans1(1,6)
     $       -trans1(3,1)*trans1(4,6)
     $       +trans1(4,1)*trans1(3,6)
        trans1(5,2)=
     $       -trans1(1,2)*trans1(2,6)
     $       +trans1(2,2)*trans1(1,6)
     $       -trans1(3,2)*trans1(4,6)
     $       +trans1(4,2)*trans1(3,6)
        trans1(5,3)=
     $       -trans1(1,3)*trans1(2,6)
     $       +trans1(2,3)*trans1(1,6)
     $       -trans1(3,3)*trans1(4,6)
     $       +trans1(4,3)*trans1(3,6)
        trans1(5,4)=
     $       -trans1(1,4)*trans1(2,6)
     $       +trans1(2,4)*trans1(1,6)
     $       -trans1(3,4)*trans1(4,6)
     $       +trans1(4,4)*trans1(3,6)
        trans1(5,6)=trans1(5,6)
     $       -(pxi*trans1(5,2)+pyi*trans1(5,4))
c        write(*,'(a/,6(1p6g11.4/))')
c     $       'tsolque ',((trans1(i,j),j=1,6),i=1,6)
        call tmultr5(trans,trans1,irad)
        if(irad .gt. 6 .or. calpol)then
          call tmulbs(beam ,trans1,.false.,.true.)
        endif
      enddo
      call tqente(trans,cod,beam,tz%aln*.5d0,bz,calpol,irad,ld)
      if(enarad)then
        bx= b1*cod(3)
        by= b1*cod(1)
        bxy= b1
        call trade(trans,beam,cod,bx,by,bz*br,bz,
     $       0.d0,bxy,0.d0,0.d0,
     $       0.d0,0.d0,0.d0,0.d0,.5d0*tz%aln)
      endif
      return
      end

      subroutine tsolquerc(trans,cod,beam,al,ak,
     $     bz,ak0x,ak0y,eps0,enarad,radcod,calpol,irad,ld)
      implicit none
      integer*4 ld,irad
      real*8 trans(6,12),cod(6),beam(42)
      real*8 al,ak,eps0,bz,ak0x,ak0y
      logical*4 enarad,radcod,calpol
      call tsolque(trans,cod,beam,al,ak,
     $     bz,ak0x,ak0y,eps0,enarad,radcod,calpol,irad,ld)
      return
      end

      integer*4 function itgetirad()
      include 'inc/TMACRO.inc'
      itgetirad=irad
      return
      end

      subroutine texchg(trans,cod,beam)
      implicit none
      integer*4 i
      real*8 trans(6,12),cod(6),beam(42),x0,px0,x
      x0=cod(1)
      cod(1)=cod(3)
      cod(3)=x0
      px0=cod(2)
      cod(2)=cod(4)
      cod(4)=px0
      do i=1,12
        x0=trans(1,i)
        trans(1,i)=trans(3,i)
        trans(3,i)=x0
        px0=trans(2,i)
        trans(2,i)=trans(4,i)
        trans(4,i)=px0
      enddo
      x=beam(1)
      beam(1)=beam(6)
      beam(6)=x
      x=beam(2)
      beam(2)=beam(9)
      beam(9)=x
      x=beam(3)
      beam(3)=beam(10)
      beam(10)=x
      x=beam(5)
      beam(5)=beam(7)
      beam(7)=x
      x=beam(11)
      beam(11)=beam(13)
      beam(13)=x
      x=beam(12)
      beam(12)=beam(14)
      beam(14)=x
      x=beam(16)
      beam(16)=beam(18)
      beam(18)=x
      x=beam(17)
      beam(17)=beam(19)
      beam(19)=x

      x=beam(21+1)
      beam(21+1)=beam(21+6)
      beam(21+6)=x
      x=beam(21+2)
      beam(21+2)=beam(21+9)
      beam(21+9)=x
      x=beam(21+3)
      beam(21+3)=beam(21+10)
      beam(21+10)=x
      x=beam(21+5)
      beam(21+5)=beam(21+7)
      beam(21+7)=x
      x=beam(21+11)
      beam(21+11)=beam(21+13)
      beam(21+13)=x
      x=beam(21+12)
      beam(21+12)=beam(21+14)
      beam(21+14)=x
      x=beam(21+16)
      beam(21+16)=beam(21+18)
      beam(21+18)=x
      x=beam(21+17)
      beam(21+17)=beam(21+19)
      beam(21+19)=x
      return
      end

      integer*4 function itgetqraddiv(cod,ak,al)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 nrad
      real*8 cod(6),ak,xd,xpd,a,b,al
      xd=max(1.d-6,abs(cod(1))+abs(cod(3)))
      xpd=max(1.d-6,abs(cod(2))+abs(cod(4)))
      a=min(1.d-2,abs(ak)*xd+xpd)
      b=p0*amass/c*a/al
      nrad=int(al*crad/epsrad*(h0*b)**2)
      itgetqraddiv=max(int(emidiv*emidiq*nrad),
     1     int(abs(a)/epsrad/1.d3*emidiv*emidiq))
      return
      end

      real*8 function tbrho()
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      tbrho=brho
      return
      end
