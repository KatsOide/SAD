      module tsolz
      type tzparam
      real*8 w1,w2,ws,w12,wd,phi1,phi2,
     $     wss,bzp,akkp,aln,csw1,csws,ca1,dcw1,dcw2,cr2,cr3,
     $     c1,s1,dc1,xs1,ch2,sh2,dch2,xsh2,pr,
     $     c1p,s1p,xs1p,ch2p,sh2p,xsh2p,g1,g2,g1p,g2p,
     $     wr1,wr2,wr1p,wr2p,
     $     w1p,w2p,wsp,w12p,wdp,phi1p,phi2p,
     $     wssip, ca1p, dcw1p, dcw2p, csw1p,cswsp,
     $     cr2p, cr3p,dxs,dxsp,
     $     aw1,aw2,aw1p,aw2p,cxs1,cxs2,cxs1p,cxs2p
      end type

      contains
        subroutine tztaf(mode, tz,
     $       zwu, pxt, pyt, zw1, wst, w12t, swss, g1t,
     $       zwup, zw1p, wstp, w12tp, g1tp, dz, dzp)
        implicit none
        type (tzparam) tz
        real*8 zwu, pxt, pyt, zw1, wst, w12t, swss, g1t,
     $       zwup, zw1p, g1tp, dz, dzp, f1,f2,f3,f1p,f2p,f3p,
     $       w12xsh,wsxs,wsm,ws1,ws2,cr1,cwsd,
     $       w12xshp,wsxsp,wsmp,ws1p,ws2p,cr1p,cwsdp,
     $       wstp,w12tp,cswa,cswap
        integer*4 mode
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
          w12xsh = w12t*xsh2
          wsxs = wst*xs1
          cwsd = c1*w12xsh - ch2*wsxs
          cswa = (csw1 + cwsd)*swss
          wsm = 4*swss*akkp
          ws1 = (-1. + wsm)*wst
          ws2 = w12t*(1. + wsm)
          cr1 = bzp**2*swss
          f1 = cswa + csws - zw1
          f2 = ca1*cr1 + g1t
          f3 = 2.*cr1*cwsd + cr3*w12xsh +
     $         dcw1*ws1 + dcw2*ws2 - cr2*wsxs
          dz=zwu*(zwu*f3/2.+(f2*pxt + bzp*f1*pyt)) - (zw1*pyt**2)/2.
          if(mode .ne. 0)then
            w12xshp = w12tp*xsh2 + w12t*xsh2p
            wsxsp = wstp*xs1 + wst*xs1p
            cwsdp = c1p*w12xsh + c1*w12xshp - ch2p*wsxs - ch2*wsxsp
            cswap = swss*(csw1p + cwsdp) + cswa*wssip
            wsmp = wsm*(-1./pr + wssip)
            ws1p = wsmp*wst + (-1. + wsm)*wstp
            ws2p = w12tp*(1. + wsm) + w12t*wsmp
            cr1p = cr1*(-2./pr + wssip)
            f1p = cswap + cswsp - zw1p
            f2p = ca1p*cr1 + ca1*cr1p + g1tp
            f3p = 2*cr1p*cwsd + 2*cr1*cwsdp +
     $           cr3p*w12xsh + cr3*w12xshp +
     $           dcw1p*ws1 + dcw1*ws1p + dcw2p*ws2 + dcw2*ws2p -
     $           cr2p*wsxs - cr2*wsxsp
            dzp=(zwu**2*f3p)/2. - (zw1p*pyt**2)/2. + 
     -           zwup*(zwu*f3 + f2*pxt + tz%bzp*f1*pyt) + 
     -           zwu*(f2p*pxt + tz%bzp*(f1p - f1/tz%pr)*pyt)
          endif
          return
        end associate
        end subroutine

        subroutine tzsetparam(tz,dp,akk,bz)
        implicit none
        type (tzparam) tz
        real*8 dp,bz,akk,xsin,xsinh,wa
        external xsin,xsinh
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

        pr=1.d0+dp
        bzp=bz/pr
        akkp=akk/pr
        wa=sqrt(bzp**4+4.d0*akkp**2)
        w1=sqrt((bzp**2+wa)*.5d0)
        w2=akkp/w1
        wss=1.d0/(w1**2+w2**2)
        ws=w1+w2
        w12=w1-w2
        wd=bzp/ws
        phi1=aln*w1
        c1=cos(phi1)
        s1=sin(phi1)
        xs1=xsin(phi1)
        if(c1 .ge. 0.d0)then
          dc1=-s1**2/(1.d0+c1)
        else
          dc1=c1-1.d0
        endif
        phi2=aln*w2
        ch2=cosh(phi2)
        sh2=sinh(phi2)
        xsh2=xsinh(phi2)
        dch2=sh2**2/(1.d0+ch2)
        g1 = (s1**2*w2**2)/akkp
        g2 = -((sh2**2*w1**2)/akkp)
        wr1 = w1/w2
        wr2 = w2/w1
        cxs1 = dc1*s1 - xs1
        cxs2 = dch2*sh2 - xsh2
        aw1 = aln + (cxs2*wr1)/ws
        aw2 = aln + (cxs1*wr2)/ws
        cr2 = c1*wr2
        cr3 = ch2*wr1
        csw1 = akkp*aln*(-dc1 + dch2)
        csws = (ch2*w1*phi1 + c1*w2*phi2)*wss
        ca1=-(ch2*dc1) - dch2 - (s1*sh2*w12)/ws
        dxs=w2*xs1-w1*xsh2
        dcw1=dc1*phi2
        dcw2=dch2*phi1
        return
        end associate
        end

        subroutine tzsetparam0(tz,dp,akk)
        implicit none
        type (tzparam) tz
        real*8 dp,akk,xsin,xsinh
        external xsin,xsinh
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

        pr=(1.d0+dp)
        akkp=akk/pr
        w1=sqrt(akkp)
        phi1=aln*w1
        c1=cos(phi1)
        s1=sin(phi1)
        xs1=xsin(phi1)
        if(c1 .ge. 0.d0)then
          dc1=-s1**2/(1.d0+c1)
        else
          dc1=c1-1.d0
        endif
        ch2=cosh(phi1)
        sh2=sinh(phi1)
        xsh2=xsinh(phi1)
        dch2=sh2**2/(1.d0+ch2)
        return
        end associate
        end
        
        subroutine tzsetparamp(tz)
        implicit none
        type (tzparam) tz
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

        w1p=-w1**3/pr*wss
        w2p=-w2**3/pr*wss
        phi1p=aln*w1p
        phi2p=aln*w2p
        wsp=w1p+w2p
        wdp=-wd*(1.d0/pr+wsp/ws)
        wssip=2.d0*(w1**4+w2**4)*wss**2/pr
        w12p=w1p-w2p
        c1p=-s1*phi1p
        s1p=c1*phi1p
        xs1p=-dc1*phi1p
        ch2p=sh2*phi2p
        sh2p=ch2*phi2p
        xsh2p=-dch2*phi2p        
        g1p = g1/pr + 2*s1*(s1p*w2 + s1*w2p)/w1
        g2p = g2/pr - 2*sh2*(sh2p*w1 + sh2*w1p)/w2
        ca1p = (-(ws*(s1p*sh2*w12 + s1*sh2p*w12 + s1*sh2*w12p +
     $       (c1p*ch2 + ch2p + ch2p*dc1)*ws)) + s1*sh2*w12*wsp)/ws**2
        wr1p = (w1p - wr1*w2p)/w2
        wr2p = (-(w1p*wr2) + w2p)/w1
        cxs1p = c1p*s1 + dc1*s1p - xs1p
        cxs2p = ch2p*sh2 + dch2*sh2p - xsh2p
        aw1p = (cxs2p*wr1 + cxs2*wr1p - cxs2*wr1*wsp/ws)/ws
        aw2p = (cxs1p*wr2 + cxs1*wr2p - cxs1*wr2*wsp/ws)/ws
        cr2p = c1p*wr2 + c1*wr2p
        cr3p = ch2p*wr1 + ch2*wr1p
        csw1p = akkp*aln*(-c1p + ch2p) - csw1/pr
        cswsp = (ch2p*w1*phi1 + 2*ch2*phi1*w1p +
     $       c1p*w2*phi2 + 2*c1*phi2*w2p)*wss + csws*wssip
        dcw1p = c1p*phi2 + dc1*phi2p
        dcw2p = ch2p*phi1 + dch2*phi1p
        dxsp = w2p*xs1 + w2*xs1p - w1p*xsh2 - w1*xsh2p
        return
        end associate
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
     $     tbrhoz,bx,by,bxy,b1,br,bz0,cw,phieps,
     $     awu,dwu,awup,dwup,
     $     dz1,dz2,dz1p,dz2p
      logical*4 enarad,calpol,radcod
      external itgetqraddiv,tbrhoz
      parameter (phieps=1.d-2)
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
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      br=tbrhoz()
      b1=br*akk
      call tinitr(trans1)
c     end   initialize for preventing compiler warning
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
     $           0.d0,0.d0,0.d0,0.d0,.5d0*aln)
          else
            call trade(trans,beam,cod,bx,by,bz*br,bz,
     $           0.d0,bxy,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,aln)
          endif
        endif
        if(n .eq. 1)then
          call tqente(trans,cod,beam,aln*.5d0,bz,calpol,irad,ld)
        else
          call tqente(trans,cod,beam,aln,bz,calpol,irad,ld)
        endif
        xi0=cod(1)
        yi0=cod(3)
        xi=xi0+dx0
        yi=yi0+dy0
        bzh=bz*.5d0
        pxi=(cod(2)+yi0*bzh)/pr
        pyi=(cod(4)-xi0*bzh)/pr
        a = (w2*ws*xi-bzp*pyi)*wss
        bw= (ws*pxi-bzp*w2*yi)*wss
        b = bw*w1
        c = (w1*wd*xi+pyi)*wss
        cw=c*w2
        dw= (-wd*pxi+w1*yi)*wss
        d = dw*w2
        u1w= a*dc1+bw*s1
        u1=u1w*w1
        u2w=-a*s1+bw*dc1
        u2=u2w*w1
        v1w= c*dch2+dw*sh2
        v1=v1w*w2
        v2w= c*sh2+dw*dch2
        v2=v2w*w2
        xf =xi0+u1w+v1w*bzp
        pxf=pxi+u2 +v2*bzp
        dy =wd*u2w+ws*v2w
        yf =yi0+dy
        dpy=-wd*u1 +ws*v1
        pyf=pyi+dpy
        adx=w2*ws*wss
        adpy=-bzp*wss
        adp=a*wssip+
     $       ((w2p*ws+w2*wsp)*xi+bzp*pyi/pr)*wss
        bwdpx=ws*wss
        bwdy=-bzp*w2*wss
        bwdp=bw*wssip+
     $       (wsp*pxi-bzp*(w2p-w2/pr)*yi)*wss
        bdpx=bwdpx*w1
        bdy=bwdy*w1
        bdp=bwdp*w1+bw*w1p
        cdx=w1*wd*wss
        cdpy=wss
        cdp=c*wssip+((w1p*wd+w1*wdp)*xi)*wss
        cwdp=cdp*w2+cw*w2p
        dwdpx=-wd*wss
        dwdy=w1*wss
        dwdp=dw*wssip+(-wdp*pxi+w1p*yi)*wss
c        write(*,'(a,1p7g15.7)')': ',dw*wssip,-wdp*pxi*wss,w1p*wss*yi,
c     $       cdp*dch2*bzp,c*ch2p*bzp,dwdp*sh2*bzp,dw*sh2p*bzp
        ddpx=dwdpx*w2
        ddy=dwdy*w2
        ddp=dwdp*w2+dw*w2p
        u1wx=adx*dc1
        u1wpx=bwdpx*s1
        u1wy=bwdy*s1
        u1wpy=adpy*dc1
        u2wx=-adx*s1
        u2wpx=bwdpx*dc1
        u2wy=bwdy*dc1
        u2wpy=-adpy*s1
        v1wx =cdx*dch2
        v1wpx=dwdpx*sh2
        v1wy =dwdy*sh2
        v1wpy=cdpy*dch2
        v2wx =cdx*sh2
        v2wpx=dwdpx*dch2
        v2wy =dwdy*dch2
        v2wpy=cdpy*sh2
        u1wp=adp*dc1+a*c1p+bwdp*s1+bw*s1p
        u1p=u1wp*w1+u1w*w1p
        u2wp=-adp*s1-a*s1p+bwdp*dc1+bw*c1p
        u2p=u2wp*w1+u2w*w1p
        v1wp=cdp*dch2+c*ch2p+dwdp*sh2+dw*sh2p
        v1p=v1wp*w2+v1w*w2p
        v2wp=cdp*sh2+c*sh2p+dwdp*dch2+dw*ch2p
        v2p=v2wp*w2+v2w*w2p
        trans1(1,1)=1.d0+u1wx+v1wx*bzp
        trans1(1,2)=u1wpx+v1wpx*bzp
        trans1(1,3)=u1wy+v1wy*bzp
        trans1(1,4)=u1wpy+v1wpy*bzp
        trans1(1,6)=u1wp+(v1wp-v1w/pr)*bzp
        trans1(2,1)=(w1*u2wx+w2*v2wx*bzp)
        trans1(2,2)=1.d0+(w1*u2wpx+w2*v2wpx*bzp)
        trans1(2,3)=(w1*u2wy+w2*v2wy*bzp)
        trans1(2,4)=(w1*u2wpy+w2*v2wpy*bzp)
        trans1(2,6)=u2p+(v2p-v2/pr)*bzp
        trans1(3,1)=wd*u2wx+ws*v2wx
        trans1(3,2)=wd*u2wpx+ws*v2wpx
        trans1(3,3)=1.d0+wd*u2wy+ws*v2wy
        trans1(3,4)=wd*u2wpy+ws*v2wpy
        trans1(3,6)=wdp*u2w+wd*u2wp+wsp*v2w+ws*v2wp
        trans1(4,1)=-wd*w1*u1wx+ws*w2*v1wx
        trans1(4,2)=-wd*w1*u1wpx+ws*w2*v1wpx
        trans1(4,3)=-wd*w1*u1wy+ws*w2*v1wy
        trans1(4,4)=1.d0-wd*w1*u1wpy+ws*w2*v1wpy
        trans1(4,6)=-wdp*u1-wd*u1p+wsp*v1+ws*v1p
        awu=a/ws*w1
        dwu=d
        awup=adp/ws*w1-a/ws*(w2p-wsp/ws*w2)
        dwup=ddp
        call tztaf(1,tz,awu,pxi,pyi,aw1,ws,w12,wss,
     $       g1,awup,aw1p,wsp,w12p,g1p,
     $       dz1,dz1p)
        call tztaf(1,tz,-dwu,-pyi,pxi,aw2,-w12,ws,-wss,g2,
     $       -dwup,aw2p,-w12p,wsp,g2p,
     $       dz2,dz2p)
        cod(5)=cod(5)+
     $       bzp*(-((awu*dwu*dxs**2)/akkp) +
     $       ca1*pxi*pyi*wss)
     $       +dz1+dz2-aln*dv
        trans1(5,6)=
     $       -((bzp*dxs*(awup*dwu*dxs +
     $       awu*dwup*dxs + 2*awu*dwu*dxsp))/
     -     akkp) + bzp*pxi*pyi*wss*(ca1p +
     $       ca1*(-1./pr + wssip))
     $       +dz1p+dz2p+dvdp*aln
        cod(1)=xf
        cod(2)=pxf*pr-bzh*yf
        cod(3)=yf
        cod(4)=pyf*pr+bzh*xf
        do i=1,4
          trans1(i,2)=trans1(i,2)/pr
          trans1(i,4)=trans1(i,4)/pr
          trans1(i,1)=trans1(i,1)-bzh*trans1(i,4)
          trans1(i,3)=trans1(i,3)+bzh*trans1(i,2)
          trans1(i,6)=trans1(i,6)
     $         -(pxi*trans1(i,2)+pyi*trans1(i,4))
        enddo
        do i=1,6
          trans1(2,i)=trans1(2,i)*pr-bzh*trans1(3,i)
          trans1(4,i)=trans1(4,i)*pr+bzh*trans1(1,i)
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
      call tqente(trans,cod,beam,aln*.5d0,bz,calpol,irad,ld)
      if(enarad)then
        bx= b1*cod(3)
        by= b1*cod(1)
        bxy= b1
        call trade(trans,beam,cod,bx,by,bz*br,bz,
     $       0.d0,bxy,0.d0,0.d0,
     $       0.d0,0.d0,0.d0,0.d0,.5d0*aln)
      endif
      return
      end associate
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
      b=brhoz*a/al
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

      real*8 function tbrhoz()
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      tbrhoz=brhoz
      return
      end
