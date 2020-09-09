      module tsolz
      type tzparam0
      sequence
      real*8 aln,bzp,pr,akkp,w1,phi1,
     $     s1,xs1,c1,dc1,sh2,xsh2,ch2,dch2
      end type

      type tzparamp
      sequence
      real*8 wr1p,wr2p,g1p,g2p,c1p,s1p,xs1p,ch2p,sh2p,xsh2p,
     $     w1p,w2p,wsp,w12p,wdp,phi1p,phi2p,
     $     wssip, ca1p, dcw1p, dcw2p, csw1p,cswsp,
     $     cr2p, cr3p,dxsp,aw1p,aw2p,cxs1p,cxs2p
      end type

      type tzparam
      sequence
      type (tzparam0) tz0
      real*8 w2,ws,w12,wd,phi2,
     $     wss,csw1,csws,ca1,dcw1,dcw2,cr2,cr3,
     $     g1,g2,wr1,wr2,dxs,aw1,aw2,cxs1,cxs2
      type (tzparamp) tzp
      end type

      contains
        pure subroutine tztaf(tz,
     $     zwu, pxt, pyt, zw1, wst, w12t, swss, g1t, dz,
     $     zwup, zw1p, wstp, w12tp, g1tp, dzp)
        implicit none
        type (tzparam) , intent(in)::tz
        real*8 , intent(in)::zwu, pxt, pyt, zw1, wst, w12t, swss, g1t
        real*8 , intent(out)::dz
        real*8 , intent(in) , optional::zwup, zw1p, wstp, w12tp, g1tp
        real*8 , intent(out), optional::dzp
        real*8 f1,f2,f3,f1p,f2p,f3p,
     $       w12xsh,wsxs,wsm,ws1,ws2,cr1,cwsd,
     $       w12xshp,wsxsp,wsmp,ws1p,ws2p,cr1p,cwsdp,
     $       cswa,cswap
        associate (
     $       w1=>tz%tz0%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%tz0%phi1,phi2=>tz%phi2,wss=>tz%wss,
     $       bzp=>tz%tz0%bzp,akkp=>tz%tz0%akkp,aln=>tz%tz0%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,g1=>tz%g1,g2=>tz%g2,
     $       wr1=>tz%wr1,wr2=>tz%wr2,
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
     $     cr2p=>tz%tzp%cr2p, cr3p=>tz%tzp%cr3p,dxs=>tz%dxs,
     $       dxsp=>tz%tzp%dxsp,aw1=>tz%aw1,aw2=>tz%aw2,
     $       aw1p=>tz%tzp%aw1p,aw2p=>tz%tzp%aw2p,cxs1=>tz%cxs1,
     $       cxs2=>tz%cxs2,cxs1p=>tz%tzp%cxs1p,cxs2p=>tz%tzp%cxs2p)
          w12xsh = w12t*xsh2
          wsxs = wst*xs1
          cwsd = c1*w12xsh - ch2*wsxs
          cswa = (csw1 + cwsd)*swss
          wsm = 4.d0*swss*akkp
          ws1 = (-1.d0 + wsm)*wst
          ws2 = w12t*(1.d0 + wsm)
          cr1 = bzp**2*swss
          f1 = cswa + csws - zw1
          f2 = ca1*cr1 + g1t
          f3 = 2.d0*cr1*cwsd + cr3*w12xsh +
     $         dcw1*ws1 + dcw2*ws2 - cr2*wsxs
          dz=zwu*(zwu*f3/2.d0+(f2*pxt + bzp*f1*pyt)) - (zw1*pyt**2)/2.d0
          if(present(dzp))then
            w12xshp = w12tp*xsh2 + w12t*xsh2p
            wsxsp = wstp*xs1 + wst*xs1p
            cwsdp = c1p*w12xsh + c1*w12xshp - ch2p*wsxs - ch2*wsxsp
            cswap = swss*(csw1p + cwsdp) + cswa*wssip
            wsmp = wsm*(-1.d0/pr + wssip)
            ws1p = wsmp*wst + (-1.d0 + wsm)*wstp
            ws2p = w12tp*(1.d0 + wsm) + w12t*wsmp
            cr1p = cr1*(-2.d0/pr + wssip)
            f1p = cswap + cswsp - zw1p
            f2p = ca1p*cr1 + ca1*cr1p + g1tp
            f3p = 2*cr1p*cwsd + 2*cr1*cwsdp +
     $           cr3p*w12xsh + cr3*w12xshp +
     $           dcw1p*ws1 + dcw1*ws1p + dcw2p*ws2 + dcw2*ws2p -
     $           cr2p*wsxs - cr2*wsxsp
            dzp=(zwu**2*f3p)/2.d0 - (zw1p*pyt**2)/2.d0 + 
     -           zwup*(zwu*f3 + f2*pxt + bzp*f1*pyt) + 
     -           zwu*(f2p*pxt + bzp*(f1p - f1/pr)*pyt)
          endif
          return
        end associate
        end subroutine

        pure elemental function tzsetparam(dp,aln0,akk,bz) result(tz)
        use mathfun
        implicit none
        type (tzparam) tz
        real*8 ,intent(in):: dp,bz,aln0,akk
        real*8 wa
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

        if(bz .eq. 0.d0)then
          tz%tz0=tzsetparam0(dp,aln0,akk)
        else
          aln=aln0
          pr=1.d0+dp
          bzp=bz/pr
          akkp=akk/pr
          wa=hypot(bzp**2,2.d0*akkp)
          w1=sqrt((bzp**2+wa)*.5d0)
          w2=akkp/w1
          wss=1.d0/(w1**2+w2**2)
          ws=w1+w2
          w12=w1-w2
          wd=bzp/ws
          phi1=aln*w1
          phi2=aln*w2
          call xsincos(phi1,s1,xs1,c1,dc1)
          call xsincosh(phi2,sh2,xsh2,ch2,dch2)
          g1 = (s1**2*w2**2)/akkp
          g2 = -((sh2**2*w1**2)/akkp)
          wr1 = w1/w2
          wr2 = w2/w1
          cxs1 = dc1*s1   - xs1
          cxs2 = dch2*sh2 - xsh2
          aw1 = aln + (cxs2*wr1)/ws
          aw2 = aln + (cxs1*wr2)/ws
          cr2 = c1 *wr2
          cr3 = ch2*wr1
          csw1 = akkp*aln*(-dc1 + dch2)
          csws = (ch2*w1*phi1 + c1*w2*phi2)*wss
          ca1=-(ch2*dc1) - dch2 - (s1*sh2*w12)/ws
          dxs=w2*xs1-w1*xsh2
          dcw1=dc1 *phi2
          dcw2=dch2*phi1
        endif
        return
        end associate
        end function

        pure elemental function tzsetparam0(dp,aln0,akk) result(tz0)
        use mathfun, only:xsincos,xsincosh
        implicit none
        type (tzparam0) tz0
        real*8 , intent(in):: dp,aln0,akk
        associate (
     $       w1=>tz0%w1,phi1=>tz0%phi1,aln=>tz0%aln,
     $       bzp=>tz0%bzp,akkp=>tz0%akkp,pr=>tz0%pr,
     $       c1=>tz0%c1,s1=>tz0%s1,dc1=>tz0%dc1,
     $       xs1=>tz0%xs1,sh2=>tz0%sh2,dch2=>tz0%dch2,
     $       xsh2=>tz0%xsh2,ch2=>tz0%ch2)
        aln=aln0
        bzp=0.d0
        pr=1.d0+dp
        akkp=akk/pr
        w1=sqrt(akkp)
        phi1=aln*w1
        call xsincos(phi1,s1,xs1,c1,dc1)
        call xsincosh(phi1,sh2,xsh2,ch2,dch2)
        return
        end associate
        end function

        pure elemental function tzsetparamp(tz) result(tzp)
        implicit none
        type (tzparamp) tzp
        type (tzparam) ,intent(in):: tz
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
     $     c1p=>tzp%c1p,s1p=>tzp%s1p,xs1p=>tzp%xs1p,
     $       ch2p=>tzp%ch2p,sh2p=>tzp%sh2p,xsh2p=>tzp%xsh2p,
     $     w1p=>tzp%w1p,w2p=>tzp%w2p,wsp=>tzp%wsp,
     $       w12p=>tzp%w12p,wdp=>tzp%wdp,phi1p=>tzp%phi1p,
     $       phi2p=>tzp%phi2p,g1p=>tzp%g1p,g2p=>tzp%g2p,
     $       wr1p=>tzp%wr1p,wr2p=>tzp%wr2p,
     $     wssip=>tzp%wssip, ca1p=>tzp%ca1p,
     $       dcw1p=>tzp%dcw1p, dcw2p=>tzp%dcw2p,
     $       csw1p=>tzp%csw1p,cswsp=>tzp%cswsp,
     $     cr2p=>tzp%cr2p, cr3p=>tzp%cr3p,
     $       dxsp=>tzp%dxsp,aw1=>tz%aw1,aw2=>tz%aw2,
     $       aw1p=>tzp%aw1p,aw2p=>tzp%aw2p,cxs1=>tz%cxs1,
     $       cxs2=>tz%cxs2,cxs1p=>tzp%cxs1p,cxs2p=>tzp%cxs2p)

        if(bzp .eq. 0.d0)then
          w1p=-.5d0*w1/pr
        else
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
     $         (c1p*ch2 + ch2p + ch2p*dc1)*ws)) + s1*sh2*w12*wsp)/ws**2
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
     $         c1p*w2*phi2 + 2*c1*phi2*w2p)*wss + csws*wssip
          dcw1p = c1p*phi2 + dc1*phi2p
          dcw2p = ch2p*phi1 + dch2*phi1p
          dxsp = w2p*xs1 + w2*xs1p - w1p*xsh2 - w1*xsh2p
        endif
        return
        end associate
        end function
      
      end module

      subroutine tsolque(trans,cod,beam,srot,al,ak,
     $     bz0,ak0x,ak0y,eps0,enarad,irad)
      use tsolz
      use tmacro, only:bradprev
      use kradlib, only:tradke
      use temw,only:tmulbs
      implicit none
      type(tzparam) tz
      integer*4 n,ndiv
      integer*4 ,intent(in):: irad
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      real*8 trans1(6,6)
      real*8 ,intent(in):: al,ak,eps0,ak0x,ak0y,bz0
      real*8 bw,dw,dx0,dy0,
     $     xi0,yi0,dy,dpy,bz,a,b,c,d,akk,eps,bzh,
     $     adx,adpy,adp,bdy,bdpx,bdp,bwdy,bwdpx,bwdp,
     $     cdx,cdpy,cdp,ddy,ddpx,ddp,dwdy,dwdpx,dwdp,cwdp,
     $     u1wx,u1wpx,u1wy,u1wpy,
     $     u2wx,u2wpx,u2wy,u2wpy,
     $     v1wx,v1wpx,v1wy,v1wpy,
     $     v2wx,v2wpx,v2wy,v2wpy,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,
     $     u1p,u1wp,u2p,u2wp,v1p,v1wp,v2p,v2wp,
     $     dv,dvdp,xi,yi,pxi,pyi,xf,yf,pxf,pyf,
     $     tbrhoz,b1,br,cw,phieps,al1,
     $     awu,dwu,awup,dwup,dz1,dz2,dz1p,dz2p
      logical*4 ,intent(in):: enarad
      external tbrhoz
      parameter (phieps=1.d-2)
        associate (
     $       w1=>tz%tz0%w1,w2=>tz%w2,ws=>tz%ws,w12=>tz%w12,wd=>tz%wd,
     $       phi1=>tz%tz0%phi1,phi2=>tz%phi2,wss=>tz%wss,
     $       bzp=>tz%tz0%bzp,akkp=>tz%tz0%akkp,aln=>tz%tz0%aln,
     $       csw1=>tz%csw1,csws=>tz%csws,ca1=>tz%ca1,dcw1=>tz%dcw1,
     $       dcw2=>tz%dcw2,cr2=>tz%cr2,cr3=>tz%cr3,g1=>tz%g1,g2=>tz%g2,
     $       wr1=>tz%wr1,wr2=>tz%wr2,
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
     $     cr2p=>tz%tzp%cr2p, cr3p=>tz%tzp%cr3p,dxs=>tz%dxs,
     $       dxsp=>tz%tzp%dxsp,aw1=>tz%aw1,aw2=>tz%aw2,
     $       aw1p=>tz%tzp%aw1p,aw2p=>tz%tzp%aw2p,cxs1=>tz%cxs1,
     $       cxs2=>tz%cxs2,cxs1p=>tz%tzp%cxs1p,cxs2p=>tz%tzp%cxs2p)

      if(ak .eq. 0.d0)then
        call tdrife(trans,cod,beam,srot,al,
     $       bz0,ak0x,ak0y,al,.true.,enarad,irad)
        return
c      elseif(ak .lt. 0.d0)then
c        write(*,*)'tsolque-implementation error ',ak
c        stop
      endif
      bz=bz0
      eps=merge(0.1d0,0.1d0*eps0,eps0 .eq. 0.d0)
      ndiv=1+int(abs(al*hypot(ak,bz))/eps)
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      br=tbrhoz()
      b1=br*akk
      call tinitr(trans1)
c     end   initialize for preventing compiler warning
      tz=tzsetparam(cod(6),aln,akk,bz)
      tz%tzp=tzsetparamp(tz)
      call tgetdv(cod(6),dv,dvdp)
      al1=aln*.5d0
      do n=1,ndiv
        call tqente(trans,cod,beam,al1,bz,irad)
c        write(*,'(a,i5,1p8g13.4)')'tsolque-n ',n,
c     $       al,ak,trans(5,1:6)
        xi0=cod(1)
        yi0=cod(3)
        xi=xi0+dx0
        yi=yi0+dy0
        bzh=bz*.5d0
        if(bzh .eq. 0.d0)then
          pxi=cod(2)/pr
          pyi=cod(4)/pr
          u1 =   xi*dc1 +pxi*s1/w1
          u2 =-xi*w1*s1 +pxi*dc1
          v1 =   yi*dch2+pyi*sh2/w1
          v2 = yi*w1*sh2+pyi*dch2
          cod(1)=xi0 +u1
          pxf=pxi +u2 
          cod(2)=pxf*pr
          cod(3)=yi0 +v1
          pyf=pyi +v2 
          cod(4)=pyf*pr
          cod(5)=cod(5)-0.25d0*(
     $         w1*(xi**2*xs1-yi**2*xsh2)
     $         +(pxi**2+pyi**2)*aln
     $         +u1*(u2+pxi)+xi*pxi*dc1
     $         +v1*(v2+pyi)+yi*pyi*dch2)
     $         -dv*aln
          trans1(1,1)=c1
          trans1(1,2)=s1/w1
          trans1(1,6)=(-aln*xi*s1+pxi*(phi1*dc1+xs1)/w1**2)*w1p
          trans1(2,1)=-w1*s1
          trans1(2,2)=c1
          trans1(2,6)=(-xi*(phi1*c1+s1)-aln*pxi*s1)*w1p
          trans1(3,3)=ch2
          trans1(3,4)=sh2/w1
          trans1(3,6)=(aln*yi*sh2+pyi*(phi1*dch2+xsh2)/w1**2)*w1p
          trans1(4,3)=w1*sh2
          trans1(4,4)=ch2
          trans1(4,6)=(yi*(phi1*ch2+sh2)+aln*pyi*sh2)*w1p
          trans1(5,6)=-0.25d0*(
     $          w1p*(xi**2*(xs1-phi1*c1)-yi**2*(xsh2-phi1*ch2)+
     $         aln*(-xi*pxi*s1+yi*pyi*sh2))
     $         +trans1(1,6)*(u2+pxi)+u1*trans1(2,6)
     $         +trans1(3,6)*(v2+pyi)+v1*trans1(4,6))
     $         +dvdp*aln
          trans1(1:4,2)=trans1(1:4,2)/pr
          trans1(1:4,4)=trans1(1:4,4)/pr
          trans1(1:4,6)=trans1(1:4,6)
     $         -(pxi*trans1(1:4,2)+pyi*trans1(1:4,4))
          trans1(2,1:6)=trans1(2,1:6)*pr
          trans1(4,1:6)=trans1(4,1:6)*pr
        else
c cod has canonical momenta!
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
     $         ((w2p*ws+w2*wsp)*xi+bzp*pyi/pr)*wss
          bwdpx=ws*wss
          bwdy=-bzp*w2*wss
          bwdp=bw*wssip+
     $         (wsp*pxi-bzp*(w2p-w2/pr)*yi)*wss
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
c     write(*,'(a,1p7g15.7)')': ',dw*wssip,-wdp*pxi*wss,w1p*wss*yi,
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
          call tztaf(tz,awu,pxi,pyi,aw1,ws,w12,wss,g1,dz1,
     $         awup,aw1p,wsp,w12p,g1p,dz1p)
          call tztaf(tz,-dwu,-pyi,pxi,aw2,-w12,ws,-wss,g2,dz2,
     $         -dwup,aw2p,-w12p,wsp,g2p,dz2p)
          cod(5)=cod(5)
     $         +bzp*(-((awu*dwu*dxs**2)/akkp) +
     $         ca1*pxi*pyi*wss)
     $         +dz1+dz2-aln*dv
          trans1(5,6)=
     $         -((bzp*dxs*(awup*dwu*dxs +
     $         awu*dwup*dxs + 2*awu*dwu*dxsp))/
     -         akkp) + bzp*pxi*pyi*wss*(ca1p +
     $         ca1*(-1./pr + wssip))
     $         +dz1p+dz2p+dvdp*aln
          cod(1)=xf
          cod(2)=pxf*pr-bzh*yf
          cod(3)=yf
          cod(4)=pyf*pr+bzh*xf
c cod has canonical momenta!
          trans1(1:4,2)=trans1(1:4,2)/pr
          trans1(1:4,4)=trans1(1:4,4)/pr
          trans1(1:4,1)=trans1(1:4,1)-bzh*trans1(1:4,4)
          trans1(1:4,3)=trans1(1:4,3)+bzh*trans1(1:4,2)
          trans1(1:4,6)=trans1(1:4,6)
     $         -(pxi*trans1(1:4,2)+pyi*trans1(1:4,4))
          trans1(2,1:6)=trans1(2,1:6)*pr-bzh*trans1(3,1:6)
          trans1(4,1:6)=trans1(4,1:6)*pr+bzh*trans1(1,1:6)
        endif
        trans1(2,6)=trans1(2,6)+pxf
        trans1(4,6)=trans1(4,6)+pyf
        trans1(5,1:4)=
     $       -trans1(1,1:4)*trans1(2,6)
     $       +trans1(2,1:4)*trans1(1,6)
     $       -trans1(3,1:4)*trans1(4,6)
     $       +trans1(4,1:4)*trans1(3,6)
        trans1(5,6)=trans1(5,6)
     $       -(pxi*trans1(5,2)+pyi*trans1(5,4))
        call tmultr5(trans,trans1,irad)
        if(irad .gt. 6)then
          call tmulbs(beam ,trans1,.true.)
        endif
        if(enarad .and. n .ne. ndiv)then
          call tradke(trans,cod,beam,srot,aln,0.d0,bzh)
        endif
        al1=aln
      enddo
      call tqente(trans,cod,beam,aln*.5d0,bz,irad)
      if(enarad)then
        call tradke(trans,cod,beam,srot,aln,0.d0,bzh)
      endif
      bradprev=0.d0
      return
      end associate
      end

      integer*4 function itgetirad()
      use tmacro
      implicit none
      itgetirad=irad
      return
      end

      integer*4 function itgetqraddiv(cod,ak,al,bzh)
      use tfstk
      use tmacro
      implicit none
      integer*4 nrad
      real*8, intent(in):: cod(6),ak,al,bzh
      real*8 xpd,a,b,xk
      real*8 ,parameter :: xmin=1.d-6,pmin=1.d-6,amax=1.d-2
      xk=abs(ak)*max(xmin,abs(cod(1))+abs(cod(3)))
      xpd=max(pmin,abs(cod(2))+abs(cod(4)))
      a=min(amax,xk+xpd)
      b=brhoz*(a*abs(bzh)+xk/abs(al))
      nrad=int(abs(al*crad/epsrad*(h0*b)**2))
      itgetqraddiv=max(int(emidiv*emidiq*nrad),
     1       int(abs(h0*b/brhoz*anrad)/epsrad/1.d6*emidiv*emidib))
c     1     int(abs(a)/epsrad/1.d3*emidiv*emidiq))
c      write(*,'(a,i5,1p7g15.7)')'itgetqraddiv ',
c     $     itgetqraddiv,al,ak,bzh,cod(1:4)
      return
      end

      real*8 function tbrho()
      use tfstk
      use tmacro
      implicit none
      tbrho=brho
      return
      end

      real*8 function tbrhoz()
      use tfstk
      use tmacro
      implicit none
      tbrhoz=brhoz
      return
      end
