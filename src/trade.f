      subroutine trade(trans,beam,cod,bx,by,bz,br,
     $     bxx,bxy,byy,dldx,al,s,ala,f1,f2,prev,next)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 i
      real*8 f1,al,dldx,bx,by,bz,
     $     bxx,bxy,byy,pr,pxi,pyi,dtdsc,ar,br,
     $     p1,brad,dtds,h1sq,alr,tlc,dp,tlc1,
     $     dp1,de,sq,fact,fact1,brad2,brad3,dbrad2,
     $     vx,vy,vz,qx,qy,qzi,pzi,htlc,htlc2,
     $     s,ala,f2,alr2,alr3,f1r,f2r,alr1
      real*8 , parameter :: pmax=0.9999d0, pmin=0.99999d0
      real*8 trans(6,12),beam(42),cod(6)
      real*8 radi(6,6)
      logical*4 prev,next
c      write(*,'(a,1p6g15.7)')'trade ',bx,by,bxx,bxy,byy,dldx
      pr=1.d0+cod(6)
      p1=p0*pr
      ar=.5d0*br
      pxi=(cod(2)+ar*cod(3))/pr
      pyi=(cod(4)-ar*cod(1))/pr
      pzi=1.d0+pxy2dpz(pxi,pyi)
      sq=pxi**2+pyi**2
c      pzi=sqrtl(1.d0-sq)
      vx=pyi*bz-pzi*by
      vy=pzi*bx-pxi*bz
      vz=pxi*by-pyi*bx
      brad2=vx**2+vy**2+vz**2
      alr=max(0.d0,al+cod(1)*dldx)
      brad=sqrt(brad2)
      brad3=brad2*brad
      if(prev .and. s .eq. 0.d0)then
        dbrad2=-(brad-bradprev)**2*f1/alr/6.d0
        brad2=brad2+dbrad2
        brad3=brad3+dbrad2*1.5d0*(brad+bradprev)
      else
        if(prev)then
          f1r=0.d0
        else
          f1r=f1
        endif
        if(next)then
          f2r=0.d0
        else
          f2r=f2
        endif
        if(f1r .ne. 0.d0 .or. f2r .ne. 0.d0)then
          call tradel(al,f1r,f2r,s,ala,alr1,alr2,alr3)
          brad2=brad2*alr2/al
          brad3=brad3*alr3/al
        endif
      endif
      bradprev=brad
      dtds=1.d0/pzi
      dtdsc=dtds*crad
      h1sq=1.d0+p1**2
      tlc=alr*dtdsc
      htlc=-h1sq*tlc
      htlc2=htlc*2.d0
      dp=max(htlc*brad2,-pmin*p0)
      u0=u0-dp
      if(irad .gt. 6)then
        qzi=-sq/pzi
        qx=pyi*bz-qzi*by
        qy=qzi*bx-pxi*bz
        radi=0.d0
        radi(6,1)=htlc2*
     $       (pzi*(-vx*bxy+vy*bxx)+vz*(pxi*bxy-pyi*bxx))
     1       -h1sq*brad2*dtdsc*dldx
        radi(6,2)=(dp*dtds**2*pxi+htlc2*
     $       ((vx*by-vy*bx)*pxi/pzi-vy*bz+vz*by))/pr
        radi(6,3)=htlc2*
     $       (pzi*(-vx*byy+vy*bxy)+vz*(pxi*byy-pyi*bxy))
        radi(6,4)=(dp*dtds**2*pyi+htlc2*
     $       ((-vx*by+vy*bx)*pyi/pzi+vx*bz-vz*bx))/pr
        radi(6,5)=0.d0
        radi(6,6)=-2.d0*p1*p0*tlc*brad2
     1       -dp*dtds**2*sq/pr
     1       -htlc2*(qx*vx+qy*vy+vz**2)/pr
c        write(*,'(1p6g15.7)')(radi(6,i),i=1,6)
        call tradp(radi,pxi,pyi,dp,pr,ar)
        call tmuld(trans,radi)
        do i=1,6
          radi(i,i)=radi(i,i)+1.d0
        enddo
        tlc1=alr*dtdsc*.5d0
        dp1=-h1sq*brad3*tlc1
        de=-erad*dp1*h1sq*sqrt(h1sq)/brhoz/pr
        beam(21)=beam(21)+de
        fact=1.d0+de/pr**2
        fact1=1.d0-de/pr**2
        beam(3)=beam(3)*fact
     $       +(ar*(2.d0*beam(5)+ar*beam(6))+pxi**2)*de
        beam(8)=beam(8)*fact
     $       +(ar*(beam(9)-beam(2)-ar*beam(4))+pxi*pyi)*de
        beam(10)=beam(10)*fact
     $       +(ar*(-2.d0*beam(7)+ar*beam(1))+pyi**2)*de
        beam(17)=beam(17)*fact1+(-ar*beam(18)/pr**2+pxi)*de
        beam(19)=beam(19)*fact1+( ar*beam(16)/pr**2+pyi)*de
        call tmulbs(beam,radi,.false.,.true.)
        beam(21)=beam(21)+de
        beam(3)=beam(3)*fact
     $       +(ar*(2.d0*beam(5)+ar*beam(6))+pxi**2)*de
        beam(8)=beam(8)*fact
     $       +(ar*(beam(9)-beam(2)-ar*beam(4))+pxi*pyi)*de
        beam(10)=beam(10)*fact
     $       +(ar*(-2.d0*beam(7)+ar*beam(1))+pyi**2)*de
        beam(17)=beam(17)*fact1+(-ar*beam(18)/pr**2+pxi)*de
        beam(19)=beam(19)*fact1+( ar*beam(16)/pr**2+pyi)*de
      endif
      if(radcod)then
        cod(6)=cod(6)+dp
        cod(2)=cod(2)+pxi*dp
        cod(4)=cod(4)+pyi*dp
        call tesetdv(cod(6))
      endif
      return
      end

      recursive subroutine tradel(al,f1,f2,s,ala,alr1,alr2,alr3)
      implicit none
      real*8 al,f1,f2,s,ala,alr1,alr2,alr3,sc
      real*8 sa,sb,s1,s2,ra,rb,r1,r2,da1,d12,d2b,r0,dl
      if(al .lt. 0.d0)then
        call tradel(-al,f1,f2,s,-ala,alr1,alr2,alr3)
        alr1=-alr1
        alr2=-alr2
        alr3=-alr3
        return
      endif
      alr1=al
      alr2=al
      alr3=al
      if(s .eq. 0.d0)then
        if(f1 .ne. 0.d0)then
          r0=min(.5d0+al/f1,1.d0)
          alr1=r0**2/2.d0*f1
          alr2=alr1*r0/1.5d0
          alr3=alr2*r0*.75d0
          dl=max(al-f1*.5d0,0.d0)
          alr1=alr1+dl
          alr2=alr2+dl
          alr3=alr3+dl
        endif
      elseif(s .eq. ala)then
        if(f2 .ne. 0.d0)then
          r0=min(.5d0+al/f2,1.d0)
          alr1=r0**2/2.d0*f1
          alr2=alr1*r0/1.5d0
          alr3=alr2*r0*.75d0
          dl=max(al-f2*.5d0,0.d0)
          alr1=alr1+dl
          alr2=alr2+dl
          alr3=alr3+dl
        endif
      else
        if(f1 .ne. 0.d0 .or. f2 .ne. 0.d0)then
          sa=s-al*.5d0
          sb=s+al*.5d0
          sc=f1*ala/(f1+f2)
          s1=max(sa,min(sc,sb,    f1*.5d0))
          s2=min(sb,max(sc,sa,ala-f2*.5d0))
          if(f1 .ne. 0.d0)then
            if(f2 .ne. 0.d0)then
              ra=min(.5d0+sa/f1,.5d0+(ala-sa)/f2,1.d0)
              rb=min(.5d0+sb/f1,.5d0+(ala-sb)/f2,1.d0)
              r1=min(.5d0+s1/f1,.5d0+(ala-s1)/f2,1.d0)
              r2=min(.5d0+s2/f1,.5d0+(ala-s2)/f2,1.d0)
            else
              ra=min(.5d0+sa/f1,1.d0)
              rb=min(.5d0+sb/f1,1.d0)
              r1=min(.5d0+s1/f1,1.d0)
              r2=min(.5d0+s2/f1,1.d0)
            endif
          else
            ra=min(.5d0+(ala-sa)/f2,1.d0)
            rb=min(.5d0+(ala-sb)/f2,1.d0)
            r1=min(.5d0+(ala-s1)/f2,1.d0)
            r2=min(.5d0+(ala-s2)/f2,1.d0)
          endif
          da1=max(s1-sa,0.d0)
          d12=max(s2-s1,0.d0)
          d2b=max(sb-s2,0.d0)
          alr1=(da1*(ra+r1)+
     $          d12*(r1+r2)+
     $          d2b*(r2+rb))/2.d0
          alr2=(da1*(ra*(ra+r1)+r1**2)+
     $         d12*(r1*(r1+r2)+r2**2)+
     $         d2b*(r2*(r2+rb)+rb**2))/3.d0
          alr3=(da1*(ra**2+r1**2)*(ra+r1)+
     $         d12*(r1**2+r2**2)*(r1+r2)+
     $         d2b*(r2**2+rb**2)*(r2+rb))/4.d0
        endif
      endif
      return
      end

      subroutine tradp(trans,pxi,pyi,dp,pr,ar)
      implicit none
      real*8 trans(6,6),pxi,pyi,dp,pr,ar
      trans(6,1)=trans(6,1)-ar*trans(6,4)
      trans(6,3)=trans(6,3)+ar*trans(6,2)
      trans(2,1)=pxi*trans(6,1)
      trans(2,2)=pxi*trans(6,2)
      trans(2,3)=pxi*trans(6,3)
      trans(2,4)=pxi*trans(6,4)
      trans(2,5)=pxi*trans(6,5)
      trans(2,6)=pxi*trans(6,6)
      trans(2,2)=trans(2,2)+dp/pr
      trans(2,6)=trans(2,6)-pxi*dp/pr
      trans(4,1)=pyi*trans(6,1)
      trans(4,2)=pyi*trans(6,2)
      trans(4,3)=pyi*trans(6,3)
      trans(4,4)=pyi*trans(6,4)
      trans(4,5)=pyi*trans(6,5)
      trans(4,6)=pyi*trans(6,6)
      trans(4,4)=trans(4,4)+dp/pr
      trans(4,6)=trans(4,6)-pyi*dp/pr
      return
      end

      subroutine trades(trans,beam,cod,bzs0,bzs1,f1,brhoz)
      use tmacro, only:bradprev
      implicit none
      real*8 trans(6,12),beam(42),cod(6),bxr,byr,bzs0,bzs1,
     $             bxx,f1,bzr,brhoz
      bzr=(bzs1-bzs0)*brhoz
      return
      if(bzr .eq. 0.d0)then
        return
      endif
      bxx=-.5d0*bzr/f1
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,bzs0,
     $     bxx,0.d0,bxx,0.d0,
     $     f1*.25d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,.5d0*(bzs0+bzs1),
     $     bxx,0.d0,bxx,0.d0,
     $     f1*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,bzs1,
     $     bxx,0.d0,bxx,0.d0,
     $     f1*.25d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
      bradprev=0.d0
      return
      end

      subroutine tradke(trans,cod,beam,srot,al,phir0,bzh)
      use tmacro
      use temw
      use tfstk, only:pxy2dpz,p2h
      use ffs_flag,only:radcod,calpol
      use tspin,only:sprotm
      implicit none
      real*8 , intent(inout)::trans(6,12),cod(6),beam(42),
     $     srot(3,9)
      real*8 , intent(in)::al,bzh,phir0
      real*8 transi(6,6),tr1(6,6),dppa(6),tr2(6,6),
     $     ddpz(6),dal(6),duc(6),dddpx(6),dddpy(6),ddg(6),
     $     dtheta(6),danp(6),dbeam(21),dpxi(6),dpyi(6),
     $     c1,dpx,dpy,ddpx,ddpy,pxr0,ct,pz00,das,
     $     pr,px,py,pz,pz0,ppx,ppy,ppz,ppa,theta,
     $     p,h1,al1,anp,uc,dg,g,pr1,pxi,pyi,
     $     p2,h2,de,cp,sp,b,pxm,pym,gi,a,dh1r,
     $     pxh,pyh,pzh,ppzb,btx,bty,btz,dct,
     $     dpxh(6),dpyh(6),dpzh(6),bp,dbp(6),dpxr0(6),dpz0(6),
     $     dppx(6),dppy(6),dppzb(6),dblx(6),dbly(6),dblz(6),
     $     dbtx(6),dbty(6),dbtz(6),dgx(6),dgy(6),dgz(6),dpz00(6)
      real*8, parameter:: gmin=-0.9999d0,
     $     cave=8.d0/15.d0/sqrt(3.d0),cuu=11.d0/27.d0,
     $     cl=1.d0+gspin
      integer*4 i
      gi=codr0(6)
      pr=1.d0+gi
      cp=cos(phir0)
      sp=sin(phir0)
      pxi=codr0(2)+bzhr0*codr0(3)
      pyi=codr0(4)-bzhr0*codr0(1)
      pz00=pr*(1.d0+pxy2dpz(pxi/pr,pyi/pr))
      pxr0= cp*pxi+sp*pz00
      pz0 =-sp*pxi+cp*pz00
      px=cod(2)+bzh*cod(3)
      py=cod(4)-bzh*cod(1)
      pz=pr*(1.d0+pxy2dpz(px/pr,py/pr))
      dpx=px-pxr0
      dpy=py-pyi
      ppx=(py*pz0 -pz*pyi)
      ppy=(pz*pxr0-px*pz0)
      ppz=(px*pyi-py*pxr0)
      ppa=abs(dcmplx(ppx,abs(dcmplx(ppy,ppz))))
      theta=asin(min(1.d0,ppa/pr**2))
      p=p0*pr
      h1=p2h(p)
      al1=al-cod(5)+codr0(5)
      anp=anrad*h1*theta
      uc=cuc*h1**3/p0*theta/al1
      dg=-cave*anp*uc
      u0=u0-dg
      g=max(gmin,gi+dg)
      pr1=1.d0+g
      ddpx=-.5d0*dpx*dg
      ddpy=-.5d0*dpy*dg
      c1=al/pr/3.d0
      if(radcod)then
        cod(1)=cod(1)+ddpx*c1
        cod(3)=cod(3)+ddpy*c1
        cod(2)=px*pr1/pr+ddpx-bzh*cod(3)
        cod(4)=py*pr1/pr+ddpy+bzh*cod(1)
        cod(6)=g
        p2=p0*pr1
        h2=p2h(p2)
        cod(5)=cod(5)*p2/h2*h1/p
        call tesetdv(g)
      else
        p2=p
        h2=h1
      endif
      if(irad .gt. 6)then
        call tinv6(transr,transi)
        call tmultr(transi,trans(:,1:6),6)
        tr2=transi
        if(bzh .ne. 0.d0)then
          tr2(2,:)=tr2(2,:)+bzh*tr2(3,:)
          tr2(4,:)=tr2(4,:)-bzh*tr2(1,:)
        endif
        ddpz=(tr2(6,:)*pr-tr2(2,:)*px-tr2(4,:)*py)/pz
        ppzb=ppz+(bsi+bzh*2.d0*al)/pr
        dpxi=(/0.d0,1.d0,bzhr0,0.d0,0.d0,0.d0/)
        dpyi=(/-bzhr0,0.d0,0.d0,1.d0,0.d0,0.d0/)
        dpz00=(-pxi*dpxi-pyi*dpyi)/pz00
        dpz00(6)=dpz00(6)+pr/pz00
        dpxr0=cp*dpxi+sp*dpz00
        dpz0=(-pxr0*dpxr0-pyi*dpyi)/pz0
        dpz0(6)=dpz0(6)+pr/pz0
        dppx=tr2(4,:)*pz0+py*dpz0-ddpz*pyi-pz*dpyi
        dppy=ddpz*pxr0+pz*dpxr0-tr2(2,:)*pz0-px*dpz0
        dppzb=tr2(2,:)*pyi+px*dpyi-tr2(4,:)*pxr0-py*dpxr0
        dh1r=p*p0/h1**2
        if(ppa .ne. 0.d0)then
          dppa=(ppx*dppx+ppy*dppy+ppz*dppzb)/ppa
          dal=-tr2(5,:)
          dal(5)=dal(5)+1.d0
          das=1.d0/sqrt(1.d0-ppa**2)
          dtheta=dppa*das
          dtheta(6)=dtheta(6)-2.d0*ppa/pr**3*das
          danp=anrad*h1*dtheta
          danp(6)=danp(6)+anp*dh1r
          duc=uc*(dtheta/theta-dal/al1)
          duc(6)=duc(6)+3.d0*uc*dh1r
          ddg=-cave*(danp*uc+anp*duc)
          dddpx=-.5d0*((tr2(2,:)-dpxr0)*dg+ddpx*ddg)
          dddpy=-.5d0*((tr2(4,:)-dpyi )*dg+ddpy*ddg)
          tr1(1,:)=c1*dddpx
          tr1(1,6)=tr1(1,6)-ddpx/pr
          tr1(3,:)=c1*dddpy
          tr1(3,6)=tr1(3,6)-ddpy/pr
          tr1(2,:)=(tr2(2,:)*dg+px*ddg)/pr+dddpx
          tr1(2,6)=tr1(2,6)-px*dg/pr
          tr1(4,:)=(tr2(4,:)*dg+py*ddg)/pr+dddpy
          tr1(4,6)=tr1(4,6)-py*dg/pr
c     derivative of dz has been ignored.
          tr1(5,:)=0.d0
          tr1(6,:)=ddg
c          write(*,'(a,1p8g15.7)')'tradke  ',tr2(2,:)
c          write(*,'(a,1p8g15.7)')' ddg    ',ddg,dg
c          write(*,'(a,1p8g15.7)')' danp   ',danp,anp
c          write(*,'(a,1p8g15.7)')' duc    ',duc,uc
c          write(*,'(a,1p8g15.7)')' dtheta ',dtheta,theta
c          write(*,'(a,1p8g15.7)')' dppy   ',dppy,ppy
c          write(*,'(a,1p8g15.7)')' dpxr0  ',dpxr0,pxr0
c          do i=1,6
c            write(*,'(1p6g15.7)')tr1(i,:)
c          enddo
          if(bzh .ne. 0.d0)then
            tr1(2,:)=tr1(2,:)-bzh  *tr1(3,:)
            tr1(4,:)=tr1(4,:)+bzh  *tr1(1,:)
          endif
          call tmuld6(trans,tr1)
          tr1(1,1)=tr1(1,1)+1.d0
          tr1(2,2)=tr1(2,2)+1.d0
          tr1(3,3)=tr1(3,3)+1.d0
          tr1(4,4)=tr1(4,4)+1.d0
          tr1(5,5)=tr1(5,5)+1.d0
          tr1(6,6)=tr1(6,6)+1.d0
          call tmulbs(beam,tr1,.false.,calint)
          de=anp*uc**2*cuu
          pxm=pxi+px
          pym=pyi+py
          pxh=pxm*.5d0
          pyh=pym*.5d0
          b=bzh*.5d0
          dbeam=0.d0
          dbeam(3)=(beam(3)+b*(2.d0*beam(5)+b*beam(6))
     $         +(pxm**2+pxi**2+px**2)/6.d0)*de
          dbeam(8) =(beam(8)-b*(beam(2)-beam(10)+b*beam(4))
     $         +(pxm*pym+pxi*pyi+px*py)/6.d0)*de
          dbeam(10)=(beam(10)+b*(-2.d0*beam(7)+b*beam(1))
     $         +(pym**2+pyi**2+py**2)/6.d0)*de
          dbeam(17)=pxh*de
          dbeam(19)=pyh*de
          dbeam(21)=de
          beam(1:21)=beam(1:21)+dbeam
          if(calint)then
            beam(22:42)=beam(22:42)+dbeam
          endif
        endif
        if(calpol)then
          pxh=(pxr0+px)/pr*.5d0
          pyh=(pyi+py)/pr*.5d0
          dpxh=(tr2(2,:)+dpxr0)/pr*.5d0
          dpxh(6)=dpxh(6)-pxh/pr
          dpyh=(tr2(4,:)+dpyi)/pr*.5d0
          dpyh(6)=dpyh(6)-pyh/pr
          pzh=1.d0+pxy2dpz(pxh,pyh)
          dpzh=-(pxh*dpxh+pyh*dpyh)/pzh
          bp=(ppx*pxh+ppy*pyh+ppzb*pzh)/pr
          dbp=(dppx*pxh+ppx*dpxh+dppy*pyh
     $         +ppy*dpyh+dppzb*pzh+ppzb*dpzh)/pr
          dbp(6)=dbp(6)-bp/pr
          btx=ppx/pr-bp*pxh
          bty=ppy/pr-bp*pyh
          btz=ppzb/pr-bp*pzh
          dblx=dbp*pxh+bp*dpxh
          dbly=dbp*pyh+bp*dpyh
          dblz=dbp*pzh+bp*dpzh
          dbtx=dppx/pr-dblx
          dbtx(6)=dbtx(6)-ppx/pr**2
          dbty=dppy/pr-dbly
          dbty(6)=dbty(6)-ppy/pr**2
          dbtz=dppzb/pr-dblz
          dbtz(6)=dbtz(6)-ppzb/pr**2
          ct=h1*gspin
          dct=ct*dh1r
          ct=ct+1.d0
          dgx=ct*dbtx+cl*dblx
          dgx(6)=dgx(6)+dct*btx
          dgy=ct*dbty+cl*dbly
          dgy(6)=dgy(6)+dct*bty
          dgz=ct*dbtz+cl*dblz
          dgz(6)=dgz(6)+dct*btz
          srot(1,4:9)=srot(1,4:9)
     $         +dgx(1)*transr(1,:)+dgx(2)*transr(2,:)
     $         +dgx(3)*transr(3,:)+dgx(4)*transr(4,:)
     $         +dgx(5)*transr(5,:)+dgx(6)*transr(6,:)
          srot(2,4:9)=srot(2,4:9)
     $         +dgy(1)*transr(1,:)+dgy(2)*transr(2,:)
     $         +dgy(3)*transr(3,:)+dgy(4)*transr(4,:)
     $         +dgy(5)*transr(5,:)+dgy(6)*transr(6,:)
          srot(3,4:9)=srot(3,4:9)
     $         +dgz(1)*transr(1,:)+dgz(2)*transr(2,:)
     $         +dgz(3)*transr(3,:)+dgz(4)*transr(4,:)
     $         +dgz(5)*transr(5,:)+dgz(6)*transr(6,:)
          if(ppa .eq. 0.d0)then
            a=1.d0
          else
            a=theta/ppa
          endif
          call sprotm(9,srot,
     $         pxh/pr1,pyh/pr1,ppx,ppy,ppz,bsi,a,p2,h2,cp,sp)
        endif
      endif
      codr0(1:6)=cod(1:6)
      transr=trans(:,1:6)
      bzhr0=bzh
      bsi=0.d0
      return
      end
