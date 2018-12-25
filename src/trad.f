      module tspin
      use macphys

      real*8, parameter :: pst=8.d0*sqrt(3.d0)/15.d0,
     $     sflc=.75d0*(elradi/finest)**2

c      type spin
c      sequence
c      real*8 sx,sy,sz
c      end type

      contains
        subroutine tradkf1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     px0,py0,zr0,bsi,al)
        use tfstk, only:pxy2dpz,p2h
        use ffs_flag
        use tmacro
        implicit none
        real*8, parameter:: gmin=-0.9999d0,
     $       cave=8.d0/15.d0/sqrt(3.d0)
        real*8 x,px,y,py,z,g,dv,px0,py0,zr0,bsi,al,
     $       dpx,dpy,dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,
     $       pxm,pym,al1,uc,ddpx,ddpy,h1,p2,h2,sx,sy,sz,
     $       ppa,an,a,dph,r1,r2
        dpx=px-px0
        dpy=py-py0
        dpz=pxy2dpz(px,py)
        dpz0=pxy2dpz(px0,py0)
        ppx=py*dpz0-dpz*py0+dpy
        ppy=dpz*px0-px*dpz0-dpx
        ppz=px*py0-py*px0
        ppa=abs(dcmplx(ppx,abs(dcmplx(ppy,ppz))))
        theta=asin(min(1.d0,max(-1.d0,ppa)))
        pr=1.d0+g
        p=p0*pr
        anp=anrad*p*theta
        call tdusrn(anp,dph,r1,r2,an)
        if(an .ne. 0.d0)then
          al1=al-z+zr0
c          al1=al*(1.d0+(dldx*pxm**2+pym**2)*.5d0)
          uc=cuc*(1.d0+p**2)*theta/al1*pr
          dg=-dph*uc
          g=max(gmin,g+dg)
          ddpx=-r1*dpx*dg
          ddpy=-r1*dpy*dg
          x=x+r2*ddpx*al
          y=y+r2*ddpy*al
          px=px+ddpx
          py=py+ddpy
          pr=1.d0+g
          p2=p0*pr
          h2=p2h(p2)
          dv=-g*(1.d0+pr)/h2/(h2+p2)+dvfs
          h1=p2h(p)
          z=z*p2/h2*h1/p
          if(calpol)then
            if(ppa .ne. 0.d0)then
              a=theta/ppa
            else
              a=0.d0
            endif
            pxm=px0+dpx*.5d0
            pym=py0+dpy*.5d0
            call sprot(sx,sy,sz,pxm,pym,
     $           ppx,ppy,ppz,bsi,a,
     $           al1,p2,h2,an)
          endif
        elseif(calpol)then
          h1=p2h(p)
          if(ppa .ne. 0.d0)then
            a=theta/ppa
          else
            a=0.d0
          endif
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi,a,
     $         al,p,h1,-1.d0)
        endif
        return
        end subroutine

        subroutine tradk1(x,px,y,py,z,g,dv,sx,sy,sz,px0,py0,zr0,bsi,al)
        use tfstk, only:pxy2dpz,p2h
        use ffs_flag
        use tmacro
        implicit none
        real*8 x,px,y,py,z,g,dv,px0,py0,zr0,bsi,al,a,
     $       dpz,dpz0,ppx,ppy,ppz,theta,pr,p,anp,dg,dpx,dpy,
     $       pxm,pym,al1,uc,ddpx,ddpy,h2,h1,sx,sy,sz,ppa,p2
        real*8, parameter:: gmin=-0.9999d0,
     $       cave=8.d0/15.d0/sqrt(3.d0)
        dpx=px-px0
        dpy=py-py0
        dpz=pxy2dpz(px,py)
        dpz0=pxy2dpz(px0,py0)
        ppx=py*dpz0-dpz*py0+dpy
        ppy=dpz*px0-px*dpz0-dpx
        ppz=px*py0-py*px0
        ppa=abs(dcmplx(ppx,abs(dcmplx(ppy,ppz))))
        theta=asin(min(1.d0,max(-1.d0,ppa)))
        pr=1.d0+g
        p=p0*pr
        h1=p2h(p)
        al1=al-z+zr0
c        al1=al*(1.d0+(dldx*pxm**2+pym**2)*.5d0)
        anp=anrad*p*theta
        uc=cuc*(1.d0+p**2)*theta/al1*pr
        dg=-cave*anp*uc
        g=max(gmin,g+dg)
        ddpx=-.5d0*dpx*dg
        ddpy=-.5d0*dpy*dg
        x=x+ddpx*al/3.d0
        y=y+ddpy*al/3.d0
        px=px+ddpx
        py=py+ddpy
        pr=1.d0+g
        p2=p0*pr
        h2=p2h(p2)
        dv=-g*(1.d0+pr)/h2/(h2+p2)+dvfs
        z=z*p2/h2*h1/p
        if(calpol)then
          if(ppa .ne. 0.d0)then
            a=theta/ppa
          else
            a=0.d0
          endif
          pxm=px0+dpx*.5d0
          pym=py0+dpy*.5d0
          call sprot(sx,sy,sz,pxm,pym,ppx,ppy,ppz,bsi,a,
     $         al1,p2,h2,anp)
        endif
        return
        end subroutine

        subroutine sprot(sx,sy,sz,pxm,pym,bx0,by0,bz0,bsi,a,
     $     al,p,h,anph)
        use tfstk,only:pxy2dpz,ktfenanq,sqrt1
        use tmacro
        implicit none
        real*8 pxm,pym,bsi,h,pzm,bx0,by0,bz0,sx,sy,sz,
     $       bx,by,bz,bp,blx,bly,blz,btx,bty,btz,ct,
     $       gx,gy,gz,g,a,p,al,
     $       gnx,gny,gnz,
     $       sux,suy,suz,
     $       tnx,tny,tnz,
     $       slx,sly,slz,
     $       bt,st,st1,r,sl,sl1,
     $       sw,anph,cosu,sinu,dcosu
        real*8 , parameter :: cl=1.d0+gspin
        pzm=1.d0+pxy2dpz(pxm,pym)
        bx=bx0*a
        by=by0*a
        bz=bz0*a+bsi
        bp=bx*pxm+by*pym+bz*pzm
        blx=bp*pxm
        bly=bp*pym
        blz=bp*pzm
        btx=bx-blx
        bty=by-bly
        btz=bz-blz
        if(anph .gt. 0.d0)then
          bt=abs(dcmplx(btx,abs(dcmplx(bty,btz))))
          if(bt .ne. 0.d0)then
            tnx=btx/bt
            tny=bty/bt
            tnz=btz/bt
            st=sx*tnx+sy*tny+sz*tnz
            slx=sx-st*tnx
            sly=sy-st*tny
            slz=sz-st*tnz
            sl=abs(dcmplx(slx,abs(dcmplx(sly,slz))))
            st1=min(1.d0,max(-1.d0,
     $           st+(pst-st)*sflc*anph*(bt*h*p/al)**2))
            sl1=1.d0+sqrt1(-st1**2)
            if(sl .ne. 0.d0)then
              r=sl1/sl
              sx=slx*r+st1*tnx
              sy=sly*r+st1*tny
              sz=slz*r+st1*tnz
            else
              sx=sl1*pxm+st1*tnx
              sy=sl1*pym+st1*tny
              sz=sl1*pzm+st1*tnz
            endif
          endif
        endif
        ct=1.d0+h*gspin
        gx=ct*btx+cl*blx
        gy=ct*bty+cl*bly
        gz=ct*btz+cl*blz
        g=abs(dcmplx(gx,abs(dcmplx(gy,gz))))
        if(g .ne. 0.d0)then
          gnx=gx/g
          gny=gy/g
          gnz=gz/g
          sinu=sin(g)
          cosu=cos(g)
          dcosu=2.d0*sin(g*.5d0)**2
          sw=(sx*gnx+sy*gny+sz*gnz)*dcosu
          sux=sy*gnz-sz*gny
          suy=sz*gnx-sx*gnz
          suz=sx*gny-sy*gnx
          sx=cosu*sx+sinu*sux+sw*gnx
          sy=cosu*sy+sinu*suy+sw*gny
          sz=cosu*sz+sinu*suz+sw*gnz
        endif
        return
        end subroutine

      end module

      subroutine trad(np,x,px,y,py,g,dv,by,bx,dbydx,dldx,dldxe,al,
     $     f1,f2,als,ala,dir)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*4 np,i
      real*8 dbydx,dldx,dldxe,al,dir,brad,bx,by,tdusr,dg,
     $     pr,p,hh,dp,h1,alc,ald,dldxx,dldpx,al1,alr1,
     $     als,ala,alr2,alr3,f1,f2,uc,rhoinv0,anp,an
      real*8, parameter:: gmin=-0.9999d0
      real*8 x(np),px(np),y(np),py(np),dv(np),g(np)
      call tradel(al,f1,f2,als,ala,alr1,alr2,alr3)
c      write(*,'(a,1p7g14.6)')'trad ',al,f1,f2,als,ala,alr,alr1
      ald=dir*al
      dldxx=dldx+dldxe
c      dldpx=dldx*ald
      dldpx=dldx*ald*.5d0
      if(rfluct .and. al .gt. 0.d0)then
c        er=c/amass*erad*alr1/alr
c        call tran_array(dv(1),np)
c        do i=1,np
c          dv(i)=(dv(i)-.5d0)*3.46410161513775461d0
c        enddo
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          pr=1.d0+g(i)
          p=p0*pr
          al1=(1.d0+x(i)*dldxx+px(i)*dldpx)*alr1
     1         *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          rhoinv0=sqrt(brad)/brhoz
          anp=anrad*p0*rhoinv0*al1
          uc=cuc*(1.d0+p**2)*rhoinv0
c          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
c     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
c          de=er*sqrt(hh*brad)/p*dp*hh
c          g(i)=max(gmin,g(i)+dp+sqrt(abs(de))*dv(i))
c          if(i .eq. 1)then
c            write(*,'(a,1p8g15.7)')'trad ',anp,uc,1.d0/rhoinv,
c     $           al1*rhoinv,dg
c          endif
          dg=-uc*tdusr(anp,an)
          if(dg .ne. 0.d0)then
            g(i)=max(gmin,g(i)+dg)
            pr=1.d0+g(i)
            h1=p2h(p0*pr)
            dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
          endif
        enddo
      else
        alc=alr2*crad
        do i=1,np
          brad=((by+dbydx*x(i))**2+dbydx*ald*
     $         (by+dbydx*(x(i)+px(i)*ald/3.d0))*px(i))*(1.d0-py(i)**2)
     $        +((bx+dbydx*y(i))**2+dbydx*ald*
     $         (bx+dbydx*(y(i)+py(i)*ald/3.d0))*py(i))*(1.d0-px(i)**2)
          hh=1.d0+(p0*(1.d0+g(i)))**2
          dp=-hh*brad*(1.d0+x(i)*dldxx+px(i)*dldpx)*alc
     1        *(1.d0+(px(i)**2+py(i)**2)*.5d0)
          g(i)=max(gmin,g(i)+dp)
          pr=1.d0+g(i)
          h1=p2h(p0*pr)
          dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
        enddo
      endif
      return
      end

      subroutine tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,px0,py0,zr0,bsi,al)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),dv(np),z(np),g(np),
     $     px0(np),py0(np),zr0(np),bsi(np),al
      real*8 sx(np),sy(np),sz(np)
      if(rfluct .and. al .ne. 0.d0)then
        do i=1,np
          call tradkf1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),
     $         sx(i),sy(i),sz(i),
     $         px0(i),py0(i),zr0(i),bsi(i),al)
        enddo
      else
        do i=1,np
          call tradk1(x(i),px(i),y(i),py(i),z(i),g(i),dv(i),
     $         sx(i),sy(i),sz(i),
     $         px0(i),py0(i),zr0(i),bsi(i),al)
        enddo
      endif
      return
      end
