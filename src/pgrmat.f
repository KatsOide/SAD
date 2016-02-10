      subroutine pgrmat(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 kx,kax,kax1,kaxfirst,kaxlast,kfirst,klast,kam,kaxi,
     $     ka1,kaa1,ka,kaa
      logical*4 normalmode,xplane,angle,periodic
      integer*4 isp1,irtc
      integer*4 narg,nc,itfloc,itfmessage,im,m
      real*8 psix,psiy,psix0,psiy0
      character*16 keyword,tfgetstr
      
      narg=isp-isp1
      if(narg .lt. 6 .or. narg .gt. 7)then
        irtc=itfmessage(9,'General::narg','"6 or 7"')
        return
      endif
      if(ktflistq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"non-List for #1"')
        return
      endif
      irtc=0
      im=itfloc(ktastk(isp1+1),irtc)
      if(irtc .ne. 0)then
        return
      endif

      if(ktfnonlistq(ktastk(isp1+4)))then
        irtc=itfmessage(9,'General::wrongtype','"List for #4"')
        return
      endif
      kax=ktfaddr(ktastk(isp1+4))
      kax1=ktfaddr(klist(kax+1))
      kaxfirst=ktfaddr(klist(kax1+1))
      kaxlast=ktfaddr(klist(kax1+nlat))
      kfirst=ktfaddr(klist(kaxfirst+1))
      klast=ktfaddr(klist(kaxlast+1))

      if(ktfnonrealq(ktastk(isp1+5)))then
        irtc=itfmessage(9,'General::wrongtype','"Real for #5"')
        return
      endif
      periodic=rtastk(isp1+5) .ne. 0.d0

      if(ktfnonstringq(ktastk(isp1+6)))then
        irtc=itfmessage(9,'General::wrongtype','"String for #6"')
        return
      endif
      keyword=tfgetstr(ktfaddr(ktastk(isp1+6)),nc)
      if(keyword(1:1) .eq. 'P')then
        normalmode=.false.
        angle=.true.  
        if(keyword .eq. 'PX')then
          xplane=.true.
        elseif(keyword .eq. 'PY')then
          xplane=.false.
        endif
      elseif(keyword(1:1) .eq. 'M')then
        normalmode=.true.
        if(keyword .eq. 'MX')then
          xplane=.true.
          angle=.false.
        elseif(keyword .eq. 'MY')then
          xplane=.false.
          angle=.false.
        elseif(keyword .eq. 'MPX')then
          xplane=.true.
          angle=.true.
        elseif(keyword .eq. 'MPY')then
          xplane=.false.
          angle=.true.
        endif
      elseif(keyword .eq. 'X')then
        normalmode=.false.
        xplane=.true.
        angle=.false.  
      elseif(keyword .eq. 'Y')then
        normalmode=.false.
        xplane=.false.
        angle=.false.  
      endif
      psix=0.d0
      psiy=0.d0
      if(narg .eq. 7)then
        if(ktfnonlistq(ktastk(isp1+7)))then
          irtc=itfmessage(9,'General::wrongtype','"List for #7"')
          return
        endif
        kax=ktfaddr(ktastk(isp1+7))
        if(ilist(2,kax-1) .ne. 2)then
          irtc=itfmessage(9,'General::wrongleng',
     $         '"#7","2"')
          return
        endif
        kax1=ktfaddr(klist(kax+1))
        psix=rlist(kax1+1)*pi2
        psiy=rlist(kax1+2)*pi2
      endif

c     dp00=rlist(ifirst+20)+1.d0
      psix0=rlist(klast+3)-rlist(kfirst+3)
      psiy0=rlist(klast+6)-rlist(kfirst+6)
      if(psix .eq. 0.d0)then
        psix=psix0
      endif
      if(psiy .eq. 0.d0)then
        psiy=psiy0
      endif

      kaxi=ilist(2,kax1+im)
      kam=ktfaddr(klist(kaxi+1))
      if(ktfrealq(ktastk(isp1+2)))then
        if(ktfnonrealq(ktastk(isp1+3)))then
          irtc=itfmessage(9,'General::wrongtype','"Real for #3"')
          return
        endif
        
        call pgrmat1(ilist(1,ilattp+1),rlist(ifgamm),im,
     $       rtastk(isp1+2),rtastk(isp1+3),1,
     $       rlist(kam+1),kax1,kx,psix,psiy,psix0,psiy0,
     $       periodic,normalmode,angle,xplane)
      elseif(ktflistq(ktastk(isp1+2)))then
        if(ktfnonlistq(ktastk(isp1+3)))then
          irtc=itfmessage(9,'General::wrongtype','"List for #3"')
          return
        endif
        ka=ktfaddr(ktastk(isp1+2))
        kaa=ktfaddr(klist(ka+1))
        m=ilist(2,ka-1)
c       print *,'length=',m
        ka1=ktfaddr(ktastk(isp1+3))
        kaa1=ktfaddr(klist(ka1+1))
        if(ilist(2,ka1-1) .lt. m)then
          irtc=itfmessage(9,'General::equalleng','"#2 and #3"')
          return
        endif
        kax=ktavaloc(-1,m)
        call pgrmat1(ilist(1,ilattp+1),rlist(ifgamm),im,rlist(kaa+1),
     $       rlist(kaa1+1),m,rlist(kam+1),kax1,rlist(kax+1),psix,psiy,
     $       psix0,psiy0,periodic,normalmode,angle,xplane)
        kx=ktflist+kax
      else
        kx=ktfoper+mtfnull
      endif
      return
      end

      subroutine pgrmat1 (latt,gammab,im,st,stt,ns,twsi,kax1,ar,psix,
     $     psiy,psix0,psiy0,periodic,normalmode,angle,xplane)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 kax1
      logical*4 periodic,normalmode,angle,xplane,entrance
      integer*4 latt(2,nlat),im,js,ns,i,i1,iaxj,ias
      real*8 gammab(nlat),twsi(27),st(ns),stt(ns),psix,psiy,ar(ns),
     $     psix0,psiy0,axi,bxi,ayi,byi,ri11,ri12,ri21,ri22,pinx,piny,
     $     fx,fy,rotation,cost,sint,r11,r12,r21,r22,det,qr,um,
     $     tx,tpx,ty,tpy,axj,bxj,dphix,ayj,byj,dphiy,phx,phy,
     $     sx,cx,sy,cy,sg,v,x,px,y,py

      axi=twsi(1)
      bxi=twsi(2)
      ayi=twsi(4)
      byi=twsi(5)
      ri11=twsi(11)
      ri12=twsi(12)
      ri21=twsi(13)
      ri22=twsi(14)

      if(.not.periodic) go to 199

c----- Ring ------------------------------------------------------------
      pinx=0.5d0*psix
      piny=0.5d0*psiy
      fx=0.5d0/sin(pinx)
      fy=0.5d0/sin(piny)
      do 101 i=1,ns
        js=st(i)
        if(idtype(latt(1,js)) .eq. icbend) then
          rotation=-rlist(idval(latt(1,js))+5)
          cost=cos(rotation)
          sint=sin(rotation)
        else
          if(xplane) then
            cost=1.d0
            sint=0.d0
          else
            cost=0.d0
            sint=1.d0
          endif
        endif
        ar(i)=0.d0
        do 100 i1=1,2
          if(i1.eq.2) js=stt(i)+1
          entrance=im.eq.js .and. i1.eq.1
          iaxj=ilist(2,kax1+js)
          ias=ilist(2,iaxj+1)
        
          r11=rlist(ias+11)
          r12=rlist(ias+12)
          r21=rlist(ias+13)
          r22=rlist(ias+14)
c     write(*,'(a,1p,4e15.7)') 'rj=',r11,r12,r21,r22
          det=r11*r22-r12*r21
          if(det.gt.1d0) then
            qr=sqrt((det-1d0)/det)
            um=sqrt(det)
            tx =  r22*qr*cost
            tpx= -r21*qr*cost +    um *sint
            ty =                r22*qr*sint
            tpy=  um    *cost + r12*qr*sint
          else
            um=sqrt(1d0-det)
            tx =            r12*sint
            tpx=  um *cost -r11*sint
            ty =  r12*cost
            tpy=  r22*cost + um*sint
          endif
c     write(*,'(a,1p,4e15.7)') 't=',tx,tpx,ty,tpy
          axj=rlist(ias+1)
          bxj=rlist(ias+2)
          dphix=twsi(3)-rlist(ias+3)

          ayj=rlist(ias+4)
          byj=rlist(ias+5)
          dphiy=twsi(6)-rlist(ias+6)

c     write(*,'(a,1p,6e15.7)')'twsj=',axj,bxj,rlist(ias+3),ayj,byj,twsj(6)
c     write(*,'(a,1p,6e15.7)')'twsi=',axi,bxi,twsi(3),ayi,byi,twsi(6)
          phx=pinx - abs(dphix)*psix/psix0
          phy=piny - abs(dphiy)*psiy/psiy0
          sx=sin(phx)
          cx=cos(phx)
          sy=sin(phy)
          cy=cos(phy)
          if( entrance )then
            sg=-1.d0
          else
            sg=sign(1,im-js)
          endif
          if( normalmode )then
            if( angle )then
              if( xplane )then
                v=fx * (
     $               -tx /sqrt(bxi*bxj) *
     $               ( (axi-axj)*sg*sx + (1.d0+axi*axj)*cx )
     $               + tpx * sqrt(bxj/bxi) *
     $               ( sg*sx - axi*cx )
     $               )
              else
                v=fy * (
     $               -ty /sqrt(byi*byj) *
     $               ( (ayi-ayj)*sg*sy + (1.d0+ayi*ayj)*cy )
     $               + tpy * sqrt(byj/byi) *
     $               ( sg*sy - ayi*cy )
     $               )
              endif                      
            else
              if( xplane )then
                v=fx * (
     $               tx * sqrt(bxi/bxj)* ( sg*sx + axj*cx )
     $               + tpx * sqrt(bxi*bxj) * cx
     $               )
              else
                v=fy * (
     $               ty * sqrt(byi/byj)* ( sg*sy + ayj*cy )
     $               + tpy * sqrt(byi*byj) * cy
     $               )
              endif
            endif
          else
            x=fx * (
     $           tx * sqrt(bxi/bxj)* ( sg*sx + axj*cx )
     $           + tpx * sqrt(bxi*bxj) * cx
     $           )
            px=fx * (
     $           -tx /sqrt(bxi*bxj) *
     $           ( (axi-axj)*sg*sx + (1.d0+axi*axj)*cx )
     $           + tpx * sqrt(bxj/bxi) *
     $           ( sg*sx - axi*cx )
     $           )
            y=fy * (
     $           ty * sqrt(byi/byj)* ( sg*sy + ayj*cy )
     $           + tpy * sqrt(byi*byj) * cy
     $           )
            py=fy * (
     $           -ty /sqrt(byi*byj) *
     $           ( (ayi-ayj)*sg*sy + (1.d0+ayi*ayj)*cy )
     $           + tpy * sqrt(byj/byi) *
     $           ( sg*sy - ayi*cy )
     $           )
c     write(*,'(a,1p,4e15.7)') 'xnormal=',x,px,y,py
c     write(*,'(a,1p,4e15.7)') 'ri=',r11,r12,r21,r22
            det=ri11*ri22-ri12*ri21
            if(det.gt.1d0) then
              qr=sqrt((det-1d0)/det)
              um=sqrt(det)
              if( angle )then
                if( xplane )then
                  v= (-ri11*x  -ri12*px)*qr +              um*py
                else
                  v=             um*px      + (-ri11*y  +ri21*py)*qr
                endif
              else
                if( xplane )then
                  v= (-ri21*x  -ri22*px)*qr +      um*y
                else
                  v=    um*x                + ( ri12*y  -ri22*py)*qr
                endif
              endif
c            x(1)= (-ri21*x  -ri22*px)*qr +      um*y
c            x(2)= (-ri11*x  -ri12*px)*qr +              um*py
c            x(3)=    um*x                + ( ri12*y  -ri22*py)*qr
c            x(4)=             um*px      + (-ri11*y  +ri21*py)*qr
            else
              um=sqrt(1d0-det)
              if( angle )then
                if( xplane )then
                  v=           um*px     -ri21*y  +ri11*py
                else
                  v= -ri21*x  -ri22*px             + um*py
                endif
              else
                if( xplane )then
                  v=   um*x              +ri22*y  -ri12*py
                else
                  v= -ri11*x  -ri12*px     + um*y
                endif
              endif
c            x(1)=   um*x              +ri22*y  -ri12*py
c            x(2)=           um*px     -ri21*y  +ri11*py
c            x(3)= -ri11*x  -ri12*px     + um*y
c            x(4)= -ri21*x  -ri22*px             + um*py
            endif
c     write(*,'(a,1p,4e15.7)') 'response=',v
          endif
          ar(i)=ar(i)+v
 100    continue
        ar(i)=-0.5d0*ar(i)
        if(trpt)then
          ar(i)=ar(i)*sqrt(gammab(js)/gammab(im))
        endif
 101  continue
      return

 199  continue
c----- Transport line --------------------------------------------------
      do 201 i=1,ns
        js=st(i)
        if(idtype(latt(1,js)) .eq. icbend) then
          rotation=-rlist(idval(latt(1,js))+5)
          cost=cos(rotation)
          sint=sin(rotation)
        else
          if(xplane) then
            cost=1.d0
            sint=0.d0
          else
            cost=0.d0
            sint=1.d0
          endif
        endif
        ar(i)=0.d0
        do 200 i1=1,2
          if(i1.eq.2) js=stt(i)+1
          entrance=im.eq.js .and. i1.eq.1
          iaxj=ilist(2,kax1+js)
          ias=ilist(2,iaxj+1)
        
          r11=rlist(ias+11)
          r12=rlist(ias+12)
          r21=rlist(ias+13)
          r22=rlist(ias+14)
c     write(*,'(a,1p,4e15.7)') 'rj=',r11,r12,r21,r22
          det=r11*r22-r12*r21
          if(det.gt.1d0) then
            qr=sqrt((det-1d0)/det)
            um=sqrt(det)
            tx =  r22*qr*cost
            tpx= -r21*qr*cost +    um *sint
            ty =                r22*qr*sint
            tpy=  um    *cost + r12*qr*sint
          else
            um=sqrt(1d0-det)
            tx =            r12*sint
            tpx=  um *cost -r11*sint
            ty =  r12*cost
            tpy=  r22*cost + um*sint
          endif
c     write(*,'(a,1p,4e15.7)') 't=',tx,tpx,ty,tpy
          axj=rlist(ias+1)
          bxj=rlist(ias+2)
          dphix=twsi(3)-rlist(ias+3)

          ayj=rlist(ias+4)
          byj=rlist(ias+5)
          dphiy=twsi(6)-rlist(ias+6)

          sx=sin(dphix)
          cx=cos(dphix)
          sy=sin(dphiy)
          cy=cos(dphiy)
          if( im .lt. js .or. entrance ) then
            v=0.d0
            go to 200
          endif
          if( normalmode )then
            if( angle )then
              if( xplane )then
                v= -tx /sqrt(bxi*bxj) *
     $               ( (axi-axj)*cx + (1.d0+axi*axj)*sx )
     $               + tpx * sqrt(bxj/bxi) *( cx - axi*sx )
              else
                v= -ty /sqrt(byi*byj) *
     $               ( (ayi-ayj)*cy + (1.d0+ayi*ayj)*sy )
     $               + tpy * sqrt(byj/byi) *( cy - ayi*sy )
              endif                      
            else
              if( xplane )then
                v= tx * sqrt(bxi/bxj) * ( cx + axj*sx )
     $               + tpx * sqrt(bxi*bxj) * sx
              else
                v= ty * sqrt(byi/byj) * ( cy + ayj*sy )
     $               + tpy * sqrt(byi*byj) * sy
              endif
            endif
          else
            x= tx * sqrt(bxi/bxj) * ( cx + axj*sx )
     $           + tpx * sqrt(bxi*bxj) * sx
            px= -tx /sqrt(bxi*bxj) *
     $           ( (axi-axj)*cx + (1.d0+axi*axj)*sx )
     $           + tpx * sqrt(bxj/bxi) *( cx - axi*sx )
            y= ty * sqrt(byi/byj) * ( cy + ayj*sy )
     $           + tpy * sqrt(byi*byj) * sy
            py= -ty /sqrt(byi*byj) *
     $           ( (ayi-ayj)*cy + (1.d0+ayi*ayj)*sy )
     $           + tpy * sqrt(byj/byi) *( cy - ayi*sy )
            det=ri11*ri22-ri12*ri21
            if(det.gt.1d0) then
              qr=sqrt((det-1d0)/det)
              um=sqrt(det)
              if( angle )then
                if( xplane )then
                  v= (-ri11*x  -ri12*px)*qr +              um*py
                else
                  v=             um*px      + (-ri11*y  +ri21*py)*qr
                endif
              else
                if( xplane )then
                  v= (-ri21*x  -ri22*px)*qr +      um*y
                else
                  v=    um*x                + ( ri12*y  -ri22*py)*qr
                endif
              endif
            else
              um=sqrt(1d0-det)
              if( angle)then
                if( xplane )then
                  v=           um*px     -ri21*y  +ri11*py
                else
                  v= -ri21*x  -ri22*px             + um*py
                endif
              else
                if( xplane )then
                  v=   um*x              +ri22*y  -ri12*py
                else
                  v= -ri11*x  -ri12*px     + um*y
                endif
              endif
            endif
          endif
          ar(i)=ar(i)+v
 200    continue
        ar(i)=-0.5d0*ar(i)
        if(trpt)then
          ar(i)=ar(i)*sqrt(gammab(js)/gammab(im))
        endif
 201  continue
      return
      end
