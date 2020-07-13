      subroutine pgrmat(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer , only:gammab
      implicit none
      type(sad_descriptor) kx
      type (sad_rlist) , pointer :: kl
      type (sad_dlist) , pointer :: klopt,klopt1,klst,klstt
      integer*8 kax,kax1,kaxfirst,kaxlast,kam,kaa1,kaa
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
      if(tflistq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"non-List for #1"')
        return
      endif
      irtc=0
      im=itfloc(ktastk(isp1+1),irtc)
      if(irtc .ne. 0)then
        return
      endif

      if(tfnonlistq(dtastk(isp1+4),klopt))then
        irtc=itfmessage(9,'General::wrongtype','"List for #4"')
        return
      endif
      kax=ktfaddr(ktastk(isp1+4))
      if(tfnonlistq(klopt%dbody(1),klopt1))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List for #4[[1]]"')
        return
      endif
      kax1=ktfaddr(klopt%dbody(1)%k)
      kaxfirst=ktfaddr(klopt1%dbody(1)%k)
      kaxlast=ktfaddr(klopt1%dbody(klopt1%nl)%k)

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
        if(tfnonreallistq(dtastk(isp1+7),kl))then
          irtc=itfmessage(9,'General::wrongtype','"Real List for #7"')
          return
        endif
        if(kl%nl .ne. 2)then
          irtc=itfmessage(9,'General::wrongleng','"#7","2"')
          return
        endif
        psix=kl%rbody(1)*pi2
        psiy=kl%rbody(2)*pi2
      endif
c      write(*,*)'pgrmat-6 ',kaxfirst,kaxlast

c     dp00=rlist(ifirst+20)+1.d0
      psix0=rlist(kaxlast+3)-rlist(kaxfirst+3)
      psiy0=rlist(kaxlast+6)-rlist(kaxfirst+6)
c      write(*,*)'pgrmat-7 ',psix0,psiy0
      if(psix .eq. 0.d0)then
        psix=psix0
      endif
      if(psiy .eq. 0.d0)then
        psiy=psiy0
      endif

      kam=ktfaddr(klopt1%dbody(im)%k)
      if(ktfrealq(ktastk(isp1+2)))then
        if(ktfnonrealq(ktastk(isp1+3)))then
          irtc=itfmessage(9,'General::wrongtype','"Real for #3"')
          return
        endif
c      write(*,*)'pgrmat-8.1 '
        
        call pgrmat1(gammab,im,
     $       rtastk(isp1+2:isp1+2),rtastk(isp1+3:isp1+3),1,
     $       rlist(kam+1:kam+ntwissfun),klopt1,kx%x(1:1),
     $       psix,psiy,psix0,psiy0,
     $       periodic,normalmode,angle,xplane)
c      write(*,*)'pgrmat-8.2 '
      elseif(tflistq(ktastk(isp1+2),klst))then
        if(tfnonlistq(ktastk(isp1+3),klstt))then
          irtc=itfmessage(9,'General::wrongtype','"List for #3"')
          return
        endif
        kaa=sad_loc(klst%head)
        m=klst%nl
c       print *,'length=',m
        kaa1=sad_loc(klstt%head)
        if(klstt%nl .lt. m)then
          irtc=itfmessage(9,'General::equalleng','"#2 and #3"')
          return
        endif
        kax=ktavaloc(-1,m)
c      write(*,*)'pgrmat-9.2 ',kaa,kaa1,m,kam
        call pgrmat1(gammab,im,rlist(kaa+1:kaa+m),
     $       rlist(kaa1+1:kaa+m),m,rlist(kam+1:kam+ntwissfun),
     $       klopt1,rlist(kax+1:kax+m),psix,psiy,
     $       psix0,psiy0,periodic,normalmode,angle,xplane)
        kx%k=ktflist+kax
c      write(*,*)'pgrmat-9.3 '
      else
        kx%k=ktfoper+mtfnull
      endif
      return
      end

      subroutine pgrmat1 (gammab,im,st,stt,ns,twsi,klopt1,ar,psix,
     $     psiy,psix0,psiy0,periodic,normalmode,angle,xplane)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec,idvalc
      use kyparam
      implicit none
      type (sad_dlist) klopt1
      type (sad_rlist), pointer :: kls
      logical*4 periodic,normalmode,angle,xplane,entrance
      integer*4 im,js,ns,i,i1
      real*8 gammab(nlat),twsi(28),st(ns),stt(ns),psix,psiy,ar(ns),
     $     psix0,psiy0,axi,bxi,ayi,byi,ri11,ri12,ri21,ri22,pinx,piny,
     $     fx,fy,rotation,cost,sint,r11,r12,r21,r22,det,um,
     $     tx,tpx,ty,tpy,axj,bxj,dphix,ayj,byj,dphiy,phx,phy,
     $     sx,cx,sy,cy,sg,v,x,px,y,py,detr,detri

      axi=twsi(mfitax)
      bxi=twsi(mfitbx)
      ayi=twsi(mfitay)
      byi=twsi(mfitby)
      ri11=twsi(mfitr1)
      ri12=twsi(mfitr2)
      ri21=twsi(mfitr3)
      ri22=twsi(mfitr4)
      detri=twsi(mfitdetr)

      if(.not.periodic) go to 199

c----- Ring ------------------------------------------------------------
      pinx=0.5d0*psix
      piny=0.5d0*psiy
      fx=0.5d0/sin(pinx)
      fy=0.5d0/sin(piny)
      do 101 i=1,ns
        js=int(st(i))
        if(idtypec(js) .eq. icbend) then
          rotation=-rlist(idvalc(js)+ky_ROT_BEND)
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
          if(i1.eq.2) js=int(stt(i))+1
          entrance=im.eq.js .and. i1.eq.1
          if(.not. tfreallistq(klopt1%dbody(js),kls))then
            cycle
          endif
c          iaxj=klist(kax1+js)
c          call loc_rlist(ktfaddr(),klas)
          r11=kls%rbody(mfitr1)
          r12=kls%rbody(mfitr2)
          r21=kls%rbody(mfitr3)
          r22=kls%rbody(mfitr4)
c     write(*,'(a,1p,4e15.7)') 'rj=',r11,r12,r21,r22
          detr=kls%rbody(mfitdetr)
          det=r11*r22-r12*r21
          um=sqrt(1.d0-det)
          if(detr.gt.1d0) then
            tx =  r12*cost
            tpx= -r11*cost +    um *sint
            ty =                r12*sint
            tpy=  um    *cost + r22*sint
          else
            tx =            r12*sint
            tpx=  um *cost -r11*sint
            ty =  r12*cost
            tpy=  r22*cost + um*sint
          endif
c     write(*,'(a,1p,4e15.7)') 't=',tx,tpx,ty,tpy
          axj=kls%rbody(mfitax)
          bxj=kls%rbody(mfitbx)
          dphix=twsi(mfitnx)-kls%rbody(mfitnx)

          ayj=kls%rbody(mfitay)
          byj=kls%rbody(mfitby)
          dphiy=twsi(mfitny)-kls%rbody(mfitny)

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
            um=sqrt(1d0-det)
            if(detri.gt.1d0) then
              if( angle )then
                if( xplane )then
                  v= (-ri21*x  +ri11*px) +              um*py
                else
                  v=          um*px      + (-ri21*y  -ri22*py)
                endif
              else
                if( xplane )then
                  v= ( ri22*x  -ri12*px) +      um*y
                else
                  v=    um*x             + ( -ri12*y  -ri22*py)
                endif
              endif
c            x(1)= (-ri21*x  -ri22*px)*qr +      um*y
c            x(2)= (-ri11*x  -ri12*px)*qr +              um*py
c            x(3)=    um*x                + ( ri12*y  -ri22*py)*qr
c            x(4)=             um*px      + (-ri11*y  +ri21*py)*qr
            else
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
        js=int(st(i))
        if(idtypec(js) .eq. icbend) then
          rotation=-rlist(idvalc(js)+5)
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
          if(i1.eq.2) js=int(stt(i))+1
          entrance=im.eq.js .and. i1.eq.1
          if(.not. tfreallistq(klopt1%dbody(js),kls))then
            cycle
          endif
c          ias=ktfaddr(klopt1%body(js))
c          iayj=klist(kay1+js)
c          ias=klist(iayj+1)
        
          r11=kls%rbody(mfitr1)
          r12=kls%rbody(mfitr2)
          r21=kls%rbody(mfitr3)
          r22=kls%rbody(mfitr4)
c     write(*,'(a,1p,4e15.7)') 'rj=',r11,r12,r21,r22
          detr=kls%rbody(mfitdetr)
          det=r11*r22-r12*r21
          um=sqrt(1d0-det)
          if(detr.gt.1d0) then
            tx =  r12*cost
            tpx= -r11*cost +    um *sint
            ty =                r12*sint
            tpy=  um    *cost + r22*sint
          else
            tx =            r12*sint
            tpx=  um *cost -r11*sint
            ty =  r12*cost
            tpy=  r22*cost + um*sint
          endif
c     write(*,'(a,1p,4e15.7)') 't=',tx,tpx,ty,tpy
          axj=kls%rbody(mfitax)
          bxj=kls%rbody(mfitbx)
          dphix=twsi(mfitnx)-kls%rbody(mfitnx)

          ayj=kls%rbody(mfitay)
          byj=kls%rbody(mfitby)
          dphiy=twsi(mfitny)-kls%rbody(mfitny)

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
            um=sqrt(1d0-det)
            if(detri.gt.1d0) then
              if( angle )then
                if( xplane )then
                  v= (-ri21*x  +ri11*px) +              um*py
                else
                  v=          um*px      + (-ri21*y  -ri22*py)
                endif
              else
                if( xplane )then
                  v= ( ri22*x  -ri12*px) +      um*y
                else
                  v=    um*x             + ( -ri12*y  -ri22*py)
                endif
              endif
c            x(1)= (-ri21*x  -ri22*px)*qr +      um*y
c            x(2)= (-ri11*x  -ri12*px)*qr +              um*py
c            x(3)=    um*x                + ( ri12*y  -ri22*py)*qr
c            x(4)=             um*px      + (-ri11*y  +ri21*py)*qr
            else
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
