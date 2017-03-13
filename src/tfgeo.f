      subroutine tfgeo(calgeo0)
      use kyparam
      use tfstk
      use ffs
      use ffs_flag
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 i,lxp,lt,nv,j,it
      real*8 geo1(3,4),vsave(256),dpos,offset,xp,fr,pos0,
     $     dox,doy,doz,g1,tffsmarkoffset,g2,rgetgl1,poso,circ,
     $     posi
      logical*4 calgeo,calgeo0,calpol0,chg
      call tfsetparam
      if(.not. geocal)then
        if(gammab(1) .ne. p0)then
          gammab(1)=p0
          call tfgeo1(1,nlat,calgeo,.true.)
          go to 9000
        endif
        return
      endif
      irad=6
      calpol0=calpol
      calpol=.false.
      calgeo=calgeo0
      if(.not. calgeo)then
        do 1030 i=1,nlat-1
          if(idtypec(i) .eq. 20)then
            calgeo=.true.
            goto 1040
          endif
1030    continue
      endif
1040  if(calgeo)then
        geo(:,:,1)=0.d0
        geo(2,1,1)=-1.d0
        geo(3,2,1)=-1.d0
        geo(1,3,1)= 1.d0
      endif
      pos(1)=0.d0
      gammab(1)=p0
      inext=0
      iprev=0
      call tfgeo1(1,nlat,calgeo,.false.)
      do i=1,nlat-1
        it=idtypec(i)
        if(it .eq. icMARK)then
          offset=tffsmarkoffset(i)
          if(offset .ne. 0.d0)then
            xp=max(1.d0,min(dble(nlat),offset+dble(i)))
            lxp=int(xp)
            fr=xp-lxp
            if(fr .eq. 0.d0)then
              geo(:,:,i)=geo(:,:,lxp)
              pos(i)=pos(lxp)
            else
              lt=idtypec(lxp)
              nv=kytbl(kwmax,lt)-1
              if(ideal)then
                call tmov(rlist(idvalc(lxp)+1),vsave,nv)
              else
                call tmov(rlist(latt(lxp)+1),vsave,nv)
              endif
              call qfraccomp(lxp,0.d0,fr,ideal,chg)
              if(chg)then
                call tmov(geo(1,1,lxp+1),geo1,12)
                pos0=pos(lxp+1)
                call tfgeo1(lxp,lxp+1,calgeo,.false.)
                if(i .ne. lxp+1)then
                  call tmov(geo(1,1,lxp+1),geo(1,1,i),12)
                  pos(i)=pos(lxp+1)
                  call tmov(geo1,geo(1,1,lxp+1),12)
                  pos(lxp+1)=pos0
                endif
              else
                call tmov(geo(1,1,lxp+1),geo(1,1,i),12)
                pos(i)=pos(lxp+1)
              endif
              if(ideal)then
                call tmov(vsave,rlist(idvalc(lxp)+1),nv)
              else
                call tmov(vsave,rlist(latt(lxp)+1),nv)
              endif
              if(i .eq. 1)then
                dpos=pos(1)
                pos=pos-dpos
              endif
              if(i .eq. nlat-1)then
                call tmov(geo(1,1,i),geo(1,1,nlat),12)
                pos(nlat)=pos(i)
              endif
            endif
          endif
        else
          j=i+1
          posi=pos(j)
          do while(.true.)
            if(pos(j) .ne. posi)then
              exit
            elseif(j .eq. nlat)then
              j=0
              posi=pos(1)
            elseif(idtypec(j) .eq. it)then
              inext(i)=j
              iprev(j)=i
              exit
            endif
            j=j+1
          enddo
        endif
      enddo
      if(sorg)then
        poso=pos(iorgr)
        pos=pos-poso
      endif
      if(calgeo)then
        dox=geo(1,4,iorgr)+ffv%geo0(1,4)*geo(1,3,iorgr)
     $       -ffv%geo0(2,4)*geo(1,1,iorgr)-ffv%geo0(3,4)*geo(1,2,iorgr)
        doy=geo(2,4,iorgr)+ffv%geo0(1,4)*geo(2,3,iorgr)
     $       -ffv%geo0(2,4)*geo(2,1,iorgr)-ffv%geo0(3,4)*geo(2,2,iorgr)
        doz=geo(3,4,iorgr)+ffv%geo0(1,4)*geo(3,3,iorgr)
     $       -ffv%geo0(2,4)*geo(3,1,iorgr)-ffv%geo0(3,4)*geo(3,2,iorgr)
        do 110 i=1,nlat
          geo(1,4,i)=geo(1,4,i)-dox
          geo(2,4,i)=geo(2,4,i)-doy
          geo(3,4,i)=geo(3,4,i)-doz
110     continue
        do 100 i=1,3
          geo1(1,i)=ffv%geo0(1,i)*geo(1,3,iorgr)
     $         -ffv%geo0(2,i)*geo(1,1,iorgr)
     1             -ffv%geo0(3,i)*geo(1,2,iorgr)
          geo1(2,i)=ffv%geo0(1,i)*geo(2,3,iorgr)
     $         -ffv%geo0(2,i)*geo(2,1,iorgr)
     1             -ffv%geo0(3,i)*geo(2,2,iorgr)
          geo1(3,i)=ffv%geo0(1,i)*geo(3,3,iorgr)
     $         -ffv%geo0(2,i)*geo(3,1,iorgr)
     1             -ffv%geo0(3,i)*geo(3,2,iorgr)
100     continue
        do 120 j=1,4
          do 130 i=1,nlat
            g1=geo(1,j,i)
            g2=geo(2,j,i)
            geo(1,j,i)=geo1(1,1)*g1+geo1(2,1)*g2+geo1(3,1)*geo(3,j,i)
            geo(2,j,i)=geo1(1,2)*g1+geo1(2,2)*g2+geo1(3,2)*geo(3,j,i)
            geo(3,j,i)=geo1(1,3)*g1+geo1(2,3)*g2+geo1(3,3)*geo(3,j,i)
130       continue
120     continue
        if(sorg)then
          poso=pos(iorgr)
          pos=pos-poso
        endif
      endif
      calpol=calpol0
 9000 pgev=rgetgl1('MOMENTUM')
      call tphyzp
      circ=pos(nlat)-pos(1)
      rlist(elatt%aux+1)=circ
      dleng=rgetgl1('FSHIFT')*circ
      if(circ .ne. 0.d0)then
        omega0=pi2*c*p0/h0/circ
        if(omega0 .eq. 0.d0)then
          write(*,*)'Design orbit length =',circ
        endif
        call rsetgl1('OMEGA0',omega0)
      endif
      return
      end

      subroutine tfgeo1(istart0,istop,calgeo,acconly)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_main
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 istart,istop,ke,ke1,i,k,i1,istart0
      integer*8 id
      real*8 p1,h1,ali,v,zetau,b,a,xiu,dchi3,coschi,sinchi,
     $     x1,x2,x3,y1,y2,y3,rho0,sp0,cp0,r1,r2,cchi1,schi1,
     $     cchi2,schi2,cchi3,schi3,dx,dy,dz,r11,r12,r13,
     $     r31,r32,r33,oneev,etau,phi,fb1,fb2,
     $     theta,cost,sint,r21,r22,r23,ald,tfacc
      parameter (oneev=1.d0+3.83d-12)
      logical*4 sol,calgeo,dir,acconly
c     begin initialize for preventing compiler warning
      y1=0.d0
      y2=0.d0
      y3=0.d0
      cost=1.d0
      sint=0.d0
c     end   initialize for preventing compiler warning
      sol=.false.
      istart=istart0
 1    ke=0
      ke1=0
      p1=gammab(istart)
      h1=p2h(p1)
c      h1=p1*sqrt(1.d0+1.d0/p1**2)
c      h1=sqrt(1.d0+p1**2)
      do i=istart,istop-1
        k=idtypec(i)
        if(trpt .and. (k .eq. icCAVI .or. k .eq. icMULT))then
          p1=tfacc(i,p1,h1,.true.)
        endif
        gammab(i+1)=p1
      enddo
      if(acconly)then
        return
      endif
      dvfs=0.d0
      call tesetdv(0.d0)
      do i=istart,istop-1
        i1=i+1
        if(sol)then
          sol=i .lt. ke
          go to 10
        endif
        k=idtypec(i)
        if(ideal)then
          id=idvalc(i)
        else
          id=elatt%comp(i)
        endif
        call loc_comp(id,cmp)
        if(kytbl(kwANGL,k) .ne. 0)then
          ali=cmp%value(kytbl(kwL,k))
          phi=cmp%value(kytbl(kwANGL,k))
          if(cmp%value(kytbl(kwFRMD,k)) .ne. 0.d0)then
            if(phi*ali .ne. 0.d0)then
              if(k .eq. icBEND)then
                fb1=cmp%value(kytbl(kwF1,icBEND))
     $               +cmp%value(kytbl(kwFB1,icBEND))
                fb2=cmp%value(kytbl(kwF1,icBEND))
     $               +cmp%value(kytbl(kwFB2,icBEND))
              else
                fb1=cmp%value(kytbl(kwFB1,k))
                fb2=cmp%value(kytbl(kwFB2,k))
              endif
              if(fb1 .ne. 0.d0 .or. fb2 .ne. 0.d0)then
                ali=ali-((phi*fb1)**2+(phi*fb2)**2)/ali/48.d0
     $               *sin(.5d0*(phi*(1.d0-cmp%value(kytbl(kwE1,k))
     $               -cmp%value(kytbl(kwE2,k)))
     $               -cmp%value(kytbl(kwAE1,k))
     $               -cmp%value(kytbl(kwAE2,k))))
     $               /sin(.5d0*phi)
              endif
            endif
          endif
        elseif(k .eq. icSOL)then
          call tsgeo(i,ke,ke1,sol)
          go to 10
        elseif(kytbl(kwL,k) .ne. 0)then
          ali=cmp%value(kytbl(kwL,k))
        else
          ali=0.d0
        endif
        pos(i1)=pos(i)+ali
        if(calgeo)then
          if(k .eq. icMARK)then
            if(ke .ne. 0 .and. ke1 .ne. 0)then
              if(rlist(idvalc(i)+ky_GEO_MARK) .ne. 0)then
                zetau=geo(1,3,ke1)*geo(1,1,i)+geo(2,3,ke1)*geo(2,1,i)
     1               +geo(3,3,ke1)*geo(3,1,i)
                b=min(1.d0,max(-1.d0,-zetau/geo(3,2,i)))
                a=1.d0+sqrt1(-b**2)
c                a=sqrt((1.d0-b)*(1.d0+b))
                 etau=geo(1,2,ke1)*geo(1,1,i)+geo(2,2,ke1)*geo(2,1,i)
     1               +geo(3,2,ke1)*geo(3,1,i)
                  xiu=geo(1,1,ke1)*geo(1,1,i)+geo(2,1,ke1)*geo(2,1,i)
     1               +geo(3,1,ke1)*geo(3,1,i)
                dchi3=atan2(a*etau-b*geo(3,3,ke1)*xiu,
     1                      a*xiu+b*geo(3,3,ke1)*etau)
                rlist(latt(ke)+ky_CHI3_SOL)=
     $               rlist(latt(ke)+ky_CHI3_SOL)+dchi3
                coschi= cos(dchi3)
                sinchi=-sin(dchi3)
                call trotg(geo(1,1,ke1),geo(1,3,ke1),coschi,sinchi)
                istart=ke1
                go to 1
              endif
            endif
          endif
          if(k .ge. 36)then
            call tmov(geo(1,1,i),geo(1,1,i1),12)
          endif
          if(kytbl(kwANGL,k) .ne. 0)then
            v=cmp%value(kytbl(kwANGL,k))
            ald=ali
            if(v .eq. 0.d0)then
              geo(1,4,i1)=geo(1,3,i)*ald+geo(1,4,i)
              geo(2,4,i1)=geo(2,3,i)*ald+geo(2,4,i)
              geo(3,4,i1)=geo(3,3,i)*ald+geo(3,4,i)
              call tmov(geo(1,1,i),geo(1,1,i1),9)
            else
              if(kytbl(kwROT,k) .ne. 0)then
                theta=cmp%value(kytbl(kwROT,k))
              else
                theta=0.d0
              endif
              if(theta .ne. 0.d0)then
                cost=cos(theta)
                sint=sin(theta)
                x1= cost*geo(1,1,i)-sint*geo(1,2,i)
                x2= cost*geo(2,1,i)-sint*geo(2,2,i)
                x3= cost*geo(3,1,i)-sint*geo(3,2,i)
                y1= sint*geo(1,1,i)+cost*geo(1,2,i)
                y2= sint*geo(2,1,i)+cost*geo(2,2,i)
                y3= sint*geo(3,1,i)+cost*geo(3,2,i)
              else
                x1= geo(1,1,i)
                x2= geo(2,1,i)
                x3= geo(3,1,i)
              endif
              rho0=ald/v
              sp0=sin(v)
              cp0=cos(v)
              r1=rho0*sp0
              if(cp0 .ge. 0.d0)then
                r2=rho0*sp0**2/(1.d0+cp0)
              else
                r2=rho0*(1.d0-cp0)
              endif
c              r2=2.d0*rho0*sin(v*.5d0)**2
              geo(1,4,i1)=geo(1,4,i)+(r1*geo(1,3,i)-r2*x1)
              geo(2,4,i1)=geo(2,4,i)+(r1*geo(2,3,i)-r2*x2)
              geo(3,4,i1)=geo(3,4,i)+(r1*geo(3,3,i)-r2*x3)
              geo(1,1,i1)= cp0*x1+sp0*geo(1,3,i)
              geo(1,3,i1)=-sp0*x1+cp0*geo(1,3,i)
              geo(2,1,i1)= cp0*x2+sp0*geo(2,3,i)
              geo(2,3,i1)=-sp0*x2+cp0*geo(2,3,i)
              geo(3,1,i1)= cp0*x3+sp0*geo(3,3,i)
              geo(3,3,i1)=-sp0*x3+cp0*geo(3,3,i)
              if(theta .ne. 0.d0)then
                geo(1,2,i1)=-sint*geo(1,1,i1)+cost*y1
                geo(1,1,i1)= cost*geo(1,1,i1)+sint*y1
                geo(2,2,i1)=-sint*geo(2,1,i1)+cost*y2
                geo(2,1,i1)= cost*geo(2,1,i1)+sint*y2
                geo(3,2,i1)=-sint*geo(3,1,i1)+cost*y3
                geo(3,1,i1)= cost*geo(3,1,i1)+sint*y3
              else
                geo(1,2,i1)=geo(1,2,i)
                geo(2,2,i1)=geo(2,2,i)
                geo(3,2,i1)=geo(3,2,i)
              endif
            endif
          elseif(k .eq. icMAP)then
            call tgmap(i)
          elseif(k .eq. icCOORD)then
            id=idvalc(i)
            cchi1=cos(cmp%value(ky_CHI1_COORD))
            schi1=sin(cmp%value(ky_CHI1_COORD))
            cchi2=cos(cmp%value(ky_CHI2_COORD))
            schi2=sin(cmp%value(ky_CHI2_COORD))
            cchi3=cos(cmp%value(ky_CHI3_COORD))
            schi3=sin(cmp%value(ky_CHI3_COORD))
            dir=cmp%value(ky_DIR_COORD) .eq. 0.d0
            if(dir)then
              dx=cmp%value(ky_DX_COORD)
              dy=cmp%value(ky_DY_COORD)
              dz=cmp%value(ky_DZ_COORD)
              r11= cchi1*cchi3+schi1*schi2*schi3
              r12=-cchi2*schi3
              r13= schi1*cchi3-cchi1*schi2*schi3
              r21=-schi1*schi2*cchi3+cchi1*schi3
              r22= cchi2*cchi3
              r23= cchi1*schi2*cchi3+schi1*schi3
              r31=-schi1*cchi2
              r32=-schi2
              r33= cchi1*cchi2
              geo(1,4,i1)=geo(1,4,i)
     1                  +dx*geo(1,1,i)+dy*geo(1,2,i)+dz*geo(1,3,i)
              geo(2,4,i1)=geo(2,4,i)
     1                  +dx*geo(2,1,i)+dy*geo(2,2,i)+dz*geo(2,3,i)
              geo(3,4,i1)=geo(3,4,i)
     1                  +dx*geo(3,1,i)+dy*geo(3,2,i)+dz*geo(3,3,i)
            else
              dx=-cmp%value(ky_DX_COORD)
              dy= cmp%value(ky_DY_COORD)
              dz= cmp%value(ky_DZ_COORD)
              r11= cchi1*cchi3+schi1*schi2*schi3
              r21= cchi2*schi3
              r31=-schi1*cchi3+cchi1*schi2*schi3
              r12= schi1*schi2*cchi3-cchi1*schi3
              r22= cchi2*cchi3
              r32= cchi1*schi2*cchi3+schi1*schi3
              r13= schi1*cchi2
              r23=-schi2
              r33= cchi1*cchi2
            endif
            geo(1,1,i1)=r11*geo(1,1,i)+r12*geo(1,2,i)+r13*geo(1,3,i)
            geo(2,1,i1)=r11*geo(2,1,i)+r12*geo(2,2,i)+r13*geo(2,3,i)
            geo(3,1,i1)=r11*geo(3,1,i)+r12*geo(3,2,i)+r13*geo(3,3,i)
            geo(1,2,i1)=r21*geo(1,1,i)+r22*geo(1,2,i)+r23*geo(1,3,i)
            geo(2,2,i1)=r21*geo(2,1,i)+r22*geo(2,2,i)+r23*geo(2,3,i)
            geo(3,2,i1)=r21*geo(3,1,i)+r22*geo(3,2,i)+r23*geo(3,3,i)
            geo(1,3,i1)=r31*geo(1,1,i)+r32*geo(1,2,i)+r33*geo(1,3,i)
            geo(2,3,i1)=r31*geo(2,1,i)+r32*geo(2,2,i)+r33*geo(2,3,i)
            geo(3,3,i1)=r31*geo(3,1,i)+r32*geo(3,2,i)+r33*geo(3,3,i)
            if(.not. dir)then
              geo(1,4,i1)=geo(1,4,i)
     1                  +dx*geo(1,1,i1)+dy*geo(1,2,i1)+dz*geo(1,3,i1)
              geo(2,4,i1)=geo(2,4,i)
     1                  +dx*geo(2,1,i1)+dy*geo(2,2,i1)+dz*geo(2,3,i1)
              geo(3,4,i1)=geo(3,4,i)
     1                  +dx*geo(3,1,i1)+dy*geo(3,2,i1)+dz*geo(3,3,i1)
            endif
          else
            geo(1,4,i1)=geo(1,3,i)*ali+geo(1,4,i)
            geo(2,4,i1)=geo(2,3,i)*ali+geo(2,4,i)
            geo(3,4,i1)=geo(3,3,i)*ali+geo(3,4,i)
            geo(1,1,i1)=geo(1,1,i)
            geo(2,1,i1)=geo(2,1,i)
            geo(3,1,i1)=geo(3,1,i)
            geo(1,2,i1)=geo(1,2,i)
            geo(2,2,i1)=geo(2,2,i)
            geo(3,2,i1)=geo(3,2,i)
            geo(1,3,i1)=geo(1,3,i)
            geo(2,3,i1)=geo(2,3,i)
            geo(3,3,i1)=geo(3,3,i)
          endif
        endif
c        write(*,*)'tfgeo1 ',i1,i1,k,ali,geo(
c     $       1,4,i1),geo(1,3,i),geo(1,4,i)
 10     continue
      enddo
      call tsetdvfs
      call tesetdv(0.d0)
      return
      end

      subroutine tfinitgeo
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      logical*4 err,geocal0
      if(geocal .or. omega0 .eq. 0.d0)then
        geocal0=geocal
        geocal=.true.
        call tffsa(0,kx,irtc)
        geocal=geocal0
      endif
      call tffssaveparams(0,-1,err)
      return
      end

      real*8 function tfacc(i,p1,h1,dir)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 i,id
      integer*8 ip
      real*8 harm,w,v,phic,dh,h2,p2,oneev,h1,p1
      logical*4 dir
      parameter (oneev=1.d0+3.83d-12)
      if(ideal)then
        ip=idvalc(i)
      else
        ip=latt(i)
      endif
      id=idtypec(i)
      if(dir)then
        v=abs(charge)*rlist(ip+kytbl(kwVOLT,id))/amass
      else
        v=-abs(charge)*rlist(ip+kytbl(kwVOLT,id))/amass
      endif
      if(v .ne. 0.d0)then
        harm=rlist(ip+kytbl(kwHARM,id))
        if(harm .eq. 0.d0)then
          w=pi2*rlist(ip+kytbl(kwFREQ,id))/c
        else
          w=omega0*harm/c
        endif
        phic=rlist(ip+kytbl(kwPHI,id))*sign(1.d0,charge)
        dh=max(oneev-h1,-v*sin(phic))
        h2=h1+dh
        p2=h2p(h2)
c        p2=h2*sqrt(1.d0-1.d0/h2**2)
c        p2=sqrt((h2-1.d0)*(h2+1.d0))
        p1=p1+dh*(h2+h1)/(p2+p1)
        h1=h2
      endif
      tfacc=p1
      return
      end

      real*8 function tfpos(k)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k
      tfpos=rlist(ifpos+k-1)
      return
      end

      subroutine tforbitgeo(og,geo,dx,dpx,dy,dpy)
      implicit none
      real*8 og(3,4),geo(3,4),g1,dx,dpx,dy,dpy
      real*8 chi2i,cchi2i,schi2i,chi1i,cchi1i,schi1i
      integer*4 j
      chi2i =-asin(min(1.d0,max(-1.d0,dpy)))
      cchi2i=cos(chi2i)
      schi2i=sin(chi2i)
      chi1i =-asin(min(1.d0,max(-1.d0,dpx/cchi2i)))
      cchi1i=cos(chi1i)
      schi1i=sin(chi1i)
      do j=1,3
        og(j,4)=geo(j,4)+geo(j,1)*dx+geo(j,2)*dy
        og(j,1)= cchi1i*geo(j,1)+schi1i*geo(j,3)
        g1     =-schi1i*geo(j,1)+cchi1i*geo(j,3)
        og(j,3)= cchi2i*g1-schi2i*geo(j,2)
        og(j,2)= schi2i*g1+cchi2i*geo(j,2)
      enddo
      return
      end

      subroutine tsetg(geo,chi)
      use ffs
      use tffitcode
      implicit none
      real*8 geo(3,3),chi(3),c1,c2,c3,s1,s2,s3
      c1=cos(pi/180.d0*chi(1))
      s1=sin(pi/180.d0*chi(1))
      c2=cos(pi/180.d0*chi(2))
      s2=sin(pi/180.d0*chi(2))
      c3=cos(pi/180.d0*chi(3))
      s3=sin(pi/180.d0*chi(3))
      geo(1,1)=c1*c2
      geo(2,1)=s1*c2
      geo(3,1)=s2
      geo(1,2)=-c3*s1+c1*s2*s3
      geo(2,2)=c1*c3+s1*s2*s3
      geo(3,2)=-c2*s3
      geo(1,3)=-c1*s2*c3-s1*s3
      geo(2,3)=-s1*s2*c3-c1*s3
      geo(3,3)= c2*c3
      return
      end
