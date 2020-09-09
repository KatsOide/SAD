      subroutine tfgeo(calgeo0)
      use kyparam
      use tfstk
      use ffs
      use ffs_flag
      use ffs_pointer
      use tffitcode
      use geolib
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 i,lxp,lt,nv,j,it
      real*8 geo1(3,4),geo2(3,3),vsave(256),dpos,offset,xp,fr,pos0,
     $     dgo(3),tffsmarkoffset,rgetgl1,poso,posi
      logical*4 calgeo,calpol0,chg
      logical*4 ,intent(in):: calgeo0
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
        do i=1,nlat-1
          if(idtypec(i) .eq. 20)then
            calgeo=.true.
            exit
          endif
        enddo
      endif
      if(calgeo)then
        geo(:,:,1)=geoini
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
              vsave(1:nv)=merge(rlist(idvalc(lxp)+1:idvalc(lxp)+nv),
     $             rlist(latt(lxp)+1:latt(lxp)+nv),ideal)
              call compelc(lxp,cmp)
              call qfraccomp(cmp,0.d0,fr,ideal,chg)
              if(chg)then
                geo1=geo(:,:,lxp+1)
                pos0=pos(lxp+1)
                call tfgeo1(lxp,lxp+1,calgeo,.false.)
                if(i .ne. lxp+1)then
                  geo(:,:,i)=geo(:,:,lxp+1)
                  pos(i)=pos(lxp+1)
                  geo(:,:,lxp+1)=geo1
                  pos(lxp+1)=pos0
                endif
              else
                geo(:,:,i)=geo(:,:,lxp+1)
                pos(i)=pos(lxp+1)
              endif
              if(ideal)then
                rlist(idvalc(lxp)+1:idvalc(lxp)+nv)=vsave(1:nv)
              else
                rlist(latt(lxp)+1:latt(lxp)+nv)=vsave(1:nv)
              endif
              if(i .eq. 1)then
                dpos=pos(1)
                pos=pos-dpos
              endif
              if(i .eq. nlat-1)then
                geo(:,:,nlat)=geo(:,:,i)
                pos(nlat)=pos(i)
              endif
            endif
          endif
        else
          j=i+1
          posi=pos(j)
          do
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
        dgo=geo(:,4,iorgr)+ffv%geo0(1,4)*geo(:,3,iorgr)
     $       -ffv%geo0(2,4)*geo(:,1,iorgr)-ffv%geo0(3,4)*geo(:,2,iorgr)
        do concurrent (i=1:nlat)
          geo(:,4,i)=geo(:,4,i)-dgo
        enddo
        geo2=matmul(ffv%geo0(:,1:3),transpose(geo(:,1:3,iorgr)))
        do concurrent (i=1:nlat)
          geo(:,1:4,i)=matmul(geo2,geo(:,1:4,i))
        enddo
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
      use tfcsi, only:icslfno
      use mathfun
      use geolib
      use kradlib, only:tallocvar,bsi
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 ,intent(in):: istart0,istop
      integer*4 istart,ke,ke1,i,k,i1
      integer*8 id
      real*8 p1,h1,ali,v,dchi3,
     $     rho0,sp0,cp0,r1,r2,cchi1,schi1,
     $     cchi2,schi2,cchi3,schi3,dx,dy,dz,r11,r12,r13,
     $     r31,r32,r33,oneev,phi,fb1,fb2,x(3),y(3),
     $     theta,cost,sint,r21,r22,r23,ald,tfacc
      parameter (oneev=1.d0+3.83d-12)
      logical*4 ,intent(in):: calgeo,acconly
      logical*4 sol,dir
      call tallocvar(bsi,1)
c     begin initialize for preventing compiler warning
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
          cycle
        endif
        k=idtypec(i)
        id=merge(idvalc(i),elatt%comp(i),ideal)
        if(k .eq. icSOL)then
          call tsgeo(i,ke,ke1,sol)
          cycle
        endif
        call compelc(i,cmp)
        ali=merge(cmp%value(kytbl(kwL,k)),0.d0,kytbl(kwL,k) .gt. 0)
        if(kytbl(kwANGL,k) .ne. 0)then
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
        endif
        pos(i1)=pos(i)+ali
        if(calgeo)then
          if(k .eq. icMARK)then
            if(ke .ne. 0 .and. ke1 .ne. 0)then
              if(rlist(idvalc(i)+ky_GEO_MARK) .ne. 0)then
                dchi3=tfchi(geo(:,:,i),3)
                geo(:,1:3,ke1)=matmul(matmul(
     $               tfderotgeo(geo(:,1:3,i),(/0.d0,0.d0,dchi3/)),
     $               transpose(geo(:,1:3,i))),geo(:,1:3,ke1))
                call compelc(ke,cmp)
                cmp%value(ky_CHI1_SOL:ky_CHI3_SOL)=
     $               tgrot(geo(:,1:3,ke),geo(:,1:3,ke1))
                istart=ke1
                go to 1
              endif
            endif
          endif
          if(k .ge. 36)then
            geo(:,:,i1)=geo(:,:,i)
          endif
          if(kytbl(kwANGL,k) .ne. 0)then
            v=cmp%value(kytbl(kwANGL,k))
            ald=ali
            if(v .eq. 0.d0)then
              geo(:,4,i1)=geo(:,3,i)*ald+geo(:,4,i)
              geo(:,1:3,i1)=geo(:,1:3,i)
            else
              theta=merge(cmp%value(kytbl(kwROT,k)),0.d0,
     $             kytbl(kwROT,k) .ne. 0)
              if(theta .ne. 0.d0)then
                cost=cos(theta)
                sint=sin(theta)
                x= cost*geo(:,1,i)-sint*geo(:,2,i)
                y= sint*geo(:,1,i)+cost*geo(:,2,i)
              else
                x= geo(:,1,i)
                y=0.d0
              endif
              rho0=ald/v
              sp0=sin(v)
              cp0=cos(v)
              r1=rho0*sp0
              r2=merge(rho0*sp0**2/(1.d0+cp0),rho0*(1.d0-cp0),
     $             cp0 .ge. 0.d0)
c              r2=2.d0*rho0*sin(v*.5d0)**2
              geo(:,4,i1)=geo(:,4,i)+(r1*geo(:,3,i)-r2*x)
              geo(:,1,i1)= cp0*x+sp0*geo(:,3,i)
              geo(:,3,i1)=-sp0*x+cp0*geo(:,3,i)
              if(theta .ne. 0.d0)then
                geo(:,2,i1)=-sint*geo(:,1,i1)+cost*y
                geo(:,1,i1)= cost*geo(:,1,i1)+sint*y
              else
                geo(:,2,i1)=geo(:,2,i)
              endif
            endif
          else
            select case (k)
            case (icMAP)
              call tgmap(i)
            case (icCOORD)
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
                geo(:,4,i1)=geo(:,4,i)
     1               +dx*geo(:,1,i)+dy*geo(:,2,i)+dz*geo(:,3,i)
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
              geo(:,1,i1)=r11*geo(:,1,i)+r12*geo(:,2,i)+r13*geo(:,3,i)
              geo(:,2,i1)=r21*geo(:,1,i)+r22*geo(:,2,i)+r23*geo(:,3,i)
              geo(:,3,i1)=r31*geo(:,1,i)+r32*geo(:,2,i)+r33*geo(:,3,i)
              if(.not. dir)then
                geo(:,4,i1)=geo(:,4,i)
     1               +dx*geo(:,1,i1)+dy*geo(:,2,i1)+dz*geo(:,3,i1)
              endif
            case default
              geo(:,4,i1)=geo(:,3,i)*ali+geo(:,4,i)
              geo(:,1:3,i1)=geo(:,1:3,i)
            end select
          endif
        endif
      enddo
      call tsetdvfs
      call tesetdv(0.d0)
      return
      end

      subroutine tfinitgeo
      use tfstk
      use ffs_flag, only:geocal
      use tmacro , only:omega0
      use tfcsi , only:lfni
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      logical*4 err,geocal0
      if(geocal .or. omega0 .eq. 0.d0)then
        geocal0=geocal
        geocal=.true.
        call tffsa(0,lfni,kx,irtc)
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
      use mathfun ,only:h2p
      implicit none
      integer*4 ,intent(in):: i
      integer*4 id
      integer*8 ip
      real*8 ,intent(inout):: h1,p1
      real*8 harm,w,v,phic,dh,h2,p2
      logical*4 ,intent(in):: dir
      real*8 , parameter :: oneev=1.d0+3.83d-12
      ip=merge(idvalc(i),latt(i),ideal)
      id=idtypec(i)
      v=merge(1.d0,-1.d0,dir)*
     $       abs(charge)*rlist(ip+kytbl(kwVOLT,id))/amass
      if(v .ne. 0.d0)then
        harm=rlist(ip+kytbl(kwHARM,id))
        w=merge(pi2*rlist(ip+kytbl(kwFREQ,id)),omega0*harm,
     $       harm .eq. 0.d0)/c
        phic=rlist(ip+kytbl(kwPHI,id))*sign(1.d0,charge)
        dh=max(oneev-h1,-v*sin(phic))
        h2=h1+dh
        p2=h2p(h2)
        p1=p1+dh*(h2+h1)/(p2+p1)
        h1=h2
      endif
      tfacc=p1
      return
      end
