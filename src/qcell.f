      subroutine qcell(ibegin,idp,
     1     hstab,vstab,tracex,tracey,fam,over)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ibegin,idp
      real*8 tracex,tracey
      logical*4 hstab,vstab,fam,over
      call qcell1(1,0.d0,nlat,0.d0,idp,
     1     hstab,vstab,tracex,tracey,fam,over,.true.,0)
      return
      end

      subroutine qcell1(ibegin,frbegin,iend,frend,
     $     idp,hstab,vstab,tracex,tracey,fam,over,
     $     chgini,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      real*8 bmin,bmax,amax
      integer*4 itmax
      parameter (bmin=1.d-16,bmax=1.d16,amax=1.d16)
      parameter (itmax=63)
      real*8 ftwiss(ntwissfun),tracex,tracey,frbegin,frend,
     $     tffselmoffset,r1,r2,r3,r4,c1,
     $     s11,s12,s13,s14,s21,s22,s23,s24,
     $     s31,s32,s33,s34,s41,s42,s43,s44,
     $     a11,a12,a21,a22,b11,b12,b21,b22,
     $     a13,a23,b13,b23,
     $     cosmux,cosmuy,sinmux,sinmuy,
     $     amux,amuy,dcosmux,dcosmuy,
     $     xb,xe,xp,fr,fra,frb,tr,
     $     dpsix,dpsiy,cosx,sinx,cosy,siny,
     $     x11,x22,y11,y22
      integer*4 ibegin,iend,idp,level,
     $     ie1,l,nm,lx,k,lfno
      logical*4 hstab,vstab,over,stab,codfnd,fam,chgini,pri
      real*8 trans(4,5),cod(6),
     $     tm11,tm12,tm13,tm14,tm15,
     $     tm21,tm22,tm23,tm24,tm25,
     $     tm31,tm32,tm33,tm34,tm35,
     $     tm41,tm42,tm43,tm44,tm45
      equivalence
     $     (trans(1,1),tm11),(trans(1,2),tm12),(trans(1,3),tm13),
     $     (trans(1,4),tm14),(trans(1,5),tm15),
     $     (trans(2,1),tm21),(trans(2,2),tm22),(trans(2,3),tm23),
     $     (trans(2,4),tm24),(trans(2,5),tm25),
     $     (trans(3,1),tm31),(trans(3,2),tm32),(trans(3,3),tm33),
     $     (trans(3,4),tm34),(trans(3,5),tm35),
     $     (trans(4,1),tm41),(trans(4,2),tm42),(trans(4,3),tm43),
     $     (trans(4,4),tm44),(trans(4,5),tm45)
      pri=.false.
      if(cell)then
        codfnd=fam
        cod=twiss(ibegin,idp,mfitdx:mfitddp)
        call qcod(idp,ibegin,frbegin,iend,frend,
     $       trans,cod,codfnd,over)
        if(orbitcal)then
          twiss(ibegin,idp,mfitdx:mfitddp )=cod
        endif
        if(over)then
          if(lfno .gt. 0)then
            if(.not. codfnd)then
              write(lfno,*)
     $             '*****qcod---> Overflow & closed orbit not found'
            else
              write(lfno,*)'*****qcod---> Overflow'
            endif
            pri=.true.
          endif
          hstab=.false.
          vstab=.false.
        elseif(.not. codfnd .and. lfno .gt. 0)then
          write(lfno,*)'*****qcod---> Closed orbit not found'
          pri=.true.
        endif
        call qmdiag(
     $       tm11,tm12,tm13,tm14,
     1       tm21,tm22,tm23,tm24,
     1       tm31,tm32,tm33,tm34,
     1       tm41,tm42,tm43,tm44,
     1       r1,r2,r3,r4,c1,level,stab,lfno)
C     ----------------------------
        if(.not. stab)then
          go to 1
        endif
        twiss(ibegin,idp,mfitr1) = r1
        twiss(ibegin,idp,mfitr2) = r2
        twiss(ibegin,idp,mfitr3) = r3
        twiss(ibegin,idp,mfitr4) = r4
        twiss(ibegin,idp,mfitdetr) = r1*r4-r2*r3
C------s : transformation matrix
C     ( c1*I  J.Transpose[R].J )
C     ( R     c1*I             )
        s11 = c1
        s12 = 0.d0
        s13 = -r4
        s14 = r2
        s21 = 0.d0
        s22 = c1
        s23 = r3
        s24 = -r1
        s31 = r1
        s32 = r2
        s33 = c1
        s34 = 0.d0
        s41 = r3
        s42 = r4
        s43 = 0.d0
        s44 = c1
C------A,B:diagonalized transfer matrix
        a11 = s11*(tm11*s11-tm13*s31-tm14*s41)
     1       + s13*(tm31*s11-tm33*s31-tm34*s41)
     1       + s14*(tm41*s11-tm43*s31-tm44*s41)
        a12 = s11*(tm12*s22-tm13*s32-tm14*s42)
     1       + s13*(tm32*s22-tm33*s32-tm34*s42)
     1       + s14*(tm42*s22-tm43*s32-tm44*s42)
        a21 = s22*(tm21*s11-tm23*s31-tm24*s41)
     1       + s23*(tm31*s11-tm33*s31-tm34*s41)
     1       + s24*(tm41*s11-tm43*s31-tm44*s41)
        a22 = s22*(tm22*s22-tm23*s32-tm24*s42)
     1       + s23*(tm32*s22-tm33*s32-tm34*s42)
     1       + s24*(tm42*s22-tm43*s32-tm44*s42)
        b11 = s31*(-tm11*s13-tm12*s23+tm13*s33)
     1       + s32*(-tm21*s13-tm22*s23+tm23*s33)
     1       + s33*(-tm31*s13-tm32*s23+tm33*s33)
        b12 = s31*(-tm11*s14-tm12*s24+tm14*s44)
     1       + s32*(-tm21*s14-tm22*s24+tm24*s44)
     1       + s33*(-tm31*s14-tm32*s24+tm34*s44)
        b21 = s41*(-tm11*s13-tm12*s23+tm13*s33)
     1       + s42*(-tm21*s13-tm22*s23+tm23*s33)
     1       + s44*(-tm41*s13-tm42*s23+tm43*s33)
        b22 = s41*(-tm11*s14-tm12*s24+tm14*s44)
     1       + s42*(-tm21*s14-tm22*s24+tm24*s44)
     1       + s44*(-tm41*s14-tm42*s24+tm44*s44)
C-----check off diagonal component
c     tmd13 = s11*(-tm11*s13-tm12*s23+tm13*s33+tm14*s43)
c     1        + s12*(-tm21*s13-tm22*s23+tm23*s33+tm24*s43)
c     1        + s13*(-tm31*s13-tm32*s23+tm33*s33+tm34*s43)
c     1        + s14*(-tm41*s13-tm42*s23+tm43*s33+tm44*s43)
c     tmd14 = s11*(-tm11*s14-tm12*s24+tm13*s34+tm14*s44)
c     1        + s12*(-tm21*s14-tm22*s24+tm23*s34+tm24*s44)
c     1        + s13*(-tm31*s14-tm32*s24+tm33*s34+tm34*s44)
c     1        + s14*(-tm41*s14-tm42*s24+tm43*s34+tm44*s44)
c     tmd23 = s21*(-tm11*s13-tm12*s23+tm13*s33+tm14*s43)
c     1        + s22*(-tm21*s13-tm22*s23+tm23*s33+tm24*s43)
c     1        + s23*(-tm31*s13-tm32*s23+tm33*s33+tm34*s43)
c     1        + s24*(-tm41*s13-tm42*s23+tm43*s33+tm44*s43)
c     tmd24 = s21*(-tm11*s14-tm12*s24+tm13*s34+tm14*s44)
c     1        + s22*(-tm21*s14-tm22*s24+tm23*s34+tm24*s44)
c     1        + s23*(-tm31*s14-tm32*s24+tm33*s34+tm34*s44)
c     1        + s24*(-tm41*s14-tm42*s24+tm43*s34+tm44*s44)
c     tmd31 = s31*(tm11*s11+tm12*s21-tm13*s31-tm14*s41)
c     1        + s32*(tm21*s11+tm22*s21-tm23*s31-tm24*s41)
c     1        + s33*(tm31*s11+tm32*s21-tm33*s31-tm34*s41)
c     1        + s34*(tm41*s11+tm42*s21-tm43*s31-tm44*s41)
c     tmd32 = s31*(tm11*s12+tm12*s22-tm13*s32-tm14*s42)
c     1        + s32*(tm21*s12+tm22*s22-tm23*s32-tm24*s42)
c     1        + s33*(tm31*s12+tm32*s22-tm33*s32-tm34*s42)
c     1        + s34*(tm41*s12+tm42*s22-tm43*s32-tm44*s42)
c     tmd41 = s41*(tm11*s11+tm12*s21-tm13*s31-tm14*s41)
c     1        + s42*(tm21*s11+tm22*s21-tm23*s31-tm24*s41)
c     1        + s43*(tm31*s11+tm32*s21-tm33*s31-tm34*s41)
c     1        + s44*(tm41*s11+tm42*s21-tm43*s31-tm44*s41)
c     tmd42 = s41*(tm11*s12+tm12*s22-tm13*s32-tm14*s42)
c     1        + s42*(tm21*s12+tm22*s22-tm23*s32-tm24*s42)
c     1        + s43*(tm31*s12+tm32*s22-tm33*s32-tm34*s42)
c     1        + s44*(tm41*s12+tm42*s22-tm43*s32-tm44*s42)
c     -deb
c     print *,'=======check off diagonal component======'
c     print *,'-------(1,2) part'
c     print *,tmd13,tmd14
c     print *,tmd23,tmd24
c     print *,'-------(2,1) part'
c     print *,tmd31,tmd32
c     print *,tmd41,tmd42
        tracex=a11+a22
        tracey=b11+b22
        hstab=tracex .ge. -2.d0 .and. tracex .le. 2.d0
        vstab=tracey .ge. -2.d0 .and. tracey .le. 2.d0
        if(hstab)then
          cosmux=.5d0*tracex
        else
          cosmux=2.d0/tracex
        endif
        if(chgini .or. hstab)then
          if(hstab)then
            sinmux=sign(sqrt(abs(-a12*a21-.25d0*(a11-a22)**2)),a12)
          else
            sinmux=sign(sqrt(1.d0-cosmux**2),a12)
          endif
          amux=atan2(sinmux,cosmux)
          dcosmux=2.d0*sin(.5d0*amux)**2
          twiss(ibegin,idp,mfitax)=.5d0*(a11-a22)/sinmux
          if(a12*a21 .lt. 0.d0)then
            twiss(ibegin,idp,mfitbx)=sqrt(-a12/a21
     $           *(1.d0+twiss(ibegin,idp,mfitax)**2))
          else
            twiss(ibegin,idp,mfitbx)=a12/sinmux
          endif
        else
          dcosmux=1.d0-cosmux
        endif
        if(vstab)then
          cosmuy=.5d0*tracey
        else
          cosmuy=2.d0/tracey
        endif
        if(chgini .or. vstab)then
          if(vstab)then
            sinmuy=sign(sqrt(abs(-b12*b21-.25d0*(b11-b22)**2)),b12)
          else
            sinmuy=sign(sqrt(1.d0-cosmuy**2),b12)
          endif
          amuy=atan2(sinmuy,cosmuy)
          dcosmuy=2.d0*sin(.5d0*amuy)**2
          twiss(ibegin,idp,mfitay)=.5d0*(b11-b22)/sinmuy
          if(b12*b21 .lt. 0.d0)then
            twiss(ibegin,idp,mfitby)=sqrt(-b12/b21
     $           *(1.d0+twiss(ibegin,idp,mfitay)**2))
          else
            twiss(ibegin,idp,mfitby)=b12/sinmuy
          endif
        else
          dcosmuy=1.d0-cosmuy
        endif
C--   deb
c     print *,'   new parameters------'
c     print *,'    axi,bxi  =',twiss(ibegin,idp,1),twiss(ibegin,idp,2)
c     print *,'    ayi,byi  =',twiss(ibegin,idp,4),twiss(ibegin,idp,5)
c     print *,'    r11,r12  =',twiss(ibegin,idp,11),twiss(ibegin,idp,12)
c     print *,'    r21,r22  =',twiss(ibegin,idp,13),twiss(ibegin,idp,14)
C========Dispersion ========
c     (Note) Disperdion is defined in 2*2 world
        a13=s11*tm15+s13*tm35+s14*tm45
        a23=s22*tm25+s23*tm35+s24*tm45
        b13=s31*tm15+s32*tm25+s33*tm35
        b23=s41*tm15+s42*tm25+s44*tm45
        if(chgini .or. hstab)then
          if( dcosmux.ne.0.d0 ) then
            twiss(ibegin,idp,mfitex) =
     $           0.5d0*(a13+a12*a23-a22*a13)/dcosmux
            twiss(ibegin,idp,mfitepx) =
     $           0.5d0*(a23+a21*a13-a11*a23)/dcosmux
          endif
        endif
        if(chgini .or. vstab)then
          if( dcosmuy.ne.0.d0 ) then
            twiss(ibegin,idp,mfitey) =
     $           0.5d0*(b13+b12*b23-b22*b13)/dcosmuy
            twiss(ibegin,idp,mfitepy)=
     $           0.5d0*(b23+b21*b13-b11*b23)/dcosmuy
          end if
        endif
      endif
 1    continue
      if(frbegin .gt. 0.d0)then
        call qtwissfrac1(ftwiss,trans,cod,idp,ibegin,frbegin,1.d0,
     $       .false.,.true.,over)
        if(.not. over)then
          call qtwiss(twiss,idp,ibegin+1,iend,over)
        endif
      else
        call qtwiss(twiss,idp,ibegin,iend,over)
      endif
      ie1=iend
      if(over)then
        if(lfno .gt. 0)then
          write(lfno,*)'***qtwiss---> Overflow'
          pri=.true.
        endif
        hstab=.false.
        vstab=.false.
      endif
      if(frend .gt. 0.d0)then
        call qtwissfrac1(ftwiss,trans,cod,idp,iend,0.d0,frend,
     $       .false.,.true.,over)
        if(over)then
          hstab=.false.
          vstab=.false.
        endif
        ie1=iend+1
      endif
      xb=dble(ibegin)
      xe=dble(min(iend+1,nlat))
      do l=ibegin+1,min(iend,nlat-1)
        nm=0
        if(idtype(latt(1,l)) .eq. icMARK)then
          xp=tffselmoffset(l)
          if(xp .ne. dble(l))then
c            write(*,*)'qcell ',l,iend,nlat,xp,xe
            if(xp .ge. xb .and. xp .lt. xe)then
              lx=int(xp)
              fr=xp-lx
 8101         if(fr .eq. 0.d0)then
                do k=1,ntwissfun
                  twiss(l,idp,k)=twiss(lx,idp,k)
                enddo
              else
                if(lx .eq. ibegin)then
                  fra=frbegin
                  frb=max(fr,fra)
                elseif(lx .eq. iend)then
                  fra=0.d0
                  frb=min(frend,fr)
                else
                  fra=0.d0
                  frb=fr
                endif
                if(fra .eq. frb)then
                  fr=0.d0
                  go to 8101
                endif
                call qtwissfrac1(ftwiss,
     $               tr,cod,idp,lx,fra,frb,
     $               .false.,.false.,over)
                twiss(l,idp,1:ntwissfun)=ftwiss
              endif
            endif
          endif
        endif
      enddo
      dpsix = twiss(ie1,idp,mfitnx) - twiss(ibegin,idp,mfitnx)
      dpsiy = twiss(ie1,idp,mfitny) - twiss(ibegin,idp,mfitny)
      cosx=cos(dpsix)
      sinx=sin(dpsix)
      cosy=cos(dpsiy)
      siny=sin(dpsiy)
      x11=sqrt(twiss(ie1,idp,mfitbx)/twiss(ibegin,idp,mfitbx))
     1     *(cosx+twiss(ibegin,idp,mfitax)*sinx)
      x22=sqrt(twiss(ibegin,idp,mfitbx)/twiss(ie1,idp,mfitbx))*
     1     (cosx-twiss(ie1,idp,mfitax)*sinx)
      y11=sqrt(twiss(ie1,idp,mfitby)/twiss(ibegin,idp,mfitby))
     1     *(cosy+twiss(ibegin,idp,mfitay)*siny)
      y22=sqrt(twiss(ibegin,idp,mfitby)/twiss(ie1,idp,mfitby))*
     1     (cosy-twiss(ie1,idp,mfitay)*siny)
      tracex=x11+x22
      tracey=y11+y22
      if(cell)then
        hstab=tracex .ge. -2.d0 .and. tracex .le. 2.d0
        vstab=tracey .ge. -2.d0 .and. tracey .le. 2.d0
        hstab=hstab .and. stab .and. (codfnd .or. fam)
        vstab=vstab .and. stab .and. (codfnd .or. fam)
      endif
      if(pri)then
        lfno=0
      endif
      return
      end
