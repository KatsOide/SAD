      subroutine qcell(idp,optstat,fam)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit, only:ffs_stat
      use tffitcode
      implicit none
      type (ffs_bound) fb
      type (ffs_stat) optstat
      integer*4 idp,lfno
      logical*4 fam
      fb%lb=1
      fb%le=nlat
      fb%fb=0.d0
      fb%fe=0.d0
      lfno=0
      call qcell1(fb,idp,optstat,fam,.true.,lfno)
      return
      end

      subroutine qcell1(fbound,idp,optstat,fam,chgini,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit , only:ffs_stat
      use tffitcode
      implicit none
      type (ffs_bound) , intent(in)::fbound
      type (ffs_stat) , intent(out)::optstat
      real*8 bmin,bmax,amax
      integer*4 itmax
      parameter (bmin=1.d-16,bmax=1.d16,amax=1.d16)
      parameter (itmax=63)
      real*8 ftwiss(ntwissfun),
     $     tffselmoffset,r1,r2,r3,r4,c1,
     $     s11,s12,s13,s14,s21,s22,s23,s24,
     $     s31,s32,s33,s34,s41,s42,s43,s44,
     $     a11,a12,a21,a22,b11,b12,b21,b22,
     $     a13,a23,b13,b23,
     $     cosmux,cosmuy,sinmux,sinmuy,
     $     amux,amuy,dcosmux,dcosmuy,
     $     xb,xe,xp,fr,fra,frb,tr(4,5),
     $     dpsix,dpsiy,cosx,sinx,cosy,siny,
     $     x11,x22,y11,y22
      integer*4 , intent(in) :: idp
      integer*4 , intent(inout) :: lfno
      integer*4 ie1,l,nm,lx
      logical*4 , intent(in)::fam,chgini
      logical*4 stab,codfnd,pri,nanq
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
      if(calc6d)then
        call qcell6d(fbound,idp,optstat,lfno)
        return
      endif
      pri=.false.
      if(cell)then
        codfnd=fam
        cod=twiss(fbound%lb,idp,mfitdx:mfitddp)
        call qcod(idp,fbound,trans,cod,codfnd,optstat%over)
        if(orbitcal)then
          twiss(fbound%lb,idp,mfitdx:mfitddp)=cod
        endif
        if(optstat%over)then
          if(lfno .gt. 0)then
            if(.not. codfnd)then
              write(lfno,*)
     $             '*****qcod---> Overflow & closed orbit not found'
            else
              write(lfno,*)'*****qcod---> Overflow'
            endif
            pri=.true.
          endif
          optstat%stabx=.false.
          optstat%staby=.false.
        elseif(.not. codfnd .and. lfno .gt. 0)then
          write(lfno,*)'*****qcod---> Closed orbit not found'
          pri=.true.
        endif
        if(.not. calopt)then
          if(codfnd)then
            optstat%stabx=.true.
            optstat%staby=.true.
          else
            optstat%stabx=.false.
            optstat%staby=.false.
          endif
          return
        endif
        call qmdiag(
     $       tm11,tm12,tm13,tm14,
     1       tm21,tm22,tm23,tm24,
     1       tm31,tm32,tm33,tm34,
     1       tm41,tm42,tm43,tm44,
     1       r1,r2,r3,r4,c1,stab,nanq,lfno)
C     ----------------------------
        if(.not. stab)then
          optstat%stabx=.false.
          optstat%staby=.false.
          if(nanq)then
            return
          endif
          go to 1
        endif
        twiss(fbound%lb,idp,mfitr1) = r1
        twiss(fbound%lb,idp,mfitr2) = r2
        twiss(fbound%lb,idp,mfitr3) = r3
        twiss(fbound%lb,idp,mfitr4) = r4
        twiss(fbound%lb,idp,mfitdetr) = r1*r4-r2*r3
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
        optstat%tracex=a11+a22
        optstat%tracey=b11+b22
        optstat%stabx=.true.
        optstat%staby=.true.
        cosmux=.5d0*optstat%tracex
        if(cosmux .gt. 1.d0)then
          cosmux=1.d0/cosmux
          if(intres)then
            optstat%stabx=.false.
          else
            optstat%tracex=2.d0*cosmux
          endif
        else
        endif
        if(cosmux .lt. -1.d0)then
          cosmux=1.d0/cosmux
          if(halfres)then
            optstat%stabx=.false.
          else
            optstat%tracex=2.d0*cosmux
          endif
        endif
        cosmuy=.5d0*optstat%tracey
        if(cosmuy .gt. 1.d0)then
          cosmuy=1.d0/cosmuy
          if(intres)then
            optstat%staby=.false.
          else
            optstat%tracey=2.d0*cosmuy
          endif
        else
        endif
        if(cosmuy .lt. -1.d0)then
          cosmuy=1.d0/cosmuy
          if(halfres)then
            optstat%staby=.false.
          else
            optstat%tracey=2.d0*cosmuy
          endif
        endif
        sinmux=sign(sqrt(1.d0-cosmux**2),a12)
        if(chgini)then
          twiss(fbound%lb,idp,mfitax)=(a11-cosmux)/sinmux
          twiss(fbound%lb,idp,mfitbx)=a12/sinmux
        endif
        amux=atan2(sinmux,cosmux)
        dcosmux=2.d0*sin(.5d0*amux)**2
        sinmuy=sign(sqrt(1.d0-cosmuy**2),b12)
        if(chgini)then
          twiss(fbound%lb,idp,mfitay)=(b11-cosmuy)/sinmuy
          twiss(fbound%lb,idp,mfitby)=b12/sinmuy
        endif
        amuy=atan2(sinmuy,cosmuy)
        dcosmuy=2.d0*sin(.5d0*amuy)**2
C--   deb
c     print *,'   new parameters------'
c     print *,'    axi,bxi  =',twiss(fbound%lb,idp,1),twiss(fbound%lb,idp,2)
c     print *,'    ayi,byi  =',twiss(fbound%lb,idp,4),twiss(fbound%lb,idp,5)
c     print *,'    r11,r12  =',twiss(fbound%lb,idp,11),twiss(fbound%lb,idp,12)
c     print *,'    r21,r22  =',twiss(fbound%lb,idp,13),twiss(fbound%lb,idp,14)
C========Dispersion ========
c     (Note) Disperdion is defined in 2*2 world
        a13=s11*tm15+s13*tm35+s14*tm45
        a23=s22*tm25+s23*tm35+s24*tm45
        b13=s31*tm15+s32*tm25+s33*tm35
        b23=s41*tm15+s42*tm25+s44*tm45
        if(chgini .or. optstat%stabx)then
          if( dcosmux.ne.0.d0 ) then
            twiss(fbound%lb,idp,mfitex) =
     $           0.5d0*(a13+a12*a23-a22*a13)/dcosmux
            twiss(fbound%lb,idp,mfitepx) =
     $           0.5d0*(a23+a21*a13-a11*a23)/dcosmux
          endif
        endif
        if(chgini .or. optstat%staby)then
          if( dcosmuy.ne.0.d0 ) then
            twiss(fbound%lb,idp,mfitey) =
     $           0.5d0*(b13+b12*b23-b22*b13)/dcosmuy
            twiss(fbound%lb,idp,mfitepy)=
     $           0.5d0*(b23+b21*b13-b11*b23)/dcosmuy
          end if
        endif
      endif
 1    continue
      if(fbound%fb .gt. 0.d0)then
        call qtwissfrac1(ftwiss,trans,cod,idp,fbound%lb,fbound%fb,1.d0,
     $       .false.,.true.,optstat%over)
        if(.not. optstat%over)then
          call qtwiss(twiss,idp,fbound%lb+1,fbound%le,optstat%over)
        endif
      else
        call qtwiss(twiss,idp,fbound%lb,fbound%le,optstat%over)
      endif
      ie1=fbound%le
      if(optstat%over)then
        if(lfno .gt. 0)then
          write(lfno,*)'***qtwiss---> Overflow'
          pri=.true.
        endif
        optstat%stabx=.false.
        optstat%staby=.false.
      endif
      if(fbound%fe .gt. 0.d0)then
        call qtwissfrac1(ftwiss,trans,cod,idp,fbound%le,0.d0,fbound%fe,
     $       .false.,.true.,optstat%over)
        if(optstat%over)then
          optstat%stabx=.false.
          optstat%staby=.false.
        endif
        ie1=fbound%le+1
      endif
      xb=dble(fbound%lb)
      xe=dble(min(fbound%le+1,nlat))
      do l=fbound%lb+1,min(fbound%le,nlat-1)
        nm=0
        if(idtypec(l) .eq. icMARK)then
          xp=tffselmoffset(l)
          if(xp .ne. dble(l))then
            if(xp .ge. xb .and. xp .lt. xe)then
              lx=int(xp)
              fr=xp-lx
 8101         if(fr .eq. 0.d0)then
                twiss(l,idp,:)=twiss(lx,idp,:)
              else
                if(lx .eq. fbound%lb)then
                  fra=fbound%fb
                  frb=max(fr,fra)
                elseif(lx .eq. fbound%le)then
                  fra=0.d0
                  frb=min(fbound%fe,fr)
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
     $               .false.,.false.,optstat%over)
                twiss(l,idp,1:ntwissfun)=ftwiss
              endif
            endif
          endif
        endif
      enddo
      dpsix = twiss(ie1,idp,mfitnx) - twiss(fbound%lb,idp,mfitnx)
      dpsiy = twiss(ie1,idp,mfitny) - twiss(fbound%lb,idp,mfitny)
      cosx=cos(dpsix)
      sinx=sin(dpsix)
      cosy=cos(dpsiy)
      siny=sin(dpsiy)
      x11=sqrt(twiss(ie1,idp,mfitbx)/twiss(fbound%lb,idp,mfitbx))
     1     *(cosx+twiss(fbound%lb,idp,mfitax)*sinx)
      x22=sqrt(twiss(fbound%lb,idp,mfitbx)/twiss(ie1,idp,mfitbx))*
     1     (cosx-twiss(ie1,idp,mfitax)*sinx)
      y11=sqrt(twiss(ie1,idp,mfitby)/twiss(fbound%lb,idp,mfitby))
     1     *(cosy+twiss(fbound%lb,idp,mfitay)*siny)
      y22=sqrt(twiss(fbound%lb,idp,mfitby)/twiss(ie1,idp,mfitby))*
     1     (cosy-twiss(ie1,idp,mfitay)*siny)
      optstat%tracex=x11+x22
      optstat%tracey=y11+y22
      if(cell)then
        if(optstat%tracex .gt. 2.d0)then
          if(intres)then
            optstat%stabx=.false.
          else
            optstat%tracex=4.d0/optstat%tracex
          endif
        endif
        if(optstat%tracex .lt. -2.d0)then
          if(halfres)then
            optstat%stabx=.false.
          else
            optstat%tracex=4.d0/optstat%tracex
          endif
        endif
        if(optstat%tracey .gt. 2.d0)then
          if(intres)then
            optstat%staby=.false.
          else
            optstat%tracey=4.d0/optstat%tracey
          endif
        endif
        if(optstat%tracey .lt. -2.d0)then
          if(halfres)then
            optstat%staby=.false.
          else
            optstat%tracey=4.d0/optstat%tracey
          endif
        endif
        optstat%stabx=optstat%stabx .and. stab .and. (codfnd .or. fam)
        optstat%staby=optstat%staby .and. stab .and. (codfnd .or. fam)
      endif
      optstat%tracez=0.d0
      if(pri)then
        lfno=0
      endif
      return
      end

      subroutine qcell6d(fbound,idp,optstat,lfno)
      use ffs
      use ffs_fit ,only:ffs_stat
      use ffs_pointer
      use temw,only:calint,normali,tinv6,etwiss2ri,tsymp,nparams,
     $     tfinibeam,iaez,tfetwiss
      use maccbk, only:i00
      implicit none
      type (ffs_bound) ,intent(in):: fbound
      type (ffs_stat) ,intent(out):: optstat
      integer*4 ,intent(in):: lfno,idp
      integer*4 ir0
      real*8 trans(6,12),cod(6),cod0(6),beam(42),srot(3,9),
     $     btr(21,21),ri0(6,6)
      logical*4 stab,cell0,codplt0,ci0,rt,wspaccheck,calint1
      real*8 params(nparams)
      ir0=irad
      cell0=cell
      cell=cell0 .and. .not. trpt
      codplt0=codplt
      ci0=calint
      calint1=wspac .or. intra
      rt=radtaper
      codplt=.true.
      ri0=etwiss2ri(twiss(fbound%lb,idp,1:ntwissfun),normali)
      cod0=twiss(1,idp,mfitdx:mfitddp)
      cod=cod0
 1    if(cell)then
        call temit(trans,cod,beam,btr,
     $     calint1,iaez,.true.,params,stab,0)
        if(.not. stab)then
          write(lfno,*)'*****qcell6d---> Unstable optics'
          optstat%stabx=.false.
          optstat%staby=.false.
          cell=.false.
          go to 1
        endif
        normali=.true.
      else
        if(wspaccheck())then
          optstat%stabx=.false.
          optstat%staby=.false.
          return
        endif
        cod=cod0
        twiss(fbound%lb,idp,1:ntwissfun)=tfetwiss(ri0,cod,.true.)
        beam(1:21)=tfinibeam(fbound%lb)
        beam(22:42)=0.d0
        call tinitr(trans)
        call tturne0(trans,cod,beam,srot,fbound,
     $       iaez,idp,.true.,rt,.true.)
      endif
      calint=ci0
      codplt=codplt0
      cell=cell0
      irad=ir0
      return
      end
