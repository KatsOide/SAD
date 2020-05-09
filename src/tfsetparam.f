      subroutine tfsetparam
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,latt,idtypec,gammab
      use mathfun
      implicit none
      integer*8 ix
      real*8 rgetgl1,df
      logical*4 calgeo
      emx   =rgetgl1('EMITX')
      emy   =rgetgl1('EMITY')
      emz   =rgetgl1('EMITZ')
      sigzs =rgetgl1('SIGZ')
      sizedp=rgetgl1('SIGE')
      dpmax =rfromd(kxsymbolv('DP',2))
      if(idtypec(1) .eq. icMARK)then
        ix=latt(1)
        rlist(ix+ky_EMIX_MARK)=emx
        rlist(ix+ky_EMIY_MARK)=emy
        rlist(ix+ky_EMIZ_MARK)=emz
        rlist(ix+ky_SIGZ_MARK)=sigzs
        rlist(ix+ky_SIGE_MARK)=sizedp
        rlist(ix+ky_DP_MARK)=dpmax
        dp0=rlist(ix+ky_DDP_MARK)
      else
        dp0=0.d0
      endif
      np0   =int(rgetgl1('NP'))
      nturn =int(rgetgl1('TURNS'))
      charge=rgetgl1('CHARGE')
      amass =rgetgl1('MASS')
      pgev  =rgetgl1('MOMENTUM')
      pbunch=rgetgl1('PBUNCH')
      anbunch=rgetgl1('NBUNCH')
      coumin=rgetgl1('MINCOUP')
      emidiv=rgetgl1('EMITDIV')
      emidib=rgetgl1('EMITDIVB')
      emidiq=rgetgl1('EMITDIVQ')
      emidis=rgetgl1('EMITDIVS')
      fridiv=rgetgl1('FRINGDIV')
      alost =rgetgl1('LOSSAMPL')
      zlost =rgetgl1('LOSSDZ')
      df    =rgetgl1('FSHIFT')
      dleng =rlist(klist(ilattp)+1)*df
      pspac_nx =max(1,int(rgetgl1('PSPACNX')))
      pspac_ny =max(1,int(rgetgl1('PSPACNY')))
      pspac_nz =max(1,int(rgetgl1('PSPACNZ')))
      pspac_dx =rgetgl1('PSPACDX')
      pspac_dy =rgetgl1('PSPACDY')
      pspac_dz =rgetgl1('PSPACDZ')
      pspac_nturn =max(1,int(rgetgl1('PSPACNTURN')))
      pspac_nturncalc =max(0,int(rgetgl1('PSPACNTURNCALC')))
      call tphyzp
c      write(*,*)'tfsetparam ',gammab(1),p0
      if(gammab(1) .ne. p0 .and.  geocal)then
        gammab(1)=p0
        call tfgeo1(1,nlat,calgeo,.true.)
c        write(*,*)'gammab(1): ',gammab(1)
      endif
      call tsetgcut
      xixf  =rfromd(kxsymbolv('XIX',3))*pi2
      xiyf  =rfromd(kxsymbolv('XIY',3))*pi2
      nparallel=max(1,int(rgetgl1('NPARA')))
c      iwakepold=ifwakep
      return
      end

      subroutine tphyzp
      use tfstk
      use tmacro
      use mathfun
      use macphys
      implicit none
      brhoz =pgev/cveloc
      brho  =brhoz/abs(charge)
      p0    =pgev/amass
      h0    =p2h(p0)
      re0   =elradi
      rclassic=charge**2*re0
      crad  =sign(rclassic*(cveloc/amass)**2/p0/1.5d0,charge)
      urad  =sign(1.5d0*hp*cveloc/p0/amass/e,charge)
      erad  =55.d0/24.d0/sqrt(3.d0)*urad
      rcratio=rclassic/(hp*cveloc/amass/e)
      cuc=1.5d0*rclassic/rcratio
      anrad =5.d0/2.d0/sqrt(3.d0)*rcratio
      ccintr=(rclassic/h0**2)**2/8.d0/pi
      if(rlist(klist(ilattp)+1) .ne. 0.d0)then
        omega0=pi2*cveloc*p0/h0/rlist(klist(ilattp)+1)
      else
        omega0=0.d0
      endif
      call rsetgl1('OMEGA0',omega0)
      return
      end

      subroutine tfresetparam
      use ffs
      implicit none
      call rsetgl1('EMITX',emx)
      call rsetgl1('EMITY',emy)
      call rsetgl1('EMITZ',emz)
      call rsetgl1('SIGZ',sigzs)
      call rsetgl1('SIGE',sizedp)
      return
      end
