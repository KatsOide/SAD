      module tmacro
        use mackw
        use macphys
        use macfile
        real*8, parameter :: c=cveloc,hp=plankr,e=elemch
        real*8
     $       amass,charge,h0,p0,omega0,trf0,crad,erad,epsrad,
     $       codin(6),dleng,anrad,urad,u0,vc0,wrfeff,dp0,brho,
     $       ccintr,cintrb,pbunch,coumin,re0,pgev,emidiv,
     $       emidib,emidiq,emidis,ctouck,dvemit,h1emit,
     $       anbunch,tdummy(6),zlost,alost,
     $       taurdx,taurdy,taurdz,fridiv,beamin(21),
     $       vccos,vcsin,vcphic,vcalpha,vceff,
     $       vcacc,dvcacc,ddvcacc,
     $       pspac_dx,pspac_dy,pspac_dz,dvfs,rcratio,rclassic,brhoz,
     $       bradprev
        integer*8 ilattp,lspect,ipoltr,ipolb,ipoll,ipolid,ipolo
        integer*4 nflag0,nlat,np0,nturn,isynch,nspect,
     $       lplot,nplot,nuse,nclas,irad,novfl,npelm,ipelm,
     $       nparallel,pspac_nx,pspac_ny,pspac_nz,
     $       pspac_nturn,pspac_nturncalc
        logical*4 oldflagsdummy,calint,caltouck,tparaed
      end module

      module sad_main
        integer*4, parameter ::expnsize=7

        type sad_el
        sequence
        integer*4 nlat1,dum1
        integer*8 aux,comp(1:2**31-2)
        end type

        type sad_comp
        sequence
        integer*4 ncomp2,id
        integer*8 param
        real*8 value(1:2**31-2)
        end type

        contains
        integer*8 function kmelaloc(n,el)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*4 n
        integer*8 ktaloc
        type (sad_el), pointer, intent(out) :: el
        kmelaloc=ktaloc(n+1)
        call c_f_pointer(c_loc(klist(kmelaloc-1)),el)
        return
        end function

        integer*8 function kmcompaloc(n,cmp)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*4 n
        integer*8 ktcaloc
        type (sad_comp), pointer, intent(out) :: cmp
        kmcompaloc=ktcaloc(n+expnsize+2)
        call c_f_pointer(c_loc(klist(kmcompaloc-1)),cmp)
        return
        end function

        subroutine loc_el(k,el)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*8 k
        type (sad_el), pointer, intent(out) ::el
        call c_f_pointer(c_loc(klist(k-1)),el)
        return
        end subroutine

        subroutine loc_comp(k,cmp)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*8 k
        type (sad_comp), pointer, intent(out) ::cmp
        call c_f_pointer(c_loc(klist(k-1)),cmp)
        return
        end subroutine

c$$$        subroutine latt_el(latt,el)
c$$$        use iso_c_binding
c$$$        implicit none
c$$$        integer*8 , target ::latt(1:)
c$$$        type (sad_el), pointer, intent(out) ::el
c$$$        call c_f_pointer(c_loc(latt(-1)),el)
c$$$        return
c$$$        end subroutine
c$$$
        integer*4 function idcomp(el,i)
        implicit none
        type (sad_el) :: el
        type (sad_comp), pointer :: cmp
        integer*4 i
        call loc_comp(el%comp(i),cmp)
        idcomp=cmp%id
        return
        end

        real*8 function dircomp(el,i)
        implicit none
        integer*4 i
        type (sad_el) :: el
        type (sad_comp), pointer :: cmp
        call loc_comp(el%comp(i),cmp)
        dircomp=cmp%value(cmp%ncomp2-2)
        return
        end function

        subroutine tclrpara(el,nl)
        use tmacro
c        use maccbk, only:idtype,pname
        use tfstk, only:itfcbk
        implicit none
        type (sad_el), pointer :: el
        type (sad_comp), pointer :: cmp
        integer*4 nl,i
        integer*8 lp
        do i=1,nl
          lp=el%comp(i)
          if(lp .gt. 0)then
            call loc_comp(lp,cmp)
            if(cmp%param .gt. 0)then
c              if(itfcbk(cmp%param) .eq. 0)then
c                write(*,*)i,idtype(cmp%id),pname(cmp%id)
c              else
                call tfree(cmp%param)
c              endif
              cmp%param=0
            endif
          endif
c$$$  if(idtype(ilist(2,lp)) .eq. 31)then
c$$$  iwpl=ilist(1,lp+kytbl(kwLWAK,icCAVI))
c$$$  if(iwpl .gt. 0)then
c$$$  if(ilist(2,iwpl) .gt. 0)then
c$$$  call tfree(int8(ilist(2,iwpl)))
c$$$  ilist(2,iwpl)=0
c$$$  endif
c$$$  endif
c$$$  iwpt=ilist(1,lp+kytbl(kwTWAK,icCAVI))
c$$$  if(iwpt .gt. 0)then
c$$$  if(ilist(2,iwpt) .gt. 0)then
c$$$  call tfree(int8(ilist(2,iwpt)))
c$$$  ilist(2,iwpt)=0
c$$$  endif
c$$$  endif
c$$$  endif
        enddo
        tparaed=.false.
        return
        end subroutine

      end module

      module tffitcode
      implicit none
      integer*4, parameter ::
     $     mfitax=1,mfitbx=mfitax+1,mfitnx=mfitbx+1,
     $     mfitay=mfitnx+1,mfitby=mfitay+1,mfitny=mfitby+1,
     $     mfitex=mfitny+1,mfitepx=mfitex+1,
     $     mfitey=mfitepx+1,mfitepy=mfitey+1,
     $     mfitr1=mfitepy+1,mfitr2=mfitr1+1,
     $     mfitr3=mfitr2+1,mfitr4=mfitr3+1,
     $     mfitdetr=mfitr4+1,
     $     mfitdx=mfitdetr+1,mfitdpx=mfitdx+1,
     $     mfitdy=mfitdpx+1,mfitdpy=mfitdy+1,
     $     mfitdz=mfitdpy+1,mfitddp=mfitdz+1,
     $     mfitpex=mfitddp+1,mfitpepx=mfitpex+1,
     $     mfitpey=mfitpepx+1,mfitpepy=mfitpey+1,
     $     mfittrx=mfitpepy+1,mfittry=mfittrx+1,mfitleng=mfittry+1,
     $     mfitgx=mfitleng+1,mfitgy=mfitgx+1,mfitgz=mfitgy+1,
     $     mfitchi1=mfitgz+1,mfitchi2=mfitchi1+1,mfitchi3=mfitchi2+1,
     $     ntwissfun=mfitddp,mfito=mfittry,mfit=mfitchi3,
     $     mfit1=mfit+12
      end module

      module ffs0
      use tffitcode
      use tmacro
      implicit none
      integer*4 ndim1
      parameter (ndim1=201)
      integer*4 nelmhash
      parameter (nelmhash=1023)
      real*8 , parameter :: xyth=0.375d0
      integer*4 , parameter :: maxcond=4001,lblname=maxcond/4

      type ffsv
        sequence
        integer*8 ifaux,ifibzl,ifmult,ifklp,ifival,iftwissp,
     $       iftwis,ifpos,ifgeo,ifsize,ifgamm ,ifele,ifcoup,
     $       iferrk,ifvarele,ifvvar,ifvalvar,ifele1,ifele2,
     $       ifmast,iftouchele,iftouchv,lfnp,iffserr,
     $       ifivcomp,ifvlim,iffssave,iut,ifiprev,ifinext,
     $       ielmhash
        real*8 emx,emy,dpmax,geo0(3,4),xixf,xiyf,sizedp,
     $       ctime0,ctime2,rsconv,fitval(maxcond)
        integer*4 mfitp(maxcond),ifitp(maxcond),ifitp1(maxcond),
     $       kdp(maxcond),kfitp(maxcond),kfit(maxcond),
     $       icalc(3,maxcond),iqcol(maxcond),lfp(2,maxcond),
     $       nvar,ntouch,itmax,measp,nfc,ncalc,
     $       blname(lblname),pading,mcommon
        integer*4 ndim,ndima,nele,nfit,marki,iorgx,iorgy,iorgr,
     $       mfpnt,mfpnt1,id1,id2,nve
        logical*4 updatesize
      end type

      type (ffsv), target, save:: ffv
      type (ffsv), pointer :: flv

      type flagset
        sequence
        logical*4 flags(1:0)
        logical*4
     $       rad,rfsw,radcod,calcod,
     $       intra,trpt,emiout,gauss,
     $       bipol,cell,ffsprmpt,dapert,
     $       fseed,ideal,codplt,canon,
     $       calpol,rfluct,cmplot,fourie,
     $       trsize,simulate,absweit,jitter,
     $       trgauss,lwake,twake,smearp,
     $       bunchsta,convgo,cellstab,spac,
     $       radlight,geocal,photons,wspac,
     $       selfcod,pspac,convcase,preservecase,
     $       lossmap,orbitcal,radtaper,sorg,
     $       intres,halfres,sumres,diffres
      end type

      integer*4 ,parameter :: nflag=48
      type (flagset), target, save :: fff
      character*8, save :: fname(1:nflag)= (/
     $     'RAD     ','RFSW    ','RADCOD  ','COD     ',
     1     'INTRA   ','TRPT    ','EMIOUT  ','GAUSS   ',
     1     'BIPOL   ','CELL    ','FFSPRMPT','DAPERT  ',
     1     'FIXSEED ','IDEAL   ','CODPLOT ','CANON   ',
     1     'POL     ','FLUC    ','CMPLOT  ','FOURIER ',
     1     'TRACKSIZ','SIMULATE','ABSW    ','JITTER  ',
     1     'TRGAUSS ','LWAKE   ','TWAKE   ','BARYCOD ',
     1     'BUNCHSTA','CONV    ','STABLE  ','SPAC    ',
     $     'RADLIGHT','GEOCAL  ','PHOTONS ','WSPAC   ',
     $     'SELFCOD ','PSPAC   ','CONVCASE','PRSVCASE',
     $     'LOSSMAP ','ORBITCAL','RADTAPER','SORG    ',
     $     'INTRES  ','HALFRES ','SUMRES  ','DIFFRES '/),
     $     sino(1:nflag)= (/
     $     '        ','        ','        ','        ',
     1     '        ','RING    ','        ','UNIFORM ',
     1     'UNIPOL  ','INS     ','        ','        ',
     1     'MOVESEED','REAL    ','        ','        ',
     1     '        ','DAMPONLY','        ','        ',
     1     '        ','OPERATE ','RELW    ','QUIET   ',
     1     'TRUNI   ','        ','        ','        ',
     1     'BATCHSTA','        ','UNSTABLE','        ',
     $     '        ','GEOFIX  ','        ','        ',
     $     '        ','        ','        ','        ',
     $     '        ','        ','        ','        ',
     $     '        ','        ','        ','        '/)

      integer*8, pointer :: ifvlim,ifibzl,ifmult,ifklp,ifival,iftwissp,
     $     iftwis,ifpos,ifgeo,ifsize,ifgamm ,ifele,ifcoup,
     $     iferrk,ifvarele,ifvvar,ifvalvar,ifele1,ifele2,
     $     ifmast,iftouchele,iftouchv,lfnp,iffserr,ifivcomp,iffssave,
     $     ifiprev,ifinext,ielmhash
      real*8, pointer :: emx,emy,dpmax,xixf,xiyf,sizedp
      real*8, pointer, dimension(:,:) :: geo0
      integer*4, pointer :: ndim,ndima,nele,nfit,marki,iorgx,iorgy,
     $     iorgr,mfpnt,mfpnt1,id1,id2,nve
      logical*4 , pointer :: updatesize

      contains
        subroutine tffsvinit
        ifvlim=>ffv%ifvlim
        ifibzl=>ffv%ifibzl
        ifmult=>ffv%ifmult
        ifklp=>ffv%ifklp
        ifival=>ffv%ifival
        iftwissp=>ffv%iftwissp
        iftwis=>ffv%iftwis
        ifpos=>ffv%ifpos
        ifgeo=>ffv%ifgeo
        ifsize=>ffv%ifsize
        ifgamm=>ffv%ifgamm
        ifele=>ffv%ifele
        ifcoup=>ffv%ifcoup
        iferrk=>ffv%iferrk
        ifvarele=>ffv%ifvarele
        ifvvar=>ffv%ifvvar
        ifvalvar=>ffv%ifvalvar
        ifele1=>ffv%ifele1
        ifele2=>ffv%ifele2
        ifmast=>ffv%ifmast
        iftouchele=>ffv%iftouchele
        iftouchv=>ffv%iftouchv
        lfnp=>ffv%lfnp
        iffserr=>ffv%iffserr
        ifivcomp=>ffv%ifivcomp
        iffssave=>ffv%iffssave
        ifiprev=>ffv%ifiprev
        ifinext=>ffv%ifinext
        emx=>ffv%emx
        emy=>ffv%emy
        dpmax=>ffv%dpmax
        sizedp=>ffv%sizedp
        xixf=>ffv%xixf
        xiyf=>ffv%xiyf
        geo0=>ffv%geo0(1:3,1:4)
        ndim=>ffv%ndim
        ndima=>ffv%ndima
        nele=>ffv%nele
        nfit=>ffv%nfit
        marki=>ffv%marki
        iorgx=>ffv%iorgx
        iorgy=>ffv%iorgy
        iorgr=>ffv%iorgr
        mfpnt=>ffv%mfpnt
        mfpnt1=>ffv%mfpnt1
        id1=>ffv%id1
        id2=>ffv%id2
        nve=>ffv%nve
        ielmhash=>ffv%ielmhash
        updatesize=>ffv%updatesize
        flv=>ffv
        return
        end subroutine

        real*8 function gettwiss(key,lp)
        use tfstk
        implicit none
        integer*4 key,lp
        gettwiss=rlist(iftwis+((key-1)*(2*ndim+1)+ndim)*nlat+lp-1)
        return
        end function

        subroutine r2twiss(r1,r2,r3,r4,rt1,rt2,rt3,rt4)
        implicit none
        real*8 , intent (in) :: r1,r2,r3,r4
        real*8 , intent (out) :: rt1,rt2,rt3,rt4
        real*8 detr,trr,detr1,alambda,sqrdet1,detr2,a
        detr=r1*r4-r2*r3
        trr=r1+r4
        detr1=2.d0*(1.d0-detr)
        if(trr .ge. 0.d0)then
          alambda=detr1/(trr+sqrt(trr**2+2.d0*detr1))
        else
          alambda=detr1/(trr-sqrt(trr**2+2.d0*detr1))
        endif
        sqrdet1=sqrt(1.d0+xyth-detr)
        rt1=(r1+alambda)*sqrdet1
        rt2=r2*sqrdet1
        rt3=r3*sqrdet1
        rt4=(r4+alambda)*sqrdet1
        detr2=rt1*rt4-rt2*rt3
        do while(detr2 .lt. 1.d0)
          a=1.d0/sqrt(detr2)
          rt1=rt1*a
          rt2=rt2*a
          rt3=rt3*a
          rt4=rt4*a
          detr2=rt1*rt4-rt2*rt3
        enddo
        return
        end subroutine

        subroutine dr2dtwiss(r1,r2,r3,r4,dr1,dr2,dr3,dr4,
     $     drt1,drt2,drt3,drt4)
        implicit none
        real*8 , intent (in) :: r1,r2,r3,r4,dr1,dr2,dr3,dr4
        real*8 , intent (out) :: drt1,drt2,drt3,drt4
        real*8 detr,trr,detr1,alambda,sqrdet1,d,
     $       ddetr,dtrr,ddetr1,dalambda,dsqrdet1,dd,e
        detr=r1*r4-r2*r3
        ddetr=dr1*r4+r1*dr4-dr2*r3-r2*dr3
        trr=r1+r4
        dtrr=dr1+dr4
        detr1=2.d0*(1.d0-detr)
        ddetr1=-2.d0*ddetr
        if(trr .ge. 0.d0)then
          d=sqrt(trr**2+2.d0*detr1)
        else
          d=-sqrt(trr**2+2.d0*detr1)
        endif
        dd=(dtrr*trr+ddetr1)/d
        e=1.d0/(trr+d)
        alambda=detr1*e
        dalambda=(ddetr1-detr1*e*(dtrr+dd))*e
        sqrdet1=sqrt(1.d0+xyth-detr)
        dsqrdet1=-.5d0*ddetr/sqrdet1
        drt1=(dr1+dalambda)*sqrdet1+(r1+alambda)*dsqrdet1
        drt2=dr2*sqrdet1+r2*dsqrdet1
        drt3=dr3*sqrdet1+r3*dsqrdet1
        drt4=(dr4+dalambda)*sqrdet1+(r4+alambda)*dsqrdet1
        return
        end subroutine

        subroutine twiss2r(rt1,rt2,rt3,rt4,r1,r2,r3,r4,cc)
        implicit none
        real*8 , intent (in)  :: rt1,rt2,rt3,rt4
        real*8 , intent (out) :: r1,r2,r3,r4,cc
        real*8 detr,trr1,detr1,alambda,sqrdet1
        detr=rt1*rt4-rt2*rt3
        sqrdet1=sqrt(detr)
        trr1=(rt1+rt4)/sqrdet1
        detr1=2.d0*(detr-xyth)
        if(trr1 .ge. 0.d0)then
          alambda=detr1/
     $         (trr1+sqrt(max(0.d0,trr1**2-2.d0*detr1)))
        else
          alambda=detr1/
     $         (trr1-sqrt(max(0.d0,trr1**2-2.d0*detr1)))
        endif
        r1=rt1/sqrdet1-alambda
        r2=rt2/sqrdet1
        r3=rt3/sqrdet1
        r4=rt4/sqrdet1-alambda
        cc=sqrt(detr-xyth)
        return
        end subroutine

        subroutine dtwiss2dr(rt1,rt2,rt3,rt4,
     $     drt1,drt2,drt3,drt4,dr1,dr2,dr3,dr4,dcc)
        implicit none
        real*8 , intent (in)  :: rt1,rt2,rt3,rt4,
     $       drt1,drt2,drt3,drt4
        real*8 , intent (out) :: dr1,dr2,dr3,dr4,dcc
        real*8 detr,trr1,detr1,alambda,sqrdet1,d,dd,e,cc,
     $       ddetr,dtrr1,ddetr1,dalambda,dsqrdet1
        detr=rt1*rt4-rt2*rt3
        ddetr=drt1*rt4+rt1*drt4-drt2*rt3-rt2*drt3
        sqrdet1=sqrt(detr)
        dsqrdet1=.5d0*ddetr/sqrdet1
        detr1=2.d0*(detr-xyth)
        ddetr1=2.d0*ddetr
        trr1=(rt1+rt4)/sqrdet1
        dtrr1=((drt1+drt4)-trr1*dsqrdet1)/sqrdet1
        if(trr1 .ge. 0.d0)then
          d=sqrt(max(0.d0,trr1**2-2.d0*detr1))
        else
          d=-sqrt(max(0.d0,trr1**2-2.d0*detr1))
        endif
        dd=(trr1*dtrr1-ddetr1)/d
        e=1.d0/(trr1+d)
        alambda=detr1*e
        dalambda=(ddetr1-detr1*(dtrr1+dd)*e)*e
        dr1=(drt1-rt1*dsqrdet1/sqrdet1)/sqrdet1-dalambda
        dr2=(drt2-rt2*dsqrdet1/sqrdet1)/sqrdet1
        dr3=(drt3-rt3*dsqrdet1/sqrdet1)/sqrdet1
        dr4=(drt4-rt4*dsqrdet1/sqrdet1)/sqrdet1-dalambda
        cc=sqrt(detr-xyth)
        dcc=.5d0*ddetr/cc
        return
        end subroutine

      end module

      module ffs_flag
        use ffs0, only:fff
        logical*4, pointer ::flags(:),
     $       rad,rfsw,radcod,calcod,
     $       intra,trpt,emiout,gauss,
     $       bipol,cell,ffsprmpt,dapert,
     $       fseed,ideal,codplt,canon,
     $       calpol,rfluct,cmplot,fourie,
     $       trsize,simulate,absweit,jitter,
     $       trgauss,lwake,twake,smearp,
     $       bunchsta,convgo,cellstab,spac,
     $       radlight,geocal,photons,wspac,
     $       selfcod,pspac,convcase,preservecase,
     $       lossmap,orbitcal,radtaper,sorg,
     $       intres,halfres,sumres,diffres
        
        contains
        subroutine ffs_init_flag
        flags=>fff%flags
        rad=>fff%rad
        rfsw=>fff%rfsw
        radcod=>fff%radcod
        calcod=>fff%calcod      
        trpt=>fff%trpt
        emiout=>fff%emiout
        gauss=>fff%gauss
        bipol=>fff%bipol
        intra=>fff%intra
        cell=>fff%cell
        fseed=>fff%fseed
        dapert=>fff%dapert
        ideal=>fff%ideal
        codplt=>fff%codplt
        convgo=>fff%convgo
        calpol=>fff%calpol
        rfluct=>fff%rfluct
        cmplot=>fff%cmplot
        fourie=>fff%fourie
        ffsprmpt=>fff%ffsprmpt
        trsize=>fff%trsize
        simulate=>fff%simulate
        absweit=>fff%absweit
        jitter=>fff%jitter
        trgauss=>fff%trgauss
        lwake=>fff%lwake
        twake=>fff%twake
        smearp=>fff%smearp
        bunchsta=>fff%bunchsta
        canon=>fff%canon
        cellstab=>fff%cellstab
        spac=>fff%spac
        radlight=>fff%radlight
        geocal=>fff%geocal
        photons=>fff%photons
        wspac=>fff%wspac
        selfcod=>fff%selfcod
        pspac=>fff%pspac
        convcase=>fff%convcase
        preservecase=>fff%preservecase
        lossmap=>fff%lossmap
        orbitcal=>fff%orbitcal
        radtaper=>fff%radtaper
        sorg=>fff%sorg
        intres=>fff%intres
        halfres=>fff%halfres
        sumres=>fff%sumres
        diffres=>fff%diffres
        return
        end subroutine
      end module

      module ffs
        use ffs0
        use ffs_flag
      end module

      module ffs_pointer
      use sad_main
      implicit none
      real*8 , pointer, contiguous :: errk(:,:),couple(:),
     $     valvar(:),valvar2(:,:)
      integer*8, pointer, dimension(:) :: iele2 
      integer*4, pointer, dimension(:) :: mult,iele,iele1,
     $     ival,klp,master,itouchele,itouchv,ivarele,ivcomp,ivvar,
     $     iprev,inext
      integer*4, pointer, contiguous, dimension(:,:) :: ibzl
      real*8 , pointer :: pos(:), gammab(:)
      real*8 , pointer , contiguous :: twiss(:,:,:),twiss2(:,:),
     $     geo(:,:,:),vlim(:,:),beamsize(:,:)
      integer*4 , pointer, contiguous :: itwissp(:)
      integer*8 , pointer :: latt(:)
      real*8 , pointer , contiguous :: utwiss(:,:,:)
      type (sad_el), pointer :: elatt

      contains
        subroutine ffs_init_pointer
        use tfstk
        use ffs
        use iso_c_binding
        implicit none
        call c_f_pointer(c_loc(ilist(1,ifibzl)),ibzl,[3,nlat])
        call c_f_pointer(c_loc(rlist(ifcoup)),couple,[nlat])
        call c_f_pointer(c_loc(rlist(iferrk)),errk,[2,nlat])
        call c_f_pointer(c_loc(ilist(1,ifmult)),mult,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifmast)),master,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifival)),ival,[nele])
        call c_f_pointer(c_loc(ilist(1,ifele)),iele,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifele1)),iele1,[nlat])
        call c_f_pointer(c_loc(klist(ifele2)),iele2,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifklp)),klp,[nele])
        call c_f_pointer(c_loc(rlist(ifvalvar)),valvar,[nve*2])
        call c_f_pointer(c_loc(rlist(ifvalvar)),valvar2,[nve,2])
        call c_f_pointer(c_loc(ilist(1,iftouchele)),itouchele,[nve*2])
        call c_f_pointer(c_loc(ilist(1,iftouchv)),itouchv,[nve*2])
        call c_f_pointer(c_loc(ilist(1,ifvvar)),ivvar,[nve*2])
        call c_f_pointer(c_loc(ilist(1,ifiprev)),iprev,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifinext)),inext,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifvarele)),ivarele,[nve*2])
        call c_f_pointer(c_loc(ilist(1,ifivcomp)),ivcomp,[nve*2])
        call c_f_pointer(c_loc(rlist(ifvlim)),vlim,[nele,2])
        return
        end subroutine

        subroutine ffs_init_sizep
        use tfstk
        use ffs
        use iso_c_binding
        implicit none
        integer*8 ktaloc
        if(ifsize .eq. 0)then
          ifsize=ktaloc(21*nlat)
          call c_f_pointer(c_loc(rlist(ifsize)),beamsize,[21,nlat])
          updatesize=.false.
        endif
        return
        end subroutine

        subroutine ffs_twiss_pointer
        use tfstk
        use ffs
        use iso_c_binding
        implicit none
        call c_f_pointer(c_loc(klist(ilattp+1)),latt,[nlat])
        call c_f_pointer(c_loc(rlist(iftwis)),twiss,
     $       [nlat,(2*ndim+1),ntwissfun])
        twiss(1:nlat,-ndim:ndim,1:ntwissfun)=>twiss
        twiss2(1:nlat*(2*ndim+1),1:ntwissfun)=>twiss
        call c_f_pointer(c_loc(rlist(ifpos)),pos,[nlat])
        call c_f_pointer(c_loc(rlist(ifgeo)),geo,[3,4,nlat])
        call c_f_pointer(c_loc(rlist(ifgamm)),gammab,[nlat])
        call c_f_pointer(c_loc(ilist(1,iftwissp)),itwissp,[nlat])
        return
        end subroutine

        integer*4 function idelc(i)
        implicit none
        integer*4 i
        idelc=idcomp(elatt,i)
        return
        end function

        integer*4 function idtypec(i)
        use maccbk, only:idtype
        implicit none
        integer*4 i
        idtypec=idtype(idcomp(elatt,i))
        return
        end function

        integer*8 function idvalc(i)
        use maccbk, only:idval
        implicit none
        integer*4 i
        idvalc=idval(idcomp(elatt,i))
        return
        end function

        integer*4 function lpnamec(i)
        implicit none
        integer*4 i,lpname
        lpnamec=lpname(idcomp(elatt,i))
        return
        end function

        subroutine compelc(i,cmp)
        implicit none
        type (sad_comp),pointer, intent(out) :: cmp
        integer*4 i
        call loc_comp(elatt%comp(i),cmp)
        return
        end subroutine

        character*(MAXPNAME) function pnamec(i)
        use maccbk, only:pname,MAXPNAME
        implicit none
        integer*4 i
        pnamec=pname(idcomp(elatt,i))
        return
        end function

        real*8 function direlc(i)
        implicit none
        integer*4 i
        type (sad_comp), pointer :: cmp
        call loc_comp(elatt%comp(i),cmp)
        direlc=cmp%value(cmp%ncomp2-2)
        return
        end function

        subroutine setdirelc(i,v)
        implicit none
        integer*4 i
        real*8 v
        type (sad_comp), pointer :: cmp
        call loc_comp(elatt%comp(i),cmp)
        cmp%value(cmp%ncomp2-2)=v
        return
        end subroutine

        subroutine tsetfringep(cmp,ic,dir,akk,table)
        use mackw
        implicit none
        type (sad_comp) cmp
        integer*4 ic
        real*8 dir,akk,table(4),f1in,f1out,f2in,f2out
        f1in =cmp%value(kytbl(kwF1,ic))
     $       +cmp%value(kytbl(kwF1K1F,ic))
        f1out=cmp%value(kytbl(kwF1,ic))
     $       +cmp%value(kytbl(kwF1K1B,ic))
        f2in =cmp%value(kytbl(kwF2,ic))
     $       +cmp%value(kytbl(kwF2K1F,ic))
        f2out=cmp%value(kytbl(kwF2,ic))
     $       +cmp%value(kytbl(kwF2K1B,ic))
        if(dir .ge. 0.d0)then
          table(1)=-akk*f1in*abs(f1in)/24.d0
          table(2)= akk*f2in
          table(3)=-akk*f1out*abs(f1out)/24.d0
          table(4)= akk*f2out
        else
          table(3)=-akk*f1in*abs(f1in)/24.d0
          table(4)= akk*f2in
          table(1)=-akk*f1out*abs(f1out)/24.d0
          table(2)= akk*f2out
        endif
        return
        end subroutine

        subroutine tsetfringepe(cmp,ic,dir,table)
        use mackw
        implicit none
        type (sad_comp) cmp
        integer*4 ic
        real*8 dir,table(4),f1in,f1out,f2in,f2out
        if(cmp%value(kytbl(kwL,ic)) .ne. 0.d0)then
          if(cmp%value(kytbl(kwK1,ic)) .ne. 0.d0 .or.
     $         kytbl(kwSK1,ic) .ne. 0 .and.
     $         cmp%value(kytbl(kwSK1,ic)) .ne. 0.d0)then
            f1in =cmp%value(kytbl(kwF1,ic))
     $           +cmp%value(kytbl(kwF1K1F,ic))
            f1out=cmp%value(kytbl(kwF1,ic))
     $           +cmp%value(kytbl(kwF1K1B,ic))
            f2in =cmp%value(kytbl(kwF2,ic))
     $           +cmp%value(kytbl(kwF2K1F,ic))
            f2out=cmp%value(kytbl(kwF2,ic))
     $           +cmp%value(kytbl(kwF2K1B,ic))
            if(dir .ge. 0.d0)then
              table(1)=f1in
              table(2)=f2in
              table(3)=f1out
              table(4)=f2out
            else
              table(3)=f1in
              table(4)=f2in
              table(1)=f1out
              table(2)=f2out
            endif
          else
            table=0.d0
          endif
        else
          table=0.d0
        endif
        return
        end subroutine
      end module

      module ffs_fit
      use ffs, only:ntwissfun,maxcond
      use tffitcode, only:mfit1
      integer*4 , parameter :: ndimmax=500
      integer*4 iuid(-ndimmax:ndimmax),
     $     jfam(-ndimmax:ndimmax),kfam(-ndimmax:ndimmax)
      real*8 dp(-ndimmax:ndimmax),scale(mfit1),
     $     tracex(-ndimmax:ndimmax),tracey(-ndimmax:ndimmax),
     $     dfam(4,-ndimmax:ndimmax),residual(-ndimmax:ndimmax),
     $     uini(ntwissfun,-ndimmax:ndimmax),wfit(maxcond),
     $     wiq(maxcond)
      logical*4 hstab(-ndimmax:ndimmax),vstab(-ndimmax:ndimmax)
      logical*4 fitflg,geomet,inicond,chgini
      integer*4 nut,nfam,nfam1,nfr,nqcol,nqcol1,nfcol,nfc0
      real*8 wexponent,offmw,etamax,avebeta,wsum
      character*8 , save :: nlist(1:mfit1)=(/
     $     'AX      ','BX      ','NX      ','AY      ',
     1     'BY      ','NY      ','EX      ','EPX     ',
     1     'EY      ','EPY     ','R1      ','R2      ',
     1     'R3      ','R4      ','DETR    ',
     $     'DX      ','DPX     ',
     1     'DY      ','DPY     ','DZ      ','DDP     ',
     1     'PEX     ','PEPX    ','PEY     ','PEPY    ',
     $     'TRX     ','TRY     ','LENG    ','GX      ',
     $     'GY      ','GZ      ','CHI1    ','CHI2    ',
     $     'CHI3    ','DEX     ','DEPX    ','DEY     ',
     $     'DEPY    ','DDX     ','DDPX    ','DDY     ',
     $     'DDPY    ','PDEX    ','PDEPX   ','PDEY    ',
     $     'PDEPY   '/)
      end module

      module ffs_wake
      integer*8 kwakep,kwakeelm
      integer*4 nwakep
      logical*4 wake
      integer*8 , pointer :: kwaketbl(:,:)
      integer*4 , pointer :: iwakeelm(:)
      end module

      subroutine tffs
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      
      logical*4 err
      call tffsa(1,kx,irtc)
      if(irtc .ne. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
      call tffssaveparams(-1,0,err)
      return
      end

      subroutine tffsinitparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_descriptor) kdp
      type (sad_symdef), pointer :: symddp
      integer*8 j
      real*8 rgetgl1,sigz0
      j=idvalc(1)
      emx=rlist(j+kytbl(kwEMIX,icMARK))
      if(emx .le. 0.d0)then
        emx=rgetgl1('EMITX')
      else
        call rsetgl1('EMITX',emx)
      endif
      emy=rlist(j+kytbl(kwEMIY,icMARK))
      if(emy .le. 0.d0)then
        emy=rgetgl1('EMITY')
      else
        call rsetgl1('EMITY',emy)
      endif
      dpmax=max(0.d0,rlist(j+kytbl(kwDP,icMARk)))
      kdp=kxsymbolz('DP',2,symddp)
      if(dpmax .le. 1.d-30)then
        dpmax=rfromd(symddp%value)
      endif
      if(dpmax .le. 1.d-30)then
        dpmax=0.01d0
      endif
      symddp%value=dfromr(dpmax)
      if(rlist(latt(1)+kytbl(kwBX,icMARK)) .le. 0.d0)then
        rlist(latt(1)+kytbl(kwBX,icMARK))=1.d0
      endif
      if(rlist(latt(1)+kytbl(kwBY,icMARK)) .le. 0.d0)then
        rlist(latt(1)+kytbl(kwBY,icMARK))=1.d0
      endif
      sigz0=max(0.d0,rlist(j+kytbl(kwSIGZ,icMARk)))
      call rsetgl1('SIGZ',sigz0)
      return
      end

      subroutine tffsalloc()
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_main
      implicit none
      integer*8 ktaloc,j
      integer*4 l,ntwis,k,i,itehash
      marki=1
      nlat=elatt%nlat1-1
      latt(1:nlat)=>elatt%comp(1:nlat)
      if(idtypec(1) .ne. icMARK)then
        write(*,*)'The first element must be a MARK element.'
        call forcesf()
      endif
      call tfhashelement
      ifele1=ktaloc(nlat/2+1)
      ifmult =ktaloc(nlat/2+1)
      nele=0
      LOOP_L: do l=1,nlat-1
        j=ielmhash+itehash(pnamec(l),MAXPNAME)*2
        LOOP_K: do k=1,ilist(1,j+1)
          i=ilist(1,k+klist(j+2)-1)
          if(i .ge. l)then
            exit LOOP_K
          endif
          if(idelc(i) .eq. idelc(l))then
            ilist(l,ifele1)=ilist(i,ifele1)
            cycle LOOP_L
          endif
        enddo LOOP_K
        nele=nele+1
        ilist(nele,ifmult)=l
        ilist(l,ifele1)=nele
      enddo LOOP_L
      nve=(nele+nlat)/2+10
      ifibzl =ktaloc(nlat*3/2+2)
      ifcoup=ktaloc(nlat)
      iferrk=ktaloc(nlat*2)
      ifmast =ktaloc(nlat/2+1)
      ifival=ktaloc(nele/2+1)
      ifele =ktaloc(nlat/2+1)
      ifele2=ktaloc(nlat)
      ifklp =ktaloc(nele/2+1)
      ifiprev =ktaloc(nlat/2+1)
      ifinext =ktaloc(nlat/2+1)
      ifvalvar=ktaloc(nve*2)
      iftouchele=ktaloc(nve)
      iftouchv=ktaloc(nve)
      ifvvar=ktaloc(nve)
      ifvarele=ktaloc(nve)
      ifivcomp=ktaloc(nve)
      ndim=1
      ndima=ndim*2+1
      ntwis =nlat*ndima
      iftwis=ktaloc(ntwis*ntwissfun)
      ifpos =ktaloc(nlat)
      ifgeo =ktaloc(nlat*12)
      ifgamm=ktaloc(nlat)
      iftwissp=ktaloc(nlat/2+1)
      ifvlim =ktaloc(nele*2)
      call ffs_init_pointer
      call ffs_twiss_pointer
      call tfinit
      ifsize=0
c      iwakepold=ktaloc(16)
c      ifwakep=iwakepold
c      call tclr(rlist(iwakepold),16)
c      ilist(1,iwakepold)=1
c      ilist(2,iwakepold)=italoc(3)
c      rlist(ilist(2,iwakepold))=0.d0
c      rlist(iwakepold+1)=0.d0
c      ilist(1,iwakepold+2)=1
c      ilist(2,iwakepold+2)=italoc(3)
c      ilist(1,iwakepold+5)=ndim
c      ilist(2,iwakepold+5)=nlat
c      ilist(1,iwakepold+6)=int(iftwis)
c      ilist(2,iwakepold+6)=int(ifsize)
      return
      end

      subroutine tffsfree
      use tfstk
      use ffs
      use tffitcode
      implicit none
      call tfresethash
c      call tfree(ilist(2,ifwakep))
c      call tfree(ilist(2,ifwakep+2))
c      call tfree(ifwakep)
      call tfree(iftouchv)
      call tfree(iftouchele)
      call tfree(ifivcomp)
      call tfree(ifvalvar)
      call tfree(ifvvar)
      call tfree(ifvarele)
      call tfree(iftwissp)
      if(ifsize .gt. 0)then
        call tfree(ifsize)
        ifsize=0
      endif
c      call tfree(ifgamm)
      call tfree(ifgeo)
      call tfree(ifpos)
      call tfree(iftwis)
c      call tfree(ifibzl)
      call tfree(ifmast)
      call tfree(ifmult)
      call tfree(ifcoup)
      call tfree(ifvlim)
      call tfree(ifival)
      call tfree(ifinext)
      call tfree(ifiprev)
      call tfree(ifklp)
      call tfree(ifele2)
      call tfree(ifele1)
      call tfree(ifele)
      call tfree(iferrk)
      call tfree(ifibzl)
      call tfree(ifgamm)
      return
      end

      subroutine tfffs(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      use tfcsi
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      integer*4 outfl1,irtc,narg,
     $     lfno1,lfni1,lfn11,lfret,lfrecl,
     $     isp1,itfmessage
      character*10 strfromis
      narg=isp-isp1
      if(narg .gt. 2)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      if(.not. ktfstringqd(dtastk(isp1+1),str))then
        irtc=itfmessage(9,'General::wrongtype','"String for #1"')
        return
      endif
      call tftclupdate(7)
      outfl1=outfl
      if(narg .eq. 2)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"File number for #2"')
          return
        endif
        outfl=int(rtastk(isp))
        if(outfl .eq. -1)then
          outfl=icslfno()
        endif
      else
        outfl=0
      endif
      levele=levele+1
      lfno1=icslfno()
      lfni1=icslfni()
      lfn11=icslfn1()
      lfret=icsmrk()
      lfrecl=icslrecl()
      call cssetp(lfrecl)
      call setbuf(str%str,str%nch)
      call cssetp(lfrecl)
      call tffsa(lfnp+1,kx,irtc)
      call tclrfpe
      call cssetp(lfret)
      call cssetl(lfrecl)
      call cssetlfno(lfno1)
      call cssetlfni(lfni1)
      call cssetlfn1(lfn11)
      outfl=outfl1
      if(irtc .eq. 0 .and. iffserr .ne. 0)then
        irtc=itfmessage(9,'FFS::error',strfromis(iffserr))
      endif
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfsetconvcase(f)
      use ffs_flag
      implicit none
      logical*4 f
      convcase=f
      return
      end
