      module sad_main
        use tfstk, only:sad_descriptor,mbody
        integer*4, parameter ::expnsize=7
        integer*4 ,parameter :: iaidx(6,6)=reshape(
     $       (/ 1, 2, 4, 7,11,16,
     $          2, 3, 5, 8,12,17,
     $          4, 5, 6, 9,13,18,
     $          7, 8, 9,10,14,19,
     $         11,12,13,14,15,20,
     $         16,17,18,19,20,21/), (/6,6/))

        type sad_el
        sequence
        integer*4 nl,nlat0
c        type (sad_descriptor) elmv
        integer*8 aux,comp(1:mbody)
        end type

        type sad_comp
        sequence
        integer*4 ncomp2,id
        integer*4 nparam
        logical*1 ori,update,updateseg,ldummy2
        integer*4 ievar1,ievar2
        type (sad_descriptor) dvalue(1:0)
        integer*8 kvalue(1:0)
        integer*4 ivalue(2,1:0)
        logical*1 lvalue(8,1:0)
        real*8 value(1:mbody)
        end type

        contains
        integer*8 function kmelaloc(n,el)
        use tfstk, only:klist,ktaloc
        use iso_c_binding
        implicit none
        integer*4 ,intent(in):: n
        type (sad_el), pointer, intent(out) :: el
        kmelaloc=ktaloc(n+1)
        call c_f_pointer(c_loc(klist(kmelaloc-1)),el)
        el%nlat0=n
c        el%elmv%k=0
        el%aux=0
        return
        end function

        integer*8 function kmcompaloc(n,cmp)
        use tfstk, only:klist,ktraaloc
        use iso_c_binding
        implicit none
        integer*4 ,intent(in)::n
        type (sad_comp), pointer, intent(out) :: cmp
        kmcompaloc=ktraaloc(0,n+expnsize+2)+2
        call c_f_pointer(c_loc(klist(kmcompaloc-2)),cmp)
        cmp%ncomp2=n+expnsize+3
        return
        end function

        subroutine loc_el(k,el)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*8 ,intent(in)::k
        type (sad_el), pointer, intent(out) ::el
        call c_f_pointer(c_loc(klist(k-1)),el)
        return
        end subroutine

        subroutine loc_comp(k,cmp)
        use tfstk, only:klist
        use iso_c_binding
        implicit none
        integer*8 , intent(in)::k
        type (sad_comp), pointer, intent(out) ::cmp
        call c_f_pointer(c_loc(klist(k-2)),cmp)
        return
        end subroutine

        integer*4 function idcomp(el,i)
        implicit none
        type (sad_el),intent(in) :: el
        type (sad_comp), pointer :: cmp
        integer*4 ,intent(in)::i
        call loc_comp(el%comp(i),cmp)
        idcomp=cmp%id
        return
        end

        real*8 function dircomp(el,i)
        implicit none
        integer*4 ,intent(in)::i
        type (sad_el),intent(in) :: el
        type (sad_comp), pointer :: cmp
        call loc_comp(el%comp(i),cmp)
        dircomp=merge(1.d0,-1.d0,cmp%ori)
        return
        end function

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
     $     mfitaz=mfitddp+1,mfitbz=mfitaz+1,mfitnz=mfitbz+1,
     $     mfitzx=mfitnz+1,mfitzpx=mfitzx+1,
     $     mfitzy=mfitzpx+1,mfitzpy=mfitzy+1,
     $     mfitpex=mfitzpy+1,
     $     mfitpepx=mfitpex+1,
     $     mfitpey=mfitpepx+1,mfitpepy=mfitpey+1,
     $     mfitpzx=mfitpepy+1,
     $     mfitpzpx=mfitpzx+1,
     $     mfitpzy=mfitpzpx+1,mfitpzpy=mfitpzy+1,
     $     mfitgmx=mfitpzpy+1,mfitgmy=mfitgmx+1,mfitgmz=mfitgmy+1,
     $     mfitbmagx=mfitgmz+1,mfitbmagy=mfitbmagx+1,mfitbmagz=mfitbmagy+1,
     $     mfittrx=mfitbmagz+1,mfittry=mfittrx+1,mfittrz=mfittry+1,
     $     mfitleng=mfittrz+1,
     $     mfitgx=mfitleng+1,mfitgy=mfitgx+1,mfitgz=mfitgy+1,
     $     mfitchi1=mfitgz+1,mfitchi2=mfitchi1+1,mfitchi3=mfitchi2+1,
     $     ntwissfun=mfitzpy,mfito=mfitbmagz,mfit=mfitchi3,
     $     mfit1=mfit+12
      logical*4 :: inittws=.true.

      contains
      subroutine tfinittws
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      character*2 strfromis
      if( .not. inittws)then
        return
      endif
      call tfevals('BeginPackage[tws$`];Begin[tws$`]',kx,irtc)
      call tfevals(
     $     ' mfax='//strfromis(mfitax)//
     $     ';mfbx='//strfromis(mfitbx)//
     $     ';mfnx='//strfromis(mfitnx)//
     $     ';mfay='//strfromis(mfitay)//
     $     ';mfby='//strfromis(mfitby)//
     $     ';mfny='//strfromis(mfitny)//
     $     ';mfex='//strfromis(mfitex)//
     $     ';mfepx='//strfromis(mfitepx)//
     $     ';mfey='//strfromis(mfitey)//
     $     ';mfepy='//strfromis(mfitepy)//
     $     ';mfr1='//strfromis(mfitr1)//
     $     ';mfr2='//strfromis(mfitr2)//
     $     ';mfr3='//strfromis(mfitr3)//
     $     ';mfr4='//strfromis(mfitr4)//
     $     ';mfdetr='//strfromis(mfitdetr)//
     $     ';mfdx='//strfromis(mfitdx)//
     $     ';mfdpx='//strfromis(mfitdpx)//
     $     ';mfdy='//strfromis(mfitdy)//
     $     ';mfdpy='//strfromis(mfitdpy)//
     $     ';mfdz='//strfromis(mfitdz)//
     $     ';mfddp='//strfromis(mfitddp)//
     $     ';mfaz='//strfromis(mfitaz)//
     $     ';mfbz='//strfromis(mfitbz)//
     $     ';mfnz='//strfromis(mfitnz)//
     $     ';mfzx='//strfromis(mfitzx)//
     $     ';mfzpx='//strfromis(mfitzpx)//
     $     ';mfzy='//strfromis(mfitzy)//
     $     ';mfzpy='//strfromis(mfitzpy)//
     $     ';mfpex='//strfromis(mfitpex)//
     $     ';mfpepx='//strfromis(mfitpepx)//
     $     ';mfpey='//strfromis(mfitpey)//
     $     ';mfpepy='//strfromis(mfitpepy)//
     $     ';mfpzx='//strfromis(mfitpzx)//
     $     ';mfpzpx='//strfromis(mfitpzpx)//
     $     ';mfpzy='//strfromis(mfitpzy)//
     $     ';mfpzpy='//strfromis(mfitpzpy)//
     $     ';mfgmx='//strfromis(mfitgmx)//
     $     ';mfgmy='//strfromis(mfitgmy)//
     $     ';mfgmz='//strfromis(mfitgmz)//
     $     ';mftrx='//strfromis(mfittrx)//
     $     ';mftry='//strfromis(mfittry)//
     $     ';mfleng='//strfromis(mfitleng)//
     $     ';mfgx='//strfromis(mfitgx)//
     $     ';mfgy='//strfromis(mfitgy)//
     $     ';mfgz='//strfromis(mfitgz)//
     $     ';mfchi1='//strfromis(mfitchi1)//
     $     ';mfchi2='//strfromis(mfitchi2)//
     $     ';mfchi3='//strfromis(mfitchi3)//
     $     ';ntwissfun='//strfromis(ntwissfun)//
     $     ";End[];EndPackage[]",kx,irtc)
      inittws=.false.
      return
      end subroutine

      end module

      module ffs0
      use tfstk
      use tffitcode
      use tmacro
      implicit none
      integer*4 ndim1
      parameter (ndim1=201)
      integer*4 nelmhash
      parameter (nelmhash=1023)
      real*8 , parameter :: xyth=0.375d0
      integer*4 , parameter :: maxcond=4001,lblname=maxcond/4

      integer*4 , parameter :: lnvev=5,lnelv=4
      type nvev
        sequence
        real*8 valvar,valvar2
        integer*4 itouchele,itouchv,ivvar,ivarele,ivcomp,idummy
      end type

      type nelv
        sequence
        type (sad_descriptor) dcomp
        real*8 vlim(2)
        integer*4 ival,klp
      end type

      type ffsv
        sequence
        integer*8 ifaux,ifibzl,ifmult,iftwissp,
     $       iftwis,ifpos,ifgeo,ifsize,ifgamm,ifcomp,ifcoup,
     $       iferrk,ifele1,ifele2,ifmast,iffserr,
     $       iffssave,iut,ifiprev,ifinext,
     $       ielmhash,ifnvev,ifnelv
        real*8 emx,emy,emz,sigzs,fshifts,
     $       dpmax,geo0(3,4),xixf,xiyf,sizedp,
     $       ctime0,ctime2,rsconv,fitval(maxcond)
        integer*4 mfitp(maxcond),ifitp(maxcond),ifitp1(maxcond),
     $       kdp(maxcond),kfitp(maxcond),kfit(maxcond),
     $       icalc(3,maxcond),iqcol(maxcond),lfp(2,maxcond),
     $       nvar,ntouch,itmax,measp,nfc,ncalc,
     $       blname(lblname),pading,mcommon
        integer*4 ndim,ndima,nele,nfit,marki,iorgx,iorgy,iorgr,
     $       mfpnt,mfpnt1,id1,id2,nve,modesize
        logical*4 setref,evarini
      end type

      type (ffsv), target, save:: ffv
      type (ffsv), pointer :: flv

c$$$ft={
c$$${"TRPT    ","RING    ","trpt"},
c$$${"CELL    ","INS     ","cell"},
c$$${"RFSW    ","        ","rfsw"},
c$$${"RAD     ","        ","rad"},
c$$${"FLUC    ","DAMPONLY","rfluct"},
c$$${"RADCOD  ","        ","radcod"},
c$$${"RADTAPER","        ","radtaper"},
c$$${"KEEPEXP ","        ","keepexp"},
c$$${"CALEXP  ","        ","calexp"},
c$$${"CALC6D  ","CALC4D  ","calc6d"},
c$$${"EMIOUT  ","        ","emiout"},
c$$${"CODPLOT ","        ","codplt"},
c$$${"COD     ","        ","calcod"},
c$$${"INTRA   ","        ","intra"},
c$$${"POL     ","        ","calpol"},
c$$${"RADPOL  ","        ","radpol"},
c$$${"LWAKE   ","        ","lwake"},
c$$${"TWAKE   ","        ","twake"},
c$$${"WSPAC   ","        ","wspac"},
c$$${"SPAC    ","        ","spac"},
c$$${"SELFCOD ","        ","selfcod"},
c$$${"CONV    ","        ","convgo"},
c$$${"STABLE  ","UNSTABLE","cellstab"},
c$$${"GEOCAL  ","GEOFIX  ","geocal"},
c$$${"RADLIGHT","        ","radlight"},
c$$${"PHOTONS ","        ","photons"},
c$$${"LOSSMAP ","        ","lossmap"},
c$$${"SORG    ","        ","sorg"},
c$$${"INTRES  ","        ","intres"},
c$$${"HALFRES ","        ","halfres"},
c$$${"SUMRES  ","        ","sumres"},
c$$${"DIFFRES ","        ","diffres"},
c$$${"FFSPRMPT","        ","ffsprmpt"},
c$$${"CONVCASE","        ","convcase"},
c$$${"PRSVCASE","        ","preservecase"},
c$$${"SUS     ","        ","suspend"},
c$$${"K64     ","LEGACY  ","k64"},
c$$${"GAUSS   ","UNIFORM ","gauss"},
c$$${"FIXSEED ","MOVESEED","fseed"},
c$$${"BIPOL   ","UNIPOL  ","bipol"},
c$$${"PSPAC   ","        ","pspac"},
c$$${"ORBITCAL","        ","orbitcal"},
c$$${"CALOPT  ","ORBONLY ","calopt"},
c$$${"DAPERT  ","        ","dapert"},
c$$${"IDEAL   ","REAL    ","ideal"},
c$$${"FOURIER ","        ","fourie"},
c$$${"TRACKSIZ","        ","trsize"},
c$$${"SIMULATE","OPERATE ","simulate"},
c$$${"ABSW    ","RELW    ","absweit"},
c$$${"JITTER  ","QUIET   ","jitter"},
c$$${"TRGAUSS ","TRUNI   ","trgauss"},
c$$${"BARYCOD ","        ","smearp"},
c$$${"WAKEOPT ","        ","wakeopt"},
c$$${"        ","        ","dummyf2"}
c$$$};
c$$$
c$$$ftt=Thread[ft];
c$$$(Print["     $  ",Null@@Table["'"//#[[k]]//"',",{k,4}]]&/@Partition[#,4])&/@ftt[[{1,2}]];
c$$$Print["     $  ",Null@@Table[#[[k]]//",",{k,4}]]&/@Partition[ftt[[3]],4];
c$$$susp;

      type flagset
        sequence
        logical*4 flags(1:0)
        logical*4
     $       trpt,cell,rfsw,rad,
     $       rfluct,radcod,radtaper,calc6d,
     $       keepexp,calexp,emiout,codplt,
     $       calcod,intra,calpol,radpol,
     $       lwake,twake,wspac,spac,
     $       selfcod,convgo,cellstab,geocal,
     $       radlight,photons,lossmap,sorg,
     $       intres,halfres,sumres,diffres,
     $       ffsprmpt,convcase,preservecase,suspend,
     $       k64,gauss,fseed,bipol,
     $       pspac,orbitcal,calopt,dapert,
     $       ideal,fourie,trsize,simulate,
     $       absweit,jitter,trgauss,smearp,
     $       wakeopt
      end type

      integer*4 ,parameter :: nflag=53
      type (flagset), target, save :: fff
      character*8, save :: fname(1:nflag)=(/
     $  'TRPT    ','CELL    ','RFSW    ','RAD     ',
     $  'FLUC    ','RADCOD  ','RADTAPER','CALC6D  ',
     $  'KEEPEXP ','CALEXP  ','EMIOUT  ','CODPLOT ',
     $  'COD     ','INTRA   ','POL     ','RADPOL  ',
     $  'LWAKE   ','TWAKE   ','WSPAC   ','SPAC    ',
     $  'SELFCOD ','CONV    ','STABLE  ','GEOCAL  ',
     $  'RADLIGHT','PHOTONS ','LOSSMAP ','SORG    ',
     $  'INTRES  ','HALFRES ','SUMRES  ','DIFFRES ',
     $  'FFSPRMPT','CONVCASE','PRSVCASE','SUS     ',
     $  'K64     ','GAUSS   ','FIXSEED ','BIPOL   ',
     $  'PSPAC   ','ORBITCAL','CALOPT  ','DAPERT  ',
     $  'IDEAL   ','FOURIER ','TRACKSIZ','SIMULATE',
     $  'ABSW    ','JITTER  ','TRGAUSS ','BARYCOD ',
     $  'WAKEOPT '/),
     $     sino(1:nflag)=(/
     $  'RING    ','INS     ','        ','        ',
     $  'DAMPONLY','        ','        ','CALC4D  ',
     $  '        ','        ','        ','        ',
     $  '        ','        ','        ','        ',
     $  '        ','        ','        ','        ',
     $  '        ','        ','UNSTABLE','GEOFIX  ',
     $  '        ','        ','        ','        ',
     $  '        ','        ','        ','        ',
     $  '        ','        ','        ','        ',
     $  'LEGACY  ','UNIFORM ','MOVESEED','UNIPOL  ',
     $  '        ','        ','ORBONLY ','        ',
     $  'REAL    ','        ','        ','OPERATE ',
     $  'RELW    ','QUIET   ','TRUNI   ','        ',
     $  '        '/)

      integer*8, pointer :: ifibzl,ifmult,iftwissp,
     $     iftwis,ifpos,ifgeo,ifsize,ifgamm,ifcomp,ifcoup,
     $     iferrk,ifele1,ifele2,ifmast,iffserr,iffssave,
     $     ifiprev,ifinext,ielmhash,ifnvev,ifnelv
      real*8, pointer :: emx,emy,emz,dpmax,xixf,xiyf,
     $     sizedp,sigzs,fshifts
      real*8, pointer, dimension(:,:) :: geo0
      integer*4, pointer :: ndim,ndima,nele,nfit,marki,iorgx,iorgy,
     $     iorgr,mfpnt,mfpnt1,id1,id2,nve,ntouch,modesize
      logical*4 , pointer :: setref,evarini
      type (nvev), pointer,dimension(:)::nvevx
      type (nelv), pointer,dimension(:)::nelvx

      type ffs_bound
      sequence
      integer*4 lb,le
      real*8 fb,fe
      end type

      contains
        subroutine tffsvinit
        ifibzl=>ffv%ifibzl
        ifmult=>ffv%ifmult
        iftwissp=>ffv%iftwissp
        iftwis=>ffv%iftwis
        ifpos=>ffv%ifpos
        ifgeo=>ffv%ifgeo
        ifsize=>ffv%ifsize
        ifgamm=>ffv%ifgamm
        ifcomp=>ffv%ifcomp
        ifcoup=>ffv%ifcoup
        iferrk=>ffv%iferrk
        ifele1=>ffv%ifele1
        ifele2=>ffv%ifele2
        ifmast=>ffv%ifmast
        iffserr=>ffv%iffserr
        iffssave=>ffv%iffssave
        ifiprev=>ffv%ifiprev
        ifinext=>ffv%ifinext
        emx=>ffv%emx
        emy=>ffv%emy
        emz=>ffv%emz
        dpmax=>ffv%dpmax
        sizedp=>ffv%sizedp
        sigzs=>ffv%sigzs
        fshifts=>ffv%fshifts
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
        ifnvev=>ffv%ifnvev
        ifnelv=>ffv%ifnelv
        ielmhash=>ffv%ielmhash
        modesize=>ffv%modesize
        ntouch=>ffv%ntouch
        setref=>ffv%setref
        evarini=>ffv%evarini
        flv=>ffv
        return
        end subroutine

        real*8 pure elemental function gettwiss(key,lp)
        implicit none
        integer*4 ,intent(in):: key,lp
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
        alambda=merge(detr1/(trr+sqrt(trr**2+2.d0*detr1)),
     $       detr1/(trr-sqrt(trr**2+2.d0*detr1)),
     $       trr >= 0.d0)
        sqrdet1=sqrt(1.d0+xyth-detr)
        rt1=(r1+alambda)*sqrdet1
        rt2=r2*sqrdet1
        rt3=r3*sqrdet1
        rt4=(r4+alambda)*sqrdet1
        detr2=rt1*rt4-rt2*rt3
        do while(detr2 < 1.d0)
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
        d=sign(sqrt(trr**2+2.d0*detr1),trr)
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
        alambda=detr1/
     $       (trr1+sign(sqrt(max(0.d0,trr1**2-2.d0*detr1)),trr1))
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
        d=sign(sqrt(max(0.d0,trr1**2-2.d0*detr1)),trr1)
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

        integer*8 function ktatwissaloc(mode,kl)
        use tffitcode
        implicit none
        type (sad_rlist), pointer,intent(out)::kl
        integer*4 ,intent(in):: mode
        integer*8 kax
        kax=ktraaloc(mode,28,kl)
        kl%rbody(mfitbz)=1.d0
        ktatwissaloc=kax
        return
        end function

      end module

      module ffs_flag
        use ffs0, only:fff
        logical*4, pointer ::flags(:),
     $       trpt,cell,rfsw,rad,
     $       rfluct,radcod,radtaper,calc6d,
     $       keepexp,calexp,emiout,codplt,
     $       calcod,intra,calpol,radpol,
     $       lwake,twake,wspac,spac,
     $       selfcod,convgo,cellstab,geocal,
     $       radlight,photons,lossmap,sorg,
     $       intres,halfres,sumres,diffres,
     $       ffsprmpt,convcase,preservecase,suspend,
     $       k64,gauss,fseed,bipol,
     $       pspac,orbitcal,calopt,dapert,
     $       ideal,fourie,trsize,simulate,
     $       absweit,jitter,trgauss,smearp,wakeopt
        
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
        k64=>fff%k64
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
        wakeopt=>fff%wakeopt
        radpol=>fff%radpol
        calc6d=>fff%calc6d
        keepexp=>fff%keepexp
        calexp=>fff%calexp
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
        calopt=>fff%calopt
        suspend=>fff%suspend
        return
        end subroutine

        integer*4 function ndivrad(phi,ak,bz,eps)
        use tmacro, only:h0
        implicit none
        real*8 , parameter :: arad=0.02d0
        real*8 , intent(in)::ak,bz,phi,eps
        real*8 aka,eps1
        eps1=merge(1.d0,eps,eps <= 0.d0)
        aka=(abs(phi)+hypot(ak,bz)*arad)/eps1
        ndivrad=merge(int(1.d0+h0/100.d0*aka*10.d0),
     $       int(1.d0+h0/100.d0*aka),trpt)
        return
        end function 

      end module

      module ffs
        use ffs0
        use ffs_flag
        use tfstk, only:resetnan

        contains
        subroutine limitcod(cod)
        real*8 , intent(inout)::cod(6)
        real*8 , parameter :: omax=1.d5,dpmin=-1.d0+1.d-5,
     $       pmax=1.d5,zmax=1.d100,codnan=1.d100
        real*8 ptmax
        call resetnan(cod,codnan)
        cod(6)=max(dpmin,min(pmax,cod(6)))
        ptmax=1.d0+cod(6)
        cod(1)=max(-omax,min(omax,cod(1)))
        cod(3)=max(-omax,min(omax,cod(3)))
        cod(2)=max(-ptmax,min(ptmax,cod(2)))
        cod(4)=max(-ptmax,min(ptmax,cod(4)))
        cod(5)=max(-zmax,min(zmax,cod(5)))
        return
        end subroutine
      end module

      module ffs_pointer
      use sad_main
      use ffs
      implicit none
      real*8 , pointer, contiguous :: errk(:,:),couple(:)
      integer*8, pointer, dimension(:) :: kele2 
      integer*4, pointer, dimension(:) :: mult,icomp,iele1,
     $     master,iprev,inext
      integer*4, pointer, contiguous, dimension(:,:) :: ibzl
      real*8 , pointer :: pos(:), gammab(:)
      real*8 , pointer , contiguous :: twiss(:,:,:),twiss2(:,:),
     $     geo(:,:,:),beamsize(:,:)
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
c        call c_f_pointer(c_loc(ilist(1,ifival)),ival,[nele])
c        call c_f_pointer(c_loc(dlist(ifdcomp)),dcomp,[nele])
        call c_f_pointer(c_loc(ilist(1,ifcomp)),icomp,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifele1)),iele1,[nlat])
        call c_f_pointer(c_loc(klist(ifele2)),kele2,[nlat])
c        call c_f_pointer(c_loc(ilist(1,ifklp)),klp,[nele])
        call c_f_pointer(c_loc(ilist(1,ifiprev)),iprev,[nlat])
        call c_f_pointer(c_loc(ilist(1,ifinext)),inext,[nlat])
        call c_f_pointer(c_loc(klist(ifnvev)),nvevx,[nve])
        call c_f_pointer(c_loc(klist(ifnelv)),nelvx,[nele])
        return
        end subroutine

        subroutine ffs_init_sizep
        use tfstk
        use ffs
        use iso_c_binding
        implicit none
        if(ifsize == 0)then
          ifsize=ktaloc(21*nlat)
          call c_f_pointer(c_loc(rlist(ifsize)),beamsize,[21,nlat])
          modesize=0
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
        integer*4 ,intent(in):: i
        idelc=idcomp(elatt,i)
        return
        end function

        integer*4 function idtypec(i)
        use maccbk, only:idtype
        implicit none
        integer*4 ,intent(in)::i
        idtypec=idtype(idcomp(elatt,i))
        return
        end function

        integer*4 function idtypecx(i)
        use maccbk, only:idtype
        use maccode
        use tmacro, only:nlat
        implicit none
        integer*4 ,intent(in):: i
        idtypecx=merge(icNull,idtype(idcomp(elatt,i)),
     $       i <= 0 .or. i .ge. nlat)
        return
        end function

        integer*8 function idvalc(i)
        use maccbk, only:idval
        implicit none
        integer*4 ,intent(in):: i
        idvalc=idval(idcomp(elatt,i))
        return
        end function

        integer*4 function lpnamec(i)
        use maccbk
        implicit none
        integer*4 ,intent(in):: i
        lpnamec=lpname(idcomp(elatt,i))
        return
        end function

        subroutine compelc(i,cmp)
        use sad_main
        implicit none
        type (sad_comp),pointer, intent(out) :: cmp
        integer*4,intent(in)::i
        call loc_comp(elatt%comp(i),cmp)
        return
        end subroutine

        character*(MAXPNAME) function pnamec(i)
        use maccbk, only:pname,MAXPNAME
        implicit none
        integer*4,intent(in)::i
        pnamec=pname(idcomp(elatt,i))
        return
        end function

        real*8 function direlc(i)
        implicit none
        integer*4,intent(in)::i
        type (sad_comp), pointer :: cmp
        call loc_comp(elatt%comp(i),cmp)
        direlc=merge(1.d0,-1.d0,cmp%ori)
        return
        end function

        subroutine setdirelc(i,v)
        implicit none
        integer*4 ,intent(in)::i
        real*8 ,intent(in)::v
        type (sad_comp), pointer :: cmp
        call loc_comp(elatt%comp(i),cmp)
        cmp%ori=v > 0.d0
        return
        end subroutine

        subroutine tffsbound(fbound)
        use tfstk
c        use ffs_pointer
        use mackw
        implicit none
        type (ffs_bound) ,intent(out):: fbound
        call tffsbound1(1,nlat,fbound)
        return
        end subroutine

        subroutine tffsbound1(lb,le,fbound)
        use tfstk
c        use ffs_pointer
        use mackw
        implicit none
        type (ffs_bound) ,intent(out):: fbound
        integer*4 ,intent(in):: lb,le
        integer*4 le1
        real*8 xnlat,offset
        xnlat=nlat
        offset=tffsmarkoffset(lb)
        if(offset /= 0.d0)then
          offset=max(1.d0,min(xnlat,lb+offset))
          fbound%lb=int(offset)
          fbound%fb=offset-fbound%lb
        else
          fbound%lb=lb
          fbound%fb=0.d0
        endif
        le1=le
        if(le == nlat)then
          if(idtypec(nlat-1) == icMARK)then
            le1=le1-1
          else
            fbound%le=le1
            fbound%fe=0.d0
            return
          endif
        endif
        offset=tffsmarkoffset(le1)
        if(offset /= 0.d0)then
          offset=max(1.d0,min(xnlat,le1+offset))
          fbound%le=int(offset)
          fbound%fe=offset-fbound%le
        else
          fbound%le=le1
          fbound%fe=0.d0
        endif
        return
        end subroutine

        real*8 function tffselmoffset(l)
        use tfstk
        use ffs
        use mackw
        use tfcsi, only:icslfnm
        implicit none
        integer*4 ,intent(in):: l
        integer*4 lm,nm,lx,nmmax
        parameter (nmmax=256)
        real*8 offset,xp,xe
        nm=0
        xp=l
        if(l /= nlat .and. kytbl(kwOFFSET,idtypec(l)) /= 0)then
          xe=nlat
          lm=l
          do
            offset=tffsmarkoffset(lm)
            if(offset /= 0.d0)then
              xp=offset+lm
              if(xp >= 1.d0 .and. xp <= xe)then
                lx=int(xp)
                if(idtypec(lx) == icMARK)then
                  nm=nm+1
                  if(nm .lt. nmmax)then
                    lm=lx
                    cycle
                  else
                    call termes('?Recursive OFFSET in',pnamec(l))
                  endif
                endif
              endif
            endif
            exit
          enddo
        endif
        tffselmoffset=xp
        return
        end function

        integer*4 function nextl(i)
        use tmacro, only:nlat
        use mackw
        implicit none
        integer*4 ,intent(in):: i
        integer*4 it,i1,k
        type (sad_comp), pointer :: cmp
        i1=i+1
 1      if(i1 .ge. nlat)then
          nextl=nlat
        else
          it=idtypec(i1)
          k=kytbl(kwOFFSET,it)
          if(k /= 0)then
            call loc_comp(elatt%comp(i1),cmp)
            if(cmp%value(k) /= 0.d0)then
              i1=i1+1
              go to 1
            endif
          endif
          nextl=i1
        endif
        return
        end function

        real*8 function tfvalvar(i,k) result(v)
        implicit none
        type (sad_comp), pointer :: cmps
        integer*4 ,intent(in):: i,k
        call compelc(i,cmps)
        v=cmps%value(k)
        return
        end function

        subroutine tsetfringep(cmp,ic,akk,table)
        use mackw
        implicit none
        type (sad_comp) ,intent(in):: cmp
        integer*4 ,intent(in):: ic
        real*8 ,intent(in):: akk
        real*8 ,intent(out):: table(4)
        real*8 f1in,f1out,f2in,f2out
        f1in =cmp%value(kytbl(kwF1,ic))
     $       +cmp%value(kytbl(kwF1K1F,ic))
        f1out=cmp%value(kytbl(kwF1,ic))
     $       +cmp%value(kytbl(kwF1K1B,ic))
        f2in =cmp%value(kytbl(kwF2,ic))
     $       +cmp%value(kytbl(kwF2K1F,ic))
        f2out=cmp%value(kytbl(kwF2,ic))
     $       +cmp%value(kytbl(kwF2K1B,ic))
        if(cmp%ori)then
          table(1)=-abs(akk*f1in*f1in)/24.d0
          table(2)= abs(akk)*f2in
          table(3)=-abs(akk*f1out*f1out)/24.d0
          table(4)= abs(akk)*f2out
        else
          table(3)=-abs(akk*f1in*f1in)/24.d0
          table(4)= abs(akk)*f2in
          table(1)=-abs(akk*f1out*f1out)/24.d0
          table(2)= abs(akk)*f2out
        endif
        return
        end subroutine

        subroutine tsetfringepe(cmp,ic,table)
        use mackw
        implicit none
        type (sad_comp) ,intent(in):: cmp
        integer*4 ,intent(in):: ic
        real*8 ,intent(out):: table(4)
        real*8 f1in,f1out,f2in,f2out
        if(cmp%value(kytbl(kwL,ic)) /= 0.d0)then
          if(cmp%value(kytbl(kwK1,ic)) /= 0.d0 .or.
     $         kytbl(kwSK1,ic) /= 0 .and.
     $         cmp%value(kytbl(kwSK1,ic)) /= 0.d0)then
            f1in =cmp%value(kytbl(kwF1,ic))
     $           +cmp%value(kytbl(kwF1K1F,ic))
            f1out=cmp%value(kytbl(kwF1,ic))
     $           +cmp%value(kytbl(kwF1K1B,ic))
            f2in =cmp%value(kytbl(kwF2,ic))
     $           +cmp%value(kytbl(kwF2K1F,ic))
            f2out=cmp%value(kytbl(kwF2,ic))
     $           +cmp%value(kytbl(kwF2K1B,ic))
            if(cmp%ori)then
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

        integer*4 pure function ielma(i)
        use ffs_flag,only:trpt
        use tmacro, only:nlat
        implicit none
        integer*4 ,intent(in):: i
        ielma=merge(merge(min(nlat,i),max(1,nlat+i+1),i > 0),
     $       merge(mod(i-1,nlat)+1,max(nlat-mod(-i,nlat)+1,1),
     $       i > 0),trpt)
        return
        end function

       real*8 function tffsmarkoffset(lp)
       use tfstk
       use sad_main
       use mackw
       use tmacro, only:nlat
       use ffs
       implicit none
       type (sad_comp), pointer :: cmp
       integer*4 ,intent(in):: lp
       integer*4 k
       if(lp >= nlat)then
         tffsmarkoffset=0.d0
         return
       endif
       call compelc(lp,cmp)
       k=kytbl(kwOFFSET,idtype(cmp%id))
       if(k > 0)then
         tffsmarkoffset=cmp%value(k)
         if(tffsmarkoffset /= 0.d0)then
           tffsmarkoffset=(tffsmarkoffset-.5d0)*direlc(lp)+.5d0
         endif
       else
         tffsmarkoffset=0
       endif
       return
       end function

      end module

      module multa
        use kyparam, only:nmult
        real*8 ,parameter :: fact(0:nmult+1)=[
     $       1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1       720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $       39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $       1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $       6402373705728000.d0,121645100408832000.d0,
     $       2432902008176640000.d0,51090942171709440000.d0,
     $       1124000727777607680000.d0],
     $       aninv(0:nmult+1)=[1.d0,1.d0,
     $       0.5d0,
     $       0.33333333333333333333d0,
     $       0.25d0,
     $       0.2d0,
     $       0.166666666666666666667d0,
     $       0.142857142857142857143d0,
     $       0.125d0,
     $       0.111111111111111111111d0,
     $       0.1d0,
     $       0.090909090909090909091d0,
     $       0.083333333333333333333d0,
     $       0.076923076923076923077d0,
     $       0.071428571428571428571d0,
     $       0.066666666666666666667d0,
     $       0.0625d0,
     $       0.058823529411764705882d0,
     $       0.055555555555555555556d0,
     $       0.052631578947368421053d0,
     $       0.05d0,
     $       0.047619047619047619048d0,
     $       0.045454545454545454545d0]
        logical*4 :: gknini=.true.
        real*8 :: gkn(0:nmult,0:nmult)=0.d0

        contains
        subroutine gkninit
        implicit none
        integer*4 n,k
        do n=0,nmult
          gkn(n,0)=1.d0
          do k=1,nmult-n
            gkn(n,k)=gkn(n,k-1)
     $           *dble((2*k-3)*(2*k+1))/8.d0/dble(k*(k+n+1))
          enddo
        enddo
        gknini=.false.
        return
        end subroutine
      end module

      module tparastat
        use sad_main
        use tmacro
        real*8, save :: rstat(2,icMARK)=0.d0
        logical*4, save :: lstat(2,icMARK)=.true.

        contains

        logical*4 function tparacheck(l,cmp)
        use ffs_flag
        implicit none
        type (sad_comp),intent(in):: cmp
        integer*4 ,intent(in):: l

        select case (l)
        case (icMULT,icCAVI)
          tparacheck=rstat(1,l) /= amass .or.
     $         rstat(2,l) /= charge .or.
     $         lstat(1,l) .neqv. trpt
          if(tparacheck)then
            rstat(1,l)=amass
            rstat(2,l)=charge
            lstat(1,l)=trpt
          else
            tparacheck=.not. cmp%update
          endif

        case default
          tparacheck=.not. cmp%update

        end select
        return
        end function

        subroutine tpara(cmp)
        use kyparam
        use tfstk
        use ffs_pointer, only:idelc,idvalc,compelc,tsetfringep
        use ffs_flag, only:trpt
        use multa, only:fact,aninv
        use mathfun, only:akang
        implicit none
        type (sad_comp) , intent(inout):: cmp
        integer*4 ltyp,n,ndiv,nmmax
        real*8 phi,al,psi1,psi2,theta,dtheta,w,akk,sk1,
     $       fb1,fb2,harm,vnominal,frmd,eps,
     $       cchi1,cchi2,cchi3,schi1,schi2,schi3
        complex*16 cr1,cr,ck
        cmp%update=.true.
        ltyp=idtype(cmp%id)
        if(kytbl(kwNPARAM,ltyp) == 0)then
          return
        endif
        select case (ltyp)
        case (icBEND)
          al=cmp%value(ky_L_BEND)
          phi=cmp%value(ky_ANGL_BEND)
          if(cmp%ori)then
            psi1=cmp%value(ky_E1_BEND)*phi+cmp%value(ky_AE1_BEND)
            psi2=cmp%value(ky_E2_BEND)*phi+cmp%value(ky_AE2_BEND)
            fb1=cmp%value(ky_F1_BEND)+cmp%value(ky_FB1_BEND)
            fb2=cmp%value(ky_F1_BEND)+cmp%value(ky_FB2_BEND)
          else
            psi1=cmp%value(ky_E2_BEND)*phi+cmp%value(ky_AE2_BEND)
            psi2=cmp%value(ky_E1_BEND)*phi+cmp%value(ky_AE1_BEND)
            fb2=cmp%value(ky_F1_BEND)+cmp%value(ky_FB1_BEND)
            fb1=cmp%value(ky_F1_BEND)+cmp%value(ky_FB2_BEND)
          endif
          if(cmp%value(ky_FRMD_BEND) <= 0.d0)then
            if(cmp%value(ky_FRMD_BEND) /= -1.d0)then
              fb1=0.d0
            endif
            if(cmp%value(ky_FRMD_BEND) /= -2.d0)then
              fb2=0.d0
            endif
          endif
          w=phi-psi1-psi2
          if((fb1 /= 0.d0 .or. fb2 /= 0.d0) .and.
     1         al /= 0.d0 .and. phi /= 0.d0)then
            al=al-((phi*fb1)**2+(phi*fb2)**2)/al/48.d0
     1           *sin(.5d0*w)/sin(.5d0*phi)
          endif
          cmp%value(p_L_BEND)=al
          dtheta=cmp%value(ky_DROT_BEND)
          theta=cmp%value(ky_ROT_BEND)
          cmp%value(p_PSI1_BEND)=psi1
          cmp%value(p_PSI2_BEND)=psi2
          cmp%value(p_COSPSI1_BEND)=cos(psi1)
          cmp%value(p_SINPSI1_BEND)=sin(psi1)
          cmp%value(p_COSPSI2_BEND)=cos(psi2)
          cmp%value(p_SINPSI2_BEND)=sin(psi2)
          cmp%value(p_COSTHETA_BEND)=cos(theta)
          cmp%value(p_SINTHETA_BEND)=sin(theta)
          cmp%value(p_COSW_BEND)=cos(w)
          cmp%value(p_SINW_BEND)=sin(w)
          cmp%value(p_SQWH_BEND)=merge(cmp%value(p_SINW_BEND)**2
     $         /(1.d0+cmp%value(p_COSW_BEND)),
     $         1.d0-cmp%value(p_COSW_BEND),
     $         cmp%value(p_COSW_BEND) .ge. 0.d0)
          cmp%value(p_SINWP1_BEND)=sin(phi-psi2)
c          cmp%value(p_DPHIX_BEND)=phi*sin(.5d0*dtheta)**2
c          cmp%value(p_DPHIY_BEND)=.5d0*phi*sin(dtheta)
          cmp%value(p_THETA_BEND)=theta
          cmp%value(p_FB1_BEND)=fb1
          cmp%value(p_FB2_BEND)=fb2
          cmp%ivalue(1,p_FRMD_BEND)=int(cmp%value(ky_FRMD_BEND))

        case (icQUAD)
          al=cmp%value(ky_L_QUAD)
          theta=cmp%value(ky_ROT_QUAD)
          cmp%value(p_THETA2_QUAD)=theta
     $         +akang(dcmplx(cmp%value(ky_K1_QUAD),0.d0),al,cr1)
          if(al /= 0.d0)then
            akk=cmp%value(ky_K1_QUAD)/al
            call tsetfringep(cmp,icQUAD,akk,
     $           cmp%value(p_AKF1F_QUAD:p_AKF2B_QUAD))
          else
            cmp%value(p_AKF1F_QUAD:p_AKF2B_QUAD)=0.d0
          endif
          frmd=cmp%value(ky_FRMD_QUAD)
          if(.not. cmp%ori)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%ivalue(1,p_FRMD_QUAD)=int(frmd)

        case (icSEXT,icOCTU,icDECA,icDODECA)
          theta=cmp%value(ky_ROT_THIN)
          cmp%value(p_COSTHETA_THIN)=cos(theta)
          cmp%value(p_SINTHETA_THIN)=sin(theta)
          cmp%value(p_THETA_THIN)=theta

        case (icUND)
          call undinit(cmp%value(1),cmp%value(p_PARAM_UND))

        case (icWIG)
          call twigp()

        case (icMULT)
          al=cmp%value(ky_L_MULT)
          phi=cmp%value(ky_ANGL_MULT)
          if(cmp%ori)then
            psi1=cmp%value(ky_E1_MULT)*phi+cmp%value(ky_AE1_MULT)
            psi2=cmp%value(ky_E2_MULT)*phi+cmp%value(ky_AE2_MULT)
            fb1=cmp%value(ky_FB1_MULT)
            fb2=cmp%value(ky_FB2_MULT)
          else
            psi1=cmp%value(ky_E2_MULT)*phi+cmp%value(ky_AE2_MULT)
            psi2=cmp%value(ky_E1_MULT)*phi+cmp%value(ky_AE1_MULT)
            fb2=cmp%value(ky_FB1_MULT)
            fb1=cmp%value(ky_FB2_MULT)
          endif
          frmd=cmp%value(ky_FRMD_MULT)
          if(.not. cmp%ori)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%ivalue(1,p_FRMD_MULT)=int(frmd)
          if(frmd /= 3.d0 .and. frmd /= 1.d0)then
            fb1=0.d0
          endif
          if(frmd /= 3.d0 .and. frmd /= 2.d0)then
            fb2=0.d0
          endif
          w=phi-psi1-psi2
          if((fb1 /= 0.d0 .or. fb2 /= 0.d0) .and.
     1         al /= 0.d0 .and. phi /= 0.d0)then
            al=al-((phi*fb1)**2+(phi*fb2)**2)/al/48.d0
     $           *sin(.5d0*w)/sin(.5d0*phi)
          endif
          cmp%value(p_L_MULT)=al
          cmp%value(p_ANGL_MULT)=phi
          cmp%value(p_PSI1_MULT)=psi1
          cmp%value(p_PSI2_MULT)=psi2
          cmp%value(p_FB1_MULT)=fb1
          cmp%value(p_FB2_MULT)=fb2
          if(al /= 0.d0)then
            sk1=cmp%value(ky_SK1_MULT)
            if(sk1 == 0.d0)then
              akk=cmp%value(ky_K1_MULT)/al
            else
              akk=hypot(cmp%value(ky_K1_MULT),sk1)/al
c              akk=sqrt(cmp%value(ky_K1_MULT)**2+sk1**2)/al
            endif
            call tsetfringep(cmp,icMULT,akk,
     $           cmp%value(p_AKF1F_MULT:p_AKF2B_MULT))
          else
            cmp%value(p_AKF1F_MULT:p_AKF2B_MULT)=0.d0
          endif
          harm=cmp%value(ky_HARM_MULT)
          w=merge(pi2*cmp%value(ky_FREQ_MULT)/c,
     $         omega0*harm/c,harm == 0.d0)
          vnominal=merge(cmp%value(ky_VOLT_MULT)/amass*abs(charge)
     $         *sin(-cmp%value(ky_PHI_MULT)*sign(1.d0,charge)),
     $         0.d0,trpt)
          cmp%value(p_W_MULT)=w
          cmp%value(p_VNOMINAL_MULT)=vnominal
          cmp%value(p_THETA2_MULT)=cmp%value(ky_ROT_MULt)+
     $         cmp%value(ky_DROT_MULT)+
     $         akang(dcmplx(cmp%value(ky_K1_MULT),
     $         cmp%value(ky_SK1_MULT)),al,cr1)
          cmp%value(p_CR1_MULT)=dble(cr1)
          cmp%value(p_CR1I_MULT)=imag(cr1)

          nmmax=-1
          do n=nmult,0,-1
            if(cmp%value(ky_K0_MULT+n*2) /= 0.d0
     $           .or. cmp%value(ky_SK0_MULT+n*2) /= 0.d0)then
              nmmax=n
              exit
            endif
          enddo
          cmp%ivalue(2,p_NM_MULT)=merge(nmmax,max(nmmax,0),al == 0.d0
     $         .and. cmp%value(ky_VOLT_MULT) == 0.d0)
          cr=cr1
          do n=0,nmmax
            ck=dcmplx(cmp%value(ky_K0_MULT+n*2),
     $           cmp%value(ky_SK0_MULT+n*2))*cr
            cmp%value(p_K0R_MULT+n*2)=dble(ck)
            cmp%value(p_K0R_MULT+n*2+1)=imag(ck)
            cr=cr*cr1
          enddo
          eps=eps00m*merge(1.d0,cmp%value(ky_EPS_MULT),
     $         cmp%value(ky_EPS_MULT) == 0.d0)
          ndiv=1
          if(al /= 0.d0)then
            do n=2,nmmax
              ndiv=max(ndiv,int(
     $             sqrt(ampmaxm**(n-1)/6.d0/fact(n-1)/eps*
     $             hypot(cmp%value(ky_K0_MULT+n*2),
     $             cmp%value(ky_SK0_MULT+n*2))*al))+1)
            enddo
          endif
          cmp%ivalue(1,p_NM_MULT)=min(ndiv,ndivmaxm)

          do n=0,nmmax
            cmp%lvalue(n+1,p_DOFR_MULT)=
     $           merge(merge(.true.,
     $           hypot(cmp%value(ky_K0_MULT+n*2),
     $           cmp%value(ky_SK0_MULT+n*2))/al*ampmaxm**(n+1)*aninv(n+1)
     $           > eps,n <= 2),.false.,
     $           cmp%value(ky_K0_MULT+n*2) /= 0.d0
     $           .or. cmp%value(ky_SK0_MULT+n*2) /= 0.d0)
          enddo

        case (icCAVI)
          frmd=cmp%value(ky_FRMD_CAVI)
          if(.not. cmp%ori)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%ivalue(1,p_FRMD_CAVI)=int(frmd)
          harm=cmp%value(ky_HARM_CAVI)
          w=merge(pi2*cmp%value(ky_FREQ_CAVI)/c,omega0*harm/c,
     $         harm == 0.d0)
          vnominal=merge(cmp%value(ky_VOLT_CAVI)/amass*abs(charge)
     $         *sin(-cmp%value(ky_PHI_CAVI)*sign(1.d0,charge)),
     $         0.d0,trpt)
          cmp%value(p_W_CAVI)=w
          cmp%value(p_VNOMINAL_CAVI)=vnominal

        case (icBEAM)
          call bbinit(cmp%value(1),cmp%value(p_PARAM_BEAM))

        case (icPROT)
          call phsinit(cmp%value(1),cmp%value(p_PARAM_Prot))

        case (icSOL)
          cchi1=cos(cmp%value(ky_CHI1_SOL))
          cchi2=cos(cmp%value(ky_CHI2_SOL))
          cchi3=cos(cmp%value(ky_CHI3_SOL))
          schi1=sin(cmp%value(ky_CHI1_SOL))
          schi2=sin(cmp%value(ky_CHI2_SOL))
          schi3=sin(cmp%value(ky_CHI3_SOL))
          cmp%value(p_R11_SOL)= cchi1*cchi3+schi1*schi2*schi3
          cmp%value(p_R12_SOL)=-cchi2*schi3
          cmp%value(p_R13_SOL)= schi1*cchi3-cchi1*schi2*schi3
          cmp%value(p_R21_SOL)=-schi1*schi2*cchi3+cchi1*schi3
          cmp%value(p_R22_SOL)= cchi2*cchi3
          cmp%value(p_R23_SOL)= cchi1*schi2*cchi3+schi1*schi3
          cmp%value(p_R31_SOL)=-schi1*cchi2
          cmp%value(p_R32_SOL)=-schi2
          cmp%value(p_R33_SOL)= cchi1*cchi2

        case default
        end select
        return
        end subroutine

        subroutine tclrpara
        use mackw
        use tmacro
        use ffs_pointer, only:compelc
        implicit none
        type (sad_comp), pointer :: cmp
        integer*4 i
        do i=1,nlat-1
          call compelc(i,cmp)
          cmp%update=cmp%nparam <= 0
        enddo
        tparaed=.false.
        return
        end subroutine

      end module

      module ffs_fit
      use ffs, only:ntwissfun,maxcond
      use tffitcode, only:mfit1
      type ffs_stat
        real*8 tracex,tracey,tracez
        logical*4 stabx,staby,stabz,over
      end type
      type ffs_res
      sequence
      real*8 r
      integer*4 nstab
      end type

      integer*4 , parameter :: ndimmax=500
      type (ffs_stat) optstat(-ndimmax:ndimmax)
      type (ffs_res) residual(-ndimmax:ndimmax)
      integer*4 iuid(-ndimmax:ndimmax),
     $     jfam(-ndimmax:ndimmax),kfam(-ndimmax:ndimmax)
      real*8 dp(-ndimmax:ndimmax),scale(mfit1),
     $     dfam(4,-ndimmax:ndimmax),
     $     uini(ntwissfun,-ndimmax:ndimmax),wfit(maxcond),
     $     wiq(maxcond)
      logical*4 fitflg,geomet,inicond,chgini
      integer*4 nut,nfam,nfam1,nfr,nqcol,nqcol1,nfcol,nfc0
      real*8 wexponent,offmw,etamax,avebeta,wsum
      character*8 , save :: nlist(1:mfit1)=(/
     $     'AX      ','BX      ','NX      ',
     $     'AY      ','BY      ','NY      ',
     $     'EX      ','EPX     ','EY      ','EPY     ',
     $     'R1      ','R2      ','R3      ','R4      ','DETR    ',
     $     'DX      ','DPX     ','DY      ','DPY     ',
     $     'DZ      ','DDP     ',
     $     'AZ      ','BZ      ','NZ      ',
     $     'ZX      ','ZPX     ','ZY      ','ZPY     ',
     1     'PEX     ','PEPX    ','PEY     ','PEPY    ',
     1     'PZX     ','PZPX    ','PZY     ','PZPY    ',
     $     'GMX     ','GMY     ','GMZ     ',
     $     'BMAGX   ','BMAGY   ','BMAGZ   ',
     $     'TRX     ','TRY     ','TRZ     ',
     $     'LENG    ',
     $     'GX      ','GY      ','GZ      ',
     $     'CHI1    ','CHI2    ','CHI3    ',
     $     'DEX     ','DEPX    ','DEY     ','DEPY    ',
     $     'DDX     ','DDPX    ','DDY     ','DDPY    ',
     $     'PDEX    ','PDEPX   ','PDEY    ','PDEPY   '/)

      contains
      logical*4 function resle(r1,r2,tol) result(f)
      type (ffs_res) ,intent(in):: r1,r2
      real*8 ,intent(in),optional::tol
      if(present(tol))then
        f=r1%nstab < r2%nstab .or. r1%nstab == r2%nstab .and. r1%r <= r2%r/tol
      else
        f=r1%nstab < r2%nstab .or. r1%nstab == r2%nstab .and. r1%r <= r2%r
      endif
      return
      end function

      real*8 function res2r(r1) result(f)
      type (ffs_res),intent(in):: r1
      real*8 rb
      rb=10.d0
      do while (rb < r1%r)
        rb=rb*10.d0
      enddo
      f=rb*r1%nstab+r1%r
      return
      end function

      end module

      module ffs_wake
      integer*8 kwakep,kwakeelm,iwpl,iwpt
      integer*4 nwakep,lwl0,lwt0
      logical*4 wake
      integer*4 nwp
      integer*8 , pointer :: kwaketbl(:,:)
      integer*4 , pointer :: iwakeelm(:)
      integer*4 ,allocatable,dimension(:)::itab,izs
      real*8 ,pointer :: wakel(:,:),waket(:,:)

      contains
      subroutine twxiwp(nwak)
      use tfstk
      use iso_c_binding
      implicit none
      integer*4 ,intent(in):: nwak
      if(nwak > 0)then
        iwpl=abs(kwaketbl(1,nwak))
        lwl0=merge((ilist(1,iwpl-1)-2)/2,0,iwpl /= 0)
        iwpt=abs(kwaketbl(2,nwak))
        lwt0=merge((ilist(1,iwpt-1)-2)/2,0,iwpt /= 0)
        if(lwl0 > 0)then
          call c_f_pointer(c_loc(rlist(iwpl)),wakel,[2,lwl0])
        endif
        if(lwt0 > 0)then
          call c_f_pointer(c_loc(rlist(iwpt)),waket,[2,lwt0])
        endif
      endif
      return
      end subroutine

      subroutine twxalloc(np)
      implicit none
      integer*4 ,intent(in):: np
      if(.not. allocated(itab))then
        nwp=np
        allocate(itab(nwp),izs(nwp))
      elseif(np > nwp)then
        deallocate(itab,izs)
        nwp=np
        allocate(itab(nwp),izs(nwp))
      endif
      itab(np)=np
      itab(1)=1
      return
      end subroutine

      subroutine twxdealloc
      implicit none
      if(allocated(itab))then
        deallocate(itab,izs)
      endif
      nwp=0
      return
      end subroutine 

      end module

      module tflinepcom
      use tfstk
      implicit none
      integer*8, save :: iflinep=0,ifinitlinep=0,ifelementp=0,
     $     ifelementkeyp=0,iftypekeyp=0
      type (sad_descriptor) ksdumm

      contains

      integer*4 function itftypekey(i,word1,lw1,kx)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) , optional,intent(out)::kx
      type (sad_descriptor) ks,ka
      integer*4 ,intent(in):: i,lw1
      integer*4 isp0,irtc
      character*(lw1) ,intent(in):: word1
      isp0=isp
      ktastk(isp0+1)=iftypekeyp
      rtastk(isp0+2)=dble(i)
      isp=isp0+3
      ks=kxsalocb(-1,word1,lw1)
      dtastk(isp)=ks
      levele=levele+1
      ka=tfefunref(isp0+1,.false.,irtc)
      call tfconnect(ka,irtc)
      isp=isp0
      if(irtc /= 0)then
        call tfreseterror
        itftypekey=0
      elseif(ktflistq(ka))then
        itftypekey=-1
        if(present(kx))then
          kx=ka
        endif
      else
        itftypekey=merge(ifromd(ka),0,ktfrealq(ka))
      endif        
      return
      end function

      integer*4 function itfelementkey(i,word1,lw1)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx,ks
      integer*4 ,intent(in):: i,lw1
      integer*4 irtc,isp0
      character*(lw1) ,intent(in):: word1
      isp0=isp
      ktastk(isp0+1)=ifelementkeyp
      rtastk(isp0+2)=dble(i)
      isp=isp0+3
      ks=kxsalocb(-1,word1,lw1)
      dtastk(isp)=ks
      levele=levele+1
      kx=tfefunref(isp0+1,.false.,irtc)
      call tfconnect(kx,irtc)
      isp=isp0
      if(irtc /= 0)then
        call tfreseterror
        itfelementkey=0
      else
        itfelementkey=merge(ifromd(kx),0,ktfrealq(kx))
      endif        
      return
      end function

      function tftypekey(isp1,irtc) result(kx)
      use tfstk
      use maccbk
      use maccode
      use mackw
      implicit none
      type (sad_descriptor) kx
      real*8 c
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 ic,isp0,itfmessage
      if(isp /= isp1+1)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp),c))then
        kx=dxnullo
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Typecode for #1"')
        return
      endif
      ic=int(c)
      if(ic < 0 .or. ic > icMXEL)then
        kx=dxnullo
        irtc=itfmessage(9,'General::wrongval',
     $       '"Typecode out of range"')
        return
      endif
      isp0=isp
      call tftypekeystk(ic,.true.)
      kx=kxmakelist(isp0)
      isp=isp0
      irtc=0
      return
      end function

      subroutine tftypekeystk(ic,all)
      use tfstk
      use mackw
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: ic
      integer*4 i,isp0
      logical*4 ,intent(in):: all
      character*(MAXPNAME) key,tfkwrd,tfkwrd1
      do i=1,kytbl(kwMAX,ic)-1
        key=tfkwrd(ic,i)
        if(all .or. key /= '-')then
          isp=isp+1
          if(key == '-')then
            dtastk(isp)=ksdumm
          else
            dtastk(isp)=kxsalocb(-1,key,len_trim(key))
            if(kyindex1(i,ic) /= 0)then
              key=tfkwrd1(ic,i)
              isp0=isp
              isp=isp+1
              dtastk(isp)=kxsalocb(-1,key,len_trim(key))
              kx=kxmakelist(isp0-1)
              isp=isp0
              dtastk(isp)=kx
            endif
          endif
        endif
      enddo
      return
      end subroutine

      type (sad_descriptor) function tfkeyv(i,keyword,ia,cmp,ref,saved)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,elatt,idtypec,idvalc,sad_comp,
     $     compelc,iele1
      implicit none
      type (sad_rlist), pointer :: kld,klr
      integer*4 ,intent(in):: i
      integer*4 it,kl,l,j,lk,lenw
      integer*8 ,intent(out):: ia
      character*(*) keyword
      character*128 key
      logical*4 plus
      logical*4 ,intent(in):: ref,saved
      real*8 s
      type (sad_comp), pointer ,intent(out):: cmp
c     begin initialize for preventing compiler warning
      kl=merge(i,nelvx(-i)%klp,i > 0)
      call compelc(kl,cmp)
      plus=.false.
      lk=lenw(keyword)
      if(lk > 4 .and. keyword(lk-3:lk) == '$SUM')then
        plus=.true.
        lk=lk-4
      endif
      key(1:lk)=keyword(1:lk)
      it=idtype(cmp%id)
      do l=1,kytbl(kwMAX,it)-1
        j=kyindex(l,it)
        if(j > 0)then
          if(pname(kytbl(j,0))(2:) == key(1:lk))then
            go to 1
          endif
          j=kyindex1(l,it)
          if(j > 0)then
            if(pname(kytbl(j,0))(2:) == key(1:lk))then
              go to 1
            endif
          endif
        endif
      enddo
      ia=0
      tfkeyv%k=0
      return
 1    if(i > 0)then
        ia=elatt%comp(i)+l
        tfkeyv=dlist(ia)
        if(.not. ref)then
          call tftouch(iele1(i),l)
        endif
      elseif(saved)then
        ia=idvalc(kl)+l
        tfkeyv=dlist(ia)
      else
        ia=elatt%comp(kl)+l
        if(l == nelvx(-i)%ival)then
          if(tfreallistq(dlist(ia),klr))then
            tfkeyv=kxavaloc(-1,klr%nl,kld)
            kld%rbody(1:klr%nl)=klr%rbody(1:klr%nl)
     $           /rlist(iferrk+(kl-1)*2)
          else
            tfkeyv=dlist(ia)
          endif
        else
          tfkeyv=dlist(ia)
        endif
        if(.not. ref)then
          call tftouch(iele1(kl),l)
        endif
      endif
      if(plus .and. tfreallistq(tfkeyv,klr))then
        s=sum(klr%rbody(1:klr%nl))
        tfkeyv=dfromr(s)
      endif
      return
      end function

      subroutine tftouch(i,iv)
      use ffs
      use ffs_pointer
      implicit none
      integer*4 ,intent(in):: iv,i
      integer*4 j
      do j=1,flv%ntouch
        if(nvevx(j)%itouchele == i .and. nvevx(j)%itouchv == iv)then
          return
        endif
      enddo
      call tffsnvealloc(flv%ntouch+1)
      flv%ntouch=flv%ntouch+1
      nvevx(flv%ntouch)%itouchele=i
      nvevx(flv%ntouch)%itouchv=iv
      return
      end subroutine

      subroutine elcompl(i,kl)
      use tfstk
      use ffs
      use ffs_pointer
      implicit none
      type (sad_rlist), pointer,intent(out) :: kl
      integer*4 ,intent(in)::i
      integer*4 l,isp1
      if(nelvx(i)%dcomp%k == 0)then
        isp1=isp
        do l=1,nlat-1
          if(iele1(l) == i)then
            isp=isp+1
            rtastk(isp)=dble(l)
          endif
        enddo
        nelvx(i)%dcomp=dtfcopy(kxmakelist(isp1,kl))
        isp=isp1
      else
        call descr_sad(nelvx(i)%dcomp,kl)
      endif
      return
      end subroutine

      end module

      module ffs_seg
      contains
        logical*4 function tcheckseg(cmp,ltyp,al,lsegp,irtc) result(seg)
        use tfstk
        use mackw
        use sad_main
        use kyparam
        implicit none
        type (sad_comp) ,intent(inout)::cmp
        type (sad_dlist) , pointer :: lprof
        type (sad_dlist) , pointer ,intent(inout):: lsegp
        real*8 , intent(out)::al
        integer*4 ,intent(out)::irtc
        integer*4 ,intent(in)::ltyp
        integer*4 kl,kprof
        irtc=0
        seg=.false.
        kl=kytbl(kwL,ltyp)
        if(kl == 0)then
          al=0.d0
          return
        else
          al=cmp%value(kl)
        endif
        kprof=kytbl(kwPROF,ltyp)
        if(kprof == 0 .or. tfnonlistq(cmp%dvalue(kprof),lprof))then
          return
        endif
        if(.not. cmp%updateseg)then
          call tsetupseg(cmp,lprof,lsegp,irtc)
          if(irtc /= 0)then
            return
          endif
        else
          call descr_sad(cmp%dvalue(p_PROF_MULT),lsegp)
        endif
        seg=.true.
        return
        end function

        subroutine tsetupseg(cmp,lprof,lsegp,irtc)
        use tfstk
        use kyparam
        use mackw
        use tflinepcom
        use sad_main
        implicit none
        type (sad_comp) ,intent(inout):: cmp
        type (sad_rlist) , pointer :: lvl
        type (sad_dlist) , pointer ,intent(in):: lprof
        type (sad_dlist) , pointer ,intent(out):: lsegp
        type (sad_dlist) , pointer :: lpi,lpi1,lk1
        type (sad_string), pointer :: stri1,stri2
        integer*4 ltyp,lls,i,nseg,irtc,itfmessage,l,isp0,i1,nk
        integer*4 ,parameter :: nc=ky_PROF_MULT-1
        logical*4 defk(nc,nc)
        irtc=0
        ltyp=idtype(cmp%id)
        isp0=isp
        select case (ltyp)
        case (icMULT)
          defk=.false.
          nseg=0
          lls=0
          levele=levele+1
          do i=1,lprof%nl
            if(tfnonlistq(lprof%dbody(i),lpi) .or. lpi%nl < 2)then
              go to 9100
            endif
            if(tflistq(lpi%dbody(1),lpi1))then
              if(lpi1%nl > 2 .or. lpi1%nl <= 0)then
                go to 9100
              endif
              if(.not. ktfstringq(lpi1%dbody(1),stri1))then
                go to 9100
              endif
              if(lpi1%nl == 1)then
                stri2=>stri1
              elseif(.not. ktfstringq(lpi1%dbody(2),stri2))then
                go to 9100
              endif
            elseif(ktfstringq(lpi%dbody(1),stri1))then
              stri2=>stri1
            else
              go to 9100
            endif
            if(tfnonreallistq(lpi%dbody(2),lvl))then
              go to 9100
            endif
            if(nseg == 0)then
              nseg=lvl%nl
            elseif(lvl%nl /= nseg)then
              irtc=itfmessage(99,"FFS::unequalkeyleng",'""')
              go to 9000
            endif
            isp=isp+1
            itastk(1,isp)=itftypekey(icMULT,stri1%str,stri1%nch)
            itastk(2,isp)=itftypekey(icMULT,stri2%str,stri2%nch)
            if(itastk(1,isp) <= 0 .or. itastk(1,isp) > nc .or.
     $           itastk(2,isp) <= 0 .or. itastk(2,isp) > nc)then
              irtc=itfmessage(99,"FFS::wrongkey",'""')
              go to 9000
            elseif(defk(itastk(1,isp),itastk(2,isp)))then
              call tfdebugprint(lprof%dbody(i),
     $             'Dupricated Key in '//
     $             pname(cmp%id)(1:lpname(cmp%id)),3)
              irtc=itfmessage(99,"FFS::dupkey",'""')
              go to 9000
            endif
            defk(itastk(1,isp),itastk(2,isp))=.true.
            dtastk2(isp)=lpi%dbody(2)
            if(stri1%str(1:stri1%nch) == 'L')then
              lls=isp
            endif
          enddo
          if(lls == 0)then
            irtc=itfmessage(99,"FFS::noLseg",'""')
            go to 9000
          endif
          nk=isp-isp0
          if(tfnonlistq(cmp%dvalue(p_PROF_MULT),lsegp) .or.
     $         lsegp%nl /= nk)then
            call tflocald(cmp%dvalue(p_PROF_MULT))
            cmp%dvalue(p_PROF_MULT)=kxadaloc(0,nk,lsegp)
            lsegp%dbody(1:nk)%k=ktfoper+mtfnull
          endif
          call tflocald(lsegp%dbody(1))
          lsegp%dbody(1)=kxadaloc(0,2,lk1)
          lk1%dbody(1)=dtastk(lls)
          lk1%dbody(2)=dtfcopy(dtastk2(lls))
          i1=1
          do i=isp0+1,isp
            if(i /= lls .and. itastk(1,i) == itastk(2,i))then
              i1=i1+1
              call tflocald(lsegp%dbody(i1))
              lsegp%dbody(i1)=kxadaloc(0,2,lk1)
              lk1%dbody(1)=dtastk(i)
              lk1%dbody(2)=dtfcopy(dtastk2(i))
            endif
          enddo
          do i=isp0+1,isp
            if(i /= lls .and. itastk(1,i) /= itastk(2,i))then
              i1=i1+1
              call tflocald(lsegp%dbody(i1))
              lsegp%dbody(i1)=kxadaloc(0,2,lk1)
              lk1%dbody(1)=dtastk(i)
              lk1%dbody(2)=dtfcopy(dtastk2(i))
            endif
          enddo
          cmp%updateseg=.true.
 9000     l=itfdownlevel()
          isp=isp0
          return
 9100     l=itfdownlevel()
          isp=isp0
          irtc=itfmessage(99,"FFS::wrongkeylist",'""')
        case default
          cmp%updateseg=.true.
        end select
        return
        end

        subroutine tfvcopycmp(cmps,cmpd,k,coeff)
        use tfstk
        use mackw
        use kyparam
        use sad_main
        implicit none
        type (sad_comp),intent(in):: cmps
        type (sad_comp),intent(out):: cmpd
        real*8 ,intent(in):: coeff
        integer*4 ,intent(in):: k
        if(ktfrealq(cmps%dvalue(k)))then
          if(ktfnonrealq(cmpd%dvalue(k)))then
            call tflocald(cmpd%dvalue(k))
          endif
          cmpd%value(k)=cmps%value(k)*coeff
        elseif(tflistq(cmps%dvalue(k)))then
          if(ktfnonrealq(cmpd%dvalue(k)))then
            call tflocald(cmpd%dvalue(k))
          endif
          cmpd%dvalue(k)=dtfcopy1(cmps%dvalue(k))
          if(idtype(cmpd%id) == icMULT .and.
     $         k == ky_PROF_MULT)then
            cmpd%updateseg=.false.
          endif
        endif
        cmpd%update=cmpd%nparam <= 0
        return
        end subroutine

        subroutine tfvcopy(is,id,k,coeff)
        use tfstk
        use mackw
        use kyparam
        use sad_main
        use ffs_pointer
        implicit none
        type (sad_comp), pointer :: cmps,cmpd
        real*8 ,intent(in):: coeff
        integer*4 ,intent(in):: is,id,k
        call compelc(id,cmpd)
        call compelc(is,cmps)
        call tfvcopycmp(cmps,cmpd,k,coeff)
        if(idtype(cmpd%id) == icMULT .and.
     $       k == ky_PROF_MULT)then
          cmpd%updateseg=.false.
        endif
        cmpd%update=cmpd%nparam <= 0
        return
        end subroutine

        subroutine tfvcopycmpall(cmps,cmpd,n)
        use tfstk
        use mackw
        use sad_main
        implicit none
        type (sad_comp),intent(in):: cmps
        type (sad_comp),intent(out):: cmpd
        integer*4 ,intent(in):: n
        integer*4 k
        do k=1,n
          call tfvcopycmp(cmps,cmpd,k,1.d0)
        enddo
        return
        end subroutine

        real*8 function tfvcmp(cmps,k) result(v)
        use mackw
        use tfstk
        use sad_main
        implicit none
        type (sad_comp) ,intent(in):: cmps
        type (sad_rlist), pointer :: las
        integer*4 ,intent(in):: k
        if(ktfnonrealq(cmps%dvalue(k),v))then
          if(tfreallistq(cmps%dvalue(k),las))then
            v=sum(las%rbody(1:las%nl))
          endif
        endif
        return
        end function

        subroutine tfsetcmp(v,cmpd,i)
        use mackw
        use tfstk
        use sad_main
        implicit none
        type (sad_comp) :: cmpd
        type (sad_rlist), pointer :: lad
        integer*4 ,intent(in):: i
        real*8 r0,v
        if(ktfrealq(cmpd%dvalue(i)))then
          cmpd%value(i)=v
        elseif(tfreallistq(cmpd%dvalue(i),lad))then
          r0=sum(lad%rbody(1:lad%nl))
          lad%rbody(1:lad%nl)=merge(v/r0*lad%rbody(1:lad%nl),
     $         v/lad%nl,r0 /= 0.d0)
        endif
        return
        end subroutine

      end module

      module temw
      use tffitcode
      implicit none

      private

      integer*4 , public, parameter ::
     $     ipdx=1,ipdpx=2,ipdy=3,ipdpy=4,ipdz=5,ipddp=6,
     $     ipnx=7,ipny=8,ipnz=9,
     $     ipu0=10,ipvceff=11,iptrf0=12,ipalphap=13,ipdleng=14,
     $     ipbh=15,ipdampx=16,ipdampy=17,ipdampz=18,
     $     ipjx=19,ipjy=20,ipjz=21,
     $     ipemx=22,ipemy=23,ipemz=24,ipsige=25,ipsigz=26,
     $     ipheff=27,ipdnux=28,ipdnuy=29,ipdnuz=30,
     $     iptwiss=31,iptws0=iptwiss-1,
     $     ipnnup=iptwiss+ntwissfun,iptaup=ipnnup+1,
     $     ippolx=iptaup+1,ippoly=ippolx+1,ippolz=ippoly+1,
     $     ipequpol=ippolz+1,ipequpol2=ipequpol+1,
     $     ipequpol4=ipequpol2+1,ipequpol6=ipequpol4+1,
     $     ipnup=ipequpol6+1,iprevf=ipnup+1,
     $     nparams=iprevf

c     inverse normal ccordinate at the origin (normal -> physical)
      real(8), public :: r(6, 6) = RESHAPE((/
     $     1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/),
     $     (/6, 6/))

c     Inverse matrix of r: NormalCoordinate (physical -> normal)
      real(8), public :: ri(6, 6) = RESHAPE((/
     $     1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0,
     $     0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/),
     $     (/6, 6/))

      real*8 , public::dsg(6)

      type iaemit
      sequence
      integer*8 iatr,iacod,iamat,iabmi
      end type

      type (iaemit) ,public,parameter::
     $     iaez=iaemit(int8(0),int8(0),int8(0),int8(0))
c     $     iaez%iatr=int8(0),iaez%iacod=int8(0),
c     $     iaez%iamat=int8(0),iaez%iabmi=int8(0)

      real*8 , public :: transr(6,6),codr0(6),bzhr0,bsir0
      real*8 , public :: gintd(3)
      real(8), public :: eemx,eemy,eemz,
     $     emxmin,emymin,emzmin,emxmax,emymax,emzmax,
     $     emxr,emyr,emzr,emx0,emy0,emz0,demin,
     $     equpol(3),spinmu,code(6)
      real(8), public :: bh,heff,phirf,omegaz,sige,sigz
      complex*16 , public :: ceig0(6)
      logical*4, public :: normali, initemip=.true.,
     $     initr,calint,caltouck,fndcod,synchm,epi,beamplt
      real*8 , parameter :: toln=0.1d0

      character*10, parameter ::
     $     label1(6)=['        X ','       Px ','        Y ',
     1                '       Py ','        Z ','       Pz '],
     $     label2(6)=['        x ','    px/p0 ','        y ',
     1                '    py/p0 ','        z ','    dp/p0 ']

      public :: tfetwiss,etwiss2ri,tfnormalcoord,toln,
     $     tfinitemip,tsetr0,
     $     label1,label2,tfinibeam,tmulbs,iaemit

      contains

      subroutine tsetr0(trans,cod,bzh,bsi0)
      implicit none
      real*8 ,intent(in)::trans(6,6),cod(6),bzh,bsi0
c     cod are canonical momenta!
      codr0=cod
      transr=trans
      bzhr0=bzh
      bsir0=bsi0
      return
      end subroutine

      subroutine tfinitemip
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc
      character*2 strfromis
      if( .not. initemip)then
        return
      endif
      call tfevals('BeginPackage[emit$`];Begin[emit$`]',kx,irtc)
      call tfevals(
     $'nparams='//strfromis(nparams)//
     $';ipdx='//strfromis(ipdx)//
     $';ipdpx='//strfromis(ipdpx)//
     $';ipdy='//strfromis(ipdy)//
     $';ipdpy='//strfromis(ipdpy)//
     $';ipdz='//strfromis(ipdz)//
     $';ipddp='//strfromis(ipddp)//
     $';ipnx='//strfromis(ipnx)//
     $';ipny='//strfromis(ipny)//
     $';ipnz='//strfromis(ipnz)//
     $';ipu0='//strfromis(ipu0)//
     $';ipvceff='//strfromis(ipvceff)//
     $';iptrf0='//strfromis(iptrf0)//
     $';ipalphap='//strfromis(ipalphap)//
     $';ipdleng='//strfromis(ipdleng)//
     $';ipbh='//strfromis(ipbh)//
     $';ipheff='//strfromis(ipheff)//
     $';iptwiss='//strfromis(iptwiss)//
     $';iptws0='//strfromis(iptws0)//
     $';ipdampx='//strfromis(ipdampx)//
     $';ipdampy='//strfromis(ipdampy)//
     $';ipdampz='//strfromis(ipdampz)//
     $';ipdnux='//strfromis(ipdnux)//
     $';ipdnuy='//strfromis(ipdnuy)//
     $';ipdnuz='//strfromis(ipdnuz)//
     $';ipjx='//strfromis(ipjx)//
     $';ipjy='//strfromis(ipjy)//
     $';ipjz='//strfromis(ipjz)//
     $';ipemx='//strfromis(ipemx)//
     $';ipemy='//strfromis(ipemy)//
     $';ipemz='//strfromis(ipemz)//
     $';ipsige='//strfromis(ipsige)//
     $';ipsigz='//strfromis(ipsigz)//
     $';ipnnup='//strfromis(ipnnup)//
     $';iptaup='//strfromis(iptaup)//
     $';ippolx='//strfromis(ippolx)//
     $';ippoly='//strfromis(ippoly)//
     $';ippolz='//strfromis(ippolz)//
     $';ipequpol='//strfromis(ipequpol)//
     $';ipequpol2='//strfromis(ipequpol2)//
     $';ipequpol4='//strfromis(ipequpol4)//
     $';ipequpol6='//strfromis(ipequpol6)//
     $';ipnup='//strfromis(ipnup)//
     $';iprevf='//strfromis(iprevf)//
     $';SetAttributes[{'//
     $ 'nparams,ipdx,ipdpx,ipdy,ipdpy,ipdz,ipddp,'//
     $ 'ipnx,ipny,ipnz,'//
     $ 'ipu0,ipvceff,iptrf0,ipalphap,ipdleng,'//
     $ 'ipbh,ipheff,iptwiss,ipdampx,ipdampy,ipdampz,'//
     $ 'ipdnux,ipdnuy,ipdnuz,ipjx,ipjy,ipjz,'//
     $ 'ipemx,ipemy,ipemz,ipsige,ipsigz,iptaup,ipnup,'//
     $ 'ippolx,ippoly,ippolz,'//
     $ 'ipequpol,ipequpol2,ipequpol4,ipequpol6,'//
     $ 'ipnup,iprevf'//
     $ '},Constant];'//
     $ 'End[];EndPackage[];',kx,irtc)
c      call tfdebugprint(kx,'initemip',1)
      initemip=.false.
      return
      end subroutine

      real*8 function tfetwiss(r,cod,normi) result(twiss)
      use ffs
      implicit none
      real*8 , intent(in):: r(6,6),cod(6)
      dimension twiss(ntwissfun)
      real*8 hi(6,6)
      real*8 ax,ay,az,axy,f,detm,his(4),
     $     uz11,uz12,uz21,uz22,
     $     hx11,hx12,hx21,hx22,
     $     hy11,hy12,hy21,hy22,
     $     r11,r12,r21,r22,
     $     crx,cry,crz,cx,cy,cz,sx,sy,sz,
     $     bx21,bx22,by21,by22,bz21,bz22
      logical*4 normal
      logical*4 , intent(in) :: normi
c      write(*,'(a,1p5g15.7)')'tfetwiss ',r(5,5),r(5,6),r(6,5),r(6,6),
c     $     r(5,5)*r(6,6)-r(6,5)*r(5,6)
      az=sqrt(r(5,5)*r(6,6)-r(6,5)*r(5,6))
      uz11=r(5,5)/az
      uz12=r(5,6)/az
      uz21=r(6,5)/az
      uz22=r(6,6)/az
      hx11= r(6,2)*uz11-r(5,2)*uz21
      hx12= r(6,2)*uz12-r(5,2)*uz22
      hx21=-r(6,1)*uz11+r(5,1)*uz21
      hx22=-r(6,1)*uz12+r(5,1)*uz22
      hy11= r(6,4)*uz11-r(5,4)*uz21
      hy12= r(6,4)*uz12-r(5,4)*uz22
      hy21=-r(6,3)*uz11+r(5,3)*uz21
      hy22=-r(6,3)*uz12+r(5,3)*uz22
      f=1.d0/(1.d0+az)
      ax=1.d0-(hx11*hx22-hx21*hx12)*f
      ay=1.d0-(hy11*hy22-hy21*hy12)*f
      hi(2,2)=ax
      hi(1,2)=0.d0
      hi(4,2)= (hx12*hy21-hx11*hy22)*f
      hi(3,2)=-(-hx12*hy11+hx11*hy12)*f
      hi(6,2)=-hx11
      hi(5,2)= hx12
      hi(2,1)= 0.d0
      hi(1,1)= ax
      hi(4,1)=-(hx22*hy21-hx21*hy22)*f
      hi(3,1)=(-hx22*hy11+hx21*hy12)*f
      hi(6,1)= hx21
      hi(5,1)=-hx22
      hi(2,4)= hi(3,1)
      hi(1,4)=-hi(3,2)
      hi(4,4)= ay
      hi(3,4)= 0.d0
      hi(6,4)=-hy11
      hi(5,4)= hy12
      hi(2,3)=-hi(4,1)
      hi(1,3)= hi(4,2)
      hi(4,3)= 0.d0
      hi(3,3)= ay
      hi(6,3)= hy21
      hi(5,3)=-hy22
      hi(2,6)= hx22
      hi(1,6)= hx21
      hi(4,6)= hy22
      hi(3,6)= hy21
      hi(6,6)= az
      hi(5,6)= 0.d0
      hi(2,5)= hx12
      hi(1,5)= hx11
      hi(4,5)= hy12
      hi(3,5)= hy11
      hi(6,5)= 0.d0
      hi(5,5)= az
      hi=matmul(r,hi)
c      call tmultr(hi,r,6)
c      write(*,*)'tfetwiss-3'
      detm=(hi(1,1)*hi(2,2)-hi(2,1)*hi(1,2)
     $     +hi(3,3)*hi(4,4)-hi(4,3)*hi(3,4))*.5d0
c      write(*,'(a,1p6g15.7)')'tfetwiss-4 ',detm,ax,ay,az,f,xyth
      normal=detm > xyth
      if(.not. normal)then
        his=hi(1,1:4)
        hi(1,1:4)=hi(3,1:4)
        hi(3,1:4)=his
        his=hi(2,1:4)
        hi(2,1:4)=hi(4,1:4)
        hi(4,1:4)=his
        detm=(hi(1,1)*hi(2,2)-hi(2,1)*hi(1,2)
     $       +hi(3,3)*hi(4,4)-hi(4,3)*hi(3,4))*.5d0
      endif
      axy=sqrt(detm)
      r11=( hi(4,4)*hi(3,1)-hi(3,4)*hi(4,1))/axy
      r12=( hi(4,4)*hi(3,2)-hi(3,4)*hi(4,2))/axy
      r21=(-hi(4,3)*hi(3,1)+hi(3,3)*hi(4,1))/axy
      r22=(-hi(4,3)*hi(3,2)+hi(3,3)*hi(4,2))/axy
      crx=hypot(hi(1,2),hi(2,2))
c      crx=sqrt(hi(1,2)**2+hi(2,2)**2)
      cx= hi(2,2)/crx
      sx=-hi(1,2)/crx
      bx21=(-sx*hi(1,1)+cx*hi(2,1))/axy
      bx22=crx/axy
      cry=hypot(hi(3,4),hi(4,4))
c      cry=sqrt(hi(3,4)**2+hi(4,4)**2)
      cy= hi(4,4)/cry
      sy=-hi(3,4)/cry
      by21=(-sy*hi(3,3)+cy*hi(4,3))/axy
      by22=cry/axy
      crz=hypot(uz12,uz22)
c      crz=sqrt(uz12**2+uz22**2)
      cz= uz22/crz
      sz=-uz12/crz
      bz21=-sz*uz11+cz*uz21
      bz22=crz
      twiss(mfitr1)=r11
      twiss(mfitr2)=r12
      twiss(mfitr3)=r21
      twiss(mfitr4)=r22
      if(.not. normi)then
        normal=.not. normal
      endif
      if(normal)then
        twiss(mfitax)= bx21*bx22
        twiss(mfitbx)= bx22**2
        twiss(mfitnx)= atan2(sx,cx)
        twiss(mfitay)= by21*by22
        twiss(mfitby)= by22**2
        twiss(mfitny)= atan2(sy,cy)
        twiss(mfitex)= axy*hx12-r22*hy12+r12*hy22
        twiss(mfitepx)= axy*hx22+r21*hy12-r11*hy22
        twiss(mfitey) = axy*hy12+r11*hx12+r12*hx22
        twiss(mfitepy)=axy*hy22+r21*hx12+r22*hx22
        twiss(mfitdetr)=r11*r22-r12*r21
        twiss(mfitzx) =axy*hx11-r22*hy11+r12*hy21
        twiss(mfitzpx)=axy*hx21+r21*hy11-r11*hy21
        twiss(mfitzy) =axy*hy11+r11*hx11+r12*hx21
        twiss(mfitzpy)=axy*hy21+r21*hx11+r22*hx21
      else
        twiss(mfitay)= bx21*bx22
        twiss(mfitby)= bx22**2
        twiss(mfitny)= atan2(sx,cx)
        twiss(mfitax)= by21*by22
        twiss(mfitbx)= by22**2
        twiss(mfitnx)= atan2(sy,cy)
        twiss(mfitey)= axy*hx12-r22*hy12+r12*hy22
        twiss(mfitepy)= axy*hx22+r21*hy12-r11*hy22
        twiss(mfitex) = axy*hy12+r11*hx12+r12*hx22
        twiss(mfitepx)=axy*hy22+r21*hx12+r22*hx22
        twiss(mfitdetr)=1.d0+xyth-r11*r22+r12*r21
        twiss(mfitzy) =axy*hx11-r22*hy11+r12*hy21
        twiss(mfitzpy)=axy*hx21+r21*hy11-r11*hy21
        twiss(mfitzx) =axy*hy11+r11*hx11+r12*hx21
        twiss(mfitzpx)=axy*hy21+r21*hx11+r22*hx21
      endif
      twiss(mfitdx:mfitddp)=cod
      twiss(mfitaz)=bz21*bz22
      twiss(mfitbz)=bz22**2
      twiss(mfitnz)=atan2(sz,cz)
      return
      end function

      real*8 function etwiss2ri(twiss1,normal) result(ria)
      use ffs
      implicit none
      real*8 , intent(in)::twiss1(ntwissfun)
      dimension ria(6,6)
      real*8 h(4,6),br(4,4),
     $     hx11,hx12,hx21,hx22,
     $     hy11,hy12,hy21,hy22,
     $     ex,epx,ey,epy,zx,zpx,zy,zpy,
     $     r1,r2,r3,r4,amu,detr,sqrbx,sqrby,sqrbz,aa,
     $     ax,ay,az,dethx,dethy,amux,amuy,azz
      logical*4 ,intent(out):: normal
      r1=twiss1(mfitr1)
      r2=twiss1(mfitr2)
      r3=twiss1(mfitr3)
      r4=twiss1(mfitr4)
      detr=twiss1(mfitdetr)
      normal=detr < 1.d0
      if(normal)then
        amu=sqrt(1.d0-detr)
        ex =twiss1(mfitex)
        epx=twiss1(mfitepx)
        ey =twiss1(mfitey)
        epy=twiss1(mfitepy)
        zx =twiss1(mfitzx)
        zpx=twiss1(mfitzpx)
        zy =twiss1(mfitzy)
        zpy=twiss1(mfitzpy)
        sqrbx=sqrt(twiss1(mfitbx))
        ax=twiss1(mfitax)
        sqrby=sqrt(twiss1(mfitby))
        ay=twiss1(mfitay)
      else
        amu=sqrt(1.d0-r1*r4+r2*r3)
        ey =twiss1(mfitex)
        epy=twiss1(mfitepx)
        ex =twiss1(mfitey)
        epx=twiss1(mfitepy)
        zy =twiss1(mfitzx)
        zpy=twiss1(mfitzpx)
        zx =twiss1(mfitzy)
        zpx=twiss1(mfitzpy)
        sqrby=sqrt(twiss1(mfitbx))
        ay=twiss1(mfitax)
        sqrbx=sqrt(twiss1(mfitby))
        ax=twiss1(mfitay)
      endif
      sqrbz=sqrt(twiss1(mfitbz))
      az=twiss1(mfitaz)
      hx11=amu*zx -r2*zpy+r4*zy
      hx12=amu*ex -r2*epy+r4*ey
      hx21=amu*zpx+r1*zpy-r3*zy
      hx22=amu*epx+r1*epy-r3*ey
      hy11=amu*zy -r2*zpx-r1*zx
      hy12=amu*ey -r2*epx-r1*ex
      hy21=amu*zpy-r4*zpx-r3*zx
      hy22=amu*epy-r4*epx-r3*ex
      dethx=hx11*hx22-hx21*hx12
      dethy=hy11*hy22-hy21*hy12
      azz=sqrt(1.d0-dethx-dethy)
      aa=1.d0/(1.d0+azz)
      amux=1.d0-dethx*aa
      amuy=1.d0-dethy*aa
      h(1,1)=amux
      h(1,2)=0.d0
      h(1,3)=( hx12*hy21-hx11*hy22)*aa
      h(1,4)=(-hx12*hy11+hx11*hy12)*aa
      h(1,5)=-hx11
      h(1,6)=-hx12
      h(2,1)=0.d0
      h(2,2)=amux
      h(2,3)=( hx22*hy21-hx21*hy22)*aa
      h(2,4)=(-hx22*hy11+hx21*hy12)*aa
      h(2,5)=-hx21
      h(2,6)=-hx22
      h(3,1)= h(2,4)
      h(3,2)=-h(1,4)
      h(3,3)=amuy
      h(3,4)=0.d0
      h(3,5)=-hy11
      h(3,6)=-hy12
      h(4,1)=-h(2,3)
      h(4,2)= h(1,3)
      h(4,3)=0.d0
      h(4,4)=amuy
      h(4,5)=-hy21
      h(4,6)=-hy22
      ria(5,1)= hx22
      ria(5,2)=-hx12
      ria(5,3)= hy22
      ria(5,4)=-hy12
      ria(5,5)=azz
      ria(5,6)=0.d0
      ria(6,1)=-hx21
      ria(6,2)= hx11
      ria(6,3)=-hy21
      ria(6,4)= hy11
      ria(6,5)=0.d0
      ria(6,6)=azz
      br(1,1)= amu/sqrbx
      br(1,2)=0.d0
      br(1,3)=-r4/sqrbx
      br(1,4)= r2/sqrbx
      br(2,1)= amu*ax/sqrbx
      br(2,2)= amu*sqrbx
      br(2,3)= r3*sqrbx-r4*ax/sqrbx
      br(2,4)= r2*ax/sqrbx-r1*sqrbx
      br(3,1)= r1/sqrby
      br(3,2)= r2/sqrby
      br(3,3)= amu/sqrby
      br(3,4)=0.d0
      br(4,1)= r1*ay/sqrby+r3*sqrby
      br(4,2)= r2*ay/sqrby+r4*sqrby
      br(4,3)= amu*ay/sqrby
      br(4,4)= amu*sqrby
      ria(1,1)=br(1,1)*h(1,1)+br(1,2)*h(2,1)
     $       +br(1,3)*h(3,1)+br(1,4)*h(4,1)
      ria(1,2)=br(1,1)*h(1,2)+br(1,2)*h(2,2)
     $       +br(1,3)*h(3,2)+br(1,4)*h(4,2)
      ria(1,3)=br(1,1)*h(1,3)+br(1,2)*h(2,3)
     $       +br(1,3)*h(3,3)+br(1,4)*h(4,3)
      ria(1,4)=br(1,1)*h(1,4)+br(1,2)*h(2,4)
     $       +br(1,3)*h(3,4)+br(1,4)*h(4,4)
      ria(1,5)=br(1,1)*h(1,5)+br(1,2)*h(2,5)
     $       +br(1,3)*h(3,5)+br(1,4)*h(4,5)
      ria(1,6)=br(1,1)*h(1,6)+br(1,2)*h(2,6)
     $       +br(1,3)*h(3,6)+br(1,4)*h(4,6)
      ria(2,1)=br(2,1)*h(1,1)+br(2,2)*h(2,1)
     $       +br(2,3)*h(3,1)+br(2,4)*h(4,1)
      ria(2,2)=br(2,1)*h(1,2)+br(2,2)*h(2,2)
     $       +br(2,3)*h(3,2)+br(2,4)*h(4,2)
      ria(2,3)=br(2,1)*h(1,3)+br(2,2)*h(2,3)
     $       +br(2,3)*h(3,3)+br(2,4)*h(4,3)
      ria(2,4)=br(2,1)*h(1,4)+br(2,2)*h(2,4)
     $       +br(2,3)*h(3,4)+br(2,4)*h(4,4)
      ria(2,5)=br(2,1)*h(1,5)+br(2,2)*h(2,5)
     $       +br(2,3)*h(3,5)+br(2,4)*h(4,5)
      ria(2,6)=br(2,1)*h(1,6)+br(2,2)*h(2,6)
     $       +br(2,3)*h(3,6)+br(2,4)*h(4,6)
      ria(3,1)=br(3,1)*h(1,1)+br(3,2)*h(2,1)
     $       +br(3,3)*h(3,1)+br(3,4)*h(4,1)
      ria(3,2)=br(3,1)*h(1,2)+br(3,2)*h(2,2)
     $       +br(3,3)*h(3,2)+br(3,4)*h(4,2)
      ria(3,3)=br(3,1)*h(1,3)+br(3,2)*h(2,3)
     $       +br(3,3)*h(3,3)+br(3,4)*h(4,3)
      ria(3,4)=br(3,1)*h(1,4)+br(3,2)*h(2,4)
     $       +br(3,3)*h(3,4)+br(3,4)*h(4,4)
      ria(3,5)=br(3,1)*h(1,5)+br(3,2)*h(2,5)
     $       +br(3,3)*h(3,5)+br(3,4)*h(4,5)
      ria(3,6)=br(3,1)*h(1,6)+br(3,2)*h(2,6)
     $       +br(3,3)*h(3,6)+br(3,4)*h(4,6)
      ria(4,1)=br(4,1)*h(1,1)+br(4,2)*h(2,1)
     $       +br(4,3)*h(3,1)+br(4,4)*h(4,1)
      ria(4,2)=br(4,1)*h(1,2)+br(4,2)*h(2,2)
     $       +br(4,3)*h(3,2)+br(4,4)*h(4,2)
      ria(4,3)=br(4,1)*h(1,3)+br(4,2)*h(2,3)
     $       +br(4,3)*h(3,3)+br(4,4)*h(4,3)
      ria(4,4)=br(4,1)*h(1,4)+br(4,2)*h(2,4)
     $       +br(4,3)*h(3,4)+br(4,4)*h(4,4)
      ria(4,5)=br(4,1)*h(1,5)+br(4,2)*h(2,5)
     $       +br(4,3)*h(3,5)+br(4,4)*h(4,5)
      ria(4,6)=br(4,1)*h(1,6)+br(4,2)*h(2,6)
     $       +br(4,3)*h(3,6)+br(4,4)*h(4,6)
      ria(5,1:6)=ria(5,1:6)/sqrbz
      ria(6,1:6)=ria(6,1:6)*sqrbz+ria(5,1:6)*az
      return
      end function

      subroutine tfnormalcoord(isp1,kx,irtc)
      use tfstk
      use ffs
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      logical*4 normal
      type (sad_rlist), pointer :: kl
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. tfnumlistqn(dtastk(isp),ntwissfun,kl))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real List of Length 28"')
        return
      endif
      kx=kxm2l(etwiss2ri(kl%rbody(1:ntwissfun),normal),6,6,6,.false.)
      irtc=0
      return
      end subroutine

      subroutine tmulbs(beam,trans1,beamr)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      real*8 ,intent(inout):: beam(:)
      real*8 ,intent(in):: trans1(6,6)
c      real*8 ,intent(in),target:: trans(:,:)
c      real*8 , pointer:: trans1(:,:)
c      real*8 ,target::transs(6,6)
      real*8 s(6,6)
      logical*4 ,intent(in):: beamr
      if(irad < 12)then
        return
      endif
        s(1,1:6)=trans1(1:6,1)*beam(1) +trans1(1:6,2)*beam(2)
     1          +trans1(1:6,3)*beam(4) +trans1(1:6,4)*beam(7)
     1          +trans1(1:6,5)*beam(11)+trans1(1:6,6)*beam(16)
        s(2,1:6)=trans1(1:6,1)*beam(2) +trans1(1:6,2)*beam(3)
     1          +trans1(1:6,3)*beam(5) +trans1(1:6,4)*beam(8)
     1          +trans1(1:6,5)*beam(12)+trans1(1:6,6)*beam(17)
        s(3,1:6)=trans1(1:6,1)*beam(4) +trans1(1:6,2)*beam(5)
     1          +trans1(1:6,3)*beam(6) +trans1(1:6,4)*beam(9)
     1          +trans1(1:6,5)*beam(13)+trans1(1:6,6)*beam(18)
        s(4,1:6)=trans1(1:6,1)*beam(7) +trans1(1:6,2)*beam(8)
     1          +trans1(1:6,3)*beam(9) +trans1(1:6,4)*beam(10)
     1          +trans1(1:6,5)*beam(14)+trans1(1:6,6)*beam(19)
        s(5,1:6)=trans1(1:6,1)*beam(11)+trans1(1:6,2)*beam(12)
     1          +trans1(1:6,3)*beam(13)+trans1(1:6,4)*beam(14)
     1          +trans1(1:6,5)*beam(15)+trans1(1:6,6)*beam(20)
        s(6,1:6)=trans1(1:6,1)*beam(16)+trans1(1:6,2)*beam(17)
     1          +trans1(1:6,3)*beam(18)+trans1(1:6,4)*beam(19)
     1          +trans1(1:6,5)*beam(20)+trans1(1:6,6)*beam(21)
c      enddo
        beam(1) =dot_product(s(:,1),trans1(1,:))
        beam(2) =dot_product(s(:,1),trans1(2,:))
        beam(3) =dot_product(s(:,2),trans1(2,:))
        beam(4) =dot_product(s(:,1),trans1(3,:))
        beam(5) =dot_product(s(:,2),trans1(3,:))
        beam(6) =dot_product(s(:,3),trans1(3,:))
        beam(7) =dot_product(s(:,1),trans1(4,:))
        beam(8) =dot_product(s(:,2),trans1(4,:))
        beam(9) =dot_product(s(:,3),trans1(4,:))
        beam(10)=dot_product(s(:,4),trans1(4,:))
        beam(11)=dot_product(s(:,1),trans1(5,:))
        beam(12)=dot_product(s(:,2),trans1(5,:))
        beam(13)=dot_product(s(:,3),trans1(5,:))
        beam(14)=dot_product(s(:,4),trans1(5,:))
        beam(15)=dot_product(s(:,5),trans1(5,:))
        beam(16)=dot_product(s(:,1),trans1(6,:))
        beam(17)=dot_product(s(:,2),trans1(6,:))
        beam(18)=dot_product(s(:,3),trans1(6,:))
        beam(19)=dot_product(s(:,4),trans1(6,:))
        beam(20)=dot_product(s(:,5),trans1(6,:))
        beam(21)=dot_product(s(:,6),trans1(6,:))
       if(calint .and. beamr)then
c         do j=1,6
           s(1,1:6)=trans1(1:6,1)*beam(21+1) +trans1(1:6,2)*beam(21+2)
     1          +trans1(1:6,3)*beam(21+4) +trans1(1:6,4)*beam(21+7)
     1          +trans1(1:6,5)*beam(21+11)+trans1(1:6,6)*beam(21+16)
           s(2,1:6)=trans1(1:6,1)*beam(21+2) +trans1(1:6,2)*beam(21+3)
     1          +trans1(1:6,3)*beam(21+5) +trans1(1:6,4)*beam(21+8)
     1          +trans1(1:6,5)*beam(21+12)+trans1(1:6,6)*beam(21+17)
           s(3,1:6)=trans1(1:6,1)*beam(21+4) +trans1(1:6,2)*beam(21+5)
     1          +trans1(1:6,3)*beam(21+6) +trans1(1:6,4)*beam(21+9)
     1          +trans1(1:6,5)*beam(21+13)+trans1(1:6,6)*beam(21+18)
           s(4,1:6)=trans1(1:6,1)*beam(21+7) +trans1(1:6,2)*beam(21+8)
     1          +trans1(1:6,3)*beam(21+9) +trans1(1:6,4)*beam(21+10)
     1          +trans1(1:6,5)*beam(21+14)+trans1(1:6,6)*beam(21+19)
           s(5,1:6)=trans1(1:6,1)*beam(21+11)+trans1(1:6,2)*beam(21+12)
     1          +trans1(1:6,3)*beam(21+13)+trans1(1:6,4)*beam(21+14)
     1          +trans1(1:6,5)*beam(21+15)+trans1(1:6,6)*beam(21+20)
           s(6,1:6)=trans1(1:6,1)*beam(21+16)+trans1(1:6,2)*beam(21+17)
     1          +trans1(1:6,3)*beam(21+18)+trans1(1:6,4)*beam(21+19)
     1          +trans1(1:6,5)*beam(21+20)+trans1(1:6,6)*beam(21+21)
c         enddo
        beam(21+1) =dot_product(s(:,1),trans1(1,:))
        beam(21+2) =dot_product(s(:,1),trans1(2,:))
        beam(21+3) =dot_product(s(:,2),trans1(2,:))
        beam(21+4) =dot_product(s(:,1),trans1(3,:))
        beam(21+5) =dot_product(s(:,2),trans1(3,:))
        beam(21+6) =dot_product(s(:,3),trans1(3,:))
        beam(21+7) =dot_product(s(:,1),trans1(4,:))
        beam(21+8) =dot_product(s(:,2),trans1(4,:))
        beam(21+9) =dot_product(s(:,3),trans1(4,:))
        beam(21+10)=dot_product(s(:,4),trans1(4,:))
        beam(21+11)=dot_product(s(:,1),trans1(5,:))
        beam(21+12)=dot_product(s(:,2),trans1(5,:))
        beam(21+13)=dot_product(s(:,3),trans1(5,:))
        beam(21+14)=dot_product(s(:,4),trans1(5,:))
        beam(21+15)=dot_product(s(:,5),trans1(5,:))
        beam(21+16)=dot_product(s(:,1),trans1(6,:))
        beam(21+17)=dot_product(s(:,2),trans1(6,:))
        beam(21+18)=dot_product(s(:,3),trans1(6,:))
        beam(21+19)=dot_product(s(:,4),trans1(6,:))
        beam(21+20)=dot_product(s(:,5),trans1(6,:))
        beam(21+21)=dot_product(s(:,6),trans1(6,:))
       endif
       return
       end subroutine

      real*8 function tfinibeam(is) result(beam)
      use sad_main, only:iaidx
      use ffs
      use ffs_pointer,only:twiss
      use tffitcode,only:ntwissfun
      use ffs_flag, only:wspac,intra
      use sad_basics
      implicit none
      dimension beam(21)
      integer*4 ,intent(in):: is
      integer*4 ir
      real*8 ris(6,6),rgetgl1,emmin,emxm,emym
      logical*4 normal
      ris=etwiss2ri(twiss(is,0,1:ntwissfun),normal)
      beam=0.d0
      if(wspac .or. INTRA)then
        emmin=(emx+emy)*rgetgl1('MINCOUP')
        emxm=max(emmin,emx)
        emym=max(emmin,emy)
      else
        emxm=emx
        emym=emy
      endif
      beam(iaidx(1,1))=emxm
      beam(iaidx(2,2))=emxm
      beam(iaidx(3,3))=emym
      beam(iaidx(4,4))=emym
      beam(iaidx(5,5))=emz
      beam(iaidx(6,6))=emz
      ir=irad
      irad=12
      call tmulbs(beam,tinv6(ris),.false.)
      irad=ir
      if(is == 1)then
        ri=ris
        r=tinv6(ri)
      endif
      return
      end function

      end module temw

      subroutine tffsinitparam
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_descriptor) kdp
      type (sad_symdef), pointer :: symddp
      integer*8 j
      real*8 rgetgl1
      j=idvalc(1)
      emx=rlist(j+ky_EMIX_MARK)
      if(emx <= 0.d0)then
        emx=rgetgl1('EMITX')
      else
        call rsetgl1('EMITX',emx)
      endif
      emy=rlist(j+ky_EMIY_MARK)
      if(emy <= 0.d0)then
        emy=rgetgl1('EMITY')
      else
        call rsetgl1('EMITY',emy)
      endif
      emz=rlist(j+ky_EMIZ_MARK)
      if(emz <= 0.d0)then
        emz=rgetgl1('EMITZ')
      else
        call rsetgl1('EMITZ',emz)
      endif
      sizedp=rlist(j+ky_SIGE_MARK)
      if(sizedp <= 0.d0)then
        sizedp=rgetgl1('SIGE')
      else
        call rsetgl1('SIGE',sizedp)
      endif
      sigzs=max(0.d0,rlist(j+ky_SIGZ_MARk))
      if(sigzs <= 0.d0)then
        sigzs=rgetgl1('SIGZ')
      else
        call rsetgl1('SIGZ',sigzs)
      endif
      dpmax=max(0.d0,rlist(j+ky_DP_MARk))
      kdp=kxsymbolz('DP',2,symddp)
      if(dpmax <= 1.d-30)then
        dpmax=rfromd(symddp%value)
      endif
      if(dpmax <= 1.d-30)then
        dpmax=0.01d0
      endif
      symddp%value=dfromr(dpmax)
      if(rlist(latt(1)+ky_BX_MARK) <= 0.d0)then
        rlist(latt(1)+ky_BX_MARK)=1.d0
      endif
      if(rlist(latt(1)+ky_BY_MARK) <= 0.d0)then
        rlist(latt(1)+ky_BY_MARK)=1.d0
      endif
      if(rlist(latt(1)+ky_BZ_MARK) <= 0.d0)then
        rlist(latt(1)+ky_BZ_MARK)=1.d0
      endif
      return
      end

      subroutine tffs_init_elmv
      use tfstk
      use ffs
      use ffs_pointer
      implicit none
      return
      end

      subroutine tffsalloc()
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_main
c      use ffsa
      implicit none
      integer*4 ntwis
      marki=1
      nlat=elatt%nlat0+1
      latt(1:nlat)=>elatt%comp(1:nlat)
      if(idtypec(1) /= icMARK)then
        write(*,*)'The first element must be a MARK element.'
        call abort
      endif
      call tfhashelement
      call tffs_init_nele
      ifibzl =ktaloc(nlat*3/2+2)
      ifcoup=ktaloc(nlat)
      iferrk=ktaloc(nlat*2)
      ifmast =ktaloc(nlat/2+1)
      ifcomp =ktaloc(nlat/2+1)
      ifele2=ktaloc(nlat)
      ifiprev =ktaloc(nlat/2+1)
      ifinext =ktaloc(nlat/2+1)
      nve=0
      call tffsnvealloc(nele)
      ndim=1
      ndima=ndim*2+1
      ntwis =nlat*ndima
      iftwis=ktaloc(ntwis*ntwissfun)
      ifpos =ktaloc(nlat)
      ifgeo =ktaloc(nlat*12)
      ifgamm=ktaloc(nlat)
      iftwissp=ktaloc(nlat/2+1)
      ifnelv=ktaloc(nele*lnelv)
      call ffs_init_pointer
      nelvx(1:nele)%dcomp%k=i00
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

      subroutine tffs_init_nele
      use tfstk
      use ffs
      use ffs_pointer
      implicit none
      integer*8 j
      integer*4 i,l,k,itehash
      ifele1=ktaloc(nlat/2+1)
      ifmult=ktaloc(nlat/2+1)
      nele=0
      LOOP_L: do l=1,nlat-1
        j=ielmhash+itehash(pnamec(l),MAXPNAME)*2
        LOOP_K: do k=1,ilist(1,j+1)
          i=ilist(1,k+klist(j+2)-1)
          if(i .ge. l)then
            exit LOOP_K
          endif
          if(idelc(i) == idelc(l))then
            ilist(l,ifele1)=ilist(i,ifele1)
            cycle LOOP_L
          endif
        enddo LOOP_K
        nele=nele+1
        ilist(nele,ifmult)=l
        ilist(l,ifele1)=nele
      enddo LOOP_L
c      if(elatt%elmv%k == 0)then
c        call tffs_init_elmv
c      endif
      return
      end subroutine

      subroutine tffsnvealloc(nv)
      use tfstk
      use ffs
      use ffs_pointer
      use iso_c_binding
      implicit none
      integer*4 , intent(in)::nv
      integer*4 nve0
      integer*8 ifnvev1
      if(nv < nve)then
        return
      endif
      nve0=nve
      nve=max(nv,nve0+128)
      ifnvev1=ktaloc(nve*lnvev)
      if(nve0 /= 0)then
        klist(ifnvev1:ifnvev1+nve0*lnvev-1)=
     $       klist(ifnvev:ifnvev+nve0*lnvev-1)
        call tfree(ifnvev)
      endif
      ifnvev=ifnvev1
      call c_f_pointer(c_loc(klist(ifnvev)),nvevx,[nve])
      return
      end subroutine

      subroutine tffsfree
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 l,i
      levele=levele+1
      call tfresethash
c      call tfree(ilist(2,ifwakep))
c      call tfree(ilist(2,ifwakep+2))
c      call tfree(ifwakep)
      call tfree(ifnvev)
      call tfree(iftwissp)
      if(ifsize > 0)then
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
      call tfree(ifinext)
      call tfree(ifiprev)
      call tfree(ifele2)
      call tfree(ifele1)
      call tfree(ifcomp)
      do i=1,nele
        call tflocald(nelvx(i)%dcomp)
      enddo
      call tfree(ifnelv)
      call tfree(iferrk)
      call tfree(ifibzl)
      call tfree(ifgamm)
      l=itfdownlevel()
      return
      end

      subroutine tfsetconvcase(f)
      use ffs_flag
      implicit none
      logical*4 f
      convcase=f
      return
      end

      subroutine tfmain(isp1,kx,irtc)
      use tfstk
      use trackbypass, only: bypasstrack
      use tfrbuf
      use tmacro
      use ffsfile
      use tfcsi,only:lfni
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,infl0,itfmessage,lfn,ierrfl,ierr,lfnb0,lfni0
      ierr=0
      if(isp == isp1+2)then
        ierr=int(rtastk(isp))
      elseif(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1 or 2"')
        return
      endif
      call tfopenread(isp1,kx,irtc)
      if(irtc /= 0)then
        return
      endif
      if(.not. ktfrealq(kx,lfn))then
        irtc=itfmessage(9,'General::wonrgtype','"Real"')
        return
      endif
      lfni0=lfni
      infl0=infl
      infl=lfn
      if(infl /= lfni)then
        call trbassign(infl)
      endif
      ierrfl=errfl
      errfl=ierr
      bypasstrack=.true.
      lfnb0=lfnbase
      lfnbase=lfnp
      call toplvl
      lfnp=lfnbase
      lfnbase=lfnb0
      bypasstrack=.false.
      errfl=ierrfl
      call trbclose(lfn)
      infl=lfni0
      if(lfni /= lfni0)then
        lfni=lfni0
        call trbassign(lfni)
      endif
      kx%k=ktfoper+mtfnull
      return
      end
