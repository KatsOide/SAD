      subroutine tffsmatch(flv,twiss,pos,geo,gammab,utwiss,
     $     latt,mult,iele,iele1,iele2,couple,itwissp,
     $     ival,klp,ivarele,ivvar,ivcomp,valvar,
     $     errk,vlim,nlat,nele,ndim,nut,nvar,
     $     dp,tracex,tracey,hstab,vstab,df,nfr,nfam,nfam1,
     $     dfam,jfam,kfam,inicond,iuid,uini,
     $     wake,iwakeelm,kwaketbl,nwakep,
     $     nqcol,nqcol1,nfcol,nfc0,maxcond,
     $     nlist,brho,
     $     emx,emy,dpmax,dp0,coumin,r,residual,absweit,
     $     cell,fitflg,geomet,cellstab,convgo,nparallel,
     $     orbitcal,lfno,irtc)
      use tfstk
      use ffslocal, only:ffslocalv
      use tffitcode
      use tfshare
      implicit none
      include 'inc/MACCODE.inc'
      type (ffslocalv) flv
c      include 'DEBUG.inc'
      integer*8 ifqu,ifqu0,itmmapp,iuta1,kqu
      real*8 flim1,flim2,aimp1,aimp2,badc1,badc2,amedc1,amedc2,alit,
     $     wlmin,eps,eps1,brho
      parameter (flim1=-4.d0,flim2=-3.d0,aimp1=-1.8d0,aimp2=-.8d0,
     $     badc1=-3.5d0,badc2=-2.5d0,amedc1=-2.3d0,amedc2=-1.3d0,
     $     alit=0.75d0,wlmin=0.009d0,eps=1.d-5,eps1=1.d-8)
      integer*4 nlat,nele,ndim,nfr,maxcond,ibegin,nqcol,nut,
     $     lfno,irtc,nqumax,nfcol,nfam,nfam1,nvar,
     $     nqcol0,nparallel,nqcol1,nqcol00,iuid(-nfam:nfam),
     $     nqcola1,nqcol1a1,nqcola2,nqcol1a2,lout,
     $     nwakep,iwakeelm(nwakep)
      integer*8 kwaketbl(2,nwakep)
      integer*4 latt(2,nlat),iele(nlat),iele1(nlat),iele2(nlat),
     $     mult(nlat),itwissp(nut),ivcomp(nvar),
     $     ival(nele),ivarele(nvar),ivvar(nvar),klp(nele),
     $     kdpa1(maxcond),kdpa2(maxcond),
     $     iqcol0(maxcond),jfam(-nfam:nfam),
     $     iqcola1(maxcond),iqcola2(maxcond),
     $     kfam(-nfam:nfam),nfc0,nretry,nstab,nstaba1
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),pos(nlat),
     $     geo(3,4,nlat),
     $     gammab(nlat),couple(nlat),utwiss(ntwissfun,-nfam:nfam,nut),
     $     vlim(nele,2),dval(nvar),errk(2,nlat),
     $     valvar(nele*2),dfam(4,-nfam:nfam),
     $     dp(-nfam:nfam),tracex(-nfam:nfam),tracey(-nfam:nfam),
     $     df(maxcond),df0(maxcond),emx,emy,wfit(maxcond),
     $     wiq(maxcond),df1(maxcond),df2(maxcond),
     $     ddf1(maxcond),ddf2(maxcond),r,wexponent,dpmax,dp0,
     $     residual(-nfam:nfam),residuala1(-nfam:nfam),
     $     offmw,coumin,wsum,v00,uini(6,-nfam:nfam)
      logical*4 free(nele),cell,zcal,fitflg,geomet,cellstab,wcal,
     $     hstab(-nfam:nfam),vstab(-nfam:nfam),error,error2,
     $     absweit,inicond,wake,limited1,chgini,wcal1,zcal1
      character*8 nlist(mfit1)
      integer*4 iter,ii,j,kc,nvara, i,
     $     lfpa1(2,maxcond),lfpa2(2,maxcond),
     $     ip,fork_worker,ipr,istep,npr(nparallel)
      real*8 rl,fuzz,a,b,x,crate,aimprv,fact,r0,r00,ra,alate,
     $     smallf,badcnv,amedcv,tweigh,rstab,vl1,vl2,
     $     avebeta,etamax,aitm1,aitm2,bestval(nvar),wvar(nvar),
     $     dg,f1,f2,g1,g2,ra1,rpa1,rstaba1,valvar0,
     $     tffsfmin,rp,rp0,wlimit(nvar),dv,vl,dvkc
      real*8 twisss(ntwissfun)
      logical*4 chgmod,newton,convgo,imprv,limited,over,
     $     parallel,nderiv,outt,orbitcal,nderiv0,dlim
      integer*4 itfgetrecl,kkk,kkkk,npa
      character ch
      integer*8 intffs,inumderiv,iexponent,inumw,iconvergence
      data intffs,inumderiv,iexponent,inumw,iconvergence /0,0,0,0,0/
c
      fuzz(x,a,b)=max(0.d0,min(1.d0,(x-a)/(b-a)))
c     
c     begin initialize for preventing compiler warning
      if(inumderiv .eq. 0)then
        inumderiv   =ktfsymbolz('FFS$NumericalDerivative',23)-4
        intffs      =ktfsymbolz('FFS$Interrupt',13)-4
        iexponent   =ktfsymbolz('ExponentOfResidual',18)-4
        inumw       =ktfsymbolz('OffMomentumWeight',17)-4
        iconvergence=ktfsymbolz('CONVERGENCE',11)-4
      endif
      nderiv0=rlist(inumderiv) .ne. 0.d0
      nvara=0
      aitm1=0
      aitm2=0
      r0=0.d0
      r00=0.d0
      ra=0.d0
      aimprv=0.d0
      crate=1.d0
      nqcol0=0
      outt=.true.
      dlim=.false.
c     end   initialize for preventing compiler warning
      parallel=nparallel .gt. 1 .and. .not. inicond
      if(itfgetrecl() .lt. 120)then
        ch=char(10)
      else
        ch=' '
      endif
      irtc=0
      nqcol=0
      nqumax=0
      wcal=.true.
      wexponent=max(.1d0,min(4.d0,rlist(iexponent)))
      offmw=rlist(inumw)
      rl=abs(rlist(iconvergence))*max(nfcol,1)
      zcal=.true.
      if(cell)then
         call twmov(latt(2,1),twisss,1,0,.true.)
      endif
      ibegin=1
      rstab=0.d0
      if(fitflg)then
        call tffssetlimit(valvar,vlim,ivvar,ival,nvar,
     $     latt,ivarele,ivcomp,
     $     klp,couple,errk,iele,iele1,dlim)
        fact=1.d0
        iter=0
        newton=.true.
        chgmod=.true.
        bestval(1:nvar)=valvar(1:nvar)
        free=.false.
        do i=1,nvar
          free(ivarele(i))=.true.
        enddo
        aitm1=flv%itmax*alit
        aitm2=flv%itmax
        nvara=nvar
      endif
      chgini=.true.
      if(cell)then
        nretry=1
      else
        nretry=0
      endif
      lout=lfno
      do 9000: do while(.true.)
        do 200: do kkk=1,1
          call tftclupdate(int(rlist(intffs)))
          dp0=rlist(latt(2,1)+mfitddp)
          call tffscalc(flv,twiss,pos,geo,gammab,utwiss,
     $         latt,mult,ivarele,ivcomp,iele,iele1,ivvar,ival,
     $         nvar,klp,itwissp,nut,nlat,nele,ndim,
     $         dp,tracex,tracey,hstab,vstab,nfr,nfam,nfam1,
     $         flv%kdp,df,
     $         wfit,wiq,wsum,flv%iqcol,flv%lfp,nqcol,nqcol1,nfcol,nfc0,
     $         maxcond,nlist,
     $         ibegin,dfam,jfam,kfam,inicond,iuid,uini,
     $         wake,iwakeelm,kwaketbl,nwakep,
     $         r,rp,rstab,nstab,residual,wexponent,
     $         offmw,etamax,avebeta,emx,emy,dpmax,coumin,
     $         cell,zcal,fitflg,geomet,cellstab,wcal,parallel,
     $         chgini,orbitcal,
     $         lout,error)
          if(error)then
            if(irtc .eq. 20001)then
              exit do9000
            else
              fitflg=.false.
              irtc=20001
              bestval(1:nvar)=valvar(1:nvar)
              if(cell)then
                call twmov(latt(2,1),twisss,1,0,.false.)
              endif
              exit do200
            endif
          endif
          convgo=r .le. rl
          fitflg=fitflg .and. nqcol .gt. 0
          if(convgo)then
            write(lfno,9501)' Matched. (',r,')',
     $           dpmax,dp0,wexponent,ch,offmw
 9501       format(a,1pG11.4,a,
     $           ' DP =',0pf8.5,'  DP0 =',f8.5,
     $           '  ExponentOfResidual =',f4.1,a,
     $           ' OffMomentumWeight =',f8.3)
            fitflg=.false.
            if(.not. geomet)then
              call tfgeo(latt,geo,pos,gammab,.true.)
            endif
            exit do9000
          elseif(.not. fitflg)then
            write(lfno,9501)' Residual =',r,' ',
     $           dpmax,dp0,wexponent,ch,offmw
            exit do9000
          else
            if(chgini .and. cell)then
              call twmov(latt(2,1),twisss,1,0,.true.)
            endif
            chgini=.true.
            do1082: do kkkk=1,1
              iter=iter+1
              if(chgmod)then
                fact=min(1.d0,fact*2.d0)
                f1=0.d0
                g1=rp
                chgmod=.false.
                aimprv=0.d0
                crate=1.d0
                r0=r
                rp0=rp
                r00=r0*1.06d0
                ra=r0*1.000000001d0
                bestval(1:nvar)=valvar(1:nvar)
                if(cell)then
                  call twmov(latt(2,1),twisss,1,0,.true.)
                endif
              else
                imprv=r .lt. r0
                if(imprv)then
                  if(r .lt. r00*.95d0)then
                    lout=lfno
                    if(outt)then
                      write(lfno,*)
     $     'Iterations  Residual    Method     Reduction  Variables'
                      outt=.false.
                    endif
                    if(newton)then
                      write(lfno,9701)iter,r,'  (NEWTON)  ',fact,nvara
c     $                     2.d0*(rp-rp0)/dg/fact-1.d0
 9701                 format(i8,3x,1pG11.4,a,1pG11.4,i7)
                    else
                      write(lfno,9701)iter,r,'  (DESCEND) ',fact,nvara
                    endif
                    r00=r
                  endif
                  aimprv=fuzz(log10((ra-r)/ra),aimp1,aimp2)
                  crate=(r0-r)/r0
                  if(newton)then
                    fact=min(1.d0,fact*4.d0)
                  else
                    fact=min(1.d0,fact*2.d0)
                  endif
                  f1=0.d0
                  g1=rp
                  bestval(1:nvar)=valvar(1:nvar)
                  if(cell)then
                    call twmov(latt(2,1),twisss,1,0,.true.)
                  endif
                  rp0=rp
                  r0=r
                endif
                alate=fuzz(dble(iter),aitm1,aitm2)
                smallf=1.d0-fuzz(log10(fact),flim1,flim2)
                badcnv=1.d0-fuzz(log10(crate),badc1,badc2)
                amedcv=1.d0-fuzz(log10(crate),amedc1,amedc2)
                if(newton)then
                  chgmod=max(smallf,badcnv,min(alate,amedcv)) .gt. .5d0
                else
                  if(aimprv .gt. .5d0)then
                    chgmod=.true.
                  elseif(max(smallf,badcnv,min(alate,amedcv))
     $                   .gt. .5d0)then
                    chgmod=.true.
                    chgini=.true.
                    if(nretry .gt. 0)then
                      nretry=nretry-1
                    else
                      fitflg=.false.
                    endif
                  else
                    chgmod=.false.
                  endif
                endif
                if(chgmod)then
                  r=r0
                  newton=.not. newton
                  fact=1.d0
                  valvar(1:nvar)=bestval(1:nvar)
                  if(cell)then
                    call twmov(latt(2,1),twisss,1,0,.false.)
                  endif
                  exit do200
                elseif(.not. imprv)then
                  if(newton)then
                    f2=f1
                    g2=g1
                    f1=fact
                    g1=rp
                    fact=tffsfmin(f1,f2,g1,g2,rp0,dg)
                  else
                    f1=fact
                    fact=fact*.5d0
                  endif
                  if(cell)then
                    call twmov(latt(2,1),twisss,1,0,.false.)
                  endif
                  if(nvara .eq. nvar)then
                    a=fact/f1
                    valvar(1:nvar)=valvar(1:nvar)*a
     $                   +bestval(1:nvar)*(1.d0-a)
c                    do i=1,nvar
c                      valvar(i)=valvar(i)*a+bestval(i)*(1.d0-a)
c                    enddo
                    exit do200
                  else
                    nqcol=nqcol0
                    flv%iqcol(1:nqcol)=iqcol0(1:nqcol)
                    df(1:nqcol)=df0(1:nqcol)
                    nvara=nvar
                    wlimit(1:nvar)=max(wlmin,wlimit(1:nvar))
c                    do i=1,nvar
c                      wlimit(i)=max(wlmin,wlimit(i))
c                    enddo
                    exit do1082
                  endif
                endif
              endif
              nderiv=nderiv0
              do kc=1,nvar
                i=ivarele(kc)
                if(ival(i) .gt. 0)then
                  v00=rlist(latt(2,klp(i))+ival(i))
                else
                  v00=0.d0
                endif
                wvar(kc)=tweigh(latt(1,klp(i)),
     $               idtype(latt(1,klp(i))),
     $               ivvar(kc),bestval(kc),v00,
     $               avebeta,etamax,dpmax,emx,emy,brho,
     $               absweit)
                if(.not. nderiv)then
                  nderiv=idtype(latt(1,klp(i))) .eq. icSOL
                endif
              enddo
              npa=min(nvar,nparallel)
              if(nderiv)then
                call tffssetupqu(ifqu,ifqu0,nqumax,nqcol,nvar,lfno)
                ipr=-1
                if(npa .gt. 1)then
                  istep=npa
                  ip=0
                  do while(ipr .ne. 0 .and. ip .lt. npa-1)
                    ip=ip+1
                    ipr=fork_worker()
                    npr(ip)=ipr
                  enddo
                  if(ipr .gt. 0)then
                    ip=ip+1
                  endif
                else
                  ip=1
                  istep=1
                endif
                if(ipr .eq. 0)then
                  iuta1=itmmapp(nut*(2*nfam+1)*ntwissfun)
                else
                  iuta1=sad_loc(utwiss(1,-nfam,1))
                endif
                wcal1=wcal
                zcal1=zcal
                do kc=ip,nvar,istep
                  dvkc=max(abs(valvar(kc))*eps,abs(eps1/wvar(kc)))
                  valvar0=valvar(kc)
                  valvar(kc)=valvar0+dvkc
                  call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,nvar,
     $                 nele,klp,ival,couple,errk,iele,iele1,nlat)
                  call tffscalc(flv,twiss,pos,geo,gammab,rlist(iuta1),
     $                 latt,mult,ivarele,ivcomp,iele,iele1,ivvar,ival,
     $                 nvar,klp,itwissp,nut,nlat,nele,ndim,
     $                 dp,tracex,tracey,hstab,vstab,nfr,nfam,nfam1,
     $                 kdpa1,df1,
     $                 wfit,wiq,wsum,iqcola1,lfpa1,nqcola1,nqcol1a1,
     $                 nfcol,nfc0,maxcond,nlist,
     $                 ibegin,dfam,jfam,kfam,inicond,iuid,uini,
     $                 wake,iwakeelm,kwaketbl,nwakep,
     $                 ra1,rpa1,rstaba1,nstaba1,residuala1,wexponent,
     $                 offmw,etamax,avebeta,emx,emy,dpmax,coumin,
     $                 cell,zcal1,fitflg,geomet,cellstab,wcal1,.false.,
     $                 chgini,orbitcal,lfno,error)
                  valvar(kc)=valvar0-dvkc
                  call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,nvar,
     $                 nele,klp,ival,couple,errk,iele,iele1,nlat)
                  call tffscalc(flv,twiss,pos,geo,gammab,rlist(iuta1),
     $                 latt,mult,ivarele,ivcomp,iele,iele1,ivvar,ival,
     $                 nvar,klp,itwissp,nut,nlat,nele,ndim,
     $                 dp,tracex,tracey,hstab,vstab,nfr,nfam,nfam1,
     $                 kdpa2,df2,
     $                 wfit,wiq,wsum,iqcola2,lfpa2,nqcola2,nqcol1a2,
     $                 nfcol,nfc0,maxcond,nlist,
     $                 ibegin,dfam,jfam,kfam,inicond,iuid,uini,
     $                 wake,iwakeelm,kwaketbl,nwakep,
     $                 ra1,rpa1,rstaba1,nstaba1,residuala1,wexponent,
     $                 offmw,etamax,avebeta,emx,emy,dpmax,coumin,
     $                 cell,zcal1,fitflg,geomet,cellstab,wcal1,.false.,
     $                 chgini,orbitcal,lfno,error2)
                  valvar(kc)=valvar0
                  if(error .or. error2)then
                    call tclr(ddf1,nqcol)
                    call tclr(ddf2,nqcol)
                  else
                    call tffsddf(ddf1,df,df1,flv%iqcol,iqcola1,flv%lfp,
     $                   lfpa1,flv%kdp,kdpa1,nqcol,nqcola1)
                    call tffsddf(ddf2,df,df2,flv%iqcol,iqcola2,flv%lfp,
     $                   lfpa2,flv%kdp,kdpa2,nqcol,nqcola2)
                  endif
c                  ddf1(1:nqcol)=(ddf1(1:nqcol)-ddf2(1:nqcol))
c     $                 /2.d0/dvkc/wvar(kc)
                  do j=1,nqcol
                    rlist((kc-1)*nqcol+j+ifqu-1)=
     $                   (ddf1(j)-ddf2(j))/2.d0/dvkc/wvar(kc)
                  enddo
                enddo
                call tffswait(ipr,npa,npr,iuta1,
     $               'tffsmatch-NumDerv',irtc)
              else
                call tffsqu(flv,nqcol,nqcol1,nvar,nqumax,ifqu0,ifqu,
     $               latt,ivarele,ivvar,ivcomp,
     $               free,ival,iele,iele1,iele2,
     $               couple,nlat,nele,
     $               nfam,nfam1,utwiss,itwissp,nut,gammab,geo,pos,
     $               nparallel,cell,lfno,irtc)
                if(irtc .ne. 0)then
                  irtc=20003
                  exit do9000
                endif
                do kc=1,nvar
                  do j=1,nqcol1
                    kqu=(kc-1)*nqcol+j+ifqu-1
                    rlist(kqu)=rlist(kqu)*wiq(j)/wvar(kc)
                  enddo
                enddo
                if(nqcol .gt. nqcol1)then
                  ipr=-1
                  if(npa .gt. 1)then
                    istep=npa
                    ip=0
                    do while(ipr .ne. 0 .and. ip .lt. npa-1)
                      ip=ip+1
                      ipr=fork_worker()
                      npr(ip)=ipr
                    enddo
                    if(ipr .gt. 0)then
                      ip=ip+1
                    endif
                  else
                    ip=1
                    istep=1
                  endif
                  do kc=ip,nvar,istep
                    i=ivarele(kc)
                    if(nqcol .gt. nqcol1)then
                      valvar(kc)=valvar(kc)+eps1/wvar(kc)
                      if(cell)then
                        call twmov(latt(2,1),twisss,1,0,.false.)
                      endif
                      call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,nvar,
     $                     nele,klp,ival,couple,errk,iele,iele1,nlat)
                      call twmov(latt(2,1),twiss,nlat,ndim,.true.)
                      if(zcal)then
                        call tfgeo(latt,geo,pos,gammab,
     $                       geomet .or. .not. fitflg)
                      endif
                      over=.false.
                      if(ibegin .ne. 1)then
                        do j=1,ntwissfun
                          twiss(ibegin,0,j)=utwiss(j,0,itwissp(ibegin))
                        enddo
                      else
                        twiss(1,0,3)=0.d0
                        twiss(1,0,6)=0.d0
                      endif
                      call qcell(latt,twiss,gammab,ibegin,0,
     1                     hstab(0),vstab(0),tracex(0),tracey(0),
     $                     .false.,over)
                      nqcol00=nqcol
                      nqcol=nqcol1
                      call tffsfitfun(nqcol,df1,flv%iqcol,flv%kdp,
     $                     maxcond,error)
                      if(error .or. nqcol .ne. nqcol00)then
                        irtc=20003
                      endif
                      valvar(kc)=valvar(kc)-eps1/wvar(kc)
                      call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,
     $                  nvar,nele,klp,ival,couple,errk,iele,iele1,nlat)
                      do j=nqcol1+1,nqcol
                        kqu=(kc-1)*nqcol+j+ifqu-1
                        rlist(kqu)=(df(j)-df1(j))/eps1
                      enddo
                    endif
                  enddo
                  call tffswait(ipr,npa,npr,int8(0),
     $                 'tffsmatch-EVDeriv',irtc)
                endif
              endif
              if(irtc .ne. 0)then
                exit do9000
              endif
              df0(1:nqcol)=df(1:nqcol)
              nqcol0=nqcol
              iqcol0(1:nqcol)=flv%iqcol(1:nqcol)
              wlimit(1:nvar)=1.d0
              nvara=nvar
            enddo do1082
 1082       if(newton)then
              call tfsolv(rlist(ifqu),rlist(ifqu0),
     $             df,dval,wlimit,nqcol,nvar,flv%iqcol,
     $             flv%kfitp,flv%mfitp,dg,wexponent,1.d-8/fact)
              if(wexponent .ne. 2.d0)then
                dg=dg*(rp0/wsum)**(1.d0-wexponent/2.d0)
              endif
              if(dg .gt. 0.d0)then
                newton=.false.
                df0(1:nqcol)=df(1:nqcol)
                wlimit(1:nvar)=1.d0
                nvara=nvar
                go to 1082
              endif
            else
              call tgrad(rlist(ifqu0),rlist(ifqu),
     $             df,dval,wlimit,wexponent,nqcol,nvar)
            endif
            limited=.false.
            do ii=1,nvar
              i=ivarele(ii)
              dv=dval(ii)*fact/wvar(ii)*wlimit(ii)
              valvar(ii)=bestval(ii)+dv
              call tffsvlimit(ii,i,latt(1,klp(i)),valvar(ii),
     $             bestval(ii),
     $             vl,vlim,vl1,vl2,ivvar,ival,nvar,limited1,dlim)
              if(limited1)then
                valvar(ii)=min(vl2,max(vl1,bestval(ii)))
                if(dv .ne. 0.d0)then
                  wlimit(ii)=wlimit(ii)*
     $                 min(abs((vl-valvar(ii))/dv),.3d0)
                  limited=.true.
                  if(wlimit(ii) .lt. wlmin)then
                    wlimit(ii)=0.d0
                    nvara=nvara-1
                  endif
                endif
              endif
            enddo
            if(limited .and. nvara .gt. 0)then
              df(1:nqcol)=df0(1:nqcol)
              go to 1082
            endif
          endif
        enddo do200
        call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,nvar,nele,
     $       klp,ival,couple,errk,iele,iele1,nlat)
        if(dlim)then
          call tffssetlimit(valvar,vlim,ivvar,ival,nvar,
     $         latt,ivarele,ivcomp,
     $         klp,couple,errk,iele,iele1,dlim)
        endif
      enddo do9000
      if(nqumax .gt. 0)then
        call tfree(ifqu0)
        call tmunmapp(ifqu)
      endif
      call tclrfpe
      return
      end

      subroutine tffsddf(ddf,df,df1,iqcol,iqcola1,lfp,lfpa1,
     $                 kdp,kdpa1,nqcol,nqcola1)
      implicit none
      integer*4 nqcol,nqcola1,iqcol(nqcol),iqcola1(nqcola1),
     $     lfp(2,nqcol),lfpa1(2,nqcola1),kdp(nqcol),kdpa1(nqcola1),
     $     i,j,j1
      real*8 ddf(nqcol),df(nqcol),df1(nqcola1)
      j=1
      do i=1,nqcol
        if(j .gt. nqcola1)then
          ddf(i)=0.d0
        elseif(iqcol(i) .eq. iqcola1(j) .and.
     $       lfp(1,i) .eq. lfpa1(1,j) .and.
     $       lfp(2,i) .eq. lfpa1(2,j) .and.
     $       kdp(i) .eq. kdpa1(j))then
          ddf(i)=df(i)-df1(j)
          j=j+1
        else
          do j1=j,nqcola1
            if(iqcol(i) .lt. iqcola1(j1))then
              ddf(i)=0.d0
              exit
            elseif(iqcol(i) .gt. iqcol(j1))then
              j=j1+1
            elseif(iqcol(i) .eq. iqcola1(j1) .and.
     $       lfp(1,i) .eq. lfpa1(1,j1) .and.
     $       lfp(2,i) .eq. lfpa1(2,j1) .and.
     $       kdp(i) .eq. kdpa1(j1))then
              ddf(i)=df(i)-df1(j1)
              if(j .eq. j1)then
                j=j+1
              endif
              exit
            endif
          enddo
          ddf(i)=0.d0
        endif
      enddo
      return
      end

      subroutine tffswait(ipr,npa,npr,kash,tag,irtc)
      use tfshare
      implicit none
      integer*4 ipr,npa,npr(npa),irtc,wait,ist, i,j
      integer*8 kash
      character*(*) tag
      if(ipr .eq. 0)then
        if(kash .ne. 0)then
          call tfreeshared(kash,-1)
        endif
        call tfresetsharedmap()
        call exit_without_hooks(0)
      elseif(ipr .gt. 0)then
        do i=1,npa-1
          dowait: do while(.true.)
            ipr=wait(ist)
            do j=1,npa-1
              if(npr(j) .eq. ipr)then
                irtc=0
                npr(j)=0
                exit dowait
              endif
            enddo
            if(ipr .ne. -1)then
              write(*,*)'???-'//tag//' Unexpected process:',ipr
            endif
          enddo dowait
          if(ist .ne. 0)then
            irtc=20003
          endif
        enddo
      endif
      return
      end

      subroutine tffssetupqu(ifqu,ifqu0,nqumax,nqcol,nvar,lfno)
      implicit none
      integer*8 ifqu,ifqu0,itmmapp,ktaloc
      integer*4 nqumax,nqu,nqcol,nvar,lfno
      nqu=nqcol*nvar
      if(nqu .gt. nqumax)then
        if(nqumax .gt. 0)then
          call tfree(ifqu0)
          call tmunmapp(ifqu)
        endif
        ifqu=itmmapp(nqu)
        if(ifqu .le. 0)then
          go to 9000
        endif
        ifqu0=ktaloc(nqu)
        if(ifqu0 .le. 0)then
          call tmunmapp(ifqu)
          go to 9000
        endif
        nqumax=nqu
      endif
      return
 9000 call termes(lfno,'?Too many conditions*variables.',' ')
      if(nqumax .gt. 0)then
        ifqu=itmmapp(nqumax)
        ifqu0=ktaloc(nqumax)
      endif
      return
      end

      subroutine tffssetlimit(valvar,vlim,ivvar,ival,nvar,
     $     latt,ivarele,ivcomp,
     $     klp,couple,errk,iele,iele1,dlim)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 ii,i,nvar,ivvar(nvar),ival(nele)
      integer*4 latt(2,nlat),iele(nlat),iele1(nlat),
     $     ivcomp(nvar),ivarele(nvar),klp(nele)
      real*8 valvar(nvar),vlim(nele,2),vl,vl1,vl2,
     $     couple(nlat),errk(2,nlat)     
      logical*4 dlim,limited,limited1
      limited=.false.
      do ii=1,nvar
        i=ivarele(ii)
        call tffsvlimit(ii,i,latt(1,klp(i)),valvar(ii),valvar(ii),
     $       vl,vlim,vl1,vl2,ivvar,ival,nvar,limited1,dlim)
        if(limited1)then
          limited=.true.
          valvar(ii)=vl
        endif
      enddo
      if(limited)then
        call tfsetv(latt,ivarele,ivvar,ivcomp,valvar,nvar,
     $       nele,klp,ival,couple,errk,iele,iele1,nlat)
      endif
      return
      end

      subroutine tffsvlimit(ii,i,ld,val,val0,vl,vlim,
     $     vl1,vl2,ivvar,ival,nvar,limited,dlim)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_list), pointer :: klx
      integer*8 kx
      integer*4 ii,i,ld,nvar,ivvar(nvar),ival(nele),ltyp,irtc
      real*8 val,val0,vlim(nele,2),vl,vl1,vl0,vl2
      logical*4 limited,dlim
      limited=.false.
      ltyp=idtype(ld)
      vl=val
      vl0=val0
      if(ivvar(ii) .eq. ival(i))then
        vl1=vlim(i,1)
        vl2=vlim(i,2)
        if(vl .lt. vl1)then
          vl=vl1
          limited=.true.
        elseif(vl .gt. vl2)then
          vl=vl2
          limited=.true.
        endif
        vl0=max(vl1,min(vl2,vl0))
        if(.not. bipol)then
          if(vl*vl0 .lt. 0.d0)then
            vl=0.d0
            limited=.true.
          endif
        endif
      else
        vl1=-1.d100
        vl2=1.d100
      endif
      if(ltyp .eq. 41)then
        if(ivvar(ii) .eq. mfitbx .or.ivvar(ii) .eq. mfitby)then
          if(vl .le. 1.d-9)then
            vl=0.d0
            limited=.true.
          endif
        endif
      endif
      call tffsvarfun(1,ld,ivvar(ii),val,kx,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ktfrealq(kx))then
        dlim=.true.
        if(kx .eq. 0)then
          vl=vl0
          go to 2009
        endif
      elseif(tfreallistq(kx,klx))then
        if(klx%nl .eq. 2)then
          dlim=.true.
          vl1=max(vl1,klx%rbody(1))
          vl2=min(vl2,klx%rbody(2))
          if(vl .lt. vl1)then
            vl=vl1
            go to 2009
          endif
          if(vl .gt. vl2)then
            vl=vl2
            go to 2009
          endif
        endif
      endif
      return
 2009 limited=.true.
      return
      end

      subroutine tffsvarfun(id,ld,k,x,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 kx
      integer*4 id,ld,irtc,isp1,level,itfuplevel,ltyp,
     $     itfdownlevel,k
      real*8 x
      character*(MAXPNAME) vn,tfkwrd
      integer*8 ifvr,ifvw,ifvvarn,ifvkey
      save ifvr,ifvw,ifvvarn,ifvkey
      data ifvr /0/
      if(ifvr .eq. 0)then
        ifvr=ktfsymbolz('`VariableRange',14)
        ifvw=ktfsymbolz('`VariableWeight',15)
        ifvvarn=ktsalocb(0,'        ',MAXPNAME+16)
        ifvkey=ktsalocb(0,'        ',MAXPNAME)
      endif
      ltyp=idtype(ld)
      ilist(1,ifvvarn)=len_trim(pname(ld))
      call tfpadstr(pname(ld),ifvvarn+1,ilist(1,ifvvarn))
      vn=tfkwrd(ltyp,k)
      ilist(1,ifvkey)=len_trim(vn)
      call tfpadstr(vn,ifvkey+1,ilist(1,ifvkey))
      isp1=isp+1
      if(id .eq. 1)then
        ktastk(isp1)=ktfsymbol+ifvr
      elseif(id .eq. 2)then
        ktastk(isp1)=ktfsymbol+ifvw
      endif
      isp=isp+2
      ktastk(isp)=ktfstring+ifvvarn
      isp=isp+1
      ktastk(isp)=ktfstring+ifvkey
      isp=isp+1
      rtastk(isp)=x
      call tclrfpe
      level=itfuplevel()
      call tfefunref(isp1,kx,.false.,irtc)
      isp=isp1-1
      if(irtc .ne. 0)then
        level=itfdownlevel()
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',2,6)
        endif
        if(id .eq. 1)then
          call termes(6,'Error in VariableRange '//
     $         pname(ld)//' '//vn,' ')
        elseif(id .eq. 2)then
          call termes(6,'Error in VariableWeight '//
     $         pname(ld)//' '//vn,' ')
        endif
      else
        call tfconnectk(kx,irtc)
      endif
      return
      end

      subroutine tffsqu(flv,nqcol,nqcol1,nvar,nqumax,ifqu0,ifqu,
     $     latt,ivarele,ivvar,ivcomp,
     $     free,ival,iele,iele1,iele2,
     $     couple,nlat,nele,
     $     nfam,nfam1,utwiss,itwissp,nut,gammab,geo,pos,
     $     nparallel,cell,lfno,irtc)
      use tfstk
      use ffslocal, only:ffslocalv
      use tffitcode
      use tfshare
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      type (ffslocalv) flv
      integer*8 itmmapp,ifqu,ifqu0,ktaloc,kcm,kkqu,kqu,ic
      integer*4 nqcol,nqcol1,nvar,nqumax,nlat,nele,
     $     irtc,lfno,nut,nfam,nfam1
      integer*4 latt(2,nlat),ivarele(nvar),ivvar(nvar),
     $     iele(nlat),iele1(nlat),
     $     iele2(nlat),ival(nele),itwissp(nlat),
     $     ivcomp(nvar),npp
      real*8 utwiss(ntwissfun,-nfam:nfam,nut),gammab(nlat),
     $     geo(3,4,nlat),
     $     couple(nlat),pos(nlat),frbegin,frend
      logical*4 free(nele),cell
      integer*4 nqu,k,kk,i,kq,j,kf,lp,kp,iv,kkf,kkq,kkk,
     $     ii,ltyp,jj,lbegin,lend,kc,ik1,nparallel,
     $     iclast(-nfam:nfam),iec,ik,nk,kk1,
     $     ip,ipr,istep,npr(nparallel),fork_worker
      real*8 s,dtwiss(mfittry),coup,posk,wk,ctrans(27,-nfam:nfam)
      logical*4 col(2,nqcol),disp,nzcod
      call tffscoupmatrix(iele2,kcm,lfno)
      irtc=0
      nqu=nqcol*nvar
      do9000: do kkk=1,1
        if(nqu .gt. nqumax)then
          if(nqumax .gt. 0)then
            call tfree(ifqu0)
            call tmunmapp(ifqu)
          endif
          ifqu=itmmapp(nqu)
          if(ifqu .le. 0)then
            exit
          endif
          ifqu0=ktaloc(nqu)
          if(ifqu0 .le. 0)then
            call tmunmapp(ifqu)
            exit
          endif
          nqumax=nqu
        endif
        call tffsbound(nlat,latt,lbegin,frbegin,lend,frend)
        call tclr(rlist(ifqu),nqu)
        npp=min(nvar,nparallel)
        ipr=-1
        if(npp .gt. 1)then
          istep=npp
          ip=0
          do while(ipr .ne. 0 .and. ip .lt. npp-1)
            ip=ip+1
            ipr=fork_worker()
            npr(ip)=ipr
          enddo
          if(ipr .ne. 0)then
            ip=ip+1
          endif
        else
          ip=1
          istep=1
        endif  
        do k=1,nlat-1
          kc=iele(k)
          kk=iele1(kc)
          kk1=iele1(k)
          iec=iele2(k)
          if(iec .eq. 0)then
            nk=0
          else
            nk=ilist(1,iec)
          endif
          if(kk .gt. 0 .and. free(kk) .or. iec .ne. 0)then
            posk=pos(k)
            wk=1.d0
            if(k .eq. lbegin)then
              wk=1.d0-frbegin
            endif
            if(k .eq. lend)then
              wk=frend
            endif
            ltyp=idtype(latt(1,k))
            do ii=ip,nvar,istep
              if(kk .le. 0 .or.
     $             (.not. free(kk) .and. .not. free(kk1)))then
                ik1=1
              else
                ik1=0
                do ik=1,nk
                  if(ivvar(ii) .eq. ilist(2,iec+ik))then
                    ik1=1
                    exit
                  endif
                enddo
              endif
              do11: do ik=ik1,nk
                if(ik .ne. 0)then
                  if(kcm .eq. 0)then
                    cycle
                  else
                    ic=kcm+(ilist(1,iec+ik)-1)*nvar-1
                    coup=rlist(ic+ii)*wk
                    if(coup .ne. 0.d0)then
                      iv=ilist(2,iec+ik)
                    else
                      cycle
                    endif
                  endif
                else
                  iv=ivvar(ii)
c     write(*,*)'tffssqu ',iv,ival(kk),kk,kk1,
c     $               ivarele(ii),k,ivcomp(ii),ii
                  if(iv .eq. ival(kk) .and. ivarele(ii) .eq. kk .and.
     $                 (ivcomp(ii) .eq. 0 .or.
     $                 ivcomp(ii) .eq. kc))then
                    coup=couple(k)*wk
                  elseif(iv .ne. ival(kk) .and.
     $                   ivarele(ii) .eq. kk1 .and.
     $                   (ivcomp(ii) .eq. 0 .or.
     $                   ivcomp(ii) .eq. k))then
                    coup=wk
                  else
                    cycle
                  endif
                endif
                iclast(nfam1:nfam)=0
                do kq=1,nqcol1
                  col(1,kq)=.true.
                  col(2,kq)=flv%lfp(2,kq) .gt. 0
                enddo
                do kq=1,nqcol1
                  kqu=(ii-1)*nqcol+kq+ifqu-1
                  s=coup
                  do i=1,2
                    if(col(i,kq))then
                      j=flv%iqcol(kq)
                      kf=flv%kfit(flv%kfitp(j))
                      lp=flv%lfp(i,kq)
                      if(kf .le. mfittry)then
                        if(k .lt. lp .or.
     $                       posk .lt. pos(lp) .or. cell)then
                          nzcod= kf .eq. mfitdz
     $                         .or. kf .eq. mfitddp
                          disp=kf .ge. mfitex .and. kf .le. mfitepy
     $                         .or. nzcod
                          if(.not. disp)then
                            do kkq=kq+1,nqcol1
                              jj=flv%iqcol(kkq)
                              kkf=flv%kfit(flv%kfitp(jj))
                              disp=kkf .ge. mfitex
     $                             .and. kkf .le. mfitepy
                              if(disp)then
                                exit
                              endif
                            enddo
                          endif
                          kp=flv%kdp(kq)
                          call qdcell(latt,
     1                         utwiss,itwissp,gammab,dtwiss,
     $                         k,lp,kp,iv,ctrans,iclast,
     $                         nfam,nut,disp,nzcod)
                          rlist(kqu)=rlist(kqu)+s*dtwiss(kf)
                          do kkq=kq+1,nqcol1
                            jj=flv%iqcol(kkq)
                            kkf=flv%kfit(flv%kfitp(jj))
                            if(kkf .le. mfito .and.
     $                           kp .eq. flv%kdp(kkq))then
                              kkqu=(ii-1)*nqcol+kkq+ifqu-1
                              if(col(1,kkq) .and.
     $                             flv%lfp(1,kkq) .eq. lp)then
                                rlist(kkqu)=rlist(kkqu)
     $                               +dtwiss(kkf)*coup
                                col(1,kkq)=.false.
                              elseif(col(2,kkq) .and.
     $                               flv%lfp(2,kkq) .eq. lp)then
                                rlist(kkqu)=rlist(kkqu)
     $                               -dtwiss(kkf)*coup
                                col(2,kkq)=.false.
                              endif
                            endif
                          enddo
                        endif
                      elseif(kf .eq. mfitleng)then
                        if(kytbl(kwL,ltyp) .eq. iv)then
                          if(posk .lt. pos(lp) .and. posk .ge. 0.d0)then
                            rlist(kqu)=rlist(kqu)+s
                          elseif(posk .ge. pos(lp) .and.
     $                           posk .lt. 0.d0)then
                            rlist(kqu)=rlist(kqu)-s
                          endif
                        endif
c                        write(*,'(a,1p3g15.7,3i5)')'tffsqu ',
c     $                       posk,pos(lp),rlist(kqu),ltyp,iv
                      else
                        call tdgeo(latt,utwiss,itwissp
     $                       ,gammab,s,rlist(kqu),
     $                       kf,lp,k,geo,ltyp,iv,nut,nfam)
                      endif
                    endif
                    s=-s
                  enddo
                enddo
              enddo do11
            enddo
          endif
        enddo
        call tffswait(ipr,npp,npr,int8(0),'tffsqu',irtc)
        if(kcm .ne. 0)then
          call tfree(kcm)
        endif
        return
      enddo do9000
      call termes(lfno,'?Too many conditions*variables.',' ')
      if(nqumax .gt. 0)then
        ifqu=itmmapp(nqumax)
        ifqu0=ktaloc(nqumax)
      endif
      irtc=20002
      if(kcm .ne. 0)then
        call tfree(kcm)
      endif
      return
      end

      real*8 function tffsfmin(f1,f2,g1,g2,g0,dg)
      implicit none
      real*8 f1,f2,g1,g2,g0,dg,a,b
      if(f2 .eq. 0.d0)then
        tffsfmin=-.5d0*f1*dg/((g1-g0)/f1-dg)
      else
        a=((g1-g0)/f1**2-(g2-g0)/f2**2)/(f1-f2)+dg/f1/f2
        b=(-f2*(g1-g0)/f1**2+f1*(g2-g0)/f2**2)-dg*(f1+f2)/f1/f2
        if(b .gt. 0.d0)then
          tffsfmin=-dg/(sqrt(max(0.d0,b**2-3.d0*a*dg))+b)
        else
          tffsfmin=(sqrt(max(0.d0,b**2-3.d0*a*dg))-b)/3.d0/a
        endif
      endif
      tffsfmin=min(.577d0*f1,max(f1/16.d0,tffsfmin))
      return
      end

      subroutine twmov(lp,twiss,n1,n2,right)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 lp,n1,n2,i,ntfun
      real*8 twiss(n1,-n2:n2,1:ntwissfun)
      logical*4 right
c
      if(orbitcal)then
        ntfun=ntwissfun
      else
        ntfun=mfitdetr
      endif
      if(right)then
        twiss(1,0,1:ntfun)=rlist(lp+1:lp+ntfun)
        if(rlist(lp+ilist(1,lp)) .lt. 0.d0)then
          twiss(1,0,mfitax)=-rlist(lp+mfitax)
          twiss(1,0,mfitay)=-rlist(lp+mfitay)
          twiss(1,0,mfitepx)=-rlist(lp+mfitepx)
          twiss(1,0,mfitepy)=-rlist(lp+mfitepy)
          twiss(1,0,mfitr2)=-rlist(lp+mfitr2)
          twiss(1,0,mfitr3)=-rlist(lp+mfitr3)
          if(orbitcal)then
            twiss(1,0,mfitdpx)=-rlist(lp+mfitdpx)
            twiss(1,0,mfitdpy)=-rlist(lp+mfitdpy)
          endif
        endif
      else
        do i=1,ntwissfun
          rlist(lp+i)=twiss(1,0,i)
        enddo
        if(rlist(lp+ilist(1,lp)) .lt. 0.d0)then
          rlist(lp+mfitax)=-twiss(1,0,mfitax)
          rlist(lp+mfitay)=-twiss(1,0,mfitay)
          rlist(lp+mfitepx)=-twiss(1,0,mfitepx)
          rlist(lp+mfitepy)=-twiss(1,0,mfitepy)
          rlist(lp+mfitr2)=-twiss(1,0,mfitr2)
          rlist(lp+mfitr3)=-twiss(1,0,mfitr3)
          rlist(lp+mfitdpx)=-twiss(1,0,mfitdpx)
          rlist(lp+mfitdpy)=-twiss(1,0,mfitdpy)
        endif
      endif
      return
      end

      subroutine tgrad(qu,qu0,df,grad,
     $wlimit,wexponent,nqcol,nvar)
      implicit none
      integer*4 nqcol,nvar,i,j
      real*8 qu(nqcol,nvar),qu0(nqcol,nvar),df(nqcol),
     $     grad(nvar),wlimit(nvar),s,sg,r,wexponent,
     $     dfw(nqcol),dfwi
      do i=1,nvar
        qu(:,i)=qu0(:,i)*wlimit(i)
c        do j=1,nqcol
c          qu(j,i)=qu0(j,i)*wlimit(i)
c        enddo
      enddo
      r=0.d0
      do i=1,nqcol
        if(df(i) .ne. 0.d0)then
          dfwi=abs(df(i))**wexponent
          r=r+dfwi
          dfw(i)=dfwi/df(i)
        else
          dfw(i)=0.d0
        endif
      enddo
      sg=0.d0
      do i=1,nvar
        s=0.d0
        do j=1,nqcol
          s=s+qu(j,i)*dfw(j)
        enddo
        grad(i)=s
        sg=sg+s**2
      enddo
      if(sg .eq. 0.d0)then
        return
      endif
      s=r/sg
      grad=grad*s
      return
      end

      real*8 function tweigh(i,ltyp,iv,val0,vk,
     $     avebeta,etamax,dpmax,emx,emy,brho,absweit)
      use tfstk
      use tfcode
      implicit none
      include 'inc/CBKMAC.inc'
      integer*8 kx
      integer*4 i,ltyp,iv,irtc
      real*8 val0,avebeta,dpmax,emx,emy,etamax,gw,vmin,brho,rfromk,vk
      logical*4 absweit
      gw=1.d0
      if(iv .eq. kytbl(kwK1,ltyp))then
      elseif(iv .eq. kytbl(kwL,ltyp))then
        gw=1.d0/avebeta**2
      elseif(iv .eq. kytbl(kwANGL,ltyp) .or.
     $       iv .eq. kytbl(kwK0,ltyp) .or.
     $       iv .eq. kytbl(kwSK0,ltyp))then
        gw=1.d0/avebeta
      elseif(iv .eq. kytbl(kwK2,ltyp))then
        gw=max(1.d-3,dpmax)*max(1.d-2,etamax)
      elseif(iv .eq. kytbl(kwK3,ltyp))then
        gw=(max(1.d-3,dpmax)*max(1.d-2,etamax))**2
      elseif(iv .eq. kytbl(kwK4,ltyp))then
        gw=(max(1.d-3,dpmax)*max(1.d-2,etamax))**3
      elseif(iv .eq. kytbl(kwK4,ltyp))then
        gw=(max(1.d-3,dpmax)*max(1.d-2,etamax))**4
      elseif(ltyp .eq. icMULT)then
        if(iv .ge. kytbl(kwK1,icMULT))then
          gw=(max(1.d-3,dpmax)*max(1.d-2,etamax))
     $         **((iv-kytbl(kwK1,icMULT))/2+1)
        endif
      elseif(iv .eq. kytbl(kwVOLT,ltyp))then
        gw=100.d0/max(1.d0,abs(vk))
      elseif(iv .eq. kytbl(kwFREQ,ltyp))then
        gw=100.d0/max(1.d0,abs(vk))
      elseif(iv .eq. kytbl(kwDX,ltyp) .or.
     $       iv .eq. kytbl(kwDY,ltyp))then
        if(kytbl(kwK1,ltyp) .ne. 0)then
          gw=0.1d0/avebeta
        endif
      elseif(iv .eq. kytbl(kwBZ,ltyp))then
        gw=1.d0/brho
      elseif(ltyp .eq. icMARK)then
        if(iv .eq. kytbl(kwAX,icMARK) .or.
     $       iv .eq. kytbl(kwAY,icMARK) .or.
     $       iv .eq. kytbl(kwR1,icMARK) .or.
     $       iv .eq. kytbl(kwR4,icMARK))then
          gw=1.d0/avebeta
        elseif(iv .eq. kytbl(kwBX,icMARK) .or.
     $         iv .eq. kytbl(kwBY,icMARK))then
          gw=1.d0/avebeta**2
        elseif(iv .eq. kytbl(kwEX,icMARK) .or.
     $         iv .eq. kytbl(kwEPX,icMARK))then
          gw=max(1.d-3,dpmax)*sqrt(avebeta/(abs(emx)+abs(emy)))
     $         /avebeta**2
        elseif(iv .eq. kytbl(kwEPX,icMARK) .or.
     $         iv .eq. kytbl(kwEPX,icMARK))then
          gw=max(1.d-3,dpmax)*sqrt(avebeta/(abs(emx)+abs(emy)))
     $         /avebeta
        elseif(iv .eq. kytbl(kwR2,icMARK))then
          gw=1.d0/avebeta**2
        endif
      endif
      if(.not. absweit)then
        vmin=1.d-6
        if(iv .eq. kytbl(kwK1,ltyp))then
          vmin=1.d-5
        elseif(iv .eq. kytbl(kwK2,ltyp))then
          vmin=1.d-3
        elseif(iv .eq. kytbl(kwK3,ltyp))then
          vmin=1.d-1
        elseif(iv .eq. kytbl(kwK4,ltyp))then
          vmin=1.d1
        elseif(iv .eq. kytbl(kwK5,ltyp))then
          vmin=1.d3
        elseif(iv .eq. kytbl(kwK6,ltyp))then
          vmin=1.d5
        elseif(ltyp .eq. icMULT)then
          if(iv .ge. kytbl(kwK1,icMULT))then
            vmin=10.d0**(((iv-kytbl(kwK1,icMULT))/2)*2-5)
          endif
        endif
        gw=sqrt(max(vmin,abs(val0))/gw)
      endif
      call tffsvarfun(2,i,iv,gw,kx,irtc)
      if(irtc .eq. 0 .and. ktfrealq(kx))then
        gw=rfromk(kx)
      endif
      tweigh=gw
      return
      end

      subroutine tfsolv(qu,qu0,df,dval,wlimit,nqcol,nvar,
     $     iqcol,kfitp,mfitp,dg,wexponent,eps)
      implicit none
      integer*4 nqcol,nvar,j
      real*8 qu(nqcol,nvar),qu0(nqcol,nvar),df(nqcol),dval(nvar)
      integer*4 iqcol(nqcol),kfitp(*),mfitp(*)
      real*8 b(nqcol),s,eps,dg,wexponent,wlimit(nvar)
      logical*4 fit(nqcol),again,allneg
      integer*4 nagain,i,nj
      allneg=.true.
      do i=1,nqcol
        fit(i)=mfitp(kfitp(iqcol(i))) .gt. 0
        allneg=allneg .and. .not. fit(i)
      enddo
      if(allneg)then
        fit=.true.
      endif
      nagain=0
1     nj=0
      do i=1,nqcol
        if(fit(i))then
          nj=nj+1
          do j=1,nvar
            qu0(nj,j)=qu(i,j)*wlimit(j)
          enddo
          b(nj)=df(i)
        endif
      enddo
      call tsolva(qu0,b,dval,nj,nvar,nqcol,eps)
      again=.false.
      dg=0.d0
      do i=1,nqcol
        s=0.d0
        do j=1,nvar
          s=s+qu(i,j)*wlimit(j)*dval(j)
        enddo
        if(df(i) .ne. 0.d0)then
          if(wexponent .eq. 2.d0)then
            dg=dg-df(i)*s
          else
            dg=dg-abs(df(i))**wexponent/df(i)*s
          endif
        endif
        if(.not. fit(i))then
          if((s-df(i))*df(i) .lt. 0.d0)then
            fit(i)=.true.
            again=.true.
          endif
        endif
      enddo
      if(again)then
        nagain=nagain+1
        if(nagain .le. nqcol)then
          go to 1
        endif
      endif
      dg=wexponent*dg
      return
      end

      integer*8 function itmmapp(n)
      use tfmem
      use tfstk
      use tfshare
      implicit none
      include 'inc/TMACRO1.inc'
      integer*8 ktaloc
      integer*4 n,irtc
      if(nparallel .gt. 1)then
        irtc=1
        itmmapp=ktfallocshared(n)
c        write(*,*)'mmapp ',itmmapp,n
      else
        itmmapp=ktaloc(n)
      endif
      lastpend=max(lastpend,itmmapp)
      return
      end

      subroutine tmunmapp(i)
      use tfstk
      use tfshare
      implicit none
      include 'inc/TMACRO1.inc'
      integer*8 i
      if(nparallel .gt. 1)then
        call tfreeshared(i)
c        if(mapfree(rlist(i)) .ne. 0)then
c          write(*,*)'???tmunmapp-error in unmap.'
c        endif
      else
        call tfree(i)
      endif
      return
      end

      subroutine tffscoupmatrix(iele2,kcm,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 km,k1,kam,kcm,ktfmaloc,k2
      integer*4 iele2(nlat),lfno,irtc,n,m
      real*8 rfromk
      integer*8 itfcoupm
      data itfcoupm /0/
      if(iele2(nlat) .eq. 0)then
        kcm=0
        return
      endif
      if(itfcoupm .eq. 0)then
        itfcoupm=ktfsymbolz('CouplingMatrix',14)
      endif
      levele=levele+1
      call tfsyeval(itfcoupm,km,irtc)
      call tfconnectk(km,irtc)
      if(irtc .ne. 0)then
        go to 9010
      endif
      if(.not. tflistqk(km))then
        go to 9000
      endif
      kam=ktfaddr(km)
      k1=klist(kam+1)
      if(ktfnonrealq(k1) .or. rfromk(k1) .le. 0.d0)then
        go to 9100
      endif
      k2=klist(kam+2)
      kcm=ktfmaloc(k2,n,m,.false.,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9010
      endif
      return
 9010 if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
 9000 call termes(lfno,
     $     'Error or Non-numeric results in coupling matrix',' ')
 9100 kcm=0
      return
      end
