      subroutine tffsmatch(df,dp0,r,nparallel,lfno,irtc)
      use kyparam
      use tfstk
      use ffs, only: flv,dpmax,nele,ndim,nlat,maxcond,nvevx,nelvx
      use ffs_pointer
      use ffs_flag
      use ffs_fit
      use tffitcode
      use tfshare
      use iso_c_binding
      use mackw
      implicit none
c      type (sad_descriptor) kx
c      include 'DEBUG.inc'
      integer*8 ifqu,ifquw,iuta1,kqu
      real*8 , parameter :: flim1=-4.d0,flim2=-3.d0,aimp1=-1.8d0,
     $     aimp2=-.8d0,badc1=-3.5d0,badc2=-2.5d0,amedc1=-2.3d0,
     $     amedc2=-1.3d0,alit=0.75d0,wlmin=0.009d0,eps=1.d-5,
     $     eps1=1.d-8,rtol=1.05d0,rtol1=1.05d0
      real*8, parameter :: aloadmax=2.d4
      integer*4 ,intent(in):: nparallel,lfno
      integer*4 ,intent(out):: irtc
      integer*4 ibegin,nqumax,nqcol0,nqcol00,
     $     nqcola1,nqcol1a1,nqcola2,nqcol1a2,lout
      integer*4 kdpa1(maxcond),kdpa2(maxcond),iqcol0(maxcond),
     $     iqcola1(maxcond),iqcola2(maxcond),nretry,nstab,nstaba1
      real*8 ,intent(out):: df(maxcond),r,dp0
      real*8 dval(flv%nvar),df0(maxcond),df1(maxcond),df2(maxcond),
     $     ddf1(maxcond),ddf2(maxcond),
     $     residuala1(-ndimmax:ndimmax),v00
      logical*4 free(nele),zcal,error,error2,
     $     limited1,wcal1,zcal1
      integer*4 iter,ii,j,kc,nvara, i,
     $     lfpa1(2,maxcond),lfpa2(2,maxcond),
     $     ip,ipr,istep,npr(nparallel)
      real*8 rl,fuzz,a,b,x,crate,aimprv,fact,r0,r00,ra,alate,
     $     smallf,badcnv,amedcv,rstab,vl1,vl2,
     $     aitm1,aitm2,bestval(flv%nvar),wvar(flv%nvar),
     $     dg,f1,f2,g1,g2,ra1,rpa1,rstaba1,valvar0,
     $     rp,rp0,wlimit(flv%nvar),dv,vl,dvkc
      real*8 twisss(ntwissfun)
      real*8 , pointer :: qu(:,:),quw(:,:)
      logical*4 chgmod,newton,imprv,limited,over,wcal,
     $     parallel,nderiv,outt,nderiv0,dlim
      integer*4 kkk,kkkk,npa
      integer*4 , external :: itfgetrecl, fork_worker
      integer*8 , external :: itmmapp
      real*8 , external :: tffsfmin, tweigh
      character ch
      character*10 sexp
      integer*8 intffs,inumderiv,iexponent,inumw,iconvergence
      data intffs,inumderiv,iexponent,inumw,iconvergence /0,0,0,0,0/
c
      fuzz(x,a,b)=max(0.d0,min(1.d0,(x-a)/(b-a)))
c     
c     begin initialize for preventing compiler warning
      associate (nvar=>flv%nvar)
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
         call twmov(1,twisss,1,0,.true.)
      endif
      ibegin=1
      rstab=0.d0
      if(fitflg)then
        call tffssetlimit(nvar,dlim)
        fact=1.d0
        iter=0
        newton=.true.
        chgmod=.true.
        bestval(1:nvar)=nvevx(1:nvar)%valvar
        free=.false.
        free(nvevx(1:nvar)%ivarele)=.true.
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
          dp0=rlist(latt(1)+mfitddp)
          call tffscalc(flv%kdp,df,flv%iqcol,flv%lfp,
     $         nqcol,nqcol1,ibegin,
     $         r,rp,rstab,nstab,residual,
     $         zcal,wcal,parallel,lout,error)
          if(error)then
            if(irtc .eq. 20001)then
              exit do9000
            else
              fitflg=.false.
              irtc=20001
              bestval(1:nvar)=nvevx(1:nvar)%valvar
              if(cell)then
                call twmov(1,twisss,1,0,.false.)
              endif
              exit do200
            endif
          endif
          convgo=r .le. rl
          fitflg=fitflg .and. nqcol .gt. 0
          if(calexp)then
            sexp='  CALEXP'
          else
            sexp='  NOCALEXP'
          endif
          if(convgo)then
            write(lfno,9501)' Matched. (',r,')',
     $           dpmax,dp0,wexponent,ch,offmw,sexp
 9501       format(a,1pG11.4,a,
     $           ' DP =',0pf8.5,'  DP0 =',f8.5,
     $           '  ExponentOfResidual =',f4.1,a,
     $           ' OffMomentumWeight =',f8.3,a)
            fitflg=.false.
            if(.not. geomet)then
              call tfgeo(latt,.true.)
            endif
            exit do9000
          elseif(.not. fitflg)then
            write(lfno,9501)' Residual =',r,' ',
     $           dpmax,dp0,wexponent,ch,offmw,sexp
            exit do9000
          else
            if(chgini .and. cell)then
              call twmov(1,twisss,1,0,.true.)
            endif
c            chgini=.true.
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
                r00=r0
                ra=r0*1.000000001d0
                if(cell)then
                  call twmov(1,twisss,1,0,.true.)
                endif
              else
                imprv=r .lt. r0
                if(imprv)then
                  if(r .lt. r00/rtol1)then
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
                  bestval(1:nvar)=nvevx(1:nvar)%valvar
                  if(cell)then
                    call twmov(1,twisss,1,0,.true.)
                  endif
                  rp0=rp
                  r0=r
                else
                  aimprv=0.d0
                endif
                alate=fuzz(dble(iter),aitm1,aitm2)
                smallf=1.d0-fuzz(log10(fact),flim1,flim2)
                badcnv=1.d0-fuzz(log10(crate),badc1,badc2)
                amedcv=1.d0-fuzz(log10(crate),amedc1,amedc2)
                if(newton)then
                  chgmod=max(smallf,badcnv,min(alate,amedcv)) .gt. .5d0
                else
                  if(iter .gt. flv%itmax*10)then
                    fitflg=.false.
                    chgmod=.true.
                  elseif(aimprv .gt. .5d0)then
                    chgmod=.true.
                  elseif(max(smallf,badcnv,min(alate,amedcv))
     $                   .gt. .5d0)then
                    chgmod=.true.
                    if(nretry .gt. 0)then
                      nretry=nretry-1
                    else
                      chgini=.true.
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
                  nvevx(1:nvar)%valvar=bestval(1:nvar)
                  if(cell)then
                    call twmov(1,twisss,1,0,.false.)
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
                    call twmov(1,twisss,1,0,.false.)
                  endif
                  if(nvara .eq. nvar)then
                    a=fact/f1
                    nvevx(1:nvar)%valvar=nvevx(1:nvar)%valvar*a
     $                   +bestval(1:nvar)*(1.d0-a)
                    exit do200
                  else
                    nqcol=nqcol0
                    flv%iqcol(1:nqcol)=iqcol0(1:nqcol)
                    df(1:nqcol)=df0(1:nqcol)
                    nvara=nvar
                    wlimit(1:nvar)=max(wlmin,wlimit(1:nvar))
                    exit do1082
                  endif
                endif
              endif
              nderiv=nderiv0
              do kc=1,nvar
                i=nvevx(kc)%ivarele
                if(nelvx(i)%ival .gt. 0)then
                  v00=rlist(latt(nelvx(i)%klp)+nelvx(i)%ival)
                else
                  v00=0.d0
                endif
                wvar(kc)=tweigh(idelc(nelvx(i)%klp),
     $               idtypec(nelvx(i)%klp),
     $               nvevx(kc)%ivvar,bestval(kc),v00,absweit)
                if(.not. nderiv)then
                  nderiv=idtypec(nelvx(i)%klp) .eq. icSOL
                endif
              enddo
              nderiv=nderiv .or.
     $             cell .and. nstab .gt. 0
     $             .and. dble(nvar*nfam*nlat) .lt. aloadmax
              npa=min(nvar,nparallel)
              chgini=(nstab .eq. 0) .or. nderiv
              if(nderiv)then
                call tffssetupqu(ifqu,ifquw,nqumax,nqcol,nvar,lfno)
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
                  call c_f_pointer(c_loc(rlist(iuta1)),utwiss,
     $                 [ntwissfun,2*nfam+1,nut])
                  utwiss(1:ntwissfun,-nfam:nfam,1:nut)=>utwiss
                endif
                wcal1=wcal
                zcal1=zcal
                do kc=ip,nvar,istep
                  dvkc=max(abs(nvevx(kc)%valvar)*eps,abs(eps1/wvar(kc)))
                  valvar0=nvevx(kc)%valvar
                  nvevx(kc)%valvar=valvar0+dvkc
                  call tfsetv(nvar)
                  call tffscalc(kdpa1,df1,iqcola1,lfpa1,
     $                 nqcola1,nqcol1a1,ibegin,
     $                 ra1,rpa1,rstaba1,nstaba1,residuala1,
     $                 zcal1,wcal1,.false.,lfno,error)
                  nvevx(kc)%valvar=valvar0-dvkc
                  call tfsetv(nvar)
                  call tffscalc(kdpa2,df2,iqcola2,lfpa2,
     $                 nqcola2,nqcol1a2,ibegin,
     $                 ra1,rpa1,rstaba1,nstaba1,residuala1,
     $                 zcal1,wcal1,.false.,lfno,error2)
                  nvevx(kc)%valvar=valvar0
                  if(error .or. error2)then
                    ddf1(1:nqcol)=0.d0
                    ddf2(1:nqcol)=0.d0
                  else
                    call tffsddf(ddf1,df,df1,flv%iqcol,iqcola1,flv%lfp,
     $                   lfpa1,flv%kdp,kdpa1,nqcol,nqcola1)
                    call tffsddf(ddf2,df,df2,flv%iqcol,iqcola2,flv%lfp,
     $                   lfpa2,flv%kdp,kdpa2,nqcol,nqcola2)
                  endif
                  rlist((kc-1)*nqcol+ifqu:kc*nqcol+ifqu-1)=
     $                 (ddf1(1:nqcol)-ddf2(1:nqcol))/2.d0/dvkc/wvar(kc)
                enddo
                call tffswait(ipr,npa,npr,iuta1,
     $               'tffsmatch-NumDerv',irtc)
              else
                call tffsqu(nqcol,nqcol1,nvar,nqumax,ifquw,ifqu,
     $               free,nlat,nele,nfam,nfam1,nut,
     $               nparallel,cell,lfno,irtc)
                if(irtc .ne. 0)then
                  irtc=20003
                  exit do9000
                endif
                do concurrent (kc=1:nvar)
                  kqu=(kc-1)*nqcol+ifqu
                  rlist(kqu:kqu+nqcol1-1)=rlist(kqu:kqu+nqcol1-1)
     $                 *wiq(1:nqcol1)/wvar(kc)
c     enddo
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
                    i=nvevx(kc)%ivarele
                    if(nqcol .gt. nqcol1)then
                      nvevx(kc)%valvar=nvevx(kc)%valvar+eps1/wvar(kc)
                      if(cell)then
                        call twmov(1,twisss,1,0,.false.)
                      endif
                      call tfsetv(nvar)
                      call twmov(1,twiss,nlat,ndim,.true.)
                      if(zcal)then
                        call tfgeo(latt,geomet .or. .not. fitflg)
                      endif
                      over=.false.
                      if(ibegin .ne. 1)then
                        twiss(ibegin,0,1:ntwissfun)=
     $                       utwiss(1:ntwissfun,0,itwissp(ibegin))
                      else
                        twiss(1,0,3)=0.d0
                        twiss(1,0,6)=0.d0
                      endif
                      call qcell(0,optstat(0),.false.)
                      nqcol00=nqcol
                      nqcol=nqcol1
                      call tffsfitfun(nqcol,df1,flv%iqcol,flv%kdp,
     $                     maxcond,error)
                      if(error .or. nqcol .ne. nqcol00)then
                        irtc=20003
                      endif
                      nvevx(kc)%valvar=nvevx(kc)%valvar-eps1/wvar(kc)
                      call tfsetv(nvar)
                      do concurrent (j=nqcol1+1:nqcol)
                        kqu=(kc-1)*nqcol+j+ifqu-1
                        rlist(kqu)=(df(j)-df1(j))/eps1
                      enddo
                    endif
                  enddo
                  call tffswait(ipr,npa,npr,i00,
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
              call c_f_pointer(c_loc(rlist(ifqu)),qu,[nqcol,nvar])
              call c_f_pointer(c_loc(rlist(ifquw)),quw,[nqcol,nvar])
              call tfsolv(qu,quw,
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
              call c_f_pointer(c_loc(rlist(ifqu)),qu,[nqcol,nvar])
              call c_f_pointer(c_loc(rlist(ifquw)),quw,[nqcol,nvar])
              call tgrad(qu,quw,df,dval,wlimit,wexponent,nqcol,nvar)
            endif
            limited=.false.
            do ii=1,nvar
              i=nvevx(ii)%ivarele
              dv=dval(ii)*fact/wvar(ii)*wlimit(ii)
              nvevx(ii)%valvar=bestval(ii)+dv
              call tffsvlimit(i,idelc(nelvx(i)%klp),nvevx(ii)%valvar,
     $             bestval(ii),vl,vl1,vl2,nvevx(ii)%ivvar,
     $             limited1,dlim)
              if(limited1)then
                nvevx(ii)%valvar=min(vl2,max(vl1,bestval(ii)))
                if(dv .ne. 0.d0)then
                  wlimit(ii)=wlimit(ii)*
     $                 min(abs((vl-nvevx(ii)%valvar)/dv),.3d0)
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
        call tfsetv(nvar)
        if(dlim)then
          call tffssetlimit(nvar,dlim)
        endif
      enddo do9000
      if(nqumax .gt. 0)then
        call tfree(ifquw)
        call tmunmapp(ifqu)
      endif
      call tclrfpe
      return
      end associate
      end

      subroutine tffsddf(ddf,df,df1,iqcol,iqcola1,lfp,lfpa1,
     $                 kdp,kdpa1,nqcol,nqcola1)
      implicit none
      integer*4 ,intent(in):: nqcol,nqcola1,iqcol(nqcol),
     $     iqcola1(nqcola1),lfp(2,nqcol),lfpa1(2,nqcola1),kdp(nqcol),
     $     kdpa1(nqcola1)
      integer*4 i,j,j1
      real*8 ,intent(out):: ddf(nqcol)
      real*8 ,intent(in):: df(nqcol),df1(nqcola1)
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
      integer*4 ,intent(out):: irtc,ipr
      integer*4 ,intent(in):: npa
      integer*4 ,intent(inout):: npr(npa)
      integer*4 ist,i,j,wait
      integer*8 ,intent(in):: kash
      character*(*) ,intent(in):: tag
      if(ipr .eq. 0)then
        if(kash .ne. 0)then
          call tfreeshared(kash,-1)
        endif
        call tfresetsharedmap()
        call exit_without_hooks(0)
      elseif(ipr .gt. 0)then
        irtc=0
        do i=1,npa-1
          dowait: do while(.true.)
            ipr=wait(ist)
            do j=1,npa-1
              if(npr(j) .eq. ipr)then
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

      subroutine tffssetupqu(ifqu,ifquw,nqumax,nqcol,nvar,lfno)
      use tfmem, only:ktaloc,tfree
      implicit none
      integer*8 ,intent(out):: ifqu,ifquw
      integer*8 itmmapp
      integer*4 ,intent(inout):: nqumax
      integer*4 ,intent(in):: nqcol,nvar,lfno
      integer*4 nqu
      integer*4 , parameter :: minnqu=512
      nqu=max(minnqu,nqcol*nvar)
      if(nqu .gt. nqumax)then
        if(nqumax .gt. 0)then
          call tfree(ifquw)
          call tmunmapp(ifqu)
        endif
        ifqu=itmmapp(nqu)
        if(ifqu .le. 0)then
          go to 9000
        endif
        ifquw=ktaloc(nqu)
        if(ifquw .le. 0)then
          call tmunmapp(ifqu)
          go to 9000
        endif
        nqumax=nqu
      endif
      return
 9000 call termes(lfno,'?Too many conditions*variables.',' ')
      if(nqumax .gt. 0)then
        ifqu=itmmapp(nqumax)
        ifquw=ktaloc(nqumax)
      endif
      return
      end

      subroutine tffssetlimit(nvar,dlim)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ii,i,nvar
      real*8 vl,vl1,vl2
      logical*4 dlim,limited,limited1
      limited=.false.
      do ii=1,nvar
        i=nvevx(ii)%ivarele
        call tffsvlimit(i,idelc(nelvx(i)%klp),
     $       nvevx(ii)%valvar,nvevx(ii)%valvar,
     $       vl,vl1,vl2,nvevx(ii)%ivvar,limited1,dlim)
        if(limited1)then
          limited=.true.
          nvevx(ii)%valvar=vl
        endif
      enddo
      if(limited)then
        call tfsetv(nvar)
      endif
      return
      end

      subroutine tffsvlimit(i,ld,val,val0,vl,
     $     vl1,vl2,ivv,limited,dlim)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_rlist), pointer :: klr
      type (sad_descriptor) kx
      integer*4 ,intent(in):: i,ld
      integer*4 ,intent(out):: ivv
      integer*4 ltyp,irtc
      real*8 ,intent(in):: val,val0
      real*8 ,intent(out):: vl,vl1,vl2
      real*8 vl0
      logical*4 ,intent(out):: limited,dlim
      limited=.false.
      ltyp=idtype(ld)
      vl=val
      vl0=val0
      if(ivv .eq. nelvx(i)%ival)then
        vl1=nelvx(i)%vlim(1)
        vl2=nelvx(i)%vlim(2)
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
      if(ltyp .eq. icMARK)then
        if(ivv .eq. mfitbx .or.ivv .eq. mfitby)then
          if(vl .le. 1.d-9)then
            vl=0.d0
            limited=.true.
          endif
        endif
      endif
      call tffsvarfun(1,ld,ivv,val,kx,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ktfrealq(kx))then
        dlim=.true.
        if(kx%k .eq. 0)then
          vl=vl0
          go to 2009
        endif
      elseif(tfreallistq(kx,klr))then
        if(klr%nl .eq. 2)then
          dlim=.true.
          vl1=max(vl1,klr%rbody(1))
          vl2=min(vl2,klr%rbody(2))
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
      use iso_c_binding
      use efun
      implicit none
      type (sad_string), pointer, save :: svarn, skey
      type (sad_descriptor) , save ::ifvr,ifvw
      data ifvr%k/0/
      type (sad_descriptor) ,intent(out):: kx
      type (sad_symdef) ,pointer:: symd
      integer*4 ,intent(in):: id,ld,k
      integer*4 ,intent(out):: irtc
      integer*4 isp1,level,itfuplevel,ltyp,itfdownlevel
      real*8 ,intent(in):: x
      character*(MAXPNAME) vn,tfkwrd
      integer*8, save :: ifvvarn,ifvkey
      if(ifvr%k .eq. 0)then
        ifvr=dtfcopy1(kxsymbolz('`VariableRange',14))
        ifvw=dtfcopy1(kxsymbolz('`VariableWeight',15))
        ifvvarn=ktsalocb(0,'        ',MAXPNAME+16)
        ifvkey=ktsalocb(0,'        ',MAXPNAME)
        call loc_string(ifvvarn,svarn)
        call loc_string(ifvkey,skey)
      endif
      isp=isp+1
      isp1=isp
      level=itfuplevel()
      if(id .eq. 1)then
c        dtastk(isp1)=ifvr
        call tfsyeval(ifvr,dtastk(isp1),irtc)
      elseif(id .eq. 2)then
c        dtastk(isp1)=ifvw
        call tfsyeval(ifvw,dtastk(isp1),irtc)
      endif
      if(irtc .eq. 0)then
        if(.not. ktfsymbolqdef(ktastk(isp1),symd) .or.
     $       symd%sym%override .eq. 0 .or. symd%downval .eq. 0)then
          isp=isp1-1
          kx=dxnull
          irtc=-1
          level=itfdownlevel()
          return
        endif
        ltyp=idtype(ld)
        svarn%nch=lpname(ld)
        svarn%str(1:svarn%nch+1)=pname(ld)(1:svarn%nch+1)//char(0)
        vn=tfkwrd(ltyp,k)
        skey%nch=len_trim(vn)
        skey%str(1:skey%nch+1)=vn(1:skey%nch+1)//char(0)
        isp=isp+1
        ktastk(isp)=ktfstring+ifvvarn
        isp=isp+1
        ktastk(isp)=ktfstring+ifvkey
        isp=isp+1
        rtastk(isp)=x
        call tclrfpe
c        write(*,'(a,1x,a,1x,a)')'vlimit:',svarn%str(1:svarn%nch),
c     $       vn(1:skey%nch)
        kx=tfefunref(isp1,.false.,irtc)
      endif
      if(irtc .ne. 0)then
        kx=dxnull
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
        call tfconnect(kx,irtc)
      endif
      isp=isp1-1
      return
      end

      subroutine tffsqu(nqcol,nqcol1,nvar,nqumax,ifquw,ifqu,
     $     free,nlat,nele,nfam,nfam1,nut,
     $     nparallel,cell,lfno,irtc)
      use tfstk
      use ffs, only:flv,ffs_bound,nvevx,nelvx
      use ffs_pointer
      use tffitcode
      use tfshare
      use mackw
      implicit none
      type (ffs_bound) fbound
      integer*8 ,intent(inout):: ifquw,ifqu
      integer*8 itmmapp,kcm,kkqu,kqu,ic,iec
      integer*4 ,intent(in):: nqcol,nqcol1,nvar,nlat,nele,
     $     lfno,nut,nfam,nfam1,nparallel
      integer*4 ,intent(out):: irtc,nqumax
      integer*4 npp
      logical*4 ,intent(in):: free(nele),cell
      integer*4 nqu,k,kk,i,kq,j,kf,lp,kp,iv,kkf,kkq,kkk,
     $     ii,ltyp,jj,kc,ik1,iclast(-nfam:nfam),ik,nk,kk1,
     $     ip,ipr,istep,npr(nparallel),fork_worker
      real*8 s,dtwiss(mfittry),coup,posk,wk,ctrans(4,7,-nfam:nfam)
      logical*4 col(2,nqcol),disp,nzcod
      integer*4 , parameter :: minnqu=512
      call tffscoupmatrix(kcm,lfno)
      irtc=0
      nqu=max(minnqu,nqcol*nvar)
      do9000: do kkk=1,1
        if(nqu .gt. nqumax)then
          if(nqumax .gt. 0)then
            call tfree(ifquw)
            call tmunmapp(ifqu)
          endif
          ifqu=itmmapp(nqu)
          if(ifqu .le. 0)then
            exit
          endif
          ifquw=ktaloc(nqu)
          if(ifquw .le. 0)then
            call tmunmapp(ifqu)
            exit
          endif
          nqumax=nqu
        endif
        call tffsbound(fbound)
        rlist(ifqu:ifqu+nqu-1)=0.d0
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
          kc=icomp(k)
          kk=iele1(kc)
          kk1=iele1(k)
          iec=kele2(k)
          if(iec .eq. 0)then
            nk=0
          else
            nk=ilist(1,iec)
          endif
          if(kk .gt. 0 .and. free(kk) .or. iec .ne. 0)then
            posk=pos(k)
            wk=1.d0
            if(k .eq. fbound%lb)then
              wk=1.d0-fbound%fb
            endif
            if(k .eq. fbound%le)then
              wk=fbound%fe
            endif
            ltyp=idtypec(k)
            do ii=ip,nvar,istep
              if(kk .le. 0 .or.
     $             (.not. free(kk) .and. .not. free(kk1)))then
                ik1=1
              else
                ik1=0
                do ik=1,nk
                  if(nvevx(ii)%ivvar .eq. ilist(2,iec+ik))then
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
                  iv=nvevx(ii)%ivvar
                  if(iv .eq. nelvx(kk)%ival .and.
     $                 nvevx(ii)%ivarele .eq. kk .and.
     $                 (nvevx(ii)%ivcomp .eq. 0 .or.
     $                 nvevx(ii)%ivcomp .eq. kc))then
                    coup=couple(k)*wk
                  elseif(iv .ne. nelvx(kk)%ival .and.
     $                   nvevx(ii)%ivarele .eq. kk1 .and.
     $                   (nvevx(ii)%ivcomp .eq. 0 .or.
     $                   nvevx(ii)%ivcomp .eq. k))then
                    coup=wk
                  else
                    cycle
                  endif
                endif
                iclast(nfam1:nfam)=0
                do concurrent (kq=1:nqcol1)
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
                          call qdcell(dtwiss,
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
                      else
                        call tdgeo(s,rlist(kqu),
     $                       kf,lp,k,ltyp,iv,nut,nfam)
                      endif
                    endif
                    s=-s
                  enddo
                enddo
              enddo do11
            enddo
          endif
        enddo
        call tffswait(ipr,npp,npr,i00,'tffsqu',irtc)
        if(kcm .ne. 0)then
          call tfree(kcm)
        endif
        return
      enddo do9000
      call termes(lfno,'?Too many conditions*variables.',' ')
      if(nqumax .gt. 0)then
        ifqu=itmmapp(nqumax)
        ifquw=ktaloc(nqumax)
      endif
      irtc=20002
      if(kcm .ne. 0)then
        call tfree(kcm)
      endif
      return
      end

      real*8 pure function tffsfmin(f1,f2,g1,g2,g0,dg)
      implicit none
      real*8 ,intent(in):: f1,f2,g1,g2,g0,dg
      real*8 a,b
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

      subroutine twmov(l,twiss,n1,n2,right)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer,only:direlc,compelc
      implicit none
      type (sad_comp), pointer::cmp
      integer*4 n1,n2,ntfun,l
      real*8 twiss(n1,-n2:n2,1:ntwissfun)
      logical*4 right
c
      call compelc(l,cmp)
      if(right)then
        if(orbitcal .or. calc6d)then
          ntfun=ntwissfun
        else
          ntfun=mfitdetr
        endif
        twiss(1,0,1:ntfun)=cmp%value(1:ntfun)
        if(direlc(l) .lt. 0.d0)then
          twiss(1,0,mfitax)=-cmp%value(mfitax)
          twiss(1,0,mfitay)=-cmp%value(mfitay)
          twiss(1,0,mfitaz)=-cmp%value(mfitaz)
          twiss(1,0,mfitepx)=-cmp%value(mfitepx)
          twiss(1,0,mfitepy)=-cmp%value(mfitepy)
          twiss(1,0,mfitzpx)=-cmp%value(mfitzpx)
          twiss(1,0,mfitzpy)=-cmp%value(mfitzpy)
          twiss(1,0,mfitr2)=-cmp%value(mfitr2)
          twiss(1,0,mfitr3)=-cmp%value(mfitr3)
          if(orbitcal)then
            twiss(1,0,mfitdpx)=-cmp%value(mfitdpx)
            twiss(1,0,mfitdpy)=-cmp%value(mfitdpy)
          endif
        endif
      else
        cmp%value(1:ntwissfun)=twiss(1,0,1:ntwissfun)
        if(direlc(l) .lt. 0.d0)then
          cmp%value(mfitax)=-twiss(1,0,mfitax)
          cmp%value(mfitay)=-twiss(1,0,mfitay)
          cmp%value(mfitaz)=-twiss(1,0,mfitaz)
          cmp%value(mfitepx)=-twiss(1,0,mfitepx)
          cmp%value(mfitepy)=-twiss(1,0,mfitepy)
          cmp%value(mfitzpx)=-twiss(1,0,mfitzpx)
          cmp%value(mfitzpy)=-twiss(1,0,mfitzpy)
          cmp%value(mfitr2)=-twiss(1,0,mfitr2)
          cmp%value(mfitr3)=-twiss(1,0,mfitr3)
          cmp%value(mfitdpx)=-twiss(1,0,mfitdpx)
          cmp%value(mfitdpy)=-twiss(1,0,mfitdpy)
        endif
      endif
      return
      end

      subroutine tgrad(qu,quw,df,grad,wlimit,wexponent,nqcol,nvar)
      implicit none
      integer*4 ,intent(in):: nqcol,nvar
      integer*4 i
      real*8 ,intent(out):: quw(nqcol,nvar),grad(nvar)
      real*8 ,intent(in):: qu(nqcol,nvar),df(nqcol),
     $     wlimit(nvar),wexponent
      real*8 dfw(nqcol),dfwi,sg,r
      do concurrent (i=1:nvar)
        quw(:,i)=qu(:,i)*wlimit(i)
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
      do concurrent (i=1:nvar)
        grad(i)=dot_product(quw(:,i),dfw)
      enddo
      sg=sum(grad**2)
      if(sg .ne. 0.d0)then
        grad=grad*r/sg
      endif
      return
      end

      real*8 function tweigh(i,ltyp,iv,val0,vk,absweit)
      use kyparam
      use tfstk
      use ffs, only:dpmax,emx,emy,brho,emminv
      use ffs_fit
      use cbkmac
      implicit none
      type (sad_descriptor) kx
      integer*4 i,ltyp,iv,irtc
      real*8 val0,gw,vmin,vk
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
        if(iv .ge. ky_K1_MULT)then
          gw=(max(1.d-3,dpmax)*max(1.d-2,etamax))
     $         **((iv-ky_K1_MULT)/2+1)
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
        if(iv .eq. ky_AX_MARK .or.
     $       iv .eq. ky_AY_MARK .or.
     $       iv .eq. ky_R1_MARK .or.
     $       iv .eq. ky_R4_MARK)then
          gw=1.d0/avebeta
        elseif(iv .eq. ky_BX_MARK .or.
     $         iv .eq. ky_BY_MARK)then
          gw=1.d0/avebeta**2
        elseif(iv .eq. ky_EX_MARK .or.
     $         iv .eq. ky_EPX_MARK)then
          gw=max(1.d-3,dpmax)*sqrt(avebeta/
     $         max(emminv,abs(emx)+abs(emy)))
     $         /avebeta**2
        elseif(iv .eq. ky_EPX_MARK .or.
     $         iv .eq. ky_EPX_MARK)then
          gw=max(1.d-3,dpmax)*sqrt(avebeta/
     $         max(emminv,abs(emx)+abs(emy)))
     $         /avebeta
        elseif(iv .eq. ky_R2_MARK)then
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
          if(iv .ge. ky_K1_MULT)then
            vmin=10.d0**(((iv-ky_K1_MULT)/2)*2-5)
          endif
        endif
        gw=sqrt(max(vmin,abs(val0))/gw)
      endif
      call tffsvarfun(2,i,iv,gw,kx,irtc)
      if(irtc .eq. 0 .and. ktfrealq(kx))then
        gw=rfromd(kx)
      endif
      tweigh=gw
      return
      end

      subroutine tfsolv(qu,quw,df,dval,wlimit,nqcol,nvar,
     $     iqcol,kfitp,mfitp,dg,wexponent,eps)
      implicit none
      integer*4 nqcol,nvar
      real*8 qu(nqcol,nvar),quw(nqcol,nvar)
      real*8 df(nqcol),dval(nvar)
      integer*4 iqcol(*),kfitp(*),mfitp(*)
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
          quw(nj,1:nvar)=qu(i,1:nvar)*wlimit(1:nvar)
          b(nj)=df(i)
        endif
      enddo
      call tsolva(quw,b,dval,nj,nvar,nqcol,eps)
      again=.false.
      dg=0.d0
      do i=1,nqcol
        s=sum(qu(i,:)*wlimit*dval)
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
      use tmacro
      implicit none
      integer*4 n,irtc
      if(nparallel .gt. 1)then
        irtc=1
        itmmapp=ktfallocshared(n)
      else
        itmmapp=ktaloc(n)
      endif
      return
      end

      subroutine tmunmapp(i)
      use tfstk
      use tfshare
      use tmacro
      implicit none
      integer*8 i
      if(nparallel .gt. 1)then
        call tfreeshared(i)
      else
        call tfree(i)
      endif
      return
      end

      subroutine tffscoupmatrix(kcm,lfno)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer,only:kele2
      implicit none
      type (sad_descriptor) km
      integer*8 k1,kam,kcm,ktfmaloc,k2
      integer*4 lfno,irtc,n,m
      real*8 rfromk
      integer*8 itfcoupm
      data itfcoupm /0/
      if(kele2(nlat) .eq. 0)then
        kcm=0
        return
      endif
      if(itfcoupm .eq. 0)then
        itfcoupm=ktfsymbolz('CouplingMatrix',14)
      endif
      levele=levele+1
      call tfsyeval(dfromk(itfcoupm),km,irtc)
      call tfconnect(km,irtc)
      if(irtc .ne. 0)then
        go to 9010
      endif
      if(.not. tflistq(km))then
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
