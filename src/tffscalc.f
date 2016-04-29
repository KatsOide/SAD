      subroutine tffscalc(flv,twiss,pos,geo,gammab,utwiss,
     $     latt,mult,ivarele,ivcomp,
     $     iele,iele1,ivvar,ival,
     $     nvar,klp,itwissp,nut,
     $     nlat,nele,ndim,
     $     dp,tracex,tracey,hstab,vstab,nfr,nfam,nfam1,
     $     kdp,df,
     $     wfit,wiq,wsum,iqcol,lfp,
     $     nqcol,nqcol1,nfcol,nfc0,maxcond,nlist,ibegin,
     $     dfam,jfam,kfam,inicond,iuid,uini,
     $     wake,iwakeelm,kwaketbl,nwakep,
     $     r,rp,rstab,nstab,residual,wexponent,
     $     offmw,etamax,avebeta,emx,emy,dpmax,coumin,
     $     cell,zcal,fitflg,geomet,cellstab,wcal,parallel,
     $     chgini,orbitcal,intres,halfres,sumres,diffres,
     $     lout,error)
      use tfstk
      use ffslocal, only: ffslocalv
      use tffitcode
      use tfshare
      implicit none
      type (ffslocalv) flv
      include 'inc/MACMATH.inc'
c      include 'DEBUG.inc'
      integer*8 kx
      integer*4 nlat,nele,ndim,nfr,maxcond,ibegin,nqcol,
     $     lfno,nfam,nfam1,nut,lbegin,lend,irtc
      integer*4 i1,i2,i3,i,ii,j,nfcol,iter,nvar,kt,iq,l,maxf,
     $     nqcol1,ie,ie1,iv,nfc0,nstab,lout,icslfno
      integer*4 latt(2,nlat),mult(nlat),itwissp(nlat),
     $     iele(nlat),iele1(nlat),ivvar(nvar),ival(nele),
     $     ivarele(nvar),klp(nele),ivcomp(nvar),
     1     kdp(maxcond),iqcol(maxcond),lfp(2,maxcond),
     $     jfam(-nfam:nfam),kfam(-nfam:nfam),iuid(-nfam:nfam)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),pos(nlat),geo(3,4,nlat),
     $     gammab(nlat),utwiss(ntwissfun,-nfam:nfam,nut),
     $     dp(-nfam:nfam),tracex(-nfam:nfam),tracey(-nfam:nfam),
     $     df(maxcond),wfit(maxcond),wiq(maxcond),
     $     r,rp,wexponent,dfam(4,-nfam:nfam),offmw,
     $     residual(-nfam:nfam),emx,emy,dpmax,coumin,wsum,wi,
     $     frbegin,frend,uini(27,-nfam:nfam)
      logical*4 cell,zcal,fitflg,geomet,cellstab,wcal,
     $     hstab(-nfam:nfam),vstab(-nfam:nfam),error,parallel
      character*8 nlist(mfit1)
      real*8 anux0,anuy0,anux0h,anuy0h,anuxi,anuyi,anuxih,anuyih,
     $     etamax,avebeta,rw,drw,rstab,
     $     anusumi,anusum0,anudiffi,anudiff0,physd(4)
      logical*4 fam,over(-nfam:nfam),beg,zerores,inicond,wake,
     $     chgini,orbitcal,intres,halfres,sumres,diffres
      integer*4 fork_worker,wait,irw,isw,ipr,ifb,ife,idir,
     $     jjfam(-nfam:nfam),ivoid,itfuplevel,itfdownlevel,
     $     nwakep,iwakeelm(nwakep),itgetfpe,ifpe,ntfun
      integer*8 iutm,jb,kwaketbl(2,nwakep)
      parameter (ivoid=9999)
      integer*8 iprolog,iepilog,imr,inr,isl
      data iprolog,iepilog,imr,inr,isl/0,0,0,0,0/
c     begin initialize for preventing compiler warning
      anux0=0.d0
      anuy0=0.d0
      anux0h=0.d0
      anuy0h=0.d0
      anusum0=0.d0
      anudiff0=0.d0
      iutm=0
      lfno=icslfno()
c     end   initialize for preventing compiler warning
      if(iprolog .eq. 0)then
        iprolog=ktfsymbolz('OpticsProlog',12)
        iepilog=ktfsymbolz('OpticsEpilog',12)
        imr=    ktfsymbolz('MatchingResidual',16)-4
        inr=    ktfsymbolz('NetResidual',11)-4
        isl=    ktfsymbolz('StabilityLevel',14)-4
      endif
      ifpe=itgetfpe()
      call tclrfpe
      l=itfuplevel()
      call tfeeval(ktfsymbol+iprolog,kx,.true.,irtc)
      l=itfdownlevel()
      if(irtc .ne. 0)then
        if(irtc .gt. 0)then
          if(ierrorprint .ne. 0)then
            call tfaddmessage(' ',0,lfno)
          endif
          call tclrfpe
          call termes(lfno,'Error in OpticsProlog.',' ')
        endif
        error=.true.
        return
      endif
      call tsetfpe(ifpe)
      if(zcal)then
        call tfgeo(latt,geo,pos,gammab,geomet .or. .not. fitflg)
      endif
      call tffsbound(nlat,latt,lbegin,frbegin,lend,frend)
      over(nfam1:nfam)=.false.
      hstab(nfam1:nfam)=.true.
      vstab(nfam1:nfam)=.true.
      jjfam(nfam1:nfam)=ivoid
      if(orbitcal)then
        ntfun=ntwissfun
      else
        ntfun=mfitdetr
      endif
      beg=ibegin .gt. lbegin
      if(beg)then
        twiss(ibegin,0,1:ntfun)
     $       =utwiss(1:ntfun,0,itwissp(ibegin))
        frbegin=0.d0
      else
        call twmov(latt(2,1),twiss,nlat,ndim,.true.)
        twiss(1,0,mfitnx)=0.d0
        twiss(1,0,mfitny)=0.d0
        if(orbitcal)then
          twiss(1,0,mfitddp)=twiss(1,0,mfitddp)+dp(0)
        endif
        if(lbegin .gt. 1)then
          twiss(lbegin,0,1:ntfun)=twiss(1,0,1:ntfun)
        endif
        ibegin=lbegin
      endif
      if(wake)then
        call tffswake(latt,twiss,gammab,utwiss,itwissp,uini,
     $       iwakeelm,kwaketbl,nwakep,
     $       ibegin,frbegin,lend,frend,
     $       nfam,nfam1,nut,hstab,vstab,tracex,tracey,over,beg)
      else
        call qcell1(latt,twiss,gammab,
     $       ibegin,frbegin,lend,frend,0,
     1       hstab(0),vstab(0),tracex(0),tracey(0),.false.,over(0),
     $       chgini,lout)
        call tffssetutwiss(latt,twiss,utwiss,itwissp,0,nlat,
     $       ndim,nfam,nut,lbegin,lend,frend,beg,.true.,.true.)
        if(cell)then
          anux0=aint(twiss(nlat,0,mfitnx)/pi2)
          anuy0=aint(twiss(nlat,0,mfitny)/pi2)
          anux0h=aint(twiss(nlat,0,mfitnx)/pi)
          anuy0h=aint(twiss(nlat,0,mfitny)/pi)
          anusum0=aint((twiss(nlat,0,mfitnx)+twiss(nlat,0,mfitny))/pi2)
          anudiff0=twiss(nlat,0,mfitnx)/pi2-
     $         aint(twiss(nlat,0,mfitnx)/pi2)-
     $         twiss(nlat,0,mfitny)/pi2+
     $         aint(twiss(nlat,0,mfitny)/pi2)
          if(hstab(0) .and. vstab(0) .or. chgini)then
            call twmov(latt(2,1),twiss,nlat,ndim,.false.)
          endif
        endif
        if(nfam .ne. 0)then
          if(parallel)then
            irtc=0
            iutm=ktfallocshared((2*nfam+1)*4)
c            iutm=mapalloc8(rlist(1),(2*nfam+1)*4,8,irtc)
            if(irtc .eq. 0)then
              ipr=fork_worker()
            else
              ipr=-1
            endif
          else
            ipr=-1
          endif
          if(ipr .gt. 0)then
            idir=-1
            ifb=-1
            ife=-nfr
          else
            idir=1
            ifb=1
            ife=nfr
          endif
 2        i1=0
          i2=0
          i3=0
          fam=.false.
 1        do ii=ifb,ife,idir
            if(fam)then
              if(kfam(ii) .eq. 0)then
                cycle
              elseif(ipr .gt. 0 .and. jfam(ii) .ge. 0 .or.
     $               ipr .eq. 0 .and. jfam(ii) .lt. 0)then
                cycle
              endif
              i3=i2
              i2=jfam(ii)
              i1=ii
            else
              i3=i2
              i2=i1
              i1=ii
            endif
            if(over(i2))then
              utwiss(1,ii,1:nut)=utwiss(1,i2,1:nut)
              utwiss(mfitddp,ii,1:nut)=utwiss(mfitddp,ii,1:nut)+dp(ii)
              over(ii)=.true.
              hstab(ii)=.false.
              vstab(ii)=.false.
              tracex(ii)=tracex(i2)
              tracey(ii)=tracey(i2)
            else
              if(beg)then
                twiss(ibegin,1,1:ntfun)
     $               =utwiss(1:ntfun,ii,itwissp(ibegin))
              else
                twiss(1,1,1:mfitdetr)=utwiss(1:mfitdetr,0,1)
                if(inicond)then
                  if(uini(1,ii) .ge. 0.d0)then
                    twiss(1,1,1:mfitdetr)=uini(1:mfitdetr,ii)
                  endif
                  twiss(1,1,mfitdx:mfitddp)=
     $                 utwiss(mfitdx:mfitddp,0,1)+
     $                 uini(mfitdx:mfitddp,ii)
                else
                  if(fam)then
                    twiss(1,1,mfitdx:mfitdpy )=
     $                   utwiss(mfitdx:mfitdpy ,i2,1)+dfam(1:4,ii)
                  else
                    if(dp(i2) .eq. dp(i3))then
                      call tgetphysdispu(utwiss(1,i2,1),physd)
                    else
                      physd=(utwiss(mfitdx:mfitdpy,i2,1)
     $                     -utwiss(mfitdx:mfitdpy,i3,1))/(dp(i2)-dp(i3))
                    endif
                    twiss(1,1,mfitdx:mfitdpy)=
     $                   utwiss(mfitdx:mfitdpy,i2,1)
     $                   +(dp(ii)-dp(i2))*physd
c                    if(ii .eq. nfr)then
c                      write(*,'(a,1p10g12.4)')'tffscalc ',
c     $                     twiss(1,1,mfitdx:mfitdpy),
c     $                     utwiss(mfitex:mfitepy,i2,1),dp(ii),dp(i2)
c                    endif
                  endif
                  twiss(1,1,mfitdz )=utwiss(mfitdz ,i2,1)
                  twiss(1,1,mfitddp)=twiss(1,0,mfitddp)+dp(ii)
                endif
                twiss(1,1,mfitnx)=0.d0
                twiss(1,1,mfitny)=0.d0
                if(lbegin .gt. 1)then
                  twiss(lbegin,1,1:ntfun)=twiss(1,1,1:ntfun)
                endif
              endif
              call qcell1(latt,twiss,gammab,
     $             ibegin,frbegin,lend,frend,
     $             1,hstab(ii),vstab(ii),tracex(ii),tracey(ii),
     $             fam,over(ii),chgini,lout)
              call tffssetutwiss(latt,twiss,utwiss,itwissp,ii,nlat,
     $             ndim,nfam,nut,lbegin,lend,frend,
     $             beg,.true.,.true.)
            endif
            if(cell)then
              anuxih=aint(twiss(nlat,1,mfitnx)/pi)
              anuyih=aint(twiss(nlat,1,mfitny)/pi)
              anuxi=aint(twiss(nlat,1,mfitnx)/pi2)
              anuyi=aint(twiss(nlat,1,mfitny)/pi2)
              anusumi=aint((twiss(nlat,1,mfitnx)
     $             +twiss(nlat,1,mfitny))/pi2)
              anudiffi=twiss(nlat,1,mfitnx)/pi2-
     $             aint(twiss(nlat,1,mfitnx)/pi2)-
     $             twiss(nlat,1,mfitny)/pi2+
     $             aint(twiss(nlat,1,mfitny)/pi2)
              hstab(ii)=hstab(ii) .and. (fam .or.
     $             (.not. intres .or. anuxi .eq. anux0) .and.
     $             (.not. halfres .or. anuxih .eq. anux0h) .and.
     $             (.not. sumres .or. anusumi .eq. anusum0) .and.
     $             (.not. diffres .or. anudiffi .eq. anudiff0))
              vstab(ii)=vstab(ii) .and. (fam .or.
     $             (.not. intres .or. anuyi .eq. anuy0) .and.
     $             (.not. halfres .or. anuyih .eq. anuy0h) .and.
     $             (.not. sumres .or. anusumi .eq. anusum0) .and.
     $             (.not. diffres .or. anudiffi .eq. anudiff0))
            endif
          enddo
          if(ipr .eq. -1 .and. .not. fam .and. idir .eq. 1)then
            idir=-1
            ifb=-1
            ife=-nfr
            go to 2
          endif
          if(.not. fam .and. nfam .gt. nfr)then
            fam=.true.
            ifb=nfam1
            ife=nfam
            idir=1
            go to 1
          endif
          if(ipr .gt. 0)then
            irw=0
            do while(irw .ne. ipr)
              irw=wait(isw)
            enddo
            if(isw .ne. 0)then
              call termes(lfno,
     1             '?Error in parallel process.',' ')
            endif
            do i=nfam1,nfam
              if(i .gt. 0 .and. i .le. nfr .or.
     $             kfam(i) .ne. 0 .and. jfam(i) .ge. 0)then
                jb=iutm+(i+nfam)*4
                call tmovi(rlist(jb+1),hstab(i),1)
                call tmovi(rlist(jb+2),vstab(i),1)
                tracex(i)=rlist(jb+3)
                tracey(i)=rlist(jb+4)
c                write(*,*)'tffscalc1 ',i,jb,tracex(i),rlist(jb+3)
              endif
            enddo
          elseif(ipr .eq. 0)then
            do i=nfam1,nfam
              if(i .gt. 0 .and. i .le. nfr .or.
     $             kfam(i) .ne. 0 .and. jfam(i) .ge. 0)then
                jb=iutm+(i+nfam)*4
                call tmovi(hstab(i),rlist(jb+1),1)
                call tmovi(vstab(i),rlist(jb+2),1)
                rlist(jb+3)=tracex(i)
                rlist(jb+4)=tracey(i)
c                write(*,*)'tffscalc2 ',i,jb,tracex(i),rlist(jb+3)
              endif
            enddo
            stop
          endif
          if(ipr .gt. 0)then
            call tfreeshared(iutm)
c            if(mapfree(rlist(iutm+1)) .ne. 0)then
c              call termes(lfno,
c     1             '?tffscalc-error in unmap.',' ')
c            endif
          endif
        endif
      endif
      call tdfun(flv,nfcol,iqcol,lfp,nqcol,nqcol1,
     1     kdp,nfr,nfam,nfam1,jfam,kfam,
     $     inicond,iuid,
     $     df,utwiss,itwissp,nut,
     1     tracex,tracey,pos,geo,
     $     latt,mult,
     $     dp,maxcond,nlist,error)
      if(error)then
        call termes(lfno,
     1         '?Too many fit conditions.',' ')
        return
      endif
      if(wcal)then
        iter=1
        zcal=.false.
        do i=1,nvar
          kt=idtype(latt(1,klp(ivarele(i))))
          if(kt .ne. 6 .and. kt .ne. 8 .and. kt .ne. 10
     $         .and. kt .ne. 12)then
            zcal=.true.
            exit
          endif
        enddo
        if(fitflg)then
          if(nvar .le. 0)then
            call termes(lfno,'?No variable.',' ')
            error=.true.
          endif
          if(.not. cell .and. .not. geomet)then
            do i=ibegin,nlat-1
              ii=iele(i)
              ie=iele1(ii)
              ie1=iele1(i)
              do j=1,nvar
                iv=ivvar(j)
                if(iv .eq. ival(ie) .and. ivarele(j) .eq. ie
     $               .and. (ivcomp(j) .eq. 0 .or.
     $               ivcomp(j) .eq. ii))then
                  ibegin=i
                  go to 1023
                elseif(iv .ne. ival(ie) .and. ivarele(j) .eq. ie1
     $                 .and. (ivcomp(j) .eq. 0 .or.
     $                 ivcomp(j) .eq. i))then
                  ibegin=i
                  go to 1023
                endif
                if(ivarele(j) .gt. ie)then
                  exit
                endif
              enddo
            enddo
            ibegin=nlat
          endif
        endif
1023    etamax=1.d-2
        if(cell)then
          maxf=lend
        else
          maxf=lbegin
          do i=1,nfcol
            maxf=max(maxf,
     $           flv%ifitp(flv%kfitp(i)),flv%ifitp1(flv%kfitp(i)))
          enddo
          if(pos(maxf) .le. pos(1))then
            maxf=lend
          endif
        endif
        do l=lbegin,maxf
          etamax=max(abs(twiss(l,0,mfitex)),
     $         abs(twiss(l,0,mfitey)),etamax)
        enddo
        avebeta=(pos(maxf)-pos(1))/
     $       max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
      endif
      call twfit(nfcol,nfc0,flv%kfit,wfit,wiq,
     1     flv%ifitp,flv%kfitp,kdp,nqcol,iqcol,twiss,pos,maxf,
     $     kfam,inicond,iuid,
     $     nlat,ndim,nfam,dp,emx,emy,dpmax,coumin,latt,mult,nlist,
     $     wcal)
      wcal=.false.
      rw=0.d0
      residual(nfam1:nfam)=0.d0
      wsum=0.d0
      do i=1,nqcol
        if(kdp(i) .ne. 0)then
          iq=iqcol(i)
          wi=(offmw/2.d0/
     $         sqrt(dble(max(1,abs(flv%mfitp(flv%kfitp(iq)))))))**2
          wsum=wsum+wi
          wiq(i)=wiq(i)*wi**(1.d0/wexponent)
        else
          wsum=wsum+1.d0
        endif
        df(i)=df(i)*wiq(i)
        if(df(i) .ne. 0.d0)then
          drw=min(1.d50,max(1.d-50,abs(df(i))))**wexponent
          rw=rw+drw
          residual(kdp(i))=residual(kdp(i))+drw
        endif
      enddo
      if(rw .gt. 0.d0)then
        r=wsum*(max(rw,1.d-50)/wsum)**(2.d0/wexponent)
      else
        r=0.d0
      endif
      rp=r
      nstab=0
      if(cell)then
        if(rstab .eq. 0.d0)then
          rstab=10.d0
          do while(rstab .lt. r)
            rstab=rstab*10.d0
          enddo
        endif
        cellstab=.true.
        zerores=.true.
        do i=nfam1,nfam
          if(residual(i) .ne. 0.d0)then
            zerores=.false.
            if(.not. hstab(i))then
              nstab=nstab+1
              residual(i)=residual(i)+rstab
              r=r+rstab
              cellstab=.false.
            endif
            if(.not. vstab(i))then
              nstab=nstab+1
              residual(i)=residual(i)+rstab
              r=r+rstab
              cellstab=.false.
            endif
          endif
        enddo
        if(zerores)then
          do i=nfam1,nfam
            cellstab=cellstab .and. hstab(i) .and. vstab(i)
          enddo
        endif
      else
        do i=nfam1,nfam
          if(.not. hstab(i))then
            nstab=nstab+2
            residual(i)=residual(i)+20.d0
          endif
        enddo
      endif
      rlist(imr)=r
      rlist(inr)=rp
      rlist(isl)=dble(nstab)
      ifpe=itgetfpe()
      call tclrfpe
      l=itfuplevel()
      call tfeeval(ktfsymbol+iepilog,kx,.true.,irtc)
      l=itfdownlevel()
      if(irtc .ne. 0)then
        if(irtc .gt. 0)then
          if(ierrorprint .ne. 0)then
            call tfaddmessage(' ',0,lfno)
          endif
          call tclrfpe
          call termes(lfno,'Error in OpticsEpilog.',' ')
        endif
        error=.true.
      endif
      call tsetfpe(ifpe)
      return
      end

      subroutine twfit(nfcol,nfc0,kfit,wfit,wiq,
     1     ifitp,kfitp,kdp,nqcol,iqcol,twiss,pos,maxf,
     $     kfam,inicond,iuid,
     $     nlat,ndim,nfam,dp,emx,emy,dpmax,coumin,latt,mult,nlist,
     $     wcal)
      use tfstk
      use tffitcode
      implicit none
c      include 'DEBUG.inc'
      integer*8 kx
      integer*4 nfcol,nlat,ndim,nfam,maxf,i,j,k,nqcol,iq
      integer*4 kfit(*),ifitp(*),kfitp(*),
     $     latt(2,nlat),mult(nlat),kdp(*),idp,iqcol(nqcol),
     $     kfam(-nfam:nfam),iuid(-nfam:nfam)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),pos(nlat),wfit(*),
     $     dp(-nfam:nfam),wiq(nqcol),
     $     emx,emy,dpmax,coum,emxx,emyy,dpm,coup,em,coumin,
     $rfromk
      integer*4 itfuplevel, level,irtc,nfc0
      character*16 name
      character*8 nlist(*)
      logical*4 wcal,inicond
      integer*8 ifv,ifvh,ifvloc,ifvfun,ifid
      save ifv,ifvh,ifvloc,ifvfun,ifid
      data ifv /0/
      if(ifv .eq. 0)then
        ifv=ktadaloc(0,4)
        ifvh=ktfsymbolz('FitWeight',9)
        ifvloc=ktsalocb(0,'                ',MAXPNAME+8)
        ifvfun=ktsalocb(0,'        ',MAXPNAME)
        ifid=ktraaloc(0,2)
        klist(ifv)=ktfsymbol+ktfcopy1(ifvh)
        klist(ifv+1)=ktfstring+ifvloc
        klist(ifv+2)=ktfstring+ifvfun
        klist(ifv+3)=ktflist+ifid
        klist(ifv+4)=0
      endif
      em=abs(emx)+abs(emy)
      coum=min(1.d0,
     $     max(coumin,0.01d0,1.d0/(abs(emx/emy)+abs(emy/emx))))
      coum=coum/(1.d0+coum)
      emxx=max(emx,coum*em)
      emyy=max(emy,coum*em)
      dpm=max(0.001d0,dpmax)
      coup=sqrt(abs(emxx/emyy)+abs(emyy/emxx))
      if(wcal)then
        do i=1,nfcol
          k=kfit(kfitp(i))
          j=ifitp(kfitp(i))
          if(k .eq. mfitex .or. k .eq. mfitpex)then
            wfit(i)=dpm/sqrt(emxx*twiss(j,0,mfitbx))
          elseif(k .eq. mfitepx .or. k .eq. mfitpepx)then
            wfit(i)=dpm*sqrt(twiss(j,0,mfitbx)
     1           /(1.d0+twiss(j,0,mfitax)**2)/emxx)
          elseif(k .eq. mfitey .or. k .eq. mfitpey)then
            wfit(i)=dpm/sqrt(emyy*twiss(j,0,mfitby))
          elseif(k .eq. mfitepy .or. k .eq. mfitpepy)then
            wfit(i)=dpm*sqrt(twiss(j,0,mfitby)
     1           /(1.d0+twiss(j,0,mfitay)**2)/emyy)
          elseif(k .eq. mfitr1)then
            wfit(i)=coup*sqrt(twiss(j,0,mfitbx)/twiss(j,0,mfitby))
          elseif(k .eq. mfitr2)then
            wfit(i)=coup/sqrt(twiss(j,0,mfitbx)*twiss(j,0,mfitby))
          elseif(k .eq. mfitr3)then
            wfit(i)=coup*sqrt(twiss(j,0,mfitbx)*twiss(j,0,mfitby))
          elseif(k .eq. mfitr4)then
            wfit(i)=coup*sqrt(twiss(j,0,mfitby)/twiss(j,0,mfitbx))
          elseif(k .eq. mfitdx)then
            wfit(i)=1.d0/sqrt(twiss(j,0,mfitbx)*em)
          elseif(k .eq. mfitdpx)then
            wfit(i)=sqrt(twiss(j,0,mfitbx)/
     $           (1.d0+twiss(j,0,mfitax)**2)/em)
          elseif(k .eq. mfitdy)then
            wfit(i)=1.d0/sqrt(twiss(j,0,mfitby)*em)
          elseif(k .eq. mfitdpy)then
            wfit(i)=sqrt(twiss(j,0,mfitby)/
     $           (1.d0+twiss(j,0,mfitay)**2)/em)
          elseif(k .eq. mfitleng)then
            wfit(i)=1.d0/(pos(maxf)-pos(1))*
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
          elseif(k .eq. mfitdz)then
            wfit(i)=0.01d0*sqrt(
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
     $           /em/(pos(maxf)-pos(1)))
          elseif(k .ge. mfitgx .and. k .le. mfitgz)then
            wfit(i)=0.01d0*sqrt(
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
     $           /em/(pos(maxf)-pos(1)))
          elseif(k .ge. mfitchi1 .and. k .le. mfitchi3)then
            wfit(i)=1.d3
          else
            wfit(i)=1.d0
          endif
          if(kfitp(i) .gt. nfc0)then
            wfit(i)=wfit(i)
          endif
        enddo
      endif
      do iq=1,nqcol
        i=iqcol(iq)
        if(i .gt. 0)then
          k=kfit(kfitp(i))
          j=ifitp(kfitp(i))
          idp=kdp(iq)
          wiq(iq)=wfit(i)
          call elname(latt,j,mult,name)
          call tfpadstr(name,ifvloc+1,len_trim(name))
          ilist(1,ifvloc)=len_trim(name)
          call tfpadstr(nlist(k),ifvfun+1,len_trim(nlist(k)))
          ilist(1,ifvfun)=len_trim(nlist(k))
          if(inicond)then
            rlist(ifid+1)=dble(iuid(idp))
          else
            rlist(ifid+1)=dble(kfam(idp))
          endif
          rlist(ifid+2)=dp(idp)
          rlist(ifv+4)=wfit(i)
          call tclrfpe
          level=itfuplevel()
          call tfleval(klist(ifv-3),kx,.true.,irtc)
          call tfconnectk(kx,irtc)
          if(irtc .ne. 0)then
            if(ierrorprint .ne. 0)then
              call tfaddmessage(' ',2,6)
            endif
            call termes(6,'Error in FitWeight '//
     $           nlist(k)//' at '//name,' ')
          elseif(ktfrealq(kx))then
c            write(*,*)'twfit ',nlist(k),vx,wfit(i)
            wiq(iq)=rfromk(kx)
          endif
        else
          wiq(iq)=1.d0
        endif
      enddo
      return
      end

      subroutine tffsbound(nlat,latt,lbegin,frbegin,lend,frend)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nlat,latt(2,nlat),lbegin,lend
      real*8 frbegin,frend,xnlat,offset,tffsmarkoffset
      xnlat=nlat
      if(idtype(latt(1,1)) .eq. icMARK)then
        offset=max(1.d0,min(xnlat,1.d0+tffsmarkoffset(latt(1,1))))
        lbegin=int(offset)
        frbegin=offset-lbegin
      else
        lbegin=1
        frbegin=0.d0
      endif
      if(idtype(latt(1,nlat-1)) .eq. icMARK)then
        offset=tffsmarkoffset(latt(1,nlat-1))
        if(offset .ne. 0.d0)then
          offset=max(dble(lbegin+1),min(xnlat,xnlat-1.d0+offset))
          lend=int(offset)
          frend=offset-lend
        else
          lend=nlat
          frend=0.d0
        endif
      else
        lend=nlat
        frend=0.d0
      endif
      return
      end

      subroutine tffsbound1(nlat,latt,la,lb,
     $     lbegin,frbegin,lend,frend)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nlat,latt(2,nlat),lbegin,lend,la,lb,lb1
      real*8 frbegin,frend,xnlat,offset,tffsmarkoffset
      xnlat=nlat
      if(idtype(latt(1,la)) .eq. icMARK)then
        offset=max(1.d0,min(xnlat,la+tffsmarkoffset(latt(1,la))))
        lbegin=int(offset)
        frbegin=offset-lbegin
      else
        lbegin=la
        frbegin=0.d0
      endif
      lb1=lb
      if(lb .eq. nlat .and. idtype(latt(1,nlat-1)) .eq. icMARK)then
        lb1=lb1-1
      endif
      if(idtype(latt(1,lb1)) .eq. icMARK)then
        offset=tffsmarkoffset(latt(1,lb1))
        if(offset .ne. 0.d0)then
          offset=max(dble(lbegin+1),min(xnlat,lb1+offset))
          lend=int(offset)
          frend=offset-lend
        else
          lend=lb1
          frend=0.d0
        endif
      else
        lend=lb1
        frend=0.d0
      endif
      return
      end

      real*8 function tffsmarkoffset(latt)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 latt(2)
      tffsmarkoffset=rlist(latt(2)+kytbl(kwOFFSET,icMARK))
      if(tffsmarkoffset .ne. 0.d0)then
        tffsmarkoffset=(tffsmarkoffset-.5d0)
     $       *rlist(latt(2)+ilist(1,latt(2)))+.5d0
      endif
      return
      end

      real*8 function tffselmoffset(latt,l,nlat)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nlat,latt(2,nlat),l,lm,nm,lx,nmmax,
     $     icslfno
      parameter (nmmax=256)
      real*8 offset,xp,xe,tffsmarkoffset
      nm=0
      xp=l
      if(idtype(latt(1,l)) .eq. icMARK)then
        xe=nlat
        lm=l
 8111   offset=tffsmarkoffset(latt(1,lm))
        if(offset .ne. 0.d0)then
          xp=offset+lm
          if(xp .ge. 1.d0 .and. xp .le. xe)then
            lx=int(xp)
            if(idtype(latt(1,lx)) .eq. icMARK)then
              nm=nm+1
              if(nm .lt. nmmax)then
                lm=lx
                go to 8111
              else
                call termes(icslfno(),
     $               '?Recursive OFFSET in',pname(latt(1,l)))
              endif
            endif
          endif
        endif
      endif
      tffselmoffset=xp
      return
      end

      subroutine tffssetutwiss(latt,twiss,utwiss,itwissp,idp,nlat,
     $     ndim,nfam,nut,lbegin,lend,frend,beg,begin,end)
      use tfstk
      use tffitcode
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nlat,ndim,nfam,nut,lbegin,lend,idp,jdp,
     $     latt(2,nlat),j,jp,le1,k,itwissp(nlat)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),
     $     utwiss(ntwissfun,-nfam:nfam,nut),frend
      logical*4 beg,begin,end
      jdp=min(1,abs(idp))
      do j=lbegin,lend
        jp=itwissp(j)
        if(jp .gt. 0)then
          do k=1,ntwissfun
            utwiss(k,idp,jp)=twiss(j,jdp,k)
          enddo
        endif
      enddo
      jp=itwissp(1)
      if(begin)then
        if(jp .gt. 0 .and. .not. beg)then
          do k=1,ntwissfun
            twiss(1,jdp,k)=twiss(lbegin,jdp,k)
            utwiss(k,idp,jp)=twiss(lbegin,jdp,k)
          enddo
        endif
      endif
      if(end)then
        if(idtype(latt(1,nlat-1)) .eq. icMARK)then
          if(frend .eq. 0.d0)then
            le1=lend
          else
            le1=lend+1
          endif
          do k=1,ntwissfun
            twiss(nlat-1,jdp,k)=twiss(le1,jdp,k)
            twiss(nlat,jdp,k)=twiss(le1,jdp,k)
          enddo
          jp=itwissp(nlat-1)
          if(jp .gt. 0)then
            do k=1,ntwissfun
              utwiss(k,idp,jp)=twiss(le1,jdp,k)
            enddo
          endif
          jp=itwissp(nlat)
          if(jp .gt. 0)then
            do k=1,ntwissfun
              utwiss(k,idp,jp)=twiss(le1,jdp,k)
            enddo
          endif
        endif
      endif
      return
      end
