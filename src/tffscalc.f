      subroutine tffscalc(kdp,df,iqcol,lfp,
     $     nqcola,nqcola1,ibegin,
     $     r,rp,rstab,nstab,residual1,
     $     zcal,wcal,parallel,lout,error)
      use tfstk
      use ffs, only:ndim,nlat,flv,maxcond,ffs_bound
      use ffs_flag
      use ffs_pointer
      use ffs_fit
      use ffs_wake
      use tffitcode
      use tfshare
      use tfcsi,only:icslfno
      use macmath
      use mathfun
      implicit none
c      include 'DEBUG.inc'
      type (ffs_bound) fbound,ibound
      integer*8 kx
      integer*4 ibegin,nqcola,lfno,irtc
      integer*4 i1,i2,i3,i,ii,j,iter,kt,iq,l,maxf,
     $     nqcola1,ie,ie1,iv,nstab,lout
      integer*4 kdp(maxcond),iqcol(maxcond),lfp(2,maxcond)
      real*8 df(maxcond),r,rp,wi,
     $     residual1(-ndimmax:ndimmax)
      logical*4 zcal,wcal,error,parallel
      real*8 anux0,anuy0,anux0h,anuy0h,anuxi,anuyi,anuxih,anuyih,
     $     rw,drw,rstab,
     $     anusumi,anusum0,anudiffi,anudiff0
      logical*4 fam,beg,zerores
      integer*4 irw,isw,ipr,ifb,ife,idir,
     $     jjfam(-nfam:nfam),ifpe,ntfun
      integer*4, external :: fork_worker,wait,itfdownlevel,itfuplevel,
     $     itgetfpe
      integer*8 iutm,jb
      integer*4 , parameter:: ivoid=9999
      integer*8 ,save :: iprolog=0,iepilog=0,imr=0,inr=0,isl=0
      associate (nvar=>flv%nvar)
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
c      call tfmemcheckprint('tffscalc-before-prolog',.true.,irtc)
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
        call tfgeo(latt,geomet .or. .not. fitflg)
      endif
      call tffsbound(fbound)
      optstat(nfam1:nfam)%over=.false.
      optstat(nfam1:nfam)%stabx=.true.
      optstat(nfam1:nfam)%staby=.true.
      optstat(nfam1:nfam)%stabz=.true.
      jjfam(nfam1:nfam)=ivoid
      if(orbitcal .or. calc6d)then
        ntfun=ntwissfun
      else
        ntfun=mfitdetr
      endif
      beg=ibegin .gt. fbound%lb
      if(beg)then
        twiss(ibegin,0,1:ntfun)
     $       =utwiss(1:ntfun,0,itwissp(ibegin))
        fbound%fb=0.d0
      else
        call twmov(1,twiss,nlat,ndim,.true.)
        twiss(1,0,mfitnx)=0.d0
        twiss(1,0,mfitny)=0.d0
        if(orbitcal)then
          twiss(1,0,mfitddp)=twiss(1,0,mfitddp)+dp(0)
        endif
        if(fbound%lb .gt. 1)then
          twiss(fbound%lb,0,1:ntfun)=twiss(1,0,1:ntfun)
        endif
        ibegin=fbound%lb
      endif
      ibound=fbound
      ibound%lb=ibegin
      if(wake)then
        call tffswake(ibound,beg)
      else
c        call tfevals('Print["PROF-1: ",LINE["PROFILE","Q1"]]',kxx,irtc)
        call qcell1(ibound,0,optstat(0),.false.,chgini,lout)
c        call tfevals('Print["PROF-2: ",LINE["PROFILE","Q1"]]',kxx,irtc)
        call tffssetutwiss(0,nlat,fbound,beg,.true.,.true.)
        if(cell)then
          anux0=aint(twiss(nlat,0,mfitnx)/pi2)
          anuy0=aint(twiss(nlat,0,mfitny)/pi2)
          anux0h=aint(twiss(nlat,0,mfitnx)/pi)
          anuy0h=aint(twiss(nlat,0,mfitny)/pi)
          anusum0=aint((twiss(nlat,0,mfitnx)+twiss(nlat,0,mfitny))/pi2)
          anudiff0=tfloor(
     $         (twiss(nlat,0,mfitnx)-twiss(nlat,0,mfitny))/pi2)
c          anudiff0=twiss(nlat,0,mfitnx)/pi2-
c     $         aint(twiss(nlat,0,mfitnx)/pi2)-
c     $         twiss(nlat,0,mfitny)/pi2+
c     $         aint(twiss(nlat,0,mfitny)/pi2)
          if(optstat(0)%stabx .and. optstat(0)%staby
     $         .and. optstat(0)%stabz .or. chgini)then
            call twmov(1,twiss,nlat,ndim,.false.)
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
            if(optstat(i2)%over)then
              utwiss(1,ii,1:nut)=utwiss(1,i2,1:nut)
              utwiss(mfitddp,ii,1:nut)=utwiss(mfitddp,ii,1:nut)+dp(ii)
              optstat(ii)%over=.true.
              optstat(ii)%stabx=.false.
              optstat(ii)%staby=.false.
              optstat(ii)%stabz=.false.
              optstat(ii)%tracex=optstat(i2)%tracex
              optstat(ii)%tracey=optstat(i2)%tracey
              optstat(ii)%tracez=optstat(i2)%tracez
            else
              if(beg)then
                twiss(ibegin,1,1:ntfun)
     $               =utwiss(1:ntfun,ii,itwissp(ibegin))
              else
                twiss(1,1,:)=utwiss(:,0,1)
                if(inicond)then
                  if(uini(mfitbx,ii) .gt. 0.d0)then
                    twiss(1,1,:)=uini(:,ii)
                  endif
                  twiss(1,1,mfitdx:mfitddp)=
     $                 utwiss(mfitdx:mfitddp,0,1)+
     $                 uini(mfitdx:mfitddp,ii)
                else
                  if(fam)then
                    twiss(1,1,mfitdx:mfitdpy )=
     $                   utwiss(mfitdx:mfitdpy ,i2,1)+dfam(1:4,ii)
                  else
c                    call tgetphysdispu(utwiss(1,i2,1),physd)
c                    physd1=(utwiss(mfitdx:mfitdpy,i2,1)
c     $                     -utwiss(mfitdx:mfitdpy,i3,1))/(dp(i2)-dp(i3))
c                    do i=1,4
c                      if(abs(physd(i)) .gt. abs(physd1(i)))then
c                        physd(i)=physd1(i)
c                      endif
c                    enddo
                    twiss(1,1,mfitdx:mfitdpy)=
     $                   utwiss(mfitdx:mfitdpy,i2,1)
c     $                   +(dp(ii)-dp(i2))*physd
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
                if(fbound%lb .gt. 1)then
                  twiss(fbound%lb,1,1:ntfun)=twiss(1,1,1:ntfun)
                endif
              endif
              call qcell1(ibound,1,optstat(ii),fam,chgini,lout)
              call tffssetutwiss(ii,nlat,fbound,beg,.true.,.true.)
            endif
            if(cell)then
              anuxih=aint(twiss(nlat,1,mfitnx)/pi)
              anuyih=aint(twiss(nlat,1,mfitny)/pi)
              anuxi=aint(twiss(nlat,1,mfitnx)/pi2)
              anuyi=aint(twiss(nlat,1,mfitny)/pi2)
              anusumi=aint((twiss(nlat,1,mfitnx)
     $             +twiss(nlat,1,mfitny))/pi2)
              anudiffi=tfloor((twiss(nlat,1,mfitnx)
     $             -twiss(nlat,1,mfitny))/pi2)
c              anudiffi=twiss(nlat,1,mfitnx)/pi2-
c     $             aint(twiss(nlat,1,mfitnx)/pi2)-
c     $             twiss(nlat,1,mfitny)/pi2+
c     $             aint(twiss(nlat,1,mfitny)/pi2)
c              write(*,*)'tffscalc ',anudiffi,anudiff0
              optstat(ii)%stabx=optstat(ii)%stabx .and. (fam .or.
     $             (.not. intres .or. anuxi .eq. anux0) .and.
     $             (.not. halfres .or. anuxih .eq. anux0h) .and.
     $             (.not. sumres .or. anusumi .eq. anusum0) .and.
     $             (.not. diffres .or. anudiffi .eq. anudiff0))
              optstat(ii)%staby=optstat(ii)%staby .and. (fam .or.
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
                call tmovi(rlist(jb+1),optstat(i)%stabx,1)
                call tmovi(rlist(jb+2),optstat(i)%staby,1)
                optstat(i)%tracex=rlist(jb+3)
                optstat(i)%tracey=rlist(jb+4)
c                write(*,*)'tffscalc1 ',i,jb,tracex(i),rlist(jb+3)
              endif
            enddo
          elseif(ipr .eq. 0)then
            do i=nfam1,nfam
              if(i .gt. 0 .and. i .le. nfr .or.
     $             kfam(i) .ne. 0 .and. jfam(i) .ge. 0)then
                jb=iutm+(i+nfam)*4
                call tmovi(optstat(i)%stabx,rlist(jb+1),1)
                call tmovi(optstat(i)%staby,rlist(jb+2),1)
                rlist(jb+3)=optstat(i)%tracex
                rlist(jb+4)=optstat(i)%tracey
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
c      call tfevals('Print["PROF-3: ",LINE["PROFILE","Q1"]]',kxx,irtc)
      call tdfun(iqcol,lfp,nqcola,nqcola1,kdp,df,error)
      if(error)then
        call termes(lfno,
     1         '?Too many fit conditions.',' ')
        return
      endif
      if(wcal)then
        iter=1
        zcal=.false.
        do i=1,nvar
          kt=idtypec(klp(ivarele(i)))
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
          maxf=fbound%le
        else
          maxf=fbound%lb
          do i=1,nfcol
            maxf=max(maxf,
     $           flv%ifitp(flv%kfitp(i)),flv%ifitp1(flv%kfitp(i)))
          enddo
          if(pos(maxf) .le. pos(1))then
            maxf=fbound%le
          endif
        endif
        do l=fbound%lb,maxf
          etamax=max(abs(twiss(l,0,mfitex)),
     $         abs(twiss(l,0,mfitey)),etamax)
        enddo
        avebeta=(pos(maxf)-pos(1))/
     $       max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
      endif
      call twfit(flv%kfit,
     1     flv%ifitp,flv%kfitp,kdp,nqcola,iqcol,maxf,wcal)
      wcal=.false.
      rw=0.d0
      residual1(nfam1:nfam)=0.d0
      wsum=0.d0
      do i=1,nqcola
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
          residual1(kdp(i))=residual1(kdp(i))+drw
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
          if(residual1(i) .ne. 0.d0)then
            zerores=.false.
            if(.not. optstat(i)%stabx)then
              nstab=nstab+1
              residual1(i)=residual1(i)+rstab
              r=r+rstab
              cellstab=.false.
            endif
            if(.not. optstat(i)%staby)then
              nstab=nstab+1
              residual1(i)=residual1(i)+rstab
              r=r+rstab
              cellstab=.false.
            endif
          endif
        enddo
        if(zerores)then
          do i=nfam1,nfam
            cellstab=cellstab .and.
     $           optstat(i)%stabx .and. optstat(i)%staby
          enddo
        endif
      else
        do i=nfam1,nfam
          if(.not. optstat(i)%stabx)then
            nstab=nstab+2
            residual1(i)=residual1(i)+20.d0
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
c      call tfevals('Print["PROF: ",LINE["PROFILE","Q1"]]',kxx,irtc)
      call tsetfpe(ifpe)
      return
      end associate
      end

      subroutine twfit(kfit,
     1     ifitp,kfitp,kdp,nqcola,iqcol,maxf,wcal)
      use tfstk
      use ffs, only:emx,emy,dpmax,coumin
      use ffs_pointer
      use ffs_fit
      use ffs_flag, only:cell
      use tffitcode
      implicit none
c      include 'DEBUG.inc'
      integer*8 kx
      integer*4 maxf,i,j,k,nqcola,iq
      integer*4 kfit(*),ifitp(*),kfitp(*),kdp(*),idp,iqcol(nqcola)
      real*8 coum,emxx,emyy,dpm,coup,em,rfromk
      integer*4 itfuplevel, level,irtc
      character*16 name
      logical*4 wcal
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
          select case(k)
          case (mfitex,mfitpex)
            wfit(i)=dpm/sqrt(emxx*twiss(j,0,mfitbx))
          case (mfitepx,mfitpepx)
            wfit(i)=dpm*sqrt(twiss(j,0,mfitbx)
     1           /(1.d0+twiss(j,0,mfitax)**2)/emxx)
          case (mfitey,mfitpey)
            wfit(i)=dpm/sqrt(emyy*twiss(j,0,mfitby))
          case (mfitepy,mfitpepy)
            wfit(i)=dpm*sqrt(twiss(j,0,mfitby)
     1           /(1.d0+twiss(j,0,mfitay)**2)/emyy)
          case (mfitr1)
            wfit(i)=coup*sqrt(twiss(j,0,mfitbx)/twiss(j,0,mfitby))
          case (mfitr2)
            wfit(i)=coup/sqrt(twiss(j,0,mfitbx)*twiss(j,0,mfitby))
          case (mfitr3)
            wfit(i)=coup*sqrt(twiss(j,0,mfitbx)*twiss(j,0,mfitby))
          case (mfitr4)
            wfit(i)=coup*sqrt(twiss(j,0,mfitby)/twiss(j,0,mfitbx))
          case (mfitdx)
            wfit(i)=1.d0/sqrt(twiss(j,0,mfitbx)*em)
          case (mfitdpx)
            wfit(i)=sqrt(twiss(j,0,mfitbx)/
     $           (1.d0+twiss(j,0,mfitax)**2)/em)
          case (mfitdy)
            wfit(i)=1.d0/sqrt(twiss(j,0,mfitby)*em)
          case (mfitdpy)
            wfit(i)=sqrt(twiss(j,0,mfitby)/
     $           (1.d0+twiss(j,0,mfitay)**2)/em)
          case (mfitleng)
            wfit(i)=1.d0/(pos(maxf)-pos(1))*
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
          case (mfitdz)
            wfit(i)=0.01d0*sqrt(
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
     $           /em/(pos(maxf)-pos(1)))
          case (mfitgx,mfitgy,mfitgz)
            wfit(i)=0.01d0*sqrt(
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
     $           /em/(pos(maxf)-pos(1)))
          case (mfitchi1,mfitchi2,mfitchi3)
            wfit(i)=1.d3
          case default
            wfit(i)=1.d0
          end select
          if(cell .and.
     $         kfitp(i) .gt. nfc0 .and. kfitp(i) .le. nfc0+4)then
            wfit(i)=wfit(i)*0.7d0
          endif
        enddo
      endif
      do iq=1,nqcola
        i=iqcol(iq)
        if(i .gt. 0)then
          k=kfit(kfitp(i))
          j=ifitp(kfitp(i))
          idp=kdp(iq)
          wiq(iq)=wfit(i)
          call elname(j,name)
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
c            write(*,*)'twfit ',nlist(k),rfromk(kx),wfit(i)
            wiq(iq)=rfromk(kx)
          endif
        else
          wiq(iq)=1.d0
        endif
      enddo
      return
      end

      subroutine tffsbound(fbound)
      use tfstk
      use ffs, only:nlat,ffs_bound
      use ffs_pointer
      use mackw
      implicit none
      type (ffs_bound) fbound
      call tffsbound1(1,nlat,fbound)
      return
      end

      subroutine tffsbound1(lb,le,fbound)
      use tfstk
      use ffs, only:nlat,ffs_bound
      use ffs_pointer
      use mackw
      implicit none
      type (ffs_bound) fbound
      integer*4 lb,le,le1
      real*8 xnlat,offset,tffsmarkoffset
      xnlat=nlat
      offset=tffsmarkoffset(lb)
      if(offset .ne. 0.d0)then
        offset=max(1.d0,min(xnlat,lb+offset))
        fbound%lb=int(offset)
        fbound%fb=offset-fbound%lb
      else
        fbound%lb=lb
        fbound%fb=0.d0
      endif
      le1=le
      if(le .eq. nlat)then
        if(idtypec(nlat-1) .eq. icMARK)then
          le1=le1-1
        else
          fbound%le=le1
          fbound%fe=0.d0
          return
        endif
      endif
      offset=tffsmarkoffset(le1)
      if(offset .ne. 0.d0)then
        offset=max(1.d0,min(xnlat,le1+offset))
        fbound%le=int(offset)
        fbound%fe=offset-fbound%le
      else
        fbound%le=le1
        fbound%fe=0.d0
      endif
      return
      end

      real*8 function tffsmarkoffset(lp)
      use tfstk
      use ffs_pointer
      use mackw
      use tmacro, only:nlat
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 lp,k
      if(lp .ge. nlat)then
        tffsmarkoffset=0.d0
        return
      endif
      call compelc(lp,cmp)
      k=kytbl(kwOFFSET,idtype(cmp%id))
      if(k .gt. 0)then
        tffsmarkoffset=cmp%value(k)
        if(tffsmarkoffset .ne. 0.d0)then
          tffsmarkoffset=(tffsmarkoffset-.5d0)*direlc(lp)+.5d0
        endif
      else
        tffsmarkoffset=0
      endif
      return
      end

      real*8 function tffselmoffset(l)
      use tfstk
      use ffs, only:nlat
      use ffs_pointer
      use mackw
      use tfcsi, only:icslfno
      implicit none
      integer*4 l,lm,nm,lx,nmmax
      parameter (nmmax=256)
      real*8 offset,xp,xe,tffsmarkoffset
      nm=0
      xp=l
      if(l .ne. nlat .and. kytbl(kwOFFSET,idtypec(l)) .ne. 0)then
        xe=nlat
        lm=l
 8111   offset=tffsmarkoffset(lm)
        if(offset .ne. 0.d0)then
          xp=offset+lm
          if(xp .ge. 1.d0 .and. xp .le. xe)then
            lx=int(xp)
            if(idtypec(lx) .eq. icMARK)then
              nm=nm+1
              if(nm .lt. nmmax)then
                lm=lx
                go to 8111
              else
                call termes(icslfno(),
     $               '?Recursive OFFSET in',pnamec(l))
              endif
            endif
          endif
        endif
      endif
      tffselmoffset=xp
      return
      end

      subroutine tffssetutwiss(idp,nlat,fbound,beg,begin,end)
      use tfstk
      use ffs, only:ffs_bound
      use ffs_pointer
      use tffitcode
      use mackw
      implicit none
      type (ffs_bound) fbound
      integer*4 nlat,idp,jdp,j,jp,le1
      logical*4 beg,begin,end
      jdp=min(1,abs(idp))
      do j=fbound%lb,fbound%le
        jp=itwissp(j)
        if(jp .gt. 0)then
c          do k=1,ntwissfun
            utwiss(:,idp,jp)=twiss(j,jdp,:)
c          enddo
        endif
      enddo
      jp=itwissp(1)
      if(begin)then
        if(jp .gt. 0 .and. .not. beg)then
c          do k=1,ntwissfun
            twiss(1,jdp,:)=twiss(fbound%lb,jdp,:)
            utwiss(:,idp,jp)=twiss(fbound%lb,jdp,:)
c          enddo
        endif
      endif
      if(end)then
        if(idtypec(nlat-1) .eq. icMARK)then
          if(fbound%fe .eq. 0.d0)then
            le1=fbound%le
          else
            le1=fbound%le+1
          endif
c          do k=1,ntwissfun
            twiss(nlat-1,jdp,:)=twiss(le1,jdp,:)
            twiss(nlat,jdp,:)=twiss(le1,jdp,:)
c          enddo
          jp=itwissp(nlat-1)
          if(jp .gt. 0)then
c            do k=1,ntwissfun
              utwiss(:,idp,jp)=twiss(le1,jdp,:)
c            enddo
          endif
          jp=itwissp(nlat)
          if(jp .gt. 0)then
c            do k=1,ntwissfun
              utwiss(:,idp,jp)=twiss(le1,jdp,:)
c            enddo
          endif
        endif
      endif
      return
      end
