      module calc

      contains
      subroutine tffscalc(kdp,df,iqcol,lfp,
     $     nqcola,nqcola1,ibegin,
     $     r,residual1,zcal,wcal,parallel,lout,error)
      use tfstk
      use maccode
      use ffs, only:ndim,nlat,flv,maxcond,ffs_bound,nvevx,nelvx,tsetintm
      use ffs_flag
      use ffs_pointer
      use ffs_fit
      use ffs_wake
      use tffitcode
      use cellm
      use dfun
      use tfshare
      use tfcsi,only:icslfnm
      use eeval
      use macmath
      use mathfun
      implicit none
c      include 'DEBUG.inc'
      type (sad_descriptor) kx
      type (ffs_bound) fbound,ibound
      integer*4 ,intent(out):: ibegin,iqcol(maxcond),lfp(2,maxcond),nqcola1,nqcola,kdp(maxcond)
      integer*4 i1,i2,i3,i,ii,j,iter,kt,iq,l,maxf,ie,ie1,iv,lout,irtc
      type (ffs_res) ,intent(out):: r,residual1(-ndimmax:ndimmax)
      real*8 ,intent(out):: df(maxcond)
      logical*4 ,intent(in):: parallel
      logical*4 ,intent(inout):: wcal,zcal
      logical*4 ,intent(out):: error
      real*8 anux0,anuy0,anux0h,anuy0h,anuxi,anuyi,anuxih,anuyih,
     $     rw,drw,wi,anusumi,anusum0,anudiffi,anudiff0
      logical*4 fam,beg,zerores
      integer*4 irw,isw,ipr,ifb,ife,idir,jjfam(-nfam:nfam),ifpe,ntfun
      integer*4, external :: waitpid,itgetfpe
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
c     end   initialize for preventing compiler warning
      if(iprolog == 0)then
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
      kx=tfeeval(dfromk(ktfsymbol+iprolog),.true.,irtc)
      l=itfdownlevel()
      if(irtc /= 0)then
        if(irtc > 0)then
          if(ierrorprint /= 0)then
            call tfaddmessage(' ',0,icslfnm())
          endif
          call tclrfpe
          call termes(icslfnm(),'Error in OpticsProlog.',' ')
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
      ntfun=merge(ntwissfun,mfitdetr,orbitcal .or. calc6d)
      beg=ibegin > fbound%lb
      if(beg)then
        twiss(ibegin,0,1:ntfun)=utwiss(1:ntfun,0,itwissp(ibegin))
        fbound%fb=0.d0
      else
        call twmov(1,twiss,nlat,ndim,.true.)
        twiss(1,0,mfitnx)=0.d0
        twiss(1,0,mfitny)=0.d0
        if(orbitcal)then
          twiss(1,0,mfitddp)=twiss(1,0,mfitddp)+dp(0)
        endif
        if(fbound%lb > 1)then
          twiss(fbound%lb,0,1:ntfun)=twiss(1,0,1:ntfun)
        endif
        ibegin=fbound%lb
      endif
      ibound=fbound
      ibound%lb=ibegin
      if(wake)then
        call tffswake(ibound,beg)
      else
        call qcell1(ibound,0,optstat(0),.false.,chgini,lout)
        call tffssetutwiss(0,fbound,beg,.true.,.true.)
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
        if(nfam /= 0)then
          if(parallel)then
            irtc=0
            iutm=ktfallocshared((2*nfam+1)*4)
c            iutm=mapalloc8(rlist(1),(2*nfam+1)*4,8,irtc)
            if(irtc == 0)then
              ipr=itffork()
            else
              ipr=-1
            endif
          else
            ipr=-1
          endif
          if(ipr > 0)then
            idir=-1
            ifb=-1
            ife=-nfr
          else
            idir=1
            ifb=1
            ife=nfr
            call tsetintm(-1.d0)
          endif
 2        i1=0
          i2=0
          i3=0
          fam=.false.
 1        do ii=ifb,ife,idir
            if(fam)then
              if(kfam(ii) == 0)then
                cycle
              elseif(ipr > 0 .and. jfam(ii) >= 0 .or.
     $               ipr == 0 .and. jfam(ii) .lt. 0)then
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
              optstat(ii)=ffs_stat(optstat(i2)%tracex,optstat(i2)%tracey,optstat(i2)%tracez,
     $             .false.,.false.,.false.,.true.)
            else
              if(beg)then
                twiss(ibegin,1,1:ntfun)
     $               =utwiss(1:ntfun,ii,itwissp(ibegin))
              else
                twiss(1,1,:)=utwiss(:,0,1)
                if(inicond)then
                  if(uini(mfitbx,ii) > 0.d0)then
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
c                      if(abs(physd(i)) > abs(physd1(i)))then
c                        physd(i)=physd1(i)
c                      endif
c                    enddo
                    twiss(1,1,mfitdx:mfitdpy)=
     $                   utwiss(mfitdx:mfitdpy,i2,1)
c     $                   +(dp(ii)-dp(i2))*physd
c                    if(ii == nfr)then
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
                if(fbound%lb > 1)then
                  twiss(fbound%lb,1,1:ntfun)=twiss(1,1,1:ntfun)
                endif
              endif
              call qcell1(ibound,1,optstat(ii),fam,chgini,lout)
              call tffssetutwiss(ii,fbound,beg,.true.,.true.)
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
              optstat(ii)%stabx=fam .or. optstat(ii)%stabx .and. (
     $             (.not. intres .or. anuxi == anux0) .and.
     $             (.not. halfres .or. anuxih == anux0h) .and.
     $             (.not. sumres .or. anusumi == anusum0) .and.
     $             (.not. diffres .or. anudiffi == anudiff0))
              optstat(ii)%staby=fam .or. optstat(ii)%staby .and. (
     $             (.not. intres .or. anuyi == anuy0) .and.
     $             (.not. halfres .or. anuyih == anuy0h) .and.
     $             (.not. sumres .or. anusumi == anusum0) .and.
     $             (.not. diffres .or. anudiffi == anudiff0))
            endif
          enddo
          if(ipr == -1 .and. .not. fam .and. idir == 1)then
            idir=-1
            ifb=-1
            ife=-nfr
            go to 2
          endif
          if(.not. fam .and. nfam > nfr)then
            fam=.true.
            ifb=nfam1
            ife=nfam
            idir=1
            go to 1
          endif
          if(ipr > 0)then
            irw=0
            do while(irw /= ipr)
              irw=waitpid(-1,isw)
            enddo
            if(isw /= 0)then
              call termes(icslfnm(),
     1             '?Error in parallel process.',' ')
            endif
            do i=nfam1,nfam
              if(i > 0 .and. i <= nfr .or.
     $             kfam(i) /= 0 .and. jfam(i) >= 0)then
                jb=iutm+(i+nfam)*4
                optstat(i)%stabx=ilist(1,jb+1) /= 0
                optstat(i)%staby=ilist(1,jb+2) /= 0
                optstat(i)%tracex=rlist(jb+3)
                optstat(i)%tracey=rlist(jb+4)
              endif
            enddo
          elseif(ipr == 0)then
            do i=nfam1,nfam
              if(i > 0 .and. i <= nfr .or.
     $             kfam(i) /= 0 .and. jfam(i) >= 0)then
                jb=iutm+(i+nfam)*4
                ilist(1,jb+1)=merge(1,0,optstat(i)%stabx)
                ilist(1,jb+2)=merge(1,0,optstat(i)%staby)
                rlist(jb+3)=optstat(i)%tracex
                rlist(jb+4)=optstat(i)%tracey
              endif
            enddo
            stop
          endif
          if(ipr > 0)then
            call tfreeshared(iutm)
          endif
        endif
      endif
      call tdfun(iqcol,lfp,nqcola,nqcola1,kdp,df,error)
      if(error)then
        call termes(icslfnm(),
     1         '?Too many fit conditions.',' ')
        return
      endif
      if(wcal)then
        iter=1
        zcal=.false.
        do i=1,nvar
          kt=idtypec(nelvx(nvevx(i)%ivarele)%klp)
          if(kt /= icSEXT .and. kt /= icOCTU .and. kt /= icDECA
     $         .and. kt /= icDODECA)then
            zcal=.true.
            exit
          endif
        enddo
        if(fitflg)then
          if(nvar <= 0)then
            call termes(icslfnm(),'?No variable.',' ')
            error=.true.
          endif
          if(.not. cell .and. .not. geomet)then
            do i=ibegin,nlat-1
              ii=icomp(i)
              ie=iele1(ii)
              ie1=iele1(i)
              do j=1,nvar
                iv=nvevx(j)%ivvar
                if(iv == nelvx(ie)%ival .and. nvevx(j)%ivarele == ie
     $               .and. (nvevx(j)%ivcomp == 0 .or.
     $               nvevx(j)%ivcomp == ii))then
                  ibegin=i
                  go to 1023
                elseif(iv /= nelvx(ie)%ival .and.
     $                 nvevx(j)%ivarele == ie1
     $                 .and. (nvevx(j)%ivcomp == 0 .or.
     $                 nvevx(j)%ivcomp == i))then
                  ibegin=i
                  go to 1023
                endif
                if(nvevx(j)%ivarele > ie)then
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
          if(pos(maxf) <= pos(1))then
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
      residual1(nfam1:nfam)=ffs_res(0.d0,0)
      wsum=0.d0
      do i=1,nqcola
        if(kdp(i) /= 0)then
          iq=iqcol(i)
          wi=(offmw/2.d0/
     $         sqrt(dble(max(1,abs(flv%mfitp(flv%kfitp(iq)))))))**2
          wsum=wsum+wi
          wiq(i)=wiq(i)*wi**(1.d0/wexponent)
        else
          wsum=wsum+1.d0
        endif
        df(i)=df(i)*wiq(i)
        if(df(i) /= 0.d0)then
          drw=min(1.d50,max(1.d-50,abs(df(i))))**wexponent
          rw=rw+drw
          residual1(kdp(i))%r=residual1(kdp(i))%r+drw
        endif
      enddo
      r=ffs_res(merge(wsum*(max(rw,1.d-50)/wsum)**(2.d0/wexponent),0.d0,rw > 0.d0),0)
      if(cell)then
        cellstab=.true.
        zerores=.true.
        do i=nfam1,nfam
          if(residual1(i)%r /= 0.d0)then
            zerores=.false.
            if(.not. optstat(i)%stabx)then
              r%nstab=r%nstab+1
              residual1(i)%nstab=residual1(i)%nstab+1
              cellstab=.false.
            endif
            if(.not. optstat(i)%staby)then
              r%nstab=r%nstab+1
              residual1(i)%nstab=residual1(i)%nstab+1
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
            r%nstab=r%nstab+2
            residual1(i)%nstab=residual1(i)%nstab+2
          endif
        enddo
      endif
      rlist(imr)=res2r(r)
      rlist(inr)=r%r
      rlist(isl)=dble(r%nstab)
      ifpe=itgetfpe()
      call tclrfpe
      l=itfuplevel()
      kx=tfeeval(dfromk(ktfsymbol+iepilog),.true.,irtc)
      l=itfdownlevel()
      if(irtc /= 0)then
        if(irtc > 0)then
          if(ierrorprint /= 0)then
            call tfaddmessage(' ',0,icslfnm())
          endif
          call tclrfpe
          call termes(icslfnm(),'Error in OpticsEpilog.',' ')
        endif
        error=.true.
      endif
c      call tfevals('Print["PROF: ",LINE["PROFILE","Q1"]]',kxx,irtc)
      call tsetfpe(ifpe)
      return
      end associate
      end
      
      subroutine twfit(kfit,ifitp,kfitp,kdp,nqcola,iqcol,maxf,wcal)
      use tfstk
      use ffs, only:emx,emy,dpmax,coumin,emminv
      use ffs_pointer
      use ffs_fit
c      use ffs_flag, only:cell
      use tffitcode
      use eeval
      use tfcsi,only:icslfnm
      implicit none
      integer*4 ,intent(in):: nqcola,maxf,
     $     kfit(*),ifitp(*),kfitp(*),kdp(*),iqcol(nqcola)
      logical*4 ,intent(in):: wcal
      integer*4 i,j,k,iq
      real*8 coum,emxx,emyy,dpm,coup,em
      integer*4 level,irtc,idp
      character*16 name
      type (sad_descriptor) kx
      type (sad_descriptor) ,save::kfv
      data kfv%k /0/
      type (sad_dlist), pointer , save::klv
      type (sad_rlist), pointer , save::klid
      integer*8 , save:: ifvloc,ifvfun
      real*8 , parameter :: almin=1.d0
      if(kfv%k == 0)then
        kfv=kxadaloc(0,4,klv)
        klv%head=dtfcopy(kxsymbolz('`FitWeight',10))
        ifvloc=ktsalocb(0,'                ',MAXPNAME+8)
        ifvfun=ktsalocb(0,'        ',MAXPNAME)
        klv%body(1)=ktfstring+ifvloc
        klv%body(2)=ktfstring+ifvfun
        klv%dbody(3)=kxraaloc(0,2,klid)
        klv%body(4)=0
      endif
      em=max(emminv,abs(emx)+abs(emy))
      coum=min(1.d0,
     $     max(coumin,0.01d0,
     $     1.d0/(abs(emx/max(emminv,emy))+abs(emy/max(emminv,emx)))))
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
     $           /em/max(almin,pos(maxf)-pos(1)))
          case (mfitbmagx,mfitbmagy,mfitbmagz)
            wfit(i)=3.d0
          case (mfitgx,mfitgy,mfitgz)
            wfit(i)=0.01d0*sqrt(
     $           max(twiss(maxf,0,mfitnx),twiss(maxf,0,mfitny))
     $           /em/max(almin,pos(maxf)-pos(1)))
          case (mfitchi1,mfitchi2,mfitchi3)
            wfit(i)=1.d3
          case default
            wfit(i)=1.d0
          end select
c          if(cell .and.
c     $         kfitp(i) > nfc0 .and. kfitp(i) <= nfc0+4)then
c            wfit(i)=wfit(i)*0.7d0
c          endif
        enddo
      endif
      do iq=1,nqcola
        i=iqcol(iq)
        if(i > 0)then
          k=kfit(kfitp(i))
          j=ifitp(kfitp(i))
          idp=kdp(iq)
          wiq(iq)=wfit(i)
          call elname(j,name)
          call tfpadstr(name,ifvloc+1,len_trim(name))
          ilist(1,ifvloc)=len_trim(name)
          call tfpadstr(nlist(k),ifvfun+1,len_trim(nlist(k)))
          ilist(1,ifvfun)=len_trim(nlist(k))
          klid%rbody(1)=merge(dble(iuid(idp)),dble(kfam(idp)),inicond)
          klid%rbody(2)=dp(idp)
          klv%rbody(4)=wfit(i)
          call tclrfpe
          level=itfuplevel()
          kx=tfleval(klv,.true.,irtc)
          call tfconnect(kx,irtc)
          if(irtc /= 0)then
            if(ierrorprint /= 0)then
              call tfaddmessage(' ',2,icslfnm())
            endif
            call termes(6,'Error in FitWeight '//
     $           nlist(k)//' at '//name,' ')
          elseif(ktfrealq(kx,wiq(iq)))then
          endif
        else
          wiq(iq)=1.d0
        endif
      enddo
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
      integer*4 ,intent(in):: n1,n2,l
      integer*4 ntfun
      real*8 ,intent(out):: twiss(n1,-n2:n2,1:ntwissfun)
      logical*4 ,intent(in):: right
c
      call compelc(l,cmp)
      if(right)then
        ntfun=merge(ntwissfun,mfitdetr,
     $       orbitcal .or. calc6d)
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

      subroutine tffssetutwiss(idp,fbound,beg,begin,end)
      use tfstk
      use ffs, only:ffs_bound,nlat
      use ffs_pointer
      use tffitcode
      use mackw
      implicit none
      type (ffs_bound) ,intent(in):: fbound
      integer*4 ,intent(in):: idp
      integer*4 jdp,j,jp,le1
      logical*4 ,intent(in):: beg,begin,end
      jdp=min(1,abs(idp))
      do concurrent (j=fbound%lb:fbound%le)
        jp=itwissp(j)
        if(jp > 0)then
          utwiss(:,idp,jp)=twiss(j,jdp,:)
        endif
      enddo
      jp=itwissp(1)
      if(begin)then
        if(jp > 0 .and. .not. beg)then
          twiss(1,jdp,:)=twiss(fbound%lb,jdp,:)
          utwiss(:,idp,jp)=twiss(fbound%lb,jdp,:)
        endif
      endif
      if(end)then
        if(idtypec(nlat-1) == icMARK)then
          le1=merge(fbound%le,fbound%le+1,fbound%fe == 0.d0)
          twiss(nlat-1,jdp,:)=twiss(le1,jdp,:)
          twiss(nlat,jdp,:)=twiss(le1,jdp,:)
          jp=itwissp(nlat-1)
          if(jp > 0)then
            utwiss(:,idp,jp)=twiss(le1,jdp,:)
          endif
          jp=itwissp(nlat)
          if(jp > 0)then
            utwiss(:,idp,jp)=twiss(le1,jdp,:)
          endif
        endif
      endif
      return
      end

      end module
