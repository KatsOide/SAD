      subroutine temits(
     $     mphi2,amus0,amus1,amusstep,
     $     emix,emiy,res,params,
     $     lfno,kx,irtc)
      use tfstk
      use ffs_pointer
      use ffs
      use tmacro
      use tffitcode, only:ntwissfun
      use temw, only:nparams
      implicit none
      integer*8 kx
      integer*4 irtc,mphi2,lfno,mphiz,ndims,
     $     ndps,nzz
      real*8 amus0,amus1,amusstep,emix,emiy,res
      real*8 params(nparams)
      logical*4 trpt0,stab,radcod0,radtaper0,rfsw0,intra0
      integer*4 , parameter :: ndpmin=15
      radcod0=radcod
      radtaper0=radtaper
      rfsw0=rfsw
      trpt0=trpt
      intra0=intra
      radcod=.false.
      radtaper=.false.
      rfsw=.false.
      trpt=.false.
      intra=.false.
      mphiz=mphi2*2-1
      ndims=mphiz
      ndps=max(ndpmin,mphiz)
      nzz=mphiz*2+1
      call temits1(
     $     .false.,stab,ndps,nzz,mphiz,mphi2,ndims,
     $     amus0,amus1,amusstep,
     $     emix,emiy,res,params,
     $     lfno,kx,irtc)
      radcod=radcod0
      radtaper=radtaper0
      rfsw=rfsw0
      trpt=trpt0
      intra=intra0
      return
      end

      subroutine temits1(
     $     plot,stab,ndp,nz,mphi,mphi2,ndims,
     $     amus0,amus1,amusstep,
     $     emix,emiy,res0,params,
     $     lfno,kx,irtc)
      use tfstk
      use ffs
      use temw
      use tmacro
      implicit none
      type (sad_dlist), pointer :: kl
      type (sad_rlist), pointer :: kvl
      integer*4 , parameter:: itmax=20,ndiffh=16
      real*8 , parameter :: resconv=1.d-4,minconv=1.d-7
      integer*8 kx,kax,kai
      integer*4 mphi,mphi2,ndims,i1,nd1,
     $     lfno,irtc,ndp,i,k,kk,ns,it,j,nz
      real*8 amus0,amus1,amusstep, emix,emiy,res0,
     $     amus,damp,dj,dpndim,emix0,emiy0,
     $     fz,phi0s,res,sigea,vx,vy,w,
     $     emixp,emiyp,cod0(6)
      real*8 beam(42),trans(6,12),cod(6),srot(3,9),
     $     beams(10,-ndims:ndims),
     $     trads(5,5,-ndims:ndims),
     $     tws(ntwissfun,-ndims:ndims),
     $     bc(10,mphi2,ndp),bs(10,mphi2,ndp),
     $     ba(10,10,mphi,ndp),bb(10,mphi2,ndp),
     $     bd(10,4,mphi,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     ha(4,4,mphi,ndp),hb(4,mphi2,ndp),
     $     bff(20*mphi2,20*mphi2),
     $     bfx(20*mphi2),bfb(20*mphi2),
     $     hff(8*mphi2,8*mphi2),
     $     hfx(8*mphi2),hfb(8*mphi2),
c     $     dbc(10,mphi2,ndp),dbs(10,mphi2,ndp),
c     $     dhc(4,mphi2,ndp),dhs(4,mphi2,ndp),
     $     amuj(ndp)
      real*8 params(nparams),btr(21,21),rm(6,6),
     $     rx(6,6),rxi(6,6),cmu(mphi2),smu(mphi2),disppi(6),
     $     beamr(42),fj(256),conv,tw0(ntwissfun),dispp(6)
      logical*4 plot,calpol0,stab
      calpol0=calpol
      calpol=.false.
      codin=0.d0
      beamin(1:21)=0.d0
      cod=0
      call temit(trans,cod,beamr,btr,
     $     .true.,iaez,
     $     plot,params,stab,lfno)
      dispp=r(:,6)
      call tinitr(rx)
      rxi=rx
      rx(1:4,6)=dispp(1:4)
      rxi(1:4,6)=-dispp(1:4)
      dispp(5)=0.d0
      dispp(6)=1.d0
      sige=params(ipsige)
      damp=params(ipdampz)
      conv=max(resconv*abs(damp),minconv)
      sigea=sqrt(12.d0)*sige
      dpndim=sigea/ndims
      tw0=params(iptwiss:iptws0+mfitzpy)
      codin=cod
      tw0(ipnx)=0.d0
      tw0(ipny)=0.d0
      tw0(ipnz)=0.d0
      calint=intra
      disppi(5)=0.d0
      disppi(6)=1.d0
      i=0
      do i1=1,2*ndims+1
        irad=6
        if(i .eq. 0)then
          cod=codin
        elseif(i .gt. 0)then
          call tgetphysdispu(tws(1,i-1),disppi)
          cod=tws(mfitdx:mfitddp,i-1)+disppi*dpndim
        else
          call tgetphysdispu(tws(1,i+1),disppi)
          cod=tws(mfitdx:mfitddp,i+1)-disppi*dpndim
        endif
        cod(5)=0.d0
        cod0=cod
        call tcod(trans,cod,beam,.false.,fndcod)
        if(.not. fndcod .and. (i .gt. 1 .or. i .lt. -1))then
          if(i .gt. 1)then
            cod(1:4)=
     $           2.d0*tws(mfitdx:mfitdpy,i-1)-tws(mfitdx:mfitdpy,i-2)
          else
            cod(1:4)=
     $           2.d0*tws(mfitdx:mfitdpy,i+1)-tws(mfitdx:mfitdpy,i+2)
          endif
          cod(5)=0.d0
          call tcod(trans,cod,beam,.false.,fndcod)
          if(.not. fndcod)then
            cod=cod0
          endif
        endif
        irad=12
        call tinitr(trans)
        trans(:,7:12)=0.d0
        beam(1:21)=0.d0
        call tturne(trans,cod,beam,srot,iaez,
     $       .false.,.false.,.false.,.false.)
        if(fndcod)then
          cod(1:4)=cod(1:4)-codin(1:4)
        else
          cod(1:4)=cod0(1:4)-codin(1:4)
c          tra=trans(1:4,1:4)
c          dcod=cod0(1:4)-cod(1:4)
c          do k=1,4
c            tra(k,k)=tra(k,k)-1.d0
c          enddo
c          call tsolva(tra,dcod,cod,4,4,4,1.d-10)
c          cod(1:4)=cod(1:4)+cod0(1:4)-codin(1:4)
c      h0+x == h1+tr1.x
c      h0-h1 == (tr1-1).x
          write(*,'(a,i5,1p6g12.4)')'temits1-no cod ',i,cod(1:6)
        endif
c        write(*,'(a,2i5,1p7g12.4)')'temits1-cod ',i,i1,cod
        rm=rx
        call tmultr(rm,trans(:,7:12),6)
        call tmultr(rm,rxi,6)
        call tmov65(rm,trads(1,1,i))
        tws(1:ntwissfun,i)=tfetwiss(matmul(ri,tinv6(trans(:,1:6))),
     $       cod,.true.)
c        call tinv6(trans,rm)
c        call tmultr(rm,ri,6)
c        call tfetwiss(rm,cod,tws(1,i),.true.)
        call tmulbs(beam,rxi,.false.)
        beams(1:10,i)=beam(1:10)
c        write(*,'(a,i5,1p10g12.4)')'te ',i,beams(1:10,i)
        if(i .le. 0)then
          i=-i+1
        else
          i=-i
        endif
      enddo
      if(amusstep .eq. 0.d0)then
        ns=1
      else
        ns=int((amus1-amus0)/amusstep+1.1d0)
      endif
      dj=sigea**2/ndp*.5d0
      w=0.d0
      do j=1,ndp
        fj(j)=exp(-(j-.5d0)*dj/sige**2)
        w=w+fj(j)*dj/sige**2
      enddo
      call tefsetup(ha,hb,ba,bb,bd,ndp,mphi,mphi2,nz,ndims,
     $     beams,trads,tws,dj,btr,dpndim,sigea,tw0)
      phi0s=params(ipheff)*omega0*(params(iptrf0)+params(ipdz))/c
      fz=cos(phi0s)/(1.d0+cos(phi0s))*.25d0
     $     *(params(ipheff)*params(ipalphap)*pi2)**2
      if(kx .ne. 0)then
        kax=ktadaloc(-1,ns,kl)
      else
        kax=0
      endif
      do kk=1,ns
        amus=amus0+amusstep*(kk-1)
        do j=1,ndp
          amuj(j)=amus-fz*(j-.5d0)*dj/amus
        enddo
        call tesolvd(hff,hfb,hfx,
     $       hb,hc,hs,ha,amuj,cmu,smu,
     $       bd,hc,hs,
     $       ndp,mphi2,mphi,4,
     $       dj,damp,sige)
        nd1=max(1,min(int(abs(0.05d0/damp)),ndiffh))
        do k=1,nd1
          call tevdif(hfb,hb,hc,hs,ha,amuj,cmu,smu,
     $         bd,hc,hs,
     $         ndp,mphi2,mphi,4,
     $         dj,damp,sige,' ')
        enddo
c        write(*,'(a,i5,1p6g15.7)')'temits-4.1 ',kk,damp,sige,hc(1,1,1:4)
        call tesolvd(bff,bfb,bfx,
     $       bb,bc,bs,ba,amuj,cmu,smu,
     $       bd,hc,hs,
     $       ndp,mphi2,mphi,10,
     $       dj,damp,sige)
        call tesumb(bc,hc,mphi2,ndp,
     $       w,dj,sige,tws,ndims,beam,fj)
        beamr(1:21)=beam(1:21)
        call tmulbs(beamr,ri,.false.)
        vx=beamr(1)*beamr(3)-beamr(2)**2
        vy=beamr(6)*beamr(10)-beamr(9)**2
        emix0=sign(sqrt(abs(vx)),vx)
        emiy0=sign(sqrt(abs(vy)),vy)
        res0=1.d100
        res=1.d0
        it=0
        do while(res .gt. conv)
          call tevdif(hfb,hb,hc,hs,ha,amuj,cmu,smu,
     $         bd,hc,hs,
     $         ndp,mphi2,mphi,4,
     $         dj,damp,sige,' ')
          call tevdif(bfb,bb,bc,bs,ba,amuj,cmu,smu,
     $         bd,hc,hs,
     $         ndp,mphi2,mphi,10,
     $         dj,damp,sige,'b')
          call tesumb(bc,hc,mphi2,ndp,
     $         w,dj,sige,tws,ndims,beam,fj)
          beamr(1:21)=beam(1:21)
          call tmulbs(beamr,ri,.false.)
          vx=beamr(1)*beamr(3)-beamr(2)**2
          vy=beamr(6)*beamr(10)-beamr(9)**2
          emix=sign(sqrt(abs(vx)),vx)
          emiy=sign(sqrt(abs(vy)),vy)
          res=((emix-emix0)**2+(emiy-emiy0)**2)/(emix+emiy)**2
          emix0=emix
          emiy0=emiy
          it=it+1
          if(it .gt. itmax
     $         .or. it .gt. 1 .and.
     $         (res .gt. resconv .or. res .gt. res0))then
c            if(lfno .gt. 0)then
c              write(lfno,'(a,1pg15.7)')'Convergence failed. ',res
c            endif
            res0=-res
            res=0.d0
          else
            res0=res  
          endif
        enddo
        if(kx .eq. 0)then
          if(lfno .gt. 0)then
            vx=beam(1)*beam(3)-beam(2)**2
            vy=beam(6)*beam(10)-beam(9)**2
            emixp=sign(sqrt(abs(vx)),vx)
            emiyp=sign(sqrt(abs(vy)),vy)
            write(lfno,'(a,f9.6,4(a,1pg13.6),a,1pg11.3)')
     $           ' Nuz =',amus/pi2,
     $           '  Emittances X:',emix,
     $           '  Y:',emiy,
     $           '  Xproj:',emixp,
     $           '  Yproj:',emiyp,
     $           '  Conv.:',res0
            if(emiout)then
              write(lfno,*)
              write(lfno,*)'   Equiliblium beam matrix:'
              call tputbs(beam,label2,lfno)
              call tputbs(beamr,label1,lfno)
              call tedrawf(8,bc,bs,hc,hs,
     $             sige,amus/pi2,mphi2,ndp)
            endif
          endif
        else
          vx=beam(1)*beam(3)-beam(2)**2
          vy=beam(6)*beam(10)-beam(9)**2
          emixp=sign(sqrt(abs(vx)),vx)
          emiyp=sign(sqrt(abs(vy)),vy)
          kai=ktavaloc(0,6,kvl)
          kvl%rbody(1)=amus/pi2
          kvl%rbody(2)=emix
          kvl%rbody(3)=emiy
          kvl%rbody(4)=emixp
          kvl%rbody(5)=emiyp
          kvl%rbody(6)=res0
          kl%body(kk)=ktflist+kai
        endif
      enddo
      if(kx .ne. 0)then
        kx=ktflist+kax
      endif
c      call tfdebugprint(kx,'temits1',1)
      calpol=calpol0
      irtc=0
      return
      end
      
      subroutine tevdif(bfb,bb,bc,bs,ba,amuj,cmu,smu,
     $     bd,hc,hs,
     $     ndp,mphi2,mphi,nd,
     $     dj,damp0,sige,tag)
      use tfstk, only:ktfenanq
      implicit none
      integer*4 ndp,mphi2,mphi,nd
      real*8 bfb(nd*2*mphi2),bb(nd,mphi2,ndp),
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     ba(nd,nd,mphi,ndp),
     $     bd(10,4,mphi,ndp),
     $     amuj(ndp),
     $     cmu(mphi2),smu(mphi2),
     $     dj,damp,sige,eps
      integer*4 j,m,mb,ms,nsb,k1,k2,m1,m2,i1,nsbp
      real*8 diff,diff1,diff2,aj,bci1,bsi1,
     $     diff3,damp0,d,dcb(nd),dsb(nd)
      character*(*) tag
      nsb=nd*mphi2
      nsbp=nsb*ndp
      eps=sige**2
      damp=damp0*.5d0
      d=2.d0*damp
      call tevdifj (bc,bs,ndp,mphi2,nd,dj,eps,damp)
      call tevdifj1(bc,bs,ndp,mphi2,nd,dj,eps,damp)
c      write(*,'(a,1p7g15.7)')'tevdif-1 ',dj,eps,damp,bc(1,1,1:4)
      do j=1,ndp
        aj=(j-.5d0)*dj
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        bc(:,m,j)=bc(:,m,j)+(diff+diff3)*bc(:,m,j)
        bs(:,m,j)=bs(:,m,j)+(diff-diff3)*bs(:,m,j)
c        do i=1,nd
c          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
c          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
c        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          bc(:,m,j)=diff*bc(:,m,j)
          bs(:,m,j)=diff*bs(:,m,j)
c          do i=1,nd
c            bc(i,m,j)=diff*bc(i,m,j)
c            bs(i,m,j)=diff*bs(i,m,j)
c          enddo
        enddo
        do m=1,mphi2-2
          m2=m+2
          diff2=-d*(m+1)*
     $         (-1.d0+eps/aj*.5d0*(m-1))*.25d0
          bc(:,m2,j)=bc(:,m2,j)+diff2*bc(:,m,j)
          bs(:,m2,j)=bs(:,m2,j)+diff2*bs(:,m,j)
c          do i=1,nd
c            bc(i,m2,j)=bc(i,m2,j)+diff2*bc(i,m,j)
c            bs(i,m2,j)=bs(i,m2,j)+diff2*bs(i,m,j)
c          enddo
        enddo
        do m=mphi2,4,-1
          m1=m-2
          diff1=-d*(m-3)*
     $         (1.d0+eps/aj*.5d0*(m-1))*.25d0
          bc(:,m1,j)=bc(:,m1,j)+diff1*bc(:,m,j)
          bs(:,m1,j)=bs(:,m1,j)+diff1*bs(:,m,j)
c          do i=1,nd
c            bc(i,m1,j)=bc(i,m1,j)+diff1*bc(i,m,j)
c            bs(i,m1,j)=bs(i,m1,j)+diff1*bs(i,m,j)
c          enddo
        enddo
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        bc(:,m,j)=bc(:,m,j)+(diff+diff3)*bc(:,m,j)
        bs(:,m,j)=bs(:,m,j)+(diff-diff3)*bs(:,m,j)
c        do i=1,nd
c          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
c          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
c        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          bc(:,m,j)=diff*bc(:,m,j)
          bs(:,m,j)=diff*bs(:,m,j)
c          do i=1,nd
c            bc(i,m,j)=diff*bc(i,m,j)
c            bs(i,m,j)=diff*bs(i,m,j)
c          enddo
        enddo
      enddo
c      write(*,*)'tevdif-2 ',bc(1,1,1:4)
      do j=1,ndp
        bfb=0.d0
c        call tclr(bfb,2*nsb)
        if(nd .eq. 10)then
          call tesetdm(j,bfb,bd,hc,hs,ndp,mphi,mphi2)
        endif
        call tadd(bb(1,1,j),bfb,bfb,nsb)
        call tetrig(amuj(j),cmu,smu,mphi2)
        do m=1,mphi2
          mb=(m-1)*nd
          ms=mb+nsb
          do i1=1,nd
            bci1=bc(i1,1,j)
            bsi1=bs(i1,1,j)
            bfb(mb+1:mb+nd)=bfb(mb+1:mb+nd)+ba(:,i1,m,j)*bci1
            bfb(ms+1:ms+nd)=bfb(ms+1:ms+nd)+ba(:,i1,m,j)*bsi1
c            do i=1,nd
c              bfb(mb+i)=bfb(mb+i)+ba(i,i1,m,j)*bci1
c              bfb(ms+i)=bfb(ms+i)+ba(i,i1,m,j)*bsi1
c            enddo
          enddo
          do m1=2,mphi2
            k1=abs(m-m1)+1
            k2=abs(m+m1-2)+1
            do i1=1,nd
              bci1=bc(i1,m1,j)
              bsi1=bs(i1,m1,j)
              bfb(mb+1:mb+nd)=bfb(mb+1:mb+nd)
     $             +(ba(:,i1,k1,j)+ba(:,i1,k2,j))*bci1
              bfb(ms+1:ms+nd)=bfb(ms+1:ms+nd)
     $             +(ba(:,i1,k1,j)-ba(:,i1,k2,j))*bsi1
c              do i=1,nd
c                bfb(mb+i)=bfb(mb+i)
c     $               +(ba(i,i1,k1,j)+ba(i,i1,k2,j))*bci1
c                bfb(ms+i)=bfb(ms+i)
c     $               +(ba(i,i1,k1,j)-ba(i,i1,k2,j))*bsi1
c              enddo
            enddo
          enddo
        enddo
        do m=1,mphi2
          mb=(m-1)*nd
          ms=mb+nsb
          dcb=bfb(mb+1:mb+nd)
          dsb=bfb(ms+1:ms+nd)
          bc(:,m,j)= cmu(m)*dcb+smu(m)*dsb
          bs(:,m,j)=-smu(m)*dcb+cmu(m)*dsb
c          do i=1,nd
c            dc=bfb(mb+i)
c            ds=bfb(ms+i)
c            bc(i,m,j)=cmu(m)*dc+smu(m)*ds
c            bs(i,m,j)=-smu(m)*dc+cmu(m)*ds
c          enddo
        enddo
      enddo
      do j=1,ndp
        aj=(j-.5d0)*dj
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        bc(:,m,j)=bc(:,m,j)+(diff+diff3)*bc(:,m,j)
        bs(:,m,j)=bs(:,m,j)+(diff-diff3)*bs(:,m,j)
c        do i=1,nd
c          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
c          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
c        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          bc(:,m,j)=diff*bc(:,m,j)
          bs(:,m,j)=diff*bs(:,m,j)
c          do i=1,nd
c            bc(i,m,j)=diff*bc(i,m,j)
c            bs(i,m,j)=diff*bs(i,m,j)
c          enddo
        enddo
        do m=mphi2,4,-1
          m1=m-2
          diff1=-d*(m-3)*
     $         (1.d0+eps/aj*.5d0*(m-1))*.25d0
          bc(:,m1,j)=bc(:,m1,j)+diff1*bc(:,m,j)
          bs(:,m1,j)=bs(:,m1,j)+diff1*bs(:,m,j)
c          do i=1,nd
c            bc(i,m1,j)=bc(i,m1,j)+diff1*bc(i,m,j)
c            bs(i,m1,j)=bs(i,m1,j)+diff1*bs(i,m,j)
c          enddo
        enddo
        do m=1,mphi2-2
          m2=m+2
          diff2=-d*(m+1)*
     $         (-1.d0+eps/aj*.5d0*(m-1))*.25d0
          bc(:,m2,j)=bc(:,m2,j)+diff2*bc(:,m,j)
          bs(:,m2,j)=bs(:,m2,j)+diff2*bs(:,m,j)
c          do i=1,nd
c            bc(i,m2,j)=bc(i,m2,j)+diff2*bc(i,m,j)
c            bs(i,m2,j)=bs(i,m2,j)+diff2*bs(i,m,j)
c          enddo
        enddo
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        bc(:,m,j)=bc(:,m,j)+(diff+diff3)*bc(:,m,j)
        bs(:,m,j)=bs(:,m,j)+(diff-diff3)*bs(:,m,j)
c        do i=1,nd
c          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
c          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
c        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          bc(:,m,j)=diff*bc(:,m,j)
          bs(:,m,j)=diff*bs(:,m,j)
c          do i=1,nd
c            bc(i,m,j)=diff*bc(i,m,j)
c            bs(i,m,j)=diff*bs(i,m,j)
c          enddo
        enddo
      enddo
c      write(*,*)'tevdif-3 ',bc(1,1,1:4)
      call tevdifj1(bc,bs,ndp,mphi2,nd,dj,eps,damp)
      call tevdifj (bc,bs,ndp,mphi2,nd,dj,eps,damp)
c      write(*,*)'tevdif-9 ',bc(1,1,1:4)
      return
      end

      subroutine tevdifj(bc,bs,ndp,mphi2,nd,dj,eps,damp)
      implicit none
      integer*4 j,m,nd,mphi2,m1,ndp
      real*8 eps,damp,dj,
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     fbc(nd),fbs(nd),dfbc(nd),dfbs(nd)
      do j=1,ndp
        do m=1,mphi2
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          if(m .eq. 1)then
            bc(:,m,j)=bc(:,m,j)+2.d0*dfbc
          else
            bc(:,m,j)=bc(:,m,j)+dfbc
            bs(:,m,j)=bs(:,m,j)+dfbs
          endif
          m1=m+2
          if(m1 .le. mphi2)then
            bc(:,m1,j)=bc(:,m1,j)+dfbc-m*fbc
            bs(:,m1,j)=bs(:,m1,j)+dfbs-m*fbs
          endif
        enddo
      enddo
      return
      end

      subroutine tevdifj1(bc,bs,ndp,mphi2,nd,dj,eps,damp)
      implicit none
      integer*4 j,m,nd,mphi2,m2,ndp
      real*8 eps,damp,dj,
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     fbc(nd),fbs(nd),dfbc(nd),dfbs(nd)
      do j=ndp,1,-1
        do m=mphi2,1,-1
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          m2=m-2
          if(m .eq. 3)then
            bc(:,m2,j)=bc(:,m2,j)+dfbc+2.d0*fbc
          elseif(m .gt. 3)then
            bc(:,m2,j)=bc(:,m2,j)+dfbc+(m-2)*fbc
            bs(:,m2,j)=bs(:,m2,j)+dfbs+(m-2)*fbs
          endif
          bc(:,m,j)=bc(:,m,j)+dfbc
          bs(:,m,j)=bs(:,m,j)+dfbs
        enddo
      enddo
      return
      end
      
      subroutine tesolvd(bff,bfb,bfx,
     $     bb,bc,bs,ba,amuj,cmu,smu,
     $     bd,hc,hs,
     $     ndp,mphi2,mphi,nd,
     $     dj,damp0,sige)
      use tfstk,only:ktfenanq,resetnan
      implicit none
      integer*4 ndp,mphi2,mphi,nd
      real*8 bfb(nd*2*mphi2),bfx(nd*2*mphi2),bb(nd,mphi2,ndp),
     $     bff(2*nd*mphi2,2*nd*mphi2),
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     ba(nd,nd,mphi,ndp),
     $     bd(10,4,mphi,ndp),
     $     amuj(ndp),
     $     cmu(mphi2),smu(mphi2),
     $     dj,damp,sige,eps
      integer*4 j,m,mb,ms,nsb,k1,k2,m1,m2,
     $     ia,ma,maci,masi,mc
      real*8 diff,diff1,diff2,aj,bci1,bsi1,
     $     diff3,dcb(nd),dsb(nd),damp0,d
      nsb=nd*mphi2
      eps=sige**2
      damp=damp0*.5d0
      d=2.d0*damp
      do j=1,ndp
        aj=(j-.5d0)*dj
        do ma=1,mphi2
          do ia=1,nd
            bc(:,:,j)=0.d0
            bs(:,:,j)=0.d0
c            call tclr(bc(1,1,j),nsb)
c            call tclr(bs(1,1,j),nsb)
            bc(ia,ma,j)=1.d0
            bs(ia,ma,j)=1.d0
            m=ma
            if(m .eq. 2)then
              diff=d*eps*(m-1)**2*.125d0/aj
              diff3=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.125d0
              bc(ia,m,j)=bc(ia,m,j)+(diff+diff3)*bc(ia,m,j)
              bs(ia,m,j)=bs(ia,m,j)+(diff-diff3)*bs(ia,m,j)
            else
              diff=1.d0+d*eps*(m-1)**2*.125d0/aj
              bc(ia,m,j)=diff*bc(ia,m,j)
              bs(ia,m,j)=diff*bs(ia,m,j)
            endif
            m2=ma+2
            if(m2 .le. mphi2)then
              diff2=-d*(m+1)*
     $             (-1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m2,j)=bc(ia,m2,j)+diff2*bc(ia,m,j)
              bs(ia,m2,j)=bs(ia,m2,j)+diff2*bs(ia,m,j)
              diff1=-d*(m2-3)*
     $             (1.d0+eps/aj*.5d0*(m2-1))*.25d0
              bc(ia,ma,j)=bc(ia,ma,j)+diff1*bc(ia,m2,j)
              bs(ia,ma,j)=bs(ia,ma,j)+diff1*bs(ia,m2,j)
            endif
            if(ma .gt. 3)then
              m1=ma-2
              diff1=-d*(ma-3)*
     $             (1.d0+eps/aj*.5d0*(ma-1))*.25d0
              bc(ia,m1,j)=bc(ia,m1,j)+diff1*bc(ia,ma,j)
              bs(ia,m1,j)=bs(ia,m1,j)+diff1*bs(ia,ma,j)
            endif
            do m=ma-2,ma+2,2
              if(m .eq. 2)then
                diff=d*eps*(m-1)**2*.125d0/aj
                diff3=-d*(m-3)*
     $               (1.d0+eps/aj*.5d0*(m-1))*.125d0
                bc(ia,m,j)=bc(ia,m,j)+(diff+diff3)*bc(ia,m,j)
                bs(ia,m,j)=bs(ia,m,j)+(diff-diff3)*bs(ia,m,j)
              elseif(m .ge. 3 .and. m .le. mphi2)then
                diff=1.d0+d*eps*(m-1)**2*.125d0/aj
                bc(ia,m,j)=diff*bc(ia,m,j)
                bs(ia,m,j)=diff*bs(ia,m,j)
              endif
            enddo
            call tetrig(amuj(j),cmu,smu,mphi2)
            bfb=0.d0
c            call tclr(bfb,2*nsb)
            do m=1,mphi2
              mb=(m-1)*nd
              ms=mb+nsb
              do m1=ma-2,ma+2,2
                if(m1 .eq. 1)then
                  bci1=bc(ia,1,j)
                  bsi1=bs(ia,1,j)
                  bfb(mb+1:mb+nd)=bfb(mb+1:mb+nd)+ba(:,ia,m,j)*bci1
                  bfb(ms+1:ms+nd)=bfb(ms+1:ms+nd)+ba(:,ia,m,j)*bsi1
                elseif(m1 .gt. 0 .and. m1 .le. mphi2)then
                  k1=abs(m-m1)+1
                  k2=abs(m+m1-2)+1
                  bci1=bc(ia,m1,j)
                  bsi1=bs(ia,m1,j)
                  bfb(mb+1:mb+nd)=bfb(mb+1:mb+nd)
     $                 +(ba(:,ia,k1,j)+ba(:,ia,k2,j))*bci1
                  bfb(ms+1:ms+nd)=bfb(ms+1:ms+nd)
     $                 +(ba(:,ia,k1,j)-ba(:,ia,k2,j))*bsi1
                endif
              enddo
            enddo
            maci=(ma-1)*nd+ia
            masi=maci+nsb
            do m=1,mphi2
              mc=(m-1)*nd
              ms=mc+nsb
              dcb=bfb(mc+1:mc+nd)
              dsb=bfb(ms+1:ms+nd)
              bff(mc+1:mc+nd,maci)=-cmu(m)*dcb
              bff(ms+1:ms+nd,maci)= smu(m)*dcb
              bff(mc+1:mc+nd,masi)=-smu(m)*dsb
              bff(ms+1:ms+nd,masi)=-cmu(m)*dsb
            enddo
            bc(:,:,j)=0.d0
            bs(:,:,j)=0.d0
            bc(ia,ma,j)=1.d0
            bs(ia,ma,j)=1.d0
            if(ma .eq. 2)then
              m=2
              diff=d*eps*(m-1)**2*.125d0/aj
              diff3=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.125d0
              bc(ia,m,j)=bc(ia,m,j)/(1.d0+(diff+diff3))
              bs(ia,m,j)=bs(ia,m,j)/(1.d0+(diff-diff3))
            elseif(ma .gt. 2)then
              m=ma
              diff=1.d0+d*eps*(m-1)**2*.125d0/aj
              bc(ia,m,j)=bc(ia,m,j)/diff
              bs(ia,m,j)=bs(ia,m,j)/diff
            endif
            if(ma .le. mphi2-2)then
              m=ma
              m2=m+2
              diff2=-d*(m+1)*
     $             (-1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m2,j)=bc(ia,m2,j)-diff2*bc(ia,m,j)
              bs(ia,m2,j)=bs(ia,m2,j)-diff2*bs(ia,m,j)
              if(m2 .gt. 3)then
                diff1=-d*(m2-3)*
     $               (1.d0+eps/aj*.5d0*(m2-1))*.25d0
                bc(ia,ma,j)=bc(ia,ma,j)-diff1*bc(ia,m2,j)
                bs(ia,ma,j)=bs(ia,ma,j)-diff1*bs(ia,m2,j)
              endif
            endif
            if(ma .gt. 3)then
              m=ma
              m1=m-2
              diff1=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m1,j)=bc(ia,m1,j)-diff1*bc(ia,m,j)
              bs(ia,m1,j)=bs(ia,m1,j)-diff1*bs(ia,m,j)
            endif
            do m=ma-2,ma+2,2
              if(m .eq. 2)then
                diff=d*eps*(m-1)**2*.125d0/aj
                diff3=-d*(m-3)*
     $               (1.d0+eps/aj*.5d0*(m-1))*.125d0
                bc(ia,m,j)=bc(ia,m,j)/(1.d0+(diff+diff3))
                bs(ia,m,j)=bs(ia,m,j)/(1.d0+(diff-diff3))
              elseif(m .ge. 3 .and. m .le. mphi2)then
                diff=1.d0+d*eps*(m-1)**2*.125d0/aj
                bc(ia,m,j)=bc(ia,m,j)/diff
                bs(ia,m,j)=bs(ia,m,j)/diff
              endif
            enddo
            call tadd(bff(1,maci),bc(1,1,j),bff(1,maci),nsb)
            call tadd(bff(nsb+1,masi),bs(1,1,j),bff(nsb+1,masi),nsb)
          enddo
        enddo
        if(nd .eq. 10)then
          call tesetdm(j,bfb,bd,hc,hs,ndp,mphi,mphi2)
        else
          bfb=0.d0
        endif
        call tadd(bb(1,1,j),bfb,bfb,nsb)
        do m=1,mphi2
          mc=(m-1)*nd
          ms=mc+nsb
          dcb=bfb(mc+1:mc+nd)
          dsb=bfb(ms+1:ms+nd)
          bfb(mc+1:mc+nd)= cmu(m)*dcb+smu(m)*dsb
          bfb(ms+1:ms+nd)=-smu(m)*dcb+cmu(m)*dsb
        enddo
        call tsolva(bff,bfb,bfx,2*nsb,2*nsb,2*nsb,1.d-20)
c        if(nd .eq. 4 .and. bfx(1) .gt. 1.d10)then
c          write(*,'(a,1p10g11.3)')'tesolvd-8 ',j,nsb,bfx(1),bfb(1)
c        endif
        call resetnan(bfx)
        call tmov(bfx,bc(1,1,j),nsb)
        call tmov(bfx(nsb+1),bs(1,1,j),nsb)
      enddo
c      write(*,'(a,1p10g11.3)')'tesolvd-9 ',bs(:,1,ndp)
      return
      end

      subroutine tesetdm(j,bfb,bd,hc,hs,ndp,mphi,mphi2)
      implicit none
      integer*4 ndp,mphi,mphi2,j,m,m1,i1,k1,k2,mb,ms,nsb
      real*8
     $     bfb(20*mphi2),bd(10,4,mphi,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     hci1,hsi1
      nsb=10*mphi2
      bfb=0.d0
c      call tclr(bfb,2*nsb)
      do m=1,mphi2
        mb=(m-1)*10
        ms=mb+nsb
        do i1=1,4
          hci1=hc(i1,1,j)
          bfb(mb+1:mb+10)=bfb(mb+1:mb+10)+bd(:,i1,m,j)*hci1
        enddo
        do m1=2,mphi2
          k1=abs(m-m1)+1
          k2=abs(m+m1-2)+1
          do i1=1,4
            hci1=hc(i1,m1,j)
            hsi1=hs(i1,m1,j)
            bfb(mb+1:mb+10)=bfb(mb+1:mb+10)
     $           +(bd(:,i1,k1,j)+bd(:,i1,k2,j))*hci1
            bfb(ms+1:ms+10)=bfb(ms+1:ms+10)
     $           +(bd(:,i1,k1,j)-bd(:,i1,k2,j))*hsi1
          enddo
        enddo
      enddo
      return
      end

      subroutine tetrig(amu,cmu,smu,mphi2)
      implicit none
      integer*4 mphi2,m
      real*8 cmu(mphi2),smu(mphi2),phim,amu
      do m=1,mphi2
        phim=(m-1)*amu
        cmu(m)=cos(phim)
        smu(m)=sin(phim)
      enddo
      return
      end

      subroutine tefsetup(ha,hb,ba,bb,bd,
     $     ndp,mphi,mphi2,nz,ndims,
     $     beams,trads,tws,dj,btr,dpndim,sigea,tw0)
      use macmath
      use tfstk,only:ktfenanq
      use tffitcode,only:ntwissfun
      use sad_main, only:iaidx
      implicit none
      integer*4 mphi,mphi2,nz,ndp,i,j,k,l,m,ndims
      real*8 
     $     beams(10,-ndims:ndims),
     $     trads(5,5,-ndims:ndims),
     $     tws(ntwissfun,-ndims:ndims),tw0(ntwissfun),
     $     ba(10,10,mphi,ndp),bb(10,mphi2,ndp),
     $     bd(10,4,mphi,ndp),
     $     ha(4,4,mphi,ndp),hb(4,mphi2,ndp)
      real*8 phi,cs,phim,csm,p,f,tr1(5,5),h1(4),
     $     beama(10),trd(5,5),dpndim,sigea,btr(10,10),
     $     aj,pj,dj
      integer*4 ip,kk,k1,ll,l1,nj,j1
      parameter (nj=2)
c      iaidx(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      bb=0.d0
      bd=0.d0
      ba=0.d0
      hb=0.d0
      ha=0.d0
      do j=1,ndp
        do j1=-nj,nj
          aj=(j-.5d0)*dj+dj/(2*nj+1)*j1
          pj=sqrt(aj*2.d0)
          do i=1,nz
            phi=(i-1)*pi/(nz-1)
            cs=cos(phi)
            p=pj*cs
            ip=int((p+sigea)/dpndim)-ndims
            f=p/dpndim-ip
            call teintp1(ip,f,tws,tr1,h1,tw0,ndims)
c            write(*,'(a,5(1p5g12.4))')'setup1  ',h1
c            call teintp(ip,f,tws,tr1,h1,1.d0,ndims)
c            write(*,'(a,5(1p5g12.4/))')'setup11 ',(tr1(k,:),k=1,5)
            call teintm(ip,f,trads,trd,ndims,1.d0)
            tr1=tr1+trd
            call teintb(ip,f,beams,beama,ndims)
            do k=1,4
              do k1=1,k
                kk=iaidx(k,k1)
                do l=1,4
                  do l1=1,l
                    ll=iaidx(l,l1)
                    if(l .eq. l1)then
                      btr(kk,ll)=tr1(k,l)*tr1(k1,l)
c     btrd(kk,ll)=tr1(k,l)*trd(k1,l)+
c     $                   trd(k,l)*tr1(k1,l)+trd(k,l)*trd(k1,l)
                    else
                      btr(kk,ll)=tr1(k,l)*tr1(k1,l1)+
     $                     tr1(k,l1)*tr1(k1,l)
c     btrd(kk,ll)=tr1(k,l)*trd(k1,l1)+
c     $                   trd(k,l)*tr1(k1,l1)+
c     $                   trd(k,l)*trd(k1,l1)+
c     $                   tr1(k,l1)*trd(k1,l)+
c     $                   trd(k,l1)*tr1(k1,l)+
c     $                   trd(k,l1)*trd(k1,l)
                    endif
                  enddo
                enddo
              enddo
            enddo
            do m=1,mphi2
              phim=phi*(m-1)
              csm=cos(phim)/(nz-1)/(2*nj+1)
              if(i .eq. 1 .or. i .eq. nz)then
                csm=csm*.5d0
              endif
              do k=1,4
                ha(1:4,k,m,j)=ha(1:4,k,m,j)+tr1(1:4,k)*csm
                hb(k,m,j)=hb(k,m,j)+h1(k)*csm
                do l=1,k
                  kk=iaidx(k,l)
                  bb(kk,m,j)=bb(kk,m,j)+h1(k)*h1(l)*csm
                  bd(kk,:,m,j)=bd(kk,:,m,j)+
     $                 (tr1(k,1:4)*h1(l)+tr1(l,1:4)*h1(k))*csm
                enddo
              enddo
              ba(:,:,m,j)=ba(:,:,m,j)+btr*csm
              bb(:,m,j)=bb(:,m,j)+beama*csm
c              do k=1,10
c                do l=1,10
c                  ba(l,k,m,j)=ba(l,k,m,j)+btr(l,k)*csm
cc     bd(l,k,m,j)=bd(l,k,m,j)+btrd(l,k)*csm
c                enddo
c                bb(k,m,j)=bb(k,m,j)+beama(k)*csm
c              enddo
            enddo
            do m=mphi2+1,mphi
              phim=phi*(m-1)
              csm=cos(phim)/(nz-1)/(2*nj+1)
              if(i .eq. 1 .or. i .eq. nz)then
                csm=csm*.5d0
              endif
              do k=1,4
                ha(:,k,m,j)=ha(:,k,m,j)+tr1(1:4,k)*csm
                do l=1,k
                  kk=iaidx(k,l)
                  bd(kk,:,m,j)=bd(kk,:,m,j)+
     $                 (tr1(k,1:4)*h1(l)+tr1(l,1:4)*h1(k))*csm
                enddo
              enddo
              ba(:,:,m,j)=ba(:,:,m,j)+btr*csm
c              do k=1,10
c                do l=1,10
c                  ba(l,k,m,j)=ba(l,k,m,j)+btr(l,k)*csm
c                enddo
c              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine tespl(bx,bdx,bddx,ndp,mphi2,nd,dj,e,damp,j0,m)
      implicit none
      integer*4 ndp,mphi2,nd,m,j0
      real*8 bx(nd,mphi2,ndp),bddx(nd),bdx(nd),dj,e,damp
      real*8 ajj1,ajj2,dy1(nd),d,ajj,f1,f2,bx1(nd),bx2(nd)
      real*8 , parameter :: c1=11.d0/60.d0
      d=abs(2.d0*damp)
      ajj=(j0-.5d0)*dj
      f2=d*.5d0
      f1=e*f2
      ajj1=ajj-dj
      ajj2=ajj+dj
      if(j0 .eq. 1)then
        bx1=bx(:,m,1)
        bx2=bx(:,m,2)
      elseif(j0 .eq. ndp)then
        bx1=bx(:,m,j0-1)
        bx2=0.d0
      else
        bx1=bx(:,m,j0-1)
        bx2=bx(:,m,j0+1)
      endif
      dy1=(bx2-bx1)/dj*.5d0
      bdx=f1*dy1
      bddx=f2*(-(e+ajj)*dy1
     $     +e*(ajj1*bx1-2.d0*ajj*bx(:,m,j0)+ajj2*bx2)/dj**2)
      if(j0 .eq. ndp)then
        bddx=max(0.d0,bddx)
      endif
c        do i=1,nd
c          dy1=(bx(i,m,j0+1)-bx(i,m,j0-1))/dj*.5d0
c          bdx(i)=f1*dy1
c          bddx(i)=f2*(-(e+ajj)*dy1
c     $         +e*(ajj1*bx(i,m,j0-1)-2.d0*ajj*bx(i,m,j0)
c     $         +ajj2*bx(i,m,j0+1))/dj**2)
c        enddo
c$$$      if(up)then
c$$$        do j=1,j0-1
c$$$          do i=1,nd
c$$$            bdx(i,m,j)=bx(i,m,j)*fj(j)
c$$$          enddo
c$$$        enddo
c$$$        do i=1,nd
c$$$          bdx(i,m,j0)=bx(i,m,j0)*fj(j0)*0.5d0
c$$$        enddo
c$$$        do j=j0+1,ndp
c$$$          do i=1,nd
c$$$            bdx(i,m,j)=0.d0
c$$$          enddo
c$$$        enddo
c$$$      else
c$$$        do j=1,j0-1
c$$$          do i=1,nd
c$$$            bdx(i,m,j)=0.d0
c$$$          enddo
c$$$        enddo
c$$$        do i=1,nd
c$$$          bdx(i,m,j0)=bx(i,m,j0)*fj(j0)*0.5d0
c$$$        enddo
c$$$        do j=j0+1,ndp
c$$$          do i=1,nd
c$$$            bdx(i,m,j)=bx(i,m,j)*fj(j)
c$$$          enddo
c$$$        enddo
c$$$      endif
c$$$      call spline1m(nd,nd*mphi2,ndp,bdx(1,m,1),bddx(1,m,1),work)
c$$$      if(j0 .eq. ndp)then
c$$$        do i=1,nd
c$$$          ddy0=bddx(i,m,ndp-1)/dj**2
c$$$          y0=bdx(i,m,ndp-1)
c$$$          ddy1=bddx(i,m,ndp)/dj**2
c$$$          y1=bdx(i,m,ndp)
c$$$          bdx(i,m,ndp)=(y1-y0)/dj+(ddy0+2.d0*ddy1)*dj
c$$$          bddx(i,m,ndp)=y1+(e+ajj)*bdx(i,m,ndp)
c$$$     $         +e*ajj*6.d0*ddy1
c$$$        enddo
c$$$      else
c$$$        do i=1,nd
c$$$          ddy0=bddx(i,m,j0)/dj**2
c$$$          y0=bdx(i,m,j0)
c$$$          ddy1=bddx(i,m,j0+1)/dj**2
c$$$          y1=bdx(i,m,j0+1)
c$$$          bdx(i,m,j0)=(y1-y0)/dj-(2.d0*ddy0+ddy1)*dj
c$$$          bddx(i,m,j0)=y0+(e+ajj)*bdx(i,m,j0)
c$$$     $         +e*ajj*6.d0*ddy0
c$$$        enddo
c$$$      endif
c$$$      call ttimes(bddx(1,m,j0),.5d0*d/fj(j0),bddx(1,m,j0),nd)
c$$$      call ttimes(bdx(1,m,j0),.5d0*d*e/fj(j0),
c$$$     $     bdx(1,m,j0),nd)
      return
      end

      subroutine tesumb(bc,hc,mphi2,ndp,
     $     w,dj,sige,tws,ndims,beam,fj)
      use tffitcode
      use sad_main, only:iaidx
      implicit none
      integer*4 mphi2,ndp,ndims
      real*8 beam(42),sige,f,dj,w,e,aj
      real*8 bc(10,mphi2,ndp),hc(4,mphi2,ndp),
     $     tws(ntwissfun,-ndims:ndims),fj(ndp),dispp(4)
      integer*4 i,j,i1,ii
c      iaidx(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      e=sige**2
      beam(1:21)=0.d0
c      call tclr(beam,21)
      do j=1,ndp
        f=fj(j)*dj/e/w
        beam(1:10)=beam(1:10)+bc(:,1,j)*f
      enddo
      do j=1,ndp
        aj=(j-.5d0)*dj
        f=fj(j)*dj/e/w*sqrt(2.d0*aj)
        beam(16:19)=beam(16:19)+hc(:,2,j)*f
c        do i=1,4
c          beam(15+i)=beam(15+i)+hc(i,2,j)*f
c        enddo
      enddo
      call tgetphysdispu(tws(1,0),dispp)
      do i=1,4
        do i1=1,i
          ii=iaidx(i,i1)
          beam(ii)=beam(ii)+dispp(i)*dispp(i1)*e
        enddo
      enddo
      beam(21)=e
      return
      end
      
      subroutine teintb(ip,f,beams,beama,ndims)
      implicit none
      integer*4 ip,ip1,ndims
      real*8 f,beams(10,-ndims:ndims),beama(10),f1
      f1=1.d0-f
      ip1=ip+1
      beama=f1*beams(:,ip)+f*beams(:,ip1)
c      do i=1,10
c        beama(i)=f1*beams(i,ip)+f*beams(i,ip1)
c      enddo
c     write(*,'(1p10g12.4)')beama
      return
      end
      
      subroutine teintm(ip,f,trads,trd,ndims,wd)
      implicit none
      integer*4 ip,ndims
      real*8 f,trads(5,5,-ndims:ndims),f1,f0,wd
      real*8 trd(5,5)
      integer*4 ip1
      f1=(1.d0-f)*wd
      f0=f*wd
      ip1=ip+1
      trd=f1*trads(:,:,ip)+f0*trads(:,:,ip1)
c      do i=1,5
c        trd(1,i)=f1*trads(1,i,ip)+f0*trads(1,i,ip1)
c        trd(2,i)=f1*trads(2,i,ip)+f0*trads(2,i,ip1)
c        trd(3,i)=f1*trads(3,i,ip)+f0*trads(3,i,ip1)
c        trd(4,i)=f1*trads(4,i,ip)+f0*trads(4,i,ip1)
c        trd(5,i)=f1*trads(5,i,ip)+f0*trads(5,i,ip1)
c      enddo
      return
      end
      
      subroutine teintp(ip,f,tws,tr1,h1,rs,ndims)
      use tfstk, only:ktfenanq
      use tffitcode
      implicit none
      integer*4 ip,ip1,ndims
      real*8 tws(ntwissfun,-ndims:ndims),tr1(5,5),h1(4)
      real*8 f,rs,f1,axi,bxi,amuxi,ayi,byi,amuyi,
     $     exi,epxi,eyi,epyi,r1i,r2i,r3i,r4i,
     $     cosmux,sinmux,cosmuy,sinmuy,amu,
     $     a11,a12,a21,a22,a33,a34,a43,a44,
     $     b11,b12,b13,b14,b21,b22,b23,b24,
     $     b31,b32,b33,b34,b41,b42,b43,b44
      ip1=ip+1
      f1=1.d0-f
      axi=f*tws(mfitax,ip1)+f1*tws(mfitax,ip)
      bxi=f*tws(mfitbx,ip1)+f1*tws(mfitbx,ip)
      amuxi=(f*tws(mfitnx,ip1)+f1*tws(mfitnx,ip))*rs
      ayi=f*tws(mfitay,ip1)+f1*tws(mfitay,ip)
      byi=f*tws(mfitby,ip1)+f1*tws(mfitby,ip)
      amuyi=(f*tws(mfitny,ip1)+f1*tws(mfitny,ip))*rs
      exi=f*tws(mfitex,ip1)+f1*tws(mfitex,ip)
      epxi=f*tws(mfitepx,ip1)+f1*tws(mfitepx,ip)
      eyi=f*tws(mfitey,ip1)+f1*tws(mfitey,ip)
      epyi=f*tws(mfitepy,ip1)+f1*tws(mfitepy,ip)
      r1i=f*tws(mfitr1,ip1)+f1*tws(mfitr1,ip)
      r2i=f*tws(mfitr2,ip1)+f1*tws(mfitr2,ip)
      r3i=f*tws(mfitr3,ip1)+f1*tws(mfitr3,ip)
      r4i=f*tws(mfitr4,ip1)+f1*tws(mfitr4,ip)
      h1(1:4)=f*tws(mfitdx:mfitdpy,ip1)+f1*tws(mfitdx:mfitdpy,ip)
      cosmux=cos(amuxi)
      sinmux=sin(amuxi)
      cosmuy=cos(amuyi)
      sinmuy=sin(amuyi)
      a11=cosmux+axi*sinmux
      a12=bxi*sinmux
      a21=-(1.d0+axi**2)/bxi*sinmux
      a22=cosmux-axi*sinmux
      a33=cosmuy+ayi*sinmuy
      a34=byi*sinmuy
      a43=-(1.d0+ayi**2)/byi*sinmuy
      a44=cosmuy-ayi*sinmuy
      amu=sqrt(1.d0-r1i*r4i+r2i*r3i)
      b11=amu*a11
      b12=amu*a12
      b13= r4i*a33-r2i*a43
      b14= r4i*a34-r2i*a44
      b21=amu*a21
      b22=amu*a22
      b23=-r3i*a33+r1i*a43
      b24=-r3i*a34+r1i*a44
      b33=amu*a33
      b34=amu*a34
      b31=-r1i*a11-r2i*a21
      b32=-r1i*a12-r2i*a22
      b43=amu*a43
      b44=amu*a44
      b41=-r3i*a11-r4i*a21
      b42=-r3i*a12-r4i*a22
      tr1(1,1)=b11*amu+b13*r1i+b14*r3i
      tr1(1,2)=b12*amu+b13*r2i+b14*r4i
      tr1(1,3)=b13*amu-b11*r4i+b12*r3i
      tr1(1,4)=b14*amu+b11*r2i-b12*r1i
      tr1(2,1)=b21*amu+b23*r1i+b24*r3i
      tr1(2,2)=b22*amu+b23*r2i+b24*r4i
      tr1(2,3)=b23*amu-b21*r4i+b22*r3i
      tr1(2,4)=b24*amu+b21*r2i-b22*r1i
      tr1(3,1)=b31*amu+b33*r1i+b34*r3i
      tr1(3,2)=b32*amu+b33*r2i+b34*r4i
      tr1(3,3)=b33*amu-b31*r4i+b32*r3i
      tr1(3,4)=b34*amu+b31*r2i-b32*r1i
      tr1(4,1)=b41*amu+b43*r1i+b44*r3i
      tr1(4,2)=b42*amu+b43*r2i+b44*r4i
      tr1(4,3)=b43*amu-b41*r4i+b42*r3i
      tr1(4,4)=b44*amu+b41*r2i-b42*r1i
      tr1(1,5)=exi
      tr1(2,5)=epxi
      tr1(3,5)=eyi
      tr1(4,5)=epyi
      tr1(5,1:4)=0.d0
      tr1(5,5)=1.d0
      return
      end
      
      subroutine teintp1(ip,f,tws,tr1,h1,tw0,ndims)
      use tfstk, only:ktfenanq
      use temw, only:etwiss2ri,ri,tinv6
      use ffs, only:xyth
      use tffitcode
      implicit none
      integer*4 ip,ip1,ndims
      real*8 tws(ntwissfun,-ndims:ndims),tr1(5,5),h1(4),
     $     rxi(6,6),rt(6,6)
      real*8 f,twf(ntwissfun),tw0(ntwissfun),
     $     dnx,dny,dnz,cx,cy,cz,sx,sy,sz
      logical*4 normal
      ip1=ip+1
      twf=f*tws(:,ip1)+(1.d0-f)*tws(:,ip)
      twf(mfitdetr)=twf(mfitr1)*twf(mfitr4)-twf(mfitr2)*twf(mfitr3)
      rxi=etwiss2ri(twf,normal)
c      call etwiss2ri(tw0,ri,normal)
      dnx=twf(mfitnx)-tw0(mfitnx)
      dny=twf(mfitny)-tw0(mfitny)
      dnz=twf(mfitnz)-tw0(mfitnz)
      cx=cos(dnx)
      sx=sin(dnx)
      cy=cos(dny)
      sy=sin(dny)
      cz=cos(dnz)
      sz=sin(dnz)
      rt(1,:)= cx*ri(1,:)+sx*ri(2,:)
      rt(2,:)=-sx*ri(1,:)+cx*ri(2,:)
      rt(3,:)= cy*ri(3,:)+sy*ri(4,:)
      rt(4,:)=-sy*ri(3,:)+cy*ri(4,:)
      rt(5,:)= cz*ri(5,:)+sz*ri(6,:)
      rt(6,:)=-sz*ri(5,:)+cz*ri(6,:)
      call tmov65(matmul(tinv6(rxi),rt),tr1)
c      call tinv6(rxi,rx)
c      call tmultr(rt,tinv6(rxi),6)
c      call tmov65(rt,tr1)
      h1=twf(mfitdx:mfitdpy)
      return
      end

      subroutine tediag(trans,cod,tws,i,ndims,stab,lfno)
      use tffitcode
      implicit none
      integer*4 i,ndims,lfno
      real*8 trans(6,6),cod(6),tws(ntwissfun,-ndims:ndims)
      real*8 r1,r2,r3,r4,amu,
     $     r1a,r2a,r3a,r4a,
     $     a11,a12,a22,a33,a34,a44,
     $     cosamux,sinamux,cosamuy,sinamuy,
     $     ax,bx,ay,by,amux,amuy
      logical*4 stab,nanq
      call qmdiag(
     $     trans(1,1),trans(1,2),trans(1,3),trans(1,4),
     $     trans(2,1),trans(2,2),trans(2,3),trans(2,4),
     $     trans(3,1),trans(3,2),trans(3,3),trans(3,4),
     $     trans(4,1),trans(4,2),trans(4,3),trans(4,4),
     $     r1,r2,r3,r4,amu,stab,nanq,lfno)
c     amu=sqrt(1.d0-r1*r4+r2*r3)
      r1a=r1/amu
      r2a=r2/amu
      r3a=r3/amu
      r4a=r4/amu
      a11=trans(1,1)-r4a*trans(3,1)+r2a*trans(4,1)
      a12=trans(1,2)-r4a*trans(3,2)+r2a*trans(4,2)
c      a21=trans(2,1)+r3a*trans(3,1)-r1a*trans(4,1)
      a22=trans(2,2)+r3a*trans(3,2)-r1a*trans(4,2)
      a33=trans(3,3)+r1a*trans(1,3)+r2a*trans(2,3)
      a34=trans(3,4)+r1a*trans(1,4)+r2a*trans(2,4)
c      a43=trans(4,3)+r3a*trans(1,3)+r4a*trans(2,3)
      a44=trans(4,4)+r3a*trans(1,4)+r4a*trans(2,4)
      cosamux=(a11+a22)*.5d0
      if(abs(cosamux) .gt. 1.d0)then
        stab=.false.
      endif
      sinamux=sign(sqrt((1.d0-cosamux)*(1.d0+cosamux)),
     $     a12)
      amux=atan2(sinamux,cosamux)
      if(sinamux .eq. 0.d0)then
        ax=0.d0
        bx=1.d0
      else
        ax=(a11-a22)*.5d0/sinamux
        bx=a12/sinamux
      endif
      tws(mfitax,i)=ax
      tws(mfitbx,i)=bx
      tws(mfitnx,i)=amux
      cosamuy=(a33+a44)*.5d0
      if(abs(cosamuy) .gt. 1.d0)then
        stab=.false.
      endif
      sinamuy=sign(sqrt((1.d0-cosamuy)*(1.d0+cosamuy)),
     $     a34)
      amuy=atan2(sinamuy,cosamuy)
      if(sinamuy .eq. 0.d0)then
        ay=0.d0
        by=1.d0
      else
        ay=(a33-a44)*.5d0/sinamuy
        by=a34/sinamuy
      endif
      tws(mfitay,i)=ay
      tws(mfitby,i)=by
      tws(mfitny,i)=amuy
      tws(mfitex:mfitepy,i)=trans(1:4,6)
      tws(mfitr1,i)=r1
      tws(mfitr2,i)=r2
      tws(mfitr3,i)=r3
      tws(mfitr4,i)=r4
      tws(mfitdetr,i)=r1*r4-r2*r3
      tws(mfitdx:mfitdpy,i)=cod(1:4)
      return
      end

      subroutine tmov65(a,b)
      implicit none
      real*8 a(6,6),b(5,5)
      b(1:4,1:4)=a(1:4,1:4)
      b(5,1:4)=0.d0
      b(1:4,5)=a(1:4,6)
      b(5,5)=a(6,6)
      return
      end

      subroutine ttimes(a,b,c,n)
      implicit none
      integer*4 n
      real*8 a(n),c(n),b
      c=a*b
      return
      end
