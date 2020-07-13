      subroutine tfbeam(k,theta,beam)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_main,only:iaidx
      implicit none
      integer*4 k,ip
      real*8 trans(4,5),trans1(6,12),cod(6),beam(21)
      real*8 theta,r,pfi,emxi,emyi
c      iaidx(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      ip=nlat*ndim+k
      beam=0.d0
      r=gammab(1)/gammab(k)
      if(gauss)then
        pfi=dpmax**2*r**2
      else
        pfi=dpmax**2/3.d0*r**2
      endif
      emxi=emx*r
      emyi=emy*r
      beam(iaidx(1,1))= twiss2(ip,2)*emxi+twiss2(ip,7)**2*pfi
c      write(*,*)'tfbeam ',emxi,pfi,twiss2(ip,2),twiss2(ip,7),
c     $     beam(1)
      beam(iaidx(1,2))=-twiss2(ip,1)*emxi+twiss2(ip,7)*twiss2(ip,8)*pfi
      beam(iaidx(1,3))= twiss2(ip,7)*twiss2(ip,9)*pfi
      beam(iaidx(1,4))= twiss2(ip,7)*twiss2(ip,10)*pfi
      beam(iaidx(2,2))=(1.d0+twiss2(ip,1)**2)/twiss2(ip,2)*emxi
     1             +twiss2(ip,8)**2*pfi
      beam(iaidx(2,3))= twiss2(ip,8)*twiss2(ip,9)*pfi
      beam(iaidx(2,4))= twiss2(ip,8)*twiss2(ip,10)*pfi
      beam(iaidx(3,3))= twiss2(ip,5)*emyi+twiss2(ip,9)**2*pfi
      beam(iaidx(3,4))=-twiss2(ip,4)*emyi+twiss2(ip,9)*twiss2(ip,10)*pfi
      beam(iaidx(4,4))=(1.d0+twiss2(ip,4)**2)/twiss2(ip,5)*emyi
     1             +twiss2(ip,10)**2*pfi
      beam(iaidx(1,6))= twiss2(ip,7)*pfi
      beam(iaidx(2,6))= twiss2(ip,8)*pfi
      beam(iaidx(3,6))= twiss2(ip,9)*pfi
      beam(iaidx(4,6))= twiss2(ip,10)*pfi
      beam(iaidx(6,6))= pfi
      cod(6)=0.d0
      call qtent(trans,cod,k,0,.false.)
      call qchg(trans,cod,0.d0,0.d0,theta,.true.)
      trans1(:,5:6)=0.d0
      trans1(5:6,:)=0.d0
      trans1(1:4,1:4)=trans(1:4,1:4)
      trans1(5,5)=1.d0
      trans1(6,6)=1.d0
      irad=12
      call tmulbs(beam,trans1,.false.)
      return
      end

      subroutine tfbeamfrac(k,f,theta,beam)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 k,kf
      real*8 beam(21),twisss(ntwissfun),ftwiss(ntwissfun)
      real*8 f,fr,theta,gb0
      logical*4 over
      kf=int(f)
      fr=f-kf
      kf=min(nlat,max(1,k+kf))
      if(kf .eq. nlat .or. kf .eq. 1)then
        fr=0.d0
      endif
      if(fr .eq. 0.d0)then
        call tfbeam(kf,theta,beam)
      else
        twisss=twiss(kf+1,0,1:ntwissfun)
        gb0=gammab(kf+1)
        call qtwissfrac(ftwiss,kf,fr,over)
        twiss(kf+1,0,1:ntwissfun)=ftwiss
        gammab(kf+1)=(gb0-gammab(kf))*fr+gammab(kf)
        call tfbeam(kf+1,theta,beam)
        gammab(kf+1)=gb0
        twiss(kf+1,0,1:ntwissfun)=twisss
      endif
      return
      end

      subroutine tfsizefrac(k,f,beam)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_comp), pointer ::cmp
      type (sad_descriptor) :: dsave(kwMAX)
      integer*4 k,nvar,le,itfdownlevel,irtc
      real*8 f,beam(42),trans(6,12),cod(6),srot(3,9)
      logical*4 chg,sol,cp0,int0
      if(.not. updatesize .or. sizedp .ne. dpmax)then
        call tfsize
        updatesize=.true.
      endif
      beam(1:21)=beamsize(:,k)
      beam(22:42)=beam(1:21)
      if(f .ne. 0.d0)then
        levele=levele+1
        call qfracsave(k,dsave,nvar,.true.)
        call compelc(k,cmp)
        call qfracseg(cmp,cmp,0.d0,f,chg,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(k,irtc)
        else
          if(.not. chg)then
            le=itfdownlevel()
            return
          endif
          cod(1)=gettwiss(mfitdx,k)
          cod(2)=gettwiss(mfitdpx,k)
          cod(3)=gettwiss(mfitdy,k)
          cod(4)=gettwiss(mfitdpy,k)
          cod(5)=gettwiss(mfitdz,k)
          cod(6)=gettwiss(mfitddp,k)
          call tinitr(trans)
          trans(:,7:12)=0.d0
          sol=.false.
          cp0=codplt
          codplt=.false.
          int0=intra
          calint=.false.
          call tturne1(trans,cod,beam,srot,
     $         i00,i00,i00,0,.false.,sol,.false.,
     $         .false.,k,k)
          codplt=cp0
          calint=int0
        endif
        call qfracsave(k,dsave,nvar,.false.)
        le=itfdownlevel()
      endif
      return
      end

      subroutine tfsize
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 i
      call ffs_init_sizep
      do i=1,nlat
        call tfbeam(i,0.d0,beamsize(:,i))
      enddo
      updatesize=.true.
      return
      end
