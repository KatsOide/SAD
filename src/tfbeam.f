      subroutine tfbeam(k,theta,beam)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_main,only:iaidx
      use temw,only:tmulbs
      implicit none
      integer*4 ,intent(in):: k
      integer*4 ip
      real*8 ,intent(inout):: beam(21)
      real*8 ,intent(in):: theta
      real*8 trans(4,5),trans1(6,6),cod(6)
      real*8 r,pfi,emxi,emyi,emmin,rgetgl1
      if(calc6d)then
        call tfsize(.false.)
        beam=beamsize(:,k)
      else
        ip=nlat*ndim+k
        beam=0.d0
        r=gammab(1)/gammab(k)
        if(gauss)then
          pfi=sizedp**2*r**2
        else
          pfi=sizedp**2/3.d0*r**2
        endif
        emxi=emx*r
        emyi=emy*r
        emmin=(emxi+emyi)*rgetgl1('MINCOUP')
        emxi=max(emmin,emxi)
        emyi=max(emmin,emyi)
        beam(iaidx(1,1))= twiss2(ip,2)*emxi+twiss2(ip,7)**2*pfi
c     write(*,*)'tfbeam ',emxi,pfi,twiss2(ip,2),twiss2(ip,7),
c     $     beam(1)
        beam(iaidx(1,2))=
     $       -twiss2(ip,1)*emxi+twiss2(ip,7)*twiss2(ip,8)*pfi
        beam(iaidx(1,3))= twiss2(ip,7)*twiss2(ip,9)*pfi
        beam(iaidx(1,4))= twiss2(ip,7)*twiss2(ip,10)*pfi
        beam(iaidx(2,2))=(1.d0+twiss2(ip,1)**2)/twiss2(ip,2)*emxi
     1       +twiss2(ip,8)**2*pfi
        beam(iaidx(2,3))= twiss2(ip,8)*twiss2(ip,9)*pfi
        beam(iaidx(2,4))= twiss2(ip,8)*twiss2(ip,10)*pfi
        beam(iaidx(3,3))= twiss2(ip,5)*emyi+twiss2(ip,9)**2*pfi
        beam(iaidx(3,4))=
     $       -twiss2(ip,4)*emyi+twiss2(ip,9)*twiss2(ip,10)*pfi
        beam(iaidx(4,4))=(1.d0+twiss2(ip,4)**2)/twiss2(ip,5)*emyi
     1       +twiss2(ip,10)**2*pfi
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
      endif
      return
      end

      subroutine tfbeamfrac(k,f,theta,beam)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use temw,only:iaez
      implicit none
      type (ffs_bound) fbound
      integer*4 k,kf
      real*8 beam(21),twisss(ntwissfun),ftwiss(ntwissfun)
      real*8 trans(6,12),cod(6),beam1(42),srot(3,9)
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
        if(calc6d)then
          fbound%lb=kf
          fbound%fb=0.d0
          fbound%le=kf
          fbound%fe=fr
          beam1(1:21)=beamsize(:,k)
          beam1(22:42)=0.d0
          call tinitr12(trans)
          cod=twiss(k,0,mfitdx:mfitddp)
          irad=12
          call tturne0(trans,cod,beam1,srot,fbound,
     $     iaez,.false.,.false.,.false.)
          beam=beam1(1:21)
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
      endif
      return
      end

      subroutine tfsize(force)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_flag
      use tmacro
      use temw, only:nparams,tfinibeam,iaez,beamplt
      use tffitcode
      implicit none
      integer*4 i
      real*8 params(nparams),trans(6,12),cod(6),beam(42),btr(21,21)
      logical*4 ,intent(in):: force
      logical*4 stab,bpl
      call ffs_init_sizep
      if(calc6d)then
        if(force .or. modesize .ne. 6)then
          beamin=0.d0
          if(trpt)then
            beamin=tfinibeam(1)
          endif
          bpl=beamplt
          beamplt=.true.
          cod=twiss(1,0,mfitdx:mfitddp)
          call temit(trans,cod,beam,btr,
     $         .true.,iaez,.true.,params,stab,0)
          beamplt=bpl
        endif
        modesize=6
      else
        if(force .or. modesize .ne. 4)then
          do i=1,nlat
            call tfbeam(i,0.d0,beamsize(:,i))
          enddo
        endif
        modesize=4
      endif
      return
      end
