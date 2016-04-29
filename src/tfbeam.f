      subroutine tfbeam(twiss,gammab,k,theta,beam)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k,ip,i,m,n,ia
      real*8 twiss(nlat*(2*ndim+1),ntwissfun),gammab(nlat)
      real*8 trans(4,5),trans1(6,12),cod(6),beam(21)
      real*8 theta,r,pfi,emxi,emyi
      ia(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
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
      beam(ia(1,1))= twiss(ip,2)*emxi+twiss(ip,7)**2*pfi
c      write(*,*)'tfbeam ',emxi,pfi,twiss(ip,2),twiss(ip,7),
c     $     beam(1)
      beam(ia(1,2))=-twiss(ip,1)*emxi+twiss(ip,7)*twiss(ip,8)*pfi
      beam(ia(1,3))= twiss(ip,7)*twiss(ip,9)*pfi
      beam(ia(1,4))= twiss(ip,7)*twiss(ip,10)*pfi
      beam(ia(2,2))=(1.d0+twiss(ip,1)**2)/twiss(ip,2)*emxi
     1             +twiss(ip,8)**2*pfi
      beam(ia(2,3))= twiss(ip,8)*twiss(ip,9)*pfi
      beam(ia(2,4))= twiss(ip,8)*twiss(ip,10)*pfi
      beam(ia(3,3))= twiss(ip,5)*emyi+twiss(ip,9)**2*pfi
      beam(ia(3,4))=-twiss(ip,4)*emyi+twiss(ip,9)*twiss(ip,10)*pfi
      beam(ia(4,4))=(1.d0+twiss(ip,4)**2)/twiss(ip,5)*emyi
     1             +twiss(ip,10)**2*pfi
      beam(ia(1,6))= twiss(ip,7)*pfi
      beam(ia(2,6))= twiss(ip,8)*pfi
      beam(ia(3,6))= twiss(ip,9)*pfi
      beam(ia(4,6))= twiss(ip,10)*pfi
      beam(ia(6,6))= pfi
      cod(6)=0.d0
      call qtent(trans,cod,twiss,k,0,.false.)
      call qchg(trans,cod,0.d0,0.d0,theta,.true.)
      call tclr(trans1,36)
      do 10 i=1,4
        trans1(1,i)=trans(1,i)
        trans1(2,i)=trans(2,i)
        trans1(3,i)=trans(3,i)
        trans1(4,i)=trans(4,i)
10    continue
      trans1(5,5)=1.d0
      trans1(6,6)=1.d0
      irad=12
      call tmulbs(beam,trans1,.false.,.false.)
      return
      end

      subroutine tfbeamfrac(twiss,gammab,k,f,theta,beam)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k,i,kf
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
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
        call tfbeam(twiss,gammab,kf,theta,beam)
      else
        do i=1,ntwissfun
          twisss(i)=twiss(kf+1,0,i)
        enddo
        gb0=gammab(kf+1)
        call qtwissfrac(ftwiss,ilist(1,ilattp+1),twiss,gammab,
     $       kf,fr,over)
        do i=1,ntwissfun
          twiss(kf+1,0,i)=ftwiss(i)
        enddo
        gammab(kf+1)=(gb0-gammab(kf))*fr+gammab(kf)
        call tfbeam(twiss,gammab,kf+1,theta,beam)
        gammab(kf+1)=gb0
        do i=1,ntwissfun
          twiss(kf+1,0,i)=twisss(i)
        enddo
      endif
      return
      end

      subroutine tfsizefrac(size,k,f,beam)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k,nvar
      real*8 size(21,*),f,beam(42),trans(6,12),cod(6),vsave(256)
      logical*4 chg,sol,cp0,int0
      call tmov(size(1,k),beam,21)
      call tmov(size(1,k),beam(22),21)
      if(f .ne. 0.d0)then
        call qfracsave(ilist(1,ilattp+k),vsave,nvar,ideal,.true.)
        call qfraccomp(ilist(1,ilattp+k),0.d0,f,ideal,chg)
        if(.not. chg)then
          return
        endif
        cod(1)=gettwiss(mfitdx,k)
        cod(2)=gettwiss(mfitdpx,k)
        cod(3)=gettwiss(mfitdy,k)
        cod(4)=gettwiss(mfitdpy,k)
        cod(5)=gettwiss(mfitdz,k)
        cod(6)=gettwiss(mfitddp,k)
        call tinitr(trans)
        call tclr(trans(1,7),36)
        sol=.false.
        cp0=codplt
        codplt=.false.
        int0=intra
        calint=.false.
        call tturne1(ilist(1,ilattp+1),trans,cod,beam,
     $       rlist(iftwis),size,rlist(ifgamm),
     $       0,0,0,ndim,.false.,sol,.false.,k,k)
        codplt=cp0
        calint=int0
        call qfracsave(ilist(1,ilattp+k),vsave,nvar,ideal,.false.)
      endif
      return
      end
