      module touschek_table
      implicit none
      private

      public :: initialize_tampl

      integer(8), public:: itoul = 0
      integer, public :: id = 0

      integer, public, parameter :: ntouckx = 34, ntouckz = 25
      integer, public, parameter :: ntouckl = 120

c     ntouckx > ntouckz

c     Table of amplitude(index, axis[1->x, 2->y, 3->z])
      real(8), public :: tampl(ntouckx, 3)

c     Table of loss-rate
      real(8), public :: touckl(ntouckl)
      real(8), public :: touckm(ntouckz, ntouckx, 3)
      real(8), public, allocatable :: toucke(:,:)

      contains

      subroutine initialize_tampl()
      implicit none
      tampl = 0.d0

      tampl( 1,1)=1.d0
      tampl( 2,1)=2.d0
      tampl( 3,1)=3.d0
      tampl( 4,1)=4.d0

      tampl( 1+4,1)=5.d0
      tampl( 2+4,1)=6.d0
      tampl( 3+4,1)=7.d0
      tampl( 4+4,1)=8.d0
      tampl( 5+4,1)=9.d0

      tampl( 1+5+4,1)=10.d0
      tampl( 2+5+4,1)=11.d0
      tampl( 3+5+4,1)=13.d0
      tampl( 4+5+4,1)=15.d0
      tampl( 5+5+4,1)=17.d0
      tampl( 6+5+4,1)=19.d0
      tampl( 7+5+4,1)=21.d0
      tampl( 8+5+4,1)=24.d0
      tampl( 9+5+4,1)=28.d0
      tampl(10+5+4,1)=32.d0
      tampl(11+5+4,1)=36.d0
      tampl(12+5+4,1)=41.d0
      tampl(13+5+4,1)=46.d0
      tampl(14+5+4,1)=53.d0
      tampl(15+5+4,1)=60.d0
      tampl(16+5+4,1)=68.d0
      tampl(17+5+4,1)=77.d0
      tampl(18+5+4,1)=88.d0

c      do i = 28, ntouckx
        tampl(28:ntouckx, 1) = tampl(10:ntouckx-18, 1) * 10.d0
c      enddo

      tampl(1:ntouckx,2) = tampl(1:ntouckx,1)

      tampl(7:ntouckx-3,3) = tampl(10:ntouckx,1)
      tampl( 1,3)= 4.6d0
      tampl( 2,3)= 5.3d0
      tampl( 3,3)= 6.d0
      tampl( 4,3)= 6.8d0
      tampl( 5,3)= 7.7d0
      tampl( 6,3)= 8.8d0

      return
      end subroutine initialize_tampl

      end module touschek_table

      subroutine temit(trans,cod,beam,btr,
     $     calem0,iae,plot,params,stab,lfno)
      use tfstk
      use temw
      use ffs_flag
      use ffs_pointer
      use tmacro
      use tffitcode
      use tspin, only:spnorm,sremit,srotinit
      use eigen
      use macmath
      implicit none
      type (iaemit),intent(in):: iae
      real*8 ,parameter ::conv=1.d-12
      integer*4 ,intent(in):: lfno
      integer*4 it,i,iret
      real*8 ,intent(out):: trans(6,12),cod(6),beam(42),
     $     btr(21,21),params(nparams)
      real*8 gr(6),dc
      complex*16 cd(6),ceig(6)
      real*8 emitn(21),emitp(21),beamn(21),beamp(21),
     $     srot(3,9),srot1(3,3),sps(3,3),spm(3,3),rx(6,6),rd(6,6)
c     $     ,tw(ntwissfun)
      logical*4 ,intent(out):: stab
      logical*4 ,intent(in):: plot,calem0
      logical*4 pri,calem,calcodr,rt,radpol0,intpri,econv,
     $     inical,postcal
      it=0
      trf0=0.d0
      vcalpha=1.d0
      demin=1.d100
      beam(1:21)=beamin
c      write(*,*)'temit ',beamin(iaidx(3,3)),
c     $     beamin(iaidx(3,4)),beamin(iaidx(4,4))
      beam(22:42)=0.d0
      params=0.d0
      ceig0=(0.d0,0.d0)
      call tsetdvfs
      gr(1:5:2)=1.d0
      gr(2:6:2)=gammab(nlat)/gammab(1)
      emx0=0.d0
      emy0=0.d0
      emz0=0.d0
      emxmin=1.d-200
      emymin=1.d-200
      emzmin=1.d-200
      emxmax=1.d200
      emymax=1.d200
      emzmax=1.d200
      irad=6
      calint=.false.
      epi=.false.
      econv=.not. (intra .or. wspac)
      calem=calem0 .and. .not. trpt
      radpol0=radpol
      intpri=.true.
      caltouck=.false.
      calcodr=.not. trpt .and. calcod
      pri=lfno .gt. 0
      rt=radcod .and. radtaper
      inical=.true.
      fndcod=.true.
      postcal=.true.
      dc=1.d100
c      tw=tfetwiss(ri,codin,.true.)
c      write(*,'(a,1p6g15.7)')'temit-etwiss-ent0 ',
c     $     tw(mfitax:mfitny)/[1d0,1d0,m_2pi,1d0,1d0,m_2pi]
      if(calcodr)then
        call tcod(trans,cod,beam,.false.,fndcod)
        call limitnan(cod,-1.d10,1.d10)
        if(.not. fndcod)then
          write(lfno,*)'???-Emittance[]-closed orbit not found.'
        endif
        codin=cod
      endif
      irad=12
      if(.not. calem .and. econv)then
        irad=6
      endif
      converge: do while(.true.)
        if(inical)then
          cod=codin
          beam(1:21)=beamin
          call tinitr12(trans)
          call srotinit(srot)
          call tsetr0(trans,cod,0.d0,0.d0)
c          write(*,*)'temit-inical ',calint,intra,econv
          if(econv)then
            call tturne(trans,cod,beam,srot,iae,
     1           plot,.false.,rt,.false.)
            postcal=.false.
          else
            call tturne(trans,cod,beam,srot,iaez,
     1           .false.,.false.,rt,.false.)
          endif
          if(trpt)then
            do i=1,12
              trans(:,i)=gr*trans(:,i)
            enddo
          endif
          call limitnan(cod,-1.d10,1.d10)
          irad=12
          call tecalc(trans,cod,beam,beamn,beamp,
     $         emitn,emitp,btr,srot,srot1,
     $         rx,rd,params,ceig,dc,cd,stab,calem)
          inical=.false.
        endif
        if(pri)then
          call temitprint(lfno,intpri,
     $         rx,trans(:,7:12),cod,rd,beamn,beamp,emitn,emitp,
     $         srot,srot1,spm,sps,ceig,cd,params,calem)
          intpri=.false.
        endif
        if(calcodr .and. .not. stab .and. intra)then
          write(lfno,*)
     $         'Skip intrabeam because of unstable.'
        elseif(.not. econv)then
          call tintraconv(lfno,it,emitn,beam,dc,pri,iret)
c          if(plot)then
c            write(*,'(a,i5,1p6g15.7)')'temit-intraconv-end ',
c     $           iret,dc,emx0,emy0,eemx,eemy
c          endif
          select case (iret)
          case (2)
            inical=.true.
            cycle converge
          case (3)
            econv=.true.
            inical=.true.
            cycle converge
          case (4)
            econv=.true.
            cycle converge
          end select
        endif
        exit converge
      enddo converge
      if(plot .and. postcal)then
        call tinitr12(trans)
        cod=codin
        if(trpt)then
          beam(1:21)=beamin
        else
          beam(1:21)=beam(22:42)
        endif
        beam(22:42)=0.d0
        call srotinit(srot)
        call tturne(trans,cod,beam,srot,iae,
     $       .true.,.false.,rt,.false.)
      endif
      if(iae%iamat .eq. 0)then
        if(beamplt .and. plot .and. charge .lt. 0.d0)then
          beamsize=-beamsize
        endif
      else
        call setiamat(iae%iamat,ri,codin,beamp,beam,trans)
      endif
      return
      end

      subroutine tecalc(trans,cod,beam,beamn,beamp,
     $     emitn,emitp,btr,srot,srot1,
     $     rx,rd,params,ceig,dc,cd,stab,calem)
      use tfstk
      use temw
      use ffs_flag
      use ffs_pointer
      use tmacro
      use tffitcode
      use tspin, only:spnorm,sremit
      use eigen
      use macmath
      implicit none
      real*8 ,parameter:: conv=1.d-12
      integer*4 lfno,i,j,k,k1,k2,k3,m,l,n
      real*8 ,intent(inout):: beam(42),srot1(3,3),srot(3,9)
      real*8 ,intent(in):: trans(6,12),cod(6)
      real*8 ,intent(out):: dc
      real*8 sdamp,sqr2,bb,bbv(21),sr,rgetgl1,
     $     tune,ab(6),emxe,emye,emze,rirx(6,6)
      complex*16 ,intent(out):: cd(6)
      complex*16 dceig(6),cc(6),ceig(6)
      real*8 ,intent(out):: btr(21,21),emitn(21),emitp(21),beamn(21),
     1     beamp(21),params(nparams),rx(6,6),rd(6,6)
      real*8 tw(ntwissfun),sps(3,3),spm(3,3),rdb(6,12)
      logical*4 ,intent(in):: calem
      logical*4 ,intent(out):: stab
      rx=trans(:,1:6)
      if(.not. rfsw)then
        rx(6,1)=0.d0
        rx(6,2)=0.d0
        rx(6,3)=0.d0
        rx(6,4)=0.d0
        rx(6,5)=0.d0
        rx(6,6)=1.d0
      endif
      emxe=rgetgl1('EMITXE')
      emye=rgetgl1('EMITYE')
      emze=rgetgl1('EMITZE')
c     write(*,'(a/,6(1p6g15.7/))')'trans: ',(trans(i,1:6),i=1,6)
      if(trpt)then
        params(iptwiss:iptwiss+ntwissfun-1)=tfetwiss(ri,codin,.true.)
        rirx=matmul(ri,tinv6(rx))
        tw=tfetwiss(rirx,cod,.true.)
c      write(*,'(a,1p6g15.7)')'temit-etwiss-exit ',
c     $       tw(mfitax:mfitny)/[1d0,1d0,m_2pi,1d0,1d0,m_2pi]
        ceig(1)=exp(dcmplx(0.d0,tw(mfitnx)))
        ceig(3)=exp(dcmplx(0.d0,tw(mfitny)))
        ceig(5)=exp(dcmplx(0.d0,tw(mfitnz)))
        ceig(2:6:2)=conjg(ceig(1:5:2))
      else
        r=rx
        call teigenc(r,ri,ceig,6,6)
        call tnorm(r,ceig,lfno)
        r=tsymp(r)
        ri=tinv6(r)
      endif
c      tw=tfetwiss(ri,codin,.true.)
c      write(*,'(a,1p6g15.7)')'temit-etwiss-ent ',
c     $       tw(mfitax:mfitny)/[1d0,1d0,m_2pi,1d0,1d0,m_2pi]
      dceig=ceig-ceig0
      dc=sum(abs(dceig))
      cc(1:5:2)=ceig(1:5:2)
      cc(2:6:2)=conjg(cc(1:5:2))
      ceig0=ceig
      cd(4:6)=log(cc(1:5:2))
      if(vceff .ne. 0.d0)then
        phirf=asin(u0*pgev/vceff)
        heff=wrfeff*cveloc/omega0
      else
        phirf=0.d0
        heff=0.d0
      endif
      synchm=.not. trpt .and. rfsw .and. imag(cd(6)) .ne. 0.d0
      if(synchm)then
        if(wrfeff .ne. 0.d0)then
          alphap=-imag(cd(6))*abs(imag(cd(6)))/(c*m_2pi/omega0)
     $         /(dvcacc/pgev)
        else
          alphap=0.d0
        endif
        omegaz=abs(imag(cd(6)))*omega0/m_2pi
      else
        alphap=-rx(5,6)/m_2pi/cveloc/p0*h0*omega0
        if(trpt)then
          omegaz=0.d0
        else
          omegaz=sqrt(abs(alphap*m_2pi*heff*vceff/pgev*cos(phirf)))
     $         *omega0/m_2pi
        endif
      endif
      if(vceff .ne. 0.d0)then
        bh=sqrt(abs(vceff/pi/abs(alphap)/heff/pgev*
     1       (2.d0*cos(phirf)-(pi-2.d0*phirf)*u0*pgev/vceff)))
      else
        bh=0.d0
      endif
      call setparams(params,cod)
      stab=(abs(dble(cd(4))) .lt. 1.d-6
     $     .and. abs(dble(cd(5))) .lt. 1.d-6
     1     .and. abs(dble(cd(6))) .lt. 1.d-6) .and. fndcod
      params(ipnx:ipnz)=imag(cd(4:6))/m_2pi
      params(iptwiss:iptwiss+ntwissfun-1)=tfetwiss(ri,cod,.true.)
      if(.not. calem .and. .not. trpt)then
        return
      endif
      rd=matmul(ri,matmul(trans(:,7:12),r))
      do i=1,5,2
        cd(int(i/2)+1)=dcmplx((rd(i,i)+rd(i+1,i+1))*.5d0,
     1       (rd(i,i+1)-rd(i+1,i))*.5d0)/cc(i)
      enddo
      if(dble(cd(1)) .ne. 0.d0)then
        taurdx=-m_2pi/omega0/dble(cd(1))
      else
        taurdx=0.d0
      endif
      if(dble(cd(2)) .ne. 0.d0)then
        taurdy=-m_2pi/omega0/dble(cd(2))
      else
        taurdy=0.d0
      endif
      if(dble(cd(3)) .ne. 0.d0)then
        taurdz=-m_2pi/omega0/dble(cd(3))
      else
        taurdz=0.d0
      endif
      params(ipdampx:ipdampz)=dble(cd(1:3))
      params(ipdnux:ipdnuz)=imag(cd(1:3))
      sr=params(ipdampx)+params(ipdampy)+params(ipdampz)
      if(sr .ne. 0.d0)then
        sr=4.d0/sr
        params(ipjx:ipjz)=params(ipdampx:ipdampz)*sr
      else
        sr=0.d0
        params(ipjx:ipjz)=0.d0
      endif
      beamp(1:21)=beam(1:21)
      call tmulbs(beam,ri,.false.)
      beamn(1:21)=beam(1:21)
      if(trpt)then
        emitn=beam(1:21)
      else
        rdb(:,1:6)=rd
        if(.not. synchm)then
c          do i=1,6
          beam(iaidx(5,1:6))=0.d0
          rdb(1:6,5)=0.d0
          rdb(5,1:6)=0.d0
c          enddo
        endif
        btr=0.d0
        rdb(:,7:12)=0.d0
        do i=1,5,2
          tune=imag(cd(int(i/2)+4))
          rdb(i  ,i+6)= cos(tune)
          rdb(i  ,i+7)= sin(tune)
          rdb(i+1,i+6)=-sin(tune)
          rdb(i+1,i+7)= cos(tune)
        enddo
        do i=1,6
          do j=1,i
            k=iaidx(i  ,j  )
            do m=1,6
              do n=1,6
                l=iaidx(m,n)
                btr(k,l)=btr(k,l)-(rdb(i,m)+rdb(i,m+6))*
     1               (rdb(j,n)+rdb(j,n+6))
              enddo
            enddo
          enddo
        enddo
        sqr2=sqrt(.5d0)
        do i=1,5,2
          k1=iaidx(i  ,i  )
          k2=iaidx(i+1,i+1)
          k3=iaidx(i  ,i+1)
c     do j=1,21
          bbv=btr(k1,1:21)
          btr(k1,1:21)=( bbv+btr(k2,1:21))*sqr2
          btr(k2,1:21)=(-bbv+btr(k2,1:21))*sqr2
          btr(k3,1:21)=btr(k3,1:21)/sqr2
c     enddo
          bb=beam(k1)
          beam(k1)=( bb+beam(k2))*sqr2
          beam(k2)=(-bb+beam(k2))*sqr2
          beam(k3)=beam(k3)/sqr2
c     do j=1,21
          bbv=btr(1:21,k1)
          btr(1:21,k1)=( bbv+btr(1:21,k2))*sqr2
          btr(1:21,k2)=(-bbv+btr(1:21,k2))*sqr2
          btr(1:21,k3)=btr(1:21,k3)*sqr2
c     enddo
        enddo
        emitn=0.d0
        do i=1,21
          btr(i,i)=btr(i,i)+1.d0
        enddo
        do i=1,5,2
          k=iaidx(i,i)
          btr(k,k)=-(rdb(i  ,i)**2+rdb(i  ,i+1)**2+
     1         rdb(i+1,i)**2+rdb(i+1,i+1)**2)*.5d0-
     1         rdb(i,i+6)*(rdb(i,i)+rdb(i+1,i+1))-
     1         rdb(i,i+7)*(rdb(i,i+1)-rdb(i+1,i))
        enddo
        if(.not. synchm)then
          btr(15,21)=0.d0
          btr(15,15)=btr(15,15)*2.d0
          btr(21,15)=-btr(21,21)
          beam(21)=0.d0
        endif
        do i=1,5,2
          k1=iaidx(i,i)
          k2=iaidx(i+1,i+1)
          if(btr(k2,k2) .ne. 0.d0 .and. btr(k1,k1) .ne. 0.d0)then
            ab(i)=sqrt(abs(btr(k1,k1)/btr(k2,k2)))
c     do j=1,21
            btr(k1,1:21)=btr(k1,1:21)/ab(i)
            btr(k2,1:21)=btr(k2,1:21)*ab(i)
c     enddo
            beam(k1)=beam(k1)/ab(i)
            beam(k2)=beam(k2)*ab(i)
          else
            ab(i)=1.d0
          endif
        enddo
        call tsolva(btr,beam,emitn,21,21,21,1d-8)
        do i=1,5,2
          k1=iaidx(i,i)
          k2=iaidx(i+1,i+1)
          bb=emitn(k1)
          emitn(k1          )=(bb-emitn(k2))*sqr2
          emitn(k2          )=(bb+emitn(k2))*sqr2
          emitn(iaidx(i  ,i+1))=emitn(iaidx(i  ,i+1))*sqr2
        enddo
        emitn(iaidx(1,1))=sign(max(abs(emitn(iaidx(1,1))),emxe),
     $       emitn(iaidx(1,1)))
        emitn(iaidx(2,2))=sign(max(abs(emitn(iaidx(2,2))),emxe),
     $       emitn(iaidx(2,2)))
        emitn(iaidx(3,3))=sign(max(abs(emitn(iaidx(3,3))),emye),
     $       emitn(iaidx(3,3)))
        emitn(iaidx(4,4))=sign(max(abs(emitn(iaidx(4,4))),emye),
     $       emitn(iaidx(4,4)))
        emitn(iaidx(5,5))=sign(max(abs(emitn(iaidx(5,5))),emze),
     $       emitn(iaidx(5,5)))
        emitn(iaidx(6,6))=sign(max(abs(emitn(iaidx(6,6))),emze),
     $       emitn(iaidx(6,6)))
      endif
      if(.not. epi)then
        eemx= sign(sqrt(abs(emitn(iaidx(1,1))*emitn(iaidx(2,2))
     $       -emitn(iaidx(1,2))**2)),emitn(iaidx(2,2))*charge)
        eemy= sign(sqrt(abs(emitn(iaidx(3,3))*emitn(iaidx(4,4))
     $       -emitn(iaidx(3,4))**2)),emitn(iaidx(4,4))*charge)
        eemz= sign(sqrt(abs(emitn(iaidx(5,5))*emitn(iaidx(6,6))
     $       -emitn(iaidx(5,6))**2)),emitn(iaidx(6,6))*charge)
      endif
      emitp(1:21)=emitn
      call tmulbs(emitp,r,.false.)
      sige=sqrt(abs(emitp(21)))
      if(synchm .or. trpt)then
        sigz=sqrt(abs(emitp(15)))
      else
        if(omegaz .ne. 0.d0)then
          sigz=abs(alphap)*sige*cveloc*p0/h0/omegaz
        else
          sigz=0.d0
        endif
        eemz=sigz*sige
      endif
      params(ipemx:ipemz)=(/eemx,eemy,eemz/)
      params(ipsige)=sige
      params(ipsigz)=sigz
      params(ipnnup)=h0*gspin
      if(calpol)then
        call spnorm(srot,sps,spinmu,sdamp)
        params(iptaup)=1.d0/sdamp/params(iprevf)
        srot1=srot(:,1:3)
        params(ipnup)=spinmu/m_2pi
        call sremit(srot,sps,params,beamn,sdamp,spm,equpol)
        params(ipequpol)=equpol(3)
        params(ipequpol2:ipequpol6)=equpol
        params(ippolx:ippolz)=sps(:,1)
      endif
      if(calem .and. .not. trpt)then
        call rsetgl1('EMITX',eemx)
        call rsetgl1('EMITY',eemy)
        call rsetgl1('EMITZ',eemz)
        call rsetgl1('SIGE',sige)
        call rsetgl1('SIGZ',sigz)
      endif
      return
      end

      subroutine tintraconv(lfno,it,emitn,beam,dc,pri,iret)
      use tfstk
      use ffs_flag
      use touschek_table
      use tmacro
      use sad_main , only:iaidx
      use temw
      implicit none
      type (sad_dlist), pointer :: klx1
      type (sad_dlist), pointer :: klx2,klx
      type (sad_rlist), pointer :: klx1d,klx1l
      integer*4 ,parameter ::itmax=100
      real*8 ,parameter ::resib=3.d-6,dcmin=1.d-6
      integer*8 kax,kax1,kax1d,kax1l,kax2
      integer*4 ,intent(in):: lfno
      integer*4 ,intent(inout):: it
      integer*4 ,intent(out):: iret
      integer*4 i,iii,k,j,m
      real*8 ,intent(in):: dc
      real*8 ,intent(inout):: emitn(21),beam(42)
      real*8 trans1(6,6),rx,ry,rz,emx1,emy1,emz1,emmin,
     $     de,tf,tt,eintrb,rr
      logical*4 ,intent(out):: pri
      character*11 autofg,vout(5)
      emx1=eemx
      emy1=eemy
      emz1=eemz
      if(trpt)then
        call rsetgl1('EMITXC',eemx)
        call rsetgl1('EMITYC',eemy)
        call rsetgl1('EMITZC',eemz)
        de=0.d0
        calint=.true.
        iret=3
        return
      else
        emmin=(eemx+eemy)*coumin
        eemx=max(emmin,eemx)
        eemy=max(emmin,eemy)
        eemz=max(emz0*0.1d0,eemz)
        if(eemx .le. 0.d0 .or. eemy .le. 0.d0 .or. eemz .le. 0.d0)then
          write(lfno,'(2a,/,a,1p3g15.7/)')
     $         ' Negative emittance, ',
     $         'No intrabeam/space charge calculation:',
     $         'eem(x,y,z) =',eemx,eemy,eemz
          it=itmax+1
        else
          eemx=min(emxmax,max(emxmin,eemx))
          eemy=min(emymax,max(emymin,eemy))
          eemz=min(emzmax,max(emzmin,eemz))
          de=(1.d0-emx0/eemx)**2+
     1         (1.d0-emy0/eemy)**2+(1.d0-emz0/eemz)**2
          demin=min(de,demin)
          if(it .ge. 20)then
            if(it .eq. 20)then
              write(*,*)' Poor convergence... '
              write(*,*)
     $     '     EMITX          EMITY          EMITZ           conv'
            endif
            write(*,'(1p4G15.7)')eemx,eemy,eemz,de
          endif
        endif
      endif
 7301 if(it .gt. 1 .and. dc .lt. dcmin
     $     .and. de .lt. resib .or. it .gt. itmax)then
c        write(*,*)'tintraconv ',it,dc,de
        pri=lfno .gt. 0
        if(.not. trpt)then
          if(de .ge. resib .or. dc .ge. dcmin)then
            write(*,*)' Intrabeam/space charge convergence failed.'
            write(*,*)
          elseif(intra .and. .not. caltouck)then
            de=resib*1.01d0
            go to 7301
          endif
          if(itoul .eq. 0)then
            itoul=ktfsymbolz('TouschekTable',13)-4
          endif
          tf=rclassic**2*pbunch*sqrt(pi)/h0*omega0/2.d0/pi/p0*h0
          if(caltouck)then
            id=id+1
            kax=ktadaloc(0,4,klx)
            kax1=ktadaloc(0,2,klx1)
            kax1d=ktavaloc(0,ntouckl,klx1d)
            kax1l=ktavaloc(0,ntouckl,klx1l)
            do i=1,ntouckl
              klx1d%rbody(i)=(i+1)*2.d-3
              klx1l%rbody(i)=touckl(i)*tf
              do j=1,nlat
c factor: tf for toucke(#dp/p0,#element)
                toucke(i,j)=toucke(i,j)*tf
              enddo
            enddo
            klx1%dbody(1)%k=ktflist+kax1d
            klx1%dbody(2)%k=ktflist+kax1l
            kax2=ktadaloc(0,3,klx2)
            klx2%dbody(1)=
     $           dtfcopy1(kxm2l(tampl,ntouckx,3,ntouckx,.true.))
            do i=1,ntouckx
c              do 1:ntouckz=1,ntouckz
              touckm(1:ntouckz,i,1)=touckm(1:ntouckz,i,1)*tf
              touckm(1:ntouckz,i,2)=touckm(1:ntouckz,i,2)*tf
c              enddo
            enddo
            klx2%dbody(2)=dtfcopy1(kxm2l(touckm(1,1,1),
     $           ntouckz,ntouckx,ntouckz,.true.))
            klx2%dbody(3)=dtfcopy1(kxm2l(touckm(1,1,2),
     $           ntouckz,ntouckx,ntouckz,.true.))
c set list for toucke
            klx%rbody(1)=dble(id)
            klx%dbody(2)%k=ktflist+kax1
            klx%dbody(3)%k=ktflist+kax2
            klx%dbody(4)=dtfcopy1(
     $           kxm2l(toucke,ntouckl,nlat,ntouckl,.true.))
            call tflocal(klist(itoul))
            klist(itoul)=ktflist+kax
c            if(tfcheckelement(ktflist+kax,.true.))then
c              write(*,*)'itoul: ',itoul
c            endif
          endif
          if(pri)then
            if(caltouck)then
              write(lfno,*)
              do iii=0,int((ntouckl-1)/5)
                write(lfno,9104)((5*iii+i+1)*0.2d0,i=1,5)
 9104           format(
     1               ' Momentum acceptance:  ',5(f8.1,2x),'  %')
                do i=1,5
                  if(5*iii+i .le. ntouckl)then
                    tt=touckl(5*iii+i)*tf
                  else
                    tt=0.d0
                  endif
                  if(tt .ne. 0.d0)then
                    vout(i)=autofg(1.d0/tt,'9.6')
                  else
                    vout(i)='   ---'
                  endif
                enddo
                write(lfno,9105)(vout(i)(1:9),i=1,5)
 9105           format(
     1               ' Touschek lifetime:    ',5(a,1x),' sec')
              enddo
              write(lfno,*)
              write(lfno,'(a)')
     1             ' Touschek lifetime/100s for aperture '//
     1             '2Jx/(Nx**2 emitx'') + 2Jz/(Nz**2 emitz) < 1:'
              write(lfno,9131)'Nx',(int(tampl(k,1)),
     1             (max(0,min(999,
     $             nint(.01d0/touckm(m,k,1)))),m=1,ntouckz),
     1             k=1,ntouckx)
 9131         format(
     1             '  ',
     1             'Nz: 4.6 5.3   6 6.8 7.7 8.8  10  11  13  15',
     $             '  17  19  21  24  28  32  36  41  46  53',
     $             '  60  68  77  88 100'/
     1             '  ',a,
     1   '!---1---1---1---1---1---1---1---1---1---1---1---1',
     1   '---1---1---1---1---1---1---1---1---1---1---1---1---1'/,
     1             34(i4,'!',25(i4)/))
              write(lfno,'(a)')
     1             ' Touschek lifetime/100s for aperture '//
     1             '2Jy/(Ny**2 emitx'') + 2Jz/(Nz**2 emitz) < 1:'
              write(lfno,9131)'Ny',(int(tampl(k,2)),
     1             (max(0,min(999,
     $             nint(0.01d0/touckm(m,k,2)))),m=1,ntouckz),
     1             k=1,ntouckx)
            endif
          endif
        endif
        iret=4
        call rsetgl1('EMITXC',eemx)
        call rsetgl1('EMITYC',eemy)
        call rsetgl1('EMITZC',eemz)
        return
      else
        caltouck=intra .and. de .lt. resib*10.d0 .and. .not. trpt
        pri=.false.
        if(.not.trpt)then
          if(calint)then
            if(intra)then
              rx=eintrb(emx0,eemx,emxr)/eemx
              ry=eintrb(emy0,eemy,emyr)/eemy
              rz=eintrb(emz0,eemz,emzr)/eemz
              rr=min(100.d0,max(0.01d0,(rx*ry*rz)**(1.d0/3.d0)))
              eemx=eemx*rr
              eemy=eemy*rr
              eemz=eemz*rr
            elseif(emx0 .ne. 0.d0)then
              if(it .ge. 20)then
                emxmax=min(max(eemx,emx0),emxmax)
                emxmin=max(min(eemx,emx0),emxmin)
                eemx=sqrt(emxmax*emxmin)
              else
                eemx=sqrt(eemx*emx0)
              endif
              if(it .gt. 30)then
                emymax=min(max(eemy,emy0),emymax)
                emymin=max(min(eemy,emy0),emymin)
                eemy=sqrt(emymax*emymin)
              else
                eemy=sqrt(eemy*emy0)
              endif
              emzmax=min(max(eemz,emz0),emzmax)
              emzmin=max(min(eemz,emz0),emzmin)
              eemz=sqrt(emzmax*emzmin)
            endif
          else
            emxr=eemx
            emyr=eemy
            emzr=eemz
          endif
          emx0=eemx
          emy0=eemy
          emz0=eemz
          calint=.true.
c     ccintr=(rclassic/h0**2)**2/8.d0/pi
c     cintrb=ccintr*pbunch/eemx/eemy/eemz
c
c     cintrb=rclassic**2/8.d0/pi
c     1           *pbunch/(eemx*h0)/(eemy*h0)/(eemz*h0)/h0
c     Here was the factor 2 difference from B-M paper.
c     Pointed out by K. Kubo on 6/18/2001.
c
          cintrb=rclassic**2/4.d0/pi*pbunch
          if(emx1 .gt. 0.01d0*eemx)then
            rx=sqrt(eemx/emx1)
          else
            emitn(iaidx(1,1))=eemx
            emitn(iaidx(2,2))=eemx
            rx=1.d0
          endif
          if(emy1 .gt. 0.01d0*eemy)then
            ry=sqrt(eemy/emy1)
          else
            emitn(iaidx(3,3))=eemy
            emitn(iaidx(4,4))=eemy
            ry=1.d0
          endif
          if(emz1 .gt. 0.01d0*eemz)then
            rz=sqrt(eemz/emz1)
          else
            emitn(iaidx(5,5))=eemz
            emitn(iaidx(6,6))=eemz
            rz=1.d0
          endif
          call tinitr(trans1)
          trans1(1,1)=rx
          trans1(2,2)=rx
          trans1(3,3)=ry
          trans1(4,4)=ry
          trans1(5,5)=rz
          trans1(6,6)=rz
          call tmulbs(emitn,trans1,.false.)
          if(.not. synchm)then
            emitn(iaidx(5,1))=0.d0
            emitn(iaidx(5,2))=0.d0
            emitn(iaidx(5,3))=0.d0
            emitn(iaidx(5,4))=0.d0
            emitn(iaidx(5,5))=sigz**2
            emitn(iaidx(5,6))=0.d0
            emitn(iaidx(6,6))=sige**2
            r(1:4,5)=0.d0
            r(5,5)=1.d0
            r(6,5)=0.d0
            r(6,1:5)=0.d0
            r(6,6)=1.d0
            r(5,1:4)=0.d0
            r(5,5)=1.d0
            r(5,6)=0.d0
          endif
          call tmulbs(emitn,r,.false.)
          beam(22:42)=emitn
c          write(*,*)'intraconv ',
c     $         sqrt(emitn(1)*emitn(3)-emitn(2)**2),
c     $         sqrt(emitn(6)*emitn(10)-emitn(9)**2),
c     $         sqrt(emitn(15)*emitn(21)-emitn(20)**2)
          it=it+1
        else
          beam(22:42)=0.d0
          beam(1:21)=beamin
          it=itmax+1
        endif
        iret=2
        return
      endif
      iret=1
      return
      end

      real*8 function eintrb(em0,y,emr)
      implicit none
      real*8 ,intent(in):: em0,y,emr
      real*8 ,parameter::eps=1.d-10
      integer*4 ,parameter ::itmax=30
      real*8 em,a,y1
      integer*4 it
      em=em0
      a=(y-emr)*em0**2
      y1=y
      do it=1,itmax
        em=em+(y1-em)/(em**3+2.d0*a)*em**3
        y1=emr+a/em**2
        if(abs(y1-em) .lt. eps*em)then
          eintrb=min(100.d0*em0,max(emr,y1,0.01d0*em0))
          return
        endif
      enddo
      eintrb=min(100.d0*em0,max(0.01d0*em0,em))
      return
      end

      subroutine tinv(r,ri,n,ndimr)
      implicit none
      integer*4 i,j,n,ndimr
      real*8 r(ndimr,n),ri(ndimr,n)
      do 10 i=1,n-1,2
        do 20 j=1,n-1,2
          ri(i  ,j  )= r(j+1,i+1)
          ri(i  ,j+1)=-r(j  ,i+1)
          ri(i+1,j  )=-r(j+1,i  )
          ri(i+1,j+1)= r(j  ,i  )
20      continue
10    continue
      return
      end

      subroutine setiamat(iamat,ri,codin,beamn,emitp,trans)
      use tfstk
      implicit none
      integer*8 , intent(in)::iamat
      real*8 , intent(in)::ri(6,6),codin(6),beamn(21),emitp(21),
     $     trans(6,12)
      if(iamat .gt. 0)then
        dlist(iamat+4)=
     $       dtfcopy1(kxm2l(ri,6,6,6,.false.))
        dlist(iamat+1)=
     $       dtfcopy1(kxm2l(codin,0,6,1,.false.))
        dlist(iamat+5)=
     $       dtfcopy1(kxm2l(beamn,0,21,1,.false.))
        dlist(iamat+6)=
     $       dtfcopy1(kxm2l(emitp,0,21,1,.false.))
        dlist(iamat+2)=
     $       dtfcopy1(kxm2l(trans,6,6,6,.false.))
        dlist(iamat+3)=
     $       dtfcopy1(kxm2l(trans(1,7),6,6,6,.false.))
      endif
      return
      end

      subroutine setparams(params,cod)
      use temw
      use tmacro
      implicit none
      real*8 , intent(out)::params(nparams)
      real*8 , intent(in)::cod(6)
      params(iprevf)=omega0/m_2pi
      params(ipdx:ipddp)=cod
      params(ipu0)=u0*pgev
      params(ipvceff)=vceff
      params(iptrf0)=trf0
      params(ipalphap)=alphap
      params(ipdleng)=dleng
      params(ipbh)=bh
      params(ipheff)=heff
      return
      end

      subroutine temitprint(lfno,intpri,
     $     rx,rxd,cod,rd,beamn,beamp,emitn,emitp,srot,srot1,spm,sps,
     $     ceig,cd,params,calem)
      use tfstk
      use temw
      use ffs_flag
      use ffs_pointer
      use tmacro
      use tffitcode
      use macmath
      implicit none
      integer*4 ,intent(in):: lfno
      logical*4 ,intent(in):: intpri,calem
      real*8 ,intent(in):: cod(6),rx(6,6),rxd(6,6),params(nparams),
     $     srot(3,9),srot1(3,3),spm(3,3),rd(6,6),sps(3,3),
     $     beamn(21),beamp(21),emitn(21),emitp(21)
      complex*16 ,intent(in):: ceig(6),cd(6)
      integer*4 i,j
      real*8 rxrxi(6,6),so,btilt,xxs,yys,sig1,sig2,sigx,sigy
      character*11 autofg,vout(nparams)
      character*13 at
      character*9 vout9(nparams)
      real*8 ,parameter :: somin=1.d-8
      if(intpri)then
        write(lfno,*)
        write(lfno,*)'   Closed orbit:'
        write(lfno,'(10X,6A)')label2
        call tput(codin,label2,' Entrance ','9.6',1,lfno)
        call tput(cod,label2,'     Exit ','9.6',1,lfno)
        write(lfno,*)
      endif
      if(intpri .or. wspac)then
        if(.not. intpri)then
          if(intra)then
            write(lfno,*)'   Parameters with intrabeam scattering ',
     $           '+ space charge:'
          else
            write(lfno,*)'   Parameters with space charge:'
          endif
          write(lfno,*)
        endif
        rxrxi=matmul(rx,tinv6(rx))
        if(emiout)then
          write(lfno,*)'   Symplectic part of the transfer matrix:'
          call tput(rx,label2,label2,'9.6',6,lfno)
          call tput(matmul(r,matmul(rx,ri)),
     $         label1,label1,'9.6',6,lfno)
          call tput(rxrxi,label2,label2,'9.6',6,lfno)
        endif
        do i=1,6
          rxrxi(i,i)=rxrxi(i,i)-1.d0
        enddo
        so=sum(abs(rxrxi))
        if(so .gt. somin)then
          write(lfno,*)' *** Deviation from symplectic matrix = ',so
        endif
        do i=1,ntwissfun
          vout9(i)=autofg(params(iptws0+i),'9.6')
        enddo
        if(trpt)then
          at=' at entrance:'
        else
          at=':'
        endif
        write(lfno,9001)at(1:len_trim(at)),
     $         vout9(mfitax),vout9(mfitbx),vout9(mfitzx),vout9(mfitex),
     $         vout9(mfitnx),vout9(mfitzpx),vout9(mfitey),
     $         vout9(mfitr1),vout9(mfitr2),vout9(mfitay),vout9(mfitby),
     $         vout9(mfitzy),vout9(mfitey),
     $         vout9(mfitr3),vout9(mfitr4),vout9(mfitny),
     $         vout9(mfitzpy),vout9(mfitepy),
     $         vout9(mfitaz),vout9(mfitbz),vout9(mfitnz)
 9001     format('    Extended Twiss Parameters',a,/,
     $         'AX:',a,' BX:',a,              26x,'  ZX:',a,'  EX:',a,/
     $         11x,'PSIX:',a,              26x,' ZPX:',a,' EPX:',a,/
     $         'R1:',a,' R2:',a,' AY:',a,' BY:',a,'  ZY:',a,'  EY:',a,/
     $         'R3:',a,' R4:',a,    12x,'PSIY:',a,' ZPY:',a,' EPY:',a,/
     $         51x,'  AZ:',a,'  BZ:',a,/
     $         65x,'PSIZ:',a,/
     $         '    Units: B(X,Y,Z), E(X,Y), R2: m ',
     $         '| PSI(X,Y,Z): radian | ZP(X,Y), R3: 1/m',/)
          vout(1) =autofg(pgev/1.d9       ,'10.7')
          vout(2) =autofg(omega0/m_2pi      ,'10.7')
          vout(3) =autofg(u0*pgev/1.d6    ,'10.7')
          vout(4) =autofg(vceff/1.d6      ,'10.7')
          vout(5) =autofg(trf0*1.d3       ,'10.7')
          vout(6) =autofg(alphap          ,'10.7')
          vout(7) =autofg(-dleng*1.d3     ,'10.7')
          vout(8) =autofg(heff            ,'10.7')
          vout(9) =autofg(bh              ,'10.7')
          vout(10)=autofg(omegaz/m_2pi      ,'10.7')
          write(lfno,9101)vout(1:10)(1:10)
 9101     format(   'Design momentum      P0 =',a,' GeV',
     1         1x,'Revolution freq.     f0 =',a,' Hz '/
     1         'Energy loss per turn U0 =',a,' MV ',
     1         1x,'Effective voltage    Vc =',a,' MV '/
     1         'Equilibrium position dz =',a,' mm ',
     1         1x,'Momentum compact. alpha =',a/
     1         'Orbit dilation       dl =',a,' mm ',
     1         1x,'Effective harmonic #  h =',a,/
     1         'Bucket height     dV/P0 =',a,'    ',
     $         1x,'Synchrotron frequency   =',a,' Hz '/)
          if(emiout)then
            if(trpt)then
              write(lfno,*)
     $    '   Phase rotations and normal coordinates at entrance:'
            else
              write(lfno,*)'   Eigen values and eigen vectors:'
            endif
            write(lfno,*)
            write(lfno,9011)'     Real:',(dble(ceig(j)),j=1,6)
            write(lfno,9011)'Imaginary:',(imag(ceig(j)),j=1,6)
          endif
          write(lfno,9012)'Imag.tune:',(dble(cd(j))/m_2pi,j=4,6)
          write(lfno,9012)'Real tune:',(imag(cd(j))/m_2pi,j=4,6)
          write(lfno,*)
 9011     format(2x,a,6f10.7)
 9012     format(2x,a,3(f10.7,10x))
          if(emiout)then
            call tput(r ,label1,label2,'9.6',6,lfno)
            call tput(ri,label2,label1,'9.6',6,lfno)
            if(.not. trpt)then
              call tput(matmul(r,matmul(rx,ri)),
     $             label1,label1,'9.6',6,lfno)
            endif
            write(lfno,*)'   Radiation part of the transfer matrix:'
            call tput(rxd,label2,label2,'9.6',6,lfno)
            call tput(rd ,label1,label1,'9.6',6,lfno)
          endif
          if(trpt)then
            write(lfno,*)'   Damping through the beam line:'
          else
            write(lfno,*)'   Damping per one revolution:'
          endif
          write(lfno,9013)
     1         'X :',dble(cd(1)),'Y :',dble(cd(2)),'Z :',dble(cd(3))
 9013     format(10x,3(a,1p,g14.6,3x))
          if(.not. trpt)then
            write(lfno,*)'   Damping time (sec):'
            write(lfno,9013)
     1           'X :',-m_2pi/omega0/dble(cd(1)),
     $           'Y :',-m_2pi/omega0/dble(cd(2)),
     $           'Z :',-m_2pi/omega0/dble(cd(3))
            write(lfno,*)'   Tune shift due to radiation:'
            write(lfno,9013)
     1           'X :',imag(cd(1))/m_2pi,
     $           'Y :',imag(cd(2))/m_2pi,'Z :',imag(cd(3))/m_2pi
          endif
          write(lfno,*)'   Damping partition number:'
          write(lfno,9014)
     1         'X :',params(ipjx),'Y :',params(ipjy),'Z :',params(ipjz)
c     1         'X :',dble(cd(1))*sr,'Y :',dble(cd(2))*sr,
c     $         'Z :',dble(cd(3))*sr
 9014     format(10x,3(a,f10.4,7x))
          write(lfno,*)
        endif
        if(.not. calem .and. .not. trpt)then
          return
        endif
        if(.not. intpri .and. intra .and. .not. wspac)then
          write(lfno,*)'   Parameters with intrabeam scattering:'
          write(lfno,*)
        endif
        if(emiout)then
          if(calint)then
            if(intra)then
              if(wspac)then
                write(lfno,*)
     1      '   Beam matrix by radiation+intrabeam fluctuation ',
     $      'with space charge:'
              else
                write(lfno,*)
     1      '   Beam matrix by radiation+intrabeam fluctuation:'
              endif
            elseif(wspac)then
              write(lfno,*)
     1             '   Beam matrix with space charge:'
            endif
          else
            write(lfno,*)'   Beam matrix by radiation fluctuation:'
          endif
          call tputbs(beamp,label2,lfno)
          call tputbs(beamn,label1,lfno)
          if(trpt)then
            write(lfno,*)'   Beam matrix at exit:'
          else
            write(lfno,*)'   Equiliblium beam matrix:'
          endif
          call tputbs(emitp,label2,lfno)
          call tputbs(emitn,label1,lfno)
        endif
        if(trpt)then
          write(lfno,*)'   Parameters at exit:'
        endif
        if(.not. intpri .and. (intra .or. wspac))then
          vout(1)=autofg(pbunch/1.d10 ,'11.8')
          vout(2)=autofg(coumin*1.d2,'11.8')
          write(lfno,9103)(vout(i)(1:11),i=1,2)
 9103     format( 'Particles/bunch/1e10   =',a,'    ',
     $         1x,'Minimum coupling       =',a,' %  ')
        endif
        vout(1)=autofg(eemx             ,'11.8')
        vout(2)=autofg(eemy             ,'11.8')
        vout(3)=autofg(eemz             ,'11.8')
        vout(4)=autofg(sige            ,'11.8')
        vout(5)=autofg(sigz*1.d3       ,'11.8')
        vout(9)=autofg(params(ipnnup)   ,'11.7')
        vout(10)=autofg(params(iptaup)/60.d0,'11.7')
c     kiku <------------------
        xxs=emitp(1)-emitp(6)
        yys=-2.d0*emitp(4)
        if(xxs .ne. 0.d0 .and. yys .ne. 0.d0)then
          btilt= atan2(yys,xxs) /2d0
        else
          btilt=0.d0
        endif
        sig1 = abs(emitp(1)+emitp(6))/2d0
c     sig2 = 0.5d0* sqrt(abs((emitp(1)-emitp(6))**2+4d0*emitp(4)**2))
        sig2=0.5d0*hypot(emitp(1)-emitp(6),2.d0*emitp(4))
        sigx = max(sqrt(sig1+sig2),sqrt(abs(sig1-sig2)))
        sigy = min(sqrt(sig1+sig2),sqrt(abs(sig1-sig2)))
        vout(6)=autofg(btilt,'11.8')
        vout(7)=autofg(sigx*1.d3  ,'11.8')
        vout(8)=autofg(sigy*1.d3  ,'11.8')
        write(lfno,9102)vout(1:10)
 9102   format( 'Emittance X            =',a,' m  ',
     1       1x,'Emittance Y            =',a,' m'/
     1       'Emittance Z            =',a,' m  ',
     1       1x,'Energy spread          =',a,/
     1       'Bunch Length           =',a,' mm ',
     1       1X,'Beam tilt              =',a,' rad'/
     1       'Beam size xi           =',a,' mm ',
     1       1X,'Beam size eta          =',a,' mm'/
     $       ,'Nominal spin tune      =',a,'    ',
     $       1x,'Polarization time      =',a,' min'/)
c     9103   format(3X,'Beam dimension along principal axis:'/
ccc   call putsti(eemx,eemy,eemz,sige,sigz,btilt,sigx,sigy,
ccc   1              calint,fndcod)
c     kiku ------------------>
        if(calpol)then
          write(lfno,*)'  Polarization vector at entrance:'
          write(lfno,9013)
     1         'x :',sps(1,1),'y :',sps(2,1),'z :',sps(3,1)
          if(emiout)then
            write(lfno,*)'\n'//'   Spin precession matrix:'
            write(lfno,*)'                sx             sy'//
     $           '             sz'
            write(*,'(a,1p3g15.7)')'        sx',srot1(1,:)
            write(*,'(a,1p3g15.7)')'        sy',srot1(2,:)
            write(*,'(a,1p3g15.7)')'        sz',srot1(3,:)
            write(lfno,*)'\n'//'   One-turn depolarization vectors:'
            write(lfno,*)'                x              px'//
     $           '             y              py'//
     $           '             z              pz'
            write(*,'(a,1p6g15.7)')'        sx',srot(1,4:9)
            write(*,'(a,1p6g15.7)')'        sy',srot(2,4:9)
            write(*,'(a,1p6g15.7)')'        sz',srot(3,4:9)
            write(lfno,*)'\n'//'   Spin depolarization matrix:'
            write(lfno,*)'                s1             s2'//
     $           '             s3'
            write(*,'(a,1p6g15.7)')'        s1',spm(1,:)
            write(*,'(a,1p6g15.7)')'        s2',spm(2,:)
            write(*,'(a,1p6g15.7)')'        s3',spm(3,:)
          endif
          vout(1)=autofg(spinmu/m_2pi,'11.8')
          vout(2)=autofg(equpol(3)*100.d0,'11.8')
          write(lfno,9110)vout(1:2)
 9110     format(/'Spin tune              =',a,'    ',
     1         1x,'Equil. polarization    =',a,' %'/)
        endif
        return
        end
