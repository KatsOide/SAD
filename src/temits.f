      subroutine temits(
     $     ndim,ntwissfun,
     $     mphi2,amus0,amus1,amusstep,
     $     emx,emy,res,params,
     $     lfni,lfno,kx,irtc)
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tmacro
      implicit none
      integer*8 kx
      integer*4 ntwissfun,irtc,ndim,mphi2,lfni,lfno,mphiz,ndims,
     $     ndps,nzz
      real*8 amus0,amus1,amusstep,emx,emy,res,size
      real*8 params(59)
      logical*4 trpt0,stab
      trpt0=trpt
      trpt=.false.
      mphiz=mphi2*2-1
      ndims=mphiz
      ndps=mphiz
      nzz=mphiz*2+1
      call temits1(twiss,size,gammab,
     $     ndim,ntwissfun,codplt,stab,ndps,nzz,mphiz,mphi2,ndims,
     $     amus0,amus1,amusstep,
     $     emx,emy,res,params,
     $     lfni,lfno,kx,irtc)
      trpt=trpt0
      return
      end

      subroutine temits1(twiss,size,gammab,
     $     ndim,ntwissfun,plot,stab,ndp,nz,mphi,mphi2,ndims,
     $     amus0,amus1,amusstep,
     $     emx,emy,res0,params,
     $     lfni,lfno,kx,irtc)
      use tfstk
      use ffs_flag
      use temw, only: r, ri
      use tmacro
      implicit none
      integer*4 , parameter:: npara=59,itmax=1000
      integer*8 kx,kax,kai
      integer*4 ndim,ntwissfun,mphi,mphi2,ndims,
     $     lfni,lfno,irtc,ndp,i,k,kk,ns,it,j,nz
      real*8 amus0,amus1,amusstep, emx,emy,res0,
     $     amus,amuss,damp,dj,dpf,dpndim,emx0,emy0,eps,
     $     fz,phi0s,pa,res,sige,sigea,vx,vy,w,dpi
      real*8 beam(42),trans(6,12),cod(6),
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
     $     dbc(10,mphi2,ndp),dbs(10,mphi2,ndp),
     $     dhc(4,mphi2,ndp),dhs(4,mphi2,ndp),
     $     amuj(ndp)
      real*8 params(npara),twiss(nlat,-ndim:ndim,ntwissfun),
     $     size(21,nlat),gammab(nlat),btr(21,21)
      real*8 rm(6,12),rx(6,6),rxi(6,6),cmu(256),smu(256),
     $     beamr(42),fj(256),conv
      character*10 label1(6),label2(6)
      logical*4 rfsw0,plot,stab
      data label1/'        X ','       Px ','        Y ',
     1     '       Py ','        Z ','       Pz '/
      data label2/'        x ','    px/p0 ','        y ',
     1     '    py/p0 ','        z ','    dp/p0 '/
      rfsw0=rfsw
      rfsw=.false.
      call tclr(codin,6)
      call tclr(beamin,21)
      call temit(trans,cod,beamr,btr,
     $     .true.,0,0,0,0,
     $     plot,params,stab,lfni,lfno)
      call tinitr(rx)
      call tmov(rx,rxi,36)
      do i=1,4
        rx(i,6)=r(i,6)
        rxi(i,6)=-r(i,6)
      enddo
      sige=params(25)
      damp=params(18)
c      conv=1.d-4*abs(damp)
      conv=1.d-7
c      write(*,*)'emits ',damp,sige
      sigea=sqrt(12.d0)*sige
      dpndim=sigea/ndims
      calint=intra
      irad=12
      do i=-ndims,ndims
        dpi=i*dpndim
        do j=1,4
          cod(j)=codin(j)+r(j,6)*dpi
        enddo
        cod(5)=0.d0
        cod(6)=codin(6)+dpi
        call tinitr(trans)
        call tclr(trans(1,7),36)
        call tclr(beam,21)
        call tturne(trans,cod,beam,int8(0),int8(0),int8(0),
     $       .false.,.false.,.false.)
        call tmultr(trans,rxi,12)
        call tmov(rx,rm,36)
        call tmultr(rm,trans,6)
        call tmov(rx,trans,36)
        call tmultr(trans,trans(1,7),6)
        dpf=cod(6)-codin(6)
        do j=1,4
          cod(j)=cod(j)-codin(j)-dpf*r(j,6)
        enddo
        call tediag(rm,cod,tws,i,ndims,ntwissfun,lfno)
        call tmulbs(beam,rxi,.false.,.false.)
        call tmov(beam,beams(1,i),10)
        call tmov65(trans,trads(1,1,i))
      enddo
      if(amusstep .eq. 0.d0)then
        ns=1
      else
        ns=int((amus1-amus0)/amusstep+1.1d0)
      endif
      if(ns .eq. 1)then
        amuss=0.d0
      else
        amuss=(amus1-amus0)/(ns-1)
      endif
      dj=sigea**2/ndp*.5d0
      w=0.d0
      do j=1,ndp
        fj(j)=exp(-(j-.5d0)*dj/sige**2)
        w=w+fj(j)*dj/sige**2
      enddo
      call tefsetup(ha,hb,ba,bb,bd,ndp,mphi,mphi2,nz,ndims,
     $     ntwissfun,
     $     beams,trads,tws,dj,btr,btr(1,6),dpndim,sigea)
      eps=sige**2
      phi0s=params(27)*omega0*params(12)/c
      fz=cos(phi0s)/(1.d0+cos(phi0s))*.25d0
     $     *(params(27)*params(13)*pi2)**2
      if(kx .ne. 0)then
        kax=ktadaloc(-1,ns)
      else
        kax=0
      endif
      do kk=1,ns
        amus=amus0+amusstep*(kk-1)
        do j=1,ndp
          amuj(j)=amus*(1.d0-fz*(j-.5d0)*dj/amus**2)
        enddo
        call tesolvd(hff,hfb,hfx,
     $       hb,hc,hs,ha,amuj,cmu,smu,
     $       bd,hc,hs,
     $       ndp,mphi2,mphi,4,
     $       dj,damp,sige)
        do k=1,16
          call tevdif(hfb,hb,hc,hs,ha,amuj,cmu,smu,
     $         dhc,dhs,
     $         bd,hc,hs,fj,
     $         ndp,mphi2,mphi,4,
     $         dj,damp,sige)
        enddo
        call tesolvd(bff,bfb,bfx,
     $       bb,bc,bs,ba,amuj,cmu,smu,
     $       bd,hc,hs,
     $       ndp,mphi2,mphi,10,
     $       dj,damp,sige)
c          write(*,'(1p 8g12.4)')((hfx(i+(k-1)* 4),i=1, 4),k=1,2*mphi2)
        emx0=0.d0
        emy0=0.d0
        call tesumb(bc,hc,mphi,mphi2,ndp,
     $       w,dj,sige,tws,ndims,ntwissfun,beam,fj)
        call tmov(beam,beamr,21)
        call tmulbs(beamr,ri,.false.,.false.)
        vx=beamr(1)*beamr(3)-beamr(2)**2
        vy=beamr(6)*beamr(10)-beamr(9)**2
        emx0=sign(sqrt(abs(vx)),vx)
        emy0=sign(sqrt(abs(vy)),vy)
c        write(*,*)emx0,emy0
        res0=1.d100
        res=1.d0
        it=0
        do while(res .gt. conv)
c          write(*,'(1p4g15.7)')(hc(i,3,10),i=1,4)
          call tevdif(hfb,hb,hc,hs,ha,amuj,cmu,smu,
     $         dhc,dhs,
     $         bd,hc,hs,fj,
     $         ndp,mphi2,mphi,4,
     $         dj,damp,sige)
          call tevdif(bfb,bb,bc,bs,ba,amuj,cmu,smu,
     $         dbc,dbs,
     $         bd,hc,hs,fj,
     $         ndp,mphi2,mphi,10,
     $         dj,damp,sige)
          call tesumb(bc,hc,mphi,mphi2,ndp,
     $         w,dj,sige,tws,ndims,ntwissfun,beam,fj)
          call tmov(beam,beamr,21)
          call tmulbs(beamr,ri,.false.,.false.)
          vx=beamr(1)*beamr(3)-beamr(2)**2
          vy=beamr(6)*beamr(10)-beamr(9)**2
          emx=sign(sqrt(abs(vx)),vx)
          emy=sign(sqrt(abs(vy)),vy)
          res=((emx-emx0)**2+(emy-emy0)**2)/(emx+emy)**2
          emx0=emx
          emy0=emy
c         write(*,*)'synchrob ',it,emx,emy,res
          it=it+1
          if((it .gt. itmax .and. res .gt. conv)
     $         .or. it .gt. 1. and.
     $         (res .gt. 1.d-4 .or. res .gt. res0))then
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
            write(lfno,'(a,f9.6,2(a,1pg13.6),a,1pg11.3)')
     $           ' Nuz =',amus/pi2,
     $           '  Emittances X:',emx,
     $           '  Y:',emy,
     $           '  Conv.:',res0
            if(emiout)then
              write(lfno,*)
              write(lfno,*)'   Equiliblium beam matrix:'
              call tputbs(beam,label2,lfno)
              call tputbs(beamr,label1,lfno)
              call tedrawf(8,bc,bs,hc,hs,pa,
     $             sige,amus/pi2,mphi,mphi2,ndp)
            endif
          endif
        else
          kai=ktavaloc(0,4)
          rlist(kai+1)=amus/pi2
          rlist(kai+2)=emx
          rlist(kai+3)=emy
          rlist(kai+4)=res0
          klist(kax+kk)=ktflist+kai
        endif
      enddo
      if(kx .ne. 0)then
        kx=ktflist+kax
      endif
      rfsw=rfsw0
      irtc=0
      return
      end
      
      subroutine tevdif(bfb,bb,bc,bs,ba,amuj,cmu,smu,
     $     dbc,dbs,
     $     bd,hc,hs,fj,
     $     ndp,mphi2,mphi,nd,
     $     dj,damp0,sige)
      implicit none
      integer*4 ndp,mphi2,mphi,nd
      real*8 bfb(nd*2*mphi2),bb(nd,mphi2,ndp),
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     dbc(nd,mphi2,ndp),dbs(nd,mphi2,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     ba(nd,nd,mphi,ndp),
     $     bd(10,4,mphi,ndp),
     $     amuj(ndp),
     $     cmu(32),smu(32),
     $     dj,damp,sige,eps
      integer*4 j,i,m,mb,ms,nsb,k1,k2,m1,m2,
     $     i1,nsbp
      real*8 diff,diff1,diff2,aj,bci1,bsi1,
     $     diff3,dc,ds,damp0,d,fj(50),
     $     fbc(10),dfbc(10),fbs(10),dfbs(10)
      nsb=nd*mphi2
      nsbp=nsb*ndp
      eps=sige**2
      damp=damp0*.5d0
      d=2.d0*damp
c      call tmov(bc,dbc,nsbp)
c      call tmov(bs,dbs,nsbp)
      do j=1,ndp
        do m=1,mphi2
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          if(m .eq. 2)then
            do i=1,nd
              bc(i,m,j)=bc(i,m,j)+2.d0*dfbc(i)
            enddo
          else
            do i=1,nd
              bc(i,m,j)=bc(i,m,j)+dfbc(i)
              bs(i,m,j)=bs(i,m,j)+dfbs(i)
            enddo
          endif
          m1=m+2
          if(m1 .le. mphi2)then
            do i=1,nd
              bc(i,m1,j)
     $             =bc(i,m1,j)+dfbc(i)-m*fbc(i)
              bs(i,m1,j)
     $             =bs(i,m1,j)+dfbs(i)-m*fbs(i)
            enddo
          endif
        enddo
      enddo
      do j=ndp,1,-1
        do m=mphi2,1,-1
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          m2=m-2
          if(m .eq. 3)then
            do i=1,nd
              bc(i,m2,j)
     $             =bc(i,m2,j)+dfbc(i)+2.d0*fbc(i)
            enddo
          elseif(m .gt. 3)then
            do i=1,nd
              bc(i,m2,j)
     $             =bc(i,m2,j)+dfbc(i)+(m-2)*fbc(i)
              bs(i,m2,j)
     $             =bs(i,m2,j)+dfbs(i)+(m-2)*fbs(i)
            enddo
          endif
          do i=1,nd
            bc(i,m,j)=bc(i,m,j)+dfbc(i)
            bs(i,m,j)=bs(i,m,j)+dfbs(i)
          enddo
        enddo
      enddo
      do j=1,ndp
        aj=(j-.5d0)*dj
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        do i=1,nd
          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          do i=1,nd
            bc(i,m,j)=diff*bc(i,m,j)
            bs(i,m,j)=diff*bs(i,m,j)
          enddo
        enddo
        do m=1,mphi2-2
          m2=m+2
          diff2=-d*(m+1)*
     $         (-1.d0+eps/aj*.5d0*(m-1))*.25d0
          do i=1,nd
            bc(i,m2,j)=bc(i,m2,j)+diff2*bc(i,m,j)
            bs(i,m2,j)=bs(i,m2,j)+diff2*bs(i,m,j)
          enddo
        enddo
        do m=mphi2,4,-1
          m1=m-2
          diff1=-d*(m-3)*
     $         (1.d0+eps/aj*.5d0*(m-1))*.25d0
          do i=1,nd
            bc(i,m1,j)=bc(i,m1,j)+diff1*bc(i,m,j)
            bs(i,m1,j)=bs(i,m1,j)+diff1*bs(i,m,j)
          enddo
        enddo
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        do i=1,nd
          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          do i=1,nd
            bc(i,m,j)=diff*bc(i,m,j)
            bs(i,m,j)=diff*bs(i,m,j)
          enddo
        enddo
      enddo
      do j=1,ndp
        call tclr(bfb,2*nsb)
        if(nd .eq. 10)then
          call tesetdm(j,bfb,bd,hc,hs,ndp,mphi,mphi2)
        endif
        call tadd(bb(1,1,j),bfb,bfb,nsb)
        call tetrig(amuj,cmu,smu,mphi2,ndp,j,1.d0)
        do m=1,mphi2
          mb=(m-1)*nd
          ms=mb+nsb
          do i1=1,nd
            bci1=bc(i1,1,j)
            bsi1=bs(i1,1,j)
            do i=1,nd
              bfb(mb+i)=bfb(mb+i)+ba(i,i1,m,j)*bci1
              bfb(ms+i)=bfb(ms+i)+ba(i,i1,m,j)*bsi1
            enddo
          enddo
          do m1=2,mphi2
            k1=abs(m-m1)+1
            k2=abs(m+m1-2)+1
            do i1=1,nd
              bci1=bc(i1,m1,j)
              bsi1=bs(i1,m1,j)
              do i=1,nd
                bfb(mb+i)=bfb(mb+i)
     $               +(ba(i,i1,k1,j)+ba(i,i1,k2,j))*bci1
                bfb(ms+i)=bfb(ms+i)
     $               +(ba(i,i1,k1,j)-ba(i,i1,k2,j))*bsi1
              enddo
            enddo
          enddo
        enddo
        do m=1,mphi2
          mb=(m-1)*nd
          ms=mb+nsb
          do i=1,nd
            dc=bfb(mb+i)
            ds=bfb(ms+i)
            bc(i,m,j)=
     $           cmu(m)*dc+smu(m)*ds
            bs(i,m,j)=
     $           -smu(m)*dc+cmu(m)*ds
          enddo
        enddo
      enddo
      do j=1,ndp
        aj=(j-.5d0)*dj
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        do i=1,nd
          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          do i=1,nd
            bc(i,m,j)=diff*bc(i,m,j)
            bs(i,m,j)=diff*bs(i,m,j)
          enddo
        enddo
        do m=mphi2,4,-1
          m1=m-2
          diff1=-d*(m-3)*
     $         (1.d0+eps/aj*.5d0*(m-1))*.25d0
          do i=1,nd
            bc(i,m1,j)=bc(i,m1,j)+diff1*bc(i,m,j)
            bs(i,m1,j)=bs(i,m1,j)+diff1*bs(i,m,j)
          enddo
        enddo
        do m=1,mphi2-2
          m2=m+2
          diff2=-d*(m+1)*
     $         (-1.d0+eps/aj*.5d0*(m-1))*.25d0
          do i=1,nd
            bc(i,m2,j)=bc(i,m2,j)+diff2*bc(i,m,j)
            bs(i,m2,j)=bs(i,m2,j)+diff2*bs(i,m,j)
          enddo
        enddo
        m=2
        diff=d*eps*(m-1)**2*.125d0/aj
        diff3=-d*(m-3)*
     $       (1.d0+eps/aj*.5d0*(m-1))*.125d0
        do i=1,nd
          bc(i,m,j)=bc(i,m,j)+(diff+diff3)*bc(i,m,j)
          bs(i,m,j)=bs(i,m,j)+(diff-diff3)*bs(i,m,j)
        enddo
        do m=3,mphi2
          diff=1.d0+d*eps*(m-1)**2*.125d0/aj
          do i=1,nd
            bc(i,m,j)=diff*bc(i,m,j)
            bs(i,m,j)=diff*bs(i,m,j)
          enddo
        enddo
      enddo
      do j=ndp,1,-1
        do m=mphi2,1,-1
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          m2=m-2
          if(m .eq. 3)then
            do i=1,nd
              bc(i,m2,j)
     $             =bc(i,m2,j)+dfbc(i)+2.d0*fbc(i)
            enddo
          elseif(m .gt. 3)then
            do i=1,nd
              bc(i,m2,j)
     $             =bc(i,m2,j)+dfbc(i)+(m-2)*fbc(i)
              bs(i,m2,j)
     $             =bs(i,m2,j)+dfbs(i)+(m-2)*fbs(i)
            enddo
          endif
          do i=1,nd
            bc(i,m,j)=bc(i,m,j)+dfbc(i)
            bs(i,m,j)=bs(i,m,j)+dfbs(i)
          enddo
        enddo
      enddo
      do j=1,ndp
        do m=1,mphi2
          call tespl(bc,fbc,dfbc,ndp,mphi2,nd,dj,eps,damp,j,m)
          call tespl(bs,fbs,dfbs,ndp,mphi2,nd,dj,eps,damp,j,m)
          if(m .eq. 2)then
            do i=1,nd
              bc(i,m,j)=bc(i,m,j)+2.d0*dfbc(i)
            enddo
          else
            do i=1,nd
              bc(i,m,j)=bc(i,m,j)+dfbc(i)
              bs(i,m,j)=bs(i,m,j)+dfbs(i)
            enddo
          endif
          m1=m+2
          if(m1 .le. mphi2)then
            do i=1,nd
              bc(i,m1,j)
     $             =bc(i,m1,j)+dfbc(i)-m*fbc(i)
              bs(i,m1,j)
     $             =bs(i,m1,j)+dfbs(i)-m*fbs(i)
            enddo
          endif
        enddo
      enddo
      return
      end
      
      subroutine tesolvd(bff,bfb,bfx,
     $     bb,bc,bs,ba,amuj,cmu,smu,
     $     bd,hc,hs,
     $     ndp,mphi2,mphi,nd,
     $     dj,damp0,sige)
      implicit none
      integer*4 ndp,mphi2,mphi,nd
      real*8 bfb(nd*2*mphi2),bfx(nd*2*mphi2),bb(nd,mphi2,ndp),
     $     bff(2*nd*mphi2,2*nd*mphi2),
     $     bc(nd,mphi2,ndp),bs(nd,mphi2,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     ba(nd,nd,mphi,ndp),
     $     bd(10,4,mphi,ndp),
     $     amuj(ndp),
     $     cmu(32),smu(32),
     $     dj,damp,sige,eps
      integer*4 j,i,m,mb,ms,nsb,k1,k2,m1,m2,
     $     ia,ma,maci,masi,mc
      real*8 diff,diff1,diff2,aj,bci1,bsi1,
     $     diff3,dc,ds,damp0,d
      nsb=nd*mphi2
      eps=sige**2
      damp=damp0*.5d0
      d=2.d0*damp
      do j=1,ndp
        aj=(j-.5d0)*dj
        do ma=1,mphi2
          do ia=1,nd
            call tclr(bc(1,1,j),nsb)
            call tclr(bs(1,1,j),nsb)
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
            m2=m+2
            if(m2 .le. mphi2)then
              diff2=-d*(m+1)*
     $             (-1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m2,j)=bc(ia,m2,j)+diff2*bc(ia,m,j)
              bs(ia,m2,j)=bs(ia,m2,j)+diff2*bs(ia,m,j)
              m=m2
              m1=m-2
              diff1=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m1,j)=bc(ia,m1,j)+diff1*bc(ia,m,j)
              bs(ia,m1,j)=bs(ia,m1,j)+diff1*bs(ia,m,j)
            endif
            m=ma
            if(m .gt. 3)then
              m1=m-2
              diff1=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m1,j)=bc(ia,m1,j)+diff1*bc(ia,m,j)
              bs(ia,m1,j)=bs(ia,m1,j)+diff1*bs(ia,m,j)
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
            call tetrig(amuj,cmu,smu,mphi2,ndp,j,1.d0)
            call tclr(bfb,2*nsb)
            do m=1,mphi2
              mb=(m-1)*nd
              ms=mb+nsb
              do m1=ma-2,ma+2,2
                if(m1 .eq. 1)then
                  bci1=bc(ia,1,j)
                  bsi1=bs(ia,1,j)
                  do i=1,nd
                    bfb(mb+i)=bfb(mb+i)+ba(i,ia,m,j)*bci1
                    bfb(ms+i)=bfb(ms+i)+ba(i,ia,m,j)*bsi1
                  enddo
                elseif(m1 .gt. 0 .and. m1 .le. mphi2)then
                  k1=abs(m-m1)+1
                  k2=abs(m+m1-2)+1
                  bci1=bc(ia,m1,j)
                  bsi1=bs(ia,m1,j)
                  do i=1,nd
                    bfb(mb+i)=bfb(mb+i)
     $                   +(ba(i,ia,k1,j)+ba(i,ia,k2,j))*bci1
                    bfb(ms+i)=bfb(ms+i)
     $                   +(ba(i,ia,k1,j)-ba(i,ia,k2,j))*bsi1
                  enddo
                endif
              enddo
            enddo
            maci=(ma-1)*nd+ia
            masi=maci+nsb
            do m=1,mphi2
              mc=(m-1)*nd
              ms=mc+nsb
              do i=1,nd
                dc=bfb(mc+i)
                ds=bfb(ms+i)
                bff(mc+i,maci)=-cmu(m)*dc
                bff(ms+i,maci)= smu(m)*dc
                bff(mc+i,masi)=-smu(m)*ds
                bff(ms+i,masi)=-cmu(m)*ds
              enddo
            enddo
            call tclr(bc(1,1,j),nsb)
            call tclr(bs(1,1,j),nsb)
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
              m=m2
              m1=m-2
              diff1=-d*(m-3)*
     $             (1.d0+eps/aj*.5d0*(m-1))*.25d0
              bc(ia,m1,j)=bc(ia,m1,j)-diff1*bc(ia,m,j)
              bs(ia,m1,j)=bs(ia,m1,j)-diff1*bs(ia,m,j)
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
          call tclr(bfb,2*nsb)
        endif
        call tadd(bb(1,1,j),bfb,bfb,nsb)
        do m=1,mphi2
          mc=(m-1)*nd
          ms=mc+nsb
          do i=1,nd
            dc=bfb(mc+i)
            ds=bfb(ms+i)
            bfb(mc+i)= cmu(m)*dc+smu(m)*ds
            bfb(ms+i)=-smu(m)*dc+cmu(m)*ds
          enddo
        enddo
        call tsolva(bff,bfb,bfx,2*nsb,2*nsb,2*nsb,1.d-20)
        call tmov(bfx,bc(1,1,j),nsb)
        call tmov(bfx(nsb+1),bs(1,1,j),nsb)
      enddo
      return
      end

      subroutine tesetdm(j,bfb,bd,hc,hs,ndp,mphi,mphi2)
      implicit none
      integer*4 ndp,mphi,mphi2,j,m,m1,i1,i,k1,k2,
     $     mb,ms,nsb
      real*8
     $     bfb(20*mphi2),bd(10,4,mphi,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     hci1,hsi1
      nsb=10*mphi2
      call tclr(bfb,2*nsb)
      do m=1,mphi2
        mb=(m-1)*10
        ms=mb+nsb
        do i1=1,4
          hci1=hc(i1,1,j)
          hsi1=hs(i1,1,j)
          do i=1,10
            bfb(mb+i)=bfb(mb+i)
     $           +(bd(i,i1,m,j))*hci1
            bfb(ms+i)=bfb(ms+i)
     $           +(bd(i,i1,m,j))*hsi1
          enddo
        enddo
        do m1=2,mphi2
          k1=abs(m-m1)+1
          k2=abs(m+m1-2)+1
          do i1=1,4
            hci1=hc(i1,m1,j)
            hsi1=hs(i1,m1,j)
            do i=1,10
              bfb(mb+i)=bfb(mb+i)
     $             +(bd(i,i1,k1,j)+bd(i,i1,k2,j))
     $             *hci1
              bfb(ms+i)=bfb(ms+i)
     $             +(bd(i,i1,k1,j)-bd(i,i1,k2,j))
     $             *hsi1
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine tetrig(amuj,cmu,smu,mphi2,ndp,j,fact)
      implicit none
      integer*4 mphi2,m,ndp,j
      real*8 amuj(ndp),cmu(mphi2),smu(mphi2),
     $     phim,fact
      do m=1,mphi2
        phim=(m-1)*amuj(j)
        cmu(m)=cos(phim)
        smu(m)=sin(phim)
      enddo
      return
      end

      subroutine tefsetup(ha,hb,ba,bb,bd,
     $     ndp,mphi,mphi2,nz,ndims,ntwissfun,
     $     beams,trads,tws,dj,btr,btrd,dpndim,sigea)
      use macmath
      implicit none
      integer*4 mphi,mphi2,nz,ndp,i,j,k,l,m,ndims,ntwissfun
      real*8 
     $     beams(10,-ndims:ndims),
     $     trads(5,5,-ndims:ndims),
     $     tws(ntwissfun,-ndims:ndims),
     $     ba(10,10,mphi,ndp),bb(10,mphi2,ndp),
     $     bd(10,4,mphi,ndp),
     $     ha(4,4,mphi,ndp),hb(4,mphi2,ndp)
      real*8 phi,cs,phim,csm,p,f,tr1(5,5),h1(5),
     $     beama(15),trd(5,5),dpndim,sigea,btr(10,10),
     $     btrd(10,10),aj,pj,dj
      integer*4 ip,kk,k1,ll,l1,ia,n,nj,j1
      parameter (nj=2)
      ia(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      call tclr(bb,10*mphi2*ndp)
      call tclr(bd,40*mphi*ndp)
      call tclr(ba,100*mphi*ndp)
      call tclr(hb,4*mphi2*ndp)
      call tclr(ha,16*mphi*ndp)
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
            call teintp(ip,f,tws,tr1,h1,1.d0,ndims,ntwissfun)
            call teintm(ip,f,trads,trd,ndims,1.d0)
            call tadd(tr1,trd,tr1,25)
            call teintb(ip,f,beams,beama,ndims)
            do k=1,4
              do k1=1,k
                kk=ia(k,k1)
                do l=1,4
                  do l1=1,l
                    ll=ia(l,l1)
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
                ha(1,k,m,j)=ha(1,k,m,j)+tr1(1,k)*csm
                ha(2,k,m,j)=ha(2,k,m,j)+tr1(2,k)*csm
                ha(3,k,m,j)=ha(3,k,m,j)+tr1(3,k)*csm
                ha(4,k,m,j)=ha(4,k,m,j)+tr1(4,k)*csm
                hb(k,m,j)=hb(k,m,j)+h1(k)*csm
                do l=1,k
                  kk=ia(k,l)
                  bb(kk,m,j)=bb(kk,m,j)+h1(k)*h1(l)*csm
                  bd(kk,1,m,j)=bd(kk,1,m,j)+
     $                 (tr1(k,1)*h1(l)+tr1(l,1)*h1(k))*csm
                  bd(kk,2,m,j)=bd(kk,2,m,j)+
     $                 (tr1(k,2)*h1(l)+tr1(l,2)*h1(k))*csm
                  bd(kk,3,m,j)=bd(kk,3,m,j)+
     $                 (tr1(k,3)*h1(l)+tr1(l,3)*h1(k))*csm
                  bd(kk,4,m,j)=bd(kk,4,m,j)+
     $                 (tr1(k,4)*h1(l)+tr1(l,4)*h1(k))*csm
                enddo
              enddo
              do k=1,10
                do l=1,10
                  ba(l,k,m,j)=ba(l,k,m,j)+btr(l,k)*csm
c     bd(l,k,m,j)=bd(l,k,m,j)+btrd(l,k)*csm
                enddo
                bb(k,m,j)=bb(k,m,j)+beama(k)*csm
              enddo
            enddo
            do m=mphi2+1,mphi
              phim=phi*(m-1)
              csm=cos(phim)/(nz-1)/(2*nj+1)
              if(i .eq. 1 .or. i .eq. nz)then
                csm=csm*.5d0
              endif
              do k=1,4
                ha(1,k,m,j)=ha(1,k,m,j)+tr1(1,k)*csm
                ha(2,k,m,j)=ha(2,k,m,j)+tr1(2,k)*csm
                ha(3,k,m,j)=ha(3,k,m,j)+tr1(3,k)*csm
                ha(4,k,m,j)=ha(4,k,m,j)+tr1(4,k)*csm
                do l=1,k
                  kk=ia(k,l)
                  bd(kk,1,m,j)=bd(kk,1,m,j)+
     $                 (tr1(k,1)*h1(l)+tr1(l,1)*h1(k))*csm
                  bd(kk,2,m,j)=bd(kk,2,m,j)+
     $                 (tr1(k,2)*h1(l)+tr1(l,2)*h1(k))*csm
                  bd(kk,3,m,j)=bd(kk,3,m,j)+
     $                 (tr1(k,3)*h1(l)+tr1(l,3)*h1(k))*csm
                  bd(kk,4,m,j)=bd(kk,4,m,j)+
     $                 (tr1(k,4)*h1(l)+tr1(l,4)*h1(k))*csm
                enddo
              enddo
              do k=1,10
                do l=1,10
                  ba(l,k,m,j)=ba(l,k,m,j)+btr(l,k)*csm
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine tespl(bx,bdx,bddx,ndp,mphi2,nd,dj,e,damp,j0,m)
      implicit none
      integer*4 ndp,mphi2,nd
      real*8 bx(nd,mphi2,ndp),
     $     bddx(nd),bdx(nd),dj,e,damp
      real*8 ajj1,ajj2
      integer*4 m,i,j0
      real*8 c1,dy1,d,ajj,f1,f2
      parameter (c1=11.d0/60.d0)
      d=-2.d0*damp
      ajj=(j0-.5d0)*dj
      f1=d*e*.5d0
      f2=d*.5d0
      if(j0 .eq. 1)then
c        ajj2=ajj+dj
        do i=1,nd
          bdx(i)=0.d0
c          bddx(i)=f2*(
c     $         e*(-ajj*bx(i,m,j0)+ajj2*bx(i,m,j0+1))/dj**2)
          bddx(i)=0.d0
        enddo
      elseif(j0 .eq. ndp)then
c        ajj1=ajj-dj
        do i=1,nd
          bdx(i)=0.d0
c          bddx(i)=f2*(
c     $         e*(ajj1*bx(i,m,j0-1)-ajj*bx(i,m,j0))/dj**2)
          bddx(i)=0.d0
        enddo
      else
        ajj1=ajj-dj
        ajj2=ajj+dj
        do i=1,nd
          dy1=(bx(i,m,j0+1)-bx(i,m,j0-1))/dj*.5d0
          bdx(i)=f1*dy1
          bddx(i)=f2*(-(e+ajj)*dy1
     $         +e*(ajj1*bx(i,m,j0-1)-2.d0*ajj*bx(i,m,j0)
     $         +ajj2*bx(i,m,j0+1))/dj**2)
        enddo
      endif
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

      subroutine tesumb(bc,hc,mphi,mphi2,ndp,
     $     w,dj,sige,tws,ndims,ntwissfun,beam,fj)
      implicit none
      integer*4 mphi,mphi2,ndp,ndims,ntwissfun
      real*8 beam(42),sige,f,dj,w,e,aj
      real*8 bc(10,mphi2,ndp),hc(4,mphi2,ndp),
     $     tws(ntwissfun,-ndims:ndims),fj(ndp)
      integer*4 i,j,m,ia,i1,ii,n
      ia(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      e=sige**2
      call tclr(beam,21)
      do j=1,ndp
        aj=(j-.5d0)*dj
        f=fj(j)*dj/e/w
        do i=1,10
          beam(i)=beam(i)+bc(i,1,j)*f
        enddo
      enddo
      do j=1,ndp
        aj=(j-.5d0)*dj
        f=fj(j)*dj/e/w*sqrt(2.d0*aj)
        do i=1,4
          beam(15+i)=beam(15+i)+hc(i,2,j)*f
        enddo
      enddo
      do i=1,4
        do i1=1,i
          ii=ia(i,i1)
          beam(ii)=beam(ii)+tws(6+i,0)*tws(6+i1,0)*e
        enddo
      enddo
      beam(16)=beam(16)+tws(7,0)*e
      beam(17)=beam(17)+tws(8,0)*e
      beam(18)=beam(18)+tws(9,0)*e
      beam(19)=beam(19)+tws(10,0)*e
      beam(21)=e
      return
      end
      
      subroutine teintb(ip,f,beams,beama,ndims)
      implicit none
      integer*4 ip,i,ip1,ndims
      real*8 f,beams(10,-ndims:ndims),beama(15),f1
      f1=1.d0-f
      ip1=ip+1
      do i=1,10
        beama(i)=f1*beams(i,ip)+f*beams(i,ip1)
      enddo
c     write(*,'(1p10g12.4)')beama
      return
      end
      
      subroutine teintm(ip,f,trads,trd,ndims,wd)
      implicit none
      integer*4 ip,ndims
      real*8 f,trads(5,5,-ndims:ndims),f1,f0,wd
      real*8 trd(5,5)
      integer*4 i,ip1
      f1=(1.d0-f)*wd
      f0=f*wd
      ip1=ip+1
      do i=1,5
        trd(1,i)=f1*trads(1,i,ip)+f0*trads(1,i,ip1)
        trd(2,i)=f1*trads(2,i,ip)+f0*trads(2,i,ip1)
        trd(3,i)=f1*trads(3,i,ip)+f0*trads(3,i,ip1)
        trd(4,i)=f1*trads(4,i,ip)+f0*trads(4,i,ip1)
        trd(5,i)=f1*trads(5,i,ip)+f0*trads(5,i,ip1)
      enddo
      return
      end
      
      subroutine teintp(ip,f,tws,tr1,h1,rs,ndims,ntwissfun)
      implicit none
      integer*4 ip,ip1,ndims,ntwissfun
      real*8 tws(ntwissfun,-ndims:ndims),tr1(5,5),h1(4)
      real*8 f,rs,f1,axi,bxi,amuxi,ayi,byi,amuyi,
     $     exi,epxi,eyi,epyi,r1i,r2i,r3i,r4i,
     $     cosmux,sinmux,cosmuy,sinmuy,amu,
     $     a11,a12,a21,a22,a33,a34,a43,a44,
     $     b11,b12,b13,b14,b21,b22,b23,b24,
     $     b31,b32,b33,b34,b41,b42,b43,b44
      ip1=ip+1
      f1=1.d0-f
      axi=f*tws(1,ip1)+f1*tws(1,ip)
      bxi=f*tws(2,ip1)+f1*tws(2,ip)
      amuxi=(f*tws(3,ip1)+f1*tws(3,ip))*rs
      ayi=f*tws(4,ip1)+f1*tws(4,ip)
      byi=f*tws(5,ip1)+f1*tws(5,ip)
      amuyi=(f*tws(6,ip1)+f1*tws(6,ip))*rs
      exi=f*tws(7,ip1)+f1*tws(7,ip)
      epxi=f*tws(8,ip1)+f1*tws(8,ip)
      eyi=f*tws(9,ip1)+f1*tws(9,ip)
      epyi=f*tws(10,ip1)+f1*tws(10,ip)
      r1i=f*tws(11,ip1)+f1*tws(11,ip)
      r2i=f*tws(12,ip1)+f1*tws(12,ip)
      r3i=f*tws(13,ip1)+f1*tws(13,ip)
      r4i=f*tws(14,ip1)+f1*tws(14,ip)
      h1(1)=f*tws(15,ip1)+f1*tws(15,ip)
      h1(2)=f*tws(16,ip1)+f1*tws(16,ip)
      h1(3)=f*tws(17,ip1)+f1*tws(17,ip)
      h1(4)=f*tws(18,ip1)+f1*tws(18,ip)
c     write(*,*)axi,bxi,amuxi,ayi,byi,amuyi,ip,f
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
      tr1(5,1)=0.d0
      tr1(5,2)=0.d0
      tr1(5,3)=0.d0
      tr1(5,4)=0.d0
      tr1(5,5)=1.d0
      return
      end
      
      subroutine tediag(trans,cod,tws,i,ndims,ntwissfun,lfno)
      implicit none
      integer*4 i,ndims,ntwissfun,lfno
      real*8 trans(6,6),cod(6),tws(ntwissfun,-ndims:ndims)
      real*8 r1,r2,r3,r4,amu,
     $     r1a,r2a,r3a,r4a,
     $     a11,a12,a21,a22,a33,a34,a43,a44,
     $     cosamux,sinamux,cosamuy,sinamuy,
     $     ax,bx,ay,by,amux,amuy
      logical stab
      call qmdiag(
     $     trans(1,1),trans(1,2),trans(1,3),trans(1,4),
     $     trans(2,1),trans(2,2),trans(2,3),trans(2,4),
     $     trans(3,1),trans(3,2),trans(3,3),trans(3,4),
     $     trans(4,1),trans(4,2),trans(4,3),trans(4,4),
     $     r1,r2,r3,r4,amu,stab,lfno)
c     amu=sqrt(1.d0-r1*r4+r2*r3)
      r1a=r1/amu
      r2a=r2/amu
      r3a=r3/amu
      r4a=r4/amu
      a11=trans(1,1)-r4a*trans(3,1)+r2a*trans(4,1)
      a12=trans(1,2)-r4a*trans(3,2)+r2a*trans(4,2)
      a21=trans(2,1)+r3a*trans(3,1)-r1a*trans(4,1)
      a22=trans(2,2)+r3a*trans(3,2)-r1a*trans(4,2)
      a33=trans(3,3)+r1a*trans(1,3)+r2a*trans(2,3)
      a34=trans(3,4)+r1a*trans(1,4)+r2a*trans(2,4)
      a43=trans(4,3)+r3a*trans(1,3)+r4a*trans(2,3)
      a44=trans(4,4)+r3a*trans(1,4)+r4a*trans(2,4)
c     write(*,*)a11,a12,a21,a22
      cosamux=(a11+a22)*.5d0
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
      tws(1,i)=ax
      tws(2,i)=bx
      tws(3,i)=amux
c     write(*,*)ax,bx,amux
      cosamuy=(a33+a44)*.5d0
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
      tws(4,i)=ay
      tws(5,i)=by
      tws(6,i)=amuy
      tws(7,i)=trans(1,6)
      tws(8,i)=trans(2,6)
      tws(9,i)=trans(3,6)
      tws(10,i)=trans(4,6)
      tws(11,i)=r1
      tws(12,i)=r2
      tws(13,i)=r3
      tws(14,i)=r4
      tws(15,i)=cod(1)
      tws(16,i)=cod(2)
      tws(17,i)=cod(3)
      tws(18,i)=cod(4)
      return
      end

      subroutine tmov65(a,b)
      implicit none
      real*8 a(6,6),b(5,5)
      integer*4 i
      do i=1,4
        b(1,i)=a(1,i)
        b(2,i)=a(2,i)
        b(3,i)=a(3,i)
        b(4,i)=a(4,i)
        b(5,i)=0.d0
      enddo
      b(1,5)=a(1,6)
      b(2,5)=a(2,6)
      b(3,5)=a(3,6)
      b(4,5)=a(4,6)
      b(5,5)=a(6,6)
      return
      end

      subroutine ttimes(a,b,c,n)
      implicit none
      integer*4 n,i
      real*8 a(n),c(n),b
      do i=1,n
        c(i)=a(i)*b
      enddo
      return
      end
