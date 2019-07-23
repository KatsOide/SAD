      module besseltab
      integer,parameter :: lbtab=1000
      real*8, parameter :: hbar=6.58211928e-16 ! eV.s
      real*8 xi(lbtab),Sw(lbtab),SwI(lbtab)
      end module besseltab

      subroutine setcolb(blist,colb)
      use wsbb
      real*8 blist(nblist)
      type(sbeam) colb

      colb%gamma=blist(1)
      colb%gambet=blist(2)
      colb%ce=blist(3)
      colb%xangle(1:8)=blist(4:11)
      colb%Benv(1:36)=blist(12:47)
      colb%Benv5(1:25)=blist(48:72)
      colb%v_cen(1:5)=blist(73:77)
      colb%cod(1:6)=blist(78:83)
      colb%Luminosity=blist(84)
      colb%bstrl=int(blist(85)+0.1)
      colb%nslice=int(blist(99)+0.1)
      colb%zslice(1:nslimax*2)=blist(100:99+nslimax*2)

      return
      end

      subroutine storecolb(blist,colb)
      use wsbb
      real*8 blist(nblist)
      type(sbeam) colb
      
!      blist=>colb%zslice
!      colb%zslice=>blist
      blist(1)=colb%gamma
      blist(2)=colb%gambet
      blist(3)=colb%ce
      blist(4:11)=colb%xangle(1:8)
      blist(12:47)=colb%Benv(1:36)
      blist(48:72)=colb%Benv5(1:25)
      blist(73:77)=colb%v_cen(1:5)
      blist(78:83)=colb%cod(1:6)
      blist(84)=colb%Luminosity
      blist(85)=dble(colb%bstrl)
      blist(99)=dble(colb%nslice)
      blist(100:99+nslimax*2)=colb%zslice(1:nslimax*2)

      return
      end

      
!     18/01/93 303061342  MEMBER NAME  BEAMBEAM *.FORT     M  E2FORT
      subroutine beambeam(np,x,px,y,py,z,g,dv,spx,spy,spz,
     $     p_in,blist,iturn)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 
!  Beam-Beam tracking subroutine with a six dimension symplectic  
!  mothod.
!  Beamstrahlung, Multi IPs 
!                          
!                   written by K.Ohmi      2015.9.2
!                          
!                                                  
!
!        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        emix  =p_in(22), emiy =(23),     dp=(24),  sigz=(27)
!        nsli  =p_in(28), rne  =(29),
!  Modified to accept spin variables spx, spy , spz by K. O., 11 Dec 2018.
!  As for blist(nblist), modified to use kfromr/rfromk to store/recall integer data, K. O. 12 Dec 2018.
      use besseltab
      use wsbb
      use tfstk
      use tmacro
      implicit none
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     spx(np),spy(np),spz(np)
      real*8 p_in(70),blist(nblist)
      real*8 sqrpi
      integer*4 np,iturn,idummy,bstrhl,nbsemit,ns1,ns2,nsli,nss
c      integer*4 n,i,is,j,jl,jm,ju,iseed,irtc
      integer*4 n,i,is,j,jl,jm,ju,irtc
      LOGICAL*4 icross,tilt,deform,cod
      real*8 a2,a2x,a2y,acx,acy,asp,c1,cdu,d1x,d1y,dux,duy,delg,dLum
      real*8 fx0,fy0,fxy,gamp,gxy,h1,hi,pn,pxy2,pz
      real*8 sigz,xx,yy,w1x,w1y
      real*8 cosx,sinx,tanx,sig11x,sig11y,sig12x,sig12y,sqrsig11,sz
      real*8 sint,cost
      real*8 cnbs,cp2bs,dp2bs,cpbs,cubs,dpzbs,dsz,rhoi,rhoi2
      real*8 rnbs,rlum_col,w2x,w2y,xij,xtmp,ucE
      real*8 tran,rgetgl
      type (sad_symdef), pointer :: symd
!      equilavence (blist(1),gamp)

!      type(sbeam), target :: colb
      type(sbeam) colb

      call setcolb(blist,colb)

      bstrhl=colb%bstrl

!      write(*,'(A,I6,1P,3E12.3)') 'Welcome to beambeam',
!     &     iturn,colb%ce,colb%gamma
!      write(*,'(4F10.5,I5)') colb%xangle(1),colb%xangle(2),
!     &     colb%xangle(3),colb%xangle(4),colb%nslice
!      write(*,'(5E12.4)') (p_in(i),i=1,30)
!      call m_print(colb%Benv,6)
!      call m_print(colb%Benv5,5)
!      write(*,'(1P,(4D12.5,A)/)') (colb%v_cen(i),i=1,4),' *z'
      
!     write(*,*) colb
!     blist=>colb
!     begin initialize for preventing compiler warning

      sqrpi=sqrt(pi)
!     end   initialize for preventing compiler warning
      nsli=colb%nslice
      if(nsli .le. 0) return
!
!   Lorentz factor gamma is stored in h0.
!   Factor gammabeta is stored in p0.
!     gamma=blist(i_p0+2)=h0
!     gammabet=blist(i_p0+1)=p0
!
!
!
      icross=.false.
      tilt=.false.
      deform=.false.
      cod=.false.
      if(colb%xangle(1).ne.0.) icross=.true.
      if(colb%xangle(5).ne.0.) tilt=.true.
      do i=1,4
         if(colb%v_cen(i).ne.0.) then
            deform=.true.
            exit
         endif
      enddo

      do i=1,4
         if(colb%cod(i).ne.0.) then
            cod=.true.
            exit
         endif
      enddo

!
!    sigz should be given dynamically for bstrhl
      sigz=p_in(27)

      gamp=colb%gamma
      cnbs=5.*sqrt(3.)*finest*gamp/6.
      cpbs=2*re*gamp*gamp*gamp/3
      cubs=1.5d0*hbar*cveloc*gamp*gamp/am_e

!      write(*,'(A,1P,4E12.4)') 'bst ',gamp,cnbs,cpbs,cubs
        

      do i=1,np

!   Get Canonical variables
!       g(i)=pn-1.d0   Not good for precision
!       g(i)=g(i)*(g(i)+2.d0)
!   pn=|p|/p0=1+delta
         pn=1.d0+g(i)
         px(i)=px(i)*pn
         py(i)=py(i)*pn
      enddo
!       write(*,*) 'Before BB interation'
!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)

! ### Crossing angle correction ######################################

      sinx=colb%xangle(2)
      cosx=colb%xangle(3)
      tanx=colb%xangle(4)
c      write(*,*) 'xangle ',sinx,cosx,tanx
c      call tfmemcheckprint('beambeam',1,.false.,irtc)

      if(icross) then
         do i=1,np
            pxy2=px(i)**2+py(i)**2
            pn=1.d0+g(i)
            pz=sqrt(pn**2-pxy2)
            h1=pxy2/(pn+pz)
            pz=1.d0/pz
            g(i)=g(i)-tanx*px(i)+tanx*tanx*h1
            px(i)=(px(i)-tanx*h1)/cosx
            py(i)=py(i)/cosx
            pxy2=px(i)**2+py(i)**2
            pz=1.d0/sqrt((1.+g(i))**2-pxy2)
            xx=x(i)
            x(i)=tanx*z(i)+(1.d0+sinx*px(i)*pz)*xx
            y(i)=y(i)+sinx*py(i)*xx*pz
            z(i)=z(i)/cosx-tanx/cosx*h1*xx*pz
         enddo
      endif
!      write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)


!  Closed orbit correction

      if(cod) then
         do i=1,np
            x(i)=x(i)-colb%cod(1)
            px(i)=px(i)+colb%cod(2)
            y(i)=y(i)-colb%cod(3)
            py(i)=py(i)+colb%cod(4)
            z(i)=z(i)+colb%cod(5)
         enddo
      endif

! t_angle rotation

      if(tilt) then
         sint=colb%xangle(5)
         cost=colb%xangle(6)
         do i=1,np
            xtmp=cost*x(i)+sint*y(i)
            y(i)=cost*y(i)-sint*x(i)
            x(i)=xtmp
            xtmp=cost*px(i)+sint*py(i)
            py(i)=cost*py(i)-sint*px(i)
            px(i)=xtmp
         enddo
      endif
!      write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)

!  Slice loop

      nss=max(0,int(rgetgl('BBCUT',idummy)))
      ns1=nss+1
      ns2=nsli-nss
      do 15 i=1,np
         do 20 is=1,nsli
!            if(is .lt. ns1 .or. is .gt. ns2)cycle
!            dsz=(colb%zslice(nslimax+is-1)-colb%zslice(nslimax+is))*
!     &           sigz*0.5
            dsz=(colb%zslice(nslimax+is-1)-colb%zslice(nslimax+is))*
     &           sigz*0.5
            if(is.eq.1) dsz=(5.-colb%zslice(nslimax+is))*sigz*0.5
            if(is.eq.colb%nslice) dsz=(colb%zslice(nslimax+is-1)+5.)*
     &           sigz*0.5
            
!       write(*,*) ' In loop ',is
            pxy2=px(i)**2+py(i)**2
            sz=(z(i)-colb%zslice(is)*sigz)/2.d0
!
!  exp(-:D:)x

            x(i)=x(i)+px(i)*sz
            y(i)=y(i)+py(i)*sz
            g(i)=g(i)-pxy2*0.25d0

!     center of mass of each slice

!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
            if(deform) then
               x(i) =x(i)-colb%zslice(is)*sigz*colb%v_cen(1)
               px(i)=px(i)-colb%zslice(is)*sigz*colb%v_cen(2)
               y(i) =y(i)-colb%zslice(is)*sigz*colb%v_cen(3)
               py(i)=py(i)-colb%zslice(is)*sigz*colb%v_cen(4)
            endif
!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
!   ###########################################################
!  Beam Beam interaction Body
!   ###########################################################
            sig11x=colb%Benv5(1)+2.d0*colb%Benv5(2)*sz 
     &           +colb%Benv5(7)*sz**2
            sig12x=colb%Benv5(2)+colb%Benv5(7)*sz
            sig11y=colb%Benv5(13)+2.d0*colb%Benv5(14)*sz
     &           +colb%Benv5(19)*sz**2
            sig12y=colb%Benv5(14)+colb%Benv5(19)*sz

            asp=sqrt(sig11y/sig11x)

!   ##### Luminosity caluculation

!      dLum=0.d0

            xx=x(i)*x(i)/sig11x
            yy=y(i)*y(i)/sig11y
            if(xx.lt.36.d0.and.yy.lt.36.d0) then 
               dLum=exp(-xx*0.5)*exp(-yy*0.5)/pi2/sqrt(sig11x*sig11y)
               colb%Luminosity=colb%Luminosity+dLum
            endif
!
!    ##### x-y kick ###########################################
!
!     case: sigma_x = sigma_y (round beam)
            if(abs((sig11x-sig11y)/sig11x).le.1.d-5) then
!      blist(1100)=gamma  blist(1101)=sli   blist(1102)=partn/sli*re
               c1=colb%ce/colb%nslice
               hi=-x(i)**2/2.d0/sig11x-y(i)**2/2.d0/sig11y
               a2=2.d0*(1.d0-exp(hi))/(x(i)**2+y(i)**2)
               fx0=c1*a2*x(i)
               fy0=c1*a2*y(i)
!   pz kick
               gxy=c1*sig12x/sig11x*exp(hi)

!    case: sigma_x > sigma_y
            elseif(sig11x.gt.sig11y) then
               sqrsig11=sqrt(2.d0*(sig11x-sig11y))
               c1=colb%ce*2.d0*sqrpi/sqrsig11/colb%nslice

               acx=abs(x(i))
               acy=abs(y(i))
               d1x=acx/sqrsig11
               d1y=acy/sqrsig11
               call cerrf(d1x,d1y,w1y,w1x)
               hi=-x(i)**2/2.d0/sig11x-y(i)**2/2.d0/sig11y
               call cerrf(d1x*asp,d1y/asp,w2y,w2x)
               a2=exp(hi)
               a2x=a2*w2x
               a2y=a2*w2y
               fx0=c1*(w1x-a2x)
               fy0=c1*(w1y-a2y)
               if(x(i) .lt. 0.d0)fx0=-fx0
               if(y(i) .lt. 0.d0)fy0=-fy0
!   pz kick
               cdu=-1.d0/(sig11x-sig11y)/2.d0
               fxy=x(i)*fx0+y(i)*fy0
               dux=cdu*(fxy+2.d0*colb%ce/colb%nslice*(a2*asp-1.d0))
               duy=-cdu*(fxy+2.d0*colb%ce/colb%nslice*(a2/asp-1.d0))
               gxy=sig12x*dux+sig12y*duy
!         write(*,*) sz,gxy,sig12x,sig12y

!     case: sigma_x < sigma_y      
            else 
               sqrsig11=sqrt(2.d0*(sig11y-sig11x))
               c1=colb%ce/colb%nslice*2.d0*sqrpi/sqrsig11

               acx=abs(x(i))
               acy=abs(y(i))
               d1x=acx/sqrsig11
               d1y=acy/sqrsig11
               call cerrf(d1y,d1x,w1x,w1y)
               hi=-x(i)**2/2.d0/sig11x-y(i)**2/2.d0/sig11y
               call cerrf(d1y/asp,d1x*asp,w2x,w2y)
               a2=exp(hi)
               a2x=a2*w2x
               a2y=a2*w2y
               fx0=c1*(w1x-a2x)
               fy0=c1*(w1y-a2y)
               if(x(i) .lt. 0.d0)fx0=-fx0
               if(y(i) .lt. 0.d0)fy0=-fy0
!     pz kick
               cdu=-1.d0/(sig11x-sig11y)/2.d0
               fxy=x(i)*fx0+y(i)*fy0
               dux=cdu*(fxy+2.d0*colb%ce/colb%nslice*(a2*asp-1.d0))
               duy=-cdu*(fxy+2.d0*colb%ce/colb%nslice*(a2/asp-1.d0))
               gxy=sig12x*dux+sig12y*duy
            endif
!
            px(i)=px(i)-fx0
            py(i)=py(i)-fy0
            g(i)=g(i)-gxy

!     Beamstrahlung

            if(bstrhl.gt.0) then
               rhoi2=(fx0*fx0+fy0*fy0)/dsz/dsz
               rhoi=sqrt(rhoi2)
               rnbs=cnbs*rhoi*dsz
               dpzbs=cpbs*rhoi2*dsz
               dp2bs=sqrt(cp2bs*rhoi*rhoi2*dsz)
!     write(*,'(I6,1P,3E12.4)') i,rnbs,dpzbs,dp2bs
               ucE=cubs*rhoi
!               write(*,'(2I3,1P,6E12.4)') is,i,rhoi,ucE,dsz,rnbs,fx0,fy0

               if(rnbs.gt.0.5) then
                  write(*,*) 'Nbs is too large',rnbs
                  stop
               endif
!               g(i)=g(i)+dpzbs  !+tgauss(iseed)*dp2bs
!  Poisson generator
c               xx=tran(iseed)
               xx=tran()
!           if(xx.gt.prb) then
!  Binomial generator
               if(xx.lt.rnbs) then
!   Quantum excitation model using 2nd order momemnt
!              g(i)=g(i)-(dpzbs+tgauss(iseed)*dp2bs*sqrt(rnbs))/rnbs
!!!!!              g(i)=g(i)-dpzbs/rnbs  !  both gave the same result
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  K bessel function table
                  xx=xx/rnbs
!     Locate j for SwI()=xx 
                  n=lbtab        ! dim of SwI
                  jl=0
                  ju=n+1
 100              if(ju-jl.gt.1) then
                     jm=(ju+jl)/2
                     if((SwI(n).gt.SwI(1)).eqv.(xx.gt.SwI(jm))) then
                        jl=jm
                     else
                        ju=jm
                     endif
                     goto 100
                  endif
                  j=jl
                  xij=xi(j)+(xi(j+1)-xi(j))*(xx-SwI(j))/
     &                 (SwI(j+1)-SwI(j))
                  delg=cubs*rhoi*exp(xij)
!                  write(*,'(I4,1P,7E12.4)') 
!     &                 j,xi(j),xi(j+1),xij,SwI(j),SwI(j+1),xx,delg
                  g(i)=g(i)-delg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                  nbsemit=nbsemit+1
               endif
!     dpzbsav=dpzbsav+dpzbs/rnbs   ! same as upper formula
!              goto 2
            endif

!     center of mass of each slice
!
!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
            if(deform) then
               x(i) =x(i)+colb%zslice(is)*sigz*colb%v_cen(1)
               px(i)=px(i)+colb%zslice(is)*sigz*colb%v_cen(2)
               y(i) =y(i)+colb%zslice(is)*sigz*colb%v_cen(3)
               py(i)=py(i)+colb%zslice(is)*sigz*colb%v_cen(4)
            endif
!
!     exp(:D:)x

            pxy2=px(i)**2+py(i)**2
            x(i)=x(i)-px(i)*sz
            y(i)=y(i)-py(i)*sz
            g(i)=g(i)+pxy2*0.25d0
!     
 20      enddo
 15   enddo
c      call tfmemcheckprint('beambeam',5,.false.,irtc)

! t_angle rotation

!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
      if(tilt) then
         do i=1,np
            xtmp=cost*x(i)-sint*y(i)
            y(i)=cost*y(i)+sint*x(i)
            x(i)=xtmp
            xtmp=cost*px(i)-sint*py(i)
            py(i)=cost*py(i)+sint*px(i)
            px(i)=xtmp
         enddo
      endif

!  Closed orbit correction

      if(cod) then
         do i=1,np
            x(i)=x(i)+colb%cod(1)
            px(i)=px(i)-colb%cod(2)
            y(i)=y(i)+colb%cod(3)
            py(i)=py(i)-colb%cod(4)
            z(i)=z(i)-colb%cod(5)
         enddo
      endif

! ###  Crossing angle correction ###################################

!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
      if(icross) then
         do i=1,np
            pxy2=px(i)**2+py(i)**2
            pn=1.d0+g(i)
            pz=sqrt(pn**2-pxy2)
            h1=pxy2/(pn+pz)
            pz=1.d0/pz
            x(i)=(x(i)-sinx*z(i))/(1.d0+(sinx*(px(i)+sinx*h1))*pz)
            y(i)=y(i)-sinx*py(i)*pz*x(i)
            z(i)=z(i)*cosx+sinx*cosx*h1*pz*x(i)
            h1=h1*cosx**2
            xx=px(i)
            px(i)=xx*cosx+tanx*h1
            py(i)=py(i)*cosx
            g(i)=g(i)+sinx*xx
         enddo
      endif

!       write(*,*) 'After BB interation'
!       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)

!   Return to SAD variables

!   pn=|p|/p0=1+delta
c      call tfmemcheckprint('beambeam',7,.false.,irtc)
      do  i=1,np
         pn=1.d0+g(i)
         px(i)=px(i)/pn
         py(i)=py(i)/pn
!     When you change momentum, change dv(i) as
         h1=sqrt(1.d0+(pn*colb%gambet)**2)
         dv(i)=-g(i)*(1.d0+pn)/h1/(h1+pn*colb%gamma)+dvfs
      enddo

!      write(*,*) 'blist(nblist)',blist(nblist)
c      call descr_sad(dfromk(int8(blist(nblist)+0.1d0)),symd)
      call descr_sad(dfromr(blist(nblist)),symd)
      rlum_col=colb%Luminosity/colb%nslice/np
c      call tfmemcheckprint('beambeam',8,.false.,irtc)
      symd%value=dfromr(rlum_col)
!     call tfsetlist(ntfreal,0,rlum_col,iax,1)
!!!!      iax=colb%iax
!!!!      rlist(iax)=rlum_col
c      call tfmemcheckprint('rlum_col',0,.false.,irtc)
      call LumRead(rlum_col)
c      call tfmemcheckprint('rlum_col',1,.false.,irtc)
!       write(*,*) rlum_col
!


      return
      end subroutine beambeam

!   18/01/93 303072023  MEMBER NAME  BBINIT   *.FORT     M  E2FORT
      subroutine bbinit(p_in,blist)
      use besseltab
      use wsbb
      use tfstk
      use tmacro
      implicit none
      integer*4 idummy,nsli,i,irtc
!  parameter (nblist=1600,tbuf0=1100,nslimax=500)
      real*8 benv_inv(36),m_sub5(25),v_sub5(5),c_sub
      real*8 m_tmp(36),m_Crs(36),eig5(10)
      real*8 m_B(36),m_R(36),m_H(36),m_HRB(36),m_Emit(36)
      real*8 p_in(70),rne,blist(nblist)
      real*8 rgetgl,gauinv,rfromk
      type (sad_descriptor) kv
      type (sad_symdef), pointer :: vsymd
!      integer*4 itlookup,itfsymbol
      real*8 sigx,sigy,sigz,az,bz,emiz,dp,pwb0,yy
      character*17 vname
      data vname/'LuminosityMonitor'/
!cccccccccc parameter input

      type(sbeam) colb

!      i8list=colb
!      write(*,*) 'Welcome to bbinit',associated(colb)
!      if(.not.associated(colb)) then
!         allocate(colb)
!         write(*,*) 'colb is allocated'
!      endif
      
      rne=p_in(29)
      nsli=int(p_in(28))
      colb%bstrl=p_in(53)
            
      colb%nslice=nsli
      if(nsli .le. 0) then
         colb%nslice=1
      endif
      if(nsli .gt. nslimax) then
         write(*,*)'*******************************'
         write(*,*)' Too many number of slices !!! '
         write(*,*)'    nslice=',nsli
         write(*,*)' Maximum allowable number is ', nslimax
         write(*,*)'*******************************'
         stop
      end if
!        alphax= p_in(1), betax= (2), alphay= (3), betay= (4)
!        r1    = p_in(5), r2   = (6),     r3= (7),   r4 = (8)
!        ex    = p_in(9), epx  =(10),     ey=(11),  epy =(12)
!        zx    =p_in(13), zpx  =(14),     zy=(15),  zpy =(16)
!        dx    =p_in(17), dpx  =(18),     dy=(19),  dpy =(20)
!        emix  =p_in(22), emiy =(23),     dp=(24),  sigz=(27)
!        nsli  =p_in(28), rne  =(29),     bstrl=(53)
!---  |----1----|----2----|----3----|----4----|----5----|----6----|----7
!       write(*,'(//A)') 
!     &   '*** Strong Beam setup *****************************'
!       WRITE(*,'(A,1PD12.6)') 'Number of Particles= ',rne
!       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
!     &    'betau      = ',p_in(2) ,'m   betav     = ',p_in(4) ,'m'
!       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
!     &    'etaX       = ',p_in(9) ,'m   etaY      = ',p_in(11),'m'
!       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
!     &    'emittance x= ',p_in(22),'m  emittance y= ',p_in(23),'m'
!       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
!     &    'sigz       = ',p_in(27),'m  sigp       = ',p_in(24),'m'
!       WRITE(*,'(A,F10.2,A)') 'Crossing angle ',p_in(21)*1.D3,' mrad'
!       if(colb%bstrl.eq.1) then
!          WRITE(*,'(A)') 'Beam strahlung ON'
!       else
!          WRITE(*,'(A)') 'Beam strahlung OFF'
!       endif
!     
!  Strong beam envelope calculation

      if(p_in(31).gt.0..and.p_in(37).gt.0..and.p_in(42).gt.0..and. 
     &     p_in(46).gt.0..and.p_in(49).gt.0..and.p_in(51).gt.0.) then
!         write(*,*) 'Strong beam is constructed from envelope matrix'
      call m_symmet_set_6d(p_in(31),colb%Benv)
      else
!         write(*,*) 'Strong beam is constructed from twiss parameters'
!      p_in(25)=<z pz>
         emiz=p_in(24)**2*p_in(27)**2-p_in(25)**2
         if(emiz.le.0.) then
            write(*,*) ' bb.f : Error in emitance z calculation'
            stop
         endif
         emiz=sqrt(emiz)
         bz=p_in(27)**2/emiz
         az=-p_in(25)/emiz
         call m_B_set_6d(p_in(1),p_in(2),p_in(3),p_in(4),az,bz,m_B)
         call m_R_set_6d(p_in(5),p_in(6),p_in(7),p_in(8),m_R)
         call m_H_set_6d(p_in(9),p_in(10),p_in(11),p_in(12),
     &        p_in(13),p_in(14),p_in(15),p_in(16),m_H)
         call m_Emit_set_6d(p_in(22),p_in(23),p_in(24)*p_in(27),m_Emit)

         call m_mul3_6d(m_H,m_R,m_B,m_HRB)
!       call is_symplectic_6d(blist(m_HRB))

!         call m_print(m_HRB,6)
         call m_sim_tr_6d(m_HRB,m_Emit,colb%Benv)
!      call m_print(blist(m_Emit),6)
      endif
!      write(*,*) ' Beam envelope matrix on the head on frame'
!      call m_print(colb%Benv,6)

      if(p_in(21).ne.0.) then
         call m_Crs_set_6d(p_in(21),m_Crs)
         call m_sim_tr_6d(m_Crs,colb%Benv,colb%Benv)
      endif

!       write(*,*) ' Envelope of Strong beam at IP (Streak image)'
!      call m_print(colb%Benv,6)
!       write(*,*) ' '
      dp=0.d0
      sigx=sqrt(colb%Benv(1))
      sigy=sqrt(colb%Benv(15))
      sigz=sqrt(colb%Benv(29))
!      WRITE(*,'(A,1PD12.6,A)') 'sigma x=',sigx ,'m',
!     &   'sigma y=',sigy,'m','sigma z=',sigz,'m'
!      com=re*partn/2.d0/pi/h0/(sigx+sigy)
!      xix=com*p_in(2)/sigx
!      xiy=com*p_in(4)/sigy
!         write(*,*) 'Convensional beam beam parameters '
!         write(*,'(A,F10.5,A,F10.5/)') 'xix=',xix,'  xiy=',xiy

!   z slicing


!   Weak beam    p0=momentum(eV)  p0+1=p0/m_e=beta gamma  p0+2=gamma

      pwb0=rgetgl('MOMENTUM',idummy)
      colb%gambet=pwb0/am_e
      colb%gamma=sqrt(colb%gambet**2+1)
!       write(*,*) ' Weak beam momentum '
!       write(*,'(3F15.3/)') (blist(i_p0+i),i=0,2)

!  Crossing angle and fundamental parameter set

      colb%xangle(1)=p_in(21)
      colb%xangle(2)=sin(p_in(21))
      colb%xangle(3)=cos(p_in(21))
      colb%xangle(4)=tan(p_in(21))
      colb%ce=rne*re/colb%gamma
!      WRITE(*,'(A,1P,E12.4,A/)') 'Ne re/gamma ',colb%ce,' mrad'

!  Closed orbit of strong beam

      colb%cod(1)=p_in(17)
      colb%cod(2)=p_in(18)/colb%xangle(3)
      colb%cod(3)=p_in(19)
      colb%cod(4)=p_in(20)/colb%xangle(3)
      colb%cod(5)=p_in(26)/colb%xangle(3)

!     write(*,*) '  Beam displacement'
!       write(*,'(A/,1P,(5D12.5)/)')
!     &    '    dx          dpx         dy         dpy        dz',
!     &            (blist(i_cod+i),i=0,4)


!-------------------------------------------------------------------
!         z slicing
!-------------------------------------------------------------------


!    get 5x5 matrix

      call m_inverse_6d(colb%Benv,Benv_inv)
!      call m_print(blist(m_benv_inv),6)
      call m_zslice_5d(Benv_inv,m_sub5,v_sub5,c_sub)
!      call m_print(blist(m_sub5),5)
      call teigen(m_sub5,m_tmp,eig5,5,5)
!      call m_print(blist(m_sub5),5)
      call m_normalize_5d(m_sub5,eig5)
!      call m_print(blist(m_sub5),5)
!      write(*,'(5D12.6)') (blist(eig5+2*i-2),i=1,5)
!      write(*,'(5D12.6)') (blist(eig5+2*i-1),i=1,5)
      do i=1,5
         if(eig5(2*(i-1)+1).lt.0.) then
            write(*,*) 'bbinit : 5dim. beam envelope has zero or',  
     &           ' negative eigen values'
             stop
          endif
          eig5(2*(i-1)+1)=1.d0/eig5(2*(i-1)+1)
       enddo
       call m_diag_5d(eig5(1),eig5(3),eig5(5),eig5(7),eig5(9),
     &      colb%Benv5)
       call m_sim_tr_5d(m_sub5,colb%Benv5,colb%Benv5)
       call mv_mul_5d(colb%Benv5,v_sub5,colb%v_cen)
       call v_chg_sign(colb%v_cen,5)
!       write(*,*) 'Strong beam envelope of each slice'
!       call m_print(colb%Benv5,5)
!       write(*,*) 'Its center of mass'
!       write(*,'(1P,(4D12.5,A)/)') (colb%v_cen(i),i=1,4),' *z'
       call benv_norm_axis(colb%Benv5,colb%v_cen,colb%xangle(5))
       colb%xangle(7)=cos(colb%xangle(5))
       colb%xangle(6)=sin(colb%xangle(5))
!      write(*,'(A,F10.2,A/)') ' Tilt angle =',
!     &     blist(i_t_angle)*1.D+3,' mrad'

!      write(*,*) 'Envelope and its center on the tilt frame'
!      call m_print(blist(m_Benv5),5)
!      write(*,'(1P,(4D12.5,A)/)') (blist(v_bcen_fac+i-1),i=1,4),' *z'

!  *****   <x_i>=v_bcen_fac(i)*z   *****
!   beambeam requires m_Benv5, v_bcen_fac, and tilt_angle


!         WRITE(*,*)'DP= ',DP
!ccccccccccc zstar(i)
! Slicing

       do i=1,nsli-1
          yy=1.d0/nsli*i
          colb%zslice(i+nslimax)=gauinv(yy)
          if (i .eq. 1)then
             colb%zslice(i)=-exp(-colb%zslice(i+nslimax)**2/2.d0)/
     &            sqrt(2.d0*pi)*nsli
          else
             colb%zslice(i)=(exp(-colb%zslice(i+nslimax-1)**2/2.d0) 
     &        -exp(-colb%zslice(i+nslimax)**2/2.d0))/sqrt(2.d0*pi)*nsli
          endif
       enddo
       if (nsli .eq. 1)then
          colb%zslice(1)=0.d0
       else
          colb%zslice(nsli)=exp(-colb%zslice(nsli+nslimax-1)**2/2.d0)
     &         /sqrt(2.d0*pi)*nsli
       endif
       colb%zslice(nslimax)=100.
       colb%zslice(nslimax+nsli)=100.
!       write(*,*)'************* Slice ***************'
!       write(*,'(A)')'Slice    zc           zboundary   (unit sigz)'
!       write(*,'(20X,F10.5)') colb%zslice(nslimax)
       do i=1,nsli
          colb%zslice(i)=-colb%zslice(i)
          colb%zslice(i+nslimax)=-colb%zslice(i+nslimax)
!          write(*,'(I3,F10.5,10X,F10.5)') i,colb%zslice(i),1./dble(nsli)
!          write(*,'(20X,F10.5)') colb%zslice(i+nslimax)
       enddo
!       write(*,*)'***********************************'

       call mkbesseltab
!
c
c If we return list the allocation should be done in each track.
c that is this if block should be remove to beambeam.
       if(itfcontext.gt.0) then
         kv=kxsymbolz(vname,len(vname),vsymd)
         call tflocald(vsymd%value)
c     In the case of only one real variable.
         blist(nblist)=rfromk(kv%k)
         vsymd%value=dfromr(0.d0)
       else
         kv%k=0
       endif

!      write(*,*) kv%k

c       call tfmemcheckprint('storecolb',0,.true.,irtc)
       call storecolb(blist,colb)
c       call tfmemcheckprint('storecolb',1,.true.,irtc)

       return
       end subroutine bbinit


      subroutine LumRead(rlum_col)
      use tfstk
c      use tfcbk
      implicit none
      type (sad_descriptor) kx,kem
      integer*4 irtc,level,itfdownlevel,isp0
      real*8 rlum_col
      save kem
      data kem%k/0/
      if(itfcontext .le. 0)then
        return
      endif
      if(kem%k .eq. 0)then
        kem=kxsymbolz('LuminosityMonitor',17)
      endif
      call tclrfpe
      levele=levele+1
      isp0=isp
      isp=isp+1
      dtastk(isp)=kem
      isp=isp+1
      rtastk(isp)=rlum_col
c      isp=isp+1
c      itastk(1,isp)=ntfreal
c      vstk(ivstkoffset+isp)=sigx
c      isp=isp+1
c      itastk(1,isp)=ntfreal
c      vstk(ivstkoffset+isp)=sigy
c      isp=isp+1
c      itastk(1,isp)=ntfreal
c      vstk(ivstkoffset+isp)=sigz

      call tfefunref(isp0+1,kx,.false.,irtc)
c      call tfdebugprint(kx,'LumRead',1)
      isp=isp0
      level=itfdownlevel()
      return
      end


      subroutine mkbesseltab
!  use physconst
      use besseltab
      implicit none
      integer i
      real*8 x,dx,xnu,ri,rk,rip,rkp,x0,x2,dx6,rk0,rk2
      real*8 Pgtot

!      write(*,*) 'Make S(x) table for radiation (beamstrahlung)'
      dx=0.1d0
      xnu=5./3.
      x=1.1**40
      Sw(1)=0.
      SwI(1)=0.
      do i=1,lbtab
      x0=x
      x=x/1.02
      x2=(x+x0)*0.5
      dx6=(x0-x)/6.
      xi(i)=log(x)
      call bessik(x2,xnu,ri,rk2,rip,rkp)
      call bessik(x,xnu,ri,rk,rip,rkp)
      if(i.gt.1) Sw(i)=Sw(i-1)+(rk0+4*rk2+rk)*dx6
      if(i.gt.1) Pgtot=x2*Sw(i)*(x0-x)
      if(i.gt.1) SwI(i)=SwI(i-1)+Sw(i)*(x0-x)
      rk0=rk
!     write(*,'(1P,6E12.4)') x,bessk0(x),rk,Sw(i),Pgtot,SwI(i)
      enddo
      do i=1,lbtab
         SwI(i)=SwI(i)/SwI(lbtab)
!     write(*,'(I4,1P,6E12.4)') i,xi(i),SwI(i)
      enddo
      return
      end subroutine mkbesseltab




