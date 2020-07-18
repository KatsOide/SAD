C   18/01/93 303061342  MEMBER NAME  BEAMBEAM *.FORT     M  E2FORT
        subroutine beambeam(np,x,px,y,py,z,g,dv,work,blist,iturn)
      use tfstk
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c  Beam-Beam tracking subroutine with a six dimension symplectic  c
c  mothod proposed by K. Hirata.                                  c
c                                                                 c
c                           written by K.Ohmi                     c
c                           1st version is written by Y.Funakoshi c
c                                                                 c
c                                        July 1994                c
c        with trkick                                              c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use tmacro
      use wsbb
      implicit real*8(a-h,o-z)
      parameter(m_Benv=1103,m_Benv5=1493,iv_bc_fac=1518)
      parameter(m_B5_11=1493,m_B5_12=1494,m_B5_22=1499)
      parameter(m_B5_33=1505,m_B5_34=1506,m_B5_44=1511)
      parameter(iv_cen=1518,i_t_angle=1523,i_cr_angle=1526,i_cod=1530)
      parameter(i_p0=1533)
      parameter(mean_n=1580,mean_x=1591,Luminosity=1596)
      parameter(sqrpi=m_sqrtpi,sqr2=m_sqrt2)
      
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),work(np)
      real*8 blist(*),x_mon1(6)
      type (sad_symdef), pointer :: symd
      integer*4 idummy
      LOGICAL*4 icross,tilt,deform,cod,lum_mon
      data mon_n/1000/
c     begin initialize for preventing compiler warning
      sinx=0.
      cosx=0.
      tanx=0.
c     end   initialize for preventing compiler warning
      nsli=int(blist(1101))
      if(nsli .le. 0) return
      mturn=int(blist(1590))
c
c   Lorentz factor gamma is stored in h0.
c   Factor gammabeta is stored in p0.
c     gamma=blist(i_p0+2)=h0
c     gammabet=blist(i_p0+1)=p0
c
c
c
      icross=.false.
      tilt=.false.
      deform=.false.
      cod=.false.
      if(iturn.eq.1) lum_mon=.false.
      if(iturn.eq.mturn) then
c        write(*,'(A,I8,A,I6)') 'Luminosity monitor start from ',iturn,
c     &     '     Nparticle is ',np
        lum_mon=.true.
      endif
      if(blist(i_cr_angle).ne.0.) icross=.true.
      if(blist(i_t_angle).ne.0.) tilt=.true.
      do 5 i=1,4
        if(blist(iv_cen+i-1).ne.0.) then
           deform=.true.
           goto 6
        endif
 5      continue
 6      continue
      do 7 i=1,4
        if(blist(i_cod+i-1).ne.0.) then
           cod=.true.
           goto 8
        endif
 7    continue
 8    continue
c        
c
      do 11 i=1,np
c
c   Get Canonical variables
c       g(i)=pn-1.d0   Not good for precision
c       g(i)=g(i)*(g(i)+2.d0)
c   pn=|p|/p0=1+delta
       pn=1.d0+g(i)
       px(i)=px(i)*pn
       py(i)=py(i)*pn
 11    continue
c       write(*,*) 'Before BB interation'
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
c
c  Beam size monitor
c
      do 21 i=1,6
 21      x_mon1(i)=0.
      do 22 i=1,np
       x_mon1(1)=x_mon1(1)+x(i)
       x_mon1(2)=x_mon1(2)+y(i)
       x_mon1(3)=x_mon1(3)+x(i)*x(i)
       x_mon1(4)=x_mon1(4)+y(i)*x(i)
 22    x_mon1(5)=x_mon1(5)+y(i)*y(i)
c
c ### Crossing angle correction ######################################
c
       if(icross) then
         sinx=blist(i_cr_angle+1)
         cosx=blist(i_cr_angle+2)
         tanx=blist(i_cr_angle+3)
         do 12 i=1,np
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
 12      continue
       endif
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
c
c
c  Closed orbit correction
c
      if(cod) then
        do 13 i=1,np
        x(i)=x(i)-blist(i_cod)
        px(i)=px(i)+blist(i_cod+1)
        y(i)=y(i)-blist(i_cod+2)
        py(i)=py(i)+blist(i_cod+3)
        z(i)=z(i)+blist(i_cod+4)
 13     continue
      endif
c
c t_angle rotation
c
      if(tilt) then
        do 14 i=1,np
        xtmp=blist(i_t_angle+1)*x(i)+blist(i_t_angle+2)*y(i)
        y(i)=blist(i_t_angle+1)*y(i)-blist(i_t_angle+2)*x(i)
        x(i)=xtmp
        xtmp=blist(i_t_angle+1)*px(i)+blist(i_t_angle+2)*py(i)
        py(i)=blist(i_t_angle+1)*py(i)-blist(i_t_angle+2)*px(i)
        px(i)=xtmp
 14     continue
      endif
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
c
c  Slice loop
c
       nss=max(0,int(rgetgl('BBCUT',idummy)))
       ns1=nss+1
       ns2=nsli-nss
       do 15 i=1,np
       do 20 is=1,nsli
        if(is .lt. ns1 .or. is .gt. ns2)cycle
c
c       write(*,*) ' In loop ',is
        pxy2=px(i)**2+py(i)**2
        sz=(z(i)-blist(is))/2.d0
c
c  exp(-:D:)x
c
        x(i)=x(i)+px(i)*sz
        y(i)=y(i)+py(i)*sz
        g(i)=g(i)-pxy2*0.25d0
c
c     center of mass of each slice
c
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
       if(deform) then
           x(i) =x(i)-blist(is)*blist(iv_cen)
           px(i)=px(i)-blist(is)*blist(iv_cen+1)
           y(i) =y(i)-blist(is)*blist(iv_cen+2)
           py(i)=py(i)-blist(is)*blist(iv_cen+3)
       endif
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
c   ###########################################################
c  Beam Beam interaction Body
c   ###########################################################
       sig11x=blist(m_B5_11)+2.d0*blist(m_B5_12)*sz
     &       +blist(m_B5_22)*sz**2
       sig12x=blist(m_B5_12)+blist(m_B5_22)*sz
       sig11y=blist(m_B5_33)+2.d0*blist(m_B5_34)*sz
     &       +blist(m_B5_44)*sz**2
       sig12y=blist(m_B5_34)+blist(m_B5_44)*sz
c
       asp=sqrt(sig11y/sig11x)
c
c   ##### Luminosity caluculation
c
c      dLum=0.d0
c
      xx=x(i)*x(i)/sig11x
      yy=y(i)*y(i)/sig11y
      if(xx.lt.36.d0.and.yy.lt.36.d0) then 
         dLum=exp(-xx*0.5)*exp(-yy*0.5)/pi2/sqrt(sig11x*sig11y)
         x_mon1(6)=x_mon1(6)+dLum
      endif
c
c    ##### x-y kick ###########################################
c
c     case: sigma_x = sigma_y (round beam)
      if(abs((sig11x-sig11y)/sig11x).le.1.d-5) then
c      blist(1100)=gamma  blist(1101)=sli   blist(1102)=partn/sli*re
         c1=blist(1102)/blist(1100)
         hi=-x(i)**2/2.d0/sig11x-y(i)**2/2.d0/sig11y
         a2=2.d0*(1.d0-exp(hi))/(x(i)**2+y(i)**2)
         fx0=c1*a2*x(i)
         fy0=c1*a2*y(i)
c   pz kick
         gxy=c1*sig12x/sig11x*exp(hi)

c    case: sigma_x > sigma_y
      elseif(sig11x.gt.sig11y) then
         sqrsig11=sqrt(2.d0*(sig11x-sig11y))
         c1=blist(1102)*2.d0*sqrpi/sqrsig11/blist(1100)

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
c   pz kick
         cdu=-1.d0/(sig11x-sig11y)/2.d0
         fxy=x(i)*fx0+y(i)*fy0
         dux=cdu*(fxy+2.d0*blist(1102)/blist(1100)*(a2*asp-1.d0))
         duy=-cdu*(fxy+2.d0*blist(1102)/blist(1100)*(a2/asp-1.d0))
         gxy=sig12x*dux+sig12y*duy
!         write(*,*) sz,gxy,sig12x,sig12y
c 
c     case: sigma_x < sigma_y      
      else 
         sqrsig11=sqrt(2.d0*(sig11y-sig11x))
         c1=blist(1102)*2.d0*sqrpi/sqrsig11/blist(1100)

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
c   pz kick
         cdu=-1.d0/(sig11x-sig11y)/2.d0
         fxy=x(i)*fx0+y(i)*fy0
         dux=cdu*(fxy+2.d0*blist(1102)/blist(1100)*(a2*asp-1.d0))
         duy=-cdu*(fxy+2.d0*blist(1102)/blist(1100)*(a2/asp-1.d0))
         gxy=sig12x*dux+sig12y*duy
      endif
c 
      px(i)=px(i)-fx0
      py(i)=py(i)-fy0
      g(i)=g(i)-gxy
c
c
c     center of mass of each slice
c
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
       if(deform) then
           x(i) =x(i)+blist(is)*blist(iv_cen)
           px(i)=px(i)+blist(is)*blist(iv_cen+1)
           y(i) =y(i)+blist(is)*blist(iv_cen+2)
           py(i)=py(i)+blist(is)*blist(iv_cen+3)
       endif
c
c     exp(:D:)x
c
        pxy2=px(i)**2+py(i)**2
        x(i)=x(i)-px(i)*sz
        y(i)=y(i)-py(i)*sz
        g(i)=g(i)+pxy2*0.25d0
c
 20   continue
 15   continue
c
c t_angle rotation
c
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
      if(tilt) then
        do 16 i=1,np
        xtmp=blist(i_t_angle+1)*x(i)-blist(i_t_angle+2)*y(i)
        y(i)=blist(i_t_angle+1)*y(i)+blist(i_t_angle+2)*x(i)
        x(i)=xtmp
        xtmp=blist(i_t_angle+1)*px(i)-blist(i_t_angle+2)*py(i)
        py(i)=blist(i_t_angle+1)*py(i)+blist(i_t_angle+2)*px(i)
        px(i)=xtmp
 16     continue
      endif
c
c  Closed orbit correction
c
      if(cod) then
        do 17 i=1,np
        x(i)=x(i)+blist(i_cod)
        px(i)=px(i)-blist(i_cod+1)
        y(i)=y(i)+blist(i_cod+2)
        py(i)=py(i)-blist(i_cod+3)
        z(i)=z(i)-blist(i_cod+4)
 17     continue
      endif
c
c ###  Crossing angle correction ###################################
c
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
      if(icross) then
        do 18 i=1,np
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
 18     continue
      endif
c
c       write(*,*) 'After BB interation'
c       write(*,'(2D25.15)') x(1),px(1),y(1),py(1),z(1),g(1)
c
c   Return to SAD variables
c
c   pn=|p|/p0=1+delta
      do 19 i=1,np
      pn=1.d0+g(i)
      px(i)=px(i)/pn
      py(i)=py(i)/pn
c   When you change momentum, change dv(i) as
      h1=sqrt(1.d0+(pn*blist(i_p0+1))**2)
      dv(i)=-g(i)*(1.d0+pn)/h1/(h1+pn*blist(i_p0+2))+dvfs
 19   continue
c
      continue
c
c    Luminosity Monitor
c
c
c       write(*,'(1P(7D11.4))') x_mon1,blist(Luminosity)
       if(mod(iturn,mon_n).eq.1) then
        do 1 i=1,6
 1         blist(mean_n+i-1)=x_mon1(i)/np
       else
        do 27 i=1,6
         blist(mean_n+i-1)=blist(mean_n+i-1)+x_mon1(i)/np
 27      continue
       endif
       call descr_sad(dfromk(int8(blist(nblist)+0.1d0)),symd)
       rlum_col=x_mon1(6)/nsli/np
c       call tfsetlist(ntfreal,0,rlum_col,iax,1)
       symd%value=dfromr(rlum_col)
c       write(*,*) 'Calling LumRead ',rlum_col
       call LumRead(rlum_col)
c
      if(mod(iturn,mon_n).eq.0) then
        do 28 i=1,6
 28        blist(mean_n+i-1)=blist(mean_n+i-1)/mon_n
c         if(mod(iturn,1000).eq.0) then
c         write(*,*) ' Lum. monitor ',
c     &      iturn-mon_n+1,iturn,'  np=',np
c         write(*,'(A,A)') '      <x>          <y>         <xx>',
c     &      '         <xy>         <yy>      Lum(m-2)'
c        write(*,'(1P,(6E13.5)/)') (blist(mean_n+i),i=0,4),
c     &      blist(mean_n+5)/nsli
c        endif
       endif
c
      if(lum_mon) then
      do 24 i=1,5
 24      blist(mean_x+i-1)=blist(mean_x+i-1)+x_mon1(i)/np
         blist(Luminosity)=blist(Luminosity)+x_mon1(6)/np/nsli
c
      if(mod(iturn,1000).eq.0) then
c       write(*,*) 'mean values ',mturn,'-',iturn
c       write(*,'(A,A)') '      <x>          <y>         <xx>',
c     &      '         <xy>         <yy>      Lum(m-2)'
c        write(*,'(1P,(6E13.5))') (x_mon1(i),i=1,5),x_mon1(6)/nsli
       do 26 i=1,5
 26       x_mon1(i)=blist(mean_x+i-1)/(iturn-mturn+1)
          x_mon1(6)=blist(Luminosity)/(iturn-mturn+1)
c        write(*,'(1P,(6E13.5)/)') (x_mon1(i),i=1,6)
      endif
c
      if(iturn.eq.nturn) then
        do 25 i=1,5
 25        blist(mean_x+i-1)=blist(mean_x+i-1)/(nturn-mturn+1)
        blist(Luminosity)=blist(Luminosity)/
     &         (nturn-mturn+1)/1.e4
c        write(*,'(/A,1PE15.7,A/)') 'Luminosity/Ns/Nw/Ncoll = ',
c     &      blist(Luminosity),' cm-2'
c       write(*,'(A,A)') '       <x>            <y>           <xx>',
c     &      '           <xy>           <yy>'
c        write(*,'(1P,(5E15.7)/)') (blist(mean_x+i),i=0,4)
      endif
      endif
c
      return
      end
c
C   18/01/93 303072023  MEMBER NAME  BBINIT   *.FORT     M  E2FORT
      subroutine bbinit(p_in,blist)
      use tfstk
      use tmacro
      use wsbb, only:nblist
      implicit real*8(a-h,o-z)
c     elradi:	Classical electron radius
c     elmass:	Electron mass energy equivalent in eV
      parameter (am_e=elmass,re=elradi)
      integer*4 tbuf0,idummy,nslimax
      parameter (tbuf0=1100,nslimax=500)
      real*8 p_in(70),blist(nblist),rfromk
      type (sad_descriptor) kv
      type (sad_symdef), pointer :: vsymd
      integer v_sub5,c_sub,eig5,v_bcen_fac
      character*17 vname
      data vname/'LuminosityMonitor'/
ccccccccccc parameter input
        partn=p_in(29)
          sli=p_in(28)
           nsli=int(sli)
        if(nsli .le. 0) then
          blist(tbuf0+1)=sli
          return
        endif
        if(nsli .gt. nslimax) then
         write(*,*)'*******************************'
         write(*,*)' Too many number of slices !!! '
         write(*,*)'    nslice=',nsli
         write(*,*)' Maximum allowable number is ', nslimax
         write(*,*)'*******************************'
         stop
        end if
c        alphax= p_in(1), betax= (2), alphay= (3), betay= (4)
c        r1    = p_in(5), r2   = (6),     r3= (7),   r4 = (8)
c        ex    = p_in(9), epx  =(10),     ey=(11),  epy =(12)
c        zx    =p_in(13), zpx  =(14),     zy=(15),  zpy =(16)
c        dx    =p_in(17), dpx  =(18),     dy=(19),  dpy =(20)
c        emix  =p_in(22), emiy =(23),     dp=(24),  sigz=(27)
c---|----1----|----2----|----3----|----4----|----5----|----6----|----7
c       write(*,'(//A)') 
c     &   '*** Strong Beam setup *****************************'
c       WRITE(*,'(A,1PD12.6)') 'Number of Particles= ',partn
c 603      format(e10.5)
c       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
c     &    'betau      = ',p_in(2) ,'m   betav     = ',p_in(4) ,'m'
c       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
c     &    'etaX       = ',p_in(9) ,'m   etaY      = ',p_in(11),'m'
c       WRITE(*,'(A,1PD12.6,A,D12.6,A)')
c     &    'emittance x= ',p_in(22),'m  emittance y= ',p_in(23),'m'
c       WRITE(*,'(A,F10.2,A/)') 'Crossing angle ',p_in(21)*1.D3,' mrad'
c
c  Strong beam envelope calculation
c
       m_Benv=tbuf0+3
       m_Emit=m_benv+36
       m_B=m_Emit+36
       m_R=m_B+36
       m_H=m_R+36
       m_HRB=m_H+36
       m_BRHt=m_HRB+36
       m_Crs=m_BRHt+36
       if(p_in(31).gt.0..and.p_in(37).gt.0..and.p_in(42).gt.0..and.
     &    p_in(46).gt.0..and.p_in(49).gt.0..and.p_in(51).gt.0.) then
          write(*,*) 'Strong beam is constructed from envelope matrix'
         call m_symmet_set_6d(p_in(31),blist(m_Benv))
       else
         write(*,*) 'Strong beam is constructed from twiss parameters'
c      p_in(25)=<z pz>
       emiz=p_in(24)**2*p_in(27)**2-p_in(25)**2
       if(emiz.lt.0.) then
          write(*,*) ' bb.f : Error in emitance z calculation'
          stop
       endif
       emiz=sqrt(emiz)
       bz=p_in(27)**2/emiz
       az=-p_in(25)/emiz
       call m_B_set_6d(p_in(1),p_in(2),p_in(3),p_in(4),az,bz,
     +      blist(m_B))
       call m_R_set_6d(p_in(5),p_in(6),p_in(7),p_in(8),
     +      blist(m_R))
       call m_H_set_6d(p_in(9),p_in(10),p_in(11),p_in(12),
     +      p_in(13),p_in(14),p_in(15),p_in(16),blist(m_H))
       call m_Emit_set_6d(p_in(22),p_in(23),p_in(24)*p_in(27),
     +      blist(m_Emit))
c
       call m_mul3_6d(blist(m_H),blist(m_R),blist(m_B),
     +      blist(m_HRB))
c       call is_symplectic_6d(blist(m_HRB))
c
c       call m_print(blist(m_HRB),6)
       call m_sim_tr_6d(blist(m_HRB),blist(m_Emit),
     +      blist(m_Benv))
c      call m_print(blist(m_Emit),6)
      endif
c      write(*,*) ' Beam envelope matrix on the head on frame'
c      call m_print(blist(m_Benv),6)
c
      if(p_in(21).ne.0.) then
        call m_Crs_set_6d(p_in(21),blist(m_Crs))
        call m_sim_tr_6d(blist(m_Crs),blist(m_Benv),blist(m_Benv))
      endif
c
c       write(*,*) ' Envelope of Strong beam at IP (Streak image)'
c       call m_print(blist(m_Benv),6)
c       write(*,*) ' '
          dp=0.d0
      sigx=sqrt(blist(m_Benv))
      sigy=sqrt(blist(m_Benv+2*6+2))
      sigz=sqrt(blist(m_Benv+4*6+4))
c      WRITE(*,'(A,1PD12.6,A)') 'sigma x=',sigx ,'m',
c     &   'sigma y=',sigy,'m','sigma z=',sigz,'m'
c      com=re*partn/2.d0/pi/h0/(sigx+sigy)
c      xix=com*p_in(2)/sigx
c      xiy=com*p_in(4)/sigy
c         write(*,*) 'Convensional beam beam parameters '
c         write(*,'(A,F10.5,A,F10.5/)') 'xix=',xix,'  xiy=',xiy
c
c   z slicing
c
      m_Benv_inv=m_Crs+36
      m_tmp=m_Benv_inv+36
      m_sub5=m_tmp+25
      v_sub5=m_sub5+25
      c_sub=v_sub5+5
      eig5=c_sub+1
      m_Benv5=eig5+10
      v_bcen_fac=m_Benv5+25
      i_t_angle=v_bcen_fac+5
      i_cr_angle=i_t_angle+3
      i_cod=i_cr_angle+4
      i_p0=i_cod+5
      i_blmax=i_p0+3
      if(i_blmax.gt.nblist) then
         write(*,*) 
     &   'bbinit : Number of variables exceed the allocation',i_blmax
         stop
      endif
c
c   Weak beam    p0=momentum(eV)  p0+1=p0/m_e=beta gamma  p0+2=gamma
c
       blist(i_p0)=rgetgl('MOMENTUM',idummy)
       blist(i_p0+1)=blist(i_p0)/amass
       blist(i_p0+2)=sqrt(blist(i_p0+1)**2+1)
c       write(*,*) ' Weak beam momentum '
c       write(*,'(3F15.3/)') (blist(i_p0+i),i=0,2)
c
c  Crossing angle and fundamental parameter set
c
        blist(i_cr_angle)=p_in(21)
        blist(i_cr_angle+1)=sin(p_in(21))
        blist(i_cr_angle+2)=cos(p_in(21))
        blist(i_cr_angle+3)=tan(p_in(21))
        blist(tbuf0)=blist(i_p0+1)
        blist(tbuf0+1)=sli
        blist(tbuf0+2)=partn/sli*re
c
c  Closed orbit of strong beam
c
       blist(i_cod)=p_in(17)
       blist(i_cod+1)=p_in(18)/blist(i_cr_angle+2)
       blist(i_cod+2)=p_in(19)
       blist(i_cod+3)=p_in(20)/blist(i_cr_angle+2)
       blist(i_cod+4)=p_in(26)/blist(i_cr_angle+2)
c       write(*,*) '  Beam displacement'
c       write(*,'(A/,1P,(5D12.5)/)')
c     &    '    dx          dpx         dy         dpy        dz',
c     &            (blist(i_cod+i),i=0,4)
c
c
c-------------------------------------------------------------------
c         z slicing
c-------------------------------------------------------------------
c
c
c    get 5x5 matrix
c
      call m_inverse_6d(blist(m_Benv),blist(m_benv_inv))
c      call m_print(blist(m_benv_inv),6)
      call m_zslice_5d(blist(m_Benv_inv),blist(m_sub5),
     +    blist(v_sub5),blist(c_sub))
c      call m_print(blist(m_sub5),5)
      call teigen(blist(m_sub5),blist(m_tmp),blist(eig5),5,5)
c      call m_print(blist(m_sub5),5)
      call m_normalize_5d(blist(m_sub5),blist(eig5))
c      call m_print(blist(m_sub5),5)
c      write(*,'(5D12.6)') (blist(eig5+2*i-2),i=1,5)
c      write(*,'(5D12.6)') (blist(eig5+2*i-1),i=1,5)
      do 15 i=1,5
         if(blist(eig5+2*(i-1)).lt.0.) then
            write(*,*) 'bbinit : 5dim. beam envelope has zero or', 
     &         ' negative eigen values'
            stop
         endif
 15      blist(eig5+2*(i-1))=1.d0/blist(eig5+2*(i-1))
      call m_diag_5d(blist(eig5),blist(eig5+2),blist(eig5+4),
     +    blist(eig5+6),blist(eig5+8),blist(m_Benv5))
      call m_sim_tr_5d(blist(m_sub5),blist(m_Benv5),
     +    blist(m_Benv5))
      call mv_mul_5d(blist(m_Benv5),blist(v_sub5),
     +    blist(v_bcen_fac))
      call v_chg_sign(blist(v_bcen_fac),5)
c      write(*,*) 'Strong beam envelope of each slice'
c      call m_print(blist(m_Benv5),5)
c      write(*,*) 'Its center of mass'
c      write(*,'(1P,(4D12.5,A)/)') (blist(v_bcen_fac+i-1),i=1,4),' *z'
      call benv_norm_axis(blist(m_Benv5),blist(v_bcen_fac),
     $    blist(i_t_angle))
      blist(i_t_angle+1)=cos(blist(i_t_angle))
      blist(i_t_angle+2)=sin(blist(i_t_angle))
c      write(*,'(A,F10.2,A/)') ' Tilt angle =',
c     &     blist(i_t_angle)*1.D+3,' mrad'
c
c      write(*,*) 'Envelope and its center on the tilt frame'
c      call m_print(blist(m_Benv5),5)
c      write(*,'(1P,(4D12.5,A)/)') (blist(v_bcen_fac+i-1),i=1,4),' *z'
c
c  *****   <x_i>=v_bcen_fac(i)*z   *****
c   beambeam requires m_Benv5, v_bcen_fac, and tilt_angle
c
c---|----1----|----2----|----3----|----4----|----5----|----6----|----7-
c         WRITE(*,*)'DP= ',DP
cccccccccccc zstar(i)
c Slicing
c
       do 10 i=1,nsli-1
         yy=1.d0/nsli*i
         blist(i+nslimax)=gauinv(yy)
        if (i .eq. 1)then
         blist(i)=-exp(-blist(i+nslimax)**2/2.d0)/sqrt(2.d0*pi)*nsli
        else
         blist(i)=(exp(-blist(i+nslimax-1)**2/2.d0)
     &            -exp(-blist(i+nslimax)**2/2.d0))/sqrt(2.d0*pi)*nsli
        endif
  10   continue
        if (nsli .eq. 1)then
         blist(1)=0.d0
        else
         blist(nsli)=exp(-blist(nsli+nslimax-1)**2/2.d0)/
     &          sqrt(2.d0*pi)*nsli
        endif
c         write(*,*)'************* Slice ***************'
       do 11 i=1,nsli
         blist(i)=-blist(i)
c         write(*,*) i,blist(i)
         blist(i)=blist(i)*sigz
c        if(i.ne.nsli) write(*,*) '   ',blist(i+nslimax)
  11   continue
c         write(*,*)'***********************************'
c
c   Luminosity monitor
c
       blist(1590)=p_in(30)
       if(blist(1590).lt.0.99) blist(1590)=1.d0
       do 12 i=1,10
 12       blist(1590+i)=0.
c
c If we return list the allocation should be done in each track.
c that is this if block should be remove to beambeam.
       if(itfcontext.gt.0) then
         kv=kxsymbolz(vname,len(vname),vsymd)
         call tflocald(vsymd%value)
c          iax=itaaloc(0,1)
c          iax=itaaloc(0,len_list) nlist=length of the list
c          blist(1600)=iax
c          ilist(1,iv-1)=ntflist
c          ilist(2,iv-1)=itfcopy(ntflist,iax)
c     In the case of only one real variable.
          blist(nblist)=rfromk(kv%k)
          vsymd%value=dfromr(0.d0)
       else
         kv%k=0
       endif
       return
         end
c
      subroutine LumRead(rlum_col)
      use tfstk
      use tfcbk
      use efun
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

      kx=tfefunref(isp0+1,.false.,irtc)
c      call tfdebugprint(kx,'LumRead',1)
      isp=isp0
      level=itfdownlevel()
      return
      end
c
C   09/02/93 305021512  MEMBER NAME  GAUINV   *.FORT     M  E2FORT
      FUNCTION GAUINV(P0)
C  INVERSE OF (INTEGRATED) NORMAL DISTRIBUTION FUNCTION
C              1         X= Y
C     P(Y)=-----------* INTEGRAL EXP(-X**2/2) DX
C          SQRT(2*PI)    X= -INF
C     IF P(Y)=P0, THEN GAUINV(P0)=Y.
C        0 < P0 < 1 ,   -INF < Y < +INF
C  IF THIS ROUTINE IS USED TO CONVERT UNIFORM RANDOM NUMBERS TO
C  GAUSSIAN, MAXIMUM RELATIVE ERROR IN THE DISTRIBUTION FUNCTION
C  DP/DX=EXP(-X**2/2)/SQRT(2*PI) IS LESS THAN 0.640E-3 EVERYWHERE
C  IN THE RANGE  2**(-31) < P0 < 1-2**31.  (MINIMAX APPROXIMATION)
C
      IMPLICIT REAL*8(A-H,O-Z)
C------------------------
      DATA PP1/0.334624883253D0/, QQ2/0.090230446775D0/,
     1     QQ3/0.049905685242D0/, QQ4/0.027852994157D0/,
     2     QQ5/0.015645650215D0/
      DATA A3/ 4.5585614D+01/, A2/ 2.1635544D+00/, A1/ 2.7724523D+00/,
     1     A0/ 2.5050240D+00/,
     2     B4/ 4.0314354D+02/, B3/-2.7713713D+02/, B2/ 7.9731883D+01/,
     3     B1/-1.4946512D+01/, B0/ 2.2157257D+00/,
     4     C4/ 4.1394487D+03/, C3/-1.5585873D+03/, C2/ 2.4648581D+02/,
     5     C1/-2.4719139D+01/, C0/ 2.4335936D+00/,
     6     D4/ 4.0895693D+04/, D3/-8.5400893D+03/, D2/ 7.4942805D+02/,
     7     D1/-4.1028898D+01/, D0/ 2.6346872D+00/,
     8     E4/ 3.9399134D+05/, E3/-4.6004775D+04/, E2/ 2.2566998D+03/,
     9     E1/-6.8317697D+01/, E0/ 2.8224654D+00/,
     O     F0/-8.1807613D-02/, F1/-2.8358733D+00/, F2/ 1.4902469D+00/
C------------------------
      GAUINV=0.D0
      P=P0-0.5D0
      P1=ABS(P)
      IF(P1.GE.PP1) GOTO 120
      P2=P**2
      GAUINV=(((A3*P2+A2)*P2+A1)*P2+A0)*P
      RETURN
  120 Q=0.5D0-P1
      IF(Q.LE.QQ2) GOTO 140
      GAUINV=(((B4*Q+B3)*Q+B2)*Q+B1)*Q+B0
      GOTO 200
  140 IF(Q.LE.QQ3) GOTO 150
      GAUINV=(((C4*Q+C3)*Q+C2)*Q+C1)*Q+C0
      GOTO 200
  150 IF(Q.LE.QQ4) GOTO 160
      GAUINV=(((D4*Q+D3)*Q+D2)*Q+D1)*Q+D0
      GOTO 200
  160 IF(Q.LE.QQ5) GOTO 170
      GAUINV=(((E4*Q+E3)*Q+E2)*Q+E1)*Q+E0
      GOTO 200
  170 IF(Q.LE.0D0) GOTO 900
      T=SQRT(-2D0*LOG(Q))
      GAUINV=T+F0+F1/(F2+T)
  200 IF(P.LT.0D0) GAUINV=-GAUINV
      RETURN
  900 WRITE(6,910) P0
  910 FORMAT(' (FUNC.GAUINV) INVALID INPUT ARGUMENT ',1PD20.13)
      RETURN
      END
C
C CERRF *************************************************
      SUBROUTINE CERRF(XX, YY, WX, WY)
*----------------------------------------------------------------------*
* Purpose:                                                             *
*   Modification of WWERF, double precision complex error function,    *
*   written at CERN by K. Koelbig.                                     *
* Input:                                                               *
*   XX, YY    (real)    Argument to CERF.                              *
* Output:                                                              *
*   WX, WY    (real)    Function result.                               *
*----------------------------------------------------------------------*

*---- Single precision version.
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (I-N)
C     REAL RWX, RWY
      PARAMETER         (MWFLT = 1, MREAL = 3)
*     IMPORTANT: MWNAM must be an integer multiple of MWFLT.
      PARAMETER         (MCWRD = 4)
      PARAMETER         (MCNAM = 16, MWNAM = MCNAM / MCWRD)

      PARAMETER         (CC     = 1.12837 91670 9551D0)
      PARAMETER         (ONE    = 1.D0)
      PARAMETER         (TWO    = 2.D0)
      PARAMETER         (XLIM   = 5.33D0)
      PARAMETER         (YLIM   = 4.29D0)
      DIMENSION         RX(33), RY(33)

c     XX=REAL(ZI)
c     YY=AIMAG(ZI)

      X = ABS(XX)
      Y = ABS(YY)

      IF (Y .LT. YLIM  .AND.  X .LT. XLIM) THEN
         Q = (ONE - Y / YLIM) * SQRT(ONE - (X/XLIM)**2)
         H = ONE / (3. 2* Q)
         NC =  7+ INT(23.0*Q)
         XL = H**( 1- NC)
         XH = Y + 0.5/H
         YH = X
         NU =  10+ INT(21.0*Q)
         RX(NU+1) = 0.
         RY(NU+1) = 0.

         DO 10 N = NU, 1, -1
            TX = XH + N * RX(N+1)
            TY = YH - N * RY(N+1)
            TN = TX*TX + TY*TY
            RX(N) = 0. 5* TX / TN
            RY(N) = 0. 5* TY / TN
   10    CONTINUE

         SX = 0.
         SY = 0.

         DO 20 N = NC, 1, -1
            SAUX = SX + XL
            SX = RX(N) * SAUX - RY(N) * SY
            SY = RX(N) * SY + RY(N) * SAUX
            XL = H * XL
   20    CONTINUE

         WX = CC * SX
         WY = CC * SY
      ELSE
         XH = Y
         YH = X
         RX(1) = 0.
         RY(1) = 0.

         DO 30 N = 9, 1, -1
            TX = XH + N * RX(1)
            TY = YH - N * RY(1)
            TN = TX*TX + TY*TY
            RX(1) = 0. 5* TX / TN
            RY(1) = 0. 5* TY / TN
   30    CONTINUE

         WX = CC * RX(1)
         WY = CC * RY(1)
      ENDIF

      IF(Y .EQ. 0.) WX = EXP(-X**2)
      IF(YY .LT. 0.) THEN
         WX = TWO * EXP(Y*Y-X*X) * COS(TWO*X*Y) - WX
         WY = - TWO * EXP(Y*Y-X*X) * SIN(TWO*X*Y) - WY
         IF(XX .GT. 0.) WY = -WY
      ELSE
         IF(XX .LT. 0.) WY = -WY
      ENDIF
*      print *,'wx= ',wx,' wy= ',wy
C     RWX = WX
C     RWY = WY
*      print *,'rwx= ',rwx,' rwy= ',rwy

      RETURN
      END
