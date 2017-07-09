      subroutine tcod(trans,cod,beam,fndcod)
      use tfstk
      use ffs_flag
      use ffs_pointer
      use ffs_fit
      use ffs, only:mfitddp
      use tmacro
      implicit none
      integer*4 im,lmax
      real*8 conv0,epsr0,epsrr,rmax,fmin,a,ddpmax,red
c      real*8 dzmax
      logical fndcod
      parameter (lmax=100,
     $     conv0=1.d-10,epsr0=1.d-6,ddpmax=3.e-5,
     1     epsrr=1.d-4,rmax=1.d200,fmin=1.d-4,a=0.25d0)
      real*8 trans(6,12),cod(6),codi(6),codf(6),dcod(6),beam(42),
     1     r0,fact,trf00,dtrf0,r,dcod1(6),codw(6),conv,trs(6,6),
     $     dcod0(6),s,trw,dz,alambdarf,trf0s,v0
      integer*4 loop,i
      logical*4 isnan,rt,rtr
      real*8 , parameter :: codw0(6) =
     $     [1.d-6,1.d-5,1.d-6,1.d-5,1.d-5,1.d-6]
      if(rfsw)then
        im=6
      else
        im=5
      endif
      rtr=radcod .and. radtaper
      trf0=0.d0
      vcalpha=1.d0
      trf0s=0.d0
      v0=p0/h0
      dp0=rlist(latt(1)+mfitddp)
      rt=.false.
 10   dcod=0.d0
      dcod0=0.d0
      epsrad=epsr0
      r0=rmax
      fact=1.d0
      loop=lmax
      codi=cod
      fndcod=.true.
      codw=codw0
      conv=conv0
      trw=codw(5)
      if(radtaper)then
        if(rtr)then
          codw(6)=codw(6)*10.d0
          conv=conv*1.d3
        else
          codw(6)=codw(6)*10.d0
          conv=conv*1.d5
        endif
        trw=codw(5)
        codi(6)=dp0
      endif
 1    loop=loop-1
      if(loop .le. 0)then
        if(rtr)then
          rtr=.false.
          go to 10
        else
          fndcod=.false.
        endif
        return
      endif
      call tinitr(trans)
      codf=codi
      trf00=trf0
      call tturne(trans,codf,beam,
     $     int8(0),int8(0),int8(0),.false.,.true.,rt)
      dz=(codi(5)+codf(5))*0.5d0
      rt=radtaper
      dcod1=codi-codf
      if(radtaper)then
        dleng=dleng+dcod1(5)
        call rsetgl1('FSHIFT',-dleng/circ)
        dcod1(5)=0.d0
      endif
      r=0.d0
      do i=1,im
        r=r+(dcod1(i)/codw(i))**2
      enddo
      dtrf0=trf0-trf00
      if(.not. radcod)then
        r=r+(dtrf0/trw)**2
      endif
      write(6,'(a,1p7g12.5)')' tcod ',r,r0,fact,trf0,dtrf0,dleng
      write(6,'(1p6g12.5)')codi,codf,dcod
      if(r .lt. conv)then
        cod=codi
        return
      endif
      if(isnan(r))then
        trf0=trf0s
        loop=0
        go to 1
      endif
      if(r .ge. r0)then
        fact=fact*.5d0
        codi=cod+fact*dcod
        if(fact .lt. fmin)then
          cod=codi
          loop=0
        endif
        go to 1
      endif
      trf0s=trf0
      r0=r
      s=0.d0
      do i=1,6
        s=s+dcod(i)*dcod0(i)/codw(i)**2
      enddo
      red=r/r0
      if(red .lt. 0.7d0)then
        fact=min(fact*2.d0,1.d0)
      elseif(s .lt. 0.d0)then
        fact=max(fact*0.5d0,fmin)
      endif
      trs(:,1:6)=trans(:,1:6)
      do i=1,6
        trs(i,i)=trs(i,i)-1.d0
      enddo
      if(.not. rfsw)then
        trs(1:6,5)=0.d0
        trs(6,6)=0.d0
        dcod1(6)=0.d0
      elseif(radtaper)then
        trs(1:6,6)=0.d0
        dcod1(5)=0.d0
      endif
      dcod0=dcod
      call tsolvg(trs,dcod1,dcod,im,6,6)
      if(radcod)then
        codi(5)=codi(5)-dz-dtrf0*v0
        alambdarf=pi2/wrfeff
        trf0=trf0+dz/v0
        if(trf0 .lt. 0.d0)then
          trf0=-mod(-trf0+0.5d0*alambdarf,alambdarf)+alambdarf*0.5d0
        else
          trf0= mod(trf0-0.5d0*alambdarf,alambdarf)+alambdarf*0.5d0
        endif
      endif
      cod=codi
      codi=codi+fact*dcod
      go to 1
      end
