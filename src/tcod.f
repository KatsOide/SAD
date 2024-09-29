      module codm
      real*8, parameter :: 
     $     conv0=1.d-4,epsr0=1.d-6,ddpmax=3.e-5,
     1     epsrr=1.d-4,rmax=1.d200,fmin=1.d-4,a=0.25d0, dpthre=3.e-4
      real*8 , parameter :: codw0(6)=[1.d-6,1.d-7,1.d-7,1.d-8,1.d-5,1.d-6]
      integer*8 ::kcodconv=0

      contains
      subroutine tgetcodw(codw)
      use tfstk
      implicit none
      real*8 ,intent(out):: codw(6)
      type (sad_rlist),pointer :: lcodw
      integer*4 nc
      real*8 w
      if(kcodconv == 0)then
        kcodconv=ktfsymbolz('`CODCONV',8)-4
      endif
      if(tfreallistq(dlist(kcodconv),lcodw) .and. lcodw%nl > 0)then
        codw=codw0
        nc=min(6,lcodw%nl)
        codw(1:nc)=codw0(1:nc)*lcodw%rbody(1:nc)
      elseif(ktfrealq(klist(kcodconv),w))then
        codw=codw0*w
      else
        codw=codw0
      endif
      return
      end subroutine

      recursive subroutine tcod(trans,cod,beam,optics,fndcod)
      use tfstk
      use ffs_flag
      use ffs_pointer
      use ffs_fit
      use ffs, only:mfitddp
      use tmacro
      use temw, only:iaez
      use sad_basics
      use mathfun,only:asinz
      implicit none
      integer*4 im
      logical*4 ,intent(out):: fndcod
      logical*4 ,intent(in):: optics
      integer*4 , parameter :: lmax=100
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42)
      real*8 codi(6),codf(6),dcod(6),r0,fact,trf00,dtrf0,r,
     $     dcod1(6),codw(6),conv,trs(6,6),dcod0(6),s,trw,dz,alambdarf,trf0s,v0,red,ddp,srot(3,9)
      integer*4 loop,i
      logical*4 rt,rtr
      vcalpha=1.d0
      trf0=0.d0
      rtr=radcod .and. radtaper
      if(rfsw)then
        rfsw=.false.
        if(rtr)then
          radcod=.false.
          radtaper=.false.
        endif
        call tcod(trans,cod,beam,optics,fndcod)
        rfsw=.true.
        if(rtr)then
          radcod=.true.
          radtaper=.true.
        endif
        if(fndcod .and. vceff /= 0.d0)then
          cod(5)=asinz((u0*pgev-vcacc)/vceff)/wrfeff-trf0
        endif
c        write(*,'(a,l2,1p7g15.7)')'tcod-0 ',
c     $       fndcod,vceff,wrfeff,trf0,vcacc,u0*pgev,cod(5)
        im=6
      else
        im=5
      endif
      trf0s=trf0
      v0=p0/h0
      rt=.false.
      dptaper=0.d0
 10   dcod=0.d0
      dcod0=0.d0
      r0=rmax
      fact=1.d0
      loop=lmax
      codi=cod
      fndcod=.true.
      call tgetcodw(codw)
      if(.not. rfsw)then
        codw=codw*0.01d0
      endif
      conv=conv0
      trw=codw(5)
      if(radcod .and. radtaper)then
        if(.not. rtr)then
          conv=conv*1.d2
        endif
        trw=codw(5)
        codi(6)=dp0
      endif
 1    loop=loop-1
      if(loop <= 0)then
        if(rtr)then
          rtr=.false.
          go to 10
        else
          fndcod=.false.
        endif
        return
      endif
      call tinitr12(trans)
      codf=codi
      trf00=trf0
      call tturne(trans,codf,beam,srot,iaez,.false.,.true.,rt,optics)
      dz=(codi(5)+codf(5))*0.5d0
      rt=radtaper
      dcod1=codi-codf
      if(rfsw)then
        if(radtaper)then
          if(.not. ktfenanq(dcod1(5)) .and. trans(5,6) /= 0.d0)then
            dptaper=dptaper-dcod1(5)/trans(5,6)
          endif
        endif
      else
        dcod1(5)=0.d0
      endif
      r=sum((dcod1(1:im)/codw(1:im))**2)
      dtrf0=trf0-trf00
      if(.not. radcod)then
        r=r+(dtrf0/trw)**2
      endif
c      write(6,'(a,2l3,1p10g12.5)')' tcod ',rtr,rfsw,r,r0,fact,trf0,dtrf0,dleng
c      write(6,'(1p6g12.5)')codi,codf,dcod1
      if(r < conv)then
        cod=codi
        return
      endif
      if(ktfenanq(r))then
        trf0=trf0s
        loop=0
        go to 1
      endif
      if(r >= r0)then
        fact=fact*.5d0
        codi=cod+fact*dcod
        if(fact < fmin)then
          cod=codi
          loop=0
        endif
        go to 1
      endif
      red=r/r0
      trf0s=trf0
      r0=r
      s=sum(dcod*dcod0/codw**2)
      if(red < 1.0d0)then
        fact=min(fact*(2.d0-max(red-0.7d0,0.d0)/0.3d0),1.d0)
      elseif(s < 0.d0)then
        fact=max(fact*0.5d0,fmin)
      endif
      trs(:,1:6)=trans(:,1:6)
      do i=1,6
        trs(i,i)=trs(i,i)-1.d0
      enddo
      if(.not. rfsw)then
        trs(1:6,5)=0.d0
        trs(1:6,6)=0.d0
        dcod1(6)=0.d0
        dcod1(5)=0.d0
      elseif(radtaper)then
        trs(1:6,6)=0.d0
c        dcod1(5)=0.d0
      endif
      dcod0=dcod
      ddp=dcod1(6)
      call tsolvg(trs,dcod1,dcod,im,6,6)
      if(radcod)then
        codi(5)=codi(5)-dz-dtrf0*v0
        alambdarf=pi2/wrfeff
        trf0=trf0+dz/v0
        trf0=sign(1.d0,trf0)*mod(
     $       sign(1.d0,trf0)*(trf0-0.5d0*alambdarf),alambdarf)
     $       +alambdarf*0.5d0
c        write(*,*)'tcod-dz ',dz,dvcacc,dz*dvcacc/pgev,ddp
        if(abs((dz-dtrf0*v0)*dvcacc/pgev)+abs(ddp) .gt. dpthre)then
          dcod(1:4)=0.d0
        endif
      endif
      cod=codi
      codi=codi+fact*dcod
      go to 1
      end subroutine

      end module
