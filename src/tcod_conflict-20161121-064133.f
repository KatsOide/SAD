      recursive subroutine tcod(trans,cod,beam,fndcod)
      use tfstk
      use ffs_flag
      use ffs_pointer
      use ffs_fit
      use tmacro
      implicit none
      integer*4 im,lmax
      real*8 conv0,epsr0,epsrr,rmax,fmin,a,ddpmax,red
      logical fndcod
      parameter (lmax=300,
     $     conv0=1.d-10,epsr0=1.d-6,ddpmax=3.e-5,
     1     epsrr=1.d-4,rmax=1.d200,fmin=1.d-4,a=0.25d0)
      real*8 trans(6,12),cod(6),codi(6),codf(6),dcod(6),beam(42),
     1     r0,fact,trf00,dtrf0,r,dcodw(6),codw(6),conv,trs(6,6),
     $     dcod0(6),s
      integer*4 loop,i
      logical*4 isnan,rt
      real*8 , parameter :: codw0(6) =
     $     [1.d-6,1.d-5,1.d-6,1.d-5,1.d-5,1.d-6]
      if(rfsw)then
        im=6
      else
        im=5
      endif
      dcod=0.d0
      dcod0=0.d0
      trf0=0.d0
      vcphic=0.d0
      vcalpha=1.d0
      epsrad=epsr0
      r0=rmax
      fact=1.d0
      loop=lmax
      codi=cod
      fndcod=.true.
      rt=.false.
      codw=codw0
      conv=conv0
      if(radtaper)then
        if(codplt)then
          codw(5)=codw(5)*100.d0
          conv=conv*1.d4
        else
          codw(5)=codw(5)*1000.d0
          conv=conv*1.d5
        endif
      endif
 1    loop=loop-1
      if(loop .le. 0)then
        if(radtaper .and. codplt)then
          codplt=.false.
          call tcod(trans,cod,beam,fndcod)
          codplt=.true.
        else
          fndcod=.false.
        endif
        return
      endif
      call tinitr(trans)
      if(radtaper .and. codplt)then
        cod(6)=0.d0
        codi(6)=0.d0
      endif
      codf=codi
      trf00=trf0
      call tturne(trans,codf,beam,
     $     int8(0),int8(0),int8(0),.false.,.true.,rt)
      rt=.true.
      dcodw=(codi-codf)/codw
      r=0.d0
      do i=1,im
        r=r+dcodw(i)**2
      enddo
      if(isnan(r))then
        trf0=trf0*a+trf00*(1.d0-a)
      endif
      dtrf0=trf0-trf00
c      write(6,'(a,1p5g12.5)')' tcod ',r,r0,fact,trf0,trf00
c      write(6,'(1p6g12.5)')codi,codf,dcod
      if(r .lt. conv)then
c        trf0=trf0+codi(5)
c        cod=codi
c        cod(5)=0.d0
        return
      endif
      if(r .ge. r0 .or. isnan(r))then
        trf0=trf00
        fact=fact*.5d0
        if(radcod .and. rfsw)then
          cod(5)=cod(5)+dtrf0
        endif
        codi=cod+fact*dcod
        if(isnan(r))then
          go to 1
        endif
        if(fact .lt. fmin)then
          loop=0
        endif
        go to 1
      endif
      r0=r
      s=0.d0
      do i=1,6
        s=s+dcod(i)*dcod0(i)/codw(i)**2
      enddo
      red=r/r0
      if(red .lt. 0.7d0 .or. s .ge. 0.d0)then
        fact=min(fact*2.d0,1.d0)
      else
        fact=max(fact*0.5d0,fmin)
      endif
      do i=1,6
        trs(:,i)=trans(:,i)/codw
      enddo
      do i=1,6
        trs(i,1:6)=trs(i,1:6)*codw
      enddo
      do i=1,6
        trs(i,i)=trs(i,i)-1.d0
      enddo
      if( .not. rfsw)then
        trs(1:6,5)=0.d0
        trs(6,6)=0.d0
        cod(6)=0.d0
      elseif(radtaper .and. codplt)then
        trs(1:6,6)=0.d0
      endif
      dcod0=dcod
      call tsolvg(trs,dcodw,dcod,im,6,6)
      dcod=dcod*codw
      if(radtaper)then
        dcod(6)=min(ddpmax,max(-ddpmax,dcod(6)))
      endif
      cod=codi
      codi=codi+fact*dcod
      go to 1
      end
