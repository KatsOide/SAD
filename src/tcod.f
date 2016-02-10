      subroutine tcod(latt,trans,cod,beam,twiss,gammab,ndim,fndcod)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 im,lmax,lthre,ndim
      real*8 codmin,conv,epsr0,epsrr,rmax,fmin,a
      logical fndcod
      parameter (lmax=200,lthre=190,
     $     codmin=1.d-4,conv=1.d-10 ,epsr0=1.d-6,
     1     epsrr=1.d-4,rmax=1.d200,fmin=3.d-5,a=0.25d0)
      real*8 trans(6,12),cod(6),codi(6),codf(6),dcod(6),beam(42),
     1     gammab(*),r0,fact,trf00,dtrf0,r,
     $     epsr2,twiss(nlat,-ndim:ndim,*)
      integer*4 latt(2,nlat),loop,i
      logical*4 isnan,rt
      if(rfsw)then
        im=6
      else
        im=5
      endif
      trf0=0.d0
      vcphic=0.d0
      vcalpha=1.d0
      epsrad=epsr0*1000.d0
      r0=rmax
      fact=1.d0
      loop=lmax
      codi=cod
      fndcod=.true.
      rt=.false.
1     loop=loop-1
      if(loop .le. 0)then
        write(*,*)'Closed orbit was not found.'
        fndcod=.false.
        return
      endif
      codf=codi
      trans=0.d0
      call tinitr(trans)
      trf00=trf0
      call tturne(latt,trans,codf,beam,twiss,0.d0,gammab,
     $     int8(0),int8(0),int8(0),1,.false.,.true.,rt)
      rt=.true.
      r=0.d0
      do i=1,im
        r=r+((codi(i)-codf(i))/codmin)**2
      enddo
      if(isnan(r))then
        trf0=trf0*a+trf00*(1.d0-a)
      endif
      if(radcod .and. rfsw)then
        dtrf0=trf0-trf00
        cod(5)=cod(5)-dtrf0
        codi(5)=codi(5)-dtrf0
        codf(5)=codf(5)-dtrf0
      endif
c      write(6,'(a,1p5g12.5)')' tcod ',r,r0,epsrad,fact,trf0
c      write(6,'(1p6g12.5)')codi,codf
      if(r .lt. conv)then
        if(epsrad .gt. epsr0)then
          if(radcod .or. calpol)then
            epsr2=min(epsrad,max(epsr0,min(1.d-4,epsrr*sqrt(r))))
            r=(epsrad/epsrr)**2
            r0=r
            go to 22
          endif
        endif
        epsrad=epsr0
        cod=codi
c        if(rfsw)then
c          trf0=trf0+cod(5)
c          cod(5)=0.d0
c        endif
        return
      endif
      if(r .ge. r0 .or. isnan(r))then
        trf0=trf00
        fact=fact*.5d0
        codi=cod+fact*dcod
        if(isnan(r))then
          go to 1
        endif
        if(fact .lt. fmin)then
          if(radcod)then
            epsrad=max(epsr0,epsrad*.5d0)
          endif
          fact=fact*2.d0
          r=r0*2.d0
          go to 21
        endif
        go to 1
      endif
 21   epsr2=min(epsrad,max(epsr0,min(1.d-4,epsrr*sqrt(r))))
 22   if(epsr2 .lt. epsrad)then
        r0=r*(epsrad/epsr2)**2
        epsrad=epsr2
      endif
      fact=min(fact*2.d0,1.d0)
      cod=codi-codf
      do i=1,6
        trans(i,i)=trans(i,i)-1.d0
      enddo
      if( .not. rfsw)then
        trans(1:6,5)=0.d0
        trans(6,6)=0.d0
        cod(6)=0.d0
c      elseif(radcod .and. radtaper .and. loop .lt. lthre)then
c        dvfs=dvfs-cod(5)/(c*p0/h0/omega0*pi2)
c        call RsetGL1('FSHIFT',dvfs)
c        cod(5)=0
      endif
      call tsolvg(trans,cod,dcod,im,6,6)
c      write(*,'(a,1p7g14.6)')'tcod-dcod ',fact,dcod
      cod=codi
      codi=codi+fact*dcod
      go to 1
      end
