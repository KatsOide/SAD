      subroutine spline(np,x,y,ddy,work,mode,dy)
C
C     Normal Spline y''(1)=y''(np)=0.
C     ddy returns y''(i)/6.
C     y=A y(i)+B y(i+1)+C ddy(i)+D ddy(i+1)
C     y'=(y(i+1)-y(i))/(x(i+1)-x(i))
C        -(3A^2-1)(x(i+1)-x(i))ddy(i)
C        +(3B^2-1)(x(i+1)-x(i))ddy(i+1)
C     A=(x(i+1)-x)/(x(i+1)-x(i))
C     B=1-A=(x-x(i))/(x(i+1)-x(i))
C     C=(A^3-A)(x(i+1)-x(i))^2
C     D=(B^3-B)(x(i+1)-x(i))^2
C
      implicit none
      integer*4 np,i,mode
      real*8 ,intent(in):: x(np),y(np),dy(2)
      real*8 ,intent(out):: ddy(np),work(np)
      real*8 dy1a,ddy1
      real*8, allocatable :: dddy(:)
      real*8 f
      if(np .le. 2)then
        ddy=0.d0
        return
      endif
      allocate (dddy(np))
      if(mode .eq. 0 .or. mode .eq. 2)then
        ddy(1)=0.d0
        ddy(2)=(y(3)-y(2))/(x(3)-x(2))-(y(2)-y(1))/(x(2)-x(1))
        work(2)=(x(3)-x(1))*2.d0
      elseif(mode .eq. 1 .or. mode .eq. 3)then
        work(1)=(x(2)-x(1))*.5d0
        work(2)=(x(3)-x(1))*2.d0-work(1)
        ddy(1)=.5d0*((y(2)-y(1))/work(1)/2.d0-dy(1))
        ddy(2)=(y(3)-y(2))/(x(3)-x(2))-(y(2)-y(1))/(x(2)-x(1))-ddy(1)
      elseif(mode .eq. 4)then
        dy1a=(y(2)-y(1)+y(np)-y(np-1))/(x(2)-x(1)+x(np)-x(np-1))
        work(1)=(x(2)-x(1))*.5d0
        work(2)=(x(3)-x(1))*2.d0-work(1)
        ddy(1)=.5d0*((y(2)-y(1))/work(1)/2.d0-dy1a)
        ddy(2)=(y(3)-y(2))/(x(3)-x(2))-(y(2)-y(1))/(x(2)-x(1))-ddy(1)
        dddy(1)=-0.5d0
        dddy(2)= 0.5d0
        do i=2,np-2
          f=x(i+1)-x(i)
          work(i+1)=(x(i+2)-x(i))*2.d0-f**2/work(i)
          ddy(i+1) =(y(i+2)-y(i+1))/(x(i+2)-x(i+1))
     1         -(y(i+1)-y(i))/f-f/work(i)*ddy(i)
          dddy(i+1)=-f/work(i)*dddy(i)
        enddo
        work(np)=x(np)-x(np-1)
        ddy(np)=dy1a-(y(np)-y(np-1))/work(np)
        ddy(np-1)=(ddy(np-1)-ddy(np)*.5d0)/
     $       (work(np-1)-.5d0*work(np))
        ddy(np)=(ddy(np)/work(np)-ddy(np-1))*.5d0
        dddy(np-1)=(dddy(np-1)-.5d0)/
     $       (work(np-1)-.5d0*work(np))
        dddy(np)=(1.d0/work(np)-dddy(np-1))*.5d0
        do i=np-2,2,-1
          ddy(i) =( ddy(i)- ddy(i+1)*(x(i+1)-x(i)))/work(i)
          dddy(i)=(dddy(i)-dddy(i+1)*(x(i+1)-x(i)))/work(i)
        enddo
        ddy(1) =(ddy(1)/work(1)- ddy(2))*.5d0
        dddy(1)=(-0.5d0/work(1)-dddy(2))*.5d0
        if(dddy(1) .ne. dddy(np))then
          ddy1=-(ddy(np)-ddy(1))/(dddy(np)-dddy(1))
          ddy=ddy+ddy1*dddy
        endif
        deallocate (dddy)
        return
      endif
      do i=2,np-2
        f=x(i+1)-x(i)
c        write(*,'(a,i4,1p5g12.5)')'spline ',
c     $       i,f,work(i),x(i),x(i+1),x(i+2)
        work(i+1)=(x(i+2)-x(i))*2.d0-f**2/work(i)
        ddy(i+1)=(y(i+2)-y(i+1))/(x(i+2)-x(i+1))
     1           -(y(i+1)-y(i))/f-f/work(i)*ddy(i)
      enddo
      if(mode .eq. 0 .or. mode .eq. 1)then
        ddy(np)=0.d0
        ddy(np-1)=ddy(np-1)/work(np-1)
      elseif(mode .eq. 2 .or. mode .eq. 3)then
        work(np)=x(np)-x(np-1)
        ddy(np)=dy(2)-(y(np)-y(np-1))/work(np)
        ddy(np-1)=(ddy(np-1)-ddy(np)*.5d0)/
     $       (work(np-1)-.5d0*work(np))
        ddy(np)=(ddy(np)/work(np)-ddy(np-1))*.5d0
      endif
      do i=np-2,2,-1
        ddy(i)=(ddy(i)-ddy(i+1)*(x(i+1)-x(i)))/work(i)
c        write(*,*)'spline ',i,work(i),ddy(i)
      enddo
      if(mode .eq. 1 .or. mode .eq. 3)then
        ddy(1)=(ddy(1)/work(1)-ddy(2))*.5d0
      endif
      deallocate (dddy)
      return
      end

      subroutine spline1(np,y,ddy,work,mode1,mode2)
C
C     Spline when dx=const.=1
C     Normal Spline y''(1)=y''(np)=0.
C     ddy returns y''(i)/6.
C     y=A y(i)+B y(i+1)+C ddy(i)+D ddy(i+1)
C     y'=(y(i+1)-y(i))
C        -(3A^2-1)ddy(i)
C        +(3B^2-1)ddy(i+1)
C     A=(x(i+1)-x)
C     B=1-A=(x-x(i))
C     C=(A^3-A)
C     D=(B^3-B)
C
      implicit none
      integer*4 ,intent(in)::np,mode1,mode2
      integer*4 i
      real*8 ,intent(in)::y(np)
      real*8 ,intent(out)::ddy(np),work(np)
      if(mode1 .eq. 0)then
        ddy(1)=0.d0
        if(np .le. 2)then
          ddy(2)=0.d0
          return
        endif
        ddy(2)=y(3)-2.d0*y(2)+y(1)
        work(2)=4.d0
      elseif(mode1 .eq. 1)then
        work(1)=.5d0
        work(2)=3.5d0
        ddy(1)=.5d0*(y(2)-y(1)-ddy(1))
        ddy(2)=y(3)-2.d0*y(2)+y(1)-ddy(1)
      elseif(mode1 .eq. 2)then
        ddy(2)=y(3)-2.d0*y(2)+y(1)-ddy(1)
        work(2)=4.d0
      endif
      do i=2,np-2
        work(i+1)=4.d0-1.d0/work(i)
        ddy(i+1)=y(i+2)-2.d0*y(i+1)+y(i)
     1           -ddy(i)/work(i)
      enddo
      if(mode2 .eq. 0)then
        ddy(np)=0.d0
        ddy(np-1)=ddy(np-1)/work(np-1)
      elseif(mode2 .eq. 1)then
        ddy(np)=ddy(np)-(y(np)-y(np-1))
        ddy(np-1)=(ddy(np-1)-ddy(np)*.5d0)/
     $       (work(np-1)-.5d0)
        ddy(np)=(ddy(np)-ddy(np-1))*.5d0
      elseif(mode2 .eq. 2)then
        ddy(np-1)=(ddy(np-1)-ddy(np))/work(np-1)
      endif
      do 20 i=np-2,2,-1
        ddy(i)=(ddy(i)-ddy(i+1))/work(i)
20    continue
      if(mode1 .eq. 1)then
        ddy(1)=ddy(1)-ddy(2)*.5d0
      endif
      return
      end

      real*8 function splint1(np,y,mode1,mode2,dy)
      implicit none
      integer*4 , intent(in)::np,mode1,mode2
c      integer*4 i
      real*8 ,intent(in)::y(np),dy(2)
c      real*8 s
      real*8,allocatable:: work(:),ddy(:)
      allocate (work(np),ddy(np))
      ddy(1)=dy(1)
      ddy(np)=dy(2)
      call spline1(np,y,ddy,work,mode1,mode2)
      splint1=(y(1)+y(np)-.5d0*(ddy(1)+ddy(np)))*.5d0
     $     +sum(y(2:np-1)-ddy(2:np-1)*.5d0)
      deallocate (work,ddy)
      return
      end

      real*8 function splint(f,x0,x1,mode,dy,eps,epsabs,n0)
      implicit none
      integer*4 nmax,mode
      parameter (nmax=8192)
      integer*4 n,k,k1,i,i1
      integer*4 , intent(in)::n0
      real*8 , external:: f
      real*8 , intent(in)::x0,x1,eps,epsabs,dy(2)
      real*8 , allocatable::x(:,:),y(:,:),ddy(:),work(:)
      real*8 xi,xstep,sddy,s,s0,dx,dx2,s20,s2,ddyi
      logical*4 first
      allocate (x(nmax,2),y(nmax,2),ddy(nmax),work(nmax))
      first=.true.
      n=max(n0,4)
      k=1
      xi=x0
      xstep=(x1-x0)/(n-1)
      s0=0.d0
      do i=1,n
        x(i,1)=xi
        y(i,1)=f(xi)
        xi=xi+xstep
        s0=s0+y(i,1)
      enddo
      s0=(s0-.5d0*(y(1,1)+y(n,1)))*xstep
      s20=s0
 1    call spline(n,x(1,k),y(1,k),ddy,work,mode,dy)
      s=0.d0
      sddy=0.d0
      do i=1,n-1
        dx=x(i+1,k)-x(i,k)
        dx2=.5d0*dx**2
        s=s+(y(i+1,k)+y(i,k)-(ddy(i+1)+ddy(i))*dx2)*dx
        ddyi=merge(-(ddy(i)-ddy(i-1))/(x(i,k)-x(i-1,k)),0.d0,
     $       i .gt. 1)
        if(i .lt. n-1)then
          ddyi=ddyi+(ddy(i+2)-ddy(i+1))/(x(i+2,k)-x(i+1,k))
        endif
        work(i)=abs(ddyi)*dx2**2
        sddy=sddy+work(i)
      enddo
      s=s*.5d0
      s2=merge(s,(4.d0*s-s0)/3.d0,first)
c      write(*,*)n,s2,s
      if(abs(s2-s20) .le. max(eps*s20,epsabs))then
        splint=s2
        deallocate (x,y,ddy,work)
        return
      endif
      first=.false.
      sddy=sddy/(n-1)
      k1=3-k
      x(1,k1)=x(1,k)
      y(1,k1)=y(1,k)
      i1=1
      do i=1,n-1
        if(work(i) .ge. sddy)then
          i1=i1+1
          x(i1,k1)=(x(i,k)+x(i+1,k))*.5d0
          y(i1,k1)=f(x(i1,k1))
        endif
        i1=i1+1
        x(i1,k1)=x(i+1,k)
        y(i1,k1)=y(i+1,k)
        if(i1 .ge. nmax-2)then
          go to 9000
        endif
      enddo
      k=k1
      n=i1
      s0=s
      s20=s2
      go to 1
 9000 write(*,*)'Splint convergence failed. ',s,s0
      splint=s2
      deallocate (x,y,ddy,work)
      return
      end

      subroutine spline1m(m,md,np,y,ddy,work)
C
C     Spline when dx=const.=1, for vector.
C     Normal Spline y''(1)=y''(np)=0.
C     ddy returns y''(i)/6.
C     y=A y(i)+B y(i+1)+C ddy(i)+D ddy(i+1)
C     y'=(y(i+1)-y(i))
C        -(3A^2-1)ddy(i)
C        +(3B^2-1)ddy(i+1)
C     A=(x(i+1)-x)
C     B=1-A=(x-x(i))
C     C=(A^3-A)
C     D=(B^3-B)
C
      implicit none
      integer*4 np,i,m,j,md
      real*8 y(md,np),ddy(md,np),work(np)
      do j=1,m
        ddy(j,1)=0.d0
        ddy(j,2)=y(j,3)-2.d0*y(j,2)+y(j,1)
      enddo
      work(2)=4.d0
      do 10 i=2,np-2
        work(i+1)=4.d0-1.d0/work(i)
        do j=1,m
          ddy(j,i+1)=y(j,i+2)-2.d0*y(j,i+1)+y(j,i)
     1           -ddy(j,i)/work(i)
        enddo
10    continue
      do j=1,m
        ddy(j,np)=0.d0
        ddy(j,np-1)=ddy(j,np-1)/work(np-1)
      enddo
      do 20 i=np-2,2,-1
        do j=1,m
          ddy(j,i)=(ddy(j,i)-ddy(j,i+1))/work(i)
        enddo
20    continue
      return
      end

      subroutine tfspline(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kd,k1,k2
      type (sad_dlist), pointer :: klx,kld
      integer*8 ka1,ka2,ix,iy,iddy,ka
      integer*4 isp1,irtc,itfmessage,n,m,mode
      real*8 dy(2)
      type (sad_descriptor), save :: kxperiodic
      data kxperiodic%k /0/
      if(kxperiodic%k .eq. 0)then
        kxperiodic=kxsymbolz('`Periodic',9)
      endif
      mode=0
      if(isp .eq. isp1+2)then
        call tfgetoption('Derivative',ktastk(isp),kd,irtc)
        if(irtc .ne. 0)then
          if(irtc .eq. -1)then
            go to 9100
          endif
          return
        endif
        if(tfsameq(kd,kxperiodic))then
          mode=4
        else
          if(ktfnonlistq(kd,kld))then
            go to 9100
          endif
          if(kld%nl .ne. 2)then
            go to 9100
          endif
          k1=kld%dbody(1)
          if(ktfoperq(k1,ka1))then
            if(ka1 .ne. mtfnull)then
              go to 9100
            endif
          elseif(ktfnonrealq(k1,dy(1)))then
            go to 9100
          else
            mode=1
          endif
          k2=kld%dbody(2)
          if(ktfoperq(k2,ka2))then
            if(ka2 .ne. mtfnull)then
              go to 9100
            endif
          elseif(ktfnonrealq(k2,dy(2)))then
            go to 9100
          else
            mode=mode+2
          endif
        endif
      elseif(isp .ne. isp1+1)then
        go to 9100
      endif
      call tfmsize(ktastk(isp1+1),n,m,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(n .ne. 2)then
        go to 9000
      endif
      ka=ktfaddr(ktastk(isp1+1))
      ix=ktfaddr(klist(ka+1))
      iy=ktfaddr(klist(ka+2))
      iddy=ktavaloc(0,max(m,3))+1
      call spline2(m,
     $     rlist(ix+1:ix+m),
     $     rlist(iy+1:iy+m),
     $     rlist(iddy:iddy+m-1),mode,dy)
      kx=kxadaloc(-1,3,klx)
      klx%dbody(1)%k=ktflist+ktfcopy1(ix)
      klx%dbody(2)%k=ktflist+ktfcopy1(iy)
      klx%dbody(3)%k=ktflist+iddy-1
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $     '"Argument","{{x1, ..}, {y1, ..}}"')
      return
 9100 irtc=itfmessage(9,'General::narg','"1 + Derivative->{dy1,dy2}"')
      return
      end

      subroutine spline2(m,x,y,ddy,mode,dy)
      implicit none
      integer*4 m,mode
      real*8 x(m),y(m),ddy(m),work(m),dy(2)
      call spline(m,x,y,ddy,work,mode,dy)
      return
      end

      subroutine tffindindex(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer ::kl
      type (sad_rlist), pointer ::klr,kll
      integer*4 isp1,irtc,itfmessage,i1,i2,ih,n,i,m
      real*8 x
      if(isp1+2 .ne. isp)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. tflistq(dtastk(isp1+1),kl))then
        irtc=itfmessage(9,'General::wrongval','"#1","List"')
        return
      endif
      if(ktfnonreallistqo(kl))then
        irtc=itfmessage(9,'General::wrongval','"#1","List of Reals"')
        return
      endif
      if(ktfrealq(ktastk(isp)))then
        i1=1
        i2=kl%nl
        x=rtastk(isp)
        irtc=0
        do while(i2 .gt. i1+1)
          ih=i1+(i2-i1)/2
          if(kl%rbody(ih) .eq. x)then
            kx=dfromr(dble(ih))
            return
          elseif(kl%rbody(ih) .gt. x)then
            i2=ih
          else
            i1=ih
          endif
        enddo
        kx=dfromr(dble(i1))
      elseif(tfreallistq(ktastk(isp),kll))then
        n=kll%nl
        m=kl%nl
        kx=kxavaloc(-1,n,klr)
        i1=1
        do i=1,n
          i2=min(i1+1,m)
          x=kll%rbody(i)
          if(x .lt. kl%rbody(i1))then
            i2=i1
            i1=1
          elseif(x .eq. kl%rbody(i2))then
            i1=i2
          elseif(x .gt. kl%rbody(i2))then
            i1=i2
            i2=m
          endif
          do while(i2 .gt. i1+1)
            ih=i1+(i2-i1)/2
            if(kl%rbody(ih) .eq. x)then
              i1=ih
              exit
            elseif(kl%rbody(ih) .gt. x)then
              i2=ih
            else
              i1=ih
            endif
          enddo
          klr%rbody(i)=min(i1,m-1)
        enddo
        irtc=0
      else
        irtc=itfmessage(9,'General::wrongval',
     $       '"#2","Real or list of Reals"')
      endif
      return
      end
