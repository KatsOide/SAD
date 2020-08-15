      real*8 function bint(f,a,b,epslon,epsabs)
      implicit none
      integer itmax,n,iter,i
      parameter (itmax=30)
      real*8 ,intent(in):: a,b,epslon,epsabs
      real*8 ba,s0,s20,s,xstep,x,s1,s2
      real*8 f
      external f

      ba=b-a
      s0=ba*(f(a)+f(b))*.5d0
c      write(*,*)'bint ',a,b,f(a),f(b)
      s20=s0
      n=1
      iter=0
1     s=0.d0
      xstep=ba/n
      x=a-xstep*.5d0
      do 10 i=1,n
        x=x+xstep
        s=s+f(x)
10    continue
      s=s*xstep
      s1=(s+s0)*.5d0
      s2=(2.d0*s+s0)/3.d0
c     write(*,*)s2
      if(abs(s2-s20) .lt. max(epslon*abs(s2),epsabs)
     1   .or. iter .gt. itmax)then
        bint=s2
        if(iter .gt. itmax)then
          write(*,*)'BINT Convergence failed. ',s2,s20
        endif
        return
      endif
      n=n*2
      s0=s1
      s20=s2
      iter=iter+1
      go to 1
      end

      real*8 function rombint(f,a,b,epslon,epsabs)
      implicit none
      integer*4 itmax,n,iter,i,k
      parameter (itmax=30,k=5)
      real*8 a,b,epslon,epsabs,ba,s0,s20,s,xstep,x,s1,s2,
     $     h(itmax),y(itmax),ds2
      real*8 f
      external f
      ba=b-a
      s0=ba*(f(a)+f(b))*.5d0
      s20=s0
      h(1)=1.d0
      y(1)=s0
      n=1
      iter=2
1     s=0.d0
      xstep=ba/n
      x=a-xstep*.5d0
      do i=1,n
        x=x+xstep
        s=s+f(x)
      enddo
      s=s*xstep
      s1=(s+s0)*.5d0
      h(iter)=h(iter-1)*.25d0
      y(iter)=s1
      if(iter .lt. k)then
        s2=(2.d0*s+s0)/3.d0
        if(abs(s2-s20) .lt. max(epslon*abs(s2),epsabs))then
          rombint=s2
          return
        endif
        s20=s2
      else
        call polint(h(iter-k+1),y(iter-k+1),k,0.d0,s2,ds2)
c        write(*,*)iter,s2,ds2
        rombint=s2
        if(abs(ds2) .lt. max(epslon*abs(s2),epsabs))then
          return
        endif
        if(iter .ge. itmax)then
          write(*,*)'ROMBINT Convergence failed. '
          return
        endif
      endif
      s0=s1
      n=n*2
      iter=iter+1
      go to 1
      end

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      implicit none
      INTEGER*4 n,NMAX
      PARAMETER (NMAX=10) 
      REAL*8 dy,x,y,xa(n),ya(n)
c Given arrays xa and ya, each of length n, and given a value x, this routine returns a
c value y, and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that
c P(xai) = yai, i = 1, . . . , n, then the returned value y = P(x).
      INTEGER*4 i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n 
        dift=abs(x-xa(i))
        if(dift .lt. dif)then
          ns=i
          dif=dift
        endif
        c(i)=ya(i) 
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den .eq. 0.d0)then
            write(*,*)'polint failure in polint.'
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if(2*ns .lt. n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
      return
      END
