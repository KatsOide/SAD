      subroutine eigs33(a,r,eig)
      implicit none
      real*8 ,intent(inout):: a(3,3)
      real*8 ,intent(out):: r(3,3),eig(3)

c     Use extended double precision real for internal variables if possible
      integer ex_real_kind
      parameter (ex_real_kind=max(
     $     selected_real_kind(precision(1.d0) * 2),
     $     selected_real_kind(precision(1.d0) + 1),
     $     selected_real_kind(precision(1.d0))))

      integer*4 ,parameter ::itmax=8,itmax1=3

      integer*4 it,i2,i1,i3
c      real(ex_real_kind) one,d,c,s,x21,x31,x22,x33,x32,
      real*8 one,d,c,s,x21,x31,x22,x33,x32,
     $     r11,r21,r31,r12,r22,r32,r13,r23,r33,rr,
     $     fact,dd,u,g,de,dde,dx22,dx32,dx21,tr
      parameter (one=1.)

      real*8 a1(3,3),a2(3,3),w(3,3)
      complex*16 ceig(3)
      integer*4 ic
      
      a(1,2)=a(2,1)
      a(1,3)=a(3,1)
      a(2,3)=a(3,2)
      if(abs(a(3,3)) >= max(abs(a(1,1)),abs(a(2,2))))then
        i3=3
        i2=merge(1,2,abs(a(2,2)) <= abs(a(1,1)))
      elseif(abs(a(1,1)) >= max(abs(a(2,2)),abs(a(3,3))))then
        i3=1
        i2=merge(3,2,abs(a(2,2)) <= abs(a(3,3)))
      else
        i3=2
        i2=merge(3,1,abs(a(1,1)) <= abs(a(3,3)))
      endif
      i1=6-i2-i3

      d=hypot(a(i2,i1),a(i3,i1))
      if(d /= 0.d0)then
        c=a(i2,i1)/d
        s=a(i3,i1)/d
      else
        c=1.d0
        s=0.d0
      endif
      de=a(i1,i1)
      x21=d
      x31=0.d0
      dx22=s*((a(i3,i3)-a(i2,i2))*s+2.d0*a(i3,i2)*c)
      dx32=s*((a(i3,i3)-a(i2,i2))*c-2.d0*a(i3,i2)*s)
      x22=a(i2,i2)+dx22-de
      x33=a(i3,i3)-dx22-de
      x32=a(i3,i2)+dx32
c      x22=a(i2,i2)*c**2+a(i3,i3)*s**2+2.d0*a(i3,i2)*c*s
c      x33=a(i2,i2)*s**2+a(i3,i3)*c**2-2.d0*a(i3,i2)*c*s-de
c      x32=a(i3,i2)*(c-s)*(c+s)+(a(i3,i3)-a(i2,i2))*c*s
      r11=1.d0
      r21=0.d0
      r31=0.d0
      r12=0.d0
      r22=c
      r32=s
      r13=0.d0
      r23=-s
      r33=c
      fact=1.d0
      do it=1,itmax
        dd=abs(x22)
        if(abs(x21)/fact+dd == dd)then
          eig(i1)=0.d0
          u=x22-x33
          d=hypot(u,2.d0*x32)
          if(d == 0.d0)then
            eig(i2)=dble(x22)
            eig(i3)=dble(x33)
          else
            s=-sign(sqrt(.5d0*(d+abs(u))/d),u)
            c=abs(x32/s)/d
c            c=sqrt(.5d0*(1.d0+abs(u)/d))
c            s=sign(sqrt(2.d0*x32**2/d/(d+abs(u))),u*x32)
c     eig(2) * eig(3) == x22 * x33 - x32^2
c     eig(2) + eig(3) == x22 + x33
c     d := Sqrt[(x22 - x33)^2 + (2 * x32)^2]
            tr=x22+x33
            eig(i2)=dble(.5d0*(tr+sign(d,tr)))
            eig(i3)=dble(.5d0*(tr-sign(d,tr)))
            rr=r12
            r12= c*rr+s*r13
            r13=-s*rr+c*r13
            rr=r22
            r22= c*rr+s*r23
            r23=-s*rr+c*r23
            rr=r32
            r32= c*rr+s*r33
            r33=-s*rr+c*r33
          endif
c          if(eig(1)*eig(2)*eig(3) <= 0.d0)then
c            write(*,'(a,3i5,1p10g12.4)')'eigs33-1 ',i1,i2,i3,eig,x22+x33,d,u
c            stop
c          endif
          ic=1
          exit
        endif
        dd=sqrt(abs(x22*x33))
        if(abs(x32)/fact+dd == dd)then
          eig(i3)=dble(x33)
          d=hypot(x22,2.d0*x21)
          if(d == 0.d0)then
            eig(i1)=0.d0
            eig(i2)=0.d0
          else
            s=sign(sqrt(.5d0*(d+abs(x22))/d),x22)
            c=abs(x21/s)/d
c            s=sign(sqrt(2.d0*x21**2/d/(d+abs(u))),u*x21)
c     eig(1) * eig(2) == x11 * x22 - x21^2
c     eig(1) + eig(2) == x11 + x22
c     d := Sqrt[(x11 - x22)^2 + (2 * x21)^2]
c            eig(i1)=dble(.5d0*(x22-sign(d,x22)))
            eig(i1)=dble(.5d0*(x22+sign(d,x22)))
            eig(i2)=dble(-x21*(x21/eig(i1)))
            rr=r11
            r11= c*rr+s*r12
            r12=-s*rr+c*r12
            rr=r21
            r21= c*rr+s*r22
            r22=-s*rr+c*r22
            rr=r31
            r31= c*rr+s*r32
            r32=-s*rr+c*r32
          endif
          ic=2
          exit
        endif
        d=hypot(x22,2.d0*x21)
        g=2.d0*x21*(x21/(x22+sign(d,x22)))
        d=hypot(x32,x33-g)
        c=(x33-g)/d
        s=-x32/d
        dx22=((x33-x22)*s+2.d0*x32*c)*s
        dx32=((x33-x22)*c-2.d0*x32*s)*s
c        x22=y22*c**2+x33*s**2+2.d0*x32*c*s
c        x33=y22*s**2+x33*c**2-2.d0*x32*c*s
c        x32=x32*(c-s)*(c+s)-(y22-x33)*c*s
        x22=x22+dx22
        x33=x33-dx22
        x32=x32+dx32
        x31=-x21*s
        x21= x21*c
        rr=r12
        r12= c*rr+s*r13
        r13=-s*rr+c*r13
        rr=r22
        r22= c*rr+s*r23
        r23=-s*rr+c*r23
        rr=r32
        r32= c*rr+s*r33
        r33=-s*rr+c*r33
        if(abs(x31)/fact+abs(x32) /= abs(x32))then
          d=hypot(x32,x31)
          c=x32/d
          s=-x31/d
          dde =(x22*s+2.d0*x21*c)*s
          dx21=(x22*c-2.d0*x21*s)*s
          de=de+dde
          x33=x33-dde
          x22=x22-2.d0*dde
          x21=x21+dx21
c          x22=y22*c**2-2.d0*x21*c*s-dde
c          x21=x21*(c-s)*(c+s)+y22*c*s
          x32=d
          x31=0.d0
          rr=r11
          r11= c*rr+s*r12
          r12=-s*rr+c*r12
          rr=r21
          r21= c*rr+s*r22
          r22=-s*rr+c*r22
          rr=r31
          r31= c*rr+s*r32
          r32=-s*rr+c*r32
        endif
        if(it >= itmax1)then
          fact=fact*2.d0
        endif
      enddo
      r(i1,i1)=dble(r11)
      r(i2,i1)=dble(r21)
      r(i3,i1)=dble(r31)
      r(i1,i2)=dble(r12)
      r(i2,i2)=dble(r22)
      r(i3,i2)=dble(r32)
      r(i1,i3)=dble(r13)
      r(i2,i3)=dble(r23)
      r(i3,i3)=dble(r33)
      eig=eig+dble(de)
      a1=matmul(transpose(r),matmul(a,r))
      if(abs(a1(2,1)) > 1e-5*sqrt(abs(a1(1,1)*a1(2,2))))then
        write(*,'(a,4i3,1p9g12.4)')'eigs33 ',it,i1,i2,i3,a1
        write(*,'(a,i3,1p9g12.4)')'                ',ic,a
        write(*,'(a,1p9g12.4)')'                   ',r
        write(*,'(a,1p9g12.4)')'                   ',eig,de,dde,.5d0*x22,0.5d0*sign(d,x22)
        a2=a
        call teigen(a2,w,ceig,3,3)
        write(*,'(a,1p9g12.4)')'                   ',a2
        write(*,'(a,1p9g12.4)')'                   ',ceig
        call refs33(a,r,eig)
        stop
      endif
      return
      end

      subroutine refs33(a,r,eig)
      implicit none
      real*8 ,intent(in):: a(3,3),r(3,3),eig(3)
      integer*4 i1,i2,i3
      real*8 rr(3,3),aa(3,3),da(3,3)
      if(abs(eig(3)) >= max(abs(eig(1)),abs(eig(2))))then
        i1=3
        i2=merge(1,2,abs(eig(1)) >= abs(eig(2)))
      elseif(abs(eig(2)) >= max(abs(eig(3)),abs(eig(1))))then
        i1=2
        i2=merge(3,1,abs(eig(3)) >= abs(eig(1)))
      else
        i1=1
        i2=merge(2,3,abs(eig(2)) >= abs(eig(3)))
      endif
      i3=6-i1-i2
      rr(:,1)=r(:,i1)/norm2(r(:,i1))
      rr(:,2)=r(:,i2)-dot_product(rr(:,1),r(:,i2))*rr(:,1)
      rr(:,2)=rr(:,2)/norm2(rr(:,2))
      rr(:,3)=r(:,i3)-dot_product(rr(:,1),r(:,i3))*rr(:,1)-dot_product(rr(:,2),r(:,i3))*rr(:,2)
      rr(:,3)=rr(:,3)/norm2(rr(:,3))
      aa(1,:)=rr(:,1)*eig(i1)
      aa(2,:)=rr(:,2)*eig(i2)
      aa(3,:)=rr(:,3)*eig(i3)
      aa=matmul(rr,aa)
      da=a-aa
      WRITE(*,'(a,1p10g12.4)')'refs33 ',aa
      WRITE(*,'(a,1p10g12.4)')'       ',a
      WRITE(*,'(a,1p10g12.4)')'       ',da
      return
      end

      real*8 function eigr33(r,u,ndim)
      use macmath
      implicit none
      integer*4 ,intent(in):: ndim
      real*8 ,intent(in):: r(ndim,3)
      real*8 ,intent(out):: u(3)
      real*8 a,b,s2,u1,u2,u3,w,x1,x2,x3,c,s
      s2=r(1,3)**2+r(2,3)**2
      if(s2 == 0.d0)then
        if(r(3,3) > 0.d0)then
          u(1)=0.d0
          u(2)=0.d0
          u(3)=1.d0
          eigr33=atan2(-r(1,2),r(1,1))
        else
          u1=r(1,2)
          u2=1.d0-r(1,1)
          w=hypot(u1,u2)
          u(1)=u1/w
          u(2)=u2/w
          u(3)=0.d0
          eigr33=pi
        endif
      else
        b=r(3,2)*r(1,3)-r(3,1)*r(2,3)
        a=-(1.d0+r(1,1))*r(2,3)**2+(r(1,2)+r(2,1))*r(1,3)*r(2,3)
     1    -(1.d0+r(2,2))*r(1,3)**2
        u1=b*r(1,3)+a*r(2,3)
        u2=b*r(2,3)-a*r(1,3)
        u3=b*(1.d0+r(3,3))
c        w=sign(hypot(u1,hypot(u2,u3)),u3*r(3,3))
c        w=sign(hypot3(u1,u2,u3),u3*r(3,3))
        w=sign(norm2([u1,u2,u3]),u3*r(3,3))
        u1=u1/w
        u2=u2/w
        u3=u3/w
        x1=r(1,1)*u2-r(1,2)*u1
        x2=r(2,1)*u2-r(2,2)*u1
        x3=r(3,1)*u2-r(3,2)*u1
        c=x1*u2-x2*u1
        s=(x1*u1+x2*u2)*u3-x3*(u1**2+u2**2)
        eigr33=atan2(s,c)
        u(1)=u1
        u(2)=u2
        u(3)=u3
      endif
      return
      end

      subroutine rot33(r,u,phi,ndim)
      implicit none
      integer*4 ,intent(in):: ndim
      real*8 ,intent(out):: r(ndim,3)
      real*8 ,intent(in):: u(3)
      real*8 phi,w,u1,u2,u3,c,s,q,th
c      w=hypot(u(1),hypot(u(2),u(3)))
      w=norm2([u(1),u(2),u(3)])
      u1=u(1)/w
      u2=u(2)/w
      u3=u(3)/w
      th=tan(.5d0*phi)
      s=2.d0*th/(1.d0+th**2)
      q=th*s
      c=1.d0-q
c      q=2.d0*sin(.5d0*phi)**2
      r(1,1)=u1**2+(u2**2+u3**2)*c
      r(2,1)=u1*u2*q+u3*s
      r(3,1)=u1*u3*q-u2*s
      r(2,2)=u2**2+(u3**2+u1**2)*c
      r(3,2)=u2*u3*q+u1*s
      r(1,2)=u2*u1*q-u3*s
      r(3,3)=u3**2+(u1**2+u2**2)*c
      r(1,3)=u3*u1*q+u2*s
      r(2,3)=u3*u2*q-u1*s
      return
      end
