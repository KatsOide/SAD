      subroutine tquad(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak0,bz,
     1                 dx,dy,theta,theta2,krad,chro,
     1                 fringe,f1in,f2in,f1out,f2out,mfring,eps0,kin)
      use ffs_flag
      use tmacro
c      use ffs_pointer, only:inext,iprev
      use photontable,only:tsetpcvt,pcvt
      use sol,only:tsolrot
      use mathfun, only:pxy2dpz,sqrt1,akang
      use kradlib
      implicit none
      integer*4 ,intent(in):: np,mfring
      logical*4 , intent(in)::krad,chro,fringe,kin
      integer*4 i
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),g(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak0,dx,dy,theta,theta2,
     $     f1in,f1out,f2in,f2out,eps0,bz
      real*8 ak,alr,p,a,ea,b,pxf,pyf,bxs,bys,bzs
      real*8, parameter :: ampmax=0.9999d0
      if(al .eq. 0.d0)then
        call tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,4,0.d0,ak0,
     $       dx,dy,theta,.false.,.false.)
        return
      elseif(ak0 .eq. 0.d0)then
        if(bz .eq. 0.d0)then
          call tdrift_free(np,x,px,y,py,z,dv,al)
        else
          call tdrift_solenoid(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,bz,krad)
        endif
        return
      endif
c      theta2=theta+akang(dcmplx(ak0,0.d0),al,cr1)
      ak=sign(ak0,al)
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.true.)
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
        if(calpol)then
          bsi=0.d0
        endif
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 2)then
        call ttfrin(np,x,px,y,py,z,g,4,ak,al,bzs)
      endif
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        do concurrent (i=1:np)
          p=1.d0+g(i)
          a=f1in/p
          ea=exp(a)
          b=f2in/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
        enddo
      endif
      if(krad)then
        if(photons)then
          call tsetpcvt(l_track,dx,dy,theta2,0.d0,0.d0,al)
          pcvt%fr0=-0.5d0*f1in/al
        endif
        if(f1in .ne. 0.d0)then
          call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,f1in,0.d0)
        endif
        pxr0=px
        pyr0=py
        zr0=z
        call tsolqur(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,
     $       bzs,0.d0,0.d0,eps0,alr)
      else
        call tsolqu(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,ak,bzs,0.d0,0.d0,0,eps0)
      endif
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        do concurrent (i=1:np)
          p=(1.d0+g(i))**2
          p=1.d0+g(i)
          a=-f1out/p
          ea=exp(a)
          b= f2out/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
        enddo
      endif
      if(fringe .and. mfring .gt. -4 .and. mfring .ne. 1)then
        call ttfrin(np,x,px,y,py,z,g,4,-ak,al,bzs)
      endif
      if(krad .and. f1out .ne. 0.d0)then
        if(photons)then
          pcvt%fr0=1.d0-.5d0*f1out/al
        endif
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,f1out,0.d0)
      endif
      call tsolrot(np,x,px,y,py,z,g,sx,sy,sz,
     $     al,bz,dx,dy,0.d0,
     $     0.d0,0.d0,theta2,bxs,bys,bzs,.false.)
      return
      end
c
      subroutine tthin(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     nord,al,ak,
     1     dx,dy,theta,krad,fringe)
      use ffs_flag
      use tmacro
      use kradlib
      use mathfun
      use multa,only:fact
      use photontable
      use kyparam, only:nmult
      implicit none
c     alpha=1/sqrt(12),beta=1/6-alpha/2,gamma=1/40-1/24/sqrt(3)
      real*8,parameter::alpha=2.88675134594812882d-1,
     1           beta =2.23290993692602255d-2,
     1           gamma=9.43738783765593145d-4,
     1           alpha1=.5d0-alpha
      integer*4 ,intent(in):: nord,np
      integer*4 i,kord
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np),
     $     dv(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,ak,dx,dy,theta
      logical*4, intent(in)::krad,fringe
      real*8 sint,cost,ala,alb,aki,akf,dpz,al1,
     $     f1,f2,f3,f4,f5,xi,yi,zi,pr,r,rcx1,pxi,rk1,rk
      complex*16 cx
c     begin initialize for preventing compiler warning
      ala=0.d0
      alb=0.d0
      aki=0.d0
c     end   initialize for preventing compiler warning
      if(ak .eq. 0.d0)then
        call tdrift_free(np,x,px,y,py,z,dv,al)
        return
      endif
      include 'inc/TENT.inc'
      if(krad)then
        pxr0=px
        pyr0=py
        zr0=z
        if(calpol)then
          bsi=0.d0
        endif
        if(photons)then
          call tsetpcvt(l_track,dx,dy,theta,0.d0,0.d0,al)
        endif
      endif
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,ak,al,0.d0)
      endif
      kord=nord/2-1
      akf=ak/fact(kord)
      if(al .ne. 0.d0)then
        ala=al*alpha1
        alb=al*alpha
        do concurrent (i=1:np)
          dpz=pxy2dpz(px(i),py(i))
          al1=ala/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*ala
        enddo
        akf=akf*.5d0
      endif
      if(kord .eq. 2)then
        do concurrent (i=1:np)
          aki=akf/(1.d0+g(i))
          px(i)=px(i)-aki*(x(i)-y(i))*(x(i)+y(i))
          py(i)=py(i)+2.d0*aki*x(i)*y(i)
        enddo
      else
        do concurrent (i=1:np)
          pr=1.d0+g(i)
          aki=akf/pr
          cx=dcmplx(x(i),-y(i))**kord
          px(i)=px(i)-aki*dble(cx)
          py(i)=py(i)-aki*imag(cx)
        enddo
      endif
      if(al .ne. 0.d0)then
        if(kord .le. 0)then
          do concurrent (i=1:np)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)*2.d0
            x(i)=x(i)+px(i)*al1
            y(i)=y(i)+py(i)*al1
            z(i)=z(i)+dpz  *al1-dv(i)*alb*2.d0
          enddo
        elseif(kord .eq. 2)then
          f1=beta*2.d0*al
          f2=4.d0*gamma*al**2
          f3=gamma*6.d0*al**2
          f4=.5d0*beta*al
          f5=f2
          do concurrent (i=1:np)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
            pr=1.d0+g(i)
            aki=akf/pr*2.d0
            r=xi  **2+yi  **2
            cx=dcmplx(xi  ,-yi  )**2
            rcx1=xi  *dble(cx)+yi  *imag(cx)
            px(i)=px(i)+aki**2*(f1*r*xi
     1        -aki*(f2*xi  *rcx1
     1             +f3*r*dble(cx)) )
            py(i)=py(i)+aki**2*(f1*r*yi
     1        -aki*(f2*yi  *rcx1
     1             +f3*r*imag(cx)) )
            zi  =zi  +aki**2*r*(f4*r-f5*aki*rcx1)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
          enddo
        elseif(kord .eq. 1)then
          f1=beta*al
          f3=gamma*2.d0*al**2
          f4=.5d0*f1
          f5=f3
          do concurrent (i=1:np)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
            r=xi**2+yi**2
            rcx1=(xi-yi)*(xi+yi)
            px(i)=px(i)+aki**2*(f1-aki*f3)*xi
            py(i)=py(i)+aki**2*(f1+aki*f3)*yi
            zi  =zi  +aki**2*(f4*r-f5*aki*rcx1)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
          enddo
        else
          f1=beta*kord*al
          f2=2.d0*gamma*kord*(kord-1)*al**2
          f3=gamma*kord*(kord+1)*al**2
          f4=.5d0*beta*al
          f5=2.d0*gamma*kord*al**2
          do concurrent (i=1:np)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            xi  =x(i)+px(i)*al1
            yi  =y(i)+py(i)*al1
            zi  =z(i)+dpz  *al1-dv(i)*alb
            pr=1.d0+g(i)
            aki=akf/pr*2.d0
            r=xi  **2+yi  **2
            rk1=r**(kord-2)
            rk=r*rk1
            cx=dcmplx(xi  ,-yi  )**kord
            rcx1=xi  *dble(cx)+yi  *imag(cx)
            px(i)=px(i)+aki**2*(f1*rk*xi
     1        -aki*(f2*rk1*xi  *rcx1
     1             +f3*rk*dble(cx)) )
            py(i)=py(i)+aki**2*(f1*rk*yi
     1        -aki*(f2*rk1*yi  *rcx1
     1             +f3*rk*imag(cx)) )
            zi  =zi  +aki**2*rk*(f4*r-f5*aki*rcx1)
            dpz=pxy2dpz(px(i),py(i))
            al1=alb/(1.d0+dpz)
            x(i)=xi  +px(i)*al1
            y(i)=yi  +py(i)*al1
            z(i)=zi  +dpz  *al1-dv(i)*alb
          enddo
        endif
        if(kord .eq. 2)then
          do concurrent (i=1:np)
            aki=akf/(1.d0+g(i))
            px(i)=px(i)-aki*(x(i)-y(i))*(x(i)+y(i))
            py(i)=py(i)+2.d0*aki*x(i)*y(i)
          enddo
        else
          do concurrent (i=1:np)
            aki=akf/(1.d0+g(i))
            cx=dcmplx(x(i),-y(i))**kord
            px(i)=px(i)-aki*dble(cx)
            py(i)=py(i)-aki*imag(cx)
          enddo
        endif
        do concurrent (i=1:np)
          dpz=pxy2dpz(px(i),py(i))
          al1=ala/(1.d0+dpz)
          x(i)=x(i)+px(i)*al1
          y(i)=y(i)+py(i)*al1
          z(i)=z(i)+dpz  *al1-dv(i)*ala
        enddo
      endif
      if(fringe)then
        call ttfrin(np,x,px,y,py,z,g,nord,-ak,al,0.d0)
      endif
      if(krad)then
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,al,0.d0)
      endif
      include 'inc/TEXIT.inc'
      return
      end
c
      subroutine ttfrin(np,x,px,y,py,z,g,nord,ak,al,bz)
      use multa,only:fact,nmult
      implicit none
      integer*4 ,intent(in):: np,nord
      integer*4 i,kord
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ,intent(in):: ak,al,bz
      real*8 akk,aki,a,b,ab,t,dx1,dy1,d,xx,yy,
     $     h,f,px1,py1,xi,yi,bzph,pr,pxa,pya
      complex*16 cx,cp,cz,cz1,cx1,cp1,ca
      real*8 ,parameter ::
     $     an(nmult+2)=[
     $     1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,
     $     11.d0,12.d0,13.d0,14.d0,15.d0,16.d0,17.d0,18.d0,19.d0,
     $     20.d0,21.d0,22.d0,23.d0]
      if(al .eq. 0.d0 .or. ak .eq. 0.d0)then
        return
      endif
      if(nord .eq. 4)then
        akk=ak/al/4.d0
        if(bz .eq. 0.d0)then
          do concurrent (i=1:np)
            aki=akk/(1.d0+g(i))
            xi=x(i)
            yi=y(i)
            a=aki*xi**2
            b=aki*yi**2
            ab=aki*(xi-yi)*(xi+yi)
            t=ab**2/6.d0
            dx1= xi*(a/3.d0+b+t)
            dy1=-yi*(a+b/3.d0-t)
            x(i)=xi+dx1
            y(i)=yi+dy1
            d=a+b
            xx=1.d0+d+ab*(5.d0*a-b)/6.d0
            yy=1.d0-d+ab*(a-5.d0*b)/6.d0
            h=2.d0*aki*xi*yi*(1.d0-ab/3.d0)
            f=xx*yy+h**2
            px1=(px(i)*yy+py(i)*h)/f
            py1=(py(i)*xx-px(i)*h)/f
            z(i)=z(i)-px1*(dx1+t*xi)-py1*(dy1+t*yi)
            px(i)=px1
            py(i)=py1
          enddo
        else
          do concurrent (i=1:np)
            pr=(1.d0+g(i))
            aki=akk/pr
            bzph=.5d0*bz/pr
            xi=x(i)
            yi=y(i)
            a=aki*xi**2
            b=aki*yi**2
            ab=aki*(xi-yi)*(xi+yi)
            t=ab**2/6.d0
            dx1= xi*(a/3.d0+b+t)
            dy1=-yi*(a+b/3.d0-t)
            x(i)=xi+dx1
            y(i)=yi+dy1
            d=a+b
            xx=1.d0+d+ab*(5.d0*a-b)/6.d0
            yy=1.d0-d+ab*(a-5.d0*b)/6.d0
            h=2.d0*aki*xi*yi*(1.d0-ab/3.d0)
            f=xx*yy+h**2
            pxa=px(i)+bzph*y(i)
            pya=py(i)-bzph*x(i)
            px(i)=(pxa*yy+pya*h)/f-bzph*yi
            py(i)=(pya*xx-pxa*h)/f+bzph*xi
            z(i)=z(i)-px(i)*(dx1+t*xi)-py(i)*(dy1+t*yi)
          enddo
        endif
      elseif(nord .eq. 2)then
        akk=ak/al
        do concurrent (i=1:np)
          aki=akk/(1.d0+g(i))
          x(i)=x(i)+.5d0*aki*y(i)**2
          py(i)=py(i)-aki*px(i)*y(i)
          z(i)=z(i)-.5d0*aki*px(i)*y(i)**2
        enddo
      elseif(nord .eq. 6)then
        akk=ak/al/24.d0
        do concurrent (i=1:np)
          aki=akk/(1.d0+g(i))
          cx=dcmplx(x(i),y(i))
          cp=dcmplx(px(i),-py(i))
          cz1=cx*cx
          cz=cz1*cx
          a=aki*imag(cz)*2.d0
          ca=-aki*cx*(cz/4.d0-conjg(cz))
          cx1=cx+ca
c          d=1.d0+a**2-(aki*4.d0)**2*
          d=1.d0+a**2-(aki*3.d0)**2*
     1         (dble(cz)**2+imag(cz)**2)
c          cp1=(dcmplx(1.d0,a)*cp-(aki*4.d0)*cz1*conjg(cx*cp))/d
          cp1=(dcmplx(1.d0,a)*cp-(aki*3.d0)*cz1*conjg(cx*cp))/d
          z(i)=z(i)-dble(ca*cp1)
          x(i)=dble(cx1)
          y(i)=imag(cx1)
          px(i)=dble(cp1)
          py(i)=-imag(cp1)
        enddo
      else
        kord=nord/2
        akk=ak/al/fact(kord)/4.d0
        do concurrent (i=1:np)
          aki=akk/(1.d0+g(i))
          cx=dcmplx(x(i),y(i))
          cp=dcmplx(px(i),-py(i))
          cz1=cx**(kord-1)
          cz=cz1*cx
          a=aki*imag(cz)*2.d0
          ca=-aki*cx*(cz/an(kord+1)-conjg(cz))
          cx1=cx+ca
          d=1.d0+a**2-(aki*an(kord))**2*
     1         (dble(cz)**2+imag(cz)**2)
          cp1=(dcmplx(1.d0,a)*cp-(aki*an(kord))*cz1*conjg(cx*cp))/d
          z(i)=z(i)-dble(ca*cp1)
          x(i)=dble(cx1)
          y(i)=imag(cx1)
          px(i)=dble(cp1)
          py(i)=-imag(cp1)
        enddo
      endif
      return
      end

      subroutine ttfrins(np,x,px,y,py,z,g,nord,ak,al,bz)
      use macmath
      implicit none
      integer*4 ,intent(in):: np,nord
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ,intent(in):: ak(2),al,bz
      integer*4 i
      real*8 x0,px0,theta,cost,sint,aka
      if((ak(1) .ne. 0.d0 .or. ak(2) .ne. 0.d0) .and. al .ne. 0.d0)then
c        theta=pi/nord
        if(ak(1) .eq. 0.d0)then
          theta=pi/nord
          aka=ak(2)
        elseif(ak(2) .eq. 0.d0)then
          theta=0.d0
          aka=ak(1)
        else
          theta=atan2(ak(2),ak(1))*2.d0/nord
          aka=hypot(ak(1),ak(2))
        endif
        if(theta .ne. 0.d0)then
          cost=cos(theta)
          sint=sin(theta)
          do concurrent (i=1:np)
            x0   =x(i)
            x(i) = cost*x0 -sint*y(i)
            y(i) = sint*x0 +cost*y(i)
            px0  =px(i)
            px(i)= cost*px0-sint*py(i)
            py(i)= sint*px0+cost*py(i)
          enddo
          call ttfrin(np,x,px,y,py,z,g,nord,aka,al,bz)
          do concurrent (i=1:np)
            x0   =x(i)
            x(i) = cost*x0 +sint*y(i)
            y(i) =-sint*x0 +cost*y(i)
            px0  =px(i)
            px(i)= cost*px0+sint*py(i)
            py(i)=-sint*px0+cost*py(i)
          enddo
        else
          call ttfrin(np,x,px,y,py,z,g,nord,aka,al,bz)
        endif
      endif
      return
      end

      subroutine tqlfri(np,x,px,y,py,z,g,f1,f2,bz)
      implicit none
      integer*4 ,intent(in):: np
      integer*4 i
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 ,intent(in):: bz,f1,f2
      real*8 xf,yf,pxf,pyf,a,b,bb,bzph,ea,f,p
      do concurrent (i=1:np)
        p=(1.d0+g(i))
        bzph=.5d0*bz/p
        a=f1/p
        ea=exp(a)
        b=f2/p
        bb=b*bzph/p
        f=1.d0/((1.d0-bb)*(1.d0+bb))
        xf=(ea*x(i)+b*(px(i)+bzph*y(i)/ea-bb*py(i)))*f
        yf=(y(i)/ea-b*(py(i)-bzph*x(i)*ea-bb*px(i)))*f
        pxf=(px(i)+bzph*yf)/ea
        pyf=(py(i)-bzph*xf)*ea
        z(i)=z(i)-(a*ea*x(i)+b*ea*(1.d0+.5d0*a)*pxf)*pxf
     $       +(a*y(i)/ea+b/ea*(1.d0-.5d0*a)*pyf)*pyf
        px(i)=pxf-bzph*y(i)
        py(i)=pyf+bzph*x(i)
        x(i)=xf
        y(i)=yf
      enddo
      return
      end

      subroutine zcheck(x,tag)
      implicit none
      real*8 x
      character*(*) tag
      write(*,*)'zcheck-0 ',tag
      if(x .ne. 0.d0)then
        write(*,*)'zcheck-x ',x
        stop
      endif
      return
      end
