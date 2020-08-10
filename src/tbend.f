      module multa
        integer*4, parameter :: nmult=21
        real*8 ,parameter :: fact(0:nmult+1)=[
     $       1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1       720.d0,     5040.d0,     40320.d0,362880.d0,3628800.d0,
     $       39916800.d0,479001600.d0,6227020800.d0,87178291200.d0,
     $       1307674368000.d0,20922789888000.d0,355687428096000.d0,
     $       6402373705728000.d0,121645100408832000.d0,
     $       2432902008176640000.d0,51090942171709440000.d0,
     $       1124000727777607680000.d0],
     $       an(0:nmult+1)=[1.d0,1.d0,
     $       0.5d0,
     $       0.33333333333333333333d0,
     $       0.25d0,
     $       0.2d0,
     $       0.166666666666666666667d0,
     $       0.142857142857142857143d0,
     $       0.125d0,
     $       0.111111111111111111111d0,
     $       0.1d0,
     $       0.090909090909090909091d0,
     $       0.083333333333333333333d0,
     $       0.076923076923076923077d0,
     $       0.071428571428571428571d0,
     $       0.066666666666666666667d0,
     $       0.0625d0,
     $       0.058823529411764705882d0,
     $       0.055555555555555555556d0,
     $       0.052631578947368421053d0,
     $       0.05d0,
     $       0.047619047619047619048d0,
     $       0.045454545454545454545d0]
        logical*4 :: gknini=.true.
        real*8 :: gkn(0:nmult,0:nmult)=0.d0

        contains
        subroutine gkninit
        implicit none
        integer*4 n,k
        do n=0,nmult
          gkn(n,0)=1.d0
          do k=1,nmult-n
            gkn(n,k)=gkn(n,k-1)
     $           *dble((2*k-3)*(2*k+1))/8.d0/dble(k*(k+n+1))
          enddo
        enddo
        gknini=.false.
        return
        end subroutine
      end module

      module tbendcom
        real*8 rho0,rhob,f1r,f2r,fb1,fb2

        contains
        pure subroutine tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
        use tfstk
        use ffs_flag, only:calpol
        use mathfun
        implicit none
        integer*4 ,intent(in):: np
        integer*4 i
        real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $       sx(np),sy(np),sz(np)
        real*8 ,intent(in):: phi0,dtheta
        real*8 r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $       pxi,pyi,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf,
     $       cphi0,sphi0,cdt,sdt,sdth2,sxf,syf
        cphi0=cos(phi0*.5d0)
        sphi0=sin(phi0*.5d0)
        sdth2=sin(dtheta*.5d0)**2
        cdt=1.d0-2.d0*sdth2
        sdt=sin(dtheta)
        r11=cdt*cphi0**2+sphi0**2
        r12=-cphi0*sdt
        r13=-2.d0*sdth2*sphi0*cphi0
        r21=-r12
        r22=cdt
        r23=sphi0*sdt
        r31=r13
        r32=-r23
        r33=cdt*sphi0**2+cphi0**2
        if(calpol)then
          do concurrent (i=1:np)
            xi=x(i)
            yi=y(i)
            pxi=px(i)
            pyi=py(i)
            pzi=1.d0+pxy2dpz(pxi,pyi)
            xf =r11*xi +r12*yi
            yf =r21*xi +r22*yi
            zf =r31*xi +r32*yi
            pxf=r11*pxi+r12*pyi+r13*pzi
            pyf=r21*pxi+r22*pyi+r23*pzi
            pzf=r31*pxi+r32*pyi+r33*pzi
            px(i)=pxf
            py(i)=pyf
            x(i)=xf-pxf/pzf*zf
            y(i)=yf-pyf/pzf*zf
            z(i)=z(i)+zf/pzf
            sxf  =r11*sx(i)+r12*sy(i)+r13*sz(i)
            syf  =r21*sx(i)+r22*sy(i)+r23*sz(i)
            sz(i)=r31*sx(i)+r32*sy(i)+r33*sz(i)
            sx(i)=sxf
            sy(i)=syf
          enddo
        else
          do concurrent (i=1:np)
            xi=x(i)
            yi=y(i)
            pxi=px(i)
            pyi=py(i)
            pzi=1.d0+pxy2dpz(pxi,pyi)
            xf =r11*xi +r12*yi
            yf =r21*xi +r22*yi
            zf =r31*xi +r32*yi
            pxf=r11*pxi+r12*pyi+r13*pzi
            pyf=r21*pxi+r22*pyi+r23*pzi
            pzf=r31*pxi+r32*pyi+r33*pzi
            px(i)=pxf
            py(i)=pyf
            x(i)=xf-pxf/pzf*zf
            y(i)=yf-pyf/pzf*zf
            z(i)=z(i)+zf/pzf
          enddo
        endif
        return
        end subroutine

        pure subroutine tbshift(np,x,px,y,py,z,dx,dy,phi0,cost,sint,ent)
        use tfstk
        use mathfun
        implicit none
        integer*4 ,intent(in):: np
        integer*4 i
        real*8 , intent(in)::dx,dy,phi0,cost,sint
        real*8 , intent (inout)::x(np),px(np),y(np),py(np),z(np)
        real*8 phih,dcph,sph,ds,dx1,dy1,dxa,st1,x1,px1,al,y1
        logical*4 , intent(in)::ent
        if(dx .ne. 0.d0 .or. dy .ne. 0.d0)then
          phih=phi0*.5d0
          sph=sin(phih)
          dcph=-2.d0*sin(phih*.5d0)**2
          if(ent)then
            st1=sint
          else
            st1=-sint
          endif
          dxa=dx*cost-dy*st1
          dx1=dx+dcph*cost*dxa
          dy1=dy-dcph*st1 *dxa
          ds=dxa*sph
          dxa=dxa*(1.d0+dcph)
          if(ent)then
            do concurrent (i=1:np)
              al=ds/(1.d0+pxy2dpz(px(i),py(i)))
              x1  =x(i)+px(i)*al-dx1
              y(i)=y(i)+py(i)*al-dy1
              z(i)=z(i)-al
              x(i)=x1*cost-y(i)*sint
              y(i)=x1*sint+y(i)*cost
              px1=px(i)
              px(i)=px1*cost-py(i)*sint
              py(i)=px1*sint+py(i)*cost
            enddo
          else
            do concurrent (i=1:np)
              x1  =x(i)*cost-y(i)*sint
              y(i)=x(i)*sint+y(i)*cost
              px1=px(i)
              px(i)=px1*cost-py(i)*sint
              py(i)=px1*sint+py(i)*cost
              al=ds/(1.d0+pxy2dpz(px(i),py(i)))
              x(i)=x1  +px(i)*al-dx1
              y(i)=y(i)+py(i)*al-dy1
              z(i)=z(i)-al
            enddo
          endif
        elseif(sint .ne. 0.d0 .or. cost .ne. 1.d0)then
          do concurrent (i=1:np)
            y1=y(i)
            y(i)=x(i)*sint+y1*cost
            x(i)=x(i)*cost-y1*sint
            px1=px(i)
            px(i)=px1*cost-py(i)*sint
            py(i)=px1*sint+py(i)*cost
          enddo
        endif
        return
        end subroutine 

      end module

      module bendib
      use tfstk
      implicit none
      real*8 xi,pxi,yi,pyi,zi,p,dp,rhoe,rho0,rhob,drhob,drhop,
     $     rhosq,
     $     akk,akxsq,akysq,akx,aky,dcx,aksx,dcy,aksy,phix,phiy,
     $     spx,spy,sxkx,syky,dcxkx,xsxkx
      real*8 , save :: akxi=0.d0,alxi=0.d0,dpxi=0.d0
      real*8 ,parameter :: a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1     a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1     a13=231.d0/13312.d0,a15=143.d0/10240.d0,
     $     psqmax=0.9999d0,epsbend=1.d-3,
     $     rbh=.5d0+1.d0/sqrt(12.d0),rbl=1.d0/6.d0/rbh

      contains

      subroutine tbendal(n,ndiv,f1r,f2r,aln,alx,alr)
      implicit none
      integer*4 , intent(in)::n,ndiv
      real*8 , intent(in) ::f1r,f2r,aln
      real*8 , intent(out)::alx,alr
      if(n .gt. 0 .and. n .le. ndiv)then
        alx=aln
        alr=aln
      elseif(n .eq. -1)then
        alx=rbl*f1r
        alr=f1r
      elseif(n .eq. 0)then
        alx=rbh*f1r
        alr=f1r
      elseif(n .eq. ndiv+1)then
        alx=rbh*f2r
        alr=f2r
      else
        alx=rbl*f2r
        alr=f2r
      endif
      return
      end subroutine

      subroutine tbendiinit(ak1,al,force)
      use mathfun
      implicit none
      real*8, intent (in):: ak1,al
      logical*4 , intent(in), optional::force
      real*8 xspx
c      if(present(force))then
c        write(*,*)'tbendiinit ',force,ak1,akxi,al,alxi,dpxi,dp
c      else
c        write(*,*)'tbendiinit-noforce ',ak1,akxi,al,alxi,dpxi,dp
c      endif
      if(ak1 .eq. akxi .and. al .eq. alxi .and. dpxi .eq. dp)then
        if(.not. present(force) .or. .not. force)then
          return
        endif
      endif
      akxi=ak1
      alxi=al
      dpxi=dp
      p=1.d0+dp
      rhosq=rho0*rhoe
      drhop=(rhoe-rho0)/rhosq
      akk=ak1/al
      akxsq=-1.d0/rhosq-akk/p
      akysq=akk/p
      if(akxsq .gt. 0.d0)then
        akx=sqrt(akxsq)
        phix=akx*al
        dcx=2.d0*sinh(.5d0*phix)**2
        call sxsinh(phix,spx,xspx)
c        spx=sinh(phix)
        aksx=akx*spx
        dcxkx=dcx/akxsq
        sxkx=spx/akx
        xsxkx=xspx/akx/akxsq
      elseif(akxsq .lt. 0.d0)then
        akx=sqrt(-akxsq)
        phix=akx*al
        dcx=-2.d0*sin(.5d0*phix)**2
        call sxsin(phix,spx,xspx)
c        spx=sin(phix)
        aksx=-akx*spx
        dcxkx=dcx/akxsq
        sxkx=spx/akx
        xsxkx=-xspx/akx/akxsq
      else
        akx=0.d0
        phix=0.d0
        dcx=0.d0
        spx=0.d0
        aksx=0.d0
        dcxkx=0.5d0*al**2
        sxkx=al
        xsxkx=-1.d0/6.d0*al**3
      endif
      if(akysq .gt. 0.d0)then
        aky=sqrt(akysq)
        phiy=aky*al
        dcy=2.d0*sinh(.5d0*phiy)**2
        spy=sinh(phiy)
        aksy=aky*spy
        syky=spy/aky
      elseif(akysq .lt. 0.d0)then
        aky=sqrt(-akysq)
        phiy=aky*al
        dcy=-2.d0*sin(.5d0*phiy)**2
        spy=sin(phiy)
        aksy=-aky*spy
        syky=spy/aky
      else
        aky=0.d0
        phiy=0.d0
        dcy=0.d0
        spy=0.d0
        aksy=0.d0
        syky=al
      endif
      return
      end subroutine

      subroutine tbendicorr(ak1,al,phi)
      use mathfun
      implicit none
      real*8 , intent(in) :: ak1,al,phi
      real*8 pzi,dpzi,xr,s
      s=min(psqmax,pxi**2+pyi**2)
      dpzi=sqrt1(-s)
      pzi=1.d0+dpzi
      xi =(-al*dpzi*pxi+pzi*xi)/(pzi-phi*pxi)
      pxi=pxi+phi*dpzi
      yi =yi+(-al*dpzi+xi*phi)*pyi/pzi
      zi =zi-phi*dpzi*xi/pzi+al*dpzi*(.5d0*pzi+.5d0-1.d0/pzi)
      xr=xi/rho0
      pxi=pxi-ak1*xi*xr*(0.5d0-xr*(2.d0-xr)/12.d0)/p
      return
      end subroutine

      subroutine tbendibody(al)
      implicit none
      real*8, intent(in):: al
      real*8 dxf,dpxf,dyf,dpyf,hi,xiksq
c      dxf =(drhopak+xi)*dcx+aksx/akxsq*pxi
c      dpxf= aksx*(drhopak+xi)+pxi*dcx
c      dyf = dcy*yi+aksy/akysq*pyi
c      dxf = drhop*dcxkx+xi*dcx+sxkx*pxi
      xiksq = drhop    +akxsq*xi
      dxf = dcxkx*xiksq+sxkx*pxi
      dpxf= sxkx*xiksq +dcx*pxi
      dyf = dcy*yi     +syky*pyi
      dpyf= aksy*yi    +dcy*pyi
      hi=.5d0*(pxi**2+pyi**2-akxsq*xi**2-akysq*yi**2)
     $     -drhop*xi
      zi =zi-.25d0*(dxf*pxi+xi*dpxf+dxf*dpxf
     $     +dyf*pyi+yi*dpyf+dyf*dpyf
     $     -(5.d0/rho0-1.d0/rhoe)*(drhop*xsxkx-sxkx*xi-dcxkx*pxi))
     $     -hi*al*.5d0
      xi =xi +dxf
      pxi=pxi+dpxf
      yi =yi +dyf
      pyi=pyi+dpyf
      return
      end subroutine

      end module

      subroutine tbend(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,phib,phi0,psi1,psi2,
     1     cosp1,sinp1,cosp2,sinp2,
     1     ak,dx,dy,theta,dtheta,cost,sint,
     1     fb10,fb20,mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     krad,eps,ini,iniph)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:inext,iprev
      use multa, only:nmult
      use tbendcom
      use kradlib
      use photontable
      use mathfun, only:akang
      implicit none
      integer*4 ,parameter::ndivmax=1024
      real*8 ,parameter:: a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0
      real*8 ,parameter::smax=0.99d0,smin=0.01d0,rphidiv=3e-3
      integer*4 ,intent(in):: np,mfring
      real*8 ,intent(in):: al,phib,phi0,cosp1,sinp1,cosp2,sinp2,ak,
     $     dx,dy,theta,cost,sint,cosw,sinw,sqwh,sinwp1,eps,
     $     psi1,psi2,fb10,fb20,dtheta
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),g(np),sx(np),sy(np),sz(np)
      integer*4 ndiv,iniph,n1,n2
      real*8 theta2,phir
      complex*16 akm(0:nmult),cr1
      logical*4 ,intent(in):: fringe,ini,krad
      if(phi0 .eq. 0.d0)then
        if(ak .eq. 0.d0)then
          call tsteer(np,x,px,y,py,z,g,dv,sx,sy,sz,al,-phib,
     1         dx,dy,theta+dtheta,
     1         cosp1,sinp1,cosp2,sinp2,
     $         fb10,fb20,fringe,eps,krad)
        elseif(phib .eq. 0.d0)then
          call tquad(np,x,px,y,py,z,g,dv,sx,sy,sz,al,ak,0.d0,
     1         dx,dy,theta+dtheta,theta+dtheta,krad,.true.,
     1         fringe,0.d0,0.d0,0,eps,.true.)
        else
          akm=0.d0
          akm(0)=phib-phi0
          akm(1)=ak
          theta2=theta+dtheta+akang(dcmplx(ak,0.d0),al,cr1)
          call tmulti(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $         al,akm,0.d0,0.d0,0.d0,0.d0,
     $         dx,dy,0.d0,0.d0,0.d0,theta+dtheta,0.d0,
     $         theta2,cr1,
     $         eps,krad,fringe,
     $         0.d0,0.d0,0.d0,0.d0,
     $         mfring,fb10,fb20,
     $         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     $         .false.,.false.,0,i00)
        endif
        return
      elseif(ak .ne. 0.d0)then
        call tbendi(np,x,px,y,py,z,g,dv,sx,sy,sz,al,phib,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       ak,dx,dy,theta,dtheta,cost,sint,
     1       fb10,fb20,mfring,krad,fringe,eps)
        return
      endif
      call tbshift(np,x,px,y,py,z,dx,dy,phi0,cost,sint,.true.)
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      if(phib .eq. 0.d0)then
        call tbdrift(np,x,px,y,py,z,dv,sx,sz,al,phi0)
        go to 9000
      elseif(al .eq. 0.d0)then
        call tbthin(np,x,px,y,py,z,g,sx,sy,sz,phib,phi0,dx,dy,
     1              dtheta,cost,sint)
        go to 9000
      endif
      rhob=al/phib
      rho0=al/phi0
      fb1=fb10
      fb2=fb20
      n1=1
      n2=0
      ndiv=1
      f1r=0.d0
      f2r=0.d0
      if(krad)then
        if(ini)then
          if(photons .and. iniph .eq. 0)then
            call tsetpcvt(l_track,dx,dy,theta,dtheta,phi0,al)
          endif
          pxr0=px
          pyr0=py
          zr0=z
        endif
        if(iprev(l_track) .eq. 0 .and. fb1 .ne. 0.d0)then
          n1=-1
          f1r=fb1*.5d0
        endif
        if(inext(l_track) .eq. 0 .and. fb2 .ne. 0.d0)then
          n2=2
          f2r=fb2*.5d0
        endif
        phir=phib*(al-f1r-f2r)/al
        ndiv=min(ndivmax,ndivrad(phir,0.d0,0.d0,eps))
c        if(photons)then
c          select case (iniph)
c            case (0)
c              call tsetphotongeo(al/ndiv,phi0/ndiv,theta,.true.)
c            case (1)
c              call tsetphotongeo(al/ndiv,phi0/ndiv,pp%theta,.true.)
c            case default
c              call tsetphotongeo(al/ndiv,phi0/ndiv,0.d0,.false.)
c          end select
c        endif
      endif
      n2=ndiv+n2
      if(n1 .eq. n2)then
        call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,phi0,
     1       cosp1,sinp1,cosp2,sinp2,
     1       mfring,fringe,
     $       cosw,sinw,sqwh,sinwp1,
     1       krad,al,1.d0,1.d0)
      else
        call tbendr(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       al,phib,phi0,psi1,psi2,
     1       cosp1,sinp1,cosp2,sinp2,f1r,f2r,
     1       mfring,fringe,n1,n2,ndiv)
      endif
 9000 if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      call tbshift(np,x,px,y,py,z,-dx,-dy,-phi0,cost,-sint,.false.)
      return
      end

      subroutine tbdrift(np,x,px,y,py,z,dv,sx,sz,al,phi0)
      use tfstk
      use ffs_flag, only:calpol
      use mathfun
      implicit none
      integer*4 np,i
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),sx(np),sz(np)
      real*8 ,intent(in):: al,phi0
      real*8 cp,sp,rho0,dx,xi,pzi,pzf,dl,dcp,xsp
c      th=tan(.5d0*phi0)
c      sp=2.d0*th/(1.d0+th**2)
c      dcp=th*sp
c      cp=1.d0-dcp
c      cp=cos(phi0)
      call xsincos(phi0,sp,xsp,cp,dcp)
c      sp=sin(phi0)
c      if(cp .ge. 0.d0)then
c        dcp=sp**2/(1.d0+cp)
c      else
c        dcp=1.d0-cp
c      endif
      rho0=al/phi0
      call tdrift_free(np,x,px,y,py,z,dv,rho0*sp)
      dx=-rho0*dcp
      dl=rho0*xsp
      if(calpol)then
        do concurrent (i=1:np)
          xi=x(i)+dx
          pzi=1.d0+pxy2dpz(px(i),py(i))
          pzf=pzi*cp-px(i)*sp
          x(i)=xi*pzi/pzf
          y(i)=y(i)+xi*sp*py(i)/pzf
          z(i)=z(i)-xi*sp/pzf+(1.d0-dv(i))*dl
          px(i)=px(i)*cp+pzi*sp
          sx(i)= cp*sx(i)+sp*sz(i)
          sz(i)=(sz(i)-sp*sx(i))/cp
        enddo
      else
        do concurrent (i=1:np)
          xi=x(i)+dx
          pzi=1.d0+pxy2dpz(px(i),py(i))
          pzf=pzi*cp-px(i)*sp
          x(i)=xi*pzi/pzf
          y(i)=y(i)+xi*sp*py(i)/pzf
          z(i)=z(i)-xi*sp/pzf+(1.d0-dv(i))*dl
          px(i)=px(i)*cp+pzi*sp
        enddo
      endif
      return
      end

      subroutine tbthin(np,x,px,y,py,z,g,sx,sy,sz,phib,phi0,dx,dy,
     1                 dtheta,cost,sint)
      use tbendcom, only:tbrot,tbshift
      implicit none
      integer*4 ,intent(in):: np
      integer*4 i
      real*8 ,intent(in):: phib,phi0,dx,dy,cost,sint,dtheta
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     g(np),sx(np),sy(np),sz(np)
      call tbshift(np,x,px,y,py,z,dx,dy,phi0,cost,sint,.true.)
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,phi0,dtheta)
      endif
      do concurrent (i=1:np)
        px(i)=px(i)+phi0-phib/(1.d0+g(i))
        z(i)=z(i)-x(i)*phi0
      enddo
      if(dtheta .ne. 0.d0)then
        call tbrot(np,x,px,y,py,z,sx,sy,sz,-phi0,-dtheta)
      endif
      call tbshift(np,x,px,y,py,z,-dx,-dy,-phi0,cost,-sint,.false.)
      return
      end

      subroutine tbendr(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,phib,phi0,psi1,psi2,
     1     cosp1,sinp1,cosp2,sinp2,f1r,f2r,
     1     mfring,fringe,n1,n2,ndiv)
      use tfstk
      use ffs_flag
      use tmacro
      use tspin
      use photontable
      use bendib, only:rbh,rbl,tbendal
      use mathfun, only:xsincos
      implicit none
      integer*4 ,intent(in):: np,mfring,ndiv,n1,n2
      integer*4 mfr1,n
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),g(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,phib,phi0,cosp1,sinp1,cosp2,sinp2,
     $     psi1,psi2,f1r,f2r
      real*8 wn,aln,phibn,phi0n,alr,
     $     coswn,sinwn,sqwhn,sinwpn,bsi1,bsi2,alx,xsinwn,
     $     cosp1n,sinp1n,cosp2n,sinp2n,psi1n,psi2n
      logical*4 ,intent(in):: fringe
      aln=(al-f1r-f2r)/ndiv
      do n=n1,n2
        mfr1=0
        psi1n=0.d0
        psi2n=0.d0
        cosp1n=1.d0
        sinp1n=0.d0
        cosp2n=1.d0
        sinp2n=0.d0
        bsi1=0.d0
        bsi2=0.d0
        call tbendal(n,ndiv,f1r,f2r,aln,alx,alr)
        if(n .eq. n1)then
          if(mfring .gt. 0 .or. mfring .eq. -1)then
            mfr1=-1
          endif
          psi1n=psi1
          cosp1n=cosp1
          sinp1n=sinp1
          bsi1=1.d0
        elseif(n .eq. n2)then
          if(mfring .gt. 0 .or. mfring .eq. -2)then
            mfr1=-2
          endif
          psi2n=psi2
          cosp2n=cosp2
          sinp2n=sinp2
          bsi2=1.d0
        endif
        phi0n=alx/al*phi0
        phibn=alx/al*phib
        if(n .le. 2 .or. n .ge. ndiv)then
          wn=phi0n-psi1n-psi2n
          call xsincos(wn,sinwn,xsinwn,coswn,sqwhn)
          sinwpn=sin(phi0n-psi2n)
        endif
        call tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       alx,phi0n,
     1       cosp1n,sinp1n,cosp2n,sinp2n,
     1       mfr1,fringe,
     $       coswn,sinwn,-sqwhn,sinwpn,
     1       .true.,alr,bsi1,bsi2)
        pcvt%fr0=pcvt%fr0+alx/pcvt%al
      enddo
      return
      end

      subroutine tbendcore(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,phi0,
     1     cosp1,sinp1,cosp2,sinp2,
     1     mfring,fringe,
     $     cosw,sinw,sqwh,sinwp1,
     1     krad,alr,bsi1,bsi2)
      use tfstk
      use ffs_flag
      use tmacro
      use multa, only:nmult
      use tbendcom
      use kradlib
      use mathfun
      implicit none
      integer*4 ,intent(in):: np,mfring
      integer*4 i
      real*8 ,intent(in):: al,phi0,cosp1,sinp1,cosp2,sinp2,
     $     cosw,sinw,sqwh,sinwp1,alr
      real*8 drhob,dp,p,
     $     pinv,rhoe,pxi,pyi,dpzi,pzi,sp1,x1,dz1,y1,z1,px1,
     $     py1,pv1sqi,f,ff,x2,py2,z2,dph2,ph2,dpx2,pz2,drho,
     $     t2,dpx3,px3,dpz3,pz3,t3,x3,da,y3,z3,pv2sqi,x4,py4,z4,dpz4,
     $     dz4,dxfr1,dyfr1,dzfr1,dxfr2,dyfr2,dzfr2,dpz32,
     $     dyfra1,dyfra2,fa,t4,dpx3a,t2t3,dcosp,px1px3,
     $     phi0a,bsi1,bsi2
      real*8, parameter :: smin=1.d-4
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),g(np),sx(np),sy(np),sz(np)
      logical*4 ,intent(in):: krad,fringe
      if((mfring .gt. 0 .or. mfring .eq. -1) .and. fb1 .ne. 0.d0)then
        dxfr1=fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        dzfr1=dxfr1*sinp1
        if(fringe)then
          dyfra1=4.d0*dyfr1/fb1**2
        else
          dyfra1=0.d0
        endif
      else
        dxfr1=0.d0
        dyfr1=0.d0
        dyfra1=0.d0
        dzfr1=0.d0
      endif
      if((mfring .gt. 0 .or. mfring .eq. -2) .and. fb2 .ne. 0.d0)then
        dxfr2=fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        dzfr2=dxfr2*sinp2
        if(fringe)then
          dyfra2=4.d0*dyfr2/fb2**2
        else
          dyfra2=0.d0
        endif
      else
        dxfr2=0.d0
        dyfr2=0.d0
        dyfra2=0.d0
        dzfr2=0.d0
      endif
      if(cosp1*cosp2 .gt. 0.d0)then
        dcosp=(sinp2-sinp1)*(sinp2+sinp1)/(cosp1+cosp2)
      else
        dcosp=cosp1-cosp2
      endif
      drhob=rhob-rho0
      do concurrent (i=1:np)
        bsi(i)=bsi(i)+bsi1*y(i)/rhob
        dp=g(i)
        p=1.d0+dp
        pinv=1.d0/p
        rhoe=rhob*p
        pxi=px(i)
        pyi=py(i)
        dpzi=pxy2dpz(pxi,pyi)
        pzi=1.d0+dpzi
        sp1=sinp1/pzi
        x1=x(i)/(cosp1-pxi*sp1)
        dz1=x1*sp1
        y1=y(i)+pyi*dz1
        z1=z(i)-dz1
        px1= pxi*cosp1+pzi*sinp1
        x1=x1+dxfr1*dp*pinv
        py1=pyi+(dyfr1-dyfra1*y1**2)*y1*pinv**2
        z1=z1+(dxfr1*px1+
     $       (.5d0*dyfr1-.25d0*dyfra1*y1**2)*y1**2*pinv)*pinv-dzfr1
        pv1sqi=1.d0/max(smin,1.d0-px1**2)
        fa=y1/rhoe*sqrt(pv1sqi)
        f=(1.d0-(y1/rhob)**2/6.d0)*fa
        ff=.25d0*(f+fa)*y1*pv1sqi
        x2=x1+ff
        py2=py1-px1*f
        z2=z1-px1*ff
        dph2=pxy2dpz(0.d0,py2)
        ph2=1.d0+dph2
        dpx2=pxi*cosp1+(dpzi-dph2)*sinp1
        pz2=1.d0+pxy2dpz(px1,py2)
        drho=drhob+rhoe*dph2+rhob*dp
        t2=(px1+ph2*sinp1)/(pz2+ph2*cosp1)
        dpx3a=(x2*sinw-drho*(sinp2+sinwp1))/rhoe
        dpx3=dpx3a-dpx2*(cosw-sinw*t2)
        px3=ph2*sinp2+dpx3
        px1px3=ph2*(sinp2+sinp1)+dpx3a+dpx2*(sqwh+sinw*t2)
        dpz3=pxy2dpz(px3,py2)
        pz3=1.d0+dpz3
        dpz32=px1px3*(px1-px3)/(pz2+pz3)
        t3=(px3+ph2*sinp2)/(pz3+ph2*cosp2)
        t2t3=(ph2*sinp2+px1px3)/(pz3+ph2*cosp2)
     $       +ph2*sinp1/(pz2+ph2*cosp1)
     $       +px1*(dpz32-ph2*dcosp)
     $       /(pz3+ph2*cosp2)/(pz2+ph2*cosp1)
        t4=(cosp2+t3*sinp2)*(pz2*cosp1+px1*sinp1)
        x3=x2*(cosw-rho0/rhoe*t3*sinw/ph2)
     1       +(rho0*(cosw*t2t3+sinw*(1.d0-t2*t3))*dpx2-
     1       drho*(-(sinp2+sinwp1)*rho0/rhoe*t3
     $       -dpz32-sqwh*pz2-sinw*px1))/ph2
        da=asin(min(1.d0,max(-1.d0,
     $       (dpx2*(
     $       sinp1*(t2*(pz3*cosp2+px3*sinp2)-px1*cosp2)
     $       -sinp2*(t3*(pz2*cosp1+px1*sinp1)-px3*cosp1)
     $       +cosp1*cosp2*dpz32
     $       +t4*(sqwh+sinw*t2))
     1       +dpx3a*t4)/ph2**2)))
        phi0a=phi0+da
        y3=y1+py2*rhoe*phi0a
        z3=z2-phi0*(dp*rhob+drhob)-da*rhoe-dv(i)*al
        pv2sqi=1.d0/max(smin,1.d0-px3**2)
        fa=y3/rhoe*sqrt(pv2sqi)
        f=(1.d0-(y3/rhob)**2/6.d0)*fa
        ff=.25d0*(f+fa)*y3*pv2sqi
        py4=py2-px3*f
        z4=z3-px3*ff
        x4=x3-ff-dxfr2*dp*pinv
        py4=py4+(dyfr2-dyfra2*y3**2)*y3*pinv**2
        z4=z4+(dxfr2*px3+
     $       (.5d0*dyfr2-.25d0*dyfra2*y3**2)*y3**2*pinv)*pinv-dzfr2
        dpz4=pxy2dpz(px3,py4)
        px(i)=-cosp2*dpx3+sinp2*(dpz4-dpz3-dpx3*t3)
        dz4=x4*sinp2/(cosp2*(1.d0+dpz4)+sinp2*px3)
        x(i)=x4*cosp2+px(i)*dz4
        py(i)=py4
        y(i)=y3+py4*dz4
        z(i)=z4-dz4
        bsi(i)=bsi(i)-bsi2*y(i)/rhob
      enddo
      if(krad)then
        call tradk(np,x,px,y,py,z,g,dv,sx,sy,sz,alr,phi0)
      endif
      return
      end
