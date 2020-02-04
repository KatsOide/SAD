c Obsolete
c
      subroutine tqrad(np,x,px,y,py,z,g,dv,sx,sy,sz,l,al,ak,
     1                 dx,dy,theta,cost,sint,radlvl,
     1                 f1in,f2in,f1out,f2out,mfring)
      use ffs_flag
      use tmacro
      implicit real*8(a-h,o-z)
      integer*4 np,l
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     sx(np),sy(np),sz(np)
      LOGICAL GAUSSR
      real*8 byx,bxa,bya,dprad,dpradx,dprady
      real*8 f1in,f2in,f1out,f2out
      include 'inc/TENT.inc'
      call ttfrin(np,x,px,y,py,z,g,4,ak,al,0.d0)
      al2=al*.5d0
      do 2010 i=1,np
        a=px(i)**2+py(i)**2
        dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        r=-dpz/(1.d0+dpz)*al2
        x(i)=x(i)+px(i)*r
        y(i)=y(i)+py(i)*r
        z(i)=z(i)-(3.d0+dpz)*a/2.d0/(2.d0+dpz)*r
2010  continue
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        do 2110 i=1,np
c          p=(1.d0+g(i))**2
          p=(1.d0+g(i))
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
2110    continue
      endif
      GAUSSR=RADLVL .EQ. 2.D0
      akk0=ak/al
      sqrtk=sqrt(abs(akk0))
      byx=akk0*brhoz
      b1=(byx)**2
      ur=urad*p0**3
      an=anrad*p0
      sp=0.d0
      do 100 i=1,np
        alsum=0.d0
        al0=al
c        p=g(i)*(2.d0+g(i))
        p=g(i)
        pr=1.d0+p
c        sqrtp=1.d0+g(i)
        sqrtp=sqrt(1.d0+g(i))
 110    brad0=x(i)**2+y(i)**2
        brad1=x(i)*px(i)+y(i)*py(i)
        brad2=(px(i)**2+py(i)**2-akk0*(x(i)-y(i))*(x(i)+y(i)))/3.d0
        brad=b1*(brad0+al0*(brad1+al0*brad2))
        if(brad .ne. 0.d0)then
          rho1=brhoz/sqrt(brad)
          alx=min(al0,.09d0*rho1/an)
111       brad=b1*(brad0+alx*(brad1+alx*brad2))
          rho=brhoz/sqrt(brad)
          alr=alx*(1.d0+(px(i)**2+py(i)**2)*.5d0)
          prob=alr*an/rho
          if(prob .gt. .1d0)then
            alx=alx*(.09d0/prob)
            go to 111
          endif
          if(tran() .gt. 1.d0-prob)then
            bxa=byx*(y(i)*(1.d0+akk0*alx**2/6.d0)+py(i)*alx*.5d0)
            bya=byx*(x(i)*(1.d0-akk0*alx**2/6.d0)+px(i)*alx*.5d0)
            call tsynchrad(pr,alr,bxa,bya,
     $           dprad,dpradx,dprady,
     $           i,l,alsum,0.d0,0.d0,theta,
     $           x(i),y(i),px(i),py(i))
            p=max(-.999d0,p-dprad*.5d0)
            px(i)=px(i)-dpradx*.5d0
            py(i)=py(i)-dprady*.5d0
            pr=1.d0+p
            sp=sp-dprad
            sqrtp=sqrt(pr)
          else
            dprad=0.d0
          endif
        else
          alx=al0
          dprad=0.d0
        endif
        alsum=alsum+alx
        if(ak .gt. 0.d0)then
          akk=sqrtk/sqrtp
          phi=akk*alx*.5d0
          s=phi**2
          t=1.d0-s/(3.d0-s/(5.d0-s/(7.d0-s/9.d0)))
          th=1.d0+s/(3.d0+s/(5.d0+s/(7.d0+s/9.d0)))
          u=t**2+phi**2
          a11=(t-phi)*(t+phi)/u
          a12=2.d0*t*phi/u
          u=(th-phi)*(th+phi)
          b11=(th**2+phi**2)/u
          b12=2.d0*th*phi/u
          a21=-a12*akk
          a12=a12/akk
          b21=b12*akk
          b12=b12/akk
c    Here I ignore the change of dv(i)
          xi=x(i)
          yi=y(i)
          x(i)=a11*xi+a12*px(i)
          y(i)=b11*yi+b12*py(i)
          z(i)=z(i)-dv(i)*alx-
     1       (((px(i)**2+(xi*akk)**2)*(alx+a11*a12)+
     1         (py(i)-yi*akk)*(py(i)+yi*akk)*(alx+b11*b12))*.5d0+
     1         xi*x(i)*a21+yi*y(i)*b21)*.5d0
          px(i)=a21*xi+a11*px(i)
          py(i)=b21*yi+b11*py(i)
        elseif(ak .lt. 0.d0)then
          akk=sqrtk/sqrtp
          phi=akk*alx
          phi=akk*alx*.5d0
          s=phi**2
          t=1.d0-s/(3.d0-s/(5.d0-s/(7.d0-s/9.d0)))
          th=1.d0+s/(3.d0+s/(5.d0+s/(7.d0+s/9.d0)))
          u=t**2+phi**2
          b11=(t-phi)*(t+phi)/u
          b12=2.d0*t*phi/u
          u=(th-phi)*(th+phi)
          a11=(th**2+phi**2)/u
          a12=2.d0*th*phi/u
          b21=-b12*akk
          b12=b12/akk
          a21=a12*akk
          a12=a12/akk
          xi=x(i)
          yi=y(i)
          x(i)=xi*a11+px(i)*a12
          y(i)=yi*b11+py(i)*b12
          z(i)=z(i)-dv(i)*alx-
     1       (((px(i)-xi*akk)*(px(i)+xi*akk)*(alx+a11*a12)+
     1         (py(i)**2+(yi*akk)**2)*(alx+b11*b12))*.5d0+
     1         xi*x(i)*a21+yi*y(i)*b21)*.5d0
          px(i)=a21*xi+a11*px(i)
          py(i)=b21*yi+b11*py(i)
        endif
        if(dprad .ne. 0.d0)then
          p=max(-.999d0,p-dprad*.5d0)
          px(i)=px(i)-dpradx*.5d0
          py(i)=py(i)-dprady*.5d0
          pr=1.d0+p
          sqrtp=sqrt(pr)
        endif
        al0=al0-alx
        if(al0 .gt. 0.d0)then
          go to 110
        endif
        g(i)=p
100   continue
      if(radcod)then
        sp=0.d0
      else
        sp=sp/np
      endif
      do 310 i=1,np
        dp=max(-.99d0,g(i)-sp)
        pr=1.d0+dp
        r=(1.d0+g(i))/pr
c       px(i)=px(i)*r
c       py(i)=py(i)*r
c        g(i)=dp/(1.d0+sqrt(pr))
        g(i)=dp
        p1=pr*p0
        h=sqrt(1.d0+p1**2)
        dv(i)=-dp*(1.d0+pr)/h/(h+pr*h0)+dvfs
310   continue
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        do 2120 i=1,np
c          p=(1.d0+g(i))**2
          p=(1.d0+g(i))
          a=f1out/p
          ea=exp(a)
          b=f2out/p
          pxf=px(i)/ea
          pyf=py(i)*ea
          z(i)=z(i)-(a*x(i)+b*(1.d0+.5d0*a)*pxf)*px(i)
     $             +(a*y(i)+b*(1.d0-.5d0*a)*pyf)*py(i)
          x(i)=ea*x(i)+b*px(i)
          y(i)=y(i)/ea-b*py(i)
          px(i)=pxf
          py(i)=pyf
2120    continue
      endif
      do 2020 i=1,np
        a=px(i)**2+py(i)**2
        dpz=a*(-.5d0-a*(.125d0+a*.0625d0))
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        dpz=(dpz**2-a)/(2.d0+2.d0*dpz)
        r=-dpz/(1.d0+dpz)*al2
        x(i)=x(i)+px(i)*r
        y(i)=y(i)+py(i)*r
        z(i)=z(i)-(3.d0+dpz)*a/2.d0/(2.d0+dpz)*r
2020  continue
      call ttfrin(np,x,px,y,py,z,g,4,-ak,al,0.d0)
      include 'inc/TEXIT.inc'
      return
      end
