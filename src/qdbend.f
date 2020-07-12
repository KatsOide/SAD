      subroutine qdbend(dtrans,dcod,al,phib0,phi0,psi1,psi2,ak,
     1     utwiss,dx,dy,theta,iv)
            use tfstk
      use ffs
      use tffitcode
      use mathfun
      implicit none
      integer*4 ,intent(in):: iv
      real*8 ,intent(in):: utwiss(ntwissfun),
     $     al,phib0,phi0,psi1,psi2,ak,dx,dy,theta
      real*8 phi,phib,rhob,rho0,aind,ax,ay,pr,
     $     rho,tanp1,tanp2,ak1,ak2,xi,pxi,yi,pyi,rhos,rhox,phix,
     $     rhoa,rhoy,phiy,bpr,dphix,csphix,snphix,sinsq,scphix,
     $     xcs,drho,dx1,dpx1,xd,t11,t12,t21,t16,t33,t34,b11,b12,
     $     b21,b16,b26,a11,a12,a16,a21,a22,a26,a33,a34,a36,a43,
     $     a44,a46,chphix,chphiy,csphiy,ddrho,dk1,dk2,dphiy,
     $     scphiy,u,shphix,shphiy,ycs,snphiy,x(5),y(5)
      real*8 ,intent(out):: dtrans(4,5),dcod(6)
      real*8 cod(6)
      if(al .eq. 0.d0)then
c        write(*,*)' QDBTHI Not installed yet.'
        dtrans=0.d0
        dcod=0.d0
        return
      endif
      if(phi0 .eq. 0.d0 .or. phib0 .eq. 0.d0)then
        phi=1.d-7
        phib=1.d-7
      else
        phi=phi0
        phib=phib0
      endif
      dx1=0.d0
      dpx1=0.d0
      a12=0.d0
      a16=0.d0
      a22=0.d0
      a26=0.d0
      a34=0.d0
      a44=0.d0
      a46=0.d0
      a36=0.d0
      call qtentu(dtrans,cod,utwiss,.true.)
      call qchg(dtrans,cod,-dx,-dy,theta,.true.)
      rhob=al/phib
      rho0=al/phi
      aind=ak*rho0/phi
      if(aind .eq. -1.d0)then
        aind=-1.d0+1.d-8
      endif
      ax=sqrt(abs(1.d0+aind))
      ay=sqrt(abs(aind))
      pr=1.d0+utwiss(mfitddp)
      rho=rhob*pr
      tanp1=tan(psi1*phi0)
      tanp2=tan(psi2*phi0)
      ak1=-tanp1/rho
      ak2=-tanp2/rho
      xi =cod(1)
      pxi=cod(2)-ak1*xi
      yi =cod(3)
      pyi=cod(4)+ak1*yi
      rhos=sqrt(rho*rho0)
      if(ax .ne. 0.d0)then
        rhox=rhos/ax
        phix=al/rhox
        rhoa=rho/(1.d0+aind)
      else
        phix=0.d0
        rhox=1.d19
        rhoa=1.d19
      endif
      if(ay .ne. 0.d0)then
        rhoy=rhos/ay
        phiy=al/rhoy
      else
        phiy=0.d0
        rhoy=1.d19
      endif
      if(iv .eq. 2)then
        bpr=rho/rho0
        if(aind .gt. -1.d0)then
          dphix=phib/pr/phix
          csphix=cos(phix)
          snphix=sin(phix)
          if(csphix .ge. 0.d0)then
            sinsq=snphix**2/(1.d0+csphix)
          else
            sinsq=1.d0-csphix
          endif
c          sinsq=2.d0*sin(phix*.5d0)**2
          scphix=sinc(phix)
          xcs =phix*csphix+snphix
          drho=(rho-rho0)/(1.d0+aind)
          dx1 =drho*((aind-1.d0)*sinsq/phix+snphix)*dphix
          dpx1=drho*(xcs+(aind-1.d0)*snphix)*dphix/al
          xd=xi-drho
          t11=csphix
          t12=rhox*snphix
          t21=-snphix/rhox
          t16=((phix*snphix*xd
     1         -rhox*scphix*pxi)*.5d0
     1         +rhoa*sinsq)/pr
          if(aind .gt. 0.d0)then
            t33=cosh(phiy)
            t34=rhoy*sinh(phiy)
          elseif(aind .eq. 0.d0)then
            t33=1.d0
            t34=al
          else
            t33=cos(phiy)
            t34=rhoy*sin(phiy)
          endif
          b11 =-snphix*dphix
          b12 =rhox/phix*scphix*dphix
          b21 =-xcs*dphix/al
          b16=( drho*(1.d0-aind)*snphix
     1         +2.d0*rhoa*((aind-1.d0)/phix*sinsq+snphix)
     1         +xd*xcs
     1         +pxi*rhox/phix*(scphix+phix**2*snphix)
     1        )*dphix*.5d0/pr
          b26=( pr*rhoa*(aind-1.d0)*snphix
     1         +(pr*(rhoa+pxi*al*.5d0)+drho*(1.d0-aind)*.5d0)*xcs
     1         +xd*(3.d0*phix*csphix+(1.d0-phix**2)*snphix)*.5d0
     1        )*dphix/al/pr
        elseif(aind .eq. -1.d0)then
          snphix=phib/pr
          dx1=(bpr-1.d0)/bpr*al
          dpx1=dx1/al
          t11=1.d0
          t12=al
          t21=0.d0
          t16=0.d0
          t33=cos(phiy)
          t34=rhoy*sin(phiy)
          b11=-snphix
          b12=0.d0
          b21=0.d0
          b16=-2.d0*bpr/pr*phi*al
          b26=b16/al
        else
          dphix=phib/pr/phix
          chphix=cosh(phix)
          shphix=sinh(phix)
          sinsq=shphix**2/(1.d0+chphix)
c          sinsq=2.d0*sin(phix*.5d0)**2
          scphix=sinhc(phix)
          xcs=phix*chphix+shphix
          drho=(rho-rho0)/(1.d0+aind)
          dx1=drho*((aind-1.d0)*sinsq/phix+shphix)*dphix
          dpx1=drho*(xcs+(aind-1.d0)*shphix)*dphix/al
          xd=xi-drho
          t11=chphix
          t12=rhox*shphix
          t21=shphix/rhox
          t16=((-phix*shphix*xd
     1         -rhox*scphix*pxi)*.5d0
     1         -rhoa*sinsq)/pr
          t33=cos(phiy)
          t34=rhoy*sin(phiy)
          b11=-shphix*dphix
          b12=-rhox/phix*scphix*dphix
          b21=-xcs*dphix/al
          b16=( drho*(1.d0-aind)*shphix
     1         +2.d0*rhoa*((aind-1.d0)/phix*sinsq+shphix)
     1         +xd*xcs
     1         -pxi*rhox/phix*(scphix-phix**2*shphix)
     1        )*dphix*.5d0/pr
          b26=( pr*rhoa*(aind-1.d0)*shphix
     1         +(pr*(rhoa+pxi*al*.5d0)+drho*(1.d0-aind)*.5d0)*xcs
     1         +xd*(3.d0*phix*chphix+(1.d0+phix**2)*shphix)*.5d0
     1        )*dphix/al/pr
        endif
        dk1=2.d0*psi1/rho
        dk2=2.d0*psi2/rho
        a11=b11+dk1*t12
        a12=b12
        a21=b21+(dk1+dk2)*t11
        a22=b11+dk2*t12
        a33=-dk1*t34
        a34=0.d0
        a43=-t33*(dk1+dk2)
        a44=-dk2*t34
        a16=b16
        a26=b26+dk2*t16
        a36=0.d0
        a46=0.d0
      elseif(iv .eq. 8)then
        bpr=rho/rho0
        if(aind .gt. 0.d0)then
          dphix= al/2.d0/phix/bpr
          dphiy=-al/2.d0/phiy/bpr
          csphix=cos(phix)
          snphix=sin(phix)
          if(csphix .ge. 0.d0)then
            sinsq=snphix**2/(1.d0+csphix)
          else
            sinsq=1.d0-csphix
          endif
c          sinsq=2.d0*sin(phix*.5d0)**2
          chphiy=cosh(phiy)
          shphiy=sinh(phiy)
          scphix=sinc(phix)
          scphiy=sinhc(phiy)
          xcs=phix*csphix+snphix
          ycs=phiy*chphiy+shphiy
          drho=(rho-rho0)/(1.d0+aind)
          ddrho=-rho0/phi*drho/(1.d0+aind)
          dx1=ddrho*sinsq+drho*snphix*dphix
          dpx1=ddrho*snphix/rhox+drho/al*xcs*dphix
          a11=-snphix*dphix
          a12= rhox/phix*scphix*dphix
          a21=-xcs/al*dphix
          a16=( drho*snphix
     1         +rhoa*(-2.d0*sinsq/phix+snphix)
     1         -(drho-xi)*xcs*.5d0
     1         +pxi*rhox/phix*(scphix+phix**2*snphix)*.5d0
     1        )*dphix/pr
          a26=(-rhoa*2.d0*snphix
     1         +(rhoa+pxi*al*.5d0+drho)*xcs
     1         -(drho-xi)*.5d0*(3.d0*phix*csphix+(1.d0-phix**2)*snphix)
     1        )/al/pr*dphix
          a33=-shphiy*dphiy
          a34=-rhoy/phiy*scphiy*dphiy
          a43=-(shphiy+phiy*chphiy)/al*dphiy
          a36=(yi*ycs
     1         -pyi*rhoy/phiy*(scphiy-phiy**2*shphiy)
     1        )*.5d0/pr*dphiy
          a46=(pyi*ycs
     1         -yi*rhoy/phiy*(3.d0*phiy*chphiy+(1.d0+phiy**2)*shphiy)
     1        )*.5d0/pr*dphiy
        elseif(aind .eq. 0.d0)then
          u=phi**2/al
          dphix= phix/2.d0/u
          csphix=cos(phix)
          snphix=sin(phix)
          if(csphix .ge. 0.d0)then
            sinsq=snphix**2/(1.d0+csphix)
          else
            sinsq=1.d0-csphix
          endif
c          sinsq=2.d0*sin(phix*.5d0)**2
          scphix=sinc(phix)
          drho=(rho-rho0)
          ddrho=-rho0/phi*drho
          dx1=ddrho*sinsq+drho*snphix*dphix
          dpx1=ddrho*snphix/rhox+drho/al*(snphix+phix*csphix)*dphix
          a11=-snphix*dphix
          a12= rhox/phix*scphix*dphix
          a21=-(snphix+phix*csphix)/al*dphix
          a16=((-2.d0*sinsq/phix+phix*(rho0-rho)/rho*csphix/2.d0
     1         +(3.d0*rho-rho0)/rho*snphix/2.d0)*rhox**2/rho0
     1         +xi/2.d0*(phix*csphix+snphix)
     1         +pxi*rhox/2.d0*(scphix/phix+phix*snphix)
     1        )*dphix/pr
          a26=((scphix/phix*(1.d0+rho0/rho)+(1.d0-rho0/rho)*phix*snphix
     1         )*rhox/rho0
     1         +xi*(3.d0*csphix+(1.d0/phix-phix)*snphix)/rhox
     1         +pxi*(phix*csphix+snphix))*dphix/2.d0/pr
          a33=-rho0*al/rho/2.d0
          a34= rho0/rho/6.d0*al**2
          a43=-rho0/rho
          a36= (yi*2.d0
     1          +pyi*al/1.5d0)*rho0/rho/4.d0/pr
          a46=-(yi*4.d0/al
     1          +pyi*2.d0)*rho0/rho/4.d0/pr
        elseif(aind .gt. -1.d0)then
          dphix= al/2.d0/phix/bpr
          dphiy=-al/2.d0/phiy/bpr
c          sinsq=2.d0*sin(phix*.5d0)**2
          csphix=cos(phix)
          snphix=sin(phix)
          if(csphix .ge. 0.d0)then
            sinsq=snphix**2/(1.d0+csphix)
          else
            sinsq=1.d0-csphix
          endif
c          sinsq=2.d0*sin(phix*.5d0)**2
          csphiy=cos(phiy)
          snphiy=sin(phiy)
          scphix=sinc(phix)
          scphiy=sinc(phiy)
          xcs=phix*csphix+snphix
          ycs=phiy*csphiy+snphiy
          drho=(rho-rho0)/(1.d0+aind)
          ddrho=-rho0/phi*drho/(1.d0+aind)
          dx1=ddrho*sinsq+drho*snphix*dphix
          dpx1=ddrho*snphix/rhox+drho/al*xcs*dphix
          a11=-snphix*dphix
          a12= rhox/phix*scphix*dphix
          a21=-xcs/al*dphix
          a16=( drho*snphix
     1         +rhoa*(-2.d0*sinsq/phix+snphix)
     1         -(drho-xi)*xcs*.5d0
     1         +pxi*rhox/phix*(scphix+phix**2*snphix)*.5d0
     1        )*dphix/pr
          a26=(-rhoa*2.d0*snphix
     1         +(rhoa+pxi*al*.5d0+drho)*xcs
     1         -(drho-xi)*.5d0*(3.d0*phix*csphix+(1.d0-phix**2)*snphix)
     1        )/al/pr*dphix
          a33=-snphiy*dphiy
          a34= rhoy/phiy*scphiy*dphiy
          a43=-(snphiy+phiy*csphiy)/al*dphiy
          a36=(yi*ycs
     1         +pyi*rhoy/phiy*(scphiy+phiy**2*snphiy)
     1        )*.5d0/pr*dphiy
          a46=(pyi*ycs
     1         +yi*rhoy/phiy*(3.d0*phiy*csphiy+(1.d0-phiy**2)*snphiy)
     1        )*.5d0/pr*dphiy
        else
          dphix= al/2.d0/phix/bpr
          dphiy=-al/2.d0/phiy/bpr
          chphix=cosh(phix)
          shphix=sinh(phix)
          sinsq=shphix**2/(1.d0+chphix)
c          sinsq=2.d0*sinh(phix*.5d0)**2
          csphiy=cos(phiy)
          snphiy=sin(phiy)
          scphix=sinhc(phix)
          scphiy=sinc(phiy)
          xcs=phix*chphix+shphix
          ycs=phiy*csphiy+snphiy
          drho=(rho-rho0)/(1.d0+aind)
          ddrho=-rho0/phi*drho/(1.d0+aind)
          dx1=-ddrho*sinsq+drho*shphix*dphix
          dpx1=-ddrho*shphix/rhox+drho/al*xcs*dphix
          a11=-shphix*dphix
          a12=-rhox/phix*scphix*dphix
          a21=-xcs/al*dphix
          a16=( drho*shphix
     1         +rhoa*(-2.d0*sinsq/phix+shphix)
     1         -(drho-xi)*xcs*.5d0
     1         -pxi*rhox/phix*(scphix-phix**2*shphix)*.5d0
     1        )*dphix/pr
          a26=(-rhoa*2.d0*shphix
     1         +(rhoa+pxi*al*.5d0+drho)*xcs
     1         -(drho-xi)*.5d0*(3.d0*phix*chphix+(1.d0+phix**2)*shphix)
     1        )/al/pr*dphix
          a33=-snphiy*dphiy
          a34= rhoy/phiy*scphiy*dphiy
          a43=-(snphiy+phiy*csphiy)/al*dphiy
          a36=(yi*ycs
     1         +pyi*rhoy/phiy*(scphiy+phiy**2*snphiy)
     1        )*.5d0/pr*dphiy
          a46=(pyi*ycs
     1         +yi*rhoy/phiy*(3.d0*phiy*csphiy+(1.d0-phiy**2)*snphiy)
     1        )*.5d0/pr*dphiy
        endif
        a22=a11
        a44=a33
      else
c     iv != kytbl(kwANGL, icBEND) && iv != kytbl(kwK1, icBEND)
c     This block is never invoked by src/qdtwiss.f
        write(*,*)'qdbend @ src/qdbend.f: ',
     $        'Invoked with invalid iv=',iv,
     $        '(FIXME)' 
        call abort
      endif
c      do 4010 i=1,5
        dtrans(2,:)=dtrans(2,:)-ak1*dtrans(1,:)
        dtrans(4,:)=dtrans(4,:)+ak1*dtrans(3,:)
c4010  continue
c      do 10 i=1,5
        x=dtrans(1,:)
        dtrans(1,:)= a11*x+a12*dtrans(2,:)
        dtrans(2,:)= a21*x+a22*dtrans(2,:)
        y=dtrans(3,:)
        dtrans(3,:)= a33*y+a34*dtrans(4,:)
        dtrans(4,:)= a43*y+a44*dtrans(4,:)
c10    continue
      dtrans(1,5)=dtrans(1,5)+a16
      dtrans(2,5)=dtrans(2,5)+a26
      dtrans(3,5)=dtrans(3,5)+a36
      dtrans(4,5)=dtrans(4,5)+a46
c     write(*,*)a16,a26,dtrans(1,5),dtrans(2,5)
      dcod(1)= a11*xi+a12*pxi+dx1
      dcod(2)= a21*xi+a22*pxi+dpx1
      dcod(3)= a33*yi+a34*pyi
      dcod(4)= a43*yi+a44*pyi
      dcod(6)=pr
c      do 4020 i=1,5
        dtrans(2,:)=dtrans(2,:)-ak2*dtrans(1,:)
        dtrans(4,:)=dtrans(4,:)+ak2*dtrans(3,:)
c4020  continue
      dcod(2)=dcod(2)-ak2*dcod(1)
      dcod(4)=dcod(4)+ak2*dcod(3)
      call qchg(dtrans,dcod,0.d0,0.d0,-theta,.false.)
      return
      end
