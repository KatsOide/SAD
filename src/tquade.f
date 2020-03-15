      subroutine tquade(trans,cod,beam,srot,al,ak,
     1     dx,dy,theta,enarad,
     $     fringe,f1in,f2in,f1out,f2out,mfring,eps0,
     $     kin,achro,next)
      use tfstk
      use ffs_flag
      use tmacro
      use temw
      use tspin, only:tradke
      use mathfun
      implicit none
      integer*4 ndiv,i,mfring,n,itgetqraddiv
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6),srot(3,9),
     $     al,ak,dx,dy,theta,f1in,f2in,f1out,f2out,
     $     eps0,f1r,f2r,eps,b1,akn,
     $     aln,pr,akk,phi,scphi,shcphi,sinc2,sinhc2,akr,
     $     xsin2,xsinh2,a11,a12,a21,b11,b12,b21,als,
     $     xi,pxi,yi,pyi,xf,pxf,yf,pyf,
     $     zx,zy,zxp,zyp,x,y,bm21
      logical*4 enarad,fringe,kin,next,prev,achro,krad
      real*8 , parameter:: pramin=1.d-4
      integer*4 , parameter :: ndivmax=512
      if(al .eq. 0.d0)then
        call tthine(trans,cod,beam,srot,4,
     $       al,ak,dx,dy,theta,.false.)
        return
      elseif(ak .eq. 0.d0)then
        call tdrife(trans,cod,beam,srot,al,
     $       0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
        return
      endif
      bm21=beam(21)
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,0.d0,0.d0,.true.)
      krad=enarad .and. al .ne. 0.d0
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      if(fringe .and. mfring .ge. 0. and. mfring .ne. 2)then
        call tqfrie(trans,cod,beam,ak,al,0.d0)
      endif
      prev=bradprev .ne. 0.d0
      if(mfring .eq. 1 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak,f1in,f2in,0.d0)
        f1r=f1in
      else
        f1r=0.d0
      endif
      if(krad .and. f1in .ne. 0.d0)then
        call tradke(trans,cod,beam,srot,f1in,0.d0,0.d0)
      else
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        f2r=f1out
      else
        f2r=0.d0
      endif
      if(eps0 .eq. 0.d0)then
        eps=.1d0
      else
        eps=.1d0*eps0
      endif
      b1=0.d0
      ndiv=1+min(10000,int(abs(ak*al)/eps))
      if(enarad)then
        b1=brhoz*ak/al
        ndiv=min(ndivmax,max(ndiv,itgetqraddiv(cod,ak,al,0.d0)))
      endif
      akn=ak/ndiv
      aln=al/ndiv
      if(irad .gt. 6)then
        call tinitr(trans1)
      endif
      if(achro)then
        pr=1.d0
      else
        pr=max(pramin,1.d0+cod(6))
      endif
      akr=ak/al/pr
      akk=sqrt(abs(akr))
      phi=aln*akk
      scphi=sinc(phi)
      shcphi=sinhc(phi)
      sinc2=sinc(2.d0*phi)
      sinhc2=sinhc(2.d0*phi)
      xsin2=xsin(2.d0*phi)
      xsinh2=xsinh(2.d0*phi)
      if(ak .gt. 0.d0)then
        a11=cos(phi)
        a12=sin(phi)/akk
        a21=-a12*akk**2
        b11=cosh(phi)
        b12=sinh(phi)/akk
        b21=b12*akk**2
      else
        b11=cos(phi)
        b12=sin(phi)/akk
        b21=-b12*akk**2
        a11=cosh(phi)
        a12=sinh(phi)/akk
        a21=a12*akk**2
      endif
      als=0.d0
      do 100 n=1,ndiv
c        if(enarad)then
c          bx= b1*cod(3)
c          by= b1*cod(1)
c          bxy= b1
c          if(n .eq. 1)then
c            call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $           0.d0,bxy,0.d0,0.d0,
c     $           .5d0*aln,als,al,f1r,f2r,prev,next)
c          else
c            call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $           0.d0,bxy,0.d0,0.d0,
c     $           aln,als,al,f1r,f2r,prev,next)
c          endif
c          If(abs(cod(2)) .gt. 1.d0)then
c            write(*,*)'tquade-trade-2 ',cod(2),cod(5),cod(6)
c          endif
c          als=als+aln
        if(kin)then
          if(n .eq. 1)then
            call tqente(trans,cod,beam,aln*.5d0,0.d0,irad)
          else
            call tqente(trans,cod,beam,aln     ,0.d0,irad)
          endif
        endif
        xi =cod(1)
        pxi=cod(2)/pr
        yi =cod(3)
        pyi=cod(4)/pr
        if(akr .gt. 0.d0)then
          xf= a11*xi+a12*pxi
          pxf=a21*xi+a11*pxi
          yf= b11*yi+b12*pyi
          pyf=b21*yi+b11*pyi
          zx=pxi**2+(xi*akk)**2
          zy=(pyi-yi*akk)*(pyi+yi*akk)
          trans1(5,1)=(-.5d0*akk*xsin2*xi-a21*a12*pxi)*.5d0
          trans1(5,3)=(.5d0*akk*xsinh2*yi-b21*b12*pyi)*.5d0
          trans1(2,6)=(akk*scphi*xi-aln*a21*pxi)*.5d0
          trans1(4,6)=(-akk*shcphi*yi-aln*b21*pyi)*.5d0
          trans1(1,6)=-(aln*pxf+a12*pxi)*.5d0/pr
          trans1(3,6)=-(aln*pyf+b12*pyi)*.5d0/pr
          zxp=-(2.d0*pxi**2+(akk*xi)**2)/pr
          zyp=-(2.d0*pyi**2-(akk*yi)**2)/pr
          trans1(5,6)=-.25d0*(zxp*(aln+a11*a12)+zyp*(aln+b11*b12)
     1           -.25*(zx*sinc2+zy*sinhc2)/akk/pr)
     1           -.5*xi*(trans1(1,6)*a21-.5d0*xf*(a21-akk*phi*a11)/pr)
     1           -.5*yi*(trans1(3,6)*b21-.5d0*yf*(b21+akk*phi*b11)/pr)
     1           +h0/h1emit**3*aln
        else
          xf= a11*xi+a12*pxi
          pxf=a21*xi+a11*pxi
          yf= b11*yi+b12*pyi
          pyf=b21*yi+b11*pyi
          zx=(pxi-xi*akk)*(pxi+xi*akk)
          zy=pyi**2+(yi*akk)**2
          trans1(5,1)=( .5d0*akk*xsinh2*xi-a21*a12*pxi)*.5d0
          trans1(5,3)=(-.5d0*akk*xsin2 *yi-b21*b12*pyi)*.5d0
          trans1(2,6)=(-akk*shcphi*xi-aln*a21*pxi)*.5d0
          trans1(4,6)=(akk*scphi*yi-aln*b21*pyi)*.5d0
          trans1(1,6)=-(aln*pxf+a12*pxi)*.5d0/pr
          trans1(3,6)=-(aln*pyf+b12*pyi)*.5d0/pr
          zxp=-(2.d0*pxi**2-(akk*xi)**2)/pr
          zyp=-(2.d0*pyi**2+(akk*yi)**2)/pr
          trans1(5,6)=-.25d0*(zxp*(aln+a11*a12)+zyp*(aln+b11*b12)
     1           -.25d0*(zx*sinhc2+zy*sinc2)/akk/pr)
     1           -.5*xi*(trans1(1,6)*a21-.5d0*xf*(a21+akk*phi*a11)/pr)
     1           -.5*yi*(trans1(3,6)*b21-.5d0*yf*(b21-akk*phi*b11)/pr)
     1           +h0/h1emit**3*aln
        endif
        trans1(1,1)=a11
        trans1(2,1)=pr*a21
        trans1(1,2)=a12/pr
        trans1(2,2)=a11
        trans1(5,2)=-(aln*pxi+a12*pxf)*.5d0/pr
        trans1(3,3)=b11
        trans1(4,3)=pr*b21
        trans1(3,4)=b12/pr
        trans1(4,4)=b11
        trans1(5,4)=-(aln*pyi+b12*pyf)*.5d0/pr
        do i=1,irad
          x=trans(1,i)
          y=trans(3,i)
          trans(5,i)=trans(5,i)
     $         +trans1(5,1)*x+trans1(5,2)*trans(2,i)
     $         +trans1(5,3)*y+trans1(5,4)*trans(4,i)
     $         +trans1(5,6)*trans(6,i)
          trans(1,i)=trans1(1,1)*x+trans1(1,2)*trans(2,i)
     $         +trans1(1,6)*trans(6,i)
          trans(2,i)=trans1(2,1)*x+trans1(2,2)*trans(2,i)
     $         +trans1(2,6)*trans(6,i)
          trans(3,i)=trans1(3,3)*y+trans1(3,4)*trans(4,i)
     $         +trans1(3,6)*trans(6,i)
          trans(4,i)=trans1(4,3)*y+trans1(4,4)*trans(4,i)
     $         +trans1(4,6)*trans(6,i)
        enddo
        call tmulbs(beam ,trans1,.true.)
        cod(1)= xf
        cod(3)= yf
        cod(5)=cod(5)-dvemit*aln-
     1     ((zx*(aln+a12*a11)+
     1       zy*(aln+b12*b11))*.5d0+
     1      xi*xf*a21+yi*yf*b21)*.5d0
        cod(2)=pxf*pr
        cod(4)=pyf*pr
        if(krad)then
          if(n .eq. 1)then
            bsir0=ak/aln*xi*yi
          endif
          call tradke(trans,cod,beam,srot,aln,0.d0,0.d0)
          if(radcod)then
            if(achro)then
              pr=1.d0
            else
              pr=max(pramin,1.d0+cod(6))
            endif
            akr=ak/al/pr
            akk=sqrt(abs(akr))
            phi=akk*aln
            scphi=sinc(phi)
            shcphi=sinhc(phi)
            sinc2=sinc(2.d0*phi)
            sinhc2=sinhc(2.d0*phi)
            xsin2=xsin(2.d0*phi)
            xsinh2=xsinh(2.d0*phi)
            if(akr .gt. 0.d0)then
              a11=cos(phi)
              a12=sin(phi)/akk
              a21=-a12*akk**2
              b11=cosh(phi)
              b12=sinh(phi)/akk
              b21=b12*akk**2
            else
              b11=cos(phi)
              b12=sin(phi)/akk
              b21=-b12*akk**2
              a11=cosh(phi)
              a12=sinh(phi)/akk
              a21=a12*akk**2
            endif
          endif
        endif
100   continue
      if(kin)then
        call tqente(trans,cod,beam,aln*.5d0,0.d0,irad)
      endif
c      if(krad)then
c        bx= b1*cod(3)
c        by= b1*cod(1)
c        bxy= b1
c        call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $       0.d0,bxy,0.d0,0.d0,
c     $       .5d0*aln,al,al,f1r,f2r,prev,next)
c      endif
      if(.not. next)then
        bradprev=0.d0
      endif
      if(mfring .eq. 2 .or. mfring .eq. 3)then
        call tqlfre(trans,cod,beam,al,ak,-f1out,f2out,0.d0)
      endif
      if(fringe .and. mfring .ge. 0 .and. mfring .ne. 1)then
        call tqfrie(trans,cod,beam,-ak,al,0.d0)
      endif
      if(krad .and. f1out .ne. 0.d0)then
        bsir0=-ak/aln*cod(1)*cod(3)
        call tradke(trans,cod,beam,srot,f1out,0.d0,0.d0)
      endif
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,0.d0,0.d0,.false.)
      return
      end
