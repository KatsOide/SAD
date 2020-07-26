      subroutine tsteee(trans,cod,beam,srot,al,phib,dx,dy,theta,enarad,
     $     apsi1,apsi2,fb1,fb2,mfring,fringe,next)
      use ffs_flag
      use tmacro
      use bendib, only:rbh,rbl,tbendal
      use temw, only:tsetr0,tmulbs
      use kradlib, only:tradke      
      use mathfun
      implicit none
      real*8 epslon,a3,a5,a7,a9,a11,a13,a15
      parameter (epslon=1.d-6)
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      integer*4 mfring,nrad,ndiv,n,n1,n2
      real*8 trans(6,12),cod(6),beam(42),srot(3,9),
     $     al,phib,dx,dy,theta,f1r,f2r,
     $     fb1,fb2,rhob,rbc,alc,phic,alx,alr,
     $     dxfr1,dyfr1,
     $     dxfr2,dyfr2,
     $     b,aln,phin,pr,pxi,pyi,rhoe,s,dpz1,
     $     pz1,dpx,pxf,d,w,u,spz,spx,phsq,dl,dpz2,pz2,
     $     dyfra1,dyfra2,apsi1,apsi2,tanp1,tanp2
      real*8 trans1(6,6)
      logical*4 enarad,fringe,next,prev,krad
      if(al .eq. 0.d0)then
        call tthine(trans,cod,beam,srot,2,al,-phib,dx,dy,theta,
     $       .false.)
        return
      elseif(phib .eq. 0.d0)then
        call tdrife(trans,cod,beam,srot,al,
     $         0.d0,0.d0,0.d0,0.d0,.true.,.false.,irad)
        return
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,theta,0.d0,0.d0,.true.)
      krad=enarad .and. al .ne. 0.d0
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      rhob=al/phib
      prev=bradprev .ne. 0.d0
      rbc=1.d0
      f1r=0.d0
      f2r=0.d0
      if(fb1 .ne. 0.d0 .and. (mfring .gt. 0 .or. mfring .eq. -1))then
        dxfr1=-fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        if(fringe)then
          dyfra1=4.d0*dyfr1/fb1**2
        else
          dyfra1=0.d0
        endif
        call tblfre(trans,cod,beam,dxfr1,dyfr1,dyfra1)
        f1r=.5d0*fb1
        rbc=rbc-f1r/al
        n1=-1
      else
        n1=1
      endif
      if(fb2 .ne. 0.d0 .and.
     $       (mfring .gt. 0 .or. mfring .eq. -2))then
        f2r=0.5d0*fb2
        rbc=rbc-f2r/al
        n2=2
      else
        n2=0
      endif
      alc=al*rbc
      phic=phib*rbc
      if(enarad)then
        tanp1=tan(apsi1)
        tanp2=tan(apsi2)
        b=-brhoz/rhob
        nrad=int(abs(alc/epsrad*crad*(h0*b)**2))
        ndiv=1+max(int(nrad*emidiv*emidib),
     1       int(abs(phib*h0*anrad)/epsrad/1.d6*emidiv*emidib))
c        write(*,*)'tsteee ',ndiv,nrad,phib,b,epsrad,crad
      else
c     begin initialize for preventing compiler warning
        tanp1=0.d0
        tanp2=0.d0
c     end   initialize for preventing compiler warning
        ndiv=1
      endif
      aln=alc/ndiv
      phin=phic/ndiv
      if(mfring .ne. -2)then
c        call tbfrie(trans,cod,beam,-rhob,0.d0,.true.)
        call tbedge(trans,cod,beam,al,-phib,apsi1,.true.)
      endif
      call tinitr(trans1)
      do 100 n=n1,ndiv+n2
        call tbendal(n,ndiv,f1r,f2r,aln,alx,alr)
        pr=1.d0+cod(6)
        pxi=cod(2)/pr
        pyi=cod(4)/pr
        rhoe=rhob*pr
        s=pxi**2+pyi**2
        dpz1=-s/(1.d0+sqrtl(1.d0-s))
        pz1=1.d0+dpz1
        dpx=alx/rhoe
        pxf=pxi+dpx
        s=pxf**2+pyi**2
        dpz2=-s/(1.d0+sqrtl(1.d0-s))
        pz2=1.d0+dpz2
        phsq=pxi**2+pz1**2
        d=pxf*pz1+pxi*pz2
        if(d .eq. 0.d0 .or. pxi*pxf .lt. 0.d0)then
          u=(dpx+pxf*dpz1-pxi*dpz2)/phsq
          w=(pxf*dpz1-pxi*dpz2+dpx*pyi**2)/phsq
        else
          u=dpx*(pxf+pxi)/d
          w=-dpx*(pxf*dpz1+pxi*dpz2)/d
        endif
        s=u**2
        if(s .gt. 2.d-2)then
          dl=rhoe*((asin(u)/u-1.d0)*u+w)
        else
          if(s .gt. 2.d-4)then
            dl=rhoe*(s*(a3+s*(a5+s*(a7+s*(a9+s*(a11+s*(a13+s*a15))))))
     $           *u+w)
          else
            dl=rhoe*(s*(a3+s*(a5+s*(a7+s*a9)))*u+w)
          endif
        endif
        spz=pz2+pz1
        spx=pxf+pxi
        phsq=(1.d0-pyi)*(1.d0+pyi)
        trans1(1,2)=(2.d0*spz+spx*(pxf/pz2+pxi/pz1))/spz**2*alx/pr
        trans1(1,6)=-spx/pz1/pz2/spz*alx/pr
        trans1(1,4)=-pyi*trans1(1,6)
        trans1(5,2)=trans1(1,6)
        trans1(5,6)= rhob*u/pz1/pz2
        trans1(5,4)=-pyi*trans1(5,6)
        trans1(3,2)=trans1(1,4)
        trans1(3,4)=(alx+dl)/pr-pyi*trans1(5,4)
        trans1(3,6)=trans1(5,4)
        trans1(5,6)=trans1(5,6)-(alx+dl)/pr+h0/h1emit**3*alx
        call tmultr5(trans,trans1,irad)
        call tmulbs(beam ,trans1,.true.)
        cod(1)=cod(1)+alx*spx/spz
        cod(2)=pxf*pr
        cod(3)=cod(3)+pyi*(alx+dl)
        cod(5)=cod(5)-(dl+dvemit*alx)
        if(n .ne. ndiv+n2 .and. krad)then
          call tradke(trans,cod,beam,srot,alr,0.d0,0.d0)
        endif
100   continue
      if(.not. next)then
        bradprev=0.d0
      endif
      if(mfring .ne. -1)then
c        call tbfrie(trans,cod,beam, rhob,0.d0,.false.)
        call tbedge(trans,cod,beam,al,-phib,apsi2,.false.)
      endif
      if(f2r .ne. 0.d0)then
        dxfr2=fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        if(fringe)then
          dyfra2=4.d0*dyfr2/fb2**2
        else
          dyfra2=0.d0
        endif
        call tblfre(trans,cod,beam,dxfr2,dyfr2,dyfra2)
      endif
      if(krad)then
        call tradke(trans,cod,beam,srot,alr,0.d0,0.d0)
      endif
      call tchge(trans,cod,beam,srot,
     $     -dx,-dy,-theta,0.d0,0.d0,.false.)
      return
      end
