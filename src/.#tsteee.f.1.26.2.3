      subroutine tsteee(trans,cod,beam,al,phib,dx,dy,theta,enarad,
     $     apsi1,apsi2,fb1,fb2,mfring,fringe,ld)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      real*8 epslon,a3,a5,a7,a9,a11,a13,a15
      parameter (epslon=1.d-6)
      parameter (a3=1.d0/6.d0,a5=3.d0/40.d0,a7=5.d0/112.d0,
     1           a9=35.d0/1152.d0,a11=63.d0/2816.d0,
     1           a13=231.d0/13312.d0,a15=143.d0/10240.d0)
      integer*4 mfring,ld,nrad,ndiv,n
      real*8 trans(6,12),cod(6),beam(42),al,phib,dx,dy,theta,
     $     fb1,fb2,
     $     rhob,
     $     dxfr1,dyfr1,
     $     dxfr2,dyfr2,
     $     f1r,b,aln,phin,pr,pxi,pyi,rhoe,s,dpz1,
     $     pz1,dpx,pxf,d,w,u,spz,spx,phsq,dl,dpz2,pz2,
     $     dyfra1,dyfra2,apsi1,apsi2,tanp1,tanp2,als
      real*8 trans1(6,6)
      logical*4 enarad,fringe
      if(al .le. 0.d0)then
        call tthine(trans,cod,beam,2,al,-phib,dx,dy,theta,
     $       .false.,ld)
        return
      elseif(phib .eq. 0.d0)then
        call tdrife(trans,cod,beam,al,
     $         0.d0,0.d0,0.d0,.true.,enarad,calpol,irad,ld)
        return
      endif
      call tchge(trans,cod,beam,-dx,-dy,theta,.true.,ld)
      rhob=al/phib
      if(fb1 .ne. 0.d0 .and. (mfring .gt. 0 .or. mfring .eq. -1))then
        dxfr1=-fb1**2/rhob/24.d0
        dyfr1=fb1/rhob**2/6.d0
        if(fringe)then
          dyfra1=4.d0*dyfr1/fb1**2
        else
          dyfra1=0.d0
        endif
        call tblfre(trans,cod,beam,dxfr1,dyfr1,dyfra1,ld)
        f1r=fb1
      else
        f1r=0.d0
      endif
      if(enarad)then
        tanp1=tan(apsi1)
        tanp2=tan(apsi2)
        b=-brhoz/rhob
        nrad=int(al/epsrad*crad*(h0*b)**2)
        ndiv=1+max(int(nrad*emidiv*emidib),
     1       int(abs(phib)/epsrad/1.d3*emidiv*emidib))
      else
c     begin initialize for preventing compiler warning
        tanp1=0.d0
        tanp2=0.d0
c     end   initialize for preventing compiler warning
        ndiv=1
      endif
      aln=al/ndiv
      phin=phib/ndiv
      if(mfring .ne. -2)then
c        call tbfrie(trans,cod,beam,-rhob,0.d0,.true.,ld)
        call tbedge(trans,cod,beam,al,-phib,apsi1,.true.,ld)
      endif
      call tinitr(trans1)
      als=0.d0
      do 100 n=1,ndiv
        if(enarad)then
          if(n .eq. 1)then
            call trade(trans,beam,cod,0.d0,b,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,-tanp1,
     $           f1r,f1r,0.d0,al,.5d0*aln)
          else
            call trade(trans,beam,cod,0.d0,b,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,0.d0,
     $           f1r,f1r,als,al,aln)
          endif
          als=als+aln
        endif
        pr=1.d0+cod(6)
        pxi=cod(2)/pr
        pyi=cod(4)/pr
        rhoe=rhob*pr
        s=min(.9d0,pxi**2+pyi**2)
        dpz1=-s/(1.d0+sqrt(1.d0-s))
        pz1=1.d0+dpz1
        dpx=aln/rhoe
        pxf=pxi+dpx
        s=min(.9d0,pxf**2+pyi**2)
        dpz2=-s/(1.d0+sqrt(1.d0-s))
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
        trans1(1,2)=(2.d0*spz+spx*(pxf/pz2+pxi/pz1))/spz**2*aln/pr
        trans1(1,6)=-spx/pz1/pz2/spz*aln/pr
        trans1(1,4)=-pyi*trans1(1,6)
        trans1(5,2)=trans1(1,6)
        trans1(5,6)= rhob*u/pz1/pz2
        trans1(5,4)=-pyi*trans1(5,6)
        trans1(3,2)=trans1(1,4)
        trans1(3,4)=(aln+dl)/pr-pyi*trans1(5,4)
        trans1(3,6)=trans1(5,4)
        trans1(5,6)=trans1(5,6)-(aln+dl)/pr+h0/h1emit**3*aln
        call tmultr5(trans,trans1,irad)
        call tmulbs(beam ,trans1,.true.,.true.)
        if(calpol)then
          if(n .eq. 1)then
            call polpar(21,ld,aln,0.d0,-phin,0.d0,0.d0,cod)
          else
            call polpar(22,ld,aln,0.d0,-phin,0.d0,0.d0,cod)
          endif
        endif
        cod(1)=cod(1)+aln*spx/spz
        cod(2)=pxf*pr
        cod(3)=cod(3)+pyi*(aln+dl)
        cod(5)=cod(5)-(dl+dvemit*aln)
        if(n .eq. ndiv)then
          if(calpol)then
            if(irad .eq. 6)then
              npelm=npelm+1
            else
              ipelm=ipelm+1
              call tinitr(trans1)
              call tmov(trans1,rlist(ipoltr+(ipelm-1)*36),36)
              call polpar(23,ld,aln,0.d0,-phin,0.d0,0.d0,cod)
            endif
          endif
          if(enarad)then
            if(fb2 .ne. 0.d0 .and.
     $           (mfring .gt. 0 .or. mfring .eq. -2))then
              f1r=fb2
            endif
            call trade(trans,beam,cod,0.d0,b,0.d0,0.d0,
     $           0.d0,0.d0,0.d0,-tanp2,
     $           f1r,f1r,al,al,.5d0*aln)
          endif
        endif
100   continue
      if(mfring .ne. -1)then
c        call tbfrie(trans,cod,beam, rhob,0.d0,.false.,ld)
        call tbedge(trans,cod,beam,al,-phib,apsi2,.false.,ld)
      endif
      if(fb2 .ne. 0.d0 .and. (mfring .gt. 0 .or. mfring .eq. -2))then
        dxfr2=fb2**2/rhob/24.d0
        dyfr2=fb2/rhob**2/6.d0
        if(fringe)then
          dyfra2=4.d0*dyfr2/fb2**2
        else
          dyfra2=0.d0
        endif
        call tblfre(trans,cod,beam,dxfr2,dyfr2,dyfra2,ld)
      endif
      call tchge(trans,cod,beam,dx,dy,-theta,.false.,ld)
      return
      end
