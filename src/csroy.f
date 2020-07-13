      subroutine tfcsroysetup(isp1,kx,irtc)
      use tfstk
      use iso_c_binding
      implicit none
      type (sad_rlist), pointer :: lazl,lal
      type (sad_descriptor) kx
      integer*8 kawi,kaaj,kaomega,kazjn,kadpjn,karhoj,karho,
     $     kax
      integer*4 isp1,irtc,itfmessage,nr,nj,nphi,mphi
      real*8 dz,aw,r56,r65,zspan,dt
      real*8 ,pointer ::zjn(:,:),dpjn(:,:)
      if(isp .ne. isp1+3)then
        go to 9000
      endif
      if(.not. tfreallistq(ktastk(isp1+1)))then
        go to 9010
      endif
      karho=ktfaddr(ktastk(isp1+1))
      nr=ilist(2,karho-1)
      if(.not. tfreallistq(dtastk(isp1+2),lazl))then
        go to 9010
      endif
      if(.not. tfreallistq(dtastk(isp1+3),lal))then
        go to 9010
      endif
      if(lal%nl .ne. 7)then
        go to 9010
      endif
      dz=lal%rbody(1)
      nj=int(lal%rbody(2))
      nphi=int(lal%rbody(3))
      mphi=nphi*16
      aw=lal%rbody(4)
      r56=lal%rbody(5)
      r65=lal%rbody(6)
      zspan=lal%rbody(7)
      kawi=ktavaloc(0,nr)
      kaaj=ktavaloc(0,nj)
      kaomega=ktavaloc(0,nj)
      kazjn=ktavaloc(0,nj*mphi)
      kadpjn=ktavaloc(0,nj*mphi)
      karhoj=ktavaloc(0,nj)
      dt=min(1.d0,0.01d0/sqrt(abs(r56*r65)))
      call c_f_pointer(c_loc(rlist(kazjn)),zjn,[mphi,nj])
      call c_f_pointer(c_loc(rlist(kazjn)),dpjn,[mphi,nj])
      call csroysetup(
     $     rlist(karho+1:karho+nr),
     $     lazl%rbody(1:nr),
     $     rlist(kawi+1:kawi+nr),
     $     rlist(kaaj+1:kaaj+nj),
     $     rlist(kaomega+1:kaomega+nj),
     $     zjn,dpjn,
     $     rlist(karhoj+1:karhoj+nj),
     $     aw,r56,r65,dz,nr,nj,mphi,zspan,dt,irtc)
      kax=ktadaloc(-1,6)
      klist(kax+1)=ktflist+kawi
      klist(kax+2)=ktflist+kaaj
      klist(kax+3)=ktflist+kaomega
      klist(kax+4)=ktflist+kazjn
      klist(kax+5)=ktflist+kadpjn
      klist(kax+6)=ktflist+karhoj
      kx%k=ktflist+kax
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'CSR::irreg','""')
      endif
      return
 9000 irtc=itfmessage(9,'General::narg','"3"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $'"rho, zl, {dz, nj, nphi, inten, r56, r65, zspan}"')
c
c   rho:     particle density in z, as the output tfcsrhaissin. Sum(rho)*dz = 1.
c   zl:      list of impedance (nr points) in the "real fft" format, as the result of rftr(a)
c   dz:      z mesh size (m)
c   nj:      number of J meshes
c   nphi:    number of phi modes
c   inten:   bunch charge/energy (pC/eV)
c   r56:     one turn matrix component r56 (m)
c   r65:     one turn matrix component r65 (1/m)
c   zspan:   span in Gaussian distribution (sigma).
c
      return
      end

      subroutine tfcsroymat(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_dlist), pointer :: las
      type (sad_rlist), pointer :: lawi,laomega,lazj,larhoj,lal
      type (sad_descriptor) kx
      integer*8 kax
      integer*4 isp1,irtc,itfmessage,nr,nj,nphi,mm
      real*8 dz,pi,bw,anus,sigz,sige
      parameter (pi=3.14159265358979324d0)
      if(isp .ne. isp1+2)then
        go to 9000
      endif
      if(.not. tflistq(dtastk(isp1+1),las))then
        go to 9010
      endif
      if(las%nl .ne. 6)then
        go to 9010
      endif
      if(.not. tfreallistq(las%dbody(1),lawi))then
        go to 9010
      endif
      if(.not. tfreallistq(las%dbody(3),laomega))then
        go to 9010
      endif
      if(.not. tfreallistq(las%dbody(4),lazj))then
        go to 9010
      endif
      if(.not. tfreallistq(las%dbody(6),larhoj))then
        go to 9010
      endif
      nr=lawi%nl
      nj=laomega%nl
      nphi=lazj%nl/nj/16
      if(.not. tfreallistq(dtastk(isp1+2),lal))then
        go to 9010
      endif
      if(lal%nl .ne. 5)then
        go to 9010
      endif
      dz=lal%rbody(1)
      anus=lal%rbody(3)
      sigz=lal%rbody(4)
      sige=lal%rbody(5)
      bw=lal%rbody(2)*4.d0/pi/(2.d0*pi*anus*sige*sigz)
      mm=nj*nphi
      call csroymat(kax,
     $     lawi   %rbody(1),
     $     laomega%rbody(1),
     $     lazj   %rbody(1),
     $     larhoj %rbody(1),
     $     bw,dz,nj,nphi,nr)
      kx%k=ktflist+kax
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"2"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $'"{wi, J, omega, zj, dpj, rhoj}, {dz, inten, nus, sigz, sige}"')
c
c   The list {wi, .., rhoj} is the output of tfcsroysetup.
c   wi:        list of primitive function of w(z)
c   J(nj):     list of actions
c   omega(nj): list of omega(J) corresponding to J(nj)
c   zj(mphi,nj):  list of z coeesponding to J and phi
c   dpj:       list of dp corresponding to J and phi
c   rhoj:      list of (-g' Delta J)^(1/2), Sum(rhoj^2) = 1
c
c   dz:        z mesh size (m)
c   inten:     bunch charge/energy (pC/eV)
c   nus:       synchrotron tune
c   sigz:      0 intensity bunch length (m)
c   sige:      0 intensity momentum spread
c
      return
      end

      subroutine csroymat(kax,wi0,omega,zjn,rhoj,bw,
     $     dz,nj,nphi,nr)
      use tfstk
      implicit none
      integer*8 kax
      integer*4 nj,nphi,nr,mphi,nrh,m,j,j1,m1,n,n1,
     $     mz,jn,jn1,jn10
      real*8 bm(nj*nphi,nj*nphi),wi0(nr),omega(nj),
     $     zjn(nphi*16,nj),ddwi(nr),work(nr),wi(nr),dz,pi,
     $     rhoj(nj),cosmn(nphi*16,nphi),dphi,bw,
     $     wj,wj1,w,a,b,c,d,f,zmj,zmj1,zd
      parameter (pi=3.14159265358979324d0)
      mphi=nphi*16
      dphi=pi/mphi
      nrh=nr/2
      wi(1:nrh)   =bw*wi0(nrh+1:nr)
      wi(nrh+1:nr)=bw*wi0(1:nrh)
      call spline1(nr,wi,ddwi,work,0,0)
      do n=1,nphi
        do m=1,mphi
          cosmn(m,n)=-n*cos(n*pi*(m-0.5)/mphi)
        enddo
      enddo
      bm=0.d0
      do j=1,nj
        wj=rhoj(j)*omega(j)*dphi
        do n=1,nphi
          jn=nphi*(j-1)+n
          do m=1,mphi
            zmj=zjn(m,j)
            do j1=1,nj
              wj1=rhoj(j1)*omega(j1)*dphi
              w=wj*wj1*cosmn(m,n)
              jn10=nphi*(j1-1)
              do m1=1,mphi
                zmj1=zjn(m1,j1)
                zd=(zmj-zmj1)/dz+nrh+1
                mz=int(zd)
                b=zd-mz
                a=1.d0-b
                c=-a*b*(a+1.d0)
                d=-a*b*(b+1.d0)
                f=(a*wi(mz)+b*wi(mz+1)+c*ddwi(mz)+d*ddwi(mz+1))*w
                do n1=1,nphi
                  jn1=jn10+n1
                  bm(jn1,jn)=bm(jn1,jn)+f*cosmn(m1,n1)
                enddo
              enddo
            enddo
          enddo
          bm(jn,jn)=bm(jn,jn)+(n*omega(j))**2
        enddo
      enddo
      kax=ktfaddr(kxm2l(bm,nj*nphi,nj*nphi,nj*nphi,.false.))
      return
      end

      subroutine csroysetup(
     $ rho,zl,wi,aj,omega,zjn,dpjn,
     $ rhoj,
     $ aw,r56,r65,dz,nr,nj,mphi,zspan,dt,irtc)
      implicit none
      integer*4 nr,i,nrh,mphi,nj,i1,irtc
      real*8 rho(nr),zl(nr),zi(nr),wi(nr),aw,dz,
     $     v(nr),ddv(nr/2),work(nr/2),ajz(nr/2),
     $     aj(nj),omega(nj),zjn(mphi,nj),f,f1,
     $     zc,dk,zspan,pi,ox,ak,dt,csrwkick,s,
     $     r56,r65,rhoth,z,dpjn(mphi,nj),rhoj(nj)
      parameter (pi=3.14159265358979324d0)
      nrh=nr/2
      dk=pi/dz/nr
      do i=3,nr-1,2
        ak=dk*(i-1)
        zi(i)  = zl(i+1)/ak
        zi(i+1)=-zl(i)/ak
      enddo
      zi(1)=0.d0
      zi(2)=0.d0
      wi=zi/nr/dz
      call tftrr(wi,nr,.false.)
      v=rho
      call csrwake(v,zl,nr)
      call spline1(nrh,v,ddv,work,0,0)
      zc=0.d0
      s=0.d0
      do i=1,nrh
        zc=zc+rho(i)*csrwkick(nrh,i*dz,0.d0,dz,v,ddv,aw,r65)
        s=s+rho(i)
      enddo
      zc=zc/s/r65
      rhoth=exp(-zspan**2/2.d0)/sqrt(2.d0*pi)/dz;
      do i=1,nrh
        if(rho(i) .ge. rhoth)then
          i1=i
          exit
        endif
      enddo
      f1=csrwkick(nrh,dz-zc,zc,dz,v,ddv,aw,r65)
      do i=1,nrh
        z=dz*i-zc
        f=csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
        if(f*f1 .ge. 0.d0)then
          call csrgetj(
     $         nrh,z,zc,dz,v,ddv,aw,r56,r65,dt,ajz(i),ox)
        else
          ajz(i)=0.d0
        endif
      enddo
      call csrjasign(nrh,nj,i1,aj,omega,
     $     zjn,dpjn,rhoj,ajz,rho,zc,dz,v,ddv,
     $     aw,r56,r65,dt,mphi,irtc)
      return
      end

      subroutine csrjasign(
     $     nrh,nj,i1,aj,omega,zjn,dpjn,rhoj,
     $     ajz,rho,zc,dz,v,ddv,
     $     aw,r56,r65,dt,mphi,irtc)
      implicit none
      integer*4 nrh,nj,i1,m,m0,k,i,mphi,iz,irtc
      real*8 aj(nj),omega(nj),zjn(mphi,nj),dpjn(mphi,nj),
     $     ajz(nrh),v(nrh),ddv(nrh),rho(nrh),rhoj(nj),zi,
     $     zc,dz,aw,r56,r65,dt,pi,dj,ajn,dt1,z,p,csrwkick,f,
     $srho
      parameter (pi=3.14159265358979324d0)
      dj=ajz(i1)/nj
      k=0
      m=1
      m0=0
      do while(m0 .ne. m)
        ajn=dj*(k+0.5d0)
        m0=m
        do i=1,nrh-1
          if(ajz(i) .gt. ajn)then
            if(ajz(i+1) .gt. 0.d0 .and. ajz(i+1) .le. ajn)then
              zjn(1,m)=i*dz-zc+(ajz(i)-ajn)/(ajz(i)-ajz(i+1))*dz
              call csrgetj(
     $             nrh,zjn(1,m),zc,dz,v,ddv,aw,r56,r65,dt,aj(m),
     $             omega(m))
              m=m+1
              if(m .gt. nj)then
                go to 10
              endif
            endif
          endif
        enddo
        k=k+1
      enddo
      irtc=-1
      return
 10   srho=0.d0
      do i=1,nj
        dt1=pi/omega(i)/mphi
        z=zjn(1,i)
        p=0.5d0*dt1*csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
        zi=(z+zc)/dz
        iz=int(zi)
        rhoj(i)=(rho(iz)+(zi-iz)*(rho(iz+1)-rho(iz)))
        srho=srho+rhoj(i)
        do k=1,mphi
          zjn(k,i)=z
          dpjn(k,i)=p
          z=z+p*r56*dt1
          f=csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
          p=p+dt1*f
        enddo
      enddo
      srho=1.d0/2.d0/pi/srho
      do i=1,nj
        rhoj(i)=sqrt(rhoj(i)*srho)
      enddo
      irtc=0
      return
      end

      real*8 function csrfindc(nrh,zc,dz,v,ddv,aw,r65)
      implicit none
      integer*4 nrh,i
      real*8 zc,dz,v(nrh),ddv(nrh),r65,aw,z,f,f1,
     $     csrwkick
      f1=csrwkick(nrh,dz-zc,zc,dz,v,ddv,aw,r65)
      do i=2,nrh
        z=dz*i-zc
        f=csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
        if(f*f1 .le. 0)then
          csrfindc=z-dz*f1/(f-f1)
          exit
        endif
        f1=f
      enddo
      csrfindc=0.d0
      return
      end

      subroutine csrgetj(nrh,z1,zc,dz,v,ddv,aw,r56,r65,dt,aj,omega)
      implicit none
      integer*4 nrh,n
      real*8 z1,zc,dz,v(nrh),ddv(nrh),aw,r56,r65,aj,omega,pi,ze,
     $     z,p,p0,f,t0,dt,csrwkick
      
      parameter (pi=3.14159265358979324d0)
      z=z1
      p=0.5d0*dt*csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
      p0=p
      aj=0.d0
      n=0
      f=1.d0
      do while(p*p0 .gt. 0.d0)
        z=z+p*r56*dt
        f=csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
        aj=aj+p*r56*(p+dt*f*.5d0)
        p=p+dt*f
        n=n+1
      enddo
      t0=2.d0*(n*dt+.5d0*dt-p/f)
      omega=2.d0*pi/t0
      aj=abs(aj*dt/2.d0/pi)
      ze=z
      return
      end

      real*8 function csrwkick(nrh,z,zc,dz,v,ddv,aw,r65)
      implicit none
      integer*4 nrh,m
      real*8 v(nrh),ddv(nrh),aw,z,dz,zc,zk,r65,
     $     a,b,c,d,w
      
      zk=(z+zc)/dz
      m=int(zk)
      if(m .gt. 0 .and. m .lt. nrh)then
        b=zk-m
        a=1.d0-b
        c=-a*b*(a+1.d0)
        d=-a*b*(b+1.d0)
        w=a*v(m)+b*v(m+1)+c*ddv(m)+d*ddv(m+1)
        csrwkick=aw*w+r65*z
      else
        csrwkick=r65*z
      endif
      return
      end
