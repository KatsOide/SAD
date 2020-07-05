      subroutine tfcsrvec(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klpipe,klr,kli,klxi
      integer*4, allocatable :: ind(:)
      integer*4 isp1,irtc,m1,n1,i,j,jm,mm,itfmessage
      real*8 an
      if(isp .ne. isp1+2)then
        go to 9000
      endif
      if(.not. tflistq(dtastk(isp1+1),kl))then
        go to 9010
      endif
      if(.not. tfreallistq(dtastk(isp1+2),klpipe))then
        go to 9010
      endif
      m1=klpipe%nl
      an=0.d0
      do i=1,m1
        an=max(an,klpipe%rbody(i))
      enddo
      n1=int(an)
      allocate(ind(1:(n1+2)*(m1+2)))
c      kaind=ktaloc(((n1+2)*(m1+2)+1)/2+1)
      call csrind(ind,klpipe%rbody(1:m1),m1,n1,mm)
      if(kl%nl .eq. m1)then
        kx=kxavaloc(-1,mm,klr)
        do i=1,m1
          if(.not. tfreallistq(kl%dbody(i),kli))then
            deallocate(ind)
            go to 9010
          endif
          do j=1,n1
            jm=ind(j*(m1+2)+i+1)
            if(jm .ne. 0)then
              klr%rbody(jm)=kli%rbody(j)
            endif
          enddo
        enddo
      elseif(kl%nl .eq. mm)then
        kx=kxadaloc(-1,m1,klx)
        do i=1,m1
          klx%dbody(i)=kxavaloc(0,n1,klxi)
          do j=1,n1
            jm=ind(j*(m1+2)+i+1)
            if(jm .ne. 0)then
              klxi%rbody(j)=kl%rbody(jm)
            else
              klxi%rbody(j)=0.d0
            endif
          enddo
        enddo
      else
        deallocate(ind)
        go to 9010
      endif
      deallocate(ind)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"5"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $       '"Vector or Matrix, Pipe"')
      return
      end

      subroutine tfcsrmat(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: klpipe,klp
      integer*4, allocatable :: ind(:)
      integer*4 isp1,irtc,i,ic,m1,n1,itfmessage,mm,mode,nm
      real*8 dx,dy,omega,rho,an,fudge0,fudge1,fudge2,fudge3,
     $     tilt
      if(isp .ne. isp1+6)then
        go to 9000
      endif
      if(.not. tfreallistq(dtastk(isp1+1),klpipe))then
        go to 9010
      endif
      if(tfreallistq(dtastk(isp1+2),klp))then
        if(klp%nl .ne. 7)then
          go to 9010
        endif
        dx    =klp%rbody(1)
        dy    =klp%rbody(2)
        fudge0=klp%rbody(3)
        fudge1=klp%rbody(4)
        fudge2=klp%rbody(5)
        fudge3=klp%rbody(6)
        tilt  =klp%rbody(7)
      else
        go to 9010
      endif
      if(ktfnonrealq(ktastk(isp1+3)))then
        go to 9010
      endif
      omega=rtastk(isp1+3)
      if(ktfnonrealq(ktastk(isp1+4)))then
        go to 9010
      endif
      rho=rtastk(isp1+4)
      if(ktfnonrealq(ktastk(isp1+5)))then
        go to 9010
      endif
      ic=int(rtastk(isp1+5))
      if(ktfnonrealq(ktastk(isp1+6)))then
        go to 9010
      endif
      mode=int(rtastk(isp1+6))
      m1=klpipe%nl
      an=0.d0
      nm=0
      do i=1,m1
        nm=nm+int(klpipe%rbody(i))
        an=max(an,klpipe%rbody(i))
      enddo
      n1=int(an)
      if(mode .eq. 1000)then
        call csrhelmat(kx,klpipe%rbody(1),dx,dy,m1,n1,ic,nm)
      else
        allocate(ind(1:(n1+2)*(m1+2)))
        call csrind(ind,klpipe%rbody(1),m1,n1,mm)
        call csrmat(kx,ind,
     $       dx,dy,omega,rho,
     $       fudge0,fudge1,fudge2,fudge3,tilt,
     $       mode,mm,m1,n1,ic)
        deallocate(ind)
      endif
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"5"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $'"Pipe, {dx, dy, fudge0, fudge1, fudge2, fudge3, tilt},'//
     $     ' omega/c, rho, center, mode"')
      return
      end

      subroutine csrmat(kx,ind,dx,dy,omega,rho,
     $     fudge0,fudge1,fudge2,fudge3,tilt,
     $     mode,mm,m,n,ic)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 m,n,mm,i,j,jm,ic,mode,i1,j1,im1,jm1
      integer*4 ind(0:m+1,0:n+1)
      real*8 acm(mm,mm),dx,dy,omega,rho,or,dx1,dy1,
     $     x,r,cfx,cfy,cg,ce,rdx,rdy,ac,
     $     cfy1,cfxa1,cfxa2,cfya,cgy,cfxy,
     $     fudge0,fudge1,fudge2,fudge3,fc,tilt
      acm=0.d0
c      parameter (fudge1=4.d0/3.d0,fudge2=2.d0/3.d0)
      if(mode .eq. 110)then
c Es resistive
      elseif(mode .eq. 220)then
c Er resistive
        ac=0.5d0/omega
        cfy=ac*2.d0/dy**3
        do i=1,m
          do j=1,n
            jm=ind(i,j)
            if(jm .eq. 0)then
              exit
            endif
            if(ind(i,j+1) .eq. 0)then
              acm(jm,1)=cfy
            endif
          enddo
        enddo
      elseif(mode .eq. 330)then
c Ey resistive
        ac=0.5d0/omega
        cfx=ac*2.d0/dx**3
        do i=1,m
          do j=1,n
            jm=ind(i,j)
            if(jm .eq. 0)then
              exit
            endif
            if(ind(i-1,j) .eq. 0)then
              acm(jm,1)=cfx
            endif
            if(ind(i+1,j) .eq. 0)then
              acm(jm,1)=acm(jm,1)+cfx
            endif
          enddo
        enddo
      elseif(mode .eq. 22 .or. mode .eq. 33)then
c 11: Es
c 22: Er
c 33: Ey
        or=omega*rho
        do i=1,m
          dx1=dx-(i-ic-0.5d0)*tilt/m*dx
          dy1=dy+(i-ic-0.5d0)*tilt/m*dy
          if(rho .eq. 0.d0)then
            rdx=1.d0/dx1
            rdy=1.d0/dy1
            ce=0.d0
            cg=0.d0
            ac=0.5d0/omega
          else
            x=(i-ic-0.5d0)*dx
            r=rho+x
            rdx=r/dx1
            rdy=r/dy1
            if(mode .eq. 22)then
              ac=0.5d0*omega/(or**2-1.d0)
              ce=(3.d0+omega**2*x*(r+rho))*ac
              cg=1.5d0*rdx*ac
            else
              ac=0.5d0/or/rho
              ce=omega**2*x*(r+rho)*ac
              cg=0.5d0*rdx*ac
            endif
          endif
          cfx=rdx**2*ac
          cfy=rdy**2*ac
          If(mode .eq. 22)then
            cfy1 = cfy
            cfya =-cfy
            cfxa1= cfx-cg
            cfxa2= cfx+cg
          else
            cfy1 =-cfy
            cfya = cfy
            cfxa1=-cfx+cg
            cfxa2=-cfx-cg
          endif
          do j=1,n
            jm=ind(i,j)
            if(jm .eq. 0)then
              exit
            endif
            j1=ind(i,j+1)
            im1=ind(i-1,j)
            i1 =ind(i+1,j)
            if((im1 .eq. 0 .or. i1 .eq. 0) .and. j1 .eq. 0)then
              fc=fudge0
            elseif(j1 .ne. 0 .and. im1 .ne. 0 .and. i1 .ne. 0 .and.
     $             (ind(i+1,j+1) .eq. 0 .or. ind(i-1,j+1) .eq. 0))then
              fc=fudge2
            else
              fc=1.d0
            endif
            acm(jm,jm)=ce-2.d0*(cfx+cfy)*fc
            if(j .eq. 1)then
              acm(jm,jm)=acm(jm,jm)+cfy1
              acm(jm,j1)=cfy
            else
              jm1=ind(i,j-1)
              acm(jm,jm1)=cfy
              if(j1 .ne. 0)then
                acm(jm,j1 )=cfy
              else
                acm(jm,jm)=acm(jm,jm)+cfya*fc
              endif
            endif
            if(im1 .ne. 0)then
              acm(jm,im1)=cfx-cg
            else
              acm(jm,jm)=acm(jm,jm)+cfxa1*fc
            endif
            if(i1 .ne. 0)then
              acm(jm,i1)=cfx+cg
            else
              acm(jm,jm)=acm(jm,jm)+cfxa2*fc
            endif
          enddo
        enddo
      elseif(mode .eq. 23 .or. mode .eq. 32)then
c 23; Ey -> Er
c 32: Er -> Ey
        or=omega*rho
        do i=1,m
          dx1=dx-(i-ic-0.5d0)*tilt/m*dx
          dy1=dy+(i-ic-0.5d0)*tilt/m*dy
          if(rho .eq. 0.d0)then
            rdx=1.d0/dx1
            rdy=1.d0/dy1
            ac=0.5d0/omega
            cgy=0.d0
          else
            x=(i-ic-0.5d0)*dx
            r=rho+x
            rdx=r/dx1
            rdy=r/dy1
            if(mode .eq. 23)then
              ac=0.5d0*omega/(or**2-1.d0)
              cgy=rdy*ac
            else
              ac=0.5d0/or/rho
              cgy=0.d0
            endif
          endif
          cfx =rdx**2*ac
          cfy =rdy**2*ac
          cfxy=rdx*rdy*ac
          do j=1,n
            jm=ind(i,j)
            if(jm .eq. 0)then
              exit
            endif
            j1 =ind(i,j+1)
            jm1=ind(i,j-1)
            im1=ind(i-1,j)
            i1 =ind(i+1,j)
            if(j .eq. 1)then
              acm(jm,jm )= cgy
            endif
            if(jm1 .ne. 0)then
              acm(jm,jm1)=-cgy
            endif
            if(j1 .eq. 0)then
              acm(jm,jm)=acm(jm,jm)+cgy
              if(im1 .eq. 0)then
                acm(jm,jm) =acm(jm,jm)-cfxy*fudge1
              endif
              if(i1 .eq. 0)then
                acm(jm,jm) =acm(jm,jm)+cfxy*fudge1
              endif
            else
              acm(jm,j1)=cgy
              if(im1 .ne. 0 .and. i1 .ne. 0)then
                if(ind(i-1,j+1) .eq. 0)then
                  acm(jm,jm) =acm(jm,jm)-cfxy*fudge3
                endif
                if(ind(i+1,j+1) .eq. 0)then
                  acm(jm,jm) =acm(jm,jm)+cfxy*fudge3
                endif
              endif
            endif
          enddo
        enddo
      endif
      if(mode .ge. 100)then
        kx=kxm2l(acm,0,mm,mm,.false.)
      else
        kx=kxm2l(acm,mm,mm,mm,.false.)
      endif
      return
      end

      subroutine csrind(ind,pipe,m,n,mm)
      implicit none
      integer*4 m,n,i,j,k,mm
      integer*4 ind(0:m+1,0:n+1)
      real*8 pipe(m)
      k=0
      do i=1,m
        do j=1,int(pipe(i))
          k=k+1
          ind(i,j)=k
        enddo
        do j=int(pipe(i))+1,n+1
          ind(i,j)=0
        enddo
        ind(i,0)=0
      enddo
      do j=0,n+1
        ind(0,j)=0
        ind(m+1,j)=0
      enddo
      mm=k
      return
      end

      subroutine tfcsrinit(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kex,key
      type (sad_dlist), pointer :: klx,kl1
      type (sad_rlist), pointer :: klprof,kl11
      integer*4 isp1,irtc,nphimax,mx,my,nphi,m,itfmessage,ic
      parameter (nphimax=16384)
      real*8 dx,dy,rho
      real*8 phix(nphimax),phiy(nphimax)
      dx=0.d0
      if(isp .ne. isp1+1)then
        go to 9030
      endif
      if(.not. tflistq(dtastk(isp1+1),kl1))then
        go to 9000
      endif
      if(kl1%nl .ne. 4)then
        go to 9000
      endif
      if(tfreallistq(kl1%dbody(1),kl11))then
        if(kl11%nl .ne. 2)then
          go to 9000
        endif
        dx=kl11%rbody(1)
        dy=kl11%rbody(2)
      elseif(ktfnonrealq(kl1%dbody(1)))then
        go to 9000
      else
        dy=dx
      endif
      if(dx .le. 0.d0 .or. dy .le. 0.d0)then
        go to 9000
      endif
      if(ktfnonrealq(kl1%dbody(2),ic))then
        go to 9000
      endif
      if(ic .le. 0)then
        go to 9000
      endif
      if(ktfnonrealq(kl1%dbody(3),rho))then
        go to 9000
      endif
      if(.not. tfreallistq(kl1%dbody(4),klprof))then
        go to 9000
      endif
      mx=klprof%nl
      if(mx .le. ic)then
        go to 9010
      endif
      call csrphimesh(klprof%rbody(1),
     $     phix,phiy,dx,dy,ic,mx,nphi,my,nphimax)
      if(nphi .ge. nphimax)then
        go to 9020
      endif
      m=mx*my
      call csr2dpoisson(phix,phiy,kex,key,
     $     dx,dy,klprof%rbody(1),nphi,ic,mx,my,rho,3)
      kx=kxadaloc(-1,2,klx)
      klx%dbody(1)=dtfcopy1(kex)
      klx%dbody(2)=dtfcopy1(key)
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $'"{dx or {dx, dy}, center, rho, pipe:{y1, y2, ...}} for #1"')
      return
 9010 irtc=itfmessage(9,'General::wrongval',
     $       '"Center must be within beam pipe"')
      return
 9020 irtc=itfmessage(9,'General::wrongval',
     $       '"Too large number of meshes"')
      return
 9030 irtc=itfmessage(9,'General::narg',
     $       '"1"')
      return
      end

      subroutine csrhelmat(kx,pipe,dx,dy,mx,my,ic,nm)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 mx,my,nphi,ic,nm,nphimax,i,j,ind,k
      parameter (nphimax=4000)
      real*8 pipe(mx),phix(0:nphimax),phiy(0:nphimax),
     $     tx(nphimax),ty(nphimax),
     $     xc(nm),yc(nm),xp1(nm),xp2(nm),yp1(nm),
     $     psi(nm*2,nm*2),dx,dy,x,y,cs,sx,sy,r
      real*8 pi
      parameter (pi=3.14159265358979324d0)
      call csrphimesh(pipe,phix,phiy,dx,dy,ic,mx,nphi,my,nphimax)
      nphi=nphi-1
      call csrdoublephimesh(phix,phiy,nphi)
      do i=1,nphi
        sx=phix(i+1)-phix(i-1)
        sy=phiy(i+1)-phiy(i-1)
        r=hypot(sx,sy)
        tx(i)=sx/r
        ty(i)=sy/r
      enddo
      ind=0
      do i=1,mx
        x=(i-ic-0.5d0)*dx
        do j=1,int(pipe(i))
          y=(j-0.5d0)*dy
          ind=ind+1
          xc(ind)=x
          yc(ind)=y
          yp1(ind)=pipe(i)*dy
          xp1(ind)=-ic*dx
          xp2(ind)=(mx-ic)*dx
          do k=j-1,1,-1
            if(int(pipe(k)) .lt. j)then
              xp1(ind)=(k-ic)*dx
              exit
            endif
          enddo
          do k=j+1,mx
            if(int(pipe(k)) .lt. j)then
              xp2(ind)=(k-ic-1)*dx
              exit
            endif
          enddo
        enddo
      enddo
      cs=dx*dy/2.d0/pi
      call csr2dpoisson1(phix,phiy,psi,xc,yc,tx,ty,
     $     xp1,xp2,yp1,
     $     cs,nphi,nm)
      kx=kxm2l(psi,nm*2,nm*2,nm*2,.false.)
      return
      end

      subroutine csrdoublephimesh(phix,phiy,nphi)
      implicit none
      integer*4 nphi,i,ii
      real*8 phix(0:nphi*2),phiy(0:nphi*2)
      do i=nphi,1,-1
        ii=i*2
        phix(ii)=phix(i)
        phiy(ii)=phiy(i)
        phix(ii-1)=(phix(i)+phix(i-1))*.5d0
        phiy(ii-1)=(phiy(i)+phiy(i-1))*.5d0
      enddo
      nphi=nphi*2-1
      return
      end

      real*8 function phi1x(x,y,sa)
      implicit none
      real*8 x,y,r,sa
      real*8 pi
      parameter (pi=3.14159265358979324d0)
      r=x**2+y**2
      if(r .eq. 0.d0)then
        phi1x=.5d0-log(sa/pi)*.5d0
      else
        phi1x=-log(r)*.5d0+sa/24.d0/r
      endif
      return
      end

      real*8 function phi1y(x,y)
      implicit none
      real*8 x,y
      if(y .eq. 0.d0)then
        phi1y=0.d0
      else
        phi1y=-atan(x/y)
      endif
      return
      end

      subroutine csr2dpoisson1(phix,phiy,psi,xc,yc,tx,ty,
     $     xp1,xp2,yp1,cs,nphi,nm)
      implicit none
      integer*4 nphi,i,j,k,nv,nm
      real*8 phix(0:nphi),phiy(0:nphi),tx(nphi),ty(nphi),
     $     psi(nm*2,nm*2),xc(nm),yc(nm),
     $     xp1(nm),xp2(nm),yp1(nm),
     $     a1(nphi-1,nphi),er1(nphi,nm*2),b1(nphi-1,nm*2),
     $     xp,yp,xf,yf,exa1(nm),eya1(nm),exa2(nm),eya2(nm),
     $     rsqi1,rsqi2,cs,eps,sa,ws2,w1,w2,
     $     phi1x,phi1y
      real*8 pi,el
      parameter (pi=3.14159265358979324d0,el=1.d-8)
      sa=cs*2.d0*pi
      ws2=sa/16.d0
      nv=nphi
      eps=1.d-4
      do i=2,nphi-1
        xf=phix(i)
        yf=phiy(i)
        do k=1,nm
          xp=xc(k)
          yp=yc(k)
          b1(i,k)=-tx(i)*( phi1x(xf-xp,yf-yp,sa)+phi1x(xf-xp,yf+yp,sa))
     $            -ty(i)*(-phi1y(yf-yp,xf-xp)-phi1y(yf+yp,xf-xp))
          b1(i,nm+k)=
     $            -tx(i)*( phi1y(yf-yp,xf-xp)-phi1y(yf+yp,xf-xp))
     $            -ty(i)*( phi1x(yf-yp,xf-xp,sa)-phi1x(yf+yp,xf-xp,sa))
        enddo
        do j=1,nphi
          xp=phix(j)
          yp=phiy(j)
          if(i .ne. j)then
            rsqi1=1.d0/((xf-xp)**2+(yf-yp)**2)
            w1=rsqi1*(1.d0+ws2*rsqi1)
          else
            w1=0.d0
          endif
          rsqi2=1.d0/((xf-xp)**2+(yf+yp)**2)
          w2=rsqi2*(1.d0+ws2*rsqi2)
          a1(i,j)= tx(i)*(xf-xp)*(w1+w2)+ty(i)*((yf-yp)*w1+(yf+yp)*w2)
        enddo
        a1(1,i)=1.d0
      enddo
      a1(1,1)=1.d0
      a1(1,nphi)=1.d0
      do i=1,nm
        b1(1,i)= -yp1(i)+yc(i)
        b1(1,nm+i)= yc(i)
      enddo
      call tsolvm(a1,b1,er1,nphi-1,nv,nm*2,nphi-1,nphi-1,nv,eps,.true.)
      do k=1,nm
        xf=xc(k)
        yf=yc(k)
        do i=1,nm
          xp=xc(i)
          yp=yc(i)
          exa1(i)= phi1x(xf-xp,yf-yp,sa)+phi1x(xf-xp,yf+yp,sa)
          eya1(i)=-phi1y(yf-yp,xf-xp)-phi1y(yf+yp,xf-xp)
          exa2(i)= phi1y(yf-yp,xf-xp)-phi1y(yf+yp,xf-xp)
          eya2(i)= phi1x(yf-yp,xf-xp,sa)-phi1x(yf+yp,xf-xp,sa)
        enddo
        do j=1,nphi
          xp=phix(j)
          yp=phiy(j)
          rsqi1=1.d0/((xf-xp)**2+(yf-yp)**2)
          rsqi2=1.d0/((xf-xp)**2+(yf+yp)**2)
          w1=rsqi1*(1.d0+ws2*rsqi1)
          w2=rsqi2*(1.d0+ws2*rsqi2)
          do i=1,nm
            exa1(i)=exa1(i)+ (xf-xp)*(w1+w2)*er1(j,i)
            eya1(i)=eya1(i)+((yf-yp)*w1+(yf+yp)*w2)*er1(j,i)
            exa2(i)=exa2(i)+ (xf-xp)*(w1+w2)*er1(j,nm+i)
            eya2(i)=eya2(i)+((yf-yp)*w1+(yf+yp)*w2)*er1(j,nm+i)
          enddo
        enddo
        do i=1,nm
          psi(k,   i   )=exa1(i)*cs
          psi(k,   nm+i)=exa2(i)*cs
          psi(nm+k,i   )=eya1(i)*cs
          psi(nm+k,nm+i)=eya2(i)*cs
        enddo
      enddo
      return
      end

      subroutine csr2dpoisson(phix,phiy,kex,key,dx,dy,
     $     pipe,nphi,ic,mx,my,rho0,iter)
      use tfstk
      implicit none
      type (sad_descriptor) kex,key
      integer*4 nphi,mx,my,i,j,k,l,ic,it
      real*8 pipe(mx),bs
      real*8 phix(nphi),phiy(nphi),ex(mx,my),ey(mx,my),
     $     a(nphi-1,2*nphi-3),dx,dy,rho(2*nphi-3),b(nphi-1),
     $     xp,yp,xf,yf,exa,eya,rsqi1,rsqi2,rsq1,x,
     $     ddx,ddy1,ddy2,
     $     rf,rho0,y,q1(mx,my),eps
      integer*4 ind(0:mx+1,0:my+1),mm,iter
      real*8 pi
      parameter (pi=3.14159265358979324d0)
      eps=min(1.d-2,max(1.d0/(nphi-1)**2,1.d-4))
      call csrind(ind,pipe,mx,my,mm)
      do it=1,iter
        if(it .gt. 1)then
          do k=1,mx
            rf=(dx*(k-ic-0.5d0)+rho0)*(2.d0*pi)/(dx*dy)
            do l=1,int(pipe(k))
              q1(k,l)=-ex(k,l)/rf
            enddo
          enddo
        endif
        do i=2,nphi-1
          xf=phix(i)
          yf=phiy(i)
          bs=-log(xf**2+yf**2)
          do j=1,nphi-1
            xp=(phix(j)+phix(j+1))*.5d0
            yp=(phiy(j)+phiy(j+1))*.5d0
            a(i-1,j)=-log(((xf-xp)**2+(yf-yp)**2)*
     $           ((xf-xp)**2+(yf+yp)**2))
          enddo
          do j=2,nphi-1
            if(j .ne. i)then
              xp=phix(j)
              yp=phiy(j)
              a(i-1,j+nphi-2)=-log(((xf-xp)**2+(yf-yp)**2)*
     $             ((xf-xp)**2+(yf+yp)**2))
            else
              a(i-1,j+nphi-2)=0.d0
            endif
          enddo
          if(it .gt. 1)then
            do k=1,mx
              xp=dx*(k-ic-0.5d0)
              do l=1,int(pipe(k))
                yp=dy*(l-0.5d0)
                bs=bs-log(((xf-xp)**2+(yf-yp)**2)*
     $               ((xf-xp)**2+(yf+yp)**2))*q1(k,l)
              enddo
            enddo
          endif
          b(i-1)=bs
        enddo
        do i=1,2*nphi-3
          a(nphi-1,i)=1.d0
        enddo
        b(nphi-1)=0.5d0
c        call tsolvg(a,b,rho,nphi-1,nphi-1,nphi-1)
        call tsolva(a,b,rho,nphi-1,2*nphi-3,nphi-1,eps)
        do k=1,mx
          xf=dx*(k-ic-0.5d0)
          do l=1,int(pipe(k))
            yf=dy*(l-0.5d0)
            exa=xf/(xf**2+yf**2)
            eya=yf/(xf**2+yf**2)
            do i=1,nphi-1
              xp=(phix(i)+phix(i+1))*.5d0
              yp=(phiy(i)+phiy(i+1))*.5d0
              ddx=xf-xp
              ddy1=yf-yp
              ddy2=yf+yp
              rsqi1=1.d0/(ddx**2+ddy1**2)
              rsqi2=1.d0/(ddx**2+ddy2**2)
              exa=exa-rho(i)*ddx*(rsqi1+rsqi2)
              eya=eya-rho(i)*(ddy1*rsqi1+ddy2*rsqi2)
            enddo
            do i=2,nphi-1
              xp=phix(i)
              yp=phiy(i)
              ddx=xf-xp
              ddy1=yf-yp
              ddy2=yf+yp
              rsqi1=1.d0/(ddx**2+ddy1**2)
              rsqi2=1.d0/(ddx**2+ddy2**2)
              exa=exa-rho(i+nphi-2)*ddx*(rsqi1+rsqi2)
              eya=eya-rho(i+nphi-2)*(ddy1*rsqi1+ddy2*rsqi2)
            enddo
            if(it .gt. 1)then
              do i=1,mx
                x=dx*(i-ic-0.5d0)
                ddx=xf-x
                do j=1,int(pipe(i))
                  y=dy*(j-0.5d0)
                  ddy1=yf-y
                  rsq1=ddx**2+ddy1**2
                  if(rsq1 .eq. 0.d0)then
                    rsqi1=0.d0
                  else
                    rsqi1=1.d0/rsq1
                  endif
                  ddy2=yf+y
                  rsqi2=1.d0/(ddx**2+ddy2**2)
                  exa=exa+q1(i,j)*ddx*(rsqi1+rsqi2)
                  eya=eya+q1(i,j)*(ddy1*rsqi1+ddy2*rsqi2)
                enddo
              enddo
            endif
            ex(k,l)=exa
            ey(k,l)=eya
          enddo
          do l=int(pipe(k))+1,my
            ex(k,l)=0.d0
            ey(k,l)=0.d0
          enddo
        enddo
      enddo
      kex=kxm2l(ex,mx,my,mx,.false.)
      key=kxm2l(ey,mx,my,mx,.false.)
      return
      end

      subroutine csrphimesh(pipe,phix,phiy,dx,dy,ic,mx,nphi,my,nphimax)
      implicit none
      integer*4 mx,ic,i,j,k,jy,nphi,my,nphimax,l
      real*8 pipe(mx)
      real*8 phix(*),phiy(*),dx,dy,x
      jy=0
      k=1
      my=0
      do i=1,mx
        x=dx*(i-ic-1)
        l=int(pipe(i))
        my=max(my,l)
        do j=jy+1,l
          k=k+1
          if(k .ge. nphimax)then
            go to 9000
          endif
          phix(k)=x
          phiy(k)=dy*(j-0.5d0)
        enddo
        do j=jy,l+1,-1
          k=k+1
          if(k .ge. nphimax)then
            go to 9000
          endif
          phix(k)=x
          phiy(k)=dy*(j-0.5d0)
        enddo
        if(l .gt. 0)then
          k=k+1
          if(k .ge. nphimax)then
            go to 9000
          endif
          phix(k)=(i-ic-0.5d0)*dx
          phiy(k)=l*dy
        endif
        jy=l
      enddo
      x=dx*(mx-ic)
      if(jy .gt. 0)then
        do j=jy,1,-1
          k=k+1
          if(k .ge. nphimax)then
            go to 9000
          endif
          phix(k)=x
          phiy(k)=dy*(j-0.5d0)
        enddo
      endif
      phix(1)= phix(2)
      phiy(1)=-phiy(2)
      k=k+1
      if(k .ge. nphimax)then
        go to 9000
      endif
      phix(k)= phix(k-1)
      phiy(k)=-phiy(k-1)
 9000 nphi=k
c      write(*,'(i6,1p2g15.7)')(k,phix(k),phiy(k),k=1,nphi)
      return
      end
