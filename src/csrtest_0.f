      subroutine tfcsrvec(isp1,itx,iax,vx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*4 isp1,itx,iax,irtc,ia,iapipe,ipipe,
     $     m1,n1,iaind,italoc,itavaloc,itadaloc,i,j,jm,
     $     iav,iaxi,mm,itfmessage,iavi
      real*8 vx,an
      logical*4 tflistq,tfnumlistq
      if(isp .ne. isp1+2)then
        go to 9000
      endif
      ia=itastk(2,isp1+1)
      if(.not. tflistq(itastk(1,isp1+1),ia))then
        go to 9010
      endif
      iapipe=itastk(2,isp1+2)
      if(.not. tfnumlistq(itastk(1,isp1+2),iapipe))then
        go to 9010
      endif
      m1=ilist(2,iapipe-1)
      ipipe=ilist(2,iapipe+1)+1
      an=0.d0
      do i=1,m1
        an=max(an,rlist(ipipe+i-1))
      enddo
      n1=an
      iaind=italoc(((n1+2)*(m1+2)+1)/2+1)
      call csrind(ilist(1,iaind),rlist(ipipe),m1,n1,mm)
      if(ilist(2,ia-1) .eq. m1)then
        iax=itavaloc(-1,mm)
        do i=1,m1
          iavi=ilist(2,ilist(2,ilist(1,ia+1)+i)+1)
          do j=1,n1
            jm=ilist(j*(m1+2)+i+1,iaind)
            if(jm .ne. 0)then
              rlist(iax+1+jm)=rlist(iavi+j)
            endif
          enddo
        enddo
      elseif(ilist(2,ia-1) .eq. mm)then
        iax=itadaloc(-1,m1)
        iav=ilist(2,ia+1)
        do i=1,m1
          iaxi=itavaloc(0,n1)
          call tfsetlist(ntflist,iaxi,0.d0,iax,i)
          do j=1,n1
            jm=ilist(j*(m1+2)+i+1,iaind)
            if(jm .ne. 0)then
              rlist(iaxi+j+1)=rlist(iav+jm)
            else
              rlist(iaxi+j+1)=0.d0
            endif
          enddo
        enddo
      else
       go to 9010
      endif
      call tfree(iaind)
      itx=ntflist
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"5"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $       '"Vector or Matrix, Pipe"')
      return
      end

      subroutine tfcsrmat(isp1,itx,iax,vx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*4 isp1,itx,iax,irtc,
     $     iapipe,i,ic,m1,n1,itfmessage,
     $     mm,iacm,ipipe,iaind,mode,
     $     italoc,itfm2l,iai,mm1
      real*8 vx,dx,dy,omega,rho,an,
     $     fudge0,fudge1,fudge2,fudge3
      logical*4 tfnumlistq
      if(isp .ne. isp1+6)then
        go to 9000
      endif
      iapipe=itastk(2,isp1+1)
      if(.not. tfnumlistq(itastk(1,isp1+1),iapipe))then
        go to 9010
      endif
      iai=itastk(2,isp1+2)
      if(tfnumlistq(itastk(1,isp1+2),iai))then
        if(ilist(2,iai-1) .ne. 6)then
          go to 9010
        endif
        dx    =rlist(ilist(2,iai+1)+1)
        dy    =rlist(ilist(2,iai+1)+2)
        fudge0=rlist(ilist(2,iai+1)+3)
        fudge1=rlist(ilist(2,iai+1)+4)
        fudge2=rlist(ilist(2,iai+1)+5)
        fudge3=rlist(ilist(2,iai+1)+6)
      else
        go to 9010
      endif
      if(itastk(1,isp1+3) .ne. ntfreal)then
        go to 9010
      endif
      omega=vstk(ivstkoffset+isp1+3)
      if(itastk(1,isp1+4) .ne. ntfreal)then
        go to 9010
      endif
      rho=vstk(ivstkoffset+isp1+4)
      if(itastk(1,isp1+5) .ne. ntfreal)then
        go to 9010
      endif
      ic=vstk(ivstkoffset+isp1+5)
      if(itastk(1,isp1+6) .ne. ntfreal)then
        go to 9010
      endif
      mode=vstk(ivstkoffset+isp1+6)
      m1=ilist(2,iapipe-1)
      ipipe=ilist(2,iapipe+1)+1
      an=0.d0
      do i=1,m1
        an=max(an,rlist(ipipe+i-1))
      enddo
      n1=an
      iaind=italoc(((n1+2)*(m1+2)+1)/2+1)
      call csrind(ilist(1,iaind),rlist(ipipe),m1,n1,mm)
      if(mode .lt. 100)then
        mm1=mm
      else
        mm1=1
      endif
      iacm=italoc(mm*mm1)
      call tclr(rlist(iacm),mm*mm1)
      call csrmat(rlist(iacm),ilist(1,iaind),
     $     dx,dy,omega,rho,
     $     fudge0,fudge1,fudge2,fudge3,
     $     mode,mm,m1,n1,ic)
      if(mm1 .eq. 1)then
        iax=itfm2l(rlist(iacm),0,mm,mm,.false.)
      else
        iax=itfm2l(rlist(iacm),mm,mm,mm,.false.)
      endif
      call tfree(iacm)
      call tfree(iaind)
      itx=ntflist
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"5"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $'"Pipe, {dx, dy, fudge0, fudge1, fudge2, fudge3},'//
     $     ' omega/c, rho, center, mode"')
      return
      end

      subroutine csrmat(acm,ind,
     $     dx,dy,omega,rho,
     $     fudge0,fudge1,fudge2,fudge3,
     $     mode,mm,m,n,ic)
      implicit none
      integer*4 m,n,mm,i,j,jm,ic,mode,i1,j1,im1,jm1
      integer*4 ind(0:m+1,0:n+1)
      real*8 acm(mm,mm),dx,dy,omega,rho,or,
     $     x,r,cfx,cfy,cg,ce,rdx,rdy,ac,
     $     cfy1,cfxa1,cfxa2,cfya,cgy,cfxy,
     $     fudge0,fudge1,fudge2,fudge3,fc
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
          if(rho .eq. 0.d0)then
            rdx=1.d0/dx
            rdy=1.d0/dy
            ce=0.d0
            cg=0.d0
            ac=0.5d0/omega
          else
            x=(i-ic-0.5d0)*dx
            r=rho+x
            rdx=r/dx
            rdy=r/dy
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
              acm(jm,jm)=acm(jm,jm)+cfy1*fc
              if(j1 .ne. 0)then
                acm(jm,j1)=cfy
              endif
            else
              jm1=ind(i,j-1)
              acm(jm,jm1)=cfy
            endif
            if(j1 .ne. 0)then
              acm(jm,j1 )=cfy
            else
              acm(jm,jm)=acm(jm,jm)+cfya*fc
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
          if(rho .eq. 0.d0)then
            rdx=1.d0/dx
            rdy=1.d0/dy
            ac=0.5d0/omega
            cgy=0.d0
          else
            x=(i-ic-0.5d0)*dx
            r=rho+x
            rdx=r/dx
            rdy=r/dy
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

      subroutine tfcsrinit(isp1,itx,iax,vx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*4 isp1,itx,iax,irtc,nphimax,
     $     itdx,iadx,itprof,iaprof,mx,my,nphi,
     $     ibm,m,iex,iey,iaex,iaey,italoca,itfm2l,itaaloc,
     $     itfmessage,ia1,itcent,iacent,ic,itfcopy1,
     $     itrho,iarho
      parameter (nphimax=16384)
      real*8 vx,dx,dy,vcent,vprof,rho
      real*8 phix(nphimax),phiy(nphimax)
      logical*4 tflistq,tfnumlistq
      if(isp .ne. isp1+1)then
        go to 9030
      endif
      ia1=itastk(2,isp1+1)
      if(.not. tflistq(itastk(1,isp1+1),ia1))then
        go to 9000
      endif
      if(ilist(2,ia1-1) .ne. 4)then
        go to 9000
      endif
      call tfgetl(ia1,1,itdx,iadx,dx)
      if(tfnumlistq(itdx,iadx))then
        if(ilist(2,iadx-1) .ne. 2)then
          go to 9000
        endif
        dx=rlist(ilist(2,iadx+1)+1)
        dy=rlist(ilist(2,iadx+1)+2)
      elseif(itdx .ne. ntfreal)then
        go to 9000
      else
        dy=dx
      endif
      if(dx .le. 0.d0 .or. dy .le. 0.d0)then
        go to 9000
      endif
      call tfgetl(ia1,2,itcent,iacent,vcent)
      if(itcent .ne. ntfreal .or. vcent .le. 0.d0)then
        go to 9000
      endif
      ic=vcent
      call tfgetl(ia1,3,itrho,iarho,rho)
      if(itrho .ne. ntfreal)then
        go to 9000
      endif
      call tfgetl(ia1,4,itprof,iaprof,vprof)
      if(.not. tfnumlistq(itprof,iaprof))then
        go to 9000
      endif
      if(ilist(2,iaprof-1) .le. ic)then
        go to 9010
      endif
      mx=ilist(2,iaprof-1)
      ibm=ilist(2,iaprof+1)+1
      call csrphimesh(rlist(ibm),
     $     phix,phiy,dx,dy,ic,mx,nphi,my,nphimax)
      if(nphi .ge. nphimax)then
        go to 9020
      endif
      m=mx*my
      iex=italoca(m)
      iey=italoca(m)
      call csr2dpoisson(phix,phiy,
     $     rlist(iex),rlist(iey),
     $     dx,dy,rlist(ibm),nphi,ic,mx,my,rho,3)
      iaex=itfm2l(rlist(iex),mx,my,mx,.false.)
      iaey=itfm2l(rlist(iey),mx,my,mx,.false.)
      call tfree(iex)
      call tfree(iey)
      itx=ntflist
      iax=itaaloc(-1,2)
      call tfsetlist(ntflist,itfcopy1(iaex),0.d0,iax,1)
      call tfsetlist(ntflist,itfcopy1(iaey),0.d0,iax,2)
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

      subroutine csr2dpoisson(phix,phiy,ex,ey,dx,dy,
     $     pipe,nphi,ic,mx,my,rho0,iter)
      implicit none
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
        l=pipe(i)
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
