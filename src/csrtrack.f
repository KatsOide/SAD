      subroutine tfcsrtrack(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_rlist), pointer :: klr,klz,klp,klzl,kll
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 np,nr,itfmessage,nt,ndiv
      real*8 dz,aw,r56,r65,damp,sigmae,frf,phirf
      if(isp /= isp1+5)then
        go to 9000
      endif
      if(.not. tfreallistq(ktastk(isp1+1),klz))then
        go to 9010
      endif
      np=klz%nl
      if(.not. tfreallistq(ktastk(isp1+2),klp))then
        go to 9010
      endif
      if(klp%nl /= np)then
        go to 9020
      endif
      if(.not. tfreallistq(ktastk(isp1+3),klzl))then
        go to 9010
      endif
      nr=klzl%nl
      if(.not. tfreallistq(ktastk(isp1+4),kll))then
        go to 9010
      endif
      if(kll%nl /= 9)then
        go to 9010
      endif
      dz=kll%rbody(1)
      ndiv=int(max(1.d0,kll%rbody(2)))
      aw=kll%rbody(3)/np
      r56=kll%rbody(4)
      r65=kll%rbody(5)
      frf=kll%rbody(6)
      phirf=kll%rbody(7)
      damp=kll%rbody(8)
      sigmae=kll%rbody(9)
      if(ktfnonrealq(ktastk(isp1+5)))then
        go to 9010
      endif
      nt=int(rtastk(isp1+5))
      kx=kxavaloc(-1,nt*2,klr)
      call csrtrack(klz%rbody(1),klp%rbody(1),np,klzl%rbody(1),nr,
     $     dz,aw,r56,r65,frf,phirf,damp,sigmae,
     $     klr%rbody(1),ndiv,nt)
c
c This function returns {sigz1, sigdp1, ... }. When diverged, zeros are returned in sigz and sigdp.
c
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::narg','"5"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $'"z, dp, zl,'//
     $'{dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmadp},'//
     $'nturn"')
      return
c
c   z(np):   list of z for np particles (m)
c   dp(np):  list of dp/p0 for np paticles
c   zl:      list of impedance (nr points) in the "real fft" format, as the result of rftr(a)
c   dz:      mesh size in z
c   ndiv:    number of divisions per turn to apply wake.
c   inten:   q_bunch/E
c   r56:     r56 (m)
c   r65:     r65 (1/m)
c   frf      rf frequency (Hz)
c   phirf    rf phase (rad): positive for alpha > 0 electron ring
c   damp:    damping/per turn
c   sigmadp: equil. sigmadp
c   nturn:   number of turns.
c
 9020 irtc=itfmessage(9,'General::wrongleng',
     $       '"z","   equal to length of dp"')
      return
      end

      subroutine tfcsrhaissin(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_rlist), pointer :: klr,klzl,kll
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 nr,itfmessage
      real*8 dz,aw,sigz,r65
      if(isp /= isp1+2)then
        go to 9000
      endif
      if(.not. tfreallistq(ktastk(isp1+1),klzl))then
        go to 9010
      endif
      nr=klzl%nl
      if(.not. tfreallistq(ktastk(isp1+2),kll))then
        go to 9010
      endif
      if(kll%nl /= 4)then
        go to 9010
      endif
      dz=kll%rbody(1)
      r65=kll%rbody(3)
      sigz=kll%rbody(4)
      aw=kll%rbody(2)/r65/sigz**2
      kx=kxavaloc(-1,nr,klr)
      call csrhaissin(klr%rbody(1),
     $     klzl%rbody(1),nr,aw,dz,sigz,irtc)
      if(irtc /= 0)then
        go to 9020
      endif
      return
 9000 irtc=itfmessage(9,'General::narg','"2"')
      return
 9010 irtc=itfmessage(9,'General::wrongtype',
     $       '"zl, {dz, inten, r65, sigmaz}"')
c
c   zl:      list of impedance (nr points) in the "real fft" format, as the result of rftr(a)
c   dz:      mesh size in z.
c   inten:   q_bunch/E
c   r65:     r65 (1/m)
c   sigmaz:  bunch length (m)
c
 9020 irtc=itfmessage(9,'General::toomany',
     $       '"iterations"')
      return
      end

      subroutine csrhaissin(rho,zl,nr,aw,dz,sigz,irtc)
      use ftr
      use solv
      implicit none
      integer*4 nr,i,k,imax,irtc,nr1,m,iter
      real*8 rho(nr),zl(nr),zi(nr),f(0:nr/2),g(0:nr/2,0:nr/2),
     $     fact,du0,thre1,ulast(0:nr/2),
     $     u(0:nr/2),u0(0:nr/2),wi(nr),du(0:nr/2),v(nr),
     $     aw,dz,sigz,pi,zc,a,sr,ak,thre,dk,sdr,sdf,
     $     sdf0
      parameter (pi=3.14159265358979324d0,imax=100,
     $     thre=1.d-5,thre1=1.d-1)
c
c     rho: physical particle density: Sum(rho)*dz = 1;
c
c      write(*,*)'CSRHi-0 ',nr,aw,dz,sigz
      iter=0
      fact=1.d0
      sdf0=1.d100
      du0=0.d0
      du(0)=0.d0
      nr1=nr/2
      zc=nr/4*dz
      a=dz/sigz/sqrt(2.d0*pi)
      u(0)=log(a)
      do i=1,nr1
        u0(i)=-((dz*i-zc)/sigz)**2/2.d0
        u(i)=u(0)+u0(i)
        rho(i)=exp(u(i))
        rho(i+nr1)=0.d0
      enddo
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
      do i=1,imax
        v=rho/dz
        call csrwake(v,zi,nr)
        sdf=0.d0
        do k=1,nr1
          f(k)=u(k)-u(0)-u0(k)+aw*v(k)
          do m=1,k
            g(k,m)=aw*wi(k-m+1)*rho(m)
          enddo
          g(k,k)=g(k,k)+1.d0
          do m=k+1,nr1
            g(k,m)=aw*wi(nr-(m-k)+1)*rho(m)
          enddo
          g(k,0)=-1.d0
          g(0,k)=rho(k)
          sdf=sdf+abs(f(k))
        enddo
c        write(*,*)'CSRHi ',i,iter,sdf,du0,du(0)
        if(sdf .gt. sdf0)then
          fact=fact*.5d0
          go to 10
        else
          if(sdf .lt. thre)then
            rho=rho/dz
            irtc=0
            return
          elseif(iter .gt. imax)then
            irtc=-1
            return
          endif
          sdf0=sdf
          fact=min(1.d0,fact*2.d0)
          f(0)=0.d0
          g(0,0)=0.d0
          call tsolvg(g,f,du)
          ulast=u
        endif
 10     iter=iter+1
        sr=0.d0
        do k=1,nr1
          u(k)=ulast(k)-du(k)*fact
          rho(k)=exp(u(k))
          rho(k+nr1)=0.d0
          sr=sr+rho(k)
        enddo
        du0=-log(sr)
        u(0)=ulast(0)-du(0)*fact+du0
        sdr=0.d0
        do k=1,nr1
          rho(k)=rho(k)/sr
          u(k)=u(k)+du0
        enddo
      enddo
      irtc=-1
      return
      end

      subroutine csrtrack(z,dp,np,zl,nr,dz,aw,r56,r65,
     $     frf,phirf,
     $     damp,sigmae,data,ndiv,nt)
      implicit none
      real*8 dpmax
      parameter (dpmax=0.1d0)
      integer*4 np,nr,nt,i,k,m,nr1,ndiv,n,nc
      real*8 z(np),dp(np),zl(nr),r56,r65,damp,sigmae,dz,aw,data(2,nt)
      real*8 rho(nr),zc,ddrho(nr),work(nr),rand(np),
     $     a,b,c,d,w,ad,bd,sz,sdp,vz,vdp,r56n,zk,awn,frf,phirf,
     $     wr,vr,dpr,cs,pi,zlim
      parameter (pi=3.14159265358979324d0,cs=299792458.d0)
      awn=aw/ndiv
      ad=1.d0-2.d0*damp/ndiv
      bd=sqrt(4.d0*damp/ndiv*(1.d0-damp/ndiv))*sigmae
      zc=dz*nr/4.d0
      nr1=nr/2
      r56n=r56/ndiv
      if(frf /= 0.d0)then
        wr=2.d0*pi*frf/cs
        vr=r65/wr
        dpr=vr*sin(phirf)
        zlim=0.5d0*cs/frf
      else
        wr=0.d0
        vr=0.d0
        dpr=0.d0
        zlim=1.d100
      endif
      nc=(ndiv-1)/2+1
      do i=1,nt
        do n=1,ndiv
          do k=1,np
            z(k)=z(k)+r56n*dp(k)
          enddo          
          call csrrho(z,rho,dz,np,nr)
          call csrwake(rho,zl,nr)
          call spline1(nr1,rho,ddrho,work,0,0)
          call tgauss_array(rand,np)
          sz=0.d0
          sdp=0.d0
          if(n .eq. nc)then
            if(n .eq. ndiv)then
              sz=0.d0
              sdp=0.d0
              do k=1,np
                zk=(z(k)+zc)/dz
                m=int(zk)
                if(m .gt. 0 .and. m .lt. nr1)then
                  b=zk-m
                  a=1.d0-b
                  c=-a*b*(a+1.d0)
                  d=-a*b*(b+1.d0)
                  w=a*rho(m)+b*rho(m+1)+c*ddrho(m)+d*ddrho(m+1)
                  if(frf .eq. 0.d0)then
                    dp(k)=(dp(k)+awn*w+r65*z(k))*ad+bd*rand(k)
                  else
                    dp(k)=(dp(k)+awn*w+vr*sin(wr*z(k)+phirf)-dpr)*ad
     $                   +bd*rand(k)
                  endif
                else
                  if(frf .eq. 0.d0)then
                    dp(k)=(dp(k)+r65*z(k))*ad+bd*rand(k)
                  else
                    dp(k)=(dp(k)+vr*sin(wr*z(k)+phirf)-dpr)*ad
     $                   +bd*rand(k)
                  endif
                endif
              enddo
            else
              do k=1,np
                zk=(z(k)+zc)/dz
                m=int(zk)
                if(m .gt. 0 .and. m .lt. nr1)then
                  b=zk-m
                  a=1.d0-b
                  c=-a*b*(a+1.d0)
                  d=-a*b*(b+1.d0)
                  w=a*rho(m)+b*rho(m+1)+c*ddrho(m)+d*ddrho(m+1)
                  if(frf .eq. 0.d0)then
                    dp(k)=(dp(k)+awn*w+r65*z(k))*ad+bd*rand(k)
                  else
                    dp(k)=(dp(k)+awn*w+vr*sin(wr*z(k)+phirf)-dpr)*ad
     $                   +bd*rand(k)
                  endif
                else
                  if(frf .eq. 0.d0)then
                    dp(k)=(dp(k)+r65*z(k))*ad+bd*rand(k)
                  else
                    dp(k)=(dp(k)+vr*sin(wr*z(k)+phirf)-dpr)*ad
     $                   +bd*rand(k)
                  endif
                endif
              enddo
            endif
          else
            if(n .eq. ndiv)then
              sz=0.d0
              sdp=0.d0
              do k=1,np
                zk=(z(k)+zc)/dz
                m=int(zk)
                if(m .gt. 0 .and. m .lt. nr1)then
                  b=zk-m
                  a=1.d0-b
                  c=-a*b*(a+1.d0)
                  d=-a*b*(b+1.d0)
                  w=a*rho(m)+b*rho(m+1)+c*ddrho(m)+d*ddrho(m+1)
                  dp(k)=(dp(k)+awn*w)*ad+bd*rand(k)
                else
                  dp(k)=dp(k)*ad+bd*rand(k)
                endif
              enddo
            else
              do k=1,np
                zk=(z(k)+zc)/dz
                m=int(zk)
                if(m .gt. 0 .and. m .lt. nr1)then
                  b=zk-m
                  a=1.d0-b
                  c=-a*b*(a+1.d0)
                  d=-a*b*(b+1.d0)
                  w=a*rho(m)+b*rho(m+1)+c*ddrho(m)+d*ddrho(m+1)
                  dp(k)=(dp(k)+awn*w)*ad+bd*rand(k)
                else
                  dp(k)=dp(k)*ad+bd*rand(k)
                endif
              enddo
            endif
          endif
        enddo
        sz=0.d0
        sdp=0.d0
        vz=0.d0
        vdp=0.d0
        if(frf .eq. 0.d0)then
          do k=1,np
            sz=sz+z(k)
            sdp=sdp+dp(k)
          enddo
          sz=sz/np
          sdp=sdp/np
          do k=1,np
            dp(k)=dp(k)-sdp
            z(k)=z(k)-sz
            vz=vz+z(k)**2
            vdp=vdp+dp(k)**2
          enddo
        else
          do k=1,np
            if(dp(k) .lt. -dpmax)then
              z(k)=-zlim*10.d0
            elseif(z(k) .gt. -zlim .and. z(k) .lt. zlim)then
              sz=sz+z(k)
              sdp=sdp+dp(k)
            endif
          enddo
          sz=sz/np
          sdp=sdp/np
          do k=1,np
            if(z(k) .gt. -zlim .and. z(k) .lt. zlim)then
              vz=vz+(z(k)-sz)**2
              vdp=vdp+(dp(k)-sdp)**2
            endif
          enddo
        endif
        data(1,i)=sqrt(vz/np)
        data(2,i)=sqrt(vdp/np)
        if(data(2,i) .gt. dpmax
     $       .or. data(1,i) .gt. zc)then
          do k=i+1,nt
            data(1,k)=0.d0
            data(2,k)=0.d0
          enddo
          return
        endif
      enddo
      return
      end

      subroutine csrwake(rho,zl,nr)
      use ftr
      implicit none
      integer*4 ,intent(in):: nr
      integer*4 i,i2
      real*8 ,intent(inout):: rho(nr)
      real*8 ,intent(in):: zl(nr)
      real*8 ri2
      call trftr(rho,nr,.true.)
      do i=2,nr/2
        i2=i*2-1
        ri2=rho(i2)
        rho(i2)  =zl(i2)*ri2      -zl(i2+1)*rho(i2+1)
        rho(i2+1)=zl(i2)*rho(i2+1)+zl(i2+1)*ri2
      enddo
      rho(1)=zl(1)*rho(1)
      rho(2)=zl(2)*rho(2)
      call tftrr(rho,nr,.false.)
      rho=rho/nr
      return
      end

      subroutine csrrho(z,rho,dz,np,nr)
      implicit none
      integer*4 ,intent(in):: np,nr
      integer*4 ic,i,l
      real*8 ,intent(in):: z(np)
      real*8 ,intent(out):: rho(nr)
      real*8 dz,a,rmin,rmax,r,zc
      ic=nr/4
      zc=dz*ic
      rmin=1.d0
      rmax=nr/2-1.d-10
      rho=0.d0
      do i=1,np
        r=(z(i)+zc)/dz
        if(r .ge. rmin .and. r .le. rmax)then
          l=int(r)
          a=r-l
          rho(l)=rho(l)+1.d0-a
          rho(l+1)=rho(l+1)+a
        endif
      enddo
      rho=rho/dz
      return
      end
