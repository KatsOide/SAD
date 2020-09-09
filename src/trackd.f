c     CAUTION: kptbl(#,3) MUST be `0' before trackd() called
      subroutine trackd(range,r1,n1,nturn,
     $     trval,phi,dampr,dampenough,ivar1,ivar2,lfno)
      use tfstk
      use ffs_flag
      use tmacro,only:np0,codin,dvfs,omega0,taurdx,taurdy,taurdz,
     $     nparallel
      use ffs_pointer, only:idelc
      use tftr
      use tfshare
      use macmath
      use iso_c_binding
      implicit none
      integer*4 ,parameter ::n1p0=256,n2p=51,maxturn=2**29,
     $     maxpara=256,nw=16,nkptbl = 6, minnp=16
      real*8, parameter :: big=1.d300
      character*15 ,parameter ::vname='`ResultOfDAPERT'
      character ,parameter ::label(3)=['X','Y','Z']
      integer*4 ,intent(in):: ivar1,ivar2,lfno,n1,nturn
      real*8 ,intent(in):: dampr,dampenough,range(3,3),r1(n1)
      type (sad_descriptor) kv,kxm
      type (sad_symdef),pointer::symd
      type (sad_dlist),pointer::kal,kal1,kal2,kal2i
      type (sad_rlist),pointer::kal11,kal12,kal13,kal2iv
      integer*4 ,pointer ::ntloss(:,:),kzx(:,:)
      integer*4 ,allocatable::kptbl(:,:),mturn(:)
      integer*4 n1p,j,n,jzout,np1,k,np,kp,kx,kz,irw,nsc,iw,
     $     jj,ip,isw,kseed,npmax,npara,nxm(n1p0),np00,
     $     muls,irtc,i,fork_worker,waitpid,ncons,nscore,ivar3
      integer*8 intlm
      real*8 ,allocatable ::x(:),px(:),y(:),py(:),z(:),g(:),dv(:),
     $     spx(:),spy(:),spz(:),aenox(:),aenoy(:),aenoz(:)
      real*8 a2min,a2max,a3min,a3max,a1min,a1max,a2step,a3step,a1step,
     $     emx,emz,rgetgl1,cx,cy,cz,sx,sy,sz,dampx,dampy,dampz,t0,
     $     trval,phi(3),a1i(n1p0)
      character*12 autos
      logical*4 damp,ini,remain,pol0
      character rad62a
c     begin initialize for preventing compiler warning
      ipr=0
c     end   initialize for preventing compiler warning
      damp=dampr .ne. 0.d0
c      write(*,*)'trackd ',damp
      if(damp)then
        t0=m_2_pi/omega0
        dampx=exp(-t0/taurdx)
        dampy=exp(-t0/taurdy)
        dampz=exp(-t0/taurdz)
      else
        dampx=1.d0
        dampy=1.d0
        dampz=1.d0
      endif
c      write(*,*)'trackd0 ',damp,dampx,t0,omega0,taurdx
      ncons=max(1,nint(rgetgl1('DAPWIDTH')))
      trpt=.false.
      pol0=calpol
      calpol=.false.
      call tsetdvfs
      trval=0.d0
      nscore=0
      muls=merge((nturn/6200+1)*100,10,nturn .gt. 600)
c      write(*,*)'trackd-muls: ',nturn,muls,nturn/6200
      a2min=range(1,1)
      a2max=range(2,1)
      a3min=range(1,2)
      a3max=range(2,2)
      a1min=minval(r1(1:2))
      a1max=maxval(r1(1:2))
      if(a1max .eq. a1min)then
        a1min=0.d0
      endif
      a2step=(a2max-a2min)/(n2p-1)
      a3step=(a3max-a3min)/(n2p-1)
      if(n1 .eq. 0)then
        a1step=abs(r1(3))
        if(a1step .eq. 0.d0 .or. a1min .eq. a1max)then
          n1p=1
        else
          n1p=min(n1p0,int(max(2.d0,(a1max-a1min)/a1step+1.1d0)))
          a1step=(a1max-a1min)/(n1p-1)
        endif
        do i=1,n1p
          a1i(i)=a1min+a1step*(i-1)
        enddo
      else
        n1p=min(n1p0,n1)
        a1i(1:n1p)=r1(1:n1p)
        a1min=minval(a1i(1:n1p))
        a1max=maxval(a1i(1:n1p))
      endif
      intlm=ktfallocshared((n1p*n2p+2)/2)
      call c_f_pointer(c_loc(klist(intlm)),ntloss,[n1p,n2p])
      ivar3=6-ivar1-ivar2
      write(lfno,'(
     $   '' Range    '',a,''min:'',f8.3,'' '',a,''max:'',f8.3,/
     $   ''         ('',a,''min:'',f8.3,'' '',a,''max:'',f8.3,'')'',/
     $   ''          '',a,''min:'',f8.3,'' '',a,''max:'',f8.3,/
     $   '' Display: '',i5,'' turns/character'')')
     $     label(ivar2),a2min,label(ivar2),a2max,
     $     label(ivar3),a3min,label(ivar3),a3max,
     $     label(ivar1),a1min,label(ivar1),a1max,muls
      write(lfno,*)
     $'     N'//label(ivar1)//
     $     '     0----|----1----|----2----|----3----|----4----|----5'
      kv=kxsymbolv(vname,len(vname),symd)
c      kv=ktfsymbolz(vname,len(vname))-4
      call tflocald(symd%value)
      symd%value=kxadaloc(0,2,kal)
      kal%dbody(1)=kxadaloc(0,3,kal1)
      kal1%dbody(1)=kxavaloc(0,2,kal11)
      kal11%rbody(1)=a2min
      kal11%rbody(2)=a2max
      kal1%dbody(2)=kxavaloc(0,2,kal12)
      kal12%rbody(1)=a3min
      kal12%rbody(2)=a3max
      kal1%dbody(3)=kxavaloc(0,2,kal13)
      kal13%rbody(1)=a1min
      kal13%rbody(2)=a1max
      kal%dbody(2)=kxadalocnull(0,n1p,kal2)
c      lp0=latt(1)+kytbl(kwmax,idtype(idelc(1)))+1
      emx=sqrt(abs(rgetgl1('EMITX'))+abs(rgetgl1('EMITY')))
      emz=merge(sqrt(abs(rgetgl1('EMITZ'))),abs(rgetgl1('SIGE')),rfsw)
      ntloss(1:n1p,1:n2p)=maxturn
      nxm(1:n1p)=n2p+1
      npara=min(nparallel,nprmax,np0/minnp)
      npr=0
      ipn=0
c      write(*,*)'trackd-npara ',npara,np0
      if(npara .gt. 1)then
        kseed=0
        do while(npr .lt. npara-1)
          npri=npr
          iprid=fork_worker()
          if(iprid .eq. 0)then
            call tfaddseed(kseed,irtc)
            if(irtc .ne. 0)then
              write(*,*)'addseed-error ',irtc
              call exit_without_hooks(0)
            endif
            exit
          endif
          npr=npr+1
          ipr(npr)=iprid
c          write(*,*)'trackd-ipr ',npara,npr,iprid,muls
        enddo
        npr=npr+1
        npri=npr
        npmax=max(1,min(n1p*ncons/npara,np0))
      else
        iprid=getpid()
        npri=0
        npmax=max(1,min(n1p*ncons,np0))
      endif
      npz=npmax
      np00=np0
      np0=npz
      allocate (x(npz))
      allocate (px(npz))
      allocate (y(npz))
      allocate (py(npz))
      allocate (z(npz))
      allocate (g(npz))
      allocate (dv(npz))
      allocate (spx(npz))
      allocate (spy(npz))
      allocate (spz(npz))
      allocate (aenox(npz))
      allocate (aenoy(npz))
      allocate (aenoz(npz))
      call tfevals('`ExtMap$@InitMap['//autos(dble(npz))//',1]',
     $     kxm,irtc)
      allocate(kptbl(npmax,6))
      allocate(mturn(npmax))
      allocate(kzx(2,npmax))
c      if(muls .ne. 10)then
c        write(*,*)'trackd-5 ',npri,muls
c      endif
      do i=1,npmax
        kptbl(i,1)=i
        kptbl(i,2)=i
      enddo
      kzx(1,1:npmax)=0
      kzx(2,1:npmax)=0
c      if(muls .ne. 10)then
c        write(*,*)'trackd-6 ',npri,muls
c      endif
      n=1
      jzout=1
      loop_1: do
        np1=npmax
        iw=nw
        remain=.true.
        cx=cos(phi(1))
        sx=sin(phi(1))
        cy=cos(phi(2))
        sy=sin(phi(2))
        cz=cos(phi(3))
        sz=sin(phi(3))
        LOOP_K: do k=1,npmax
          if(kzx(1,k) .eq. 0)then
            ip=kptbl(k,1)
            if(ip .le. np1)then
              if(remain)then
                do i=1,n1p
                  do j=nxm(i)-1,1,-1
                    if(j .lt. nxm(i)-ncons)then
                      exit
                    elseif(ntloss(i,j) .lt. nturn)then
                      nxm(i)=j
                    elseif(ntloss(i,j) .eq. maxturn)then
                      ntloss(i,j)=nturn
                      kzx(1,k)=i
                      kzx(2,k)=j
                      mturn(k)=0
                      select case(ivar1)
                      case (1)
                        x(ip)=a1i(i)
                        px(ip)=0.d0
                      case (2)
                        y(ip)=a1i(i)
                        py(ip)=0.d0
                      case default
                        z(ip)=0.d0
                        g(ip)=a1i(i)
                      end select
                      select case(ivar2)
                      case (1)
                        x(ip)=a2step*(j-1)+a2min
                        px(ip)=-x(ip)*sx
                        x(ip)=x(ip)*cx
                      case (2)
                        y(ip)=a2step*(j-1)+a2min
                        py(ip)=-y(ip)*sy
                        y(ip)=y(ip)*cy
                      case default
                        z(ip)=a2step*(j-1)+a2min
                        g(ip)=-z(ip)*sz
                        z(ip)=z(ip)*cz
                      end select
                      select case(ivar3)
                      case (1)
                        x(ip)=a3step*(j-1)+a3min
                        px(ip)=-x(ip)*sx
                        x(ip)=cx*x(ip)
                      case (2)
                        y(ip)=a3step*(j-1)+a3min
                        py(ip)=-y(ip)*sy
                        y(ip)=cy*y(ip)
                      case default
                        z(ip)=a3step*(j-1)+a3min
                        g(ip)=-z(ip)*sz
                        z(ip)=cz*z(ip)
                      end select
                      aenox(k)=(x(ip)**2+px(ip)**2)*emx**2*dampenough
                      aenoy(k)=(y(ip)**2+py(ip)**2)*emx**2*dampenough
                      aenoz(k)=(z(ip)**2+ g(ip)**2)*emz**2*dampenough
                      select case (ivar3)
                      case (1)
                        aenox(k)=big
                      case (2)
                        aenoy(k)=big
                      case default
                        aenoz(k)=big
                      end select
c     Reinit kptbl(ip,4) to reuse particle array slot `ip'
                      kptbl(ip,4)=0
                      call tinip1(x(ip),px(ip),y(ip),py(ip),
     $                     z(ip),g(ip),dv(ip),
     $                     emx,emz,codin,dvfs)
c     write(*,'(a,7i5,1p3g15.7)')
c     $                   ' trackd-Launch ',npr1,i,j,ivar1,ivar2,ivar3,
c     $                     nxm(i),x(ip),y(ip),py(ip)
                      cycle LOOP_K
                    endif
                  enddo
                enddo
              endif
              remain=.false.
c     Swap particle k <-> j[array index ip <-> np1]
              j=kptbl(np1,2)
c     - Update maps between partice ID and array index
              kptbl(j,  1)=ip
              kptbl(k,  1)=np1
              kptbl(ip, 2)=j
              kptbl(np1,2)=k
c     - Overwrite slot[np1] to slot[ip](Drop particle[k] information)
              kptbl(ip,3:nkptbl) = kptbl(np1,3:nkptbl)
              x(ip)=x(np1)
              px(ip)=px(np1)
              y(ip)=y(np1)
              py(ip)=py(np1)
              z(ip)=z(np1)
              g(ip)=g(np1)
              dv(ip)=dv(np1)
              np1=np1-1
              if(np1 .le. 0)then
                exit LOOP_K
              endif
            endif
          endif
        enddo LOOP_K
        np=np1
        do while(np .gt. 0)
c     write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-1 ',n,np,(kptbl(i,1),y(i),i=1,14)
          call tturn(np,x,px,y,py,z,g,dv,spx,spy,spz,kptbl,n)
c     if(muls .ne. 10)then
c     write(*,*)'trackd-tturn ',npri,muls
c     endif
c     write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-2 ',n,np,(kptbl(i,1),y(i),i=1,14)
          if(damp .or. dampenough .ne. 0.d0)then
            call tpdamp(np,x,px,y,py,z,g,dv,dampx,dampy,dampz,damp,
     $           aenox,aenoy,aenoz,kptbl(1,2),mturn)
          endif
          n=n+1
          ini=.false.
          do i=1,npmax
            if(kzx(1,i) .le. 0)then
              cycle
            endif
            kp=kptbl(i,1)
            if(kp .le. np)then
              mturn(i)=mturn(i)+1
              if(mturn(i) .ge. nturn)then
                kz=kzx(1,i)
                kx=kzx(2,i)
                ntloss(kz,kx)=nturn
                kzx(1,i)=0
                ini=.true.
              endif
            else
              kz=kzx(1,i)
              kx=kzx(2,i)
c     if(kz .eq. 1)then
c     write(*,'(a,1x,8i10)')'trackd-Lost: ',
c     $           kz,kx,mturn(i),i,np,np1,kp,npri
c     endif
              ntloss(kz,kx)=mturn(i)
              kzx(1,i)=0
              ini=.true.
            endif
          enddo
          if(ini)then
            cycle loop_1
          elseif(np1 .lt. npmax)then
            iw=iw-1
            if(iw .le. 0)then
              cycle loop_1
            endif
          endif
        enddo
        exit
      enddo loop_1
      call tfevals('`ExtMap$@ResetMap[]',kxm,irtc)
      if(iprid .eq. 0)then
c        write(*,*)'trackd-stop ',npr1
c        stop
        call tfresetsharedmap()
        call exit_without_hooks(0)
      endif
      np0=np00
c      write(*,*)'trackd-wait ',npr,n1p,n2p,np0
      loop_j: do j=1,npr-1
        irw=waitpid(-1,isw)
c        write(*,*)'trackd-wait-j ',j,irw
        do k=1,npr-1
          if(irw .eq. ipr(k))then
c            write(*,*)'trackd-wait-k ',k
            n=k
            ipr(k)=0
            cycle loop_j
          endif
        enddo
        write(*,*)'???trackd-Strange child: ',irw,isw
      enddo loop_j
      do i=1,n1p
        nsc=n2p
        do k=1,n2p
          if(ntloss(i,k) .lt. nturn)then
            nsc=k-1
            exit
          endif
        enddo
        nscore=nscore+nsc
c        write(*,*)'trackd: ',lfno,i,muls,
c     $       rad62a(ntloss(i,1)/muls,1)
        write(lfno,'(1x,f8.2,i3,1x,51a1)')a1i(i),nsc,
     $       (rad62a(ntloss(i,jj)/muls,jj),jj=1,n2p)
        kal2%dbody(i)=kxadaloc(0,3,kal2i)
        kal2i%rbody(1)=a1i(i)
        kal2i%rbody(2)=dble(nsc)
        kal2i%dbody(3)=kxavaloc(0,n2p,kal2iv)
        kal2iv%rbody(1:n2p)=dble(ntloss(i,1:n2p))
      enddo
      call tfreeshared(intlm)
      write(lfno,*)
     $'     N'//label(ivar1)//
     $     '     0----|----1----|----2----|----3----|----4----|----5'
      write(lfno,'(a,i5)')'    Score: ',nscore
c      call tfdebugprint(symd%value,'res of DA',10)
c      call tfevals('Print['//vname//']',kx,irtc)
      trval=nscore
      calpol=pol0
      return
      end

      character function rad62a(m,j)
      implicit none
      integer*4 ,intent(in):: m,j
      character rad62
      if(m .eq. 0)then
        rad62a=' '
        if(mod(j,10) .eq. 1 .and. j .ne. 1)then
          rad62a='.'
        endif
      else
        rad62a=rad62(m)
      endif
      return
      end

      subroutine tinip1(x,px,y,py,z,g,dv,emx,emz,codin,dvfs)
      implicit none
      real*8 ,intent(inout):: x,px,y,py,z,dv,g
      real*8 ,intent(in):: codin(6),dvfs,emx,emz
      real*8 xa(8)
      xa(1)=x *emx
      xa(2)=px*emx
      xa(3)=y *emx
      xa(4)=py*emx
      xa(5)=z *emz
      xa(6)=g *emz
      call tmap(xa,xa,-1)
      xa(1:6)=xa(1:6)+codin
      call tconv(xa,xa,-1)
      x =xa(1)
      px=xa(2)
      y =xa(3)
      py=xa(4)
      z =xa(5)
      g =xa(6)
      dv=xa(7)+dvfs
      return
      end

      subroutine tinip(np,x,px,y,py,z,g,dv,emx,emz,codin,dvfs,cmplot)
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np)
      real*8 xa(8),tgauss,dxi,dyi,dzi,emx,emz,codin(6),dvfs,rgetgl1
      logical*4 cmplot
      dxi=rgetgl1('TDXI')
      dyi=rgetgl1('TDYI')
      dzi=rgetgl1('TDZI')
      do i=1,np
        xa(1)=x(i) *emx
        xa(2)=px(i)*emx
        xa(3)=y(i) *emx
        xa(4)=py(i)*emx
        xa(5)=z(i) *emz
        xa(6)=g(i) *emz
        if(cmplot)then
          xa(2)=xa(1)*tgauss()
          xa(1)=xa(1)*tgauss()
          xa(4)=xa(3)*tgauss()
          xa(3)=xa(3)*tgauss()
          xa(5)=xa(6)*tgauss()
          xa(6)=xa(6)*tgauss()
        endif
        call tmap(xa,xa,-1)
        xa(1:6)=xa(1:6)+codin
c        call tadd(xa,codin,xa,6)
        call tconv(xa,xa,-1)
        x(i) =xa(1)+dxi
        px(i)=xa(2)
        y(i) =xa(3)+dyi
        py(i)=xa(4)
        z(i) =xa(5)+dzi
        g(i) =xa(6)
        dv(i)=xa(7)+dvfs
      enddo
c      write(*,'(a,1p6g15.7)')'tinip ',xa(6),emx,emz
      return
      end

      subroutine tpdamp(np,x,px,y,py,z,g,dv,dampx,dampy,dampz,
     $     damp,aenox,aenoy,aenoz,kptbl,mturn)
      use tfstk
      use tmacro
      implicit none
      integer*4 ,intent(in):: np,kptbl(np)
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),
     $     dv(np),g(np)
      real*8 xa(7)
      real*8 ,intent(in):: dampx,dampy,dampz,
     $     aenox(np),aenoy(np),aenoz(np)
      integer*4 ,intent(inout):: mturn(np)
      integer*4 i,j
      logical*4 damp
      call tconvm(np,px,py,g,dv,1)
      do i=1,np
        xa(1)=x(i)-codin(1)
        xa(2)=px(i)-codin(2)
        xa(3)=y(i)-codin(3)
        xa(4)=py(i)-codin(4)
        xa(5)=z(i)-codin(5)
        xa(6)=g(i)-codin(6)
        call tmap(xa,xa,1)
        j=kptbl(i)
        if(  xa(1)**2+xa(2)**2 .lt. aenox(j) .and.
     $       xa(3)**2+xa(4)**2 .lt. aenoy(j) .and.
     $       xa(5)**2+xa(6)**2 .lt. aenoz(j))then
c          write(*,'(a,2i5,1p6g15.7)')'tpdamp ',j,i,
c     $xa(1),xa(2),xa(5),xa(6),aenox(j),aenoz(j)
          mturn(j)=nturn-1
c          write(*,*)'tpdamp ',j,nturn-1
        endif
        if(damp)then
          xa(1)=xa(1)*dampx
          xa(2)=xa(2)*dampx
          xa(3)=xa(3)*dampy
          xa(4)=xa(4)*dampy
          xa(5)=xa(5)*dampz
          xa(6)=xa(6)*dampz
          call tmap(xa,xa,-1)
          x(i) =xa(1)+codin(1)
          px(i)=xa(2)+codin(2)
          y(i) =xa(3)+codin(3)
          py(i)=xa(4)+codin(4)
          z(i) =xa(5)+codin(5)
          g(i) =xa(6)+codin(6)
        endif
      enddo
      call tconvm(np,px,py,g,dv,-1)
      return
      end
