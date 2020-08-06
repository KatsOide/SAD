c     CAUTION: kptbl(#,3) MUST be `0' before trackd() called
      subroutine trackd(latt,kptbl,x,px,y,py,z,g,dv,pz,
     $     kzx,mturn,trval,phi,damp,dampenough,ivar1,ivar2,lfno)
      use iso_c_binding
      use tfstk
      use ffs_flag
      use tfshare
      use tmacro
      implicit none
      integer*4 n1p0,nxp
      parameter (n1p0=256,nxp=51)
      integer*4, contiguous, pointer, dimension(:,:) :: ntloss
      integer*4 kptbl(np0,6),mturn(np0),kzx(2,np0)
      integer*8 intlm,latt(nlat)
      integer*4 ivar1,ivar2,lfno
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 trval,phi(3),damp,dampenough
      intlm=ktfallocshared((n1p0*nxp+2)/2)
      call c_f_pointer(c_loc(klist(intlm)),ntloss,[n1p0,nxp])
      call trackd0(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     mturn,kzx,trval,phi,
     $     damp .ne. 0.d0,dampenough,ntloss,ivar1,ivar2,lfno)
      call tfreeshared(intlm)
c      if(mapfree(ntloss(intlm+1)) .ne. 0)then
c        write(*,*)'???trackd-error in munmap.'
c      endif
      return
      end

      subroutine trackd0(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     mturn,kzx,trval,phi,
     $     damp,dampenough,ntloss,ivar1,ivar2,lfno)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:idelc
      implicit none
      integer*4 n1p0,n2p,maxturn,maxpara,nw,lfno,ncons,nscore,
     $     ivar1,ivar2,ivar3
      parameter (n1p0=200,n2p=51,maxturn=2**29,maxpara=256,nw=16)
      integer*4, parameter :: nkptbl = 6, minnp=16
      integer*8 kv,kax,kax11,kax12,kax13,kax2,
     $     kaxi,kaxi3,kax1,latt(nlat)
      integer*4 kptbl(np0,nkptbl),n1p,npr1,
     $     ipr,j,n,jzout,np1,k,np,kp,kz,kx,irw,nsc,iw,
     $     jj,ip,isw,kseed,npmax,npara,nxm(n1p0)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     1     pz(np0),dampenough,aenox(np0),aenoy(np0),aenoz(np0),
     $     spx(np0),spy(np0),spz(np0)
      real*8 a2min,a2max,a3min,a3max,a1min,a1max,a2step,a3step,a1step,
     $     emx,emz,rgetgl1,cx,cy,cz,sx,sy,sz,
     $     dampx,dampy,dampz,t0
      real*8 trval,phi(3)
      real*8 a1i(n1p0)
      real*8, parameter :: big=1.d300
      character*14 vname
      data vname /'ResultOfDAPERT'/
      integer*4 ntloss(n1p0,n2p)
      integer*4 mturn(np0),kzx(2,np0),muls,irtc
      integer*4 i,fork_worker,wait,ichpid(maxpara)
      logical*4 ini,remain,damp,pol0
      character label(3)
      data label /'X','Y','Z'/
      character rad62a
c     begin initialize for preventing compiler warning
      ipr=0
c     end   initialize for preventing compiler warning
      if(damp)then
        t0=pi2/omega0
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
      if(nturn .gt. 600)then
        muls=(nturn/6200+1)*100
      else
        muls=10
      endif
      a2min=x(1)
      a2max=x(2)
      a3min=y(1)
      a3max=y(2)
      a1min=min(g(1),g(2))
      a1max=max(g(1),g(2))
      if(a1max .eq. a1min)then
        a1min=0.d0
      endif
      a2step=(a2max-a2min)/(n2p-1)
      a3step=(a3max-a3min)/(n2p-1)
      if(pz(1) .eq. 0.d0)then
        a1step=abs(g(3))
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
        n1p=min(n1p0,int(pz(1)))
        a1i(1:n1p)=g(1:n1p)
        a1min=minval(a1i(1:n1p))
        a1max=maxval(a1i(1:n1p))
      endif
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
      kv=ktfsymbolz(vname,len(vname))-4
      call tflocal(klist(kv))
      kax=ktadaloc(0,2)
      klist(kv)=ktflist+kax
      kax1=ktadaloc(0,3)
      klist(kax+1)=ktflist+kax1
      kax11=ktavaloc(0,2)
      klist(kax1+1)=ktflist+kax11
      rlist(kax11+1)=a2min
      rlist(kax11+2)=a2max
      kax12=ktavaloc(0,2)
      klist(kax1+2)=ktflist+kax12
      rlist(kax12+1)=a3min
      rlist(kax12+2)=a3max
      kax13=ktavaloc(0,2)
      klist(kax1+3)=ktflist+kax13
      rlist(kax13+1)=a1min
      rlist(kax13+2)=a1max
      kax2=ktadalocnull(0,n1p)
      klist(kax+2)=ktflist+kax2
c      lp0=latt(1)+kytbl(kwmax,idtype(idelc(1)))+1
      emx=sqrt(abs(rgetgl1('EMITX'))+abs(rgetgl1('EMITY')))
      if(rfsw)then
        emz=sqrt(abs(rgetgl1('EMITZ')))
      else
        emz=abs(rgetgl1('SIGE'))
      endif
      ntloss(1:n1p,1:n2p)=maxturn
      nxm(1:n1p)=n2p+1
      npara=min(nparallel,maxpara,np0/minnp)
      ipr=1
      if(npara .gt. 1)then
        kseed=0
        npr1=1
        do while(npr1 .lt. npara .and. ipr .gt. 0)
          kseed=kseed+2
          ipr=fork_worker()
          if(ipr .gt. 0)then
            ichpid(npr1)=ipr
            npr1=npr1+1
            call tfaddseed(kseed,irtc)
            if(irtc .ne. 0)then
              go to 3000
            endif
          endif
        enddo
        npmax=max(1,min(n1p*ncons/npara,np0))
      else
        npr1=1
        npmax=np0
      endif
      kptbl=0
      do i=1,npmax
        kptbl(i,1)=i
        kptbl(i,2)=i
      enddo
      kzx(1,1:npmax)=0
      kzx(2,1:npmax)=0
      n=1
      jzout=1
 1    np1=npmax
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
                    if(ivar1 .eq. 1)then
                      x(ip)=a1i(i)
                      px(ip)=0.d0
                    elseif(ivar1 .eq. 2)then
                      y(ip)=a1i(i)
                      py(ip)=0.d0
                    else
                      z(ip)=0.d0
                      g(ip)=a1i(i)
                    endif
                    if(ivar2 .eq. 1)then
                      x(ip)=a2step*(j-1)+a2min
                      px(ip)=-x(ip)*sx
                      x(ip)=x(ip)*cx
                    elseif(ivar2 .eq. 2)then
                      y(ip)=a2step*(j-1)+a2min
                      py(ip)=-y(ip)*sy
                      y(ip)=y(ip)*cy
                    else
                      z(ip)=a2step*(j-1)+a2min
                      g(ip)=-z(ip)*sz
                      z(ip)=z(ip)*cz
                    endif
                    if(ivar3 .eq. 1)then
                      x(ip)=a3step*(j-1)+a3min
                      px(ip)=-x(ip)*sx
                      x(ip)=cx*x(ip)
                    elseif(ivar3 .eq. 2)then
                      y(ip)=a3step*(j-1)+a3min
                      py(ip)=-y(ip)*sy
                      y(ip)=cy*y(ip)
                    else
                      z(ip)=a3step*(j-1)+a3min
                      g(ip)=-z(ip)*sz
                      z(ip)=cz*z(ip)
                    endif
                    aenox(k)=(x(ip)**2+px(ip)**2)*emx**2*dampenough
                    aenoy(k)=(y(ip)**2+py(ip)**2)*emx**2*dampenough
                    aenoz(k)=(z(ip)**2+ g(ip)**2)*emz**2*dampenough
                    if(ivar3 .eq. 1)then
                      aenox(k)=big
                    elseif(ivar3 .eq. 2)then
                      aenoy(k)=big
                    else
                      aenoz(k)=big
                    endif
c     Reinit kptbl(ip,4) to reuse particle array slot `ip'
                    kptbl(ip,4)=0
                    call tinip(1,
     $                   x(ip),px(ip),y(ip),py(ip),
     $                   z(ip),g(ip),dv(ip),
     $                   emx,emz,codin,dvfs,.false.)
c                    write(*,'(a,7i5,1p3g15.7)')
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
              go to 3000
            endif
          endif
        endif
      enddo LOOP_K
      np=np1
 101  continue
c      write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-1 ',n,np,(kptbl(i,1),y(i),i=1,14)
      call tturn(np,x,px,y,py,z,g,dv,spx,spy,spz,kptbl,n)
c      write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-2 ',n,np,(kptbl(i,1),y(i),i=1,14)
      if(damp .or. dampenough .ne. 0.d0)then
        call tpdamp(np,x,px,y,py,z,g,dv,dampx,dampy,dampz,damp,
     $       aenox,aenoy,aenoz,kptbl(1,2),mturn)
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
c          write(*,'(a,1x,7i10)')'trackd-Lost ',
c     $         kz,kx,mturn(i),i,np,np1,kp
          ntloss(kz,kx)=mturn(i)
          kzx(1,i)=0
          ini=.true.
        endif
      enddo
      if(ini)then
        go to 1
      elseif(np1 .lt. npmax)then
        iw=iw-1
        if(iw .le. 0)then
          go to 1
        endif
      endif
      go to 101
 3000 continue
      if(ipr .eq. 0)then
c        write(*,*)'trackd-stop ',npr1
c        stop
        call tfresetsharedmap()
        call exit_without_hooks(0)
      endif
      do j=1,npr1-1
 3010   irw=wait(isw)
        do k=1,npr1-1
          if(irw .eq. ichpid(k))then
            n=k
            ichpid(k)=0
            go to 3020
          endif
        enddo
        write(*,*)'???trackd-Strange child: ',irw,isw
        go to 3010
 3020   continue
      enddo
      do i=1,n1p
        nsc=0
        do k=1,n2p
          if(ntloss(i,k) .lt. nturn)then
            nsc=k-1
            go to 431
          endif
        enddo
        nsc=n2p
 431    nscore=nscore+nsc
        write(lfno,'(1x,f8.2,i3,1x,51a1)')a1i(i),nsc,
     $       (rad62a(ntloss(i,jj)/muls,jj),jj=1,n2p)
        kaxi=ktadaloc(0,3)
        klist(kax2+i)=ktflist+kaxi
        rlist(kaxi+1)=a1i(i)
        rlist(kaxi+2)=dble(nsc)
        kaxi3=ktraaloc(0,n2p)
        klist(kaxi+3)=ktflist+kaxi3
        rlist(kaxi3+1:kaxi3+n2p)=ntloss(i,1:n2p)
c        do jj=1,n2p
c          rlist(kaxi3+jj)=ntloss(i,jj)
c        enddo
      enddo
      write(lfno,*)
     $'     N'//label(ivar1)//
     $     '     0----|----1----|----2----|----3----|----4----|----5'
      write(lfno,'(a,i5)')'    Score: ',nscore
      trval=nscore
      call tltrm(kptbl)
      calpol=pol0
      return
      end

      character function rad62a(m,j)
      implicit none
      integer*4 m,j
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
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     xa(7),dampx,dampy,dampz,aenox(np0),aenoy(np0),aenoz(np0)
      integer*4 kptbl(np),mturn(np0),j
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
        if(
     $       xa(1)**2+xa(2)**2 .lt. aenox(j) .and.
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
