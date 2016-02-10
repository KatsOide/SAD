c     CAUTION: kptbl(#,3) MUST be `0' before trackd() called
      subroutine trackd(latt,kptbl,x,px,y,py,z,g,dv,pz,
     $     mturn,kzx,trval,phix,phiy,damp,lfno)
      use tfstk
      use tfshare
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 nzp0,nxp
      parameter (nzp0=200,nxp=51)
      integer*4 latt(2,nlat),kptbl(np0,6),lfno,
     $     mturn(np0),kzx(2,np0)
      integer*8 intlm
      integer*4 irtc
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 trval,phix,phiy,damp
      irtc=1
      intlm=ktfallocshared((nzp0*nxp+1)/2)
      call trackd0(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     mturn,kzx,trval,phix,phiy,damp,klist(intlm),lfno)
      call tfreeshared(intlm)
c      if(mapfree(ntloss(intlm+1)) .ne. 0)then
c        write(*,*)'???trackd-error in munmap.'
c      endif
      return
      end

      subroutine trackd0(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     mturn,kzx,trval,phix,phiy,damp,ntloss,lfno)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'

      integer*4 nzp0,nxp,maxturn,maxpara,nw,lfno,ncons,nscore
      parameter (nzp0=200,nxp=51,maxturn=2**30,maxpara=256,nw=16)
      integer, parameter :: nkptbl = 6
      integer*8 kv,kax,kax11,kax12,kax13,kax2,
     $     kaxi,kaxi3,kax1
      integer*4 latt(2,nlat),kptbl(np0,nkptbl),nzp,lp0,nzp1,
     $     ipr,j,n,jzout,np1,k,np,kp,kz,kx,irw,nsc,iw,
     $     jj,ip,isw,kseed,npmax,npara,nxm(nzp0)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     1     pz(np0)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep,zstep,
     $     emx,emz,rgetgl1,cx,cy,sx,sy,dampx,dampy,dampz,t0
      real*8 trval,phix,phiy,damp
      real*8 zi(nzp0)
      character*14 vname
      data vname /'ResultOfDAPERT'/
      integer*4 ntloss(nzp0,nxp),mturn(np0),kzx(2,np0),muls,irtc
      integer*4 i,fork_worker,wait,ichpid(maxpara)
      logical ini,remain
      character rad62a
c     begin initialize for preventing compiler warning
      ipr=0
c     end   initialize for preventing compiler warning
      if(damp .eq. 0.d0)then
        dampx=1.d0
        dampy=1.d0
        dampz=1.d0
      else
        t0=pi2/omega0
        dampx=exp(-t0/taurdx)
        dampy=exp(-t0/taurdy)
        dampz=exp(-t0/taurdz)
      endif
      ncons=max(1,nint(rgetgl1('DAPWIDTH')))
      trpt=.false.
      call tsetdvfs
      trval=0.d0
      nscore=0
      if(nturn .gt. 600)then
        muls=(nturn/6200+1)*100
      else
        muls=10
      endif
      xmin=x(1)
      xmax=x(2)
      ymin=y(1)
      ymax=y(2)
      zmin=min(g(1),g(2))
      zmax=max(g(1),g(2))
      if(zmax .eq. zmin)then
        zmin=0.d0
      endif
      xstep=(xmax-xmin)/(nxp-1)
      ystep=(ymax-ymin)/(nxp-1)
      if(pz(1) .eq. 0.d0)then
        zstep=abs(g(3))
        if(zstep .eq. 0.d0 .or. zmin .eq. zmax)then
          nzp=1
        else
          nzp=min(nzp0,int(max(2.d0,(zmax-zmin)/zstep+1.1d0)))
          zstep=(zmax-zmin)/(nzp-1)
        endif
        do i=1,nzp
          zi(i)=zmin+zstep*(i-1)
        enddo
      else
        nzp=int(pz(1))
        zi(1:nzp)=g(1:nzp)
        zmin=minval(zi(1:nzp))
        zmax=maxval(zi(1:nzp))
      endif
      write(lfno,'(
     $     '' Range    Xmin:'',f8.3,'' Xmax:'',f8.3,/
     $     ''         (Ymin:'',f8.3,'' Ymax:'',f8.3,'')'',/
     $     ''          Zmin:'',f8.3,'' Zmax:'',f8.3,/
     $     '' Display: '',i5,'' turns/character'')')
     $     xmin,xmax,ymin,ymax,zmin,zmax,muls
      write(lfno,*)
     $'     NZ     0----|----1----|----2----|----3----|----4----|----5'
      kv=ktfsymbolz(vname,len(vname))-4
      call tflocal(klist(kv))
      kax=ktadaloc(0,2)
      klist(kv)=ktflist+kax
      kax1=ktadaloc(0,3)
      klist(kax+1)=ktflist+kax1
      kax11=ktavaloc(0,2)
      klist(kax1+1)=ktflist+kax11
      rlist(kax11+1)=xmin
      rlist(kax11+2)=xmax
      kax12=ktavaloc(0,2)
      klist(kax1+2)=ktflist+kax12
      rlist(kax12+1)=ymin
      rlist(kax12+2)=ymax
      kax13=ktavaloc(0,2)
      klist(kax1+3)=ktflist+kax13
      rlist(kax13+1)=zmin
      rlist(kax13+2)=zmax
      kax2=ktadalocnull(0,nzp)
      klist(kax+2)=ktflist+kax2
      lp0=latt(2,1)+kytbl(kwmax,idtype(latt(1,1)))+1
      emx=sqrt(abs(rgetgl1('EMITX'))+abs(rgetgl1('EMITY')))
      if(rfsw)then
        emz=sqrt(abs(rgetgl1('EMITZ')))
      else
        emz=abs(rgetgl1('SIGE'))
      endif
      call tpara(latt)
      ntloss(1:nzp,1:nxp)=maxturn
      nxm(1:nzp)=nxp+1
      npara=min(nparallel,maxpara)
      ipr=1
      if(npara .gt. 1)then
        kseed=0
        nzp1=1
        do while(nzp1 .lt. npara .and. ipr .gt. 0)
          kseed=kseed+2
          ipr=fork_worker()
          if(ipr .gt. 0)then
            ichpid(nzp1)=ipr
            nzp1=nzp1+1
            call tfaddseed(kseed,irtc)
            if(irtc .ne. 0)then
              go to 3000
            endif
          endif
        enddo
        npmax=max(1,min(nzp*ncons/npara,np0))
      else
        nzp1=1
        npmax=np0
      endif
      kzx(1,1:npmax)=0
      kzx(2,1:npmax)=0
      n=1
      jzout=1
 1    np1=npmax
      iw=nw
      remain=.true.
      cx=cos(phix)
      sx=sin(phix)
      cy=cos(phiy)
      sy=sin(phiy)
      LOOP_K: do k=1,npmax
        if(kzx(1,k) .eq. 0)then
          ip=kptbl(k,1)
          if(ip .le. np1)then
            if(remain)then
              do i=1,nzp
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
                    x(ip)=xstep*(j-1)+xmin
                    px(ip)=-x(ip)*sx
                    x(ip)=x(ip)*cx
                    y(ip)=ystep*(j-1)+ymin
                    py(ip)=-y(ip)*sy
                    y(ip)=y(ip)*cy
                    z(ip)=0.d0
                    g(ip)=zi(i)
c     Reinit kptbl(ip,4) to reuse particle array slot `ip'
                    kptbl(ip,4)=0
                    call tinip(1,
     $                   x(ip),px(ip),y(ip),py(ip),
     $                   z(ip),g(ip),dv(ip),
     $                   emx,emz,codin,dvfs,.false.)
c                    write(*,'(a,4i5,1p3g15.7)')
c     $                   ' trackd-Launch ',nzp1,i,j,
c     $                     nxm(i),x(ip),px(ip),g(ip)
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
            pz(ip)=pz(np1)
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
      call tturn(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n)
c      write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-2 ',n,np,(kptbl(i,1),y(i),i=1,14)
      if(damp .ne. 0.d0)then
        call tpdamp(np,x,px,y,py,z,g,dv,dampx,dampy,dampz)
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
          ntloss(kz,kx)=mturn(i)
          kzx(1,i)=0
c          write(*,'(a,1x,9i6)')'trackd-Lost ',
c     $         kz,kx,nxr(kz),nxm(kz),mturn(i),i,np,np1,kp
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
c        write(*,*)'trackd-stop ',nzp1
c        stop
        call tfresetsharedmap()
        call exit_without_hooks(0)
      endif
      do j=1,nzp1-1
 3010   irw=wait(isw)
        do k=1,nzp1-1
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
      do i=1,nzp
        nsc=0
        do k=1,nxp
          if(ntloss(i,k) .lt. nturn)then
            nsc=k-1
            go to 431
          endif
        enddo
        nsc=nxp
 431    nscore=nscore+nsc
        write(lfno,'(1x,f8.2,i3,1x,51a1)')zi(i),nsc,
     $       (rad62a(ntloss(i,jj)/muls,jj),jj=1,nxp)
        kaxi=ktadaloc(0,3)
        klist(kax2+i)=ktflist+kaxi
        rlist(kaxi+1)=zi(i)
        rlist(kaxi+2)=dble(nsc)
        kaxi3=ktraaloc(0,nxp)
        klist(kaxi+3)=ktflist+kaxi3
        rlist(kaxi3+1:kaxi3+nxp)=ntloss(i,1:nxp)
c        do jj=1,nxp
c          rlist(kaxi3+jj)=ntloss(i,jj)
c        enddo
      enddo
      write(lfno,*)
     $'     NZ     0----|----1----|----2----|----3----|----4----|----5'
      write(lfno,'(a,i5)')'    Score: ',nscore
      trval=nscore
      call tltrm(latt,kptbl)
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
        call tadd(xa,codin,xa,6)
        call tconv(xa,xa,-1)
        x(i) =xa(1)+dxi
        px(i)=xa(2)
        y(i) =xa(3)+dyi
        py(i)=xa(4)
        z(i) =xa(5)+dzi
        g(i) =xa(6)
        dv(i)=xa(7)+dvfs
      enddo
c      write(*,*)'tinip ',np,(xa(i),i=1,6)
      return
      end

      subroutine tpdamp(np,x,px,y,py,z,g,dv,dampx,dampy,dampz)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),
     $     xa(7),dampx,dampy,dampz
      call tconvm(np,px,py,g,dv,1)
      do i=1,np
        xa(1)=x(i)-codin(1)
        xa(2)=px(i)-codin(2)
        xa(3)=y(i)-codin(3)
        xa(4)=py(i)-codin(4)
        xa(5)=z(i)-codin(5)
        xa(6)=g(i)-codin(6)
        call tmap(xa,xa,1)
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
      enddo
      call tconvm(np,px,py,g,dv,-1)
      return
      end
