c$Header: /SAD/cvsroot/oldsad/src/trackd.f,v 1.30 2011/08/03 23:29:54 oide Exp $
      subroutine trackd(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     mturn,kzx,trval,phix,phiy,lfno)
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 nzp0,nxp,lfno,ncons,nscore
      parameter (nzp0=200,nxp=51)
      integer*4 latt(2,nlat),kptbl(np0,4),nzp,iv,
     $     iax,iax1,iax11,iax12,iax13,iax2,lp0,nzp1,
     $     nzstep,ipr,j,n,jzout,np1,k,np,kp,kz,kx,irw,nsc,
     $     iaxi,iaxi3,jj,ip,isw,kseed
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),
     1     pz(np0)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep,zstep,
     $     emx,emz,trval,rgetgl1,phix,phiy,cx,cy,sx,sy
      character rad62
      character*14 vname
      data vname /'ResultOfDAPERT'/
      integer*4 ntloss(nzp0,nxp),mturn(np0),nxm(nzp0),nxr(nzp0),
     $     kzx(2,np0),muls,ntlm(1),irtc
      integer*8 intlm,kb,mapalloc8
      integer*4 mapfree
      real*8 zi(nzp0)
      integer*4 i,itraaloc,itaaloc,itfsymbol,fork,wait,ichpid(100)
      logical ini,remain
c     begin initialize for preventing compiler warning
      ipr=0
c     end   initialize for preventing compiler warning
      ncons=max(1,nint(rgetgl1('DAPWIDTH')))
      trpt=.false.
      trval=0.d0
      nscore=0
      muls=(nturn/6200+1)*100
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
        zmin=1.d38
        zmax=-1.d38
        nzp=pz(1)
        do i=1,nzp
          zi(i)=g(i)
          zmin=min(zmin,zi(i))
          zmax=max(zmax,zi(i))
        enddo
      endif
      write(lfno,9101)xmin,xmax,ymin,ymax,zmin,zmax,muls
 9101 format(
     $     ' Range    Xmin:',f8.3,' Xmax:',f8.3,/
     1     '         (Ymin:',f8.3,' Ymax:',f8.3,')',/
     1     '          Zmin:',f8.3,' Zmax:',f8.3,/
     1     ' Display: ',i5,' turns/character')
      iv=itfsymbol(vname,len(vname),.false.)-4
      call tfreev(ilist(1,iv-1),ilist(2,iv-1))
      iax=itaaloc(0,2)
      ilist(1,iv-1)=ntflist
      ilist(2,iv-1)=iax
      iax1=itaaloc(0,3)
      call tfsetlist(ntflist,iax1,0.d0,iax,1)
      iax11=itaaloc(0,2)
      call tfsetlist(ntflist,iax11,0.d0,iax1,1)
      call tfsetlist(ntfreal,0,xmin,iax11,1)
      call tfsetlist(ntfreal,0,xmax,iax11,2)
      iax12=itaaloc(0,2)
      call tfsetlist(ntflist,iax12,0.d0,iax1,2)
      call tfsetlist(ntfreal,0,ymin,iax12,1)
      call tfsetlist(ntfreal,0,ymax,iax12,2)
      iax13=itaaloc(0,2)
      call tfsetlist(ntflist,iax13,0.d0,iax1,3)
      call tfsetlist(ntfreal,0,zmin,iax13,1)
      call tfsetlist(ntfreal,0,zmax,iax13,2)
      iax2=itaaloc(0,nzp)
      call tfsetlist(ntflist,iax2,0.d0,iax,2)
      lp0=latt(2,1)+kytbl(kwmax,idtype(latt(1,1)))+1
      np1=np0
      emx=sqrt(abs(rgetgl1('EMITX'))+abs(rgetgl1('EMITY')))
      if(rfsw)then
        emz=sqrt(abs(rgetgl1('EMITZ')))
      else
        emz=abs(rgetgl1('SIGE'))
      endif
      call tpara(latt)
      write(lfno,*)
     1'     NZ     0----!----1----!----2----!----3----!----4----!----5'
      if(nparallel .gt. 1 .and. nzp .ge. nparallel)then
        kseed=0
        irtc=1
        intlm=mapalloc8(ntlm(1), nzp0*nxp, 4, irtc)
        nzp1=1
        nzstep=nparallel
        ipr=1
        do while(nzp1 .lt. nparallel .and. ipr .gt. 0)
          kseed=kseed+2
          ipr=fork()
          if(ipr .gt. 0)then
            call tfaddseed(kseed,irtc)
            if(irtc .ne. 0)then
              return
            endif
            ichpid(nzp1)=ipr
            nzp1=nzp1+1
          endif
        enddo
      else
        nzp1=1
        nzstep=1
        intlm=0
      endif
      do i=nzp1,nzp,nzstep
        do j=1,nxp
          ntloss(i,j)=0
        enddo
        nxm(i)=nxp
        nxr(i)=nxp
      enddo
      do i=1,np0
        kzx(1,i)=0
        kzx(2,i)=0
      enddo
      n=1
      jzout=1
 1    np1=np0
      remain=.true.
      cx=cos(phix)
      sx=sin(phix)
      cy=cos(phiy)
      sy=sin(phiy)
      do 100 k=1,np0
        if(kzx(1,k) .eq. 0)then
          ip=kptbl(k,1)
          if(ip .le. np1)then
            if(remain)then
              do 110 i=nzp1,nzp,nzstep
                if(nxm(i) .ne. 0)then
                  if(nxm(i) .ge. nxr(i)-ncons+1)then
                    j=nxm(i)
                    if(ntloss(i,j) .eq. 0)then
                      kzx(1,k)=i
                      kzx(2,k)=j
                      mturn(k)=0
                      ntloss(i,j)=-1
                      nxm(i)=min(nxm(i),j-1)
                      x(ip)=xstep*(j-1)+xmin
                      px(ip)=-x(ip)*sx
                      x(ip)=x(ip)*cx
                      y(ip)=ystep*(j-1)+ymin
                      py(ip)=-y(ip)*sy
                      y(ip)=y(ip)*cy
                      z(ip)=0.d0
                      g(ip)=zi(i)
                      call tinip(1,
     1                     x(ip),px(ip),y(ip),py(ip),
     1                     z(ip),g(ip),dv(ip),pz(ip),
     1                     emx,emz,codin,.false.)
c       write(*,*)' trackd-Launch ',i,j,
c     $                     nxr(i),nxm(i),x(ip),px(ip),g(ip)
                      go to 100
                    endif
                  endif
                endif
 110          continue
            endif
            remain=.false.
            kzx(1,k)=0
            j=kptbl(np1,3)
            kptbl(j,1)=ip
            kptbl(ip,3)=j
            kptbl(k,1)=np1
            kptbl(np1,3)=k
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
 100  continue
      np=np1
 101  continue
c      write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-1 ',n,np,(kptbl(i,1),y(i),i=1,14)
      call tturn(np,latt,x,px,y,py,z,g,dv,pz,kptbl,n)
c      write(*,'(a,2i5,14(i5,1pg12.5))')
c     $     'trackd-tturn-2 ',n,np,(kptbl(i,1),y(i),i=1,14)
      n=n+1
      ini=.false.
      do 210 i=1,np0
        if(kzx(1,i) .le. 0)then
          go to 210
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
            do 220 j=min(nxp-ncons+1,kx),min(kx+ncons-1,nxp)
              if(ntloss(kz,j) .ne. nturn)then
                go to 210
              endif
 220        continue
            do 230 j=1,nxm(kz)
              ntloss(kz,j)=1000000
 230        continue
            nxm(kz)=0
          endif
        else
          kz=kzx(1,i)
          kx=kzx(2,i)
          ntloss(kz,kx)=mturn(i)
          kzx(1,i)=0
c          write(*,'(a,1x,9i6)')'trackd-Lost ',
c     1         kz,kx,nxr(kz),nxm(kz),mturn(i),i,np,np1,kp
          nxr(kz)=min(nxr(kz),kx-1)
          ini=.true.
        endif
 210  continue
      if(ini)then
        do 410 i=jzout,nzp,nzstep
          if(nxm(i) .ne. 0)then
            go to 1
          endif
          do 420 k=1,np0
            if(kzx(1,k) .eq. i)then
              go to 1
            endif
 420      continue
          jzout=i+nzstep
 410    continue
        go to 1
      else
        go to 101
      endif
 3000 continue
      if(nparallel .gt. 1 .and. nzp .ge. nparallel)then
        if(ipr .eq. 0)then
          do k=1,nxp
            kb=intlm+(k-1)*nzp0
            do i=nzp1,nzp,nzstep
              ntlm(kb+i)=ntloss(i,k)
            enddo
          enddo
          stop
        else
          do j=1,nparallel-1
 3010       irw=wait(isw)
            do k=1,nparallel-1
              if(irw .eq. ichpid(k))then
                n=k
                ichpid(k)=0
                go to 3020
              endif
            enddo
            write(*,*)'???trackd-Strange child: ',irw,isw
            go to 3010
 3020       do k=1,nxp
              kb=intlm+(k-1)*nzp0
              do i=n,nzp,nzstep
                ntloss(i,k)=ntlm(kb+i)
              enddo
            enddo
          enddo
        endif
      endif
      if(intlm .ne. 0)then
        if(mapfree(ntlm(intlm+1)) .ne. 0)then
          write(*,*)'???trackd-error in munmap.'
        endif
      endif
      do i=1,nzp
        nsc=0
        do 430 k=1,nxp
          if(ntloss(i,k) .lt. nturn)then
            nsc=k-1
            go to 431
          endif
 430    continue
        nsc=nxp
 431    nscore=nscore+nsc
        write(lfno,9001)zi(i),nsc,
     $       (rad62(ntloss(i,jj)/muls),jj=1,nxp)
 9001   format(1x,f8.2,i3,1x,51a1,i3)
c     write(*,*)'trackd-end '
        iaxi=itaaloc(0,3)
        call tfsetlist(ntflist,iaxi,0.d0,iax2,i)
        call tfsetlist(ntfreal,0,zi(i),iaxi,1)
        call tfsetlist(ntfreal,0,dble(nsc),iaxi,2)
        iaxi3=itraaloc(0,nxp)
        call tfsetlist(ntflist,iaxi3,0.d0,iaxi,3)
        do jj=1,nxp
          rlist(ilist(2,iaxi3+1)+jj)=ntloss(i,jj)
        enddo
      enddo
      write(lfno,*)
     1'     NZ     0----!----1----!----2----!----3----!----4----!----5'
      write(lfno,'(a,i5)')'    Score: ',nscore
      trval=nscore
      call tltrm(latt,kptbl)
      return
      end

      subroutine tinip(np,x,px,y,py,z,g,dv,pz,emx,emz,codin,cmplot)
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),g(np),pz(np)
      real*8 xa(8),rgetgl1,tgauss,dxi,dyi,dzi,emx,emz,codin(6)
      logical*4 cmplot
      dxi=rgetgl1('TDXI')
      dyi=rgetgl1('TDYI')
      dzi=rgetgl1('TDZI')
      do 10 i=1,np
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
        dv(i)=xa(7)
10    continue
c      write(*,*)'tinip ',np,(xa(i),i=1,6)
      return
      end
