      subroutine tffswake(fbound,beg)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use ffs_wake
      use tffitcode
      use iso_c_binding
      implicit none
      type (ffs_bound) ,intent(in):: fbound
      type (ffs_bound) ibound,ibound1
      integer*8 iwbufxy,iwbufxyl,iwl,iwt,iwsl,iwst
      integer*4 itp1,i,j,ii,iutp,iwp,iutp1,np,npf,idx,idy,
     $     nlw,ntw
      real*8 rgetgl1,utwiss1(ntwissfun,-nfam:nfam),sigz,
     $     wbuf(nfam1:nfam),wzl(nfam1:nfam),wzt(nfam1:nfam),dx,dy
      real*8 , pointer::wakel(:,:),waket(:,:),wbufxy(:,:),wbufxyl(:,:)
      logical*4 ,intent(in):: beg
      sigz=rgetgl1('SIGZ')
      ibound1%fb=0.d0
      ibound1%fe=0.d0
      ibound=fbound
      itp1=itwissp(fbound%lb)
      if(.not. beg)then
        do j=nfam1,nfam
          if(uini(mfitbx,j) .le. 0.d0)then
            utwiss(:,j,itp1)=twiss(fbound%lb,0,:)
          else
            utwiss(:,j,itp1)=uini(:,j)
          endif
          utwiss(mfitdx:mfitddp,j,itp1)=
     $         twiss(fbound%lb,0,mfitdx:mfitddp)+uini(mfitdx:mfitddp,j)
        enddo
      endif
      iwsl=0
      iwst=0
      np=nfam-nfam1+1
      iwbufxy=ktaloc(np**2)
      iwbufxyl=ktaloc(np**2)
      call c_f_pointer(c_loc(rlist(iwbufxy)),wbufxy,[np,np])
      call c_f_pointer(c_loc(rlist(iwbufxyl)),wbufxyl,[np,np])
      do iwp=1,nwakep+1
        if(iwp .eq. nwakep+1)then
          ibound%le=fbound%le
          ibound%fe=fbound%fe
        else
          ibound%le=iwakeelm(iwp)
          if(ibound%le .gt. fbound%le)then
            cycle
          endif
          if(ibound%le .eq. fbound%le)then
            ibound%fe=fbound%fe
          else
            ibound%fe=0.d0
          endif
        endif
        if(ibound%le .lt. ibound%lb)then
          cycle
        elseif(ibound%le .gt. ibound%lb .or.
     $         ibound%le .eq. ibound%lb .and. ibound%fe .ne. 0.d0)then
          do i=nfam1,nfam
            ii=min(1,abs(i))
            twiss(ibound%lb,ii,:)=utwiss(:,i,itp1)
            call qcell1(ibound,ii,optstat(i),i .ne. 0,.true.,0)
            call tffssetutwiss(i,nlat,ibound,beg,
     $           ibound%lb .eq. fbound%lb,ibound%le .eq. fbound%le)
          enddo
        elseif(ibound%lb .eq. fbound%le)then
          npf=np*ntwissfun
          if(idtypec(nlat-1) .eq. icMARK)then
            iutp=itwissp(nlat-1)
            utwiss(:,nfam1:nfam,iutp)=utwiss(:,nfam1:nfam,itp1)
c            call tmov(utwiss(1,nfam1,itp1),
c     $           utwiss(1,nfam1,iutp),npf)
          endif
          if(ibound%lb .lt. nlat)then
            iutp=itwissp(nlat)
            utwiss(:,nfam1:nfam,iutp)=utwiss(:,nfam1:nfam,itp1)
c            call tmov(utwiss(1,nfam1,itp1),
c     $           utwiss(1,nfam1,iutp),npf)
          endif
          return
        endif
        if(iwp .le. nwakep)then
          iutp=itwissp(ibound%le)
          idx=kytbl(kwDX,idtypec(ibound%le))
          if(idx .ne. 0)then
            dx=rlist(latt(ibound%le)+idx)
          else
            dx=0.d0
          endif
          idy=kytbl(kwDY,idtypec(ibound%le))
          if(idy .ne. 0)then
            dy=rlist(latt(ibound%le)+idy)
          else
            dy=0.d0
          endif
          iwl=abs(kwaketbl(1,iwp))
          iwt=abs(kwaketbl(2,iwp))
          nlw=(ilist(1,iwl-1)-1)/2
          ntw=(ilist(1,iwt-1)-1)/2
          call c_f_pointer(c_loc(rlist(iwl)),wakel,[2,nlw])
          call c_f_pointer(c_loc(rlist(iwt)),waket,[2,ntw])
          call tffswakekick(utwiss(1:ntwissfun,-nfam,iutp),utwiss1,
     $         iwl,iwt,wakel,waket,
     $         sigz,2.d0*gammab(ibound%le)*amass,nfam,nfam1,
     $         dx,dy,wbuf,wbufxy,wbufxyl,
     $         wzl,wzt,iwsl,iwst)
          ibound1%lb=ibound%le
          ibound1%le=ibound%le+1
          do i=nfam1,nfam
            ii=min(1,abs(i))
            twiss(ibound1%lb,ii,:)=utwiss1(:,i)
            call qcell1(ibound1,ii,optstat(i),i .ne. 0,.true.,0)
            utwiss1(:,i)=twiss(ibound1%le,ii,:)
          enddo
          iutp1=itwissp(ibound1%le)
          call tffswakekick(utwiss1,
     $         utwiss(1:ntwissfun,-nfam,iutp1),
     $         iwl,iwt,wakel,waket,
     $         sigz,2.d0*gammab(ibound%le+1)*amass,nfam,nfam1,
     $         dx,dy,wbuf,wbufxy,wbufxyl,
     $         wzl,wzt,iwsl,iwst)
          ibound%lb=ibound%le+1
          ibound%fb=0.d0
          itp1=itwissp(ibound%lb)
        endif
      enddo
      call tfree(iwbufxy)
      call tfree(iwbufxyl)
      return
      end

      subroutine tffswakekick(ut0,ut1,
     $     iwl,iwt,wakel,waket,sigz,p,
     $     nfam,nfam1,dx,dy,
     $     wbuf,wbufxy,wbufxyl,wzl,wzt,iwsl,iwst)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 iwl,iwt,iwsl,iwst
      integer*4 nfam,nfam1,i,j,lw,i1,i2,ic,np
      real*8 ut0(ntwissfun,-nfam:nfam),
     $     ut1(ntwissfun,-nfam:nfam),
     $     wbuf(nfam1:nfam),
     $     wbufx(nfam1:nfam),wbufy(nfam1:nfam),
     $     wbufxl(nfam1:nfam),wbufyl(nfam1:nfam),
     $     wbufxy(nfam1:nfam,nfam1:nfam),
     $     wbufxyl(nfam1:nfam,nfam1:nfam),
     $     wzl(nfam1:nfam),wzt(nfam1:nfam),
     $     wakel(2,*),waket(2,*),sigz,
     $     zi,dz,w,sigg,a,p,dz0,xi,yi,cp0,dsig,ddz0,
     $     dx,dy,z1,z2,d
      logical*4 re
      real*8 , parameter :: fact=.398942280401433d0
      np=nfam-nfam1+1
      ut1(:,1:np)=ut0(:,1:np)
c      call tmov(ut0(1,nfam1),ut1(1,nfam1),ntwissfun*np)
      cp0=abs(charge*e*pbunch/np/p)
      sigg=max(1.d-5,sigz/np/anbunch)
      dsig=.1d0*sigg
      if(iwl .ne. 0)then
        if(iwl .eq. iwsl)then
          ddz0=ut0(mfitdz,nfam1)-wzl(nfam1)
          re=.false.
          do i=nfam1+1,nfam
            if(abs(ut0(mfitdz,i)-wzl(i)-ddz0) .gt. dsig)then
              re=.true.
              exit
            endif
          enddo
        else
          re=.true.
        endif
        if(re)then
          iwsl=iwl
          lw=(ilist(1,iwl-1)-1)/2
          z1=-1d100
          i1=1
          z2=1d100
          i2=lw
          call tclr(wbuf,np)
          do i=nfam1,nfam
            zi=ut0(mfitdz,i)
            wzl(i)=zi
            do j=nfam1,nfam
              dz=zi-ut0(mfitdz,j)
              if(dz .ge. -sigg)then
                w=min(1.d0,(dz+sigg)/sigg*.5d0)
                if(dz .le. wakel(1,1))then
                  wbuf(j)=wbuf(j)+w*wakel(2,1)
                elseif(dz .lt. wakel(1,lw))then
                  if(dz .ge. z2)then
                    i1=i2
                    i2=lw
                  elseif(dz .lt. z1)then
                    i2=i1
                    i1=1
                  endif
                  do while(i2 .gt. i1+1)
                    ic=i1+(i2-i1)/2
                    if(wakel(1,ic) .gt. dz)then
                      i2=ic
                    else
                      i1=ic
                    endif
                  enddo
                  z1=wakel(1,i1)
                  z2=wakel(1,i2)
                  a=(dz-z1)/(z2-z1)
                  wbuf(j)=wbuf(j)
     $                 +w*((1.d0-a)*wakel(2,i1)+a*wakel(2,i2))
                endif
              endif
            enddo
          enddo
        endif
        ut1(mfitddp,:)=ut0(mfitddp,:)-cp0*wbuf
      endif
      if(iwt .ne. 0)then
        re=.false.
        if(iwt .eq. iwst)then
          ddz0=ut0(mfitdz,nfam1)-wzt(nfam1)
          re=.true.
          do i=nfam1+1,nfam
            if(abs(ut0(mfitdz,i)-wzt(i)-ddz0) .gt. dsig)then
              re=.false.
              exit
            endif
          enddo
        endif
        wbufx=0.d0
        wbufy=0.d0
        wbufxl=0.d0
        wbufyl=0.d0
        if(re)then
          do i=nfam1,nfam
            xi=ut0(mfitdx,i)-dx
            yi=ut0(mfitdy,i)-dy
            wbufx=wbufx+xi*wbufxy(i,:)
            wbufy=wbufy+yi*wbufxy(i,:)
            wbufxl=wbufxl-xi*wbufxyl(i,:)
            wbufyl=wbufyl-yi*wbufxyl(i,:)
          enddo
        else
          iwst=iwt
          lw=(ilist(1,iwt-1)-1)/2
          i1=1
          z1=-1d100
          i2=lw
          z2=1d100
          dz0=sigg*fact
          do i=nfam1,nfam
            zi=ut0(mfitdz,i)
            xi=ut0(mfitdx,i)-dx
            yi=ut0(mfitdy,i)-dy
            wzt(i)=zi
            do j=nfam1,nfam
              dz=zi-ut0(mfitdz,j)+dz0
              if(dz .gt. waket(1,1) .and. dz .lt. waket(1,lw))then
                if(dz .ge. z2)then
                  i1=i2
                  i2=lw
                elseif(dz .lt. z1)then
                  i2=i1
                  i1=1
                endif
                do while(i2 .gt. i1+1)
                  ic=i1+(i2-i1)/2
                  if(waket(1,ic) .gt. dz)then
                    i2=ic
                  else
                    i1=ic
                  endif
                enddo
                z1=waket(1,i1)
                z2=waket(1,i2)
                a=(dz-z1)/(z2-z1)
                d=(waket(2,i2)-waket(2,i1))/(z2-z1)
                w=waket(2,i1)+d*(dz-z1)
                wbufx(j)=wbufx(j)+xi*w
                wbufy(j)=wbufy(j)+yi*w
                wbufxl(j)=wbufxl(j)-xi*d
                wbufyl(j)=wbufyl(j)-yi*d
                wbufxy(i,j)=w
                wbufxyl(i,j)=d
              else
                wbufxy(i,j)=0.d0
                wbufxyl(i,j)=0.d0
              endif
            enddo
          enddo
        endif
        do i=nfam1,nfam
          ut1(mfitdpx,i)=ut0(mfitdpx,i)+cp0*wbufx(i)
          ut1(mfitdpy,i)=ut0(mfitdpy,i)+cp0*wbufy(i)
          ut1(mfitddp,i)=ut1(mfitddp,i)+
     $         cp0*(wbufxl(i)*(ut0(mfitdx,i)-dx)
     $         +wbufyl(i)*(ut0(mfitdy,i)-dy))
        enddo
      endif
      return
      end

      subroutine tffssetupwake(lfno,irtc)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use ffs_wake
      use eeval
      use iso_c_binding
      implicit none
      type (sad_dlist) ,pointer :: kwll,kwtl
      type (sad_descriptor) kx
      type (ffs_bound) fbound
      integer*8 kal,kalj,ktfmalocp,kat,katj
      integer*4 irtc,isp0,isp1,lfno,n,m,lenw,l,itfdownlevel,isp2,
     $     i,j,k
      character*(MAXPNAME+10) name
      integer*8 ifname,ifwfunl,ifwfunt
      save ifname,ifwfunl,ifwfunt
      data ifwfunl,ifwfunt/0,0/
      if(ifwfunl .eq. 0)then
        ifwfunl=ktadaloc(0,2)
        ifwfunt=ktadaloc(0,2)
        dlist(ifwfunl)=dtfcopy1(kxsymbolz('WakeFunction',12))
        klist(ifwfunt)=ktfcopy1(klist(ifwfunl))
        dlist(ifwfunl+1)=kxsymbolz('Longitudinal',12)
        dlist(ifwfunt+1)=kxsymbolz('Transverse',10)
        ifname=ktsalocb(1,name,len(name))
        klist(ifwfunl+2)=ktfstring+ifname
        klist(ifwfunt+2)=ktfstring+ifname
      endif
      call tffsbound(fbound)
      isp0=isp
      do i=fbound%lb,fbound%le-1
        isp1=isp+2
        kal=0
        kat=0
        call elname(i,name)
        ilist(1,ifname)=lenw(name)
        call tfpadstr(name,ifname+1,ilist(1,ifname))
        levele=levele+1
        call descr_sad(dlist(ifwfunl),kwll)
        kx=tfleval(kwll,.true.,irtc)
        if(irtc .ne. 0)then
          if(ierrorprint .ne. 0)then
            call tfaddmessage('WakeFunction Longitudinal',0,lfno)
          endif
          go to 9000
        endif
        if(ktflistq(kx))then
          kal=ktfmalocp(kx,n,m,.false.,.true.,.false.,.false.,irtc)
          if(irtc .eq. 0 .and. m .eq. 2)then
            isp=isp1
            LOOP_J_1: do j=isp-3,isp0+1,-2
              kalj=abs(ktastk(j))
              if(kalj .ne. 0)then
                if(ilist(1,kal-1) .eq. ilist(1,kalj-1))then
                  do k=0,ilist(1,kal-1)-2
                    if(rlist(kal+k) .ne. rlist(kalj+k))then
                      cycle LOOP_J_1
                    endif
                  enddo
                  call tfree(kal)
                  kal=-kalj
                  go to 20
                endif
              endif
            enddo LOOP_J_1
          else
            kal=0
          endif
        endif
 20     isp2=isp
        call descr_sad(dlist(ifwfunt),kwtl)
        kx=tfleval(kwtl,.true.,irtc)
        isp=isp2
        if(irtc .ne. 0)then
          if(ierrorprint .ne. 0)then
            call tfaddmessage('WakeFunction Transverse',0,lfno)
          endif
          go to 9000
        endif
        if(ktflistq(kx))then
          kat=ktfmalocp(kx,n,m,.false.,.true.,.false.,.false.,irtc)
          if(irtc .eq. 0 .and. m .eq. 2)then
            isp=isp1
            LOOP_J_2: do j=isp-2,isp0+2,-2
              katj=abs(ktastk(j))
              if(katj .ne. 0)then
                if(ilist(1,kat-1) .eq. ilist(1,katj-1))then
                  do k=0,ilist(1,kat-1)-2
                    if(rlist(kat+k) .ne. rlist(katj+k))then
                      cycle LOOP_J_2
                    endif
                  enddo
                  call tfree(kat)
                  kat=-katj
                  go to 120
                endif
              endif
            enddo LOOP_J_2
          else
            kat=0
          endif
        endif
 120    l=itfdownlevel()
        ktastk(isp1-1)=kal
        ktastk(isp1  )=kat
        ivstk2(1,isp1)=i
      enddo
      nwakep=(isp-isp0)/2
      if(nwakep .gt. 0)then
        kwakep=ktaloc(nwakep*2)
        kwakeelm=ktaloc(nwakep)
        call tmov(ktastk(isp0+1),klist(kwakep),nwakep*2)
        do i=1,nwakep
          ilist(i,kwakeelm)=ivstk2(1,isp0+i*2)
        enddo
      endif
      call c_f_pointer(c_loc(ilist(1,kwakeelm)),iwakeelm,[nwakep])
      call c_f_pointer(c_loc(klist(kwakep)),kwaketbl,[2,nwakep])
      irtc=0
      return
 9000 l=itfdownlevel()
      do i=isp0+1,isp
        if(ktastk(i) .gt. 0)then
          call tfree(ktastk(i))
        endif
      enddo
      isp=isp0
      return
      end

      subroutine tffsclearwake
      use tfstk
      use ffs_wake
      implicit none
      integer*8 i
      if(nwakep .gt. 0)then
        do i=kwakep,kwakep+2*nwakep-1
          if(klist(i) .gt. 0)then
            call tfree(klist(i))
          endif
        enddo
        call tfree(kwakep)
        call tfree(kwakeelm)
      endif
      return
      end
