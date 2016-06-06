      subroutine tffswake(ibegin,frbegin,lend,frend,over,beg)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use ffs_wake
      use tffitcode
      implicit none
      integer*8 ktaloc,iwbufxy,iwbufxyl,iwl,iwt,iwsl,iwst,i1,i2
      integer*4 ibegin,lend,itp1,i,j,ii,iutp,iwp,iutp1,np,npf,idx,idy
      real*8 frbegin,frend,
     $     rgetgl1,utwiss1(ntwissfun,-nfam:nfam),fr1,fr2,sigz,
     $     wbuf(nfam1:nfam),wzl(nfam1:nfam),wzt(nfam1:nfam),dx,dy
      logical*4 over(-nfam:nfam),beg
      sigz=rgetgl1('SIGZ')
      i1=ibegin
      fr1=frbegin
      itp1=itwissp(ibegin)
      if(.not. beg)then
        do j=nfam1,nfam
          if(uini(1,j) .le. 0.d0)then
            do i=1,mfitdetr
              utwiss(i,j,itp1)=twiss(ibegin,0,i)
            enddo
          else
            do i=1,mfitdetr
              utwiss(i,j,itp1)=uini(i,j)
            enddo
          endif
          do i=mfitdx,mfitddp
            utwiss(i,j,itp1)=twiss(ibegin,0,i)+uini(i,j)
          enddo
        enddo
      endif
      iwsl=0
      iwst=0
      np=nfam-nfam1+1
      iwbufxy=ktaloc(np**2)
      iwbufxyl=ktaloc(np**2)
      do iwp=1,nwakep+1
        if(iwp .eq. nwakep+1)then
          i2=lend
          fr2=frend
        else
          i2=iwakeelm(iwp)
          if(i2 .gt. lend)then
            cycle
          endif
          if(i2 .eq. lend)then
            fr2=frend
          else
            fr2=0.d0
          endif
        endif
        if(i2 .lt. i1)then
          cycle
        elseif(i2 .gt. i1 .or.
     $         i2 .eq. i1 .and. fr2 .ne. 0.d0)then
          do i=nfam1,nfam
            ii=min(1,abs(i))
            do j=1,ntwissfun
              twiss(i1,ii,j)=utwiss(j,i,itp1)
            enddo
            call qcell1(i1,fr1,i2,fr2,ii,
     1           hstab(i),vstab(i),tracex(i),tracey(i),
     $           i .ne. 0,over(i),.true.,0)
            call tffssetutwiss(i,nlat,
     $           i1,i2,fr2,beg,i1 .eq. ibegin,i2 .eq. lend)
          enddo
        elseif(i1 .eq. lend)then
          npf=np*ntwissfun
          if(idtype(latt(1,nlat-1)) .eq. icMARK)then
            iutp=itwissp(nlat-1)
            call tmov(utwiss(1,nfam1,itp1),
     $           utwiss(1,nfam1,iutp),np)
          endif
          if(i1 .lt. nlat)then
            iutp=itwissp(nlat)
            call tmov(utwiss(1,nfam1,itp1),
     $           utwiss(1,nfam1,iutp),np)
          endif
          return
        endif
        if(iwp .le. nwakep)then
          iutp=itwissp(i2)
          idx=kytbl(kwDX,idtype(latt(1,i2)))
          if(idx .ne. 0)then
            dx=rlist(latt(2,i2)+idx)
          else
            dx=0.d0
          endif
          idy=kytbl(kwDY,idtype(latt(1,i2)))
          if(idy .ne. 0)then
            dy=rlist(latt(2,i2)+idy)
          else
            dy=0.d0
          endif
          iwl=abs(kwaketbl(1,iwp))
          iwt=abs(kwaketbl(2,iwp))
          call tffswakekick(utwiss(1:ntwissfun,-nfam,iutp),utwiss1,
     $         iwl,iwt,rlist(iwl),rlist(iwt),
     $         sigz,2.d0*gammab(i2)*amass,nfam,nfam1,
     $         dx,dy,wbuf,rlist(iwbufxy),rlist(iwbufxyl),
     $         wzl,wzt,iwsl,iwst)
          do i=nfam1,nfam
            ii=min(1,abs(i))
            twiss(i2,ii,1:ntwissfun)=utwiss1(1:ntwissfun,i)
            call qcell1(i2,0.d0,i2+1,0.d0,ii,
     1           hstab(i),vstab(i),tracex(i),tracey(i),
     $           i .ne. 0,over(i),.true.,0)
            utwiss1(1:ntwissfun,i)=twiss(i2+1,ii,1:ntwissfun)
          enddo
          iutp1=itwissp(i2+1)
          call tffswakekick(utwiss1,
     $         utwiss(1:ntwissfun,-nfam,iutp1),
     $         iwl,iwt,rlist(iwl),rlist(iwt),
     $         sigz,2.d0*gammab(i2+1)*amass,nfam,nfam1,
     $         dx,dy,wbuf,rlist(iwbufxy),rlist(iwbufxyl),
     $         wzl,wzt,iwsl,iwst)
          i1=i2+1
          fr1=0.d0
          itp1=itwissp(i1)
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
     $     zi,dz,w,sigg,a,p,fact,dz0,xi,yi,cp0,dsig,ddz0,
     $     dx,dy,z1,z2,d
      logical*4 re
c      parameter (fact=1.d0/sqrt(pi2))
      parameter (fact=.398942280401433d0)
      np=nfam-nfam1+1
      call tmov(ut0(1,nfam1),ut1(1,nfam1),ntwissfun*np)
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
        do i=nfam1,nfam
          ut1(mfitddp,i)=ut0(mfitddp,i)-cp0*wbuf(i)
        enddo
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
        call tclr(wbufx,np)
        call tclr(wbufy,np)
        call tclr(wbufxl,np)
        call tclr(wbufyl,np)
        if(re)then
          do i=nfam1,nfam
            xi=ut0(mfitdx,i)-dx
            yi=ut0(mfitdy,i)-dy
            do j=nfam1,nfam
              wbufx(j)=wbufx(j)+xi*wbufxy(i,j)
              wbufy(j)=wbufy(j)+yi*wbufxy(i,j)
              wbufxl(j)=wbufxl(j)-xi*wbufxyl(i,j)
              wbufyl(j)=wbufyl(j)-yi*wbufxyl(i,j)
            enddo
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
      use iso_c_binding
      implicit none
      integer*8 kx,ktaloc,kal,kalj,ktfmalocp,kat,katj
      integer*4 irtc,isp0,isp1,lfno,n,m,lenw,l,itfdownlevel,isp2,
     $     i,lbegin,lend,j,k
      real*8 frbegin,frend
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
      call tffsbound(lbegin,frbegin,lend,frend)
      isp0=isp
      do i=lbegin,lend-1
        isp1=isp+2
        kal=0
        kat=0
        call elname(i,name)
        ilist(1,ifname)=lenw(name)
        call tfpadstr(name,ifname+1,ilist(1,ifname))
        levele=levele+1
        call tfleval(klist(ifwfunl-3),kx,.true.,irtc)
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
        call tfleval(klist(ifwfunt-3),kx,.true.,irtc)
        isp=isp2
        if(irtc .ne. 0)then
          if(ierrorprint .ne. 0)then
            call tfaddmessage(kerror,'WakeFunction Transverse',0,lfno)
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
