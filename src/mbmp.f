      subroutine mbmp(latt,twiss,mult,master,isb,nbump,istr,nstr,nbcor,
     $     yplane,iret,lfno)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      implicit real*8 (a-h,o-z)
      parameter (nconj=3, keymax=50)
      logical tmatch,yplane,mhogal,pat,enome
      character*8 wkey,name
      character ermes*30
      integer*8 latt(*)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),mult(*)
      integer*4 master(nlat),id
      integer*8 itemp
      dimension isb(nbcor+4,nbump),istr(nstra,4)
      include 'inc/common.inc'
      common/mbumpf/wkey(keymax),ifte(nconj,keymax),nw
c
      save
c
      nbp=nbump
      iret=0
      m=0
      do 11 j=1,nw
        l=ifte(3,j)-1
        do 10 i=1,nlat-1
          call elnameK(i,name)
          enome=name.eq.wkey(j)
          id=idcomp(elatt,i)
          if( tmatch(pname(id),wkey(j)) .or. enome)then
            pat=idtype(id).eq.icmark .or. master(i).gt.0
c           print *,pnamec(i),name,enome
            if(pat .or. enome) then
              if(mhogal(ifte(1,j),ifte(2,j),i)) then
                l=l+1
                if(l.eq.ifte(3,j)) then
                  m=m+1
                  isb(1,m)=i
                  l=0
                  if(m.eq.nbp) goto 12
                endif
              endif
            endif
          endif
   10   continue
   11 continue
   12 if(nw.gt.1) then
        itemp=ktaloc((nbp+1)/2)
        do 13 i=1,nbp
          ilist(mod(i-1,2)+1,itemp+(i-1)/2)=isb(1,i)
   13   continue
c       ---- sort in ascending order and eliminate identical elements.
c            nbp might be changed. ----
        call msort(ilist(1,itemp),nbp)
        do 14 i=1,nbp
          isb(1,i)=ilist(mod(i-1,2)+1,itemp+(i-1)/2)
   14   continue
        call tfree(itemp)
      endif
      nbump=nbp
      if(yplane) then
        write(lfno,'(2(A,I4))')
     1       ' No.of y-bump=',nbump,' / Correctors per bump=',nbcor
      else
        write(lfno,'(2(A,I4))')
     1       ' No.of x-bump=',nbump,' / Correctors per bump=',nbcor
      endif
c     ---- find corrector around j-th target
      do 20 j=1,nbp
        call mfnst(latt,twiss,isb(1,j),nbcor,istr,nstr,yplane,iret,
     1             lfno)
        if(iret.ne.0) then
          write(ermes,'(A,I3)') '!! Error in MBMP: iret=',iret
          call permes(' ',ermes,' ',lfno)
          return
        endif
c       print *,isb(1,j)
c       call mbufw(pnamec(isb(1,j)),.false.,lfno)
c       write(ermes,'(i3)') isb(2,j)
c       call mbufw(ermes,.false.,lfno)
c       do i=3,isb(2,j)+2
c         call elname(latt,istr(istr(isb(i,j),2),1),mult,name)
c         call mbufw(name,.false.,lfno)
c         write(ermes,'(i3,a,i4,a)') isb(i,j),'(',
c    $         istr(istr(isb(i,j),2),1),')'
c         call mbufw(ermes,.false.,lfno)
c       enddo
c       call mbufw(' ',.true.,lfno)
   20 continue
      return
      end

      subroutine mbmpf(word,wordp,latt,mult,master,nb)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      implicit real*8 (a-h,o-z)
      parameter (nconj=3, keymax=50)
      logical tmatch,exist,abbrev,fcon,mhogal,pat,enome
      character*(*) word,wordp
      character*8 wkey,conj(nconj),name
      integer*4 master(nlat),mult(*),id
      integer*8 latt(*)
      common/mbumpf/wkey(keymax),ifte(nconj,keymax),nw
      data conj/'F_ROM   ','T_O     ','E_VERY  '/
c
      save
c      
      do 50 i=1,keymax
   50   wkey(i)=' '
      do 51 j=1,keymax
        ifte(1,j)=1
        ifte(2,j)=nlat
   51   ifte(3,j)=1
      nw=0
      nwc=0
      fcon=.true.
      goto 53
   52 call getwdl2(word,wordp)
   53 do 61 i=1,nconj
        if(abbrev(word,conj(i),'_')) then
          if(nw.eq.0) return
          if(i.ne.3) then
            call getwdl2(word,wordp)
            m=ielm(wordp,exist)
            if(exist) ifte(i,nw)=m
          else
            m=int(getva(exist))
            if(exist) then
              ifte(i,nw)=m
              call getwdl2(word,wordp)
            endif
          endif
          if(fcon) then
            nwc=nwc+1
            fcon=.false.
          endif
          if(nwc.ne.nw) then
            do 60 j=nwc,nw-1
              ifte(i,j)=ifte(i,nw)
   60       continue
          endif
          goto 52
        endif
   61 continue
      do 62 i=1,nlat-1
        call elnameK(i,name)
        enome=name.eq.wordp
        if( tmatch(pname(idcomp(elatt,i)),word) .or. enome) then
          if(.not.fcon) then
            fcon=.true.
            nwc=nw
          endif
          nw=nw+1
          wkey(nw)=wordp
          goto 52
        endif
   62 continue
c     write(*,'(a,I2,10A8)') 'nw=',nw,(wkey(i),i=1,nw)
c     write(*,*) (ifte(1,i),ifte(2,i),ifte(3,i),i=1,nw)
      m=0
      do 71 j=1,nw
        l=ifte(3,j)-1
        do 70 i=1,nlat-1
          call elnameK(i,name)
          enome=name.eq.wkey(j)
          id=idcomp(elatt,i)
          if( tmatch(pname(id),wkey(j)) .or. enome)then
            pat=idtype(id).eq.icmark .or. master(i).gt.0
            if(pat .or. enome) then
              if(mhogal(ifte(1,j),ifte(2,j),i)) then
c                print *,pnamec(i),i,ifte(1,j),ifte(2,j),
c    $               master(i)
                l=l+1
                if(l.eq.ifte(3,j)) then
                  m=m+1
                  l=0
                endif
              endif
            endif
          endif
   70   continue
   71 continue
      nb=m
c     print *,'nb=',nb
      return
      end

      subroutine mbstr(latt,twiss,gammab,mult,isb,xs,nbp,ncb,istr,estr,
     $     cor)
c     ---- solve bump equation ---------
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical cor(*),xplane,yplane
      integer*8 latt(nlat)
      dimension twiss(*),gammab(nlat),mult(*)
     z,         isb(ncb+4,*),xs(ncb+2,2,*),istr(*),estr(*)
      real*8 wk(ncb+2,2)
c
      save
c
      xplane=cor(1).or.cor(3).or.cor(5).or.cor(7)
      yplane=cor(2).or.cor(4).or.cor(6).or.cor(8)
      if(xplane) then
        do 120 j=1,nbp
          call msolb(latt,twiss,gammab,mult,isb(1,j),xs(1,1,j),ncb,istr,
     1               estr,wk,wk(1,2),.false.)
  120   continue
        if(yplane)then
          do 121 j=nbp+1,2*nbp
            call msolb(latt,twiss,gammab,mult,isb(1,j),xs(1,1,j),ncb,
     1                 istr,estr,wk,wk(1,2),.true.)
  121     continue
        endif
      elseif(yplane) then
        do 122 j=1,nbp
          call msolb(latt,twiss,gammab,mult,isb(1,j),xs(1,1,j),ncb,istr,
     1               estr,wk,wk(1,2),.true.)
c         write(*,'(1p,6e11.2)')(xs(i,1,j),i=1,isb(2,j))
c         write(*,'(1p,6e11.2)')(xs(i,2,j),i=1,isb(2,j))
  122   continue
      endif
c ------ debug ------------
c     do 500 i=1,isb(2,5)
c       n=istr(isb(i+2,5))
c       rlist(latt(2,n)+11)=rlist(latt(2,n)+11) + xs(i,2,5)*1d-3
c500  continue
c -------------------------
      return
      end

      subroutine mbufg(char,out,lfno)
cHITACHI   parameter (llen=80)
      parameter (llen=131)
      logical out
      character*(*) char
      character cline*(llen)
      data ip/1/,cline/' '/
      save ip,cline
      le=lene(char)
      if(ip+le-1.gt.llen .or. out) then
        if(cline.ne.' ') write(lfno,'(a)') cline
        cline=' '
        ip=1
      endif
      if(.not.out)then
        cline(ip:ip+le-1)=char(1:le)
c       ip=ip+le+1
        ip=ip+le
      endif
      return
      end

      subroutine mbufw(char,out,lfno)
      logical out
      character*(*) char
c
      call mbufg(' ',out,lfno)
      if(.not.out) call mbufg(char,out,lfno)
      return
      end

      subroutine mc2to4(twiss,ip,n,u,x)
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),u(4),x(4)
      r11=twiss(n,ip,mfitr1)
      r12=twiss(n,ip,mfitr2)
      r21=twiss(n,ip,mfitr3)
      r22=twiss(n,ip,mfitr4)
      det=r11*r22-r12*r21
      if(det.gt.1d0) then
        qr=sqrt((det-1d0)/det)
        un=sqrt(det)
        x(1)= (-r21*u(1)  -r22*u(2))*qr +    un*u(3)
        x(2)= (-r11*u(1)  -r12*u(2))*qr +               un*u(4)
        x(3)=    un*u(1)                + ( r12*u(3)  -r22*u(4))*qr
        x(4)=               un*u(2)     + (-r11*u(3)  +r21*u(4))*qr
      else
        un=sqrt(1d0-det)
        x(1)=   un*u(1)             +r22*u(3)  -r12*u(4)
        x(2)=              un*u(2)  -r21*u(3)  +r11*u(4)
        x(3)= -r11*u(1)  -r12*u(2)  + un*u(3)
        x(4)= -r21*u(1)  -r22*u(2)             + un*u(4)
      endif
      return
      end

      subroutine mc4to2(twiss,ip,n,c4,c2)
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),c4(4),c2(4)
      r11=twiss(n,ip,11)
      r12=twiss(n,ip,12)
      r21=twiss(n,ip,13)
      r22=twiss(n,ip,14)
      det=r11*r22-r12*r21
      if(det.gt.1d0) then
        qr=sqrt((det-1d0)/det)
        un=sqrt(det)
        c2(1)=(-r12*c4(1) + r22*c4(2))*qr +  un*c4(3)
        c2(2)=( r11*c4(1) - r21*c4(2))*qr              + un*c4(4)
        c2(3)=   un*c4(1)                 +(r21*c4(3)  +r22*c4(4))*qr
        c2(4)=               un*c4(2)     +(r11*c4(3)  +r12*c4(4))*qr
      else
        un=sqrt(1d0-det)
        c2(1)=  un*c4(1)               -r22*c4(3)  +r12*c4(4)
        c2(2)=            + un*c4(2)   +r21*c4(3)  -r11*c4(4)
        c2(3)= r11*c4(1)  +r12*c4(2)   + un*c4(3)
        c2(4)= r21*c4(1)  +r22*c4(2)               + un*c4(4)
      endif
      return
      end

      subroutine mcchar(ia,ib,n)
      integer*4 ia(n),ib(n)
      do 10 i=1,n
        call mcchar1(ia(i),ib(i))
   10 continue
      return
      end

      subroutine mcchar1(ia,ib)
      integer*4 ia,ib
      ib=ia
      return
      end

      subroutine mccor(word,wordp,
     $     latt,pos,twiss,
     $     gammab,mult,master,kfit,ifitp,
     $     mfitp,fitval,nfc,istr,estr,nstr,
     $     imon,emon,nmon,newcor,lfno)
      use tfstk
c -**** Main routine for orbit correction ****--------------------------
c  Note
c   . Subroutines for correction use the array TWISS as a buffer area:
c     twiss(*,ndim,*) is used as a 'base optics' for dp=0
c     twiss(*,ndim-1,*)      "        "          for dp=+ddp/2
c     twiss(*,1-ndim,*)      "        "          for dp=-ddp/2
c     twiss(*,-ndim,*)       "     'model optics' for dp=0 (in pbump,
c                     pasex, and palgn)
c     twiss(*,ndim-2,*)      "     back-up for dx,dy,ex,ey (if operate)
c   . imon(*,3) ;area for flag (>0:error 0:normal)
c     imon(*,4) ;element number which BPM is attached to.
c     istr(*,3) ;area for flag (>0:error 0:normal)
c ----------------------------------------------------------------------
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      use iso_c_binding
      implicit real*8 (a-h,o-z)
      parameter (nobj=8,nmeth=5,nsolvr=2,nbound=2,nother=nsolvr+nbound,
     1           nkey=nobj+nmeth+nother,ddp=1d-6, nbc0=4)
      parameter (mfitc1=32,mfitc2=28)
      real*8 ,pointer :: amat(:,:)
      character*(*) word,wordp,keywrd(nkey)*8,xy*4,line*79
      integer*4 ,pointer :: isbp(:,:)
      integer*4 coex(nmeth+nother,nmeth+nother)
      integer*8 isex,istb,istb1,ia,ib,ix,iz,iqc,iqd,
     $     itemp,iu,iv,ixb,iaad,iadd,iac,ibc,iaa,iqq
      logical stab,cod,dsp,normal,xplane,yplane,both,bump,minusi,micado,
     1        cond,zsum,exist,exec,coup,operate
      logical corc(nkey),corcn(nkey),corca(nobj),method(nobj,nmeth)
      integer*8 latt(nlat),l,iofs
      dimension pos(nlat),twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat),mult(nlat),master(nlat)
      dimension dp1(-ndim:ndim)
      dimension istr(nstra,4),imon(nmona,4),emon(*),estr(*)
      dimension kfit(nfc),ifitp(nfc),mfitp(nfc),fitval(nfc)
      real*8 , allocatable :: xs(:,:,:),wk(:,:)
      include 'inc/common.inc'
      data line/' '/
      data keywrd/'MX      ','MY      ','X       ','Y       ',
     1            'EX      ','EY      ','PEX     ','PEY     ',
     1            'SK      ','BUMP    ','SEX     ','ALIGN   ',
     1            '-I      ',
     1            'SVD     ','MICADO  ','COND    ','ZSUM    '/
      data corc  /.false.,   .false.,   .false.,   .false.,
     1            .false.,   .false.,   .false.,   .false.,
     1            .true.,    .false.,   .false.,   .false.,
     1            .false.,
     1            .true.,    .false.,   .false.,   .false./
      data coex  /1,0,0,0,1,1,1,1,1,
     1            0,1,0,0,1,1,1,1,1,
     1            0,0,1,0,1,1,1,1,1,
     1            0,0,0,1,1,1,1,1,1,
     1            1,1,1,1,1,1,1,1,1,
     1            1,1,1,1,1,1,0,1,1,
     1            1,1,1,1,1,0,1,1,1,
     1            1,1,1,1,1,1,1,1,1,
     &            1,1,1,1,1,1,1,1,1/
      data micado/.false./,nmicad/0/
      save corc
c     begin initialize for preventing compiler warning
      nbc=0
      iqc=0
      iqd=0
c     end   initialize for preventing compiler warning
      call cputime(ctime0,irtc)
      operate=.not.simulate
      do 2 i=1,nobj
      do 2 j=1,nmeth
        method(i,j)=.false.
2     continue
      do 5 i=1,nkey
    5   corcn(i)=.false.
      exec=.false.
1     call getwdl2(word,wordp)
      do 10 i=1,nkey
        if(word.eq.keywrd(i)) then
          corcn(i)=.true.
          if(word.eq.'MICADO') then
            nmicad=int(getva(exist))
c7/30-------------->
          elseif(word.eq.'ALIGN')then
            lstack=int(getva(exist))
            if(.not.exist) lstack=0
c7/30<-------------
          elseif(word.eq.'BUMP')then
            nbc=int(getva(exist))
            if(.not.exist) nbc=nbc0
          endif
          go to 1
        elseif(word.eq.'NO'//keywrd(i).or.word.eq.'NON'//keywrd(i)) then
          corc(i)=.false.
          goto 1
        endif
   10 continue
      do 12 i=1,nobj
        if(corcn(i)) then
           exec=.true.
           do 11 j=1,nobj
             corc(j)=corcn(j)
   11      continue
           exit
        endif
   12 continue
      do 15 i=1,nmeth+nother
        if(corcn(nobj+i)) then
          corc(nobj+i)=.true.
          do 14 j=1,nmeth+nother
            if(coex(i,j).eq.0) corc(nobj+j)=.false.
   14     continue
        endif
   15 continue
      cod=corc(1).or.corc(2).or.corc(3).or.corc(4)
      dsp=corc(5).or.corc(6).or.corc(7).or.corc(8)
      normal=corc(1).or.corc(2).or.corc(5).or.corc(6)
      yplane=corc(2).or.corc(4).or.corc(6).or.corc(8)
      xplane=corc(1).or.corc(3).or.corc(5).or.corc(7)
      both=xplane.and.yplane
      bump=corc(10)
      minusi=corc(13)
      micado=corc(nobj+nmeth+2)
      cond=corc(nobj+nmeth+3)
      zsum=corc(nobj+nmeth+4)
      eptsol0=eptsol
      stab=.true.
c     ---Read fowllowing words if word='BUMP' or 'SEX' ---
      if(corcn(nobj+2) .or. corcn(nobj+3)) then
        call mbmpf(word,wordp,latt,mult,master,nbp)
      endif
      do 21 i=1,nkey
        j=min(5,lene(keywrd(i)))
        call mbufw(keywrd(i)(1:j),.false.,lfno)
   21 continue
      call mbufw(' ',.true.,lfno)
      k=1
      do 22 i=1,nkey
        j=min(5,lene(keywrd(i)))
        write(line(k+(j+1)/2:),'(L1)') corc(i)
        k=k+j+1
   22 continue
      write(lfno,'(a)') line
      if(.not.exec) return
c***** Execute correction ****
c 1001 continue
      do 24 i=1,nobj
        corca(i)=corc(i)
        do 25 j=1,nmeth
   25     method(i,j)=corc(j+nobj)
   24 continue
      if( minusi ) then
        call ppair1(latt,1,nlat,icsext,npair)
        write(lfno,'(A,I4)')' No. of sext. pairs (total ring) =',npair
        isex=ktaloc(npair)
        call ppair(latt,1,nlat,icsext,rlist(isex))
        call packpi(rlist(isex),npair,istr,istr(1,2),istr(1,3),nstr,
     1               nstr,.true.)
      endif
      if( corc(5).and.corc(6) .or. corc(7).and.corc(8) ) then
        ida=nmon*2
        xy='XY'
      elseif( corc(5).or.corc(7) ) then
        ida=nmon
        xy='X'
      elseif( corc(6).or.corc(8) ) then
        ida=nmon
        xy='Y'
      else
        ida=0
      endif
      if( corc(1).and.corc(2) .or. corc(3).and.corc(4) ) then
        ica=nmon*2
        xy='XY'
      elseif( corc(1).or.corc(3) ) then
        ica=nmon
        xy='X'
      elseif( corc(2).or.corc(4) ) then
        ica=nmon
        xy='Y'
      else
        ica=0
      endif
      mdim1=ica+ida
      if(operate) then
        if(ica.ne.0) then
          ica=2*nmon
        endif
        if(ida.ne.0) then
          ida=2*nmon
        endif
        xy='XY'
      endif
      nc=0
      if(cond) then
        do 26 k=1,nfc
          if(mfitp(k).ne.0) then
            if((kfit(k) .ge. mfitc1 .and. kfit(k) .le. mfitc1+3)  .or.
     1         (kfit(k) .ge. mfitc2 .and. kfit(k) .le. mfitc2+3))then
              nc=nc+1
            endif
          endif
   26   continue
      endif
      if(nc.eq.0) cond=.false.
      i=mcoupsten()
      coup=i.ne.0
      nc=nc+i
      if(zsum) then
        do i=1,nstr
          zsum=.false.
          l=idval(idcomp(elatt,istr(istr(i,2),1)))
          if(rlist(l+1).ne.0d0 .and. rlist(l+2).ne.0d0) then
            zsum=.true.
            exit
          endif
        enddo
      endif
      if(zsum) nc=nc+1
      mdim=ica+ida+nc
      if(mdim.eq.0) return
      mdim1=mdim1+nc
      do kkk=1,1
        if( bump ) then
          nbco=nbc
          if(nbp.eq.0) exit
          if(both) then
            nbump=nbp*2
          else
            nbump=nbp
          endif
          istb=ktaloc(((nbco+4)*nbump+1)/2)
c          ixs =ktaloc(2*(nbco+2)*nbump)
          allocate(xs(nbco+2,2,nbump))
        elseif( corc(nobj+3) ) then
          nbco=4
          if(nbp.eq.0) exit
          istb=ktaloc(((nbco+4)*nbp+1)/2)
          istb1=ktaloc(((nbco+4)*nbp+1)/2)
c          ixs =ktaloc(4*nbco*nbp)
          allocate(xs(nbco,4,nbp))
        endif
c     7/30--------------->
        if(corc(12)) then
c     'ALIGN'
          ia=0
          ib=0
          ix=0
          iz=0
c     7/30<--------------
        else
          ia=ktaloc(mdim*nstr)
          ib=ktaloc(mdim)
          ix=ktaloc(nstr)
          iz=ktaloc(4*nmon+nc)
        endif
c     ----- obtain optics
        dp1(ndim-1)= ddp/2d0 + dp0 +1d0
        dp1(0 )=               dp0 +1d0
        dp1(1-ndim)=-ddp/2d0 + dp0 +1d0
        if(newcor.ne.-1) then
          do 30 j=1,ntwissfun
            call tmov(twiss(1,0,j),twiss(1,ndim,j),nlat)
            twiss(1,ndim-1,j)=twiss(1,0,j)
            twiss(1,1-ndim,j)=twiss(1,0,j)
 30       continue
          call pqcell(latt,ndim-1,dp1(ndim-1),stab)
          call pqcell(latt,1-ndim,dp1(1-ndim),stab)
          newcor=1
        endif
        call pcbak(latt,twiss)
c     call prkick(latt,mult,istr,nstr,.false.,.false.,lfno)
        call mcrcod(latt,twiss,mult,imon,emon,nmon,.true.,.false.,
     1       .false.,lfno)
        irow=0
        do 1010 ic=1,nobj
          if( .not.corca(ic) ) cycle
          if( method(ic,1) ) then
c     ---- SK ----
            if(cod) then
              psix=twiss(nlat,0,3)-twiss(1,0,3)
              psiy=twiss(nlat,0,6)-twiss(1,0,6)
              call c_f_pointer(c_loc(rlist(ia+irow)),amat,[mdim,nstr])
              call mcrmat(latt,twiss,gammab,ndim,psix,psiy,amat,
c     $             rlist(ia+irow),
     1             mdim,.false.,normal,istr,nstr,imon,nmon,xy)
            elseif(dsp) then
              psix=0d0
              psiy=0d0
              iadd=ia+irow-mdim
              k=istr(1,2)
              do 31 j=1,nstr
                istr(1,2)=istr(j,2)
                iadd=iadd+mdim
                call pcrmat(latt,twiss,gammab,psix,psiy,
     $               rlist(iadd),ida,
     1               istr,1,imon,nmon,.false.,normal,xy,
     1               rlist(iz),ddp)
 31           continue
              istr(1,2)=k
            endif
          elseif(method(ic,2)) then
c     ---- BUMP ----
            if(cod) then
              psix=twiss(nlat,0,3)-twiss(1,0,3)
              psiy=twiss(nlat,0,6)-twiss(1,0,6)
              call c_f_pointer(c_loc(rlist(ia+irow)),amat,[mdim,nstr])
              call mcrmat(latt,twiss,gammab,ndim,psix,psiy,amat,
c     $             rlist(ia+irow),
     1             mdim,.false.,normal,istr,nstr,imon,nmon,xy)
c     ---- Calc steerings of unit bumps.
              if(both) then
                call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp])
                call mbmp(latt,twiss,mult,master,isbp,nbp,istr,
     1               nstr,nbco,.false.,iret,lfno)
                k=nbp*(nbco+4)
                call c_f_pointer(c_loc(ilist(mod(k,2)+1,
     $               istb+k/2)),isbp,[nbco,nbp])
                call mbmp(latt,twiss,mult,master,isbp,
     $               nbp,istr,nstr,nbco,.true.,iret,lfno)
              else
                call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp])
                call mbmp(latt,twiss,mult,master,isbp,nbp,istr,
     1               nstr,nbco,corc(2).or.corc(4),iret,lfno)
              endif
c              iwk =ktaloc(2*(nbco+2))
c     print *,'free before mbstr=',mfalloc(-1)
              call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp])
              call mbstr(latt,twiss,gammab,mult,isbp,
     $             xs,nbp,nbco,istr,estr,corc)
c     ---- calc bump response
              iqc =ktaloc(mdim*2*nbump)
              call mrmb(ilist(1,istb),xs,rlist(iqc),
     z             nbco,nbump,rlist(ia+irow),mdim)
            elseif(dsp) then
              psix=0d0
              psiy=0d0
              iadd=ia+irow-mdim
              k=istr(1,2)
              do 41 j=1,nstr
                istr(1,2)=istr(j,2)
                iadd=iadd+mdim
                call pcrmat(latt,twiss,gammab,psix,psiy,rlist(iadd),ida,
     1               istr,1,imon,nmon,.false.,normal,xy,
     1               rlist(iz),ddp)
 41           continue
              istr(1,2)=k
c     ---- calc steerings for unit bumps
              if(both) then
                call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp])
                call mbmp(latt,twiss,mult,master,isbp,nbp,istr,
     1               nstr,nbco,.false.,iret,lfno)
                k=nbp*(nbco+4)
                call mbmp(latt,twiss,mult,master,isbp,nbp,istr,
     $               nstr,nbco,.true.,iret,lfno)
              else
                call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp])
                call mbmp(latt,twiss,mult,master,isbp,nbp,istr,
     1               nstr,nbco,corc(6).or.corc(8),iret,lfno)
              endif
              call mbstr(latt,twiss,gammab,mult,isbp,
     $             xs,nbp,nbco,istr,estr,corc)
c     ---- calc bump response
              iqd =ktaloc(mdim*2*nbump)
              call mrmb(isbp,xs,rlist(iqd),nbco,nbump,
     z             rlist(ia+irow),mdim)
            endif
          elseif(method(ic,3)) then
c     ---- SEX ----
            if(cod) then
              nbp1=nbp
              call c_f_pointer(c_loc(ilist(1,istb)),isbp,[nbco+4,nbp1])
              call mbmp(latt,twiss,mult,master,isbp,nbp1,istr,nstr,
     z             nbco,.false.,iret,lfno)
              nbp1=nbp
              call mbmp(latt,twiss,mult,master,rlist(istb1),nbp1,istr,
     1             nstr,nbco,.true.,iret,lfno)
              nbp=nbp1
              call tfree(istb1)
c     elseif(dsp) then
            endif
          elseif(method(ic,4)) then
c     ---- ALIGN ----
            if(cod) then
              if( minusi ) then
c     call palgn(word,wordp,latt,twiss,gammab,mult,rlist(isex),npair,
c     z                   istr,istr(1,2),estr,nstr,imon,
c     z                   nmon,corc,lfno)
              else
                call pasex(word,wordp,latt,pos,twiss,mult,master,
     $               gammab,' ',' ',istr,estr,nstr,imon,emon,nmon,
     $               lfno)
              endif 
              exit
            endif
          endif
c     93/10/28
          if(both .or. operate) then
            corca(ic+1)=.false.
            irow=irow+2*nmon
          else
            irow=irow+nmon
          endif
 1010   continue
c     .... Take into account conditions in the matching conditions buffer ..
        jc=0
        if(cond) then
          itemp=ktaloc(nmona)
          ilist(mod(nmona,2)+1,itemp+nmona/2)=1
          iac=ia+mdim-nc
          ibc=ib+mdim-nc
          do 50 k=1,nfc
            if(mfitp(k).ne.0) then
              ilist(1,itemp)=ifitp(k)
              if(kfit(k) .eq. mfitc1) then
c     -- 'DX' ---
                call mcrmat(latt,twiss,gammab,ndim,psix,psiy,
     1               rlist(iac+jc),mdim,
     1               .false.,.false.,istr,nstr,rlist(itemp),1,'X')
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc1+1) then
c     -- 'DPX' ---
                call mcrmat(latt,twiss,gammab,ndim,psix,psiy,
     1               rlist(iac+jc),mdim,
     z               .true.,.false.,istr,nstr,rlist(itemp),1,'X')
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc1+2) then
c     -- 'DY' ---
                call mcrmat(latt,twiss,gammab,ndim,psix,psiy,
     1               rlist(iac+jc),mdim,
     z               .false.,.false.,istr,nstr,rlist(itemp),1,'Y')
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc1+3) then
c     -- 'DPY' ---
                call mcrmat(latt,twiss,gammab,ndim,psix,psiy,
     1               rlist(iac+jc),mdim,
     z               .true.,.false.,istr,nstr,rlist(itemp),1,'Y')
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc2) then
c     -- 'DEX' ---
                call pcrmat(latt,twiss,gammab,psix,psiy,
     1               rlist(iac+jc),mdim,istr,nstr,rlist(itemp),1,
     1               .false.,.false.,'X',rlist(iz),ddp)
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc2+1) then
c     -- 'DEPX' ---
                call pcrmat(latt,twiss,gammab,psix,psiy,
     1               rlist(iac+jc),mdim,istr,nstr,rlist(itemp),1,
     1               .true.,.false.,'X',rlist(iz),ddp)
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc2+2) then
c     -- 'DEY' ---
                call pcrmat(latt,twiss,gammab,psix,psiy,
     1               rlist(iac+jc),mdim,istr,nstr,rlist(itemp),1,
     1               .false.,.false.,'Y',rlist(iz),ddp)
                rlist(ibc+jc)=fitval(k)
              elseif(kfit(k) .eq. mfitc2+3) then
c     -- 'DEPY' ---
                call pcrmat(latt,twiss,gammab,psix,psiy,
     1               rlist(iac+jc),mdim,istr,nstr,rlist(itemp),1,
     1               .true.,.false.,'Y',rlist(iz),ddp)
                rlist(ibc+jc)=fitval(k)
c     elseif(kfit(k) .eq. 3) then
c     -- 'NX' ---
c     elseif(kfit(k) .eq. 6) then
c     -- 'NY' ---
              else
                cycle
              endif
              jc=jc+1
            endif
 50       continue
          call tfree(itemp)
c     print *,' nc=',nc,' mdim=',mdim,' mdim1=',mdim1
c     write(*,'(1p,5d11.3)')(rlist(ibc+i),i=0,j)
c     <------Check of Condition Matrix-------
c     print *,' Response matrix: steering number:'
c     read(*,*) nofste
c     write(*,'(1p,5e11.3)')
c     1       (rlist(iac+k),k=mdim*(nofste-1),mdim*(nofste-1)+nc-1)
c     -------Check of Condition Matrix------>
        endif
        iac=ia+mdim-nc+jc
        ibc=ib+mdim-nc+jc
        if(zsum) then
          do 51 i=1,nstr  
            l=idval(idcomp(elatt,istr(istr(i,2),1))-1)
            rlist(iac+(i-1)*mdim)=rlist(l+1)*rlist(l+2)
 51       continue
          rlist(ibc)=0d0
          iac=iac+1
          ibc=ibc+1
        endif
        if(coup) then
          n=mcoupsten()
          itemp=ktaloc(n*(nstr+1))
          call mcoupstea(rlist(itemp),n,latt,mult,istr,nstr,lfno)
          do j=1,n
            do i=1,nstr
              rlist(iac+(i-1)*mdim)=rlist(itemp+(j-1)*(nstr+1)+i-1)
            enddo
            rlist(ibc)=rlist(itemp+(j-1)*(nstr+1)+nstr)
            iac=iac+1
            ibc=ibc+1
          enddo
          call tfree(itemp)
        endif
c     ----- solve equation ---------
        if(psix.eq.0d0) then
          psix=twiss(nlat,ndim,3)-twiss(1,ndim,3)
        endif
        if(psiy.eq.0d0) then
          psiy=twiss(nlat,ndim,6)-twiss(1,ndim,6)
        endif
        write(lfno,'(2(A,F12.5))')
     1       '  Tune used (x/y) =',psix/pi2,' /',psiy/pi2
c     <----Check of Response Matrix-------
c     print *,' nc=',nc,' mdim=',mdim,' mdim1=',mdim1
c     print *,' Response matrix: steering number:'
c     read(*,*) nofste
c     write(*,'(1p,10e11.3)')
c     $     (rlist(ia+k),k=mdim*(nofste-1),mdim*nofste-1)
c     (rlist(iqd+k),k=mdim*(nofste-1),mdim*nofste-1)
c     write(*,'(2(1p,2e15.7/))') rlist(ia),rlist(ia+2),rlist(ia+1),
c     $     rlist(ia+3)
c     -----Check of Response Matrix------>
        call mcnrmc(latt,twiss,rlist(ib),
     $       normal,imon,emon,nmon,corc,lfno)
c     print *,' b:'
c     write(*,'(1p,10e11.3)')(rlist(ib+k),k=0,mdim-1)
c     write(*,'(1p,2e15.7)')rlist(ib),rlist(ib+1)
        if(bump) then
          nvar=2*nbump
          if(cod) then
            iaa=iqc
          else
            iaa=iqd
          endif
        else
          nvar=nstr
          iaa=ia
        endif
        if(operate) then
          iu=ktaloc(mdim)
          iv=ktaloc(mdim*nvar)
          call tmov(rlist(iaa),rlist(iv),mdim*nvar)
          call preduce(rlist(iaa),mdim1,rlist(iv),mdim,nc,nmon,nvar,
     $         corc,nobj)
c     <----Check of Response Matrix-------
c     print *,' nc=',nc,' mdim=',mdim,' mdim1=',mdim1
c     print *,' Response matrix: steering number:'
c     read(*,*) nofste
c     write(*,'(1p,10e11.3)')
c     1       (rlist(iaa+k),k=mdim1*(nofste-1),mdim1*nofste-1)
c     -----Check of Response Matrix------>
        endif
        write(lfno,'(2(a,i5))')
     1       ' BPMs using =',nmon,' Steerings using =',nstr
        call cputime(ctime1,irtc)
        if( bump ) then
          ixb=ktaloc(nvar)
          if( cod ) then
            iqq=iqc
          else
            iqq=iqd
          endif
          call mwght(twiss,ndim,rlist(iqq),rlist(ib),mdim1,nvar,mdim1,
     1         imon,nmon,corc)
          call msolvg(rlist(iqq),rlist(ib),rlist(ixb),mdim1,nvar,mdim1,
     $         nc,cond.or.zsum.or.coup,micado,nmicad,.true.,
     $         .false.)
          write(*,'(1p,10d12.4)') (rlist(ixb-1+i),i=1,nvar)
          call tfree(iqq)
          call pbset(latt,mult,rlist(ixb),ilist(1,istb),xs,
     1         istr,estr,nstr,nbco,nbump,.true.,lfno)
          if(operate) then
            call mbexpect(twiss,imon,nmon,
     $           ilist(1,istb),xs,rlist(ixb),nbco,nbump,
     $           rlist(ia),mdim,corc,nobj)
          endif
          call tfree(ixb)
          call tfree(istb)
        else
          call mwght(twiss,ndim,rlist(ia),rlist(ib),mdim1,nstr,mdim1,
     1         imon,nmon,corc)
          call msolvg(rlist(ia),rlist(ib),rlist(ix),mdim1,nstr,mdim1,nc,
     $         cond.or.zsum.or.coup,micado,nmicad,.true.,
     $         .false.)
c     print *,'x:'
c     write(*,'(1p,6e11.3)')(rlist(ix+k),k=0,nstr-1)
c     call pcbak(latt,twiss)
c     -------------------
          call pcset(latt,mult,rlist(ix),istr,estr,nstr,
     1         1d0,.true.,lfno)
c     -------------------
          if(operate)then
            call pinner(rlist(iv),rlist(ix),rlist(iu),mdim,nvar,mdim)
            iofs=iu-1
            if(cod) then
              do 901 i=1,nmon
                j=imon(i,2)
                twiss(imon(j,1),0,15)=
     $               twiss(imon(j,1),0,15)-rlist(iofs+i)
                twiss(imon(j,1),0,17)=twiss(imon(j,1),0,17)
     1               - rlist(iofs+nmon+i)
 901          continue
              iofs=iofs+2*nmon
            endif
            if(dsp) then
              do 902 i=1,nmon
                j=imon(i,2)
                twiss(imon(j,1),0,7)=twiss(imon(j,1),0,7)-rlist(iofs+i)
                twiss(imon(j,1),0,9)=twiss(imon(j,1),0,9)
     1               - rlist(iofs+nmon+i)
 902          continue
            endif
          endif
        endif
        if(operate)then
          call tfree(iu)
          call tfree(iv)
        endif
      enddo
      call cputime(ctime2,irtc)
      if(simulate) then
        if(.not.corc(12)) then
c       --'not ALIGN'---
          if(cell) then
            do 910 j=15,18
              twiss(1,0,j)=optiv(j)
 910        continue
          endif
          call pqcell(latt,twiss,gammab,0,dp1(0),stab)
        endif
      endif
c ----- unpack -------
      if( minusi ) then
        call packpi(rlist(isex),npair,istr,istr(1,2),istr(1,3),nstr,
     1              nstra,.false.)
        call tfree(isex)
      endif
c --------------------
      if(iz.ne.0)call tfree(iz)
      if(ix.ne.0)call tfree(ix)
      if(ib.ne.0)call tfree(ib)
      if(ia.ne.0)call tfree(ia)
      if(.not.stab) then
c       eptsol=eptsol*10
c       if(eptsol.lt.1.d0) then
          call pundo(latt,twiss,gammab,0,dp1(0))
c         print *,' EPSilon=',sngl(eptsol)
c         goto 1001
c       endif
      endif
      eptsol=eptsol0
      call prkick(latt,mult,istr,nstr,.true.,.true.,lfno)
      call mcrcod(latt,twiss,mult,imon,emon,nmon,.true.,.true.,
     1            .true.,lfno)
      call cputime(ctime3,irtc)
      write(lfno,'(t10,2(A,F10.6))')
     1                           ' ctime(solver)=',1d-6*(ctime2-ctime1),
     1                           ' ctime(mccor)=',1d-6*(ctime3-ctime0)
      return
      end
c

      subroutine preduce(a,ma,b,mb,nc,nmon,nvar,corc,nobj)
      implicit none
      logical corc
      real*8 a,b
      integer ma,mb,nc,nmon,nvar,nobj,j,k,ia,ib
      dimension a(ma,nvar),b(mb,nvar),corc(nobj)
c
      do j=1,nvar
        ib=-1
        ia=-2
        do k=1,nobj
          if(corc(k)) then
            if(mod(k,2).eq.1) then
              ia=ia+2
            elseif(.not.corc(k-1)) then
              ia=ia+2
            endif
            ib=ib+1
            call tmov(b(1+(ia+mod(k+1,2))*nmon,j),a(1+ib*nmon,j),nmon)
          endif
        enddo
        do k=nc,1,-1
          a(ma+1-k,j)=b(mb+1-k,j)
        enddo
      enddo
c
      return
      end
c

      subroutine mbexpect(twiss,imon,nmon,is,xs,x,ncb,nbp,a,ia,corc,
     $     nobje)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical corc(nobje)
      integer*8 iw
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      dimension imon(nmona,4)
      dimension is(ncb+4,nbp),xs(ncb+2,2,nbp),x(2,nbp),a(ia,*)
      include 'inc/common.inc'
c
      iw=ktaloc(ia)
      call pclr(rlist(iw),ia)
      do j=1,nbp
        ncor=is(2,j)
        do i=1,ia
          s=0d0
          t=0d0
          do k=1,ncor
            s=s+a(i,is(k+2,j))*xs(k,1,j)
            t=t+a(i,is(k+2,j))*xs(k,2,j)
          enddo
          rlist(iw-1+i)=rlist(iw-1+i) + s*x(1,j) + t*x(2,j)
        enddo
      enddo
      j=-2
      do i=1,nobje,2
        if(corc(i).or.corc(i+1)) then
          j=j+2
          if(i.le.4) then
            do k=1,nmon
              twiss(imon(imon(k,2),1),0,15)=
     $             twiss(imon(imon(k,2),1),0,15)
     $             - rlist(iw-1+j*nmon+k)
              twiss(imon(imon(k,2),1),0,17)=
     $             twiss(imon(imon(k,2),1),0,17)
     $             - rlist(iw-1+(j+1)*nmon+k)
            enddo
          else
            do k=1,nmon
              twiss(imon(imon(k,2),1),0,7)=
     $             twiss(imon(imon(k,2),1),0,7)
     $             - rlist(iw-1+j*nmon+k)
              twiss(imon(imon(k,2),1),0,9)=
     $             twiss(imon(imon(k,2),1),0,9)
     $             - rlist(iw-1+(j+1)*nmon+k)
            enddo
          endif
        endif
      enddo
      call tfree(iw)
      return
      end

      subroutine  mcepst(word,lfno)
      implicit real*8 (a-h,o-z)
      logical exist
      character*(*) word
      include 'inc/common.inc'
c     precision of TSOLVA
      eps=getva(exist)
      if(exist) then
        eptsol=eps
      endif
      call getwdl(word)
      return
      end

      subroutine mcfix(latt,twiss,gammab,newcor)
c              0: correction is initiated (first FIX or first COR)
c              1: current optics used
c             -1: optics is fixed by FIX
      use ffs
      use tffitcode
      use ffs_fit ,only: ffs_stat
      implicit real*8 (a-h,o-z)
      parameter (ddp=1d-6,epsdlt=1d-6)
c BOTH HP and DEC compiler allows dynamic array size if these arrays are
c  automatic (NY )
      type (ffs_stat) optstat(-1:1)
      integer*8 latt(nlat)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
      real*8 dp1(-ndim:ndim)
c
c     write(*,*) ' newcor in mcfix',newcor
      dp1( 1)= ddp/2d0 + dp0 +1d0
      dp1( 0)=           dp0 +1d0
      dp1(-1)=-ddp/2d0 + dp0 +1d0
      do  ip=-1,1,2
        do  j=1,ntwissfun
          twiss(1,ip,j)=twiss(1,0,j)
        enddo
      enddo
c93/11/01
      do ip=-1,1,1
        if( cell ) then
          optstat(ip)%stabx=.true.
          optstat(ip)%staby=.true.
          call qcell(ip,optstat(ip),.false.)
          if( .not. optstat(ip)%stabx .or.
     $         .not. optstat(ip)%staby ) then
            write(*,'(A,1PD12.5,A/2(A,F9.4))')
     z           '  Unstable orbit  dp=',dp1(ip)-1d0,' :',
     z           '          tracex=',optstat(ip)%tracex,
     $           ' tracey=',optstat(ip)%tracey
          endif
        else
          twiss(1,0,3)=0d0
          twiss(1,0,6)=0d0
          call qtwiss(twiss,ip,1,nlat,over)
        endif
      enddo
c     ------- save current cod & dispersion as standards to which the
c             cod and the dispersion should be corrected.
      do j=1,ntwissfun
        call tmov(twiss(1,0,j),twiss(1,ndim,j),nlat)
        call tmov(twiss(1,1,j),twiss(1,ndim-1,j),nlat)
        call tmov(twiss(1,-1,j),twiss(1,1-ndim,j),nlat)
      enddo
c       ------
      newcor=-1
      return
      end
c   --------------------------

      subroutine mcfre(newcor)
      integer*4 newcor
      if(newcor.ne.0) then
        newcor=1
      endif
      return
c   --------------------------
      end

      subroutine mckick(word,wordp,latt,mult)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
      character*(*) word,wordp
      character*(MAXPNAME) word1
c     character*255 tfconvstr
      logical exist
      integer*8 latt(nlat)
      dimension mult(nlat)
c
      save
c
c
c     itype=itfpeeko(ia,x,next)
c     if(itype.eq.1) then
c       call cssetp(next)
c       xk=x
c       print *,xk
c     else
c       call getwdl(word)
c       return
c     endif
c     itype=itfpeeko(ia,x,next)
c     print *,itype
c     if(itype.eq.1) then
c       print *,itype
c       call gtewdl(word)
c       return
c     elseif(itype.eq.101) then
c       print *,itype
c       call cssetp(next)
c       word=tfconvstr(101,ia,x,nc,'*')
c       if(word.eq.' ') then
c         call permes('?Missing element for KICK','.',' ',6)
c         call getwdl(word)
c         return
c       endif
c     endif
      xk=getva(exist)
      if( .not.exist ) then
        call getwdl2(word,wordp)
        return
      endif
      call getwdl2(word,wordp)
      npnt=ielm(wordp,exist)
      if( exist ) then
        if(idtypec(npnt).eq.2) then
          rlist(latt(npnt)+11)=rlist(latt(npnt)+11)+xk
        else
          word1=wordp
          call permes(' ???',' '//word1,' is not a BEND.',6)
        endif
      endif
      call getwdl2(word,wordp)
      return
      end

      subroutine mclear(word,latt,istr,nstr)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      logical tmatch
      character*(*) word
      integer*8 latt(nlat)
      dimension istr(nstra,4)
      include 'inc/common.inc'
c
      save
c
c
    1 call getwdl(word)
      do 10 i=1,nstr
        if(tmatch(pname(idcomp(elatt,istr(istr(i,2),1))),word)) then
          goto 11
        endif
   10 continue
      return
   11 do 12 i=1,nstr
        if(tmatch(pname(idcomp(elatt,istr(istr(i,2),1))),word)) then
          rlist(latt(istr(istr(i,2),1))+11)=0d0
        endif
   12 continue
      goto 1
      end

      subroutine mclin(a,b,xa,xb,n,iop)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      goto (1,2,3,4,5),iop
c   1 if(n.gt.16) then
    1   do 10 i=1,n
          b(i)=xb*b(i) + a(i)
   10   continue
c     else
c*VOPTION NOVEC
c       do 12 i=1,n
c         b(i)=xb*b(i) + a(i)
c  12   continue
c     endif
      return
c   2 if(n.gt.16) then
    2   do 14 i=1,n
          b(i)=xb*b(i) - a(i)
   14   continue
c     else
c*VOPTION NOVEC
c       do 16 i=1,n
c         b(i)=xb*b(i) - a(i)
c  16   continue
c     endif
      return
c   3 if(n.gt.16) then
    3   do 20 i=1,n
          b(i)=b(i) + xa*a(i)
   20   continue
c     else
c*VOPTION NOVEC
c       do 22 i=1,n
c         b(i)=b(i) + xa*a(i)
c  22   continue
c     endif
      return
c   4 if(n.gt.16) then
    4   do 24 i=1,n
          b(i)=-b(i) + xa*a(i)
   24   continue
c     else
c*VOPTION NOVEC
c       do 26 i=1,n
c         b(i)=-b(i) + xa*a(i)
c  26   continue
c     endif
      return
c   5 if(n.gt.16) then
    5   do 30 i=1,n
          b(i)=xb*b(i) + xa*a(i)
   30   continue
c     else
c*VOPTION NOVEC
c       do 32 i=1,n
c         b(i)=xb*b(i) + xa*a(i)
c  32   continue
c     endif
      return
      end

      subroutine mcmon(word,latt,mult,pos,master,ipmon,ipemon,nmon,
     &                 initial,lfno)
      use tfstk
      use ffs
      use tffitcode
      logical initial
      character*(*) word
      integer*8 latt(nlat)
      dimension pos(nlat)
      integer*4 master(nlat),mult(*)
      include 'inc/common.inc'
      if(itmon.eq.0) then
        k=0
        do 10 i=1,nlat-1
          if(master(i).ge.0) then
            k=k+1
          endif
   10   continue
        nmona=k
        ipmon =ktaloc(2*nmona)
        ipemon=ktaloc(4*nmona)
        itmon=ipmon
        itemon=ipemon
        j=0
        do 11 i=1,nlat-1
          if(master(i).ge.0) then
c           master(i)=i+nslave; if master (that means it is a head AND
c                                 idtype=2,4,6,...(i.e., magnet) )
c                    =-i+1; if slave
c                    =0; if not magnet
            ilist(mod(j,2)+1,ipmon+int(j/2))=i
            ilist(mod(j+nmona,2)+1,ipmon+int((nmona+j)/2))=j+1
            ilist(mod(j+3*nmona,2)+1,ipmon+int((3*nmona+j)/2))=
     1            monqu(latt,pos,i)
            j=j+1
          endif
   11   continue
        nmon=0
        call pclr(rlist(ipemon),4*nmona)
      endif
      if(initial) return
c     .. mcmon is called by preadmon with initial=.true. .......
      call mcmon1(word,latt,mult,rlist(itmon),rlist(itemon),nmon,lfno)
      return
      end

      subroutine  mcmon1(word,latt,mult,imon,emon,nmon,lfno)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      implicit real*8 (a-h,o-z)
      character*(*) word
      character*16 en
      logical tmatch, exist
      integer*8 latt(nlat)
      dimension mult(*)
      dimension imon(nmona,4),emon(nmona,4)
      external pack
      parameter (nvpar=4)
      include 'inc/common.inc'
c  ---------------------------------------------------------
c             ( nvpar: no.of vaules expected before Mname
c              MON V1 V2 ... Vn Mname : assumed format
c                  if Vn=@ previous errors are kept, and no
c                          errors are generated.
c       V1=rms of positioning error for x.
c       V2=rms of precision for x.
c       V3=rms of positioning error for y.
c       V4=rms of precision for y.
c           frandom errors are generated in each measurment.
c  ---------------------------------------------------------
      nm=0
      na=0
      do 7 i=1,nmona
        imon(i,3)=1
    7 continue
      do 10 k=1,nvpar
        errv=getva(exist)
        if(exist) then
          errval(k)=errv
          errflg(k)=.true.
        else
          call getwdl(word)
          if(word.eq.'@') then
            errflg(k)=.false.
          else
            do 8 i=k,nvpar
              errflg(i)=.false.
              errval(i)=0d0
    8       continue
            do 9 i=k,nvpar
              if(i.ge.3)then
                errflg(i)=errflg(i-2)
                errval(i)=errval(i-2)
              endif
    9       continue
            goto 12
          endif
       endif
   10 continue
   11 call getwdl(word)
   12 do 13 i=1,nmona
        if( tmatch(pname(idcomp(elatt,imon(imon(i,2),1))),word) )then
          imon(imon(i,2),3)=0
          nm=nm+1
        endif
   13 continue
      if(nm.ne.na) then
        na=nm
        goto 11
      else
        do 14 i=1,nmona
          call elnameK(imon(imon(i,2),1),en)
          if(en.eq.word) then
            imon(imon(i,2),3)=0
            nm=nm+1
            na=nm
            goto 11
          endif
 14     continue
      endif
      call pack(imon(1,2),imon(1,3),nmon,nmona)
      nmonact=nmon
      if(nmon.eq.0)then
        write(lfno,'(2(A,I4))')
     1                   '  BPM not available: ',nmon,' in ',nmona
      else
        write(lfno,'(2(A,I4))') '  BPM available: ',nmon,' in ',nmona
      endif
c
      do 22 j=1,4,2
        if(errflg(j))then
          do 20 k=1,nmon
            emon(imon(k,2),(j+1)/2)=errval(j) *tgauss()
   20     continue
        elseif(errval(j).eq.0d0) then
          do 21 k=1,nmon
            emon(imon(k,2),(j+1)/2)=0d0
   21     continue
        endif
   22 continue
      if(errval(1).ne.0d0 .or. errval(2).ne.0d0)then
        write(lfno,*) '  -- Statistics of BPM alignment errors --'
        call mstat2(emon(1,1),imon(1,2),nmon,tmax,tmin,trms,tave)
        tmax=max(tmax,-tmin)
        write(lfno,'(3(A,1PD12.4))') '  dxm (max/rms/mean) =',
     1                          tmax,' /',trms,' /',tave
        call mstat2(emon(1,2),imon(1,2),nmon,tmax,tmin,trms,tave)
        tmax=max(tmax,-tmin)
        write(lfno,'(3(A,1PD12.4))') '  dym (max/rms/mean) =',
     1                          tmax,' /',trms,' /',tave
      endif
      return
      end

      subroutine  mcnrmc(latt,twiss,b,normal,imon,emon,nmon,corc,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      parameter (nkey=8, ddp=1d-6)
      parameter (sqr2=1.41421356d0)
      logical corc(nkey),cod,dsp,normal
      integer*8 latt(nlat)
      dimension b(*)
     1,         imon(nmona,4),emon(nmona,4)
     1,         twiss(nlat,-ndim:ndim,ntwissfun)
     1,         u(4),xx(4)
      include 'inc/common.inc'
c
      cod=corc(1).or.corc(2).or.corc(3).or.corc(4)
      dsp=corc(5).or.corc(6).or.corc(7).or.corc(8)
      if(dsp) then
        esclx=errval(3)/dpshft
        escly=errval(4)/dpshft
      endif
      esclx2=sqr2*esclx
      escly2=sqr2*escly
c     _________
      if( cod ) then
c     :::::::::
        do 100 i=1,nmon
          j=imon(i,2)
          n=imon(j,1)
c          ------ find nearest Quad
          nq=imon(j,4)
c         ____________
          if( normal ) then
c         ::::::::::::       : COD in normal_mode
            r11=twiss(n,ndim,11)
            r12=twiss(n,ndim,12)
            r21=twiss(n,ndim,13)
            r22=twiss(n,ndim,14)
            det=r11*r22-r12*r21
            if(det.gt.1d0) then
              qr=sqrt((det-1d0)/det)
              un=sqrt(det)
              if( corc(1).and.corc(2) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                xp=twiss(n,0,16)-twiss(n,ndim,16)
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                yp=twiss(n,0,18)-twiss(n,ndim,18)
                b(i)=-r12*qr*x + r22*qr*xp + un    *y
                b(i+nmon)=un*x             + r21*qr*y + r22*qr*yp
              elseif( corc(1) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                xp=twiss(n,0,16)-twiss(n,ndim,16)
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                b(i)=-r12*qr*x + r22*qr*xp + un    *y
              elseif( corc(2) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                yp=twiss(n,0,18)-twiss(n,ndim,18)
                b(i) =    un*x             + r21*qr*y + r22*qr*yp
              endif
            else
              un=sqrt(1d0-det)
              if( corc(1).and.corc(2) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                xp=twiss(n,0,16)-twiss(n,ndim,16)
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                yp=twiss(n,0,18)-twiss(n,ndim,18)
                b(i)=      un*x          - r22*y + r12*yp
                b(i+nmon)=r11*x + r12*xp + un *y
              elseif( corc(1) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                yp=twiss(n,0,18)-twiss(n,ndim,18)
                b(i)=      un*x          - r22*y + r12*yp
              elseif( corc(2) ) then
                x=twiss(n,0,15)-twiss(n,ndim,15)
     z           -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
                xp=twiss(n,0,16)-twiss(n,ndim,16)
                y=twiss(n,0,17)-twiss(n,ndim,17)
     z           -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
                b(i) =    r11*x + r12*xp + un *y
              endif
            endif
c         ____
          else
c         ::::               : COD in real_space
            if( corc(3).and.corc(4) ) then
              b(i)=twiss(n,0,15)-twiss(n,ndim,15)
     z            -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
              b(i+nmon)=twiss(n,0,17)-twiss(n,ndim,17)
     z            -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
            elseif( corc(3) ) then
              b(i)=twiss(n,0,15)-twiss(n,ndim,15)
     z            -rlist(latt(nq)+5)-emon(j,1)+errval(2)*tgauss()
            elseif( corc(4) ) then
              b(i)=twiss(n,0,17)-twiss(n,ndim,17)
     z            -rlist(latt(nq)+6)-emon(j,2)+errval(4)*tgauss()
            endif
c         _____
          endif
c         :::::
  100   continue
c     _____
      endif
c     :::::
      nc=0
      do 101 i=1,4
        if(corc(i)) nc=nc+1
  101 continue
      nc=nc*nmon
c     _________________________
      if( corc(5).and.corc(6) ) then
c     ::::::::::::::::::::::: if dispersion x & y in normal mode
        nc1=nc+nmon
        do 110 i=1,nmon
c         ---- measurement of dispersion in real space -----
c         do 111 j=15,18
c           xx(j-14)=(twiss(imon(i),1,j)-twiss(imon(i),-1,j))/ddp
c 111     continue
c         +---Note.----------------------------------------------------+
c         |    Observed   dispersion in real space is defined as the   |
c         |    dispersion in normal mode converted by the inverse of R |
c         |    matrix of real machine. Reconverting this by the R mat- |
c         |    rix of model machine, one obtaines the dispersion in    |
c         |    normal mode inferred from the measured dipersion.       |
c         +------------------------------------------------------------+
          k=imon(i,2)
          do 111 j=7,10
            u(j-6)=twiss(imon(k,1),0,j)
  111     continue
          call mc2to4(twiss,0,imon(k,1),u,xx)
          call mc4to2(twiss,ndim,imon(k,1),xx,u)
          call metaer(twiss,ndim,esclx,escly,i,imon,emon)
          b(nc+i) =u(1)+emon(k,3) -twiss(imon(k,1),ndim,7)
          b(nc1+i)=u(3)+emon(k,4) -twiss(imon(k,1),ndim,9)
  110   continue
c     ________________
      elseif( corc(5) ) then
c     ::::::::::::::::      : if dispersion x in normal mode
        do 112 i=1,nmon
          k=imon(i,2)
          do 113 j=7,10
            u(j-6)=twiss(imon(k,1),0,j)
  113     continue
          call mc2to4(twiss,0,imon(k,1),u,xx)
          call mc4to2(twiss,ndim,imon(k,1),xx,u)
          call metaer(twiss,ndim,esclx,escly,i,imon,emon)
          b(nc+i) =u(1)+emon(k,3) -twiss(imon(k,1),ndim,7)
  112   continue
c     ________________
      elseif( corc(6) ) then
c     ::::::::::::::::      : if dispersion y in normal mode
        do 114 i=1,nmon
          k=imon(i,2)
          do 115 j=7,10
            u(j-6)=twiss(imon(k,1),0,j)
  115     continue
          call mc2to4(twiss,0,imon(k,1),u,xx)
          call mc4to2(twiss,ndim,imon(k,1),xx,u)
          call metaer(twiss,ndim,esclx,escly,i,imon,emon)
          b(nc+i) =u(3)+emon(k,4) -twiss(imon(k,1),ndim,9)
  114   continue
c     _____________________________
      elseif( corc(7).and.corc(8) ) then
c     :::::::::::::::::::::::::::::: if dispersion x & y in real space
        nc1=nc+nmon
        do 120 i=1,nmon
          j=imon(i,2)
          if(simulate) then
            u(1)=twiss(imon(j,1),0,7)
            u(2)=twiss(imon(j,1),0,8)
            u(3)=twiss(imon(j,1),0,9)
            u(4)=twiss(imon(j,1),0,10)
            call mc2to4(twiss,0,imon(j,1),u,xx)
            if(errval(3).ne.0d0) then
              emon(j,3)=tgauss()*esclx2
            else
              emon(j,3)=0d0
            endif
            if(errval(4).ne.0d0) then
              emon(j,4)=tgauss()*escly2
            else
              emon(j,4)=0d0
            endif
            b(nc+i) =xx(1)+emon(j,3)
            b(nc1+i)=xx(3)+emon(j,4)
          else
            b(nc+i) =twiss(imon(j,1),0,7)
            b(nc1+i)=twiss(imon(j,1),0,9)
          endif
          u(1)=twiss(imon(j,1),ndim,7)
          u(2)=twiss(imon(j,1),ndim,8)
          u(3)=twiss(imon(j,1),ndim,9)
          u(4)=twiss(imon(j,1),ndim,10)
          call mc2to4(twiss,ndim,imon(j,1),u,xx)
          b(nc+i) =b(nc+i) -xx(1)
          b(nc1+i)=b(nc1+i)-xx(3)
 120    continue
c     _________________
      elseif( corc(7) ) then
c     :::::::::::::::::     : if dispersion x in real space
        do 130 i=1,nmon
          j=imon(i,2)
          if(simulate) then
            u(1)=twiss(imon(j,1),0,7)
            u(2)=twiss(imon(j,1),0,8)
            u(3)=twiss(imon(j,1),0,9)
            u(4)=twiss(imon(j,1),0,10)
            call mc2to4(twiss,0,imon(j,1),u,xx)
            if(errval(3).ne.0d0) then
              emon(j,3)=tgauss()*esclx2
            else
              emon(j,3)=0d0
            endif
            b(nc+i) =xx(1)+emon(j,3)
          else
            b(nc+i)=twiss(imon(j,1),0,7)
          endif
          u(1)=twiss(imon(j,1),ndim,7)
          u(2)=twiss(imon(j,1),ndim,8)
          u(3)=twiss(imon(j,1),ndim,9)
          u(4)=twiss(imon(j,1),ndim,10)
          call mc2to4(twiss,ndim,imon(j,1),u,xx)
          b(nc+i)=b(nc+i)-xx(1)
 130    continue
c     _________________
      elseif( corc(8) ) then
c     :::::::::::::::::     : if dispersion y in real space
        do 140 i=1,nmon
          j=imon(i,2)
          if(simulate) then
            u(1)=twiss(imon(j,1),0,7)
            u(2)=twiss(imon(j,1),0,8)
            u(3)=twiss(imon(j,1),0,9)
            u(4)=twiss(imon(j,1),0,10)
            call mc2to4(twiss,0,imon(j,1),u,xx)
            if(errval(4).ne.0d0) then
              emon(j,4)=tgauss()*escly2
            else
              emon(j,4)=0d0
            endif
            b(nc+i) =xx(3)+emon(j,4)
          else
            b(nc+i)=twiss(imon(j,1),0,9)
          endif
          u(1)=twiss(imon(j,1),ndim,7)
          u(2)=twiss(imon(j,1),ndim,8)
          u(3)=twiss(imon(j,1),ndim,9)
          u(4)=twiss(imon(j,1),ndim,10)
          call mc2to4(twiss,ndim,imon(j,1),u,xx)
          b(nc+i)=b(nc+i)-xx(3)
 140    continue
c     _____
      endif
c     :::::
      if(corc(5).or.corc(7)) then
        if(corc(7).and.errval(3).ne.0d0) then
          write(lfno,*)
     $         '  -- Statistics of dispersion mesurement-errors --'
          call mstat(emon(1,3),nmon,tmax,tmin,trms,tave)
          tmax=max(tmax,-tmin)
          write(lfno,'(3(A,1PD12.4))') '  dEx (max/rms/mean) =',
     z         tmax,' /',trms,' /',tave
        endif
        if(corc(6).or.(corc(8).and.errval(4).ne.0d0))then
          call mstat(emon(1,4),nmon,tmax,tmin,trms,tave)
          tmax=max(tmax,-tmin)
          write(lfno,'(3(A,1PD12.4))') '  dEy (max/rms/mean) =',
     z         tmax,' /',trms,' /',tave
        endif
      elseif(corc(6).or.(corc(8).and.errval(4).ne.0d0))then
        write(lfno,*)
     $       '  -- Statistics of dispersion mesurement-errors --'
        call mstat(emon(1,4),nmon,tmax,tmin,trms,tave)
        tmax=max(tmax,-tmin)
        write(lfno,'(3(A,1PD12.4))') '  dEy (max/rms/mean) =',
     z                          tmax,' /',trms,' /',tave
      endif
c     write(lfno,*) ' r.h.s vector'
c     write(*,'(1P,10d12.4)') (b(i),i=1,nmon)
      return
      end

      subroutine mcrcod(latt,twiss,mult,imon,emon,nmon,mcod,push,
     1                  print,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical push,print,mcod
      integer*8 latt(*)
      dimension twiss(*),imon(*),emon(*),mult(*)
      real*8 x(nmon*2)
      call mcrcod1(latt,twiss,mult,imon,emon,nmon,
     1             x,mcod,push,print,lfno)
      return
      end

      subroutine mcrcod1(latt,twiss,mult,imon,emon,nmon,temp,mcod,
     1                   push,print,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical push,print,mcod
cHP   character*8 vout(4,2)*80,autofg,name(4,2)
      character*8 vout(4,2)*79,autofg,name(4,2)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),mult(*),
     1          imon(nmona,4),emon(nmona,4),temp(nmon,2),rc(12,2)
      include 'inc/common.inc'
      data rc /24*0d0/
      data vout(3,1)/' max(mm)'/,vout(3,2)/' max(mm)'/,
     1     vout(4,1)/'     at '/,vout(4,2)/'     at '/,
     1     vout(1,1)/' rms(mm)'/,vout(1,2)/' rms(mm)'/,
     1     vout(2,1)/' wghted.'/,vout(2,2)/' wghted.'/
      data name/8*' '/
      save vout,rc,name
c
      if(push) then
        do 10 i=1,2
          vout(3,i)(41:48)=autofg(rc(1,i)*1d3,'S8.5')
          vout(3,i)(49:56)=autofg(rc(4,i)*1d3,'S8.5')
          vout(3,i)(57:64)=autofg(rc(7,i)*1d3,'S8.5')
          vout(3,i)(65:72)=autofg(rc(10,i)*1d3,'S8.5')
          vout(4,i)(42:49)=name(1,i)
          vout(4,i)(50:57)=name(2,i)
          vout(4,i)(58:65)=name(3,i)
          vout(4,i)(66:73)=name(4,i)
          vout(1,i)(41:48)=autofg(rc(2,i)*1d3,'S8.5')
          vout(1,i)(49:56)=autofg(rc(5,i)*1d3,'S8.5')
          vout(1,i)(57:64)=autofg(rc(8,i)*1d3,'S8.5')
          vout(1,i)(65:72)=autofg(rc(11,i)*1d3,'S8.5')
          vout(2,i)(41:48)=autofg(rc(3,i),'S8.5')
          vout(2,i)(49:56)=autofg(rc(6,i),'S8.5')
          vout(2,i)(57:64)=autofg(rc(9,i),'S8.5')
          vout(2,i)(65:72)=autofg(rc(12,i),'S8.5')
   10   continue
      endif
      if(nmon .eq. 0) return
      do 11 i=1,nmon
        j=imon(i,2)
        nq=imon(j,4)
        temp(i,1)=twiss(imon(j,1),0,15)-twiss(imon(j,1),ndim,15)
     1         -rlist(latt(nq)+5)-emon(j,1)
   11 continue
      do 12 i=1,nmon
        j=imon(i,2)
        temp(i,2)=twiss(imon(j,1),0,15)-twiss(imon(j,1),ndim,15)
   12 continue
      do 14 j=1,2
        rc(3,j)=0d0
        do 13 i=1,nmon
          rc(3,j)=rc(3,j) + twiss(imon(imon(i,2),1),0,2)*temp(i,j)**2
   13   continue
        rc(3,j)=sqrt(rc(3,j)/max(nmon,1))
        call mstatp(temp(1,j),nmon,ceil,floor,rc(2,j),xmean,imax)
        rc(1,j)=max(ceil,-floor)*sign(1d0,ceil+floor)
        call elname(latt,imon(imon(imax,2),1),mult,name(1,j))
   14 continue
      do 15 i=1,nmon
        j=imon(i,2)
        nq=imon(j,4)
        temp(i,1)=twiss(imon(j,1),0,17)-twiss(imon(j,1),ndim,17)
     1         -rlist(latt(nq)+6)-emon(j,2)
   15 continue
      do 16 i=1,nmon
        j=imon(i,2)
        temp(i,2)=twiss(imon(j,1),0,17)-twiss(imon(j,1),ndim,17)
   16 continue
      do 18 j=1,2
        rc(6,j)=0d0
        do 17 i=1,nmon
          rc(6,j)=rc(6,j) + twiss(imon(imon(i,2),1),0,5)*temp(i,j)**2
   17   continue
        rc(6,j)=sqrt(rc(6,j)/max(nmon,1))
        call mstatp(temp(1,j),nmon,ceil,floor,rc(5,j),xmean,imax)
        rc(4,j)=max(ceil,-floor)*sign(1d0,ceil+floor)
        call elname(latt,imon(imon(imax,2),1),mult,name(2,j))
   18 continue
      do 20 i=1,nmon
        j=imon(i,2)
   20   temp(i,1)=twiss(imon(j,1),0,7)-twiss(imon(j,1),ndim,7)
      rc(9,1)=0d0
      do 21 i=1,nmon
        j=imon(i,2)
   21   rc(9,1)=rc(9,1) + twiss(imon(j,1),0,2)*temp(i,1)**2
      rc(9,1)=sqrt(rc(9,1)/max(nmon,1))
      call mstatp(temp,nmon,ceil,floor,rc(8,1),xmean,imax)
      rc(7,1)=max(ceil,-floor)*sign(1d0,ceil+floor)
      call elname(latt,imon(imon(imax,2),1),mult,name(3,1))
      do 22 i=1,nmon
        j=imon(i,2)
   22   temp(i,1)=twiss(imon(j,1),0,9)-twiss(imon(j,1),ndim,9)
      rc(12,1)=0d0
      do 23 i=1,nmon
   23   rc(12,1)=rc(12,1) + twiss(imon(imon(i,2),1),0,5)*temp(i,1)**2
      rc(12,1)=sqrt(rc(12,1)/max(nmon,1))
      call mstatp(temp,nmon,ceil,floor,rc(11,1),xmean,imax)
      rc(10,1)=max(ceil,-floor)*sign(1d0,ceil+floor)
      call elname(latt,imon(imon(imax,2),1),mult,name(4,1))
      do 24 i=7,12
        rc(i,2)=rc(i,1)
   24 continue
      name(3,2)=name(3,1)
      name(4,2)=name(4,1)
      if(print) then
        do 25 i=1,2
          vout(3,i)( 9:16)=autofg(rc(1,i)*1d3,'S8.5')
          vout(3,i)(17:24)=autofg(rc(4,i)*1d3,'S8.5')
          vout(3,i)(25:32)=autofg(rc(7,i)*1d3,'S8.5')
          vout(3,i)(33:40)=autofg(rc(10,i)*1d3,'S8.5')
          vout(4,i)(10:17)=name(1,i)
          vout(4,i)(18:25)=name(2,i)
          vout(4,i)(26:33)=name(3,i)
          vout(4,i)(34:41)=name(4,i)
          vout(1,i)( 9:16)=autofg(rc(2,i)*1d3,'S8.5')
          vout(1,i)(17:24)=autofg(rc(5,i)*1d3,'S8.5')
          vout(1,i)(25:32)=autofg(rc(8,i)*1d3,'S8.5')
          vout(1,i)(33:40)=autofg(rc(11,i)*1d3,'S8.5')
          vout(2,i)( 9:16)=autofg(rc(3,i),'S8.5')
          vout(2,i)(17:24)=autofg(rc(6,i),'S8.5')
          vout(2,i)(25:32)=autofg(rc(9,i),'S8.5')
          vout(2,i)(33:40)=autofg(rc(12,i),'S8.5')
   25   continue
        if(mcod) then
          j=1
          write(lfno,
     &         '(t2,a,t9,a,t17,a,t25,a,t33,a,t41,a,t49,a,t57,a,t65,a)')
     1         'MOD    ',' dx     ',' dy     ',' dEX    ',' dEY    ',
     1               ' dx_old ',' dy_old ',' dEX_old',' dEY_old'
        else
          j=2
          write(lfno,
     &         '(t2,a,t9,a,t17,a,t25,a,t33,a,t41,a,t49,a,t57,a,t65,a)')
     1         'TOD    ',' dx     ',' dy     ',' dEX    ',' dEY    ',
     1               ' dx_old ',' dy_old ',' dEX_old',' dEY_old'
        endif
        write(lfno,'(4(a:/))') (vout(i,j),i=1,4)
      endif
      return
      end

      subroutine mcrmat(latt,twiss,gammab,ip,psix,psiy,a,ida,prime,
     $     normal,istr,nstr,imon,nmon,xy)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
      implicit real*8 (a-h,o-z)
      logical cxy(3),prime,normal
      character*(*) xy
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat)
     1,         a(ida,nstr)
     1,         istr(nstra,4),imon(nmona,4)
      include 'inc/common.inc'
c
      cxy(1)=.false.
      cxy(2)=.false.
      cxy(3)=.false.
      if( xy.eq.'X') then
        cxy(1)=.true.
      elseif(xy.eq.'Y') then
        cxy(2)=.true.
      elseif(xy.eq.'XY') then
        cxy(3)=.true.
      endif
      if(psix.eq.0d0) then
        px2=0.5d0*(twiss(nlat,ip,3)-twiss(1,ip,3))
      else
        px2=0.5d0*psix
      endif
      if(psiy.eq.0d0) then
        py2=0.5d0*(twiss(nlat,ip,6)-twiss(1,ip,6))
      else
        py2=0.5d0*psiy
      endif
c     cofx=0.5/sin(px2)
c     cofy=0.5/sin(py2)
      do 300 j=1,nstr
        jp=istr(istr(j,2),1)
        jt=idtypec(jp)
        if(jt.eq.icbend) then
          call mcrmatb(latt,twiss,gammab,ip,px2,py2,a(1,j),ida,prime,
     $     normal,jp,imon,nmon,cxy)
        elseif(jt.eq.icquad) then
          call mcrmatb(latt,twiss,gammab,ip,px2,py2,a(1,j),ida,prime,
     $     normal,jp,imon,nmon,cxy)
        elseif(jt.eq.icsext) then
          call mcrmatb(latt,twiss,gammab,ip,px2,py2,a(1,j),ida,prime,
     $     normal,jp,imon,nmon,cxy)
        endif
  300 continue
      return
      end

      subroutine mcrmatb(latt,twiss,gammab,ip,px2,py2,a,ida,prime,
     1                   normal,jp,imon,nmon,cxy)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      implicit real*8 (a-h,o-z)
c     parameter (nvec=16)
c     logical*4 cxy(3),prime,normal,vec1,vec2,edge
      logical*4 cxy(3),prime,normal,          edge,enter,exit
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat),
     $     a(ida),
     $     imon(nmona,4),
     $     u(4),x(4)
      include 'inc/common.inc'
c
c     cxy(1)=.false.
c     cxy(2)=.false.
c     cxy(3)=.false.
c     if( xy.eq.'X') then
c       cxy(1)=.true.
c     elseif(xy.eq.'Y') then
c       cxy(2)=.true.
c     elseif(xy.eq.'XY') then
c       cxy(3)=.true.
c     endif
      if( .not.cell ) goto 400
c ****** periodic lattice *************************************
c     if(psix.eq.0d0) then
c       px2=0.5d0*(twiss(nlat,ip,3)-twiss(1,ip,3))
c     else
c       px2=0.5d0*psix
c     endif
c     if(psiy.eq.0d0) then
c       py2=0.5d0*(twiss(nlat,ip,6)-twiss(1,ip,6))
c     else
c       py2=0.5d0*psiy
c     endif
      cofx=0.5/sin(px2)
      cofy=0.5/sin(py2)
c     do 300 j=1,nstr
c       jp=istr(j,2)
c       cost=cos(-rlist(idvalc(istr(jp,1))+5))
c       sint=sin(-rlist(idval(idelc(istr(jp,1)))+5))
        rotation=-rlist(idval(idcomp(elatt,jp))+5)
        cost=cos(rotation)
        sint=sin(rotation)
        if(cxy(3)) then
          do 110 k=1,2*nmon
 110        a(k)=0d0
        else
          do 120 k=1,nmon
 120        a(k)=0d0
        endif
        fedge=-1d0
        exit=.false.
        enter=.false.
c      .... calc response matrix at entrance & exit
cm      nslave=master(istr(jp,1))
c       do 130 m=istr(jp,1),istr(jp,1)+1
        do 130 m=jp,jp+1
          fedge=-fedge
          lmo1=nmon
          lmo2=nmon+1
          ient=0
          edge=.false.
          do 132 i=1,nmon
            if(imon(imon(i,2),1).eq.m) then
              lmo1=i-1
              lmo2=i+1
              ient=i
              edge=.true.
              goto 133
            elseif(imon(imon(i,2),1).gt.m) then
              lmo1=i-1
              lmo2=i
              goto 133
            endif
  132     continue
  133     continue
          axj=twiss(m,ip,1)
          sqrbxj=sqrt(twiss(m,ip,2))
          pxj=twiss(m,ip,3)
          ayj=twiss(m,ip,4)
          sqrbyj=sqrt(twiss(m,ip,5))
          pyj=twiss(m,ip,6)
          r11=twiss(m,ip,11)
          r12=twiss(m,ip,12)
          r21=twiss(m,ip,13)
          r22=twiss(m,ip,14)
c     write(*,'(a,1p,4e15.7)') 'rj=',r11,r12,r21,r22
c     write(*,'(a,1p,6e15.7)')'twsj=',axj,twiss(m,ip,2),pxj,ayj,
c    $         twiss(m,ip,5),pyj
          det=r11*r22-r12*r21
          if(det.gt.1d0) then
            qr=sqrt((det-1d0)/det)
            um=sqrt(det)
            tx =  r22*qr*cost
            txp= -r21*qr*cost +    um *sint
            ty =                r22*qr*sint
            typ=  um    *cost + r12*qr*sint
          else
            um=sqrt(1d0-det)
            tx =            r12*sint
            txp=  um *cost -r11*sint
            ty =  r12*cost
            typ=  r22*cost + um*sint
          endif
c     write(*,'(a,1p,4e15.7)') 't=',tx,txp,ty,typ
c         ____________
          if( normal ) then
c         ::::::::::::       : Green's function in normal_mode
            if( cxy(3) ) then
              fx =cofx/sqrbxj *tx
              fxp=cofx*sqrbxj *txp
              fy =cofy/sqrbyj *ty
              fyp=cofy*sqrbyj *typ
              if(.not.prime) then
                do 140 i=1,lmo1
                  n=imon(imon(i,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i)     =a(i)      + sqrbxi *((sx+axj*cx)*fx
     z                 +cx*fxp)
                  a(i+nmon)=a(i+nmon) + sqrbyi *((sy+ayj*cy)*fy
     z                 +cy*fyp)
 140            continue
                do 144 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i)     =a(i)      + sqrbxi *((sx+axj*cx)*fx
     z                 +cx*fxp)
                  a(i+nmon)=a(i+nmon) + sqrbyi *((sy+ayj*cy)*fy
     z                 +cy*fyp)
 144            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-fedge*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-fedge*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(ient)  =a(ient)   + sqrbxi *((sx+axj*cx)*fx
     z                                                +cx*fxp)
                  a(ient+nmon)=a(ient+nmon)+sqrbyi *((sy+ayj*cy)*fy
     z                                                +cy*fyp)
                endif
              else
                do 150 i=1,lmo1
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i)     =a(i)      - fx /sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                 + fxp /sqrbxi * (sx - alfaxi*cx)
                  a(i+nmon)=a(i+nmon) - fy  /sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                 + fyp /sqrbyi*(sy - alfayi*cy)
 150            continue
                do 154 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i)     =a(i)      - fx /sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                 + fxp /sqrbxi * (sx - alfaxi*cx)
                  a(i+nmon)=a(i+nmon) - fy  /sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                 + fyp /sqrbyi*(sy - alfayi*cy)
 154            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sx =-fedge*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  sy =-fedge*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(ient)  =a(ient)    - fx /sqrbxi *
     z                     (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                               + fxp /sqrbxi * (sx - alfaxi*cx)
                  a(ient+nmon)=a(ient+nmon) - fy  /sqrbyi *
     z                     (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                                 + fyp /sqrbyi*(sy - alfayi*cy)
                endif
              endif
            elseif( cxy(1) ) then
              fx =cofx/sqrbxj *tx
              fxp=cofx*sqrbxj *txp
              if(.not.prime) then
                do 160 i=1,lmo1
                  n=imon(imon(i,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(i)= a(i) + sqrbxi *((sx + axj*cx)  *fx
     z                 +cx             *fxp)
 160            continue
                do 164 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(i)= a(i) + sqrbxi *((sx + axj*cx)  *fx
     z                 +cx             *fxp)
 164            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-fedge*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(ient)= a(ient) + sqrbxi *((sx + axj*cx)  *fx
     z                                     +cx             *fxp)
                endif
              else
                do 170 i=1,lmo1
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(i)= a(i) - fx  /sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                 + fxp /sqrbxi * (sx - alfaxi*cx)
 170            continue
                do 174 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(i)= a(i) - fx  /sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                 + fxp /sqrbxi * (sx - alfaxi*cx)
 174            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxi=twiss(n,ip,3)
                  sx =-fedge*sin(px2-abs(pxi-pxj))
                  cx = cos(px2-abs(pxi-pxj))
                  a(ient)= a(ient) - fx  /sqrbxi *
     z                  (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z                       + fxp /sqrbxi * (sx - alfaxi*cx)
                endif
              endif
            elseif( cxy(2) ) then
              fy =cofy/sqrbyj *ty
              fyp=cofy*sqrbyj *typ
              if(.not.prime) then
                do 180 i=1,lmo1
                  n=imon(imon(i,2),1)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i) =a(i) + sqrbyi *(( sy + ayj*cy )*fy
     z                 +cy*fyp)
 180            continue
                do 184 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i) =a(i) + sqrbyi *(( sy + ayj*cy )*fy
     z                 +cy*fyp)
 184            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-fedge*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(ient) =a(ient) + sqrbyi *(( sy + ayj*cy )*fy
     z                                     +cy*fyp)
                endif
              else
                do 190 i=1,lmo1
                  n=imon(imon(i,2),1)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i) = a(i) - fy  /sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                 + fyp /sqrbyi * (sy - alfayi*cy)
 190            continue
                do 194 i=lmo2,nmon
                  n=imon(imon(i,2),1)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(i) = a(i) - fy  /sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                 + fyp /sqrbyi * (sy - alfayi*cy)
 194            continue
                if(edge) then
                  n=imon(imon(ient,2),1)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyi=twiss(n,ip,6)
                  sy =-fedge*sin(py2-abs(pyi-pyj))
                  cy = cos(py2-abs(pyi-pyj))
                  a(ient) = a(ient) - fy  /sqrbyi *
     z                  (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z                        + fyp /sqrbyi * (sy - alfayi*cy)
                endif
              endif
            endif
c         ____
          else
c         ::::               : Green's function in real_space
            fx =cofx/sqrbxj *tx
            fxp=cofx*sqrbxj *txp
            fy =cofy/sqrbyj *ty
            fyp=cofy*sqrbyj *typ
            do 200 i=1,lmo1
              n=imon(imon(i,2),1)
              alfaxi=twiss(n,ip,1)
              sqrbxi=sqrt(twiss(n,ip,2))
              pxi=twiss(n,ip,3)
              alfayi=twiss(n,ip,4)
              sqrbyi=sqrt(twiss(n,ip,5))
              pyi=twiss(n,ip,6)
c     write(*,'(a,1p,6e15.7)')'twsi=',alfaxi,twiss(n,ip,2),pxi,alfayi,
c    $         twiss(n,ip,5),pyi
c             sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
c             cx = cos(px2-abs(pxi-pxj))
              sx =-sin(px2+pxi-pxj )
              cx = cos(px2+pxi-pxj )
c             sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
c             cy = cos(py2-abs(pyi-pyj))
              sy =-sin(py2+pyi-pyj )
              cy = cos(py2+pyi-pyj )
              u(1)= sqrbxi *( (sx+axj*cx)*fx +cx*fxp )
              u(2)= - fx/sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z            + fxp/sqrbxi * (sx - alfaxi*cx)
              u(3)= sqrbyi *( (sy+ayj*cy)*fy +cy*fyp )
c Oide modified following one line. 4/28/1993
c  fy had been fx.   There were three lines more with the same change.
              u(4)= - fy/sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z            + fyp /sqrbyi*(sy - alfayi*cy)
c     write(*,'(a,1p,4e15.7)') 'xnormal=',u(1),u(2),u(3),u(4)

              call mc2to4(twiss,ip,n,u,x)
c     write(*,'(a,1p,4e15.7)') 'xreal=',x(1),x(2),x(3),x(4)
              if( cxy(3) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(1)
                  a(i+nmon) = a(i+nmon) + x(3)
                else
                  a(i)      = a(i) + x(2)
                  a(i+nmon) = a(i+nmon) + x(4)
                endif
              elseif( cxy(1) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(1)
                else
                  a(i)      = a(i) + x(2)
                endif
c     write(*,'(a,1p,4e15.7)') 'response=',x(1)
              elseif( cxy(2) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(3)
                else
                  a(i)      = a(i) + x(4)
                endif
              endif
  200       continue
            do 202 i=lmo2,nmon
              n=imon(imon(i,2),1)
              alfaxi=twiss(n,ip,1)
              sqrbxi=sqrt(twiss(n,ip,2))
              pxi=twiss(n,ip,3)
              alfayi=twiss(n,ip,4)
              sqrbyi=sqrt(twiss(n,ip,5))
              pyi=twiss(n,ip,6)
c             sx =-sign(1d0,pxj-pxi)*sin(px2-abs(pxi-pxj))
c             cx = cos(px2-abs(pxi-pxj))
              sx = sin(px2-pxi+pxj )
              cx = cos(px2-pxi+pxj )
c             sy =-sign(1d0,pyj-pyi)*sin(py2-abs(pyi-pyj))
c             cy = cos(py2-abs(pyi-pyj))
              sy = sin(py2-pyi+pyj )
              cy = cos(py2-pyi+pyj )
              u(1)= sqrbxi *( (sx+axj*cx)*fx +cx*fxp )
              u(2)= - fx/sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z            + fxp/sqrbxi * (sx - alfaxi*cx)
              u(3)= sqrbyi *( (sy+ayj*cy)*fy +cy*fyp )
              u(4)= - fy/sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z            + fyp /sqrbyi*(sy - alfayi*cy)
              call mc2to4(twiss,ip,n,u,x)
              if( cxy(3) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(1)
                  a(i+nmon) = a(i+nmon) + x(3)
                else
                  a(i)      = a(i) + x(2)
                  a(i+nmon) = a(i+nmon) + x(4)
                endif
              elseif( cxy(1) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(1)
                else
                  a(i)      = a(i) + x(2)
                endif
              elseif( cxy(2) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(3)
                else
                  a(i)      = a(i) + x(4)
                endif
              endif
  202       continue
            if(edge) then
              n=imon(imon(ient,2),1)
              alfaxi=twiss(n,ip,1)
              sqrbxi=sqrt(twiss(n,ip,2))
              pxi=twiss(n,ip,3)
              alfayi=twiss(n,ip,4)
              sqrbyi=sqrt(twiss(n,ip,5))
              pyi=twiss(n,ip,6)
c             sx =-fedge*sin(px2-abs(pxi-pxj))
c             cx = cos(px2-abs(pxi-pxj))
              sx =-fedge*sin(px2)
              cx = cos(px2)
c             sy =-fedge*sin(py2-abs(pyi-pyj))
c             cy = cos(py2-abs(pyi-pyj))
              sy =-fedge*sin(py2)
              cy = cos(py2)
              u(1)= sqrbxi *( (sx+axj*cx)*fx +cx*fxp )
              u(2)= - fx/sqrbxi *
     z                 (  (alfaxi-axj)*sx +(1d0+alfaxi*axj)*cx  )
     z            + fxp/sqrbxi * (sx - alfaxi*cx)
              u(3)= sqrbyi *( (sy+ayj*cy)*fy +cy*fyp )
              u(4)= - fy/sqrbyi *
     z                 (  (alfayi-ayj)*sy +(1d0+alfayi*ayj)*cy  )
     z            + fyp /sqrbyi*(sy - alfayi*cy)
              call mc2to4(twiss,ip,n,u,x)
              if( cxy(3) ) then
                if(.not.prime) then
                  a(ient)      = a(ient) + x(1)
                  a(ient+nmon) = a(ient+nmon) + x(3)
                else
                  a(ient)      = a(ient) + x(2)
                  a(ient+nmon) = a(ient+nmon) + x(4)
                endif
              elseif( cxy(1) ) then
                if(.not.prime) then
                  a(ient)      = a(ient) + x(1)
                else
                  a(ient)      = a(ient) + x(2)
                endif
              elseif( cxy(2) ) then
                if(.not.prime) then
                  a(ient)      = a(ient) + x(3)
                else
                  a(ient)      = a(ient) + x(4)
                endif
              endif
            endif
c         _____
          endif
c         -----
  130   continue
c       ..... mean value of entrance & exit.
        if( cxy(3) ) then
          if(trpt) then
            do 210 i=1,nmon
              k=imon(imon(i,2),1)
              g=sqrt(gammab(m)/gammab(k))*0.5d0
              a(i)=a(i)*g
 210          a(i+nmon)=a(i+nmon)*g
          else
            do 211 i=1,2*nmon
 211          a(i)=a(i)*0.5d0
          endif
        else
          if(trpt) then
            do 220 i=1,nmon
              k=imon(imon(i,2),1)
 220          a(i)=a(i)*0.5d0*sqrt(gammab(m)/gammab(k))
          else
            do 221 i=1,nmon
 221          a(i)=a(i)*0.5d0
          endif
        endif
c 300 continue
      return
c
c ****** non-periodic lattice *********************************
  400 continue
      idwn=1
c     do 410 j=1,nstr
c       jp=istr(j,2)
c       m=istr(jp,1)
c       cost=cos(-rlist(idval(idelc(m))+5))
c       sint=sin(-rlist(idval(idelc(m))+5))
        rotation=-rlist(idval(idcomp(elatt,jp))+5)
        cost=cos(rotation)
        sint=sin(rotation)
        if(cxy(3)) then
          do 420 k=1,nmon
            a(k)=0d0
            a(k+nmon)=0d0
 420      continue
        else
          do 430 k=1,nmon
            a(k)=0d0
 430      continue
        endif
        do 440 i=idwn,nmon
          if(imon(imon(i,2),1).gt.jp) then
            idwn=i
            goto 442
          endif
  440   continue
        idwn=nmon+1
  442   continue
c      .... calc response matrix at entrance & exit
c       ms=istr(jp,1)
c       if(nslave.eq.0) then
c         mf=ms+1
c       else
c         do 444 i=ms+1,nlat-1
c           if(master(i).eq.-ms) then
c             mf=i+1
c             goto 445
c           endif
c444      continue
c       endif
c445    continue
c       do 450 m=ms,mf,mf-ms
        do 450 m=jp,jp+1
          axj=twiss(m,ip,1)
          sqrbxj=sqrt(twiss(m,ip,2))
          pxj=twiss(m,ip,3)
          ayj=twiss(m,ip,4)
          sqrbyj=sqrt(twiss(m,ip,5))
          pyj=twiss(m,ip,6)
          r11=twiss(m,ip,11)
          r12=twiss(m,ip,12)
          r21=twiss(m,ip,13)
          r22=twiss(m,ip,14)
          det=r11*r22-r12*r21
          if(det.gt.1d0) then
            qr=sqrt((det-1d0)/det)
            um=sqrt(det)
            tx =  r22*qr*cost
            txp= -r21*qr*cost +    um *sint
            ty =                r22*qr*sint
            typ=  um    *cost + r12*qr*sint
          else
            um=sqrt(1d0-det)
            tx =            r12*sint
            txp=  um *cost -r11*sint
            ty =  r12*cost
            typ=  r22*cost + um*sint
          endif
c         ____________
          if( normal ) then
c         ------------       : Green's function in normal_mode
            fx = tx /sqrbxj
            fxp=sqrbxj *txp
            fy = ty /sqrbyj
            fyp=sqrbyj *typ
            if( cxy(3) ) then
              if(.not.prime) then
                do 460 i=idwn,nmon
                  n=imon(imon(i,2),1)
c                 alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxij=twiss(n,ip,3)-pxj
c                 alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyij=twiss(n,ip,6)-pyj
                  sx = sin(pxij)
                  cx = cos(pxij)
                  sy = sin(pyij)
                  cy = cos(pyij)
                  a(i)=a(i) + sqrbxi*(
     z                fx * (cx+axj*sx) + fxp * sx )
                  a(i+nmon)=a(i+nmon) + sqrbyi*(
     z                fy * (cy+ayj*sy) + fyp * sy )
  460           continue
              else
                do 470 i=idwn,nmon
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxij=twiss(n,ip,3)-pxj
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyij=twiss(n,ip,6)-pyj
                  sx = sin(pxij)
                  cx = cos(pxij)
                  sy = sin(pyij)
                  cy = cos(pyij)
                  a(i)=a(i) +(-fx *
     z               (  (1d0+alfaxi*axj)*sx +(alfaxi-axj)*cx  )
     z               + fxp *(cx-alfaxi*sx) )/sqrbxi
                  a(i+nmon)=a(i+nmon) +(-fy *
     z               (  (1d0+alfayi*ayj)*sy +(alfayi-ayj)*cy  )
     z               + fyp *(cy-alfayi*sy) )/sqrbyi
  470           continue
              endif
            elseif( cxy(1) ) then
              if(.not.prime) then
                do 480 i=idwn,nmon
                  n=imon(imon(i,2),1)
c                 alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxij=twiss(n,ip,3)-pxj
                  sx = sin(pxij)
                  cx = cos(pxij)
                  a(i)=a(i) + sqrbxi*(
     z                 fx* (cx+axj*sx) + fxp * sx )
 480            continue
              else
                do 490 i=idwn,nmon
                  n=imon(imon(i,2),1)
                  alfaxi=twiss(n,ip,1)
                  sqrbxi=sqrt(twiss(n,ip,2))
                  pxij=twiss(n,ip,3)-pxj
                  sx = sin(pxij)
                  cx = cos(pxij)
                  a(i)=a(i) +(-fx *
     z                 (  (1d0+alfaxi*axj)*sx +(alfaxi-axj)*cx  )
     z                 + fxp *(cx-alfaxi*sx) )/sqrbxi
 490            continue
              endif
            elseif( cxy(2) ) then
              if(.not.prime) then
                do 500 i=idwn,nmon
                  n=imon(imon(i,2),1)
c                 alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyij=twiss(n,ip,6)-pyj
                  sy = sin(pyij)
                  cy = cos(pyij)
                  a(i)=a(i) + sqrbyi*(
     z                 fy * (cy+ayj*sy) + fyp * sy )
 500            continue
              else
                do 510 i=idwn,nmon
                  n=imon(imon(i,2),1)
                  alfayi=twiss(n,ip,4)
                  sqrbyi=sqrt(twiss(n,ip,5))
                  pyij=twiss(n,ip,6)-pyj
                  sy = sin(pyij)
                  cy = cos(pyij)
                  a(i)=a(i) +(-fy *
     z                 (  (1d0+alfayi*ayj)*sy +(alfayi-ayj)*cy  )
     z                 + fyp *(cy-alfayi*sy) )/sqrbyi
 510            continue
              endif
            endif
c         ____
          else
c         ::::               : Green's function in real_space
            fx = tx /sqrbxj
            fxp=sqrbxj *txp
            fy = ty /sqrbyj
            fyp=sqrbyj *typ
            do 520 i=idwn,nmon
              n=imon(imon(i,2),1)
              alfaxi=twiss(n,ip,1)
              sqrbxi=sqrt(twiss(n,ip,2))
              pxij=twiss(n,ip,3)-pxj
              alfayi=twiss(n,ip,4)
              sqrbyi=sqrt(twiss(n,ip,5))
              pyij=twiss(n,ip,6)-pyj
              sx = sin(pxij)
              cx = cos(pxij)
              sy = sin(pyij)
              cy = cos(pyij)
c             u(1)= fx*sqrbxi * (cx+axj*sx) + fxp*sqrbxi * sx
              u(1)= sqrbxi* (fx * (cx+axj*sx) + fxp * sx)
c             u(2)=-fx/sqrbxi *
c    z               (  (1d0+alfaxi*axj)*sx +(alfaxi-axj)*cx  )
c    z           + fxp/sqrbxi *(cx-alfaxi*sx)
              u(2)=(-fx *
     z               (  (1d0+alfaxi*axj)*sx +(alfaxi-axj)*cx  )
     z           + fxp *(cx-alfaxi*sx) )/sqrbxi
c             u(3)= fy*sqrbyi * (cy+ayj*sy) + fyp*sqrbyi * sy
              u(3)= sqrbyi* (fy * (cy+ayj*sy) + fyp * sy )
c             u(4)=-fy/sqrbyi *
c    z               (  (1d0+alfayi*ayj)*sy +(alfayi-ayj)*cy  )
c    z           + fyp/sqrbyi *(cy-alfayi*sy)
              u(4)=(-fy *
     z               (  (1d0+alfayi*ayj)*sy +(alfayi-ayj)*cy  )
     z           + fyp *(cy-alfayi*sy) )/sqrbyi
              call mc2to4(twiss,ip,n,u,x)
              if( cxy(3) ) then
                if(.not.prime) then
                  a(i)      = a(i) + x(1)
                  a(i+nmon) = a(i+nmon) + x(3)
                else
                  a(i)      = a(i) + x(2)
                  a(i+nmon) = a(i+nmon) + x(4)
                endif
              elseif( cxy(1) ) then
                if(.not.prime) then
                  a(i) = a(i) + x(1)
                else
                  a(i) = a(i) + x(2)
                endif
              elseif( cxy(2) ) then
                if(.not.prime) then
                  a(i) = a(i) + x(3)
                else
                  a(i) = a(i) + x(4)
                endif
              endif
  520       continue
c         _____
          endif
c         :::::
  450   continue
c       ..... mean value of entrance & exit.
        if( cxy(3) ) then
          do 530 i=1,nmon
            k=imon(imon(i,2),1)
            g=sqrt(gammab(m)/gammab(k))*0.5d0
            a(i)=a(i)*g
 530        a(i+nmon)=a(i+nmon)*g
        else
          do 540 i=1,nmon
            k=imon(imon(i,2),1)
 540        a(i)=a(i)*0.5d0*sqrt(gammab(m)/gammab(k))
        endif
c 410  continue
      return
      end

      subroutine  mcrmda(twiss,a,ap,d,ida,
     1                   nstr,imon,nmon,prime,cor)
c     ---- calc  dS/dp * S^(-1) * u(0) ----
c          a(i,j) = u ;  a(i+nmon,j) = v
c          ap(i,j)= up;  ap(i+nmon,j)= vp
      use tffitcode
      use ffs
      implicit real*8 (a-h,o-z)
      logical prime,cor(2)
      dimension a(2*nmon,nstr),ap(2*nmon,nstr),d(ida,nstr)
     z,         imon(nmona,4)
     z,         twiss(nlat,-ndim:ndim,ntwissfun)
     z,         p(4,4,2),s(4,4),sp(4,4),sum(4)
      equivalence (p(1,1,1),sp(1,1)),(p(1,1,2),s(1,1))
      logical errflg
      common /codcor/eptsol,errflg(4),errval(4),dpshft,optiv(18),
     1               ibckup,nmona,nmonact,nstra,nstract,
     1               itmon,itemon,itstr,itestr
c          mcepst msolvg mcmon mcnrmc mdpmax palgn pcbak pundo prkick
c eptsol     o      o                          o
c errflg                   o
c errval                   o     o             o
c dpshft                         o      o
c ibckup                                             o      o     o
c optiv                                              o
c nmona                    o     o             o
c nmonact
c nstra                                        o                  o
c nstract
c itmon
c itemon
c itstr
c itestr
c          pbump  corinit mcstr mccor mcrcod mcrmda monel pwrite mrecal
c eptsol     o      o             o
c errflg            o
c errval            o
c dpshft            o
c ibckup     o      o
c optiv
c nmona      o                    o     o      o      o     o      o
c nmonact
c nstra                     o     o                         o      o
c nstract
c itmon             o
c itemon            o
c itstr             o       o
c itestr            o       o
c          pkill pstati mstore ptrim petcod mclear twsdrw pbumps monact
c eptsol
c errflg
c errval
c dpshft
c ibckup
c optiv
c nmona      o     o      o      o      o            o
c nmonact                                                          o
c nstra      o     o      o                   o              o
c nstract
c itmon                                                            o
c itemon
c itstr
c itestr
c          mcrmat mbmp msolb pcrmat ptol pcset mstack mweght
c eptsol
c errflg
c errval
c dpshft
c ibckup
c optiv
c nmona                                          o      o
c nmonact
c nstra      o     o     o     o     o     o     o
c nstract
c itmon
c itemon
c itstr
c itestr
      do 100 i=1,nmon
        n=imon(imon(i,2),1)
        do 120 ip=ndim-1,1-ndim,2*(1-ndim)
          if(ip.eq.ndim-1) then
            iq=2
          else
            iq=1
          endif
          r11=twiss(n,ip,11)
          r12=twiss(n,ip,12)
          r21=twiss(n,ip,13)
          r22=twiss(n,ip,14)
          det=r11*r22-r12*r21
          if(det.gt.1d0) then
            qr=sqrt((det-1d0)/det)
            un=sqrt(det)
            p(1,1,iq)= -r12 *qr
            p(1,2,iq)=  r22 *qr
            p(2,1,iq)=  r11 *qr
            p(2,2,iq)= -r21 *qr
            p(1,3,iq)=  un
            p(1,4,iq)=  0d0
            p(2,3,iq)=  0d0
            p(2,4,iq)=  un
            p(3,1,iq)=  un
            p(3,2,iq)=  0d0
            p(4,1,iq)=  0d0
            p(4,2,iq)=  un
            p(3,3,iq)=  r21 *qr
            p(3,4,iq)=  r22 *qr
            p(4,3,iq)=  r11 *qr
            p(4,4,iq)=  r12 *qr
          else
            un=sqrt(1d0-det)
            p(1,1,iq)=  un
            p(1,2,iq)=  0d0
            p(2,1,iq)=  0d0
            p(2,2,iq)=  un
            p(1,3,iq)= -r22
            p(1,4,iq)=  r12
            p(2,3,iq)=  r21
            p(2,4,iq)= -r11
            p(3,1,iq)=  r11
            p(3,2,iq)=  r12
            p(4,1,iq)=  r21
            p(4,2,iq)=  r22
            p(3,3,iq)=  un
            p(3,4,iq)=  0d0
            p(4,3,iq)=  0d0
            p(4,4,iq)=  un
          endif
  120   continue
        do 122 k=1,4
          do 122 l=1,4
            sp(l,k)=s(l,k)-sp(l,k)
  122   continue
        r11=twiss(n,ndim,11)
        r12=twiss(n,ndim,12)
        r21=twiss(n,ndim,13)
        r22=twiss(n,ndim,14)
        det=r11*r22-r12*r21
        if(det.gt.1d0) then
          qr=sqrt((det-1d0)/det)
          un=sqrt(det)
          s(1,1)= -r21 *qr
          s(1,2)= -r22 *qr
          s(2,1)= -r11 *qr
          s(2,2)= -r12 *qr
          s(1,3)=  un
          s(1,4)=  0d0
          s(2,3)=  0d0
          s(2,4)=  un
          s(3,1)=  un
          s(3,2)=  0d0
          s(4,1)=  0d0
          s(4,2)=  un
          s(3,3)=  r12 *qr
          s(3,4)= -r22 *qr
          s(4,3)= -r11 *qr
          s(4,4)=  r21 *qr
        else
          un=sqrt(1d0-det)
          s(1,1)=  un
          s(1,2)=  0d0
          s(2,1)=  0d0
          s(2,2)=  un
          s(1,3)=  r22
          s(1,4)= -r12
          s(2,3)= -r21
          s(2,4)=  r11
          s(3,1)= -r11
          s(3,2)= -r12
          s(4,1)= -r21
          s(4,2)= -r22
          s(3,3)=  un
          s(3,4)=  0d0
          s(4,3)=  0d0
          s(4,4)=  un
        endif
        if(prime) then
          do 131 j=2,4,2
            j1=j-1
            do 130 k=1,4
              sum(k)=0d0
              sum(k)=sum(k) + sp(j,1)*s(1,k)
     1                      + sp(j,2)*s(2,k)
     1                      + sp(j,3)*s(3,k)
     1                      + sp(j,4)*s(4,k)
  130       continue
            sp(j1,1)=sum(1)
            sp(j1,2)=sum(2)
            sp(j1,3)=sum(3)
            sp(j1,4)=sum(4)
  131     continue
        else
          do 133 j=1,4,2
            do 132 k=1,4
              sum(k)=0d0
              sum(k)=sum(k) + sp(j,1)*s(1,k)
     1                      + sp(j,2)*s(2,k)
     1                      + sp(j,3)*s(3,k)
     1                      + sp(j,4)*s(4,k)
  132       continue
            sp(j,1)=sum(1)
            sp(j,2)=sum(2)
            sp(j,3)=sum(3)
            sp(j,4)=sum(4)
  133     continue
        endif
        in=i+nmon
        if( cor(1).and.cor(2) ) then
          do 200 j=1,nstr
            d(i,j)=sp(1,1)*a(i,j) + sp(1,2)*ap(i,j)
     1              + sp(1,3)*a(in,j) + sp(1,4)*ap(in,j)
            d(in,j)=sp(3,1)*a(i,j) + sp(3,2)*ap(i,j)
     1              + sp(3,3)*a(in,j) + sp(3,4)*ap(in,j)
  200     continue
        elseif( cor(1) ) then
          do 210 j=1,nstr
            d(i,j)=sp(1,1)*a(i,j) + sp(1,2)*ap(i,j)
     1              + sp(1,3)*a(in,j) + sp(1,4)*ap(in,j)
  210     continue
        elseif( cor(2) ) then
          do 220 j=1,nstr
            d(i,j)=sp(3,1)*a(i,j) + sp(3,2)*ap(i,j)
     1              + sp(3,3)*a(in,j) + sp(3,4)*ap(in,j)
  220     continue
        endif
  100 continue
      return
      end

      subroutine mcstr(word,latt,mult,master,ipstr,ipestr,nstr,initial,
     &                 lfno)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc
      logical initial
      character*(*) word
      integer*8 latt(nlat),ipstr,ipestr,itstr,itestr
      integer*4 mult(*),master(nlat)
      include 'inc/common.inc'
      include 'inc/coroper.inc'
c
      if(itstr.eq.0) then
        ns=0
        do 20 i=1,nlat-1
          if(idvalc(i).eq.icbend) then
            if(master(i).gt.0) then
              ns=ns+1
            endif
          endif
   20    continue
        if(ns.eq.0) return
        ipstr =ktaloc(2*ns)
        ipestr=ktaloc(ns)
        itstr=ipstr
        itestr=ipestr
        nstra=ns
        ns=0
        do 21 i=1,nlat-1
          if(idvalc(i).eq.icbend) then
            if(master(i).gt.0) then
              ilist(mod(ns,2)+1,ipstr+int(ns/2))=i
              ilist(mod(ns+nstra,2)+1,ipstr+int((nstra+ns)/2))=ns+1
              ns=ns+1
            endif
          endif
   21    continue
        nstr=0
        call pclr(rlist(ipestr),ns)
      endif
      if(istope.ne.0) then
        call tfree(istope)
        istope=0
        nstope=0
      endif
      call mcoupste('REMOVE',latt,mult,rlist(itstr),nstr,lfno)
      if(initial) return
c     ... mcstr is called by preadstr with initial=.true. .....
      call mcstr1(word,latt,mult,rlist(itstr),rlist(itestr),nstr,lfno)
      return
      end
c      
      subroutine mcoupste(word,latt,mult,istr,nstr,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical exist
      character*(*) word
      character*16 en
      character*11 autofg
      integer*8 latt(2,nlat),mcaloc,it,it1,last,iofs,next,l,
     $     itemp,ia
      integer*4 mult(*)
      dimension istr(nstra,4)
      include 'inc/common.inc'
      integer*8 itcoupste
      common /corcoup/ itcoupste
c
      if(word.eq.'REMOVE') then
c.......remove
        if(itcoupste.ne.0) then
          next=0
          last=itcoupste
          do while(next.ne.itcoupste)
            next=klist(last+1)
c           print *,last
            call tfree(last)
            last=next
          enddo
          itcoupste=0
        endif
        return
      endif
      itemp=ktaloc(nstr)
      ia=ktaloc(nstr)
      js=0
 1    a=getva(exist)
      if(.not.exist) a=1d0
      call getwdl(word)
      if(word.eq.' ') then
        if(js.eq.0) then
c.........write accumulated couple condition
          if(itcoupste.ne.0) then
            next=0
            last=itcoupste
            do while(next.ne.itcoupste)
              next=klist(last+1)
c.............print out
              n=ilist(1,last+1)
              l=last+3+n
              do i=1,n
                a=rlist(last+2+i)
                if(a.ne.1d0) then
                  en=autofg(a,'S11.8')
                  call mbufw(en,.false.,lfno)
                endif
                k=ilist(mod(i-1,2)+1,l+(i-1)/2)
                call elname(latt,k,mult,en)
                call mbufw(en,.false.,lfno)
              enddo
              a=rlist(last+2+n+1)
              en=autofg(a,'S11.8')
              call mbufw(en,.false.,lfno)
              call mbufw(' ',.true.,lfno)
              last=next
            enddo
          endif
        else
          it=mcaloc(itcoupste,js+3+(js+1)/2)+1
          it1=it+1
          ilist(1,it1)=js
          call tmov(rlist(ia),rlist(it1+1),js)
          rlist(it1+1+js)=a
          call pmovi(rlist(itemp),rlist(it1+1+js+1),js)
        endif
        go to 99
      endif
      if(word(lene(word)-1:).eq.'.1')word=word(1:lene(word)-2)
      do i=1,nstr
        l=istr(istr(i,2),1)
        call elname(latt,int(l),mult,en)
        if(word.eq.en) then
          js=js+1
          ilist(mod(js-1,2)+1,itemp+(js-1)/2)=l
          rlist(ia-1+js)=a
          go to 1
        endif
      enddo
 99   call tfree(itemp)
      call tfree(ia)
      return
      end
c
      integer*8 function mcaloc(itroot,n)
      use tfstk
      use ffs
      use tffitcode
      integer*8 next,itroot,last
      if(itroot.eq.0) then
        itroot=ktaloc(n+1)
        klist(itroot)=itroot
        klist(itroot+1)=itroot
        mcaloc=itroot
      else
        next=0
        last=itroot
        do while(next.ne.itroot)
          next=klist(last+1)
c         print *,last,next
          last=next
        enddo
        last=klist(itroot)
        next=ktaloc(n+1)
        klist(last+1)=next
        klist(next)=last
        klist(next+1)=itroot
        klist(itroot)=next
        mcaloc=next
      endif
c     print *,mcaloc
      return
      end
c
      integer*4 function mcoupsten()
      use tfstk
      use ffs
      use tffitcode
      common /corcoup/ itcoupste
      if(itcoupste.eq.0) then
        mcoupsten=0
      else
        n=0
        next=0
        last=itcoupste
        do while(next.ne.itcoupste)
          n=n+1
          next=ilist(2,last)
          last=next
        enddo
        mcoupsten=n
      endif
      return
      end
c
      subroutine mcoupstea(ar,nc,latt,mult,istr,nstr,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      integer*8 latt(nlat)
      integer*4 mult(*)
      dimension istr(nstra,4)
      dimension ar(nstr+1,nc)
      include 'inc/common.inc'
      common /corcoup/ itcoupste
c
      ar=0.d0
c      call pclr(ar,(nstr+1)*nc)
      icoup=0
      next=0
      last=itcoupste
      do while(next.ne.itcoupste)
        icoup=icoup+1
        next=ilist(2,last)
c.......make array     
        n=ilist(1,last+1)
        l=last+3+n
        cr=rlist(last+1+n+1)
c       print *,n,cr
        do i=1,n
          a=rlist(last+1+i)
          j=ilist(mod(i-1,2)+1,l+(i-1)/2)
c         print *,j,a
          do k=1,nstr
            if(j.eq.istr(istr(k,2),1)) then
              ar(k,icoup)=a
              go to 10
            endif
          enddo
          cr=cr-a*rlist(latt(j)+11)
 10       continue
        enddo
        ar(nstr+1,icoup)=cr
        last=next
      enddo
c     do i=1,nc
c       write(lfno,'(1p,10d11.3)')(ar(j,i),j=1,nstr+1)
c     enddo
      return
      end

      subroutine mcstr1(word,latt,mult,istr,estr,nstr,lfno)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      implicit real*8 (a-h,o-z)
      character*(*) word
      character*16 en
      logical tmatch, exist, err
      integer*8 latt(nlat)
      dimension mult(*)
      dimension istr(nstra,4), estr(nstra)
      external pack
      include 'inc/common.inc'
      errst=getva(exist)
      if(exist) then
        err=.true.
      else
        err=.false.
        errst=0d0
      endif
      istr(:,4)=istr(:,3)
c      call pmovi(istr(1,3),istr(1,4),nstra)
      do 9 i=1,nstra
        istr(i,3)=1
    9 continue
      na=0
      ns=0
    1 continue
      call getwdl(word)
      do 10 i=1,nstra
        if( tmatch(pname(idcomp(elatt,istr(istr(i,2),1))),word) )then
          istr(istr(i,2),3)=0
          ns=ns+1
        endif
   10 continue
      if(ns.ne.na) then
        na=ns
        goto 1
      else
        do 11 i=1,nstra
          call elnameK(istr(istr(i,2),1),en)
          if(en.eq.word) then
            istr(istr(i,2),3)=0
            ns=ns+1
            na=ns
            goto 1
          endif
 11     continue
      endif
      if(ns.eq.0) then
        istr(:,3)=istr(:,4)
c        call pmovi(istr(1,4),istr(1,3),nstra)
        return
      endif
      call pack(istr(1,2),istr(1,3),nstr,nstra)
c
      if( err ) then
        do 30 i=1,nstr
          j=istr(i,2)
          estr(j)=errst * tgauss()
   30   continue
      else
        do 31 i=1,nstr
          j=istr(i,2)
          estr(j)=0d0
   31   continue
      endif
      nstract=nstr
      write(lfno,'(2(A,I4))')
     1  '  Correction dipole available :',nstr,' in ',nstra
      return
      end
