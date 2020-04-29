      subroutine palgn(word,wordp,latt,twiss,gammab,mult,pairs,npair,
     $     istr,jst,estr,nstr,imon,nmo,cor,lfno)
      use ffs
      use tffitcode
      use tfstk
      use ffs_pointer, only:idelc,idtypec
      implicit real*8(a-h,o-z)
      logical cor(*),exist,bump,pvert
      character*(*) word,wordp
      integer*8 latt(nlat)
      dimension twiss(*),gammab(nlat)
     z,         mult(*),pairs(*)
     z,         istr(nstra,4),jst(*),estr(*),imon(*)
      external pack
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
      npair1=npair
      npnti=ielm(wordp,exist)
      if(exist) then
        call getwdl2(word,wordp)
        npntf=ielm(wordp,exist)
        if(exist) then
          bump=.true.
        else
          bump=.false.
        endif
      else
        return
      endif
      xamp=getva(exist)
      if(.not.exist) then
        xamp=5d-3
      endif
      call getwdl2(word,wordp)
      if(.not.bump)then
        npntk=npnti
        npnti=1
      endif
      nquad=0
      do 20 i=npnti,npntf-1
        if(idtypec(i).eq.4) then
          do 21 j=i-1,-nlat,-1
            k=j+nlat*((1-sign(1,j-1))/2)
            if(idtypec(k).ne.41.and.k.ne.nlat)then
              if(k.le.npntf-1.and.idtypec(k).eq.4) then
              else
                nquad=nquad+1
              endif
              goto 20
            endif
   21     continue
        endif
   20 continue
      print *,'nquad=',nquad,' npnti=',npnti,' npntf=',npntf
      if(.not.bump)then
        npnti=npntk
      endif
c ----- save current steers. -----
      call pcbak(latt,twiss)
c --------------------------------
      if( bump ) then
        call ppair1(latt,npnti,npntf,6,npair)
        print *,' no. of pairs (specified) =',npair
        nfc=npair+4
      else
        nfc=1
      endif
c     ..... exclude vertical steerings
      do 10 i=1,nstr
        j=istr(i,2)
        if(pvert(latt,istr(j,1))) then
          istr(j,3)=istr(j,3)+1
        endif
   10 continue
      call pack(istr(1,2),istr(1,3),ns,nstr)
      ikf=italoc((nfc+1)/2)
      ifp=italoc((nfc+1)/2)
      imp=italoc((nfc+1)/2)
      ifv=italoc(nfc)
      ipa=italoc(npair)
      iqu=italoc((nquad+1)/2+1)
      call palgn1(latt,twiss,mult,pairs,rlist(ipa),npair,
     1            rlist(iqu),ist,jst,estr,ns,imon,nmo,cor,
     1            npnti,npntf,bump,nfc,rlist(ikf),rlist(ifp),rlist(imp),
     1            rlist(ifv),xamp,lfno)
      npair=npair1
c     ..... include vertical steerings
      do 11 i=1,nstr
        j=istr(i,2)
        if(pvert(latt,istr(j,1))) then
          istr(j,3)=min(0,istr(j,3)-1)
        endif
   11 continue
      call pack(istr(1,2),istr(1,3),nstr,nstr)
      call tfree(int8(iqu))
      call tfree(int8(ipa))
      call tfree(int8(ifv))
      call tfree(int8(imp))
      call tfree(int8(ifp))
      call tfree(int8(ikf))
      return
      end

      subroutine palgn1(latt,twiss,mult,ipair,ipairt,npair,iquad,
     1                  istr,jst,estr,nstr,imon,nmo,cor,npnti,
     1                  npntf,bump,nfc,kfit,ifitp,mfitp,fitval,xamp,
     1                  lfno)
      use ffs
      use tffitcode
      use tfstk
      use sad_main
      use ffs_pointer, only:idelc,idvalc,idtypec
      implicit real*8(a-h,o-z)
      parameter (mfitc1=32,floor=1d-15,ceil=1d15,halfpi=pi*0.5d0)
      parameter (qupi=pi*0.25d0,qupi3=pi*0.75d0,qupi5=pi*1.25d0
     z,                                         qupi7=pi*1.75d0)
      parameter (cutl=1d-5)
      logical stab,xplane,yplane,cor(*),mhogan,pvert,bump
      character name*8
      integer*8 latt(nlat),id,itw,isaw,ivar,ix,idx,ia,ib,itr,isav
      integer*4 i
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
     z,         mult(*),ipair(2,npair),ipairt(2,npair),iquad(*)
     z,         istr(nstra,4),jst(*),estr(*),imon(nmona,4)
     z,         kfit(nfc),ifitp(nfc),mfitp(nfc),fitval(nfc)
     z,         trans(4,5)
      logical errflg
      common /codcor/eptsol,errflg(4),errval(4),dpshft,optiv(18),
     1               ibckup,nmona,nmonact,nstra,nstract,
     1               itmon,itemon,itstr,itestr
      save
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
      xplane=cor(1).or.cor(3)
      yplane=cor(2).or.cor(4)
      psiy=twiss(nlat,ndim,6)-twiss(1,ndim,6)
      if( bump ) then
        call ppair(latt,npnti,npntf,6,ipair)
        nmosav=nmo
        call packpi(ipair,npair,imon,imon(1,2),imon(1,3),nmo,nmosav,
     1              .true.)
      else
        npntk=npnti
        npnti=1
        do 36 i=1,nstr
          j=istr(i,2)
          if(istr(j,1).eq.npntk) then
            ipntk=i
            goto 38
          endif
   36   continue
        write(lfno,*) ' Inappropriate steering specified'
        if(lfno.ne.6) print *,' Inappropriate steering specified'
        goto 991
      endif
      call pmovi(ipair,ipairt,2*npair)
      npair0=npair
      print *,' no. of steers/monitors=',nstr,' /',nmo
c     print *,' npnti=', npnti,' npntf=',npntf
      nquad=0
      do 20 i=npnti,npntf-1
        if(idtypec(i).eq.4) then
          do 21 j=i-1,-nlat,-1
            k=j+nlat*((1-sign(1,j-1))/2)
            if(idtypec(k).ne.41.and.k.ne.nlat)then
              if(k.le.npntf-1.and.idtypec(k).eq.4) then
              else
                nquad=nquad+1
                iquad(nquad)=i
              endif
              goto 20
            endif
   21     continue
        endif
   20 continue
c     write(*,'(a,i4/(15i4))')' iquad',nquad,(iquad(i),i=1,nquad)
      nquad=0
   38 do 40 i=1,nfc
   40   mfitp(i)=1
      call pclr(fitval,nfc)
c     __________
      if(yplane) then
c     ::::::::::
        ia=italoc(nmo*(npair+nquad))
        ib=italoc(nmo)
        ix=italoc(npair+nquad)
        idx=italoc(npair+nquad)
        ivar=italoc(npair+nquad)
        isav=italoc(nstr)
        itr=italoc(npair)
        do 41 i=1,npair
          rlist(itr-1+i)=
     1         -rlist(idvalc(ipair(1,i))+2)*
     1         (twiss(ipair(1,i),0,17)+twiss(ipair(2,i),0,17)
     1         -rlist(latt(ipair(1,i))+6)
     1         -rlist(latt(ipair(2,i))+6))
   41   continue
1       continue
        do 50 i=1,nstr
c   50     rlist(isav-1+i)=rlist(latt(ist(i))+11)
   50     rlist(isav-1+i)=rlist(latt(jst(i))+11)
        do 52 i=1,nmo
          j=imon(i,2)
   52     rlist(ib-1+i)=twiss(imon(j,1),0,17)
        do 60 i=1,npair
          print *,ipair(1,i),ipair(2,i)
   60     rlist(ivar-1+i)=twiss(ipair(1,i),0,15)
        do 62 i=1,nquad
          rlist(ivar-1+npair+i)=twiss(iquad(i),0,15)
   62   continue
        if( bump ) then
          psix0=twiss(nlat,ndim,3)-twiss(1,ndim,3)
          ifitp(1)=npnti
          ifitp(2)=npnti
          kfit(1)=mfitc1
          kfit(2)=mfitc1+1
c         psix=twiss(ipair(1,1),ndim,3)
c         ifitp(3)=ipair(1,1)
          psix=twiss(ipairt(1,1),ndim,3)
          ifitp(3)=ipairt(1,1)
          kfit(3)=mfitc1
          fitval(3)=xamp
          lc=3
          do 64 i=2,npair0
c           dpsix=mod(twiss(ipair(1,i),ndim,3)-psix+psix0,psix0)
            dpsix=mod(twiss(ipairt(1,i),ndim,3)-psix+psix0,psix0)
            dpsix=mod(dpsix,pi2)
            if(dpsix.lt.qupi .or. dpsix.gt.qupi7) then
              lc=lc+1
              fitval(lc)=fitval(lc-1)
            elseif(dpsix.gt.qupi3 .and. dpsix.lt.qupi5) then
              lc=lc+1
              fitval(lc)=-fitval(lc-1)
            else
              goto 64
            endif
            ifitp(lc)=ipairt(1,i)
            kfit(lc)=mfitc1
            psix=twiss(ipairt(1,i),ndim,3)
   64     continue
          lc=lc+1
          ifitp(lc)=npntf
          kfit(lc)=mfitc1
          lc=lc+1
          ifitp(lc)=npntf
          kfit(lc)=mfitc1+1
c         print *,' bump condition, nstr=',nstr
c         write(*,'(2I4,1PD12.3)') (ifitp(i),kfit(i),fitval(i),i=1,lc)
          call pbump(latt,twiss,gammab,ndim,mult,istr,estr,nstr,kfit,
     1               ifitp,mfitp,fitval,lc,.true.,.false.,.true.,lfno)
cdebug------>
c         if(npair.eq.npair0-1) goto 199
cdebug------>
        else
          kk=istr(1,1)
          kkk=istr(1,2)
          istr(1,1)=npntk
          istr(1,2)=1
          ll=imon(1,1)
          lll=imon(1,2)
          imon(1,1)=ipairt(1,1)
          imon(1,2)=1
          call mcrmat(latt,twiss,gammab,ndim,0d0,0d0,trans,1,
     z                .false.,.false.,istr,1,ipairt(1,1),1,'X')
          istr(1,1)=kk
          istr(1,2)=kkk
          imon(1,1)=ll
          imon(1,2)=lll
          rlist(latt(npntk)+11)=rlist(latt(npntk)+11)
     z                            - xamp/trans(1,1)
        endif
        call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
        if(.not.stab) then
          call pundo(latt,twiss,gammab,0,dp0+1d0)
          goto 199
        endif
        goto 199
        psix=twiss(nlat,0,3)-twiss(1,0,3)
        psiy=twiss(nlat,0,6)-twiss(1,0,6)
        do 110 jp=1,npair
          id=idvalc(ipair(1,jp))+5
          dx=rlist(id)
          rlist(id)=halfpi
c         call mcrmat(latt,twiss,gammab,-ndim,0d0,0d0,
c    1                rlist(ia+nmo*(jp-1)),nmo,
          kk=istr(1,1)
          kkk=istr(1,2)
          istr(1,1)=ipair(1,jp)
          istr(1,2)=1
          call mcrmat(latt,twiss,gammab,0    ,0d0,0d0,
     1                rlist(ia+nmo*(jp-1)),nmo,
     1                .false.,.false.,istr,1,imon,nmo,'Y')
          istr(1,1)=kk
          istr(1,2)=kkk
          rlist(id)=dx
  110   continue
        j=ia+nmo*(npair-1)
        do 112 l=1,nquad
          id=idvalc(iquad(l))+5
          j=j+nmo
          dx=rlist(id)
          rlist(id)=halfpi
          kk=istr(1,1)
          kkk=istr(1,2)
          istr(1,1)=iquad(l)
          istr(1,2)=1
          call mcrmat(latt,twiss,gammab,-ndim,0d0,0d0,rlist(j),nmo,
     z                .false.,.false.,istr,1,imon,nmo,'Y')
          istr(1,1)=kk
          istr(1,2)=kkk
          rlist(id)=dx
  112   continue
c       print *,' response matrix'
c       write(*,'(5e12.4)')((rlist(ia+nmo*(jp-1)+j),j=0,nmo-1),
c    *                                        jp=1,npair+nquad)
        write(lfno,*) ' ddx at sext & quad'
        do 120 i=1,npair
          l=ia+(i-1)*nmo-1
          w=twiss(ipair(1,i),0,15)-rlist(ivar-1+i)
          do 121 j=1,nmo
  121       rlist(l+j)=rlist(l+j)*w
          call mwbufr(w,'(1PD12.4)',5,lfno)
  120   continue
        j=ivar-1+npair
        do 122 i=1,nquad
            j=j+1
            l=l+nmo
            w=twiss(iquad(i),0,15)-rlist(j)
            do 123 k=1,nmo
  123         rlist(l+k)=rlist(l+k)*w
            call mwbufr(w,'(1PD12.4)',5,lfno)
  122   continue
        call mwbufr(0.d0,' ',-5,lfno)
c        ... measure vertical COD change
        if(errval(4).eq.0d0) then
          do 124 i=1,nmo
            j=imon(i,2)
            rlist(ib-1+i)=twiss(imon(j,1),0,17)-rlist(ib-1+i)
  124     continue
        else
          do 125 i=1,nmo
            j=imon(i,2)
            rlist(ib-1+i)=twiss(imon(j,1),0,17)-rlist(ib-1+i)
     z                   + errval(4)*tgauss()
  125       continue
        endif
        write(lfno,*)' ddy at monitors'
        do 126 i=1,nmo
  126     call mwbufr(rlist(ib-1+i),'(1PD12.4)',5,lfno)
        call mwbufr(0.d0,' ',-5,6)
        call mstat(rlist(ib),nmo,ddymax,ddymin,ddyrms,ddyave)
        write(lfno,'(a,1pd12.4)')'   rms=',ddyrms
c       ---- solve ----
c       ... weight factor
        idb=italoc(nmo)
        do 128 i=1,nmo
          j=istr(i,2)
          fact=sqrt(twiss(imon(j,1),0,5))
          rlist(ib-1+i)=rlist(ib-1+i)*fact
          rlist(idb-1+i)=errval(4)*fact
          do 128 jp=1,npair+nquad
            l=ia+i-1+nmo*(jp-1)
            rlist(l)=rlist(l)*fact
  128   continue
        if(nmo.lt.npair+nquad) then
          call tsolvg(rlist(ia),rlist(ib),rlist(ix),nmo,npair+nquad,
     *                                                       nmo)
        else
          iu=italoc(nmo*(npair+nquad))
          iw=italoc(npair+nquad)
          iv=italoc((npair+nquad)**2)
          call plssvd(rlist(ia),rlist(ib),rlist(ix),rlist(idb),
     *                rlist(idx),chi,nmo,npair+nquad,nmo,npair+nquad,
     *                rlist(iu),rlist(iw),rlist(iv),eptsol)
          call tfree(int8(idb))
          call tfree(int8(iv))
          call tfree(int8(iw))
          call tfree(int8(iu))
        endif
c       --- reset horizontal bump
        if( bump ) then
          do 130 i=1,nstr
c  130       rlist(latt(ist(i))+11)=rlist(isav-1+i)
  130       rlist(latt(jst(i))+11)=rlist(isav-1+i)
        else
c          rlist(latt(ist(ipntk))+11)=rlist(isav-1+ipntk)
          rlist(latt(jst(ipntk))+11)=rlist(isav-1+ipntk)
        endif
c       ------------------------
        if(nmo.ge.npair+nquad) then
          write(lfno,'(a,1pg10.3,a)')
     *         ' Guess of coupling strength: (chi=',chi,')'
        else
          write(lfno,'(a)')
     *         ' Guess of coupling strength:'
        endif
        strgmx=-1d15
        do 132 i=1,npair+nquad
          if(i.le.npair.and.abs(rlist(ix-1+i)).gt.strgmx) then
            strgmx=abs(rlist(ix-1+i))
            lpmx=i
          endif
          call mwbufr(rlist(ix-1+i),'(1PD11.3)',6,lfno)
          if(nmo.ge.npair+nquad) then
            call mwbufr(rlist(idx+i-1),'(1PD11.3)',6,lfno)
          endif
  132   continue
        call mwbufr(0.d0,' ',-6,lfno)
        write(*,'(a/('' '',1p,6e11.4))') ' Guess/True :',
     1   (rlist(ix+i)/
     1     (sign( max(-rlist(itr+i),max(rlist(itr+i),1d-31))
     1            ,rlist(itr+i) )
     1     ),i=0,npair-1)
        if(strgmx.lt.cutl) then
          goto 199
        endif
        write(lfno,*) ' Kicks: pair#=',lpmx
        do 140 jp=1,npair
          if(jp.ne.lpmx) goto 140
          if(abs(rlist(ix-1+jp)).gt.cutl) then
            ns=0
            dpsi=ceil
            do 142 i=1,nstra
              if(jst(i).gt.0.and. pvert(latt,jst(i)) ) then
                if(mhogan(ipair(1,jp),ipair(2,jp)+1,jst(i))) then
                  ns=ns+1
                  dmu=mod(
     z                   twiss(ipair(2,jp)+1,ndim,6)
     z                    - twiss(jst(i),ndim,6)+psiy,psiy)
                  if(abs(dmu-halfpi).le.dpsi) then
                     js=i
                     dpsi=abs(dmu-halfpi)
                  endif
                endif
              endif
  142       continue
            if(ns.gt.0) then
              call pfmat(twiss,ndim,jst(js),ipair(2,jp),trans)
              trans(3,4)=trans(3,4)
c    z                   + 0.5d0*rlist(latt(ipair(2,jp))+1)*trans(4,4)
c             akick=rlist(ix-1+jp)/rlist(idval(ilist(2,latt(ipair(1,jp))))+2)
              akick=rlist(ix-1+jp)/rlist(latt(ipair(1,jp))+2)
     z                           /trans(3,4)
     z                    /sin(rlist(idvalc(jst(js))+5))
c             print *,' jst(js)=',jst(js)
              call elname(jst(js),name)
              print *,name
              call mwbufr(akick,'(1PD12.4)',5,lfno)
              rlist(latt(jst(js))+11)=rlist(latt(jst(js))+11)
     z                                + akick*(1d0+estr(istr(js,2)))
            else
            endif
          else
            call mwbufr(0d0,'(1PD12.4)',5,lfno)
          endif
  140   continue
        call mwbufr(0.d0,' ',-5,lfno)
        if( npair.ge.2) then
          do 150 i=lpmx+1,npair
            ipair(1,i-1)=ipair(1,i)
  150       ipair(2,i-1)=ipair(2,i)
          npair=npair-1
cdebug---->
          if(npair.eq.npair0-1) goto 199
cdebug<----
          call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
          if(stab) then
          else
            call pundo(latt,twiss,gammab,0,dp0+1d0)
            goto 199
          endif
          goto 1
        endif
  199   npair=npair0
        call pmovi(ipairt,ipair,2*npair)
c     _____
      endif
c     :::::
c     __________
      if(xplane) then
c     ::::::::::
      elseif(yplane) then
        goto 999
      else
        goto 991
c     _____
      endif
c     :::::
  999 continue
      call tfree(itr)
      call tfree(isav)
      call tfree(ivar)
      call tfree(ix)
      call tfree(idx)
      call tfree(ib)
      call tfree(ia)
      npair=npair0
      if( bump ) then
        call packpi(ipair,npair,imon,imon(1,2),imon(1,3),nmo,nmosav,
     1               .false.)
        nmo=nmosav
      endif
  991 continue
      return
      end

      subroutine palocx(ia,na,nna)
      use tfstk
      use ffs
      use tffitcode
      integer*8 ia,ib
      ib=ktaloc(nna)
      call tmov(rlist(ia),rlist(ib),min(na,nna))
      call tfree(ia)
      ia=ib
      return
      end
C   18/03/92 203182111  MEMBER NAME  PAMOEB   *.FORT     M  E2FORT

      subroutine pamoeb(p,
c    <---- funk was replaced by -------------------------------------
     &     funk,id,latt,twiss,mult,gammab,size,istr,estr,nstr,observ,
     $     iobs,synchb,dely,lfno,
c    --------------------------------------------------------------->
     &     y,ndim,pr,prr,pbar,ftol,itmax)
      implicit real*8 (a-h,o-z)
      character shape*8,autofg*11
      parameter (alpha=1.d0,beta=0.5d0,gamma=2.d0)
      character*(*) observ
      dimension p(ndim,ndim+1),id(ndim),y(ndim+1),pr(ndim),prr(ndim),
     1          pbar(ndim)
      dimension latt(*),twiss(*),mult(*),gammab(*),size(*),istr(*),
     $     estr(*),iobs(*)
      logical synchb
c     <----- output status -------------------
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
c     --------------------------------------->
      save ical
      external funk
c
      mpts=ndim+1
c     <------ Count Ncall -------
      iterint=0
      ncall=0
      ncalla=0
c     <-------- Shape -----------
      shape=' '
c     <----- Output status ------
      ylow=0d0
c     -------------------------->
    1 ilo=1
      if(y(1).gt.y(2))then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 11 i=1,mpts
        if(y(i).lt.y(ilo)) ilo=i
        if(y(i).gt.y(ihi))then
          inhi=ihi
          ihi=i
        elseif(y(i).gt.y(inhi))then
          if(i.ne.ihi) inhi=i
        endif
   11 continue
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
c     <------------ output status --------------------------------------
      if(iter.eq.0.or.iterint.ne.0) then
c.......write status if iter=0 (petune is called with INIT) or
c                    if iterint.ne.0 (iterated process in the Pamoeb)
        write(lfno,'(2(a,i4),a,1pg11.4,2(a,1pd11.4),a)')
     &       ' Iter=',iter,' Ncall=',nemitd,' rtol=',rtol,' y(hi)=',
     $       y(ihi),' y(lo)=',y(ilo),' '//shape
        do j=1,mpts
          if(j.eq.ihi .or. j.eq.ilo) then
            do i=1,ndim
              call mbufg(autofg(p(i,j),'11.4'),.false.,lfno)
            enddo
            call mbufw(autofg(y(j),'11.4'),.false.,lfno)
            if(j.eq.ihi) call mbufw('<Hi>',.false.,lfno)
            if(j.eq.ilo) call mbufw('<Lo>',.false.,lfno)
            call mbufw(' ',.true.,lfno)
          endif
        enddo
      endif
      do l=0,ncall-ncalla-1
        call pmbdata(latt,mult,iobs,iobs(1),observ,nemitd-l,y(ihi),
     $       .false.,lfno)
      enddo
      if(ical.ne.0) then
        call pmbdata(latt,mult,iobs,iobs(1),observ,ical,0d0,.true.,lfno)
      endif
      ncalla=ncall
      ical=0
      if(iter.eq.0) then
c       .. This line is needed to print always the first-step value.:
c          The first step after initialization is expected normally to
c          get LOWER than the minimum in the initialization step. When
c          it is EQUAL to the minimum, we put y(ihi) into ylow in order
c          to print out the first step value in any case
        ylow=y(ihi)
      else
        ylow=y(ilo)
      endif
c     ----------------------------------------------------------------->
      if(rtol.lt.ftol) return
      if(iterint.eq.itmax) then
c       print *,' Pameob exceeding maximum iterations',itmax
        return
      endif
c     <------ Count Ncall -----------
      iterint=iterint+1
      iter=iter+1
c     -------------------------->
      do 12 j=1,ndim
        pbar(j)=0d0
   12 continue
      do 14 i=1,mpts
        if(i.ne.ihi)then
          do 13 j=1,ndim
            pbar(j)=pbar(j)+p(j,i)
   13     continue
        endif
   14 continue
      do 15 j=1,ndim
        pbar(j)=pbar(j)/ndim
        pr(j)=(1.d0+alpha)*pbar(j)-alpha*p(j,ihi)
   15 continue
c     ypr=funk(pr)
      ypr=funk(pr,id,ndim,latt,twiss,mult,gammab,size,istr,estr,nstr,
     $     observ,iobs,.false.,synchb,dely,lfno)
c     <------- Count Ncall -------
      ncall=ncall+1
c     --------------------------->
      if(ypr.le.y(ilo))then
        do 16 j=1,ndim
          prr(j)=gamma*pr(j)+(1.d0-gamma)*pbar(j)
   16   continue
c       yprr=funk(prr)
        yprr=funk(prr,id,ndim,latt,twiss,mult,gammab,size,istr,estr,
     $       nstr,observ,iobs,.false.,synchb,dely,lfno)
c       <------- Count Ncall -------
        ncall=ncall+1
c       --------------------------->
        if(yprr.lt.y(ilo))then
          do 17 j=1,ndim
            p(j,ihi)=prr(j)
   17     continue
          y(ihi)=yprr
c         <-------- Shape ------------
          shape='extend'
          ical=nemitd
c         --------------------------->
        else
          do 18 j=1,ndim
            p(j,ihi)=pr(j)
   18     continue
          y(ihi)=ypr
c         <-------- Shape ------------
          shape='reverse'
          ical=nemitd-1
c         --------------------------->
        endif
      elseif(ypr.ge.y(inhi))then
        if(ypr.lt.y(ihi))then
          do 19 j=1,ndim
            p(j,ihi)=pr(j)
   19     continue
          y(ihi)=ypr
        endif
        do 21 j=1,ndim
          prr(j)=beta*p(j,ihi)+(1.d0-beta)*pbar(j)
   21   continue
c       yprr=funk(prr)
        yprr=funk(prr,id,ndim,latt,twiss,mult,gammab,size,istr,estr,
     $       nstr,observ,iobs,.false.,synchb,dely,lfno)
c       <------- Count Ncall -------
        ncall=ncall+1
c       --------------------------->
        if(yprr.lt.y(ihi))then
          do 22 j=1,ndim
            p(j,ihi)=prr(j)
   22     continue
          y(ihi)=yprr
c         <-------- Shape ------------
          shape='collapse'
          if(yprr.lt.y(ilo)) ical=nemitd
c         --------------------------->
        else
          do 24 i=1,mpts
            if(i.ne.ilo)then
              do 23 j=1,ndim
                pr(j)=0.5d0*(p(j,i)+p(j,ilo))
                p(j,i)=pr(j)
   23         continue
c             y(i)=funk(pr)
              y(i)=funk(pr,id,ndim,latt,twiss,mult,gammab,size,istr,
     $             estr,nstr,observ,iobs,.false.,synchb,dely,lfno)
c             <------ Count Ncall ---------
              ncall=ncall+1
              if(y(i).lt.ylow) then
                ylow=y(i)
                ical=nemitd
              endif
c             ---------------------------->
            endif
   24     continue
c         <-------- Shape ------------
          shape='shrink'
c         --------------------------->
        endif
      else
        do 25 j=1,ndim
          p(j,ihi)=pr(j)
   25   continue
        y(ihi)=ypr
c       <-------- Shape ------------
        shape='reverse'
c       --------------------------->
      endif
      goto 1
      end

      subroutine pasex(word,wordp,latt,pos,twiss,mult,master,
     $     gammab,title,case,istr,estr,nstr,imon,emon,nmon,
     $     lfno)
      use tfstk
c----- Beam-based alignment of sextupoles. -----------------------------
c     Principle: Large orbit excursion at sextupoles in the horizontal
c     plane makes orbit distorsion in the vertical plane. One can
c     estimate the vertical misalignments of sextupoles by analyzing
c     the vertical orbit:
c       yb''-K*yb=K2(xb*yb + xb*(y0-dy) + (x0-dx)*yb)
c     yb is the orbit change in the vertical plane due to a large
c     orbit in the horizontal plane. y0 is the original orbit. 
c     Assumptions: The base orbit is stored in the current buffer.
c     The orbit difference from the base is stored in the 'M'
c     stacks and the corresponding corrector is stored in the 'C' stacks
c-----------------------------------------------------------------------
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec,pnamec
      implicit real*8 (a-h,o-z)
      parameter (kstack=3)
      parameter (lioo=50,lstdo=6,nminc=20)
      logical pass,tmatch,lod,gout,result,exist,mhogal,subbase,hpos,
     $     micado,vonly
      character*(*) word,wordp,title,case
      character*16 patt
      integer*8 latt(nlat)
      dimension pos(nlat),twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat)
      dimension istr(nstra,4),imon(nmona,4),emon(*),estr(*),mult(*),
     $     master(nlat)
      dimension vscale(2)
      external pack
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      include 'inc/common.inc'
      data lio/lioo/
      data isex/0/
c
      save
c
c
      gout=.false.
      subbase=.false.
      hpos=.false.
      micado=.false.
      vonly=.false.
      vscale(1)=1d20
      vscale(2)=-1d20
      ndat=1
      lod=.false.
c
 1    lod=.false.
 2    call getwdl2(word,wordp)
      if(word.eq.'RESULT') then
        result=.true.
        go to 1
      elseif(word.eq.'SB') then
c       display a solution subtracted the base-orbit (i.e., dx and dy)
        subbase=.true.
        go to 1
      elseif(word.eq.'HPOS') then
c       horizontal axis is position, otherwise sextupole number.
        hpos=.true.
        go to 1
      elseif(word.eq.'VSCALE') then
c       vertical scale in the graphics.
        vscale(1)=getva(exist)
        if(exist) then
          vscale(2)=getva(exist)
          if(min(vscale(1),vscale(2)).eq.vscale(1)) then
            r=vscale(1)
            vscale(1)=vscale(2)
            vscale(2)=r
          endif
        endif
        go to 1
      elseif(word.eq.'MICADO') then
        micado=.true.
        nmicad=getva(exist)
        if(.not.exist) nmicad=nsex
        go to 1
      elseif(word.eq.'VONLY') then
        vonly=.true.
        go to 1
      elseif(word.eq.'DATA') then
        ndat=getva(exist)*2
        if(.not.exist) ndat=2
        goto 1
      elseif(word.eq.'GOUT') then
        lio=getva(exist)
        if(exist) then
          gout=.true.
          go to 1
        else
          lio=lioo
          gout=.true.
          lod=.true.
          go to 2
        endif
      elseif(lod) then
        if(word.eq.' ') go to 991
        inquire(unit=lio,opened=lod,iostat=ios)
        if(lod) then
          if(lio.ne.lstdo) close(lio)
          lod=.not.lod
        endif
        call texpfn(word)
        open (lio,file=word,iostat=ios,status='UNKNOWN',err=991)
        if(ios .ne. 0) then
          call permes(' Cannot open',' ',word,lfno)
          return
        endif
        go to 1
 991    call permes(' Error in Opening',' ',word,lfno)
        return
      else
        ls=1
        lf=nlat
        exist=.true.
        if(word.eq.' ') then
          patt='*'
        else
          patt=wordp
          call getwdl2(word,wordp)
          if(word.eq.' ') then
            exist=.false.
          else
            ls=ielm(wordp,exist)
            call getwdl2(word,wordp)
            if(word.eq.' ') then
              exist=.false.
            else
              lf=ielm(wordp,exist)
            endif
          endif
        endif
        if(isex.eq.0) then
 10       nsex=0
          do i=1,nlat-1
            if(idtypec(i).eq.icsext
     $           .and. master(i).gt.0) then
              if(tmatch(pnamec(i),patt)) then
                if(mhogal(ls,lf,i)) then
                  nsex=nsex+1
                  if(isex.ne.0) then
                    ilist(mod(nsex-1,2)+1,isex+(nsex-1)/2)=i
                  endif
                endif
              endif
            endif
          enddo
          if(isex.eq.0) then
            isex=italoc((nsex+1)/2)
            go to 10
          endif
          if(.not.exist) go to 20
          if(nsex.ne.0) go to 1
        endif
      endif
c      
 20   continue
      if(result) go to 100
c
      it=italoc(2*nlat)
      call tmov(twiss(1,0,15),rlist(it),nlat)
      call tmov(twiss(1,0,17),rlist(it+nlat),nlat)
      call pcbak(latt,twiss)
      im=italoc((nmona+1)/2)
      call pmovi(imon(1,3),rlist(im),nmona)
      ix=italoc(2*nsex)
      ixs=italoc(4*nsex)
      ixs0=italoc(4*nsex)
      ixsb=italoc(2*nsex)
      iex=italoc(2*nsex)
c
c     ndata=1
c     if(iud.eq.0) then
c       iud=italoc(2*nsex*nminc)
c       nmax=nminc
c     elseif(ndata.ge.nmax) then
c       call palocx(iud,2*nsex*nmax,2*nsex*(nmax+nminc))
c       nmax=nmax+nminc
c     endif
      ip=ipnt(3)
      idima=0
      do i=1,ndat,2
        lp=ilist(mod(ip-ldep-1,2)+1,iistck(3)+(ip-ldep-1)/2)
        idima=idima+ilist(1,lp)
      enddo
      ia=italoc(idima*nsex*2)
      ib=italoc(idima)
c     pass=.true.
c     ia=italoc(1)
c     ib=italoc(1)
c     idima=1
c     do i=1,2
        call pasex1(latt,twiss,mult,gammab,rlist(ia),
     $       idima,rlist(ib),rlist(ix),rlist(iex),rlist(isex),nsex,
     $       rlist(ixs),rlist(ixs0),rlist(ixsb),istr,estr,nstr,
     $       rlist(iistck(3)),rlist(iistck(1)),imon,emon,nmon,ndat,
     $       pass,vonly,micado,nmicad,subbase,lio,lfno)
c       if(pass) then
c         call tfree(int8(ia))
c         ia=italoc(idima*nsex*2)
c         print *,'idima=',idima
c         call tfree(int8(ib))
c         ib=italoc(idima)
c       endif 
c       pass=.not.pass
c     enddo
      call tfree(int8(ia))
      call tfree(int8(ib))
      if(gout) then
        call pasexg(latt,pos,mult,title,case,rlist(ix),rlist(iex),
     $       rlist(isex),nsex,vscale,hpos,vonly,lio)
      endif
c     
c      call tmov(rlist(it),twiss(1,0,15),nlat)
c      call tmov(rlist(it+nlat),twiss(1,0,17),nlat)
c      call pundo(latt,twiss,gammab,0,dp0+1d0)
      do i=1,nlat-1
        if(idtypec(i).eq.icbend) then
          rlist(latt(i)+11)=rlist(ibckup-1+i)
        endif
      enddo
      call pmovi(rlist(im),imon(1,3),nmona)
      call pack(imon(1,2),imon(1,3),nmon,nmona)
      call tfree(int8(iex))
      call tfree(int8(ix))
      call tfree(int8(ixs))
      call tfree(int8(ixs0))
      call tfree(int8(ixsb))
      call tfree(int8(im))
      call tfree(int8(it))
 100  continue

c     if(iud.ne.0) then
c       call tfree(int8(iud))
c       iud=0
c     endif
      if(isex.ne.0) then
        call tfree(int8(isex))
        isex=0
      endif
      if(gout.and.lio.ne.lstdo) then
        close(lio)
      endif
      return
      end
c
      subroutine pasex1(latt,twiss,mult,
     $     gammab,a,ia,b,x,ex,isex,nsex,xs,xs0,xsb,
     $     istr,estr,nstr,id,idc,imon,emon,nmon,ndat,pass,vonly,micado,
     $     nmicad,subbase,lio,lfno)
      use tfstk
c***** K'*xb*yb term is included in the response matrix ****************
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc
      implicit real*8 (a-h,o-z)
      parameter (kstack=3,halfpi=0.5d0*pi,tiny=1d-20)
      logical pass,corc(8),vonly,micado,stab,subbase
      character name*8,autofg*11
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
      dimension a(ia,*),b(*),x(2*nsex),ex(2*nsex),mult(*)
      dimension xs(4,nsex),xs0(4,nsex),xsb(2,nsex),isex(nsex)
      dimension istr(nstra,4),imon(nmona,4),emon(*),estr(*)
      dimension id(*),idc(*)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      include 'inc/common.inc'
      data corc/.false.,.true.,6*.false./
c
c      print *,'mfaloc=',mfalloc(-1)
c      print *,'nsex=',nsex
c     ip=ipnt(3)
c     ic=ipnt(1)
c     if(pass) then
c       ia=0
c       do ldep=1,ndat,2
c         lp=id(ip-ldep)
c         ia=ia+ilist(1,lp)
c       enddo
c       return
c     endif
c.....Get base orbit from current buffer......
      if(simulate) then
        do i=1,nsex
          xsb(1,i)=twiss(isex(i),0,15)
          xsb(2,i)=twiss(isex(i),0,17)
        enddo
      else
        call pclr(xs0,4*nsex)
        call pasex2(latt,twiss,gammab,isex,nsex,xs0,istr,estr,
     $       nc,imon,emon,no,1,lfno)
        do i=1,nsex
c     ++ assume sextupole is a single piece
          xsb(1,i)=xs0(1,i)+
     $         0.5d0*rlist(idvalc(isex(i)-1)+1)*xs0(2,i)
          xsb(2,i)=xs0(3,i)+
     $         0.5d0*rlist(idvalc(isex(i)-1)+1)*xs0(4,i)
        enddo
      endif
c.....Get data from stacks .......
      la=1
      do 100 ldep=1,ndat,2
c...... estimate positions at sextupoles from BPM data for base orbit.....
        call mfetch(3,ldep,latt,twiss,istr,nc,imon,no)
        call mfetch(1,ldep,latt,twiss,istr,nc,imon,no)
        call pclr(xs0,4*nsex)
        call pasex2(latt,twiss,gammab,isex,nsex,xs0,istr,estr,
     $       nc,imon,emon,no,1,lfno)
        do i=1,nsex
c         ++ assume sextupole is a single piece
          xs0(1,i)=xs0(1,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs0(2,i)
          xs0(3,i)=xs0(3,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs0(4,i)
        enddo
c...... fetch stack-'M' to current buffer. ......
        call mfetch(3,ldep+1,latt,twiss,istr,nc,imon,no)
        do i=1,no
          b(la-1+i)=twiss(imon(imon(i,2),1),0,17)
        enddo
c...... fetch stack-'C' to current buffer. .....
        call mfetch(1,ldep+1,latt,twiss,istr,nc,imon,no)
c...... estimate positions at sextupoles from BPM data. ...
        call pclr(xs,4*nsex)
        call pasex2(latt,twiss,gammab,isex,nsex,xs,istr,estr,nc,
     $       imon,emon,no,2,lfno)
        do 16 i=1,nsex
c         ++ assume sextupole is a single piece
          xs(1,i)=xs(1,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs(2,i)
          xs(3,i)=xs(3,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs(4,i)
 16     continue
c        print *,'to check pasex2'
c        write(*,'(1p,10d11.3)')(xs(1,i)-twiss(isex(i),0,15),xs(3,i)-
c    $        twiss(isex(i),0,17),i=1,nsex)
c...... construct response matrix ......
c   --- make optics that corresponds to large horizontal orbit . ---
        print *,'no=',no,' nc=',nc,' ia=',ia
c.......Include K'*xb*yb to the reaponse matrix. ........
c       assume thin quads are attached to both ends of sexts.
        call perr(latt,.true.)
        do i=1,nsex
          ls=isex(i)
          alh=0.5d0*rlist(idvalc(ls)+1)
          akh=0.5d0*rlist(idvalc(ls)+2)
c         rlist(latt(ls-1)+2)=akh*(xs(1,i)-alh*xs(2,i))
c         rlist(latt(ls+1)+2)=akh*(xs(1,i)+alh*xs(2,i))
          rlist(latt(ls-2)+2)=akh*(xs(1,i)-alh*xs(2,i))
          rlist(latt(ls+2)+2)=akh*(xs(1,i)+alh*xs(2,i))
          rlist(latt(ls  )+2)=0d0
        enddo
        do i=1,ntwissfun
          twiss(1,-ndim,i)=twiss(1,ndim,i)
        enddo
        call pqcell(latt,twiss,gammab,-ndim,dp0+1d0,stab)
        call perr(latt,.false.)
        do 20 i=1,nsex
c         ++ assume sextupole is a single piece
          ls=isex(i)
          r=rlist(idvalc(ls)+5)
          rlist(idvalc(ls)+5)=
     $         -halfpi+rlist(idvalc(ls)+4)
          k=istr(istr(1,2),1)
          istr(istr(1,2),1)=ls
          call mcrmat(latt,twiss,gammab,-ndim,0d0,0d0,
     $         a(la,i),ia,.false.,.false.,istr,1,imon,no,'Y')
          istr(istr(1,2),1)=k
          rlist(idvalc(ls)+5)=r
          call tmov(a(la,i),a(la,i+nsex),no)
          ak=rlist(idvalc(ls)+2)*
c    $         (xs(1,i)*(xs(3,i)+xs0(3,i))+xs(3,i)*xs0(1,i))
     $         (xs(1,i)*xs0(3,i)+xs(3,i)*xs0(1,i))
          do j=1,no
            b(la-1+j)=b(la-1+j)-ak*a(la-1+j,i)
          enddo
c         ++ assume sextupole is a single piece
          ak=rlist(idvalc(ls)+2)
c         ak=ak*xs(3,i)*(xs(1,i)+ak*al/12d0*(xs(1,i)**2+xs(3,i)**2))
          call ptimes(ak*xs(1,i),a(la,i),no)
          call ptimes(ak*xs(3,i),a(la,i+nsex),no)
 20     continue
        call mwght(twiss,ndim,a(la,1),b(la),no,2*nsex,ia,imon,no,corc)
        la=la+no
 100  continue
c.... solve .....
      if(vonly) then
        nvari=nsex
      else
        nvari=2*nsex
      endif
c     call tsolva(a,b,x,ia,nvari,ia,eptsol)
      if(micado) then
        call msolvg(a,b,x,ia,nvari,ia,0,
     &       .false.,.true.,nmicad,.true.,.true.)
        do i=1,nvari
          ex(i)=b(i)
        enddo
      else
        call tsvd(a,b,x,ia,nvari,ia,eptsol,.true.)
        do i=1,nvari
          err=0d0
          do j=1,nvari
            err=err+(b(j)*a(j,i))**2
          enddo
          err=sqrt(err)
          ex(i)=err
        enddo
      endif
c     call psolva(a,b,x,ia,nvari,ia,eptsol)
c     print *,'Input number of variable:'
      la=1
      do 220 ldep=1,ndat,2
        call mfetch(3,ldep,latt,twiss,istr,nc,imon,no)
        call mfetch(1,ldep,latt,twiss,istr,nc,imon,no)
        call pclr(xs0,4*nsex)
        call pasex2(latt,twiss,gammab,isex,nsex,xs0,istr,estr,
     $       nc,imon,emon,no,1,lfno)
        do i=1,nsex
c         ++ assume sextupole is a single piece
          xs0(1,i)=xs0(1,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs0(2,i)
          xs0(3,i)=xs0(3,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs0(4,i)
        enddo
        call mfetch(3,ldep+1,latt,twiss,istr,nc,imon,no)
        do i=1,no
          b(la-1+i)=twiss(imon(imon(i,2),1),0,17)
        enddo
        call mfetch(1,ldep+1,latt,twiss,istr,nc,imon,no)
        call pclr(xs,4*nsex)
        call pasex2(latt,twiss,gammab,isex,nsex,xs,istr,estr,nc,
     $       imon,emon,no,2,lfno)
        do i=1,nsex
          xs(1,i)=xs(1,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs(2,i)
          xs(3,i)=xs(3,i)+
     $         0.5d0*rlist(idvalc(isex(i))+1)*xs(4,i)
        enddo
        call perr(latt,.true.)
        do i=1,nsex
          ls=isex(i)
          alh=0.5d0*rlist(idvalc(ls)+1)
          akh=0.5d0*rlist(idvalc(ls)+2)
c         rlist(latt(ls-1)+2)=akh*(xs(1,i)-alh*xs(2,i))
c         rlist(latt(ls+1)+2)=akh*(xs(1,i)+alh*xs(2,i))
          rlist(latt(ls-2)+2)=akh*(xs(1,i)-alh*xs(2,i))
          rlist(latt(ls+2)+2)=akh*(xs(1,i)+alh*xs(2,i))
          rlist(latt(ls  )+2)=0d0
        enddo
        do i=1,ntwissfun
          twiss(1,-ndim,i)=twiss(1,ndim,i)
        enddo
        call pqcell(latt,twiss,gammab,-ndim,dp0+1d0,stab)
        call perr(latt,.false.)
        do 200 i=1,nsex
          ls=isex(i)
          r=rlist(idvalc(ls)+5)
          rlist(idvalc(ls)+5)=
     $         -halfpi+rlist(idvalc(ls)+4)
          k=istr(istr(1,2),1)
          istr(istr(1,2),1)=ls
          call mcrmat(latt,twiss,gammab,-ndim,0d0,0d0,
     $         a(la,i),ia,.false.,.false.,istr,1,imon,no,'Y')
          istr(istr(1,2),1)=k
          rlist(idvalc(ls)+5)=r
          call tmov(a(la,i),a(la,i+nsex),no)
          ak=rlist(idvalc(ls)+2)*
c    $         (xs(1,i)*(xs(3,i)+xs0(3,i))+xs(3,i)*xs0(1,i))
     $         (xs(1,i)*xs0(3,i)+xs(3,i)*xs0(1,i))
          do j=1,no
            b(la-1+j)=b(la-1+j)-ak*a(la-1+j,i)
          enddo
          ak=rlist(idvalc(ls)+2)
          call ptimes(ak*xs(1,i),a(la,i),no)
          call ptimes(ak*xs(3,i),a(la,i+nsex),no)
 200    continue
        la=la+no
 220  continue
      chisq=0d0
      do i=1,ia
        r=0d0
        do j=1,nvari
          r=r+a(i,j)*x(j)
        enddo
        chisq=chisq+(b(i)-r)**2
      enddo
      chisq=chisq/dble(ia-nvari)
      write(*,'(a,1pd11.3)')' Chisq=',chisq
      chisq=sqrt(chisq)
      do i=1,nvari
        ex(i)=chisq*ex(i)
      enddo
      do i=1,nsex
        call elnameK(isex(i),name)
        if(simulate) then
          dx=rlist(latt(isex(i))+5)
          dy=rlist(latt(isex(i))+6)
          x(i)=x(i)-(xsb(2,i)-dy)
          x(i+nsex)=x(i+nsex)-(xsb(1,i)-dx)
        elseif(subbase) then
          x(i)=xsb(2,i)-x(i)
          x(i+nsex)=xsb(1,i)-x(i+nsex)
        endif
        if(vonly) then
          call mbufw(name//autofg(x(i),'11.3')//
     $         autofg(ex(i),'11.3'),
     $         .false.,lfno)
        else
          call mbufw(name//autofg(x(i+nsex),'11.3')//
     $         autofg(ex(i+nsex),'11.3')//
     $         autofg(x(i),'11.3')//
     $         autofg(ex(i),'11.3'),.false.,lfno)
        endif
        twiss(isex(i),0,15)=x(i+nsex)
        twiss(isex(i),0,17)=x(i)
        twiss(isex(i),0,7)=ex(i+nsex)
        twiss(isex(i),0,9)=ex(i)
c      <<< check >>>
        rlist(latt(isex(i))+6)=-x(i)
        if(.not.vonly) rlist(latt(isex(i))+5)=-x(i+nsex)
c      <<<       >>>
      enddo
      call mbufw(' ',.true.,lfno)
      return
      end
c     
      subroutine pasex2(latt,twiss,gammab,isex,nsex,xs,istr,estr,
     $     nstr,imon,emon,no,iter,lfno)
      use tfstk
c.... estimate positions at sextupoles from BPM data. ...
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc
      implicit real*8 (a-h,o-z)
      logical exta,mhogal
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
      dimension isex(nsex),xs(4,nsex)
      dimension istr(*),estr(*),imon(nmona,4),emon(*)
      dimension tr(4,5),ra(4)
      include 'inc/common.inc'
c     begin initialize for preventing compiler warning
      px1=0.d0
      py1=0.d0
c     end   initialize for preventing compiler warning
      do 900 loop=1,iter
        exta=.false.
        ibs=1
        do 100 isx=1,nsex
          ls=isex(isx)
          al=rlist(idvalc(ls)+1)
          ak=rlist(idvalc(ls)+2)
          do 10 i=1,no-1
            if(mhogal(imon(imon(i,2),1),imon(imon(i+1,2),1),ls)) then
              l1=imon(imon(i,2),1)
              l2=imon(imon(i+1,2),1)
              go to 13
            endif
 10       continue
          if(.not.trpt) then
            l1=imon(imon(no,2),1)
            l2=imon(imon(1,2),1)
          else
            if(.not.exta .and. isx.ne.1) then
c.............extend right. ....
              do 12 i=isx,nsex
                call pasex3(twiss,latt,gammab,ra,l2,isex(i),isex(i),ibs,
     $               nstr,istr)
                call pftmat(twiss,gammab,tr,l2,isex(i),ndim)
                do 11 k=1,4
                  xs(k,i)=tr(k,1)*twiss(l2,0,15) + tr(k,2)*px1 +
c8/03$                 tr(k,3)*twiss(l2,0,17) + tr(k,4)*py1 + ra(k)
     $                 tr(k,3)*twiss(l2,0,17) + tr(k,4)*py1 +
     $                 tr(k,5)*dp0 + ra(k)
 11             continue
 12           continue
              goto 101
            endif
            exta=.true.
            go to 100
          endif
 13       continue
          call pftmat(twiss,gammab,tr,l1,l2,ndim)
          call pasex3(twiss,latt,gammab,ra,l1,l2,l2,ibs,nstr,istr)
          ibsa=ibs
c........ solve x2=m11*x1 + m12*px1 + m13*y1 + m14*py1 + m15*dp + ra1
c               y2=m31*x1 + m32*px1 + m33*y1 + m34*py1 + m35*dp + ra2
c         w.r.t. (px1,py1). ra1 and ra2 come from corrector kicks.
          a11=tr(1,2)
          a12=tr(1,4)
          a21=tr(3,2)
          a22=tr(3,4)
          c1=twiss(l2,0,15)-tr(1,1)*twiss(l1,0,15)-
     $         tr(1,3)*twiss(l1,0,17)
c8/03$         -ra(1)
     $         -tr(1,5)*dp0 -ra(1)
          c2=twiss(l2,0,17)-tr(3,1)*twiss(l1,0,15)-
     $         tr(3,3)*twiss(l1,0,17)
c8/03$         -ra(3)
     $         -tr(3,5)*dp0 -ra(3)
c       ++ include nonlinear kick of sexts. ++
          call pftmat(twiss,gammab,tr,ls,l2,ndim)
          do k=1,4
            tr(k,2)=tr(k,2)-al*tr(k,1)
            tr(k,4)=tr(k,4)-al*tr(k,3)
          enddo
          tx=-0.5d0*ak*((xs(1,isx)+al*xs(2,isx))**2-
     $         (xs(3,isx)+al*xs(4,isx))**2)
          ty=ak*(xs(1,isx)+al*xs(2,isx))*(xs(3,isx)+al*xs(4,isx))
          c1=c1-tr(1,2)*tx-tr(1,4)*ty
          c2=c2-tr(3,2)*tx-tr(3,4)*ty
c       ++
          denomi=a11*a22-a12*a21
          px1=( a22*c1-a12*c2)/denomi
          py1=(-a21*c1+a11*c2)/denomi
c
          call pftmat(twiss,gammab,tr,l1,ls,ndim)
          call pasex3(twiss,latt,gammab,ra,l1,ls,ls,ibs,nstr,istr)
          ibs=ibsa
          do 20 k=1,4
            xs(k,isx)=tr(k,1)*twiss(l1,0,15) + tr(k,2)*px1 +
c8/03$           tr(k,3)*twiss(l1,0,17) + tr(k,4)*py1 + ra(k)
     $           tr(k,3)*twiss(l1,0,17) + tr(k,4)*py1 + tr(k,5)*dp0
     $           + ra(k)
 20       continue
          if(exta) then
c ......... extend left ....
            do 22 i=isx-1,1,-1
              k=1
              call pasex3(twiss,latt,gammab,ra,isex(i),l1,isex(i),k,ibs,
     $             istr)
              call pftmat(twiss,gammab,tr,l1,isex(i),ndim)
              do 21 k=1,4
                xs(k,i)=tr(k,1)*twiss(l1,0,15) + tr(k,2)*px1 +
c8/03$               tr(k,3)*twiss(l1,0,17) + tr(k,4)*py1 + ra(k)
     $               tr(k,3)*twiss(l1,0,17) + tr(k,4)*py1 + tr(k,5)*dp0
     $               + ra(k)
 21           continue
 22         continue
            exta=.false.
          endif
 100    continue
 101    continue
 900  continue 
      return
      end
c
      subroutine pasex3(twiss,latt,gammab,ra,l1,l2,ls,is1,is2,istr)
      use tfstk
c     effect of steerings that place in between l1 and l2
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc,idtypec
      implicit real*8 (a-h,o-z)
      logical mhogan,first
      integer*8 latt(nlat)
      dimension twiss(*),gammab(*),ra(4)
      dimension istr(nstra,4)
      dimension tr(4,5),tr1(4,5)
      include 'inc/common.inc'
c
      iss=is1
      do 10 i=1,4
        ra(i)=0d0
 10   continue
      first=.true.
      do 20 i=is1,is2
        j=istr(istr(i,2),1)
        if(mhogan(l1,l2,j)) then
          if(idtypec(j).eq.icbend) then
            if(first) then
              iss=i
              first=.false.
            endif
            call pftmat(twiss,gammab,tr,j,ls,ndim)
            call pftmat(twiss,gammab,tr1,j+1,ls,ndim)
            cost=cos(-rlist(idvalc(j)+5))
            sint=sin(-rlist(idvalc(j)+5))
            ddk=rlist(latt(j)+11)
            call padd(tr1,tr,20)
            call ptimes(0.5d0,tr,20)
            do 15 k=1,4
c8/03         ra(k)=ra(k)+ddk*(tr(k,2)*cost+tr(k,4)*sint)
              ra(k)=ra(k)-ddk*(tr(k,2)*cost+tr(k,4)*sint)/(1d0+dp0)
 15         continue
          endif
        endif
 20   continue
      is1=iss
      return
      end
c
      subroutine pftmat(twiss,gammab,tr,l1,l2,ip)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(*),tr(4,5)
c   Next line is changed by Oide 10/6/1994
      call tftmat(tr,l1,l2,ip,.false.)
      if(trpt) then
        r=sqrt(gammab(l1)/gammab(l2))
        call ptimes(r,tr,20)
      endif
      return
      end
c
c      subroutine psolva(a,b,x,m,n,mp,eps)
c      include 'inc/TFMACRO.inc'
c      dimension a(mp,n),b(m),x(n)
c      iu=italoc(mp*n)
c      iw=italoc(n)
c      iv=italoc(n*n)
c      call psolvg(a,b,x,m,n,mp,n,rlist(iu),rlist(iw),rlist(iv),eps)
c      call tfree(int8(iv))
c      call tfree(int8(iw))
c      call tfree(int8(iu))
c      return
c      end
c
      subroutine perr(latt,rst)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical rst
      integer*8 latt(nlat)
      save ie
c
      if(rst) then
        ie=italoc(4*nlat)
      endif
      call perr1(latt,rst,rlist(ie))
      if(.not.rst) then
        call tfree(int8(ie))
      endif
      return
      end
c
      subroutine perr1(latt,rst,errs)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc,idtypec
      implicit real*8 (a-h,o-z)
      logical rst
      integer*8 latt(nlat),lp
      dimension errs(4,nlat)
c
      if(rst) then
        do i=1,nlat-1
          ld=idtypec(i)
          lp=latt(i)
          if(ld.eq.icbend) then
            errs(1,i)=rlist(lp+11)
            errs(2,i)=rlist(lp+5)
            rlist(lp+11)=0d0
            rlist(lp+5)=rlist(idvalc(i)+5)
          elseif(ld.eq.icquad.or.ld.eq.icsext) then
            errs(1,i)=rlist(lp+5)
            errs(2,i)=rlist(lp+6)
            errs(3,i)=rlist(lp+4)
            errs(4,i)=rlist(lp+2)
            rlist(lp+5)=0d0
            rlist(lp+6)=0d0
            rlist(lp+4)=rlist(idvalc(i)+4)
            rlist(lp+2)=rlist(idvalc(i)+2)
          endif
        enddo
      else
        do i=1,nlat-1
          ld=idtypec(i)
          lp=latt(i)
          if(ld.eq.icbend) then
            rlist(lp+11)=errs(1,i)
            rlist(lp+5)=errs(2,i)
          elseif(ld.eq.icquad.or.ld.eq.icsext) then
            rlist(lp+5)=errs(1,i)
            rlist(lp+6)=errs(2,i)
            rlist(lp+4)=errs(3,i)
            rlist(lp+2)=errs(4,i)
          endif
        enddo
      endif
      return
      end
c
      subroutine mfetch(icom,ndepth,latt,twiss,istr,nc,imon,no)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      parameter (kstack=3)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      dimension istr(nstra,4),imon(nmona,4)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      include 'inc/common.inc'
c
      call mfetch1(icom,ndepth,latt,twiss,istr,nc,imon,no,
     $     rlist(iistck(1)),rlist(iistck(2)),rlist(iistck(3)))
      return
      end
c
      subroutine mfetch1(icom,ndepth,latt,twiss,istr,nc,imon,no,idc,ido,
     $     idm)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      parameter (kstack=3)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      dimension istr(nstra,4),imon(nmona,4)
      dimension idc(*),ido(*),idm(*)
      external pack
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      include 'inc/common.inc'
c
      ip=ipnt(icom)
      if(icom.eq.1) then
c...... copy stack 'C' to current buffer. .....
        lc=idc(ip+1-ndepth)
        nc=ilist(1,lc)
        do i=1,nstra
          istr(istr(i,2),3)=1
        enddo
        do i=1,nc
          istr(ilist(mod(i-1,2)+1,lc+nc+1+(i-1)/2),3)=0
        enddo
        call pack(istr(1,2),istr(1,3),nc,nstra)
        do 14 i=1,nc
          j=istr(istr(i,2),1)
          rlist(latt(j)+11)=rlist(lc+i)
 14     continue
        do i=nc+1,nstra
          rlist(latt(istr(istr(i,2),1))+11)=0d0
        enddo
      elseif(icom.eq.2) then
      elseif(icom.eq.3) then
c...... copy stack 'M' to current buffer. ......
        lp=idm(ip+1-ndepth)
        no=ilist(1,lp)
        do i=1,nmona
          imon(imon(i,2),3)=1
        enddo
        do i=1,no
          imon(ilist(mod(i-1,2)+1,lp+2*no+1+(i-1)/2),3)=0
        enddo
        call pack(imon(1,2),imon(1,3),no,nmona)
        do 10 i=1,no
          j=imon(imon(i,2),1)
          twiss(j,0,15)=rlist(lp+i)
          twiss(j,0,17)=rlist(lp+no+i)
 10     continue
      endif
      return
      end
c
      subroutine pasexg(latt,pos,mult,title,case,x,ex,isex,nsex,vscale,
     $     hpos,vonly,lio)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical hpos,vonly
      character*(*) title,case
c      character*8 name
      character*80 line
      character*30 dat
      integer*8 latt(nlat)
      dimension pos(nlat),mult(*)
      dimension x(2*nsex),ex(2*nsex),isex(nsex),vscale(2)
      data line/' '/
c
      call fdate1(dat)
      write(lio,*)
     $     'NEWFRANE;SET FONT DUPLEX;SET CARD 80;SET TITLE SIZE -3'
      write(lio,*)'SET WINDOW X 3.3 12.1 Y 2.5 8.35'
      write(lio,*)'TITLE 7.5 9.1 CENTER '' '''
      write(lio,*)'MORE  ''',title(1:min(59,lene(title))),''''
      write(lio,*)'CASE '' ',case(1:min(59,lene(case))) ,''''
      write(lio,*)'SET TITLE SIZE -1.6'
      write(lio,*)'TITLE 8 8.5 ''',dat,  ''''
c
      write(lio,*)'SET SYMBOL 9O SIZE .5'
      if(vonly) go to 10
      write(lio,*)'SET WINDOW X 2.5 11.0 Y 5.4 8.3'
      write(lio,*)'TITLE LEFT ''DX (MM)'' SIZE 2'
      write(lio,*)'CASE       ''FL  LL '''
      write(lio,*)'SET LABEL BOTTOM OFF'
      if(hpos) then
        write(lio,*)
     $     'SET LIMITS X ',sngl(pos(1)),' TO ',sngl(pos(nlat))
      else
        write(lio,*)'SET LIMITS X 0 ',nsex+1
      endif
      if(vscale(1).ne.1d20) write(lio,*)
     $     'SET LIMITS Y ',sngl(vscale(2)),' TO ',sngl(vscale(1))
      write(lio,*)'SET ORDER X Y DY'
      do i=1,nsex
        y=x(i+nsex)*1d3
        if(abs(y).lt.1d-20) y=0d0
        if(abs(y).gt.1d20) y=1d20
        dy=ex(i+nsex)*1d3
        if(abs(dy).lt.1d-20) dy=0d0
        if(abs(dy).gt.1d20) dy=1d20
        if(hpos) then
          write(line(1+34*mod(i-1,2):),'(1p,3e11.3,a)')
     $         pos(isex(i)),y,dy,';'
        else
          write(line(1+27*mod(i-1,2):),'(i4,1p,2e11.3,a)') i,y,dy,';'
        endif
        if(mod(i,2).eq.0 .or.i.eq.nsex) then
          write(lio,'(a)') line(1:lene(line))
          line=' '
        endif
      enddo
      write(lio,*)'PLOT'
 10   write(lio,*)'SET WINDOW X 2.5 11.0 Y 2.5 5.4'
      write(lio,*)'TITLE LEFT ''DY (MM)'' SIZE 2'
      write(lio,*)'CASE       ''FL  LL '''
      write(lio,*)'SET LABEL BOTTOM ON'
      if(hpos) then
        write(lio,*)
     $     'SET LIMITS X ',sngl(pos(1)),' TO ',sngl(pos(nlat))
      else
        write(lio,*)'SET LIMITS X 0 ',nsex+1
      endif
      if(vscale(1).ne.1d20) write(lio,*)
     $     'SET LIMITS Y ',sngl(vscale(2)),' TO ',sngl(vscale(1))
      write(lio,*)'SET ORDER X Y DY'
      do i=1,nsex
        y=x(i)*1d3
        if(abs(y).lt.1d-20) y=0d0
        if(abs(y).gt.1d20) y=1d20
        dy=ex(i)*1d3
        if(abs(dy).lt.1d-20) dy=0d0
        if(abs(dy).gt.1d20) dy=1d20
        if(hpos) then
          write(line(1+34*mod(i-1,2):),'(1p,3e11.3,a)')
     $         pos(isex(i)),y,dy,';'
        else
          write(line(1+27*mod(i-1,2):),'(i4,1p,2e11.3,a)') i,y,dy,';'
        endif
        if(mod(i,2).eq.0 .or.i.eq.nsex) then
          write(lio,'(a)') line(1:lene(line))
          line=' '
        endif
      enddo
      write(lio,*)'PLOT'
c     write(lio,*)'SET WINDOW X 2.5 11.0 Y 1.30 2.2'
c     write(lio,*)'SET LIMITS X 0 ',nsex+1,' Y -1 1'
c     write(lio,*)
c    $     'SET OUTLINE ALL OFF;SET TICK ALL OFF;SET LABEL ALL OFF'
c     write(lio,*)'SET TITLE SIZE -.5'
c     do i=1,nsex
c       call elnameK(isex(i),name)
c       write(lio,*)'TITLE ',i,sngl(1.5-0.7*mod(i,2)),
c    $       ' XDATA ANGLE 90 ''',name,''''
c     enddo
      return
      end

      subroutine  pbset(latt,mult,x,isb,xs,
     z                  istr,estr,nstr,ncb,nbump,disp,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      logical disp
      integer*8 latt(nlat)
      dimensionmult(*),x(2*nbump),isb(ncb+4,nbump),xs(ncb+2,2,nbump)
     z,         istr(*),estr(*)
      iw=italoc(nstr)
      call pclr(rlist(iw),nstr)
      do 11 i=1,nbump
        i1=2*i-1
        do 10 j=1,isb(2,i)
          rlist(iw-1+isb(j+2,i))=rlist(iw-1+isb(j+2,i))
     1                  +xs(j,1,i)*x(i1)+xs(j,2,i)*x(i1+1)
   10   continue
   11 continue
      call pcset(latt,mult,rlist(iw),istr,estr,nstr,1d0,disp,lfno)
      call tfree(int8(iw))
      return
      end

      subroutine pbump(latt,twiss,gammab,ip,mult,istr,estr,nstr,kfit,
     1                 ifitp,mfitp,fitval,nfc,iter,lin,disp,lfno)
      use tfstk
      use tmacro
      implicit real*8(a-h,o-z)
      logical iter,lin,disp
      integer*8 latt(*)
      dimension twiss(*),gammab(*),mult(*),istr(*),estr(*),
     1          kfit(*),ifitp(*),mfitp(*),fitval(*)
      if(nstr.eq.0 .or. nfc.eq.0)then
        return
      endif
      iqu=italoc(nfc*nstr)
      idf=italoc(nfc)
      idv=italoc(nstr)
      iwk=italoc(nfc)
      call pbump1(latt,twiss,gammab,ip,mult,istr,estr,nstr,kfit,ifitp,
     1            mfitp,fitval,nfc,rlist(iqu),rlist(idf),rlist(idv),
     z            rlist(iwk),iter,lin,disp,lfno)
      call tfree(int8(iwk))
      call tfree(int8(idv))
      call tfree(int8(idf))
      call tfree(int8(iqu))
      return
      end

      subroutine pbump1(latt,twiss,gammab,ip,mult,istr,estr,nstr,kfit,
     1                  ifitp,mfitp,fitval,nfc,qu,df,dval,iwk,iter,lin,
     1                  disp,lfno)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
      implicit real*8 (a-h,o-z)
      parameter (mfitc1=32,mfitc2=28,epsi=1d-10,epse=1d-10,loopmx=19)
      parameter (ddp=1.d-6)
      logical idealz,stab,iter,lin,conv,disp
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat),
     1          mult(*),istr(nstra,4),estr(nstra),qu(nfc,nstr),
     1          dval(nstr),df(nfc),iwk(2,nfc),wrk(4),imw(1,4)
      dimension kfit(nfc),ifitp(nfc),mfitp(nfc),fitval(nfc)
      include 'inc/common.inc'
c
c     PBUMP makes a bump orbit on the 'backuped' orbit which is made
c     by PCBAK. It calculates corrector strength in two satges:
c     In the first stage, calculation is based on the 'model'
c     optics (twiss(*,-ndim,*)). Exception is for lin=T, in that
c     case caluculation is based on the 'fixed' optics(twiss(*,ndim,*)).
c     If lin=F, the process of the first stage is iterated.
c     In the second stage, calculated strengths are added to the
c     'backuped' correctors and calculate the orbit in the current 
c     optics buffer(twiss(*,0,*)). If iter=T, the process of the second
c     satge is iterated.
c
      idealz=ideal
      itgt=0
      iqu=0
c
      conv=.false.
      ideal=.true.
      do 22 i=1,nlat-1
        if(idtypec(i).eq.icbend) then
          rlist(latt(i)+11)=0d0
        endif
   22 continue
      idd=ip
      if(lin) then
c       idd=ip
      else
        idd=-ndim
        do 23 i=1,ntwissfun
          twiss(1,-ndim,i)=twiss(1,ip,i)
   23   continue
        call pqcell(latt,twiss,gammab,-ndim,dp0+1d0,stab)
      endif
      psix=twiss(nlat,idd,3)-twiss(1,idd,3)
      psiy=twiss(nlat,idd,6)-twiss(1,idd,6)
c Oide modified following four lines. 4/28/1993
c This change looks necessary to make a bump when nmona is zero.
c Also see the change around line label 999.
      nmonaa=nmona
      nmona=1
      imw(1,2)=1
      do 100 loop=0,loopmx
        nc=0
        do 40 k=1,nfc
          if(mfitp(k).ne.0) then
            imw(1,1)=ifitp(k)
            if(kfit(k) .eq. mfitc1) then
c             -- 'DDX' ---
              call mcrmat(latt,twiss,gammab,idd,psix,psiy,dval,1,
     1                    .false.,.false.,istr,nstr,imw,1,'X')
            elseif(kfit(k) .eq. mfitc1+1) then
c             -- 'DDPX' ---
              call mcrmat(latt,twiss,gammab,idd,psix,psiy,dval,1,
     z                    .true.,.false.,istr,nstr,imw,1,'X')
            elseif(kfit(k) .eq. mfitc1+2) then
c             -- 'DDY' ---
              call mcrmat(latt,twiss,gammab,idd,psix,psiy,dval,1,
     z                    .false.,.false.,istr,nstr,imw,1,'Y')
            elseif(kfit(k) .eq. mfitc1+3) then
c             -- 'DDPY' ---
              call mcrmat(latt,twiss,gammab,idd,psix,psiy,dval,1,
     z                    .true.,.false.,istr,nstr,imw,1,'Y')
            elseif(kfit(k) .eq. mfitc2) then
c             -- 'DEX' ---
              call pcrmat(latt,twiss,gammab,psix,psiy,dval,1,istr,nstr,
     1                    imw,1,.false.,.true.,'X',wrk,ddp)
            elseif(kfit(k) .eq. mfitc2+1) then
c             -- 'DEPX' ---
              call pcrmat(latt,twiss,gammab,psix,psiy,dval,1,istr,nstr,
     1                    imw,1,.true.,.true.,'X',wrk,ddp)
            elseif(kfit(k) .eq. mfitc2+2) then
c             -- 'DEY' ---
              call pcrmat(latt,twiss,gammab,psix,psiy,dval,1,istr,nstr,
     &                    imw,1,.false.,.true.,'Y',wrk,ddp)
            elseif(kfit(k) .eq. mfitc2+3) then
c             -- 'DEPY' ---
              call pcrmat(latt,twiss,gammab,psix,psiy,dval,1,istr,nstr,
     &                    imw,1,.true.,.true.,'Y',wrk,ddp)
            else
              goto 40
            endif
            nc=nc+1
            if(loop.eq.0) then
              df(nc)=fitval(k)
              if(kfit(k).ge.mfitc2.and.kfit(k).le.mfitc2+3) then
                iwk(1,nc)=kfit(k)-mfitc2+7
              else
                iwk(1,nc)=kfit(k)-mfitc1+15
              endif
              iwk(2,nc)=ifitp(k)
c             print *,iwk(1,nc),iwk(2,nc),sngl(fitval(k))
            endif
c           write(*,'(i3,10d10.3)') kfit(k),dval
            do 30 i=1,nstr
              qu(nc,i)=dval(i)
   30       continue
          endif
   40   continue
        if(conv) then
          iqu=italoc(nfc*nstr)
          call tmov(qu,rlist(iqu),nfc*nstr)
          goto 101
        endif
        if(loop.eq.0) then
          if(nc.eq.0) goto 999
          itgt=italoc(nc)
          do 50 i=1,nc
            rlist(itgt-1+i)=twiss(iwk(2,i),-ndim,iwk(1,i))+df(i)
   50     continue
          idd=-ndim
          psix=twiss(nlat,-ndim,3)-twiss(1,-ndim,3)
          psiy=twiss(nlat,-ndim,6)-twiss(1,-ndim,6)
        endif
c       write(*,'(a/(1p,5e11.4))')' qu',((qu(i,j),i=1,nc),j=1,nstr)
        call tsolvg(qu,df,dval,nc,nstr,nfc)
c       write(*,'(a/(1p,8e11.4))') ' dval',(dval(i),i=1,nstr)
        do 60 i=1,nstr
          j=istr(i,2)
          rlist(latt(istr(j,1))+11) = rlist(latt(istr(j,1))+11)
     1                               - dval(i)
   60   continue
        if(iter) then
c        if(.not.lin) then
          call pqcell(latt,twiss,gammab,-ndim,dp0+1d0,stab)
          dmax=-1d31
          do 70 i=1,nc
            df(i)=rlist(itgt-1+i)-twiss(iwk(2,i),-ndim,iwk(1,i))
            if(iwk(1,i).eq.16.or.iwk(1,i).eq.8) then
              dmax=max(dmax,abs(df(i))*twiss(iwk(2,i),-ndim,2))
            elseif(iwk(1,i).eq.18.or.iwk(1,i).eq.10) then
              dmax=max(dmax,abs(df(i))*twiss(iwk(2,i),-ndim,5))
            else
              dmax=max(dmax,df(i),-df(i))
            endif
   70     continue
c         write(*,'(a/(1p,6e11.4))') ' df',(df(i),i=1,nc)
c         print *,' loop=',loop,' dmax=',sngl(dmax),' (model)'
          if(dmax.lt.epsi) then
            conv=.true.
          endif
        else
          goto 101
        endif
  100 continue
      call permes(' PBUMP -->',' Iteration not converged (model)',' ',
     &            lfno)
  101 continue
      nmona=nmonaa
      do 102 i=1,nstr
        j=istr(i,2)
        dval(i)=rlist(latt(istr(j,1))+11)
  102 continue
      if(idealz) then
        ideal=.true.
        call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
        goto 999
      endif
      ideal=.false.
      do 103 i=1,nlat-1
        if(idtypec(i).eq.icbend) then
c          write(*,*)i,rlist(latt(i)+11),rlist(ibckup-1+i)
          rlist(latt(i)+11)=rlist(ibckup-1+i)
        endif
  103 continue
      do 104 i=1,18
        twiss(1,0,i)=optiv(i)
  104 continue
      if(iter .and. .not.lin .and. simulate) then
        nc=0
        do 105 k=1,nfc
          if(mfitp(k).ne.0) then
            if(kfit(k).ge.mfitc1.and.kfit(k).le.mfitc1+3 .or.
     1         kfit(k).ge.mfitc2.and.kfit(k).le.mfitc2+3) then
              nc=nc+1
c             df(nc)=fitval(k)
              rlist(itgt-1+nc)=twiss(iwk(2,nc),0,iwk(1,nc)) + fitval(k)
            endif
          endif
  105   continue
        call pcset(latt,mult,dval,istr,estr,nstr,1d0,.false.,lfno)
        call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
      else
        call pcset(latt,mult,dval,istr,estr,nstr,1d0,disp,lfno)
c        do  i=1,nlat-1
c          if(idtypec(i).eq.icbend) then
c            if(abs(rlist(latt(i)+11)).gt.1d-10) then
c              print *,'pbump1',i,rlist(latt(i)+11),rlist(ibckup-1+i)
c            endif
c          endif
c        enddo
        if(simulate) then
          call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
        endif
        goto 999
      endif
      do 200 loop=0,loopmx
        dmax=-1d31
        do 110 i=1,nc
          df(i)=rlist(itgt-1+i)-twiss(iwk(2,i),0,iwk(1,i))
          if(iwk(1,i).eq.16.or.iwk(1,i).eq.8) then
            dmax=max(dmax,abs(df(i))*twiss(iwk(2,i),-ndim,2))
          elseif(iwk(1,i).eq.18.or.iwk(1,i).eq.10) then
            dmax=max(dmax,abs(df(i))*twiss(iwk(2,i),-ndim,5))
          else
            dmax=max(dmax,df(i),-df(i))
          endif
  110   continue
c       write(*,'(a/(1p,6e11.4))') ' df',(df(i),i=1,nc)
        print *,' loop=',loop,' dmax=',sngl(dmax)
        if(dmax.lt.epse) then
          goto 201
        endif
        call tmov(rlist(iqu),qu,nfc*nstr)
c       write(*,'(a/(1p,5e11.4))')' qu',((qu(i,j),i=1,nc),j=1,nstr)
        call tsolvg(qu,df,dval,nc,nstr,nfc)
        call pcset(latt,mult,dval,istr,estr,nstr,-1d0,.false.,lfno)
        call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
  200 continue
      call permes(' PBUMP -->',' Iteration not converged',' ',lfno)
      goto 999
  201 continue
      do 210 i=1,nstr
        j=istr(i,2)
        dval(i)=rlist(latt(istr(j,1))+11)-rlist(ibckup-1+istr(j,1))
        rlist(latt(istr(j,1))+11)=rlist(ibckup-1+istr(j,1))
  210 continue
      call pcset(latt,mult,dval,istr,estr,nstr,1d0,disp,lfno)
      call pqcell(latt,twiss,gammab,0,dp0+1d0,stab)
c Oide modified following one line. 4/28/1993
  999 nmona=nmonaa
      if(iqu.ne.0)  call tfree(int8(iqu))
      if(itgt.ne.0) call tfree(int8(itgt))
      ideal=idealz
      return
      end

      subroutine pbumps(istr,nstr,kfit,ifitp,mfitp,nfc,packP)
      use ffs
      use tffitcode
      parameter (mfitc1=32,mfitc2=28)
      logical packP,mhogal
      dimension istr(nstra,4)
      dimension kfit(nfc),ifitp(nfc),mfitp(nfc)
      external pack
      include 'inc/common.inc'
c
      lmin=nlat
      lmax=1
      if(.not.cell) then
        if(packP) then
          call pmovi(istr(1,3),istr(1,4),nstra)
          do 10 k=1,nfc
            if(mfitp(k).ne.0) then
              if(kfit(k).ge.mfitc1.and.kfit(k).le.mfitc1+3 .or.
     1           kfit(k).ge.mfitc2.and.kfit(k).le.mfitc2+3) then
                lmin=min(lmin,ifitp(k))
                lmax=max(lmax,ifitp(k))
              endif
            endif
   10     continue
c         do 11 k=1,nfc
c           if(mfitp(k).ne.0) then
c             if(kfit(k).ge.mfitc1.and.kfit(k).le.mfitc1+3 .or.
c    1           kfit(k).ge.mfitc2.and.kfit(k).le.mfitc2+3) then
c               if(ifitp(k).eq.lmin) then
c                 if(fitval(k).ne.0d0) then
c                   lmin=1
c                   exit
c                 endif
c               endif
c             endif
c           endif
c11       enddo
          do 13 i=1,nstr
            j=istr(i,2)
            if(mhogal(lmin,lmax,istr(j,1))) then
              istr(j,3)=0
c Oide modified as above. 4/28/1993
c              istr(j,3)=min(0,istr(j,3)-1)
            else
              istr(j,3)=1
c Oide modified as above. 4/28/1993
c              istr(j,3)=istr(j,3)+1
            endif
 13       continue
          call pack(istr(1,2),istr(1,3),nstr,nstr)
        else
          call pmovi(istr(1,4),istr(1,3),nstra)
          call pack(istr(1,2),istr(1,3),nstr,nstra)
        endif
      endif
      return
      end

      subroutine  pcbak(latt,twiss)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
c ----- Save current steering & initial values of optics
c       If operate, save dx dy ex ey on twiss(*,ndim-2,*) -----
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      include 'inc/common.inc'
c     print *,'ibckup=',ibckup,ilist(2,latt(1))
      if(ibckup.eq.0) then
        ibckup=italoc(nlat)
      endif
      do 100 i=1,nlat-1
        if(idtypec(i).eq.icbend) then
          rlist(ibckup-1+i)=rlist(latt(i)+11)
        endif
  100 continue
      do 130 i=1,18
  130   optiv(i)=twiss(1,0,i)
      if(.not.simulate) then
        call tmov(twiss(1,0,15),twiss(1,ndim-2,15),nlat)
        call tmov(twiss(1,0,17),twiss(1,ndim-2,17),nlat)
        call tmov(twiss(1,0,7),twiss(1,ndim-2,7),nlat)
        call tmov(twiss(1,0,9),twiss(1,ndim-2,9),nlat)
      endif
      return
      end

      subroutine pclr(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
c     if(n.ge.16) then
        do 10 i=1,n
          a(i)=0d0
   10   continue
c     else
cVOPTION NOVEC
c       do 20 i=1,n
c         a(i)=0d0
c  20   continue
c     endif
      return
      end

      subroutine pclri(ia,n)
      dimension ia(n)
      if(n.ge.16) then
        do 10 i=1,n
          ia(i)=0
   10   continue
      else
*VOPTION NOVEC
        do 20 i=1,n
          ia(i)=0
   20   continue
      endif
      return
      end

      subroutine  pcrmat(latt,twiss,gammab,psix,psiy,
     1                   a,ia,istr,nstr,imon,nmon,prime,normal,xy,
     1                   wrk,ddp)
      use ffs
      use tffitcode
      logical prime,normal,dcon
      character*(*) xy
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat)
      dimension istr(nstra,4),imon(*)
      dimension a(ia,nstr),wrk(2*nmon,2),dcon(2)
      include 'inc/common.inc'
c
c     c1= 1d0/ddp -0.5d0
c     c2= 1d0/ddp +0.5d0
      c1= 1d0/ddp/(1d0+0.5d0*ddp)
      c2= 1d0/ddp/(1d0-0.5d0*ddp)
      jj=istr(1,2)
      if(xy.eq.'XY') then
        nm=2*nmon
      else
        nm=nmon
      endif
      if(normal) then
        if(xy.eq.'Y') then
          dcon(1)=.false.
          dcon(2)=.true.
        elseif(xy.eq.'X') then
          dcon(1)=.true.
          dcon(2)=.false.
        elseif(xy.eq.'XY') then
          dcon(1)=.true.
          dcon(2)=.true.
        endif
        do 110 j=1,nstr
          istr(1,2)=istr(j,2)
c         ----- sgnl kick response of x & y
          call mcrmat(latt,twiss,gammab,ndim,psix,psiy,wrk,2*nmon,
     1                .false.,normal,istr,1,imon,nmon,'XY')
c         ----- sgnl kick response of x' & y'
          call mcrmat(latt,twiss,gammab,ndim,psix,psiy,wrk(1,2),2*nmon,
     1                .true.,normal,istr,1,imon,nmon,'XY')
c         ----- calc dR/dp*R^(-1)*(x x' y y')^t
          call mcrmda(twiss,wrk,wrk(1,2),a(1,j),ia,1,imon,nmon,prime,
     1                dcon)
          call mcrmat(latt,twiss,gammab,ndim-1,0d0,0d0,wrk,ia,
     &         prime,normal,istr,1,imon,nmon,xy)
          call mcrmat(latt,twiss,gammab,1-ndim,0d0,0d0,wrk(1,2),ia,
     &         prime,normal,istr,1,imon,nmon,xy)
          call mclin(wrk,wrk(1,2),c1,-c2,nmon,5)
c         ----- add correction -dR/dp*R^(-1)*(x)
          call mclin(wrk(1,2),a(1,j),1d0,-1d0/ddp,nm,1)
110     continue
      else
        do 120 j=1,nstr
          istr(1,2)=istr(j,2)
          call mcrmat(latt,twiss,gammab,ndim-1,0d0,0d0,wrk,ia,prime,
     1                normal,istr,1,imon,nmon,xy)
          call mcrmat(latt,twiss,gammab,1-ndim,0d0,0d0,a(1,j),ia,prime,
     1                normal,istr,1,imon,nmon,xy)
          call mclin(wrk,a(1,j),c1,-c2,nm,5)
120     continue
      endif
      istr(1,2)=jj
      return
      end

      subroutine  pcset(latt,mult,x,istr,estr,nstr,pm,disp,lfno)
      use tfstk
      use ffs
      use tffitcode
      parameter (bigg=pi*0.5d0)
      logical disp
      character*8 vout*80,autofg,name
      integer*8 latt(nlat)
      dimension mult(*),
     &          istr(nstra,4),estr(nstra),
     &          x(nstr),rc(4)
      data vout/' kick calc (max/@/rms/ave)'/
      include 'inc/common.inc'
c
      it=italoc(nstr)
      big=-1d29
      if(pm.eq.-1d0) then
        do 20 i=1,nstr
          j=istr(i,2)
          rlist(it-1+i)=rlist(latt(istr(j,1))+11)-x(i)*(1d0+estr(j))
          big=max(big,rlist(it-1+i),-rlist(it-1+i))
   20   continue
      elseif(pm.eq.1d0) then
        do 21 i=1,nstr
          j=istr(i,2)
          rlist(it-1+i)=rlist(latt(istr(j,1))+11)+x(i)*(1d0+estr(j))
          big=max(big,rlist(it-1+i),-rlist(it-1+i))
   21   continue
      else
        do 22 i=1,nstr
          j=istr(i,2)
          rlist(it-1+i)=rlist(latt(istr(j,1))+11)
     1                 + pm*x(i)*(1d0+estr(j))
          big=max(big,rlist(it-1+i),-rlist(it-1+i))
   22   continue
      endif
      if(big.gt.bigg) then
        call mstatp(rlist(it),nstr,rc(1),tmin,rc(2),rc(3),imax)
        rc(1)=max(rc(1),-tmin)
        call permes(' !!!',' Too large kick angle: --> canceled.',' ',
     &              lfno)
        if(disp) then
          call elname(istr(istr(imax,2),1),name)
          vout(27:34)=autofg(rc(1)*1d3,'8.5')
          vout(36:43)=name
          vout(43:50)=autofg(rc(2)*1d3,'8.5')
          vout(51:58)=autofg(rc(3)*1d3,'8.5')
          write(lfno,'(a)') vout
        endif
      else
        do 30 i=1,nstr
          j=istr(i,2)
   30     rlist(latt(istr(j,1))+11)= rlist(it-1+i)
      endif
      call tfree(int8(it))
      return
      end

      subroutine pecorg(ft,p,nbin,dds,lab,lca,pab,pca,title,case,lfno)
      implicit real*8 (a-h,o-z)
      complex*16 ft(nbin)
      character*(*) lab,lca,pab,pca,title,case
cslac character*40 dat
c     character*30 dat
c     call fdate1(dat)
      write(lfno,*)'NEWFRAME;SET FONT DUPLEX;SET TITLE SIZE -3.5'
      write(lfno,*)'SET WINDOW X 3.3 12.1 Y 2.5  8.35'
      write(lfno,*)
     1            'TITLE TOP ''',title(1:min(59,lene(title))),''''
      write(lfno,*)'CASE      ''',case(1:min(59,lene(case))) ,''''
      write(lfno,*)'SET TITLE SIZE -1.6'
cslac write(lfno,*)'TITLE 7 8.5 ''',dat,''''
c     write(lfno,*)'TITLE 8 8.5 ''',dat,  ''''
      write(lfno,*)'SET TITLE SIZE -3'
      write(lfno,*)'TITLE LEFT ''',lab(1:lene(lab)),''''
      write(lfno,*)'CASE       ''',lca(1:lene(lca)),''''
      write(lfno,*)'TITLE BOTTOM ''',pab(1:lene(pab)),''''
      write(lfno,*)'CASE         ''',pca(1:lene(pca)),''''
      write(lfno,*)'SET LIMIT X ',sngl(-dds*float(nbin/2)),
     1                            sngl(dds*float(nbin/2-1))
      call tdinit(lfno,'JOIN 1','10 ')
      do 10 i=nbin/2+1,nbin
        call tdput(float(i-1-nbin)*dds,dble(ft(i))/p)
10    continue
      do 20 i=1,nbin/2
        call tdput(float(i-1)*dds,dble(ft(i))/p)
20    continue
      call tdterm
      return
      end

      subroutine pecorr(word,latt,pos,errk,title,case,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      parameter (ier=6)
      integer*8 latt(nlat)
      dimension pos(nlat),errk(2,nlat)
      character*(*) word,title,case
      logical exist
      dds=getva(exist)
      if(.not.exist) then
        call getwdl(word)
        return
      endif
      call getwdl(word)
      nbin=int(1d0+log((pos(nlat)-pos(1))/dds)/log(2d0))
      nbin=2**nbin
      ift=italoc(nbin*2*ier)
      ifl=italoc((ier*nbin+1)/2)
      call pecorr1(word,latt,pos,errk,rlist(ift),rlist(ifl),nbin,
     1             title,case,lfno,exist)
      call tfree(int8(ifl))
      call tfree(int8(ift))
      if(exist) call getwdl(word)
      return
      end

      subroutine pecorr1(word,latt,pos,errk,ft,ift,nbin,title,case,
     1                   lfno,exist)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit real*8 (a-h,o-z)
      parameter (ier=6)
      integer*8 jd
      complex*16 ft(nbin,ier)
      integer*8 latt(nlat),k
      dimension pos(nlat),errk(2,nlat),ift(nbin,ier),err(ier)
      equivalence (dl,err(1)),(dk,err(2)),(ddk,err(3)),(dtheta,err(4)),
     1            (dx,err(5)),(dy,err(6))
      character*(*) word,title,case
      character*20  lab(ier),lca(ier),pab(ier),pca(ier)
c     character*20  name
      logical tmatch,exist
      data (lab(i),i=1,2)/'<DL(S)DL(S+L)>0S1   ','<DK(S)DK(S+L)>0S1   '/
      data (lca(i),i=1,2)/' GL L GL L L  XLX   ',' GL L GL L L  XLX   '/
      data (lab(i),i=3,4)/'<DDK(S)DDK(S+L)>0S1 ','<DQ(S)DQ(S+L)>0S1   '/
      data (lca(i),i=3,4)/' GLL L GLL L L  XLX ',' GG L GG L L  XLX   '/
      data (lab(i),i=5,6)/'<DX(S)DX(S+L)>0S1   ','<DY(S)DY(S+L)>0S1   '/
      data (lca(i),i=5,6)/' GL L GL L L  XLX   ',' GL L GL L L  XLX   '/
      data (pab(i),i=1,2)/'VDL(n)V223          ','VDK(N)V223          '/
      data (pca(i),i=1,2)/'SGL L SX X          ','VGL L SX X          '/
      data (pab(i),i=3,4)/'VDDK(N)V223         ','VDQ(N)V223          '/
      data (pca(i),i=3,4)/'SGLL L SX X         ','SGG L SX X          '/
      data (pab(i),i=5,6)/'VDX(N)V223          ','VDY(N)V223          '/
      data (pca(i),i=5,6)/'SGL L SX X          ','SGL L SX X          '/
      exist=.false.
      dds=(pos(nlat)-pos(1))/float(nbin)
      print *,'nbin=',nbin,' dds=',sngl(dds)
        call pclr(ft,2*ier*nbin)
        call pclri(ift,ier*nbin)
        do 2010 i=1,nlat-1
          if( .not.tmatch(pnamec(i),word) )goto 2010
          exist=.true.
          j=ilist(2,latt(i))
          k=latt(i)
          id=idtype(j)
          jd=idval(j)
          go to (2110,2120,2010,2140,2010,2140,2010,2140),id
          go to 2010
2110      dl=rlist(k+1)-rlist(jd+1)
          dk=0.d0
          ddk=0.d0
          dtheta=0.d0
          dx=0.d0
          dy=0.d0
          go to 2200
2120      dl=rlist(k+1)-rlist(jd+1)
          if(rlist(jd+2) .ne. 0.d0)then
            dk=rlist(k+2)/rlist(jd+2)-1.d0
          else
            dk=errk(1,i)-1.d0
          endif
          ddk=rlist(k+11)-rlist(jd+11)
          dtheta=rlist(k+5)-rlist(jd+5)
          dx=rlist(k+9)
          dy=rlist(k+10)
          go to 2200
2140      dl=rlist(k+1)-rlist(jd+1)
          if(rlist(jd+2) .ne. 0.d0)then
            dk=rlist(k+2)/rlist(jd+2)-1.d0
          else
            dk=errk(1,i)-1.d0
          endif
          ddk=errk(2,i)-1.d0
          dtheta=rlist(k+4)-rlist(jd+4)
          dx=rlist(k+5)
          dy=rlist(k+6)
2200      continue
          do 2210 lk=1,6
            ip=min(int(pos(i)/dds)+1,nbin)
            ft(ip,lk)=ft(ip,lk)+err(lk)
            ift(ip,lk)=ift(ip,lk)+1
2210      continue
2010    continue
      do 2300 lk=ier,1,-1
        do 2310 i=1,nbin
          if(ift(i,lk).ne.0d0) then
            ft(i,lk)=ft(i,lk)/float(ift(i,lk))
          endif
2310    continue
        p=0d0
        do 2320 i=1,nbin
          p=p+dble(ft(i,lk))**2
2320    continue
        if(p.gt.1d-18)then
          call tcftr(ft(1,lk),nbin,1)
          do 2330 i=1,nbin
            ft(i,lk)=ft(i,lk)*conjg(ft(i,lk))
2330      continue
          p=p*nbin
          call pecorg(ft(1,lk),p,nbin,1d0,pab(lk),pca(lk),
     1                'MODE NUMBER N',' LLL  LLLLL L',title,case,lfno)
          call tcftr(ft(1,lk),nbin,-1)
c         ff=1d0/dble(ft(1,lk))
c         do 2340 i=1,nbin
c           ft(i,lk)=ft(i,lk)*ff
c2340     continue
c         print *,'rms=',sngl(sqrt(dble(ft(1,lk))/nbin))
          call pecorg(ft(1,lk),dble(ft(1,lk)),nbin,dds,lab(lk),lca(lk),
     $         'L(M)','L L ',title,case,lfno)
        endif
2300  continue
      return
      end

      subroutine permes(pre,mes,post,lfno)
      implicit real*8(a-h,o-z)
      parameter (iterm=6)
      character*(*) pre,mes,post
      l=lfno
 1    write(l,'(A,'' '',A,'' '',A)') pre(1:lene(pre)),mes(1:lene(mes)),
     &                      post(1:lene(post))
      if(lfno.ne.iterm) then
        l=iterm
        goto 1
      endif
      return
      end

      subroutine petcod(word,twiss,imon,nmon)
c      Read/write Tristan MR cod data in KAMADA format.
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),imon(nmona,4)
      dimension cod1(400,2),is1(400,2),rms(2)
      external pack
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
      logical exist
      character word1*8
      character*(*) word
      word1=word
      call getwdl(word)
      if(word.eq.' ') return
      if(word1.eq.'PETIN') then
        call petget(word,cod1,is1)
        cut=getva(exist)
        if(.not.exist) then
          cut=1d3
        endif
        call petfil(cod1,is1,rms,cut)
        do 110 j=1,400
          if(j.gt.nmon) goto 100
          if(is1(j,1).ne.0.or.is1(j,2).ne.0)then
            k=imon(j,2)
            imon(k,3)=imon(k,3)+1
            twiss(imon(k,1),0,15)=0d0
            twiss(imon(k,1),0,17)=0d0
          else
            k=imon(j,2)
            imon(k,3)=min(0,imon(k,3)-1)
            twiss(imon(k,1),0,15)=cod1(j,1)*1d-3
            twiss(imon(k,1),0,17)=cod1(j,2)*1d-3
          endif
110     continue
100     call pack(imon(1,2),imon(1,3),nmon,nmona)
        write(*,'(2(A,I4))')'  BPM available =',nmon,' in ',nmonact
        nmonact=nmon
      elseif(word1.eq.'PETOUT') then
        call pclr(cod1,800)
        do 120 j=1,nmon
          k=imon(j,2)
          cod1(j,1)=twiss(imon(k,1),0,15)
          cod1(j,2)=twiss(imon(k,1),0,17)
120     continue
        call pclri(is1,800)
        call petput(word,cod1,is1)
      endif
      call getwdl(word)
      return
      end

      subroutine petfil(cod1,is1,rms,cut)
      implicit real*8 (a-h,o-z)
      DIMENSION COD1(400,2),is1(400,2),rms(2),sig(2),inum(2)
c
      do 4  i=1,2
        sig(i)=0d0
        inum(i)=0
        do 5  j=1,400
          if(is1(j,i).eq.0) then
            inum(i)=inum(i)+1
            sig(i)=sig(i)+cod1(j,i)**2
          endif
    5   continue
        rms(i)=sqrt(sig(i)/real(inum(i)))
c       print *,'Rmsx(mm)=',sngl(rms(1)),' Rmsy(mm)=',sngl(rms(2)),
c    1          'before filtering'
    4 continue
      do 10 i=1,2
        gcut=cut*rms(i)
        inum(i)=0
        sig(i)=0d0
        do 12 j=1,400
          if(abs(cod1(j,i)).gt.gcut) then
            is1(j,i)=1
          else
            inum(i)=inum(i)+1
            sig(i)=sig(i)+cod1(j,i)**2
          endif
   12   continue
        rms(i)=sqrt(sig(i)/real(inum(i)))
   10 continue
      return
      END

      subroutine petget(fnam,cod1,istat)
      implicit real*8 (a-h,o-z)
      character*(*) fnam
C*DEC
      character*256 fnam2
C*DEC
c
      DIMENSION COD1(400,2),S1(400),X1(400),Y1(400),istat(400,2)
      logical lod,lex
      data in/50/
c
      inquire(unit=in,opened=lod,iostat=ios)
c     if(ios .eq. 0) then
c     print *,ios,lod
      if(lod) close(in)
C*DEC
      fnam2 = fnam
      inquire(iostat=ios,file=fnam2,exist=lex)
C*DEC
C*HP
C     inquire(iostat=ios,file=fnam,exist=lex)
C*HP
c     print *,ios,lex
      if(lex) then
        open (in,file=fnam,iostat=ios,status='OLD',
     &   access='SEQUENTIAL')
        if(ios .ne. 0) then
          print *,'Cannot open ',fnam
        else
          READ(in,5000)X1
          READ(in,5000)Y1
          READ(in,5000)S1
 5000     FORMAT (200A4,200A4)
C
          DO 110 J=1,400
            COD1(J,1)=X1(J)
            COD1(J,2)=Y1(J)
  110     CONTINUE
C         reject wrong data
          DO 4  I=1,2
            INUM=0
            DO 5  J=1,400
C             REJECT MONITORS MISSING IN THE SKELETON WITH FLAG 2
              IF(NINT(S1(J)).EQ.2) then
                istat(j,i)=2
                GOTO 5
              endif
              INUM=INUM+1
C             SET FLAG 1 AND COD 0.0 FOR MONITORS WITH
c                                       HARDWARE FAILURE FLAG 1
              IF(NINT(S1(J)).EQ.1) THEN
                istat(J,I)=1
                cod1(J,I)=0.0
              ENDIF
C SET FLAG; 0 H&V NORMAL, 3 REJECT H, 4 REJECT V, 5 REJECT BOTH
              IF(    NINT(S1(J)).EQ.0
     #          .OR. NINT(S1(J)).EQ.3.AND.I.EQ.2
     #          .OR. NINT(S1(J)).EQ.4.AND.I.EQ.1) THEN
                istat(J,I)=0
              ELSE
                istat(J,I)=1
              ENDIF
    5       CONTINUE
c           print *,inum
    4     CONTINUE
          close(in)
        endif
      else
        print *,fnam,'does not exist.'
      endif
      return
      END

      subroutine petput(fnam,cod1,istat)
      implicit real*8 (a-h,o-z)
      character*(*) fnam
c
      real*4 COD1(400,2),S1(400),X1(400),Y1(400),istat(400,2)
      logical lod
      data io/52/
c
      inquire(unit=io,opened=lod,iostat=ios)
      if(lod) close(io)
        open (io,file=fnam,iostat=ios,status='UNKNOWN',
     &   access='SEQUENTIAL',form='FORMATTED')
        if(ios .ne. 0) then
          print *,'Cannot open ',fnam,' ios=',ios
        else
          do 110 j=1,400
            x1(J)=cod1(j,1)
            y1(j)=cod1(j,2)
            s1(j)=istat(j,1)+istat(j,2)
  110     CONTINUE
          write(io,5000)x1
          write(io,5000)y1
          write(io,5000)s1
 5000     FORMAT (2(200A))
          close(io)
        endif
      return
      END

      subroutine petune(word,latt,twiss,mult,gammab,size,nlist,
     $     istr,estr,nstr,title,case,lfno)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      use tfcsi, only:cssetp
      implicit real*8(a-h,o-z)
      parameter (ftol0=1d-6,itmax0=100,demi0=0d0,nrep=10,loop0=1)
      parameter (meminc0=1800,meminc=720,item=20)
      parameter (lioo=50,stndout=6)
      logical exist,abbrev,tmatch,nepas,sequent,history,lod,init,inited,
     $     cont,result,synchb
      character*(*) word,title,case
      character name*12,observ*4
      character*8 nlist(mfit1)
      character*256 tfconvstr
      integer*8 latt(nlat)
      dimension twiss(*),mult(*),gammab(*),size(*),istr(*),
     $     estr(*)
      include 'inc/common.inc'
      common /cordefv/ ipvbmp,nvbmp
c     ..... (iter is used in pamoeb.)
      common /corbtune/iy,ip,iob,iiv,ibckup1,ipemibuf,nemitd,nememax,
     $     iter
      external pmbump
      data itmax/itmax0/,demi/demi0/,loop/loop0/,lio/lioo/
      data observ/'EMIY'/
      data sequent,history,init,inited,cont,result,synchb/7*.false./
      save itmax,demi,loop,lio,
     $     observ,
     $     sequent,history,init,inited,cont,result,synchb
c
      if(iob.ne.0) then
        nobs=ilist(1,iob)
      else
        nobs=0
      endif 
c
 1    call getwdl(word)
      if(word.eq.'INIT') then
        init=.true.
c7/30------>
        if(iob.ne.0) then
          call tfree(int8(iob))
          iob=0
          nobs=0
        endif 
c7/30<------
        go to 1
      elseif(word.eq.'RESULT') then
        result=.true.
        go to 1
      elseif(word.eq.'SYNCHB') then
        synchb=.true.
        go to 1
      elseif(abbrev(word,'N_ITER','_')) then
        itmax=getva(exist)
        if(.not.exist) itmax=itmax0
        goto 1
      elseif(abbrev(word,'D_EMI','_')) then
        demi=getva(exist)
        if(.not.exist) demi=demi0
        goto 1
      elseif(abbrev(word,'S_EQUENTIAL','_')) then
        sequent=.true.
        loop=getva(exist)
        if(.not.exist) loop=loop0
        goto 1
      elseif(word.eq.'EMIY' .or. word.eq.'SIGY' .or. word.eq.'LUMI')
     $       then
        observ=word
        goto 1
      elseif(abbrev(word,'H_ISTORY','_')) then
        itype=itfpeeko(ia,x,next)
        lfni1=x+.5d0
        if(itype.eq.1) then
          call cssetp(next)
          write(word,'(''ftn'',i2.2)') lfni1
        elseif(itype.eq.101) then
          call cssetp(next)
          word=tfconvstr(101,ia,x,nc,'*')
          if(word.eq.' ') then
            call permes('?Missing filename for BTUNE HISTORY','.',' ',
     $           lfno)
            call getwdl(word)
            return
          endif
          call texpfn(word)
        endif
        inquire(unit=lio,opened=lod,iostat=ios)
        if(lod) close(lio)

c       inquire(iostat=ios,file=word,exist=lod)
c       if(.not.lod) then
c         call permes('?File',word,' not found.',lfno)
c         call getwdl(word)
c         return
c       endif
        open (lio,file=word,iostat=ios,status='UNKNOWN',err=991)
        if(ios .ne. 0) then
          call permes('?Cannot open',word,'.',lfno)
          call getwdl(word)
          return
        endif
        history=.true.
        go to 1
 991    call permes('! Error in opening',' ',word,lfno)
        return
      else
        exist=.false.
 10     j=0
        nepas=.false.
        do 12 i=1,nlat-1
          if(tmatch(pnamec(i),word)) then
            do 11 l=1,nobs
              if(ilist(mod(l,2)+1,iob+l/2).eq.i) then
                nepas=.true.
                go to 12
              endif
 11         continue
            j=j+1
            if(exist) then
              ilist(mod(j+nobs,2)+1,iob+(j+nobs)/2)=i
            endif
          endif
 12     continue
        nepas=j.eq.0 .and. nepas
        if(j.eq.0) then
          do 14 i=1,nlat-1
            call elnameK(i,name)
            if(name.eq.word) then
              do 13 l=1,nobs
                if(ilist(mod(l,2)+1,iob+l/2).eq.i) then
                  nepas=.true.
                  go to 15
                endif
 13           continue
              j=1
              if(exist) ilist(mod(j+nobs,2)+1,iob+(j+nobs)/2)=i
              go to 16
            endif
 14       continue
 15       if(nobs.eq.0) then
            return
          elseif(nepas) then
            go to 1
          else
            go to 3
          endif
        endif
 16     exist=.not.exist
        if(exist) then
          if(nobs.eq.0) then
            iob=italoc(j/2+1)
          else
            call palocx(iob,nobs/2+1,(nobs+j)/2+1)
          endif
          ilist(1,iob)=j+nobs
          goto 10
        else
          nobs=ilist(1,iob)
        endif
        goto 1
      endif
c
 3    continue
      if(iob.eq.0) then
        iob=italoc(1)
        ilist(1,iob)=1
        ilist(2,iob)=nlat
      endif
      ftol=max(ftol0,demi)
      nd=nvbmp
      if(init) then
        inited=.false.
        cont=.false.
        result=.false.
      endif
      if(result) then
c7/13   inited=.false.
        init=.false.
        cont=.false.
      endif
      if(inited) then
        if(.not.result) then
          cont=.true.
        endif
      elseif(.not.init .and. .not.result) then
        init=.true.
        cont=.true.
        result=.true.
      endif
c     print *,'init=',init,' inited=',inited,' cont=',cont,' result=',
c    $     result
c*****Initialize********************************************************
      if(init) then
      print *,'INIT=',init,' inited=',inited,' cont=',cont,' result=',
     $     result
        print *,'ibckup1=',ibckup1,' iiv=',iiv,' iy=',iy,' ip=',ip,
     $       ' ipemibuf=',ipemibuf
c.......save corrector to backup storage ....
        call pcbak(latt,twiss)
        if(ibckup1.ne.0) call tfree(int8(ibckup1))
        ibckup1=italoc(nlat)
        call tmov(rlist(ibckup),rlist(ibckup1),nlat)
c.......save initial value.
        if(iiv.ne.0) call tfree(int8(iiv))
        iiv=italoc(nd)
        do 18 l=1,nd
          lp=ilist(mod(l-1,2)+1,ipvbmp+(l-1)/2)
          nc=ilist(1,lp)
          lv=ilist(2,lp)
          lp3=lp+3+2*((nc+1)/2)
          rlist(iiv-1+l)=rlist(lp3-1+lv)
c         print *,iiv-1+l
 18     continue
        if(iy.ne.0) call tfree(int8(iy))
        iy=italoc(nd+1)
        if(ip.ne.0) call tfree(int8(ip))
        ip=italoc(nd*(nd+1))
        iter=0
c7/13<---
        if(ipemibuf.ne.0) then
          call tfree(int8(ipemibuf))
          ipemibuf=0
          nemitd=0
        endif
c7/13--->   
c       write(*,'(1p,10d11.4)')(rlist(iiv-1+i),i=1,nd)
        if(.not.sequent .and. nd.ne.1) then
          if(itmax.ne.0) then
c...........set initial value.
            call petune1(rlist(iy),rlist(ip),rlist(ipvbmp),nd,latt,
     $           twiss,mult,gammab,size,istr,estr,nstr,observ,
     $           rlist(iob),demi,synchb,lfno)
          endif

c     
        elseif(loop.ne.0) then
          call pclr(rlist(ip),nd)
        endif
        inited=init
        init=.false.
      endif
c*****TUNE**************************************************************
      if(cont) then
      print *,'init=',init,' inited=',inited,' CONT=',cont,' result=',
     $     result
        if((sequent.or.nd.eq.1).and.loop.eq.0 .or.
     $     (.not.sequent.and.nd.ne.1).and.itmax.eq.0) then
c.........if loop.eq.0(sequential) or itmax=0(simplex)........
          itemno=ilist(1,iob)*item
          if(nemitd.eq.0) then
            ipemibuf=italoc(meminc0)
            nememax=meminc0
          elseif((nemitd+1)*itemno.gt.nememax) then
            call palocx(ipemibuf,nemitd*itemno,
     $           nemitd*itemno+max(meminc,itemno))
            nememax=nemitd*itemno+max(meminc,itemno)
          endif
          nemitd=nemitd+1
          ymin=pmeas(latt,twiss,gammab,size,observ,rlist(iob),
     $         rlist(ipemibuf+(nemitd-1)*itemno),.true.,synchb,lfno)
          call pmbdata(latt,mult,rlist(iob),ilist(1,iob),observ,nemitd,
     $         0d0,.true.,lfno)
          call pmbdata(latt,mult,rlist(iob),ilist(1,iob),observ,nemitd,
     $         ymin,.false.,lfno)
        elseif(sequent.or.nd.eq.1) then
c.........Sequential.....
          if(itmax.ge.itmax0) itmax=3
          do 17 i=1,loop
            call petunes(latt,twiss,mult,gammab,size,nlist,istr,
     $           estr,nstr,observ,rlist(iob),itmax,ftol,demi,iter,
     $           synchb,lfno)
 17       continue
        else
c.........simplex.....
          ipr=italoc(nd)
          iprr=italoc(nd)
          ipbar=italoc(nd)
          call pamoeb(rlist(ip),
c     ----- pass arguments of pmbump ----
     &         pmbump,rlist(ipvbmp),latt,twiss,mult,gammab,size,istr,
     $         estr,nstr,observ,rlist(iob),synchb,demi,lfno,
c     --------------------
     &         rlist(iy),nd,rlist(ipr),rlist(iprr),rlist(ipbar),ftol,
     $         itmax)
          call tfree(int8(ipbar))
          call tfree(int8(iprr))
          call tfree(int8(ipr))
          print *,'ibckup1=',ibckup1,' iiv=',iiv,' iy=',iy,' ip=',ip,
     $         ' ipemibuf=',ipemibuf
        endif
      endif
c*****Result************************************************************
      if(result) then
        print *,'init=',init,' inited=',inited,' cont=',cont,' RESULT=',
     $       result
        ilo=1
        if(.not.sequent.and.nd.ne.1) then
          if(itmax.ne.0) then
            do 20 i=1,nd+1
              if(rlist(iy-1+i).lt.rlist(iy-1+ilo)) ilo=i
 20         continue
            ymin=pmbump(rlist(ip+nd*(ilo-1)),rlist(ipvbmp),nd,latt,
     $           twiss,mult,gammab,size,istr,estr,nstr,observ,
     $           rlist(iob),.true.,synchb,0d0,lfno)
            call pmbdata(latt,mult,rlist(iob),ilist(1,iob),observ,
     $           nemitd,0d0,.false.,lfno)
          endif
        else
          if(loop.ne.0) then
            ymin=rlist(ip+nd)
          endif
        endif
        write(lfno,'(a,1pd11.4)')' Minimized. emiy/emix =',ymin
        do 22 l=1,nd
          lp=ilist(mod(l-1,2)+1,ipvbmp+(l-1)/2)
          nc=ilist(1,lp)
          lv=ilist(2,lp)
          kv=ilist(1,lp+1)
          lp3=lp+3+2*((nc+1)/2)
          call mcchar(ilist(2,lp+1),name,3)
          write(lfno,'(a,a,1pd11.4)') ' '//name//' ',nlist(kv),
     &         rlist(ip-1+nd*(ilo-1)+l)
c.........recover initial value ...
          rlist(lp3-1+lv)=rlist(iiv-1+l)
 22     continue
c.......write beam parameter ............
        call pmbdata(latt,mult,rlist(iob),ilist(1,iob),observ,nemitd,
     $       0d0,.true.,lfno)
c       write(*,'(1p,10d11.4)')(rlist(iiv-1+i),i=1,nd)
c7/16   if(ip.ne.0) call tfree(int8(ip))
c7/16   if(iy.ne.0) call tfree(int8(iy))
c7/16   if(iiv.ne.0) call tfree(int8(iiv))
c7/16   ip=0
c7/16   iy=0
c7/16   iiv=0
c.......output history......
        if(history) then
c          call pmbdrw(rlist(ipemibuf),rlist(iob),ilist(1,iob),nemitd,
c     $         title,case,lio)
c          if(lio.ne.lfno .and. lio.ne. stndout) close(lio)
        endif
c7/13   if(ipemibuf.ne.0) then
c7/13     call tfree(int8(ipemibuf))
c7/13     ipemibuf=0
c7/13     nemitd=0
c7/13   endif
c.......recover backup data, enabling UNDO. ....     
        call tmov(rlist(ibckup1),rlist(ibckup),nlat)
c7/16   call tfree(int8(ibckup1))
c7/16   call tfree(int8(iob))
c7/16   ibckup1=0
c7/16   iob=0
c.......recover parameters......
c7/13   iter=0
        itmax=itmax0
        demi=demi0
        loop=loop0
        lio=lioo
        sequent=.false.
        history=.false.
        result=.false.
        init=.false.
c7/13   inited=.false.
        cont=.false.
        print *,'ibckup1=',ibckup1,' iiv=',iiv,' iy=',iy,' ip=',ip,
     $       ' ipemibuf=',ipemibuf
      endif
      return
      end

      subroutine petune1(y,p,id,nd,latt,twiss,mult,gammab,size,istr,
     &     estr,nstr,observ,iobs,demi,synchb,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      character*(*) observ
      logical synchb
      integer*8 latt(*)
      dimension y(0:nd),p(nd,0:nd),id(nd)
      dimension twiss(*),mult(*),gammab(*),istr(*),estr(*),
     $     iobs(*),size(*)
c     <----- output status -------------------
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
c     --------------------------------------->
c     set initial values.
      do 10 j=1,nd
        p(j,0)=0d0
 10   continue
      do 12 l=1,nd
        do 11 j=1,nd
          p(j,l)=0d0
 11     continue
        ip=id(l)
        nc=ilist(1,ip)
        lv=ilist(2,ip)
        ip3=ip+3+2*((nc+1)/2)
        p(l,l)=rlist(ip3-1+lv)
c       print *,'ip=',ip,' nc=',nc,' lv=',lv,' p(l,l)=',p(l,l)
 12   continue
c     make mutiple bumps and get resulting emittance.
      do 13 l=0,nd
        y(l)=pmbump(p(1,l),id,nd,latt,twiss,mult,gammab,size,istr,estr,
     &       nstr,observ,iobs,.false.,synchb,demi,lfno)
 13   continue
      ymax=y(0)
      do 14 l=0,nd
        if(y(l).gt.ymax) ymax=y(l)
 14   continue
      do 15 l=0,nd
        call pmbdata(latt,mult,iobs,iobs(1),observ,l+1,
     $       ymax,.false.,lfno)
 15   continue
      return
      end

      subroutine petunes(latt,twiss,mult,gammab,size,nlist,istr,estr,
     &     nstr,observ,iobs,itmax,ftol,demi,ncallt,synchb,lfno)
      use tfstk
      use ffs
      use tffitcode
c     emittance tuning by a single bump.
      implicit real*8(a-h,o-z)
      parameter (norder=2,nrepmax=10,fext=3.d0)
      logical repeat,synchb
      character*(*) observ
      character name*12
      character*8 nlist(mfit1)
      integer*8 latt(*),ipw,iyw
      dimension twiss(*),mult(*),gammab(*),size(*),istr(*),
     $     estr(*),iobs(*)
      dimension a(0:norder),wb(norder+1),wc(norder+1,norder+1)
      common /cordefv/ ipvbmp,nvbmp
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
c
      if(itmax.le.1) return
      ipw=italoc(itmax)
      iyw=italoc(itmax)
      do 100 loop=1,nvbmp
c.......save initial value.
        lp=ilist(mod(loop-1,2)+1,ipvbmp+(loop-1)/2)
        nc=ilist(1,lp)
        lv=ilist(2,lp)
        kv=ilist(1,lp+1)
        lp3=lp+3+2*((nc+1)/2)
        call mcchar(ilist(2,lp+1),name,3)
        a0=rlist(lp3-1+lv)
c
        ncall=0
        ncalla=0
        nrep=0
c.......Repeat from here
 1      continue
        repeat=.false.
c.......eval at itmax points.
        do 10 i=1,itmax
c         rlist(ipw-1+i)=a0-(i-1)*2.d0*a0/(itmax-1)
          l=mod(itmax/2-1+i,itmax)+1
          rlist(ipw-1+i)=a0-(l-1)*2.d0*a0/(itmax-1)
          rlist(iyw-1+i)=pmbump(rlist(ipw-1+i),
     $         ilist(mod(loop-1,2)+1,ipvbmp+(loop-1)/2),1,latt,twiss,
     $         mult,gammab,size,istr,estr,nstr,observ,iobs,.false.,
     $         synchb,demi,lfno)
          ncall=ncall+1
          ncallt=ncallt+1
c         print *,rlist(ipw-1+i),rlist(iyw-1+i)
 10     continue
        ymin=rlist(iyw)
        ymax=rlist(iyw)
        do 11 i=1,itmax
          ymin=min(ymin,rlist(iyw-1+i))
          ymax=max(ymax,rlist(iyw-1+i))
 11     continue
c........record ymax ......         
        do 12 i=1,ncall-ncalla
          call pmbdata(latt,mult,iobs,iobs(1),observ,nemitd-i+1,ymax,
     $         .false.,lfno)
 12     continue
        ncalla=ncall
c..........................
        if((ymax-ymin)/(ymax+ymin).le.ftol) then
          rlist(lp3-1+lv)=a0
          xmin=0d0
          write(lfno,'(a,2(a,1pd11.4),a,i2,a,i4,a,1pd11.4)')
     $         ' '//name//' ',nlist(kv)(1:max(5,lene(nlist(kv)))),xmin,
     $         '(',rlist(ippv-1+loop),') Ncall=',ncall,'(',ncallt,
     $         ') Insignificant change. y=',rlist(iyw-1+(itmax+1)/2)
          go to 100
        endif 
c.......fit norder-th order polinomial
        call plfit(rlist(ipw),rlist(iyw),a,itmax,norder,wc,wb)
c.......from here norder=2 is assumed
        if(a(2).le.0d0) then
          rlist(lp3-1+lv)=a0
          xmin=0d0
          write(lfno,'(a,2(a,1pd11.4),a,i2,a,i4,a,1pd11.4)')
     $         ' '//name//' ',nlist(kv)(1:max(5,lene(nlist(kv)))),xmin,
     $         '(',rlist(ippv-1+loop),') Ncall=',ncall,'(',ncallt,
     $         ') Curvature negative. y=',rlist(iyw-1+(itmax+1)/2)
          go to 100
        else
          xmax=a0*fext
          xmin=-a(1)/a(2)*0.5d0
          if(abs(xmin).gt.abs(xmax)) then
            xmin=sign(xmax,xmin)
            repeat=.true.
          endif
        endif
        rlist(ipw)=xmin
        ymin=pmbump(rlist(ipw),ilist(mod(loop-1,2)+1,ipvbmp+(loop-1)/2),
     $       1,latt,twiss,mult,gammab,size,istr,estr,nstr,observ,iobs,
     $       .true.,synchb,demi,lfno)
        call pmbdata(latt,mult,iobs,iobs(1),observ,nemitd,0d0,.false.,
     $       lfno)
c.......accumulate parameter change and record ymin...............
        rlist(ippv-1+loop)=rlist(ippv-1+loop)+xmin
        rlist(ippv+nvbmp)=ymin
        ncall=ncall+1
        ncallt=ncallt+1
        rlist(lp3-1+lv)=a0
        if(repeat) then
          if(nrep.lt.nrepmax) then
            nrep=nrep+1
            call pmbdata(latt,mult,iobs,iobs(1),observ,nemitd,0d0,
     $           .true.,lfno)
            go to 1
          endif 
        endif 
        write(lfno,'(a,2(a,1pd11.4),a,i2,a,i4,a,1pd11.4)')
     $       ' '//name//' ',nlist(kv)(1:max(5,lene(nlist(kv)))),xmin,
     $       '(',rlist(ippv-1+loop),') Ncall=',ncall,
     $       '(',ncallt,') Minimized. y =',ymin
        call pmbdata(latt,mult,iobs,iobs(1),observ,nemitd,0d0,.true.,
     $       lfno)
c     
 100  continue 
c     if(history) call tdterm
      call tfree(ipw)
      call tfree(iyw)
      return
      end
