      subroutine  msolb(latt,twiss,gammab,mult,isb,x,ncb,istr,estr,wk,
     1                  iwk,yplane)
c     ----- msolb is called by mbstr ------
      use tfstk
      use ffs
      use tffitcode
      include 'inc/TFMACRO.inc'
      parameter (mfitc1=32,unit=1d-4,nfcmax=6)
      logical yplane,simulz
      dimension wk(ncb+2),iwk(ncb+2,2),isb(ncb+4),istr(nstra,4),estr(*),
     1          x(ncb+2,2)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),
     1          mult(*),kfit(nfcmax),ifitp(nfcmax),mfitp(nfcmax),
     1          fitval(nfcmax)
      data lfno/6/
      include 'inc/common.inc'
c
      nt=isb(1)
      ncor=isb(2)
      ne=istr(isb(2+ncor),2)
      ne=istr(ne,1)+1
      nfc=ncor
      do 10 i=1,nfc
   10   mfitp(i)=1
      call pclr(fitval,nfc)
      ifitp(1)=nt
      ifitp(2)=nt
      do 11 i=3,nfc
        ifitp(i)=ne
   11 continue
      if(yplane) then
        kfit(1)=mfitc1+2
        kfit(2)=mfitc1+3
        kfit(3)=mfitc1+2
        kfit(4)=mfitc1+3
        if(ncor.gt.4) then
          kfit(5)=mfitc1
          kfit(6)=mfitc1+1
        endif
      else
        kfit(1)=mfitc1
        kfit(2)=mfitc1+1
        kfit(3)=mfitc1
        kfit(4)=mfitc1+1
        if(ncor.gt.4) then
          kfit(5)=mfitc1+2
          kfit(6)=mfitc1+3
        endif
      endif
      do 20 i=1,ncor
        j=istr(isb(i+2),2)
        wk(i)=rlist(latt(istr(j,1))+11)
        iwk(i,1)=j
        iwk(i,2)=istr(i,2)
   20 continue
      do 21 i=1,ncor
        istr(i,2)=iwk(i,1)
   21 continue
c     If simulate=F and lin=T, pbump does not call PQCELL and remains only
c     correctors to make abump. This fudge is necessary for speed.
      simulz=simulate
      simulate=.false.
c     ---- make cos bump ----
      fitval(1)=unit
      call pbump(latt,twiss,gammab,ndim,mult,istr,estr,ncor,kfit,
     1           ifitp,mfitp,fitval,nfc,.false.,.true.,.false.,lfno)
      do 30 i=1,ncor
   30   x(i,1)=(rlist(latt(istr(iwk(i,1),1))+11)-wk(i))/unit
c     ---- make sin bump ----
      fitval(1)=0d0
      fitval(2)=unit
c     call cputime(ctime0,irtc)
      call pbump(latt,twiss,gammab,ndim,mult,istr,estr,ncor,kfit,
     1           ifitp,mfitp,fitval,nfc,.false.,.true.,.false.,lfno)
c     call cputime(ctime1,irtc)
c     write(*,'(A,F10.6)')' ctime(pbump)=',1d-6*(ctime1-ctime0)
      do 50 i=1,ncor
        x(i,2)=(rlist(latt(istr(iwk(i,1),1))+11)-wk(i))/unit
   50 continue
      do 60 i=1,ncor
        rlist(latt(istr(iwk(i,1),1))+11)=wk(i)
        istr(i,2)=iwk(i,2)
   60 continue
c     write(*,'(1p,6e11.2)')(x(i,1),i=1,isb(2))
c     write(*,'(1p,6e11.2)')(x(i,2),i=1,isb(2))
      simulate=simulz
      return
      end
