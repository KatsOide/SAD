C     This rotine is used in TRCOD functions in FFS ( N.yamamoto)
c     function : description
c     trcod    : calculate cod by tracking and set it in twiss.
c     trcodf   : fixcod, move cod in twiss to basic cod and store it in 
c                a private  area (itw)
c     trcods   : recover basic cod[twiss(*,ndim,*)] from itw.
c     trcodm   : replace basic cod and saved cod at one point[ nlat1 ] 
c                by current cod.
c     
C   29/06/92 207140451  MEMBER NAME  TRCODA   *.FORT     M  E2FORT
      subroutine trcoda(latt,nlat1,name,mul, print,
     &            word,title,case,exist,twiss,xa,ya,xxa,xya,yya,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
ckikuchi ... 1 line modified
      logical print, exist
ckikuchi ... 1 line added
      real*8 twiss(nlat,-ndim:ndim,ntwissfun)
      character*(*) word,title,case
      character*8 name,tname
      dimension latt(2,nlat) , mul(nlat)
      integer latt,mul
      dimension sv(5)
      real*8 sv
      logical*4 ojitter
      integer*4 nl,i,olfno
      integer*8 itw
      data itw /0/
      save itw
c
      olfno=lfno
      lfno= 0
c
c     do 1 i=1,nlat1
c         write(lfno,'(1h ,a,i3,4G15.6)')'trcod0'
c    &    ,0   , twiss(i,0   ,15),
c    &    twiss(i,  0 ,17),twiss(i,0   ,16),twiss(i, 0  ,18)
c1        write(lfno,'(1h ,a,i3,4G15.6)')'trcod0'
c    &    ,ndim, twiss(i,ndim,15),
c    &    twiss(i,ndim,17),twiss(i,ndim,16),twiss(i,ndim,18)
c
      ojitter=jitter
      jitter=.false.
      call tracke(latt,nlat1,sv,np,'STANDBY',name,lfno)
      jitter=ojitter
      nl=1
      call tracke(latt,nl   ,sv,np,'TRACK',name,lfno)
      call tracke(latt,nl   ,sv,np,'POS',name,lfno)
          twiss(1,0,mfitdx) =sv(1)
          twiss(1,0,mfitdpx)=sv(2)
          twiss(1,0,mfitdy) =sv(3)
          twiss(1,0,mfitdpy)=sv(4)
      call tracke(latt,nl   ,sv,np,'DISP',name,lfno)
          if(sv(5) .ne. 0d0) then
            twiss(1,0, 7)=sv(1)
            twiss(1,0, 8)=sv(2)
            twiss(1,0, 9)=sv(3)
            twiss(1,0,10)=sv(4)
          end if
c     write(lfno,*)'trcod(start)',(sv(j),j=1,4)
c find next monitor elemetn
      do 1000 i=2,nlat1
          nl=i
          call elname(latt,nl   ,mul,tname)
          CALL TRACKE(LATT,NL,SV,NP,'CONT',tNAME,LFNO)
          CALL TRACKE(LATT,NL,SV,NP,'POS', tNAME,LFNO)
          twiss(i,0,mfitdx) =sv(1)
          twiss(i,0,mfitdpx)=sv(2)
          twiss(i,0,mfitdy) =sv(3)
          twiss(i,0,mfitdpy)=sv(4)
          call tracke(latt,nl   ,sv,np,'DISP',name,lfno)
          if(sv(5) .ne. 0d0) then
             twiss(i,0, 7)=sv(1)
             twiss(i,0, 8)=sv(2)
             twiss(i,0, 9)=sv(3)
             twiss(i,0,10)=sv(4)
          end if
c
c         write(lfno,'(1h ,a,a,i3,4G15.6)')'trcod:',tname,i,
c    &                                     (sv(j),j=1,4)
c         write(lfno,'(1h ,a,i3,4G15.6)')'trcod0'
c    &    ,ndim, twiss(i,ndim,15),
c    &    twiss(i,ndim,17),twiss(i,ndim,16),twiss(i,ndim,18)
c
1000  continue
c
      lfno=olfno
      call tracke(latt,nl,sv,np,'CONT',name,lfno)
c      if(print)then
c        call phdrw(ix,np,word,title,case,exist,lfno)
c      endif
      call tracke(latt,nlat1,sv,np,'RESET',name,lfno)
      xa=sv(1)
      ya=sv(2)
      xxa=sv(3)
      xya=sv(5)
      yya=sv(4)
      return
c
      entry  trcodf(latt,nlat1,name,mul,print,
     &               word,title,case,exist,twiss,lfno)
      call elname(latt,nlat1,mul,tname)
      write (lfno,*) ' trcodf',nlat1 ,tname
      if(itw .ne. 0)then
         if(ilist(1,itw-1)-1 .lt. nlat*ntwissfun)then
            print *, ' size of lattice has been changed. ',
     $        nlat,'!=',(ilist(1,itw-1)-1)/ntwissfun
            call tfree(int8(itw))
            itw=0
         end if
      end if
      if(itw .eq. 0) then
         itw=ktaloc(nlat*ntwissfun)
      endif
      do 2000 i=7,10
        call tmov(twiss(1,0,i),twiss(1,ndim,i),nlat1)
        call tmov(twiss(1,0,i),rlist(itw+nlat*(i-1)),nlat1)
2000  continue
      do 2100 i=15,18
        call tmov(twiss(1,0,i),twiss(1,ndim,i),nlat1)
        call tmov(twiss(1,0,i),rlist(itw+nlat*(i-1)),nlat1)
2100  continue
      return
c
      entry  trcods(latt,nlat1,name,mul,print,
     &               word,title,case,exist,twiss,lfno)
       if(itw .eq. 0) then
         return
       end if
      do 3000 i=7,10
        call tmov(rlist(itw+nlat*(i-1)),twiss(1,ndim,i),nlat1)
3000  continue
      do 3100 i=15,18
        call tmov(rlist(itw+nlat*(i-1)),twiss(1,ndim,i),nlat1)
3100  continue
      return
c
      entry  trcodm(latt,nlat1,name,mul,print,
     &               word,title,case,exist,twiss,lfno)
          twiss(nlat1,ndim,mfitdx) =twiss(nlat1,0,mfitdx)
          twiss(nlat1,ndim,mfitdpx)=twiss(nlat1,0,mfitdpx)
          twiss(nlat1,ndim,mfitdy) =twiss(nlat1,0,mfitdy)
          twiss(nlat1,ndim,mfitdpy)=twiss(nlat1,0,mfitdpy)
       if( itw .ne. 0)then
          rlist(itw+nlat*(mfitdx-1)+nlat1-1)
     $        =twiss(nlat1,0,mfitdx)
          rlist(itw+nlat*(mfitdpx-1)+nlat1-1)
     $         =twiss(nlat1,0,mfitdpx)
          rlist(itw+nlat*(mfitdy-1)+nlat1-1)
     $         =twiss(nlat1,0,mfitdy)
          rlist(itw+nlat*(mfitdpy-1)+nlat1-1)
     $         =twiss(nlat1,0,mfitdpy)
       end if
      return
c
      end
