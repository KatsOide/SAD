      subroutine  mstore(word,latt,twiss,istr,nstr,imon,emon,nmon)
      use tfstk
      use ffs
      use tffitcode
      parameter (kfiles=3)
      parameter (meminc=16,lreal8=8,ichr=lreal8*10)
      logical extend,exist
      character*(*) word
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
      dimension istr(*),imon(*),emon(*)
      call getwdl(word)
      exist=.true.
      if(word.eq.'C')then
        icom=1
      elseif(word.eq.'O') then
        icom=2
      elseif(word.eq.'M') then
        icom=3
      else
        if(icomf.eq.0) return
        icom=icomf
        exist=.false.
      endif
      icomf=icom
      if(exist) call getwdl(word)
      if(nof(icom).eq.0) then
        iifnam(icom)=italoc(ichr/lreal8*meminc)
        iidata(icom)=italoc((meminc+1)/2)
        memax(icom)=meminc
      endif
    1 call mstor1(word,latt,twiss,istr,nstr,imon,emon,nmon,
     z            rlist(iifnam(icom)),rlist(iidata(icom)),extend)
      if(extend) then
        call palocx(iifnam(icom),memax(icom),memax(icom)+meminc)
        call palocx(iidata(icom),memax(icom),memax(icom)+meminc)
        memax(icom)=memax(icom)+meminc
        goto 1
      endif
      if(nof(icom).eq.0) then
        call tfree(int8(iidata(icom)))
        call tfree(int8(iifnam(icom)))
      endif
      call getwdl(word)
      return
      end
