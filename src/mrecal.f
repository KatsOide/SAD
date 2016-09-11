      subroutine mrecal(word,latt,twiss,istr,nstr,imon,emon,nmon,
     1                  lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (kfiles=3)
      logical exist
      character*(*) word
      integer*8 latt(nlat)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun)
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
        return
      endif
      if(word.eq.' ') then
        call permes(' ???',' No such file.',' ',lfno)
        return
      endif
      call mrecal1(word,latt,twiss,istr,nstr,imon,emon,nmon,
     1             rlist(iifnam(icom)),rlist(iidata(icom)))
      return
      end
