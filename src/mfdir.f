      subroutine  mfdir(word,lfno)
      use tfstk
      use ffs
      use tffitcode
      parameter (kfiles=3)
      parameter (lreal8=8,integr=4,ichr=lreal8*10,namel=ichr/integr)
      logical exist
      character*(*) word
      character fn(1)*(ichr)
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
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
      ia=iifnam(icom)
      do 10 i=1,nof(icom)
        ib=ilist(mod(i-1,2)+1,ia+(i-1)/2)
        call mcchar(rlist(ib),fn(1),namel)
        call mbufw(fn(1),.false.,lfno)
   10 continue
      call mbufw(' ',.true.,lfno)
      return
      end
