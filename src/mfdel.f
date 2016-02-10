      subroutine mfdel(word)
      use tfstk
      use ffs
      use tffitcode
      parameter (kfiles=3)
      logical exist
      character*(*) word
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
    1 if(exist) call getwdl(word)
      call mfdel1(word,icom,rlist(iifnam(icom)),rlist(iidata(icom)),
     1            exist)
      if(exist) goto 1
      if(nof(icom).eq.0) then
        if(iidata(icom).ne.0) then
          call tfree(int8(iidata(icom)))
          iidata(icom)=0
          call tfree(int8(iifnam(icom)))
          iifnam(icom)=0
        endif
      endif
      return
      end
