      subroutine mfdel1(word,icom,ifname,id,exist)
      use tfstk
      use ffs
      use tffitcode
      parameter (kfiles=3)
      parameter (integr=4,lreal8=8,ichr=lreal8*10,namel=ichr/integr)
      logical exist
      character*(*) word
      character fnn(1)*(ichr)
      dimension id(*),ifname(*)
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
c
      exist=.false.
      icom=icomf
      nf=nof(icom)
      do 10 i=1,nf
        call mcchar(rlist(ifname(i)),fnn,namel)
        if(word.eq.fnn(1)) then
          np=i
          call tfree(int8(id(i)))
          call tfree(int8(ifname(i)))
          goto 11
        endif
   10 continue
      return
c
   11 continue
      exist=.true.
      do 12 i=np,nf-1
        ifname(i)=ifname(i+1)
        id(i)=id(i+1)
   12 continue
      nof(icom)=nf-1
      return
      end
