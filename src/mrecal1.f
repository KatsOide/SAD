      subroutine mrecal1(word,latt,twiss,istr,nstr,imon,emon,nmon,
     1                   ifname,id)
      use tfstk
      use ffs
      use tffitcode
      parameter (kfiles=3)
      parameter (kstack=3)
      parameter (integr=4,lreal8=8,ichr=lreal8*10,namel=ichr/integr)
      character*(*) word
      character fnn(1)*(ichr)
      integer*8 latt(nlat)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun)
      dimension id(*),ifname(*)
      dimension istr(nstra,4),imon(nmona,4),emon(*)
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      include 'inc/common.inc'
c
      icom=icomf
      nf=nof(icom)
      do 10 i=1,nf
        call mcchar(rlist(ifname(i)),fnn,namel)
        if(word.eq.fnn(1)) then
          np=i
          goto 11
        endif
   10 continue
      return
c
   11 continue
      icomo=icoma
      icoma=icom
      call mstack1(1,1d0,latt,twiss,istr,nstr,imon,emon,nmon)
      icoma=icomo
      print *,np,id(np)
      if(icom.eq.1) then
        no=ilist(1,id(np))
        do 21 i=1,no
          j=ilist(mod(i-1,2)+1,id(np)+no+1+(i-1)/2)
          if(istr(j,3).eq.0) then
            rlist(latt(istr(j,1))+11)=rlist(id(np)+i)
          endif
   21   continue
      elseif(icom.eq.2) then
        call tmov(rlist(id(np)),twiss(1,0,15),nlat)
        call tmov(rlist(id(np)+nlat),twiss(1,0,16),nlat)
        call tmov(rlist(id(np)+2*nlat),twiss(1,0,17),nlat)
        call tmov(rlist(id(np)+3*nlat),twiss(1,0,18),nlat)
      elseif(icom.eq.3) then
        no=ilist(1,id(np))
        do 22 i=1,no
          j=ilist(mod(i-1,2)+1,id(np)+2*no+1+(i-1)/2)
          if(imon(j,3).eq.0)then
            twiss(imon(j,1),0,15)=rlist(id(np)+i)
            twiss(imon(j,1),0,17)=rlist(id(np)+no+i)
          endif
   22   continue
      endif
      call getwdl(word)
      return
      end
