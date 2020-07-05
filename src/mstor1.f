      subroutine mstor1(word,latt,twiss,istr,nstr,imon,emon,nmon,
     1                  ifname,id,extend)
      use tfstk
      use ffs
      use tffitcode
      use iso_c_binding
      parameter (kfiles=3)
      parameter (integr=4,lreal8=8,ichr=lreal8*10,namel=ichr/integr)
      logical extend
      character*(*) word
      character*30 dat
      character fn(1)*(ichr)
      character,pointer :: fnn (:)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun)
      dimension id(*),ifname(*)
      dimension istr(nstra,4),imon(nmona,4),emon(nmona,4)
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
      include 'inc/common.inc'
c
      if(word.eq.' ') then
        call fdate1(dat)
        fn(1)=dat(1:11)
      else
        fn(1)=word
      endif
      extend=.false.
      icom=icomf
      nf=nof(icom)
      do 10 i=1,nf
        call c_f_pointer(c_loc(rlist(ifname(i))),fnn,[ichr])
c        call mcchar(rlist(ifname(i)),fnn,namel)
        if(fn(1).eq.fnn(1)) then
          np=i
          call tfree(int8(id(np)))
          goto 11
        endif
   10 continue
      if(nf.lt.memax(icom)) then
        nf=nf+1
        np=nf
        ifname(np)=italoc(ichr/lreal8)
      else
        extend=.true.
        return
      endif
   11 continue
      nof(icom)=nf
      call mcchar(fn,rlist(ifname(np)),namel)
c
      if(icom.eq.1) then
        id(np)=italoc(1+(3*nstr+1)/2)
        ilist(1,id(np))=nstr
        do 21 i=1,nstr
          rlist(id(np)+i)=rlist(latt(istr(istr(i,2),1))+11)
          ilist(mod(i-1,2)+1,id(np)+nstr+1+(i-1)/2)=istr(i,2)
   21   continue
      elseif(icom.eq.2) then
        id(np)=italoc(4*nlat)
c       print *,np,id(np)
        call tmov(twiss(1,0,15),rlist(id(np)),nlat)
        call tmov(twiss(1,0,16),rlist(id(np)+nlat),nlat)
        call tmov(twiss(1,0,17),rlist(id(np)+2*nlat),nlat)
        call tmov(twiss(1,0,18),rlist(id(np)+3*nlat),nlat)
      elseif(icom.eq.3) then
        id(np)=italoc(1+(5*nmon+1)/2)
        ilist(1,id(np))=nmon
        do 22 i=1,nmon
          j=imon(imon(i,2),1)
          nq=imon(imon(i,2),4)
          rlist(id(np)+i)=twiss(j,0,15)-twiss(j,ndim,15)
     1                   -rlist(latt(nq)+5)+emon(imon(i,2),1)
          rlist(id(np)+nmon+i)=twiss(j,0,17)-twiss(j,ndim,17)
     1                   -rlist(latt(nq)+6)+emon(imon(i,2),2)
          ilist(mod(i-1,2)+1,id(np)+2*nmon+1+(i-1)/2)=imon(i,2)
   22   continue
      endif
      return
      end
