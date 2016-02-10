      subroutine pstati1(latt,twiss,imon,emon,nmon,istr,nstr,x,temp)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      dimension latt(2,nlat),twiss(nlat,-ndim:ndim,ntwissfun),
     1          imon(nmona,4),emon(nmona,4),istr(nstra,4),x(4,4),
     1          temp(*)
      include 'inc/common.inc'
c
      do 10 i=1,nmon
        j=imon(i,2)
        nq=imon(j,4)
   10   temp(i)=twiss(imon(j,1),0,15)-twiss(imon(j,1),ndim,15)
     1                         -rlist(latt(2,nq)+5)+emon(j,1)
      call mstat(temp,nmon,x(1,1),x(2,1),x(3,1),x(4,1))
      do 12 i=1,nmon
        j=imon(i,2)
        nq=imon(j,4)
   12   temp(i)=twiss(imon(j,1),0,17)-twiss(imon(j,1),ndim,17)
     1                         -rlist(latt(2,nq)+6)+emon(j,2)
      call mstat(temp,nmon,x(1,2),x(2,2),x(3,2),x(4,2))
      do 20 i=1,nmon
        j=imon(i,2)
   20   temp(i)=twiss(imon(j,1),0,7)-twiss(imon(j,1),ndim,7)
      call mstat(temp,nmon,x(1,3),x(2,3),x(3,3),x(4,3))
      do 22 i=1,nmon
        j=imon(i,2)
   22   temp(i)=twiss(imon(j,1),0,9)-twiss(imon(j,1),ndim,9)
      call mstat(temp,nmon,x(1,4),x(2,4),x(3,4),x(4,4))
      do 24 i=1,4
   24   x(1,i)=max(x(1,i),-x(2,i))
      do 26 i=1,nstr
        j=istr(i,2)
   26   temp(i)=rlist(latt(2,istr(j,1))+11)
      call mstat(temp,nstr,x(2,1),x(2,2),x(2,3),x(2,4))
      x(2,1)=max(x(2,1),-x(2,2))
      return
      end
