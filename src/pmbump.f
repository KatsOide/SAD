      real*8 function pmbump(p,id,nd,latt,twiss,mult,gammab,size,istr,
     &     estr,nstr,observ,iobs,leave,synchb,demi,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      real*8 tgauss
      parameter (meminc0=1800,meminc=720,item=20)
      logical leave,synchb
      character*(*) observ
      dimension p(nd),id(nd)
      dimension latt(2,*),twiss(nlat,-ndim:ndim,ntwissfun),
     $     mult(*),gammab(*),size(21,*)
      dimension istr(*),estr(*),iobs(*)
      include 'inc/common.inc'
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
c
c     print *,mfalloc(-1),' (pmbump-start)'
      call pcbak(latt,twiss)
      ibckup1=italoc(nlat)
      call tmov(rlist(ibckup),rlist(ibckup1),nlat)
c     set multiple bump.
      do 100 l=1,nd
        ip=id(l)
        nc=ilist(1,ip)
        lv=ilist(2,ip)
        ip1=ip+3
        ip2=ip+3+(nc+1)/2
        ip3=ip+3+2*((nc+1)/2)
        rlist(ip3-1+lv)=p(l)
        im=italoc((nc+1)/2)
c       print *,'pmbump: im=',im
        do 10 i=1,nc
          ilist(mod(i-1,2)+1,im+(i-1)/2)=1
 10     continue
c       if(leave) then
c         print *,'ip=',ip,' nc=',nc,' lv=',lv,' p(',l,')=',p(l)
c         write(*,*)(ilist(mod(i-1,2)+1,ip1+(i-1)/2),i=1,nc)
c         write(*,*)(ilist(mod(i-1,2)+1,ip2+(i-1)/2),i=1,nc)
c         write(*,*)(rlist(ip3+i-1),i=1,nc)
c       endif
        call pcbak(latt,twiss)
        call pbumps(istr,nstr,rlist(ip1),rlist(ip2),rlist(im),nc,
     &       .true.)
c       print *,mfalloc(-1),' (avans-pbump)',l
        call pbump(latt,twiss,gammab,ndim,mult,istr,
     1       estr,nstr,rlist(ip1),rlist(ip2),rlist(im),rlist(ip3),nc,
     1       .false.,.false.,.false.,lfno)
        call pbumps(istr,nstr,rlist(ip1),rlist(ip2),rlist(im),nc,
     &       .false.)
c       print *,'im=',im
        call tfree(int8(im))
c        do i=1,nlat-1
c          if(idtype(latt(1,i)).eq.icbend) then
c            if(abs(rlist(latt(2,i)+11)).gt.1d-10) then
c              print *,' kick at',i,rlist(latt(2,i)+11)
c            endif
c          endif
c        enddo
 100  continue
      itemno=iobs(1)*item
      if(nemitd.eq.0) then
        ipemibuf=italoc(meminc0)
        nememax=meminc0
      elseif((nemitd+1)*itemno.gt.nememax) then
        call palocx(ipemibuf,nemitd*itemno,
     $       nemitd*itemno+max(meminc,itemno))
        nememax=nemitd*itemno+max(meminc,itemno)
      endif
c     print *,mfalloc(-1),ipemibuf,nememax
      nemitd=nemitd+1
      pmbump=pmeas(latt,twiss,gammab,size,observ,iobs,
     $     rlist(ipemibuf+(nemitd-1)*itemno),.true.,synchb,lfno)
      if(demi.ne.0d0) then
C*DEC changed by Y.Tange 11-Jan-1995
        pmbump=pmbump*(1d0+tgauss()*demi)
C*DEC End
C*HP
C       pmbump=pmbump*(1d0+gauss()*demi)
C*HP
      endif
c     print *,mfalloc(-1),ipemibuf,nememax
      if(.not.leave) then 
        call tmov(rlist(ibckup1),rlist(ibckup),nlat)
        call pundo(latt,twiss,gammab,0,1d0+dp0)
      endif
      call tfree(int8(ibckup1))
      return
      end
