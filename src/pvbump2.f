      subroutine pvbump2(isb,ncor,kv,vi,vname,id,nvbmp,istr,yplane)
      use tfstk
      use ffs
c     set bump conditions to the 'vbump' buffer.
      use tffitcode
      parameter (mfitc1=32,mfitc2=28)
      logical coup,yplane
      character*12 vname
      dimension istr(nstra,4)
      dimension isb(ncor+4),id(nvbmp)
      include 'inc/common.inc'
c
      ncorb=isb(2)
      lt=isb(1)
      coup=ncorb.ne.ncor
      if(coup) then
        nc=9
      else
        nc=5
      endif 
c     write(*,'(10i4,a)')(isb(i),i=1,10),' pvbump2'
      ip=italoc(3+2*((nc+1)/2)+nc)
      ip1=ip+3
      ip2=ip+3+(nc+1)/2
      ip3=ip+3+2*((nc+1)/2)
      id(nvbmp)=ip
      ilist(1,ip)=nc
      ilist(2,ip)=1
      ilist(1,ip+1)=kv
      call mcchar(vname,ilist(2,ip+1),3)
      if(yplane) then
        m=mfitc1+2
        mo=mfitc1
      else
        m=mfitc1
        mo=mfitc1+2
      endif 
      ilist(1,ip1)=kv
      ilist(2,ip1)=m
      ilist(1,ip1+1)=m+1
      ilist(2,ip1+1)=m
      ilist(1,ip1+2)=m+1
      ilist(1,ip2)=lt
      ilist(2,ip2)=istr(istr(isb(3),2),1)
      ilist(1,ip2+1)=istr(istr(isb(3),2),1)
      ilist(2,ip2+1)=istr(istr(isb(2+ncorb),2),1)+1
      ilist(1,ip2+2)=istr(istr(isb(2+ncorb),2),1)+1
      rlist(ip3)=vi
      rlist(ip3+1)=0d0
      rlist(ip3+2)=0d0
      rlist(ip3+3)=0d0
      rlist(ip3+4)=0d0
      if(coup) then
        ilist(2,ip1+2)=mo
        ilist(1,ip1+3)=mo+1
        ilist(2,ip1+3)=mo
        ilist(1,ip1+4)=mo+1
        ilist(2,ip2+2)=istr(istr(isb(3),2),1)
        ilist(1,ip2+3)=istr(istr(isb(3),2),1)
        ilist(2,ip2+3)=istr(istr(isb(2+ncorb),2),1)+1
        ilist(1,ip2+4)=istr(istr(isb(2+ncorb),2),1)+1
        rlist(ip3+5)=0d0
        rlist(ip3+6)=0d0
        rlist(ip3+7)=0d0
        rlist(ip3+8)=0d0
      endif
      return
      end
