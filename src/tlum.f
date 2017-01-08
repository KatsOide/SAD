      real*8 function tlum(np,x,y,xm,ym,sigx,sigy,nmesh)
      use tfstk
      implicit none
      integer*8 ibin
      integer*4 np,nmesh
      real*8 x(2,np),y(2,np),xm,ym,tluma,sigx,sigy
      nmesh=max(int(sqrt(real(np))),6)
      ibin=ktaloc(2*(nmesh+1)**2)
      if(ibin .le. 0)then
        write(*,*)'Insufficient memory in tlum, mesh size =',nmesh
        tlum=0.d0
        return
      endif
      tlum=tluma(np,x,y,xm,ym,sigx,sigy,rlist(ibin),nmesh)
      call tfree(ibin)
      return
      end
