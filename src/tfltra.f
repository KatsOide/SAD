      subroutine tfltra(word,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4, parameter :: ndp=30,maxdp=ndp*2
      real*8, parameter :: dpstep=.001d0
      integer*4 lfno,ns,i,n1
      integer*8 ios,ics,isn,ic,is,itm,its,ix,imx,imy
      real*8 os1,os2,oss
      character*(*) word
      logical exist
      call tfgetr(os1,pi2,word,lfno,exist)
      if(.not. exist)then
        call termes('?Missing initial nuz for','LTR_ACK')
        go to 8001
      endif
      call tfgetr(os2,pi2,word,lfno,exist)
      if(.not. exist)then
        call termes('?Missing final nuz for','LTR_ACK')
        go to 8001
      endif
      call tfgetr(oss,pi2,word,lfno,exist)
      if(.not. exist)then
        call termes('?Missing nus step for','LTR_ACK')
        go to 8001
      endif
      if(oss .ne. 0.d0)then
        ns=min(300,nint(abs((os2-os1)/oss))+1)
      else
        ns=50
      endif
      ios=ktaloc(ns)
      if(ns .gt. 1)then
        oss=(os2-os1)/(ns-1)
      else
        oss=0.d0
      endif
      do 10 i=0,ns-1
        rlist(ios+i)=os1+i*oss
10    continue
      n1=ndp*2+1
      ics=ktaloc(ns)
      isn=ktaloc(ns)
      ic =ktaloc(ns)
      is =ktaloc(ns)
      itm=ktaloc(n1*ntwissfun)
      its=ktaloc(n1*ntwissfun)
      ix =ktaloc(maxdp*4*ns)
      imx=ktaloc(maxdp*ns/2+1)
      imy=ktaloc(maxdp*ns/2+1)
      call tfltr1(rlist(itm),rlist(its),rlist(ix),
     1            rlist(imx),rlist(imy),rlist(ios),
     1            rlist(ics),rlist(isn),rlist(ic),rlist(is),
     1            ns,lfno)
      call tfree(imy)
      call tfree(imx)
      call tfree(ix)
      call tfree(its)
      call tfree(itm)
      call tfree(is )
      call tfree(ic )
      call tfree(isn)
      call tfree(ics)
      call tfree(ios)
      return
8001  call termes('Syntax: LTR_ACK nusi nusf nusstep',' ')
      return
      end
