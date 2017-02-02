      subroutine trackb(latt,nlat1,name,print,
     z     word,title,case,exist,
     $     kx,irtc,ret,
     $     xa,ya,xxa,xya,yya,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
ckikuchi ... 1 line modified
      logical print, exist
ckikuchi ... 1 line added
      type (sad_descriptor) kx
      type (sad_list), pointer :: klx,klx1,klx2,klxi
      integer*4 irtc,lfno,nlat1,i,np
      real*8 result(6,7),xa,ya,xxa,xya,yya
      character*(*) word,title,case
      character*(MAXPNAME) name
      integer*8 latt(nlat)
      real*8 sv(5)
      logical*4 ret
      integer*4 ix0,ix,ix1,ix2,ix3,ix4
      common /tt/ ix0,ix,ix1,ix2,ix3,ix4
      call tracker(latt,nlat1,sv,result,np,'STANDBY',name,lfno)
      call tracker(latt,nlat1,sv,result,np,'CONT',name,lfno)
      xa=sv(1)
      ya=sv(2)
      xxa=sv(3)
      xya=sv(5)
      yya=sv(4)
      if(print)then
        call phdrw(ix,np,word,title,case,exist,lfno)
      endif
      if(ret)then
        kx=kxadaloc(-1,3,klx)
        klx%rbody(1)=dble(np)
        klx%dbody(2)=kxavaloc(0,6,klx1)
        klx1%rbody(1:6)=result(:,7)
        klx%dbody(3)=kxadaloc(0,6,klx2)
        do i=1,6
          klx2%dbody(i)=kxavaloc(0,6,klxi)
          klxi%rbody(1:6)=result(:,i)
        enddo
        irtc=0
      endif
      call tracker(latt,nlat1,sv,result,np,'RESET',name,lfno)
      return
      end
