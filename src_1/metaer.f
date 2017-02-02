      subroutine metaer(twiss,ip,esclx,escly,im,imon,emon)
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     imon(nmona,4),emon(nmona,*)
      include 'inc/common.inc'
c
      j=imon(im,2)
      n=imon(j,1)
      r11=twiss(n,ip,11)
      r12=twiss(n,ip,12)
      r21=twiss(n,ip,13)
      r22=twiss(n,ip,14)
      det=r11*r22-r12*r21
      if(det.gt.1d0) then
c       qr=sqrt((det-1d0)/det)
        qr=(det-1d0)/det
        un=sqrt(det)
        emon(j,3)=sqrt(
     z          (r12**2 + (r22/twiss(n,ip,2))**2)*qr + un**2
     z                  ) *tgauss() *esclx
        emon(j,4)=sqrt(
     z           un**2 + (r21**2 + (r22/twiss(n,ip,5))**2)*qr
     z                  ) *tgauss() *escly
c       c2(1)=(-r12*c4(1) + r22*c4(2))*qr +  un*c4(3)
c       c2(2)=( r11*c4(1) - r21*c4(2))*qr              + un*c4(4)
c       c2(3)=   un*c4(1)                 +(r21*c4(3)  +r22*c4(4))*qr
c       c2(4)=               un*c4(2)     +(r11*c4(3)  +r12*c4(4))*qr
      else
        un=sqrt(1d0-det)
        emon(j,3)=sqrt(
     z          un**2 +  r22**2 + (r12/twiss(n,ip,5))**2
     z                  ) *tgauss() *esclx
        emon(j,4)=sqrt(
     z          r11**2 + (r12/twiss(n,ip,2))**2+ un**2
     z                  ) *tgauss() *escly
c       c2(1)=  un*c4(1)               -r22*c4(3)  +r12*c4(4)
c       c2(2)=            + un*c4(2)   +r21*c4(3)  -r11*c4(4)
c       c2(3)= r11*c4(1)  +r12*c4(2)   + un*c4(3)
c       c2(4)= r21*c4(1)  +r22*c4(2)               + un*c4(4)
      endif
      return
      end
