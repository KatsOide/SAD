      subroutine ttinit(latt,x,px,y,py,z,g,dv)
      use tfstk
      use ffs_flag
      use tffitcode
      use tmacro
      implicit real*8 (a-h,o-z)
      dimension latt(2,nlat)
      dimension x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      trf0=0.d0
      vcphic=0.d0
      vcalpha=1.d0
      if(idtype(latt(1,1)) .ne. icMARK)then
        write(*,*)'Missing entrance condition of beam.'
        nturn=0
        return
      endif
      kx=latt(2,1)
      axi=rlist(kx+1)
      ayi=rlist(kx+4)
      bxi=rlist(kx+2)
      byi=rlist(kx+5)
      emx=rlist(kx+kytbl(kwEMIX,icMARK))
      emy=rlist(kx+kytbl(kwEMIY,icMARK))
      dpmax=rlist(kx+kytbl(kwDP,icMARK))
      exi=rlist(kx+ 7)
      epxi=rlist(kx+ 8)
      eyi=rlist(kx+ 9)
      epyi=rlist(kx+10)
      r1i=rlist(kx+11)
      r2i=rlist(kx+12)
      r3i=rlist(kx+13)
      r4i=rlist(kx+14)
      detr=rlist(kx+mfitdetr)
      amu=sqrt(1.d0-(r1i*r4i-r2i*r3i))
      dxi=rlist(kx+mfitdx)
      dpxi=rlist(kx+mfitdpx)
      dyi=rlist(kx+mfitdy)
      dpyi=rlist(kx+mfitdpy)
      sigz=rlist(kx+kytbl(kwSIGZ,icMARK))
      azi=rlist(kx+kytbl(kwAZ,icMARK))
      dxj=rlist(kx+26)
      dpxj=rlist(kx+27)
      dyj=rlist(kx+28)
      dpyj=rlist(kx+29)
      dxi=dxi+dxj*tgauss()
      dpxi=dpxi+dpxj*tgauss()
      dyi=dyi+dyj*tgauss()
      dpyi=dpyi+dpyj*tgauss()
      sigx=sqrt(bxi*emx)
      sigpx=sqrt(emx/bxi)
      sigy=sqrt(byi*emy)
      sigpy=sqrt(emy/byi)
      rz=sqrt(1.d0+azi**2)
c      if(twake .or. lwake)then
c        nb=ilist(1,iwakepold)
c        ns=ilist(1,iwakepold+2)
c        sb=rlist(iwakepold+1)
c        n=np0/nb/ns
c        ndp=n/8
c        if(ns .eq. 1)then
c          ds=0.d0
c        else
c          sp=min(3.d0,sqrt(ns/2.d0))
c          ds=sigz/(ns-1)*2.d0*sp
c        endif
c        if(sigz .eq. 0.d0)then
c          ddp=0.d0
c        else
c          ddp=-ds/sigz*azi/rz*dpmax
c        endif
c        rlist(iwakepold+3)=ds
c        call tfree(int8(ilist(2,iwakepold+2)))
c        iwsp=italoc(ns)
c        ilist(2,iwakepold+2)=iwsp
c        call tfree(int8(ilist(2,iwakepold+4)))
c        iwwp=italoc(np0)
c        ilist(2,iwakepold+4)=iwwp
c        if(ds .ne. 0.d0)then
c          s=0.d0
c          do 250 i=1,ns
c            rlist(iwsp+i-1)=exp(-((i-(ns+1)*.5d0)*ds/sigz)**2*.5d0)
c            s=s+rlist(iwsp+i-1)
c250       continue
c          do 251 i=1,ns
c            rlist(iwsp+i-1)=rlist(iwsp+i-1)/s
c251       continue
c        else
c          do 260 i=1,ns
c            rlist(iwsp+i-1)=1.d0
c260       continue
c        endif
c$$$        do 210 i=1,nb
c$$$          do 220 j=1,ns
c$$$            dz=-ds*(j-(ns+1)*.5d0)
c$$$            zij=-(i-1)*sb+dz
c$$$            dp=ddp*(j-(ns+1)*.5d0)
c$$$            l=((i-1)*ns+j-1)*n
c$$$            if(ndp .eq. 1)then
c$$$              dpc=0.d0
c$$$            else
c$$$              dpc=dpmax*2.d0/(ndp-1)
c$$$            endif
c$$$            do 230 k=1,ndp
c$$$              pc=(k-1)*dpc-dpmax
c$$$              pj=pc+dp+dp0
c$$$              koff=l+(k-1)*8
c$$$              do 240 m=koff+1,koff+8
c$$$                x(m)=dxi+exi*pj
c$$$                px(m)=dpxi-axi/bxi*x(m)+epxi*pj
c$$$                y(m)=dyi+eyi*pj
c$$$                py(m)=dpyi-ayi/byi*y(m)+epyi*pj
c$$$                z(m)=zij
c$$$                g(m)=pj
c$$$                rlist(iwwp+m-1)=rlist(iwsp+j-1)
c$$$240           continue
c$$$              x(koff+1)=x(koff+1)+sigx*2.d0
c$$$              px(koff+1)=px(koff+1)-axi/bxi*sigx*2.d0
c$$$              px(koff+3)=px(koff+3)+sigpx*2.d0
c$$$              y(koff+5)=y(koff+5)+sigy*2.d0
c$$$              py(koff+5)=py(koff+5)-ayi/byi*sigy*2.d0
c$$$              py(koff+7)=py(koff+7)+sigpy*2.d0
c$$$              x(koff+2)=x(koff+2)-sigx*2.d0
c$$$              px(koff+2)=px(koff+2)+axi/bxi*sigx*2.d0
c$$$              px(koff+4)=px(koff+4)-sigpx*2.d0
c$$$              y(koff+6)=y(koff+6)-sigy*2.d0
c$$$              py(koff+6)=py(koff+6)+ayi/byi*sigy*2.d0
c$$$              py(koff+8)=py(koff+8)-sigpy*2.d0
c$$$230         continue
c$$$220       continue
c$$$210     continue
c      else
        do 10 i=1,np0
          pk=tgauss()
          if(twake .or. lwake)then
            pk=0.d0
          endif
          if(gauss)then
            pj=dp0+dpmax*(tgauss()-azi*pk)/rz
          else
            pj=dp0+dpmax*((tran()-.5d0)*2.d0-azi*pk)/rz
          endif
          if(trgauss)then
            xi  =sigx*tgauss()
            xi0=xi+exi*pj
            pxi0=sigpx*tgauss()-axi/bxi*xi+epxi*pj
            yi  =sigy*tgauss()
            yi0=yi+eyi*pj
            pyi0=sigpy*tgauss()-ayi/byi*yi+epyi*pj
          else
1001        r1=2.d0*tran()-1.d0
            r2=2.d0*tran()-1.d0
            if(r1**2+r2**2 .gt. 1.d0)then
              go to 1001
            endif
            xi  =sigx*r1
            xi0=xi+exi*pj
            pxi0=sigpx*r2-axi/bxi*xi+epxi*pj
1002        r1=2.d0*tran()-1.d0
            r2=2.d0*tran()-1.d0
            if(r1**2+r2**2 .gt. 1.d0)then
              go to 1002
            endif
            yi  =sigy*r1
            yi0=yi+eyi*pj
            pyi0=sigpy*r2-ayi/byi*yi+epyi*pj
          endif
          if(detr .le. 1.d0)then
            x(i) =amu*xi0 +r4i*yi0-r2i*pyi0+dxi
            px(i)=amu*pxi0-r3i*yi0+r1i*pyi0+dpxi
            y(i) =amu*yi0 -r1i*xi0-r2i*pxi0+dyi
            py(i)=amu*pyi0-r3i*xi0-r4i*pxi0+dpyi
          else
            y(i) =amu*xi0 +r4i*yi0-r2i*pyi0+dyi
            py(i)=amu*pxi0-r3i*yi0+r1i*pyi0+dpyi
            x(i) =amu*yi0 -r1i*xi0-r2i*pxi0+dxi
            px(i)=amu*pyi0-r3i*xi0-r4i*pxi0+dpxi
          endif
          z(i)=sigz*pk
          g(i)=pj
10      continue
        if(.not. jitter)then
          sp=0.d0
          sz=0.d0
          sx=0.d0
          spx=0.d0
          sy=0.d0
          spy=0.d0
          do 110 i=1,np0
            sp=sp+g(i)
            sz=sz+z(i)
            sx=sx+x(i)
            spx=spx+px(i)
            sy=sy+y(i)
            spy=spy+py(i)
110       continue
          sp=sp/np0-dp0
          sz=sz/np0
          sx=sx/np0-exi*dp0-dxi
          spx=spx/np0-(epxi*dp0-axi/bxi*exi*dp0)-dpxi
          sy=sy/np0-eyi*dp0-dyi
          spy=spy/np0-(epyi*dp0-ayi/byi*eyi*dp0)-dpyi
          do 120 i=1,np0
            g(i)=g(i)-sp
            z(i)=z(i)-sz
            x(i)=x(i)-sx
            px(i)=px(i)-spx
            y(i)=y(i)-sy
            py(i)=py(i)-spy
120       continue
        endif
c      endif
      call tconvm(np0,px,py,g,dv,-1)
      return
      end
