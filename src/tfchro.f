      subroutine tfchro(latt,mult,
     1              alphax,betax,psix,dx,
     1              alphay,betay,psiy,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      real*8 alphax(nlat),betax(nlat),psix(nlat)
      real*8  alphay(nlat),betay(nlat),psiy(nlat)
      real*8 dx(nlat)
      integer*4 latt(2,nlat),mult(nlat)
      character*12 name
      write(lfno,*)'Element           X Chromaticity   Y Chromaticity'
      gx=0.d0
      gy=0.d0
      do 10 i=1,nlat-1
        k=idtype(latt(1,i))
        ali=rlist(idval(latt(1,i))+1)
        if(k .eq. 4 .or. k .eq. 6 .or. k .eq. icmult)then
          if(k .eq. 4)then
            v=rlist(latt(2,i)+2)
            xix=-( v*betax(i)+(1.d0+alphax(i)**2)/betax(i)*ali+
     1           alphax(i+1)-alphax(i))*.5d0
            xiy=-(-v*betay(i)+(1.d0+alphay(i)**2)/betay(i)*ali+
     1           alphay(i+1)-alphay(i))*.5d0
          elseif(k .eq. 6)then
            v=rlist(latt(2,i)+2)
            xix= (betax(i)*dx(i)+betax(i+1)*dx(i+1))*v*.5d0
            xiy=-(betay(i)*dx(i)+betay(i+1)*dx(i+1))*v*.5d0
          elseif(k .eq. icmult)then
            v=rlist(latt(2,i)+kytbl(kwK1,icMULT))
            xix=-( v*betax(i)+(1.d0+alphax(i)**2)/betax(i)*ali+
     1           alphax(i+1)-alphax(i))*.5d0
            xiy=-(-v*betay(i)+(1.d0+alphay(i)**2)/betay(i)*ali+
     1           alphay(i+1)-alphay(i))*.5d0
            v=rlist(latt(2,i)+kytbl(kwK2,icMULT))
            xix=xix+(betax(i)*dx(i)+betax(i+1)*dx(i+1))*v*.5d0
            xiy=xiy-(betay(i)*dx(i)+betay(i+1)*dx(i+1))*v*.5d0
          endif
          call elname(latt,i,mult,name)
          write(lfno,9001)name,xix,xiy
9001      format(1x,a,1x,2f17.3)
          gx=gx+xix
          gy=gy+xiy
        endif
10    continue
      write(lfno,9001)'TOTAL       ',gx,gy
      return
      end
