      real*8 function tluma(np,x,y,xm,ym,sigx,sigy0,bin,nmesh)
      use tmacro
      implicit real*8(a-h,o-z)
      dimension x(2,np),y(2,np)
      real*4 bin(2,0:nmesh,0:nmesh),alum
      s=nmesh/6.d0
      sigy=sigy0
      alum0=0.d0
1     call tclr(bin,2*(nmesh+1)**2)
      do 10 j=1,2
        do 20 i=1,np
          xa=((x(j,i)-xm)/sigx+3.d0)*s
          ya=((y(j,i)-ym)/sigy+3.d0)*s
          if(xa .gt. 0.d0 .and. ya .gt. 0.d0)then
            ix=int(xa)
            iy=int(ya)
            if(ix .lt. nmesh .and. iy .lt. nmesh)then
              fx=xa-ix
              fy=ya-iy
              bin(j,ix  ,iy  )=bin(j,ix  ,iy  )+(1.d0-fx)*(1.d0-fy)
              bin(j,ix+1,iy  )=bin(j,ix+1,iy  )+      fx *(1.d0-fy)
              bin(j,ix  ,iy+1)=bin(j,ix  ,iy+1)+(1.d0-fx)*      fy
              bin(j,ix+1,iy+1)=bin(j,ix+1,iy+1)+      fx *      fy
            endif
          endif
20      continue
10    continue
      alum=alum0
      do 30 i=0,nmesh
        do 40 j=0,nmesh
          alum=alum+bin(1,j,i)*bin(2,j,i)
c         alum=alum+(bin(1,j,i)**2+bin(2,j,i)**2)*.5d0
40      continue
30    continue
      tluma=np**2/alum*sigx*sigy*(6.d0/nmesh)**2/4.d0/pi
      if(tluma .lt. .8d0*sigx*sigy)then
        sigy1=tluma/sigx
        scale=max(.8d0,sigy1/sigy)
        j0=(-scale*3.d0+3.d0)*s
        j1=( scale*3.d0+3.d0)*s
        sigy1=(sigy*(j1-j0))/nmesh
        alum=alum0
        do 110 i=0,nmesh
          do 120 j=0,j0-1
            alum=alum+bin(1,i,j)*bin(2,i,j)
c           alum=alum+(bin(1,i,j)**2+bin(2,i,j)**2)*.5d0
120       continue
          do 130 j=j1+1,nmesh
            alum=alum+bin(1,i,j)*bin(2,i,j)
c           alum=alum+(bin(1,i,j)**2+bin(2,i,j)**2)*.5d0
130       continue
110     continue
        alum0=alum*sigy1/sigy
        sigy=sigy1
        go to 1
      endif
      return
      end
