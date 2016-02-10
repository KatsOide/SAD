      subroutine qdquad(dtrans,dcod,al,ak,utwiss,itwissp,
     $gammab,k1,idp,dx,dy,theta,iv,nfam,nut)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nfam,nut,k1,idp,iv,itwissp(nlat),i
      real*8  utwiss(ntwissfun,-nfam:nfam,nut),gammab(nlat),
     $     al,ak,theta,
     $     dtrans(4,5),dcod(6),cod(6),dx,dy,pr,xi,pxi,yi,pyi,
     $     akk,phi,dphi,sinphi,cosphi,a1,a2,a11,a12,a21,shphi,chphi,
     $     b1,b2,b11,b12,b21,b16,b26,a16,a26,x,y
      if(iv .eq. 4)then
        call qdrotate(dtrans,dcod,utwiss,gammab,
     $       k1,itwissp(k1),idp,dx,dy,iv,
     $       ndim,nlat,nfam,nut,ntwissfun)
      elseif(iv .eq. 2)then
        if(al .le. 0.d0)then
          call qdthin(dtrans,dcod,4,al,ak,
     $         utwiss,itwissp,gammab,k1,idp,
     $         dx,dy,theta,iv,nfam,nut)
          return
        endif
        call qtentu(dtrans,cod,utwiss(1,idp,itwissp(k1)),.true.)
        call qchg(dtrans,cod,-dx,-dy,theta,.true.)
        pr=1.d0+cod(6)
        xi=cod(1)
        pxi=cod(2)
        yi=cod(3)
        pyi=cod(4)
        akk=sqrt(abs(ak/pr)/al)
        phi=akk*al
        if(ak*pr .gt. 0.d0)then
          dphi=.5d0*phi/ak
          sinphi=sin(phi)
          cosphi=cos(phi)
          a1=cosphi+sinphi/phi
          a2=(1.d0-phi)*(1.d0+phi)/phi*sinphi
          a11=-sinphi*dphi
          a12=(cosphi-sinphi/phi)*dphi/akk
          a21=-a1*dphi*akk
          shphi=sinh(phi)
          chphi=cosh(phi)
          b1=chphi+shphi/phi
          b2=(1.d0+phi**2)*shphi/phi
          b11= shphi*dphi
          b12= (chphi-shphi/phi)*dphi/akk
          b21= b1*dphi*akk
          a16= .25d0*al/pr*(a1*xi/pr+(-a2+cosphi)/ak*pxi)
          a26= .25d0/pr*((a2+3.d0*cosphi)*xi/pr+a1*al*pxi)
          b16=-.25d0*al/pr*(b1*yi/pr-(-b2+chphi)/ak*pyi)
          b26=-.25d0/pr*((b2+3.d0*chphi)*yi/pr+b1*al*pyi)
        elseif(ak*pr .lt. 0.d0)then
          dphi=.5d0*phi/ak
          sinphi=sin(phi)
          cosphi=cos(phi)
          a1=cosphi+sinphi/phi
          a2=(1.d0-phi)*(1.d0+phi)/phi*sinphi
          b11=-sinphi*dphi
          b12=(cosphi-sinphi/phi)*dphi/akk
          b21=-a1*dphi*akk
          shphi=sinh(phi)
          chphi=cosh(phi)
          b1=chphi+shphi/phi
          b2=(1.d0+phi**2)*shphi/phi
          a11= shphi*dphi
          a12= (chphi-shphi/phi)*dphi/akk
          a21= b1*dphi*akk
          b16=-.25d0*al/pr*(a1*yi/pr-(-a2+cosphi)/ak*pyi)
          b26=-.25d0/pr*((a2+3.d0*cosphi)*yi/pr+a1*al*pyi)
          a16= .25d0*al*(b1*xi/pr+(-b2+chphi)/ak*pxi)
          a26= .25d0/pr*((b2+3.d0*chphi)*xi/pr+b1*al*pxi)
        else
          a11=al*.5d0
          a12=-al**2/6.d0
          a21=-1.d0
          b11=-al*.5d0
          b12= al**2/6.d0
          b21=-1.d0
          a16=-(al*a21*xi+(-a12+al*a11)*pxi)*.5d0/pr
          a26=-((a21-ak*a11-1.d0)*xi+al*a21*pxi)*.5d0/pr
          b16=-(al*b21*yi+(-b12+al*b11)*pyi)*.5d0/pr
          b26=-((b21+ak*b11+1.d0)*yi+al*b21*pyi)*.5d0/pr
        endif
        do 10 i=1,5
          x=dtrans(1,i)
          dtrans(1,i)= a11*x+a12*dtrans(2,i)
          dtrans(2,i)= a21*x+a11*dtrans(2,i)
          y=dtrans(3,i)
          dtrans(3,i)= b11*y+b12*dtrans(4,i)
          dtrans(4,i)= b21*y+b11*dtrans(4,i)
 10     continue
        dtrans(1,5)=dtrans(1,5)+a16
        dtrans(2,5)=dtrans(2,5)+a26
        dtrans(3,5)=dtrans(3,5)+b16
        dtrans(4,5)=dtrans(4,5)+b26
        dcod(1)= a11*xi+a12*pxi
        dcod(2)= a21*xi+a11*pxi
        dcod(3)= b11*yi+b12*pyi
        dcod(4)= b21*yi+b11*pyi
        dcod(5)=0.d0
        call qchg(dtrans,dcod,0.d0,0.d0,-theta,.false.)
      else
        call tclr(dtrans,20)
        call tclr(dcod,4)
      endif
      return
      end

      subroutine qdrotate(dtrans,dcod,utwiss,gammab,
     $k1,kk1,idp,dx,dy,iv,ndim,nlat,nfam,nut,ntwissfun)
      implicit none
      integer*4 ndim,nlat,idp,k1,iv,nfam,nut,kk1,ntwissfun
      real*8  utwiss(ntwissfun,-nfam:nfam,nut),gammab(nlat),
     $     dtrans(4,5),dcod(6),cod(6),trans1(4,5),trans2(4,5),
     $     dcod1(4),dx,dy,cod2(6)
      integer*4 i
      real*8 x,px,y,py
      call qtentu(trans2,cod,utwiss(1,idp,kk1),.true.)
      call qtentu(trans1,cod2,utwiss(1,idp,kk1+1),.true.)
c      if(canon)then
c        call qangletocanon(trans2,cod)
c        call qangletocanon(trans1,cod2)
c      endif
      cod(1)=cod(1)-dx
      cod(3)=cod(3)-dy
      call tftmatu(utwiss(1,idp,kk1),utwiss(1,idp,kk1+1),
     $     utwiss(3,idp,nut),utwiss(6,idp,nut),
     $     gammab,trans1,k1,k1+1,.false.,.true.)
      do i=1,4
        dtrans(i,3)=-trans1(i,1)
        dtrans(i,1)= trans1(i,3)
        dtrans(i,4)=-trans1(i,2)
        dtrans(i,2)= trans1(i,4)
        dtrans(i,5)=0.d0
        dcod1(i)=cod2(i)
     $       -(trans1(i,1)*cod(1)+trans1(i,2)*cod(2)
     $       +trans1(i,3)*cod(3)+trans1(i,4)*cod(4))
      enddo
      dcod1(1)=dcod1(1)-dx
      dcod1(3)=dcod1(3)-dy
      do i=1,5
        dtrans(1,i)=dtrans(1,i)+trans1(3,i)
        dtrans(3,i)=dtrans(3,i)-trans1(1,i)
        dtrans(2,i)=dtrans(2,i)+trans1(4,i)
        dtrans(4,i)=dtrans(4,i)-trans1(2,i)
      enddo
      do i=1,4
        dcod(i)=
     $       dtrans(i,1)*cod(1)+dtrans(i,2)*cod(2)
     $       +dtrans(i,3)*cod(3)+dtrans(i,4)*cod(4)
      enddo
      dcod(1)=dcod(1)+dcod1(3)
      dcod(2)=dcod(2)+dcod1(4)
      dcod(3)=dcod(3)-dcod1(1)
      dcod(4)=dcod(4)-dcod1(2)
      do i=1,4
        x =dtrans(i,1)
        px=dtrans(i,2)
        y =dtrans(i,3)
        py=dtrans(i,4)
        dtrans(i,1)=x*trans2(1,1)+px*trans2(2,1)
     $       +y*trans2(3,1)+py*trans2(4,1)
        dtrans(i,2)=x*trans2(1,2)+px*trans2(2,2)
     $       +y*trans2(3,2)+py*trans2(4,2)
        dtrans(i,3)=x*trans2(1,3)+px*trans2(2,3)
     $       +y*trans2(3,3)+py*trans2(4,3)
        dtrans(i,4)=x*trans2(1,4)+px*trans2(2,4)
     $       +y*trans2(3,4)+py*trans2(4,4)
        dtrans(i,5)=x*trans2(1,5)+px*trans2(2,5)
     $       +y*trans2(3,5)+py*trans2(4,5)+dtrans(i,5)
      enddo
c      if(canon)then
c        call qcanontoangle(dtrans,dcod)
c      endif
      return
      end
