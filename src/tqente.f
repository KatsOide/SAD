      subroutine tqente(trans,cod,beam,al,bz,irad)
      use mathfun, only: sqrtl
      use sol, only:tsoldz
      use temw,only:tmulbs
      implicit none
      integer*4 irad,i
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6),trans2(6,6),
     $     al,bz,pr,a,pz,r,pxi,pyi,bzh,z0,f,g,fx,fy,dpz,rpz,ptmin
      parameter (ptmin=0.9999d0)
      pr=1.d0+cod(6)
      if(bz .eq. 0.d0)then
        pxi=cod(2)
        pyi=cod(4)
        a=pxi**2+pyi**2
        dpz=-a/pr/(1.d0+sqrtl(1.d0-a/pr**2))
        pz=pr+dpz
        r=-dpz/pz/pr*al
        rpz=pz/pr
        cod(1)=cod(1)+r*pxi
        cod(3)=cod(3)+r*pyi
        cod(5)=cod(5)-(2.d0+rpz)*a/2.d0/(pr+pz)*r
        f=al/pz**3
        fx=f*pxi
        fy=f*pyi
        g=dpz*(1.d0+rpz+rpz**2)
        if(irad .gt. 6)then
          call tinitr(trans2)
        endif
        trans2(1,2)=r+pxi*fx
        trans2(1,4)=pyi*fx
        trans2(1,6)=g*fx
        trans2(3,2)=pxi*fy
        trans2(3,4)=r+pyi*fy
        trans2(3,6)=g*fy
        trans2(5,2)=trans2(1,6)
        trans2(5,4)=trans2(3,6)
        trans2(5,6)=-a*g/pr*f
        trans(1,1)=trans(1,1)+trans2(1,2)*trans(2,1)
     $       +trans2(1,4)*trans(4,1)+trans2(1,6)*trans(6,1)
        trans(3,1)=trans(3,1)+trans2(3,2)*trans(2,1)
     $       +trans2(3,4)*trans(4,1)+trans2(3,6)*trans(6,1)
        trans(5,1)=trans(5,1)+trans2(5,2)*trans(2,1)
     $       +trans2(5,4)*trans(4,1)+trans2(5,6)*trans(6,1)
        trans(1,2)=trans(1,2)+trans2(1,2)*trans(2,2)
     $       +trans2(1,4)*trans(4,2)+trans2(1,6)*trans(6,2)
        trans(3,2)=trans(3,2)+trans2(3,2)*trans(2,2)
     $       +trans2(3,4)*trans(4,2)+trans2(3,6)*trans(6,2)
        trans(5,2)=trans(5,2)+trans2(5,2)*trans(2,2)
     $       +trans2(5,4)*trans(4,2)+trans2(5,6)*trans(6,2)
        trans(1,3)=trans(1,3)+trans2(1,2)*trans(2,3)
     $       +trans2(1,4)*trans(4,3)+trans2(1,6)*trans(6,3)
        trans(3,3)=trans(3,3)+trans2(3,2)*trans(2,3)
     $       +trans2(3,4)*trans(4,3)+trans2(3,6)*trans(6,3)
        trans(5,3)=trans(5,3)+trans2(5,2)*trans(2,3)
     $       +trans2(5,4)*trans(4,3)+trans2(5,6)*trans(6,3)
        trans(1,4)=trans(1,4)+trans2(1,2)*trans(2,4)
     $       +trans2(1,4)*trans(4,4)+trans2(1,6)*trans(6,4)
        trans(3,4)=trans(3,4)+trans2(3,2)*trans(2,4)
     $       +trans2(3,4)*trans(4,4)+trans2(3,6)*trans(6,4)
        trans(5,4)=trans(5,4)+trans2(5,2)*trans(2,4)
     $       +trans2(5,4)*trans(4,4)+trans2(5,6)*trans(6,4)
        trans(1,5)=trans(1,5)+trans2(1,2)*trans(2,5)
     $       +trans2(1,4)*trans(4,5)+trans2(1,6)*trans(6,5)
        trans(3,5)=trans(3,5)+trans2(3,2)*trans(2,5)
     $       +trans2(3,4)*trans(4,5)+trans2(3,6)*trans(6,5)
        trans(5,5)=trans(5,5)+trans2(5,2)*trans(2,5)
     $       +trans2(5,4)*trans(4,5)+trans2(5,6)*trans(6,5)
        trans(1,6)=trans(1,6)+trans2(1,2)*trans(2,6)
     $       +trans2(1,4)*trans(4,6)+trans2(1,6)*trans(6,6)
        trans(3,6)=trans(3,6)+trans2(3,2)*trans(2,6)
     $       +trans2(3,4)*trans(4,6)+trans2(3,6)*trans(6,6)
        trans(5,6)=trans(5,6)+trans2(5,2)*trans(2,6)
     $       +trans2(5,4)*trans(4,6)+trans2(5,6)*trans(6,6)
        do i=7,irad
          trans(1,i)=trans(1,i)+trans2(1,2)*trans(2,i)
     $         +trans2(1,4)*trans(4,i)+trans2(1,6)*trans(6,i)
          trans(3,i)=trans(3,i)+trans2(3,2)*trans(2,i)
     $         +trans2(3,4)*trans(4,i)+trans2(3,6)*trans(6,i)
          trans(5,i)=trans(5,i)+trans2(5,2)*trans(2,i)
     $         +trans2(5,4)*trans(4,i)+trans2(5,6)*trans(6,i)
        enddo
      else
        bzh=bz*.5d0
c cod are canonical!
        pxi=cod(2)+bzh*cod(3)
        pyi=cod(4)-bzh*cod(1)
        a=pxi**2+pyi**2
        pz=pr*sqrtl(1.d0-a/pr**2)
        r=a/pr/(pr+pz)*al
        call tinitr(trans2)
        trans2(2,3)= bzh
        trans2(4,1)=-bzh
        trans2(5,2)= pxi/pr/pz*al
        trans2(5,3)= pxi*bzh/pr/pz*al
        trans2(5,4)= pyi/pr/pz*al
        trans2(5,1)=-pyi*bzh/pr/pz*al
        trans2(5,5)=0.d0
        trans2(5,6)=-a/pr**2/pz*al
        cod(2)=pxi
        cod(4)=pyi
        z0=cod(5)
        call tsoldz(trans1,cod,r,0.d0,0.d0,bz,.false.)
        trans2=matmul(trans1,trans2)
c        call tmultr(trans2,trans1,6)
        cod(5)=z0-(a/pr/(pr+pz))**2*(2.d0*pr+pz)/pz/2.d0*al
        trans2(5,5)=1.d0
        trans2(5,6)=a**2*(pr**2+pr*pz+pz**2)/(pr*pz)**3/(pr+pz)*al
        cod(2)=cod(2)-bzh*cod(3)
        cod(4)=cod(4)+bzh*cod(1)
        trans2(2,1)=trans2(2,1)-bzh*trans2(3,1)
        trans2(4,1)=trans2(4,1)+bzh*trans2(1,1)
        trans2(2,2)=trans2(2,2)-bzh*trans2(3,2)
        trans2(4,2)=trans2(4,2)+bzh*trans2(1,2)
        trans2(2,3)=trans2(2,3)-bzh*trans2(3,3)
        trans2(4,3)=trans2(4,3)+bzh*trans2(1,3)
        trans2(2,4)=trans2(2,4)-bzh*trans2(3,4)
        trans2(4,4)=trans2(4,4)+bzh*trans2(1,4)
        trans2(2,6)=trans2(2,6)-bzh*trans2(3,6)
        trans2(4,6)=trans2(4,6)+bzh*trans2(1,6)
        trans2(5,1)=
     $       -trans2(1,1)*trans2(2,6)
     $       +trans2(2,1)*trans2(1,6)
     $       -trans2(3,1)*trans2(4,6)
     $       +trans2(4,1)*trans2(3,6)
        trans2(5,2)=
     $       -trans2(1,2)*trans2(2,6)
     $       +trans2(2,2)*trans2(1,6)
     $       -trans2(3,2)*trans2(4,6)
     $       +trans2(4,2)*trans2(3,6)
        trans2(5,3)=
     $       -trans2(1,3)*trans2(2,6)
     $       +trans2(2,3)*trans2(1,6)
     $       -trans2(3,3)*trans2(4,6)
     $       +trans2(4,3)*trans2(3,6)
        trans2(5,4)=
     $       -trans2(1,4)*trans2(2,6)
     $       +trans2(2,4)*trans2(1,6)
     $       -trans2(3,4)*trans2(4,6)
     $       +trans2(4,4)*trans2(3,6)
        call tmultr5(trans,trans2,irad)
c      write(*,'(a/,6(1p6g11.4/))')
c     $     'tqente-2 ',((trans(i,j),j=1,6),i=1,6)
c      write(*,*)'with irad: ',irad
      endif
      if(irad .gt. 6)then
        call tmulbs(beam,trans2,.true.)
      endif
      return
      end
