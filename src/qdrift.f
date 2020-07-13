      subroutine qdrift(trans,cod,al,coup)
      use mathfun, only:sqrtl
      implicit none
      real*8 ,intent(out):: trans(4,5)
      real*8 ,intent(inout):: cod(6)
      real*8 ,intent(in):: al
      real*8 pxi,pyi,pxisq,pyisq,pzi,ale,alz
      logical*4 ,intent(out):: coup
      pxi=cod(2)
      pyi=cod(4)
      pxisq=pxi**2
      pyisq=pyi**2
      pzi=sqrtl(1.d0-pxisq-pyisq)
      ale=al/pzi
      alz=ale/pzi**2
      cod(1)=cod(1)+pxi*ale
      cod(3)=cod(3)+pyi*ale
      trans(1,1)=1.d0
      trans(1,2)=ale+pxisq*alz
      trans(1,3)=0.d0
      trans(1,4)=pxi*pyi*alz
      trans(1,5)=0.d0
      trans(2,1)=0.d0
      trans(2,2)=1.d0
      trans(2,3)=0.d0
      trans(2,4)=0.d0
      trans(2,5)=0.d0
      trans(3,1)=0.d0
      trans(3,2)=trans(1,4)
      trans(3,3)=1.d0
      trans(3,4)=ale+pyisq*alz
      trans(3,5)=0.d0
      trans(4,1)=0.d0
      trans(4,2)=0.d0
      trans(4,3)=0.d0
      trans(4,4)=1.d0
      trans(4,5)=0.d0
      coup=.true.
      return
      end

      subroutine qtest(trans,cod,al,angle,coup)
      implicit none
      real*8 trans(4,5),cod(6),al
      real*8 angle,rho
      logical coup

      rho    = al/angle

      cod(1)=cod(1)*cos(angle)+cod(2)*rho*cos(angle)
     $ +rho*(1-cos(angle))*cod(6)
      cod(2)=-sin(angle)/rho*cod(1)+cod(2)*cos(angle)
     $ +sin(angle)*cod(6)
      cod(3)=cod(3)+cod(4)*al
      cod(4)=cod(4)

      write(*,*)'qtest ',al,angle,rho,cod

      trans(1,1)=cos(angle)
      trans(1,2)=rho*sin(angle)
      trans(1,3)=0.d0
      trans(1,4)=0.d0
      trans(1,5)=rho*(1-cos(angle))
      trans(2,1)=-sin(angle)/rho
      trans(2,2)=cos(angle)
      trans(2,3)=0.d0
      trans(2,4)=0.d0
      trans(2,5)=sin(angle)
      trans(3,1)=0.d0
      trans(3,2)=0.d0
      trans(3,3)=1.d0
      trans(3,4)=al
      trans(3,5)=0.d0
      trans(4,1)=0.d0
      trans(4,2)=0.d0
      trans(4,3)=0.d0
      trans(4,4)=1.d0
      trans(4,5)=0.d0
      coup=.true.
      return
      end
