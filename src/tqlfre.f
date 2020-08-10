      subroutine tqlfre(trans,cod,beam,al,ak,f1,f2,bz)
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      implicit none
      real*8 ,intent(in):: al,ak,f1,f2,bz
      real*8 af1,af2,p,a,b,ea,bzph,bp,xf,yf,pxf,pyf,f,fdp,bb
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42)
      real*8 trans1(6,6)
c      write(*,*)'tflfre ',f1,f2
      af1=-ak/al*f1*abs(f1)/24.d0
      af2=ak/al*f2
      p=1.d0+cod(6)
      a=af1/p
      b=af2/p
      ea=exp(a)
      if(bz .eq. 0.d0)then
        bzph=0.d0
        bp=b/p
        xf=ea*cod(1)+bp*cod(2)
        yf=cod(3)/ea-bp*cod(4)
        pxf=cod(2)/ea
        pyf=cod(4)*ea
        trans1(1,1)= ea
        trans1(1,2)= bp
        trans1(1,6)= 
     $       (-a*ea*cod(1)-bp*2.d0*cod(2))/p
        trans1(3,3)= 1.d0/ea
        trans1(3,4)=-bp
        trans1(3,6)= 
     $        ( a/ea*cod(3)+bp*2.d0*cod(4))/p
        trans1(2,2)= 1.d0/ea
        trans1(2,6)= a/p*pxf
        trans1(4,4)= ea
        trans1(4,6)=-a/p*pyf
        trans1(5,1)=
     $       -trans1(1,1)*trans1(2,6)
        trans1(5,2)=
     $       -trans1(1,2)*trans1(2,6)+trans1(2,2)*trans1(1,6)
        trans1(5,3)=
     $       -trans1(3,3)*trans1(4,6)
        trans1(5,4)=
     $       -trans1(3,4)*trans1(4,6)+trans1(4,4)*trans1(3,6)
        trans1(5,6)=
     $       (((2.d0+a)*a*ea*cod(1)-bp*ea*(2.d0+.5d0*a)*pxf)/p*pxf
     $       -((2.d0-a)*a/ea*cod(3)-bp/ea*(2.d0-.5d0*a)*pyf)/p*pyf
     $       -(a*ea*cod(1)+bp*ea*(2.d0+a)*pxf)*trans1(2,6)
     $       +(a/ea*cod(3)+bp/ea*(2.d0-a)*pyf)*trans1(4,6))/p
c        do i=1,irad
          trans(5,1:irad)=trans1(5,1)*trans(1,1:irad)
     $       +trans1(5,2)*trans(2,1:irad)
     $       +trans1(5,3)*trans(3,1:irad)+trans1(5,4)*trans(4,1:irad)
     $       +            trans(5,1:irad)+trans1(5,6)*trans(6,1:irad)
          trans(1,1:irad)=trans1(1,1)*trans(1,1:irad)
     $         +trans1(1,2)*trans(2,1:irad)
     $         +trans1(1,6)*trans(6,1:irad)
          trans(2,1:irad)=trans1(2,2)*trans(2,1:irad)
     $         +trans1(2,6)*trans(6,1:irad)
          trans(3,1:irad)=trans1(3,3)*trans(3,1:irad)
     $         +trans1(3,4)*trans(4,1:irad)
     $         +trans1(3,6)*trans(6,1:irad)
          trans(4,1:irad)=trans1(4,4)*trans(4,1:irad)
     $         +trans1(4,6)*trans(6,1:irad)
c        enddo
        if(irad .gt. 6)then
          trans1(1,3)= 0.d0
          trans1(1,4)= 0.d0
          trans1(1,5)=0.d0
          trans1(3,1)= 0.d0
          trans1(3,2)= 0.d0
          trans1(3,5)=0.d0
          trans1(2,1)= 0.d0
          trans1(2,3)= 0.d0
          trans1(2,4)= 0.d0
          trans1(2,5)=0.d0
          trans1(4,1)= 0.d0
          trans1(4,2)= 0.d0
          trans1(4,3)= 0.d0
          trans1(4,5)=0.d0
          trans1(5,5)=1.d0
        endif
      else
        bzph=.5d0*bz
        bp=b/p
        bb=bp*bzph
        f=1.d0/((1.d0-bb)*(1.d0+bb))
        fdp=2.d0*f*bb**2/p
        xf=(ea*cod(1)+bp*(cod(2)+bzph*cod(3)/ea-bb*cod(4)))*f
        yf=(cod(3)/ea-bp*(cod(4)-bzph*cod(1)*ea-bb*cod(2)))*f
c cod has canonical momenta!
        pxf=(cod(2)+bzph*yf)/ea
        pyf=(cod(4)-bzph*xf)*ea
        trans1(1,1)= ea*f
        trans1(1,2)= bp*f
        trans1(1,3)= bp*bzph/ea*f
        trans1(1,4)=-bp*bb*f
        trans1(1,5)=0.d0
        trans1(1,6)= xf*fdp
     $       +(-a*ea*cod(1)-bp*(2.d0*cod(2)+bzph/ea*cod(3)*(2.d0+a)
     $       -3.d0*bb*cod(4)))*f/p
        trans1(3,1)= bp*bzph*ea*f
        trans1(3,2)= bp*bb*f
        trans1(3,3)= f/ea
        trans1(3,4)=-bp*f
        trans1(3,5)=0.d0
        trans1(3,6)= yf*fdp
     $       +( a/ea*cod(3)+bp*(2.d0*cod(4)-bzph*ea*cod(1)*(2.d0-a)
     $       -3.d0*bb*cod(2)))*f/p
        trans1(2,1)= bzph/ea*trans1(3,1)
        trans1(2,2)= (1.d0+bzph*trans1(3,2))/ea
        trans1(2,3)= bzph/ea*trans1(3,3)-bzph
        trans1(2,4)= bzph/ea*trans1(3,4)
        trans1(2,5)=0.d0
        trans1(2,6)= a/p*pxf+bzph/ea*trans1(3,6)
        trans1(4,1)=-bzph*ea*trans1(1,1)+bzph
        trans1(4,2)=-bzph*ea*trans1(1,2)
        trans1(4,3)=-bzph*ea*trans1(1,3)
        trans1(4,4)= (1.d0-bzph*trans1(1,4))*ea
        trans1(4,5)=0.d0
        trans1(4,6)=-a/p*pyf-bzph*ea*trans1(1,6)
        trans1(5,1)=
     $       -trans1(1,1)*trans1(2,6)+trans1(2,1)*trans1(1,6)
     $       -trans1(3,1)*trans1(4,6)+trans1(4,1)*trans1(3,6)
        trans1(5,2)=
     $       -trans1(1,2)*trans1(2,6)+trans1(2,2)*trans1(1,6)
     $       -trans1(3,2)*trans1(4,6)+trans1(4,2)*trans1(3,6)
        trans1(5,3)=
     $       -trans1(1,3)*trans1(2,6)+trans1(2,3)*trans1(1,6)
     $       -trans1(3,3)*trans1(4,6)+trans1(4,3)*trans1(3,6)
        trans1(5,4)=
     $       -trans1(1,4)*trans1(2,6)+trans1(2,4)*trans1(1,6)
     $       -trans1(3,4)*trans1(4,6)+trans1(4,4)*trans1(3,6)
        trans1(5,5)=1.d0
        trans1(5,6)=
     $       (((2.d0+a)*a*ea*cod(1)-bp*ea*(2.d0+.5d0*a)*pxf)/p*pxf
     $       -((2.d0-a)*a/ea*cod(3)-bp/ea*(2.d0-.5d0*a)*pyf)/p*pyf
     $       -(a*ea*cod(1)+bp*ea*(2.d0+a)*pxf)*trans1(2,6)
     $       +(a/ea*cod(3)+bp/ea*(2.d0-a)*pyf)*trans1(4,6))/p
        call tmultr5(trans,trans1,irad)
      endif
      if(irad .gt. 6)then
        trans1(6,1:5)=0.d0
        trans1(6,6)=1.d0
        call tmulbs(beam ,trans1,.true.)
      endif
      cod(5)=cod(5)-((a*ea*cod(1)+bp*ea*(1.d0+.5d0*a)*pxf)*pxf
     $     -(a*cod(3)/ea+bp/ea*(1.d0-.5d0*a)*pyf)*pyf)/p
      cod(2)=pxf-bzph*cod(3)
      cod(4)=pyf+bzph*cod(1)
c cod has canonical momenta!
      cod(1)=xf
      cod(3)=yf
      return
      end
