      subroutine tbdecoup(beam,param,param1,
     1                    emix,emiy,emix1,emiy1,sigp)
      implicit none
      integer*4 i,j
      real*8 beam(6,6),param(21),param1(21),emix,emiy,emix1,emiy1,sigp
      real*8 trans(6,6)
      real*8 a,c,d,p1,p2,p4,q1,q2,q4,m1,m2,m3,m4,r1,r2,r3,r4,f,g
      real*8 detp,detq,detr
      sigp=beam(6,6)
      if(sigp .ne. 0)then
        param(7 )=beam(6,1)/sigp
        param(8 )=beam(6,2)/sigp
        param(9 )=beam(6,3)/sigp
        param(10)=beam(6,4)/sigp
      else
        param( 7)=0.d0
        param( 8)=0.d0
        param( 9)=0.d0
        param(10)=0.d0
      endif
      do 10 i=1,4
        do 20 j=1,i
          beam(i,j)=beam(i,j)-param(6+i)*param(6+j)*sigp
          beam(j,i)=beam(i,j)
20      continue
10    continue
      d=beam(1,1)*beam(2,2)-beam(2,1)**2
      emix=sign(sqrt(abs(d)),d)
      if(emix .ne. 0)then
        param(1)=-beam(2,1)/emix
        param(2)= beam(1,1)/emix
      else
        param(1)=0.d0
        param(2)=0.d0
      endif
      d=beam(3,3)*beam(4,4)-beam(4,3)**2
      emiy=sign(sqrt(abs(d)),d)
      if(emiy .ne. 0)then
        param(4)=-beam(4,3)/emiy
        param(5)=beam(3,3)/emiy
      else
        param(4)=0.d0
        param(5)=0.d0
      endif
      p1= beam(2,2)
      p2= beam(2,1)
      p4= beam(1,1)
      q1= beam(4,4)
      q2= beam(4,3)
      q4= beam(3,3)
      m1= beam(4,2)
      m2= beam(4,1)
      m3= beam(3,2)
      m4= beam(3,1)
      r1=-m4*p1+m3*p2+m3*q2-m1*q4
      r2= m4*p2-m3*p4-m4*q2+m2*q4
      r3=-m2*p1+m1*p2+m3*q1-m1*q2
      r4= m2*p2-m1*p4-m4*q1+m2*q2
      detr=r1*r4-r2*r3
      detp=p1*p4-p2**2
      detq=q1*q4-q2**2
      if(detp .eq. detq)then
        call tclr(param1,14)
        emix1=0.d0
        emiy1=0.d0
        return
      endif
      f=detr/(detp-detq)**2
      g=4.d0*f+1.d0
      if(g .lt. 0.d0)then
        call tclr(param1,14)
        emix1=0.d0
        emiy1=0.d0
        return
      endif
      a=sqrt((g+sqrt(g))/g*.5d0)
      c=(2.d0*a**2-1.d0)/a/(detp-detq)
      r1=c*r1
      r2=c*r2
      r3=c*r3
      r4=c*r4
      param1(11)=r1
      param1(12)=r2
      param1(13)=r3
      param1(14)=r4
      call tinitr(trans)
      trans(1,1)=a
      trans(2,2)=a
      trans(3,3)=a
      trans(4,4)=a
      trans(1,3)=-r4
      trans(1,4)=r2
      trans(2,3)=r3
      trans(2,4)=-r1
      trans(3,1)=r1
      trans(3,2)=r2
      trans(4,1)=r3
      trans(4,2)=r4
      call tmultr(beam,trans,6)
      trans(1,3)=r1
      trans(1,4)=r3
      trans(2,3)=r2
      trans(2,4)=r4
      trans(3,1)=-r4
      trans(3,2)=r3
      trans(4,1)=r2
      trans(4,2)=-r1
      call tmultr(trans,beam,6)
      call tmov(trans,beam,36)
c      write(*,'(1x,:1p6g12.4)')beam
      d=beam(1,1)*beam(2,2)-beam(2,1)**2
      emix1=sign(sqrt(abs(d)),d)
      if(emix1 .ne. 0)then
        param1(1)=-beam(2,1)/emix1
        param1(2)=beam(1,1)/emix1
      else
        param1(1)=0.d0
        param1(2)=0.d0
      endif
      d=beam(3,3)*beam(4,4)-beam(4,3)**2
      emiy1=sign(sqrt(abs(d)),d)
      if(emiy1 .ne. 0)then
        param1(4)=-beam(4,3)/emiy1
        param1(5)=beam(3,3)/emiy1
      else
        param1(4)=0.d0
        param1(5)=0.d0
      endif
      param1( 7)=a*param( 7)-r4*param(9)+r2*param(10)
      param1( 8)=a*param( 8)+r3*param(9)-r1*param(10)
      param1( 9)=a*param( 9)+r1*param(7)+r2*param( 8)
      param1(10)=a*param(10)+r3*param(7)+r4*param( 8)
      return
      end
