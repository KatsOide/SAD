      subroutine sols33(a,b)
      implicit none
      real*8 ,intent(in):: a(3,3)
      real*8 ,intent(out):: b(3,3)
      real*8 r,c,s,a11,a12,a32,a13,a33,a22,a23,rr
      r=hypot(a(1,1),a(3,1))
      c=a(1,1)/r
      s=a(3,1)/r
      a11=r
      a12= a(2,1)*c+a(3,2)*s
      a32=-a(2,1)*s+a(3,2)*c
      a13= a(3,1)*c+a(3,3)*s
      a33=-a(3,1)*s+a(3,3)*c
      rr=b(1,1)
      b(1,1)= rr*c+b(3,1)*s
      b(3,1)=-rr*s+b(3,1)*c
      rr=b(1,2)
      b(1,2)= rr*c+b(3,2)*s
      b(3,2)=-rr*s+b(3,2)*c
      rr=b(1,3)
      b(1,3)= rr*c+b(3,3)*s
      b(3,3)=-rr*s+b(3,3)*c
      r=hypot(a11,a(2,1))
      c=a11/r
      s=a(2,1)/r
      a11=r
      rr=a12
      a12= rr*c+a(2,2)*s
      a22=-rr*s+a(2,2)*c
      rr=a13
      a13= rr*c+a(3,2)*s
      a23=-rr*s+a(3,2)*c
      rr=b(1,1)
      b(1,1)= rr*c+b(2,1)*s
      b(2,1)=-rr*s+b(2,1)*c
      rr=b(1,2)
      b(1,2)= rr*c+b(2,2)*s
      b(2,2)=-rr*s+b(2,2)*c
      rr=b(1,3)
      b(1,3)= rr*c+b(2,3)*s
      b(2,3)=-rr*s+b(2,3)*c
      r=hypot(a22,a32)
      c=a22/r
      s=a32/r
      a22=r
      rr=a23
      a23= rr*c+a33*s
      a33=-rr*s+a33*c
      rr=b(2,1)
      b(2,1)= rr*c+b(3,1)*s
      b(3,1)=-rr*s+b(3,1)*c
      rr=b(2,2)
      b(2,2)= rr*c+b(3,2)*s
      b(3,2)=-rr*s+b(3,2)*c
      rr=b(2,3)
      b(2,3)= rr*c+b(3,3)*s
      b(3,3)=-rr*s+b(3,3)*c
      b(3,1)=b(3,1)/a33
      b(3,2)=b(3,2)/a33
      b(3,3)=b(3,3)/a33
      b(2,1)=(b(2,1)-a23*b(3,1))/a22
      b(2,2)=(b(2,2)-a23*b(3,2))/a22
      b(2,3)=(b(2,3)-a23*b(3,3))/a22
      b(1,1)=(b(1,1)-a13*b(3,1)-a12*b(2,1))/a11
      b(1,2)=(b(1,2)-a13*b(3,2)-a12*b(2,2))/a11
      b(1,3)=(b(1,3)-a13*b(3,3)-a12*b(2,3))/a11
      return
      end
