      subroutine thess(a,w,v,n,ndim)
      implicit none
      integer*4 n,ndim,i,i1,j
      real*8 a(n,n),w(ndim,n),v(n),aa(n)
      real*8 h1,p,q
      v(2:n)=v(2:n)**2
      do 10 i=1,n-2
        i1=i+1
        do 20 j=i1+1,n
          if(a(j,i) .ne. 0.d0)then
            if(abs(a(i1,i)) .gt. abs(a(j,i)))then
              p=a(j,i)/a(i1,i)
              h1=v(i1)+v(j)*p**2
              q=v(j)*p/h1
              v(j)=v(i1)*v(j)/h1
              v(i1)=h1
              a(j,i)=0.d0
              a(j ,i1:n)=a(j ,i1:n)-p*a(i1,i1:n)
              a(i1,i1:n)=a(i1,i1:n)+q*a(j ,i1:n)
              a(:,i1)=a(:,i1)+p*a(:,j )
              a(:,j )=a(:,j )-q*a(:,i1)
              w(1:n,i1)=w(1:n,i1)+p*w(1:n,j )
              w(1:n,j )=w(1:n,j )-q*w(1:n,i1)
            else
              p=a(i1,i)/a(j,i)
              h1=v(j)+v(i1)*p**2
              q=v(i1)*p/h1
              v(j)=v(i1)*v(j)/h1
              v(i1)=h1
              aa(i1:n)=a(j,i1:n)
              a(j ,i1:n)=p*aa(i1:n)-a(i1,i1:n)
              a(i1,i1:n)=aa(i1:n)-q*a(j ,i1:n)
              a(i1,i)=a(j,i)
              a(j ,i)=0.d0
              aa=a(:,i1)
              a(:,i1)=p*aa+a(:,j )
              a(:,j )=q*a(:,i1)-aa
              aa=w(1:n,i1)
              w(1:n,i1)=p*aa+w(1:n,j )
              w(1:n,j )=q*w(1:n,i1)-aa
            endif
          endif
 20     continue
 10   continue
      v(2:n)=sqrt(v(2:n))
      return
      end
