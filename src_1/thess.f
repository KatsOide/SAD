      subroutine thess(a,w,v,n,ndim)
      implicit none
      integer*4 n,ndim,i,i1,j,k
      real*8 a(n,n),w(ndim,n),v(n)
      real*8 h1,aa,ww,p,q
      do 1 i=2,n
        v(i)=v(i)**2
1     continue
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
              do 30 k=i1,n
                a(j ,k)=a(j ,k)-p*a(i1,k)
                a(i1,k)=a(i1,k)+q*a(j ,k)
 30           continue
              do 40 k=1,n
                a(k,i1)=a(k,i1)+p*a(k,j )
                a(k,j )=a(k,j )-q*a(k,i1)
                w(k,i1)=w(k,i1)+p*w(k,j )
                w(k,j )=w(k,j )-q*w(k,i1)
 40           continue
            else
              p=a(i1,i)/a(j,i)
              h1=v(j)+v(i1)*p**2
              q=v(i1)*p/h1
              v(j)=v(i1)*v(j)/h1
              v(i1)=h1
              do 50 k=i1,n
                aa=a(j,k)
                a(j ,k)=p*aa-a(i1,k)
                a(i1,k)=aa-q*a(j ,k)
 50           continue
              a(i1,i)=a(j,i)
              a(j ,i)=0.d0
              do 60 k=1,n
                aa=a(k,i1)
                a(k,i1)=p*aa+a(k,j )
                a(k,j )=q*a(k,i1)-aa
                ww=w(k,i1)
                w(k,i1)=p*ww+w(k,j )
                w(k,j )=q*w(k,i1)-ww
 60           continue
            endif
          endif
 20     continue
 10   continue
      do 110 i=2,n
        v(i)=sqrt(v(i))
 110  continue
      return
      end
