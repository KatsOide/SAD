      subroutine tsvdm(a,b,x,n,m,ndim,ndimb,epslon,inv)
c
c   Singular Value Decomposition
c
c   a(ndim,m)   array contains the matrix a(n,m).
c   b(ndimb,n)  array to contain U .
c   x(m)        vector to contain (if inv W^-1 else W)
c   epslon      threshold.
c               V is returned in a .
c               a=U^T.W.V .
c    
c
      use tfstk, only:ktfenanq
      use mathfun,only:nmaxsvd,nomp,sqrt1
c      use omp_lib
      implicit none
      integer*4 ,parameter ::itmax=256
      integer*4 ,intent(in):: n,m,ndim,ndimb
      real*8 ,intent(inout):: a(ndim,m),b(ndimb,n)
      real*8 ,intent(out):: x(m)
      real*8 ,intent(in):: epslon
      real*8 ,allocatable::v(:),aa(:),bb(:)
      integer*4 ,allocatable::lsep(:)
      real*8 anorm
      real*8 f,g,s,r,w,h,xmin,z,vv,c,p,y,an
      real*8 q,h1,h2,t,r1,r2
      integer*4 i,j,mn,it,isep,ibegin,iend,i1,i1mn,n1,kkk
      logical*4 inv
      mn=min(n,m)
      n1=min(ndimb,n)
      if(dble(mn)**2*dble(max(n,m)) > dble(nmaxsvd))then
        write(*,*)' TSVDM Too large matrix. ',n,m
        return
      endif
c      call omp_set_dynamic(.true.)
      allocate(v(0:2*max(m,n)),lsep(0:n),aa(m),bb(n))
      v(1:n)=1.d0
      x(1:m)=1.d0
      b(1:n1,1:n1)=0.d0
      do j=1,n1
        b(j,j)=1.d0
      enddo
      do i=1,mn
        i1=i+1
        if(i < n)then
          do j=i1,n
            if(abs(a(i,i)) > abs(a(j,i)))then
              p=a(j,i)/a(i,i)
              h1=v(i)+v(j)*p**2
              q=v(j)*p/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
c              call omp_set_num_threads(int(1+(m-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,b,i,j,i1,m,n1,p,q)
              a(j,i1:m)=a(j,i1:m)-p*a(i,i1:m)
              a(i,i1:m)=a(i,i1:m)+q*a(j,i1:m)
c!$OMP END PARALLEL WORKSHARE
              a(j,i)=0.d0
              b(j,1:n1)=b(j,1:n1)-p*b(i,1:n1)
              b(i,1:n1)=b(i,1:n1)+q*b(j,1:n1)
            elseif(a(j,i) /= 0.d0)then
              p=a(i,i)/a(j,i)
              h1=v(j)+v(i)*p**2
              q=v(i)*p/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
c              call omp_set_num_threads(int(1+(m-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,b,aa,bb,i,j,i1,m,n1,p,q)
              aa(i1:m)=a(j,i1:m)
              a(j,i1:m)=p*aa(i1:m)-a(i,i1:m)
              a(i,i1:m)=aa(i1:m)-q*a(j,i1:m)
c!$OMP END PARALLEL WORKSHARE
              a(i,i)=a(j,i)
              a(j,i)=0.d0
              bb(1:n1)=b(j,1:n1)
              b(j,1:n1)=p*bb(1:n1)-b(i,1:n1)
              b(i,1:n1)=bb(1:n1)-q*b(j,1:n1)
            endif
          enddo
        endif
        if(i1 < m)then
          do j=i1+1,m
            r1=x(i1)*a(i,i1)
            r2=x(j )*a(i,j )
            r=sign(hypot(r1,r2),r1)
            if(r /= 0.d0)then
              c=r1/r
              s=r2/r
              if(abs(c) > abs(s))then
                a(i,j)=0.d0
                h1=x(i1)/c
                h2=x(j )*c
                p=s*x(i1)/h2
                q=s*x(j )/h1
                x(i1)=h1
                x(j )=h2
c                call omp_set_num_threads(int(1+(n-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,p,q,i1,j,n)
                a(i1:n,j )=a(i1:n,j )-p*a(i1:n,i1)
                a(i1:n,i1)=a(i1:n,i1)+q*a(i1:n,j )
c!$OMP END PARALLEL WORKSHARE
              else
                a(i,i1)=a(i,j)
                a(i,j)=0.d0
                h1=x(j )/s
                h2=x(i1)*s
                p=c*x(j )/h2
                q=c*x(i1)/h1
                x(i1)=h1
                x(j )=h2
c                call omp_set_num_threads(int(1+(n-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,bb,p,q,i1,j,n)
                bb(i1:n)=a(i1:n,j)
                a(i1:n,j )=p*bb(i1:n)-a(i1:n,i1)
                a(i1:n,i1)=bb(i1:n)-q*a(i1:n,j )
c!$OMP END PARALLEL WORKSHARE
              endif
            else
              c=1.d0
              s=0.d0
            endif
            a(i,j)=s
            if(j < ndim+2)then
c            if(j < n+2)then
              a(j-1,i)=c
            else
              if(c <= abs(s) .and. c /= 0.d0)then
                a(i,j)=sign(1.d0/c,s)
              endif
            endif
          enddo
        endif
      enddo
      anorm=0.d0
      xmin=1.d38
      do i=1,mn
        p=sqrt(v(i))
        v(i+mn)=x(i)
        a(i,i)=a(i,i)*p
        x(i)=a(i,i)*x(i)
        b(i,1:n1)=b(i,1:n1)*p
        if(i < m)then
          a(i,i+1)=a(i,i+1)*p
          v(i)=a(i,i+1)*x(i+1)
        else
          v(i)=0.d0
        endif
        an=abs(x(i))+abs(v(i))
        anorm=max(anorm,an)
        if(an /= 0.d0)then
          xmin=min(xmin,abs(x(i)))
        endif
      enddo
      v(2*mn+1:m)=x(mn+1:m)
      do i=min(mn,m-2),1,-1
        i1=i+1
        i1mn=i1+mn
        do j=m,i+2,-1
          s=a(i,j)
          if(j < ndim+2)then
c          if(j < n+2)then
            c=a(j-1,i)
            a(j-1,i)=0.d0
          else
            if(abs(s) > 1.d0)then
              c=abs(1.d0/s)
              s=sign(1.d0+sqrt1(-c**2),s)
c              s=sign(1.d0-c**2/(1.d0+sqrt((1.d0-c)*(1.d0+c))),s)
            else
              c=1.d0+sqrt1(-s**2)
c              c=1.d0-s**2/(1.d0+sqrt((1.d0-s)*(1.d0+s)))
            endif
          endif
          if(abs(c) > abs(s))then
            h1=v(i1mn)*c
            h2=v(j +mn)/c
            p=s*v(j +mn)/h1
            q=s*v(i1mn)/h2
            v(i1mn)=h1
            v(j +mn)=h2
            a(i,j)=q*a(i,i1)
c            call omp_set_num_threads(int(1+(mn-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,p,q,i1,mn,j)
            a(i1:mn,i1)=a(i1:mn,i1)-p*a(i1:mn,j )
            a(i1:mn,j )=a(i1:mn,j )+q*a(i1:mn,i1)
c!$OMP END PARALLEL WORKSHARE 
          else
            h1=v(j +mn)/s
            h2=v(i1mn)*s
            p=c*v(j +mn)/h2
            q=c*v(i1mn)/h1
            v(i1mn)=h1
            v(j +mn)=h2
            a(i,j )=a(i,i1)
            a(i,i1)=a(i,i1)*q
c            call omp_set_num_threads(int(1+(mn-i)/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,aa,p,q,i1,mn,j)
            aa(i1:mn)=a(i1:mn,j)
            a(i1:mn,j )=p*aa(i1:mn)+a(i1:mn,i1)
            a(i1:mn,i1)=q*a(i1:mn,j )-aa(i1:mn)
c!$OMP END PARALLEL WORKSHARE
          endif
        enddo
      enddo
      lsep(0)=0
      lsep(1)=1
      isep=1
      ibegin=1
      iend=mn
      v(mn*2+1:mn)=1.d0
      v(0)=0.d0
 1002 if(v(iend) /= 0.d0)then
        f=v(iend)
        v(iend)=0.d0
        do i=iend,ibegin,-1
          if(abs(f)+abs(x(i)) /= abs(x(i)))then
            p=hypot(x(i),f)
            vv=v(i-1)/p
            v(i-1)=vv*x(i)
            x(i  )=p
            f=-vv*f
          else
            exit
          endif
        enddo
      endif
      do3001: do
c      write(*,'(a,2i5,1p10g12.4)')'tsvdm-5 ',ibegin,iend,a(1,11:20)
        do while(ibegin >= iend)
          isep=isep-1
          if(isep <= 0)then
            exit do3001
          endif
          iend=ibegin-1
          ibegin=lsep(isep)
        enddo
        do it=1,itmax
          if(x(iend) == 0.d0)then
            iend=iend-1
            go to 1002
          endif
          do i=iend-1,ibegin,-1
c            an=max(abs(x(i)),abs(x(i+1)))
            an=abs(x(i))+abs(x(i+1))
            if(abs(v(i)) <= an*1.d-16)then
              v(i)=0.d0
              if(i == iend-1)then
                iend=iend-1
              else
                isep=isep+1
                ibegin=i+1
                lsep(isep)=ibegin
              endif
              cycle do3001
            endif
          enddo
          w=x(ibegin)
          z=x(iend)
          y=x(iend-1)
          if(ibegin < iend-1)then
            g=v(iend-2)
          else
            g=0.d0
          endif
          h=v(iend-1)
          f=((y-z)*(y+z)+(g-h)*(g+h))*.5d0
          if(w == 0.d0)then
            f=0.d0
          else
            g=h*y
            y=f+sign(hypot(f,g),f)
            if(y /= 0.d0)then
              f=((w-z)*(w+z)-h**2+g**2/y)/w
            else
              f=((w-z)*(w+z)-h**2)/w
            endif
          endif
          g=v(ibegin)
          h=g
          do1221: do kkk=1,1
            do i=ibegin,iend-1
              i1=i+1
              z=hypot(f,g)
              v(i-1)=z
              if(z /= 0.d0)then
                c=f/z
                s=g/z
              else
                x(i)=w
                v(i)=h
                exit do1221
              endif
              f= w*c+h*s
              g=-w*s+h*c
              h= x(i1)*s
              y= x(i1)*c
              z=hypot(f,h)
              x(i)=z
              if(z /= 0.d0)then
                c=f/z
                s=h/z
              else
                v(i)=g
                x(i1)=y
                exit do1221
              endif
              f= c*g+s*y
              w=-s*g+c*y
              g=v(i1)*s
              h=v(i1)*c
              if(abs(c) > abs(s))then
                h1=v(i+mn)/c
                h2=v(i1+mn)*c
                r=s*v(i+mn)/h2
                t=s*v(i1+mn)/h1
                v(i+mn)=h1
                v(i1+mn)=h2
c                call omp_set_num_threads(int(1+mn/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,b,r,t,m,n1,i,i1)
                a(i1,1:m)=a(i1,1:m)-r*a(i ,1:m)
                a(i ,1:m)=a(i ,1:m)+t*a(i1,1:m)
c!$OMP END PARALLEL WORKSHARE
                b(i1,1:n1)=b(i1,1:n1)-r*b(i,1:n1)
                b(i ,1:n1)=b(i ,1:n1)+t*b(i1,1:n1)
              else
                h1=v(i1+mn)/s
                h2=v(i+mn)*s
                r=c*v(i1+mn)/h2
                t=c*v(i+mn)/h1
                v(i+mn)=h1
                v(i1+mn)=h2
c                call omp_set_num_threads(int(1+mn/nomp))
c!$OMP PARALLEL WORKSHARE shared(aa,bb,a,b,r,t,m,n1,i,i1)
                aa(1:m)=a(i1,1:m)
                a(i1,1:m)=r*aa(1:m)-a(i ,1:m)
                a(i ,1:m)=aa(1:m)-t*a(i1,1:m)
c!$OMP END PARALLEL WORKSHARE
                bb(1:n1)=b(i1,1:n1)
                b(i1,1:n1)=r*bb(1:n1)-b(i,1:n1)
                b(i ,1:n1)=bb(1:n1)-t*b(i1,1:n1)
              endif
            enddo
            v(iend-1)=f
            x(iend)=w
          enddo do1221
          v(ibegin-1)=0.d0
        enddo
        write(*,*)' TSVDM convergence fail. ',iend,v(iend-1),
     $       x(iend-1),x(iend)
        v(iend)=0.d0
        if(ktfenanq(x(iend)))then
          x(iend)=0.d0
        endif
        iend=iend-1
        if(ktfenanq(x(iend)))then
          x(iend)=0.d0
        endif
        if(ktfenanq(v(iend)))then
          v(iend)=0.d0
        endif
      enddo do3001
      anorm=abs(x(1))
      do i=2,mn
        anorm=max(anorm,abs(x(i)))
      enddo
      anorm=anorm*epslon
      if(inv)then
        do i=1,mn
          s=abs(x(i))
          f=v(i+mn)
          if(s > anorm)then
            x(i)=1.d0/s
            w=f*x(i)
          else
            x(i)=s/anorm**2
            w=f
          endif
c          call omp_set_num_threads(int(1+mn/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,b,m,w,f)
          a(i,1:m)=a(i,1:m)*w
          b(i,1:n1)=b(i,1:n1)*f
c!$OMP END PARALLEL WORKSHARE
        enddo
      else
        do i=1,mn
          s=abs(x(i))
          f=v(i+mn)
          if(s > anorm)then
            w=f/s
          else
            s=0.d0
            w=f
          endif
          x(i)=s
c          call omp_set_num_threads(int(1+mn/nomp))
c!$OMP PARALLEL WORKSHARE shared(a,b,m,w,f)
          a(i,1:m)=a(i,1:m)*w
          b(i,1:n1)=b(i,1:n1)*f
c!$OMP END PARALLEL WORKSHARE
        enddo
      endif
      x(mn+1:m)=0.d0
      a(mn+1:n,1:m)=0.d0
      return
      end
