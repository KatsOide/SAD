c      complex*16 a(3,3),b(3,3)
c      real*8 x(3)
c      a(1,1)=1.d0
c      a(1,2)=(0.d0,1.d0)
c      a(1,3)=(1.d0,1.d0)
c      a(2,1)=(0.d0,1.d0)
c      a(2,2)=(0.d0,2.d0)
c      a(2,3)=(1.d0,0.d0)
c      a(3,1)=(1.d0,1.d0)
c      a(3,2)=(1.d0,2.d0)
c      a(3,3)=(1.d0,1.d0)
c      call tcsvdm(a,b,x,3,3,3,3,1.d-10,.false.)
c      write(*,*)a,b,x
c      stop
c      end
c
c


      subroutine tcsvdm(a,b,x,n,m,ndim,ndimb,epslon,inv)
c
c   Singular Value Decomposition for Complex Matrix
c
c   a(ndim,m)   array contains the matrix a(n,m).
c   b(ndimb,n)  array to contain U .
c   x(m)        vector to contain (if inv, W^-1 else W)
c   epslon      threshold.
c               V is returned in a .
c               a=Conjugate[U^T].W.V .
c    
c
      use tfstk, only : forcesf
      implicit none
      integer*4 nmax,itmax
      integer*4 ,intent(in)::n,m,ndim,ndimb
      parameter (nmax=100000,itmax=256)
      complex*16 ,intent(inout)::a(ndim,m),b(ndimb,n)
      complex*16 p,cp,q,zc,zs,a1,aa,bb,r1,r2
      real*8 ,intent(out)::x(m)
      real*8 ,allocatable::v(:)
      real*8 anorm
      real*8 ,intent(in)::epslon
      real*8 g,s,r,w,h,xmin,z,c,y,an,aa1,r1a,r2a,ap,f
      real*8 h1,h2,t,vv
      integer*4 i,j,k,mn,it,isep,ibegin,iend,i1,i1mn,n1,kkk
      integer*4 ,allocatable::lsep(:)
      logical*4 ,intent(in)::inv
      mn=min(n,m)
      if(max(mn+m,n) .gt. nmax)then
        write(*,*)' TCSVDM Too large matrix. ',n,m
        return
      endif
      allocate (v(0:nmax),lsep(nmax))
      n1=min(ndimb,n)
      do 1 i=1,n
        v(i)=1.d0
 1    continue
      do 2 i=1,m
        x(i)=1.d0
 2    continue
      do j=1,n1
        do i=1,n1
          b(i,j)=0.d0
        enddo
        b(j,j)=1.d0
      enddo
      do 10 i=1,mn
        i1=i+1
        if(i .lt. n)then
          do 5110 j=i1,n
            if(abs(a(i,i)) .gt. abs(a(j,i)))then
              p=a(j,i)/a(i,i)
              cp=conjg(p)
              h1=v(i)+v(j)*dble(p*cp)
              q=v(j)*cp/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
              do 5120 k=i1,m
                a(j,k)=a(j,k)-p*a(i,k)
                a(i,k)=a(i,k)+q*a(j,k)
 5120         continue
              a(j,i)=0.d0
              do k=1,n1
                b(j,k)=b(j,k)-p*b(i,k)
                b(i,k)=b(i,k)+q*b(j,k)
              enddo
            elseif(a(j,i) .ne. 0.d0)then
              p=a(i,i)/a(j,i)
              cp=conjg(p)
              h1=v(j)+v(i)*dble(p*cp)
              q=v(i)*cp/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
              do 5130 k=i1,m
                aa=a(j,k)
                a(j,k)=p*aa-a(i,k)
                a(i,k)=aa-q*a(j,k)
 5130         continue
              a(i,i)=a(j,i)
              a(j,i)=0.d0
              do k=1,n1
                bb=b(j,k)
                b(j,k)=p*bb-b(i,k)
                b(i,k)=bb-q*b(j,k)
              enddo
            endif
 5110     continue
        endif
        if(i1 .lt. m)then
          do 5210 j=i1+1,m
            r1=x(i1)*a(i,i1)
            r2=x(j )*a(i,j )
            r1a=abs(r1)
            r2a=abs(r2)
            r=hypot(r1a,r2a)
c            r=sqrt(r1a**2+r2a**2)
            if(r .ne. 0.d0)then
              if(r1a .gt. r2a)then
                c=r1a/r
                zc=c
                zs=r2*c/r1
                a(i,j)=0.d0
                h1=x(i1)/c
                h2=x(j )*c
                p=zs*x(i1)/h2
                q=conjg(zs)*x(j )/h1
                x(i1)=h1
                x(j )=h2
                do 5220 k=i1,n
                  a(k,j )=a(k,j )-p*a(k,i1)
                  a(k,i1)=a(k,i1)+q*a(k,j )
 5220           continue
              else
                s=r2a/r
                zc=r1*s/r2
                zs=s
                a(i,i1)=a(i,j)
                a(i,j)=0.d0
                h1=x(j )/s
                h2=x(i1)*s
                p=zc*x(j )/h2
                q=conjg(zc)*x(i1)/h1
                x(i1)=h1
                x(j )=h2
                do 5221 k=i1,n
                  aa=a(k,j)
                  a(k,j )=p*aa-a(k,i1)
                  a(k,i1)=aa-q*a(k,j )
 5221           continue
              endif
            else
              zc=1.d0
              zs=0.d0
            endif
            a(i,j)=zs
            if(j .lt. n+2)then
              a(j-1,i)=zc
            elseif(r .ne. 0.d0 .and. r1a .le. r2a .and.
     $             zc .ne. 0.d0)then
              a(i,j)=1.d0/zc
            endif
 5210     continue
        endif
 10   continue
      anorm=0.d0
      xmin=1.d38
      a1=a(1,1)
      do 20 i=1,mn
        ap=sqrt(v(i))
        v(i+mn)=x(i)
        aa1=abs(a1)
        if(aa1 .ne. 0.d0)then
          p=ap*conjg(a1)/aa1
        else
          p=1.d0
        endif
        a(i,i)=a(i,i)*p
        x(i)=abs(a(i,i))*x(i)
c        x(i)=aa1*ap*x(i)
        do k=1,n1
          b(i,k)=b(i,k)*p
        enddo
        if(i .lt. m)then
          a(i,i+1)=a(i,i+1)*p
          aa1=abs(a(i,i+1))
          v(i)=aa1*x(i+1)
          if(aa1 .ne. 0.d0)then
            a1=a(i+1,i+1)*conjg(a(i,i+1))/aa1
          endif
        else
          v(i)=0.d0
        endif
        an=abs(x(i))+abs(v(i))
        anorm=max(anorm,an)
        if(an .ne. 0.d0)then
          xmin=min(xmin,abs(x(i)))
        endif
 20   continue
      do 5301 i=mn+1,m
        v(i+mn)=x(i)
 5301 continue
      do 5310 i=min(mn,m-2),1,-1
        i1=i+1
        i1mn=i1+mn
        do 5320 j=m,i+2,-1
          zs=a(i,j)
          s=abs(zs)
          if(j .lt. n+2)then
            zc=a(j-1,i)
            c=abs(zc)
            a(j-1,i)=0.d0
          else
            if(s .gt. 1.d0)then
              zc=1.d0/zs
              c=abs(zc)
              s=1.d0-c**2/(1.d0+sqrt((1.d0-c)*(1.d0+c)))
            else
              c=1.d0-s**2/(1.d0+sqrt((1.d0-s)*(1.d0+s)))
c     This code path MIGHT satisfy c > s
              if(.not. (c .gt. s))then
                write(*,*)'tcsvdm @ src/tcsvdm.f: ',
     $                'Identity equation c > s ',
     $                'for j >= n+2 && 1 >= s is broken!',
     $                '(FIXME)'
                call abort
              endif
c     begin initialize for preventing compiler warning
              zc=0.d0
c     end   initialize for preventing compiler warning
            endif
          endif
          if(c .gt. s)then
            h1=v(i1mn)*c
            h2=v(j +mn)/c
            p=conjg(zs)*v(j +mn)/h1
            q=zs*v(i1mn)/h2
            v(i1mn)=h1
            v(j +mn)=h2
            a(i,j)=q*a(i,i1)
            do 5330 k=i1,mn
              a(k,i1)=a(k,i1)-p*a(k,j )
              a(k,j )=a(k,j )+q*a(k,i1)
 5330       continue
          else
            h1=v(j +mn)/s
            h2=v(i1mn)*s
            p=conjg(zc)*v(j +mn)/h2
            q=zc*v(i1mn)/h1
            v(i1mn)=h1
            v(j +mn)=h2
            a(i,j )=a(i,i1)
            a(i,i1)=a(i,i1)*q
            do 5331 k=i1,mn
              aa=a(k,j)
              a(k,j )=p*aa+a(k,i1)
              a(k,i1)=q*a(k,j )-aa
 5331       continue
          endif
 5320   continue
 5310 continue
c     write(*,9710)(v(i+mn)-1.d0,i=1,m)
c9710 format(1x,:1p11g11.3)
      lsep(1)=1
      isep=1
      ibegin=1
      iend=mn
      do 1510 i=1,mn
        v(mn+i)=1.d0
 1510 continue
      v(0)=0.d0
 1002 if(v(iend) .ne. 0.d0)then
        f=v(iend)
        v(iend)=0.d0
        do 1110 i=iend,ibegin,-1
          if(abs(f)+abs(x(i)) .ne. abs(x(i)))then
            ap=hypot(abs(f),x(i))
c            ap=sqrt(abs(f)**2+x(i)**2)
            vv=v(i-1)/ap
            v(i-1)=vv*x(i)
            x(i  )=ap
            f=-vv*f
          else
            exit
          endif
 1110   continue
      endif
      do3001: do while(.true.)
 1001   if(ibegin .ge. iend)then
          isep=isep-1
          if(isep .le. 0)then
            exit
          endif
          iend=ibegin-1
          ibegin=lsep(isep)
          go to 1001
        endif
        do 1210 it=1,itmax
c     write(*,*)it,ibegin,iend,v(iend-1)
          if(x(iend) .eq. 0.d0)then
            iend=iend-1
            go to 1002
          endif
          do 1010 i=iend-1,ibegin,-1
            an=max(abs(x(i)),abs(x(i+1)))
            if(abs(v(i))+an .eq. an)then
              v(i)=0.d0
              if(i .eq. iend-1)then
                iend=iend-1
              else
                isep=isep+1
                ibegin=i+1
                lsep(isep)=ibegin
              endif
              go to 1001
            endif
 1010     continue
          w=x(ibegin)
          z=x(iend)
          y=x(iend-1)
          if(ibegin .lt. iend-1)then
            g=v(iend-2)
          else
            g=0.d0
          endif
          h=v(iend-1)
          f=((y-z)*(y+z)+(g-h)*(g+h))*.5d0
          if(w .eq. 0.d0)then
            f=0.d0
          else
            g=h*y
            y=f+sign(abs(dcmplx(f,g)),f)
            if(y .ne. 0.d0)then
              f=((w-z)*(w+z)-h**2+g**2/y)/w
            else
              f=((w-z)*(w+z)-h**2)/w
            endif
          endif
          g=v(ibegin)
          h=g
          do1221: do kkk=1,1
            do 1140 i=ibegin,iend-1
              i1=i+1
              z=abs(dcmplx(f,g))
              v(i-1)=z
              if(z .ne. 0.d0)then
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
              z=abs(dcmplx(f,h))
              x(i)=z
              if(z .ne. 0.d0)then
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
              if(abs(c) .gt. abs(s))then
                h1=v(i+mn)/c
                h2=v(i1+mn)*c
                r=s*v(i+mn)/h2
                t=s*v(i1+mn)/h1
                v(i+mn)=h1
                v(i1+mn)=h2
                do 1150 k=1,m
                  a(i1,k)=a(i1,k)-r*a(i ,k)
                  a(i ,k)=a(i ,k)+t*a(i1,k)
 1150           continue
                do k=1,n1
                  b(i1,k)=b(i1,k)-r*b(i,k)
                  b(i ,k)=b(i ,k)+t*b(i1,k)
                enddo
              else
                h1=v(i1+mn)/s
                h2=v(i+mn)*s
                r=c*v(i1+mn)/h2
                t=c*v(i+mn)/h1
                v(i+mn)=h1
                v(i1+mn)=h2
                do 1151 k=1,m
                  aa=a(i1,k)
                  a(i1,k)=r*aa-a(i ,k)
                  a(i ,k)=aa-t*a(i1,k)
 1151           continue
                do k=1,n1
                  bb=b(i1,k)
                  b(i1,k)=r*bb-b(i,k)
                  b(i ,k)=bb-t*b(i1,k)
                enddo
              endif
 1140       continue
            v(iend-1)=f
            x(iend)=w
          enddo do1221
          v(ibegin-1)=0.d0
 1210   continue
        write(*,*)' TCSVDM convergence fail. ',iend
        iend=iend-1
        v(iend)=0.d0
        go to 1001
      enddo do3001
      anorm=abs(x(1))
      do 3010 i=2,mn
        anorm=max(anorm,abs(x(i)))
 3010 continue
      anorm=anorm*epslon
      if(inv)then
        do i=1,mn
          s=abs(x(i))
          if(s .gt. anorm)then
            x(i)=1.d0/s
            w=v(i+mn)*x(i)
          else
            x(i)=s/anorm**2
            w=v(i+mn)
          endif
          do j=1,m
            a(i,j)=a(i,j)*w
          enddo
          do k=1,n1
            b(i,k)=b(i,k)*v(i+mn)
          enddo
        enddo
      else
        do i=1,mn
          s=abs(x(i))
          if(s .gt. anorm)then
            w=v(i+mn)/s
          else
            s=0.d0
            w=v(i+mn)
          endif
          x(i)=s
          do j=1,m
            a(i,j)=a(i,j)*w
          enddo
          do k=1,n1
            b(i,k)=b(i,k)*v(i+mn)
          enddo
        enddo
      endif
      do j=mn+1,m
        x(i)=0.d0
      enddo
      do i=1,m
        do j=mn+1,n
          a(j,i)=0.d0
        enddo
      enddo
      deallocate (v,lsep)
      return
      end
