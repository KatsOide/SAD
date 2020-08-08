      subroutine tsolvg(a,b,x,n,m,ndim)
      implicit none
      integer*4 n,m,ndim
      real*8 a(ndim,m),b(n),x(m)
      call tsvd(a,b,x,n,m,ndim,1.d-8,.false.)
      return
      end

      subroutine tsolva(a,b,x,n,m,ndim,epslon)
      implicit none
      integer*4 ,intent(in):: n,m,ndim
      real*8 ,intent(inout):: a(ndim,m),b(n)
      real*8 ,intent(out):: x(m)
      real*8 ,intent(in):: epslon
      call tsvd(a,b,x,n,m,ndim,epslon,.false.)
      return
      end

      subroutine tsvd(a,b,x,n,m,ndim,epslon,svd)
c
c   Singular Value Decomposition
c
c   a(ndim,m)   array to contain the matrix a(n,m).
c   b(n)        vector of r.h.s.
c   x(m)        unknown variable.
c   epslon      threshold.
c   svd         when .true., does complete svd and returns
c                 V in a, W^-1 in b, while a is decomposed
c                 a=U^T.W.V .
c               when .false. a quicker routine is used to obtain x.
c                 a and b do not return meaningful things.
c
      implicit none
      integer*4 ,intent(in):: n,m,ndim
      integer*4 , parameter ::nmax=100000,itmax=256
      real*8 ,intent(inout):: a(ndim,m),b(n)
      real*8 ,intent(out):: x(m)
      real*8 ,intent(in):: epslon
      real*8 v(0:nmax),anorm,enorm
      real*8 aa,f,g,s,r,w,u,h,xmin,z,vv,d,c,p,bb,y,an
      real*8 q,h1,h2,t,r1,r2,ra,aam(m),bbm(n)
      integer*4 lsep(0:nmax),i,j,k,kkk,nfail,
     $     mn,it,isep,ibegin,iend,ma,i1,i1mn
      logical*4 svd
c     begin initialize for preventing compiler warning
      enorm=0.d0
c     end   initialize for preventing compiler warning

      nfail=4
      mn=min(n,m)
      if(max(mn+m,n) .gt. nmax)then
        write(*,*)' TSVD Too large matrix. ',n,m
        return
      endif
      v(1:n)=1.d0
      x=1.d0
      do 10 i=1,mn
        i1=i+1
        if(i .lt. n)then
          do 5110 j=i1,n
            if(abs(a(i,i)) .gt. abs(a(j,i)))then
              p=a(j,i)/a(i,i)
              h1=v(i)+v(j)*p**2
              q=v(j)*p/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
c              do 5120 k=i1,m
                a(j,i1:m)=a(j,i1:m)-p*a(i,i1:m)
                a(i,i1:m)=a(i,i1:m)+q*a(j,i1:m)
c 5120         continue
              a(j,i)=0.d0
              b(j)=b(j)-p*b(i)
              b(i)=b(i)+q*b(j)
            elseif(a(j,i) .ne. 0.d0)then
              p=a(i,i)/a(j,i)
              h1=v(j)+v(i)*p**2
              q=v(i)*p/h1
              v(j)=v(i)*v(j)/h1
              v(i)=h1
c              do 5130 k=i1,m
                aam(i1:m)=a(j,i1:m)
                a(j,i1:m)=p*aam(i1:m)-a(i,i1:m)
                a(i,i1:m)=aam(i1:m)-q*a(j,i1:m)
c 5130         continue
              a(i,i)=a(j,i)
              a(j,i)=0.d0
              bb=b(j)
              b(j)=p*bb-b(i)
              b(i)=bb-q*b(j)
            endif
 5110     continue
        endif
        if(i1 .lt. m)then
          do 5210 j=i1+1,m
            r1=x(i1)*a(i,i1)
            r2=x(j )*a(i,j )
            r=sign(abs(dcmplx(r1,r2)),r1)
            if(r .ne. 0.d0)then
              c=r1/r
              s=r2/r
              if(abs(c) .gt. abs(s))then
                a(i,j)=0.d0
                h1=x(i1)/c
                h2=x(j )*c
                p=s*x(i1)/h2
                q=s*x(j )/h1
                x(i1)=h1
                x(j )=h2
c                do 5220 k=i1,n
                  a(i1:n,j )=a(i1:n,j )-p*a(i1:n,i1)
                  a(i1:n,i1)=a(i1:n,i1)+q*a(i1:n,j )
c 5220           continue
              else
                a(i,i1)=a(i,j)
                a(i,j)=0.d0
                h1=x(j )/s
                h2=x(i1)*s
                p=c*x(j )/h2
                q=c*x(i1)/h1
                x(i1)=h1
                x(j )=h2
c                do 5221 k=i1,n
                  bbm(i1:n)=a(i1:n,j)
                  a(i1:n,j )=p*bbm(i1:n)-a(i1:n,i1)
                  a(i1:n,i1)=bbm(i1:n)-q*a(i1:n,j )
c 5221           continue
              endif
            else
              c=1.d0
              s=0.d0
            endif
            a(i,j)=s
C     if(j .lt. ndim+2)then
            if(j .lt. n+2)then
              a(j-1,i)=c
            else
              if(c .le. abs(s) .and. c .ne. 0.d0)then
                a(i,j)=sign(1.d0/c,s)
              endif
            endif
 5210     continue
        endif
 10   continue
      anorm=0.d0
      xmin=1.d300
      do 20 i=1,mn
        p=sqrt(v(i))
        v(i+mn)=x(i)
        a(i,i)=a(i,i)*p
        x(i)=a(i,i)*x(i)
        b(i)=b(i)*p
        if(i .lt. m)then
          a(i,i+1)=a(i,i+1)*p
          v(i)=a(i,i+1)*x(i+1)
        else
          v(i)=0.d0
        endif
        an=abs(x(i))+abs(v(i))
        anorm=max(anorm,an)
        if(an .ne. 0.d0)then
          xmin=min(xmin,abs(x(i)))
        endif
 20   continue
      if(.not. svd)then
        enorm=anorm*epslon
        if(xmin .gt. enorm)then
          x(mn+2:m)=0.d0
          r=abs(dcmplx(x(mn),v(mn)))
          if(r .ne. 0.d0)then
            w=1.d0/r
            d=x(mn)*w
            b(mn)=b(mn)*w
            v(mn)=v(mn)*w
            if(mn .lt. m)then
              x(mn+1)=v(mn)*b(mn)
            endif
            x(mn)=d*b(mn)
            ma=m
          else
            ma=mn
            d=0.d0
            v(mn)=0.d0
            b(mn)=0.d0
            x(mn)=0.d0
            if(mn .lt. m)then
              x(mn+1)=0.d0
            endif
          endif
          do 4020 i=mn-1,1,-1
            if(x(i) .ne. 0.d0)then
              p=v(i)*d
              ra=abs(x(i))+abs(v(i))
              r=(x(i)/ra)**2+((v(i)-p)/ra)*((v(i)+p)/ra)
              w=1.d0/sqrt(r)/ra
              b(i)=(b(i)-p*b(i+1))*w
              u=-p*w
              v(i+1:min(mn,ma-1))=u*v(i+1:min(mn,ma-1))
              x(i+2:min(mn+1,ma))=x(i+2:min(mn+1,ma))
     $             +b(i)*v(i+1:min(mn,ma-1))
c              do 4030 k=i+1,min(mn,ma-1)
c                v(k)=u*v(k)
c                x(k+1)=x(k+1)+b(i)*v(k)
c 4030         continue
              v(i)=(v(i)-p*d)*w
              x(i+1)=x(i+1)+b(i)*v(i)
              d=x(i)*w
              x(i)=d*b(i)
            else
              ma=i
              d=0.d0
              v(i)=0.d0
              b(i)=0.d0
              x(i)=0.d0
            endif
 4020     continue
          do 4040 i=min(mn,m-2),1,-1
            do 4050 j=m,i+2,-1
              s=a(i,j)
C     if(j .lt. ndim+2)then
              if(j .lt. n+2)then
                c=a(j-1,i)
              else
                if(abs(s) .gt. 1.d0)then
                  c=abs(1.d0/s)
                  s=sign(1.d0-c**2/(1.d0+sqrt((1.d0-c)*(1.d0+c))),s)
                else
                  c=1.d0-s**2/(1.d0+sqrt((1.d0-s)*(1.d0+s)))
                endif
              endif
              aa=x(i+1)
              x(i+1)= c*aa-s*x(j)
              x(j  )= s*aa+c*x(j)
 4050       continue
 4040     continue
          return
        endif
      endif
      v(mn*2+1:mn+m)=x(mn+1:m)
c      do 5301 i=mn+1,m
c        v(i+mn)=x(i)
c 5301 continue
      do 5310 i=min(mn,m-2),1,-1
        i1=i+1
        i1mn=i1+mn
        do 5320 j=m,i+2,-1
          s=a(i,j)
C     if(j .lt. ndim+2)then
          if(j .lt. n+2)then
            c=a(j-1,i)
            a(j-1,i)=0.d0
          else
            if(abs(s) .gt. 1.d0)then
              c=abs(1.d0/s)
              s=sign(1.d0-c**2/(1.d0+sqrt((1.d0-c)*(1.d0+c))),s)
            else
              c=1.d0-s**2/(1.d0+sqrt((1.d0-s)*(1.d0+s)))
            endif
          endif
          if(abs(c) .gt. abs(s))then
            h1=v(i1mn)*c
            h2=v(j +mn)/c
            p=s*v(j +mn)/h1
            q=s*v(i1mn)/h2
            v(i1mn)=h1
            v(j +mn)=h2
            a(i,j)=q*a(i,i1)
            a(i1:mn,i1)=a(i1:mn,i1)-p*a(i1:mn,j )
            a(i1:mn,j )=a(i1:mn,j )+q*a(i1:mn,i1)
c            do 5330 k=i1,mn
c              a(k,i1)=a(k,i1)-p*a(k,j )
c              a(k,j )=a(k,j )+q*a(k,i1)
c 5330       continue
          else
            h1=v(j +mn)/s
            h2=v(i1mn)*s
            p=c*v(j +mn)/h2
            q=c*v(i1mn)/h1
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
      lsep(0)=0
      lsep(1)=1
      isep=1
      ibegin=1
      iend=mn
      v(mn+1:mn*2)=1.d0
c      do 1510 i=1,mn
c        v(mn+i)=1.d0
c 1510 continue
      v(0)=0.d0
 1002 if(v(iend) .ne. 0.d0)then
        f=v(iend)
        v(iend)=0.d0
        do 1110 i=iend,ibegin,-1
          if(abs(f)+abs(x(i)) .ne. abs(x(i)))then
            p=abs(dcmplx(x(i),f))
            vv=v(i-1)/p
            v(i-1)=vv*x(i)
            x(i  )=p
            f=-vv*f
          else
            exit
          endif
 1110   continue
      endif
      do3001: do while(.true.)
        if(ibegin .ge. iend)then
          isep=isep-1
          if(isep .le. 0)then
            exit
          endif
          iend=ibegin-1
          ibegin=lsep(isep)
          cycle
        endif
        do 1210 it=1,itmax
c     write(*,*)it,ibegin,iend,v(iend-1)
          if(x(iend) .eq. 0.d0)then
            iend=iend-1
            go to 1002
          endif
          do 1010 i=iend-1,ibegin,-1
c            an=max(abs(x(i)),abs(x(i+1)))
            an=abs(x(i))+abs(x(i+1))
            if(abs(v(i)) .le. an*1.d-16)then
              v(i)=0.d0
              if(i .eq. iend-1)then
                iend=iend-1
              else
                isep=isep+1
                ibegin=i+1
                lsep(isep)=ibegin
              endif
              cycle do3001
            endif
 1010     continue
          if(.not. svd)then
            do1611: do kkk=1,1
              do 1610 i=1,mn
                if(abs(x(i)) .lt. enorm)then
                  exit do1611
                endif
 1610         continue
              do 1710 i=iend,ibegin+1,-1
                p=v(i+mn)*v(i-1)/v(i-1+mn)/x(i)
                a(i-1,1:m)=a(i-1,1:m)-p*a(i,1:m)
                b(i-1)=b(i-1)-p*b(i)
                v(i-1)=0.d0
 1710         continue
              iend=ibegin-1
              cycle do3001
            enddo do1611
          endif
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
            do i=ibegin,iend-1
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
c                do 1150 k=1,m
                  a(i1,1:m)=a(i1,1:m)-r*a(i ,1:m)
                  a(i ,1:m)=a(i ,1:m)+t*a(i1,1:m)
c 1150           continue
                b(i1)=b(i1)-r*b(i)
                b(i )=b(i )+t*b(i1)
              else
                h1=v(i1+mn)/s
                h2=v(i+mn)*s
                r=c*v(i1+mn)/h2
                t=c*v(i+mn)/h1
                v(i+mn)=h1
                v(i1+mn)=h2
c                do 1151 k=1,m
                  aam(1:m)=a(i1,1:m)
                  a(i1,1:m)=r*aam(1:m)-a(i ,1:m)
                  a(i ,1:m)=aam(1:m)-t*a(i1,1:m)
c 1151           continue
                bb=b(i1)
                b(i1)=r*bb-b(i)
                b(i )=bb-t*b(i1)
              endif
            enddo
            v(iend-1)=f
            x(iend)=w
          enddo do1221
          v(ibegin-1)=0.d0
 1210   continue
        if(nfail .ge. 0)then
          write(*,'(a,i10,1p3g15.7)')' TSVD convergence failed: ',
     $         iend,v(iend-1),x(iend-1),x(iend)
          if(nfail .eq. 0)then
            write(*,*)' -- further message will be suppressed.'
          endif
          nfail=nfail-1
        endif
        iend=iend-1
        v(iend)=0.d0
      enddo do3001
      anorm=maxval(abs(x(1:mn)))
c      anorm=abs(x(1))
c      do 3010 i=2,mn
c        anorm=max(anorm,abs(x(i)))
c 3010 continue
      anorm=anorm*epslon
      if(svd)then
        do concurrent (i=1:mn)
          s=x(i)
          if(abs(s) .gt. anorm)then
            w=v(i+mn)/s
          else
            w=s*v(i+mn)/anorm**2
          endif
          v(i)=w/v(i+mn)
          b(i)=b(i)*w
c          do j=1,m
            a(i,1:m)=a(i,1:m)*w
c          enddo
        enddo
        do concurrent (i=1:m)
          s=sum(a(1:mn,i)*b(1:mn))
c          s=a(1,i)*b(1)
c          do j=2,mn
c            s=s+a(j,i)*b(j)
c          enddo
          x(i)=s
        enddo
        b(1:mn)=v(1:mn)
        b(mn+1:m)=0.d0
c        do i=1,m
        a(mn+1:n,1:m)=0.d0
c        enddo
      else
        do concurrent (i=1:mn)
          s=x(i)
          if(abs(s) .gt. anorm)then
            w=(v(i+mn)/s)**2
          else
            w=(s*v(i+mn)/anorm**2)**2
          endif
          b(i)=b(i)*w
        enddo
        do concurrent (i=1:m)
c          s=a(1,i)*b(1)
c          do 3040 j=2,mn
c            s=s+a(j,i)*b(j)
c 3040     continue
          x(i)=sum(a(1:mn,i)*b(1:mn))
        enddo
      endif
      return
      end
