      module eigen
      real*8 , allocatable::vx(:)

      contains

      subroutine tbal(w,a,n,ndim)
      implicit none
c   radix is the machine radix of floating.
      real*8 , parameter ::radix=2.d0,sqrdx=radix**2,rth=0.95d0
      integer*4 n,ndim,i,j
      real*8 a(n,n),w(ndim,n)
      real*8 c,r,g,f,s,u
      logical*4 last
      a=w(1:n,:)
      w(1:n,:)=0.d0
      vx=1.d0
      do 1 i=1,n
        w(i,i)=1.d0
1     continue
      last=.true.
      do while(last)
        last=.false.
        do i=1,n
          c=0.d0
          r=0.d0
          do j=1,n
            if(j .ne. i)then
              u=vx(j)/vx(i)
              c=c+abs(a(j,i)*u)
              r=r+abs(a(i,j)/u)
            endif
          enddo
          if(c .ne. 0.d0 .and. r .ne. 0.d0)then
            g=r/radix
            f=1.d0
            s=c+r
            do while(c .lt. g)
              f=f*radix
              c=c*sqrdx
            enddo
            g=r*radix
            do while(c .gt. g)
              f=f/radix
              c=c/sqrdx
            enddo
            if((c+r)/f .lt. rth*s)then
              last=.true.
              vx(i)=vx(i)/f
            endif
          endif
        enddo
      enddo
      return
      end subroutine

      subroutine thess(a,w,n,ndim)
      implicit none
      integer*4 n,ndim,i,i1,j
      real*8 a(n,n),w(ndim,n),aa(n)
      real*8 h1,p,q
      vx(2:n)=vx(2:n)**2
      do 10 i=1,n-2
        i1=i+1
        do 20 j=i1+1,n
          if(a(j,i) .ne. 0.d0)then
            if(abs(a(i1,i)) .gt. abs(a(j,i)))then
              p=a(j,i)/a(i1,i)
              h1=vx(i1)+vx(j)*p**2
              q=vx(j)*p/h1
              vx(j)=vx(i1)*vx(j)/h1
              vx(i1)=h1
              a(j,i)=0.d0
              a(j ,i1:n)=a(j ,i1:n)-p*a(i1,i1:n)
              a(i1,i1:n)=a(i1,i1:n)+q*a(j ,i1:n)
              a(:,i1)=a(:,i1)+p*a(:,j )
              a(:,j )=a(:,j )-q*a(:,i1)
              w(1:n,i1)=w(1:n,i1)+p*w(1:n,j )
              w(1:n,j )=w(1:n,j )-q*w(1:n,i1)
            else
              p=a(i1,i)/a(j,i)
              h1=vx(j)+vx(i1)*p**2
              q=vx(i1)*p/h1
              vx(j)=vx(i1)*vx(j)/h1
              vx(i1)=h1
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
      vx(2:n)=sqrt(vx(2:n))
      return
      end subroutine

      subroutine tqr(a,w,eig,n,ndim)
      use tfstk, only: ktfenanq
      implicit none
      integer*4 , parameter ::itmax=100
      real*8 ,parameter ::vmax=1.d10,vmin=1.d0/vmax,
     $     decth=1.d-16,threj=1.d-8,alpha=1.d-6,
     $     amth=1.d-30
      integer*4 , intent(in)::n,ndim
      integer*4 ibp,ib,ie,i,j,ie1,iter,is,itm,
     $     is1,is2,i1,i2,i3,jm,j1,ii
      integer*4 ibtab(n)
      real*8 , intent(inout)::a(n,n)
      real*8 , intent(out):: w(ndim,n),eig(2,n)
      real*8 aa(n),aa1(n)
      real*8 px(2),py(2)
      real*8 anorm,c,s,p,q,v1,v2,am,am1,w1,w2,w3,
     $     dec,sqrd,r,da,pk,pl,u,v,pm,pn,a1,a2,x,y,ee,sa,
     $     r1,r2,s1
      logical*4 paired,jordan
      complex*16 ca,cc,cu,cr,ck,cm
      equivalence (ck,px),(cm,py)
      equivalence (px,pk),(px(2),pl),(py,pm),(py(2),pn)
c     begin initialize for preventing compiler warning
      is2=0
      w1=0.d0
      w2=0.d0
      w3=0.d0
c     end   initialize for preventing compiler warning
      jordan=.false.
      ibp=1
      ibtab(1)=1
      ib=1
      ie=n
      do i=1,ie-2
        a(i+2,i)=0.d0
        a(min(i+3,ie),i)=0.d0
      enddo
      anorm=0.d0
      do i=1,n
        do j=max(1,i-1),n
          anorm=anorm+abs(vx(i)/vx(j)*a(i,j))
        enddo
      enddo
      anorm=anorm*2/n**2
      do1:do while(.true.)
        if(ie .le. ib+1)then
          if(ib .eq. 1)then
            exit
          else
            ie=ib-1
            ibp=ibp-1
            ib=ibtab(ibp)
            cycle do1
          endif
        endif
        ie1=ie-1
        iter=0
        itm=min(itmax,10*max(ie-ib+1,2))
        do2:do while(.true.)
c     write(*,'(1P4G15.7)')((vx(i)/vx(j)*a(i,j),j=1,4),i=1,4)
          do i=ie-1,ib,-1
            i1=i+1
            s=max(anorm,abs(a(i1,i1))+abs(a(i,i)))
            if(abs(vx(i1)/vx(i)*a(i1,i))+s .eq. s)then
              a(i1,i)=0.d0
              if(i1 .eq. ie)then
                ie=ie-1
                cycle do1
              endif
              ibp=ibp+1
              ibtab(ibp)=ib
              ib=i1
              cycle do1
            endif
          enddo
          do i=ib,ie
            if(abs(vx(i)) .lt. vmin  .or. abs(vx(i)) .gt. vmax)then
              j1=max(i-1,ib)
              a(i,j1:n)=a(i,j1:n)*vx(i)
              j1=min(ie,i+1)
              a(1:j1,i)=a(1:j1,i)/vx(i)
              w(1:n,i)=w(1:n,i)/vx(i)
              vx(i)=1.d0
            endif
          enddo
          do is=ie-2,ib,-1
            is1=is+1
            is2=is1+1
            a1=a(ie,ie)-a(is,is)
            a2=a(ie1,ie1)-a(is,is)
            w1=((a1*a2-a(ie1,ie)*a(ie,ie1))/a(is1,is)+a(is,is1))
     1           *vx(is)/vx(is1)
            if(w1 .eq. 0.d0)then
              w1=-alpha*a(ie1,ie)*a(ie,ie1)/a(is1,is)*vx(is)/vx(is1)
            endif
            w2=a(is1,is1)-a(is,is)-a1-a2
            w3=a(is2,is1)*vx(is2)/vx(is1)
            if(is .eq. ib)then
              exit
            endif
            u=abs(a(is,is-1)*vx(is)/vx(is-1))*(abs(w2)+abs(w3))
            v=abs(w1)*(abs(a(is-1,is-1))+abs(a(is,is))+abs(a(is1,is1)))
            if(u+v .eq. v)then
              a(is1,is-1)=0.d0
              a(is2,is-1)=0.d0
              exit
            endif
          enddo
          a(is2,is)=0.d0
          a(min(is2+1,ie),is)=0.d0
          am =hypot(w1,w2)
          am1=hypot(am,w3)
          if(is .gt. ib)then
            if(am1 .ne. 0.d0)then
              a(is,is-1)=a(is,is-1)*w1/am1*vx(is)
            endif
          endif
          do i=is-1,ie-2
            i1=i+1
            i2=i1+1
            i3=i2+1
            if(i1 .gt. is)then
              w1=vx(i1)*a(i1,i)
              w2=vx(i2)*a(i2,i)
              a(i2,i)=0.d0
              am =hypot(w1,w2)
              if(i3 .le. ie)then
                w3=vx(i3)*a(i3,i)
                a(i3,i)=0.d0
                am1=hypot(am,w3)
              else
                w3=0.d0
                am1=0.d0
              endif
            endif
            if(am .ne. 0.d0)then
              c=w1/am
              s=w2/am
              if(abs(c) .gt. abs(s))then
                v1=vx(i1)/c
                v2=c*vx(i2)
                p=s*vx(i1)/v2
                q=s*vx(i2)/v1
                vx(i1)=v1
                if(abs(v2) .gt. vmin .and. abs(v2) .lt. vmax)then
                  vx(i2)=v2
                  a(i2,i1:n)=a(i2,i1:n)-p*a(i1,i1:n)
                  a(i1,i1:n)=a(i1,i1:n)+q*a(i2,i1:n)
                  j1=min(ie,i3)
                  a(1:j1,i1)=a(1:j1,i1)+p*a(1:j1,i2)
                  a(1:j1,i2)=a(1:j1,i2)-q*a(1:j1,i1)
                  w(1:n,i1)=w(1:n,i1)+p*w(1:n,i2)
                  w(1:n,i2)=w(1:n,i2)-q*w(1:n,i1)
                else
                  vx(i2)=1.d0
                  aa(i1:n)=a(i2,i1:n)-p*a(i1,i1:n)
                  a(i1,i1:n)=a(i1,i1:n)+q*aa(i1:n)
                  a(i2,i1:n)=v2*aa(i1:n)
                  j1=min(ie,i3)
                  a(1:j1,i1)=a(1:j1,i1)+p*a(1:j1,i2)
                  a(1:j1,i2)=(a(1:j1,i2)-q*a(1:j1,i1))/v2
                  w(1:n,i1)=w(1:n,i1)+p*w(1:n,i2)
                  w(1:n,i2)=(w(1:n,i2)-q*w(1:n,i1))/v2
                endif
              else
                v1=vx(i2)/s
                v2=vx(i1)*s
                if(i1 .gt. is)then
                  a(i1,i)=w2/vx(i2)
                endif
                p=c*vx(i2)/v2
                q=c*vx(i1)/v1
                vx(i1)=v1
                if(abs(v2) .gt. vmin .and. abs(v2) .lt. vmax)then
                  vx(i2)=v2
                  aa(i1:n)=a(i2,i1:n)
                  a(i2,i1:n)=p*aa(i1:n)-a(i1,i1:n)
                  a(i1,i1:n)=aa(i1:n)-q*a(i2,i1:n)
                  j1=min(ie,i3)
                  aa(1:j1)=a(1:j1,i1)
                  a(1:j1,i1)= p*aa(1:j1)+a(1:j1,i2)
                  a(1:j1,i2)= q*a(1:j1,i1)-aa(1:j1)
                  aa=w(1:n,i1)
                  w(1:n,i1)= p*aa+w(1:n,i2)
                  w(1:n,i2)= q*w(1:n,i1)-aa
                else
                  vx(i2)=1.d0
                  aa(i1:n)=a(i2,i1:n)
                  aa1(i1:n)=p*aa(i1:n)-a(i1,i1:n)
                  a(i1,i1:n)=aa(i1:n)-q*aa1(i1:n)
                  a(i2,i1:n)=aa1(i1:n)*v2
                  j1=min(ie,i3)
                  aa(1:j1)=a(1:j1,i1)
                  a(1:j1,i1)= p*aa(1:j1)+a(1:j1,i2)
                  a(1:j1,i2)= (q*a(1:j1,i1)-aa(1:j1))/v2
                  aa=w(1:n,i1)
                  w(1:n,i1)= p*aa+w(1:n,i2)
                  w(1:n,i2)= (q*w(1:n,i1)-aa)/v2
                endif
              endif
            endif
            if(am1 .ne. 0.d0)then
              c=am/am1
              s=w3/am1
              if(abs(c) .gt. abs(s))then
                v1=vx(i1)/c
                v2=c*vx(i3)
                p=s*vx(i1)/v2
                q=s*vx(i3)/v1
                vx(i1)=v1
                vx(i3)=v2
                a(i3,i1:n)=a(i3,i1:n)-p*a(i1,i1:n)
                a(i1,i1:n)=a(i1,i1:n)+q*a(i3,i1:n)
                j1=min(ie,i3+1)
                a(1:j1,i1)=a(1:j1,i1)+p*a(1:j1,i3)
                a(1:j1,i3)=a(1:j1,i3)-q*a(1:j1,i1)
                w(1:n,i1)=w(1:n,i1)+p*w(1:n,i3)
                w(1:n,i3)=w(1:n,i3)-q*w(1:n,i1)
              else
                v1=vx(i3)/s
                v2=vx(i1)*s
                if(i1 .gt. is)then
                  a(i1,i)=w3/vx(i3)
                endif
                p=c*vx(i3)/v2
                q=c*vx(i1)/v1
                vx(i1)=v1
                vx(i3)=v2
                aa(i1:n)=a(i3,i1:n)
                a(i3,i1:n)=p*aa(i1:n)-a(i1,i1:n)
                a(i1,i1:n)=aa(i1:n)-q*a(i3,i1:n)
                j1=min(ie,i3+1)
                aa(1:j1)=a(1:j1,i1)
                a(1:j1,i1)= p*aa(1:j1)+a(1:j1,i3)
                a(1:j1,i3)= q*a(1:j1,i1)-aa(1:j1)
                aa=w(1:n,i1)
                w(1:n,i1)= p*aa+w(1:n,i3)
                w(1:n,i3)= q*w(1:n,i1)-aa
              endif
c     if(abs(vx(i1)) .gt. vmax .or. abs(vx(i1)) .lt. vmin .or.
c     $         abs(vx(i2)) .gt. vmax .or. abs(vx(i2)) .lt. vmin)then
c     write(*,*)'tqs-10 ',i1,vx(i1),vx(i2),v1,v2,c,s,
c     $           vx(8437),a(8437,8435)
c     endif
            endif
c     write(*,'(1P4G15.7)')((vx(ii)/vx(j)*a(ii,j),j=1,4),ii=1,4)
          enddo
          if(is .gt. ib)then
            a(is,is-1)=a(is,is-1)/vx(is)
          endif
          iter=iter+1
          if(iter .gt. itm)then
            write(*,*)' TEIGEN convergence failed. Range =',ib,ie
            write(*,*)
     1           '        Lower right corner =',
     1           a(ie-1,ie-1),vx(ie)/vx(ie-1)*a(ie,ie-1),a(ie,ie)
            a(ie,ie-1)=0.d0
            if(ktfenanq(a(ie,ie)))then
              a(ie,ie)=0.d0
              a(ie-1,ie-1)=0.d0
            endif
            if(ktfenanq(vx(ie)))then
              vx(ie)=1.d0
            endif
            if(ktfenanq(vx(ie-1)))then
              vx(ie-1)=1.d0
            endif
            iter=0
            ie=ie-1
            cycle do1
          endif
        enddo do2
      enddo do1
      do i=1,n
        j1=max(i-1,1)
        a(i,j1:n)=a(i,j1:n)*vx(i)
        j1=min(n,i+1)
        a(1:j1,i)=a(1:j1,i)/vx(i)
        w(1:n,i)=w(1:n,i)/vx(i)
      enddo
      paired=.false.
      do i=1,n-1
        if(paired)then
          paired=.false.
          cycle
        endif
        i1=i+1
        if(a(i1,i) .ne. 0.d0)then
          if(abs(a(i1,i)) .gt. abs(a(i,i1)))then
            aa(i:n)=a(i,i:n)
            a(i,i:n)=a(i1,i:n)
            a(i1,i:n)=aa(i:n)
            aa(1:i1)=a(1:i1,i)
            a(1:i1,i)=a(1:i1,i1)
            a(1:i1,i1)=aa(1:i1)
            aa=w(1:n,i)
            w(1:n,i)=w(1:n,i1)
            w(1:n,i1)=aa
          endif
          if(a(i1,i) .eq. 0.d0)then
            eig(1,i)=a(i,i)
            eig(2,i)=0.d0
            cycle
          endif
          s=a(i,i)-a(i1,i1)
          dec=s**2+4.d0*a(i1,i)*a(i,i1)
          s1=a(i,i)**2+a(i1,i1)**2
          if(dec .ge. -decth*s1)then
            sqrd=sqrt(abs(dec))
            p=2.d0*a(i1,i)/(s+sign(sqrd,s))
            a(i1,i:n)=a(i1,i:n)-p*a(i,i:n)
            a(1:i1,i)=a(1:i1,i)+p*a(1:i1,i1)
            w(1:n,i)=w(1:n,i)+p*w(1:n,i1)
            eig(1,i)=a(i,i)
            eig(1,i1)=a(i1,i1)
            eig(2,i)=0.d0
            eig(2,i1)=0.d0
            a(i1,i)=0.d0
          else
            p=-s/a(i,i1)*.5d0
            r=sign(
     1           sqrt(abs((a(i1,i)-p*(a(i,i)-a(i1,i1)+p*a(i,i1)))
     $           /a(i,i1))),a(i,i1))
            a(i1,i:n)=a(i1,i:n)-p*a(i,i:n)
            a(i ,i:n)=a(i ,i:n)*r
            a(1:i1,i)=(a(1:i1,i)+p*a(1:i1,i1))/r
            w(1:n,i)=(w(1:n,i)+p*w(1:n,i1))/r
            eig(1,i )=a(i ,i )
            eig(1,i1)=a(i1,i1)
            eig(2,i )=a(i ,i1)
            eig(2,i1)=a(i1,i )
          endif
          paired=.true.
        else
          eig(1,i)=a(i,i)
          eig(2,i)=0.d0
        endif
      enddo
      if(.not. paired)then
        eig(1,n)=a(n,n)
        eig(2,n)=0.d0
      endif
      do i=2,n
        if(a(i,i-1) .ne. 0.d0)then
          cycle
        endif
        i1=i+1
        if(eig(2,i) .eq. 0.d0)then
          jm=0
          do  j=i-1,1,-1
            if(a(j+1,j) .ne. 0.d0)then
              cycle
            endif
            if(eig(2,j) .eq. 0.d0)then
              da=a(j,j)-a(i,i)
              sa=abs(a(j,j))+abs(a(i,i))
              if(abs(da) .gt. threj*sa)then
                p=a(j,i)/da
                a(j,i:n)=a(j,i:n)+p*a(i,i:n)
                a(1:j,i)=a(1:j,i)-p*a(1:j,j)
                w(1:n,i)=w(1:n,i)-p*w(1:n,j)
              else
                jordan=.true.
              endif
            else
              j1=j-1
              u=a(i,i)-a(j1,j1)
              v=a(j1,j)
              r=u**2+v**2
              pk=(a(j1,i)*u+a(j ,i)*v)/r
              pl=(a(j ,i)*u-a(j1,i)*v)/r
              a(j1,i:n)=a(j1,i:n)-pk*a(i,i:n)
              a(j ,i:n)=a(j ,i:n)-pl*a(i,i:n)
              a(1:j,i)=a(1:j,i)+pk*a(1:j,j1)+pl*a(1:j,j)
              w(1:n,i)=w(1:n,i)+pk*w(1:n,j1)+pl*w(1:n,j)
            endif
          enddo
        else
          do j=i-1,1,-1
            if(a(j+1,j) .ne. 0.d0)then
              cycle
            endif
            if(eig(2,j) .eq. 0.d0)then
              x=a(i,i)-a(j,j)
              y=a(i,i1)
              r=x**2+y**2
              pk=(a(j,i )*x+a(j,i1)*y)/r
              pl=(a(j,i1)*x-a(j,i )*y)/r
              a(j,i:n)=a(j,i:n)-pk*a(i,i:n)-pl*a(i1,i:n)
              a(1:j,i )=a(1:j,i )+pk*a(1:j,j)
              a(1:j,i1)=a(1:j,i1)+pl*a(1:j,j)
              w(1:n,i )=w(1:n,i )+pk*w(1:n,j)
              w(1:n,i1)=w(1:n,i1)+pl*w(1:n,j)
            else
              j1=j-1
              cu=dcmplx(a(j1,j1)-a(i,i),-a(i,i1))
              v =a(j1,j)
              cr=cu**2+v**2
              if(cr .ne. (0.d0,0.d0))then
                ca=dcmplx(a(j1,i),a(j1,i1))
                cc=dcmplx(a(j ,i),a(j ,i1))
                ck=(cu*ca-v*cc)/cr
                cm=(cu*cc+v*ca)/cr
                a(j1,i:n)=a(j1,i:n)+pk*a(i,i:n)+pl*a(i1,i:n)
                a(j,i:n )=a(j ,i:n)+pm*a(i,i:n)+pn*a(i1,i:n)
                a(1:j,i )=a(1:j,i )-pk*a(1:j,j1)-pm*a(1:j,j)
                a(1:j,i1)=a(1:j,i1)-pl*a(1:j,j1)-pn*a(1:j,j)
                w(1:n,i )=w(1:n,i )-pk*w(1:n,j1)-pm*w(1:n,j)
                w(1:n,i1)=w(1:n,i1)-pl*w(1:n,j1)-pn*w(1:n,j)
              endif
            endif
          enddo
        endif
      enddo
      if(jordan)then
        do i=n,2,-1
          if(eig(2,i) .eq. 0.d0)then
            jm=0
            do j=i-1,1,-1
              if(eig(1,j) .eq. eig(1,i) .and. eig(2,j) .eq. 0.d0)then
                r=a(j,i)
                if(abs(r) .gt. amth)then
                  if(jm .eq. 0)then
                    r1=sqrt(abs(r))
                    r2=r/r1
                    a(j,j:n)=a(j,j:n)/r1
                    a(1:j,j)=a(1:j,j)*r1
                    w(1:n,j)=w(1:n,j)*r1
                    a(1:i,i)=a(1:i,i)/r2
                    a(i,i:n)=a(i,i:n)*r2
                    w(1:n,i)=w(1:n,i)/r2
                    jm=j
                  else
                    a(j,jm:n)=a(j,jm:n)-r*a(jm,jm:n)
                    a(1:j,jm)=a(1:j,jm)+r*a(1:j,j)
                    w(1:n,jm)=w(1:n,jm)+r*w(1:n,j)
                  endif
                else
                  a(j,i)=0.d0
                endif
              endif
            enddo
          endif
        enddo
      endif
c     do 1000 k=1,n
c       write(*,'(1X,:1P10G13.6)')(a(k,j),j=1,n)
c1000  continue
c     do 1010 k=1,n
c       write(*,'(1X,:1P10G13.6)')(w(k,j),j=1,n)
c1010  continue
c     write(*,*)
      if(mod(n,2) .eq. 0)then
        ii=0
        do i=1,n
          if(eig(2,i) .gt. 0.d0 .and. ii .ne. 0)then
            aa=w(1:n,i-1)
            w(1:n,i-1)=w(1:n,i  )
            w(1:n,i  )=w(1:n,i+1)
            w(1:n,i+1)=aa
            ee=eig(1,i-1)
            eig(1,i-1)=eig(1,i  )
            eig(1,i  )=eig(1,i+1)
            eig(1,i+1)=ee
            ee=eig(2,i-1)
            eig(2,i-1)=eig(2,i  )
            eig(2,i  )=eig(2,i+1)
            eig(2,i+1)=ee
          endif
          ii=1-ii
        enddo
      endif
      return
      end subroutine

      subroutine teigenc(a,w,ceig,n,ndim)
      use iso_c_binding
      implicit none
      integer*4 ,intent(in):: n,ndim
      real*8 ,intent(inout):: a(ndim,n),w(n,n)
      real*8 ,pointer :: eig(:,:)
      complex*16 ,target ::ceig(n)
      call c_f_pointer(c_loc(ceig),eig,[2,n])
      call teigen(a,w,eig,n,ndim)
      return
      end subroutine

      end module

      subroutine teigen(a,w,eig,n,ndim)
c
c   Subroutine to obtain eigen values and eigen vectors of a real
c   matrix a.
c                                7-Dec-1987   K. Oide
c   Modified                    27-Oct-1989   K. Oide
c   Modified                    17-Nov-1989   K. Oide
c   Added Balancing              1-Nov-1990   K. Oide
c   Modified Convergence        20-Nov-1990   K. Oide
c   Modified Shifting            4-Dec-1990   K. Oide
c   Fast Givens' Transformation 21-Jan-1992   K. Oide
c   Maximum size 14,000         29-May-2010   K. Oide
c   Maximum size 32768           5-Oct-2010   K. Oide
c
c   Usage:       call teigen(a,w,eig,n,ndim)
c
c                where
c
c                a     is the input n*n matrix, and the eigen vectors
c                      are returned here.
c                w     is a work area of n**2 words.
c                eig   is a (2,n) real array where the eigen values
c                      are returned.   eig(1,i) and eig(2,i) contain
c                      the real and imaginary part of the i-th eigen-
c                      valeu respectively.   When the eigenvalue is
c                      a complex, the output are stored as
c
c                      eig(1,i  ) =      er
c                      eig(2,i  ) =     !ei!
c                      eig(1,i+1) =      er
c                      eig(2,i+1) =    -!ei!
c
c                      where er and ei are the real and imaginary
c                      part of the i-th eigenvalue respectively.
c                n     is the matrix dimension.
c                ndim  is the size of the first dimension of a.
c
c   Restriction: When the matrix cannot be diagonalized (Jordan
c                standard form), the returned eigen vectors are not
c                correct.
c
      use eigen
      implicit none
      integer*4 n,ndim,nmax
      parameter (nmax=32768)
      real*8 a(ndim,n),w(n,n),eig(2,n)
      if(n .gt. nmax)then
        write(*,*)'TEIGEN-Too large matrix.',n
      else
        allocate (vx(n))
        call tbal(a,w,n,ndim)
        call thess(w,a,n,ndim)
        call tqr(w,a,eig,n,ndim)
        deallocate (vx)
      endif
      return
      end subroutine
