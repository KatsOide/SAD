      subroutine tceigen(ca,cw,ceig,n,ndim)
      implicit none
      integer*4 ,intent(in):: n,ndim
      integer*4 ,parameter ::itmax=30
      integer*4 i,j,k,i1,i2,it
      complex*16 ,intent(inout):: ca(ndim,n)
      complex*16 ,intent(out):: cw(n,n),ceig(n)
      complex*16 cs,cc,ccs,ccc,
     $     cpk,cr,cp,cu,cd,ctr,cdet,csqr,c1,c2
      real*8 a,w,xc,xs,xp
      cw=(0.d0,0.d0)
      do i=1,n
        cw(i,i)=(1.d0,0.d0)
      enddo
      do i=1,n-1
        cp=ca(i+1,i)
        if(imag(cp) .ne. 0.d0)then
          a=abs(cp)
          cu=conjg(cp)/a
          do j=i+1,n
            ca(i+1,j)=ca(i+1,j)*cu
          enddo
          cu=conjg(cu)
          ca(1:n,i+1)=ca(1:n,i+1)*cu
          cw(1:n,i+1)=cw(1:n,i+1)*cu
        else
          a=dble(cp)
        endif
        xp=a
        do j=i+2,n
          a=abs(dcmplx(xp,abs(ca(j,i))))
          if(a .ne. 0.d0)then
            xc=xp/a
            cs=ca(j,i)/a
            ccs=conjg(cs)
            do k=i+1,n
              cpk=ca(i+1,k)
              ca(i+1,k)= cpk*xc +ca(j,k)*ccs
              ca(j  ,k)=-cpk*cs +ca(j,k)*xc
            enddo
            do k=1,n
              cpk=ca(k,i+1)
              ca(k,i+1)= cpk*xc +ca(k,j)*cs
              ca(k,j  )=-cpk*ccs+ca(k,j)*xc
              cpk=cw(k,i+1)
              cw(k,i+1)= cpk*xc +cw(k,j)*cs
              cw(k,j  )=-cpk*ccs+cw(k,j)*xc
            enddo
            xp=a
          endif
        enddo
        ca(i+1,i)=xp
      enddo
      i1=1
      i2=n
 1    if(i2 .eq. 1)then
        go to 3000
      endif
      if(i1 .eq. i2)then
        i2=i1-1
        i1=1
        go to 1
      endif
      it=0
 10   do i=i2-1,i1,-1
        w=abs(ca(i,i))+abs(ca(i+1,i+1))
        if(abs(dble(ca(i+1,i)))+w .eq. w)then
          ca(i+1,i)=0.d0
          if(i .eq. i2-1)then
            i2=i2-1
          else
            i1=i+1
          endif
          go to 1
        endif
      enddo
      ctr=ca(i2-1,i2-1)+ca(i2,i2)
      cdet=2.d0*(ca(i2-1,i2-1)*ca(i2,i2)-
     $     dble(ca(i2,i2-1))*ca(i2-1,i2))
      csqr=sqrt(ctr**2-2.d0*cdet)
      c1=ctr+csqr
      c2=ctr-csqr
      if(abs(c1) .gt. abs(c2))then
        cr=cdet/c1
      else
        cr=cdet/c2
      endif
      cp=ca(i1,i1)-cr
      a=abs(dcmplx(abs(cp),dble(ca(i1+1,i1))))
      cc=cp/a
      xs=dble(ca(i1+1,i1))/a
      ccc=conjg(cc)
      do k=i1,n
        cpk=ca(i1,k)
        ca(i1  ,k)= cpk*ccc+ca(i1+1,k)*xs
        ca(i1+1,k)=-cpk*xs +ca(i1+1,k)*cc
      enddo
      do k=1,i1+1
        cpk=ca(k,i1)
        ca(k,i1  )= cpk*cc +ca(k,i1+1)*xs
        ca(k,i1+1)=-cpk*xs +ca(k,i1+1)*ccc
      enddo
      if(i1+2 .le. i2)then
        ca(i1+2,i1  )= dble(ca(i1+2,i1+1))*xs
        ca(i1+2,i1+1)= dble(ca(i1+2,i1+1))*ccc
      endif
      do k=1,n
        cpk=cw(k,i1)
        cw(k,i1  )= cpk*cc +cw(k,i1+1)*xs
        cw(k,i1+1)=-cpk*xs +cw(k,i1+1)*ccc
      enddo
      do i=i1,i2-2
        cp=ca(i+1,i)
        a=abs(dcmplx(abs(cp),dble(ca(i+2,i))))
        if(a .ne. 0.d0)then
          cc=cp/a
          xs=dble(ca(i+2,i))/a
          ccc=conjg(cc)
          do k=i+1,n
            cpk=ca(i+1,k)
            ca(i+1,k)= cpk*ccc+ca(i+2,k)*xs
            ca(i+2,k)=-cpk*xs +ca(i+2,k)*cc
          enddo
          ca(i+1,i)=a
          do k=1,i+2
            cpk=ca(k,i+1)
            ca(k,i+1)= cpk*cc +ca(k,i+2)*xs
            ca(k,i+2)=-cpk*xs +ca(k,i+2)*ccc
          enddo
          if(i+3 .le. i2)then
            ca(i+3,i+1)= dble(ca(i+3,i+2))*xs
            ca(i+3,i+2)= dble(ca(i+3,i+2))*ccc
          endif
          do k=1,n
            cpk=cw(k,i+1)
            cw(k,i+1)= cpk*cc +cw(k,i+2)*xs
            cw(k,i+2)=-cpk*xs +cw(k,i+2)*ccc
          enddo
        endif
      enddo
      if(imag(ca(i2,i2-1)) .ne. 0.d0)then
        a=abs(ca(i2,i2-1))
        cu=conjg(ca(i2,i2-1))/a
        ca(i2,i2-1)=a
c        do i=i2+1,n
        ca(i2,i2+1:n)=ca(i2,i2+1:n)*cu
c        enddo
        cu=conjg(cu)
c        do i=1,i2-1
        ca(1:i2-1,i2)=ca(1:i2-1,i2)*cu
c        enddo
c        do i=1,n
        cw(:,i2)=cw(:,i2)*cu
c        enddo
      endif
      it=it+1
      if(it .gt. itmax)then
        write(*,*)'TCEIGEN Convergence fail.'
        ca(i2,i2-1)=0.d0
      endif
      go to 10
 3000 do i=2,n
        do j=i-1,1,-1
          cd=ca(j,j)-ca(i,i)
          if(cd .ne. 0.d0)then
            cu=ca(j,i)/cd
c            do k=i,n
            ca(j,i:n)=ca(j,i:n)+cu*ca(i,i:n)
c            enddo
c            do k=1,j-1
            ca(1:j-1,i)=ca(1:j-1,i)-cu*ca(1:j-1,j)
c            enddo
            ca(j,i)=0.d0
c            do k=1,n
            cw(:,i)=cw(:,i)-cu*cw(:,j)
c            enddo
          endif
        enddo
      enddo
      do i=1,n
        ceig(i)=ca(i,i)
      enddo
      ca(1:n,:)=cw
      return
      end
