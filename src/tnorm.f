      subroutine tnorm(r,ceig,lfno)
      implicit none
      integer*4 i,j
      integer*4 ,intent(in):: lfno
      real*8 ,intent(inout):: r(6,6)
      real*8 sa,sb,a,s,cost,sint,s1,s2,s3,b,v(6,2)
      complex*16 ,intent(out):: ceig(6)
      complex*16 cc
      do i=1,5,2
        if(imag(ceig(i)) .eq. 0.d0)then
          sa=sum(r(:,i)**2)
          sb=sum(r(:,i+1)**2)
          a=sqrt(sqrt(sa/sb))
          r(:,i  )=r(:,i  )/a
          r(:,i+1)=r(:,i+1)*a
        endif
        s=0.d0
        do j=1,5,2
          s=s+r(j,i)*r(j+1,i+1)-r(j,i+1)*r(j+1,i)
        enddo
        sa=sqrt(abs(s))
        if(s .gt. 0.d0)then
          sb=sa
        elseif(s .eq. 0.d0)then
          if(lfno .ne. 0)then
            write(lfno,*)'???-tnorm-Unstable transfer matrix.'
          endif
          sa=1.d0
          sb=1.d0
        else
          sb=-sa
          ceig(i  )=conjg(ceig(i  ))
          ceig(i+1)=conjg(ceig(i+1))
        endif
        r(:,i  )=r(:,i  )/sa
        r(:,i+1)=r(:,i+1)/sb
      enddo
      s1=(r(5,1)*r(6,2)-r(5,2)*r(6,1))**2
      s2=(r(5,3)*r(6,4)-r(5,4)*r(6,3))**2
      s3=(r(5,5)*r(6,6)-r(5,6)*r(6,5))**2
      j=maxloc((/s1,s2,s3/),1)*2-1
      if(j .ne. 5)then
        v=r(:,5:6)
        r(:,5:6)=r(:,j:j+1)
        r(:,j:j+1)=v
        cc=ceig(5)
        ceig(5)=ceig(j)
        ceig(j)=cc
        cc=ceig(6)
        ceig(6)=ceig(j+1)
        ceig(j+1)=cc
      endif
      s1=(r(1,1)*r(2,2)-r(1,2)*r(2,1))**2
      s2=(r(1,3)*r(2,4)-r(1,4)*r(2,3))**2
      if(s2 .gt. s1)then
        v=r(:,1:2)
        r(:,1:2)=r(:,3:4)
        r(:,3:4)=v
        cc=ceig(1)
        ceig(1)=ceig(3)
        ceig(3)=cc
        cc=ceig(2)
        ceig(2)=ceig(4)
        ceig(4)=cc
      endif
      do i=1,5,2
        if(imag(ceig(i)) .ne. 0.d0)then
          a=hypot(r(i,i),r(i,i+1))
          cost=r(i,i  )/a
          sint=r(i,i+1)/a
          v(:,1)=r(:,i)
          r(:,i  )= v(:,1)*cost+r(:,i+1)*sint
          r(:,i+1)=-v(:,1)*sint+r(:,i+1)*cost
        else
          a=r(i,i)
          if(a .ne. 0.d0)then
            b=r(i,i+1)
            r(:,i+1)=r(:,i+1)*a-b*r(:,i)
            r(:,i  )=r(:,i  )/a
          endif
        endif
      enddo
c     call tsymp(r)
      return
      end
