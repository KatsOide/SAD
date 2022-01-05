      subroutine temp(pexln)
      use macmath
      use cbkmac
      implicit real*8 (a-h,o-z)
      integer pexln,status
      real*8 tm(6,6),codin(6),codout(6)
      real*8 tmw(6,6),codw(6)
      real*8 beta(2),alpha(2),det(2),c(2),s(2),mu(2)
      real*8 beta2(2),alpha2(2),dmu(2)
c for debug
      call ptrace('temp',+1)
c end debug
      status = 0
      do 1100 i=1,6
        codin(i)=0.0d0
        codout(i)=0.0d0
 1100 continue
c
      nel=ilist(1,pexln)
      call m_unit(tm,6)
      do 2100 i=1,nel
        call optelm(pexln+i
     &        ,codin,tm,codw,tmw,status)
        call m_copy_6d(tmw,tm)
        call v_copy(codw,codin,6)
 2100 continue
c
      do 3100 i=1,2
        i1=2*i-1
        i2=2*i
        c(i)=(tm(i1,i1)+tm(i2,i2))/2.0d0
        if (abs(c(i)).ge. 1.d0) then
          call errmsg('tmp',
     &                'unstable optics',0,0)
          status=-1
        else
          if (tm(i1,i2) .gt. 0.0d0) then
            s(i)=sqrt(1.d0-c(i)**2)
          else
            s(i)=-sqrt(1.d0-c(i)**2)
          endif
          beta(i)=tm(i1,i2)/s(i)
          alpha(i)=(tm(i1,i1)-tm(i2,i2))/(2.d0*s(i))
          det(i)=tm(i1,i1)*tm(i2,i2)-tm(i1,i2)*tm(i2,i1)
          mu(i)=atan2(s(i),c(i))/(2.d0*pi)
        endif
 3100 continue
      if (status .eq. -1) then
        write(outfl,
     &        '((1H ,8X,6G10.4))'
     &       ) ((tm(i,j),j=1,6),i=1,6)
      else
        call prntws('        ',alpha,beta,mu,codin,det)
        call v_clear(dmu,2)
        do 3200 i=1,nel
           call m_unit(tm,6)
           call optelm(pexln+i
     &                ,codin,tm,codw,tmw,status)
           call m_copy_6d(tmw,tm)
           call v_copy(codw,codin,6)
          do 3210 j=1,2
            j1=2*j-1
            j2=2*j
            g0=(1.d0+alpha(j)**2)/beta(j)
            alpha2(j)=(tm(j1,j1)*tm(j2,j2)+tm(j1,j2)*tm(j2,j1))*alpha(j)
     &              -tm(j1,j1)*tm(j2,j1)*beta(j)-tm(j1,j2)*tm(j2,j2)*g0
            beta2(j)=tm(j1,j1)**2*beta(j)
     &              -2.d0*tm(j1,j1)*tm(j1,j2)*alpha(j)
     &              +tm(j1,j2)**2*g0
            det(j)=tm(j1,j1)*tm(j2,j2)-tm(j1,j2)*tm(j2,j1)
            s(j)=tm(j1,j2)/sqrt(beta(j)*beta2(j))
            c(j)=sqrt(beta(j)/beta2(j))*tm(j1,j1)
     &                     -alpha(j)*s(j)
            dmu(j)=dmu(j)+0.125d0*atan2(s(j),c(j))/atan(1.0d0)
 3210     continue
          call prntws(pname(ilist(1,pexln+i)),
     &                alpha2,beta2,dmu,codin,det)
          do 3220 j=1,2
            beta(j)=beta2(j)
            alpha(j)=alpha2(j)
 3220     continue
c         call v_copy(beta2,beta,2)
c         call v_copy(alpha2,alpha,2)
 3200   continue
      endif
c for debug
      call ptrace('temp',-1)
c end debug
      return
      end
