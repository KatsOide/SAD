c     Don't confuse, Emacs. This is -*- fortran -*- mode!
C TENT
      if(theta .ne. 0.d0)then
        cost=cos(theta)
        sint=sin(theta)
      else
        cost=1.d0
        sint=0.d0
      endif
      if(dx .ne. 0.d0)then
        x(1:np)=x(1:np)-dx
      endif
      if(dy .ne. 0.d0)then
        y(1:np)=y(1:np)-dy
      endif
      if(theta .ne. 0.d0)then
        if(abs(cost) .ge. abs(sint))then
          x(1:np)=cost*x(1:np)-sint*y(1:np)
          y(1:np)=(sint*x(1:np)+y(1:np))/cost
          px(1:np)=cost*px(1:np)-sint*py(1:np)
          py(1:np)=(sint*px(1:np)+py(1:np))/cost
        else
          do i=1,np
            xi=x(i)
            x(i)=xi*cost-y(i)*sint
            y(i)=xi*sint+y(i)*cost
            pxi=px(i)
            px(i)=pxi*cost-py(i)*sint
            py(i)=pxi*sint+py(i)*cost
          enddo
        endif
        if(calpol)then
          do i=1,np
            xi=sx(i)
            sx(i)=cost*xi-sint*sy(i)
            sy(i)=sint*xi+cost*sy(i)
          enddo
        endif
      endif
