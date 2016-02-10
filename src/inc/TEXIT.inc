c     Don't confuse, Emacs. This is -*- fortran -*- mode!
        if(theta .ne. 0.d0)then
          do i=1,np
            xi=x(i)
            x(i)= xi*cost+y(i)*sint
            y(i)=-xi*sint+y(i)*cost
            pxi=px(i)
            px(i)= pxi*cost+py(i)*sint
            py(i)=-pxi*sint+py(i)*cost
          enddo
        endif
        if(dx .ne. 0.d0)then
          do i=1,np
            x(i)=x(i)+dx
          enddo
        endif
        if(dy .ne. 0.d0)then
          do i=1,np
            y(i)=y(i)+dy
          enddo
        endif
