ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c  Phase space rotation including 6 dimensional twiss parameters  c
c                                                                 c
c                           written by K.Ohmi                     c
c                                                                 c
c                                        July 1994                c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine phsinit(p_in,blist)
      use ffs_flag
      use tmacro
      implicit real*8(a-h,o-z)
      real*8 blist(*),p_in(*)
      real*8 tmp(36),eig(12)
c      real*8 tmp2(36)
      integer m_PtoN,m_prot,m_NtoP,v_eig,v_difu
      m_PtoN=1
      m_prot=m_PtoN+36
      m_NtoP=m_prot+36
      m_damp=m_NtoP+36
      m_difu=m_damp+36
      m_prot_d=m_difu+36
      v_eig=m_prot_d+36
      v_difu=v_eig+6
c
      call m_PtoN_set_6d(p_in,blist(m_PtoN))
      write(*,*) 'Transfer matrix from physical to normal coordinate'
      call m_print(blist(m_PtoN),6)
      call m_prot_set_6d(p_in,blist(m_prot))
      call m_NtoP_set_6d(p_in,blist(m_NtoP))
      write(*,'(/A)') 'Matrices on the normal mode'
      write(*,*) 'Symplectic revolution matrix'
      call m_print(blist(m_prot),6)
c
      if(rad) then
      call m_transpose_6d(p_in(30),blist(m_damp))
      call m_symmet_set_6d(p_in(66),blist(m_difu))
      call m_mul_6d(blist(m_prot),blist(m_damp),tmp)
c  M=M_0 (1-D)
      call m_sub_6d(blist(m_prot),tmp,blist(m_prot_d))
      write(*,*) 'Damping matrix'
      call m_print(blist(m_damp),6)
      write(*,*) 'Diffusion matrix'
      call m_print(blist(m_difu),6)
      else
      call m_copy_6d(blist(m_prot),blist(m_prot_d))
      endif
      call m_mul3_6d(blist(m_NtoP),blist(m_prot_d),blist(m_PtoN),
     &               blist(m_prot_d))
      write(*,*) 'Revolution matrix (including radiation damping)'
      call m_print(blist(m_prot_d),6)
c
      if(rfluct) then
      call teigen(blist(m_difu),tmp,eig,6,6)
c m_difu becomes diagonal tranfformation to get the eigen values.
      write(*,'(/A)') 'U which diagonalize B'
c      call m_print(blist(m_difu),6)
      call m_normalize_6d(blist(m_difu),eig)
      call m_print(blist(m_difu),6)
      write(*,*) 'Eigen values of B'
      write(*,'(6E12.6)') (eig(2*i-1),i=1,6)
      write(*,'(6E12.6)') (eig(2*i),i=1,6)
c      write(*,*) 'check for the diagonalization'
c      call m_transpose_6d(blist(m_difu),tmp)
c      call m_symmet_set_6d(p_in(66),tmp2)
c      call m_mul3_6d(tmp,tmp2,blist(m_difu),tmp)
c      call m_print(tmp,6)
c         check OK 94.9.13
      do 15 i=1,6
         if(eig(2*i-1).lt.0.) then
            write(*,*) 'phsinit : diffusion matrix has zero or', 
     &         ' negative eigen values'
            stop
         endif
         blist(v_eig+i-1)=sqrt(eig(2*i-1))
 15      continue
      call mv_mul_6d(blist(m_difu),blist(v_eig),blist(v_difu))
      write(*,'(/A)') 'Diffusion vector'
      write(*,'(6E12.6/)') (blist(v_difu+i),i=0,5)
      endif
c      stop
      return
      end
cc
c ####### Tracking routine ##################################
c   
c  Modified args include spins without transformation, 12 Jul 2020. 
      subroutine phsrot(np,x,px,y,py,z,g,dv,sx,sy,sz,blist)
      use ffs_flag
      use tmacro
      use wsbb, only:nblist
      implicit real*8(a-h,o-z)
      integer*4 ,intent(in):: np
      integer v_eig,v_difu
      parameter(m_PtoN=1,m_prot=37,m_NtoP=73,
     &          m_damp=109,m_difu=145,m_prot_d=181,
     &          v_eig=217,v_difu=223)
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np),
     $     dv(np), sx(np),sy(np),sz(np)
      real*8 blist(*)
      real*8 xdifu(6),xo(6)
c
cc
      do i=1,np
c
c   Get Canonical variables
c       g(i)=pn-1.d0   Not good for precision
c       g(i)=g(i)*(g(i)+2.d0)
c   pn=|p|/p0=1+delta
       pn=1.d0+g(i)
       px(i)=px(i)*pn
       py(i)=py(i)*pn
      enddo
c
      do i=1,np
c
       do j=1,6
         xo(j)=blist(j-1+m_prot_d)*x(i)
     &        +blist(j+5+m_prot_d)*px(i)
     &        +blist(j+11+m_prot_d)*y(i)
     &        +blist(j+17+m_prot_d)*py(i)
     &        +blist(j+23+m_prot_d)*z(i)
     &        +blist(j+29+m_prot_d)*g(i)
       enddo
c
         x(i)=xo(1)
         px(i)=xo(2)
         y(i)=xo(3)
         py(i)=xo(4)
         z(i)=xo(5)
         g(i)=xo(6)
c
c  quantum diffusion
       if(rfluct) then
       call tgauss_array(xo(1),6)
       do j=1,6
         xo(j)=blist(v_difu+j-1)*xo(j)
       enddo
       do j=1,6
         xdifu(j)=0.
         do k=1,6
           xdifu(j)=xdifu(j)+blist(j-1+6*(k-1)+m_NtoP)*xo(k)
         enddo
       enddo
       x(i)=x(i)+xdifu(1)
       px(i)=px(i)+xdifu(2)
       y(i)=y(i)+xdifu(3)
       py(i)=py(i)+xdifu(4)
       z(i)=z(i)+xdifu(5)
       g(i)=g(i)+xdifu(6)
       endif
c
      enddo
c
c   Return to SAD variables
c
c   pn=|p|/p0=1+delta
      do i=1,np
      pn=1.d0+g(i)
      px(i)=px(i)/pn
      py(i)=py(i)/pn
c   When you change momentum, change dv(i) as
      h1=sqrt(1.d0+(pn*blist(i_p0+1))**2)
      dv(i)=-g(i)*(1.d0+pn)/h1/(h1+pn*blist(i_p0+2))+dvfs
c   g(i)=sqrt(pn)-1.d0  Not good for precision
c      g(i)=g(i)/(1.d0+sqrt(pn))
      enddo
c
c     number=number+1
c     write(*,*)'number=',number
      return
      end
