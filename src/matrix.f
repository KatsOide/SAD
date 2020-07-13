c
c     Matrix operation library
c
      subroutine v_clear(a,n)
      implicit real*8(a-h,o-z)
      real*8 a(n)
      a=0.d0
c      do 10 i=1,n
c 10      a(i)=0.d0
      return
      end
c ------------------------------------------------------------------
      subroutine m_clear(a,n)
      implicit real*8(a-h,o-z)
      real*8 a(*)
      a(1:n*n)=0.d0
c      do 10 i=1,n*n
c 10      a(i)=0.d0
      return
      end
c ------------------------------------------------------------------
      subroutine m_unit(a,n)
      implicit real*8(a-h,o-z)
      real*8 a(n,n)
      a=0.d0
c      do 20 i=1,n
c         do 20 j=1,n
c            a(j,i)=0.
c 20   continue
      do i=1,n
        a(i,i)=1.d0
      enddo
      return
      end
c ------------------------------------------------------------------
      subroutine v_copy(a,b,n)
      implicit real*8(a-h,o-z)
      real*8 a(n),b(n)
      b=a
c      do 10 i=1,n
c 10      b(i)=a(i)
      return
      end
c ------------------------------------------------------------------
      subroutine m_copy_6d(a,b)
      real*8 a(36),b(36)
      b=a
c      do 10 i=1,36
c 10      b(i)=a(i)
      return
      end
c ------------------------------------------------------------------
      subroutine v_chg_sign(a,n)
      real*8 a(n)
      a=-a
c      do 10 i=1,n
c 10      a(i)=-a(i)
      return
      end
c ------------------------------------------------------------------
      subroutine m_print(a,n)
      integer*4 ,intent(in):: n
      real*8 ,intent(in):: a(n,n)
      integer*4 i
c      write(*,*) 'Matrix print out : Ndim=',n
      do i=1,n
        write(*,'(1P,(6D12.5))') (a(i,j),j=1,n)
      enddo
      return
      end
c ------------------------------------------------------------------
      subroutine m_add_6d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(36),b(36),c(36)
      c=a+b
c      do 10 i=1,36
c      c(i)=a(i)+b(i)
c 10   continue
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_sub_6d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(36),b(36),c(36)
      c=a-b
c      do 10 i=1,36
c      c(i)=a(i)-b(i)
c 10   continue
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_mul_6d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(6,6),b(6,6),c(6,6)
c      real*8 t(6,6)
c
      c=matmul(a,b)
c      call m_clear(t,6)
c      do 10 i=1,6
c      do 20 j=1,6
c      t(i,j)=0.d0
c      do 30 k=1,6
c         t(i,j)=t(i,j)+a(i,k)*b(k,j)
c 30   continue
c 20   continue
c 10   continue
c
c     do 40 i=1,6
c      do 40 j=1,6
c      c(i,j)=t(i,j)
c 40   continue
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_mul3_6d(a,b,c,d)
      implicit real*8(a-h,o-z)
      real*8 a(*),b(*),c(*),d(*)
      real*8 t(6,6)
c
      call m_mul_6d(a,b,t)
      call m_mul_6d(t,c,d)
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_transpose_6d(a,b)
      implicit real*8(a-h,o-z)
      real*8 a(6,6),b(6,6)
c      real*8 t(6,6)
      b=transpose(a)
c
c      do 10 i=1,6
c      do 10 j=1,6
c 10      t(i,j)=a(j,i)
c      do 20 i=1,6
c      do 20 j=1,6
c 20      b(j,i)=t(j,i)
c      return
      end
c ------------------------------------------------------------------
      subroutine m_sim_tr_6d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(*),b(*),c(*)
      real*8 t(6,6)
c
      call m_transpose_6d(a,t)
      call m_mul3_6d(a,b,t,c)
      return
      end
c ------------------------------------------------------------------
      subroutine m_mul_5d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(5,5),b(5,5),c(5,5)
c      real*8 t(5,5)
c
      c=matmul(a,b)
c$$$      call m_clear(t,5)
c$$$      do 10 i=1,5
c$$$      do 20 j=1,5
c$$$      t(i,j)=0.d0
c$$$      do 30 k=1,5
c$$$         t(i,j)=t(i,j)+a(i,k)*b(k,j)
c$$$ 30   continue
c$$$ 20   continue
c$$$ 10   continue
c$$$c
c$$$      do 40 i=1,5
c$$$      do 40 j=1,5
c$$$      c(i,j)=t(i,j)
c$$$ 40   continue
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_mul3_5d(a,b,c,d)
      implicit real*8(a-h,o-z)
      real*8 a(*),b(*),c(*),d(*)
      real*8 t(5,5)
c
      call m_mul_5d(a,b,t)
      call m_mul_5d(t,c,d)
      return
      end
c
c ------------------------------------------------------------------
      subroutine m_transpose_5d(a,b)
      implicit real*8(a-h,o-z)
      real*8 a(5,5),b(5,5)
c      real*8 t(5,5)
c
      b=transpose(a)
c$$$      do 10 i=1,5
c$$$      do 10 j=1,5
c$$$ 10      t(i,j)=a(j,i)
c$$$      do 20 i=1,5
c$$$      do 20 j=1,5
c$$$ 20      b(j,i)=t(j,i)
      return
      end
c ------------------------------------------------------------------
      subroutine m_sim_tr_5d(a,b,c)
      implicit real*8(a-h,o-z)
      real*8 a(*),b(*),c(*)
      real*8 t(6,6)
c
      call m_transpose_5d(a,t)
      call m_mul3_5d(a,b,t,c)
      return
      end
c ------------------------------------------------------------------
      subroutine m_symp_set_6d(S)
      real*8 S(6,6)
c
c      call m_clear(S,6)
      s=0.d0
      S(1,2)=1.d0
      S(2,1)=-1.d0
      S(3,4)=1.d0
      S(4,3)=-1.d0
      S(5,6)=1.d0
      S(6,5)=-1.d0
      return
      end
c-----------------------------------------------------------------
      subroutine m_symmet_set_6d(p_in,blist)
      real*8 p_in(*),blist(*)
      k=1
      do i=1,6
       do j=1,6
         if(i.le.j) then
            blist(i+(j-1)*6)=p_in(k)
            k=k+1
         else
            blist(i+(j-1)*6)=blist(j+(i-1)*6)
         endif
       enddo
      enddo

      return
      end
c ------------------------------------------------------------------
      subroutine is_symplectic_6d(a)
      real*8 a(6,6),t(6,6),S(6,6)
      call m_symp_set_6d(S)
      call m_transpose_6d(a,t)
      call m_mul3_6d(t,S,a,t)
      write(*,*) ' Symplectic check'
      call m_print(t,6)
      return
      end
c ------------------------------------------------------------------
      subroutine m_B_set_6d(ax,bx,ay,by,az,bz,m_B)
      implicit real*8(a-h,o-z)
      real*8 ax,bx,ay,by,az,bz,m_B(6,6)
c
      call m_clear(m_B,6)
      sqrb=sqrt(bx)
      m_B(1,1)=sqrb
      m_B(2,1)=-ax/sqrb
      m_B(2,2)=1.d0/sqrb
      sqrb=sqrt(by)
      m_B(3,3)=sqrb
      m_B(4,3)=-ay/sqrb
      m_B(4,4)=1.d0/sqrb
      sqrb=sqrt(bz)
      m_B(5,5)=sqrb
      m_B(6,5)=-az/sqrb
      m_B(6,6)=1.d0/sqrb
      return
      end
c ------------------------------------------------------------------
      subroutine m_R_set_6d(r1,r2,r3,r4,m_R)
      implicit real*8(a-h,o-z)
      real*8 r1,r2,r3,r4,m_R(6,6)
      real*8 rmu
      rmu=sqrt(1.d0-r1*r4+r2*r3)
c
      call m_clear(m_R,6)
      m_R(1,1)=rmu
      m_R(2,2)=rmu
      m_R(3,3)=rmu
      m_R(4,4)=rmu
      m_R(5,5)=1.d0
      m_R(6,6)=1.d0
c
      m_R(3,1)=-r1
      m_R(3,2)=-r2
      m_R(4,1)=-r3
      m_R(4,2)=-r4
      m_R(1,3)=r4
      m_R(1,4)=-r2
      m_R(2,3)=-r3
      m_R(2,4)=r1
c
      return
      end
c ------------------------------------------------------------------
      subroutine m_H_set_6d(ex,epx,ey,epy,zx,zpx,zy,zpy,m_H)
      implicit real*8(a-h,o-z)
      real*8 ex,epx,ey,epy,zx,zpx,zy,zpy,m_H(6,6)
c
      call m_clear(m_H,6)
      detHx=zx*epx-zpx*ex
      detHy=zy*epy-zpy*ey
      a2=1.d0-detHx-detHy
      if(a2.gt.0.) then 
        a=sqrt(a2)
      else
        write(*,*) ' det Hx +det Hy >1 ',detHx,detHy
        stop
      endif
      ai=1.d0/(1.d0+a)
c 
      m_H(1,1)=1.d0-detHx*ai
      m_H(2,2)=m_H(1,1)
      m_H(3,3)=1.d0-detHy*ai
      m_H(4,4)=m_H(3,3)
      m_H(5,5)=a
      m_H(6,6)=m_H(5,5)
c
      m_H(1,5)=zx
      m_H(1,6)=ex
      m_H(2,5)=zpx
      m_H(2,6)=epx
      m_H(3,5)=zy
      m_H(3,6)=ey
      m_H(4,5)=zpy
      m_H(4,6)=epy
c
      m_H(5,1)=-epx
      m_H(5,2)=ex
      m_H(6,1)=zpx
      m_H(6,2)=-zx
      m_H(5,3)=-epy
      m_H(5,4)=ey
      m_H(6,3)=zpy
      m_H(6,4)=-zy
c
      m_H(1,3)=(ex*zpy-epy*zx)*ai
      m_H(1,4)=(ey*zx-ex*zy)*ai
      m_H(2,3)=(epx*zpy-epy*zpx)*ai
      m_H(2,4)=(ey*zpx-epx*zy)*ai
      m_H(3,1)=(ey*zpx-epx*zy)*ai
      m_H(3,2)=(-ey*zx+ex*zy)*ai
      m_H(4,1)=(epy*zpx-epx*zpy)*ai
      m_H(4,2)=(ex*zpy-epy*zx)*ai
      return
      end
c ------------------------------------------------------------------
      subroutine m_NtoP_set_6d(a,m_NtoP)
      implicit real*8(a-h,o-z)
      real*8 a(*),m_NtoP(6,6)
      real*8 t(6,6)
c
c      call m_clear(m_NtoP,6)
       if(a(2).le.0.) then
          write(*,*) ' Finite BETAX is required.' 
          write(*,*) ' BETAX is set to be 1..' 
          a(2)=1.d0
       endif
       if(a(5).le.0.) then
          write(*,*) ' Finite BETAY is required.' 
          write(*,*) ' BETAY is set to be 1..' 
          a(5)=1.d0
       endif
       if(a(27).le.0.) then
          write(*,*) ' Finite BETAZ is required.' 
          write(*,*) ' BETAZ is set to be 1..' 
          a(27)=1.d0
       endif
c
       call m_B_set_6d(a(1),a(2),a(4),a(5),a(22),a(27),m_NtoP)
       call m_R_set_6d(a(7),a(8),a(9),a(10),t)
       call m_mul_6d(t,m_NtoP,m_NtoP)
       call m_H_set_6d(a(11),a(12),a(13),a(14),
     +      a(15),a(16),a(17),a(18),t)
       call m_mul_6d(t,m_NtoP,m_NtoP)
c       write(*,*) ' m_NtoP_set_6d'
c       call m_print(m_NtoP,6)
      return
      end
c ------------------------------------------------------------------
      subroutine m_PtoN_set_6d(a,m_PtoN)
      implicit real*8(a-h,o-z)
      real*8 a(*),m_PtoN(6,6)
c
      call m_NtoP_set_6d(a,m_PtoN)
      call m_sym_inv_6d(m_PtoN,m_PtoN)
      return
      end
c ------------------------------------------------------------------
      subroutine m_sym_inv_6d(a,b)
      implicit real*8(a-h,o-z)
      real*8 a(6,6),b(6,6)
      real*8 S6(6,6)
c
      call m_symp_set_6d(S6)
      call m_sim_tr_6d(S6,a,b)
      return
      end 
c ------------------------------------------------------------------
      subroutine m_prot_set_6d(p,a)
      use macmath
      real*8 px,py,pz,p(*),a(6,6)
      call m_clear(a,6)
      px=p(3)*2.d0*pi
      py=p(6)*2.d0*pi
      pz=p(25)*2.d0*pi
      a(1,1)=cos(px)
      a(1,2)=sin(px)
      a(2,1)=-sin(px)
      a(2,2)=cos(px)
      a(3,3)=cos(py)
      a(3,4)=sin(py)
      a(4,3)=-sin(py)
      a(4,4)=cos(py)
      a(5,5)=cos(pz)
      a(5,6)=-sin(pz)
      a(6,5)=sin(pz)
      a(6,6)=cos(pz)
      return
      end
c ------------------------------------------------------------------
      subroutine m_Emit_set_6d(emitx,emity,emitz,m_Emit)
      implicit real*8(a-h,o-z)
      real*8 emitx,emity,emitz,m_Emit(6,6)
c
      call m_clear(m_Emit,6)
      m_Emit(1,1)=emitx
      m_Emit(2,2)=emitx
      m_Emit(3,3)=emity
      m_Emit(4,4)=emity
      m_Emit(5,5)=emitz
      m_Emit(6,6)=emitz
      return
      end
c ------------------------------------------------------------------
      subroutine m_Crs_set_6d(x_angle,m_Crs)
      implicit real*8(a-h,o-z)
      real*8 x_angle,m_Crs(6,6)
c
      call m_clear(m_Crs,6)
      cosxi=1.d0/cos(x_angle)
      tanx=tan(x_angle)
      m_Crs(1,1)=1.d0
      m_Crs(2,2)=cosxi
      m_Crs(3,3)=1.d0
      m_Crs(4,4)=cosxi
      m_Crs(5,5)=cosxi
      m_Crs(6,6)=1.d0
c
      m_Crs(1,5)=tanx
      m_Crs(6,2)=-tanx
      return
      end
c ------------------------------------------------------------------
      subroutine m_inverse_6d(a,b)
      real*8 a(6,6),b(6,6),t(6,6)
      integer indx(6)
c
      t=a
c      do 10 i=1,6
c      do 10 j=1,6
c 10     t(j,i)=a(j,i)
      call m_inverse(t,6,6,indx,b)
      return
      end
c ------------------------------------------------------------------
      subroutine m_zslice_5d(x,a,b,c)
      real*8 x(6,6),a(5,5),b(5),c
c
      a(1:4,1:4)=x(1:4,1:4)
c      do 10 j=1,4
c      do 10 i=1,4
c 10      a(i,j)=x(i,j)
      do i=1,4
        b(i)=x(5,i)
        a(5,i)=x(6,i)
        a(i,5)=x(i,6)
      enddo
      c=x(5,5)
      b(5)=x(5,6)
      a(5,5)=x(6,6)
      return
      end
c ------------------------------------------------------------------
      subroutine m_normalize_5d(a,b)
      real*8 a(5,5),b(10)
      real*8 ta(5,5),tb(10),rnorm,remax,absa
      integer naxis(5),ikeep(5)
      ikeep=0
c      do 5 i=1,5
c 5       ikeep(i)=0
c
      do i=1,5
        rnorm=0.d0
        remax=0.
        naxis(i)=0
        do j=1,5
          rnorm=rnorm+a(j,i)**2
        enddo
        rnorm=1.d0/sqrt(rnorm)
        do j=1,5
          ta(j,i)=a(j,i)*rnorm
          absa=abs(ta(j,i))
          if(absa.gt.remax.and.ikeep(j).eq.0) then
            remax=absa
            naxis(i)=j
          endif
        enddo
        ikeep(naxis(i))=1
      enddo
      do i=1,5
        if(ikeep(i).eq.0) then
          write(*,*) 'matrix.f : Normal axis is determined'
          write(*,*) naxis
          stop
        endif
      enddo
      do i=1,5
        k=naxis(i)
        do j=1,5
          a(j,k)=ta(j,i)
        enddo
        tb(2*k-1)=b(2*i-1)
        tb(2*k)=b(2*i)
      enddo
      b=tb
c      do 60 i=1,10
c 60     b(i)=tb(i)
      return
      end
c ------------------------------------------------------------------
      subroutine m_normalize_6d(a,b)
      real*8 a(6,6),b(12)
      real*8 ta(6,6),tb(12),rnorm,remax,absa
      integer naxis(6),ikeep(6)
      ikeep=0
c      do 5 i=1,6
c 5       ikeep(i)=0
c
      do i=1,6
        rnorm=0.d0
        remax=0.
        naxis(i)=0
        do j=1,6
          rnorm=rnorm+a(j,i)**2
        enddo
        rnorm=1.d0/sqrt(rnorm)
        do j=1,6
          ta(j,i)=a(j,i)*rnorm
          absa=abs(ta(j,i))
          if(absa.gt.remax.and.ikeep(j).eq.0) then
            remax=absa
            naxis(i)=j
          endif
        enddo
        ikeep(naxis(i))=1
      enddo
      do i=1,6
        if(ikeep(i).eq.0) then
          write(*,*) 'matrix.f : Normal axis is determined'
          write(*,*) naxis
          stop
        endif
      enddo
      do i=1,6
        k=naxis(i)
        do j=1,6
          a(j,k)=ta(j,i)
        enddo
        tb(2*k-1)=b(2*i-1)
        tb(2*k)=b(2*i)
      enddo
      b=tb
c      do 60 i=1,12
c 60     b(i)=tb(i)
        return
        end
c ------------------------------------------------------------------
      subroutine m_diag_5d(a1,a2,a3,a4,a5,b)
      real*8 a1,a2,a3,a4,a5,b(5,5)
      call m_clear(b,5)
      b(1,1)=a1
      b(2,2)=a2
      b(3,3)=a3
      b(4,4)=a4
      b(5,5)=a5
      return
      end
c ------------------------------------------------------------------
      subroutine mv_mul_5d(a,b,c)
      real*8 a(5,5),b(5),c(5)
      c=matmul(a,b)
c      do 10 i=1,5
c       t(i)=0.d0
c      do 10 j=1,5
c 10    t(i)=t(i)+a(i,j)*b(j)
c      do 20 i=1,5
c 20      c(i)=t(i)
      return
      end
c ------------------------------------------------------------------
      subroutine mv_mul_6d(a,b,c)
      real*8 a(6,6),b(6),c(6)
      c=matmul(a,b)
c$$$      do 10 i=1,6
c$$$       t(i)=0.d0
c$$$      do 10 j=1,6
c$$$ 10    t(i)=t(i)+a(i,j)*b(j)
c$$$      do 20 i=1,6
c$$$ 20      c(i)=t(i)
      return
      end
c ------------------------------------------------------------------
      subroutine benv_norm_axis(env,vec,t_angle)
      implicit real*8(a-h,o-z)
      real*8 env(5,5),vec(5),t_angle
      real*8 u(5,5),tmp(5,5),xx(5)
      t_angle=atan2(2.d0*env(1,3),env(1,1)-env(3,3))*0.5d0
      cost=cos(t_angle)
      sint=sin(t_angle)
      call m_clear(u,5)
      u(1,1)=cost
      u(2,2)=cost
      u(3,3)=cost
      u(4,4)=cost
      u(5,5)=1.d0
      u(1,3)=sint
      u(2,4)=sint
      u(3,1)=-sint
      u(4,2)=-sint
c
      call m_clear(tmp,5)
      do i=1,5
        do j=1,5
          do k=1,5
            tmp(i,j)=tmp(i,j)+env(i,k)*u(j,k)
          enddo
        enddo
      enddo
      call m_clear(env,5)
      do i=1,5
        do j=1,5
          do k=1,5
            env(i,j)=env(i,j)+u(i,k)*tmp(k,j)
          enddo
        enddo
      enddo
      do i=1,5
        xx(i)=0.
        do j=1,5
          xx(i)=xx(i)+u(i,j)*vec(j)
        enddo
      enddo
      vec=xx
c      do 40 i=1,5
c 40      vec(i)=xx(i)
      return
      end
c ------------------------------------------------------------------
      subroutine m_inverse(a,n,np,indx,y)
      integer n,np,indx(np)
      real*8 a(np,np),y(np,np),d
      call m_unit(y,np)
      call ludcmp(a,n,np,indx,d)
      do 10 j=1,n
      call lubksb(a,n,np,indx,y(1,j))
  10  continue
      return
      end
c ------------------------------------------------------------------
      subroutine ludcmp(a,n,np,indx,d)
      implicit real*8(a-h,o-z)
      parameter (nmax=100,tiny=1.e-20)
      integer n,np,indx(np)
      real*8 a(np,np),d,vv(nmax)
c     begin initialize for preventing compiler warning
      imax=0
c     end   initialize for preventing compiler warning
      d=1.d0
      do 10 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
  11    continue
        if(aamax.eq.0.d0) WRITE(*,*) ' Singular matrix '
        vv(i)=1./aamax
  10  continue
      do 19 j=1,n
        if(j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if(i.gt.1) then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
  13          continue
              a(i,j)=sum
            endif
  14      continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if(j.gt.1) then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
  15        continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if(dum.gt.aamax) then
            imax=i
            aamax=dum
          endif
  16    continue
        if(j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
  17      continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n) then
          if(a(j,j).eq.0.d0) a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
  18      continue
        endif
  19  continue
      if(a(n,n).eq.0.d0) a(n,n)=tiny
      return
      end
c ------------------------------------------------------------------
      subroutine lubksb(a,n,np,indx,b)
      implicit real*8(a-h,o-z)
      integer n,np,indx(np)
      real*8 a(np,np),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
  11      continue
        else if(sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
  12  continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n) then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
  13      continue
        endif
        b(i)=sum/a(i,i)
  14  continue
      return
      end
