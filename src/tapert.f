      subroutine tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl,np,kturn,
     $     ax,ay,dx,dy,
     $     xl,yl,xh,yh,pxj,pyj,dpj,theta)
      use kyparam
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:idelc,idtypec,latt
      use tracklim
      implicit none
      integer, parameter :: nkptbl = 6
      real*8 ,intent(inout):: x(np0),px(np0),y(np0),py(np0),z(np0),
     $     g(np0),dv(np0),sx(np0),sy(np0),sz(np0)
      integer*4 ,intent(inout):: kptbl(np0,nkptbl),np
      integer*4 ,intent(in):: kturn
      real*8 ,intent(in):: ax,ay,dx,dy,xl,yl,xh,yh,pxj,pyj,dpj,theta
      real*8 xh1,xl1,yh1,yl1,ax1,ay1,xa,ya,cost,sint,
     $     x1,px1,y1,py1,z1,g1,dv1,phi(2),sx1,sy1,sz1
      integer*4 i,j,k,kptmp(nkptbl)
      logical eli, dodrop

      call limitnan(x(1:np),xlimit)
      call limitnan(px(1:np),plimit)
      call limitnan(y(1:np),xlimit)
      call limitnan(py(1:np),plimit)
      call limitnan(z(1:np),zlimit)
c     Shortcut case: pxj != 0.0 || pyj != 0.0
      if(pxj .ne. 0.d0 .or. pyj .ne. 0.d0)then
        xa=min(abs(xh),abs(xl))
        ya=min(abs(yh),abs(yl))
        if(xa .ne. 0.d0)then
          ax1 = 1.d0 / xa
        else
          ax1 = 0.d0
        endif
        if(ya .ne. 0.d0)then
          ay1 = 1.d0 / ya
        else
          ay1 = 0.d0
        endif
        do i=1,np
          if(((ax1*x(i))**2 + (ay1*y(i))**2) .gt. 1.d0)then
            call tran_array(phi,2)
            phi(1) = pi2 * phi(1)
            phi(2) = pi2 * phi(2)
            x(i)   = xa  * cos(phi(1))
            y(i)   = ya  * sin(phi(1))
            px(i)  = pxj * cos(phi(2))
            py(i)  = pyj * sin(phi(2))
c            dp=g(i)*(2.d0+g(i))+dpj
            g(i)=g(i)+dpj
c            g(i)=dp/(1.d0+sqrt(1.d0+dp))
          endif
        enddo
        return
      endif

c     Shortcut case: twake || lwake
c      if(twake .or. lwake)then
c        do i=1,np
c          x(i)=min(xh,max(xl,x(i)))
c          y(i)=min(yh,max(yl,y(i)))
c        enddo
c        return
c      endif

c     General aperture case:
c
c     Don't flip aperture conditional!!
c     This conditional is designed to drop NaN particle.
c     In IEEE754 standard, any comparision operation
c     with NaN operand is defined as ``False''.
c
      eli = (ax .ne. 0.d0) .and. (ay .ne. 0.d0)
      if(eli .and. ((xh .eq. xl) .or. (yh .eq. yl)))then
        xl1 = -1.d100
        xh1 =  1.d100
        yl1 = -1.d100
        yh1 =  1.d100
      else
        xl1 = xl
        xh1 = xh
        yl1 = yl
        yh1 = yh
      endif
c      write(*,'(a,i5,1p8g14.7)')'tapert-a ',l,
c     $     xl,xh,yl,yh,xl1,xh1,yl1,yh1

      ax1 = 0.d0
      ay1 = 0.d0
      if(ax .ne. 0.d0)then
        ax1 = 1.d0 / ax
      endif
      if(ay .ne. 0.d0)then
        ay1 = 1.d0 / ay
      endif

      if(theta .ne. 0) then
        cost = cos(theta)
        sint = sin(theta)
      else
        cost = 1.d0
        sint = 0.d0
      endif

c     Marking drop particles
      dodrop = .false.
      if(eli)then
        if(theta .ne. 0.d0)then
c     Case: eli && theta != 0.0
          do i=1,np
            xa = cost * (x(i) - dx) - sint * (y(i) - dy)
            ya = sint * (x(i) - dx) + cost * (y(i) - dy)
            if(.not. (
     $           ((ax1 * xa)**2 + (ay1 * ya)**2 .le. 1.d0) .and.
     $           (      xl1 .lt. xa .and. xa .lt. xh1
     $           .and.  yl1 .lt. ya .and. ya .lt. yh1) .and.
     $           abs(z(i)) .le. zlost))then
              dodrop = .true.
              kptbl(i,4) = l_track
              kptbl(i,5) = kturn
            endif
          enddo
        else
c     Case: eli && !(theta != 0.0)
          do i=1,np
            xa = (x(i) - dx)
            ya = (y(i) - dy)
            if(.not. (
     $           ((ax1 * xa)**2 + (ay1 * ya)**2 .le. 1.d0) .and.
     $           (      xl1 .lt. xa .and. xa .lt. xh1
     $           .and.  yl1 .lt. ya .and. ya .lt. yh1) .and.
     $           abs(z(i)) .le. zlost))then
c              write(*,'(a,i5,1p8g14.7)')'tapert ',i,
c     $             xa,ya,ax1,ay1,xl1,xh1,yl1,yh1
              dodrop = .true.
              kptbl(i,4) = l_track
              kptbl(i,5) = kturn
            endif
          enddo
        endif
      else
        if(theta .ne. 0.d0)then
c     Case: !eli && (theta != 0.0)
          do i=1,np
            xa = cost * (x(i) - dx) - sint * (y(i) - dy)
            ya = sint * (x(i) - dx) + cost * (y(i) - dy)
            if(.not. (
     $           (      xl1 .lt. xa .and. xa .lt. xh1
     $           .and.  yl1 .lt. ya .and. ya .lt. yh1) .and.
     $           abs(z(i)) .le. zlost))then
              dodrop = .true.
              kptbl(i,4) = l_track
              kptbl(i,5) = kturn
            endif
          enddo
        else
c     Case: !eli && !(theta != 0.0)
          do i=1,np
            xa = (x(i) - dx)
            ya = (y(i) - dy)
            if(.not. (
     $           (      xl1 .lt. xa .and. xa .lt. xh1
     $           .and.  yl1 .lt. ya .and. ya .lt. yh1) .and.
     $           abs(z(i)) .le. zlost))then
c              write(*,'(a,i5,1p8g14.7)')'tapert ',i,
c     $             xa,ya,ax1,ay1,xl1,xh1,yl1,yh1
              dodrop = .true.
              kptbl(i,4) = l_track
              kptbl(i,5) = kturn
            endif
          enddo
        endif
      endif

c     Shortcut case: Lossless
      if(.not. dodrop)then
        return
      endif

c     Reporting drop particles
      if( (.not. dapert)
     $     .and. (.not. trpt .or. idtypec(l_track) .eq. icAprt)
     $     .and. (outfl .ne. 0))then
        call tapert_report_dropped(outfl,kturn,l_track,
     $       latt,np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl)
      endif

c     Sweaping drop particles
      i=1
      do while(i .le. np)
        if(kptbl(i,4) .ne. 0)then
c     Search alive particle from tail: (i, np]
          do while(i .lt. np .and. (kptbl(np,4) .ne. 0))
            np=np-1
          enddo
          if(kptbl(np,4) .eq. 0)then
c     Swap drop particlen slot[i] with tail alive particle slot[np]
            j=kptbl(np,2)
            k=kptbl(i, 2)
c      - Update maps between partice ID and array index
            kptbl(k, 1)=np
            kptbl(j, 1)=i
            kptbl(np,2)=k
            kptbl(i, 2)=j
c      - Swap kptbl except forward/backward[kptbl(*,1)/kptbl(*,2)]
            kptmp(   3:nkptbl) = kptbl(np,3:nkptbl)
            kptbl(np,3:nkptbl) = kptbl(i, 3:nkptbl)
            kptbl(i, 3:nkptbl) = kptmp(   3:nkptbl)
c      - Swap particle coordinates
            x1  = x(i)
            px1 = px(i)
            y1  = y(i)
            py1 = py(i)
            z1  = z(i)
            g1  = g(i)
            dv1 = dv(i)
            sx1 = sx(i)
            sy1 = sy(i)
            sz1 = sz(i)

            x(i)  = x(np)
            px(i) = px(np)
            y(i)  = y(np)
            py(i) = py(np)
            z(i)  = z(np)
            g(i)  = g(np)
            dv(i) = dv(np)
            sx(i) = sx(np)
            sy(i) = sy(np)
            sz(i) = sz(np)

            x(np)  = x1
            px(np) = px1
            y(np)  = y1
            py(np) = py1
            z(np)  = z1
            g(np)  = g1
            dv(np) = dv1
            sx(np) = sx1
            sy(np) = sy1
            sz(np) = sz1
c            if(ktfenanq(sy(np)))then
c              write(*,'(a,2i5,1p8g12.4)')'tapert-i ',i,np,
c     $           sx(i),sy(i),sz(i),sx1,sy1,sz1
c            endif
          endif
          np=np-1
        endif
        i=i+1
      enddo
      return
      end

      subroutine tapert1(x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl,np,kturn)
      use kyparam
      use tfstk
      use tmacro
      use ffs_pointer, only:latt
      implicit none
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      real*8 sx(np0),sy(np0),sz(np0)
      integer*4 kptbl(np0,6),np,kturn
      integer*8 lp
      real*8 dpxj,dpyj,ddp,dx1,dy1,dx2,dy2
      lp=latt(l_track)
      dpxj=rlist(lp+ky_JDPX_Aprt)
      dpyj=rlist(lp+ky_JDPY_Aprt)
      ddp=rlist(lp+ky_DP_Aprt)
      dx1=min(rlist(lp+ky_DX1_Aprt),
     $        rlist(lp+ky_DX2_Aprt))
      dx2=max(rlist(lp+ky_DX1_Aprt),
     $        rlist(lp+ky_DX2_Aprt))
      dy1=min(rlist(lp+ky_DY1_Aprt),
     $        rlist(lp+ky_DY2_Aprt))
      dy2=max(rlist(lp+ky_DY1_Aprt),
     $        rlist(lp+ky_DY2_Aprt))
c      write(*,*)'tapert1 ',l,
c     $     rlist(lp+ky_AX_Aprt),
c     $     rlist(lp+ky_AY_Aprt),
c     $     rlist(lp+ky_DX_Aprt),
c     $     rlist(lp+ky_DY_Aprt)
      call tapert(x,px,y,py,z,g,dv,sx,sy,sz,
     1     kptbl,np,kturn,
     $     rlist(lp+ky_AX_Aprt),
     $     rlist(lp+ky_AY_Aprt),
     $     rlist(lp+ky_DX_Aprt),
     $     rlist(lp+ky_DY_Aprt),
     $     dx1,dy1,dx2,dy2,dpxj,dpyj,ddp,
     $     rlist(lp+ky_ROT_Aprt))
      return
      end

c     Helper functions for Aperture Handling in Tracking Modules

c     Report new drop marked particles in alive area [1, np]
      subroutine tapert_report_dropped(outfd,kturn,lbegin,
     $     latt,np,x,px,y,py,z,g,dv,sx,sy,sz,kptbl)
      use tfstk
      use tmacro
      use ffs_pointer, only:idvalc,idtypec,idelc
      implicit none
      integer*4 outfd,kturn,lbegin
      integer*8 latt(nlat)
      integer*4 np
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      real*8 sx(np0),sy(np0),sz(np0)
      integer*4 kptbl(np0,6)
      integer*4 i,l,t
      character*2 ord

      do i=1,np
         l = kptbl(i,4)
         t = kptbl(i,5)
         if((l .ge. lbegin) .and. (l .gt. 0))then
           if(l .le. nlat)then
             write(outfd,'(1x,''P. '',i5,'' lost in '',i5,a,'' turn'',
     $'' at '',i5,''('',a,''), amplitudes:'',3(1x,g13.7))')
     $            kptbl(i,2), t, ord(t),l,
     $            pname(idelc(l))(1:lpname(idelc(l))),
     $            x(i),y(i),z(i)
           else
             write(outfd,'(1x,''P. '',i5,'' lost in '',i5,a,'' turn'',
     $ '' at '',i5,'', amplitudes:'',3(1x,g13.7))')
     $            kptbl(i,2), t, ord(t), l,
     $            x(i),y(i),z(i)
           endif
         endif
      enddo
      return
      end subroutine tapert_report_dropped

c     Sweep new drop marked particles from alive area [1, np]
c     This subroutine need exclusive access to given arguments
      subroutine tapert_sweep_dropped(np0,np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl)
      implicit none
      integer, parameter :: nkptbl = 6
      integer(4), intent(in)    :: np0
      integer(4), intent(inout) :: np
      real(8),    intent(inout) :: x(np0), px(np0), y(np0), py(np0),
     $     z(np0), g(np0), dv(np0), sx(np0),sy(np0),sz(np0)
      integer(4), intent(inout) :: kptbl(np0,nkptbl)
      integer(4) :: i, j, k, kptmp(nkptbl)
      real(8) :: x1, px1, y1, py1, z1, g1, dv1, sx1, sy1, sz1

c     Scan new drop marked particles from alive area [1, np]
      i = 1
      do while(i .le. np)
         if(kptbl(i,4) .ne. 0)then
c           Search alive particle from tail: (i, np]
            do while((i .lt. np) .and. (kptbl(np,4) .ne. 0))
               np = np - 1
            enddo
            if(kptbl(np,4) .eq. 0)then
c              Swap dropped particlen slot[i] with tail alive particle slot[np]
               j = kptbl(np,2)
               k = kptbl(i, 2)
c              Update maps between partice ID and array index
               kptbl(k, 1) = np
               kptbl(j, 1) = i
               kptbl(np,2) = k
               kptbl(i, 2) = j
c              Swap kptbl except forward/backward[kptbl(*,1)/kptbl(*,2)]
               kptmp(   3:nkptbl) = kptbl(np,3:nkptbl)
               kptbl(np,3:nkptbl) = kptbl(i, 3:nkptbl)
               kptbl(i, 3:nkptbl) = kptmp(   3:nkptbl)
c              Swap particle coordinates
               x1  = x (i)
               px1 = px(i)
               y1  = y (i)
               py1 = py(i)
               z1  = z (i)
               g1  = g (i)
               dv1 = dv(i)
               sx1 = sx(i)
               sy1 = sy(i)
               sz1 = sz(i)

               x (i) = x (np)
               px(i) = px(np)
               y (i) = y (np)
               py(i) = py(np)
               z (i) = z (np)
               g (i) = g (np)
               dv(i) = dv(np)
               sx(i) = sx(np)
               sy(i) = sy(np)
               sz(i) = sz(np)

               x (np) = x1
               px(np) = px1
               y (np) = y1
               py(np) = py1
               z (np) = z1
               g (np) = g1
               dv(np) = dv1
               sx(np) = sx1
               sy(np) = sy1
               sz(np) = sz1
            endif
            np = np - 1
         endif
         i = i + 1
      enddo
      return
      end subroutine tapert_sweep_dropped

c     Sweep new inject marked particles from dead area (np, np0]
c     This subroutine need exclusive access to given arguments
      subroutine tapert_sweep_injected(np0,np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     kptbl)
      implicit none
      integer, parameter :: nkptbl = 6
      integer(4), intent(in)    :: np0
      integer(4), intent(inout) :: np
      real(8),    intent(inout) :: x(np0), px(np0), y(np0), py(np0),
     $     z(np0), g(np0), dv(np0), sx(np0),sy(np0),sz(np0)
      integer(4), intent(inout) :: kptbl(np0,nkptbl)
      integer(4) :: i, j, k, m, kptmp(nkptbl)
      real(8) :: x1, px1, y1, py1, z1, g1, dv1, sx1, sy1, sz1

c     Scan dead particles from dead area (np, m = np0]
      i = np + 1
      m = np0
      do while(i .le. m)
         if(kptbl(i,4) .ne. 0)then
c           Search injected particle from tail: (i, m]
            do while((i .lt. m) .and. (kptbl(m,4) .ne. 0))
               m = m - 1
            enddo
            if(kptbl(m,4) .eq. 0)then
c              Swap dead particlen slot[i] with tail injected particle slot[m]
               j = kptbl(m,2)
               k = kptbl(i,2)
c              Update maps between partice ID and array index
               kptbl(k,1) = m
               kptbl(j,1) = i
               kptbl(m,2) = k
               kptbl(i,2) = j
c              Swap kptbl except forward/backward[kptbl(*,1)/kptbl(*,2)]
               kptmp(  3:nkptbl) = kptbl(m,3:nkptbl)
               kptbl(m,3:nkptbl) = kptbl(i,3:nkptbl)
               kptbl(i,3:nkptbl) = kptmp(  3:nkptbl)
c              Swap particle coordinates
               x1  = x (i)
               px1 = px(i)
               y1  = y (i)
               py1 = py(i)
               z1  = z (i)
               g1  = g (i)
               dv1 = dv(i)
               sx1 = sx(i)
               sy1 = sy(i)
               sz1 = sz(i)

               x (i) = x (m)
               px(i) = px(m)
               y (i) = y (m)
               py(i) = py(m)
               z (i) = z (m)
               g (i) = g (m)
               dv(i) = dv(m)
               sx(i) = sx(m)
               sy(i) = sy(m)
               sz(i) = sz(m)

               x (m) = x1
               px(m) = px1
               y (m) = y1
               py(m) = py1
               z (m) = z1
               g (m) = g1
               dv(m) = dv1
               sx(m) = sx1
               sy(m) = sy1
               sz(m) = sz1
            endif
            m = m - 1
         endif
         i = i + 1
      enddo
      np = m
      return
      end subroutine tapert_sweep_injected
