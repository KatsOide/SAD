      module tfcc
      use tfstk
      real*8 ,private::xmin,xmax,ymin,ymax
      integer*4 ,private::isp0

      contains

      subroutine tfcanvasclip(isp1,kx,irtc)
      use mathfun, only:outer2
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_rlist), pointer :: klr
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 ka,ka1,ka2,kal,kadash,kaxi
      integer*4 i,j,na,np,n,
     $     itfmessage,narg,ispa,ndash,k,isp2,isp3,ignore
      real*8 x0,x1,y0,y1,xa,ya,xb,yb,
     $     ux,uy,t(4),dmax,dash,s,s1,v,xm,ym,x,y,ds,
     $     z(2),z0(2),z1(2),dz(2)
      real*8 ,parameter::eps=1.d-6
      logical*4 pin,in0,in1
      pin(x,y)=x >=xmin .and. x <= xmax .and. y >= ymin .and. y <= ymax
      narg=isp-isp1
      if(narg /= 2 .and. narg /= 3)then
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
      endif
      ndash=0
      kadash=0
      xm=0.d0
      ym=0.d0
      if(narg == 3)then
        ka=ktfaddr(ktastk(isp))
        if(tfreallistq(ktastk(isp)))then
          ndash=ilist(2,ka-1)
          if(ndash > 1)then
            kadash=ktfaddr(klist(ka+1))
          endif
        elseif(ktfrealq(ktastk(isp)))then
          v=rtastk(isp)
          if(v > 0.d0)then
            xm=v
          elseif(v < 0.d0)then
            ym=-v
          endif
        endif
      endif
      if(.not. tflistq(ktastk(isp1+2)))then
c        write(*,*)'canvasclip-1'
        irtc=itfmessage(9,'General::wrongtype','"List for #2"')
        return
      endif
      kal=ktfaddr(ktastk(isp1+2))
      if(ilist(2,kal-1) /= 2)then
c        write(*,*)'canvasclip-2'
        irtc=itfmessage(9,'General::wrongleng','"#2","2"')
        return
      endif
      if(ktfreallistq(kal))then
c        write(*,*)'canvasclip-3'
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List for #2[[1]]"')
        return
      endif
      ka1=ktfaddr(klist(kal+1))
      ka2=ktfaddr(klist(kal+2))
      if(ktfnonreallistq(ka1) .or. ktfnonreallistq(ka2))then
c        write(*,*)'canvasclip-4'
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of Reals for #2[[1]] and #2[[2]]"')
        return
      endif
      xmin=anint(rlist(ka1+1))
      ymin=anint(rlist(ka1+2))
      xmax=anint(rlist(ka2+1))
      ymax=anint(rlist(ka2+2))
      isp0=isp
      ka=ktfaddr(ktastk(isp1+1))
      if(.not. tflistq(ktastk(isp1+1)))then
        go to 9100
      endif
      if(ilist(2,ka-1) /= 2)then
        go to 9100
      endif
      if(ktfreallistq(ka))then
        go to 9100
      endif
      ka1=ktfaddr(klist(ka+1))
      ka2=ktfaddr(klist(ka+2))
      if(.not. tfreallistq(klist(ka+1)) .or.
     $     .not. tfreallistq(klist(ka+2)))then
        go to 9100
      endif
      np=ilist(2,ka1-1)
      if(ilist(2,ka2-1) /= np)then
        go to 9100
      endif
      isp0=isp
      x0=anint(rlist(ka1+1))
      y0=anint(rlist(ka2+1))
      in0=pin(x0,y0)
      call tfccput(x0,y0)
      do n=2,np
        x1=anint(rlist(ka1+n))
        y1=anint(rlist(ka2+n))
        in1=pin(x1,y1)
        ux=x1-x0
        uy=y1-y0
        t(1:2)=merge([-1.d0,-1.d0],([xmin,xmax]-x0)/ux,ux == 0.d0)
        t(3:4)=merge([-1.d0,-1.d0],([ymin,ymax]-y0)/uy,uy == 0.d0)
        do i=1,4
          if(t(i) >= 0.d0 .and. t(i) <= 1.d0)then
            xa=anint(x0+t(i)*ux)
            ya=anint(y0+t(i)*uy)
            if(pin(xa,ya))then
              do j=i+1,4
                if(t(j) >= 0.d0 .and. t(j) <= 1.d0)then
                  xb=anint(x0+t(j)*ux)
                  yb=anint(y0+t(j)*uy)
                  if(pin(xb,yb))then
                    if(t(i) < t(j))then
                      call tfccput(xa,ya)
                      call tfccput(xb,yb)
                    else
                      call tfccput(xb,yb)
                      call tfccput(xa,ya)
                    endif
                    go to 1
                  endif
                endif
              enddo
              if(.not. in0)then
                call tfccput(x0,y0)
                call tfccputa(x0,y0,x1,y1)
              endif
              isp=isp+2
              rtastk(isp-1:isp)=[xa,ya]
              exit
            endif
          endif
        enddo
        if(in1)then
          call tfccput(x1,y1)
        elseif(.not. in0)then
          call tfccput(x0,y0)
          call tfccputa(x0,y0,x1,y1)
        endif
 1      x0=x1
        y0=y1
        in0=in1
      enddo
      call tfccput(x0,y0)
      na=isp-isp0
      if(na == 2)then
        isp=isp+2
        rtastk(isp-1:isp)=rtastk(isp-3:isp-2)+[xm,ym]
        na=4
      elseif(na > 4)then
        ignore=0
        ispa=isp0+3
        z0=rtastk(isp0+1:isp0+2)
        z1=rtastk(isp0+3:isp0+4)
        do i=isp0+5,isp-1,2
          z=rtastk(i:i+1)
          ds=1.d0
          if(ignore < 4)then
            if(dot_product(z-z1,z1-z0) >= 0)then
              ds=outer2(z-z0,z1-z0)**2/(1.d-4+dot_product(z-z0,z-z0))
            endif
          endif
          if(ds < 0.16d0)then
            rtastk(ispa:ispa+1)=z
            ignore=ignore+1
          else
            rtastk(ispa:ispa+1)=z1
            z0=z1
            ispa=ispa+2
            ignore=0
          endif
          z1=z
        enddo
        rtastk(ispa:ispa+1)=z1
        isp=ispa+1
c        write(*,*)na,isp-isp0
        na=isp-isp0
      endif
      if(kadash /= 0)then
        isp2=isp
        if(na > 1)then
          isp3=isp
          dash=0.d0
          k=1
          dmax=max(abs(rlist(kadash+1)),eps)
          z0=rtastk(isp0+1:isp0+2)
          isp=isp+2
          rtastk(isp-1:isp)=z0
          do i=isp0+3,isp2-1,2
            z1=rtastk(i:i+1)
 11         dz=z1-z0
            s=hypot(dz(1),dz(2))
            if(s > 0.d0)then
              s1=s+dash
              if(dash >= 0.d0)then
                if(s1 <= dmax)then
                  isp=isp+2
                  rtastk(isp-1:isp)=z1
                  dash=s1
                  z0=z1
                else
                  z0=z0+dz*(dmax-dash)/s
                  isp=isp+2
                  rtastk(isp-1:isp)=anint(z0)
                  kaxi=ktavaloc(-1,isp-isp3)
                  klist(kaxi+1:kaxi+isp-isp3)=ktastk(isp3+1:isp)
                  isp3=isp3+1
                  ktastk(isp3)=ktflist+kaxi
                  isp=isp3
                  k=mod(k,ndash)+1
                  dash=min(-abs(rlist(kadash+k)),-eps)
                  go to 11
                endif
              else
                if(s1 > 0.d0)then
                  z0=z0-dz*dash/s
                  isp=isp+2
                  rtastk(isp-1:isp)=anint(z0)
                  dash=0.d0
                  k=mod(k,ndash)+1
                  dmax=max(abs(rlist(kadash+k)),eps)
                  go to 11
                elseif(s1 == 0.d0)then
                  isp=isp+2
                  rtastk(isp-1:isp)=anint(z1)
                  dash=0.d0
                  k=mod(k,ndash)+1
                  dmax=max(abs(rlist(kadash+k)),eps)
                  z0=z1
                else
                  dash=s1
                  z0=z1
                endif
              endif
            endif
          enddo
          if(dash > 0.d0)then
            if(isp > isp3)then
              kaxi=ktavaloc(-1,isp-isp3)
              klist(kaxi+1:kaxi+isp-isp3)=ktastk(isp3+1:isp)
              isp3=isp3+1
              ktastk(isp3)=ktflist+kaxi
            endif
          endif
          isp=isp3
        endif
        if(isp > isp2)then
          kx=kxmakelist(isp2)
        else
          kx=dxnulll
        endif
      else
        if(na <= 1)then
          kx=dxnulll
        else
          kx=kxavaloc(-1,na,klr)
          klr%rbody(1:na)=rtastk(isp0+1:isp0+na)
        endif
      endif
      isp=isp0
      irtc=0
      return
 9100 irtc=itfmessage(9,'General::wrongtype',
     $     '"{{x1, x2, ..}, {y1, y2, ..}}"')
c      write(*,*)'canvasclip-4'
      isp=isp0
      return
      end subroutine

      subroutine tfccput(x,y)
      implicit none
      real*8 ,intent(in):: x,y
      real*8 xl,yl
      xl=max(xmin,min(xmax,x))
      yl=max(ymin,min(ymax,y))
      if(isp > isp0 .and. rtastk(isp-1) == xl .and. rtastk(isp) == yl)then
        return
      endif
      isp=isp+2
      rtastk(isp-1:isp)=[xl,yl]
      return
      end subroutine

      subroutine tfccputa(x0,y0,x1,y1)
      implicit none
      real*8 ,intent(in):: x0,y0,x1,y1
      real*8 xm(2),ym(2)
      integer*4 i
      xm=merge([xmin,xmax],[xmax,xmin],x0<=x1)
      ym=merge([ymin,ymax],[ymax,ymin],y0<=y1)
      do i=1,2
        if((y0-ym(i))*(y1-ym(i)) <= 0.d0)then
          if(x0 <= xmin)then
            call tfccput(xmin,ym(i))
          elseif(x0 >= xmax)then
            call tfccput(xmax,ym(i))
          endif
        endif
        if((x0-xm(i))*(x1-xm(i)) <= 0.d0)then
          if(y0 <= ymin)then
            call tfccput(xm(i),ymin)
          elseif(y0 >= ymax)then
            call tfccput(xm(i),ymax)
          endif
        endif
      enddo
      return
      end subroutine

      end module

      subroutine tfcanvas3dcliptriangle(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,klxi
      type (sad_rlist), pointer :: klxj
      integer*8 ka,kas,kac,kai,kavc,kavs
      integer*4 isp1,irtc,iav(3),nt,m,itfmessage,isp0,j,i,ii,m3
      real*8 xs,ys,zs,xc,yc,zc
      if(isp /= isp1+3)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      kas=ktfaddr(ktastk(isp1+2))
      if(.not. tfnumlistqn(dtastk(isp1+2),3))then
        go to 9000
      endif
      kavs=ktfaddr(klist(kas+1))
      kac=ktfaddr(ktastk(isp1+3))
      if(.not. tfnumlistqn(dtastk(isp1+3),3))then
        go to 9000
      endif
      kavc=ktfaddr(klist(kac+1))
      ka=ktfaddr(ktastk(isp1+1))
      if(.not. tflistq(ktastk(isp1+1)))then
        go to 9000
      endif
      if(ilist(2,ka-1) /= 3)then
        go to 9000
      endif
      if(ktfreallistq(ka))then
        go to 9000
      endif
      m=-1
      do i=1,3
        if(.not. tfreallistq(klist(ka+i)))then
          go to 9000
        endif
        kai=ilist(2,ka+i)
        if(m < 0)then
          m=ilist(2,kai-1)
        elseif(m /= ilist(2,kai-1))then
          go to 9000
        endif
        iav(i)=ilist(2,kai+1)
      enddo
      isp0=isp
      m3=m/3
      call tfcliptriangle1(m3,
     $     rlist(iav(1)+1:iav(1)+m3),
     $     rlist(iav(2)+1:iav(2)+m3),
     $     rlist(iav(3)+1:iav(3)+m3))
      nt=(isp-isp0)/9
      if(nt <= 0)then
        kx=dxnulll
      else
        xs=rlist(kavs+1)
        ys=rlist(kavs+2)
        zs=rlist(kavs+3)
        xc=rlist(kavc+1)
        yc=rlist(kavc+2)
        zc=rlist(kavc+3)
        kx=kxadaloc(-1,nt,klx)
        do i=1,nt
          ii=isp0+(i-1)*9
          klx%dbody(i)=kxadaloc(0,3,klxi)
          do j=1,3
            klxi%dbody(j)=kxavaloc(0,3,klxj)
            klxj%rbody(2)=vstk(ivstkoffset+ii+j)*xs+xc
            klxj%rbody(3)=vstk(ivstkoffset+ii+j+3)*ys+yc
            klxj%rbody(4)=vstk(ivstkoffset+ii+j+6)*zs+zc
          enddo
        enddo
      endif
      isp=isp0
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $  '"{{x, ..},{y, ..},{z, ..}}, scale, center"')
      return
      end

      subroutine tfcliptriangle1(m,x,y,z)
      use tfstk
      implicit none
      integer*4 isp0,m,i,isp1
      real*8 x(3,m),y(3,m),z(3,m)
      isp0=isp
      do i=1,m
        rtastk(isp+1)=x(1,i)
        rtastk(isp+2)=x(2,i)
        rtastk(isp+3)=x(3,i)
        rtastk(isp+4)=y(1,i)
        rtastk(isp+5)=y(2,i)
        rtastk(isp+6)=y(3,i)
        rtastk(isp+7)=z(1,i)
        rtastk(isp+8)=z(2,i)
        rtastk(isp+9)=z(3,i)
        isp=isp+9
        if(
     $       x(1,i) > 1.d0 .or. x(1,i) < -1.d0
     $       .or. x(2,i) > 1.d0 .or. x(2,i) < -1.d0
     $       .or. x(3,i) > 1.d0 .or. x(3,i) < -1.d0
     $       .or. y(1,i) > 1.d0 .or. y(1,i) < -1.d0
     $       .or. y(2,i) > 1.d0 .or. y(2,i) < -1.d0
     $       .or. y(3,i) > 1.d0 .or. y(3,i) < -1.d0
     $       .or. z(1,i) > 1.d0 .or. z(1,i) < -1.d0
     $       .or. z(2,i) > 1.d0 .or. z(2,i) < -1.d0
     $       .or. z(3,i) > 1.d0 .or. z(3,i) < -1.d0)then
          isp1=isp
          call tfcliptriangle2(isp1,1,.true.)
          call tfcliptriangle2(isp1,1,.true.)
          call tfcliptriangle2(isp1,1,.true.)
        endif
      enddo
      return
      end

      recursive subroutine tfcliptriangle2(isp1,idir,chg)
      use tfstk
      implicit none
      integer*4 isp1,idir,i,j
      real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3,r1,r2,r3,x,y,z
      logical*4 chg
      if(isp1 >= isp)then
        return
      endif
      do i=isp1+1,isp,9
 1      x1=rtastk(i)
        x2=rtastk(i+1)
        x3=rtastk(i+2)
        y1=rtastk(i+3)
        y2=rtastk(i+4)
        y3=rtastk(i+5)
        z1=rtastk(i+6)
        z2=rtastk(i+7)
        z3=rtastk(i+8)
        if(x1 > 1.d0)then
          if(x2 > 1.d0)then
            if(x3 > 1.d0)then
              do j=i+10,isp
                rtastk(j-9)=rtastk(j)
              enddo
              isp=isp-9
              call tfcliptriangle2(i-1,idir,.false.)
              if(chg)then
                go to 100
              endif
              return
            else
              r1=(1.d0-x3)/(x1-x3)
              r2=(1.d0-x3)/(x2-x3)
              rtastk(i)=1.d0
              rtastk(i+1)=1.d0
              rtastk(i+3)=r1*(y1-y3)+y3
              rtastk(i+4)=r2*(y2-y3)+y3
              rtastk(i+6)=r1*(z1-z3)+z3
              rtastk(i+7)=r2*(z2-z3)+z3
            endif
          else
            if(x3 > 1.d0)then
              r1=(1.d0-x2)/(x1-x2)
              r3=(1.d0-x2)/(x3-x2)
              rtastk(i)=1.d0
              rtastk(i+2)=1.d0
              rtastk(i+3)=r1*(y1-y2)+y2
              rtastk(i+5)=r3*(y3-y2)+y2
              rtastk(i+6)=r1*(z1-z2)+z2
              rtastk(i+8)=r3*(z3-z2)+z2
            else
              r2=(1.d0-x2)/(x1-x2)
              r3=(1.d0-x3)/(x1-x3)
              rtastk(isp+3)=1.d0
              rtastk(isp+6)=r3*(y1-y3)+y3
              rtastk(isp+9)=r3*(z1-z3)+z3
              rtastk(i)=1.d0
              rtastk(i+3)=r2*(y1-y2)+y2
              rtastk(i+6)=r2*(z1-z2)+z2
              rtastk(isp+1)=1.d0
              rtastk(isp+4)=rtastk(i+3)
              rtastk(isp+7)=rtastk(i+6)
              rtastk(isp+2)=x3
              rtastk(isp+5)=y3
              rtastk(isp+8)=z3
              isp=isp+9
            endif
          endif
        else
          if(x2 > 1.d0)then
            rtastk(i)=x2
            rtastk(i+1)=x3
            rtastk(i+2)=x1
            rtastk(i+3)=y2
            rtastk(i+4)=y3
            rtastk(i+5)=y1
            rtastk(i+6)=z2
            rtastk(i+7)=z3
            rtastk(i+8)=z1
            go to 1
          elseif(x3 > 1.d0)then
            rtastk(i)=x3
            rtastk(i+1)=x1
            rtastk(i+2)=x2
            rtastk(i+3)=y3
            rtastk(i+4)=y1
            rtastk(i+5)=y2
            rtastk(i+6)=z3
            rtastk(i+7)=z1
            rtastk(i+8)=z2
            go to 1
          endif
        endif
      enddo
 100  if(isp1 < isp)then
        do i=isp1+1,isp,9
          rtastk(i)=-rtastk(i)
          rtastk(i+1)=-rtastk(i+1)
          rtastk(i+2)=-rtastk(i+2)
        enddo
        if(idir > 0)then
          call tfcliptriangle2(isp1,-1,.true.)
        else
          do i=isp1+1,isp,9
            x=rtastk(i)
            rtastk(i)=rtastk(i+3)
            rtastk(i+3)=rtastk(i+6)
            rtastk(i+6)=x
            y=rtastk(i+1)
            rtastk(i+1)=rtastk(i+4)
            rtastk(i+4)=rtastk(i+7)
            rtastk(i+7)=y
            z=rtastk(i+2)
            rtastk(i+2)=rtastk(i+5)
            rtastk(i+5)=rtastk(i+8)
            rtastk(i+8)=z
          enddo
        endif
      endif
      return
      end

      subroutine tfcanvas3dlighttriangle(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kat,kac,kaci,kaci1,kaci2
      integer*4 isp1,irtc,itfmessage,nt,nc,isp0,i
      if(isp1+2 /= isp)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      kat=ktfaddr(ktastk(isp1+1))
      if(.not. tflistq(ktastk(isp1+1)))then
        go to 9000
      endif
      nt=ilist(2,kat-1)
      if(nt <= 0)then
        kx=kxnulll
        irtc=0
        return
      endif
      if(ktfreallistq(kat))then
        go to 9000
      endif
      kac=ktfaddr(ktastk(isp))
      if(.not. tflistq(ktastk(isp1+1)))then
        go to 9100
      endif
      if(ktfreallistq(kac))then
        go to 9100
      endif
      nc=ilist(2,kac-1)
      if(nc <= 0)then
        go to 9100
      endif
      isp0=isp
      do i=1,nc
        if(.not. tflistq(klist(kac+i)))then
          go to 9110
        endif
        kaci=ktfaddr(klist(kac+i))
        if(ktfreallistq(kaci) .or. ilist(2,kaci-1) /= 2)then
          go to 9110
        endif
        if(.not. tfnumlistqn(dlist(kaci+1),3) .or.
     $       .not. tfnumlistqn(dlist(kaci+2),3))then
          go to 9110
        endif
        kaci1=ktfaddr(klist(kaci+1))
        kaci2=ktfaddr(klist(kaci+2))
        rtastk(isp+1)=rlist(kaci1+1)
        rtastk(isp+2)=rlist(kaci1+2)
        rtastk(isp+3)=rlist(kaci1+3)
        rtastk(isp+4)=rlist(kaci2+1)
        rtastk(isp+5)=rlist(kaci2+2)
        rtastk(isp+6)=rlist(kaci2+3)
        isp=isp+6
      enddo
      call tfcanvas3dlighttriangle1(nt,kat,nc,
     $     rtastk(isp0+1:isp),kx,irtc)
      isp=isp0
      if(irtc /= 0)then
        return
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"List of Triangles"')
      return
 9110 isp=isp0
 9100 irtc=itfmessage(9,'General::wrongtype',
     $ '"List of light sources: {{{x1, y1, z1}, {r1, g1, b1}}, ...}"')
      return
      end

      subroutine tfcanvas3dlighttriangle1(nt,kat,nc,slight,kx,irtc)
      use tfstk
      implicit none
      integer*8 kat,kx,kaci,kati,kati1,kati2,kati3
      integer*4 nt,nc,irtc,i,j,isp0,itfmessage
      real*8 slight(3,2,nc),rgb(3,nt),u0(nc),
     $     x1,x2,x3,y1,y2,y3,cx,cy,cz,anx,any,anz,an,vx,vy,vz,u,
     $     z1,z2,z3
      do j=1,nc
        u0(j)=slight(1,1,j)**2+slight(2,1,j)**2+slight(3,1,j)**2
      enddo
      do9000: do
        do i=1,nt
          rgb(1,i)=0.d0
          rgb(2,i)=0.d0
          rgb(3,i)=0.d0
          if(.not. tflistq(klist(kat+i)))then
            exit do9000
          endif
          kati=ktfaddr(klist(kat+i))
          if(ilist(2,kati-1) /= 3 .or. ktfreallistq(kati))then
            exit do9000
          endif
          if(.not. tfnumlistqn(dlist(kati+1),3)
     $         .or. .not. tfnumlistqn(dlist(kati+2),3)
     $         .or. .not. tfnumlistqn(dlist(kati+3),3))then
            exit do9000
          endif
          kati1=ktfaddr(klist(kati+1))
          kati2=ktfaddr(klist(kati+2))
          kati3=ktfaddr(klist(kati+3))
          x1=rlist(kati1+1)
          y1=rlist(kati1+2)
          z1=rlist(kati1+3)
          x2=rlist(kati2+1)
          y2=rlist(kati2+2)
          z2=rlist(kati2+3)
          x3=rlist(kati3+1)
          y3=rlist(kati3+2)
          z3=rlist(kati3+3)
          cx=(x1+x2+x3)/3.d0
          cy=(y1+y2+y3)/3.d0
          cz=(z1+z2+z3)/3.d0
          anx=(y2-y1)*(z3-z2)-(y3-y2)*(z2-z1)
          any=(z2-z1)*(x3-x2)-(z3-z2)*(x2-x1)
          anz=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1)
          an=sqrt(anx**2+any**2+anz**2)
          anx=anx/an
          any=any/an
          anz=anz/an
          do j=1,nc
            vx=slight(1,1,j)-cx
            vy=slight(2,1,j)-cy
            vz=slight(3,1,j)-cz
            u=(vx*anx+vy*any+vz*anz)/sqrt(vx**2+vy**2+vz**2)**3*u0(j)
            if(u > 0.d0)then
              rgb(1,i)=rgb(1,i)+slight(1,2,j)*u
              rgb(2,i)=rgb(2,i)+slight(2,2,j)*u
              rgb(3,i)=rgb(3,i)+slight(3,2,j)*u
            endif
          enddo
          rgb(1,i)=min(1.d0,max(0.d0,rgb(1,i)))
          rgb(2,i)=min(1.d0,max(0.d0,rgb(2,i)))
          rgb(3,i)=min(1.d0,max(0.d0,rgb(3,i)))
        enddo
        isp0=isp
        do i=1,nt
          kaci=ktavaloc(-1,3)
          rlist(kaci+2)=rgb(1,i)
          rlist(kaci+3)=rgb(2,i)
          rlist(kaci+4)=rgb(3,i)
          isp=isp+1
          ktastk(isp)=ktflist+kaci
        enddo
        kx=ktflist+ktfmakelist(isp0)
        isp=isp0
        irtc=0
        return
      enddo do9000
      irtc=itfmessage(9,'General::wrongtype',
     $     '"List of Triangles"')
      return
      end

      subroutine tfcanvas3dprojection(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx,kae,kapx,kapy,kaoff,kave,kavx,kavy,kavoff
      integer*4 isp1,irtc,itfmessage
      real*8 d,e(3)
      if(isp /= isp1+5)then
        irtc=itfmessage(9,'General::narg','"5"')
        return
      endif
      if(.not. tfnumlistqn(dtastk(isp1+2),3))then
        go to 9000
      endif
      kae=ktfaddr(ktastk(isp1+2))
      if(.not. tfnumlistqn(dtastk(isp1+3),3))then
        go to 9000
      endif
      kapx=ktfaddr(ktastk(isp1+3))
      if(.not. tfnumlistqn(dtastk(isp1+4),3))then
        go to 9000
      endif
      kapy=ktfaddr(ktastk(isp1+4))
      if(.not. tfnumlistqn(dtastk(isp),2))then
        go to 9000
      endif
      kaoff=ktfaddr(ktastk(isp))
      kave=ktfaddr(klist(kae+1))
      kavx=ktfaddr(klist(kapx+1))
      kavy=ktfaddr(klist(kapy+1))
      kavoff=ktfaddr(klist(kaoff+1))
      d=sqrt(rlist(kave+1)**2+rlist(kave+2)**2+rlist(kave+3)**2)
      e(1)=rlist(kave+1)/d
      e(2)=rlist(kave+2)/d
      e(3)=rlist(kave+3)/d
      call tfcanvas3dprojection1(ktastk(isp1+1),e,
     $     rlist(kavx+1:kavx+3),rlist(kavy+1:kavy+3),
     $     d,rlist(kavoff+1:kavoff+2),kx,irtc)
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"[data, Eye, Xproj, Yproj, Offset]"')
      return
      end

      recursive subroutine
     $     tfcanvas3dprojection1(k,e,px,py,d,offset,kx,irtc)
      use tfstk
      implicit none
      integer*8 k,kx,kxi,ki,kax,ka1,kai,ka,kac,kaxi
      integer*4 irtc,itfmessage,n,isp0,i
      real*8 e(3),px(3),py(3),d,x,y,z,p,a,ed,ux,uy,uz,xp,yp,
     $     offset(2)
      if(.not. tflistq(k))then
        go to 9000
      endif
      ka=ktfaddr(k)
      if(ktfreallistq(ka))then
        go to 9000
      endif
      n=ilist(2,ka-1)
      ka1=ktfaddr(klist(ka+1))
      isp0=isp
      if(tfreallistq(klist(ka+1)))then
        ed=1.d300
        do i=1,n
          kai=ktfaddr(klist(ka+i))
          if(.not. tfnumlistqn(dlist(ka+i),3))then
            go to 9100
          endif
          x=rlist(kai+1)
          y=rlist(kai+2)
          z=rlist(kai+3)
          p=x*e(1)+y*e(2)+z*e(3)
          a=1.d0-p/d
          ux=(x-p*e(1))/a
          uy=(y-p*e(2))/a
          uz=(z-p*e(3))/a
          xp=ux*px(1)+uy*px(2)+uz*px(3)
          yp=ux*py(1)+uy*py(2)+uz*py(3)
          kaxi=ktavaloc(-1,2)
          rlist(kaxi+2)=anint(xp+offset(1))
          rlist(kaxi+3)=anint(yp+offset(2))
          isp=isp+1
          ktastk(isp)=ktflist+kaxi
          ed=min(ed,-((x-e(1))**2+(y-e(2))**2+(z-e(3))**2))
        enddo
        kac=ktfmakelist(isp0)
        isp=isp0
        kax=ktadaloc(-1,2)
        rlist(kax+1)=ed
        klist(kax+2)=ktflist+ktfcopy1(kac)
      else
        do i=1,n
          ki=klist(ka+i)
          call tfcanvas3dprojection1(
     $         ki,e,px,py,d,offset,kxi,irtc)
          if(irtc /= 0)then
            isp=isp0
            return
          endif
          isp=isp+1
          ktastk(isp)=kxi
        enddo
        kax=ktfmakelist(isp0)
        isp=isp0
      endif
      kx=ktflist+kax
      irtc=0
      return
 9100 isp=isp0
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"{{x1, y1, z1}, ...} or List of them"')
      return
      end

      subroutine tfcanvassymbol(isp1,kx,irtc)
      use tfstk
      use convstr,only:tfconvround
      use strbuf
      use macmath
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_strbuf), pointer :: strb
      type (sad_string), pointer :: str
      integer*8 ka,kad,kad1,kad2,ks
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i,np
      real*8 x,y,s,a,xmin,ymin,xmax,ymax,yoff,sa,sh,ar(32),s1
c      parameter (a=sqrt(0.75d0))
      parameter (a=.866025403784439d0)
      character*2 sym
      logical*4 full
      if(isp /= isp1+8)then
        irtc=itfmessage(9,'General::narg','"10"')
        return
      endif
      if(iand(ktfmask,ktastk(isp1+8)) /= ktflist)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Canvas$Range for #8"')
        return
      endif
      if(iand(ktfmask,ktastk(isp1+9)) /= ktflist)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Canvas$Offset for #9"')
        return
      endif
      kad=ktfaddr(ktastk(isp1+8))
      yoff=rlist(ktfaddr(ktastk(isp1+9))+2)
      call getstringbuf(strb,0,.true.)
      sym="1O"
      if(ktfstringq(ktastk(isp1+4)))then
        ks=ktfaddr(ktastk(isp1+4))
        if(ilist(1,ks) == 2)then
          call tmovb(ilist(1,ks+1),sym,2)
        endif
      endif
      if(ktfnonlistq(ktastk(isp1+3)))then
        go to 9000
      endif
      ka=ktfaddr(ktastk(isp1+3))
      if(klist(ka) /= ktfoper+mtflist
     $     .or. ilist(2,ka-1) /= 2 .or.
     $     ktfreallistq(ka))then
        go to 9000
      endif
      x=rlist(ka+1)
      y=rlist(ka+2)
      kad1=ktfaddr(klist(kad+1))
      kad2=ktfaddr(klist(kad+2))
      xmin=rlist(kad1+1)
      ymin=rlist(kad1+2)
      xmax=rlist(kad2+1)
      ymax=rlist(kad2+2)
      if(x < xmin .or. x > xmax .or.
     $     y < ymin .or. y > ymax)then
        kx%k=ktfoper+mtfnull
        call tfreestringbuf(strb)
        return
      endif
      if(ktfnonrealq(ktastk(isp1+5)))then
      irtc=itfmessage(9,'General::wrongtype','"size for #5"')
      return
      endif
      s=rtastk(isp1+5)
      if(ktfnonstringq(ktastk(isp1+6)))then
        irtc=itfmessage(9,'General::wrongtype','"fill color for #6"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+7)))then
        irtc=itfmessage(9,'General::wrongtype','"outline for #7"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+1)))then
        irtc=itfmessage(9,'General::wrongtype','"pathname for #1"')
        return
      endif
      call loc_string(ktfaddr(ktastk(isp1+1)),str)
      call putstringbufb(strb,str%str,str%nch,full)
      np=0
      if(sym == "1O")then
        call putstringbufb(strb,' create oval ',13,full)
        call tfconvround(strb,x-s)
        call tfconvround(strb,y-s)
        call tfconvround(strb,x+s)
        call tfconvround(strb,y+s)
        call putstringbufb(strb,' -outline {',11,full)
        call loc_string(ktfaddr(ktastk(isp1+7)),str)
        call putstringbufb(strb,str%str,str%nch,full)
        call putstringbufb1(strb,'}')
      else
        if(sym(2:2) == "O")then
          call putstringbufb(strb,' create polygon ',16,full)
          s1=s*sqrt(m_pi/a/1.5d0)
          sa=s1*a
          sh=s1*.5d0
          if(sym(1:1) == "6")then
            ar(1:6)=(/x-s1,y,x+sh,y-sa,x+sh,y+sa/)
          elseif(sym(1:1) == "7")then
            ar(1:6)=(/x+s1,y,x-sh,y-sa,x-sh,y+sa/)
          elseif(sym(1:1) == "8")then
            ar(1:6)=(/x,y+s1,x-sa,y-sh,x+sa,y-sh/)
          elseif(sym(1:1) == "9")then
            ar(1:6)=(/x,y-s1,x-sa,y+sh,x+sa,y+sh/)
          endif
          np=6
        elseif(sym == "BX")then
          s1=s*sqrt(m_pi_4)
          ar(1:8)=(/x+s1,y-s1,x+s1,y+s1,x-s1,y+s1,x-s1,y-s1/)
          np=8
        elseif(sym == "RH")then
          s1=s*sqrt(m_pi_2)
          ar(1:8)=(/x+s1,y,x,y+s1,x-s1,y,x,y-s1/)
          np=8
        elseif(sym == "PL")then
          ar(1:24)=(/x+s,y-1,x+s,y+1,x+1,y+1,x+1,y+s,
     $         x-1,y+s,x-1,y+1,x-s,y+1,x-s,y-1,x-1,y-1,
     $         x-1,y-s,x+1,y-s,x+1,y-1/)
          np=24
        elseif(sym == "TI")then
          ar(1:24)=(/x+s+1,y+s-1,x+s-1,y+s+1,
     $         x,y+1,x-s+1,y+s+1,x-s-1,y+s-1,x-1,y,
     $         x-s-1,y-s+1,x-s+1,y-s-1,x,y-1,x+s-1,y-s-1,
     $         x+s+1,y-s+1,x+1,y/)
          np=24
        endif
        do i=1,np
          call tfconvround(strb,ar(i))
        enddo
      endif
      call putstringbufb(strb,' -fill {',8,full)
      call loc_string(ktfaddr(ktastk(isp1+6)),str)
      call putstringbufb(strb,str%str,str%nch,full)
      call putstringbufb1(strb,'}')
      if(ktfstringq(ktastk(isp1+2)))then
        call putstringbufb(strb,' -tags ',7,full)
        call loc_string(ktfaddr(ktastk(isp1+2)),str)
        call putstringbufb(strb,str%str,str%nch,full)
      endif
      kx=kxstringbuftostring(strb)
      irtc=0
      return
 9000 call tfreestringbuf(strb)
      irtc=itfmessage(9,'General::wrongtype','"{x, y} for #2"')
      return
      end

      subroutine tfcanvassymboldirect(isp1,kx,irtc)
      use tfstk
      use convstr,only:ktrsaloc
      use macmath
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 ks,ka,kad,kad1,kad2,kat,kavx,kavy,ka2
      integer*4 itfmessage,isp0,m,i,isym
      real*8 x,y,s,a,xmin,xmax,ymin,ymax,yoff,wmin,sa,sh,s1
      logical*4 ol
      type (sad_descriptor), save ::
     $     ioutline,ifill,itag,iwidth,ibar
      data ioutline%k,ifill%k,itag%k,iwidth%k,ibar%k /0,0,0,0,0/
c      parameter (a=sqrt(0.75d0))
      parameter (a=.866025403784439d0,wmin=7.d0)
      character*2 sym
      if(ioutline%k == 0)then
        ioutline=kxsalocb(0,'-outline',8)
        ifill=kxsalocb(0,'-fill',5)
        itag=kxsalocb(0,'-tags',5)
        iwidth=kxsalocb(0,'-width',6)
        ibar=kxsymbolf('Bar',3,.true.)
      endif
      if(isp /= isp1+13)then
        irtc=itfmessage(9,'General::narg','"13"')
        return
      endif
      ka2=0
      sym="1O"
      if(ktfstringq(ktastk(isp1+4)))then
        ks=ktfaddr(ktastk(isp1+4))
        if(ilist(1,ks) == 2)then
          call tmovb(ilist(1,ks+1),sym,2)
        endif
      elseif(ktfsymbolq(ktastk(isp1+4)))then
        if(tfsamesymbolq(ibar,dtastk(isp1+4)))then
          sym="BA"
        endif
      endif
      if(ktfnonlistq(ktastk(isp1+8)))then
        irtc=itfmessage(9,'General::wrongtype','"Canvas$Range for #8"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+9)))then
        irtc=itfmessage(9,'General::wrongtype','"Canvas$Offset for #9"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+3)))then
        go to 9000
      endif
      ka=ktfaddr(ktastk(isp1+3))
      if(klist(ka) /= ktfoper+mtflist .or. ilist(2,ka-1) < 2 .or. ilist(2,ka-1) > 4)then
        go to 9000
      endif
      if(ktfreallistq(ka))then
        m=1
        kavx=ka
        kavy=kavx+1
      elseif(ktfnonlistq(klist(ka+1)) .or. ktfnonlistq(klist(ka+2)))then
        go to 9000
      else
        kavx=ktfaddr(klist(ka+1))
        kavy=ktfaddr(klist(ka+2))
        if(klist(kavx) /= ktfoper+mtflist .or. klist(kavy) /= ktfoper+mtflist .or.
     $       ktfnonreallistq(kavx) .or. ktfnonreallistq(kavy))then
          go to 9000
        else
          m=ilist(2,kavx-1)
          if(m /= ilist(2,kavy-1))then
            go to 9000
          endif
        endif
      endif
      kad=ktfaddr(ktastk(isp1+8))
      kad1=ktfaddr(klist(kad+1))
      kad2=ktfaddr(klist(kad+2))
      xmin=rlist(kad1+1)
      ymin=rlist(kad1+2)
      xmax=rlist(kad2+1)
      ymax=rlist(kad2+2)
      yoff=rlist(ktfaddr(ktastk(isp1+9)+2))
      if(ktfnonrealq(ktastk(isp1+5),s))then
        irtc=itfmessage(9,'General::wrongtype','"size for #5"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+6)))then
        irtc=itfmessage(9,'General::wrongtype','"fill color for #6"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+7)))then
        irtc=itfmessage(9,'General::wrongtype','"outline for #7"')
        return
      endif
c      if(itastk(1,isp1+1) /= ntfreal)then
c        irtc=itfmessage(9,'General::wrongtype',
c     $       '"TkCanvasPointer for #1"')
c        return
c      endif
      if(ktfstringq(ktastk(isp1+2)))then
        kat=ktfaddr(ktastk(isp1+2))
      elseif(ktflistq(ktastk(isp1+2)))then
        kat=ktfaddr(ktastk(isp1+2))
        if(ilist(2,kat-1) /= m)then
          irtc=itfmessage(9,'General::equalleng',
     $         '"tags (#2) and points (#3)"')
          return
        endif
        if(ktfreallistq(kat))then
          kat=0
        endif
      else
        kat=0
      endif
      ol=.true.
      isp0=isp
      if((sym == "BA" .or. sym == "BR") .and. s < wmin)then
        ka2=ktrsaloc(-1,anint(max(1.d0,s/wmin*2)))
        ol=.false.
      endif
      isp=isp0
      isym=1
      if(sym == "1O" .or. sym == "CI" .or. sym == "CL")then
        isym=1
      elseif(sym == "6O" .or. sym == "LT")then
        isym=2
      elseif(sym == "7O" .or. sym == "RT")then
        isym=3
      elseif(sym == "8O" .or. sym == "DT")then
        isym=4
      elseif(sym == "9O" .or. sym == "UT")then
        isym=5
      elseif(sym == "BX" .or. sym == "SQ")then
        isym=6
      elseif(sym == "RH" .or.  sym == "DM" .or. sym == "DI")then
        isym=7
      elseif(sym == "PL")then
        isym=8
      elseif(sym == "TI" .or. sym == "ML" .or. sym == "MU")then
        isym=9
      elseif(sym == "BA" .or. sym == "BR")then
        isym=99
      endif
      do i=1,m
        x=rlist(kavx+i)
        y=rlist(kavy+i)
        if(x < xmin .or. x > xmax)then
          cycle
        elseif(isym /= 99)then
          if(y < ymin .or. y > ymax)then
            cycle
          endif
        else
          if(max(y,yoff) < ymin .or. min(y,yoff) > ymax)then
            cycle
          endif
        endif
        isp=isp+1
        rtastk(isp)=rtastk(isp1+1)
        rtastk(isp)=rtastk(isp1+1)
        isp=isp+1
        select case (isym)
        case (1)
          rtastk(isp)=rtastk(isp1+10)
          rtastk(isp+1:isp+4)=(/x-s,y-s,x+s,y+s/)
          isp=isp+5
        case (2,3,4,5)
          s1=s*sqrt(m_pi/a/1.5d0)
          sa=s1*a
          sh=s1*.5d0
          rtastk(isp)=rtastk(isp1+11)
          select case (isym)
          case (2)
            rtastk(isp+1:isp+6)=(/x-s1,y,x+sh,y-sa,x+sh,y+sa/)
          case (3)
            rtastk(isp+1:isp+6)=(/x+s1,y,x-sh,y-sa,x-sh,y+sa/)
          case (4)
            rtastk(isp+1:isp+6)=(/x,y+s1,x-sa,y-sh,x+sa,y-sh/)
          case (5)
            rtastk(isp+1:isp+6)=(/x,y-s1,x-sa,y+sh,x+sa,y+sh/)
          end select
          isp=isp+7
        case (6)
          s1=s*sqrt(m_pi_4)
          rtastk(isp)=rtastk(isp1+11)
          rtastk(isp+1:isp+8)=
     $         (/x+s1,y-s1,x+s1,y+s1,x-s1,y+s1,x-s1,y-s1/)
          isp=isp+9
        case (7)
          s1=s*sqrt(m_pi_2)
          rtastk(isp)=rtastk(isp1+11)
          rtastk(isp+1:isp+8)=
     $         (/x+s1,y,x,y+s1,x-s1,y,x,y-s1/)
          isp=isp+9
        case (8)
          s1=s*1.2d0
          rtastk(isp)=rtastk(isp1+11)
          rtastk(isp+1:isp+24)=(/x+s1,y-1,x+s1,y+1,x+1,y+1,x+1,y+s1,
     $         x-1,y+s1,x-1,y+1,x-s1,y+1,x-s1,y-1,x-1,y-1,
     $         x-1,y-s1,x+1,y-s1,x+1,y-1/)
          isp=isp+25
        case(9)
          rtastk(isp)=rtastk(isp1+11)
          rtastk(isp+1:isp+24)=(/x+s+1,y+s-1,x+s-1,y+s+1,
     $         x,y+1,x-s+1,y+s+1,x-s-1,y+s-1,x-1,y,
     $         x-s-1,y-s+1,x-s+1,y-s-1,x,y-1,x+s-1,y-s-1,
     $         x+s+1,y-s+1,x+1,y/)
          isp=isp+25
        case (99)
          if(ol)then
            rtastk(isp)=rtastk(isp1+12)
            rtastk(isp+1)=x-s/wmin
            rtastk(isp+2)=min(max(ymin,y),ymax)
            rtastk(isp+3)=x+s/wmin
            rtastk(isp+4)=min(max(ymin,yoff),ymax)
          else
            rtastk(isp)=rtastk(isp1+13)
            rtastk(isp+1)=x
            rtastk(isp+2)=min(max(ymin,y),ymax)
            rtastk(isp+3)=x
            rtastk(isp+4)=min(max(ymin,yoff),ymax)
          endif
          isp=isp+5
        end select
        if(ol)then
          dtastk(isp)=ioutline
          isp=isp+1
          rtastk(isp)=rtastk(isp1+7)
        else
          dtastk(isp)=iwidth
          isp=isp+1
          ktastk(isp)=ktfstring+ka2
        endif
        isp=isp+1
        dtastk(isp)=ifill
        isp=isp+1
        rtastk(isp)=rtastk(isp1+6)
        if(kat /= 0)then
          if(ktfstringq(ktastk(isp1+2)))then
            isp=isp+1
            dtastk(isp)=itag
            isp=isp+1
            rtastk(isp)=rtastk(isp1+2)
          elseif(ktfstringq(klist(kat+i)))then
            isp=isp+1
            dtastk(isp)=itag
            isp=isp+1
            ktastk(isp)=klist(kat+i)
          endif
        endif
        call tftkcanvcreateitem(isp0,kx,irtc)
        isp=isp0
      enddo
      kx%k=ktfoper+mtfnull
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"{x, y} or {{x1,..,xn}, {y1,..,yn}} for #3"')
      return
      end
