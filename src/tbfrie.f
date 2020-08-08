      subroutine tbfrie(trans,cod,beam,rhob,ak,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      implicit none
      integer*4 i
      real*8 trans(6,12),cod(6),beam(42),trans1(6,6),
     $     yi,pr,rho,ak,rhob,pxi
      logical*4 ent
      yi=cod(3)
      pr=1.d0+cod(6)
      rho=rhob*pr
      if(ent)then
        trans1(1,1)=1.d0
        trans1(1,2)=0.d0
        trans1(1,3)=yi/rho
        trans1(1,4)=0.d0
        trans1(1,5)=0.d0
        trans1(1,6)=-(yi/pr)**2/rhob*.5d0
        trans1(2,1)=-ak
        trans1(2,2)=1.d0
        trans1(2,3)=-ak*trans1(1,3)
        trans1(2,4)=0.d0
        trans1(2,5)=0.d0
        trans1(2,6)=-ak*trans1(1,6)
        trans1(3,1)=0.d0
        trans1(3,2)=0.d0
        trans1(3,3)=1.d0
        trans1(3,4)=0.d0
        trans1(3,5)=0.d0
        trans1(3,6)=0.d0
        trans1(4,1)=0.d0
        trans1(4,2)=-trans1(1,3)
        trans1(4,3)=ak-cod(2)/rho
        trans1(4,4)=1.d0
        trans1(4,5)=0.d0
        trans1(4,6)=yi*cod(2)/rho/pr
        trans1(5,1)=0.d0
        trans1(5,2)=trans1(1,6)
        trans1(5,3)=-trans1(4,6)
        trans1(5,4)=0.d0
        trans1(5,5)=1.d0
        trans1(5,6)=(yi/pr)**2*cod(2)/rho
        trans1(6,1)=0.d0
        trans1(6,2)=0.d0
        trans1(6,3)=0.d0
        trans1(6,4)=0.d0
        trans1(6,5)=0.d0
        trans1(6,6)=1.d0
        do i=1,irad
          trans(4,i)=trans1(4,2)*trans(2,i)+trans1(4,3)*trans(3,i)
     $         +trans(4,i)+trans1(4,6)*trans(6,i)
          trans(5,i)=trans1(5,2)*trans(2,i)+trans1(5,3)*trans(3,i)
     $         +trans(5,i)+trans1(5,6)*trans(6,i)
          trans(2,i)=trans1(2,1)*trans(1,i)+trans(2,i)
     $         +trans1(2,3)*trans(3,i)+trans1(2,6)*trans(6,i)
          trans(1,i)=trans(1,i)+trans1(1,3)*trans(3,i)
     $         +trans1(1,6)*trans(6,i)
        enddo
c        write(*,'(a/,6(1p6g12.5/))')'tbfrie-1 ',trans1
        call tmulbs(beam,trans1,.true.)
        cod(1)=cod(1)+.5d0*yi**2/rho
        cod(4)=min(pr,max(-pr,cod(4)+(ak-cod(2)/rho)*yi))
        cod(5)=cod(5)-.5d0*(yi/pr)**2*cod(2)/rhob
        cod(2)=min(pr,max(-pr,cod(2)-ak*cod(1)))
      else
        pxi=min(pr,max(-pr,cod(2)-ak*cod(1)))
        trans1(1,1)=1.d0
        trans1(1,2)=0.d0
        trans1(1,3)=yi/rho
        trans1(1,4)=0.d0
        trans1(1,5)=0.d0
        trans1(1,6)=-(yi/pr)**2/rhob*.5d0
        trans1(2,1)=-ak
        trans1(2,2)=1.d0
        trans1(2,3)=0.d0
        trans1(2,4)=0.d0
        trans1(2,5)=0.d0
        trans1(2,6)=0.d0
        trans1(3,1)=0.d0
        trans1(3,2)=0.d0
        trans1(3,3)=1.d0
        trans1(3,4)=0.d0
        trans1(3,5)=0.d0
        trans1(3,6)=0.d0
        trans1(4,2)=-trans1(1,3)
        trans1(4,1)=-ak*trans1(4,2)
        trans1(4,3)=ak-pxi/rho
        trans1(4,4)=1.d0
        trans1(4,5)=0.d0
        trans1(4,6)=yi*pxi/rho/pr
        trans1(5,2)=trans1(1,6)
        trans1(5,1)=-ak*trans1(5,2)
        trans1(5,3)=-trans1(4,6)
        trans1(5,4)=0.d0
        trans1(5,5)=1.d0
        trans1(5,6)=(yi/pr)**2*pxi/rho
        trans1(6,1)=0.d0
        trans1(6,2)=0.d0
        trans1(6,3)=0.d0
        trans1(6,4)=0.d0
        trans1(6,5)=0.d0
        trans1(6,6)=1.d0
        do i=1,irad
          trans(4,i)=trans1(4,1)*trans(1,i)+trans1(4,2)*trans(2,i)
     $         +trans1(4,3)*trans(3,i)+trans(4,i)+trans1(4,6)*trans(6,i)
          trans(5,i)=trans1(5,1)*trans(1,i)+trans1(5,2)*trans(2,i)
     $         +trans1(5,3)*trans(3,i)+trans(5,i)+trans1(5,6)*trans(6,i)
          trans(2,i)=trans1(2,1)*trans(1,i)+trans(2,i)
          trans(1,i)=trans(1,i)+trans1(1,3)*trans(3,i)
     $         +trans1(1,6)*trans(6,i)
        enddo
c        write(*,'(a/,6(1p6g12.5/))')'tbfrie-2 ',trans1
        call tmulbs(beam,trans1,.true.)
        cod(2)=pxi
        cod(1)=cod(1)+.5d0*yi**2/rho
        cod(4)=min(pr,max(-pr,cod(4)+(ak-pxi/rho)*yi))
        cod(5)=cod(5)-.5d0*(yi/pr)**2*pxi/rhob
      endif
      return
      end

      subroutine tbedge(trans0,cod,beam,al,phib,psi,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      implicit none
      real*8 trans0(6,12),cod(6),beam(42),phib,psi,al,
     $     trans(6,6),trans1(6,6),rhob
      logical*4 ent
      rhob=al/phib
      if(ent)then
        if(psi .eq. 0.d0)then
          call tbedgemaxwell(trans,cod,rhob)
        else
          call tbedgedrift(trans,cod,psi)
          call tbedgemaxwell(trans1,cod,rhob)
          call tmultr5(trans,trans1,6)
          call tbedgebody(trans1,cod,rhob,psi)
          call tmultr5(trans,trans1,6)
        endif
      else
        if(psi .eq. 0.d0)then
          call tbedgemaxwell(trans,cod,-rhob)
        else
          call tbedgebody(trans,cod,rhob,psi)
          call tbedgemaxwell(trans1,cod,-rhob)
          call tmultr5(trans,trans1,6)
          call tbedgedrift(trans1,cod,psi)
          call tmultr5(trans,trans1,6)
        endif
      endif
      call tmultr5(trans0,trans,irad)
      call tmulbs(beam,trans,.true.)
      return
      end

      subroutine tbedgebody(trans,cod,rhob,psi)
      implicit none
      real*8 trans(6,6),rhob,cod(6),pr,psi,cosp,sinp,
     $     sinsq,pxi,pyi,s,dpzi,pzi,phsq,px0,dpz0,pz0,
     $     xi,pxf,pzf,dpzf,a,b,theta,aa,c,ddpz,s1,
     $     sqrs,sqrs1
      real*8, parameter :: ampmin=1.d-6,ampmax=1.d0-ampmin
      pr=1.d0+cod(6)
      cosp=cos(psi)
      sinp=sin(psi)
      sinsq=2.d0*sin(.5d0*psi)**2
      pxi=min(pr,max(-pr,cod(2)))
      pyi=min(pr,max(-pr,cod(4)))
      s=min(ampmax,pxi**2+pyi**2)
      sqrs=pr*sqrt(max(ampmin,1.d0-s/pr**2))
      dpzi=-s/(pr+sqrs)
      pzi=pr+dpzi
      phsq=pzi**2+pxi**2
      px0=pxi*cosp-pzi*sinp
      dpz0=pxi*sinp+dpzi*cosp-pr*sinsq
      pz0=pr+dpz0
      xi=cod(1)
      pxf=px0+xi*sinp/rhob
      cod(2)=min(pr,max(-pr,pxf))
      s1=min(ampmax,pyi**2+pxf**2)
      sqrs1=pr*sqrt(max(ampmin,1.d0-s1/pr**2))
      dpzf=-s1/(pr+sqrs1)
      pzf=pr+dpzf
      ddpz=(s-s1)/(sqrs+sqrs1)+(dpzi+pr)*sinsq-pxi*sinp
      cod(1)=rhob*ddpz+xi*cosp
      a=pr*xi*sinp/rhob+dpz0*pxf-dpzf*px0
      aa=a*cosp-(px0*pxf+pz0*pzf)*sinp
      c=(aa/phsq+sinp)/pzi/pzf
      b=dpzi+pr*sinsq-dpzf*cosp+pxf*sinp
      theta=asin(min(ampmax,max(-ampmax,a/phsq)))
      cod(3)=cod(3)-rhob*theta*pyi
      cod(5)=cod(5)+rhob*theta*pr
      trans(1,1)=cosp-pxf/pzf*sinp
      trans(1,2)=-rhob/pzi/pzf*a
      trans(1,3)=0.d0
      trans(1,4)=-pyi*rhob/pzi/pzf*b
      trans(1,5)=0.d0
      trans(1,6)= pr*rhob/pzi/pzf*b
      trans(2,1)=sinp/rhob
      trans(2,2)=pz0/pzi
      trans(2,3)=0.d0
      trans(2,4)=pyi/pzi*sinp
      trans(2,5)=0.d0
      trans(2,6)=-pr/pzi*sinp
      trans(3,1)=-pyi/pzf*sinp
      trans(3,2)=pyi*rhob/pzi/pzf*(dpzf-dpz0)
      trans(3,3)=1.d0
      trans(3,4)=-rhob*(pyi**2*c+theta)
      trans(3,5)=0.d0
      trans(3,6)=pr*pyi*rhob*c
      trans(4,1)=0.d0
      trans(4,2)=0.d0
      trans(4,3)=0.d0
      trans(4,4)=1.d0
      trans(4,5)=0.d0
      trans(4,6)=0.d0
      trans(5,1)=trans(2,1)*trans(1,6)-trans(1,1)*trans(2,6)
      trans(5,2)=trans(2,2)*trans(1,6)-trans(1,2)*trans(2,6)
      trans(5,3)=0.d0
      trans(5,4)=trans(2,4)*trans(1,6)-trans(1,4)*trans(2,6)
     $     +trans(3,6)
      trans(5,5)=1.d0
      trans(5,6)=rhob*(theta-pr**2*c)
      trans(6,1)=0.d0
      trans(6,2)=0.d0
      trans(6,3)=0.d0
      trans(6,4)=0.d0
      trans(6,5)=0.d0
      trans(6,6)=1.d0
      return
      end

      subroutine tbedgemaxwell(trans,cod,rhob)
      use mathfun, only: sqrtl
      implicit none
      real*8 trans(6,6),rhob,cod(6),pr,a,b,pxi,pyi,yi,pvi,
     $     c,af,bf,pvi2,yr
      real*8, parameter :: ampmin=1.d-6
      pr=1.d0+cod(6)
      pxi=min(pr,max(-pr,cod(2)))
      pyi=min(pr,max(-pr,cod(4)))
      pvi2=(pr-pxi)*(pr+pxi)
      pvi=pr*sqrtl(pvi2/pr**2)
      yi=cod(3)
      yr=(yi/rhob)**2
      a=(1.d0-yr/6.d0)
      af=a*yi/rhob/pvi
      b=(.5d0-yr/24.d0)
      bf=b*yi**2/rhob/pvi2/pvi
      c=(1.d0-yr/2.d0)
      cod(1)=cod(1)+bf*pr
      cod(4)=min(pr,max(-pr,cod(4)-af*pxi))
      cod(5)=cod(5)-bf*pxi
      trans(1,1)=1.d0
      trans(1,2)=3.d0*bf*pr**2*pxi/pvi2
      trans(1,3)=af*pr**2/pvi2
      trans(1,4)=0.d0
      trans(1,5)=0.d0
      trans(1,6)=-bf*pr*(pr**2+2.d0*pxi**2)/pvi2
      trans(2,1)=0.d0
      trans(2,2)=1.d0
      trans(2,3)=0.d0
      trans(2,4)=0.d0
      trans(2,5)=0.d0
      trans(2,6)=0.d0
      trans(3,1)=0.d0
      trans(3,2)=0.d0
      trans(3,3)=1.d0
      trans(3,4)=0.d0
      trans(3,5)=0.d0
      trans(3,6)=0.d0
      trans(4,1)=0.d0
      trans(4,2)=-trans(1,3)
      trans(4,3)=-c*pxi/rhob/pvi
      trans(4,4)=1.d0
      trans(4,5)=0.d0
      trans(4,6)=af*pr*pxi/pvi2
      trans(5,1)=0.d0
      trans(5,2)=trans(1,6)
      trans(5,3)=-trans(4,6)
      trans(5,4)=0.d0
      trans(5,5)=1.d0
      trans(5,6)=bf*pxi*(2.d0*pr**2+pxi**2)/pvi2
      trans(6,1)=0.d0
      trans(6,2)=0.d0
      trans(6,3)=0.d0
      trans(6,4)=0.d0
      trans(6,5)=0.d0
      trans(6,6)=1.d0
      return
      end

      subroutine tbedgedrift(trans,cod,psi)
      implicit none
      real*8 trans(6,6),psi,cod(6),cosp,sinp,
     $     s,pr,f,dpzi,pzi,sx,phsq,sxa,pxi,pyi
      real*8, parameter ::ampmin=1.d-6
      pr=1.d0+cod(6)
      cosp=cos(psi)
      sinp=sin(psi)
      pxi=min(pr,max(-pr,cod(2)))
      pyi=min(pr,max(-pr,cod(4)))
      s=pxi**2+pyi**2
      dpzi=-s/(pr+pr*sqrt(max(ampmin,1.d0-s/pr**2)))
      pzi=pr+dpzi
      f=cosp-pxi/pzi*sinp
      cod(1)=cod(1)/f
      sx=cod(1)*sinp
      sxa=sx/f/pzi**3
      cod(2)=pxi*cosp+pzi*sinp
      cod(3)=cod(3)+pyi/pzi*sx
      cod(5)=cod(5)-pr/pzi*sx
c      write(*,'(a,1p6g15.7)')'tbed ',psi,cod(5),-pr/pzi*sx
      phsq=pzi**2+pxi**2
      trans(1,1)=1.d0/f
      trans(1,2)=phsq*sxa
      trans(1,3)=0.d0
      trans(1,4)=pxi*pyi*sxa
      trans(1,5)=0.d0
      trans(1,6)=-pxi*pr*sxa
      trans(2,1)=0.d0
      trans(2,2)=f
      trans(2,3)=0.d0
      trans(2,4)=-pyi/pzi*sinp
      trans(2,5)=0.d0
      trans(2,6)=pr/pzi*sinp
      trans(3,1)=pyi/f/pzi*sinp
c      trans(3,2)=pyi*sxa*(f*pxi+phsq*sinp/pzi)
      trans(3,2)=pyi*sxa*(pxi*cosp+pzi*sinp)
      trans(3,3)=1.d0
c      trans(3,4)=sxa*(f*phsq+pxi*pyi**2*sinp/pzi)
      trans(3,4)=sxa*(f*pzi**2+pyi**2*cosp)
      trans(3,5)=0.d0
c      trans(3,6)=-pr*pyi*sxa*(f+pxi*sinp/pzi)
      trans(3,6)=-pr*pyi*sxa*(cosp/pzi)
      trans(4,1)=0.d0
      trans(4,2)=0.d0
      trans(4,3)=0.d0
      trans(4,4)=1.d0
      trans(4,5)=0.d0
      trans(4,6)=0.d0
      trans(5,1)=-trans(1,1)*trans(2,6)
      trans(5,2)= trans(2,2)*trans(1,6)-trans(1,2)*trans(2,6)
      trans(5,3)=0.d0
      trans(5,4)= trans(2,4)*trans(1,6)-trans(1,4)*trans(2,6)
     $     +trans(3,6)
      trans(5,5)=1.d0
      trans(5,6)=sxa*(f*pzi*s+pr**2*pxi*sinp)/pzi
      trans(6,1)=0.d0
      trans(6,2)=0.d0
      trans(6,3)=0.d0
      trans(6,4)=0.d0
      trans(6,5)=0.d0
      trans(6,6)=1.d0
      return
      end

      subroutine tbfrme(trans,cod,beam,ak,fb1,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      implicit none
      real*8 trans(6,6),trans1(6,6),cod(6),beam(42),ak(2),pr,
     $     akx,aky,xi,pxi,yi,pyi,y1,px1,a,dx,dpx,dy,dpy,dz,fb1,
     $     dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray,rhob
      logical*4 ent
      pr=1.d0+cod(6)
      akx= ak(1)
      aky=-ak(2)
      if(akx .ne. 0.d0 .and. fb1 .ne. 0.d0)then
        rhob=1.d0/akx
        dxfrx=fb1**2/rhob/24.d0
        dyfrx=fb1/rhob**2/6.d0
        dyfrax=4.d0*dyfrx/fb1**2
      else
        dxfrx=0.d0
        dyfrx=0.d0
        dyfrax=0.d0
      endif
      if(aky .ne. 0.d0 .and. fb1 .ne. 0.d0)then
        rhob=-1.d0/aky
        dyfry=fb1**2/rhob/24.d0
        dxfry=fb1/rhob**2/6.d0
        dxfray=4.d0*dxfry/fb1**2
      else
        dyfry=0.d0
        dxfry=0.d0
        dxfray=0.d0
      endif
      if((dxfrx .ne. 0.d0 .or. dyfry .ne. 0.d0) .and. ent)then
        call tbfrmle(trans,cod,beam,
     $       dxfrx,dyfrx,dyfrax,
     $       dxfry,dyfry,dxfray)
      endif
      trans1=0.d0
      xi =cod(1)
      yi =cod(3)
      pxi=min(pr,max(-pr,cod(2)))
      pyi=min(pr,max(-pr,cod(4)))
      y1 =-aky*xi +akx*yi
      px1= akx*pxi+aky*pyi
      a  =1.d0/(akx**2+aky**2)/pr
      dx = .5d0*a*akx*y1**2
      dpx= a*aky*px1*y1
      dy = .5d0*a*aky*y1**2
      dpy=-a*akx*px1*y1
      dz =-.5d0*a*px1*y1**2/pr
      trans1(1,1)=1.d0-a*akx*aky*y1
      trans1(1,3)= a*akx**2*y1
      trans1(1,6)=-dx/pr
      trans1(2,1)=-a*aky**2*px1
      trans1(2,2)=1.d0+a*akx*aky*y1
      trans1(2,3)= a*akx*aky*px1
      trans1(2,4)= a*aky**2*y1
      trans1(2,6)=-dpx/pr
      trans1(3,1)=-a*aky**2*y1
      trans1(3,3)=1.d0+a*akx*aky*y1
      trans1(3,6)=-dy/pr
      trans1(4,1)= a*akx*aky*px1
      trans1(4,2)=-a*akx**2*y1
      trans1(4,3)=-a*akx**2*px1
      trans1(4,4)=1.d0-a*akx*aky*y1
      trans1(4,6)=-dpy/pr
      trans1(5,1)= a*aky*px1*y1/pr
      trans1(5,2)=-.5d0*a*akx*y1**2/pr
      trans1(5,3)=-a*akx*px1*y1/pr
      trans1(5,4)=-.5d0*a*aky*y1**2/pr
      trans1(5,5)=1.d0
      trans1(5,6)=-2.d0*dz/pr
      trans1(6,6)=1.d0
      cod(1)=xi    +dx
      cod(2)=pxi   +dpx
      cod(3)=yi    +dy
      cod(4)=pyi   +dpy
      cod(5)=cod(5)+dz
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam,trans1,.true.)
      if((dxfrx .ne. 0.d0 .or. dyfry .ne. 0.d0) .and. .not. ent)then
        call tbfrmle(trans,cod,beam,
     $       dxfrx,dyfrx,dyfrax,
     $       dxfry,dyfry,dxfray)
      endif
      return
      end

      subroutine tbfrmle(trans,cod,beam,
     $     dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray)
      use tmacro
      use temw,only:tmulbs
      implicit none
      real*8 trans(6,6),cod(6),beam(42),trans1(6,6),
     $     dxfrx,dyfrx,dyfrax,
     $     dxfry,dyfry,dxfray,dpx,dpy,pr,xi,yi,dz
      trans1=0.d0
      pr=1.d0+cod(6)
      xi=cod(1)+dxfrx*cod(6)/pr
      yi=cod(3)-dyfry*cod(6)/pr
      dpx=(dxfry-dxfray*xi**2)/pr*xi
      dpy=(dyfrx-dyfrax*yi**2)/pr*yi
      cod(2)=cod(2)+dpx
      cod(4)=cod(4)+dpy
      dz=(dxfrx*cod(2)-dyfry*cod(4)+
     $     (.5d0*dyfrx-.25d0*dyfrax*yi**2)*yi**2
     $     +(.5d0*dxfry-.25d0*dxfray*xi**2)*xi**2)/pr**2
      cod(1)=xi
      cod(3)=yi
      cod(5)=cod(5)+dz
      trans1(1,1)=1.d0
      trans1(1,6)=dxfrx/pr**2
      trans1(2,1)=(dxfry-3.d0*dxfry*xi**2)/pr
      trans1(2,2)=1.d0
      trans1(2,6)=-dpx/pr-3.d0*xi**2/pr*trans1(1,6)
      trans1(3,3)=1.d0
      trans1(3,6)=-dyfry/pr**2
      trans1(4,3)=(dyfrx-3.d0*dyfrx*yi**2)/pr
      trans1(4,4)=1.d0
      trans1(4,6)=-dpy/pr-3.d0*yi**2/pr*trans1(3,6)
      trans1(5,1)=-trans1(2,6)+trans1(2,1)*trans1(1,6)
      trans1(5,2)= trans1(1,6)
      trans1(5,3)=-trans1(4,6)+trans1(4,3)*trans1(3,6)
      trans1(5,4)= trans1(3,6)
      trans1(5,5)= 1.d0
      trans1(5,6)= (dxfrx*trans1(2,6)-dyfry*trans1(4,6)+
     $     (dyfrx*yi-dyfrax*yi**2)*yi*trans1(3,6)+
     $     (dxfry*xi-dxfray*xi**2)*xi*trans1(1,6))/pr**2
     $     -2.d0*dz/pr
      trans1(6,6)=1.d0
      call tmultr5(trans,trans1,irad)
      call tmulbs(beam,trans1,.true.)
      return
      end
