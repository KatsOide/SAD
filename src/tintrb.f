      module intrb
      implicit none
      real*8 a1,a2
      end module

      subroutine tintrb(trans,cod,beam,bmi,al,al1,optics,ll)
      use intrb
      use temw, only:diagr=>r, diagri=>ri,tinv6,eemx,eemy,eemz,caltouck,
     $     tmulbs
      use touschek_table
      use tfstk
      use ffs_flag
      use ffs_pointer ,only:beamsize
      use tmacro
      use mathfun
      use sad_main, ia=>iaidx
      implicit none
      real*8 fintrb
      external fintrb
      integer*4 ,intent(in):: ll
      integer*4 i,j
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42)
      real*8 ,intent(out):: bmi(21)
      real*8 ,intent(in):: al,al1
      real*8 pl(3,3),r(3,3),eig(3),xx(3,3),xxs(3,3),
     $     touckf,bint,e1,e2,e3
      real*8 xp(3,3),transw(6,6),
     $     pxi,pyi,s,pr,pzi,alx,ale,alz,hi,a,b,d,vol,
     $     bm,ptrans,extrans,eytrans,eztrans,tf,aez,aex0,aey0,
     $     aez0,aexz,aeyz,f1,f2,f3,bn,bmax,bmin,ci,pvol,vol1,
     $     transsp(6,6)
      real*8 trans1(6,6),trans2(6,6)
      logical*4 ,intent(in):: optics
c     real*8  vmin/0.d0/
      if(al .eq. 0.d0)then
        bmi=0.d0
        return
      endif
      pxi=cod(2)
      pyi=cod(4)
      s=pxi**2+pyi**2
      pr=1.d0+cod(6)
      pzi=pr*(1.d0+sqrt1(-s/pr**2))
      call tinitr(trans1)
      alx=.5d0*al-al1
      ale=alx/pzi
      alz=ale/pzi**2
      trans1(1,2)=ale+pxi**2*alz
      trans1(1,4)=pxi*pyi*alz
      trans1(1,6)=-pxi*alz*pr
      trans1(3,2)=trans1(1,4)
      trans1(3,4)=ale+pyi**2*alz
      trans1(3,6)=-pyi*alz*pr
      trans1(5,2)=trans1(1,6)
      trans1(5,4)=trans1(3,6)
      hi=p2h(p0*pr)
      trans1(5,6)=h0/hi**3*alx+s*alz
      if(wspac)then
        if(optics)then
          bmi=beamsize(:,ll)
        else
          bmi=beam(22:42)+beam(1:21)
        endif
        call tmulbs(bmi,trans1,.false.)
        call twspace(transsp,cod,al,bmi,ll)
c        write(*,*)'tintrab-wspac-twspace-end '
c        write(*,'(1p6g16.7)')(transsp(i,:),i=1,6)
        trans2=matmul(tinv6(trans1),matmul(transsp,trans1))
        trans(:,1:irad)=matmul(trans2,trans(:,1:irad))
        if(optics)then
          return
        endif
        call tmulbs(beam,trans2,.false.)
      endif
      if(intra)then
        bmi=beam(22:42)+beam(1:21)
c        call tmov(beam(22),bmi,21)
c        call tmov(trans,transa,36)
c      call tadd(transa,trans(1,7),transa,36)
c        call tmulbs(bmi,transa,.false.)
c        call tadd(bmi,beam,bmi,21)
        if(caltouck)then
          transw=matmul(trans(:,1:6),diagr)
c          call tmov(diagr,transw,36)
c          call tmultr(transw,trans,6)
        endif
        a=p0**2/(hi+1.d0)
        b=a/hi
        d=1.d0/hi
        call tinitr(trans2)
        trans2(1,1)=1.d0+a*pxi**2
        trans2(1,3)=a*pxi*pyi
        trans2(1,5)=a*pxi*pzi
        trans2(3,1)=trans2(1,3)
        trans2(3,3)=1.d0+a*pyi**2
        trans2(3,5)=a*pyi*pzi
        trans2(5,1)=trans2(1,5)
        trans2(5,3)=trans2(3,5)
        trans2(5,5)=1.d0+a*pzi**2
        trans2(2,2)=d+b*(pzi**2+pyi**2)
        trans2(2,4)=-b*pxi*pyi
        trans2(2,6)=-b*pxi*pzi
        trans2(4,2)=trans2(2,4)
        trans2(4,4)=d+b*(pzi**2+pxi**2)
        trans2(4,6)=-b*pyi*pzi
        trans2(6,2)=trans2(2,6)
        trans2(6,4)=trans2(4,6)
        trans2(6,6)=d+b*s
        do 3010 i=1,6
          trans2(i,2)=trans2(i,2)-pxi/pzi*trans2(i,6)
          trans2(i,4)=trans2(i,4)-pyi/pzi*trans2(i,6)
          trans2(i,5)=(pxi*trans2(i,1)+pyi*trans2(i,3)
     $         +pzi*trans2(i,5))/pr
          trans2(i,6)=pr/pzi*trans2(i,6)
 3010   continue
        trans1=matmul(trans2,trans1)
        call tmulbs(bmi,trans1,.false.)
        xx(1,1)=bmi(ia(1,1))
        xx(2,1)=bmi(ia(3,1))
        xx(3,1)=bmi(ia(5,1))
        xx(2,2)=bmi(ia(3,3))
        xx(3,2)=bmi(ia(5,3))
        xx(3,3)=bmi(ia(5,5))
        call eigs33(xx,r,eig)
        vol1=sqrt(max(1.d-80,eig(1)*eig(2)*eig(3)))
        vol=sqrt((4.d0*pi)**3)*vol1
        bm=sqrt(min(abs(eig(1)),abs(eig(2)),abs(eig(3))))
        xxs=xx
c        call tmov(xx,xxs,9)
        xp(1,1)=bmi(ia(1,2))
        xp(1,2)=bmi(ia(1,4))
        xp(1,3)=bmi(ia(1,6))
        xp(2,1)=bmi(ia(3,2))
        xp(2,2)=bmi(ia(3,4))
        xp(2,3)=bmi(ia(3,6))
        xp(3,1)=bmi(ia(5,2))
        xp(3,2)=bmi(ia(5,4))
        xp(3,3)=bmi(ia(5,6))
        call sols33(xxs,xp)
        pl(1,1)=bmi(ia(2,2))
     1       -bmi(ia(2,1))*xp(1,1)-bmi(ia(2,3))*xp(2,1)
     $       -bmi(ia(2,5))*xp(3,1)
        pl(2,1)=bmi(ia(4,2))
     1       -bmi(ia(4,1))*xp(1,1)-bmi(ia(4,3))*xp(2,1)
     $       -bmi(ia(4,5))*xp(3,1)
        pl(3,1)=bmi(ia(6,2))
     1       -bmi(ia(6,1))*xp(1,1)-bmi(ia(6,3))*xp(2,1)
     $       -bmi(ia(6,5))*xp(3,1)
        pl(2,2)=bmi(ia(4,4))
     1       -bmi(ia(4,1))*xp(1,2)-bmi(ia(4,3))*xp(2,2)
     $       -bmi(ia(4,5))*xp(3,2)
        pl(3,2)=bmi(ia(6,4))
     1       -bmi(ia(6,1))*xp(1,2)-bmi(ia(6,3))*xp(2,2)
     $       -bmi(ia(6,5))*xp(3,2)
        pl(3,3)=bmi(ia(6,6))
     1       -bmi(ia(6,1))*xp(1,3)-bmi(ia(6,3))*xp(2,3)
     $       -bmi(ia(6,5))*xp(3,3)
        call eigs33(pl,r,eig)
        ptrans=sqrt(abs(eig(1)+eig(2)+eig(3)))
        pvol=sqrt(max(1.d-80,eig(1)*eig(2)*eig(3)))
        if(vol .ne. 0.d0 .and. caltouck)then
          if(ptrans .ne. 0.d0)then
            trans2=tinv6(matmul(trans1,transw))
            extrans=(trans2(1,6)**2+trans2(2,6)**2)*ptrans**2
            eytrans=(trans2(3,6)**2+trans2(4,6)**2)*ptrans**2
            eztrans=(trans2(5,6)**2+trans2(6,6)**2)*ptrans**2
            tf=al/vol/(ptrans*p0)**3
            if(eztrans .ne. 0.d0)then
              do i=1,ntouckl
                aez=(.002d0+i*.002d0)**2
     $               *(diagri(5,6)**2+diagri(6,6)**2)
                touckl(i)=touckl(i)+touckf(aez/eztrans)*tf
                toucke(i,ll)=toucke(i,ll)+touckf(aez/eztrans)*tf
              enddo
              do i=1,ntouckx
                aex0=(tampl(i,1)**2)*(abs(eemx)+abs(eemy))
                aey0=(tampl(i,2)**2)*(abs(eemx)+abs(eemy))
                do  j=1,ntouckz
                  aez0=(tampl(j,3)**2)*abs(eemz)
                  aexz=1.d0/(extrans/aex0+eztrans/aez0)
                  touckm(j,i,1)=touckm(j,i,1)+touckf(aexz)*tf
                  aeyz=1.d0/(eytrans/aey0+eztrans/aez0)
                  touckm(j,i,2)=touckm(j,i,2)+touckf(aeyz)*tf
                enddo
              enddo
            endif
          endif
        else
          touckl=1.d20
          touckm=1.d20
        endif
        a1=eig(1)/eig(2)
        a2=eig(1)/eig(3)
        f1=2.d0*bint(fintrb,0.d0,hpi,1.d-3,1.d-19)*eig(1)
        a1=eig(2)/eig(3)
        a2=eig(2)/eig(1)
        f2=2.d0*bint(fintrb,0.d0,hpi,1.d-3,1.d-19)*eig(2)
        a1=eig(3)/eig(1)
        a2=eig(3)/eig(2)
        f3=2.d0*bint(fintrb,0.d0,hpi,1.d-3,1.d-19)*eig(3)
        e1=f2+f3-2.d0*f1
        e2=f3+f1-2.d0*f2
        e3=f1+f2-2.d0*f3
        pl(1,1)=r(1,1)*r(1,1)*e1+r(1,2)*r(1,2)*e2+r(1,3)*r(1,3)*e3
        pl(2,1)=r(2,1)*r(1,1)*e1+r(2,2)*r(1,2)*e2+r(2,3)*r(1,3)*e3
        pl(3,1)=r(3,1)*r(1,1)*e1+r(3,2)*r(1,2)*e2+r(3,3)*r(1,3)*e3
        pl(2,2)=r(2,1)*r(2,1)*e1+r(2,2)*r(2,2)*e2+r(2,3)*r(2,3)*e3
        pl(3,2)=r(3,1)*r(2,1)*e1+r(3,2)*r(2,2)*e2+r(3,3)*r(2,3)*e3
        pl(3,3)=r(3,1)*r(3,1)*e1+r(3,2)*r(3,2)*e2+r(3,3)*r(3,3)*e3
        bn=abs(vol/pbunch)**(1.d0/3.d0)
        bmax=max(1.d-80,min(bm,bn))
        bmin=max(rclassic/(ptrans*p0)**2,
     1       sqrt(abs(
     1       vol/pi/(ptrans*p0/h0*c)/pbunch
     $       /max(taurdx,taurdy,taurdz)))
     1       )
c     write(*,*)'bmin,bmax,vol ',bmin,bmax,vol
c     bmin=rclassic
        ci=cintrb*al*log(2.d0*bmax/bmin)/vol1/pvol/h0**4
c        write(*,*)'tintrb ',ll,vol1*pvol,ci,pl
c     write(*,*)log(2.d0*bmax/bmin)
c     if(ci .gt. vmin)then
c     vmin=ci
c     write(*,*)vmin,cintrb,bmin
c     endif
        bmi=0.d0
        bmi(ia(2,2))=ci*pl(1,1)
        bmi(ia(4,2))=ci*pl(2,1)
        bmi(ia(6,2))=ci*pl(3,1)
        bmi(ia(4,4))=ci*pl(2,2)
        bmi(ia(6,4))=ci*pl(3,2)
        bmi(ia(6,6))=ci*pl(3,3)
        call tmulbs(bmi,tinv6(trans1),.false.)
        beam(1:21)=bmi+beam(1:21)
c        if(mod(ll,10) .eq. 0)then
c          write(*,'(a,i5,1p8g14.6)')'tintrb ',ll,ci,al,
c     $         beam(3),beam(10),beam(21),bmi(3),bmi(10),bmi(21)
c        endif
c        call tadd(beam,bmi,beam,21)
c        call tinv6(trans,transa)
        call tmulbs(bmi,tinv6(trans(:,1:6)),.false.)
      endif
      return
      end

      real*8 pure function fintrb(t)
      use intrb
      implicit none
      real*8 , intent(in)::t
      real*8 cost,sqsint
      cost=cos(t)
      sqsint=(1.d0-cost)*(1.d0+cost)
      fintrb=sqsint*cost/
     1       sqrt((sqsint+a1*cost**2)*(sqsint+a2*cost**2))
c     fintrb=t/sqrt(((1.d0-a1)*t+a1)*((1.d0-a2)*t+a2))
      return
      end

      real*8 function touckf(x)
c
c Approximation of C[x]/x \propto 1/tau .
c
      implicit none
      real*8 x,eeuler,a,b
c      parameter (eeuler=7.98221278918726d0,a=5.5077d0,b=1.1274d0)
      parameter (eeuler=7.982212789187259d0,a=5.62966d0,b=0.75159d0)
      if(x .eq. 0.d0)then
        touckf=1.d200
      else
        touckf=(log(1.d0/x/eeuler+1.d0)*exp(-x)
     1          *(b+eeuler*x)/(b+x*(a+2.d0*x)))/x
      endif
      return
      end

      subroutine twspace(trans,cod,al,beam,l)
      use tfstk
      use tmacro
      use mathfun
      use sad_main, ia=>iaidx
      use ffs_pointer,only:gammab
      implicit none
      integer*4 ,intent(in):: l
      real*8 ,intent(in):: al
      real*8 ,intent(out):: trans(6,6)
      real*8 ,intent(in):: cod(6),beam(21)
      real*8 xx1,yy1,xy1,u,v,a,c2,s2,sx,sy,p1,h1,f,akx,aky,
     $     aks,akd,sigzsq
c      real*8 fx,fy,fu,fxx,fyy,fxy
      sigzsq=beam(ia(5,5))
      xx1=beam(ia(1,1))-beam(ia(5,1))**2/sigzsq
      yy1=beam(ia(3,3))-beam(ia(5,3))**2/sigzsq
      xy1=beam(ia(3,1))-beam(ia(5,1))*beam(ia(5,3))/sigzsq
      u=xx1-yy1
      v=2.d0*xy1
      a=hypot(u,v)
      if(a .eq. 0.d0)then
        c2=1.d0
        s2=0.d0
      else
        c2=u/a
        s2=v/a
      endif
      sx=(xx1+yy1)*.5d0
      sy=sqrt(sx-a*.5d0)
      sx=sqrt(sx+a*.5d0)
      p1=(1.d0+cod(6))*(gammab(l)+gammab(l+1))*.5d0
      h1=p2h(p1)
c      call twspfu(0.d0,0.d0,sx,sy,fx,fy,fu,fxx,fyy,fxy)
c      f1=pbunch*rclassic*al/(p1**2*h1*sqrt(2.d0*pi*sigzsq))
c      akx=-f*fxx
c      aky=-f*fyy
      f=2.d0*pbunch*rclassic*al/(p1**2*h1*(sx+sy)
     $     *sqrt(2.d0*pi*sigzsq))
c      write(*,'(a,1p8g15.7)')'twspace ',f,sy,a,u,v,xx1,yy1
      akx=f/sx
      aky=f/sy
c      write(*,'(a,i5,1p6g15.7)')'twspace ',
c     $     l,akx,aky,-.5d0*fxx*f*(sx+sy),-.5d0*fyy*f*(sx+sy)
      aks=(akx+aky)*.5d0
      akd=(akx-aky)*.5d0
      call tinitr(trans)
      trans(2,1)=aks+akd*c2
      trans(2,3)=akd*s2
      trans(4,1)=trans(2,3)
      trans(4,3)=aks-akd*c2
c      if(l .lt. 10)then
c        write(*,*)'tspace ',l,gammab(l),gammab(l+1)
c      endif
      return
      end

      subroutine qwspac(trans,cod,al,beam,coup,l)
      use tfstk
      use tmacro
      implicit none
      integer*4 ,intent(in):: l
      real*8 ,intent(in):: al
      real*8 ,intent(inout)::  trans(4,5),cod(6),beam(42)
      real*8 trs(6,6)
      logical*4 ,intent(in):: coup
      if(al .eq. 0.d0)then
        return
      endif
      call twspace(trs,cod,al,beam,l)
      if(coup)then
        trans(1,1)=trans(1,1)+trans(1,2)*trs(2,1)+trans(1,4)*trs(4,1)
        trans(2,1)=trans(2,1)+trans(2,2)*trs(2,1)+trans(2,4)*trs(4,1)
        trans(3,1)=trans(3,1)+trans(3,2)*trs(2,1)+trans(3,4)*trs(4,1)
        trans(4,1)=trans(4,1)+trans(4,2)*trs(2,1)+trans(4,4)*trs(4,1)
        trans(1,3)=trans(1,3)+trans(1,2)*trs(2,3)+trans(1,4)*trs(4,3)
        trans(2,3)=trans(2,3)+trans(2,2)*trs(2,3)+trans(2,4)*trs(4,3)
        trans(3,3)=trans(3,3)+trans(3,2)*trs(2,3)+trans(3,4)*trs(4,3)
        trans(4,3)=trans(4,3)+trans(4,2)*trs(2,3)+trans(4,4)*trs(4,3)
      else
        trans(1,1)=trans(1,1)+trans(1,2)*trs(2,1)
        trans(2,1)=trans(2,1)+trans(2,2)*trs(2,1)
        trans(3,1)=                               trans(3,4)*trs(4,1)
        trans(4,1)=                               trans(4,4)*trs(4,1)
        trans(3,2)=0.d0
        trans(4,2)=0.d0
        trans(1,3)=           trans(1,2)*trs(2,3)
        trans(2,3)=           trans(2,2)*trs(2,3)
        trans(3,3)=trans(3,3)                    +trans(3,4)*trs(4,3)
        trans(4,3)=trans(4,3)                    +trans(4,4)*trs(4,3)
        trans(1,4)=0.d0
        trans(2,4)=0.d0
      endif
      return
      end

      subroutine twspac(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     al,cod,beam,l)
      use tfstk
      use ffs_flag
      use tmacro
      use mathfun
      use sad_main, ia=>iaidx
      use ffs_pointer,only:gammab
      implicit none
      integer*4 ,intent(in):: l
      integer*4 np,i
      real*8 ,intent(inout):: x(np),px(np),y(np),py(np),
     $     z(np),g(np),dv(np),sx(np),sy(np),sz(np)
      real*8 ,intent(in):: al,cod(6),beam(21)
      real*8 xx1,yy1,xy1,a,c1,s1,sigx,sigy,p1,h1,f,dx,dy,
     $     dx1,dy1,dpx,dpy,pr,u,v,theta,sigzsq,
     $     az,dg,dpr,pr1,fx,fy,fu,xc,yc,zc,fxx,fyy,fxy
      sigzsq=beam(ia(5,5))
      xx1=beam(ia(1,1))-beam(ia(5,1))**2/sigzsq
      yy1=beam(ia(3,3))-beam(ia(5,3))**2/sigzsq
      xy1=beam(ia(3,1))-beam(ia(5,1))*beam(ia(5,3))/sigzsq
      u=xx1-yy1
      v=2.d0*xy1
      a=hypot(u,v)
      theta=atan2(v,u)*.5d0
      c1=cos(theta)
      s1=sin(theta)
      sigx=(xx1+yy1)*.5d0
      sigy=sqrt(sigx-a*.5d0)
      sigx=sqrt(sigx+a*.5d0)
      if(selfcod)then
        xc=sum(x)/np
        yc=sum(y)/np
        zc=sum(z)/np
      else
        xc=cod(1)
        yc=cod(3)
        zc=cod(5)
      endif
      p1=(1.d0+cod(6))*(gammab(l)+gammab(l+1))*.5d0
      h1=p2h(p1)
c      h1=p1*sqrt(1.d0+1.d0/p1**2)
      f=pbunch*rclassic*al/(p1**2*h1*sqrt(2.d0*pi*sigzsq))
      do i=1,np
        dx=x(i)-xc
        dy=y(i)-yc
        dx1= dx*c1+dy*s1
        dy1=-dx*s1+dy*c1
        az=f*exp(-.5d0*(z(i)-zc)**2/sigzsq)
        call twspfu(dx1,dy1,sigx,sigy,fx,fy,fu,fxx,fyy,fxy)
c        bb=bbkick1(dx1,dy1,sigx,sigy)
c        call bbkick(dcmplx(sigx,sigy),dcmplx(dx1,dy1),
c     $       bb,1,tr)
c        fx=az*dble(bb)
c        fy=az*dimag(bb)
c        dg=-(z(i)-zc)/sigzsq*az*
c     $       twspu(dx1,dy1,sigx,sigy,4.d-2,1.d-4)
c        if(l .lt. 10 .and. i .eq. 1)then
c          write(*,*)l,fx,fxx*dx1,fy,fyy*dy1
c        endif
        fx=-az*fx
        fy=-az*fy
        dg=-(z(i)-zc)/sigzsq*az*fu
        dpx= fx*c1-fy*s1
        dpy= fx*s1+fy*c1
        dpr=g(i)
        pr=1.d0+dpr
        dpr=dpr+dg
        pr1=1.d0+dpr
        px(i)=(px(i)*pr+dpx)/pr1
        py(i)=(py(i)*pr+dpy)/pr1
        g(i)=dpr
        h1=p2h(p0*pr1)
        dv(i)=-dpr*(1.d0+pr1)/h1/(h1+pr1*h0)+dvfs
      enddo
      return
      end

      module wspf
      implicit none
      real*8 a,b,r,eps
      end module

      real*8 function twspu(x,y,sigx,sigy,epslon,epsabs)
      use wspf
      use gammaf, only:gamma0log
      implicit none
      real*8 x,y,sigx,sigy,u
      real*8 twspurec,epslon,epsabs,
     $     c,d,f0,f1,f2,f3,d2,
     $     ax,ay,r1,sig2,twspf1,rombint
      external twspf1
      if(sigx .lt. sigy)then
        twspu=twspurec(y,x,sigy,sigx,epslon,epsabs)
        return
      endif
      sig2=sigx*(sigx+sigy)
      ax=x/sig2
      ay=y/sig2
      a=.5d0*x*ax
      b=.5d0*y*ay
      r=sigy/sigx
      if(r .eq. 1.d0)then
        u=2.d0*(a+b)
        if(u .eq. 0.d0)then
          twspu=0.d0
        else
          twspu=gamma0log(u)
        endif
        return
      endif
      twspu=2.d0*(a+b/r)
      if(twspu .le. min(epslon,epsabs*10d0))then
        twspu=2.d0*(a+b/r)+log(.5d0+.5d0*r)
      elseif(a+b .lt. 30.d0)then
        twspu=2.d0*rombint(twspf1,0.d0,1.d0,epslon,epsabs)
     $       +log(.5d0+.5d0*r)
      else
        r1=1.d0-r
        u=a+b
        c=r1/u
        d=(a-b)/u
        d2=d**2
        f0=1.d0-d2
        f1=f0-d2
        f2=3.d0-4.d0*d2
        f3=4.d0*f0*f1-f2
        twspu=gamma0log(2.d0*u)
     $       -c*(64.d0*d+c*(-48.d0*f1
     $       +c*(-80.d0*d*f2+c*(210.d0*f3
     $       +c*(756.d0*d*(f2+2.d0*f3)
     $       +c*(3465.d0*(2.d0*(d*f2)**2-1.d0)
     $       +c*(19305.d0*d*(1.d0-4.d0*f1*(f3+f1)))))))))/128.d0
c        dy(1)=-(a+b)**2+2.d0*b*(1.d0-r)
c        r1=1.d0+r
c        ev=exp(-(a+b/r**2)*r1)
c        dy(2)=2.d0*(ev*(b/r**3+r*(a+1.d0/r1))-r/r1)/r1
c        twspu=2.d0*splint(twspf1,0.d0,1.d0,3,dy,
c     $       epslon,epsabs,8)+rlog
      endif
      return
      end

      real*8 function twspurec(x,y,sigx,sigy,epslon,epsabs)
      implicit none
      real*8 x,y,sigy,sigx,twspu,epslon,epsabs
      twspurec=twspu(x,y,sigx,sigy,epslon,epsabs)
      return
      end

      real*8 pure elemental function twspf1(v)
      use wspf
      implicit none
      real*8 ,intent(in):: v
      real*8 f,t,x,w,tsq
      t=(1.d0-r)*(1.d0-v)+r
      tsq=t**2
      f=a+b/tsq
      w=v*(1.d0+t)
      x=w*f
      if(x .lt. 0.001d0)then
        twspf1=f*(1.d0-x/2.d0*(1.d0-x/3.d0*
     $       (1.d0-x/4.d0*(1.d0-x/5.d0*
     $       (1.d0-x/6.d0)))))
      else
        twspf1=(1.d0-exp(-x))/w
      endif
      return
      end

      subroutine twspfu(x,y,sigx,sigy,fx,fy,fu,fxx,fyy,fxy)
      use tfstk
      use iso_c_binding
c      use tmacro, only:l_track
      implicit none
      integer*8 ,save:: iu=0
      integer*4 ,parameter::nr=20,nx=60,ny=60,m=(nr+1)*(nx+1)*(ny+1)
      real*8 ,intent(in):: x,y,sigy,sigx
      real*8 ,intent(out):: fx,fy,fu,fxx,fyy,fxy
      real*8 ,pointer,save::u(:,:,:),uxx(:,:,:),uyy(:,:,:),uxxyy(:,:,:),
     $     urr(:,:,:),uxxrr(:,:,:),uyyrr(:,:,:),uxxyyrr(:,:,:)     
      if(iu .eq. 0)then
        iu=ktaloc(8*m)
        call c_f_pointer(c_loc(rlist(iu)),u,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+m)),uxx,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+2*m)),uyy,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+3*m)),uxxyy,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+4*m)),urr,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+5*m)),uxxrr,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+6*m)),uyyrr,[nx+1,ny+1,nr+1])
        call c_f_pointer(c_loc(rlist(iu+7*m)),uxxyyrr,[nx+1,ny+1,nr+1])
        call twspfuinit(
     $     u,uxx,uyy,uxxyy,urr,uxxrr,uyyrr,uxxyyrr)
      endif
      if(sigx .lt. sigy)then
        call twspfu0(y,x,sigy,sigx,fy,fx,fu,fxx,fyy,fxy,
     $     u,uxx,uyy,uxxyy,urr,uxxrr,uyyrr,uxxyyrr)
      else
        call twspfu0(x,y,sigx,sigy,fx,fy,fu,fxx,fyy,fxy,
     $     u,uxx,uyy,uxxyy,urr,uxxrr,uyyrr,uxxyyrr)
      endif
      return
      end

      subroutine twspfu0(x,y,sigx,sigy,fx,fy,fu,fxx,fyy,fxy,
     $     u,uxx,uyy,uxxyy,urr,uxxrr,uyyrr,uxxyyrr)
c      use tmacro, only:l_track
      implicit none
      integer*4 nr,nx,i,j,n,ny
      real*8 ,intent(in):: x,y,sigx,sigy
      real*8 ,intent(out):: fx,fy,fu,fxx,fyy,fxy
      real*8 xm,xstep,r,rstep,rm,ym,ystep
      parameter (nr=20,nx=60,ny=60,xm=15.d0,ym=30.d0,rm=5.d0)
      parameter (xstep=xm/nx,rstep=rm/nr**2,ystep=ym/ny)
      real*8 ,intent(in)::u(0:nx,0:ny,0:nr),uxx(0:nx,0:ny,0:nr),
     $     uyy(0:nx,0:ny,0:nr),uxxyy(0:nx,0:ny,0:nr)
      real*8 ,intent(in):: urr(0:nx,0:ny,0:nr),uxxrr(0:nx,0:ny,0:nr),
     $     uyyrr(0:nx,0:ny,0:nr),uxxyyrr(0:nx,0:ny,0:nr)
      real*8 rl,ax,ay,aax,aay,ar,bax,bay,aax2,aay2,ar2,
     $     br,br2,bax2,bay2,u0,u1,u2,u3,twspu,dxs,dys,
     $     aax23,bax23,aay23,bay23,
     $     uxx0,uxx1,uxx2,uxx3,uyy0,uyy1,uyy2,uyy3,
     $     uxxyy0,uxxyy1,uxxyy2,uxxyy3,
     $     up,uq,ur,us,uyyp,uyyq,uxxr,uxxs
      complex*16 z,bbkick1
      integer*4 i1,j1,n1
      r=sigy/sigx
      rl=sqrt(-log10(r)/rstep)
      n=int(rl)
      if(n .ge. nr)then
        go to 9000
      endif
      dxs=xstep*sigx
      dys=ystep*sigy
      ax=min(abs(x)/dxs,xm*100.d0)
      i=int(ax)
      if(i .ge. nx)then
        go to 9000
      endif
      ay=min(abs(y)/dys,ym*100.d0)
      j=int(ay)
      if(j .ge. ny)then
        go to 9000
      endif
      i1=i+1
      j1=j+1
      n1=n+1
      bax=ax-i
      bay=ay-j
      br=rl-n
      aax=1.d0-bax
      aay=1.d0-bay
      ar=1.d0-br
      ar2=-br*(ar+1.d0)
      br2=-ar*(br+1.d0)
      aax2=-bax*(aax+1.d0)
      bax2=-aax*(bax+1.d0)
      aay2=-bay*(aay+1.d0)
      bay2=-aay*(bay+1.d0)
c      daay2=-(aay+1.d0)+bay=-aay-aay=-2*aay
c      dbay2=-(bay+1)+aay=-bay-bay=-2*bay
c      write(*,*)'twspfu ',sigx,sigy
      u0=ar*(u(i,j,n)+ar2*urr(i,j,n))+
     $     br*(u(i,j,n1)+br2*urr(i,j,n1))
      u1=ar*(u(i1,j,n)+ar2*urr(i1,j,n))+
     $     br*(u(i1,j,n1)+br2*urr(i1,j,n1))
      u2=ar*(u(i,j1,n)+ar2*urr(i,j1,n))+
     $     br*(u(i,j1,n1)+br2*urr(i,j1,n1))
      u3=ar*(u(i1,j1,n)+ar2*urr(i1,j1,n))+
     $     br*(u(i1,j1,n1)+br2*urr(i1,j1,n1))
      uxx0=ar*(uxx(i,j,n)+ar2*uxxrr(i,j,n))+
     $     br*(uxx(i,j,n1)+br2*uxxrr(i,j,n1))
      uxx1=ar*(uxx(i1,j,n)+ar2*uxxrr(i1,j,n))+
     $     br*(uxx(i1,j,n1)+br2*uxxrr(i1,j,n1))
      uxx2=ar*(uxx(i,j1,n)+ar2*uxxrr(i,j1,n))+
     $     br*(uxx(i,j1,n1)+br2*uxxrr(i,j1,n1))
      uxx3=ar*(uxx(i1,j1,n)+ar2*uxxrr(i1,j1,n))+
     $     br*(uxx(i1,j1,n1)+br2*uxxrr(i1,j1,n1))
      uyy0=ar*(uyy(i,j,n)+ar2*uyyrr(i,j,n))+
     $     br*(uyy(i,j,n1)+br2*uyyrr(i,j,n1))
      uyy1=ar*(uyy(i1,j,n)+ar2*uyyrr(i1,j,n))+
     $     br*(uyy(i1,j,n1)+br2*uyyrr(i1,j,n1))
      uyy2=ar*(uyy(i,j1,n)+ar2*uyyrr(i,j1,n))+
     $     br*(uyy(i,j1,n1)+br2*uyyrr(i,j1,n1))
      uyy3=ar*(uyy(i1,j1,n)+ar2*uyyrr(i1,j1,n))+
     $     br*(uyy(i1,j1,n1)+br2*uyyrr(i1,j1,n1))
      uxxyy0=ar*(uxxyy(i,j,n)+ar2*uxxyyrr(i,j,n))+
     $     br*(uxxyy(i,j,n1)+br2*uxxyyrr(i,j,n1))
      uxxyy1=ar*(uxxyy(i1,j,n)+ar2*uxxyyrr(i1,j,n))+
     $     br*(uxxyy(i1,j,n1)+br2*uxxyyrr(i1,j,n1))
      uxxyy2=ar*(uxxyy(i,j1,n)+ar2*uxxyyrr(i,j1,n))+
     $     br*(uxxyy(i,j1,n1)+br2*uxxyyrr(i,j1,n1))
      uxxyy3=ar*(uxxyy(i1,j1,n)+ar2*uxxyyrr(i1,j1,n))+
     $     br*(uxxyy(i1,j1,n1)+br2*uxxyyrr(i1,j1,n1))
      up=aax*(u0+aax2*uxx0)+bax*(u1+bax2*uxx1)
      uq=aax*(u2+aax2*uxx2)+bax*(u3+bax2*uxx3)
      uyyp=aax*(uyy0+aax2*uxxyy0)+bax*(uyy1+bax2*uxxyy1)
      uyyq=aax*(uyy2+aax2*uxxyy2)+bax*(uyy3+bax2*uxxyy3)
      fu=aay*(up+aay2*uyyp)+bay*(uq+bay2*uyyq)
      aax23=3.d0*aax2+2.d0
      bax23=3.d0*bax2+2.d0
      aay23=3.d0*aay2+2.d0
      bay23=3.d0*bay2+2.d0
      fy=((uq-up)-aay23*uyyp+bay23*uyyq)/dys
      fyy=-6.d0*(aay*uyyp+bay*uyyq)/dys**2
c      daay2=-(aay+1.d0)+bay=-aay-aay=-2*aay
c      write(*,*)'twspfu ',ay,j,fy,uyy0,uyy2,uyyp,uyyq
      ur=aay*(u0+aay2*uyy0)+bay*(u2+bay2*uyy2)
      us=aay*(u1+aay2*uyy1)+bay*(u3+bay2*uyy3)
      uxxr=aay*(uxx0+aay2*uxxyy0)+bay*(uxx2+bay2*uxxyy2)
      uxxs=aay*(uxx1+aay2*uxxyy1)+bay*(uxx3+bay2*uxxyy3)
      fx=((us-ur)-aax23*uxxr+bax23*uxxs)/dxs
      fxx=-6.d0*(aax*uxxr+bax*uxxs)/dxs**2
      fxy=(aax23* (uxx2 - uxx0 + bay23*uxxyy2 - aay23*uxxyy0)
     $     - bax23* (uxx3 - uxx1 + bay23*uxxyy3 - aay23*uxxyy1)
     $     + ((u2 - u0 + bay23*uyy2 - aay23*uyy0)
     $     - (u3 - u1 + bay23*uyy3 - aay23*uyy1)))/dxs/dys
      if(x .gt. 0.d0)then
        fx=-fx
        fxy=-fxy
      endif
      if(y .gt. 0.d0)then
        fy=-fy
        fxy=-fxy
      endif
      return
 9000 z=bbkick1(x,y,sigx,sigy)
      fx=dble(z)
      fy=dimag(z)
      fu=twspu(x,y,sigx,sigy,1.d-8,1.d-11)
      fxx=0.d0
      fyy=0.d0
      fxy=0.d0
c      write(*,*)'twspfu ',x,y,sigx,sigy
      return
      end
      
      subroutine twspfuinit(u,uxx,uyy,uxxyy,urr,uxxrr,uyyrr,uxxyyrr)
      implicit none
      integer*4 nr,nx,ny,i,j,n
      real*8 xm,xstep,ym,ystep,x,y,r,twspu,
     $     rstep,rm,ystepr
      parameter (nr=20,nx=60,ny=60,xm=15.d0,ym=30.d0,rm=5.d0)
      parameter (xstep=xm/nx,rstep=rm/nr**2,ystep=ym/ny)
      complex*16 bbkick1
      real*8 u(0:nx,0:ny,0:nr),uxx(0:nx,0:ny,0:nr),
     $     uyy(0:nx,0:ny,0:nr),uxxyy(0:nx,0:ny,0:nr)
      real*8 urr(0:nx,0:ny,0:nr),uxxrr(0:nx,0:ny,0:nr),
     $     uyyrr(0:nx,0:ny,0:nr),uxxyyrr(0:nx,0:ny,0:nr)
      real*8 s(0:ny),dds(0:ny),work(0:ny)
      real*8 sxx(0:nr),ddsxx(0:nr)
      real*8 syy(0:nr),ddsyy(0:nr)
      real*8 sxxyy(0:nr),ddsxxyy(0:nr)
      do n=0,nr
        r=10.d0**(-n**2*rstep)
        ystepr=r*ystep
        do i=0,nx
          x=xstep*i
          do j=0,ny
            y=j*ystepr
            s(j)=twspu(x,y,1.d0,r,1.d-8,1.d-11)
            u(i,j,n)=s(j)
          enddo
          dds(0)=0.d0
          dds(ny)=-dimag(bbkick1(x,ny*ystepr,1.d0,r))*ystepr
          call spline1(ny+1,s,dds,work,1,1)
          do j=0,ny
            uyy(i,j,n)=dds(j)
          enddo
        enddo
        do j=0,ny
          dds(0)=0.d0
          dds(nx)=-dble(bbkick1(xm,j*ystepr,1.d0,r))*xstep
          call spline1(nx+1,u(0,j,n),uxx(0,j,n),work,1,1)
          s(j)=uxx(0,j,n)
        enddo
        dds(0)=0.d0
        call spline1(ny+1,s,dds,work,1,0)
        do j=0,ny
          uxxyy(0,j,n)=dds(j)
          s(j)=uxx(nx,j,n)
        enddo
        call spline1(ny+1,s,dds,work,0,0)
        do j=0,ny
          uxxyy(nx,j,n)=dds(j)
          call spline1(nx+1,uyy(0,j,n),uxxyy(0,j,n),work,2,2)
        enddo
      enddo
      do i=0,nx
        do j=0,ny
          do n=0,nr
            s(n)=u(i,j,n)
            sxx(n)=uxx(i,j,n)
            syy(n)=uyy(i,j,n)
            sxxyy(n)=uxxyy(i,j,n)
          enddo
          dds(0)=0.d0
          ddsxx(0)=0.d0
          ddsyy(0)=0.d0
          ddsxxyy(0)=0.d0
          call spline1(nr+1,s,dds,work,1,0)
          call spline1(nr+1,sxx,ddsxx,work,1,0)
          call spline1(nr+1,syy,ddsyy,work,1,0)
          call spline1(nr+1,sxxyy,ddsxxyy,work,1,0)
          do n=0,nr
            urr(i,j,n)=dds(n)
            uxxrr(i,j,n)=ddsxx(n)
            uyyrr(i,j,n)=ddsyy(n)
            uxxyyrr(i,j,n)=ddsxxyy(n)
          enddo
        enddo
      enddo
      return
      end

      logical*4 function wspaccheck() result(v)
      use tfstk
      use ffs
      implicit none
      v=wspac .and. (ifsize .eq. 0 .or. modesize .ne. 6)
      if(v)then
        write(*,*)'WSPAC without beam matrix!  ',
     $       'You need EMIT with CODPLOT.'
      endif
      return
      end
