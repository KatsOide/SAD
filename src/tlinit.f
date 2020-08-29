      subroutine tlinit(np,h0,geo)
      use tfmem, only:ktaloc
      implicit none
      integer*4 np
      real*8 h0,geo(3,4)
      integer*8 katbl
      integer*4 np0,ltbl,lpoint
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
      np0=np
      ltbl=int(h0*10.d0)+1
      katbl=ktaloc(np*4*ltbl)
      lpoint=0
      gx=geo(1,4)
      gy=geo(2,4)
      gz=geo(3,4)
      gt=0.d0
      geo0(:,3)=geo(:,3)
      geo0(1,4)=gx
      geo0(2,4)=gy
      geo0(3,4)=gz
      t0=gt
      return
      end

      subroutine tlstore(np,x,y,z,dv,theta,dl,dphi,rho0,v0,dvfs,keep)
      use tfstk
      implicit none
      integer*8 katbl,new,kai
      integer*4 np,np0,ltbl,lpoint,m,i
      real*8 x(np),y(np),z(np),dv(np),theta,dphi,rho0,v0,
     $     x1,x2,x3,y1,y2,y3,sp0,cp0,r1,r2,
     $     geoi(3,3),cost,sint,dl
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt,dvfs
      logical*4 keep
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
c     begin initialize for preventing compiler warning
      cost=1.d0
      sint=0.d0
      y1=0.d0
      y2=0.d0
      y3=0.d0
c     end   initialize for preventing compiler warning
c      write(*,*)'tlstore ',lpoint,dphi,dl,theta,rho0,v0
      if(dphi .eq. 0.d0)then
        gx=geo0(1,4)+dl*geo0(1,3)
        gy=geo0(2,4)+dl*geo0(2,3)
        gz=geo0(3,4)+dl*geo0(3,3)
        geoi=geo0(:,1:3)
        gt=t0+dl/v0
      else
        if(theta .ne. 0.d0)then
          cost=cos(theta)
          sint=sin(theta)
          x1= cost*geo0(1,1)-sint*geo0(1,2)
          x2= cost*geo0(2,1)-sint*geo0(2,2)
          x3= cost*geo0(3,1)-sint*geo0(3,2)
          y1= sint*geo0(1,1)+cost*geo0(1,2)
          y2= sint*geo0(2,1)+cost*geo0(2,2)
          y3= sint*geo0(3,1)+cost*geo0(3,2)
        else
          x1= geo0(1,1)
          x2= geo0(2,1)
          x3= geo0(3,1)
        endif
        sp0=sin(dphi)
        cp0=cos(dphi)
        r1=rho0*sp0
        if(cp0 .ge. 0.d0)then
          r2=rho0*sp0**2/(1.d0+cp0)
        else
          r2=rho0*(1.d0-cp0)
        endif
        gx=geo0(1,4)+(r1*geo0(1,3)-r2*x1)
        gy=geo0(2,4)+(r1*geo0(2,3)-r2*x2)
        gz=geo0(3,4)+(r1*geo0(3,3)-r2*x3)
        geoi(1,1)= cp0*x1+sp0*geo0(1,3)
        geoi(1,3)=-sp0*x1+cp0*geo0(1,3)
        geoi(2,1)= cp0*x2+sp0*geo0(2,3)
        geoi(2,3)=-sp0*x2+cp0*geo0(2,3)
        geoi(3,1)= cp0*x3+sp0*geo0(3,3)
        geoi(3,3)=-sp0*x3+cp0*geo0(3,3)
        geoi(1,2)=geo0(1,2)
        geoi(2,2)=geo0(2,2)
        geoi(3,2)=geo0(3,2)
        gt=t0+abs(rho0*dphi)/v0
      endif
      lpoint=lpoint+1
      if(lpoint .gt. ltbl)then
c        write(*,*)'tlspect ',ltbl,np0,ltbl*8*np0
        new=ktaloc(ltbl*8*np0)
        if(new .gt. 0)then
          call tmov(rlist(katbl),rlist(new),ltbl*np0)
          call tmov(rlist(katbl+ltbl),rlist(new+ltbl*2),ltbl*np0)
          call tmov(rlist(katbl+ltbl*2),rlist(new+ltbl*4),ltbl*np0)
          call tmov(rlist(katbl+ltbl*3),rlist(new+ltbl*6),ltbl*np0)
          ltbl=ltbl*2
          call tfree(katbl)
          katbl=new
        else
          write(*,*)'Too long ltbl',ltbl
          call abort
        endif
      endif
      m=ltbl*4
      do i=1,np
        kai=(i-1)*4*ltbl+katbl+lpoint-1
        rlist(kai       )=gt-z(i)/(v0*(1.d0-dv(i)+dvfs))
        rlist(kai+ltbl  )=gx+x(i)*geoi(1,1)+y(i)*geoi(1,2)
        rlist(kai+ltbl*2)=gy+x(i)*geoi(2,1)+y(i)*geoi(2,2)
        rlist(kai+ltbl*3)=gz+x(i)*geoi(3,1)+y(i)*geoi(3,2)
      enddo
      if(keep)then
        geo0(:,3)=geoi(:,3)
        if(dphi .ne. 0 .and. theta .ne. 0.d0)then
          geo0(1,2)=-sint*geo0(1,1)+cost*y1
          geo0(1,1)= cost*geo0(1,1)+sint*y1
          geo0(2,2)=-sint*geo0(2,1)+cost*y2
          geo0(2,1)= cost*geo0(2,1)+sint*y2
          geo0(3,2)=-sint*geo0(3,1)+cost*y3
          geo0(3,1)= cost*geo0(3,1)+sint*y3
        endif
        geo0(1,4)=gx
        geo0(2,4)=gy
        geo0(3,4)=gz
        t0=gt
      endif
      return
      end

      subroutine tlresult(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,kli
      type (sad_rlist), pointer :: klrij
      integer*8 katbl,kati
      integer*4 np0,ltbl,lpoint,i,j
      real*8 geo0(3,4),tax,tay,taz,t0,gx,gy,gz,gt
      common /radlcomr/geo0,t0,tax,tay,taz,gx,gy,gz,gt
      common /radlcomi/katbl,np0,ltbl,lpoint
c      write(*,*)ltbl,lpoint,katbl,ilist(1,katbl-1)
      kx=kxadaloc(-1,np0,klx)
      do i=1,np0
        kati=katbl+ltbl*4*(i-1)
        klx%dbody(i)=kxadaloc(0,4,kli)
        do j=1,4
          kli%dbody(j)=kxavaloc(0,lpoint,klrij)
          call tmov(rlist(kati+(j-1)*ltbl),klrij%rbody(1),lpoint)
        enddo
      enddo
      call tfree(katbl)
      return
      end

      subroutine tlfield(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_rlist), pointer :: klt
      integer*8 kx, ke
      integer*4 isp1,irtc,n,m,itfmessage
      real*8 tx,ty,tz
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      call tfmsize(ktastk(isp1+1),n,m,irtc)
      if(irtc .ne. 0)then
        return
      elseif(n .ne. 4)then
        irtc=itfmessage(9,'General::wrongleng','"#1","4"')
        return
      endif
      if(.not. tfreallistq(ktastk(isp),klt))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of Reals for #2"')
        return
      endif
      if(klt%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongleng','"#2","3"')
        return
      endif
      tx=klt%rbody(1)
      ty=klt%rbody(2)
      tz=klt%rbody(3)
      call tlcalc(ktastk(isp1+1),ke,m,tx,ty,tz)
      kx=ktflist+ke
      return
      end

      subroutine tlcalc(ka,ke,lpoint,tx,ty,tz)
      use tfstk
      use tmacro
      implicit none
      integer*8 ka,ke
      integer*4 lpoint,i,j
      real*8 xtbl(lpoint,4),f(lpoint,10),
     $     tx,ty,tz,dxi,dyi,dzi,ri,r0,dri,coeff,dt,f1,z0
      call tfl2m(klist(ktfaddr(ka)-3),xtbl,lpoint,4,.true.)
      z0=4.d-7*2.d0*asin(1.d0)*c
      r0=sqrt(tx**2+ty**2+tz**2)
      do i=1,lpoint
        dxi=tx-xtbl(i,2)
        dyi=ty-xtbl(i,3)
        dzi=tz-xtbl(i,4)
        ri=sqrt(dxi**2+dyi**2+dzi**2)
        dri=-(xtbl(i,2)*(tx+dxi)+xtbl(i,3)*(ty+dyi)
     $       +xtbl(i,4)*(tz+dzi))/(r0+ri)
        xtbl(i,1)=xtbl(i,1)+dri/c
        xtbl(i,2)=dxi/ri
        xtbl(i,3)=dyi/ri
        xtbl(i,4)=dzi/ri
        f(i,8)=ri/c
      enddo
      do i=2,4
        call spline(lpoint,xtbl(1,1),xtbl(1,i),f(1,i),f(1,1),0,0.d0)
      enddo
      do i=5,7
        do j=1,lpoint
          f(j,10)=xtbl(j,i-3)/f(j,8)**2
        enddo
        call spline(lpoint,xtbl(1,1),f(1,10),f(1,i),f(1,1),0,0.d0)
        f1=f(lpoint-1,i)
        do j=1,lpoint-1
          dt=xtbl(j+1,1)-xtbl(j,1)
          f(j,i)=f(j,10)+f(j,8)*
     $         ((f(j+1,10)-f(j,10))/dt-dt*(2.d0*f(j,i)+f(j+1,i)))
        enddo
        dt=xtbl(lpoint,1)-xtbl(lpoint-1,1)
        f(lpoint,i)=f(lpoint,10)+f(lpoint,8)*
     $       ((f(lpoint,10)-f(lpoint-1,10))/dt
     $       +dt*(f1+2.d0*f(lpoint,i)))
      enddo
      coeff=1.d-7*e*charge
      do i=1,lpoint
        f(i,1)=xtbl(i,1)
        f(i,2)=coeff*(6.d0*f(i,2)+f(i,5))
        f(i,3)=coeff*(6.d0*f(i,3)+f(i,6))
        f(i,4)=coeff*(6.d0*f(i,4)+f(i,7))
        f(i,5) =(xtbl(i,3)*f(i,4)-xtbl(i,4)*f(i,3))/z0
        f(i,6) =(xtbl(i,4)*f(i,2)-xtbl(i,2)*f(i,4))/z0
        f(i,7) =(xtbl(i,2)*f(i,3)-xtbl(i,3)*f(i,2))/z0
        f(i,8) =f(i,3)*f(i,7)-f(i,4)*f(i,6)
        f(i,9) =f(i,4)*f(i,5)-f(i,2)*f(i,7)
        f(i,10)=f(i,2)*f(i,6)-f(i,3)*f(i,5)
      enddo
      ke=ktfaddr(kxm2l(f,lpoint,10,lpoint,.true.))
      return
      end

      subroutine tlspect(isp1,kx,nparallel,irtc)
      use tfstk
      use macphys
      implicit none
      type (sad_dlist), pointer :: klp,klx
      type (sad_rlist), pointer :: kl1,kl2,klt
      real*8 speedoflight
      parameter (speedoflight=cveloc)
      integer*8 kx,k1,k2,kac,kas,kal,kack,kad
      integer*4 isp1,irtc,m,na,i,k,
     $     na1,nastep,nparallel,ipr,fork_worker,getpid,j,
     $     waitpid,irw,naa,isw,ipid,itfmessage
      real*8 al0,al1,dal,w,c,s,a0,a1,x,ak,ak0,ak1,dak
      integer*4 ichpid(100)
      character*64 fn
c     begin initialize for preventing compiler warning
      ipid=0
      naa=0
      kal=0
      kas=0
c     end   initialize for preventing compiler warning
      irtc=0
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. tflistq(ktastk(isp1+1),klp))then
        irtc=itfmessage(9,'General::wrongtype','"List for #1"')
        return
      endif
      if(klp%nl .ne. 2)then
        irtc=itfmessage(9,'General::wrongleng','"#1","2"')
        return
      endif
      k1=klp%dbody(1)%k
      k2=klp%dbody(2)%k
      if(.not. tfreallistq(k1,kl1) .or. .not. tfreallistq(k2,kl2))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of List of Reals for #1"')
        return
      endif
      m=kl1%nl
      if(m .ne. kl2%nl)then
        irtc=itfmessage(9,'General::equalleng','"#1[[1]] and #[[2]]"')
        return
      endif
      if(.not. tfreallistq(ktastk(isp),klt))then
        irtc=itfmessage(9,'General::wrongtype','"List of Reals for #2"')
        return
      endif
      if(klt%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongleng','"#","3"')
        return
      endif
      al0=klt%rbody(1)
      al1=klt%rbody(2)
      dal=klt%rbody(3)
      na=int((al1-al0)/dal+2.5d0)
      ak0=2.d0*pi/al1
      ak1=2.d0*pi/al0
      dak=(ak1-ak0)/(na-1)
      na1=1
      if(nparallel .gt. 1 .and. na .ge. nparallel)then
        ipid=getpid()
        nastep=nparallel
        ipr=1
        do while(na1 .lt. nparallel .and. ipr .gt. 0)
          ipr=fork_worker()
          if(ipr .gt. 0)then
            ichpid(na1)=ipr
            na1=na1+1
          endif
        enddo
      else
        ipr=1
        nastep=1
      endif
      if(ipr .eq. 0)then
        naa=(na-na1)/nastep+1
        kac=ktaloc(naa*3)
      else
        kac=ktavaloc(0,na)
        kas=ktavaloc(0,na)
        kal=ktavaloc(0,na)
      endif
      do k=na1,na,nastep
        ak=ak0+dak*(k-1)
        w=ak*speedoflight
        c=0.d0
        s=0.d0
        do i=2,m-1
          x=w*kl1%rbody(i)
          a0=(kl2%rbody(i)-kl2%rbody(i-1))/
     $         (kl1%rbody(i)-kl1%rbody(i-1))
          a1=(kl2%rbody(i+1)-kl2%rbody(i))/
     $         (kl1%rbody(i+1)-kl1%rbody(i))
          c=c+cos(x)*(a0-a1)
          s=s+sin(x)*(a0-a1)
        enddo
        a0=(kl2%rbody(m)-kl2%rbody(m-1))/
     $       (kl1%rbody(m)-kl1%rbody(m-1))
        x=w*kl1%rbody(m)
        c=c/w+cos(x)*a0/w+kl2%rbody(m)*sin(x)
        s=s/w+sin(x)*a0/w-kl2%rbody(m)*cos(x)
        a1=(kl2%rbody(2)-kl2%rbody(1))/
     $       (kl1%rbody(2)-kl1%rbody(1))
        x=w*kl1%rbody(1)
        c=(c-cos(x)*a1/w-kl2%rbody(1)*sin(x))/w
        s=(s-sin(x)*a1/w+kl2%rbody(1)*cos(x))/w
        if(ipr .eq. 0)then
          kack=kac+(k-na1)/nastep*3
          rlist(kack  )=ak
          rlist(kack+1)=c
          rlist(kack+2)=s
        else
          rlist(kal+k)=ak
          rlist(kac+k)=c
          rlist(kas+k)=s
        endif
      enddo
      if(nparallel .gt. 1 .and. na .ge. nparallel)then
        if(ipr .eq. 0)then
          write(fn,'(a,i12.12)')'/tmp/',(ipid*nparallel)+na1
          open(999,file=fn,form='unformatted',status='unknown')
          write(999)(rlist(kac+i),i=0,naa*3-1)
          stop
        else
          kad=ktaloc(((na-1)/nastep+1)*3)
          do j=1,nparallel-1
 3001       irw=waitpid(-1,isw)
            do k=1,nparallel-1
              if(irw .eq. ichpid(k))then
                na1=k
                ichpid(k)=0
                go to 3010
              endif
            enddo
            write(*,*)'RadiationSpectrum Strange child: ',irw,isw
            go to 3001
 3010       write(fn,'(a,i12.12)')'/tmp/',(ipid*nparallel)+na1
            open(999,file=fn,form='unformatted',status='old')
            naa=(na-na1)/nastep+1
            read(999)(rlist(kad+i),i=0,naa*3-1)
            close(999,status='DELETE')
            do k=na1,na,nastep
              kack=kad+(k-na1)/nastep*3
              rlist(kal+k)=rlist(kack  )
              rlist(kac+k)=rlist(kack+1)
              rlist(kas+k)=rlist(kack+2)
            enddo
          enddo
          call tfree(kad)
        endif
      endif
      kx=ktflist+ktadaloc(-1,3,klx)
      klx%dbody(1)%k=ktflist+kal
      klx%dbody(2)%k=ktflist+kac
      klx%dbody(3)%k=ktflist+kas
      return
      end

c$$$      subroutine tsynchrad(p,alr,bxa,bya,
c$$$     $     dprad,dpradx,dprady,
c$$$     $     mp,l,al,s0,phi,theta,xi,yi,pxi,pyi)
c$$$      use tfstk
c$$$      use ffs
c$$$      use ffs_pointer
c$$$      use tffitcode
c$$$      implicit none
c$$$      integer*4 irtn,ng,l,mp
c$$$      real*8 p,alr,bxa,bya,dprad,dpradx,dprady,tran,
c$$$     $     b0,e0,eg,thx,thy,xi30,thu,thv,
c$$$     $     c1,c2,s1,s2,xi3,xi2,xi1,s0,phi1,
c$$$     $     thx1,thy1,ds,thu1,xi,yi,pxi,pyi,al,al1,
c$$$     $     phi,theta,gx,gy,gz,dpgx,dpgy,dpgz
c$$$      e0=amass*sqrt((p*p0)**2+1.d0)
c$$$      b0=sqrt(bxa**2+bya**2)
c$$$      call SYNRADCL(E0,B0,alr,2,NG,EG,thu,thv,XI30,XI2,IRTN)
c$$$      if(irtn .gt. 100 .or. ng .eq. 0)then
c$$$        dprad=0.d0
c$$$        dpradx=0.d0
c$$$        dprady=0.d0
c$$$      else
c$$$        dprad=eg/amass/p0
c$$$        c1=bya/b0
c$$$        s1=bxa/b0
c$$$        thx=thu*c1+thv*s1
c$$$        thy=-thu*s1+thv*c1
c$$$        dpradx=thx*dprad
c$$$        dprady=thy*dprad
c$$$        if(photons)then
c$$$          c2=(c1-s1)*(c1+s1)
c$$$          s2=2.d0*c1*s1
c$$$          xi3=c2*xi30
c$$$          xi1=-s2*xi30
c$$$          ds=tran()*alr
c$$$          thu1=thu-(s0+ds)*b0/brho
c$$$          thx1=thu1*c1+thv*s1
c$$$          thy1=-thu1*s1+thv*c1
c$$$          al1=al+ds+s0
c$$$          if(al .eq. 0.d0)then
c$$$            phi1=0.d0
c$$$          else
c$$$            phi1=phi*al1/al
c$$$          endif
c$$$          call tphotonconv(al1,phi1,theta,
c$$$     $         geo(:,:,l),xi,yi,
c$$$     $         dprad,pxi-thx1,pyi-thy1,
c$$$     $         gx,gy,gz,dpgx,dpgy,dpgz,xi1,xi3)
c$$$          call tphotonstore(mp,l,gx,gy,gz,
c$$$     $         dpgx*p0,dpgy*p0,dpgz*p0,
c$$$     $         xi1,xi2,xi3)
c$$$        endif
c$$$      endif
c$$$      return
c$$$      end
