      module wakez
        real*8 dzwr
        real*8 ,parameter :: dzlim=1.d-5

        contains
        subroutine txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $       dx,dy,theta,nb,
     $       fw,nwak,p0,h0,init)
        use tfstk, only: ktfenanq
        use ffs_flag,only:calpol,twake,lwake
        use ffs_wake
        implicit none
        integer*4 ,intent(in):: np,nb,nwak
        real*8 ,intent(inout):: x(np),px(np),y(np),py(np),z(np),g(np),
     $       dv(np),sx(np),sy(np),sz(np)
        real*8 ,intent(in):: dx,dy,theta
        logical*4 ,intent(in):: init
        integer*4 ns,i,n,k,l,m,l1,lwl,lwt
        real*8 ,dimension(:), allocatable:: xs,ys,zs,wx,wy,wz,ws
        real*8 dz,w,zk,dzk,pa,pb,h1,fw,dwx,dwy,dwz,p0,h0,fwp,
     $       cost,sint,xi,pxi,pmin,zmin,zw0l,zw0t
        parameter (pmin=1.d-10,zmin=-1.d30)
        call twxiwp(nwp)
        if(lwake)then
          lwl=lwl0
        else
          lwl=0
        endif
        if(twake)then
          lwt=lwt0
        else
          lwt=0
        endif
        if(lwl == 0 .and. lwt == 0)then
          return
        endif
        include 'inc/TENT.inc'
        if(init)then
          call twxalloc(np)
        endif
        if(init .or.
     $       abs(z(itab(np))-z(itab(1))-dzwr) >= dzlim*dzwr)then
          do i=1,np
            itab(i)=i
            if(ktfenanq(z(i)))then
              z(i)=zmin
            endif
          enddo
          call spsort(np,itab,z)
          call txwdefslice(np,z,nb,itab,izs)
          dzwr=z(itab(np))-z(itab(1))
        endif
        fwp=fw/p0
        ns=izs(np)
        allocate(ws(ns),xs(ns),ys(ns),zs(ns),wx(ns),wy(ns),wz(ns))
        ws=0.d0
        xs=0.d0
        ys=0.d0
        zs=0.d0
        wx=0.d0
        wy=0.d0
        wz=0.d0
        do i=1,np
          n=izs(i)
          ws(n)=ws(n)+1.d0
          k=itab(i)
          xs(n)=xs(n)+x(k)
          ys(n)=ys(n)+y(k)
          zs(n)=zs(n)+z(k)
        enddo
        do n=1,ns
          if(ws(n) /= 0.d0)then
            zs(n)=zs(n)/ws(n)
          endif
        enddo
        if(lwl > 0)then
          zw0l=min(0.d0,wakel(1,1))
        else
          zw0l=1.d300
        endif
        if(lwt > 0)then
          zw0t=min(0.d0,waket(1,1))
        else
          zw0t=1.d300
        endif
        do n=1,ns
          do m=1,ns
            dz=zs(m)-zs(n)
            if(dz >= zw0l)then
              do l=1,lwl
                if(wakel(1,l) >= dz)then
                  if(n /= m .or. zw0l < 0.d0)then
                    l1=min(max(l-1,1),lwl-1)
                    w=(dz-wakel(1,l1))/(wakel(1,l1+1)-wakel(1,l1))
                    wz(n)=wz(n)+(wakel(2,l1)*(1.d0-w)+wakel(2,l1+1)*w)*ws(m)
                  else
                    wz(n)=wz(n)+wakel(2,1)*ws(n)*.5d0
                  endif
                  exit
                endif
              enddo
            endif
            if(dz >= zw0t)then
              do l=1,lwt
                if(waket(1,l) >= dz)then
                  if(n /= m .or. zw0t < 0.d0)then
                    l1=min(max(l-1,1),lwt-1)
                    w=(dz-waket(1,l1))/(waket(1,l1+1)-waket(1,l1))
                    wx(n)=wx(n)+(waket(2,l1)*(1.d0-w)+waket(2,l1+1)*w)*xs(m)
                    wy(n)=wy(n)+(waket(2,l1)*(1.d0-w)+waket(2,l1+1)*w)*ys(m)
                  else
                    wx(n)=wx(n)+waket(2,1)*xs(n)*.5d0
                    wy(n)=wy(n)+waket(2,1)*ys(n)*.5d0
                  endif
                  exit
                endif
              enddo
            endif
          enddo
        enddo
        do i=1,np
          n=izs(i)
          k=itab(i)
          zk=z(k)
          dzk=zk-zs(n)
          if(dzk < 0.d0)then
            if(n /= 1)then
              n=n-1
              dzk=zk-zs(n)
            endif
          elseif(n == ns)then
            n=n-1
            dzk=zk-zs(n)
          endif
          w=dzk/(zs(n+1)-zs(n))
          dwx=(1.d0-w)*wx(n)+w*wx(n+1)
          dwy=(1.d0-w)*wy(n)+w*wy(n+1)
          dwz=(1.d0-w)*wz(n)+w*wz(n+1)
          pa=1.d0+g(k)
          g(k)=max(g(k)-fwp*dwz,pmin-1.d0)
          pb=1.d0+g(k)
          px(k)=(pa*px(k)+fwp*dwx)/pb
          py(k)=(pa*py(k)+fwp*dwy)/pb
          h1=sqrt(1.d0+(pa*p0)**2)
          dv(k)=-g(k)*pa/h1/(h1+pa*h0)
        enddo
        include 'inc/TEXIT.inc'
        return
        end subroutine 

        subroutine txwdefslice(np,z,nb,itab,izs)
        implicit none
        integer*4 ,intent(inout):: itab(np),izs(np)
        integer*4 ,intent(in):: np,nb
        real*8 ,intent(in):: z(np)
        integer*4 ,allocatable,dimension(:)::mb
        integer*4 i,m,m0,nm,ks,m1,j
        real*8 za,sb
        za=z(itab(np))-z(itab(1))
        if(nb > 1)then
          sb=za/(nb-1)*.5d0
        else
          sb=za
        endif
        allocate(mb(0:nb))
        mb=0
        m=1
        izs(1)=1
        do i=2,np
          if(z(itab(i))-z(itab(i-1)) > sb)then
            mb(m)=i-1
            m=m+1
            if(m > nb)then
              do j=i,np
                izs(j)=m
              enddo
              mb(m)=np
              exit
            endif
          endif
          izs(i)=m
        enddo
        if(m == 1)then
          mb(1)=np
        endif
        do m=nb,1,-1
          if(mb(m) /= 0)then
            nm=mb(m)-mb(m-1)
            m1=int(sqrt(dble(nm)))
            mb(m)=(nm+m1-1)/m1
          endif
        enddo
        m0=1
        nm=mb(1)-1
        ks=1
        do i=2,np
          if(izs(i) /= m0)then
            m0=m0+1
            ks=ks+1
            izs(i)=ks
            nm=mb(m0)-1
          else
            izs(i)=ks
            nm=nm-1
            if(nm <= 0)then
              nm=mb(m0)
              ks=ks+1
            endif
          endif
        enddo
        return
        end subroutine

      end module
