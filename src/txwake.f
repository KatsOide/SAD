      module wakez
        real*8 ,parameter :: dzlim=1.d-5
        real*8 ::dzwr=0.d0

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
        real*8 ,intent(in):: dx,dy,theta,fw,p0,h0
        logical*4 ,intent(in):: init
        integer*4 ns,i,n,k,l,m,l1,lwl,lwt,ns1
        real*8 ,dimension(:), allocatable:: xs,ys,zs,wx,wy,wz,ws
        real*8 dz,w,zk,dzk,pa,pb,h1,dwx,dwy,dwz,fwp,
     $       cost,sint,xi,pxi,pmin,zmin,zw0l,zw0t
        parameter (pmin=1.d-10,zmin=-1.d30)
        call twxiwp(nwak)
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
c        write(*,'(a,l2,i5,1p8g14.6)')'txwake ',init,nwak,p0,z(itab(np)),z(itab(1)),z(itab(np))-z(itab(1))-dzwr,dzlim*dzwr
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
        ns1=izs(np)
        ns=ns1+1
        allocate(ws(0:ns),xs(0:ns),ys(0:ns),zs(0:ns),wx(0:ns),wy(0:ns),wz(0:ns))
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
        do n=1,ns1
          if(ws(n) /= 0.d0)then
            zs(n)=zs(n)/ws(n)
          endif
        enddo
        zs(0)=z(itab(1))-dzwr/np
        zs(ns)=z(itab(np))+dzwr/np
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
        do n=0,ns
          do m=0,ns
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
            n=n-1
            dzk=zk-zs(n)
          endif
c          elseif(n == ns)then
c            n=n-1
c            dzk=zk-zs(n)
c          endif
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

        logical*4 function ewak(l,nextwake,lele,cmp,nwak,nwak1) result(a)
        use sad_main
        use maccode
        use kyparam
        implicit none
        integer*4 ,intent(in):: l,nextwake,lele,nwak
        type (sad_comp),intent(in):: cmp
        integer*4 ,intent(out):: nwak1
        if(l /= nextwake .or. nwak == 0)then
          a=.false.
          nwak1=0
        elseif(lele == icCAVI)then
          a=.false.
          nwak1=nwak
        elseif(lele == icMULT)then
          nwak1=nwak
          if(cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT) == 0.d0)then
            a=.true.
          else
            a=.false.
          endif
        else
          a=.true.
          nwak1=nwak
        endif
        return
        end function

      end module
