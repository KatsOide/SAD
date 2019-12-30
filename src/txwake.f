      subroutine txwake(np,x,px,y,py,z,g,dv,sx,sy,sz,
     $     dx,dy,theta,nb,
     $     fw,lwl,wakel,lwt,waket,p0,h0,itab,izs,init)
      use tfstk, only: ktfenanq
      use ffs_flag,only:calpol
      implicit none
      integer*4 np,ns,lwl,lwt,i,n,k,l,m,nb
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     sx(np),sy(np),sz(np),
     $     wakel(2,lwl),waket(2,lwt),xs(np),ys(np),zs(np),
     $     dz,w,wx(np),wy(np),wz(np),ws(np),zk,dzk,
     $     pa,pb,h1,fw,dwx,dwy,dwz,p0,h0,fwp,dx,dy,theta,
     $     cost,sint,xi,pxi,pmin,zmin
      parameter (pmin=1.d-10,zmin=-1.d30)
      integer*4 itab(np),izs(np)
      logical*4 init
      if(lwl .eq. 0 .and. lwt .eq. 0)then
        return
      endif
      include 'inc/TENT.inc'
      if(init)then
        do i=1,np
          itab(i)=i
          if(ktfenanq(z(i)))then
            z(i)=zmin
          endif
        enddo
        call spsort(np,itab,z)
        call txwdefslice(np,z,nb,itab,izs)
      endif
      fwp=fw/p0
      ns=izs(np)
      do n=1,ns
        ws(n)=0.d0
        xs(n)=0.d0
        ys(n)=0.d0
        zs(n)=0.d0
        wx(n)=0.d0
        wy(n)=0.d0
        wz(n)=0.d0
      enddo
      do i=1,np
        n=izs(i)
        ws(n)=ws(n)+1.d0
        k=itab(i)
        xs(n)=xs(n)+x(k)
        ys(n)=ys(n)+y(k)
        zs(n)=zs(n)+z(k)
      enddo
      do n=1,ns
        if(ws(n) .ne. 0.d0)then
          zs(n)=zs(n)/ws(n)
        endif
      enddo
      do n=1,ns
        do m=n,ns
          dz=zs(m)-zs(n)
          do l=1,lwl
            if(wakel(1,l) .ge. dz)then
              if(l .eq. 1)then
                wz(n)=wz(n)+wakel(2,l)*ws(m)*.5d0
              else
                w=(dz-wakel(1,l-1))/(wakel(1,l)-wakel(1,l-1))
                wz(n)=wz(n)+(wakel(2,l-1)*(1.d0-w)+wakel(2,l)*w)*ws(m)
              endif
              exit
            endif
          enddo
          do l=1,lwt
            if(waket(1,l) .ge. dz)then
              if(l .eq. 1)then
                wx(n)=wx(n)+waket(2,l)*xs(m)*.5d0
                wy(n)=wy(n)+waket(2,l)*ys(m)*.5d0
              else
                w=(dz-waket(1,l-1))/(waket(1,l)-waket(1,l-1))
                wx(n)=wx(n)+(waket(2,l-1)*(1.d0-w)+waket(2,l)*w)*xs(m)
                wy(n)=wy(n)+(waket(2,l-1)*(1.d0-w)+waket(2,l)*w)*ys(m)
              endif
              exit
            endif
          enddo
        enddo
      enddo
      do i=1,np
        n=izs(i)
        k=itab(i)
        zk=z(k)
        dzk=zk-zs(n)
        if(dzk .lt. 0.d0)then
          if(n .eq. 1)then
            dwx=wx(1)
            dwy=wy(1)
            dwz=wz(1)
            go to 10
          else
            n=n-1
            dzk=zk-zs(n)
          endif
        endif
        if(n .eq. ns)then
          dwx=wx(n)
          dwy=wy(n)
          dwz=wz(n)
        else
          w=dzk/(zs(n+1)-zs(n))
          dwx=(1.d0-w)*wx(n)+w*wx(n+1)
          dwy=(1.d0-w)*wy(n)+w*wy(n+1)
          dwz=(1.d0-w)*wz(n)+w*wz(n+1)
        endif
 10     pa=1.d0+g(k)
        g(k)=max(g(k)-fwp*dwz,pmin-1.d0)
        pb=1.d0+g(k)
        px(k)=(pa*px(k)+fwp*dwx)/pb
        py(k)=(pa*py(k)+fwp*dwy)/pb
        h1=sqrt(1.d0+(pa*p0)**2)
        dv(k)=-g(k)*pa/h1/(h1+pa*h0)
c        if(mod(i,200) .eq. 0)then
c          write(*,*)'txwake ',i,n,dzk,dwz,wz(n),fw
c        endif
      enddo
      include 'inc/TEXIT.inc'
      return
      end

      subroutine txwdefslice(np,z,nb,itab,izs)
      implicit none
      integer*4 np,nb,itab(np),izs(np),mb(0:nb),
     $     i,m,m0,nm,ks,m1,j
      real*8 z(np),za,sb
      za=z(itab(np))-z(itab(1))
      if(nb .gt. 1)then
        sb=za/(nb-1)*.5d0
      else
        sb=za
      endif
      do m=0,nb
        mb(m)=0
      enddo
      m=1
      izs(1)=1
      do i=2,np
        if(z(itab(i))-z(itab(i-1)) .gt. sb)then
          mb(m)=i-1
          m=m+1
          if(m .gt. nb)then
            do j=i,np
              izs(j)=m
            enddo
            mb(m)=np
            exit
          endif
        endif
        izs(i)=m
      enddo
      if(m .eq. 1)then
        mb(1)=np
      endif
      do m=nb,1,-1
        if(mb(m) .ne. 0)then
          nm=mb(m)-mb(m-1)
          m1=int(sqrt(dble(nm)))
          mb(m)=(nm+m1-1)/m1
        endif
      enddo
      m0=1
      nm=mb(1)-1
      ks=1
      do i=2,np
        if(izs(i) .ne. m0)then
          m0=m0+1
          ks=ks+1
          izs(i)=ks
          nm=mb(m0)-1
        else
          izs(i)=ks
          nm=nm-1
          if(nm .le. 0)then
            nm=mb(m0)
            ks=ks+1
          endif
        endif
c        write(*,*)'txwds ',i,izs(i)
      enddo
      return
      end
