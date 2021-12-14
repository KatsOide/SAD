      subroutine tphplt(np,zbuf,lsp,x0,y0,xj,xnu,xphi,zwork,
     1                  nturn,np0,kptbl,lfnplt)
      implicit none
      integer*4 np,lsp,i,k,n,np0,kptbl(np0,6),lfnplt,nturn
      integer*4 iyaxis(6),ixaxis(6),ix,iy,iyk
      real*8 zbuf(2,lsp,3,np0),x0(3,np0),y0(3,np0),xj(3,np0),xnu(3,np0)
      real*8 xphi(3,np0)
      real*8 zwork(2,lsp),phi1,phi2
c      real*8 zra,zia,phi,c,s,wtune
      real*8 xmax,xmin,ymax,ymin
      real*8 xsp,ysp,x,y,xa,ya,xa1,ya1,xj1,xj2,xxa,yya,xa2,ya2
      character*3 titlel(6),titleb(6)
      data titlel/'PX ','PY ','PZ ','dJX','dJY','dJY'/
      data titleb/'X  ','Y  ','Z  ','dJZ','dJZ','dJX'/
      data iyaxis/   1,   2,   3,  1,  2,  2/
      data ixaxis/   1,   2,   3,  3,  3,  1/
      do 10 k=1,3
        xmax=-1.d30
        xmin=1.d30
        ymax=-1.d30
        ymin=1.d30
        ix=ixaxis(k)
        iy=iyaxis(k)
        iyk=2-(k-1)/3
        do 110 i=1,np0
          if(kptbl(i,1) .le. np)then
            do 120 n=1,nturn
              x=zbuf(1,n,ix,i)
              y=zbuf(iyk,n,iy,i)
              xmax=max(xmax,x)
              xmin=min(xmin,x)
              ymax=max(ymax,y)
              ymin=min(ymin,y)
120         continue
          endif
110     continue
        xsp=xmax-xmin
        ysp=ymax-ymin
        if(k .lt. 4)then
          xsp=max(xsp,ysp)*1.3d0
          ysp=xsp
        else
          xsp=xsp*1.3d0
          ysp=ysp*1.3d0
        endif
        xa=xmax+xmin
        xmax=(xa+xsp)*.5d0
        xmin=(xa-xsp)*.5d0
        ya=ymax+ymin
        ymax=(ya+ysp)*.5d0
        ymin=(ya-ysp)*.5d0
        write(lfnplt,*)'NEWFRAME;SET FONT DUPLEX;SET SYMBOL .M SIZE 1.4'
        write(lfnplt,*)'SET WINDOW X 4 11 Y 2 9'
        write(lfnplt,*)'TITLE LEFT   ''',titlel(k),''''
        write(lfnplt,*)'TITLE BOTTOM ''',titleb(k),''''
        write(lfnplt,9001)sngl(xmin),sngl(xmax),sngl(ymin),sngl(ymax)
9001    format(' SET LIMIT X ',1p,2g12.4,' Y ',2g12.4)
        call tdinit(lfnplt,'PLOT',' ')
        do 130 i=1,np0
          do 140 n=1,nturn
            x=zbuf(1,n,ix,i)
            y=zbuf(iyk,n,iy,i)
            if(x .ne. 0.d0 .or. y .ne. 0.d0)then
              call tdput(x,y)
            endif
 140      continue
 130    continue
        call tdterm
10    continue
      do 1110 i=1,np0
        if(kptbl(i,1) .le. np)then
          do 1120 k=1,3
c           xnu(k,i)=wtune(zbuf(1,1,k,i),zwork,lsp)*
c    1               2.d0*3.14159265358939324d0
            xa=0.d0
            xxa=0.d0
            ya=0.d0
            yya=0.d0
            do 1130 n=1,nturn
              x=zbuf(1,n,k,i)
              y=zbuf(2,n,k,i)
              xa=xa+x
              xxa=xxa+x**2
              ya=ya+y
              yya=yya+y**2
1130        continue
            x0(k,i)=xa/nturn
            y0(k,i)=ya/nturn
            xj(k,i)=sqrt(xxa/nturn-x0(k,i)**2+yya/nturn-y0(k,i)**2)
c           zra=0.d0
c           zia=0.d0
c           do 1140 n=1,nturn
c             x=zbuf(1,n,k,i)-x0(k,i)
c             y=zbuf(2,n,k,i)-y0(k,i)
c             phi=(n-1)*xnu(k,i)
c             c=cos(phi)
c             s=sin(phi)
c             zra=zra+x*c-y*s
c             zia=zia+x*c+y*s
c140        continue
c           xphi(k,i)=atan2(zia,zra)
c           write(*,'(2i5,1p5g12.4)')
c    1                k,i,x0(k,i),y0(k,i),xj(k,i),xphi(k,i),xnu(k,i)
1120      continue
        endif
1110  continue
      do 1010 k=1,3
        xmax=-1.d30
        xmin=1.d30
        ymax=-1.d30
        ymin=1.d30
        ix=ixaxis(k+3)
        iy=iyaxis(k+3)
        do 1020 i=1,np0
          if(kptbl(i,1) .le. np)then
            xa1=x0(ix,i)
            ya1=y0(ix,i)
            xj1=xj(ix,i)
            phi1=xphi(ix,i)
            xa2=x0(iy,i)
            ya2=y0(iy,i)
            xj2=xj(iy,i)
            phi2=xphi(iy,i)
            do 1030 n=1,nturn
c             x=(zbuf(1,n,ix,i)-xa1)-xj1*cos(-(n-1)*xnu(ix,i)+phi1)
c             y=(zbuf(1,n,iy,i)-xa2)-xj2*cos(-(n-1)*xnu(iy,i)+phi2)
              x=(zbuf(1,n,ix,i)-xa1)**2+(zbuf(2,n,ix,i)-ya1)**2-xj1**2
              y=(zbuf(1,n,iy,i)-xa2)**2+(zbuf(2,n,iy,i)-ya2)**2-xj2**2
              xmax=max(xmax,x)
              xmin=min(xmin,x)
              ymax=max(ymax,y)
              ymin=min(ymin,y)
1030        continue
          endif
1020    continue
        xsp=xmax-xmin
        ysp=ymax-ymin
        xsp=xsp*1.3d0
        ysp=ysp*1.3d0
        xa=xmax+xmin
        xmax=(xa+xsp)*.5d0
        xmin=(xa-xsp)*.5d0
        ya=ymax+ymin
        ymax=(ya+ysp)*.5d0
        ymin=(ya-ysp)*.5d0
        write(lfnplt,*)'NEWFRAME;SET FONT DUPLEX;SET SYMBOL .M SIZE 1.4'
        write(lfnplt,*)'SET WINDOW X 4 11 Y 2 9'
        write(lfnplt,*)'TITLE LEFT   ''',titlel(k+3),''''
        write(lfnplt,*)'TITLE BOTTOM ''',titleb(k+3),''''
        write(lfnplt,9001)sngl(xmin),sngl(xmax),sngl(ymin),sngl(ymax)
        call tdinit(lfnplt,'PLOT',' ')
        do 1230 i=1,np0
          if(kptbl(i,1) .le. np)then
            xa1=x0(ix,i)
            ya1=y0(ix,i)
            xj1=xj(ix,i)
            phi1=xphi(ix,i)
            xa2=x0(iy,i)
            ya2=y0(iy,i)
            xj2=xj(iy,i)
            phi2=xphi(iy,i)
            do 1220 n=1,nturn
c             x=(zbuf(1,n,ix,i)-xa1)-xj1*cos(-(n-1)*xnu(ix,i)+phi1)
c             y=(zbuf(1,n,iy,i)-xa2)-xj2*cos(-(n-1)*xnu(iy,i)+phi2)
              x=(zbuf(1,n,ix,i)-xa1)**2+(zbuf(2,n,ix,i)-ya1)**2-xj1**2
              y=(zbuf(1,n,iy,i)-xa2)**2+(zbuf(2,n,iy,i)-ya2)**2-xj2**2
              call tdput(x,y)
1220        continue
          endif
1230    continue
        call tdterm
1010  continue
      return
      end
