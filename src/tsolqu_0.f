c$Header: /SAD/cvsroot/oldsad/src/tsolqu.f,v 1.15 2011/04/07 00:27:11 oide Exp $
      subroutine tsolqu(np,x,px,y,py,z,gp,dv,pz,al,ak,bz0,
     $     ak0x,ak0y,eps0)
      implicit none
      integer*4 np,i,n,ndiv
      real*8 smax
      parameter (smax=0.99d0)
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),gp(np),pz(np)
      real*8 al,ak,eps0,bz,a,b,c,d,akk,aln,eps,pr,bzp,
     $     akkp,w1,w2,ws,wd,wss,c1,s1,ch2,sh2,
     $     cbzp,dwbzp,v1wbzp,v2wbzp,vbzpsqw,dzf,
     $     bw,dw,r,ap,dpz,ak0x,ak0y,bz0,
     $     u1,u1w,u2,u2w,v1,v1w,v2,v2w,usqw,vsqw,dy,dpy,
     $     dx0,dy0,dc1,dch2,xi,yi,a12,a14,a22,a24,ra,phi,pxi
      if(ak .lt. 0.d0)then
        call tsolqurc(np,y,py,x,px,z,gp,dv,pz,al,-ak,
     $       -bz0,-ak0y,-ak0x,eps0)
        return
      endif
      bz=bz0
      if(ak .eq. 0.d0)then
        call tdrift(np,x,px,y,py,z,gp,dv,pz,al,bz,ak0x,ak0y)
        return
      endif
      if(eps0 .eq. 0.d0)then
        eps=0.1d0
      else
        eps=0.1d0*eps0
      endif
      ndiv=1+int(sqrt((ak*al)**2+(bz*al)**2)/eps)
      aln=al/ndiv
      dx0=ak0x/ak
      dy0=ak0y/ak
      akk=ak/al
      do i=1,np
c        pr=(1.d0+gp(i))**2
        pr=(1.d0+gp(i))
        bzp=bz/pr
        akkp=akk/pr
        w1=sqrt((bzp**2+sqrt(bzp**4+4.d0*akkp**2))*.5d0)
        w2=akkp/w1
        wss=w1**2+w2**2
        ws=w1+w2
        wd=bzp**2/ws
        c1=cos(aln*w1)
        s1=sin(aln*w1)
        if(c1 .ge. 0.d0)then
          dc1=-s1**2/(1.d0+c1)
        else
          dc1=c1-1.d0
        endif
        ch2=cosh(aln*w2)
        sh2=sinh(aln*w2)
        dch2=sh2**2/(1.d0+ch2)
        px(i)=px(i)+bzp*y(i)*.5d0
        py(i)=py(i)-bzp*x(i)*.5d0
        ra=aln*0.5d0
        do n=1,ndiv
          ap=min(smax,px(i)**2+py(i)**2)
          dpz=-ap/(1.d0+sqrt(1.d0-ap))
          r=-dpz/(1.d0+dpz)*ra
          ra=aln
          if(bzp .eq. 0.d0)then
            x(i)=x(i)+px(i)*r
            y(i)=y(i)+py(i)*r
          else
            phi=r*bzp
            a24=sin(phi)
            a12=a24/bzp
            a22=cos(phi)
            if(a22 .ge. 0.d0)then
              a14=a12*a24/(1.d0+a22)
            else
              a14=(1.d0-a22)/bzp
            endif
            pxi=px(i)
            x(i) =x(i)+a12*pxi+a14*py(i)
            y(i) =y(i)-a14*pxi+a12*py(i)
            px(i)=     a22*pxi+a24*py(i)
            py(i)=    -a24*pxi+a22*py(i)
          endif
          z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
          xi=x(i)+dx0
          yi=y(i)+dy0
          a=  (w2*ws*xi-bzp*py(i))/wss
          bw= (ws*px(i)-bzp*w2*yi)/wss
          b=bw*w1
          c=  (w1*wd*xi +bzp*py(i))/wss
          dw= (-wd*px(i)+w1*bzp*yi)/wss
          d=dw*w2
          usqw=(bw**2+a**2)*w1**3
          vsqw=(c-dw)*(c+dw)*w2**3
          u1w= a*dc1+bw*s1
          u1=u1w*w1
          u2w=-a*s1 +bw*dc1
          u2=u2w*w1
          v1w= c*dch2+dw*sh2
          v1=v1w*w2
          v2w= c*sh2 +dw*dch2
          v2=v2w*w2
          x(i) =x(i) +u1w+v1w
          px(i)=px(i)+u2+v2
          if(bzp .ne. 0.d0)then
            dy =( wd*u2w+ws*v2w)/bzp
            dpy=(-wd*u1 +ws*v1)/bzp
            dzf=-0.5d0*((usqw/ws+vsqw/wd)*aln
     $           +2.d0*(((u1w+a)*v2+u1w*d)*w1*wd
     $           +((u2+b)*v1w+u2*c)*w2*ws)/wss
     $           +((u1w+a)*u2+u1w*b)*w2/ws
     $           +((v1w+c)*v2+v1w*d)*w1/wd)
     $           -dv(i)*aln
          else
            cbzp=py(i)/wss
            dwbzp=w1*yi/wss
            v1wbzp=cbzp*dch2+dwbzp*sh2
            v2wbzp=cbzp*sh2 +dwbzp*dch2
            vbzpsqw=(cbzp**2-dwbzp**2)*w2**3
            dy =py(i)*sh2/w1+yi*dch2
            dpy=py(i)*dch2+yi*sh2*w1
            dzf=-0.5d0*((usqw/ws+vbzpsqw*ws)*aln
     $           +((u1w+a)*u2+u1w*b)*.5d0
     $           +((v1wbzp+cbzp)*v2wbzp+v1wbzp*dwbzp)*akkp*ws)
     $           -dv(i)*aln
          endif
          y(i)=y(i)+dy
          py(i)=py(i)+dpy
          z(i)=z(i)+dzf
        enddo
        ap=min(smax,px(i)**2+py(i)**2)
        dpz=-ap/(1.d0+sqrt(1.d0-ap))
        r=-dpz/(1.d0+dpz)*aln*.5d0
        if(bzp .eq. 0.d0)then
          x(i)=x(i)+px(i)*r
          y(i)=y(i)+py(i)*r
        else
          phi=r*bzp
          a24=sin(phi)
          a12=a24/bzp
          a22=cos(phi)
          if(a22 .ge. 0.d0)then
            a14=a12*a24/(1.d0+a22)
          else
            a14=(1.d0-a22)/bzp
          endif
          pxi=px(i)
          x(i) =x(i)+a12*pxi+a14*py(i)
          y(i) =y(i)-a14*pxi+a12*py(i)
          px(i)=     a22*pxi+a24*py(i)
          py(i)=    -a24*pxi+a22*py(i)
        endif
        z(i)=z(i)-(3.d0+dpz)*ap/2.d0/(2.d0+dpz)*r
        px(i)=px(i)-bzp*y(i)*.5d0
        py(i)=py(i)+bzp*x(i)*.5d0
      enddo
      return
      end

      subroutine tsolqurc(np,x,px,y,py,z,gp,dv,pz,al,ak,bz,
     $     ak0x,ak0y,eps0)
      implicit none
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),dv(np),gp(np),pz(np),
     $     ak0x,ak0y,eps0,al,ak,bz
      call tsolqu(np,x,px,y,py,z,gp,dv,pz,al,ak,bz,
     $     ak0x,ak0y,eps0)
      return
      end
