c obsolete 3/18/2020 V1.1.8.1k64
c
      subroutine tserad(np,x,px,y,py,g,dv,l1,rho)
      use kyparam
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      integer*8 l1
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),g(np),dv(np),
     $     rho,alf1,alc,b,f,dpr,pr,pxi,pyi,s,pzi,u,brad,p,hh,dp,de,
     $     er,h1
      alf1=rlist(l1+ky_F1_SOL)
      if(alf1 .le. 0.d0 .or. rho .eq. 0.d0)then
        return
      endif
      alc=alf1*crad
      b=(brho/rho/alf1)**2/16.d0
      f=1.d0/rho*.25d0
      if(rfluct)then
        er=c/amass*erad
        if(trpt)then
          call tgauss_array(dv(1),np)
        else
          call tran_array(dv(1),np)
          do i=1,np
            dv(i)=(dv(i)-.5d0)*3.46410161513775461d0
          enddo
        endif
        do i=1,np
c          dpr=g(i)*(2.d0+g(i))
          dpr=g(i)
          pr=1.d0+dpr
          pxi=px(i)
          pyi=py(i)
          s=pxi**2+pyi**2
          pzi=1.d0-s/(1.d0+sqrt(1.d0-s))
          u=x(i)*pxi+y(i)*pyi
          brad=b*(4.d0*(x(i)**2+y(i)**2-u**2)+
     1            alf1*(7.d0*alf1*s+8.d0*u))
          p=pr*p0
          hh=1.d0+p**2
          dp=-hh*brad*alc/pzi
          de=er*sqrt(hh*brad)/p*dp*hh
          g(i)=dpr+dp+sqrt(abs(de))*dv(i)
        enddo
      else
        do i=1,np
c          dpr=g(i)*(2.d0+g(i))
          dpr=g(i)
          pr=1.d0+dpr
          pxi=px(i)
          pyi=py(i)
          s=pxi**2+pyi**2
          pzi=1.d0-s/(1.d0+sqrt(1.d0-s))
          u=x(i)*pxi+y(i)*pyi
          brad=b*(4.d0*(x(i)**2+y(i)**2-u**2)+
     1            alf1*(7.d0*alf1*s+8.d0*u))
          p=pr*p0
          hh=1.d0+p**2
          dp=-hh*brad*alc/pzi
          g(i)=dpr+dp
        enddo
c        write(*,*)'tserad ',dp,g(np),b,alf1
      endif
      do i=1,np
        pr=1.d0+g(i)
        h1=sqrt(1.d0+(p0*pr)**2)
        dv(i)=-g(i)*(1.d0+pr)/h1/(h1+pr*h0)+dvfs
      enddo
      return
      end
