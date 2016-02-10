      subroutine ttcave(trans,cod,beam,al,ak,harm,phi,freq,
     $     dx,dy,theta,ld)
      include 'inc/TMACRO.inc'
      dimension trans(6,12),cod(6),beam(42)
      common /tem/ trans1(6,13)
      if(.not. rfsw)then
        call tdrife(trans,cod,beam,al,0.d0,0.d0,0.d0,
     $       .true.,.false.,calpol,irad,ld)
        return
      endif
      if(al .ne. 0.d0)then
        call tdrife(trans,cod,beam,al*.5d0,0.d0,0.d0,0.d0,
     $       .true.,.false.,calpol,irad,ld)
      endif
      call tchge(trans,cod,beam,-dx,-dy,theta,.true.,ld)
      if(harm .eq. 0.d0)then
        w=pi2*freq/c
      else
        w=omega0*harm/c
      endif
      v=ak*p0*w
      p1=p0*(1.d0+cod(6))
      h1=p1*sqrt(1.d0+1.d0/p1**2)
c      h1=sqrt(1.d0+p1**2)
      h1=p1+1.d0/(h1+p1)
      dh1=p0*(p0+p1)/(h1+h0)*cod(6)
      phic=phi*sign(1.d0,charge)
      v1=p1/h1
      if(rfsw)then
        t=-cod(5)/v1
      else
        t=0.d0
      endif
      phii=w*t+phic
      cosp=cos(phii)
      sinp=sin(phii)
      dh=v*cod(1)*cosp
      dh2=dh1+dh
      h2=h0+dh2
      p2=h2*sqrt(1.d0-1.d0/h2**2)
c      p2=sqrt((h2-1.d0)*(h2+1.d0))
      pf    =(h2+h0)/(p0+p2)*dh2/p0
      p2=p0*(1.d0+pf)
      v2=p2/h2
      call tinitr(trans1)
      a=-v*w*cod(1)*sinp
      trans1(2,5)= ak*w/v1*cosp
      trans1(2,6)= ak*w*cosp*t*p0/p1/h1**2
      trans1(5,1)=-t*p0/p2/h2**2*ak*w*cosp
      trans1(5,5)=(p2+a*t/p2/h2)/h2/v1
      trans1(5,6)=t*p0/h1**2/h2**2/p1/p2*(dh*(h2*(h2+h1)+p1**2)+a*t)
      trans1(6,1)= ak*w*cosp/v2
      trans1(6,5)=-a/v1/v2/p0
      trans1(6,6)=(p1-a*t/p1/h1)/h1/v2
      call tmultr(trans,trans1,irad)
      call tmulbs(beam,trans1,.true.,.true.)
      if(calpol)then
        call polpar(320,ld,0.d0,0.d0,0.d0,0.d0,0.d0,cod)
      endif
      cod(6)=pf
      cod(5)=-t*v2
      cod(2)=cod(2)-ak*sinp
      call tesetdv(cod(6))
      call tchge(trans,cod,beam, dx, dy,-theta,.false.,ld)
      if(al .ne. 0.d0)then
        call tdrife(trans,cod,beam,al*.5d0,0.d0,0.d0,0.d0,
     $       .true.,.false.,calpol,irad,ld)
      endif
      return
      end
