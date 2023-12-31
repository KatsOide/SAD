      subroutine tthine(trans,cod,beam,srot,nord,al,ak,dx,dy,theta,enarad)
      use ffs_flag
      use tmacro
      use drife
      use temw,only:bsir0,tsetr0,tmulbs
      use multa, only:fact
      use chg,only:tchge
      use sad_basics
      implicit none
      integer*4 kord,nord
      real*8 ,intent(in):: al,ak,dx,dy,theta
      real*8 b1,aki,ala,alb
      real*8 ,intent(inout):: trans(6,12),cod(6),beam(42),srot(3,9)
      complex*16 cx,cx1
      logical*4 ,intent(in):: enarad
      logical*4 krad
      real*8 trans1(6,13)
      if(ak .eq. 0.d0)then
        call tdrife0(trans,cod,beam,srot,al,0.d0,0.d0,.true.,.false.,irad)
        return
      endif
      call tchge(trans,cod,beam,srot,
     $     dx,dy,0.d0,theta,0.d0,0.d0,0.d0,0.d0,.true.)
      krad=enarad .and. al .ne. 0.d0
      if(krad)then
        call tsetr0(trans(:,1:6),cod(1:6),0.d0,0.d0)
      endif
      kord=nord/2-1
      b1=0.d0
      aki=ak/fact(kord)
      if(al .ne. 0.d0)then
        ala=al/6.d0
        alb=al/1.5d0
        call tdrife0(trans,cod,beam,srot,ala,0.d0,0.d0,.true.,.false.,irad)
        b1=brhoz/al*aki
        aki=aki*.5d0
      endif
      call tinitr(trans1)
      if(kord .eq. 0)then
        cx=(1.d0,0.d0)
        cx1=(0.d0,0.d0)
      else
        if(kord == 1)then
          cx1=(1.d0,0.d0)
        else
          cx1=dcmplx(cod(1),-cod(3))**(kord-1)
        endif
        cx=dcmplx(cod(1),-cod(3))*cx1
        cx1=cx1*kord
        trans1(2,1)=-aki*dble(cx1)
        trans1(2,3)=-aki*imag(cx1)
        trans1(4,1)=-aki*imag(cx1)
        trans1(4,3)= aki*dble(cx1)
        call tmultr5(trans,trans1,irad)
      endif
      call tmulbs(beam ,trans1,.true.)
c      if(enarad .and. al .ne. 0.d0)then
c        bx=-b1*imag(cx)
c        by= b1*dble(cx)
c        bxx=-b1*imag(cx1)
c        bxy= b1*dble(cx1)
c        call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $       bxx,bxy,0.d0,0.d0,0.d0,
c     $       al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
c      endif
      cod(2)=cod(2)-aki*dble(cx)
      cod(4)=cod(4)-aki*imag(cx)
      if(al .ne. 0.d0)then
        bsir0=bsir0-aki*imag(cx)/al
        call tdrife0(trans,cod,beam,srot,alb,0.d0,al*.5d0,.true.,krad,irad)
        call tinitr(trans1)
        if(kord .eq. 0)then
          cx=(1.d0,0.d0)
          cx1=(0.d0,0.d0)
        else
          if(kord == 1)then
            cx1=(1.d0,0.d0)
          else
            cx1=dcmplx(cod(1),-cod(3))**(kord-1)
          endif
          cx=dcmplx(cod(1),-cod(3))*cx1
          cx1=cx1*kord
          trans1(2,1)=-aki*dble(cx1)
          trans1(2,3)=-aki*imag(cx1)
          trans1(4,1)=-aki*imag(cx1)
          trans1(4,3)= aki*dble(cx1)
          call tmultr5(trans,trans1,irad)
        endif
        call tmulbs(beam ,trans1,.true.)
c        if(enarad)then
c          bx=-b1*imag(cx)
c          by= b1*dble(cx)
c          bxx=-b1*imag(cx1)
c          bxy= b1*dble(cx1)
c          call trade(trans,beam,cod,bx,by,0.d0,0.d0,
c     $         bxx,bxy,0.d0,0.d0,0.d0,
c     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
c        endif
        cod(2)=cod(2)-aki*dble(cx)
        cod(4)=cod(4)-aki*imag(cx)
        bsir0=bsir0+aki*imag(cx)/al
        call tdrife0(trans,cod,beam,srot,ala,0.d0,al*.5d0,.true.,krad,irad)
      endif
      bradprev=0.d0
      call tchge(trans,cod,beam,srot,
     $     dx,dy,0.d0,theta,0.d0,0.d0,0.d0,0.d0,.false.)
      return
      end
