      subroutine tthine(trans,cod,beam,nord,al,ak,
     1                 dx,dy,theta,enarad,ld)
      use ffs_flag
      use tmacro
      use temw
      implicit none
      integer*4 ld,kord,nord
      real*8 al,ak,dx,dy,theta,b1,aki,ala,alb,bx,by,bxx,bxy
      real*8 trans(6,12),cod(6),beam(42)
      complex*16 cx,cx1
      real*8 fact(0:10)
      logical*4 enarad
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1          720.d0,5040.d0,40320.d0,362880.d0,3628800.d0 /
      real*8 trans1(6,13)
      if(ak .eq. 0.d0)then
        call tdrife(trans,cod,beam,al,
     $       0.d0,0.d0,0.d0,.true.,enarad,calpol,irad,ld)
        return
      endif
      call tchge(trans,cod,beam,-dx,-dy,theta,0.d0,0.d0,.true.,ld)
      kord=nord/2-1
      b1=0.d0
      aki=ak/fact(kord)
      if(al .ne. 0.d0)then
        ala=al/6.d0
        alb=al/1.5d0
        call tdrife(trans,cod,beam,ala,
     $       0.d0,0.d0,0.d0,.true.,enarad,calpol,irad,ld)
        b1=brhoz/al*aki
        aki=aki*.5d0
      endif
      call tinitr(trans1)
      if(kord .eq. 0)then
        cx=(1.d0,0.d0)
        cx1=(0.d0,0.d0)
      else
        if(kord .eq. 1)then
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
      call tmulbs(beam ,trans1,.true.,.true.)
      if(enarad .and. al .ne. 0.d0)then
        bx=-b1*imag(cx)
        by= b1*dble(cx)
        bxx=-b1*imag(cx1)
        bxy= b1*dble(cx1)
        call trade(trans,beam,cod,bx,by,0.d0,0.d0,
     $       bxx,bxy,0.d0,0.d0,0.d0,
     $       al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
      endif
      cod(2)=cod(2)-aki*dble(cx)
      cod(4)=cod(4)-aki*imag(cx)
      if(al .ne. 0.d0)then
        call tdrife(trans,cod,beam,alb,
     $       0.d0,0.d0,0.d0,.true.,enarad,calpol,irad,ld)
        call tinitr(trans1)
        if(kord .eq. 0)then
          cx=(1.d0,0.d0)
          cx1=(0.d0,0.d0)
        else
          if(kord .eq. 1)then
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
        call tmulbs(beam ,trans1,.true.,.true.)
        if(enarad)then
          bx=-b1*imag(cx)
          by= b1*dble(cx)
          bxx=-b1*imag(cx1)
          bxy= b1*dble(cx1)
          call trade(trans,beam,cod,bx,by,0.d0,0.d0,
     $         bxx,bxy,0.d0,0.d0,0.d0,
     $         al*.5d0,0.d0,0.d0,0.d0,0.d0,.false.,.false.)
        endif
        cod(2)=cod(2)-aki*dble(cx)
        cod(4)=cod(4)-aki*imag(cx)
        call tdrife(trans,cod,beam,ala,
     $       0.d0,0.d0,0.d0,.true.,enarad,calpol,irad,ld)
      endif
      bradprev=0.d0
      call tchge(trans,cod,beam,dx,dy,-theta,0.d0,0.d0,.false.,ld)
      return
      end
