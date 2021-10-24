      subroutine tchge(trans,cod,beam,srot,
     $     dx0,dy,theta,dtheta,dchi2,alg,phig,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      use bendeb, only:tbrote
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 , parameter :: ms=9
      real*8, intent(inout):: trans(6,12),cod(6),beam(42),srot(3,ms)
      real*8 ,intent(in):: dx0,dy,theta,dtheta,alg,dchi2
      real*8 , target ::trans1(6,6),trans2(6,6),trans3(6,6)
      real*8 sx(ms),cost,sint,phig,th,xi,pxi,dx
      real*8 ,pointer :: transa(:,:)
      logical*4 ,intent(in):: ent
      logical*4 mcal
      mcal=irad .ge. 6
c      tb=phig /= 0.d0 .and. (dtheta /= 0.d0 .or. dchi2 /= 0.d0)
      if(phig /= 0.d0)then
        th=theta
      else
        th=theta+dtheta
      endif
      cost=cos(th)
      sint=sin(th)
      nullify(transa)
      if(ent)then
        if(th /= 0.d0)then
          if(mcal)then
            call tinitr(trans1)
            trans1(1,1)= cost
            trans1(1,3)=-sint
            trans1(3,1)= sint
            trans1(3,3)= cost
            trans1(2,2)= cost
            trans1(2,4)=-sint
            trans1(4,2)= sint
            trans1(4,4)= cost
            transa=>trans1
          endif
          xi=cod(1)
          cod(1)= cost*xi-sint*cod(3)
          cod(3)= sint*xi+cost*cod(3)
          pxi=cod(2)
          cod(2)= cost*pxi-sint*cod(4)
          cod(4)= sint*pxi+cost*cod(4)
          if(calpol .and. irad .gt. 6)then
            sx=srot(1,:)
            srot(1,:)= cost*sx-sint*srot(2,:)
            srot(2,:)= sint*sx+cost*srot(2,:)
          endif
          dx= cost*dx0-sint*dy
          cod(3)= cod(3)-(sint*dx0+cost*dy)
          transa=>trans1
        else
          dx=dx0
          cod(3)=cod(3)-dy
        endif
        if(phig /= 0.d0)then
c          write(*,'(a,l2,1p10g12.4)')'tchge-2 ',ent,dx0,dx,alg,phig,cod(1),cod(3),cod(5)
          call tbrote(trans3,cod,srot,alg,phig,dx,dtheta,dchi2,.true.)
c          write(*,'(a,l2,1p10g12.4)')'tchge-3 ',ent,dx0,dx,alg,phig,cod(1),cod(3),cod(5)
          if(mcal)then
            if(associated(transa))then
              call tmultr5(transa,trans3,6)
            else
              transa=>trans3
            endif
          endif
        endif
      else
        if(th /= 0.d0)then
          dx= cost*dx0-sint*dy
        else
          dx=dx0
        endif
        if(phig /= 0.d0)then
c          write(*,'(a,l2,1p10g12.4)')'tchge-6 ',ent,dx0,dx,alg,phig,cod(1),cod(3),cod(5)
          call tbrote(trans2,cod,srot,alg,phig,dx,dtheta,dchi2,.false.)
c          write(*,'(a,l2,1p10g12.4)')'tchge-7 ',ent,dx0,dx,alg,phig,cod(1),cod(3),cod(5)
          if(mcal)then
            transa=>trans2
          endif
        endif
        if(th /= 0.d0)then
          cod(3)= cod(3)+(sint*dx0+cost*dy)
          if(mcal)then
            call tinitr(trans1)
            trans1(1,1)= cost
            trans1(1,3)= sint
            trans1(3,1)=-sint
            trans1(3,3)= cost
            trans1(2,2)= cost
            trans1(2,4)= sint
            trans1(4,2)=-sint
            trans1(4,4)= cost
          endif
          xi=cod(1)
          cod(1)= cost*xi+sint*cod(3)
          cod(3)=-sint*xi+cost*cod(3)
          pxi=cod(2)
          cod(2)= cost*pxi+sint*cod(4)
          cod(4)=-sint*pxi+cost*cod(4)
          if(mcal)then
            if(calpol .and. irad .gt. 6)then
              sx=srot(1,:)
              srot(1,:)= cost*sx+sint*srot(2,:)
              srot(2,:)=-sint*sx+cost*srot(2,:)
            endif
            if(associated(transa))then
              call tmultr5(transa,trans1,6)
            else
              transa=>trans1
            endif
          endif
        else
          cod(3)=cod(3)+dy
        endif
      endif
      if(associated(transa))then
        trans(:,1:irad)=matmul(transa,trans(:,1:irad))
        if(irad .gt. 6)then
          call tmulbs(beam,transa,.true.)
        endif
      endif
      return
      end
