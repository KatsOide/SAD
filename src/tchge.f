      subroutine tchge(trans,cod,beam,dx,dy,theta,dtheta,phi0,ent,ld)
      use tfstk
      use ffs_flag
      use tmacro
      use bendeb, only:tbrote
      implicit none
      integer*4 ld
      real*8 trans(6,12),cod(6),beam(42),
     $     trans1(6,6),trans2(6,6),
     $     dx,dy,theta,cost,sint,xi,pxi,dtheta,phi0,th,dth
      logical ent
      if(phi0 .eq. 0.d0)then
        th=theta+dtheta
        dth=0.d0
      else
        th=theta
        dth=dtheta
      endif
      if(ent)then
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
        call tinitr(trans1)
        if(th .ne. 0.d0)then
          cost=cos(th)
          sint=sin(th)
          trans1(1,1)= cost
          trans1(1,3)=-sint
          trans1(3,1)= sint
          trans1(3,3)= cost
          trans1(2,2)= cost
          trans1(2,4)=-sint
          trans1(4,2)= sint
          trans1(4,4)= cost
          xi=cod(1)
          cod(1)= cost*xi-sint*cod(3)
          cod(3)= sint*xi+cost*cod(3)
          pxi=cod(2)
          cod(2)= cost*pxi-sint*cod(4)
          cod(4)= sint*pxi+cost*cod(4)
        endif
        if(dth .ne. 0.d0)then
          call tbrote(trans2,cod,phi0,dth)
          call tmultr5(trans1,trans2,6)
        endif
        call tmultr(trans,trans1,irad)
        call tmulbs(beam,trans1,.true.,.true.)
      else
        if(th .ne. 0.d0 .or. dth .ne. 0.d0)then
          if(dth .ne. 0.d0)then
            call tbrote(trans2,cod,phi0,dth)
          endif
          if(th .ne. 0.d0)then
            cost=cos(th)
            sint=sin(th)
            trans1(1,1)= cost
            trans1(1,3)=-sint
            trans1(3,1)= sint
            trans1(3,3)= cost
            trans1(2,2)= cost
            trans1(2,4)=-sint
            trans1(4,2)= sint
            trans1(4,4)= cost
            xi=cod(1)
            cod(1)= cost*xi-sint*cod(3)
            cod(3)= sint*xi+cost*cod(3)
            pxi=cod(2)
            cod(2)= cost*pxi-sint*cod(4)
            cod(4)= sint*pxi+cost*cod(4)
            if(dth .ne. 0.d0)then
              call tmultr5(trans2,trans1,6)
              call tmultr(trans,trans2,irad)
              call tmulbs(beam,trans2,.true.,.true.)
            else
              call tmultr(trans,trans1,irad)
              call tmulbs(beam,trans1,.true.,.true.)
            endif
          else
            call tmultr(trans,trans2,irad)
            call tmulbs(beam,trans2,.true.,.true.)
          endif
        endif
        cod(1)=cod(1)+dx
        cod(3)=cod(3)+dy
      endif
      return
      end
