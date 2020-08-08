      subroutine tchge(trans,cod,beam,srot,
     $     dx,dy,theta,dtheta,phi0,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use temw,only:tmulbs
      use bendeb, only:tbrote
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 , parameter :: ms=9
      real*8, intent(inout):: trans(6,12),cod(6),beam(42),srot(3,ms)
      real*8 , target ::trans1(6,6),trans2(6,6),trans3(6,6)
      real*8 sx(ms),dx,dy,theta,cost,sint,dtheta,phi0,th,dth,
     $     pr,dcph,sph,ds,dx1,al,phih,xi,pxi,pz,dy1,dxa,st1
      real*8 ,pointer :: transa(:,:)
      logical ent,mat,mcal
      mcal=irad .ge. 6
      if(phi0 .eq. 0.d0)then
        th=theta+dtheta
        cost=cos(th)
        sint=sin(th)
        dth=0.d0
        ds=0.d0
        dx1=dx
        dy1=dy
        pr=1.d0
        pz=1.d0
        al=0.d0
      else
        th=theta
        cost=cos(th)
        sint=sin(th)
        dth=dtheta
        phih=phi0*.5d0
        sph=sin(phih)
        dcph=-2.d0*sin(phih*.5d0)**2
        if(ent)then
          st1=sint
        else
          st1=-sint
        endif
        dxa=dx*cost-dy*st1
        dx1=dx+dcph*cost*dxa
        dy1=dy-dcph*st1 *dxa
        ds=dxa*sph
        dxa=dxa*(1.d0+dcph)
c        write(*,'(a,1p6g15.7)')'tchge ',dx,dy,dx1,dy1,ds
        pr=1.d0+cod(6)
        pz=pr*(1.d0+pxy2dpz(cod(2)/pr,cod(4)/pr))
        al=ds/pz
      endif      
      mat=ds .ne. 0.d0 .or. th .ne. 0.d0 .or. dth .ne. 0.d0
      nullify(transa)
      if(ent)then
        cod(1)=cod(1)+cod(2)*al-dx1
        cod(3)=cod(3)+cod(4)*al-dy1
        cod(5)=cod(5)-al
        if(mat)then
          if(al .ne. 0.d0 .and. mcal)then
            call tinitr(trans2)
            trans2(1,2)=al
            trans2(1,6)=-al*cod(2)*pr/pz**2
            trans2(3,4)=al
            trans2(3,6)=-al*cod(4)*pr/pz**2
            trans2(5,2)=trans2(1,6)
            trans2(5,4)=trans2(3,6)
            trans2(5,6)=al*pr/pz**2
            transa=>trans2
          endif
          if(th .ne. 0.d0)then
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
            if(associated(transa))then
              call tmultr5(transa,trans1,6)
            else
              transa=>trans1
            endif
          endif
          if(mcal)then
            if(dth .ne. 0.d0)then
              call tbrote(trans3,cod,srot,phi0,dth)
              if(associated(transa))then
                call tmultr5(transa,trans3,6)
              else
                transa=>trans3
              endif
            endif
            if(associated(transa))then
              trans(:,1:irad)=matmul(transa,trans(:,1:irad))
c     call tmultr(trans,transa,irad)
              if(irad .gt. 6)then
                call tmulbs(beam,transa,.true.)
              endif
            endif
          endif
        endif
      else
        if(mat)then
          if(dth .ne. 0.d0)then
            call tbrote(trans2,cod,srot,phi0,dth)
            transa=>trans2
          endif
          if(th .ne. 0.d0)then
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
            endif
            xi=cod(1)
            cod(1)= cost*xi-sint*cod(3)
            cod(3)= sint*xi+cost*cod(3)
            pxi=cod(2)
            cod(2)= cost*pxi-sint*cod(4)
            cod(4)= sint*pxi+cost*cod(4)
            if(mcal)then
              if(calpol .and. irad .gt. 6)then
                sx=srot(1,:)
                srot(1,:)= cost*sx-sint*srot(2,:)
                srot(2,:)= sint*sx+cost*srot(2,:)
              endif
              if(associated(transa))then
                call tmultr5(transa,trans1,6)
              else
                transa=>trans1
              endif
            endif
          endif
          if(mcal)then
            if(al .ne. 0.d0)then
              call tinitr(trans3)
              trans3(1,2)=al
              trans3(1,6)=-al*cod(2)*pr/pz**2
              trans3(3,4)=al
              trans3(3,6)=-al*cod(4)*pr/pz**2
              trans3(5,2)=trans3(1,6)
              trans3(5,4)=trans3(3,6)
              trans3(5,6)=al*pr/pz**2
              if(associated(transa))then
                call tmultr5(transa,trans3,6)
              else
                transa=>trans3
              endif
            endif
            if(associated(transa))then
              trans(:,1:irad)=matmul(transa,trans(:,1:irad))
c            call tmultr(trans,transa,irad)
              if(irad .gt. 6)then
                call tmulbs(beam,transa,.true.)
              endif
            endif
          endif
        endif
        cod(1)=cod(1)+cod(2)*al-dx1
        cod(3)=cod(3)+cod(4)*al-dy1
        cod(5)=cod(5)-al
      endif
      return
      end
