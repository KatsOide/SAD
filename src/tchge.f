      subroutine tchge(trans,cod,beam,srot,
     $     dx,dy,theta,dtheta,phi0,ent)
      use tfstk
      use ffs_flag
      use tmacro
      use bendeb, only:tbrote
      use mathfun, only:pxy2dpz
      implicit none
      integer*4 , parameter :: ms=9
      real*8, intent(inout):: trans(6,12),cod(6),beam(42),srot(3,ms)
      real*8 trans1(6,6),trans2(6,6),sx(ms),
     $     dx,dy,theta,cost,sint,dtheta,phi0,th,dth,
     $     pr,cph,sph,ds,dx1,al,phih,xi,pxi,pz
      logical ent,mat
      if(phi0 .eq. 0.d0)then
        th=theta+dtheta
        dth=0.d0
        ds=0.d0
        dx1=dx
      else
        th=theta
        dth=dtheta
        phih=phi0*.5d0
        cph=cos(phih)
        sph=sin(phih)
        ds=dx*sph
        dx1=dx*cph
      endif      
      mat=ds .ne. 0.d0 .or. th .ne. 0.d0 .or. dth .ne. 0.d0
c      write(*,*)'tchge ',mat,ds,dx,dx1,dy
      pr=1.d0+cod(6)
      pz=pr*(1.d0+pxy2dpz(cod(2)/pr,cod(4)/pr))
      al=ds/pz
      if(ent)then
        cod(1)=cod(1)+cod(2)*al-dx1
        cod(3)=cod(3)+cod(4)*al-dy
        cod(5)=cod(5)-al
        if(mat)then
          if(ds .ne. 0.d0)then
            call tinitr(trans2)
            trans2(1,2)=al
            trans2(1,6)=-al*cod(2)*pr/pz**2
            trans2(3,4)=al
            trans2(3,6)=-al*cod(4)*pr/pz**2
            trans2(5,2)=trans2(1,6)
            trans2(5,4)=trans2(3,6)
            trans2(5,6)=al*pr/pz**2
          endif
          if(th .ne. 0.d0)then
            call tinitr(trans1)
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
            if(irad .gt. 6)then
              sx=srot(1,:)
              srot(1,:)= cost*sx-sint*srot(2,:)
              srot(2,:)= sint*sx+cost*srot(2,:)
            endif
            if(ds .ne. 0.d0)then
              call tmultr5(trans2,trans1,6)
            else
              trans2=trans1
            endif
          endif
          if(dth .ne. 0.d0)then
            call tbrote(trans1,cod,srot,phi0,dth)
            if(th .ne. 0.d0 .or. dth .ne. 0.d0)then
              call tmultr5(trans2,trans1,6)
            else
              trans2=trans1
            endif
          endif
          call tmultr(trans,trans2,irad)
          if(irad .gt. 6)then
            call tmulbs(beam,trans2,.true.,.true.)
          endif
        endif
      else
        if(mat)then
          if(dth .ne. 0.d0)then
            call tbrote(trans2,cod,srot,phi0,dth)
          endif
          if(th .ne. 0.d0)then
            call tinitr(trans1)
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
            if(irad .gt. 6)then
              sx=srot(1,:)
              srot(1,:)= cost*sx-sint*srot(2,:)
              srot(2,:)= sint*sx+cost*srot(2,:)
            endif
            if(dth .ne. 0.d0)then
              call tmultr5(trans2,trans1,6)
            else
              trans2=trans1
            endif
          endif
          if(ds .ne. 0.d0)then
            call tinitr(trans1)
            trans1(1,2)=al
            trans1(1,6)=-al*cod(2)*pr/pz**2
            trans1(3,4)=al
            trans1(3,6)=-al*cod(4)*pr/pz**2
            trans1(5,2)=trans1(1,6)
            trans1(5,4)=trans1(3,6)
            trans1(5,6)=al*pr/pz**2
            if(th .ne. 0.d0 .or. dth .ne. 0.d0)then
              call tmultr5(trans2,trans1,6)
            else
              trans2=trans1
            endif
          endif
          call tmultr(trans,trans2,irad)
          if(irad .gt. 6)then
            call tmulbs(beam,trans2,.true.,.true.)
          endif
        endif
        cod(1)=cod(1)+cod(2)*al-dx1
        cod(3)=cod(3)+cod(4)*al-dy
        cod(5)=cod(5)-al
      endif
      return
      end
