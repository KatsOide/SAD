      subroutine qins(trans,cod,l1,idp,
     $         enter,param,trx,coup,mat,insmat)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use sad_basics
      implicit none
      integer*4 i,l2,k
      integer*4 l1,idp
      real*8 dp,x,px,y,py,dpsix,dpsiy,x1,px1,y1,py1
      real*8 trans(4,5),param(ntwissfun),trx(6,7),cod(6)
      logical*4 enter,coup,mat,insmat
      trans=0.d0
      trans(1,1)=1.d0
      trans(2,2)=1.d0
      trans(3,3)=1.d0
      trans(4,4)=1.d0
      coup=idp .ne. 0 .or. enter .and. (cod(1) .ne. 0.d0 .or.
     $     cod(2) .ne. 0.d0 .or. cod(3) .ne. 0.d0 .or.
     $     cod(4) .ne. 0.d0)
      if(coup)then
        dp=cod(6)
        x =cod(1)
        px=cod(2)
        y =cod(3)
        py=cod(4)
        cod(1)= trx(1,1)*x+trx(1,2)*px+trx(1,3)*y+trx(1,4)*py+
     $          trx(1,6)*dp+trx(1,7)
        cod(2)= trx(2,1)*x+trx(2,2)*px+trx(2,3)*y+trx(2,4)*py+
     $          trx(2,6)*dp+trx(2,7)
        cod(3)= trx(3,1)*x+trx(3,2)*px+trx(3,3)*y+trx(3,4)*py+
     $          trx(3,6)*dp+trx(3,7)
        cod(4)= trx(4,1)*x+trx(4,2)*px+trx(4,3)*y+trx(4,4)*py+
     $          trx(4,6)*dp+trx(4,7)
        do i=1,5
          x =trans(1,i)
          px=trans(2,i)
          y =trans(3,i)
          py=trans(4,i)
          trans(1,i)= trx(1,1)*x+trx(1,2)*px+trx(1,3)*y+trx(1,4)*py
          trans(2,i)= trx(2,1)*x+trx(2,2)*px+trx(2,3)*y+trx(2,4)*py
          trans(3,i)= trx(3,1)*x+trx(3,2)*px+trx(3,3)*y+trx(3,4)*py
          trans(4,i)= trx(4,1)*x+trx(4,2)*px+trx(4,3)*y+trx(4,4)*py
        enddo
        trans(1,5)=trans(1,5)+trx(1,6)
        trans(2,5)=trans(2,5)+trx(2,6)
        trans(3,5)=trans(3,5)+trx(3,6)
        trans(4,5)=trans(4,5)+trx(4,6)
        return
      endif
      l2=l1+1
      if(enter)then
        if(mat)then
          insmat=.true.
          return
        endif
        do k=1,ntwissfun
          twiss(l2,idp,k)=param(k)
        enddo
        dpsix=param(3)/pi2
        dpsiy=param(6)/pi2
      else
        do i=l1-1,1,-1
          if(idtypec(i) == 34)then
            go to 10
          endif
        enddo
        i=l1
 10     do k=1,ntwissfun
          twiss(l2,idp,k)=twiss(i,idp,k)
        enddo
        dpsix=(param(3)-(twiss(l1,idp,3)-twiss(i,idp,3)))/pi2
        dpsiy=(param(6)-(twiss(l1,idp,6)-twiss(i,idp,6)))/pi2
      endif
      dpsix=(dpsix-nint(dpsix))*pi2
      dpsiy=(dpsiy-nint(dpsiy))*pi2
      if(mat)then
      else
        twiss(l2,idp,3)=twiss(l1,idp,3)+dpsix
        twiss(l2,idp,6)=twiss(l1,idp,6)+dpsiy
        call tftmat(trans,l1,l2,idp,.false.)
      endif
      call tinitr(trx)
      do i=1,4
        trx(5,i)=-trx(1,i)*trx(2,6)+trx(2,i)*trx(1,6)-
     $            trx(3,i)*trx(4,6)+trx(4,i)*trx(3,6)
      enddo
      dp=cod(6)
      x =cod(1)
      px=cod(2)
      y =cod(3)
      py=cod(4)
      x1 = trx(1,1)*x+trx(1,2)*px+trx(1,3)*y+trx(1,4)*py+
     $     trx(1,6)*dp
      px1= trx(2,1)*x+trx(2,2)*px+trx(2,3)*y+trx(2,4)*py+
     $     trx(2,6)*dp
      y1 = trx(3,1)*x+trx(3,2)*px+trx(3,3)*y+trx(3,4)*py+
     $     trx(3,6)*dp
      py1= trx(4,1)*x+trx(4,2)*px+trx(4,3)*y+trx(4,4)*py+
     $     trx(4,6)*dp
      cod(1)=twiss(l2,idp,mfitdx)
      cod(2)=twiss(l2,idp,mfitdpx)
      cod(3)=twiss(l2,idp,mfitdy)
      cod(4)=twiss(l2,idp,mfitdpy)
      trx(1,7)=cod(1)-x1
      trx(2,7)=cod(2)-px1
      trx(3,7)=cod(3)-y1
      trx(4,7)=cod(4)-py1
      trx(5,7)=0.d0
      trx(6,7)=0.d0
      return
      end
