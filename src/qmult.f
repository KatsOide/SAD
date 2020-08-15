      subroutine qmult(trans,cod,k,al,ak,bz,
     $     phi0,psi1,psi2,apsi1,apsi2,
     1     dx,dy,dz,chi1,chi2,theta,dtheta,
     $     eps0,fringe,f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,bfrm,
     $     vc,harm,phi,freq,wakew1,autophi,ini,coup)
      use tfstk
      use ffs
      use tffitcode
      use kyparam, only:nmult
      implicit none
      integer*4 ,intent(in):: mfring,k
      real*8 ,intent(inout):: trans(4,5),cod(6)
      real*8 ,intent(in)::dx,dy,theta,f1in,f2in,f1out,f2out,
     $     al,bz,eps0,chi1,chi2,dz,vc,harm,phi,freq,wakew1,
     $     psi1,psi2,phi0,dtheta,apsi1,apsi2,fb1,fb2
      real*8 transe(6,12),beam(42),srot(3,9)
      complex*16 ,intent(in):: ak(0:nmult)
      logical*4 ,intent(out):: coup
      logical*4 ,intent(in):: fringe,bfrm,autophi,ini
      logical*4 rfsw0
      rfsw0=rfsw
      rfsw=rfsw .and. trpt
      if(ini)then
        call tinitr(transe)
      endif
      call tmulte(transe,cod,beam,srot,k,al,ak,bz,
     $     phi0,psi1,psi2,apsi1,apsi2,
     1     dx,dy,dz,chi1,chi2,theta,dtheta,
     $     eps0,.false.,fringe,
     $     f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,bfrm,vc,harm,phi,freq,wakew1,
     $     1.d0,autophi,0)
      call qcopymatg(trans,transe,k)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      rfsw=rfsw0
      return
      end

      subroutine qcopymat(trans,transe,dir)
      implicit none
      real*8 trans(4,5),transe(6,6)
      logical*4 dir
      if(dir)then
        transe(1,1)=trans(1,1)
        transe(1,2)=trans(1,2)
        transe(1,3)=trans(1,3)
        transe(1,4)=trans(1,4)
        transe(1,5)=0.d0
        transe(1,6)=trans(1,5)
        transe(2,1)=trans(2,1)
        transe(2,2)=trans(2,2)
        transe(2,3)=trans(2,3)
        transe(2,4)=trans(2,4)
        transe(2,5)=0.d0
        transe(2,6)=trans(2,5)
        transe(3,1)=trans(3,1)
        transe(3,2)=trans(3,2)
        transe(3,3)=trans(3,3)
        transe(3,4)=trans(3,4)
        transe(3,5)=0.d0
        transe(3,6)=trans(3,5)
        transe(4,1)=trans(4,1)
        transe(4,2)=trans(4,2)
        transe(4,3)=trans(4,3)
        transe(4,4)=trans(4,4)
        transe(4,5)=0.d0
        transe(4,6)=trans(4,5)
        transe(5,1)=0.d0
        transe(5,2)=0.d0
        transe(5,3)=0.d0
        transe(5,4)=0.d0
        transe(5,5)=1.d0
        transe(5,6)=0.d0
        transe(6,1)=0.d0
        transe(6,2)=0.d0
        transe(6,3)=0.d0
        transe(6,4)=0.d0
        transe(6,5)=0.d0
        transe(6,6)=1.d0
      else
        trans(1,1)=transe(1,1)
        trans(1,2)=transe(1,2)
        trans(1,3)=transe(1,3)
        trans(1,4)=transe(1,4)
        trans(1,5)=transe(1,6)
        trans(2,1)=transe(2,1)
        trans(2,2)=transe(2,2)
        trans(2,3)=transe(2,3)
        trans(2,4)=transe(2,4)
        trans(2,5)=transe(2,6)
        trans(3,1)=transe(3,1)
        trans(3,2)=transe(3,2)
        trans(3,3)=transe(3,3)
        trans(3,4)=transe(3,4)
        trans(3,5)=transe(3,6)
        trans(4,1)=transe(4,1)
        trans(4,2)=transe(4,2)
        trans(4,3)=transe(4,3)
        trans(4,4)=transe(4,4)
        trans(4,5)=transe(4,6)
      endif
      return
      end

      subroutine qcopymatg(trans,transe,k)
      use ffs_pointer, only:gammab
      implicit none
      real*8 trans(4,5),transe(6,6),rg
      integer*4 k
      if(gammab(k) .eq. gammab(k+1))then
        trans(1,1)=transe(1,1)
        trans(1,2)=transe(1,2)
        trans(1,3)=transe(1,3)
        trans(1,4)=transe(1,4)
        trans(1,5)=transe(1,6)
        trans(2,1)=transe(2,1)
        trans(2,2)=transe(2,2)
        trans(2,3)=transe(2,3)
        trans(2,4)=transe(2,4)
        trans(2,5)=transe(2,6)
        trans(3,1)=transe(3,1)
        trans(3,2)=transe(3,2)
        trans(3,3)=transe(3,3)
        trans(3,4)=transe(3,4)
        trans(3,5)=transe(3,6)
        trans(4,1)=transe(4,1)
        trans(4,2)=transe(4,2)
        trans(4,3)=transe(4,3)
        trans(4,4)=transe(4,4)
        trans(4,5)=transe(4,6)
      else
        rg=sqrt(gammab(k+1)/gammab(k))
        trans(1,1)=transe(1,1)*rg
        trans(1,2)=transe(1,2)*rg
        trans(1,3)=transe(1,3)*rg
        trans(1,4)=transe(1,4)*rg
        trans(1,5)=transe(1,6)*rg
        trans(2,1)=transe(2,1)*rg
        trans(2,2)=transe(2,2)*rg
        trans(2,3)=transe(2,3)*rg
        trans(2,4)=transe(2,4)*rg
        trans(2,5)=transe(2,6)*rg
        trans(3,1)=transe(3,1)*rg
        trans(3,2)=transe(3,2)*rg
        trans(3,3)=transe(3,3)*rg
        trans(3,4)=transe(3,4)*rg
        trans(3,5)=transe(3,6)*rg
        trans(4,1)=transe(4,1)*rg
        trans(4,2)=transe(4,2)*rg
        trans(4,3)=transe(4,3)*rg
        trans(4,4)=transe(4,4)*rg
        trans(4,5)=transe(4,6)*rg
      endif
      return
      end
