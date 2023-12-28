      subroutine qmult(trans,cod,k,al,ak,bz,
     $     phi0,psi1,psi2,apsi1,apsi2,
     1     dx,dy,dz,chi1,chi2,theta,dtheta,alg,phig,
     $     eps0,fringe,f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,bfrm,
     $     vc,harm,phi,freq,wakew1,autophi,ini,coup)
      use tfstk
      use ffs
      use tffitcode
      use kyparam, only:nmult
      use sad_basics
      implicit none
      integer*4 ,intent(in):: mfring,k
      real*8 ,intent(inout):: trans(4,5),cod(6)
      real*8 ,intent(in)::dx,dy,theta,f1in,f2in,f1out,f2out,
     $     al,bz,eps0,chi1,chi2,dz,vc,harm,phi,freq,wakew1,
     $     psi1,psi2,phi0,dtheta,apsi1,apsi2,fb1,fb2,alg,phig
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
     1     dx,dy,dz,chi1,chi2,theta,dtheta,alg,phig,
     $     eps0,.false.,fringe,
     $     f1in,f2in,f1out,f2out,
     $     mfring,fb1,fb2,bfrm,vc,harm,phi,freq,wakew1,
     $     1.d0,autophi,0)
c      write(*,'(a,1p10g12.4)')'qmult ',chi1,chi2,cod(1:5)
      call qcopymatg(trans,transe,k)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      rfsw=rfsw0
      return
      end
