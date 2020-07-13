      subroutine qcav(trans,cod,k,
     1     al,vc,harm,phi,freq,dx,dy,theta,
     $     v10,v20,v11,v02,fringe,mfring,autophi,coup)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 k,mfring
      real*8 trans(4,5),cod(6),srot(3,9),al,vc,harm,phi,
     $     freq,dx,dy,theta,transe(6,12),beam(42),v10,v20,v11,v02
      logical*4 coup,fringe,autophi
      logical*4 rfsw0
      rfsw0=rfsw
      rfsw=rfsw .and. trpt
      call tinitr(transe)
      call tcave(transe,cod,beam,srot,k,al,vc,harm,phi,freq,
     $     dx,dy,theta,v10,v20,v11,v02,fringe,mfring,autophi,1)
      call qcopymatg(trans,transe,k)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      rfsw=rfsw0
      return
      end

      subroutine qtcav(trans,cod,
     1     al,vc,harm,phi,freq,dx,dy,theta,coup)
      use tfstk
      use tffitcode
      use ffs
      implicit none
      real*8 trans(4,5),cod(6),srot(3,9),al,vc,harm,phi,
     $     freq,dx,dy,transe(6,12),beam(42),theta
      logical*4 coup
      call tinitr(transe)
      call ttcave(transe,cod,beam,srot,al,vc,harm,phi,freq,
     $     dx,dy,theta,1,.false.)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
