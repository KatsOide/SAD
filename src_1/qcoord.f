      subroutine qcoord(trans,cod,dx,dy,dz,chi1,chi2,chi3,dir,coup)
            use tfstk
      use ffs
      use tffitcode
      implicit none
      real*8 trans(4,5),cod(6),transe(6,12),dx,dy,dz,chi1,chi2,chi3,
     $     beam(42)
      logical coup,dir
      call tinitr(transe)
      call tcoorde(transe,cod,beam,dx,dy,dz,chi1,chi2,chi3,dir,1)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
