      subroutine qquad(trans,cod,al,ak,
     1dx,dy,theta,fringe,f1,f2,mfring,eps0,kin,coup)
      implicit none
      integer*4 mfring
      real*8 trans(4,5),cod(6),transe(6,12),beam(42),
     $     dx,dy,theta,ak,eps0,al,f1,f2
      logical*4 fringe,coup,kin
      call tinitr(transe)
      call tquade(transe,cod,beam,al,ak,
     1     dx,dy,theta,.false.,fringe,f1,f2,mfring,eps0,
     $     kin,.false.,1)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
