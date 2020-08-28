      subroutine qthin(trans,cod,nord,al,ak,
     1                 dx,dy,theta,coup)
      implicit none
      integer*4 nord
      real*8 ,intent(inout):: trans(4,5),cod(6)
      real*8 transe(6,12),beam(42),srot(3,9),
      real*8 ,intent(in):: dx,dy,theta,ak,al
      logical*4 ,intent(inout):: coup
      call tinitr(transe)
      call tthine(transe,cod,beam,srot,nord,al,ak,
     1     dx,dy,theta,.false.,1)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
