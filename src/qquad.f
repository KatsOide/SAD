      subroutine qquad(trans,cod,al,ak,
     1dx,dy,theta,fringe,f1in,f2in,f1out,f2out,mfring,eps0,
     $     kin,achro,coup)
      implicit none
      integer*4 mfring
      real*8 ,intent(inout):: trans(4,5),cod(6),beam(42),srot(3,9)
      real*8 ,intent(in):: dx,dy,theta,ak,eps0,al,f1in,f2in,f1out,f2out
      real*8 transe(6,12)
      logical*4 ,intent(out):: coup
      logical*4 ,intent(in):: fringe,kin,achro
      call tinitr(transe)
      call tquade(transe,cod,beam,srot,al,ak,0.d0,
     1     dx,dy,theta,.false.,fringe,f1in,f2in,f1out,f2out,mfring,eps0,
     $     kin,achro)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
