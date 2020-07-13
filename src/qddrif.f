      subroutine qddrif(dtrans,dcod,utwiss)
      use tffitcode
      use ffs
      implicit none
      real*8 ,intent(out)::  dtrans(4,5),dcod(6)
      real*8 ,intent(in):: utwiss(ntwissfun)
      real*8 cod(6)
      call qtentu(dtrans,cod,utwiss,.true.)
c      do i=1,5
        dtrans(1,:)= dtrans(2,:)
        dtrans(2,:)= 0.d0
        dtrans(3,:)= dtrans(4,:)
        dtrans(4,:)= 0.d0
c      enddo
      dcod(1)=cod(2)
      dcod(2)=0.d0
      dcod(3)=cod(4)
      dcod(4)=0.d0
      return
      end
