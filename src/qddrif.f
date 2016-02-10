      subroutine qddrif(dtrans,dcod,utwiss)
      use tffitcode
      use ffs
      implicit none
      integer*4 i
      real*8 dtrans(4,5),utwiss(ntwissfun),cod(6),dcod(6)
      call qtentu(dtrans,cod,utwiss,.true.)
      do i=1,5
        dtrans(1,i)= dtrans(2,i)
        dtrans(2,i)= 0.d0
        dtrans(3,i)= dtrans(4,i)
        dtrans(4,i)= 0.d0
      enddo
      dcod(1)=cod(2)
      dcod(2)=0.d0
      dcod(3)=cod(4)
      dcod(4)=0.d0
      return
      end
